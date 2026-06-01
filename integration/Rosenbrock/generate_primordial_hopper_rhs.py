#!/usr/bin/env python3

import argparse
import re
from pathlib import Path


FUNCTIONS = ("rhs_specie", "rhs_eint", "jac_nuc")


def find_matching_brace(text, open_index):
    depth = 0
    for index in range(open_index, len(text)):
        char = text[index]
        if char == "{":
            depth += 1
        elif char == "}":
            depth -= 1
            if depth == 0:
                return index
    raise ValueError("unmatched brace")


def extract_function_body(text, name):
    match = re.search(rf"\b{name}\s*\(", text)
    if match is None:
        raise ValueError(f"could not find {name}")

    open_brace = text.find("{", match.end())
    if open_brace < 0:
        raise ValueError(f"could not find opening brace for {name}")

    close_brace = find_matching_brace(text, open_brace)
    return text[open_brace + 1:close_brace]


def split_top_level_statements(body):
    statements = []
    start = 0
    paren_depth = 0
    bracket_depth = 0
    brace_depth = 0

    for index, char in enumerate(body):
        if char == "(":
            paren_depth += 1
        elif char == ")":
            paren_depth -= 1
        elif char == "[":
            bracket_depth += 1
        elif char == "]":
            bracket_depth -= 1
        elif char == "{":
            brace_depth += 1
        elif char == "}":
            brace_depth -= 1
        elif char == ";" and paren_depth == 0 and bracket_depth == 0 and brace_depth == 0:
            statement = body[start:index + 1].strip()
            if statement:
                statements.append(statement)
            start = index + 1

    tail = body[start:].strip()
    if tail:
        raise ValueError(f"unsplit function tail starts with: {tail[:80]!r}")

    return statements


def temp_indices(statements):
    indices = set()
    for statement in statements:
        indices.update(int(index) for index in re.findall(r"\bx(\d+)\b", statement))
    return indices


def replace_temporaries(expression):
    return re.sub(r"\bx(\d+)\b", r"x(\1)", expression)


def transform_statement(statement):
    stripped = statement.strip()

    if stripped == "using namespace Rates;":
        return None

    if stripped == "Real T = state.T;":
        return "    Real T = state.T;"

    decl = re.match(r"Real\s+x(\d+)\s*=\s*(.*);$", stripped, re.DOTALL)
    if decl:
        return f"    x({decl.group(1)}) = {replace_temporaries(decl.group(2).strip())};"

    assign = re.match(r"x(\d+)\s*=\s*(.*);$", stripped, re.DOTALL)
    if assign:
        return f"    x({assign.group(1)}) = {replace_temporaries(assign.group(2).strip())};"

    return "    " + replace_temporaries(stripped)


def transform_statements(statements):
    output = []
    for statement in statements:
        transformed = transform_statement(statement)
        if transformed is not None:
            output.append(transformed)
    return output


def phase_ranges(statement_count, phase_count):
    chunk = (statement_count + phase_count - 1) // phase_count
    ranges = []
    for begin in range(0, statement_count, chunk):
        ranges.append((begin, min(begin + chunk, statement_count)))
    return ranges


def emit_phase_table(name, ranges):
    lines = [
        f"constexpr PhaseRange {name}_phases[] = {{"
    ]
    for begin, end in ranges:
        lines.append(f"    {{{begin}, {end}}},")
    lines.append("};")
    return "\n".join(lines)


def emit_function(name, statements, temp_count):
    transformed = transform_statements(statements)

    if name == "rhs_specie":
        signature = [
            "template <typename YdotType>",
            "AMREX_GPU_HOST_DEVICE AMREX_INLINE",
            "void rhs_specie_scratch(const burn_t& state,",
            "                         YdotType& ydot,",
            "                         const Array1D<Real, 0, NumSpec-1>& X,",
            "                         Real const /*z*/,",
            "                         Real* x_scratch,",
            "                         const int zone,",
            "                         const int stride)"
        ]
    elif name == "rhs_eint":
        signature = [
            "AMREX_GPU_HOST_DEVICE AMREX_INLINE",
            "Real rhs_eint_scratch(const burn_t& state,",
            "                      const Array1D<Real, 0, NumSpec-1>& X,",
            "                      Real const z,",
            "                      Real* x_scratch,",
            "                      const int zone,",
            "                      const int stride)"
        ]
    elif name == "jac_nuc":
        signature = [
            "template <typename MatrixType>",
            "AMREX_GPU_HOST_DEVICE AMREX_INLINE",
            "void jac_nuc_scratch(const burn_t& state,",
            "                     MatrixType& jac,",
            "                     const Array1D<Real, 0, NumSpec-1>& X,",
            "                     Real const z,",
            "                     Real* x_scratch,",
            "                     const int zone,",
            "                     const int stride)"
        ]
    else:
        raise ValueError(name)

    lines = []
    lines.extend(signature)
    lines.append("{")
    lines.append("    using namespace Rates;")
    lines.append("    ScratchView x{x_scratch, zone, stride};")
    lines.extend(transformed)
    lines.append("}")
    return "\n".join(lines)


def emit_header(source_path, functions):
    all_indices = set()
    for statements in functions.values():
        all_indices.update(temp_indices(statements))

    if not all_indices:
        raise ValueError("no generated temporaries found")

    temp_count = max(all_indices) + 1

    lines = [
        "/* Do not edit directly.",
        "   Generated by integration/Rosenbrock/generate_primordial_hopper_rhs.py",
        f"   from {source_path}. */",
        "",
        "#ifndef PRIMORDIAL_RODAS5P_HOPPER_GENERATED_H",
        "#define PRIMORDIAL_RODAS5P_HOPPER_GENERATED_H",
        "",
        "#include <cmath>",
        "",
        "#include <AMReX_Array.H>",
        "#include <AMReX_REAL.H>",
        "",
        "#include <ArrayUtilities.H>",
        "#include <actual_network.H>",
        "#include <burn_type.H>",
        "#include <extern_parameters.H>",
        "#include <primordial_rodas5p_hopper.H>",
        "",
        "namespace primordial_rodas5p_hopper_generated {",
        "",
        "using namespace amrex;",
        "using namespace ArrayUtil;",
        "using namespace network_rp;",
        "",
        "struct PhaseRange {",
        "    int begin;",
        "    int end;",
        "};",
        "",
        "struct ScratchView {",
        "    Real* data;",
        "    int zone;",
        "    int stride;",
        "",
        "    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE",
        "    Real& operator() (const int slot) const",
        "    {",
        "        return data[slot * stride + zone];",
        "    }",
        "};",
        "",
        "struct SharedVector {",
        "    Real* data;",
        "    int zone;",
        "    int stride;",
        "",
        "    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE",
        "    Real& operator() (const int n) const",
        "    {",
        "        return data[(n - 1) * stride + zone];",
        "    }",
        "};",
        "",
        "struct SharedMatrix {",
        "    Real* data;",
        "    int zone;",
        "    int stride;",
        "",
        "    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE",
        "    Real& operator() (const int i, const int j) const",
        "    {",
        "        return data[((i - 1) * neqs + (j - 1)) * stride + zone];",
        "    }",
        "};",
        "",
        f"constexpr int generated_scratch_count = {temp_count};",
        "static_assert(generated_scratch_count <= primordial_rodas5p_hopper::scratch_count,",
        "              \"generated scratch exceeds Hopper scratch layout\");",
        f"constexpr int rhs_specie_statement_count = {len(functions['rhs_specie'])};",
        f"constexpr int rhs_eint_statement_count = {len(functions['rhs_eint'])};",
        f"constexpr int jac_nuc_statement_count = {len(functions['jac_nuc'])};",
        "",
        emit_phase_table("rhs_specie", phase_ranges(len(functions["rhs_specie"]), 3)),
        "",
        emit_phase_table("rhs_eint", phase_ranges(len(functions["rhs_eint"]), 1)),
        "",
        emit_phase_table("jac_nuc", phase_ranges(len(functions["jac_nuc"]), 3)),
        "",
    ]

    for name in FUNCTIONS:
        lines.append(emit_function(name, functions[name], temp_count))
        lines.append("")

    lines.extend([
        "AMREX_GPU_HOST_DEVICE AMREX_INLINE",
        "void actual_rhs_scratch(const burn_t& state,",
        "                        Real* rhs_tmp,",
        "                        Real* x_scratch,",
        "                        const int zone,",
        "                        const int stride)",
        "{",
        "    Real z = redshift;",
        "",
        "    Array1D<Real, 0, NumSpec-1> X;",
        "    for (int i = 0; i < NumSpec; ++i) {",
        "        X(i) = state.xn[i];",
        "    }",
        "",
        "    SharedVector ydot{rhs_tmp, zone, stride};",
        "    rhs_specie_scratch(state, ydot, X, z, x_scratch, zone, stride);",
        "    ydot(net_ienuc) = rhs_eint_scratch(state, X, z, x_scratch, zone, stride);",
        "}",
        "",
        "AMREX_GPU_HOST_DEVICE AMREX_INLINE",
        "void actual_jac_scratch(const burn_t& state,",
        "                        Real* jac_tmp,",
        "                        Real* x_scratch,",
        "                        const int zone,",
        "                        const int stride)",
        "{",
        "    Real z = redshift;",
        "",
        "    Array1D<Real, 0, NumSpec-1> X;",
        "    for (int i = 0; i < NumSpec; ++i) {",
        "        X(i) = state.xn[i];",
        "    }",
        "",
        "    SharedMatrix jac{jac_tmp, zone, stride};",
        "    jac_nuc_scratch(state, jac, X, z, x_scratch, zone, stride);",
        "}",
        "",
        "} // namespace primordial_rodas5p_hopper_generated",
        "",
        "#endif",
        "",
    ])

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=Path("networks/primordial_chem/actual_rhs.H"))
    parser.add_argument("--output", type=Path, default=Path("integration/Rosenbrock/primordial_rodas5p_hopper_generated.H"))
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    text = args.input.read_text()
    functions = {}
    for name in FUNCTIONS:
        functions[name] = split_top_level_statements(extract_function_body(text, name))

    header = emit_header(args.input, functions)

    if args.check:
        existing = args.output.read_text() if args.output.exists() else None
        if existing != header:
            raise SystemExit(f"{args.output} is out of date")
        return

    args.output.write_text(header)


if __name__ == "__main__":
    main()
