#!/usr/bin/env python3
"""Translate a pynucastro-generated Microphysics network into Loom kernel headers.

Reads the network's generated headers (tfactors.H, rate_type.H,
partition_functions.H, table_rates.H, reaclib_rates.H, derived_rates.H,
approximate_rates.H, modified_rates.H, branched_rates.H,
temperature_table_rates.H, tabular_weak_rates.H, actual_rhs.H) and writes

  pyn_rates.H    namespace mpl: tf_t / evaluate_tfactors, rate_t /
                 rate_derivs_t, the partition-function and weak-rate tables
                 (namespace-scope constexpr arrays behind compile-time table
                 ids), the interpolation helpers, and every rate function and
                 fill_rates of every rate namespace, each specialized for
                 rate_t (suffix _r) and rate_derivs_t (suffix _d).
  pyn_network.H  rhs_nuc, jac_nuc (bodies from actual_rhs.H) and the
                 per-pair blocks of evaluate_screening_{r,d}.

network.H (hand-written per framework, see assets/pynucastro/) includes both.

The translation is statement for statement. What changes, and why:
  * template <typename T> rate functions (T = rate_t or rate_derivs_t) are
    instantiated by hand: `if constexpr (std::is_same_v<T, rate_derivs_t>)`
    blocks are kept or dropped, and every call is renamed to the _r / _d copy;
  * tables: inline AMREX_GPU_DEVICE Array1D / Array3D globals become
    namespace-scope constexpr C arrays, read through `template <int ID>`
    accessors (a namespace-scope constexpr object passed by reference does
    not import: probes t01-t03; the id accessor does: t07, t08);
  * uninitialized Real locals are zero-initialized (they are passed as Real&
    out-arguments; probe r25), uninitialized record locals get `{}` (d14);
  * pf_cache_t's default constructor becomes pf_cache_init() (records must be
    aggregates: probe t04); amrex/std math -> mpl::math; C:: constants ->
    network_data.H values; `const auto& screened_rates = rate_eval.screened_rates`
    (a reference to a struct field, r10) -> direct field accesses;
    MathArray2D set/add -> operator() assignments.
The generic pynucastro library code (interp_tools.H, the evaluation part of
table_rates.H, the interpolation part of partition_functions.H, tfactors.H)
is translated by hand below; its source text is digest-checked, so a
pynucastro update that changes it stops the generator instead of silently
diverging.

usage (from the Microphysics root):
  uv run --no-project gen_pynucastro_net.py --net networks/he-burn/ase-iron \\
      --emit-dir unit_test/test_react/loom/ase_iron_vode
"""
import argparse
import hashlib
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def die(msg):
    sys.exit(f'gen_pynucastro_net: {msg}')


# ---------------------------------------------------------------------------
# digests of the hand-translated generic code (pynucastro library parts)
# ---------------------------------------------------------------------------

def norm_digest(text):
    """md5 of the text with whitespace runs collapsed (formatting-insensitive)."""
    return hashlib.md5(re.sub(r'\s+', ' ', text).strip().encode()).hexdigest()


GENERIC_DIGESTS = {
    # filled from the pinned sources (pynucastro 3.1.0-48-gad2501ae9, Microphysics d73a4c76)
    'tfactors.H': 'de8d7b5e13c469541b3a3b3464aea170',
    'interp_tools.H': 'b27fb0f6c96b3b7090a85309f4d7e086',
    'table_rates.H:code': 'bf0c6edc0e2a40eeb28f7ade31cba45e',
    'partition_functions.H:code': '5e99589541065bf89bdfc080cf16270c',
}


# ---------------------------------------------------------------------------
# C++ text helpers
# ---------------------------------------------------------------------------

def match_brace(text, i):
    """index just past the brace/paren that closes the one at text[i]."""
    open_c = text[i]
    close_c = {'{': '}', '(': ')', '[': ']'}[open_c]
    depth = 0
    j = i
    while j < len(text):
        c = text[j]
        if c == '/' and text.startswith('//', j):
            j = text.index('\n', j)
            continue
        if c == '/' and text.startswith('/*', j):
            j = text.index('*/', j) + 2
            continue
        if c == open_c:
            depth += 1
        elif c == close_c:
            depth -= 1
            if depth == 0:
                return j + 1
        j += 1
    die(f'unbalanced {open_c} at offset {i}')


def resolve_if_constexpr(text, conds):
    """Resolve `if constexpr (C) {...} [else {...}]` for every C in conds
    (condition text -> bool). The kept branch stays a `{ ... }` block."""
    out = []
    pos = 0
    pat = re.compile(r'if constexpr\s*\(')
    while True:
        m = pat.search(text, pos)
        if not m:
            out.append(text[pos:])
            break
        p0 = m.end() - 1
        p1 = match_brace(text, p0)
        cond = re.sub(r'\s+', ' ', text[p0 + 1:p1 - 1]).strip()
        if cond not in conds:
            out.append(text[pos:p1])
            pos = p1
            continue
        b0 = p1
        while text[b0].isspace():
            b0 += 1
        if text[b0] != '{':
            die(f'if constexpr ({cond}) without a braced body')
        b1 = match_brace(text, b0)
        then_blk = text[b0:b1]
        k = b1
        while k < len(text) and text[k].isspace():
            k += 1
        else_blk = None
        end = b1
        if text.startswith('else', k) and not re.match(r'else\w', text[k:k + 5]):
            e0 = k + 4
            while text[e0].isspace():
                e0 += 1
            if text[e0] != '{':
                die(f'if constexpr ({cond}) ... else without braces')
            e1 = match_brace(text, e0)
            else_blk = text[e0:e1]
            end = e1
        out.append(text[pos:m.start()])
        keep = then_blk if conds[cond] else (else_blk or '')
        # resolve nested ones inside the kept block too
        out.append(resolve_if_constexpr(keep, conds) if keep else '')
        pos = end
    return ''.join(out)


def functions(text):
    """(start, end, header, name, is_T_template) of top-level function
    definitions inside a namespace body: a definition starts at an optional
    `template <...>` line and ends at its closing brace."""
    res = []
    pat = re.compile(r'(template\s*<[^>]*>\s*)?MPL_DEVICE\s+(?:constexpr\s+)?[\w:<>, ]+?\s+(\w+)\s*\(', re.S)
    pos = 0
    while True:
        m = pat.search(text, pos)
        if not m:
            break
        p0 = m.end() - 1
        p1 = match_brace(text, p0)
        b0 = p1
        while text[b0] in ' \t\n':
            b0 += 1
        if text[b0] == ';':
            pos = b0 + 1
            continue
        if text[b0] != '{':
            die(f'function {m.group(2)}: expected a body')
        b1 = match_brace(text, b0)
        tmpl = m.group(1) or ''
        res.append((m.start(), b1, tmpl, m.group(2), 'typename T' in tmpl))
        pos = b1
    return res


# ---------------------------------------------------------------------------
# token-level rewrites shared by every file
# ---------------------------------------------------------------------------

MATH = {
    'std::exp(': 'math::exp(', 'std::log(': 'math::log(', 'std::log10(': 'math::log10(',
    'std::cbrt(': 'math::cbrt(', 'std::sqrt(': 'math::sqrt(', 'std::pow(': 'math::pow(',
    'std::abs(': 'math::abs(', 'std::max(': 'math::max(', 'std::min(': 'math::min(',
    'amrex::Math::exp10(': 'math::exp10(', 'amrex::Clamp(': 'math::clamp(',
    'amrex::Math::powi<': 'math::powi<',
    'C::n_A': 'net_n_A', 'C::enuc_conv2': 'net_enuc_conv2',
    'NetworkProperties::zion(': 'net_zion(',
}


def common(text):
    t = text
    t = re.sub(r'AMREX_GPU_(?:HOST_)?DEVICE\s+AMREX_(?:FORCE_)?INLINE', 'MPL_DEVICE', t)
    t = re.sub(r'AMREX_(?:FORCE_)?INLINE\s+AMREX_GPU_(?:HOST_)?DEVICE', 'MPL_DEVICE', t)
    t = t.replace('[[maybe_unused]] ', '')
    t = re.sub(r'\b(\d[\d.]*(?:[eE][+-]?\d+)?)_rt\b', r'\1', t)
    t = re.sub(r'amrex::Array1D<amrex::Real,', 'Array1D<Real,', t)
    t = t.replace('amrex::Real', 'Real')
    for a, b in MATH.items():
        t = t.replace(a, b)
    t = t.replace('#ifdef SCREENING', '#if MPL_GEN_SCREENING')
    # counted loops with strict upper bounds (`i < n + 1`, not `i <= n`):
    # value-identical, and Loom proves those indices in range (SUBRANGE/026)
    t = re.sub(r'for \(int (\w+) = ([^;]+); \1 <= ([^;]+); (?:\+\+\1|\1\+\+)\)',
               lambda m: f'for (int {m.group(1)} = {m.group(2)}; {m.group(1)} < {m.group(3)} + 1; ++{m.group(1)})', t)
    t = re.sub(r'^\s*AMREX_ASSERT\(.*?\);\n', '', t, flags=re.M)
    t = re.sub(r'^\s*amrex::ignore_unused\(.*?\);\n', '', t, flags=re.M)
    # uninitialized Real locals: `Real a, b;` / `Real a;` -> zero-initialized
    def zero_init(m):
        names = [n.strip() for n in m.group(2).split(',')]
        return m.group(1) + 'Real ' + ', '.join(n + '{}' for n in names) + ';'
    t = re.sub(r'^(\s*)Real ((?:\w+\s*,\s*)*\w+);', zero_init, t, flags=re.M)
    return t


def check_clean(text, where):
    code = re.sub(r'//[^\n]*', '', text)
    bad = re.findall(r'\b(?:amrex|std|autodiff|gcem|C)::\w+', code)
    if bad:
        die(f'{where}: untranslated names remain: {sorted(set(bad))[:12]}')


def namespace_body(text, ns):
    m = re.search(r'namespace\s+' + ns + r'\s*\{', text)
    if not m:
        die(f'namespace {ns} not found')
    b0 = m.end() - 1
    b1 = match_brace(text, b0)
    return text[b0 + 1:b1 - 1]


# ---------------------------------------------------------------------------
# rate namespaces: specialize template <typename T> functions
# ---------------------------------------------------------------------------

CONDS_R = {'std::is_same_v<T, rate_derivs_t>': False, 'do_T_derivatives': False}
CONDS_D = {'std::is_same_v<T, rate_derivs_t>': True, 'do_T_derivatives': True}


def specialize(body, tnames, suffix):
    """one copy of the T-template functions of a namespace body for
    T = rate_t (suffix _r) or rate_derivs_t (_d)."""
    T = 'rate_t' if suffix == '_r' else 'rate_derivs_t'
    conds = CONDS_R if suffix == '_r' else CONDS_D
    out = []
    for s, e, tmpl, name, is_t in functions(body):
        if not is_t:
            continue
        f = body[s:e]
        f = f[len(tmpl):] if f.startswith(tmpl) else re.sub(r'^template\s*<[^>]*>\s*', '', f)
        f = re.sub(r'\bT&', T + '&', f)
        f = f.replace('constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;',
                      f'constexpr int do_T_derivatives = {1 if suffix == "_d" else 0};')
        f = resolve_if_constexpr(f, conds)
        code = re.sub(r'//[^\n]*', '', f)
        if 'std::is_same_v' in code or re.search(r'\bT\b', code):
            die(f'{name}: template parameter T left after specialization')
        for n in tnames:
            f = re.sub(r'\b' + n + r'\(', n + suffix + '(', f)
        if BATCH_SETS and suffix == '_r':
            f = batch_sets(f)
        out.append('    ' + f.strip() + '\n')
    return '\n'.join(out)


BATCH_SETS = False


def batch_sets(f):
    """_r rate function: every ReacLib set gets its own ln_set_rate_k /
    set_rate_k, and the `rate += set_rate;` accumulations move, in their order,
    to just before the result is stored: the sets' exp chains are independent
    (same values, same sums), so they can be scheduled together."""
    parts = f.split('rate += set_rate;')
    if len(parts) < 3:
        return f   # 0 or 1 set: nothing to batch
    out, accs = [], []
    for k, part in enumerate(parts[:-1]):
        part = re.sub(r'\bln_set_rate\b', f'ln_set_rate_{k}', part)
        part = re.sub(r'\bset_rate\b', f'set_rate_{k}', part)
        if k == 0:
            # the declarations of ln_set_rate / set_rate are in the first part
            part = re.sub(r'Real ln_set_rate_0\{0\.0\};', '\n'.join(
                f'Real ln_set_rate_{j}{{0.0}};' for j in range(len(parts) - 1)), part)
            part = re.sub(r'Real set_rate_0\{0\.0\};', '\n'.join(
                f'Real set_rate_{j}{{0.0}};' for j in range(len(parts) - 1)), part)
        out.append(part)
        accs.append(f'rate += set_rate_{k};')
    tail = parts[-1]
    if re.search(r'\b(ln_)?set_rate\b', tail):
        die('batch_sets: set_rate used after the last accumulation')
    m = re.search(r'\n(\s*)rate_eval\.screened_rates\(', tail)
    if not m:
        die('batch_sets: result store not found')
    ind = m.group(1)
    tail = tail[:m.start()] + '\n' + ''.join(f'{ind}{x}\n' for x in accs) + tail[m.start():]
    return ''.join(out) + tail


def non_template_parts(body):
    """the namespace body without its T-template functions."""
    parts, pos = [], 0
    for s, e, tmpl, name, is_t in functions(body):
        if is_t:
            parts.append(body[pos:s])
            pos = e
    parts.append(body[pos:])
    return ''.join(parts)


# ---------------------------------------------------------------------------
# tables
# ---------------------------------------------------------------------------

def floats(init):
    vals = [v.strip() for v in re.sub(r'//[^\n]*', '', init).replace('\n', ' ').split(',')]
    return [v for v in vals if v]


def fmt_array(name, vals, per=6, indent='    '):
    lines = [f'{indent}constexpr Real {name}[{len(vals)}] = {{']
    for i in range(0, len(vals), per):
        lines.append(indent + '    ' + ', '.join(vals[i:i + per]) + ',')
    lines.append(indent + '};')
    return '\n'.join(lines)


class Tables:
    """1-D arrays (pf temperatures, pf values, weak-rate rhoy/temp) by id, and
    the weak-rate 3-D data tables by table id."""

    def __init__(self):
        self.arr1 = []   # (id, qualified C-array name, lo, hi)
        self.wt = []     # (table name j_X_Y, ntemp, nrhoy, nvars, rhoy id, temp id)

    def add1(self, cname, lo, hi):
        i = len(self.arr1) + 1
        self.arr1.append((i, cname, lo, hi))
        return i

    def accessors(self):
        L = ['// 1-D tables by compile-time id (probe t07): value, lower and upper index',
             '// of the stock amrex::Array1D (its lo() / hi()). Ids are 1..num_tab1d.',
             f'constexpr int num_tab1d = {len(self.arr1)};', '']
        L.append('template <int ID> MPL_DEVICE constexpr int tab_lo()\n{')
        L += [f'    {"if" if k == 0 else "else if"} constexpr (ID == {i}) {{ return {lo}; }}'
              for k, (i, _, lo, _) in enumerate(self.arr1)]
        L.append('    else { static_assert(ID < 0, "tab_lo: unknown table id"); return 0; }\n}\n')
        L.append('template <int ID> MPL_DEVICE constexpr int tab_hi()\n{')
        L += [f'    {"if" if k == 0 else "else if"} constexpr (ID == {i}) {{ return {hi}; }}'
              for k, (i, _, _, hi) in enumerate(self.arr1)]
        L.append('    else { static_assert(ID < 0, "tab_hi: unknown table id"); return 0; }\n}\n')
        L.append('template <int ID> MPL_DEVICE Real tab(const int i)\n{')
        L += [f'    {"if" if k == 0 else "else if"} constexpr (ID == {i}) {{ return {c}[i - {lo}]; }}'
              for k, (i, c, lo, _) in enumerate(self.arr1)]
        L.append('    else { static_assert(ID < 0, "tab: unknown table id"); return 0.0; }\n}\n')
        return '\n'.join(L)

    def wt_accessors(self):
        L = ['// weak-rate tables by table id (the Rates:: index of the rate):',
             '// amrex::Array3D<Real, 1, ntemp, 1, nrhoy, 1, num_vars> (Fortran order, T',
             '// fastest) as a flat constexpr array (probe t08); rhoy / temp as 1-D ids.']
        def chain(sig, f, dflt):
            out = [f'template <int TAB> MPL_DEVICE {sig}', '{']
            out += [f'    {"if" if k == 0 else "else if"} constexpr (TAB == Rates::{w[6]}) {{ {f(w)} }}'
                    for k, w in enumerate(self.wt)]
            out += [f'    else {{ static_assert(TAB < 0, "unknown weak-rate table"); {dflt} }}', '}', '']
            return out
        L += chain('constexpr int wt_rhoy_id()', lambda w: f'return {w[4]};', 'return 0;')
        L += chain('constexpr int wt_temp_id()', lambda w: f'return {w[5]};', 'return 0;')
        L += chain('Real wt_data(const int i, const int j, const int k)',
                   lambda w: f'return rate_tables::{w[0]}_data[(i - 1) + {w[1]} * ((j - 1) + {w[2]} * (k - 1))];',
                   'return 0.0;')
        return '\n'.join(L)


def translate_weak_tables(text, tabs, rate_of_table):
    """table_rates.H's namespace rate_tables -> constexpr arrays."""
    body = namespace_body(text, 'rate_tables')
    out = ['namespace rate_tables {', '']
    pat = re.compile(r'inline AMREX_GPU_DEVICE table_t (\w+)_meta\{\.ntemp=(\d+), \.nrhoy=(\d+), \.nvars=(\d+), \.nheader=(\d+)\};')
    for m in pat.finditer(body):
        name, nt, nr, nv = m.group(1), int(m.group(2)), int(m.group(3)), int(m.group(4))
        def arr(suffix):
            mm = re.search(r'inline AMREX_GPU_DEVICE amrex::Array(\d)D<amrex::Real,([^>]*)>\s*' + name + '_' + suffix + r'\s*\{(.*?)\};',
                           body, re.S)
            if not mm:
                die(f'table {name}_{suffix} not found')
            return mm.group(1), [x.strip() for x in mm.group(2).split(',')], floats(mm.group(3))
        d1, b1, rhoy = arr('rhoy')
        d2, b2, temp = arr('temp')
        d3, b3, data = arr('data')
        if (d1, b1) != ('1', ['1', str(nr)]) or (d2, b2) != ('1', ['1', str(nt)]):
            die(f'{name}: unexpected rhoy/temp bounds {b1} {b2}')
        if d3 != '3' or b3 != ['1', str(nt), '1', str(nr), '1', 'num_vars'] or len(data) != nt * nr * nv:
            die(f'{name}: unexpected data shape {b3}, {len(data)} values')
        if len(rhoy) != nr or len(temp) != nt:
            die(f'{name}: {len(rhoy)} rhoy / {len(temp)} temp values')
        out.append(f'    // {name}: ntemp={nt}, nrhoy={nr}, nvars={nv}, nheader={m.group(5)}')
        out.append(fmt_array(name + '_rhoy', rhoy))
        out.append(fmt_array(name + '_temp', temp))
        out.append(fmt_array(name + '_data', data, per=nt))
        out.append('')
        rid = tabs.add1(f'rate_tables::{name}_rhoy', 1, nr)
        tid = tabs.add1(f'rate_tables::{name}_temp', 1, nt)
        if name not in rate_of_table:
            die(f'no weak rate uses table {name}')
        tabs.wt.append((name, nt, nr, nv, rid, tid, rate_of_table[name]))
    out.append('} // namespace rate_tables')
    if not tabs.wt:
        die('no weak-rate tables found in table_rates.H')
    if 'num_vars = jtab_gamma' not in text or not re.search(r'jtab_rate\s*=\s*1,\s*jtab_nuloss\s*=\s*2,\s*jtab_gamma\s*=\s*3', text):
        die('table_rates.H TableVars changed')
    return '\n'.join(out)


def table_rates_code(text):
    """the generic evaluation part of table_rates.H (after the tables)."""
    i = text.index('} // namespace rate_tables') if '} // namespace rate_tables' in text else None
    if i is None:
        m = re.search(r'namespace\s+rate_tables\s*\{', text)
        i = match_brace(text, m.end() - 1)
    return text[i:]


def translate_pf(text, tabs):
    """partition_functions.H tables -> constexpr arrays; returns
    (tables text, {species: (pf id, temp id, threshold literal)})."""
    body = namespace_body(text, 'part_fun')
    out = ['namespace part_fun {', '']
    npts = dict((m.group(1), int(m.group(2))) for m in re.finditer(r'constexpr int (npts_\d+) = (\d+);', body))
    temp_ids = {}
    for m in re.finditer(r'inline AMREX_GPU_DEVICE amrex::Array1D<amrex::Real, 0, (npts_\d+)-1> (temp_array_\d+) = \{(.*?)\};', body, re.S):
        vals = floats(m.group(3))
        if len(vals) != npts[m.group(1)]:
            die(f'{m.group(2)}: {len(vals)} values, {m.group(1)} = {npts[m.group(1)]}')
        out.append(fmt_array(m.group(2), vals, per=5))
        temp_ids[m.group(2)] = (tabs.add1(f'part_fun::{m.group(2)}', 0, len(vals) - 1), m.group(1))
    pf = {}
    for m in re.finditer(r'inline AMREX_GPU_DEVICE amrex::Array1D<amrex::Real, 0, (npts_\d+)-1> (\w+)_pf_array = \{(.*?)\};', body, re.S):
        vals = floats(m.group(3))
        if len(vals) != npts[m.group(1)]:
            die(f'{m.group(2)}_pf_array: {len(vals)} values')
        out.append(fmt_array(m.group(2) + '_pf_array', vals, per=5))
        thr = re.search(r'constexpr amrex::Real ' + m.group(2) + r'_pf_threshold_T9 = ([^;]+);', body)
        if not thr:
            die(f'{m.group(2)}: no pf threshold')
        out.append(f'    constexpr Real {m.group(2)}_pf_threshold_T9 = {thr.group(1).replace("_rt", "")};')
        out.append('')
        pf[m.group(2)] = tabs.add1(f'part_fun::{m.group(2)}_pf_array', 0, len(vals) - 1)
    out.append('} // namespace part_fun')
    # the switch: species -> temp array used
    sw = re.search(r'switch \(inuc\) \{(.*?)\n    default:', text, re.S)
    if not sw:
        die('get_partition_function switch not found')
    cases = []
    for m in re.finditer(r'case (\w+):\s*if \(T9 > part_fun::(\w+)_pf_threshold_T9\) \{\s*'
                         r'part_fun::interpolate_pf\(T9, pf_cache\.index_(temp_array_\d+), part_fun::(temp_array_\d+), '
                         r'part_fun::(\w+)_pf_array, logpf, dlogpf_dT9\);\s*\}\s*break;', sw.group(1)):
        sp, a, ix, ta, b = m.groups()
        if not (sp == a == b and ix == ta):
            die(f'unexpected partition-function case for {sp}')
        cases.append((sp, ta))
    if len(cases) != sw.group(1).count('case '):
        die('could not parse every partition-function case')
    if sorted(c for c, _ in cases) != sorted(pf):
        die('partition-function cases and arrays differ')
    return '\n'.join(out), cases, pf, temp_ids


# ---------------------------------------------------------------------------
# hand-translated generic code
# ---------------------------------------------------------------------------

def interp_code(temp_ids):
    L = ['// from interp_tools.H interp_net::find_index / find_index_extrap, for a',
         '// 1-D table id (x_array(i) -> tab<ID>(i), lo() / hi() -> tab_lo / tab_hi).',
         'namespace interp_net {', '',
         'template <int ID>',
         'MPL_DEVICE int find_index(const Real x0)',
         '{',
         '    int left = tab_lo<ID>();',
         '    int right = tab_hi<ID>();',
         '',
         '    int idx = -1;',
         '',
         '    if (x0 >= tab<ID>(left) && x0 < tab<ID>(right)) {',
         '',
         '        while (left < right) {',
         '            int mid = (left + right) / 2;',
         '            if (tab<ID>(mid) > x0) {',
         '                right = mid;',
         '            } else {',
         '                left = mid + 1;',
         '            }',
         '        }',
         '',
         '        idx = right - 1;',
         '    }',
         '',
         '    return idx;',
         '}',
         '',
         'template <int ID>',
         'MPL_DEVICE int find_index_extrap(const Real x0)',
         '{',
         '    int left = tab_lo<ID>();',
         '    int right = tab_hi<ID>();',
         '',
         '    int idx{};',
         '',
         '    if (x0 < tab<ID>(left)) {',
         '        idx = left;',
         '    } else if (x0 >= tab<ID>(right)) {',
         '        idx = right - 1;',
         '    } else {',
         '        idx = find_index<ID>(x0);',
         '    }',
         '',
         '    return idx;',
         '}',
         '',
         '} // namespace interp_net', '']
    return '\n'.join(L)


TABLE_CODE = r'''
// from table_rates.H (generic part): tables by id instead of by reference;
// table_meta is not read by these functions and is dropped.
constexpr int jtab_rate = 1;
constexpr int jtab_nuloss = 2;
constexpr int jtab_gamma = 3;
constexpr int num_vars = jtab_gamma;

struct weak_rate_t
{
    Real rate;
    Real drate_dT;
    Real drate_drhoye;
    Real enu;
    Real denu_dT;
    Real denu_drhoye;
    Real gamma;
};

MPL_DEVICE Real evaluate_linear_1d(const Real fhi, const Real flo, const Real xhi, const Real xlo, const Real x)
{
    Real xx = math::clamp(x, xlo, xhi);
    Real f = flo + (fhi - flo) * (xx - xlo) / (xhi - xlo);

    return f;
}

MPL_DEVICE Real evaluate_linear_2d(const Real fip1jp1, const Real fip1j, const Real fijp1, const Real fij,
                                   const Real xhi, const Real xlo, const Real yhi, const Real ylo,
                                   const Real x, const Real y)
{
    const Real u = (math::clamp(x, xlo, xhi) - xlo) / (xhi - xlo);
    const Real v = (math::clamp(y, ylo, yhi) - ylo) / (yhi - ylo);

    const Real E = fij;
    const Real C = fijp1 - fij;
    const Real B = fip1j - fij;
    const Real A = fip1jp1 - fip1j - fijp1 + fij;

    return A * u * v + B * u + C * v + E;
}

template <int TAB>
MPL_DEVICE Real evaluate_vars(const int irhoy_lo, const int jtemp_lo,
                              const Real log_rhoy, const Real log_temp, const int component)
{
    constexpr int RHOY = wt_rhoy_id<TAB>();
    constexpr int TEMP = wt_temp_id<TAB>();

    int irhoy_hi = irhoy_lo + 1;
    int jtemp_hi = jtemp_lo + 1;

    Real rhoy_lo = tab<RHOY>(irhoy_lo);
    Real rhoy_hi = tab<RHOY>(irhoy_hi);

    Real t_lo = tab<TEMP>(jtemp_lo);
    Real t_hi = tab<TEMP>(jtemp_hi);

    Real fij     = wt_data<TAB>(jtemp_lo, irhoy_lo, component);
    Real fip1j   = wt_data<TAB>(jtemp_lo, irhoy_hi, component);
    Real fijp1   = wt_data<TAB>(jtemp_hi, irhoy_lo, component);
    Real fip1jp1 = wt_data<TAB>(jtemp_hi, irhoy_hi, component);

    Real r = evaluate_linear_2d(fip1jp1, fip1j, fijp1, fij,
                                rhoy_hi, rhoy_lo, t_hi, t_lo, log_rhoy, log_temp);

    return r;
}

template <int TAB>
MPL_DEVICE Real evaluate_dr_dtemp(const int irhoy_lo, const int jtemp_lo,
                                  const Real log_rhoy, const Real log_temp, const int component)
{
    constexpr int RHOY = wt_rhoy_id<TAB>();
    constexpr int TEMP = wt_temp_id<TAB>();

    const int irhoy_hi = irhoy_lo + 1;
    const int jtemp_hi = jtemp_lo + 1;

    const Real log_temp_lo = tab<TEMP>(jtemp_lo);
    const Real log_temp_hi = tab<TEMP>(jtemp_hi);

    Real result = 0.0;
    // The interpolated rate is constant outside the temperature table.
    if (!(log_temp < log_temp_lo || log_temp > log_temp_hi)) {

        const Real log_rhoy_lo = tab<RHOY>(irhoy_lo);
        const Real log_rhoy_hi = tab<RHOY>(irhoy_hi);

        const Real fij     = wt_data<TAB>(jtemp_lo, irhoy_lo, component);
        const Real fip1j   = wt_data<TAB>(jtemp_lo, irhoy_hi, component);
        const Real fijp1   = wt_data<TAB>(jtemp_hi, irhoy_lo, component);
        const Real fip1jp1 = wt_data<TAB>(jtemp_hi, irhoy_hi, component);

        const Real dlogT = log_temp_hi - log_temp_lo;
        const Real dlogr_dlogT_i = (fijp1 - fij) / dlogT;
        const Real dlogr_dlogT_ip1 = (fip1jp1 - fip1j) / dlogT;

        result = evaluate_linear_1d(dlogr_dlogT_ip1, dlogr_dlogT_i,
                                    log_rhoy_hi, log_rhoy_lo, log_rhoy);
    }
    return result;
}

template <int TAB>
MPL_DEVICE Real evaluate_dr_drhoYe(const int irhoy_lo, const int jtemp_lo,
                                   const Real log_rhoy, const Real log_temp, const int component)
{
    constexpr int RHOY = wt_rhoy_id<TAB>();
    constexpr int TEMP = wt_temp_id<TAB>();

    const int irhoy_hi = irhoy_lo + 1;
    const int jtemp_hi = jtemp_lo + 1;

    const Real log_rhoy_lo = tab<RHOY>(irhoy_lo);
    const Real log_rhoy_hi = tab<RHOY>(irhoy_hi);

    Real result = 0.0;
    // The interpolated rate is constant outside the density table.
    if (!(log_rhoy < log_rhoy_lo || log_rhoy > log_rhoy_hi)) {

        const Real log_temp_lo = tab<TEMP>(jtemp_lo);
        const Real log_temp_hi = tab<TEMP>(jtemp_hi);

        const Real fij     = wt_data<TAB>(jtemp_lo, irhoy_lo, component);
        const Real fip1j   = wt_data<TAB>(jtemp_lo, irhoy_hi, component);
        const Real fijp1   = wt_data<TAB>(jtemp_hi, irhoy_lo, component);
        const Real fip1jp1 = wt_data<TAB>(jtemp_hi, irhoy_hi, component);

        const Real dlogrhoy = log_rhoy_hi - log_rhoy_lo;
        const Real dlogr_dlogrhoy_j = (fip1j - fij) / dlogrhoy;
        const Real dlogr_dlogrhoy_jp1 = (fip1jp1 - fijp1) / dlogrhoy;

        result = evaluate_linear_1d(dlogr_dlogrhoy_jp1, dlogr_dlogrhoy_j,
                                    log_temp_hi, log_temp_lo, log_temp);
    }
    return result;
}

template <int do_T_derivatives, int TAB>
MPL_DEVICE void get_entries(const int irhoy_lo, const int jtemp_lo,
                            const Real log_rhoy, const Real log_temp,
                            weak_rate_t& entries)
{
    entries.rate = evaluate_vars<TAB>(irhoy_lo, jtemp_lo,
                                      log_rhoy, log_temp, jtab_rate);

    entries.enu = evaluate_vars<TAB>(irhoy_lo, jtemp_lo,
                                     log_rhoy, log_temp, jtab_nuloss);

    entries.gamma = evaluate_vars<TAB>(irhoy_lo, jtemp_lo,
                                       log_rhoy, log_temp, jtab_gamma);

    if constexpr (do_T_derivatives) {
        entries.drate_dT =
            evaluate_dr_dtemp<TAB>(irhoy_lo, jtemp_lo,
                                   log_rhoy, log_temp, jtab_rate);

        entries.drate_drhoye =
            evaluate_dr_drhoYe<TAB>(irhoy_lo, jtemp_lo,
                                    log_rhoy, log_temp, jtab_rate);

        entries.denu_dT =
            evaluate_dr_dtemp<TAB>(irhoy_lo, jtemp_lo,
                                   log_rhoy, log_temp, jtab_nuloss);

        entries.denu_drhoye =
            evaluate_dr_drhoYe<TAB>(irhoy_lo, jtemp_lo,
                                    log_rhoy, log_temp, jtab_nuloss);
    }
}

template <int do_T_derivatives, int TAB>
MPL_DEVICE void tabular_evaluate(const Real log_rhoy, const Real rhoy,
                                 const Real log_temp, const Real temp,
                                 weak_rate_t& entries)
{
    constexpr int RHOY = wt_rhoy_id<TAB>();
    constexpr int TEMP = wt_temp_id<TAB>();

    const int irhoy_lo = interp_net::find_index_extrap<RHOY>(log_rhoy);
    const int jtemp_lo = interp_net::find_index_extrap<TEMP>(log_temp);
    // find_index_extrap returns lo .. hi - 1 (MPL_ASSUME operands must be
    // locals or constants, not calls)
    constexpr int rhoy_lo_i = tab_lo<RHOY>();
    constexpr int rhoy_hi_i = tab_hi<RHOY>();
    constexpr int temp_lo_i = tab_lo<TEMP>();
    constexpr int temp_hi_i = tab_hi<TEMP>();
    MPL_ASSUME(rhoy_lo_i <= irhoy_lo && irhoy_lo < rhoy_hi_i);
    MPL_ASSUME(temp_lo_i <= jtemp_lo && jtemp_lo < temp_hi_i);

    get_entries<do_T_derivatives, TAB>(irhoy_lo, jtemp_lo,
                                       log_rhoy, log_temp, entries);

    entries.rate = math::exp10(entries.rate);
    if constexpr (do_T_derivatives) {
        entries.drate_dT = entries.rate * entries.drate_dT / temp;
        entries.drate_drhoye = entries.rate * entries.drate_drhoye / rhoy;
    } else {
        entries.drate_dT = 0.0;
        entries.drate_drhoye = 0.0;
    }
    entries.enu = -math::exp10(entries.enu);
    if constexpr (do_T_derivatives) {
        entries.denu_dT = entries.enu * entries.denu_dT / temp;
        entries.denu_drhoye = entries.enu * entries.denu_drhoye / rhoy;
    } else {
        entries.denu_dT = 0.0;
        entries.denu_drhoye = 0.0;
    }
    entries.gamma = math::exp10(entries.gamma);
}
'''


def pf_code(cases, pf_ids, temp_ids, nspec_total):
    L = ['// from partition_functions.H (generic part), arrays by id.',
         'namespace part_fun {', '',
         '// interpolate_pf: temp_array / pf_array by id',
         'template <int TA, int PF>',
         'MPL_DEVICE void interpolate_pf(const Real t9, int idx, Real& logpf, Real& dlogpf_dT9)',
         '{',
         '    // if the index is -1, we are either outside of the bounds or',
         '    // we never computed it.  try recomputing here',
         '    if (idx == -1) {',
         '        idx = interp_net::find_index<TA>(t9);',
         '    }',
         '',
         '    if (idx != -1) {',
         '        // find_index returns lo .. hi - 1 (or -1)',
         '        constexpr int ta_lo = tab_lo<TA>();',
         '        constexpr int ta_hi = tab_hi<TA>();',
         '        MPL_ASSUME(ta_lo <= idx && idx < ta_hi);',
         '',
         '        dlogpf_dT9 = (tab<PF>(idx+1) - tab<PF>(idx)) /',
         '                       (tab<TA>(idx+1) - tab<TA>(idx));',
         '',
         '        logpf = tab<PF>(idx) + dlogpf_dT9 * (t9 - tab<TA>(idx));',
         '',
         '    } else {',
         '',
         '        logpf = 0.0;',
         '        dlogpf_dT9 = 0.0;',
         '',
         '    }',
         '}', '',
         '// pf_cache_t: an aggregate (probe t04); its stock default constructor is',
         '// pf_cache_init(), called right where the stock code constructs one.',
         'struct pf_cache_t {',
         '    Array2D<Real, 1, Species::NumSpecTotal, 1, 2> data;']
    L += [f'    int index_{t};' for t in temp_ids]
    L += ['};', '',
          'MPL_DEVICE void pf_cache_init(pf_cache_t& c)',
          '{',
          '    for (int i = 1; i <= Species::NumSpecTotal; ++i) {',
          '        c.data(i, 1) = -900.0;',
          '        c.data(i, 2) = -900.0;',
          '    }']
    L += [f'    c.index_{t} = -1;' for t in temp_ids]
    L += ['}', '', '} // namespace part_fun', '',
          '// get_partition_function (cached temperature index)',
          'MPL_DEVICE void get_partition_function(const int inuc, const Real T9,',
          '                                       const part_fun::pf_cache_t& pf_cache,',
          '                                       Real& logpf, Real& dlogpf_dT9)',
          '{',
          '    logpf = 0.0;',
          '    dlogpf_dT9 = 0.0;',
          '',
          '    switch (inuc) {', '']
    for sp, ta in cases:
        L += [f'    case Species::{sp}:',
              f'        if (T9 > part_fun::{sp}_pf_threshold_T9) {{',
              f'            part_fun::interpolate_pf<{temp_ids[ta][0]}, {pf_ids[sp]}>(T9, pf_cache.index_{ta}, logpf, dlogpf_dT9);',
              '        }',
              '        break;', '']
    L += ['    default:', '',
          '        logpf = 0.0;',
          '        dlogpf_dT9 = 0.0;', '',
          '    }', '}', '',
          '// The derived rates call get_partition_function_cached(SPECIES, T9, pf_cache, ...)',
          '// with a compile-time species and tfactors.T9, after fill_rates has set',
          '// pf_cache.index_* = find_index(tfactors.T9, ...). Kernel form (value-identical):',
          '//  * the species is a template parameter: `if constexpr` picks the case the',
          '//    stock switch takes;',
          '//  * no memoization: a cache hit returns what get_partition_function',
          '//    computed for the same species and T9, i.e. the same bits, so every call',
          '//    recomputes instead (no cache check, no cache writes);',
          '//  * interpolate_pf_indexed: the index is the one fill_rates found for this',
          '//    T9, so the stock `if (idx == -1) idx = find_index(t9, ...)` would',
          '//    recompute the same value (-1) and is dropped.',
          '// The stock shape (a runtime switch, the cache check and a binary search',
          '// inlined at all 242 call sites of ase-iron, times every inlined RHS) overflowed',
          '// Loom\'s 16-bit block count (RESOURCE_EXHAUSTED / TARGET/034 block_count_u16).',
          '// The runtime-inuc versions above stay for the validators.',
          'namespace part_fun {',
          'template <int TA, int PF>',
          'MPL_DEVICE void interpolate_pf_indexed(const Real t9, const int idx, Real& logpf, Real& dlogpf_dT9)',
          '{',
          '    if (idx != -1) {',
          '        constexpr int ta_lo = tab_lo<TA>();',
          '        constexpr int ta_hi = tab_hi<TA>();',
          '        MPL_ASSUME(ta_lo <= idx && idx < ta_hi);',
          '',
          '        dlogpf_dT9 = (tab<PF>(idx+1) - tab<PF>(idx)) /',
          '                       (tab<TA>(idx+1) - tab<TA>(idx));',
          '',
          '        logpf = tab<PF>(idx) + dlogpf_dT9 * (t9 - tab<TA>(idx));',
          '    } else {',
          '        logpf = 0.0;',
          '        dlogpf_dT9 = 0.0;',
          '    }',
          '}',
          '} // namespace part_fun', '',
          'template <int inuc>',
          'MPL_DEVICE void get_partition_function_indexed(const Real T9,',
          '                                               const part_fun::pf_cache_t& pf_cache,',
          '                                               Real& logpf, Real& dlogpf_dT9)',
          '{',
          '    logpf = 0.0;',
          '    dlogpf_dT9 = 0.0;', '']
    for k, (sp, ta) in enumerate(cases):
        L += [f'    {"if" if k == 0 else "else if"} constexpr (inuc == Species::{sp}) {{',
              f'        if (T9 > part_fun::{sp}_pf_threshold_T9) {{',
              f'            part_fun::interpolate_pf_indexed<{temp_ids[ta][0]}, {pf_ids[sp]}>(T9, pf_cache.index_{ta}, logpf, dlogpf_dT9);',
              '        }',
              '    }']
    L += ['    else {',
          '        logpf = 0.0;',
          '        dlogpf_dT9 = 0.0;',
          '    }',
          '}', '']
    L += [
          'MPL_DEVICE void get_partition_function_cached(const int inuc, const Real T9,',
          '                                              part_fun::pf_cache_t& pf_cache,',
          '                                              Real& logpf, Real& dlogpf_dT9)',
          '{',
          '    if (pf_cache.data(inuc, 1) != -900.0) {',
          '        // present in cache',
          '        logpf = pf_cache.data(inuc, 1);',
          '        dlogpf_dT9 = pf_cache.data(inuc, 2);',
          '    } else {',
          '        get_partition_function(inuc, T9, pf_cache, logpf, dlogpf_dT9);',
          '        pf_cache.data(inuc, 1) = logpf;',
          '        pf_cache.data(inuc, 2) = dlogpf_dT9;',
          '    }',
          '}', '']
    return '\n'.join(L)


TFACTORS_CODE = r'''
// from tfactors.H
struct tf_t {
    Real T9;
    Real T9i;
    Real T943i;
    Real T923i;
    Real T913i;
    Real T913;
    Real T923;
    Real T953;
    Real lnT9;
};

MPL_DEVICE tf_t evaluate_tfactors(const Real T)
{
    tf_t tf{};
    tf.T9 = T * 1.e-9;
    tf.T9i = 1.0 / tf.T9;
    tf.T913 = math::cbrt(tf.T9);
    tf.T913i = 1.0 / tf.T913;
    tf.T923i = tf.T913i * tf.T913i;
    tf.T943i = tf.T9i * tf.T913i;
    tf.T923 = tf.T913 * tf.T913;
    tf.T953 = tf.T9 * tf.T923;
    tf.lnT9 = math::log(tf.T9);
    return tf;
}
'''


# ---------------------------------------------------------------------------
# rate_type.H
# ---------------------------------------------------------------------------

def rate_type(text):
    t = common(text)
    t = re.sub(r'#ifndef NET_RATE_TYPE_H.*?#include <actual_network.H>\s*', '', t, flags=re.S)
    t = re.sub(r'#endif\s*$', '', t.strip())
    t = t.replace('Rates::NumScreenPairs', 'ScreenPairs::NumScreenPairs')
    # every field zero-initialized (all are written before they are read in the
    # stock burn path; lets `rate_t r{}` replace the uninitialized local)
    t = re.sub(r'^(\s+(?:Array1D<[^;]*>|Real)\s+\w+)\s*;', r'\1{};', t, flags=re.M)
    check_clean(t, 'rate_type.H')
    return '// from rate_type.H\n' + t.strip() + '\n'


# ---------------------------------------------------------------------------
# rate namespaces
# ---------------------------------------------------------------------------

PF_PRECOMPUTE = False

RATE_FILES = [
    # (file, namespace)
    ('reaclib_rates.H', 'reaclib_rates'),
    ('temperature_table_rates.H', 'temp_tabular'),
    ('modified_rates.H', 'modified_rates'),
    ('branched_rates.H', 'branched_rates'),
    ('derived_rates.H', 'derived_rates'),
    ('approximate_rates.H', 'approximate_rates'),
    ('tabular_weak_rates.H', 'tabular_weak_rates'),
]


def rate_namespace(net, fname, ns, tnames):
    raw = open(os.path.join(net, fname)).read()
    body = namespace_body(raw, ns)
    body = common(body)
    body = re.sub(r'(\n\s*)using namespace Rates;', r'\1using namespace Rates;\1using namespace ScreenPairs;', body)
    body = body.replace('using namespace rate_tables;', '')
    # pf_cache construction
    body = re.sub(r'part_fun::pf_cache_t pf_cache\{\};', 'part_fun::pf_cache_t pf_cache{};\n        part_fun::pf_cache_init(pf_cache);', body)
    body = re.sub(r'pf_cache\.index_(temp_array_\d+) = interp_net::find_index\(tfactors\.T9, part_fun::(temp_array_\d+)\);',
                  lambda m: f'pf_cache.index_{m.group(1)} = interp_net::find_index<TA_{m.group(2)}>(tfactors.T9);', body)
    # partition functions: compile-time species (see pf_code)
    body = re.sub(r'\bget_partition_function_cached\((\w+),\s*',
                  lambda m: f'get_partition_function_indexed<Species::{m.group(1)}>(', body)
    if PF_PRECOMPUTE and ns == 'derived_rates':
        # every species' partition function once per fill_rates (into pf_cache),
        # read at the call sites: the stock cache memoizes the same function of
        # (species, T9), so the values are identical; no per-call check, no
        # recomputation (242 lookups of 30 species in ase-iron)
        sps = sorted(set(re.findall(r'get_partition_function_indexed<Species::(\w+)>', body)))
        pre = ''.join(f'''        {{
            Real pf_l{{}}, pf_d{{}};
            get_partition_function_indexed<Species::{sp}>(tfactors.T9, pf_cache, pf_l, pf_d);
            pf_cache.data(Species::{sp}, 1) = pf_l;
            pf_cache.data(Species::{sp}, 2) = pf_d;
        }}
''' for sp in sps)
        # call sites first (the precompute block below calls the same function)
        body = re.sub(r'get_partition_function_indexed<Species::(\w+)>\(tfactors\.T9, pf_cache, (\w+), (\w+)\);',
                      lambda m: f'{m.group(2)} = pf_cache.data(Species::{m.group(1)}, 1);\n        {m.group(3)} = pf_cache.data(Species::{m.group(1)}, 2);',
                      body)
        body = re.sub(r'(pf_cache\.index_temp_array_1 = interp_net::find_index<TA_temp_array_1>\(tfactors\.T9\);\n)',
                      lambda m: m.group(1) + '\n        // partition functions of every species the derived rates use, once\n' + pre, body)
    # tabular weak rates: tables by id
    body = re.sub(r'tabular_evaluate<do_T_derivatives>\((j_\w+)_meta, \1_rhoy, \1_temp, \1_data,\s*',
                  lambda m: f'tabular_evaluate<do_T_derivatives, TAB_{m.group(1)}>(', body)
    nt = non_template_parts(body)
    if re.search(r'MPL_DEVICE', nt):
        die(f'{fname}: non-template functions in namespace {ns} (not expected)')
    text = (f'namespace {ns} {{\n' + nt.strip('\n') + '\n\n' +
            specialize(body, tnames, '_r') + '\n' + specialize(body, tnames, '_d') +
            f'\n}} // namespace {ns}\n')
    return text


def template_names(net):
    names = []
    for fname, ns in RATE_FILES:
        body = common(namespace_body(open(os.path.join(net, fname)).read(), ns))
        for s, e, tmpl, name, is_t in functions(body):
            if is_t:
                names.append(name)
    # fill_rates exists in every namespace; calls are qualified (ns::fill_rates)
    return sorted(set(names), key=lambda n: (-len(n), n))


# ---------------------------------------------------------------------------
# actual_rhs.H: rhs_nuc, jac_nuc, evaluate_screening pairs
# ---------------------------------------------------------------------------

def rhs_jac(net):
    raw = open(os.path.join(net, 'actual_rhs.H')).read()
    out = []
    for name in ('rhs_nuc', 'jac_nuc'):
        m = re.search(r'(template<class MatrixType>\s*)?AMREX_GPU_DEVICE AMREX_INLINE\s*void ' + name + r'\s*\(', raw)
        if not m:
            die(f'{name} not found in actual_rhs.H')
        p1 = match_brace(raw, m.end() - 1)
        b0 = raw.index('{', p1)
        b1 = match_brace(raw, b0)
        body = common(raw[b0:b1])
        if name == 'jac_nuc':
            if 'const auto& screened_rates = rate_eval.screened_rates;' not in body:
                die('jac_nuc: screened_rates alias not found')
            body = body.replace('    const auto& screened_rates = rate_eval.screened_rates;\n', '')
            body = re.sub(r'(?<![\w.])screened_rates\(', 'rate_eval.screened_rates(', body)
            body = re.sub(r'jac\.set\(([^;]*?),\s*([^,;]*?),\s*([^;]*)\);',
                          lambda mm: f'jac({mm.group(1)}, {mm.group(2)}) = {mm.group(3)};', body)
            body = re.sub(r'jac\.add\(([^,;]*?),\s*([^,;]*?),\s*([^;]*)\);',
                          lambda mm: f'jac({mm.group(1)}, {mm.group(2)}) += {mm.group(3)};', body)
            if re.search(r'jac\.(set|add|get)\(', body):
                die('jac_nuc: MathArray2D call left')
            sig = ('MPL_DEVICE void jac_nuc(const burn_t& state,\n'
                   '                        JacNetArray2D& jac,\n'
                   '                        const Array1D<Real, 1, NumSpec>& Y,\n'
                   '                        const rate_derivs_t& rate_eval)\n')
        else:
            sig = ('MPL_DEVICE void rhs_nuc(const burn_t& state,\n'
                   '                        Array1D<Real, 1, neqs>& ydot_nuc,\n'
                   '                        const Array1D<Real, 1, NumSpec>& Y,\n'
                   '                        const Array1D<Real, 1, Rates::NumRates>& screened_rates)\n')
        body = body.replace('using namespace Rates;', 'using namespace Rates;\n    using namespace Species;')
        check_clean(body, name)
        out.append(f'// from actual_rhs.H {name}\n' + sig + body + '\n')
    # screening pairs, in evaluate_screening order
    ev = re.search(r'void evaluate_screening\(.*?\n\}\n', raw, re.S)
    pairs = re.findall(r'calculate_screen_factor\([^)]*\);\s*static_assert\(scn_fac\.z1 == [^)]*\);\s*'
                       r'actual_log_screen\(pstate, scn_fac, log_scor, dlog_scor_dT\);\s*'
                       r'rate_eval\.log_screen\((\w+)\) = log_scor;\s*'
                       r'if constexpr \(std::is_same_v<T, rate_derivs_t>\) \{\s*'
                       r'rate_eval\.dlog_screen_dT\(\1\) = dlog_scor_dT;\s*\}', ev.group(0))
    if len(pairs) != ev.group(0).count('calculate_screen_factor'):
        die('evaluate_screening: unexpected block shape')
    return '\n'.join(out), pairs


def screen_factor_values(emit_dir):
    """{pair name: {field: value}} for z1_53, z2_53, zs53 from the generated
    network_data.H (ScreenPairs enumerators and net_screen_<field>(pair))."""
    t = open(os.path.join(emit_dir, 'network_data.H')).read()
    m = re.search(r'namespace ScreenPairs \{(.*?)\}', t, re.S)
    idx = dict((n, int(v)) for n, v in re.findall(r'constexpr int (k_\w+) = (\d+);', m.group(1)))
    out = {n: {} for n in idx}
    for f in ('z1_53', 'z2_53', 'zs53'):
        mm = re.search(r'net_screen_' + f + r'\(int pair\)\n\{(.*?)\n\}', t, re.S)
        if not mm:
            return None
        byp = dict((int(a), float(b)) for a, b in re.findall(r'if \(pair == (\d+)\) \{ return ([-0-9.e+]+); \}', mm.group(1)))
        for n, k in idx.items():
            out[n][f] = byp[k]
    return out


def evaluate_screening(pairs, suffix, fvals=None):
    """actual_rhs.H evaluate_screening<do_T_derivatives, T> for one number type.
    With fvals (chabrier1998): chabrier1998_helmholtz_F is computed once per
    distinct Gamma_e * z^(5/3) argument of the network's pairs (z1_53, z2_53,
    zs53; the same doubles the stock per-pair code multiplies) and passed to
    actual_log_screen_pre_*: same arguments, same values, fewer evaluations."""
    real = suffix == '_r'
    T, num = ('rate_t', 'Real') if real else ('rate_derivs_t', 'Dual')
    L = [f'// from actual_rhs.H evaluate_screening (do_T_derivatives = {0 if real else 1}, '
         f'number_t = {"amrex::Real" if real else "autodiff::dual"})',
         f'MPL_DEVICE void evaluate_screening{suffix}(const burn_t& state,',
         '                                    const Array1D<Real, 1, NumSpec>& Y,',
         '                                    const ScreenParams& params,',
         f'                                    {T}& rate_eval)',
         '{',
         '    using namespace Rates;',
         '    using namespace ScreenPairs;', '',
         '    // Evaluate screening terms', '']
    if real:
        L += ['    Real temp = state.T;']
    else:
        L += ['    // seed the dual number for temperature before calculating anything with it',
              '    Dual temp = dual_seed(state.T);']
    L += ['', f'    plasma_state_t<{num}> pstate{{}};',
          f'    fill_plasma_state{suffix}(&pstate, temp, state.rho, Y);', '',
          '    // stock: uninitialized; Real& out-arguments (probe r25)',
          '    Real log_scor{}, dlog_scor_dT{};', '']
    fidx = {}
    if fvals is not None:
        # distinct arguments in first-use order (pair order, then z1, z2, zs)
        L += ['    // chabrier1998_helmholtz_F once per distinct Gamma_e * z^(5/3) of the',
              '    // pairs below (stock: three per pair, recomputed for repeated z)']
        L += ['    const Real cs_Gamma_e = pstate.gamma_e_fac / pstate.temp;' if real else
              '    const Dual cs_Gamma_e = dual_div(pstate.gamma_e_fac, pstate.temp);']
        for p in pairs:
            for f in ('z1_53', 'z2_53', 'zs53'):
                v = fvals[p][f]
                if v not in fidx:
                    k = len(fidx)
                    fidx[v] = k
                    arg = (f'cs_Gamma_e * cs_z{k}' if real else f'dual_mul_num(cs_Gamma_e, cs_z{k})')
                    L.append(f'    constexpr Real cs_z{k} = net_screen_{f}({p});')
                    L.append(f'    const {num} cs_F{k} = chabrier1998_helmholtz_F{suffix}({arg});')
        L.append('')
    for p in pairs:
        L += ['    {',
              f'        const screen_factors_t scn_fac = screen_factor_pair<{p}>();']
        if fvals is not None:
            a, b, c = (fidx[fvals[p][f]] for f in ('z1_53', 'z2_53', 'zs53'))
            L.append(f'        actual_log_screen_pre{suffix}(pstate, scn_fac, params, cs_F{a}, cs_F{b}, cs_F{c}, log_scor, dlog_scor_dT);')
        else:
            L.append(f'        actual_log_screen{suffix}(pstate, scn_fac, params, log_scor, dlog_scor_dT);')
        L.append(f'        rate_eval.log_screen({p}) = log_scor;')
        if not real:
            L.append(f'        rate_eval.dlog_screen_dT({p}) = dlog_scor_dT;')
        L += ['    }', '']
    L += ['}', '']
    return '\n'.join(L)


# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--net', required=True, help='network directory (networks/<...>)')
    ap.add_argument('--emit-dir', required=True)
    ap.add_argument('--batch-sets', action='store_true',
                    help='_r rate functions: exp of every ReacLib set before the accumulation')
    ap.add_argument('--pf-precompute', action='store_true',
                    help='derived rates: every partition function once per fill_rates')
    ap.add_argument('--screen-cse', action='store_true',
                    help='chabrier1998: compute helmholtz_F once per distinct argument in '
                         'evaluate_screening (needs screen.H\'s actual_log_screen_pre_*)')
    ap.add_argument('--print-digests', action='store_true',
                    help='print the digests of the hand-translated generic parts and exit')
    a = ap.parse_args()
    net = a.net

    def rd(f):
        return open(os.path.join(net, f)).read()

    for f in ['tfactors.H', 'rate_type.H', 'interp_tools.H', 'partition_functions.H', 'table_rates.H',
              'actual_rhs.H'] + [f for f, _ in RATE_FILES]:
        if not os.path.isfile(os.path.join(net, f)):
            die(f'{net}/{f} missing: not a pynucastro network?')

    pf_text = rd('partition_functions.H')
    pf_code_src = pf_text[pf_text.index('    // interpolation routine'):]
    digests = {
        'tfactors.H': norm_digest(rd('tfactors.H')),
        'interp_tools.H': norm_digest(rd('interp_tools.H')),
        'table_rates.H:code': norm_digest(table_rates_code(rd('table_rates.H'))),
        'partition_functions.H:code': norm_digest(re.sub(
            r'case \w+:.*?break;|constexpr amrex::Real \w+_pf_threshold_T9 = [^;]*;|'
            r'if constexpr \(\s*spec ==.*?\}', '', pf_code_src, flags=re.S)),
    }
    if a.print_digests:
        for k, v in digests.items():
            print(f"    '{k}': '{v}',")
        return
    for k, v in digests.items():
        if GENERIC_DIGESTS.get(k) != v:
            die(f'{k}: generic pynucastro code changed (digest {v}, translated for '
                f'{GENERIC_DIGESTS.get(k)}); update the hand translation in this script, '
                'then its digest (--print-digests)')

    global PF_PRECOMPUTE, BATCH_SETS
    PF_PRECOMPUTE = a.pf_precompute
    BATCH_SETS = a.batch_sets
    tnames = template_names(net)
    tabs = Tables()

    # weak-rate table ids: the rate that reads each table
    weak = rd('tabular_weak_rates.H')
    rate_of_table = {}
    for m in re.finditer(r'void rate_(\w+)\(T& rate_eval.*?tabular_evaluate<do_T_derivatives>\((j_\w+)_meta', weak, re.S):
        rate_of_table[m.group(2)] = 'k_' + m.group(1)
    pf_tables, pf_cases, pf_ids, temp_ids = translate_pf(pf_text, tabs)
    wt_tables = translate_weak_tables(rd('table_rates.H'), tabs, rate_of_table)

    src_md5 = {f: hashlib.md5(rd(f).encode()).hexdigest() for f in
               ['tfactors.H', 'rate_type.H', 'interp_tools.H', 'partition_functions.H', 'table_rates.H',
                'actual_rhs.H'] + [f for f, _ in RATE_FILES]}
    prov = ['// GENERATED by the loom-burn-kernel skill, scripts/gen_pynucastro_net.py, from',
            f'// {net} (pynucastro: {open(os.path.join(net, "pynucastro-info.txt")).read().strip() if os.path.isfile(os.path.join(net, "pynucastro-info.txt")) else "?"}).',
            '// Do not edit; rerun the generator. Source md5:']
    prov += [f'//   {v}  {k}' for k, v in sorted(src_md5.items())]

    L = ['// Rates of a pynucastro network for the Loom burn kernel (see README.md).'] + prov + [
        '#ifndef MPL_PYN_RATES_H', '#define MPL_PYN_RATES_H', '',
        '#include "support.H"', '#include "network_data.H"', '#include "burn.H"', '',
        '#if !MPL_GEN_PYNUCASTRO',
        '#error "pyn_rates.H needs network_data.H generated for a pynucastro network"',
        '#endif', '',
        'namespace mpl {', '',
        'using namespace Species;', '',
        TFACTORS_CODE, rate_type(rd('rate_type.H')), '',
        '// ---------------------------------------------------------------------------',
        '// tables (partition_functions.H, table_rates.H) as namespace-scope constexpr arrays',
        '// ---------------------------------------------------------------------------',
        pf_tables, '', wt_tables, '']
    L += [f'constexpr int TA_{t} = {i};' for t, (i, _) in temp_ids.items()]
    L += [f'constexpr int TAB_{w[0]} = Rates::{w[6]};' for w in tabs.wt]
    L += ['', '// X-macro lists for the validators: weak-rate tables, species with a',
          '// partition function (in get_partition_function case order)',
          '#define MPL_PYN_WEAK_TABLES(X) ' + ' '.join(f'X({w[0]})' for w in tabs.wt),
          '#define MPL_PYN_PF_SPECIES(X) ' + ' '.join(f'X({sp})' for sp, _ in pf_cases)]
    L += ['', tabs.accessors(), tabs.wt_accessors(), interp_code(temp_ids), TABLE_CODE,
          pf_code(pf_cases, pf_ids, temp_ids, None), '',
          '// ---------------------------------------------------------------------------',
          '// rate functions: each template <typename T> function as _r (T = rate_t) and',
          '// _d (T = rate_derivs_t)',
          '// ---------------------------------------------------------------------------', '']
    for fname, ns in RATE_FILES:
        L.append(f'// from {fname}')
        L.append(rate_namespace(net, fname, ns, tnames))
    L += ['} // namespace mpl', '', '#endif', '']
    rates = '\n'.join(L)
    check_clean(rates, 'pyn_rates.H')

    rj, pairs = rhs_jac(net)
    fvals = None
    if a.screen_cse:
        fvals = screen_factor_values(a.emit_dir)
        if fvals is None:
            die('--screen-cse needs network_data.H with chabrier1998 screen factors in --emit-dir')
    N = ['// RHS / Jacobian bodies and screening blocks of a pynucastro network for the',
         '// Loom burn kernel; included by network.H (see README.md).'] + prov + [
        '#ifndef MPL_PYN_NETWORK_H', '#define MPL_PYN_NETWORK_H', '',
        '// network.H includes this after pyn_rates.H, screen.H and its own',
        '// definitions; it is not self-contained.', '',
        'namespace mpl {', '',
        '#if MPL_GEN_SCREENING',
        '// X-macro list for the validators: screened pairs in evaluate_screening order',
        '#define MPL_PYN_SCREEN_PAIRS(X) ' + ' '.join(f'X({p})' for p in pairs), '',
        evaluate_screening(pairs, '_r', fvals), evaluate_screening(pairs, '_d', fvals),
        '#endif', '', rj, '} // namespace mpl', '', '#endif', '']
    netw = '\n'.join(N)
    os.makedirs(a.emit_dir, exist_ok=True)
    open(os.path.join(a.emit_dir, 'pyn_rates.H'), 'w').write(rates)
    open(os.path.join(a.emit_dir, 'pyn_network.H'), 'w').write(netw)
    print(f'gen_pynucastro_net: {len(tnames)} template rate functions (x2), '
          f'{len(pf_cases)} partition functions, {len(tabs.wt)} weak-rate tables, '
          f'{len(pairs)} screened pairs -> {a.emit_dir}/pyn_rates.H, pyn_network.H')


if __name__ == '__main__':
    main()
