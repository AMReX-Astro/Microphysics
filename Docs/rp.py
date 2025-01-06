#!/usr/bin/env python3

import os
import sys
import textwrap
import itertools

MAIN_HEADER = """
+---------------------------------------+---------------------------------------------------------+------------------------------+
| parameter                             | description                                             | default value                |
+=======================================+=========================================================+==============================+
"""

SEPARATOR = """
+---------------------------------------+---------------------------------------------------------+------------------------------+
"""

ENTRY = """
| {:37} | {:55} | {:28} |
"""

WRAP_LEN = 55

class Parameter:
    # container class for the parameters

    def __init__(self):
        self.var = ""
        self.default = ""
        self.description = []
        self.category = ""
        self.namespace = ""


def pretty_category(path):
    if "/" not in path:
        # global parameters for a given top-level directory
        return ""
    # remove the top-most directory
    subdir = path.partition('/')[2]
    if path.startswith("networks/"):
        return f"NETWORK_DIR={subdir}"
    if path.startswith("EOS/"):
        return f"EOS_DIR={subdir}"
    if path.startswith("conductivity/"):
        return f"CONDUCTIVITY_DIR={subdir}"
    if path.startswith("integration/"):
        return f"INTEGRATOR_DIR={subdir}"
    if path.startswith("opacity/"):
        return f"OPACITY_DIR={subdir}"
    return path


def make_rest_table(param_files):

    params_list = []

    microphysics_home = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    for pf in param_files:

        # each file is a category
        category = pretty_category(os.path.relpath(os.path.dirname(os.path.abspath(pf)), microphysics_home))

        # default namespace is empty
        namespace = ""

        # open the file
        try:
            f = open(pf)
        except OSError:
            sys.exit(f"ERROR: {pf} does not exist")

        descr = r""

        # read in the file
        line = f.readline()
        while line:

            # we assume that parameters have an optional descriptive
            # heading before them without any blank line between the
            # description and the parameter definition.  Therefore,
            # if we encounter a blank line, zero out the description.
            if line.strip() == "":
                descr = r""
                line = f.readline()
                continue

            if line.startswith("#------"):
                line = f.readline()
                continue

            if line.startswith("@"):
                # this is a command -- we only know namespace
                fields = line.split()
                if fields[0].startswith("@namespace"):
                    namespace = fields[1].strip()
                    line = f.readline()
                    continue

            # find the description
            if line.startswith("#"):
                # handle descriptions here
                descr += line[1:].rstrip().replace("@@", r"\newline")
                line = f.readline()
                continue

            else:
                current_param = Parameter()
                line_list = line.split()

                current_param.var = line_list[0]
                current_param.default = line_list[2]
                current_param.description = descr
                current_param.category = category
                current_param.namespace = namespace

                descr = r""

                # store the current parameter in the list
                params_list.append(current_param)

            line = f.readline()


    def by_namespace(p):
        return p.namespace

    def by_category(p):
        return p.category

    print("Parameters by Namespace")
    print("=======================")

    # group by namespace first (roughly corresponds to top-level directory)

    for nm, group in itertools.groupby(sorted(params_list, key=by_namespace), key=by_namespace):

        # print the namespace

        if nm:
            nm_formatted = f"namespace: ``{nm}``"
        else:
            nm_formatted = "namespace: none"
        nmlen = len(nm_formatted)
        print(nm_formatted)
        print(nmlen*"-" + "\n")

        for c, params in itertools.groupby(sorted(group, key=by_category), key=by_category):

            # print the heading

            if c:
                print(f"**{c}:**\n")

            print(MAIN_HEADER.strip())

            for p in params:
                desc = list(textwrap.wrap(p.description.strip(), WRAP_LEN))
                if not desc:
                    desc = [""]

                for n, d in enumerate(desc):
                    if n == 0:
                        print(ENTRY.format("``"+p.var+"``", d, p.default).strip())
                    else:
                        print(ENTRY.format(" ", d, " ").strip())

                print(SEPARATOR.strip())

            print("\n\n")

def main():

    # find all of the _parameter files
    top_dir = "../"

    param_files = []
    for root, _, files in os.walk(top_dir):
        for f in files:
            if f == "_parameters":
                param_files.append(os.path.normpath("/".join([root, f])))

    make_rest_table(param_files)

if __name__ == "__main__":
    main()
