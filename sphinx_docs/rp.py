#!/usr/bin/env python3

import os
import sys
import textwrap

MAIN_HEADER = """
+----------------------------------+---------------------------------------------------------+--------------------+
| parameter                        | description                                             | default value      |
+==================================+=========================================================+====================+
"""

SEPARATOR = """
+----------------------------------+---------------------------------------------------------+--------------------+
"""

ENTRY = """
| {:32} | {:55} | {:18} |
"""

WRAP_LEN = 55

class Parameter:
    # container class for the parameters

    def __init__(self):
        self.var = ""
        self.default = ""
        self.description = []
        self.category = ""

    def value(self):
        """ the value is what we sort based on """
        return self.category + "." + self.var

    def __lt__(self, other):
        return self.value() < other.value()


def make_rest_table(param_files):

    params_list = []

    for pf in param_files:

        # each file is a category
        category = os.path.basename(os.path.dirname(pf))

        # open the file
        try:
            f = open(pf, "r")
        except IOError:
            sys.exit("ERROR: {} does not exist".format(pf))

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

                descr = r""

                # store the current parameter in the list
                params_list.append(current_param)

            line = f.readline()


    categories = {q.category for q in params_list}

    for c in sorted(categories):

        # print the heading

        params = [q for q in params_list if q.category == c]

        clen = len(c)
        print(c)
        print(clen*"=" + "\n")

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
    for root, dirs, files in os.walk(top_dir):
        for f in files:
            if f == "_parameters":
                param_files.append(os.path.normpath("/".join([root, f])))

    make_rest_table(param_files)

if __name__ == "__main__":
    main()
