#!/usr/bin/env python
import sys
import os

# tex format stuff
Mheader=r"""
\label{ch:parameters}


%%%%%%%%%%%%%%%%
% symbol table
%%%%%%%%%%%%%%%%

\begin{landscape}
"""

header=r"""
{\small

\renewcommand{\arraystretch}{1.5}
%
\begin{center}
\begin{longtable}{|l|p{5.25in}|l|}
\caption[@@catname@@]{@@catname@@} \label{table: @@sanitizedcatname@@ runtime} \\
%
\hline \multicolumn{1}{|c|}{\textbf{parameter}} & 
       \multicolumn{1}{ c|}{\textbf{description}} & 
       \multicolumn{1}{ c|}{\textbf{default value}} \\ \hline 
\endfirsthead

\multicolumn{3}{c}%
{{\tablename\ \thetable{}---continued}} \\
\hline \multicolumn{1}{|c|}{\textbf{parameter}} & 
       \multicolumn{1}{ c|}{\textbf{description}} & 
       \multicolumn{1}{ c|}{\textbf{default value}} \\ \hline 
\endhead

\multicolumn{3}{|r|}{{\em continued on next page}} \\ \hline
\endfoot

\hline 
\endlastfoot

"""

footer=r"""

\end{longtable}
\end{center}

} % ends \small
"""

Mfooter=r"""
\end{landscape}

%

"""


class Parameter(object):
    # container class for the parameters

    def __init__(self):
        self.var=""
        self.default=""
        self.description=[]
        self.category=""

    def value(self):
        """ the value is what we sort based on """
        return self.category + "." + self.var

    def __cmp__(self, other):
        return cmp(self.value(), other.value())


def make_tex_table(param_files):

    params_list=[]

    for pf in param_files:

        # each file is a category
        category = os.path.basename(os.path.dirname(pf)).replace("_","\_")

        # open the file
        try: f = open(pf, "r")
        except IOError:
            sys.exit("ERROR: {} does not exist".format(pf))

        descr=r""


        # read in the file 

        # sometimes we have a descriptive header -- skip all lines
        # before the first empty line
        found_first_param = False

        line = f.readline()
        while line:

            if not found_first_param:
                if line.strip() == "":
                    found_first_param = True
                line = f.readline()
                continue

            if line.strip() == "": 
                line = f.readline()
                continue

            if line.startswith("#------"):
                line = f.readline()
                continue

        
            # find the description
            if line.startswith("#"):
                # handle descriptions here
                descr += line[1:].rstrip().replace("@@",r"\newline")
                line = f.readline()
                continue

            else:
                current_param = Parameter()
                line_list = line.split()

                current_param.var = line_list[0]
                current_param.default = line_list[2].replace("_","\_")
                current_param.description = descr
                current_param.category = category

                descr=r""

        
                # store the current parameter in the list
                params_list.append(current_param)
                
            line = f.readline()

    
    # dump the main header
    print Mheader

    # sort the parameters and dump them in latex-fashion.  Group things by category
    current_category = ""
    start = 1

    for param in sorted(params_list):

        if not param.category == current_category:
            if not start == 1:
                print footer

            current_category = param.category
            odd = 1
            sanitized_cat_name = param.category.replace("\\", "")
            cat_header = header.replace("@@catname@@", param.category + " parameters.")
            cat_header = cat_header.replace("@@sanitizedcatname@@", sanitized_cat_name)
            print cat_header
            start = 0

        if odd == 1:
            print "\\rowcolor{tableShade}"
            odd = 0
        else:
            odd = 1

        print "\\verb= ", \
            param.var, \
            " = & ", \
            param.description, \
            " & ", \
            param.default, \
            r"\\"

    # dump the footer
    print footer
    print Mfooter

if __name__ == "__main__":


    # find all of the _parameter files 
    top_dir = "../../"

    param_files = []
    for root, dirs, files in os.walk(top_dir):
        for f in files:
            if f == "_parameters":
                param_files.append(os.path.normpath("/".join([root, f])))

    make_tex_table(param_files)
