#!/usr/bin/env python3

"""
A class that manages the runtime parameters.  The idea is that this
can be used by all of the parsing scripts to write the C++ and
Fortran files that manage runtime parameters
"""
import sys

class Param:
    """ the basic parameter class.  For each parameter, we hold the name,
        type, and default.  For some parameters, we also take a second
        value of the default, for use in debug mode (delimited via
        #ifdef AMREX_DEBUG)

    """

    def __init__(self, name, dtype, default,
                 cpp_var_name=None,
                 namespace=None,
                 debug_default=None,
                 in_fortran=0,
                 priority=0,
                 ifdef=None):

        self.name = name
        self.dtype = dtype
        self.default = default
        self.cpp_var_name = cpp_var_name
        self.priority = priority

        self.namespace = namespace

        self.debug_default = debug_default
        self.in_fortran = in_fortran

        if ifdef == "None":
            self.ifdef = None
        else:
            self.ifdef = ifdef

        if self.namespace is None or self.namespace == "":
            self.nm_pre = ""
        else:
            self.nm_pre = f"{self.namespace}::"

    def get_declare_string(self):
        """this is the line that goes into, e.g., castro_declares.H included
        into Castro.cpp"""

        if self.dtype == "int":
            tstr = f"AMREX_GPU_MANAGED int         {self.nm_pre}{self.cpp_var_name}"
        elif self.dtype == "bool":
            tstr = f"AMREX_GPU_MANAGED bool        {self.nm_pre}{self.cpp_var_name}"
        elif self.dtype == "real":
            tstr = f"AMREX_GPU_MANAGED amrex::Real {self.nm_pre}{self.cpp_var_name}"
        elif self.dtype == "string":
            tstr = f"std::string {self.nm_pre}{self.cpp_var_name}"
        else:
            sys.exit(f"invalid data type for parameter {self.name}")

        return rf"{tstr};\n"

    def get_default_string(self):
        """this is the line that goes into, e.g., castro_declares.H included
        into Castro.cpp"""

        ostr = ""

        if not self.debug_default is None:
            ostr += "#ifdef AMREX_DEBUG\n"
            ostr += f"{self.nm_pre}{self.cpp_var_name} = {self.debug_default};\n"
            ostr += "#else\n"
            ostr += f"{self.nm_pre}{self.cpp_var_name} = {self.default};\n"
            ostr += "#endif\n"
        else:
            ostr += f"{self.nm_pre}{self.cpp_var_name} = {self.default};\n"

        return ostr

    def get_f90_default_string(self):
        """this is the line that goes into, e.g., set_castro_method_params()
        to set the default value of the variable"""

        ostr = ""

        # convert to the double precision notation Fortran knows
        # if the parameter is already of the form "#.e###" then
        # it is easy as swapping out "e" for "d"; if it is a number
        # like 0.1 without a format specifier, then add a d0 to it
        # because the C++ will read it in that way and we want to
        # give identical results (at least to within roundoff)

        if self.debug_default is not None:
            debug_default = self.debug_default
            if self.dtype == "real":
                if "d" in debug_default:
                    debug_default = debug_default.replace("d", "e")
                debug_default += "_rt"

        default = self.default
        if self.dtype == "real":
            if "d" in default:
                default = default.replace("d", "e")
            default += "_rt"

        name = self.name

        # for a character, we need to allocate its length.  We allocate
        # to 1, and the Fortran parmparse will resize
        if self.dtype == "string":
            ostr += "    allocate(character(len=1)::{})\n".format(name)
        else:
            ostr += "    allocate({})\n".format(name)

        if not self.debug_default is None:
            ostr += "#ifdef AMREX_DEBUG\n"
            ostr += f"    {name} = {debug_default};\n"
            ostr += "#else\n"
            ostr += f"    {name} = {default};\n"
            ostr += "#endif\n"
        else:
            ostr += f"    {name} = {default};\n"

        return ostr

    def get_query_string(self, language):
        """this is the line that queries the ParmParse object to get
        the value of the runtime parameter from the inputs file.
        This goes into, e.g., castro_queries.H included into Castro.cpp"""

        ostr = ""
        if language == "C++":
            ostr += f"pp.query(\"{self.name}\", {self.nm_pre}{self.cpp_var_name});\n"
        elif language == "F90":
            ostr += f"    call pp%query(\"{self.name}\", {self.name})\n"
        else:
            sys.exit("invalid language choice in get_query_string")

        return ostr

    def default_format(self):
        """return the variable in a format that it can be recognized in C++
        code--in particular, preserve the quotes for strings"""
        if self.dtype == "string":
            return '{}'.format(self.default)

        return self.default

    def get_job_info_test(self):
        """this is the output in C++ in the job_info writing"""

        ostr = (
            f'jobInfoFile << ({self.nm_pre}{self.cpp_var_name} == {self.default_format()} ? "    "' +
            ': "[*] ") << "{self.namespace}.{self.cpp_var_name} = "' +
            '<< {self.np_pre}{self.cpp_var_name} << std::endl;\n')

        return ostr

    def get_decl_string(self):
        """this is the line that goes into, e.g., castro_params.H included
        into Castro.H"""

        tstr = ""

        if self.dtype == "int":
            tstr = f"extern AMREX_GPU_MANAGED int {self.cpp_var_name};\n"
        elif self.dtype == "bool":
            tstr = f"extern AMREX_GPU_MANAGED bool {self.cpp_var_name};\n"
        elif self.dtype == "real":
            tstr = f"extern AMREX_GPU_MANAGED amrex::Real {self.cpp_var_name};\n"
        elif self.dtype == "string":
            tstr = f"extern std::string {self.cpp_var_name};\n"
        else:
            sys.exit(f"invalid data type for parameter {self.name}")

        return tstr

    def get_f90_decl_string(self):
        """this is the line that goes into, e.g., meth_params_nd.F90"""

        if not self.in_fortran:
            return None

        if self.dtype == "int":
            tstr = f"integer,  allocatable, save :: {self.name}\n"
        elif self.dtype == "real":
            tstr = f"real(rt), allocatable, save :: {self.name}\n"
        elif self.dtype == "logical":
            tstr = f"logical,  allocatable, save :: {self.name}\n"
        elif self.dtype == "string":
            tstr = f"character (len=:), allocatable, save :: {self.name}\n"
            print(f"warning: string parameter {self.name} will not be available on the GPU")
        else:
            sys.exit(f"unsupported datatype for Fortran: {self.name}")

        return tstr

    def __lt__(self, other):
        return self.priority < other.priority
