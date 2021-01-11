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
                 skip_namespace_in_declare=False,
                 debug_default=None,
                 in_fortran=0,
                 priority=0,
                 size=1,
                 in_namelist = False,
                 ifdef=None):

        self.name = name

        self.dtype = dtype
        if self.dtype == "character":
            self.dtype = "string"

        self.default = default
        self.size = size
        self.cpp_var_name = cpp_var_name
        if self.cpp_var_name is None:
            self.cpp_var_name = self.name

        self.priority = priority

        self.in_namelist = in_namelist

        if namespace is not None:
            self.namespace = namespace.strip()
        else:
            self.namespace = namespace

        # if this is true, then we use the namespace when we read the var
        # (e.g., via ParmParse), but we do not declare the C++
        # parameter to be in a namespace
        self.skip_namespace_in_declare = skip_namespace_in_declare

        self.debug_default = debug_default
        self.in_fortran = in_fortran

        if ifdef == "None":
            self.ifdef = None
        else:
            self.ifdef = ifdef

        if self.namespace is None or self.namespace == "" or self.skip_namespace_in_declare:
            self.nm_pre = ""
        else:
            self.nm_pre = f"{self.namespace}::"

    def get_cxx_decl(self):
        """ get the Fortran 90 declaration """
        if self.dtype == "real":
            return "amrex::Real"
        elif self.dtype == "string":
            return "std::string"
        elif self.dtype == "bool":
            return "bool"

        return "int"

    def get_declare_string(self):
        """this is the line that goes into, e.g., castro_declares.H included
        into Castro.cpp"""

        if self.dtype != "string":
            tstr = f"AMREX_GPU_MANAGED {self.get_cxx_decl()} {self.nm_pre}{self.cpp_var_name}"
        elif self.dtype == "string":
            tstr = f"std::string {self.nm_pre}{self.cpp_var_name}"
        else:
            sys.exit(f"invalid data type for parameter {self.name}")

        return f"{tstr};\n"

    def get_decl_string(self):
        """this is the line that goes into, e.g., castro_params.H included
        into Castro.H"""

        tstr = ""

        if self.dtype != "string":
            tstr = f"extern AMREX_GPU_MANAGED {self.get_cxx_decl()} {self.cpp_var_name};\n"
        elif self.dtype == "string":
            tstr = f"extern std::string {self.cpp_var_name};\n"
        else:
            sys.exit(f"invalid data type for parameter {self.name}")

        return tstr

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
        """For Fortran, allocate the variable and set the default.  This is typically
        triggered by @@allocation@@ in the template"""

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
                if not debug_default.endswith("_rt"):
                    debug_default += "_rt"

        default = self.default
        if self.dtype == "real":
            if "d" in default:
                default = default.replace("d", "e")
            if not default.endswith("_rt"):
                default += "_rt"
        elif self.dtype == "bool":
            if default == "true":
                default = ".true."
            elif default == "false":
                default = ".false."
        name = self.name

        if self.dtype != "string":
            if self.is_array():
                ostr += f"    allocate({name}({self.size}))\n"
            else:
                ostr += f"    allocate({name})\n"

        if self.is_array():
            name_set = f"{name}(:)"
        else:
            name_set = f"{name}"

        if not self.debug_default is None:
            ostr += "#ifdef AMREX_DEBUG\n"
            ostr += f"    {name_set} = {debug_default}\n"
            ostr += "#else\n"
            ostr += f"    {name_set} = {default}\n"
            ostr += "#endif\n"
        else:
            ostr += f"    {name_set} = {default}\n"

        return ostr

    def get_query_string(self, language):
        """this is the line that queries the ParmParse object to get
        the value of the runtime parameter from the inputs file.
        This goes into, e.g., castro_queries.H included into Castro.cpp"""

        ostr = ""
        if language == "C++":
            ostr += f"pp.query(\"{self.name}\", {self.nm_pre}{self.cpp_var_name});\n"
        elif language == "F90":
            if self.dtype == "string":
                ostr += "    allocate(character(len=1) :: dummy_string_param)\n"
                ostr += "    dummy_string_param = \"\"\n"
                ostr += f"    call pp%query(\"{self.name}\", dummy_string_param)\n"
                ostr += f"    if (dummy_string_param /= \"\") {self.name} = dummy_string_param\n"
                ostr += "    deallocate(dummy_string_param)\n"
            else:
                ostr += f"    call pp%query(\"{self.name}\", {self.name})\n"
        else:
            sys.exit("invalid language choice in get_query_string")

        return ostr

    def default_format(self):
        """return the variable in a format that it can be recognized in C++
        code--in particular, preserve the quotes for strings"""
        if self.dtype == "string":
            return f'{self.default}'

        return self.default

    def get_job_info_test(self):
        """this is the output in C++ in the job_info writing"""

        ostr = (
            f'jobInfoFile << ({self.nm_pre}{self.cpp_var_name} == {self.default_format()} ? "    "' +
            f': "[*] ") << "{self.namespace}.{self.cpp_var_name} = "' +
            f'<< {self.nm_pre}{self.cpp_var_name} << std::endl;\n')

        return ostr

    def get_f90_decl(self):
        """get the Fortran 90 declaration.  This is intended for the case
        where the parameters are managed in Fortran"""
        if self.dtype == "real":
            return "real (kind=rt)"
        elif self.dtype == "string":
            return "character (len=256)"
        elif self.dtype == "int":
            return "integer"
        elif self.dtype == "bool":
            return "logical"
        return self.dtype

    def get_f90_decl_string(self):
        """this is the line that goes into, e.g., meth_params_nd.F90"""

        if not self.in_fortran:
            return None

        if self.dtype != "string":
            if self.is_array():
                tstr = f"{self.get_f90_decl()},  allocatable, save :: {self.name}(:)\n"
            else:
                tstr = f"{self.get_f90_decl()},  allocatable, save :: {self.name}\n"
        elif self.dtype == "string":
            if self.is_array():
                sys.exit("error: cannot have a character array")
            else:
                tstr = f"character (len=256) :: {self.name}\n"
            print(f"warning: string parameter {self.name} will not be available on the GPU")
        else:
            sys.exit(f"unsupported datatype for Fortran: {self.name}")

        return tstr

    def get_f90_get_function(self):
        """this returns the "getter" function in Fortran that is called from C++"""
        ostr = ""
        if self.dtype == "string":
            ostr += f"  subroutine get_f90_{self.name}_len(slen) bind(C, name=\"get_f90_{self.name}_len\")\n"
            ostr +=  "     integer, intent(inout) :: slen\n"
            ostr += f"     slen = len(trim({self.name}))\n"
            ostr += f"  end subroutine get_f90_{self.name}_len\n\n"

            ostr += f"  subroutine get_f90_{self.name}({self.name}_in) bind(C, name=\"get_f90_{self.name}\")\n"
            ostr += f"     character(kind=c_char) :: {self.name}_in(*)\n"
            ostr +=  "     integer :: n\n"
            ostr += f"     do n = 1, len(trim({self.name}))\n"
            ostr += f"        {self.name}_in(n:n) = {self.name}(n:n)\n"
            ostr +=  "     end do\n"
            ostr += f"     {self.name}_in(len(trim({self.name}))+1) = char(0)\n"
            ostr += f"  end subroutine get_f90_{self.name}\n\n"

        elif self.dtype == "bool":
            # F90 logicals are integers in C++
            ostr += f"  subroutine get_f90_{self.name}({self.name}_in) bind(C, name=\"get_f90_{self.name}\")\n"
            ostr += f"     integer, intent(inout) :: {self.name}_in\n"
            ostr += f"     {self.name}_in = 0\n"
            ostr += f"     if ({self.name}) then\n"
            ostr += f"        {self.name}_in = 1\n"
            ostr +=  "     endif\n"
            ostr += f"  end subroutine get_f90_{self.name}\n\n"

        else:
            ostr += f"  subroutine get_f90_{self.name}({self.name}_in) bind(C, name=\"get_f90_{self.name}\")\n"
            if self.is_array():
                ostr += f"     {self.get_f90_decl()}, intent(inout) :: {self.name}_in({self.size})\n"
            else:
                ostr += f"     {self.get_f90_decl()}, intent(inout) :: {self.name}_in\n"
            ostr += f"     {self.name}_in = {self.name}\n"
            ostr += f"  end subroutine get_f90_{self.name}\n\n"

        return ostr

    def is_array(self):
        """return true if the parameter is an array"""
        try:
            isize = int(self.size)
        except ValueError:
            return True
        else:
            if isize == 1:
                return False
            return True

    def __lt__(self, other):
        return self.priority < other.priority

    def __str__(self):
        return self.name
