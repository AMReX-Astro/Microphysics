"""
A code checker to make sure that there are no simple errors.
"""

import re
import pytest
import os
from pathlib import Path

def find_fortran_files():
    # find Microphysics Fortran source files 
    try:
        microphysics_home = os.environ['MICROPHYSICS_HOME']
    except KeyError:
        # this assumes this file is run from the Microphysics/util/code_checker directory
        microphysics_home = '../..'

    p = Path(microphysics_home)
    return p.glob(r'**/*.[fF]90')

def idfunc(argvalue):
    return '/'.join(str(argvalue).split('/')[-2:])

def pytest_generate_tests(metafunc):
    files = find_fortran_files()

    metafunc.parametrize('filename', files, ids=idfunc, scope="module")


def test_check_rt(filename):
    """
    make sure that all of the numerical constants use _rt and are defined as real(rt)
    """
    ignore_dirs = ['util', 'tmp_build_dir', 't', 'python_library']
    if any([f'/{s}/' in str(filename) for s in ignore_dirs]):
        return 

    with open(filename, 'r') as file_dat:

        r = re.compile(r'\W(\d*\.\d*[de]?-?\+?\d+(?:_rt)?)')
        rt = re.compile(r'(\d*\.\d*e?-?\+?\d+_rt)')

        double_prec = re.compile(r'(double precision)')

        try:
            for l in file_dat:
                for m in re.finditer(r, l.split('!')[0]):
                    assert re.fullmatch(rt, m.group(1)) is not None

                assert re.search(double_prec, l.split('!')[0]) is None
        except UnicodeDecodeError:
            return 