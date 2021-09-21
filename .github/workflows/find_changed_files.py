import subprocess
import sys
import argparse
from contextlib import contextmanager
import os

@contextmanager
def cd(newdir):
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)


def find_files(SHAs=None):
    diff_command = ['git', 'diff', '--name-only']

    if SHAs is not None:
        diff_command += SHAs

    stdout, stderr = subprocess.Popen(diff_command,
                                      stdout=subprocess.PIPE,
                                      stderr=subprocess.STDOUT).communicate()

    if stderr is not None:
        raise Exception('git diff encountered an error')

    files = [f for f in stdout.decode('utf-8').strip().split('\n') 
             if 'networks' in f]
    print(files)

    # see which directories contain changed files
    changed_networks = set()
    for f in files:
        d = f.split('/')[1]
        changed_networks.add(d)

    return changed_networks

def run(SHAs=None, make_options=''):

    networks = find_files(SHAs)

    if len(networks) == 0:
        networks = ['aprox13']

    GITHUB_WORKSPACE = os.environ.get('GITHUB_WORKSPACE')

    for network in networks:
        make_command = f'make {make_options} USE_MPI=FALSE USE_OMP=FALSE USE_CUDA=FALSE USE_NETWORK={network}'

        print(f'make command = {make_command}')

        with cd(f'unit_test/burn_cell'):

            print(f'making unit_test/burn_cell')

            process = subprocess.run(make_command,
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.STDOUT,
                                     shell=True)
            print(process.stdout.decode('utf-8'))
            if process.stderr is not None:
                raise Exception('make encountered an error')

    # compile test_eos as well
    make_command = f'make {make_options} USE_MPI=FALSE USE_OMP=FALSE USE_CUDA=FALSE

    with cd(f'unit_test/test_eos'):
        print(f'making unit_test/test_eos')

        process = subprocess.run(make_command,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT,
                                 shell=True)
        print(process.stdout.decode('utf-8'))
        if process.stderr is not None:
            raise Exception('make encountered an error')


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-make-options',
                        default='-j 2',
                        help='make options')
    parser.add_argument('SHAs', nargs='*', default=None,
                        help='SHAs to be compared')

    args = parser.parse_args()

    run(SHAs=args.SHAs, make_options=args.make_options)
