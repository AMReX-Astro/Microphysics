import subprocess
import sys
import os
import re


def pow_to_powi(text):
    # Finds all possible std::pow(x, n) where n is a potential integer
    # to amrex::Math::powi<n>(x)

    # Check for all positive and negative integer, whole float numbers
    # with and without _rt
    pattern = r"std::pow\(([^,]+),\s*(-?(?:\d+\.0*_rt?|\d))\)"

    def replacement(match):
        x = match.group(1)
        n = match.group(2)

        # Only extracts out the integer part before decimal point
        n = n.split('.')[0] if '.' in n else n
        return f"amrex::Math::powi<{n}>({x})"

    return re.sub(pattern, replacement, text)

def process_content(dir_path):
    # This function processes all text in the given directory
    for root, dirs, filenames in os.walk(dir_path):
        for filename in filenames:
            if filename.endswith(".H") or filename.endswith(".cpp"):
                filepath = os.path.join(root, filename)

                with open(filepath, 'r') as file:
                    content = file.read()

                updated_content = pow_to_powi(content)

                with open(filepath, 'w') as file:
                    file.write(updated_content)

def git_diff():
    # Run git diff to see if there are any changes made

    git_diff_output = subprocess.run(['git', 'diff', '--color=always'],
                                     capture_output=True, text=True)

    # Print out suggested change and raise error after detecting modification
    if git_diff_output.stdout:
        print("Detected potential usage to replace std::pow" +
              "with integer powers via amrex::Math::powi\n")
        print("Below are the suggested change:\n")
        print(git_diff_output.stdout)

        raise RuntimeError("Changes detected after modification")

if __name__ == '__main__':

    # Get directory paths
    directory_paths = sys.argv[1:]

    # Modify the std::pow -> amrex::Math::powi if needed
    for dir_path in directory_paths:
        process_content(dir_path)

    # Give suggested change if there are any modifications made
    git_diff()
