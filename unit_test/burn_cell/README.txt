~~~~~~~~~~~~~~~~~~~~~~~~~

1.  Pick your network in GNUmakefile on the line NETWORK_DIR :=
2.  Run make clean
3.  Run make
4.  Set up your inputs* file. Don't forget to note down your run prefix.
    - File naming convention: inputs_runprefix
5.  Run: ./main.Linux.gfortran.exe <your inputs file here>
6.  Run: python3 burn_cell_testing.py <your run prefix here>
7.  Input a testing prefix when prompted (this is similar to run_prefix)
8.  Copy your inputs file into a file with this naming convention
    'runprefix_testprefix_inputs.txt'
9.  Create a file named 'runprefix_testprefixes.txt' and list the test prefix 
    you chose on the first line.
10. Open your inputs file, change what you need to, and repeat steps 5 - 8. 
    The run prefix must remain the same and the test prefix must change. 
    For each iteration, add the new test prefix to a new line in 
    'runprefix_testprefixes.txt'.
    - If you decide to change networks, or anything else in the GNUmakefile,
      you must rerun steps 1-3 as well. Be sure to be consistent with 
      run prefixes and numsteps (see note below if you want to change the
      value of numsteps).
    - If you change the numsteps value in the inputs file, be sure to clear the 
      directory titled runprefix_output of all files before running step 5.
11. Now you can run burn_cell_compare_[rel, abs, 2networks].py as follows
    python3 burn_cell_compare_[rel, abs, 2networks].py <runprefix> 
    and view graphs of mass and moller fraction, temperature, and energy 
    generation rate
    - To view the created files run the following code in your terminal
      ls -ht | head 

~~~~~~~~~~~~~~~~~~~~~~~~

The runprefix_testprefixes.txt is designed so that the first test listed is the
  base for the analysis in burn_cell_compare_[rel, abs].py. This means all 
  subsequent error calculations will use the first test as a common base to 
  compare to. 

burn_cell_compare_2tests.py only compares the first two tests listed in 
  runprefix_testprefixes.txt.

~~~~~~~~~~~~~~~~~~~~~~~~

burn_cell_compare_2tests.py is the only code equiped to generate graphs for 
tests involving 2 different networks. burn_cell_compare_[rel, abs].py are not.

~~~~~~~~~~~~~~~~~~~~~~~~
