# NSE check tables

This test is used to have a table that has a variety of differen rho,
T, and Y_e, to test the valitidy of the function `in_nse`.

We solve for the NSE state for variety of different conditions, and update
the current state mass fraction to the NSE state to make sure we're in
NSE. Then `in_nse` is called to see whether we're in NSE or not.
