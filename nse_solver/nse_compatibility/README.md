# NSE Network Checking Script

This is a script for checking whether integrated mass fractions of
a reaction network will eventually match with the mass fractions
calculated via NSE equations. Currently, this is only valid for networks
that do not involve weak rates.

* Note that the integration doesn't achieve NSE even though we have
USE_NSE_NET=TRUE is due to the reference size of the cell. By setting
nse.nse_dx_independet=1 will fix this.

Script will print out all burn_state cases when nse_dx_independent is
not enabled. It will only print out burn_state cases that didn't enter NSE
when nse_dx_independent=1