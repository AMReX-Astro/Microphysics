#!/bin/sh

RUNPARAMS="
unit_test.tmax=1.e-4
unit_test.nsteps=1000
unit_test.density=1.e8
unit_test.temperature=3.e9"

./main3d.gnu.ex inputs_aprox13 ${RUNPARAMS} integrator.rtol_spec=1.0e-3 integrator.rtol_enuc=1.0e-3 integrator.atol_spec=1.0e-3 integrator.atol_enuc=1.0e-3
mv state_over_time.txt state_over_time_tol_a.txt

./main3d.gnu.ex inputs_aprox13 ${RUNPARAMS} integrator.rtol_spec=1.0e-5 integrator.rtol_enuc=1.0e-5 integrator.atol_spec=1.0e-5 integrator.atol_enuc=1.0e-5
mv state_over_time.txt state_over_time_tol_b.txt

./main3d.gnu.ex inputs_aprox13 ${RUNPARAMS} integrator.rtol_spec=1.0e-8 integrator.rtol_enuc=1.0e-8 integrator.atol_spec=1.0e-8 integrator.atol_enuc=1.0e-8
mv state_over_time.txt state_over_time_tol_c.txt

./main3d.gnu.ex inputs_aprox13 ${RUNPARAMS} integrator.rtol_spec=1.0e-12 integrator.rtol_enuc=1.0e-12 integrator.atol_spec=1.0e-12 integrator.atol_enuc=1.0e-12
mv state_over_time.txt state_over_time_tol_d.txt

