#!/bin/bash

function check_okay {
   if [ $? -ne 0 ]
   then
      exit 1
   fi
}

function do_one {
 
   for dir in "$@"
   do
     cd $dir
     check_okay
     echo
     echo 'Begin: '$dir
     echo '$$$$$$$$$$$$$$$$$  CLEAN  $$$$$$$$$$$$$$$$$'
     ./clean
     check_okay
     echo
     echo '$$$$$$$$$$$$$$$$$  MAKE  $$$$$$$$$$$$$$$$$'
     ./mk
     check_okay
     echo
     echo '$$$$$$$$$$$$$$$$$  EXPORT  $$$$$$$$$$$$$$$$$'
     ./export
     check_okay

     echo
     echo "%%%%%%%%%%%%%%%%%%%%%%%"
     echo "%%%                 %%%"
     echo "   SUCCESS:" $dir
     echo "%%%                 %%%"
     echo "%%%%%%%%%%%%%%%%%%%%%%%"
     echo
     cd ..
   done
}

#----------------------------------------------------------
# version 4088 order:
#----------------------------------------------------------
#do_one adipls alert atm chem colors const diffusion eos
#do_one interp_1d interp_2d ionization kap mlt mtx net neu
#do_one num rates reaclib screen star

# utils does not have a 'mk' script
#if [ ! -e utils/mk ]; then
 # cp star/mk utils/
#fi

#do_one utils weaklib
#----------------------------------------------------------


#----------------------------------------------------------
# version 4411 order:
#----------------------------------------------------------
# utils does not have a 'mk' script
#if [ ! -e utils/mk ]; then
#  cp star/mk utils/
#fi
#do_one const alert utils mtx num interp_1d interp_2d chem
#do_one colors eos kap weaklib screen reaclib rates neu net
#do_one mlt ionization diffusion atm
# these are unecessary:
#do_one sample star adipls
#----------------------------------------------------------


#----------------------------------------------------------
# version 4849 order:
#----------------------------------------------------------
# utils does not have a 'mk' script
#if [ ! -e utils/mk ]; then
#  cp star/mk utils/
#fi
#
#do_one const alert utils mtx num interp_1d interp_2d chem
#do_one colors eos kap weaklib screen reaclib rates neu net
#do_one nse mlt ionization diffusion atm
# these are unecessary:
#do_one sample
#do_one star
#
#do_one adipls
#----------------------------------------------------------


#----------------------------------------------------------
# version 4942 order:
#----------------------------------------------------------
# utils does not have a 'mk' script
if [ ! -e utils/mk ]; then
  cp star/mk utils/
fi
do_one const alert utils mtx num interp_1d interp_2d chem
do_one colors eos kap weaklib screen reaclib rates neu net
do_one nse mlt ionization atm
# these are unecessary:
#do_one sample
#do_one star

do_one adipls
#----------------------------------------------------------


#----------------------------------------------------------
# version 5118 order
#----------------------------------------------------------
# utils does not have a 'mk' script
#if [ ! -e utils/mk ]; then
#  cp star/mk utils/
#fi
#do_one const utils mtx num interp_1d interp_2d chem colors
#do_one eos kap weaklib ecapture screen reaclib rates neu
#do_one net nse mlt ionization atm 
# these are unecessary:
#do_one sample star
#
#do_one gyre adipls
#----------------------------------------------------------


# must explicitly copy two modules used for SDC algorithm
cp net/make/net_burn.mod include
cp net/make/net_initialize.mod include

echo
echo
echo
echo "###############################################"
echo
echo "                  COMPLETE"
echo
echo "###############################################"
echo
