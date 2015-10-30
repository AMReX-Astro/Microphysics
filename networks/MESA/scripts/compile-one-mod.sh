#!/bin/bash
#
# Simple Script to build a single MESA module
#
#    -does not need to be run from $(MESA_DIR)
#    -can specify relative path to mesa module dir
#

function check_okay {
	if [ $? -ne 0 ]
	then
		exit 1
	fi
}

function do_one {
        echo
        echo '-----------------------------------------'
        echo cd $1
        echo '-----------------------------------------'
	cd $1
	check_okay

        echo
        echo 'Begin: '`basename $1`
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

        # if making the net module, we must copy the net_burn.mod and 
        # net_initialize.mod files to the include directory explicitly:
        mod=`basename $1`
        if [ $mod == 'net' ]; then
           echo
           echo "copying net_burn and net_initialize"
           echo
           cp make/net_burn.mod ../include
           cp make/net_initialize.mod ../include
        fi

        # return to starting directory
        echo
        echo '-----------------------------------------'
        echo 'Return to calling directory:'
        echo cd $2
        echo '-----------------------------------------'
        cd $2

        echo
        echo "###############################################"
        echo
        echo "        SUCCESS:" $1
        echo
        echo "###############################################"
        echo
}

if [ $# -lt 1 ]; then
   echo
   echo "Compile only one MESA Module"
   echo
   echo "Usage:"
   echo
   echo "   ./compile-one-mod.sh <path-to-module-directory>"
   echo
else
   # get starting directory
   cwd=`pwd`
   do_one $1 $cwd
fi

