#!/bin/sh

test ! -z $2 || { echo "Usage: $0 <test no.> <total no. of tests>"; exit 255; }

testnr=$1; 
total=$2;
command=`grep -m 2 --after-context=1 "${testnr}/${total}" Testing/Temporary/LastTest.log | grep -v "${testnr}/${total}" | sed -e "s/Command: //" | xargs`; $command
