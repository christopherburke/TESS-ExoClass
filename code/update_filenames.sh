#!/bin/bash

OLD_NAME="sector30_20201117"
NEW_NAME="sector31_20201218"

OLD1="SECTOR = 30"
NEW1="SECTOR = 31"

OLD2="SECTOR1 = 30"
NEW2="SECTOR1 = 31"

OLD3="SECTOR2 = 30"
NEW3="SECTOR2 = 31"

# This replaces my local directory
OLD4="sector30"
NEW4="sector31"

# This replaces SPOC data directory
OLD5="sector-030"
NEW5="sector-031"

# DV report prefix
OLD6="tess2020267090513-"
NEW6="tess2020296001112-"
# DV report postfix
OLD7="-00394"
NEW7="-00411"

#LC prefix
OLD8="tess2020266004630-s0030-"
NEW8="tess2020294194027-s0031-"
OLD9="-0195-s"
NEW9="-0198-s"

for name in `ls *py`; do
  echo $name
  cp -f ${name} temp.foo
  sed -e "s/${OLD_NAME}/${NEW_NAME}/g; s/${OLD1}/${NEW1}/g; s/${OLD2}/${NEW2}/g; s/${OLD3}/${NEW3}/g; s/${OLD4}/${NEW4}/g; s/${OLD5}/${NEW5}/g; s/${OLD6}/${NEW6}/g; s/${OLD7}/${NEW7}/g;  s/${OLD8}/${NEW8}/g;  s/${OLD9}/${NEW9}/g" temp.foo > ${name}
done

