#!/bin/bash

OLD_NAME="sector31_20201218"
NEW_NAME="sector32_20200125"

OLD1="SECTOR = 31"
NEW1="SECTOR = 32"

OLD2="SECTOR1 = 31"
NEW2="SECTOR1 = 32"

OLD3="SECTOR2 = 31"
NEW3="SECTOR2 = 32"

# This replaces my local directory
OLD4="sector31"
NEW4="sector32"

# This replaces SPOC data directory
OLD5="sector-031"
NEW5="sector-032"

# DV report prefix
OLD6="tess2020296001112-"
NEW6="tess2020325171311-"
# DV report postfix
OLD7="-00411"
NEW7="-00419"

#LC prefix
OLD8="tess2020294194027-s0031-"
NEW8="tess2020324010417-s0032-"
OLD9="-0198-s"
NEW9="-0200-s"

for name in `ls *py`; do
  echo $name
  cp -f ${name} temp.foo
  sed -e "s/${OLD_NAME}/${NEW_NAME}/g; s/${OLD1}/${NEW1}/g; s/${OLD2}/${NEW2}/g; s/${OLD3}/${NEW3}/g; s/${OLD4}/${NEW4}/g; s/${OLD5}/${NEW5}/g; s/${OLD6}/${NEW6}/g; s/${OLD7}/${NEW7}/g;  s/${OLD8}/${NEW8}/g;  s/${OLD9}/${NEW9}/g" temp.foo > ${name}
done

