#!/bin/bash

OLD_NAME="sector32_20200125"
NEW_NAME="sector33_20200208"

OLD1="SECTOR = 32"
NEW1="SECTOR = 33"

OLD2="SECTOR1 = 32"
NEW2="SECTOR1 = 33"

OLD3="SECTOR2 = 32"
NEW3="SECTOR2 = 33"

# This replaces my local directory
OLD4="sector32"
NEW4="sector33"

# This replaces SPOC data directory
OLD5="sector-032"
NEW5="sector-033"

# DV report prefix
OLD6="tess2020325171311-"
NEW6="tess2020353052510-"
# DV report postfix
OLD7="-00419"
NEW7="-00430"

#LC prefix
OLD8="tess2020324010417-s0032-"
NEW8="tess2020351194500-s0033-"
OLD9="-0200-s"
NEW9="-0203-s"

for name in `ls *py`; do
  echo $name
  cp -f ${name} temp.foo
  sed -e "s/${OLD_NAME}/${NEW_NAME}/g; s/${OLD1}/${NEW1}/g; s/${OLD2}/${NEW2}/g; s/${OLD3}/${NEW3}/g; s/${OLD4}/${NEW4}/g; s/${OLD5}/${NEW5}/g; s/${OLD6}/${NEW6}/g; s/${OLD7}/${NEW7}/g;  s/${OLD8}/${NEW8}/g;  s/${OLD9}/${NEW9}/g" temp.foo > ${name}
done

