#!/bin/bash

OLD_NAME="sector33_20200208"
NEW_NAME="sector34_20210303"

OLD1="SECTOR = 33"
NEW1="SECTOR = 34"

OLD2="SECTOR1 = 33"
NEW2="SECTOR1 = 34"

OLD3="SECTOR2 = 33"
NEW3="SECTOR2 = 34"

# This replaces my local directory
OLD4="sector33"
NEW4="sector34"

# This replaces SPOC data directory
OLD5="sector-033"
NEW5="sector-034"

# DV report prefix
OLD6="tess2020353052510-"
NEW6="tess2021014055109-"
# DV report postfix
OLD7="-00430"
NEW7="-00444"

#LC prefix
OLD8="tess2020351194500-s0033-"
NEW8="tess2021014023720-s0034-"
OLD9="-0203-s"
NEW9="-0204-s"

for name in `ls *py`; do
  echo $name
  cp -f ${name} temp.foo
  sed -e "s/${OLD_NAME}/${NEW_NAME}/g; s/${OLD1}/${NEW1}/g; s/${OLD2}/${NEW2}/g; s/${OLD3}/${NEW3}/g; s/${OLD4}/${NEW4}/g; s/${OLD5}/${NEW5}/g; s/${OLD6}/${NEW6}/g; s/${OLD7}/${NEW7}/g;  s/${OLD8}/${NEW8}/g;  s/${OLD9}/${NEW9}/g" temp.foo > ${name}
done

