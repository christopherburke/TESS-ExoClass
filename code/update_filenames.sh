#!/bin/bash

OLD_NAME="sector28_20200925"
NEW_NAME="sector29_20201014"

OLD1="SECTOR = 28"
NEW1="SECTOR = 29"

OLD2="SECTOR1 = 28"
NEW2="SECTOR1 = 29"

OLD3="SECTOR2 = 28"
NEW3="SECTOR2 = 29"

# This replaces my local directory
OLD4="sector28"
NEW4="sector29"

# This replaces SPOC data directory
OLD5="sector-028"
NEW5="sector-029"

# DV report prefix
OLD6="tess2020213081515-"
NEW6="tess2020239173514-"
# DV report postfix
OLD7="-00364"
NEW7="-00382"

#LC prefix
OLD8="tess2020212050318-s0028-"
NEW8="tess2020238165205-s0029-"
OLD9="-0190-s"
NEW9="-0193-s"

for name in `ls *py`; do
  echo $name
  cp -f ${name} temp.foo
  sed -e "s/${OLD_NAME}/${NEW_NAME}/g; s/${OLD1}/${NEW1}/g; s/${OLD2}/${NEW2}/g; s/${OLD3}/${NEW3}/g; s/${OLD4}/${NEW4}/g; s/${OLD5}/${NEW5}/g; s/${OLD6}/${NEW6}/g; s/${OLD7}/${NEW7}/g;  s/${OLD8}/${NEW8}/g;  s/${OLD9}/${NEW9}/g" temp.foo > ${name}
done

