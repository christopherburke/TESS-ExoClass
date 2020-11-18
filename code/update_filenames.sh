#!/bin/bash

OLD_NAME="sector29_20201014"
NEW_NAME="sector30_20201117"

OLD1="SECTOR = 29"
NEW1="SECTOR = 30"

OLD2="SECTOR1 = 29"
NEW2="SECTOR1 = 30"

OLD3="SECTOR2 = 29"
NEW3="SECTOR2 = 30"

# This replaces my local directory
OLD4="sector29"
NEW4="sector30"

# This replaces SPOC data directory
OLD5="sector-029"
NEW5="sector-030"

# DV report prefix
OLD6="tess2020239173514-"
NEW6="tess2020267090513-"
# DV report postfix
OLD7="-00382"
NEW7="-00394"

#LC prefix
OLD8="tess2020238165205-s0029-"
NEW8="tess2020266004630-s0030-"
OLD9="-0193-s"
NEW9="-0195-s"

for name in `ls *py`; do
  echo $name
  cp -f ${name} temp.foo
  sed -e "s/${OLD_NAME}/${NEW_NAME}/g; s/${OLD1}/${NEW1}/g; s/${OLD2}/${NEW2}/g; s/${OLD3}/${NEW3}/g; s/${OLD4}/${NEW4}/g; s/${OLD5}/${NEW5}/g; s/${OLD6}/${NEW6}/g; s/${OLD7}/${NEW7}/g;  s/${OLD8}/${NEW8}/g;  s/${OLD9}/${NEW9}/g" temp.foo > ${name}
done

