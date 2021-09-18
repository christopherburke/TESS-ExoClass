#!/bin/bash

OLD_NAME="sector1-39_20210827"
NEW_NAME="sector41_20210917"

OLD1="SECTOR = -1"
NEW1="SECTOR = 41"

OLD2="SECTOR1 = 1"
NEW2="SECTOR1 = 41"

OLD3="SECTOR2 = 39"
NEW3="SECTOR2 = 41"

# This replaces my local directory
OLD4="sector1-39"
NEW4="sector41"

# This replaces SPOC data directory
OLD5="sector-001-013+027-039"
NEW5="sector-041"

# DV report prefix
OLD6="tess2018206190142-"
NEW6="tess2021205113501-"
# DV report postfix
OLD7="-00493"
NEW7="-00511"

#LC prefix
OLD8="tess2021146024351-s0039-"
NEW8="tess2021204101404-s0041-"
OLD9="-0210-s"
NEW9="-0212-s"

for name in `ls *py`; do
  echo $name
  cp -f ${name} temp.foo
  sed -e "s/${OLD_NAME}/${NEW_NAME}/g; s/${OLD1}/${NEW1}/g; s/${OLD2}/${NEW2}/g; s/${OLD3}/${NEW3}/g; s/${OLD4}/${NEW4}/g; s/${OLD5}/${NEW5}/g; s/${OLD6}/${NEW6}/g; s/${OLD7}/${NEW7}/g;  s/${OLD8}/${NEW8}/g;  s/${OLD9}/${NEW9}/g" temp.foo > ${name}
done

