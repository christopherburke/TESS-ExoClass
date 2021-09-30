#!/bin/bash

OLD_NAME="sector41_20210917"
NEW_NAME="sector42_20210930"

OLD1="SECTOR = 41"
NEW1="SECTOR = 42"

OLD2="SECTOR1 = 41"
NEW2="SECTOR1 = 42"

OLD3="SECTOR2 = 41"
NEW3="SECTOR2 = 42"

# This replaces my local directory
OLD4="sector41"
NEW4="sector42"

# This replaces SPOC data directory
OLD5="sector-041"
NEW5="sector-042"

# DV report prefix
OLD6="tess2021205113501-"
NEW6="tess2021233042500-"
# DV report postfix
OLD7="-00511"
NEW7="-00517"

#LC prefix
OLD8="tess2021204101404-s0041-"
NEW8="tess2021232031932-s0042-"
OLD9="-0212-s"
NEW9="-0213-s"

for name in `ls *py`; do
  echo $name
  cp -f ${name} temp.foo
  sed -e "s/${OLD_NAME}/${NEW_NAME}/g; s/${OLD1}/${NEW1}/g; s/${OLD2}/${NEW2}/g; s/${OLD3}/${NEW3}/g; s/${OLD4}/${NEW4}/g; s/${OLD5}/${NEW5}/g; s/${OLD6}/${NEW6}/g; s/${OLD7}/${NEW7}/g;  s/${OLD8}/${NEW8}/g;  s/${OLD9}/${NEW9}/g" temp.foo > ${name}
done

