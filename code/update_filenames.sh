#!/bin/bash

OLD_NAME="sector1-46_20220328"
NEW_NAME="sector-49_20220424"

OLD1="SECTOR = -1"
NEW1="SECTOR = 49"

OLD2="SECTOR1 = 1"
NEW2="SECTOR1 = 49"

OLD3="SECTOR2 = 46"
NEW3="SECTOR2 = 49"

# This replaces my local directory
OLD4="sector1-46"
NEW4="sector49"

# This replaces SPOC data directory
OLD5="sector-001-046"
NEW5="sector-049"

# DV report prefix
OLD6="tess2018206190142-"
NEW6="tess2022057231053-"
# DV report postfix
OLD7="-00555"
NEW7="-00585"

#LC prefix
OLD8="tess2021336043614-s0046-"
NEW8="tess2022057073128-s0049-"
OLD9="-0217-s"
NEW9="-0221-s"

# TOI Federation file
OLD10="FIXED-20220318"
NEW10="FIXED-20220424"

for name in `ls *py`; do
  echo $name
  cp -f ${name} temp.foo
  sed -e "s/${OLD_NAME}/${NEW_NAME}/g; s/${OLD1}/${NEW1}/g; s/${OLD2}/${NEW2}/g; s/${OLD3}/${NEW3}/g; s/${OLD4}/${NEW4}/g; s/${OLD5}/${NEW5}/g; s/${OLD6}/${NEW6}/g; s/${OLD7}/${NEW7}/g;  s/${OLD8}/${NEW8}/g;  s/${OLD9}/${NEW9}/g; s/${OLD10}/${NEW10}/g" temp.foo > ${name}
done

