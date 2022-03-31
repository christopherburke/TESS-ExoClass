#!/bin/bash

OLD_NAME="sector48_20220314"
NEW_NAME="sector1-46_20220328"

OLD1="SECTOR = 48"
NEW1="SECTOR = -1"

OLD2="SECTOR1 = 48"
NEW2="SECTOR1 = 1"

OLD3="SECTOR2 = 48"
NEW3="SECTOR2 = 46"

# This replaces my local directory
OLD4="sector48"
NEW4="sector1-46"

# This replaces SPOC data directory
OLD5="sector-001-048"
NEW5="sector-001-046"

# DV report prefix
OLD6="tess2022028101454-"
NEW6="tess2018206190142-"
# DV report postfix
OLD7="-00580"
NEW7="-00555"

#LC prefix
OLD8="tess2022027120115-s0048-"
NEW8="tess2021336043614-s0046-"
OLD9="-0219-s"
NEW9="-0217-s"

# TOI Federation file
OLD10="FIXED-20220314"
NEW10="FIXED-20220328"

for name in `ls *py`; do
  echo $name
  cp -f ${name} temp.foo
  sed -e "s/${OLD_NAME}/${NEW_NAME}/g; s/${OLD1}/${NEW1}/g; s/${OLD2}/${NEW2}/g; s/${OLD3}/${NEW3}/g; s/${OLD4}/${NEW4}/g; s/${OLD5}/${NEW5}/g; s/${OLD6}/${NEW6}/g; s/${OLD7}/${NEW7}/g;  s/${OLD8}/${NEW8}/g;  s/${OLD9}/${NEW9}/g; s/${OLD10}/${NEW10}/g" temp.foo > ${name}
done

