#!/bin/bash

OLD_NAME="sector44_20211122"
NEW_NAME="sector45_20211220"

OLD1="SECTOR = 44"
NEW1="SECTOR = 45"

OLD2="SECTOR1 = 44"
NEW2="SECTOR1 = 45"

OLD3="SECTOR2 = 44"
NEW3="SECTOR2 = 45"

# This replaces my local directory
OLD4="sector44"
NEW4="sector45"

# This replaces SPOC data directory
OLD5="sector-044"
NEW5="sector-045"

# DV report prefix
OLD6="tess2021285162058-"
NEW6="tess2021311000057-"
# DV report postfix
OLD7="-00532"
NEW7="-00542"

#LC prefix
OLD8="tess2021284114741-s0044-"
NEW8="tess2021310001228-s0045-"
OLD9="-0215-s"
NEW9="-0216-s"

for name in `ls *py`; do
  echo $name
  cp -f ${name} temp.foo
  sed -e "s/${OLD_NAME}/${NEW_NAME}/g; s/${OLD1}/${NEW1}/g; s/${OLD2}/${NEW2}/g; s/${OLD3}/${NEW3}/g; s/${OLD4}/${NEW4}/g; s/${OLD5}/${NEW5}/g; s/${OLD6}/${NEW6}/g; s/${OLD7}/${NEW7}/g;  s/${OLD8}/${NEW8}/g;  s/${OLD9}/${NEW9}/g" temp.foo > ${name}
done

