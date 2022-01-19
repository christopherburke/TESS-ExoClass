#!/bin/bash

OLD_NAME="sector45_20211220"
NEW_NAME="sector46_20220118"

OLD1="SECTOR = 45"
NEW1="SECTOR = 46"

OLD2="SECTOR1 = 45"
NEW2="SECTOR1 = 46"

OLD3="SECTOR2 = 45"
NEW3="SECTOR2 = 46"

# This replaces my local directory
OLD4="sector45"
NEW4="sector46"

# This replaces SPOC data directory
OLD5="sector-045"
NEW5="sector-046"

# DV report prefix
OLD6="tess2021311000057-"
NEW6="tess2021337012456-"
# DV report postfix
OLD7="-00542"
NEW7="-00547"

#LC prefix
OLD8="tess2021310001228-s0045-"
NEW8="tess2021336043614-s0046-"
OLD9="-0216-s"
NEW9="-0217-s"

for name in `ls *py`; do
  echo $name
  cp -f ${name} temp.foo
  sed -e "s/${OLD_NAME}/${NEW_NAME}/g; s/${OLD1}/${NEW1}/g; s/${OLD2}/${NEW2}/g; s/${OLD3}/${NEW3}/g; s/${OLD4}/${NEW4}/g; s/${OLD5}/${NEW5}/g; s/${OLD6}/${NEW6}/g; s/${OLD7}/${NEW7}/g;  s/${OLD8}/${NEW8}/g;  s/${OLD9}/${NEW9}/g" temp.foo > ${name}
done

