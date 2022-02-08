#!/bin/bash

OLD_NAME="sector46_20220118"
NEW_NAME="sector47_20220207"

OLD1="SECTOR = 46"
NEW1="SECTOR = 47"

OLD2="SECTOR1 = 46"
NEW2="SECTOR1 = 47"

OLD3="SECTOR2 = 46"
NEW3="SECTOR2 = 47"

# This replaces my local directory
OLD4="sector46"
NEW4="sector47"

# This replaces SPOC data directory
OLD5="sector-046"
NEW5="sector-047"

# DV report prefix
OLD6="tess2021337012456-"
NEW6="tess2021365070455-"
# DV report postfix
OLD7="-00547"
NEW7="-00560"

#LC prefix
OLD8="tess2021336043614-s0046-"
NEW8="tess2021364111932-s0047-"
OLD9="-0217-s"
NEW9="-0218-s"

# TOI Federation file
OLD10="FIXED-20220118"
NEW10="FIXED-20220207"

for name in `ls *py`; do
  echo $name
  cp -f ${name} temp.foo
  sed -e "s/${OLD_NAME}/${NEW_NAME}/g; s/${OLD1}/${NEW1}/g; s/${OLD2}/${NEW2}/g; s/${OLD3}/${NEW3}/g; s/${OLD4}/${NEW4}/g; s/${OLD5}/${NEW5}/g; s/${OLD6}/${NEW6}/g; s/${OLD7}/${NEW7}/g;  s/${OLD8}/${NEW8}/g;  s/${OLD9}/${NEW9}/g; s/${OLD10}/${NEW10}/g" temp.foo > ${name}
done

