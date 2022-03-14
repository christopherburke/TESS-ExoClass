#!/bin/bash

OLD_NAME="sector47_20220207"
NEW_NAME="sector48_20220314"

OLD1="SECTOR = 47"
NEW1="SECTOR = 48"

OLD2="SECTOR1 = 47"
NEW2="SECTOR1 = 48"

OLD3="SECTOR2 = 47"
NEW3="SECTOR2 = 48"

# This replaces my local directory
OLD4="sector47"
NEW4="sector48"

# This replaces SPOC data directory
OLD5="sector-047"
NEW5="sector-048"

# DV report prefix
OLD6="tess2021365070455-"
NEW6="tess2022028101454-"
# DV report postfix
OLD7="-00560"
NEW7="-00580"

#LC prefix
OLD8="tess2021364111932-s0047-"
NEW8="tess2022027120115-s0048-"
OLD9="-0218-s"
NEW9="-0219-s"

# TOI Federation file
OLD10="FIXED-20220207"
NEW10="FIXED-20220314"

for name in `ls *py`; do
  echo $name
  cp -f ${name} temp.foo
  sed -e "s/${OLD_NAME}/${NEW_NAME}/g; s/${OLD1}/${NEW1}/g; s/${OLD2}/${NEW2}/g; s/${OLD3}/${NEW3}/g; s/${OLD4}/${NEW4}/g; s/${OLD5}/${NEW5}/g; s/${OLD6}/${NEW6}/g; s/${OLD7}/${NEW7}/g;  s/${OLD8}/${NEW8}/g;  s/${OLD9}/${NEW9}/g; s/${OLD10}/${NEW10}/g" temp.foo > ${name}
done

