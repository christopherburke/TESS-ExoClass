#!/bin/bash

OLD_NAME="sector43_20211101"
NEW_NAME="sector44_20211122"

OLD1="SECTOR = 43"
NEW1="SECTOR = 44"

OLD2="SECTOR1 = 43"
NEW2="SECTOR1 = 44"

OLD3="SECTOR2 = 43"
NEW3="SECTOR2 = 44"

# This replaces my local directory
OLD4="sector43"
NEW4="sector44"

# This replaces SPOC data directory
OLD5="sector-043"
NEW5="sector-044"

# DV report prefix
OLD6="tess2021259155059-"
NEW6="tess2021285162058-"
# DV report postfix
OLD7="-00522"
NEW7="-00532"

#LC prefix
OLD8="tess2021258175143-s0043-"
NEW8="tess2021284114741-s0044-"
OLD9="-0214-s"
NEW9="-0215-s"

for name in `ls *py`; do
  echo $name
  cp -f ${name} temp.foo
  sed -e "s/${OLD_NAME}/${NEW_NAME}/g; s/${OLD1}/${NEW1}/g; s/${OLD2}/${NEW2}/g; s/${OLD3}/${NEW3}/g; s/${OLD4}/${NEW4}/g; s/${OLD5}/${NEW5}/g; s/${OLD6}/${NEW6}/g; s/${OLD7}/${NEW7}/g;  s/${OLD8}/${NEW8}/g;  s/${OLD9}/${NEW9}/g" temp.foo > ${name}
done

