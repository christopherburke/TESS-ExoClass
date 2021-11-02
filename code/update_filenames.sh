#!/bin/bash

OLD_NAME="sector42_20210930"
NEW_NAME="sector43_20211101"

OLD1="SECTOR = 42"
NEW1="SECTOR = 43"

OLD2="SECTOR1 = 42"
NEW2="SECTOR1 = 43"

OLD3="SECTOR2 = 42"
NEW3="SECTOR2 = 43"

# This replaces my local directory
OLD4="sector42"
NEW4="sector43"

# This replaces SPOC data directory
OLD5="sector-042"
NEW5="sector-043"

# DV report prefix
OLD6="tess2021233042500-"
NEW6="tess2021259155059-"
# DV report postfix
OLD7="-00517"
NEW7="-00522"

#LC prefix
OLD8="tess2021232031932-s0042-"
NEW8="tess2021258175143-s0043-"
OLD9="-0213-s"
NEW9="-0214-s"

for name in `ls *py`; do
  echo $name
  cp -f ${name} temp.foo
  sed -e "s/${OLD_NAME}/${NEW_NAME}/g; s/${OLD1}/${NEW1}/g; s/${OLD2}/${NEW2}/g; s/${OLD3}/${NEW3}/g; s/${OLD4}/${NEW4}/g; s/${OLD5}/${NEW5}/g; s/${OLD6}/${NEW6}/g; s/${OLD7}/${NEW7}/g;  s/${OLD8}/${NEW8}/g;  s/${OLD9}/${NEW9}/g" temp.foo > ${name}
done

