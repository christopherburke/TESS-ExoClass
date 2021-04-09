#!/bin/bash

OLD_NAME="sector34_20210303"
NEW_NAME="sector35_20210407"

OLD1="SECTOR = 34"
NEW1="SECTOR = 35"

OLD2="SECTOR1 = 34"
NEW2="SECTOR1 = 35"

OLD3="SECTOR2 = 34"
NEW3="SECTOR2 = 35"

# This replaces my local directory
OLD4="sector34"
NEW4="sector35"

# This replaces SPOC data directory
OLD5="sector-034"
NEW5="sector-035"

# DV report prefix
OLD6="tess2021014055109-"
NEW6="tess2021040113508-"
# DV report postfix
OLD7="-00444"
NEW7="-00453"

#LC prefix
OLD8="tess2021014023720-s0034-"
NEW8="tess2021039152502-s0035-"
OLD9="-0204-s"
NEW9="-0205-s"

for name in `ls *py`; do
  echo $name
  cp -f ${name} temp.foo
  sed -e "s/${OLD_NAME}/${NEW_NAME}/g; s/${OLD1}/${NEW1}/g; s/${OLD2}/${NEW2}/g; s/${OLD3}/${NEW3}/g; s/${OLD4}/${NEW4}/g; s/${OLD5}/${NEW5}/g; s/${OLD6}/${NEW6}/g; s/${OLD7}/${NEW7}/g;  s/${OLD8}/${NEW8}/g;  s/${OLD9}/${NEW9}/g" temp.foo > ${name}
done

