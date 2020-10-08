#!/bin/bash

OLD_NAME="sector27_20200905"
NEW_NAME="sector28_20200925"

OLD1="SECTOR = 27"
NEW1="SECTOR = 28"

OLD2="SECTOR1 = 27"
NEW2="SECTOR1 = 28"

OLD3="SECTOR2 = 27"
NEW3="SECTOR2 = 28"

# This replaces my local directory
OLD4="sector27"
NEW4="sector28"

# This replaces SPOC data directory
OLD5="sector-027"
NEW5="sector-028"

# DV report prefix
OLD6="tess2020187183116-"
NEW6="tess2020213081515-"
# DV report postfix
OLD7="-00362"
NEW7="-00364"

#LC prefix
OLD8="tess2020186164531-s0027-"
NEW8="tess2020212050318-s0028-"
OLD9="-0189-s"
NEW9="-0190-s"

for name in `ls *py`; do
  echo $name
  cp -f ${name} temp.foo
  sed -e "s/${OLD_NAME}/${NEW_NAME}/g; s/${OLD1}/${NEW1}/g; s/${OLD2}/${NEW2}/g; s/${OLD3}/${NEW3}/g; s/${OLD4}/${NEW4}/g; s/${OLD5}/${NEW5}/g; s/${OLD6}/${NEW6}/g; s/${OLD7}/${NEW7}/g;  s/${OLD8}/${NEW8}/g;  s/${OLD9}/${NEW9}/g" temp.foo > ${name}
done

