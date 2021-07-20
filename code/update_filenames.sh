#!/bin/bash

OLD_NAME="sector1-36_20210615"
NEW_NAME="sector38_20210719"

OLD1="SECTOR = -1"
NEW1="SECTOR = 38"

OLD2="SECTOR1 = 1"
NEW2="SECTOR1 = 38"

OLD3="SECTOR2 = 36"
NEW3="SECTOR2 = 38"

# This replaces my local directory
OLD4="sector1-36"
NEW4="sector38"

# This replaces SPOC data directory
OLD5="sector-001-036"
NEW5="sector-038"

# DV report prefix
OLD6="tess2018206190142-"
NEW6="tess2021119082105-"
# DV report postfix
OLD7="-00471"
NEW7="-00488"

#LC prefix
OLD8="tess2021065132309-s0036-"
NEW8="tess2021118034608-s0038-"
OLD9="-0207-s"
NEW9="-0209-s"

for name in `ls *py`; do
  echo $name
  cp -f ${name} temp.foo
  sed -e "s/${OLD_NAME}/${NEW_NAME}/g; s/${OLD1}/${NEW1}/g; s/${OLD2}/${NEW2}/g; s/${OLD3}/${NEW3}/g; s/${OLD4}/${NEW4}/g; s/${OLD5}/${NEW5}/g; s/${OLD6}/${NEW6}/g; s/${OLD7}/${NEW7}/g;  s/${OLD8}/${NEW8}/g;  s/${OLD9}/${NEW9}/g" temp.foo > ${name}
done

