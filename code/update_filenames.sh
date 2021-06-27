#!/bin/bash

OLD_NAME="sector37_20210614"
NEW_NAME="sector1-36_20210615"

OLD1="SECTOR = 37"
NEW1="SECTOR = -1"

OLD2="SECTOR1 = 37"
NEW2="SECTOR1 = 1"

OLD3="SECTOR2 = 37"
NEW3="SECTOR2 = 36"

# This replaces my local directory
OLD4="sector37"
NEW4="sector1-36"

# This replaces SPOC data directory
OLD5="sector-037"
NEW5="sector-001-036"

# DV report prefix
OLD6="tess2021092173506-"
NEW6="tess2018206190142-"
# DV report postfix
OLD7="-00478"
NEW7="-00471"

#LC prefix
OLD8="tess2021091135823-s0037-"
NEW8="tess2021065132309-s0036-"
OLD9="-0208-s"
NEW9="-0207-s"

for name in `ls *py`; do
  echo $name
  cp -f ${name} temp.foo
  sed -e "s/${OLD_NAME}/${NEW_NAME}/g; s/${OLD1}/${NEW1}/g; s/${OLD2}/${NEW2}/g; s/${OLD3}/${NEW3}/g; s/${OLD4}/${NEW4}/g; s/${OLD5}/${NEW5}/g; s/${OLD6}/${NEW6}/g; s/${OLD7}/${NEW7}/g;  s/${OLD8}/${NEW8}/g;  s/${OLD9}/${NEW9}/g" temp.foo > ${name}
done

