#!/bin/bash

OLD_NAME="sector36_20210421"
NEW_NAME="sector37_20210614"

OLD1="SECTOR = 36"
NEW1="SECTOR = 37"

OLD2="SECTOR1 = 36"
NEW2="SECTOR1 = 37"

OLD3="SECTOR2 = 36"
NEW3="SECTOR2 = 37"

# This replaces my local directory
OLD4="sector36"
NEW4="sector37"

# This replaces SPOC data directory
OLD5="sector-036"
NEW5="sector-037"

# DV report prefix
OLD6="tess2021066093107-"
NEW6="tess2021092173506-"
# DV report postfix
OLD7="-00460"
NEW7="-00478"

#LC prefix
OLD8="tess2021065132309-s0036-"
NEW8="tess2021091135823-s0037-"
OLD9="-0207-s"
NEW9="-0208-s"

for name in `ls *py`; do
  echo $name
  cp -f ${name} temp.foo
  sed -e "s/${OLD_NAME}/${NEW_NAME}/g; s/${OLD1}/${NEW1}/g; s/${OLD2}/${NEW2}/g; s/${OLD3}/${NEW3}/g; s/${OLD4}/${NEW4}/g; s/${OLD5}/${NEW5}/g; s/${OLD6}/${NEW6}/g; s/${OLD7}/${NEW7}/g;  s/${OLD8}/${NEW8}/g;  s/${OLD9}/${NEW9}/g" temp.foo > ${name}
done

