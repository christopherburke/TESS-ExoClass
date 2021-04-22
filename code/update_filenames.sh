#!/bin/bash

OLD_NAME="sector35_20210407"
NEW_NAME="sector36_20210421"

OLD1="SECTOR = 35"
NEW1="SECTOR = 36"

OLD2="SECTOR1 = 35"
NEW2="SECTOR1 = 36"

OLD3="SECTOR2 = 35"
NEW3="SECTOR2 = 36"

# This replaces my local directory
OLD4="sector35"
NEW4="sector36"

# This replaces SPOC data directory
OLD5="sector-035"
NEW5="sector-036"

# DV report prefix
OLD6="tess2021040113508-"
NEW6="tess2021066093107-"
# DV report postfix
OLD7="-00453"
NEW7="-00460"

#LC prefix
OLD8="tess2021039152502-s0035-"
NEW8="tess2021065132309-s0036-"
OLD9="-0205-s"
NEW9="-0207-s"

for name in `ls *py`; do
  echo $name
  cp -f ${name} temp.foo
  sed -e "s/${OLD_NAME}/${NEW_NAME}/g; s/${OLD1}/${NEW1}/g; s/${OLD2}/${NEW2}/g; s/${OLD3}/${NEW3}/g; s/${OLD4}/${NEW4}/g; s/${OLD5}/${NEW5}/g; s/${OLD6}/${NEW6}/g; s/${OLD7}/${NEW7}/g;  s/${OLD8}/${NEW8}/g;  s/${OLD9}/${NEW9}/g" temp.foo > ${name}
done

