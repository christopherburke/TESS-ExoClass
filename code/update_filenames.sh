#!/bin/bash

OLD_NAME="sector38_20210719"
NEW_NAME="sector39_20210720"

OLD1="SECTOR = 38"
NEW1="SECTOR = 39"

OLD2="SECTOR1 = 38"
NEW2="SECTOR1 = 39"

OLD3="SECTOR2 = 38"
NEW3="SECTOR2 = 39"

# This replaces my local directory
OLD4="sector38"
NEW4="sector39"

# This replaces SPOC data directory
OLD5="sector-038"
NEW5="sector-039"

# DV report prefix
OLD6="tess2021119082105-"
NEW6="tess2021147062104-"
# DV report postfix
OLD7="-00488"
NEW7="-00491"

#LC prefix
OLD8="tess2021118034608-s0038-"
NEW8="tess2021146024351-s0039-"
OLD9="-0209-s"
NEW9="-0210-s"

for name in `ls *py`; do
  echo $name
  cp -f ${name} temp.foo
  sed -e "s/${OLD_NAME}/${NEW_NAME}/g; s/${OLD1}/${NEW1}/g; s/${OLD2}/${NEW2}/g; s/${OLD3}/${NEW3}/g; s/${OLD4}/${NEW4}/g; s/${OLD5}/${NEW5}/g; s/${OLD6}/${NEW6}/g; s/${OLD7}/${NEW7}/g;  s/${OLD8}/${NEW8}/g;  s/${OLD9}/${NEW9}/g" temp.foo > ${name}
done

