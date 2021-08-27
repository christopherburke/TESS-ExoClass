#!/bin/bash

OLD_NAME="sector39_20210720"
NEW_NAME="sector40_20210826"

OLD1="SECTOR = 39"
NEW1="SECTOR = 40"

OLD2="SECTOR1 = 39"
NEW2="SECTOR1 = 40"

OLD3="SECTOR2 = 39"
NEW3="SECTOR2 = 40"

# This replaces my local directory
OLD4="sector39"
NEW4="sector40"

# This replaces SPOC data directory
OLD5="sector-039"
NEW5="sector-040"

# DV report prefix
OLD6="tess2021147062104-"
NEW6="tess2021176033103-"
# DV report postfix
OLD7="-00491"
NEW7="-00503"

#LC prefix
OLD8="tess2021146024351-s0039-"
NEW8="tess2021175071901-s0040-"
OLD9="-0210-s"
NEW9="-0211-s"

for name in `ls *py`; do
  echo $name
  cp -f ${name} temp.foo
  sed -e "s/${OLD_NAME}/${NEW_NAME}/g; s/${OLD1}/${NEW1}/g; s/${OLD2}/${NEW2}/g; s/${OLD3}/${NEW3}/g; s/${OLD4}/${NEW4}/g; s/${OLD5}/${NEW5}/g; s/${OLD6}/${NEW6}/g; s/${OLD7}/${NEW7}/g;  s/${OLD8}/${NEW8}/g;  s/${OLD9}/${NEW9}/g" temp.foo > ${name}
done

