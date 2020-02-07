#!/bin/bash

OLD_NAME="sector18_20191227"
NEW_NAME="sector19_20200124"

OLD1="SECTOR = 18"
NEW1="SECTOR = 19"

OLD2="SECTOR1 = 18"
NEW2="SECTOR1 = 19"

OLD3="SECTOR2 = 18"
NEW3="SECTOR2 = 19"

# This replaces my local directory
OLD4="sector18"
NEW4="sector19"

# This replaces SPOC data directory
OLD5="sector-18"
NEW5="sector-19"


for name in `ls *py`; do
  echo $name
  cp -f ${name} temp.foo
  sed -e "s/${OLD_NAME}/${NEW_NAME}/g; s/${OLD1}/${NEW1}/g; s/${OLD2}/${NEW2}/g; s/${OLD3}/${NEW3}/g; s/${OLD4}/${NEW4}/g; s/${OLD5}/${NEW5}/g" temp.foo > ${name}
done

