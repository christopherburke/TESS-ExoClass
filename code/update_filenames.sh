#!/bin/bash

OLD_NAME="sector12_20190812"
NEW_NAME="sector1-13_20190812"

OLD1="SECTOR = 13"
NEW1="SECTOR = 13"

OLD2="SECTOR1 = 13"
NEW2="SECTOR1 = 1"

OLD3="SECTOR2 = 13"
NEW3="SECTOR2 = 13"

# This replaces my local directory
OLD4="sector13"
NEW4="sector1-13"

# This replaces SPOC data directory
OLD5="sector-13"
NEW5="sector-01-13"


for name in `ls *py`; do
  echo $name
  cp -f ${name} temp.foo
  sed -e "s/${OLD_NAME}/${NEW_NAME}/g; s/${OLD1}/${NEW1}/g; s/${OLD2}/${NEW2}/g; s/${OLD3}/${NEW3}/g; s/${OLD4}/${NEW4}/g; s/${OLD5}/${NEW5}/g" temp.foo > ${name}
done

