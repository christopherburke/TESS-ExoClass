#!/bin/bash

OLD_NAME="sector1-13_20190812"
NEW_NAME="sector14_20190918"

OLD1="SECTOR = 13"
NEW1="SECTOR = 14"

OLD2="SECTOR1 = 1"
NEW2="SECTOR1 = 14"

OLD3="SECTOR2 = 13"
NEW3="SECTOR2 = 14"

# This replaces my local directory
OLD4="sector1-13"
NEW4="sector14"

# This replaces SPOC data directory
OLD5="sector-01-13"
NEW5="sector-14"


for name in `ls *py`; do
  echo $name
  cp -f ${name} temp.foo
  sed -e "s/${OLD_NAME}/${NEW_NAME}/g; s/${OLD1}/${NEW1}/g; s/${OLD2}/${NEW2}/g; s/${OLD3}/${NEW3}/g; s/${OLD4}/${NEW4}/g; s/${OLD5}/${NEW5}/g" temp.foo > ${name}
done

