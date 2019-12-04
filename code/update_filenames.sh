#!/bin/bash

OLD_NAME="sector14-16_20191104"
NEW_NAME="sector17_20191127"

OLD1="SECTOR = -1"
NEW1="SECTOR = 17"

OLD2="SECTOR1 = 14"
NEW2="SECTOR1 = 17"

OLD3="SECTOR2 = 16"
NEW3="SECTOR2 = 17"

# This replaces my local directory
OLD4="sector14-16"
NEW4="sector17"

# This replaces SPOC data directory
OLD5="sector-14-16"
NEW5="sector-17"


for name in `ls *py`; do
  echo $name
  cp -f ${name} temp.foo
  sed -e "s/${OLD_NAME}/${NEW_NAME}/g; s/${OLD1}/${NEW1}/g; s/${OLD2}/${NEW2}/g; s/${OLD3}/${NEW3}/g; s/${OLD4}/${NEW4}/g; s/${OLD5}/${NEW5}/g" temp.foo > ${name}
done

