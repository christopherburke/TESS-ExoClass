#!/bin/bash

OLD_NAME="sector15_20190927"
NEW_NAME="sector16_20191029"

OLD1="SECTOR = 15"
NEW1="SECTOR = 16"

OLD2="SECTOR1 = 15"
NEW2="SECTOR1 = 16"

OLD3="SECTOR2 = 15"
NEW3="SECTOR2 = 16"

# This replaces my local directory
OLD4="sector15"
NEW4="sector16"

# This replaces SPOC data directory
OLD5="sector-15"
NEW5="sector-16"


for name in `ls *py`; do
  echo $name
  cp -f ${name} temp.foo
  sed -e "s/${OLD_NAME}/${NEW_NAME}/g; s/${OLD1}/${NEW1}/g; s/${OLD2}/${NEW2}/g; s/${OLD3}/${NEW3}/g; s/${OLD4}/${NEW4}/g; s/${OLD5}/${NEW5}/g" temp.foo > ${name}
done

