#!/bin/bash

OLD_NAME="sector14_20190918"
NEW_NAME="sector15_20190927"

OLD1="SECTOR = 14"
NEW1="SECTOR = 15"

OLD2="SECTOR1 = 14"
NEW2="SECTOR1 = 15"

OLD3="SECTOR2 = 14"
NEW3="SECTOR2 = 15"

# This replaces my local directory
OLD4="sector14"
NEW4="sector15"

# This replaces SPOC data directory
OLD5="sector-14"
NEW5="sector-15"


for name in `ls *py`; do
  echo $name
  cp -f ${name} temp.foo
  sed -e "s/${OLD_NAME}/${NEW_NAME}/g; s/${OLD1}/${NEW1}/g; s/${OLD2}/${NEW2}/g; s/${OLD3}/${NEW3}/g; s/${OLD4}/${NEW4}/g; s/${OLD5}/${NEW5}/g" temp.foo > ${name}
done

