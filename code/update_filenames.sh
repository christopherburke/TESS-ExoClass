#!/bin/bash

OLD_NAME="sector16_20191029"
NEW_NAME="sector14-16_20191104"

OLD1="SECTOR = 16"
NEW1="SECTOR = -1"

OLD2="SECTOR1 = 16"
NEW2="SECTOR1 = 14"

OLD3="SECTOR2 = 16"
NEW3="SECTOR2 = 16"

# This replaces my local directory
OLD4="sector16"
NEW4="sector14-16"

# This replaces SPOC data directory
OLD5="sector-16"
NEW5="sector-14-16"


for name in `ls *py`; do
  echo $name
  cp -f ${name} temp.foo
  sed -e "s/${OLD_NAME}/${NEW_NAME}/g; s/${OLD1}/${NEW1}/g; s/${OLD2}/${NEW2}/g; s/${OLD3}/${NEW3}/g; s/${OLD4}/${NEW4}/g; s/${OLD5}/${NEW5}/g" temp.foo > ${name}
done

