#!/bin/bash

OLD_NAME="sector21_20200309"
NEW_NAME="sector22_20200403"

OLD1="SECTOR = 21"
NEW1="SECTOR = 22"

OLD2="SECTOR1 = 21"
NEW2="SECTOR1 = 22"

OLD3="SECTOR2 = 21"
NEW3="SECTOR2 = 22"

# This replaces my local directory
OLD4="sector21"
NEW4="sector22"

# This replaces SPOC data directory
OLD5="sector-21"
NEW5="sector-22"


for name in `ls *py`; do
  echo $name
  cp -f ${name} temp.foo
  sed -e "s/${OLD_NAME}/${NEW_NAME}/g; s/${OLD1}/${NEW1}/g; s/${OLD2}/${NEW2}/g; s/${OLD3}/${NEW3}/g; s/${OLD4}/${NEW4}/g; s/${OLD5}/${NEW5}/g" temp.foo > ${name}
done

