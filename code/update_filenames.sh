#!/bin/bash

OLD_NAME="sector22_20200403"
NEW_NAME="sector24_20200604"

OLD1="SECTOR = 22"
NEW1="SECTOR = 24"

OLD2="SECTOR1 = 22"
NEW2="SECTOR1 = 24"

OLD3="SECTOR2 = 22"
NEW3="SECTOR2 = 24"

# This replaces my local directory
OLD4="sector22"
NEW4="sector24"

# This replaces SPOC data directory
OLD5="sector-22"
NEW5="sector-024"


for name in `ls *py`; do
  echo $name
  cp -f ${name} temp.foo
  sed -e "s/${OLD_NAME}/${NEW_NAME}/g; s/${OLD1}/${NEW1}/g; s/${OLD2}/${NEW2}/g; s/${OLD3}/${NEW3}/g; s/${OLD4}/${NEW4}/g; s/${OLD5}/${NEW5}/g" temp.foo > ${name}
done

