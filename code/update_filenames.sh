#!/bin/bash

OLD_NAME="sector25_20200717"
NEW_NAME="sector26_20200730"

OLD1="SECTOR = 25"
NEW1="SECTOR = 26"

OLD2="SECTOR1 = 25"
NEW2="SECTOR1 = 26"

OLD3="SECTOR2 = 25"
NEW3="SECTOR2 = 26"

# This replaces my local directory
OLD4="sector25"
NEW4="sector26"

# This replaces SPOC data directory
OLD5="sector-025"
NEW5="sector-026"


for name in `ls *py`; do
  echo $name
  cp -f ${name} temp.foo
  sed -e "s/${OLD_NAME}/${NEW_NAME}/g; s/${OLD1}/${NEW1}/g; s/${OLD2}/${NEW2}/g; s/${OLD3}/${NEW3}/g; s/${OLD4}/${NEW4}/g; s/${OLD5}/${NEW5}/g" temp.foo > ${name}
done

