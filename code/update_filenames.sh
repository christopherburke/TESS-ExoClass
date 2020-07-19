#!/bin/bash

OLD_NAME="sector24_20200604"
NEW_NAME="sector25_20200717"

OLD1="SECTOR = 24"
NEW1="SECTOR = 25"

OLD2="SECTOR1 = 24"
NEW2="SECTOR1 = 25"

OLD3="SECTOR2 = 24"
NEW3="SECTOR2 = 25"

# This replaces my local directory
OLD4="sector24"
NEW4="sector25"

# This replaces SPOC data directory
OLD5="sector-024"
NEW5="sector-025"


for name in `ls *py`; do
  echo $name
  cp -f ${name} temp.foo
  sed -e "s/${OLD_NAME}/${NEW_NAME}/g; s/${OLD1}/${NEW1}/g; s/${OLD2}/${NEW2}/g; s/${OLD3}/${NEW3}/g; s/${OLD4}/${NEW4}/g; s/${OLD5}/${NEW5}/g" temp.foo > ${name}
done

