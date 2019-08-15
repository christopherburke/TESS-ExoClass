#!/bin/bash

OLD_NAME="sector12_20190712"
NEW_NAME="sector13_20190812"

OLD1="SECTOR = 12"
NEW1="SECTOR = 13"

OLD2="SECTOR1 = 12"
NEW2="SECTOR1 = 13"

OLD3="SECTOR2 = 12"
NEW3="SECTOR2 = 13"

OLD4="sector12"
NEW4="sector13"

OLD5="sector-12"
NEW5="sector-13"


for name in `ls *py`; do
  echo $name
  cp -f ${name} temp.foo
  sed -e "s/${OLD_NAME}/${NEW_NAME}/g; s/${OLD1}/${NEW1}/g; s/${OLD2}/${NEW2}/g; s/${OLD3}/${NEW3}/g; s/${OLD4}/${NEW4}/g; s/${OLD5}/${NEW5}/g" temp.foo > ${name}
done

