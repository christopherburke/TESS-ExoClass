#!/bin/bash

OLD_NAME="sector8_20190405"
NEW_NAME="sector1-6_20190428"

OLD1="SECTOR = 8"
NEW1="SECTOR = 6"

OLD2="SECTOR1 = 8"
NEW2="SECTOR1 = 1"

OLD3="SECTOR2 = 8"
NEW3="SECTOR2 = 6"

OLD4="sector8"
NEW4="sector1-6"

OLD5="sector-08"
NEW5="sector-01-06"


for name in `ls *py`; do
  echo $name
  cp -f ${name} temp.foo
  sed -e "s/${OLD_NAME}/${NEW_NAME}/g; s/${OLD1}/${NEW1}/g; s/${OLD2}/${NEW2}/g; s/${OLD3}/${NEW3}/g; s/${OLD4}/${NEW4}/g; s/${OLD5}/${NEW5}/g" temp.foo > ${name}
done

