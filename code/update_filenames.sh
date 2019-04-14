#!/bin/bash

OLD_NAME="sector7_20190305"
NEW_NAME="sector8_20190405"

OLD1="SECTOR = 7"
NEW1="SECTOR = 8"

OLD2="SECTOR1 = 7"
NEW2="SECTOR1 = 8"

OLD3="SECTOR2 = 7"
NEW3="SECTOR2 = 8"

OLD4="sector7"
NEW4="sector8"

OLD5="sector-07"
NEW5="sector-08"


for name in `ls *py`; do
  echo $name
  cp -f ${name} temp.foo
  sed -e "s/${OLD_NAME}/${NEW_NAME}/g; s/${OLD1}/${NEW1}/g; s/${OLD2}/${NEW2}/g; s/${OLD3}/${NEW3}/g; s/${OLD4}/${NEW4}/g; s/${OLD5}/${NEW5}/g" temp.foo > ${name}
done

