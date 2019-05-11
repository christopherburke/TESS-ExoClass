#!/bin/bash

OLD_NAME="sector1-6_20190428"
NEW_NAME="sector9_20190505"

OLD1="SECTOR = -1"
NEW1="SECTOR = 9"

OLD2="SECTOR1 = 1"
NEW2="SECTOR1 = 9"

OLD3="SECTOR2 = 6"
NEW3="SECTOR2 = 9"

OLD4="sector1-6"
NEW4="sector9"

OLD5="sector-01-06"
NEW5="sector-09"


for name in `ls *py`; do
  echo $name
  cp -f ${name} temp.foo
  sed -e "s/${OLD_NAME}/${NEW_NAME}/g; s/${OLD1}/${NEW1}/g; s/${OLD2}/${NEW2}/g; s/${OLD3}/${NEW3}/g; s/${OLD4}/${NEW4}/g; s/${OLD5}/${NEW5}/g" temp.foo > ${name}
done

