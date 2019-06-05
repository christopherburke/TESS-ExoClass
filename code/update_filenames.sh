#!/bin/bash

OLD_NAME="sector1-9_20190517"
NEW_NAME="sector10_20190604"

OLD1="SECTOR = -1"
NEW1="SECTOR = 10"

OLD2="SECTOR1 = 1"
NEW2="SECTOR1 = 10"

OLD3="SECTOR2 = 9"
NEW3="SECTOR2 = 10"

OLD4="sector1-9"
NEW4="sector10"

OLD5="sector-01-09"
NEW5="sector-10"


for name in `ls *py`; do
  echo $name
  cp -f ${name} temp.foo
  sed -e "s/${OLD_NAME}/${NEW_NAME}/g; s/${OLD1}/${NEW1}/g; s/${OLD2}/${NEW2}/g; s/${OLD3}/${NEW3}/g; s/${OLD4}/${NEW4}/g; s/${OLD5}/${NEW5}/g" temp.foo > ${name}
done

