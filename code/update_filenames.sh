#!/bin/bash

OLD_NAME="sector11_20190615"
NEW_NAME="sector12_20190712"

OLD1="SECTOR = 11"
NEW1="SECTOR = 12"

OLD2="SECTOR1 = 11"
NEW2="SECTOR1 = 12"

OLD3="SECTOR2 = 11"
NEW3="SECTOR2 = 12"

OLD4="sector11"
NEW4="sector12"

OLD5="sector-11"
NEW5="sector-12"


for name in `ls *py`; do
  echo $name
  cp -f ${name} temp.foo
  sed -e "s/${OLD_NAME}/${NEW_NAME}/g; s/${OLD1}/${NEW1}/g; s/${OLD2}/${NEW2}/g; s/${OLD3}/${NEW3}/g; s/${OLD4}/${NEW4}/g; s/${OLD5}/${NEW5}/g" temp.foo > ${name}
done

