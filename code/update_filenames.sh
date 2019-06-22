#!/bin/bash

OLD_NAME="sector10_20190604"
NEW_NAME="sector11_20190615"

OLD1="SECTOR = 10"
NEW1="SECTOR = 11"

OLD2="SECTOR1 = 10"
NEW2="SECTOR1 = 11"

OLD3="SECTOR2 = 10"
NEW3="SECTOR2 = 11"

OLD4="sector10"
NEW4="sector11"

OLD5="sector-10"
NEW5="sector-11"


for name in `ls *py`; do
  echo $name
  cp -f ${name} temp.foo
  sed -e "s/${OLD_NAME}/${NEW_NAME}/g; s/${OLD1}/${NEW1}/g; s/${OLD2}/${NEW2}/g; s/${OLD3}/${NEW3}/g; s/${OLD4}/${NEW4}/g; s/${OLD5}/${NEW5}/g" temp.foo > ${name}
done

