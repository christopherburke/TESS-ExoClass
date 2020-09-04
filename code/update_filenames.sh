#!/bin/bash

OLD_NAME="sector26_20200730"
NEW_NAME="sector14-26_20200825"

OLD1="SECTOR = 26"
NEW1="SECTOR = -1"

OLD2="SECTOR1 = 26"
NEW2="SECTOR1 = 14"

OLD3="SECTOR2 = 26"
NEW3="SECTOR2 = 26"

# This replaces my local directory
OLD4="sector26"
NEW4="sector14-26"

# This replaces SPOC data directory
OLD5="sector-026"
NEW5="sector-014-026"

# DV report prefix
OLD6="tess2020161181517-"
NEW6="tess2019199201929-"
# DV report postfix
OLD7="-00350_dvr.pdf"
NEW7="-00353_dvr.pdf"

for name in `ls *py`; do
  echo $name
  cp -f ${name} temp.foo
  sed -e "s/${OLD_NAME}/${NEW_NAME}/g; s/${OLD1}/${NEW1}/g; s/${OLD2}/${NEW2}/g; s/${OLD3}/${NEW3}/g; s/${OLD4}/${NEW4}/g; s/${OLD5}/${NEW5}/g; s/${OLD6}/${NEW6}/g; s/${OLD7}/${NEW7}/g " temp.foo > ${name}
done

