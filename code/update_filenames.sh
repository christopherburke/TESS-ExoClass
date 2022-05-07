#!/bin/bash

OLD_NAME="sector-49_20220424"
NEW_NAME="sector-50_20220506"

OLD1="SECTOR = 49"
NEW1="SECTOR = 50"

OLD2="SECTOR1 = 49"
NEW2="SECTOR1 = 50"

OLD3="SECTOR2 = 49"
NEW3="SECTOR2 = 50"

# This replaces my local directory
OLD4="sector49"
NEW4="sector50"

# This replaces SPOC data directory
OLD5="sector-049"
NEW5="sector-050"

# DV report prefix
OLD6="tess2022057231053-"
NEW6="tess2022085182052-"
# DV report postfix
OLD7="-00585"
NEW7="-00597"

#LC prefix
OLD8="tess2022057073128-s0049-"
NEW8="tess2022085151738-s0050-"
OLD9="-0221-s"
NEW9="-0222-s"

# TOI Federation file
OLD10="FIXED-20220424"
NEW10="FIXED-20220506"

for name in `ls *py`; do
  echo $name
  cp -f ${name} temp.foo
  sed -e "s/${OLD_NAME}/${NEW_NAME}/g; s/${OLD1}/${NEW1}/g; s/${OLD2}/${NEW2}/g; s/${OLD3}/${NEW3}/g; s/${OLD4}/${NEW4}/g; s/${OLD5}/${NEW5}/g; s/${OLD6}/${NEW6}/g; s/${OLD7}/${NEW7}/g;  s/${OLD8}/${NEW8}/g;  s/${OLD9}/${NEW9}/g; s/${OLD10}/${NEW10}/g" temp.foo > ${name}
done

echo Retrieving TOI catalog file
curl https://tev.mit.edu/data/collection/193/csv/6/ > csv-file-toi-catalog.csv
sed -e 's/\"\"//g' -e 's/,\"[^\"]*/,\"NOCOMMENT/g' csv-file-toi-catalog.csv  > csv-file-toi-catalog-${NEW10}.csv
