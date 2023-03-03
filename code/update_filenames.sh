#!/bin/bash

OLD_NAME="sector-14-60_20230219"
NEW_NAME="sector-61_20230302"

# If multisector use SECTOR = -1
OLD1="SECTOR = -1"
NEW1="SECTOR = 61"

OLD2="SECTOR1 = 14"
NEW2="SECTOR1 = 61"

OLD3="SECTOR2 = 60"
NEW3="SECTOR2 = 61"

# This replaces my local directory
OLD4="sector14-60"
NEW4="sector61"

# This replaces SPOC data directory
OLD5="sector-014-026+040-060"
NEW5="sector-061"

# DV report prefix
OLD6="tess2019199201929-"
NEW6="tess2023018070040-"
# DV report postfix
OLD7="-00692"
NEW7="-00699"

#LC prefix
OLD8="tess2022357055054-s0060-"
NEW8="tess2023018032328-s0061-"
OLD9="-0249-s"
NEW9="-0250-s"

# TOI Federation file
OLD10="FIXED-20230219"
NEW10="FIXED-20230302"

for name in `ls *py`; do
  echo $name
  cp -f ${name} temp.foo
  sed -e "s/${OLD_NAME}/${NEW_NAME}/g; s/${OLD1}/${NEW1}/g; s/${OLD2}/${NEW2}/g; s/${OLD3}/${NEW3}/g; s/${OLD4}/${NEW4}/g; s/${OLD5}/${NEW5}/g; s/${OLD6}/${NEW6}/g; s/${OLD7}/${NEW7}/g;  s/${OLD8}/${NEW8}/g;  s/${OLD9}/${NEW9}/g; s/${OLD10}/${NEW10}/g" temp.foo > ${name}
done

# IF multisector then in dvts_bulk_resamp.py SECTOR_OVRRIDE = -1; otherwise
# SECTOR_OVRRIDE = None
# Get the sector numbers from the NEW2 and NEW3 variables
secarr=(${NEW2}) # this splits variable on space character into array
secstrt=${secarr[2]} # this is element of array
secarr=(${NEW3})
secend=${secarr[2]}
echo Sectors: $secstrt - $secend
if [ $secstrt -eq $secend ]
then
 echo Single Sector
 cp -f dvts_bulk_resamp.py temp.foo
 sed -e "s/SECTOR_OVRRIDE = -1/SECTOR_OVRRIDE = None/" temp.foo > dvts_bulk_resamp.py
else
 echo Multisector
 cp -f dvts_bulk_resamp.py temp.foo
 sed -e "s/SECTOR_OVRRIDE = None/SECTOR_OVRRIDE = -1/" temp.foo > dvts_bulk_resamp.py
fi

echo Retrieving TOI catalog file
curl https://tev.mit.edu/data/collection/193/csv/6/ > csv-file-toi-catalog.csv
sed -e 's/\"\"//g' -e 's/,\"[^\"]*/,\"NOCOMMENT/g' csv-file-toi-catalog.csv  > csv-file-toi-catalog-${NEW10}.csv
