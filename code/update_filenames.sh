#!/bin/bash

OLD_NAME="sector-61_20230302"
NEW_NAME="sector-62_20230404"

# If multisector use SECTOR = -1
OLD1="SECTOR = 61"
NEW1="SECTOR = 62"

OLD2="SECTOR1 = 61"
NEW2="SECTOR1 = 62"

OLD3="SECTOR2 = 61"
NEW3="SECTOR2 = 62"

# This replaces my local directory
OLD4="sector61"
NEW4="sector62"

# This replaces SPOC data directory
OLD5="sector-061"
NEW5="sector-062"

# DV report prefix
OLD6="tess2023018070040-"
NEW6="tess2023043222439-"
# DV report postfix
OLD7="-00699"
NEW7="-00706"

#LC prefix
OLD8="tess2023018032328-s0061-"
NEW8="tess2023043185947-s0062-"
OLD9="-0250-s"
NEW9="-0254-s"

# TOI Federation file
OLD10="FIXED-20230302"
NEW10="FIXED-20230404"

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
