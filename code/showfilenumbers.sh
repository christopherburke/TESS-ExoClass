#!/bin/bash

# convert first argument from sector number to
# 3 characters with leading zeros
printf -v SECTORN "%03d" $1
echo "Sector ${SECTORN}"
DATA_DIR="/pdo/spoc-data/sector-${SECTORN}"

ls ${DATA_DIR}/dv-reports | head -4
ls ${DATA_DIR}/dv-reports/*dvm.pdf | head -2
ls ${DATA_DIR}/light-curve | head -4
ls ${DATA_DIR}/target-pixel | head -4
#ls ${DATA_DIR}/mini-reports | head -4

