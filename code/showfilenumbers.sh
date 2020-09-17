#!/bin/bash

DATA_DIR='/pdo/spoc-data/sector-027'

ls ${DATA_DIR}/dv-reports | head -4
ls ${DATA_DIR}/dv-reports/*dvm.pdf | head -2
ls ${DATA_DIR}/light-curve | head -4
ls ${DATA_DIR}/target-pixel | head -4
#ls ${DATA_DIR}/mini-reports | head -4

