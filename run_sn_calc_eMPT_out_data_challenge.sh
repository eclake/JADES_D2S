#!/bin/bash

#First run S/N calculations for each pointing - point source
mkdir /appch/data/ec296/JWST_data_challenge_1/pointing_2/R1000_lines/
python sn_calc_eMPT_based_parallel.py --eMPT-out /appch/data/ec296/JWST_data_challenge_1/pointing_2/pointing_summary.dat -c /appch/data/ec296/JWST_data_challenge_1/catalogue_IDs.fits -B -s /appch/data/ec296/JWST_data_challenge_1/catalogue_sizes.fits -o /appch/data/ec296/JWST_data_challenge_1/sn_R1000_lines_pointing_2.p -r /appch/data/ec296/JWST_data_challenge_1/pointing_2/R1000_lines_es/ --n-group 19 --n-int 2 --n-exp 3  --R1000-lines --n-proc 30 --ps-off
mkdir /appch/data/ec296/JWST_data_challenge_1/pointing_7/R1000_lines/
python sn_calc_eMPT_based_parallel.py --eMPT-out /appch/data/ec296/JWST_data_challenge_1/pointing_7/pointing_summary.dat -c /appch/data/ec296/JWST_data_challenge_1/catalogue_IDs.fits -B -s /appch/data/ec296/JWST_data_challenge_1/catalogue_sizes.fits -o /appch/data/ec296/JWST_data_challenge_1/sn_R1000_lines_pointing_7.p -r /appch/data/ec296/JWST_data_challenge_1/pointing_7/R1000_lines_es/ --n-group 19 --n-int 2 --n-exp 3  --R1000-lines --n-proc 30 --ps-off
mkdir /appch/data/ec296/JWST_data_challenge_1/pointing_8/R1000_lines/
python sn_calc_eMPT_based_parallel.py --eMPT-out /appch/data/ec296/JWST_data_challenge_1/pointing_8/pointing_summary.dat -c /appch/data/ec296/JWST_data_challenge_1/catalogue_IDs.fits -B -s /appch/data/ec296/JWST_data_challenge_1/catalogue_sizes.fits -o /appch/data/ec296/JWST_data_challenge_1/sn_R1000_lines_pointing_8.p -r /appch/data/ec296/JWST_data_challenge_1/pointing_8/R1000_lines/ --n-group 19 --n-int 2 --n-exp 3  --R1000-lines --n-proc 30 --ps-off

#Then combine them
python combine_files.py -f1 /appch/data/ec296/JWST_data_challenge_1/sn_R1000_lines_pointing_2.p -f2 /appch/data/ec296/JWST_data_challenge_1/sn_R1000_lines_pointing_7.p -f3 /appch/data/ec296/JWST_data_challenge_1/sn_R1000_lines_pointing_8.p -o /appch/data/ec296/JWST_data_challenge_1/sn_R1000_lines_combined.fits --R1000-lines -i /appch/data/ec296/JWST_data_challenge_1/catalogue_IDs.fits --1910

#now run S/N calculations for each pointing - extended source
mkdir /appch/data/ec296/JWST_data_challenge_1/pointing_2/R1000_lines_es/
python sn_calc_eMPT_based_parallel.py --eMPT-out /appch/data/ec296/JWST_data_challenge_1/pointing_2/pointing_summary.dat -c /appch/data/ec296/JWST_data_challenge_1/catalogue_IDs.fits -B -s /appch/data/ec296/JWST_data_challenge_1/catalogue_sizes.fits -o /appch/data/ec296/JWST_data_challenge_1/sn_R1000_lines_pointing_2_es.p -r /appch/data/ec296/JWST_data_challenge_1/pointing_2/R1000_lines_es/ --n-group 19 --n-int 2 --n-exp 3  --R1000-lines --n-proc 30 --es-off
mkdir /appch/data/ec296/JWST_data_challenge_1/pointing_7/R1000_lines_es/
python sn_calc_eMPT_based_parallel.py --eMPT-out /appch/data/ec296/JWST_data_challenge_1/pointing_7/pointing_summary.dat -c /appch/data/ec296/JWST_data_challenge_1/catalogue_IDs.fits -B -s /appch/data/ec296/JWST_data_challenge_1/catalogue_sizes.fits -o /appch/data/ec296/JWST_data_challenge_1/sn_R1000_lines_pointing_7_es.p -r /appch/data/ec296/JWST_data_challenge_1/pointing_7/R1000_lines_es/ --n-group 19 --n-int 2 --n-exp 3  --R1000-lines --n-proc 30 --es-off
mkdir /appch/data/ec296/JWST_data_challenge_1/pointing_8/R1000_lines_es/
python sn_calc_eMPT_based_parallel.py --eMPT-out /appch/data/ec296/JWST_data_challenge_1/pointing_8/pointing_summary.dat -c /appch/data/ec296/JWST_data_challenge_1/catalogue_IDs.fits -B -s /appch/data/ec296/JWST_data_challenge_1/catalogue_sizes.fits -o /appch/data/ec296/JWST_data_challenge_1/sn_R1000_lines_pointing_8_es.p -r /appch/data/ec296/JWST_data_challenge_1/pointing_8/R1000_lines/ --n-group 19 --n-int 2 --n-exp 3  --R1000-lines --n-proc 30 --es-off

#Then combine them
python combine_files.py -f1 /appch/data/ec296/JWST_data_challenge_1/sn_R1000_lines_pointing_2_es.p -f2 /appch/data/ec296/JWST_data_challenge_1/sn_R1000_lines_pointing_7_es.p -f3 /appch/data/ec296/JWST_data_challenge_1/sn_R1000_lines_pointing_8_es.p -o /appch/data/ec296/JWST_data_challenge_1/sn_R1000_lines_combined_es.fits --R1000-lines -i /appch/data/ec296/JWST_data_challenge_1/catalogue_IDs.fits --es --1910
