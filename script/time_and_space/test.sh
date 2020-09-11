#!/usr/bin/bash
start=$(date +%s.%N)

echo "---START---" >> space.csv
echo "megahit" >> space.csv

for i in {1..5000000}; do
    x=2
done
echo "---END---" >> space.csv
# end=$(date +%s.%N)
# echo "${start} ${end}"
# runtime=$( echo "${end} - ${start}" | bc -l )
# echo $runtime