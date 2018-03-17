#!/bin/bash

source PATHS

cd ${DOSAGEDIR}/combined

RANDSTRING=`cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 4 | head -n 1`

while read LOCUSPREFIX; do
    wc -l ${LOCUSPREFIX}.map
done < ${LOCUSNAMES} | head -n 200 | sort -n | tail -n 25 | awk '{print $2}' | sort > tmp_${RANDSTRING}

while read LOCUSPREFIX; do
    echo ${LOCUSPREFIX%.map}
done < tmp_${RANDSTRING} > ${USELOCI}

rm -rf tmp_${RANDSTRING}

cd ${CURDIR}
