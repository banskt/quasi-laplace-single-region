#!/bin/bash

source PATHS
LOCUSNAME="BigLocus"
echo ${LOCUSNAME} > ${BIGLOCUS}

for STUDY in ${STUDYNAMES[@]} "combined"; do
    OUTDIR="${DOSAGEDIR}/${STUDY}"

    echo "Creating big locus genotype and map file for ${STUDY}"
    while read LOCUSPREFIX; do
        cat ${OUTDIR}/${LOCUSPREFIX}.gen >> ${OUTDIR}/${LOCUSNAME}.gen
        cat ${OUTDIR}/${LOCUSPREFIX}.map >> ${OUTDIR}/${LOCUSNAME}.map
    done < ${USELOCI}

    echo "Creating bgen file"
    ${QCTOOL} -g ${OUTDIR}/${LOCUSNAME}.gen -og ${OUTDIR}/${LOCUSNAME}.bgen
    
    echo "Creating map file in BIMBAM format for ${STUDY}"
    cat ${OUTDIR}/${LOCUSNAME}.map | awk '{printf "%s, %s, 1\n", $2, $4}' > ${OUTDIR}/${LOCUSNAME}.matmap
    echo "Creating dosage file in BIMBAM format for ${STUDY}"
    cat ${OUTDIR}/${LOCUSNAME}.gen | awk '{printf "%s, %s, %s", $2, $5, $4; for(i=6; i<NF; i+=3) {printf ", "$(i+0)*0+$(i+1)*1+$(i+2)*2}; printf "\n" }' > ${OUTDIR}/${LOCUSNAME}.matgen

done
