#!/bin/bash

source PATHS
BIGLOCUSNAME="BigLocus"
echo ${BIGLOCUSNAME} > ${BIGLOCUS}

for STUDY in ${STUDYNAMES[@]} "combined"; do
    OUTDIR="${DOSAGEDIR}/${STUDY}"
    for EXT in gen map bgen matgen matmap; do
        if [ -f ${OUTDIR}/${BIGLOCUSNAME}.${EXT} ]; then rm -rf ${OUTDIR}/${BIGLOCUSNAME}.${EXT}; fi
    done

    echo "Creating big locus genotype and map file for ${STUDY}"
    while read LOCUSPREFIX; do
        cat ${OUTDIR}/${LOCUSPREFIX}.gen >> ${OUTDIR}/${BIGLOCUSNAME}.gen
        cat ${OUTDIR}/${LOCUSPREFIX}.map >> ${OUTDIR}/${BIGLOCUSNAME}.map
    done < ${USELOCI}

    echo "Creating bgen file"
    ${QCTOOL} -g ${OUTDIR}/${BIGLOCUSNAME}.gen -og ${OUTDIR}/${BIGLOCUSNAME}.bgen
    
    echo "Creating map file in BIMBAM format for ${STUDY}"
    cat ${OUTDIR}/${BIGLOCUSNAME}.map | awk '{printf "%s, %s, 1\n", $2, $4}' > ${OUTDIR}/${BIGLOCUSNAME}.matmap
    echo "Creating dosage file in BIMBAM format for ${STUDY}"
    cat ${OUTDIR}/${BIGLOCUSNAME}.gen | awk '{printf "%s, %s, %s", $2, $5, $4; for(i=6; i<NF; i+=3) {printf ", "$(i+0)*0+$(i+1)*1+$(i+2)*2}; printf "\n" }' > ${OUTDIR}/${BIGLOCUSNAME}.matgen

done
