#!/bin/bash

source PATHS

THIS_JOBSUBDIR="${JOBSUBDIR}/create_loci"
if [ ! -d ${THIS_JOBSUBDIR} ]; then mkdir -p ${THIS_JOBSUBDIR}; fi
if [ ! -d ${DOSAGEDIR}      ]; then mkdir -p ${DOSAGEDIR}; fi
cp ${REF_LOCUSNAMES} ${LOCUSNAMES}

for STUDY in ${STUDYNAMES[@]}; do
    if [ ! -d ${DOSAGEDIR}/${STUDY} ]; then mkdir -p ${DOSAGEDIR}/${STUDY}; fi
    cp ${REF_DOSAGEDIR}/${STUDY}/*.sample ${DOSAGEDIR}/${STUDY}/
done

## ======= DO NOT CHANGE BELOW =======================

cd ${THIS_JOBSUBDIR}

while read LOCUSPREFIX; do
    JOBNAME="common_SNPs_${LOCUSPREFIX}"
    sed -e "s|_JOBNAME|${JOBNAME}|g;
            s|_SCRIPT_|${CREATELOCI}|g;
            s|_QCTOOL_|${QCTOOL}|g;
            s|_LOCUSP_|${LOCUSPREFIX}|g;
            s|_STUDYN_|\"${STUDYNAMES[*]}\"|g;
            s|_LOCIDO_|${REF_DOSAGEDIR}|g;
            s|_LOCIDN_|${DOSAGEDIR}|g;" ${MASTER_BSUBDIR}/create_loci.bsub > ${JOBNAME}.bsub
    bsub < ${JOBNAME}.bsub
done < ${LOCUSNAMES}

cd ${CURDIR}
