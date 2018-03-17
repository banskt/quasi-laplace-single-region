#!/bin/bash

source PATHS

RANDSTRING=`cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 4 | head -n 1`
THIS_SIMDIR="${BASEDIR}/metaanalysis"
SNPTEST_JOBSUBDIR="${JOBSUBDIR}/metaanalysis"

for STUDY in ${STUDYNAMES[@]}; do
    SAMPLEDIR="${THIS_SIMDIR}/samples/${STUDY}"
    if [ ! -d ${SAMPLEDIR} ]; then mkdir -p ${SAMPLEDIR}; fi
    cp ${DOSAGEDIR}/${STUDY}/${STUDY}${SAMPLEPREFIX}.sample ${SAMPLEDIR}/phenotypes.sample
done

if [ ! -d ${SNPTEST_JOBSUBDIR} ]; then mkdir -p ${SNPTEST_JOBSUBDIR}; fi
cd ${SNPTEST_JOBSUBDIR}

SNPTEST_JOBNAME="snptest_${RANDSTRING}"
for STUDY in ${STUDYNAMES[@]}; do
    JOBNAME="${SNPTEST_JOBNAME}_${STUDY}"
    sed "s|_JOBNAME|${JOBNAME}|g;
         s|_SIMDIR_|${THIS_SIMDIR}|g;
         s|_GSTUDY_|${STUDY}|g;
         s|_SNPTEST|${SNPTEST}|g;
         s|_LOCIDIR|${DOSAGEDIR}|g;
         s|_USELOCI|${BIGLOCUS}|g;
        " ${MASTER_BSUBDIR}/snptest.bsub > ${JOBNAME}.bsub
    bsub < ${JOBNAME}.bsub
done

META_JOBNAME="meta_${RANDSTRING}"
sed "s|_JOBNAME|${META_JOBNAME}|g;
     s|_SIMDIR_|${THIS_SIMDIR}|g;
     s|_STUDYN_|\"${STUDYNAMES[*]}\"|g;
     s|_SAMPLES|\"${STUDYSAMPLES[*]}\"|g;
     s|_LOCUSN_|${BIGLOCUS}|g;
     s|_SCRIPT_|${GENINF}|g;
     s|__META__|${META}|g;
    " ${MASTER_BSUBDIR}/meta.bsub > ${META_JOBNAME}.bsub
bsub -w "done(${SNPTEST_JOBNAME}*)" < ${META_JOBNAME}.bsub

cd ${CURDIR}
