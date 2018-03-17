#!/bin/bash

CONFIGFILE=$1

if [ -z ${CONFIGFILE} ] || [ ! -f ${CONFIGFILE} ]; then 
    echo "Fatal! No configuration file found.";
    echo "Use this script as: ./05_simulate.sh CONFIGFILE"
    exit 1
fi

source ${CONFIGFILE}
source PATHS

SUBDIR="${CURDIR}/utils"
JOBSUBDIR="${JOBSUBDIR}/${PHENO_SIM_TYPE}_${HERITABILITY}"

for (( SIM=$START; SIM<=$END; SIM++ )); do

    INDEX=`echo $SIM | awk '{printf "%03d", $1}'`
    SIMFOLDER="sim${INDEX}"
    RANDSTRING=`cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 4 | head -n 1`

    THIS_JOBSUBDIR="${JOBSUBDIR}/${SIMFOLDER}"
    THIS_SIMDIR="${SIMDIR}/${SIMFOLDER}"

    if [ "${bDeleteOld}" = "true" ]; then
        if [   -d ${THIS_JOBSUBDIR} ]; then rm -rf ${THIS_JOBSUBDIR}; fi
        if [   -d ${THIS_SIMDIR} ];    then rm -rf ${THIS_SIMDIR};    fi
    fi

    if [ ! -d ${THIS_JOBSUBDIR} ]; then mkdir -p ${THIS_JOBSUBDIR}; fi
    if [ ! -d ${THIS_SIMDIR} ];    then mkdir -p ${THIS_SIMDIR};    fi

    cd ${THIS_JOBSUBDIR}

    BLORE_JOBSUBDIR="blore"
    SNPTEST_JOBSUBDIR="snptest"
    FINEMAP_JOBSUBDIR="finemap"
    PHENO_JOBSUBDIR="makepheno"
    PIMASS_JOBSUBDIR="pimass"
    GEMMA_JOBSUBDIR="gemma"

    BLORE_RESDIR="blore"

    if [ "${bMakePheno}" = "true" ]; then source ${SUBDIR}/create_phenotype; fi
    if [ "${bBloreSumm}" = "true" ]; then source ${SUBDIR}/blore_summary; fi
    if [ "${bBloreMeta}" = "true" ]; then source ${SUBDIR}/blore_meta; fi
    if [ "${bSnptest}"   = "true" ]; then source ${SUBDIR}/snptest; fi
    if [ "${bMeta}"      = "true" ]; then source ${SUBDIR}/meta; fi
    if [ "${bFinemap}"   = "true" ]; then source ${SUBDIR}/finemap; fi
    if [ "${bRegCov}"    = "true" ]; then source ${SUBDIR}/regresscov; fi
    if [ "${bPimass}"    = "true" ]; then source ${SUBDIR}/pimass; fi
    if [ "${bGemma}"     = "true" ]; then source ${SUBDIR}/gemma; fi

    cd ${CURDIR}

done
