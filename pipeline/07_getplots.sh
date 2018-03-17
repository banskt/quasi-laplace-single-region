#!/bin/bash

CONFIGFILE=$1

if [ -z ${CONFIGFILE} ] || [ ! -f ${CONFIGFILE} ]; then 
    echo "Fatal! No configuration file found.";
    echo "Use this script as: ./06_getplots.sh CONFIGFILE"
    exit 1
fi

source ${CONFIGFILE}
source PATHS

SUBDIR="${CURDIR}/utils"
PLOT_JOBSUBDIR="${JOBSUBDIR}/getplots"
CREDIBLE="true false"
WHICH=""
if [ "${bBloreMeta}" = "true" ]; then WHICH+="blore "; fi
if [ "${bFinemap}"   = "true" ]; then WHICH+="finemap "; fi
if [ "${bPimass}"    = "true" ]; then WHICH+="${MODEL_PIMASS}"; fi

OUTDIR="${POSTPROBDIR}/plots"

if [ ! -d ${PLOT_JOBSUBDIR} ]; then mkdir -p ${PLOT_JOBSUBDIR}; fi
cd ${PLOT_JOBSUBDIR}

for CRED in ${CREDIBLE}; do
    CREDFLAG=""
    THISH2=`echo ${HERITABILITY} | tr -d .`
    THISL=`echo ${CASE_CONTROL_RATIO} | tr -d . | awk '{printf "%-3s", $1}' | tr ' ' '0'`
    OUTPREFIX="pip_prc_${PHENO_SIM_TYPE}_h${THISH2}_c${NCAUSAL}_l${THISL}"
    if [ "${USE_AGESEX}" = "true" ]; then
        OUTPREFIX="${OUTPREFIX}_cov"
    fi
    if [ "${CRED}" = "true" ]; then 
        CREDFLAG="--credible"
        OUTPREFIX="${OUTPREFIX}_cred"
    fi
    
    PLOT_JOBNAME="plot_${CONFIGFILE}_c${CRED}"
    sed "s|_JOBNAME|${PLOT_JOBNAME}|g;
         s|_START__|${START}|g;
         s|_E_N_D__|${END}|g;
         s|__OUTF__|${OUTDIR}/${OUTPREFIX}.pdf|g;
         s|__CRED__|\"${CREDFLAG}\"|g;
         s|_WHICH__|\"${WHICH}\"|g;
         s|_SIMDIR_|${SIMDIR}|g;
         s|_LOCUSF_|${USELOCI}|g;
         s|__CMAX__|${NCAUSAL}|g;
         s|_WORKDIR|${POSTPROBDIR}|g;
         s|_PLT_PIP|${PLOTPIP}|g;
        " ${MASTER_BSUBDIR}/getplots.bsub > ${PLOT_JOBNAME}.bsub
    bsub < ${PLOT_JOBNAME}.bsub
done

cd ${CURDIR}
