#!/bin/sh 
#
# gtblatpipeline.sh
###########################################################################
#
#  Purpose:  
#
#  Usage:
#
#
#  Env Vars:
#
#      See the configuration file
#
#  Inputs:
#
#      - FASTA file of  gene trap sequences
#      - gtblatpipeline.config
#  Outputs:
#
#      - An archive file
#      - Log files defined by the environment variables ${LOG_PROC},
#        ${LOG_DIAG}, ${LOG_CUR} and ${LOG_VAL}
#      Output from blat
#      - all_blat_hits.psl (used by pslReps)
#      Output from pslReps
#      - best_blat_hits.psr
#      - best_blat_hits.psl (used by gene trap alo load and Bob Script)
#      Output from "Bobs Script" - originally designed for GBrowse, it creates
#	the singel hits file we need, see below.
#      - best_blat_hits_multi.gff
#      - best_blat_hits_multi_Gbrowse.gff
#      - best_blat_hits_notes.txt
#      - best_blat_hits_seqIDs.txt
#      - best_blat_hits_seqIDs_multi.txt
#      - best_blat_hits_seqIDs_single.txt
#      - best_blat_hits_single.gff
#      - best_blat_hits_single_Gbrowse.gff - this is best_blat_hits_single.gff
#		formatted for use by Gbrowse, we use it because it determines
#		the coordinate range for each sequence (see line specified by
#		'gene' in column three used by gene trap coordload and alo load)
#      - Exceptions written to standard error
#      - Configuration and initialization errors are written to a log file
#        for the shell script
#
#  Exit Codes:
#
#      0:  Successful completion
#      1:  Fatal error occurred
#      2:  Non-fatal error occurred
#
#  Assumes: 
#
#  Implementation:  This script is a wrapper over:
#	1) Fetch BLAT input files from RADAR (default)
#	2) Create single input file for BLAT from files fetched from RADAR
#	3) BLAT - run blat on a fasta file
#	4) pslReps - run on output from BLAT (all_blat_hits.psl)
#	5) Bobs Script - run on output from pslReps (best_blat_hits.psl)
#       6) Log all_blat_hits.psl in RADAR for use by gene trap alo load
#	7) Log best_blat_hits_single_Gbrowse.gff in RADAR for use by gene
#		trap alo load
#	8) Log best_blat_hits_single_Gbrowse.gff in RADAR for use by gene
#		trap coordinate load
#  Notes:  None
#
###########################################################################
#
#  Set up a log file for the shell script in case there is an error
#  during configuration and initialization.
#
cd `dirname $0`/..

SCRIPT_NAME=`basename $0`

LOG=`pwd`/gtblatpipeline.log
rm -f ${LOG}

#
#  Verify the argument(s) to the shell script.
#
if [ $# -ne 0 ]
then
    echo "Usage: $0" | tee -a ${LOG}
    exit 1
fi

#
#  Establish the configuration file name 
#
CONFIG_LOAD=`pwd`/gtblatpipeline.config

#
#  Make sure the configuration file is readable.
#
if [ ! -r ${CONFIG_LOAD} ]
then
    echo "Cannot read configuration file: ${CONFIG_LOAD}" | tee -a ${LOG}
    exit 1
fi

# source config file
echo "Sourcing Configuration"
. ${CONFIG_LOAD}

#
#  Source the DLA library functions.
#
if [ "${DLAJOBSTREAMFUNC}" != "" ]
then
    if [ -r ${DLAJOBSTREAMFUNC} ]
    then
        . ${DLAJOBSTREAMFUNC}
    else
        echo "Cannot source DLA functions script: ${DLAJOBSTREAMFUNC}"
        exit 1
    fi
else
    echo "Environment variable DLAJOBSTREAMFUNC has not been defined."
    exit 1
fi

##################################################################
##################################################################
#
# main
#
##################################################################
##################################################################

#
# createArchive including OUTPUTDIR, startLog, getConfigEnv, get job key
#
preload 

#
# rm all files/dirs from OUTPUTDIR
#
cleanDir ${OUTPUTDIR}

#
# Wait for the "GT Filter Done" flag to be set. Stop waiting if the number
# of retries expires or the abort flag is found.
#
date | tee -a ${LOG_DIAG}
echo 'Wait for the "GT Filter Done" flag to be set' | tee -a ${LOG_DIAG}

RETRY=${PROC_CTRL_RETRIES}
while [ ${RETRY} -gt 0 ]
do
    READY=`${PROC_CTRL_CMD_PROD}/getFlag ${NS_DATA_LOADS} ${FLAG_GTFILTER}`
    ABORT=`${PROC_CTRL_CMD_PROD}/getFlag ${NS_DATA_LOADS} ${FLAG_ABORT}`

    if [ ${READY} -eq 1 -o ${ABORT} -eq 1 ] 
    then
        break
    else
        sleep ${PROC_CTRL_WAIT_TIME}
    fi

    RETRY=`expr ${RETRY} - 1`
done

#
# Terminate the script if the number of retries expired or the abort flag
# was found.
#
if [ ${RETRY} -eq 0 ]
then
   echo "${SCRIPT_NAME} timed out" | tee -a ${LOG_DIAG}
   date | tee -a ${LOG_DIAG}
   exit 1
elif [ ${ABORT} -eq 1 ]
then
   echo "${SCRIPT_NAME} aborted by process controller" | tee -a ${LOG_DIAG}
   date | tee -a ${LOG_DIAG}
   exit 1
fi

#
# Don't clear "GT Filter Done" flag. saturdaytasks.csh is also waiting for it
#

#
# select the Gene Trap input files that are ready to be processed
#
if [ ${APP_RADAR_INPUT} = true ]
then
     echo "Getting Gene Trap input files" | tee -a ${LOG_PROC} ${LOG_DIAG}
     APP_INFILES=`${RADAR_DBUTILS}/bin/getFilesToProcess.csh \
        ${RADAR_DBSCHEMADIR} ${JOBSTREAM} ${FILETYPE} 0`
     STAT=$?
     checkStatus ${STAT} "getFilesToProcess.csh"
fi

#
#  Make sure there is at least one Gene Trap input file to process
#
if [ "${APP_INFILES}" = "" ]
then
    echo "There are no Gene Trap input files to process" | \
	tee -a ${LOG_PROC} ${LOG_DIAG}
    shutDown

    echo `date` | tee -a ${LOG_PROC} ${LOG_DIAG}
    echo 'Set process control flag: GT Blat Done' | tee -a ${LOG_DIAG}
    ${PROC_CTRL_CMD_PROD}/setFlag ${NS_DATA_LOADS} ${FLAG_GTBLAT} ${SCRIPT_NAME}

    exit 0
fi

#
# Blat requires one input file
#
${APP_CAT_METHOD} ${APP_INFILES} > ${WORK_FASTA_FILE}
STAT=$?
checkStatus ${STAT} "${APP_CAT_METHOD} ${APP_INFILES} > ${WORK_FASTA_FILE}"

echo `date` | tee -a ${LOG_PROC} ${LOG_DIAG}
echo "Running Gene Trap Blat Pipeline" | tee -a ${LOG_PROC} ${LOG_DIAG}
${GFCLIENT} ${GFHOST} ${GFPORT} ${GFROOT} ${GFARGS} ${WORK_FASTA_FILE} \
    ${ALL_HITS_FILE} >> ${LOG_DIAG}
STAT=$?
checkStatus  ${STAT} "${GFCLIENT} ${GFHOST} ${GFPORT} ${GFROOT} ${GFARGS} ${WORK_FASTA_FILE}"

echo `date` | tee -a ${LOG_PROC} ${LOG_DIAG}
echo "Running PSL Filter" | tee -a ${LOG_PROC} ${LOG_DIAG}
${PSLFILTER} ${PSLARGS} ${ALL_HITS_FILE} ${BEST_HITS_FILE} ${BEST_HITS_RPT} \
	>> ${LOG_DIAG}
STAT=$?
checkStatus  ${STAT} "${PSLFILTER} ${PSLARGS} ${ALL_HITS_FILE} ${BEST_HITS_FILE} ${BEST_HITS_RPT}"

echo `date` | tee -a ${LOG_PROC} ${LOG_DIAG}
echo "Running Get Single Best Hits" | tee -a ${LOG_PROC} ${LOG_DIAG}
${PERL} ${BOBSCRIPT} ${BEST_HITS_FILE} ${BOBLABELARG} >> ${LOG_DIAG}
STAT=$?
checkStatus  ${STAT} "${BOBSCRIPT} ${BEST_HITS_FILE} ${BOBLABELARG}"

# add best single hits to master file, truncating if configured to do so
if [ ${TRUNCATE_MASTER} = true ] 
then
    echo "Truncating best single hits master file"
    rm ${MASTER_GFF_FILE}
fi

echo `date` | tee -a ${LOG_PROC} ${LOG_DIAG}
echo "Appending best single hits to master file"
${APP_CAT_METHOD} ${BEST_HITS_SGL_FILE} >> ${MASTER_GFF_FILE}
STAT=$?
checkStatus  ${STAT} "Appending best single hits to master file"

#
# create a new, unique file name for output we will log in radar
#
if [ ! -f ${FILECOUNTER} ]
then
    echo "Creating new file counter" 
    fileCounter=1
    echo ${fileCounter} > ${FILECOUNTER}
else
    echo "Incrementing file counter" 
    fileCounter=`cat ${FILECOUNTER}`
    fileCounter=`expr ${fileCounter} + 1`
    echo ${fileCounter} > ${FILECOUNTER}
fi

#
# move the original file names to the new file names
#
mv ${BEST_HITS_FILE} ${BEST_HITS_FILE}.${fileCounter}
mv ${BEST_HITS_SGL_FILE} ${BEST_HITS_SGL_FILE}.${fileCounter}

#
# reset the variable names
#
BEST_HITS_FILE="${BEST_HITS_FILE}.${fileCounter}"
BEST_HITS_SGL_FILE="${BEST_HITS_SGL_FILE}.${fileCounter}"

#
# log files to process
#

#  create path to /data/downloads output directory and log in radar
BEST_HITS_BASENAME=`basename ${BEST_HITS_FILE}`
BEST_HITS_SGL_BASENAME=`basename ${BEST_HITS_SGL_FILE}`

echo "Logging file to process ${BEST_HITS_FILE} for file type \
    ${FILETYPE_GTLOAD_BEST}"
${RADAR_DBUTILS}/bin/logFileToProcess.csh ${RADAR_DBSCHEMADIR} \
    ${BEST_HITS_FILE} ${DOWNLOADDIR} ${FILETYPE_GTLOAD_BEST}
STAT=$?
checkStatus ${STAT} "logFileToProcess.csh ${BEST_HITS_FILE} \
    for file type ${FILETYPE_GTLOAD_BEST}"

echo "Logging file to process ${BEST_HITS_SGL_FILE} for file type \
    ${FILETYPE_GTLOAD_SGL}"
${RADAR_DBUTILS}/bin/logFileToProcess.csh ${RADAR_DBSCHEMADIR} \
    ${BEST_HITS_SGL_FILE} ${DOWNLOADDIR} ${FILETYPE_GTLOAD_SGL}
STAT=$?
checkStatus ${STAT} "logFileToProcess.csh ${BEST_HITS_SGL_FILE} \
    for file type ${FILETYPE_GTLOAD_SGL}"

echo "Logging file to process ${BEST_HITS_SGL_FILE} for file type \
    ${FILETYPE_GTLOAD_SGL}"
${RADAR_DBUTILS}/bin/logFileToProcess.csh ${RADAR_DBSCHEMADIR} \
    ${BEST_HITS_SGL_FILE} ${DOWNLOADDIR} ${FILETYPE_COORDLOAD}
STAT=$?
checkStatus ${STAT} "logFileToProcess.csh ${BEST_HITS_SGL_FILE} \
    for file type ${FILETYPE_COORDLOAD}"

#
# move the files from "output" directory to the downloads
# subdirectory (ex. "/data/downloads")
# so the downstream processes can find them
#
mv -f ${BEST_HITS_FILE} ${DOWNLOADDIR}
mv -f ${BEST_HITS_SGL_FILE} ${DOWNLOADDIR}

#
# log the processed files
#
if [ ${APP_RADAR_INPUT} = true ]
then
    echo "Logging processed files"
    for file in ${APP_INFILES}
    do
        echo ${file} ${FILETYPE}
        ${RADAR_DBUTILS}/bin/logProcessedFile.csh ${RADAR_DBSCHEMADIR} \
            ${JOBKEY} ${file} ${FILETYPE}
        STAT=$?
        checkStatus ${STAT} "logProcessedFile.csh ${file}"
    done
fi

echo `date` | tee -a ${LOG_PROC} ${LOG_DIAG}
echo 'Set process control flag: GT Blat Done' | tee -a ${LOG_DIAG}
${PROC_CTRL_CMD_PROD}/setFlag ${NS_DATA_LOADS} ${FLAG_GTBLAT} ${SCRIPT_NAME}

echo "Running postload" | tee -a ${LOG_DIAG}
shutDown

exit 0
