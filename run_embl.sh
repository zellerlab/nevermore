#!/bin/bash
#SBATCH -t 08:00:00
#SBATCH -e NM_%a_%j_err.txt
#SBATCH -o NM_%a_%j_out.txt
#SBATCH --qos=normal
#SBATCH --mem 8G
#SBATCH --no-requeue                                                                                                                                                                                                                                                                                                              
#SBATCH --job-name="nevermore"

# set bash strict mode http://redsymbol.net/articles/unofficial-bash-strict-mode/
set -euo pipefail
IFS=$'\n\t'

module load Nextflow/22.10.6

TIMESTAMP=$(date "+%Y%m%d_%H%M%S")

# later when the pipeline is more stable, we can switch to using a remote version:
#WORKFLOW=zellerlab/nevermore
#WORKFLOW_VERSION=test1

# or rely on your own local nevermore repo
WORKFLOW="/g/scb/zeller/$USER/software/nevermore"

# example workdir
WORKDIR=/scratch/$USER/WORK
mkdir -p ${WORKDIR}
ln -sf $WORKDIR work

# optionally you can pull specific versions with nextflow like this:
# nextflow pull ${WORKFLOW} -r ${WORKFLOW_VERSION}

nextflow run ${WORKFLOW} \
    -c /g/scb/zeller/$USER/software/nevermore/config/run_embl.config \
    -params-file /g/scb/zeller/$USER/software/nevermore/config/params_embl.yml \
    -work-dir ${WORKDIR} \
    -with-trace trace.${TIMESTAMP}.txt \
    -with-report report.${TIMESTAMP}.html \
    -resume

if [[ "$?" == "0" ]]; then
    status="finished"
else
    status="failed"
fi

echo Pipeline ${status} for ${SLURM_JOB_NAME} in ${PWD}.

# or send yourself an email
# echo Pipeline ${status} for ${SLURM_JOB_NAME} in ${PWD}. | mail -s "${WORKFLOW}::${SLURM_JOB_NAME} ${status}" $USERl@embl.de