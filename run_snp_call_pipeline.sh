#! /bin/bash

# Usage: copy this script to the working directory.  And prepare a
# fastq_map.txt file in the working directory.  The details about the
# fastq_map.txt format see the comments in pipeline file:
# /local/data/Malaria/Projects/Takala-Harrison/Cambodia_Bing/snp_call_wdl/pipelines.wdl

prefix=/local/data/Malaria/Projects/Takala-Harrison/Cambodia_Bing/snp_call_wdl_git

# prepare folders and clean previous runs

if [ -d output ]; then rm -rf output; fi
mkdir output
mkdir output/coverage output/bam output/flagstat output/gvcf output/db output/vcf 

if [ -d "cromwell-executions" ]; then rm -rf cromwell-executions; fi
if [ -d "cromwell-workflow-logs" ]; then rm -rf cromwell-workflow-logs ; fi


# modify inputs.json
cat "$prefix"/input/inputs.json  | tr  '"' '|' \
    | awk -v FS='|' -v pwd=`pwd` -v OFS='|' '\
        $2 == "out_dir" {$4 = pwd "/output" };
        $2 == "mapping_and_snp_calling.fastq_map" {$4 = pwd "/fastq_map.txt"};
        {print }' | tr '|' '"' > inputs.json

# run pipeline (using taskset to limit cpu usage)

taskset -a -c 1-24\
	java \
	-Dthread.pool.size=4 \
	-Dconfig.file="$prefix"/local.conf \
	-jar "$prefix"/bin/cromwell-59.jar \
	run "$prefix"/pipeline.wdl \
	-i inputs.json
