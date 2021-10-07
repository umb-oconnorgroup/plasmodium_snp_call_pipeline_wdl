# using WDL version 1.0
version 1.0  

# =============== WORKFLOW ================================
# About inputs:
#
# 1.  Ref genome of parasite species and host species.  The reason why host
# species reference genome is included here is that MalariaGen pipeline
# suggests deplete raw reads mapped to the host species.  The raws reads are
# first mapped to the host genome.  The unmapped reads are then converted back to
# fastq files and then mapped to the parasite genome.
# 
# 2.  Mapping reads from a sample of other parasite species.  This is mainly for
# discover ancestral allele status.  In the path_arr array, the first Path
# struct has the ref genome of studied parasite and its host species.  The 2nd
# and onwards Path structs have the information for ref genome of the
# closely-related species and their host species.  If ancestral allele status
# is not needed, just specify the a single Path object in the Path_arr for the
# studied species and its host.
# 
# 3. About the fastq map file format
# Each line in the fastq_map file represents fastq files for a single sample.
#  1) The first column is the id for host species, id=0 is for the host species of the
#  studied parasite species being studies. id=1 ...2 are for the host spcies of the
#  closely related parasites, which are only used for ancestral allele inference.
#  2) The 2nd column is the unique sample name
#  3) From the3rd column onward, each two columns represent fastq files for a single
#  run. If it is single-end sequencing, the 2nd of each pair is set to be a place holder
#  'single', the first one of each pair is the actually fastq file; if it is paired-end
#  sequencing, both two of each pair are fastq files and should not contain the place
#  holder 'single'.
# 
#    Example 1: Single-run pair-end sequencing data
#
#		 0	s1	s1r1_1.fastq.gz	s1r1_2.fastq.gz
#
#    Example 2: Pair-end sequencing data with two runs
#
#		 0	s1	s1r1_1.fastq.gz	s1r1_2.fastq.gz	s1r2_1.fastq.gz	s1r2_2.fastq.gz
#
#    Example 3: Single-end sequencing data
#    
#    0	s1 s1r1.fastq.gq	single

workflow mapping_and_snp_calling {
	input {
		File fastq_map
		File chr_list
		Array[Path] path_arr 
		Int fastq2bam_nthreads
	}

	# check shared files and folders in each Path struct
	scatter (path in path_arr) {
		call check_path{input: path=path}
	}

 	# Sample level -- scatter
  	scatter (row in read_tsv(fastq_map)) {
  		Int num_pairs = (length(row) - 2)/2 
  		Int host_id = row[0]
  		String sample = row[1]
		Boolean path_ok = check_path.path_ok[host_id]
  
  		# run level -- scatter 
  		# for the case with multiple pairs of runs for a sample
  		scatter (i in range(num_pairs)) {
  			String f1 = row[(2*i+2)]
  			String f2 = row[(2*i+3)]
  			String run_name = sample + "_run_~{i}"
  			call fastq2bam { input: f1=f1, f2=f2, path=path_arr[host_id], 
			  	path_ok = path_ok,
			  	run_name=run_name, 
			  	nthread=fastq2bam_nthreads }
  		}
  		call bam2gvcf as parasite_bam2gvcf{ input: sample_i = sample, bam_files = fastq2bam.out_bam, path = path_arr[host_id] }
  
  		# Run this task to understand how raw reads are mapped to host genome.
  		# call bam2gvcf as host_bam2gvcf { input: sample_i = sample, bam_files = fastq2bam.host_bam, 
  		# 	path = path_arr[0], use_host = true }
  	}
  	# per chromosome joint calling	
  	     
  	call make_gvcf_map{ input: gvcfs = parasite_bam2gvcf.gvcf_1f,
  		outdir = path_arr[0].out_dir}
  
  	scatter (chrom in read_lines(chr_list)) {
  		call import_gvcf{ input: gvcf_map = make_gvcf_map.gvcf_map, chrom=chrom, path = path_arr[0]}
  		call genotype_gvcf{ input: db_dir=import_gvcf.db_dir, chrom=chrom, path = path_arr[0]}
		  
  		# VCF per chr filter
  		call filter_per_chr{input: vcf_chr=genotype_gvcf.chr_vcf, path=path_arr[0]}
  	}
  
	# Merge vcf and filter together
	call filter_genome_wide{input: filter1_vcfs=filter_per_chr.filter1_vcf, 
			path=path_arr[0]}

}

# =============== Definitions ================================

# Define a struct to simplify passing ref or exe files paths to tasks
struct Path {
		String java_exe 
		String samtools_exe
		String bedtools_exe
		String gatk_exe
		String bcftools_exe
		String bowtie2_exe
		String ref_fasta
		String ref_dict
		String ref_fai
		String ref_idx
		Array[String] known_snp_vcf
		Array[String] know_snp_idx
		String host_ref_fasta
		String host_ref_dict
		String host_ref_fai
		String host_ref_idx
		Array[String] host_known_snp_vcf
		Array[String] host_know_snp_idx
		String out_dir
}

# TASK: check folders and files
task check_path {
	input {Path path}
	command <<<
		check_files (){
			for file in "${@:1}";
			do
				if [ ! -f "$file" ]; then 
				echo ERROR: File "$file" NOT exist 1>&2;
				exit 1
				fi
			done
		}
		check_dirs (){
			for dir in "${@:1}";
			do
				if [ ! -d "$dir" ]; then 
				echo ERROR: Folder "$dir" NOT exist 1>&2;
				exit 1
				fi
			done
		}

		check_files \
			~{path.java_exe} \
			~{path.samtools_exe} \
			~{path.bedtools_exe} \
			~{path.gatk_exe} \
			~{path.bcftools_exe} \
			~{path.bowtie2_exe} \
			~{path.ref_fasta} \
			~{path.ref_dict} \
			~{path.ref_fai} \
			~{path.ref_idx} \
			~{sep=" " path.known_snp_vcf} \
			~{sep=" " path.know_snp_idx} \
			~{path.host_ref_fasta} \
			~{path.host_ref_dict} \
			~{path.host_ref_fai} \
			~{path.host_ref_idx} \
			~{sep=" " path.host_known_snp_vcf} \
			~{sep=" " path.host_know_snp_idx}

		check_dirs ~{path.out_dir}

	>>>
	output {
		Boolean path_ok = true
	}
}


# TASK: deplete reads aligned to host and map the unmapped reads for
# a pair of fastq files from a single run
# REF:
# https://sites.google.com/site/wiki4metagenomics/tools/short-read/remove-host-sequences
# http://broadinstitute.github.io/picard/explain-flags.html
task fastq2bam {
	input {
		String f1
		String f2
		Path path
		Boolean path_ok
		String run_name
		Int nthread = 1
	}
	String ref_prefix=sub(path.ref_fasta,".fasta", "")
	String host_ref_prefix=sub(path.host_ref_fasta,".fasta.*$", "")
	String parasite_raw_bam = path.out_dir + "/bam/" + run_name + '.bam'
	String host_raw_bam = path.out_dir + "/bam/host_" + run_name + '.bam'

	# if f2 == single, then f1 is the single end sequence reads
	# if f2 != single, then f1, f2 are the two paired-end read files.
	command <<<
		set -xeo pipefail

		if ! ~{path_ok}; then 
			echo CHECK_PATH Failed 1>&2;
			exit 1; 
		fi


		# Pair-end
		if [ ~{f2} != "single" ]; then
			if [ ! -f ~{f1} ]; then echo FILE NOT EXIST: ~{f1} 1>&2; exit 1; fi
			if [ ! -f ~{f2} ]; then echo FILE NOT EXIST: ~{f2} 1>&2; exit 1; fi

			# 1. Map to human genome
			~{path.bowtie2_exe} -x ~{host_ref_prefix} -1 ~{f1} -2 ~{f2} | \
				~{path.samtools_exe} view -q 0 -bS > host_raw.bam

			# 2. Keep (both reads) unmapped and sort by name
			# 20210426: change the filter flag from -f4 to -f 12 -F 256 to be more conserved
			# ref: https://sites.google.com/site/wiki4metagenomics/tools/short-read/remove-host-sequences
			~{path.samtools_exe} view -b -f 12 -F 256 host_raw.bam | \
				~{path.samtools_exe} sort -n -o unmapped_sorted.bam

			# 3. Convert unmapped to paired end fastq files.
			~{path.samtools_exe} fastq -@ ~{nthread} unmapped_sorted.bam \
				-1 f1.fastq.gz -2 f2.fastq.gz \
				-0 /dev/null -s /dev/null -n	

			# 4. Map unmapped reads to plasmodium genome
			# 	add read group id here. Sample name is set to run_name temporarily
			~{path.bowtie2_exe} -x ~{ref_prefix} -1 f1.fastq.gz -2 f2.fastq.gz \
				--rg-id ~{run_name} --rg SM:~{run_name} | \
				~{path.samtools_exe} view -q 0 -bS > raw.bam

			# 5. Clean up
			rm unmapped_sorted.bam f1.fastq.gz f2.fastq.gz

		# Unpaired
		else
			if [ ! -f ~{f1} ]; then echo FILE NOT EXIST: ~{f1} 1>&2; exit 1; fi
			# 1. Map to human genome
			~{path.bowtie2_exe} -x ~{host_ref_prefix} -U ~{f1} | \
				~{path.samtools_exe} view -q 0 -bS > host_raw.bam

			# 2. Filter out unmapped and sort by name
			~{path.samtools_exe} view -b -f 4 host_raw.bam | \
				~{path.samtools_exe} sort -n -o unmapped_sorted.bam

			# 3. Convert unmapped to single end to fastq files.
			~{path.samtools_exe} fastq -@ ~{nthread} unmapped_sorted.bam \
				-o /dev/null -s /dev/null \
				-0 f.fastq.gz -n  

			# 4. Map unmapped reads to plasmodium genome
			~{path.bowtie2_exe} -x ~{ref_prefix} -U f.fastq.gz \
				--rg-id ~{run_name} --rg SM:~{run_name} | \
				~{path.samtools_exe} view -q 0 -bS > raw.bam

			# 5. Clean up
			rm unmapped_sorted.bam f.fastq.gz

		fi

		mv raw.bam ~{parasite_raw_bam}
		mv host_raw.bam ~{host_raw_bam}
	>>>
	output {
		String out_bam = parasite_raw_bam
		String host_bam = host_raw_bam
	}
}

# TASK: recalibrate and run hyplotypeCaller
#
task bam2gvcf{
	input {
		String sample_i
		Array[String] bam_files
		Path path
		Boolean use_host = false
	}
	String barcode="unknown"

	# if mapped to host, use host reference and add prefix to sample name
	String ref_fasta = if use_host then path.host_ref_fasta else path.ref_fasta
	String sample = (if use_host then "host_" else "") + sample_i
	Array[String] known_snp_vcf = if use_host then path.host_known_snp_vcf else path.known_snp_vcf
	command <<<
		set -xeo pipefail
		~{path.gatk_exe} --java-options "-Xmx10G" \
		MergeSamFiles \
			-I ~{sep=" -I " bam_files} -O raw.bam
		
		~{path.gatk_exe} --java-options "-Xmx10G" \
		AddOrReplaceReadGroups \
			-I raw.bam -O fixed.bam \
			--SORT_ORDER coordinate --VALIDATION_STRINGENCY SILENT \
			--RGLB LIB --RGPL Illumina --RGPU ~{barcode} \
			--RGSM ~{sample_i} 
		rm raw.bam
		# Already added --RGID in bowtie2

		~{path.gatk_exe} --java-options "-Xmx10G" \
		SortSam \
			-I fixed.bam -O sorted.bam --SORT_ORDER coordinate
		rm fixed.bam
		
		# MarkDuplicates: PCR duplicates and Optical duplicates
		# Ref1: https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-
		# Ref2: https://www.cureffi.org/2012/12/11/how-pcr-duplicates-arise-in-next-generation-sequencing/
		~{path.gatk_exe} --java-options "-Xmx10G" \
		MarkDuplicates \
			--REMOVE_DUPLICATES \
			-I sorted.bam -O dedup.bam --METRICS_FILE dedup_metrics.txt
		rm sorted.bam dedup_metrics.txt 

		~{path.samtools_exe} index dedup.bam

		# Base quality score recalibration: use machine learning method to
		# correct some systematic error of base score generated by sequencing
		# machine
		#
		# Importance of ReadGroup: 
		# The recalibration system is read-group aware, meaning it uses @RG
		# tags to partition the data by read group.  This allows it to perform
		# the recalibration per read group, which reflects which library a read
		# belongs to and what lane it was sequenced in on the flowcell.
		#
		# Ref: https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR- 

		~{path.gatk_exe} --java-options "-Xmx10G" \
		BaseRecalibrator \
		-I dedup.bam -O recal_data.table \
		-R ~{ref_fasta} --known-sites ~{sep=" --known-sites " known_snp_vcf}

		~{path.gatk_exe} --java-options "-Xmx10G" \
		ApplyBQSR \
		-I dedup.bam -O recalibrated.bam \
		-R ~{ref_fasta} --bqsr-recal-file recal_data.table
		rm dedup.* recal_data.table

		# Call hyplotypeCaller:
		# Ref: https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller
		# Ref: https://www.youtube.com/watch?v=vocBk2MHp1A
		#
		# skip these if not running haplotyecaller for parasite 
		if ! ~{use_host}; then
			~{path.gatk_exe} --java-options "-Xmx10G" \
			HaplotypeCaller \
			-I recalibrated.bam -O ~{sample}.g.vcf \
			-R ~{ref_fasta} -ERC GVCF

			mv ~{sample}.g.vcf ~{path.out_dir}/gvcf/
			mv ~{sample}.g.vcf.idx ~{path.out_dir}/gvcf/
		fi

		~{path.bedtools_exe} genomecov -d -ibam recalibrated.bam \
			-g ~{ref_fasta} | gzip -v > recalibrated.coverage.gz
		mv recalibrated.coverage.gz ~{path.out_dir}/coverage/~{sample}_coverage.gz

		~{path.samtools_exe} flagstat recalibrated.bam >recalibrated.flagstat
		mv recalibrated.flagstat ~{path.out_dir}/flagstat/~{sample}_flagstat

		rm recalibrated.bam

	>>>
	output {
		String gvcf_1f= if use_host then "" else (path.out_dir + "/gvcf/" + sample + ".g.vcf") # abs path
	}
}

task make_gvcf_map {
	input {
		Array[String] gvcfs
		String outdir
	}
	command <<<
		set -xeo pipefail
		echo -e '~{sep="\n" gvcfs}'  > gvcf_list
		sed -e 's/^.*\///;s/.g.vcf$//' gvcf_list > sample_list
		paste sample_list gvcf_list > map.txt
		mv map.txt ~{outdir}/gvcf_map.txt
	>>>

	output {
		String gvcf_map = outdir + "/gvcf_map.txt"
	}	
}

# TASK: import all g.vcf to database for each chromosome

task import_gvcf{
	input {
		String gvcf_map
		String chrom
		Path path
	}
		String fname = sub(chrom, ":|-", "_") 
	
	command <<<
		set -xeo pipefail

		~{path.gatk_exe} --java-options "-Xmx20G" \
		GenomicsDBImport \
			--batch-size 500 --reader-threads 5 --consolidate \
			--sample-name-map ~{gvcf_map} --genomicsdb-workspace-path my_database_~{fname}  -L ~{chrom}

		mv my_database_~{fname} ~{path.out_dir}/db/
	>>>
	runtime {
		memory: "100G"
		cpu: 10
		sge_queue: "threaded.q"
	}
	output {
		String db_dir = path.out_dir + "/db/my_database_" + fname
	}	
}

# TASK: Joint genotype calling by chromosome using the database

task genotype_gvcf{
	input {
		String db_dir
		String chrom
		Path path
	}
	String fname = sub(chrom, ":|-", "_") 
	
	command <<<
		# 20210422: Added the `--genomicsdb-use-vcf-codec` option:  
		# BCFCodec performs slightly better but currently does not support 64-bit width 
		# positions and INFO fields and for computed annotation sizes to exceed 32-bit 
		# integer space. The new versions may use vcf codes by default but not for 4.1.7.

		# Without this option we can get the error:
		# java.lang.NullPointerException
		#         at htsjdk.variant.bcf2.BCF2Decoder.decodeInt(BCF2Decoder.java:226)
		# Ref:
		# 4.1.7: URL: https://gatk.broadinstitute.org/hc/en-us/articles/360042914991-GenotypeGVCFs
		# 4.2.0: URL: https://gatk.broadinstitute.org/hc/en-us/articles/360056970432-GenotypeGVCFs
		set -xeo pipefail
		ln -s ~{db_dir} my_database

		~{path.gatk_exe} --java-options "-Xmx20G" \
		GenotypeGVCFs \
			--genomicsdb-use-vcf-codec \
			-V gendb://my_database -O ~{fname}.vcf -R ~{path.ref_fasta}

		mv ~{fname}.vcf* ~{path.out_dir}/vcf/
	>>>
	runtime {
		memory: "100G"
		cpu: 10
		sge_queue: "threaded.q"
	}
	output {
		String chr_vcf = path.out_dir + "/vcf/" + fname + ".vcf"
	}	
}

task filter_per_chr {
	input {
		String vcf_chr
		Path path
	}
	command <<<
		# Only include SNPs (no indels)
		~{path.gatk_exe} --java-options "-Xmx24G" SelectVariants \
			-R ~{path.ref_fasta} \
			-V ~{vcf_chr} \
			--select-type-to-include SNP \
			-O 1.vcf

		# Hard filtering:
		# 	The GATK does not recommend use of compound filtering
		# expressions, e.g.  the logical || "OR".  For such expressions, if
		# a record is null for or missing a particular annotation in the
		# expression, the tool negates the entire compound expression and
		# so automatically passes the variant record even if it fails on
		# one of the expressions.
		# Resource: https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering
		~{path.gatk_exe} --java-options "-Xmx24G" VariantFiltration \
			-R ~{path.ref_fasta} \
			-V 1.vcf \
			-filter "QD<2.0" --filter-name "QD2" \
			-filter "FS>60.0" --filter-name "FS60" \
			-filter "MQ<40.0" --filter-name "MQ40" \
			-filter "MQRankSum<-12.5" --filter-name "MQRankSum-12.5" \
			-filter "ReadPosRankSum<-8.0" --filter-name "ReadPosRankSum-8" \
			-O 2.vcf

		# Keep SNPs that pass the above GATK filteration criteria
		~{path.bcftools_exe} view -f PASS 2.vcf \
			-Oz -o out.vcf.gz

		rm ./*.vcf
	>>>

	runtime {
		memory: "50G"
		cpu: 10
		sge_queue: "threaded.q"
	}

	output {
		File filter1_vcf = "out.vcf.gz"
	}

}

task filter_genome_wide {
	input {
		Array[File] filter1_vcfs
		Path path
	}
	command <<<
		# merge vcf file for all chromosomes
		~{path.bcftools_exe} concat ~{sep=" " filter1_vcfs} -o a1.vcf.gz

		# Only keep sites with missing < 0.7
		~{path.bcftools_exe} view -i 'F_MISSING<0.7' a1.vcf.gz -Oz -o a2.vcf.gz 

		# find sample with missing genotype < 0.7 (to keep)
		~{path.bcftools_exe} stats -s - a2.vcf.gz \
			| awk -v FS='\t' '
				$1=="SN" && $3=="number of records:" {
					num_rec=$4
				}; 
				$1=="PSC" {
					n_miss = $(14); 
					sample=$3; 
					if(n_miss/num_rec < 0.7) print sample;
				}' > sample_to_keep.txt

		~{path.bcftools_exe} view -S sample_to_keep.txt a2.vcf.gz -Oz -o a3.vcf.gz 

		# Remove singletons
		~{path.bcftools_exe} view  -i 'MAC>1' a3.vcf.gz -Oz -o a4.vcf.gz

		# Only keep sites with missing < 0.8 after removing singletons
		~{path.bcftools_exe} view -i 'F_MISSING<0.8' a4.vcf.gz -Oz -o a5.vcf.gz 

		mv a5.vcf.gz ~{path.out_dir}/filtered.vcf.gz
		rm  ./*.vcf.gz
	>>>

	runtime {
		memory: "50G"
		cpu: 10
		sge_queue: "threaded.q"
	}

	
	output {
		String filtered_vcf = path.out_dir + '/filtered.vcf.gz'
	}
}
