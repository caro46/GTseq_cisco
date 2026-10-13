#!/bin/env python3

import os, sys, argparse, subprocess

def main():

	## arguments with help messages
	parser = argparse.ArgumentParser(description="Processing RADseq data in order to develop microhaplotypes GT-seq panels: read mapping against reference genome (bwa), genotypes calling (GATK and stacks) and filtering.", allow_abbrev=False)
	parser.add_argument("--fastq", nargs = "*", help = "Fastq file(s) to use. Default: No", default = []) #default = empty list
	parser.add_argument("--genome", type = str, help = "Reference genome to use. Default: No")
	parser.add_argument("--output_directory", type = str, default = "", help = "Ouptut directory for bam files. Default: No")
	parser.add_argument("--output_directory_vcf", type = str, default = "", help = "Output directory for vcf files. Default: No")
	parser.add_argument("--vcf_files", nargs = "*", default = [], help = "vcf files to filter individuals from. Default: No")
	parser.add_argument("--threads", type = int, help = "Number of threads. Default: 2", default = 2)
	parser.add_argument("--prefix_out", type = str, help = "Output prefix. Default: No")
	parser.add_argument("--step", type = str, help = "Step of the analysis: GenomeIndexing, ReadsMapping, GenotypeCalling. Default: ReadsMapping", default = "ReadsMapping")
	parser.add_argument("--Read1_info", type = str, help = "String distinguishing Read1 and Read2 files for Paired-end reads. Ex: _R1. Default: _R1", default = "_R1")
	parser.add_argument("--Read2_info", type = str, help = "String distinguishing Read1 and Read2 files for Paired-end reads. Ex: _R2. Default: _R2", default = "_R2")
	parser.add_argument("--gdb_path", type = str, help = "Path to GATK database for joining multi-individuals gvcf. Default: No")
	parser.add_argument("--temp_dir", type = str, default = "", help = "Temporary directory for GATK gvcf joining. Default: No")
	parser.add_argument("--pop_map", type = str, help = "Pop file for stacks population. Default: No")	
	parser.add_argument("--output_directory_stacks", type = str, default = "", help = "Ouptut directory for stacks files. Default: No")
	parser.add_argument("--numb_int", type = int, help = "Number of intervals for GATK joining. Default: 5", default = 5)
	parser.add_argument("--list_interval", type = str, help = "File containing the location to keep when calling with GATK. Default: No")
	parser.add_argument("--dp_filt_geno", type = int, help = "Filter on genotype depth for bcftools filter. Default: 4", default = 4)
	parser.add_argument("--perc_ind_filt", type = float, help = "Filter on percentage of individuals with missing data for bcftools filter. Default: 0.1", default = 0.1)
	parser.add_argument("--Geno_qual_filt", type = int, help = "Filter on Genotype quality for bcftools filter. Default: 20", default = 20)
	parser.add_argument("--maf", type = float, help = "Filter on minimum allele frequency. Default: 0.05", default = 0.05)

	## Parsing the arguments
	args = parser.parse_args()
#	print("Arguments: ",args,"\n")
	
	## Indexing genome for mapping
	if args.step == "GenomeIndexing":
		cmd = f"bwa index -a bwtsw {args.genome}"
		print(f"Command: {cmd}")
		process = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines = True, shell = True)
		process.wait()

	## Read mapping (paired-end)
	elif args.step == "ReadsMappingPE_RAD": #no mark duplicate for RAD data otherwise looses data
		for filename in args.fastq:
			if filename.endswith(".fastq.gz") | filename.endswith(".fq.gz"):
				if args.Read1_info in filename:
				## Extract individual IDs
					file_base=os.path.basename(filename)
					file_no_gz = os.path.splitext(file_base)[0]
					file_no_ext = os.path.splitext(file_no_gz)[0]
					#pair info is after a '.' so let's replace that with '_'
					file_no_ext = file_no_ext.replace('.', '_')
					Species, PopID, PopNum, Ind, pair  = file_no_ext.split("_")
					R1_filename = filename
					R2_filename = R1_filename.replace(args.Read1_info, args.Read2_info)
					cmd = f"""bwa mem -M -t {args.threads} -R \'@RG\\tID:{Species}_{PopID}_{PopNum}_{Ind}\\tSM:{Species}_{PopID}_{PopNum}_{Ind}\\tLB:library1\\tPL:illumina\' {args.genome} {R1_filename} {R2_filename} > {args.output_directory}/{Species}_{PopID}_{PopNum}_{Ind}.sam"""
					print(f"Command: {cmd}")
					process = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines = True, shell = True) #execute the command 
					process.wait()
					## Converting sam to bam
					cmd = f"samtools view --threads {args.threads} -bt {args.genome} -o {args.output_directory}/{Species}_{PopID}_{PopNum}_{Ind}.bam {args.output_directory}/{Species}_{PopID}_{PopNum}_{Ind}.sam"
					print(f"Command: {cmd}")
					process = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines = True, shell = True) #execute the command 
					process.wait()
					## Sorting bam (by position)
					cmd = f"samtools sort --threads {args.threads} {args.output_directory}/{Species}_{PopID}_{PopNum}_{Ind}.bam -o {args.output_directory}/{Species}_{PopID}_{PopNum}_{Ind}_sorted.bam"
					print(f"Command: {cmd}")
					process = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines = True, shell = True) #execute the command 
					process.wait()
					## Indexing bam file
					cmd = f"samtools index -@ {args.threads} {args.output_directory}/{Species}_{PopID}_{PopNum}_{Ind}_sorted.bam"
					print(f"Command: {cmd}")
					process = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines = True, shell = True) #execute the command 
					process.wait()
					## Deleting unecessary files
					cmd = f"rm -f {args.output_directory}/{Species}_{PopID}_{PopNum}_{Ind}.bam {args.output_directory}/{Species}_{PopID}_{PopNum}_{Ind}.sam"
					print(f"Command: {cmd}")
					process = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines = True, shell = True) #execute the command 
					process.wait()
					print(f"Done with Individual: {Ind}")

	## GATK pipeline

	## 1- GATK indexing
	elif args.step == "GATK_index":
		cmd_ref = f"samtools faidx {args.genome}"
		print(f"Command: {cmd_ref}")
		process = subprocess.Popen(cmd_ref, stdout=subprocess.PIPE, universal_newlines = True, shell = True) #execute the command
		process.wait()
		cmd_dic = f"/programs/gatk4/gatk CreateSequenceDictionary -R {args.genome}"
		print(f"Command: {cmd_dic}")
		process = subprocess.Popen(cmd_dic, stdout=subprocess.PIPE, universal_newlines = True, shell = True) #execute the command
		process.wait()

	## 2- Individual genotype calling (HaplotypeCaller)
	elif args.step == "GATK_HaplotypeCaller_gvcf":
		for filename in os.listdir(args.output_directory): 
			if filename.endswith("_sorted.bam"):
				file_base=os.path.basename(filename)
				file_no_bam = os.path.splitext(file_base)[0]
				cmd = f"/programs/gatk4/gatk --java-options \"-Xmx20G\" HaplotypeCaller -R {args.genome} -I {args.output_directory}/{filename} -O {args.output_directory_vcf}/{file_no_bam}.g.vcf.gz -ERC GVCF"
				print(f"Command: {cmd}")
				process = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines = True, shell = True) #execute the command
				process.wait()

	## 3- GATK database
	elif args.step == "GATK_db":
		db_cmd = f"/programs/gatk4/gatk --java-options \"-Xmx20G\" GenomicsDBImport --genomicsdb-workspace-path {args.gdb_path}" 
		cmd_files = ' '
		for filename in os.listdir(args.output_directory_vcf):
			if filename.endswith(".g.vcf.gz"):
				#cmd_files += filename + " "
				cmd_files += f"-V {args.output_directory_vcf}/{filename} "
		cmd_db_after = f"--tmp-dir \"{args.temp_dir}\" --max-num-intervals-to-import-in-parallel {args.numb_int} --intervals {args.list_interval}"
		cmd = db_cmd + cmd_files + cmd_db_after 
		print(f"Command: {cmd}")
		process = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines = True, shell = True)
		process.wait()

	## 4- Joining individual genotype calls
	elif args.step == "GATK_joint":
		cmd = f"/programs/gatk4/gatk --java-options \"-Xmx20G\" GenotypeGVCFs -R {args.genome} -V gendb://{args.gdb_path} -O {args.output_directory_vcf}/{args.prefix_out}.vcf.gz"
		print(f"Command: {cmd}")
		process = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines = True, shell = True)
		process.wait()

	## 5- GATK best practice with hard filtering: 
	elif args.step == "GATK_hard_filtering_BP":
		for filename in os.listdir(args.output_directory_vcf):
			if filename.endswith(".vcf.gz") and not filename.endswith("g.vcf.gz"): #add and if not g.vcf.gz
				file_base=os.path.basename(filename)
				file_no_gz = os.path.splitext(file_base)[0]
				file_no_ext = os.path.splitext(file_no_gz)[0] 
				cmd = f"/programs/gatk4/gatk VariantFiltration \
				-V {args.output_directory_vcf}/{filename} \
				-filter \"QD < 2.0\" --filter-name \"QD2\" \
				-filter \"QUAL < 30.0\" --filter-name \"QUAL30\" \
				-filter \"SOR > 3.0\" --filter-name \"SOR3\" \
				-filter \"FS > 60.0\" --filter-name \"FS60\" \
				-filter \"MQ < 40.0\" --filter-name \"MQ40\" \
				-filter \"MQRankSum < -12.5\" --filter-name \"MQRankSum-12.5\" \
				-filter \"ReadPosRankSum < -8.0\" --filter-name \"ReadPosRankSum-8\" \
				-O {args.output_directory_vcf}/{file_no_ext}_GATK_hard_filtered.vcf.gz"
				print(f"Command: {cmd}")
				process = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines = True, shell = True) #execute the command
				process.wait()
				cmd_exclude = f"/programs/gatk4/gatk --java-options -Xmx8G SelectVariants \
				--exclude-filtered \
				-V {args.output_directory_vcf}/{file_no_ext}_GATK_hard_filtered.vcf.gz \
				-O {args.output_directory_vcf}/{file_no_ext}_GATK_hard_filtered_removed.vcf.gz"
				print(f"Command: {cmd_exclude}")
				process = subprocess.Popen(cmd_exclude, stdout=subprocess.PIPE, universal_newlines = True, shell = True) #execute the command
				process.wait()

	## Running stacks (gstacks + populations) to out a vcf file
	elif args.step == "stacks_vcf":
		cmd = f"/programs/stacks-2.67/bin/gstacks -I {args.output_directory} -M {args.pop_map} -O {args.output_directory_stacks} -t {args.threads}"
		print(f"Command: {cmd}")
		process = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines = True, shell = True) #execute the command
		process.wait()
		cmd_pop = f"/programs/stacks-2.67/bin/populations -P {args.output_directory_stacks} -M {args.pop_map} --vcf -t {args.threads}"
		print(f"Command: {cmd_pop}")
		process = subprocess.Popen(cmd_pop, stdout=subprocess.PIPE, universal_newlines = True, shell = True) #execute the command
		process.wait()

	## Quality filtering of vcf (maf, geno qual, ind dp, missingness)
	elif args.step == "SNP_vcftools_SNPfiltering":
		#for filename in os.listdir(args.output_directory_vcf):
		for filename in args.vcf_files:
			if filename.endswith(".vcf.gz"):
				file_base=os.path.basename(filename)
				file_no_gz = os.path.splitext(file_base)[0]
				file_no_ext = os.path.splitext(file_no_gz)[0]
				path = os.path.dirname(filename)
				#cmd = f"vcftools --gzvcf {args.output_directory_vcf}/{filename} --maf {args.maf} --minGQ {args.Geno_qual_filt} --minDP {args.dp_filt_geno} --max-missing {args.perc_ind_filt} --recode --stdout | gzip -c > {args.output_directory_vcf}/{file_no_ext}_vcftools_filtered.vcf.gz" 
				cmd = f"vcftools --gzvcf {filename} --maf {args.maf} --minGQ {args.Geno_qual_filt} --minDP {args.dp_filt_geno} --max-missing {args.perc_ind_filt} --recode --stdout | gzip -c > {path}/{file_no_ext}_vcftools_filtered.vcf.gz" 
				print(f"Command: {cmd}")
				process = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines = True, shell = True) #execute the command
				process.wait()
			
if __name__ == "__main__":
	main()
