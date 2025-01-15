# To exclude issues related to server and/or software - rerun the problematic population

```bash
python3 fasta_to_vcf.py --step ReadsMappingPE_RAD --genome /home2/SardLab/cauret/reference_genomes/Coregonus_artedi_v1.fa.gz --threads 8 --fastq /home2/SardLab/cauret/Ackiss_RADseq_Data/cisco/rerun_A_MER/*.fq.gz --Read1_info .1 --Read2_info .2 --output_directory  /home2/SardLab/cauret/Ackiss_RADseq_Data/cisco/rerun_A_MER/bam_files
```

# Extract reads from one problematic individual

```bash
samtools view -u -f 12 -F 256 A_MER_XX_T5786_sorted.bam > A_MER_XX_T5786_sorted_unmapped.bam
bamToFastq -i A_MER_XX_T5786_sorted_unmapped.bam -fq A_MER_XX_T5786_unmapped.fastq
```

# Blast 

The 1st seven sequences on NCBI database, 1 has not hit but all the others map to pseudomonas.

# Download the bacteria genome that had good blast hits

```bash
#folder:/home2/SardLab/cauret/Ackiss_RADseq_Data/cisco/checking_contamination
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/834/565/GCF_009834565.1_ASM983456v1/GCF_009834565.1_ASM983456v1_genomic.fna.gz
```

# Mapping the reads from a problematic individual

```bash
python3 fasta_to_vcf.py --GenomeIndexing --genome GCF_009834565.1_ASM983456v1_genomic.fna.gz

python3 fasta_to_vcf.py --step ReadsMappingPE_RAD --genome GCF_009834565.1_ASM983456v1_genomic.fna.gz --threads 8 --fastq /home2/SardLab/cauret/Ackiss_RADseq_Data/cisco/rerun_A_MER/A_MER_XX_T5786*.fq.gz --Read1_info .1 --Read2_info .2 --output_directory /home2/SardLab/cauret/Ackiss_RADseq_Data/cisco/checking_contamination
```
