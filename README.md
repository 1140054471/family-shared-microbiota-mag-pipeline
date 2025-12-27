# family-shared-microbiota-mag-pipeline
This repository contains the complete computational workflow and scripts for reproducing the analysis in ‘Types and characteristics of shared microbiota within families of ethnic minorities in Inner Mongolia’, from raw sequencing reads (64 samples) to metagenome-assembled genomes (MAGs) and downstream analyses.

The raw sequencing data generated in this study have been deposited in the NCBI SRA database under BioProject accession number PRJNA919082 and are accessible via the following URL: https://www.ncbi.nlm.nih.gov/sra/PRJNA919082.

1./#数据质控
klab_metaqc list -i Rawdata -s list_sample  /#生成样品列表
klab_metaqc qc -s list_sample -o qc_family_result -t human -j 3  /#质控 

2./#数据组装（MEGAHIT：ver.1.2.9）
mkdir 01megahitAssembly
for i in E*r1* ; do echo ${i%%.rmhost.r1.fq.gz} >> sample.list ; done
cat sample.list | parallel --dryrun megahit -t 10 --out-dir {}_Output_ass --out-prefix {} -1 {}.rmhost.r1.fq.gz -2 {}.rmhost.r2.fq.gz

3./#Metabat2分箱（MetaBAT 2：ver.2.12.1）
mkdir 02metabat2
for i in E* ; do echo ${i%%_Output_ass} >> sample.list ; done
cat sample.list | parallel --j 40 cp -p {}_Output_ass/{}.contigs.fa ../02metabat2
cat ../01megahitAssembly/sample.list | parallel --plus --j 40 190910_rename_genome.py -i {}.contigs.fa -m {} -o {}_rename.fa
cat ../01megahitAssembly/sample.list | parallel --plus --j 20 bwa index {}_rename.fa
cat ../01megahitAssembly/sample.list | parallel --plus --j 35 "bwa mem -t 32 {}_rename.fa ../../qc/{}.rmhost.r1.fq.gz ../../qc/{}.rmhost.r2.fq.gz | samtools view -bS - | samtools sort -@ 16 - -o {}_rename.fa.sort.bam"
cat ../01megahitAssembly/sample.list | parallel --plus --j 5 metabat2 -i {}/{}_rename.fa -a {}/{}.depth -o ./metabat2/{}/{} --minContig 2000 --unbinned
cat ../01megahitAssembly/sample.list | parallel --plus --dryrun checkm lineage_wf -x fa -t 40 --pplacer_threads 40 metabat2/{} metabat2/{}.metabat2.checkm
cat ../01megahitAssembly/sample.list | parallel --plus --dryrun summarize_checkm.py metabat2/{}.metabat2.checkm '>' metabat2/{}.metabat2.checkm.summary

4./#Vamb分箱（VAMB：ver1.0.1）
mkdir 03vamb
cat sample.list| parallel --plus --dryrun formatting_jgi.sh {}/{}.depth '>' {}/{}.depth.mo
cat sample.list | parallel --plus --j 10 vamb --fasta {}/{}_rename.fa --jgi {}/{}.depth.mo --outdir vamb/{}.t64.min2k 
-m 2000 -p 16 -t 64 --minfasta  200000 -o _

5./#采用DAS_Tool挑选vamb和metabat2组装结果中好的基因组
mkdir 04DAS_tools
ls -d E* | parallel -j 5 DAS_Tool -i {}/{}/metabat2_bins.scaffolds2bin.tsv,{}/{}.t64.min2k/bins/vamb_bins.scaffolds2bin.tsv -l metabat2,vamb -c ../02megahitbwadepth/{}/{}_rename.2k.fa -o {}/DAS_Tool -t 30 --search_engine diamond --write_bins
for i in `ls -d E*`; do Fasta_to_Scaffolds2Bin.sh -i ${i}/DAS_Tool_DASTool_bins/ -e fa > ${i}/DAS_Tool_DASTool_bins/${i}_DAS_Tools_bins.scaffolds2bin.tsv; done
cat sample.list |parallel --plus --j 5 checkm lineage_wf -x fa -t 40 --pplacer_threads 40 {}/DAS_Tool_DASTool_bins {}.checkm
cat ../sample.list | rush -j 10 'summarize_checkm.py {}.checkm > {}.checkm.s'

6./#挑选质量较高的基因组
awk '{print $9}' all_80_5_bins_checkm > all_80_5_bins_names 
cat all_80_5_bins_names | parallel -j 10 cp -r all.bins/{}.fa all_80_5_bins /#挑选出完整度大于等于80，污染率小于等于5的bins

7./#基因组聚类（drep：(ver.2.2.2）
cat all_80_5_bins_checkm | awk 'BEGIN{OFS=",";}{print $9 ".fa", $2, $3}' > ./bins_80_5_checkm_info
sed 's/rename_bins.fa,/genome,/g' bins_80_5_checkm_info > bins_80_5_checkm_info.csv /#这里第一列表头一定要是genome
rm -rf bins_80_5_checkm_info
dRep dereplicate dereplicate_result -pa 0.95 -sa 0.95 -g ./bins_80_5/* -p 32 --genomeInfo bins_80_5_checkm_info.csv
dRep cluster drep_out95 -pa 0.95 -sa 0.95 -g raw_hi_folder/*.fasta -p 30  /#去除MAG的冗余基因基因组
dRep dereplicate dRep_derep_out -g ./raw_hi_folder/*.fasta -p 30 -pa 0.95 -sa 0.95 

8./#NR数据注释
mkdir all_rename_SGBs_faa
ls -d all_rename_SGBs/*fa | parallel --plus --dryrun prodigal -i {} -a all_rename_SGBs_faa/{/.}.faa -f gff -o all_rename_SGBs_faa/{/.}.gff -q -d all_rename_SGBs_faa/{/.}.ffn    
cat all_rename_SGBs/*.fa > all_rename_bins.fa    
cat all_rename_SGBs_faa/*faa > all_combine.faa
diamond blastp --threads 32 --max-target-seqs 10 --db  /nvmessdnode3/opt/database/uniport/uniprot_trembl_sport.dmnd --query all_combine.faa --outfmt 6qseqid sseqid stitle pident qlen slen length mismatch gapopen qstart qend sstart send evalue bitscore --out all_combine.dia
mkdir all_rename_SGBs_speci_out
ls -d all_rename_SGBs_faa/*faa | parallel -j 30 specI.py {.}.ffn {} 2 {/.} all_rename_SGBs_speci_out
cat all_rename_SGBs_speci_out/S*/*results | cut -f 1,4,5 > all_rename_SGBs_speci_out.summary
ls -d all_rename_SGBs/*fa | parallel -j 20 so {} '>' {}.stat

9./#计算SGB含量（BWA-MEM：ver.0.7.17-r1188 ；Samtools：ver.1.6）
cat ../sample.list | parallel --dryrun '/opt/software/bwa-mem2-2.1_x64-linux/bwa-mem2 mem -t 16 bwa2_index/all_rename_bins.fa qc/{}.rmhost_1.fastq.gz qc/{}.rmhost_2.fastq.gz | samtools view -bS - -@ 5 | samtools sort -@ 5 - -o all_sample_bam_results/{}.sort.bam' | parallel -j 10 {}
mkdir all_bins_abundance
ls -d all_sample_bam_results/*bam | parallel --plus --dryrun 'source activate /nvmessdnode4/opt/conda/envs/CoverM; coverm genome --bam-files {} -d all_rename_SGBs -x fa -t 10 --min-read-percent-identity 0.95 --min-covered-fraction 0.4 -m covered_fraction mean relative_abundance rpkm > all_bins_abundance/{/..}.coverm' | parallel -j 10 {}
