import os
import pandas as pd

workdir: config['work_dir']

samples = pd.read_csv(config['samples'], sep='\t', index_col='Sample')

rule all:
	input:
		'summary/junction_counts.txt',
		'summary/feature_counts.txt'

rule download_genome:
	output:
		'%s.gz' % config['genome']
	threads: 1
	params:
		'--silent --create-dirs'
	shell:
		'curl %s {params} --output {output}' % config['genome_url']

rule unzip_genome:
	input:
		'%s.gz' % config['genome']
	output:
		config['genome']
	threads: 1
	shell:
		'gunzip {input}'

rule download_annotation:
	output:
		config['features']
	threads: 1
	params:
		'--silent --create-dirs'
	shell:
		'curl %s {params} --output {output}' % config['features_url']

rule trim:
	input:
		R1= lambda wildcards: '%s/%s' % (config['data_dir'], samples.loc[wildcards.sample, 'R1']),
		R2= lambda wildcards: '%s/%s' % (config['data_dir'], samples.loc[wildcards.sample, 'R2'])
	output:
		R1='trim/{sample}_R1.fastq.gz',
		R2='trim/{sample}_R2.fastq.gz',
		html='logs/trim_report/{sample}.html',
		json='logs/trim_report/{sample}.json'
	threads: config['trimming']['threads']
	log:
		'logs/trim/{sample}.log'
	params:
		'--umi --umi_loc=read1 --umi_len=10 --adapter_sequence CTGTCTCTTATACACATCT --adapter_sequence_r2 CTGTCTCTTATACACATCT'
	shell:
		'fastp {params} -w {threads} -i {input.R1} -I {input.R2} -o {output.R1} -O {output.R2} --html {output.html} --json {output.json} 2> {log}'

rule prepare_annotations:
	input:
		feature = config['features'],
		target = config['targets']
	output:
		genes = 'genome/genes.bed',
		exons = 'genome/exons.bed',
		introns = 'genome/introns.bed',
		target_introns = 'genome/target_introns.bed',
		alias = 'genome/alias.txt'
	threads: 1
	script:
		'prepare_annotations.py'

rule prepare_hisat_ss:
	input:
		exon = 'genome/exons.bed',
		intron = 'genome/introns.bed'
	output:
		exon = 'genome/hisat/hisat_exons.txt',
		ss = 'genome/hisat/hisat_ss.txt'
	threads: 1
	script:
		'prepare_hisat.py' 

rule build_hisat_index:
	input:
		genome = config['genome'],
		exon = 'genome/hisat/hisat_exons.txt',
		ss = 'genome/hisat/hisat_ss.txt'
	output:
		expand('genome/hisat/index.{number}.ht2', number=[1,2,3,4,5,6,7,8])
	threads: 1
	shell:
		'hisat2-build -q --ss {input.ss} --exon {input.exon} {input.genome} genome/hisat/index'

rule align:
	input:
		expand('genome/hisat/index.{number}.ht2', number=[1,2,3,4,5,6,7,8]),
		R1='trim/{sample}_R1.fastq.gz',
		R2='trim/{sample}_R2.fastq.gz'
	output:
		pipe('{sample}_unfiltered.pipe')
	threads: config['alignment']['threads']
	log:
		'logs/align/{sample}_alignment.log'
	params:
		'--new-summary --max-intronlen 2000 --no-unal --rna-strandness RF -x genome/hisat/index'
	shell:
		'hisat2 -p {threads} {params} -1 {input.R1} -2 {input.R2} > {output} 2> {log}'

rule filter_bam:
	input:
		'{sample}_unfiltered.pipe'
	output:
		pipe('{sample}_filtered.pipe')
	threads: 1
	params: 
		'-bh -q %d' % config['alignment']['mapQ_filter'] 
	shell:
		'samtools view {params} {input} > {output}'

rule sort_bam:
	input:
		'{sample}_filtered.pipe'
	output:
		'align/{sample}.bam'
	log:
		'logs/align/{sample}_sort.log'
	threads: 1
	shell:
		'samtools sort {input} -o {output} 2> {log}'

rule index_alignment:
	input:
		'align/{sample}.bam'
	output:
		'align/{sample}.bam.bai'
	threads: 1
	shell:
		'samtools index {input}'

rule junction_count:
	input:
		alignment = 'align/{sample}.bam',
		index = 'align/{sample}.bam.bai',
		genes = 'genome/genes.bed',
		exons = 'genome/exons.bed'
	output:
		'junctions/{sample}.txt'
	log:
		'logs/junctions/{sample}.log'
	threads: 1
	script:
		'count_junctions.py'

rule combine_junction_counts:
	input:
		expand('junctions/{sample}.txt', sample=samples.index.values)
	output:
		'summary/junction_counts.txt'
	threads: 1
	script:
		'combine.py'

rule feature_counts:
	input:
		alignment = 'align/{sample}.bam',
		index = 'align/{sample}.bam.bai',
		introns = 'genome/target_introns.bed'
	output:
		'features/{sample}.txt'
	threads: 1
	script:
		'feature_count.py'

rule combine_feature_count:
	input:
		expand('features/{sample}.txt', sample=samples.index.values)
	output:
		'summary/feature_counts.txt'
	threads: 1
	script:
		'combine_features.py' 
