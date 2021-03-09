# VirPy
VirPy is an all-in-one workflow for the detection, quantification, and visualization of viral species from RNA-seq data acquired from clinical samples. The objective of ViriPy is to provide an accessible and easy-to-use bioinformatics tool for investigation of viruses in clinically-derived RNA-seq data.

## Dependencies

Quality check and read trimming
* FASTQC (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* Trimmomatic (http://www.usadellab.org/cms/?page=trimmomatic)

Third-party software used in ViriPy
* STAR (https://github.com/alexdobin/STAR)
* HISAT2 (http://daehwankimlab.github.io/hisat2/)
* samtools (http://samtools.sourceforge.net/)
* eXpress (https://pachterlab.github.io/eXpress/index.html)
* Subread (http://subread.sourceforge.net/)
* iVar (https://andersen-lab.github.io/ivar/html/index.html)
* IGV (http://software.broadinstitute.org/software/igv/)

## Installation
```
git clone https://github.com/Vahidnezhad-Lab/ViriPy.git
cd ViriPy
```
## Building STAR Human Genome Index Files
Human genome assembly can be obtained from: http://hgdownload.soe.ucsc.edu/downloads.html

To build STAR genome index:
```
star --runMode genomeGenerate --runThreadN NumberOfThreads --genomeDir /path/to/genomeDir --genomeFastaFiles /path/to/genome/fasta
```
## Options
Required
* -1 read1.fastq, --fq1 read1.fastq
	- Path to Read 1 of the paired-end RNA-seq
* -2 read2.fastq, --fq2 read2.fastq
	- Path to Read 2 of the paired-end RNA-seq
* -o outputDir, --out outputDir	
	- Path to the output directory
* -index human_reference, --index_dir human_reference
	- Path to directory containing index files for human genome
* -index_vir virus_reference, --index_vir virus_reference
	- Path to directory containing index files and fasta for virus genomes

Optional
* -h, --help
	- Show this help message
* -t N_THREAD, --n_thread N_THREAD
	- Number of threads, default 4
* -g True/False, --gzip True/False
	- Toggle '-g True' for using gunzipped FASTQ input, default: False

## Example Usage
```
python VirPy.py -t 4 -1 Test.1.fastq -2 Test.2.fastq -o outputDir -index human_reference/ -index_vir viruses_reference/
```

## Output
* Unmapped_aln_Coord_sorted.bam
	- Coordinate-sorted, indexed BAM file aligned to virus genomes
* viruses_detected.txt
	- Tab-delimited text file containing potential viruses detected in the sample
* eXpress/results.xprs
	- Tab-delimited file with viral abundance quantification
* featureCounts/virus_counts.txt
	- Text file for each virus detected with viral feature quantification
* VariantCalling directory
	- Variant calling *.vcf and *.tsv files
