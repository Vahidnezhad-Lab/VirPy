# ViriPy
ViriPy is an all-in-one workflow for the detection, quantification, and visualization of viral species from RNA-seq data acquired from clinical samples. The objective of ViriPy is to provide an accessible and easy-to-use bioinformatics tool for investigation of viruses in clinically-derived RNA-seq data.

## Dependencies

* FASTQC (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* Trimmomatic (http://www.usadellab.org/cms/?page=trimmomatic)
* STAR (https://github.com/alexdobin/STAR)
* HISAT2 (http://daehwankimlab.github.io/hisat2/)
* samtools (http://samtools.sourceforge.net/)
* eXpress (https://pachterlab.github.io/eXpress/index.html)
* Subread (http://subread.sourceforge.net/)
* IGV (http://software.broadinstitute.org/software/igv/)

## Installation
```
git clone https://github.com/Vahidnezhad-Lab/ViriPy.git
cd ViriPy
```
## Building STAR Human Genome Index Files
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
* -index human_references, --index_dir human_references
	- Path to directory containing index files for human genome
* -index_vir virus_references, --index_vir virus_references
	- Path to directory containing index files and fasta for virus genomes

Optional
* -h, --help
	- Show this help message
* -t N_THREAD, --n_thread N_THREAD
	- Number of threads, default 4
* -gzip True/False, --gzip True/False
	- Toggle -g 'True' for using gunzipped FASTQ input

## Example Usage
```
python ViriFy.py -t 8 -1 Test1.fastq -2 Test2.fastq -o outputDir -index human_reference/ -index_vir viruses_reference/ --gzip True
```

## Output
* Unmapped_aln_Coord_sorted.bam
	- Coordinate-sorted, indexed BAM file aligned to virus genomes
* viruses_detected.txt
	- Tab-delimited text file containing potential viruses detected in the sample
* eXpress/results.xprs
	- Tab-delimited file with viral abundance quantification
* featureCounts/virus_counts.txt
	- Tect file for each virus detected with viral feature quantification
