# VirPy
VirPy is an all-in-one workflow for the detection, quantification, and visualization of viral species from RNA-seq data acquired from clinical samples. The objective of VirPy is to provide an accessible and easy-to-use bioinformatics tool for investigation of viruses in clinically-derived RNA-seq data.

## Dependencies

Quality check and read trimming
* FASTQC (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* Trimmomatic (http://www.usadellab.org/cms/?page=trimmomatic)

Third-party software used in VirPy
* STAR (https://github.com/alexdobin/STAR)
* HISAT2 (http://daehwankimlab.github.io/hisat2/)
* samtools (http://samtools.sourceforge.net/)
* eXpress (https://pachterlab.github.io/eXpress/index.html)
* Subread (http://subread.sourceforge.net/)
* Freebayes (https://github.com/freebayes/freebayes)
* SnpEff (https://pcingola.github.io/SnpEff/)
* GATK (https://gatk.broadinstitute.org/hc/en-us)
* Picard (https://broadinstitute.github.io/picard/)
* IGV (http://software.broadinstitute.org/software/igv/)

Packages can be installed using [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)
```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

conda install -y fastqc trimmomatic star hisat2 samtools subread freebayes
```

## Installation
```
git clone https://github.com/Vahidnezhad-Lab/VirPy.git
cd VirPy
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
	- Number of threads, default: 4
* -g True/False, --gzip True/False
	- Toggle '-g True' for using gunzipped FASTQ input, default: False
* -c MIN_COVERAGE, --min_coverage MIN_COVERAGE
	- Minimum read coverage, default: 12
* -s True/False, --consensus True/False
	- Toggle '-s True' to generate consensus sequence from VCF file, default: False

## Example Usage
```
ulimit -n 2048
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
