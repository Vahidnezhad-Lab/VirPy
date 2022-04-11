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
## Building HISAT Viral Genome Index Files
The viral genome index files are included within the data directory. If you are interested in adding custom viruses of interest to the directory, please recreate the viral genome index files using the following steps:
```
hisat2-build -p [THREADS] viruses.fasta viruses
```
More information can be found on the HISAT2 webpage: http://daehwankimlab.github.io/hisat2/howto/#building-indexes
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

## Setting Up Custom Viruses of Interest
Custom viruses can be added for investigators interested in viruses not already included within the data files. 
1. Append the genome file (.fasta) for the virus of interest to the "viruses.fasta" file following the same formatting of other entries. Many can be found on https://www.ncbi.nlm.nih.gov/. Re-build the HISAT2 Index file as described in the "Building HISAT Viral Genome Index Files" section above.
2. On the NCBI page, go to "Send to:" -> "File" -> "GFF3" -> "Create File".
3. Convert the GFF3 file to a GTF file. This can be accomplished by using genometools using:
```
conda install -c bioconda genometools
gt gff3_to_gtf [example.gff3]
```
*Note: Often you will have to add "##gff-version 3" to the header at the top of the gff3 file
4. Make sure to check that the ".gtf" file is named the same as it is in "viruses.fasta" file and is in the annotationFiles directory
5. Alternatively, a ".saf" file can be generated manually, using the following header:
```
GeneID	Chr	Start	END	Strand
```
