##################################################################################
# Author: Charles Huang (cyh002@students.jefferson.edu)
# Created Time: 2/12/2021
# Description: python script for virus detection and quantification from human RNA-seq samples
##################################################################################

import sys
import argparse
import os
import os.path
import csv
from operator import itemgetter

prog = "VirPy_v2.py"

version = """%prog
VirPy is free for non-commercial use without warranty.
============================================================================
"""

usage = """Usage: %prog[-h] -1 Read1.fastq -2 Read2.fastq -o outputDir -index human_reference -index_vir
                  virus_reference [-t # threads, default: 4] [-g Using gzip input files, default: False] [-c Minimum coverage, default: 12]"""


def main():
    parser = argparse.ArgumentParser(description='VirPy: An All-In-One Workflow for Virus Detection, Quantification, and Visualization')

    parser.add_argument('-t', '--n_thread', required=False, default='4',
                        type=str, help='Number of threads')

    parser.add_argument('-1', '--fq1', required=True, metavar='read1.fastq', type=str,
                        help='Path to Read 1 of the paired-end RNA-seq')

    parser.add_argument('-2', '--fq2', required=True, metavar='read2.fastq', type=str,
                        help='Path to Read 2 of the paired-end RNA-seq')

    parser.add_argument('-o', '--out', required=True, metavar='outputDir', type=str,
                        help='Path to the output directory')

    parser.add_argument('-index', '--index_dir', required=True, metavar='human_reference', type=str,
                        help='Path to directory containing index files for human genome')

    parser.add_argument('-index_vir', '--index_vir', required=True, metavar='virus_reference', type=str,
                        help='Path to directory containing index files and fasta for virus genomes')

    parser.add_argument('-c', '--min_coverage', required=False, metavar='min_coverage', default='12', type=str,
                        help='Minimum read coverage')

    parser.add_argument('-g', '--gzip', required=False,
                        metavar='True/False', default='False', type=str,
                        help="""Toggle "-g True" for using gunzipped FASTQ input""")

    parser.add_argument('-s', '--consensus', required=False,
                        metavar='True/False', default='False', type=str,
                        help="""Toggle "-s True" to generate consensus fasta file""")

    args = parser.parse_args()

    fq1 = os.path.abspath(args.fq1)
    fq2 = os.path.abspath(args.fq2)

    try:
        os.path.isfile(fq1)
    except IOError:
        print('Error: Unable to locate Read 1 FASTQ file. Please check your path and try again.')
        sys.exit()

    try:
        os.path.isfile(fq2)
    except IOError:
        print('Error: Unable to locate Read 2 FASTQ file. Please check your path and try again.')
        sys.exit()

    index_dir = os.path.abspath(args.index_dir)
    index_vir = os.path.abspath(args.index_vir)

    try:
        os.path.isfile(index_dir)
    except IOError:
        print('Error: Unable to locate human genome reference files. Please check your path and try again.')
        sys.exit()

    try:
        os.path.isfile(index_vir)
    except IOError:
        print('Error: Unable to locate viral genome reference files. Please check your path and try again.')
        sys.exit()

    out = os.path.abspath(args.out)
    n_thread = args.n_thread
    gzip = args.gzip
    min_coverage = args.min_coverage
    consensus = args.consensus

    os.system('ulimit -n 2048')

    print("Aligning to human reference using STAR")

    def alignment():
        #cmd1='hisat2 -x '+index_dir+' -1 '+fq1+' -2'+fq2+' -S '+out+'/accepted_hits.sam -p '+n_thread+' --quiet'
        if gzip == "False":
            cmd1 = 'STAR --runThreadN ' + n_thread + ' --genomeDir ' + index_dir + ' --readFilesIn ' + fq1 + ' ' + fq2 + ' --outFileNamePrefix ' + out + '/acceptedHits. --outSAMtype BAM SortedByCoordinate'
        elif gzip =="True":
            cmd1 = 'STAR --runThreadN ' + n_thread + ' --genomeDir ' + index_dir + ' --readFilesIn ' + fq1 + ' ' + fq2 + ' --outFileNamePrefix ' + out + '/acceptedHits. --outSAMtype BAM SortedByCoordinate --readFilesCommand gunzip -c'
        else:
            print("Error: Invalid input for -g/--gzip. Input should be 'True' if using gunzipped files")
            sys.exit()
        print('Running ', cmd1)
        os.system(cmd1)

    alignment()

    print("Resorting BAM")

    def sort():
        cmd5 = 'samtools sort -n ' + out + '/acceptedHits.Aligned.sortedByCoord.out.bam -o' + out + '/acceptedHits.Aligned.sorted.out.bam'
        print('Running ', cmd5)
        os.system(cmd5)
    sort()

    print("Indexing BAM")

    def index():
        cmd6 = 'samtools index ' + out + "/acceptedHits.Aligned.sortedByCoord.out.bam"
        print('Running ', cmd6)
        os.system(cmd6)

    index()

    print("Processing Viral Counts using eXpress")

    def express():
        cmd7 = 'express ' + index_vir + '/viruses.fasta ' + out + '/acceptedHits.Aligned.sorted.out.bam -o ' + out + '/eXpress/'
        print('Running ', cmd7)
        os.system(cmd7)

    express()

    def kallisto():
        cmd6 = 'kallisto quant -i ' + index_vir + '/viruses.idx -o ' + out + '/kallisto/ ' + out + '/accepted_hits.Unmapped.out.mate1 ' + out + '/accepted_hits.Unmapped.out.mate2'
        print('Running ', cmd6)
        os.system(cmd6)

    kallisto()

    print("Getting Top Viruses")

    virus = []
    vir = []

    def getTopVirus():
        os.system('mkdir ' + out + '/featurecounts')
        missingAnnot = open(out + '/featurecounts/missingAnnots.txt', 'w')
        topList = open(out + '/viruses_detected.txt', 'w')
        reader = csv.reader(open(out + "/eXpress/results.xprs"), delimiter="\t")

        topList.write('\t'.join(next(reader)))
        topList.write('\n')

        for l in sorted(reader, key=itemgetter(5), reverse=True):
            if int(float(l[5])) > 50:
                post = l[1]
                virus.append(post)
                topList.write('\t'.join(l))
                topList.write('\n')

        print("Quantifying Viral Features using featureCounts")

        for v in virus:
            i = v.split("|")[-1]
            if os.path.isfile(index_vir + "/annotationFiles/" + i + ".gtf"):
                cmd8 = 'featurecounts -p -a ' + index_vir + '/annotationFiles/' + i + '.gtf ' + '-o ' + out + '/featureCounts/' + i + '_counts.txt ' + out + '/unmapped_aln_sorted.bam'
                os.system(cmd8)
                vir.append(v)
            elif os.path.isfile(index_vir + "/annotationFiles/" + i + ".saf"):
                cmd8 = 'featurecounts -p -F SAF -a ' + index_vir + '/annotationFiles/' + i + '.saf ' + '-o ' + out + '/featureCounts/' + i + '_counts.txt ' + out + '/unmapped_aln_sorted.bam'
                os.system(cmd8)
                vir.append(v)
            else:
                print('> There was no gene annotation file (.gtf/.saf) found for ' + i)
                missingAnnot.write('> There was no gene annotation file (.gtf/.saf) found for ' + i)
                missingAnnot.write('\n')

        missingAnnot.close()
        topList.close()

    getTopVirus()

    def call_variants():
        os.system('mkdir ' + out + '/VariantCalling')
        cmd9 = 'freebayes -f ' + index_vir + '/viruses.fasta -b ' + out + '/unmapped_aln_Coord_sorted.bam --min-coverage ' + min_coverage + ' > ' + out + '/VariantCalling/variants.vcf'
        print('Running ', cmd9)
        os.system(cmd9)
        cmd10 = 'java -jar snpEff/snpEff.jar viruses ' + out + '/VariantCalling/variants.vcf > ' + out + '/VariantCalling/variants.ann.vcf'
        print('Running ', cmd10)
        os.system(cmd10)
        os.system('mv snpEff_genes.txt ' + out + '/VariantCalling')
        os.system('mv snpEff_summary.html ' + out + '/VariantCalling')
    call_variants()

    def consensus():
        cmd11 = './gatk IndexFeatureFile ' + out + '/VariantCalling/variants.ann.vcf'
        cmd12 = './gatk FastaAlternateReferenceMaker -R ' + index_vir + '/viruses.fasta -O ' + out + '/VariantCalling/consensus.fasta -V ' + out + '/VariantCalling/variants.ann.vcf'
        print('Running ', cmd11)
        os.system(cmd11)
        print('Running ', cmd12)
        os.system(cmd12)

    if consensus == 'True':
        consensus()

    print("Outputs are stored in " + out + ". Thank you for using VirPy!")

if __name__ == '__main__':
    main()
