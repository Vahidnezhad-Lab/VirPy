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

prog = "VirFy.py"

version = """%prog
ViriPy is free for non-commercial use without warranty.
============================================================================
"""

usage = """Usage: %prog[-h] -1 Read1.fastq -2 Read2.fastq -o outputDir -index human_reference -index_vir
                  virus_reference [-t # threads, default: 4] [-g Using gzip input files, default: False]"""


def main():
    parser = argparse.ArgumentParser(description='ViriPy: An All-In-One Workflow for Virus Detection, Quantification, and Visualization')

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

    parser.add_argument('-g', '--gzip', required=False,
                        metavar='True/False', default='False', type=str,
                        help="""Toggle "-g 'True'" for using gunzipped FASTQ input""")

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

    os.system('ulimit -n 2048')

    print("Aligning to human reference using STAR")

    def alignment():
        #cmd1='hisat2 -x '+index_dir+' -1 '+fq1+' -2'+fq2+' -S '+out+'/accepted_hits.sam -p '+n_thread+' --quiet'
        if gzip == "False":
            cmd1 = 'STAR --runThreadN ' + n_thread + ' --genomeDir ' + index_dir + ' --readFilesIn ' + fq1 + ' ' + fq2 + ' --outFileNamePrefix ' + out + '/accepted_hits. --outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate'
        elif gzip =="True":
            cmd1 = 'STAR --runThreadN ' + n_thread + ' --genomeDir ' + index_dir + ' --readFilesIn ' + fq1 + ' ' + fq2 + ' --outFileNamePrefix ' + out + '/accepted_hits. --outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate --readFilesCommand gunzip -c'
        else:
            print("Error: Invalid input for -g/--gzip. Input should be 'True' if using gunzipped files")
            sys.exit()
        print('Running ', cmd1)
        os.system(cmd1)

    alignment()

    print("Aligning to virus reference using HISAT2")

    def virus_alignment():
        cmd2 = 'hisat2 -x ' + index_vir + '/virus -1 ' + out + '/accepted_hits.Unmapped.out.mate1 -2 ' + out + '/accepted_hits.Unmapped.out.mate2 -S ' + out + '/unmapped_aln.sam -p ' + n_thread + ' --quiet'
        print('Running ', cmd2)
        os.system(cmd2)

    virus_alignment()

    print("Converting SAM to BAM")

    def sam_to_bam():
        cmd3 = 'samtools view -Sb -h ' + out + '/unmapped_aln.sam > ' + out + '/unmapped_aln.bam'
        print('Running ', cmd3)
        os.system(cmd3)

    sam_to_bam()

    print("Sorting BAM")

    def sort():
        cmd4 = 'samtools sort ' + out + '/unmapped_aln.bam -o' + out + '/unmapped_aln_Coord_sorted.bam'
        print('Running ', cmd4)
        os.system(cmd4)

        cmd5 = 'samtools sort -n ' + out + '/unmapped_aln.bam -o' + out + '/unmapped_aln_sorted.bam'
        print('Running ', cmd5)
        os.system(cmd5)
    sort()

    print("Indexing BAM")

    def index():
        cmd6 = '''samtools index ''' + out + "/unmapped_aln_Coord_sorted.bam" + ''' '''
        print('Running ', cmd6)
        os.system(cmd6)

    index()

    print("Processing Viral Counts using eXpress")

    def express():
        cmd7 = 'express ' + index_vir + '/viruses.fasta ' + out + '/unmapped_aln_sorted.bam -o ' + out + '/eXpress/'
        print('Running ', cmd7)
        os.system(cmd7)

    express()

    print("Getting Top Viruses")

    def getTopVirus():
        virus = []
        missingAnnot = open(out + '/featurecounts/missingAnnots.txt', 'w')
        topList = open(out + '/viruses_detected.txt', 'w')
        reader = csv.reader(open(out + "/eXpress/results.xprs"), delimiter="\t")

        topList.write('\t'.join(next(reader)))

        for l in sorted(reader, key=itemgetter(5), reverse=True):
            if int(float(l[5])) > 50:
                post = l[1].split("|")[-1]
                virus.append(post)
                topList.write('\t'.join(l))
                topList.write('\n')

        print("Quantifying Viral Features using featureCounts")
        os.system('mkdir ' + out + '/featurecounts')
        for v in virus:
            if os.path.isfile(index_vir + "/annotationFiles/" + v + ".gtf"):
                cmd8 = 'featurecounts -p -a ' + index_vir + '/annotationFiles/' + v + '.gtf ' + '-o ' + out + '/featureCounts/' + v + '_counts.txt ' + out + '/unmapped_aln_sorted.bam'
                os.system(cmd8)
            elif os.path.isfile(index_vir + "/annotationFiles/" + v + ".saf"):
                cmd8 = 'featurecounts -p -F SAF -a ' + index_vir + '/annotationFiles/' + v + '.saf ' + '-o ' + out + '/featureCounts/' + v + '_counts.txt ' + out + '/unmapped_aln_sorted.bam'
                os.system(cmd8)
            else:
                print('> There was no gene annotation file (.gtf/.saf) found for ' + v)
                missingAnnot.write('> There was no gene annotation file (.gtf/.saf) found for ' + v)
                missingAnnot.write('\n')
        missingAnnot.close()
        topList.close()

    getTopVirus()
    print("Outputs are stored in " + out + ". Thank you for using ViriPy!")

if __name__ == '__main__':
    main()
