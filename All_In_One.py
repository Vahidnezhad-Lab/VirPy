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

prog = "VirPy_v3.py"

version = """%prog
VirPy is free for non-commercial use without warranty.
============================================================================
"""

usage = """Usage: %prog[-h] -1 Read1.fastq -2 Read2.fastq -o outputDir -index human_reference -index_vir
                  virus_reference [-t # threads, default: 4] [-g Using gzip input files, default: False] [-c Minimum coverage, default: 12]"""


def main():
    parser = argparse.ArgumentParser(
        description='VirPy: An All-In-One Workflow for Virus Detection, Quantification, and Visualization')

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

    parser.add_argument('-transcriptome', '--transcriptome', required=True, metavar=' reference_transcriptome',
                        type=str,
                        help='Path to directory containing index files and fasta for reference transcriptome')

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
    transcriptome = os.path.abspath(args.transcriptome)

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

    try:
        os.path.isfile(transcriptome)
    except IOError:
        print('Error: Unable to locate viral genome reference files. Please check your path and try again.')
        sys.exit()

    out = os.path.abspath(args.out)
    n_thread = args.n_thread
    gzip = args.gzip
    min_coverage = args.min_coverage
    consensus = args.consensus

    human_out = out + '/Human/'
    virus_out = out + '/Virus/'

    os.system('mkdir ' + human_out)
    os.system('mkdir ' + virus_out)

    print("Aligning to human reference using STAR")

    #   STAR --runThreadN 10 --runMode genomeGenerate --genomeDir /path/to/genomeDir --genomeFastaFiles /path/to/genome/fasta1 /path/to/genome/fasta2 --sjdbGTFfile /path/to/annotations.gtf
    '''
    def alignment():
        #cmd1='hisat2 -x '+index_dir+' -1 '+fq1+' -2'+fq2+' -S '+out+'/accepted_hits.sam -p '+n_thread+' --quiet'
        if gzip == "False":
            cmd = 'STAR --runThreadN ' + n_thread + ' --genomeDir ' + index_dir + ' --readFilesIn ' + fq1 + ' ' + fq2 + ' --outFileNamePrefix ' + out + '/acceptedHits. --outSAMtype BAM SortedByCoordinate'
        elif gzip =="True":
            cmd = 'STAR --runThreadN ' + n_thread + ' --genomeDir ' + index_dir + ' --readFilesIn ' + fq1 + ' ' + fq2 + ' --outFileNamePrefix ' + out + '/acceptedHits. --outSAMtype BAM SortedByCoordinate --readFilesCommand gunzip -c'
        else:
            print("Error: Invalid input for -g/--gzip. Input should be 'True' if using gunzipped files")
            sys.exit()
        print('Running ', cmd)
        os.system(cmd)
    '''

    def alignment():
        # cmd1='hisat2 -x '+index_dir+' -1 '+fq1+' -2'+fq2+' -S '+out+'/accepted_hits.sam -p '+n_thread+' --quiet'
        if gzip == "False":
            cmd = 'STAR --runThreadN ' + n_thread + ' --genomeDir ' + index_dir + ' --readFilesIn ' + fq1 + ' ' + fq2 + ' --outFileNamePrefix ' + out + '/accepted_hits. --outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate'
        elif gzip == "True":
            cmd = 'STAR --runThreadN ' + n_thread + ' --genomeDir ' + index_dir + ' --readFilesIn ' + fq1 + ' ' + fq2 + ' --outFileNamePrefix ' + out + '/accepted_hits. --outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate --readFilesCommand gunzip -c'
        else:
            print("Error: Invalid input for -g/--gzip. Input should be 'True' if using gunzipped files")
            sys.exit()
        print('Running ', cmd)
        os.system(cmd)
    alignment()

    ###########################################################################################################################################################
    # RNA-SeqVariantCalling
    ###########################################################################################################################################################

    print("Aligning to reference transcriptome using BWA")
    def transcriptome_alignment():
        cmd = 'bwa mem ' + fq1 + ' ' + fq2 + ' ' + transcriptome + '/*.fasta' + ' -t ' + n_thread + ' > ' + human_out + 'trans.sam'
        print('Running ', cmd)
        os.system(cmd)
    transcriptome_alignment()

    print("Converting transcriptome SAM to BAM")
    def trans_sam_to_bam():
        cmd = 'samtools view -Sb -h ' + human_out + 'trans.sam > ' + human_out + 'trans.bam'
        print('Running ', cmd)
        os.system(cmd)
    trans_sam_to_bam()

    print("Sorting transcriptome BAM")
    def trans_sort():
        cmd = 'samtools sort -@ ' + n_thread + ' ' + human_out + 'trans.bam -o ' + human_out + 'trans_Coord_sorted.bam'
        print('Running ', cmd)
        os.system(cmd)
    trans_sort()

    print("Adding Read Groups")
    def addreadgroup():
        cmd = 'picard AddOrReplaceReadGroups ' + 'I=' + human_out + 'accepted_hits.Aligned.sortedByCoord.out.bam ' + 'O=' + human_out + 'RG.bam' + ' RGID=4 RGLB=twist RGPL=illumina RGPU=unit1 RGSM=' + human_out
        print('Running ', cmd)
        os.system(cmd)

        cmd1 = 'picard AddOrReplaceReadGroups ' + 'I=' + human_out + 'trans_Coord_sorted.bam' + 'O=' + human_out + 'Trans_RG.bam' + 'RGID=4 RGLB=twist RGPL=illumina RGPU=unit1 RGSM=' + human_out
        print('Running ', cmd1)
        os.system(cmd1)
    addreadgroup()

    print("Marking Duplicates")
    def markduplicate():
        cmd = 'picard MarkDuplicates ' + 'I=' + human_out + 'RG.bam' + 'O=' + human_out + 'MarkDup_RG.Bam' + 'M=' + human_out + 'marked_dup_metrics.tx'
        os.system(cmd)
        cmd1 = 'picard MarkDuplicates ' + 'I=' + human_out + 'Trans_RG.bam' + 'O=' + human_out + 'MarkDup_Trans_RG.bam ' + 'M=' + human_out + 'trans_marked_dup_metrics.tx'
        os.system(cmd1)
    markduplicate()

    print("Indexing marked bam")
    def Indexmarkedbam():
        cmd = 'picard BuildBamIndex' + 'I=' + human_out + 'MarkDup_RG.Bam'
        print('Running ', cmd)
        os.system(cmd)
        cmd1 = 'picard BuildBamIndex' + 'I=' + human_out + 'MarkDup_Trans_RG.bam '
        print('Running ', cmd1)
        os.system(cmd1)
    Indexmarkedbam()

    print("Splitting Cigar Reads")
    def splitNcigar():
        cmd = './gatk SplitNCigarReads -R ' + index_dir + '/*.fasta -I ' + human_out + 'MarkDup_RG.Bam -O ' + human_out + 'splitNcigar.bam '
        print('Running ', cmd)
        os.system(cmd)
    splitNcigar()

    print("Sammtools MPileUp")
    def variantCalling():
        cmd = 'samtools mpileup -f ' + transcriptome + '/*.fasta  -g ' + human_out + 'MarkDup_Trans_RG.bam ' + '|bcftools call -m > ' + human_out + 'trans.vcf'
        print('Running ', cmd)
        os.system(cmd)

        cmd1 = './gatk HaplotypeCaller -R ' + index_dir + '/*.fasta -I ' + human_out + 'splitNcigar.bam -O ' + 'gen.vcf'
        print('Running ', cmd1)
        os.system(cmd1)
    variantCalling()

    print("Merging VCFs")
    def vcfmerg():
        cmd = 'bcftools merg ' + human_out + '/*.vcf -Oz -o ' + human_out + 'Merged.vcf'
        print('Running ', cmd)
        os.system(cmd)
    vcfmerg()

    print("Homozygosity Mapping")
    def ROH():
        cmd = 'vcftools --vcf ' + '' + human_out + 'Merged.vcf --plink --out' + human_out + 'plink'
        print('Running ', cmd)
        os.system(cmd)
        cmd1 = 'plink --file ' + human_out + 'plink' + ' –homozyg '
        print('Running ', cmd1)
        os.system(cmd1)
    ROH()

    ###########################################################################################################################################################
    # Viral Variant Calling
    ###########################################################################################################################################################

    print("Aligning to virus reference using HISAT2")
    def virus_alignment():
        cmd = 'hisat2 -x ' + index_vir + '/viruses -1 ' + out + '/accepted_hits.Unmapped.out.mate1 -2 ' + out + '/accepted_hits.Unmapped.out.mate2 -S ' + out + '/unmapped_aln.sam -p ' + n_thread + ' --quiet'
        print('Running ', cmd)
        os.system(cmd)
    virus_alignment()

    print("Converting viral SAM to BAM")
    def sam_to_bam():
        cmd = 'samtools view -Sb -h ' + out + '/unmapped_aln.sam > ' + out + '/unmapped_aln.bam'
        print('Running ', cmd)
        os.system(cmd)
    sam_to_bam()

    print("Sorting BAM")
    def sort():
        cmd = 'samtools sort -@ ' + n_thread + ' ' + out + '/unmapped_aln.bam -o' + out + '/unmapped_aln_Coord_sorted.bam'
        print('Running ', cmd)
        os.system(cmd)

        cmd1 = 'samtools sort -n -@ ' + n_thread + ' ' + out + '/unmapped_aln.bam -o' + out + '/unmapped_aln_sorted.bam'
        print('Running ', cmd1)
        os.system(cmd1)
    sort()

    print("Indexing viral BAM")
    def viral_index():
        cmd = '''samtools index ''' + out + "/unmapped_aln_Coord_sorted.bam" + ''' '''
        print('Running ', cmd)
        os.system(cmd)
    viral_index()

    print("Processing Viral Counts using eXpress")
    def express():
        cmd = 'express ' + index_vir + '/viruses.fasta ' + out + '/unmapped_aln_sorted.bam --no-bias-correct -o ' + virus_out + '/eXpress/'
        print('Running ', cmd)
        os.system(cmd)
    express()

    print("Processing Viral Counts using StringTie")
    def stringtie():
        cmd = 'stringtie ' + out + '/unmapped_aln_Coord_sorted.bam -A -l' + virus_out + '/stringtie/viral_abundance.txt -o ' + virus_out + '/stringtie/stringtie.gtf -p ' + n_thread
        print('Running ', cmd)
        os.system(cmd)
    stringtie()

    print("Getting Top Viruses")
    virus = []
    vir = []
    def getTopVirus():
        os.system('mkdir ' + virus_out + '/featurecounts')
        missingAnnot = open(virus_out + '/featurecounts/missingAnnots.txt', 'w')
        topList = open(virus_out + '/viruses_detected.txt', 'w')
        reader_e = csv.reader(open(virus_out + "/eXpress/results.xprs"), delimiter="\t")

        topList.write('\t'.join(next(reader_e)))
        topList.write('\n')

        for l in sorted(reader_e, key=itemgetter(5), reverse=True):
            if int(float(l[5])) > 50:
                post = l[1]
                virus.append(post)
                topList.write('\t'.join(l))
                topList.write('\n')

        reader_s = csv.reader(open(virus_out + "/stringtie/viral_abundance.txt"), delimiter="\t")

        for l in sorted(reader_s, key=itemgetter(5), reverse=True)[1:]:
            if int(float(l[6])) > int(min_coverage):
                post = l[2]
                if post not in virus:
                    virus.append(post)

        print("Quantifying Viral Features using featureCounts")
        for v in virus:
            i = v.split("|")[-1]
            if os.path.isfile(index_vir + "/annotationFiles/" + i + ".gtf"):
                cmd = 'featurecounts -M -O -T ' + n_thread + ' -a ' + index_vir + '/annotationFiles/' + i + '.gtf ' + '-o ' + virus_out + '/featureCounts/' + i + '_counts.txt ' + out + '/unmapped_aln_sorted.bam'
                os.system(cmd)
                vir.append(v)
            elif os.path.isfile(index_vir + "/annotationFiles/" + i + ".saf"):
                cmd = 'featurecounts -M -O -T ' + n_thread + ' -F SAF -a ' + index_vir + '/annotationFiles/' + i + '.saf ' + '-o ' + virus_out + '/featureCounts/' + i + '_counts.txt ' + out + '/unmapped_aln_sorted.bam'
                os.system(cmd)
                vir.append(v)
            else:
                print('> There was no gene annotation file (.gtf/.saf) found for ' + i)
                missingAnnot.write('> There was no gene annotation file (.gtf/.saf) found for ' + i)
                missingAnnot.write('\n')

        missingAnnot.close()
        topList.close()
    getTopVirus()

    print("Calling Viral Variants")
    def call_variants():
        os.system('mkdir ' + virus_out + '/VariantCalling')
        cmd = 'freebayes -f ' + index_vir + '/viruses.fasta -b ' + out + '/unmapped_aln_Coord_sorted.bam --min-coverage ' + min_coverage + ' > ' + virus_out + '/VariantCalling/variants.vcf'
        print('Running ', cmd)
        os.system(cmd)

        cmd1 = 'java -jar snpEff/snpEff.jar viruses ' + out + '/VariantCalling/variants.vcf > ' + virus_out + '/VariantCalling/variants.ann.vcf'
        print('Running ', cmd1)
        os.system(cmd1)

        os.system('mv snpEff_genes.txt ' + virus_out + '/VariantCalling')
        os.system('mv snpEff_summary.html ' + virus_out + '/VariantCalling')
    call_variants()

    def generate_consensus():
        cmd = './gatk IndexFeatureFile ' + virus_out + '/VariantCalling/variants.ann.vcf'
        print('Running ', cmd)
        os.system(cmd)

        cmd1 = './gatk FastaAlternateReferenceMaker -R ' + index_vir + '/viruses.fasta -O ' + virus_out + '/VariantCalling/consensus.fasta -V ' + virus_out + '/VariantCalling/variants.ann.vcf'
        print('Running ', cmd1)
        os.system(cmd1)

    if consensus == 'True':
        print("Creating viral consensus sequence")
        generate_consensus()

if __name__ == '__main__':
    main()
