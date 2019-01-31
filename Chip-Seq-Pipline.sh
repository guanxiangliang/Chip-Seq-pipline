# All Unix-based command (Mac terminal)

# Enter into sample folder
    cd $sample/

# Quality check
    ##FASTQC

# Sequencing reads trimming (Optional for low quality reads)
    ##Trimmomatic
    java -jar /Users/guanxiang/Desktop/SW/Trimmomatic-0.35/trimmomatic-0.35.jar SE -threads 1 -phred64 -trimlog trimlog.txt  ../Mapping/NEZH2.fq out.fastq  ILLUMINACLIP:/Users/guanxiang/Desktop/SW/Trimmomatic-0.35/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20

# Mapping sequence to mm9 genome

    gunzip *.fq.gz
    bowtie -q -m 1 -v 0 --sam --best --strata ~/Desktop/RD/mm9/mm9 *.fq > good.sam
    samtools view -Sb good.sam > nonSorted.bam
    samtools sort nonSorted.bam -o Sorted.bam
    samtools rmdup -s Sorted.bam Sorted-rmdup.bam
    samtools index Sorted-rmdup.bam
    /Users/guanxiang/Desktop/SW/deepTools-2.1.1/bin/bamCoverage -b Sorted-rmdup.bam --normalizeTo1x 2150570000 --binSize 10 -e=200  -o good.bw
        #bowtie2 --local -q -x  ~/Desktop/RD/mm92/mm92 -U ../Mapping/NEZH2.fq -S teat-bowtie2.sam
        #bamToBed -i $sample.bam > $sample.bed
        #gzip $sample.sam
        #rm $sample.fq
        #rm $sample.sam



# Peaking callling using MACS
    #macs2 with broad peak calling
        macs2 callpeak -t $sample.bam -c control.bam -f BAM -g 1.87e9 -n AYY1 -B -q 0.1 -s 50 cutoff 0.1 --outdir AYY1_MACS2 --keep-dup 1 --SPMR
        macs2 callpeak -t $sample.bam -c control.bam -f BAM -g 1.87e9 -n AYY1 -B -q 0.1 -s 50 --broad --broad-cutoff 0.1 --outdir AYY1_MACS2 --keep-dup 1 --SPMR
    #macs2 callpeak -t $sample.bam -c control.bam -f BAM -g 1.87e9 -n AYY1 -B -q 0.05 -s 50 --outdir AYY1_MACS2 --keep-dup 1 --SPMR
    #macs14 not used any more
        #macs14 -t $sample.bam -c ../Input/control.bam -n $sample -g 1.87e9 > MACS.log
        #macs14 -t $sample.bam -c ../Input/control.bam -g 1.87e9 -n $sample -w --call-subpeaks

    # Visualization of track
        # MACS14
        cd $sample_MACS_wiggle/
        cd treat/ && cat * > $sample.wig.gz
        gunzip $sample.wig.gz
        sed '/track/d' $sample.wig | awk 'BEGIN{print "track type=wiggle_0 name="$sample_treat" description="Extended tag pileup from MACS version 1.4.2 20120305 every 10 bp""}{print $0}' > $sample-tem.wig && /Users/guanxiang/Desktop/SW/wigToBigWig $sample-tem.wig /Users/guanxiang/Desktop/RD/mm9m.size.txt $sample.bw
        # bedtools
        bedtools genomecov -ibam input.bam -bg -scale X -g genome.chrom.sizes > normalised.bg
        wigToBigWig -clip normalised.bg genome.chrom.sizes normalised.bw
            # Track file generation (Optional, since next step will generate wig files)
            genomeCoverageBed -ibam $sample.sorted.bam -bg -trackline -split -g mm9 > $sample.bedGraph
        # deeptools (Prefer)
        /Users/guanxiang/Desktop/SW/deepTools-2.1.1/bin/bamCoverage -b AYY1.bam --normalizeTo1x 2150570000 --binSize 10 -e=200  -o AYY1_normalized_2.bw
        # MACS2
        Convert bdg to bw
        # $sample.bw load into UCSC or IGV
        # Usiing Homer
            makeTagDirectory pre-YY1.fastq.sam  -format sam
            makeTagDirectory ./pre-YY1.fastq.sam  -format sam
            makeTagDirectory ./pre-YY1.fastq.sam
            makeTagDirectory homer-analysis  pre-YY1.fastq.sam
            makeUCSCfile homer-analysis/ -o pre-yy1
            makeUCSCfile homer-analysis/ -norm 1e6  -o auto

# Merge peaks
cat *bed | cut -f1-3 | sort -k1,1 -k2,2n | awk '{print $1"\t"$2"\t"$3}'| bedtools merge > merge.bed
    #filter blacklist
        bedtools intersect -a merge.bed -b mm9-blacklist.bed -v > all_merge.bed

#bedtools intersect -a Mcf7H3k27acUcdAlnRep1_peaks.filtered.bed -b Mcf7H3k27acUcdAlnRep2_peaks.filtered.bed -wa | cut -f1-3 | sort | uniq > Mcf7Rep1_peaks.bed
#bedtools intersect -a Mcf7H3k27acUcdAlnRep1_peaks.filtered.bed -b Mcf7H3k27acUcdAlnRep2_peaks.filtered.bed -wb | cut -f1-3 | sort | uniq > Mcf7Rep2_peaks.bed
#bedtools intersect -a Panc1H3k27acUcdAlnRep1_peaks.filtered.bed -b Panc1H3k27acUcdAlnRep2_peaks.filtered.bed -wa | cut -f1-3 | sort | uniq > Panc1Rep1_peaks.bed
#bedtools intersect -a Panc1H3k27acUcdAlnRep1_peaks.filtered.bed -b Panc1H3k27acUcdAlnRep2_peaks.filtered.bed -wb | cut -f1-3 | sort | uniq > Panc1Rep2_peaks.bed
#wc -l *
#rm *filtered*
#cat *bed | sort -k1,1 -k2,2n | bedtools merge | tee merge.bed | wc -l


# Count each sample
awk -F "\t" '{$1="peak_"NR FS$1;$4=$4FS"+"}1' all_merge.bed | awk 'BEGIN{print "GeneID""\t""Chr""\t""Start""\t""End""\t""Strand"}{print $1"\t"$2"\t"$3"\t"$4"\t"$5}' > merge.saf
/Users/guanxiang/Desktop/SW/subread-1.5.0-p1-source/bin/featureCounts -a merge.saf -F SAF -o $sample-counts.txt $sample.bam -T 1

# Heatmap
/Users/guanxiang/Desktop/SW/deepTools-2.1.1/bin/computeMatrix reference-point  -S /Volumes/Atchison-Lab-Storage/Chip-Seq-2015/Chip-Seq/*/Mapping/*_normalized.bw   -R /Users/guanxiang/Desktop/RD/mm9.TSS.bed -b 1000 -o test

/Users/guanxiang/Desktop/SW/deepTools-2.1.1/bin/plotHeatmap -m test  -out ExampleHeatmap1.png      --plotTitle "Test data with default settings"




# Annotation
    # Using CEAS
    ceas -g /Users/guanxiang/Desktop/SW/CEAS-Package-1.0.2/bin/mm9.refGene -b ../MAnome-Compare/sample1_peaks.bed
    # Usng Homer (Alternative)
    ## Format BED files (From MACS bed file and save as windows txt)
    annotatePeaks.pl AYY1-Homer.txt mm9 -annStats Stats  > outputfile.txt
        #(Alternative: annotated by Ensmbl)
        annotatePeaks.pl AYY1-Homer.txt mm9 -gtf /Users/guanxiang/Desktop/RD/mm9.gtf -annStats Stats  > outputfile.txt



## Converting bedGraph to BigWig
    sort -k1,1 -k2,2n $sample.bedGraph > sorted.bedGraph
    ### Remeber to move the title row from the end to the top using TXT editor
    /Users/guanxiang/Desktop/SW/bedGraphToBigWig sorted.bedGraph /Users/guanxiang/Desktop/RD/mm9.size.txt $sample.bw
## Converting Wig to BigWig
/Users/guanxiang/Desktop/SW/wigToBigWig $sample.wig /Users/guanxiang/Desktop/RD/mm9.size.txt $sample.bw

## Finding common peak
    intersectBed -a peak1.bed -b peak2.bed -u | sort -k1,1 -k2,2n -k3,3n > common_peak1.bed
    intersectBed -a peak2.bed -b peak1.bed -u | sort -k1,1 -k2,2n -k3,3n > common_peak2.bed
    intersectBed -a peak1.bed -b peak2.bed -v | sort -k1,1 -k2,2n -k3,3n > unique_peak1.bed
    intersectBed -a peak2.bed -b peak1.bed -v | sort -k1,1 -k2,2n -k3,3n > unique_peak2.bed
    cat common_peak1.bed common_peak2.bed | sort -k1,1 -k2,2n -k3,3n | mergeBed -i - > common_peak.bed

##
for i in *fold; do base=$(echo $i | awk -F. '{print $1}'); cut -f1-3 ${i}/${base}.fastq.gz.bam_peaks.gappedPeak > ${base}.bed;intersectBed -a yy1.bed -b ${base}.bed -u | sort -k1,1 -k2,2n -k3,3n > c1.bed;intersectBed -b yy1.bed -a ${base}.bed -u | sort -k1,1 -k2,2n -k3,3n > c2.bed;cat c1.bed c2.bed | sort -k1,1 -k2,2n -k3,3n | mergeBed -i - > ${base}.over.bed; done
for i in *.bed; do base=$(echo $i | awk -F. '{print $1}');intersectBed -a fyy1.bed -b ${base}.bed -u | sort -k1,1 -k2,2n -k3,3n > c1.bed;intersectBed -b fyy1.bed -a ${base}.bed -u | sort -k1,1 -k2,2n -k3,3n > c2.bed;cat c1.bed c2.bed | sort -k1,1 -k2,2n -k3,3n | mergeBed -i - > ${base}.over.bed;done




# Differentail peak analysis using MAnorm (example: compare AYY1 vs CYY1, put all four files in one folder; NOTE: all bed files should follow a unique format that MAnorm required:(1) The first 2 files have ONLY 3 columns: chromosome, start, end. (2)The next 2 files should have 4 columns: chromosome, start, end, strand (+/-))
cp $sample_peaks.bed ./ && cp $sample.bed
awk '{print $1,$2,$3}' $sample_peaks.bed > $sample_peaks-tem.bed && rm $sample_peaks.bed && mv $sample_peaks-tem.bed $sample_peaks.bed
awk '{print $1,$2,$3,$6}' $sample.bed > $sample-tem.bed && rm $sample.bed && mv $sample-tem.bed $sample.bed
./MAnorm.sh AYY1_peaks.bed CYY1_peaks.bed AYY1.bed CYY1.bed 78  157

#genomic location to sequence
http://genome.ucsc.edu/cgi-bin/das/mm9/dna?segment=chr12:114466450,114466850

#from bam to sequence
samtools view ASMC-bowtie2-local.bam chr12:114466450-114466850 > ASMC4-local.alg

#from name and fq to fq, extract name from fastq
seqtk subseq in.fq name.lst > out.fq

#length distribution in fastq
awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' file.fastq

#sam to fa
awk '{OFS="\t"; print ">"$1"\n"$10}' read.sam > read.fa

# SRR download
/Users/guanxiang/Desktop/SW/sratoolkit.2.5.7-mac64/bin/fastq-dump SRR2532636
/Users/guanxiang/Desktop/SW/sratoolkit.2.5.7-mac64/bin/fastq-dump --split-files SRR443885

# genome fa
../SW/twoBitToFa mm9.2bit mm9.fa

# from bed to sequence
samtools faidx mm9.fa
bedtools getfasta -fi Your_fasta_File.fa -bed Coordinates.bed -fo Output_Sequences.fa

# remove "^M"
perl -p -i -e "s/\r/\n/g" SE-pro-B.txt

# check folder size
du -sh ./*

#merge bed file
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' ncrna.bed | sort -k1,1 -k2,2n | bedtools merge -s -c 4 -o collapse -delim "|" > merge.bed
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' ncrna.bed | sort -k1,1 -k2,2n | bedtools merge -s -d 200 -c 4 -o collapse -delim "|" > merge2.bed

# bulk search bioproject
esearch -db bioproject -query "(B[All Fields] AND cell[All Fields] AND chip-seq[All Fields]) AND (mouse[Organism])" | efetch --format runinfor > chip-seq-B-cell-list
esearch -db pubmed -query "(((viromes[Text Word]) OR virome[Text Word]) AND (mice[text word] or mouse[Text Word] or human[text word] or humans[text word]))  " | efetch --format runinfor

#separate read1 and read2
cat H1.fastq | perl -ne '$o=$_; $w=<STDIN>; $r=<STDIN>; $_=<STDIN>; print "$o$w$r$_" if ($o =~ /\/1/); '

#python install
python setup.py install  --user


##change the maximal open files number
#Mac
sudo sysctl -w kern.maxfilesperproc=40000
sudo sysctl -w kern.maxfiles=40000
ulimit -S -n 400000
#Linux
set rlim_fd_max = 1663840
set rlim_fd_cur = 819200
ulimit -s -n 400000
set max_nprocs=6554600
set pidmax=10000000
set maxuprc=3000000
echo 200000 > /proc/sys/fs/file-max


# mount docker image
docker run -v /Volumes/Atchison-Lab-Storage/Juan/:/notebooks  -ti mgmarin/juncbase

# delete line contain pattern
sed '/pattern to match/d'
sed '/chrUn/d' sheep-intron.txt > sheep-intron.txt.2


cufflinks -g your-gtf-file -library-type  fr-unstranded(##note: this is for unstranded library or you can use "fr-firststrand" for stranded library)  your-bam-file.bam

genes.fpkm_tracking file  is for fpkm gene expression


find ./ -type f -exec cp '{}' ../upload/ \;

