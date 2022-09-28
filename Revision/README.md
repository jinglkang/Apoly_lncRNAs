## overlap at least two program
```bash
# Kang@fishlab3 Mon Oct 18 17:08:52 ~/Desktop/LncRNA
vi Overlap_CPAT_CPC_FEELnc.lncRNA.gene.txt    # 9209 genes
perl temp10.pl Overlap_CPAT_CPC_merged_all_ind_CPAT.xls CPC_all_trans_nc.txt CPAT_all_trans_nc.txt FEELnc_all_trans_nc.txtNA.gene.txt >Candidate_lncRNA.gtf
FEELnc_classifier.pl --lncrna=Candidate_lncRNA.gtf --mrna=apoly_primary_geneannotation_v1.gtf --log=classifier/Apoly.classifier.log > classifier/Candidate_lncRNA.classifier.txt
# Kang@fishlab3 Mon Oct 18 19:59:47 ~/Desktop/LncRNA/classifier
less Candidate_lncRNA.classifier.txt|perl -alne 'my $info;for (my $i = 1; $i < @F; $i++){$info.=$F[$i]."\t"};$info=~s/\s+$//;print $info if /^isBest/ || /^1/' >Candidate_lncRNA.classifier.best.txt
perl temp1.pl Candidate_lncRNA.classifier.txt > Candidate_lncRNA.classifier.type.txt
```
## the location information of lncRNA
```bash
# Kang@fishlab3 Wed Nov 17 22:14:22 ~/Desktop/LncRNA/classifier
perl temp2.pl Candidate_lncRNA.integenic.txt >Candidate_lncRNA.integenic.location.txt
# summary the information of LncRNA genes
perl temp4.pl >Candidate_lncRNA_genes.classifier.txt
# add the location information
perl temp5.pl Candidate_lncRNA_genes.classifier.txt >Candidate_lncRNA_genes.classifier.final.txt
## add the LncRNA without interaction
cat Candidate_lncRNA_genes.classifier.final.txt Candidate_lncRNA.no_interaction.location.txt >Candidate_lncRNA_genes.classifier.final.1.txt
mv Candidate_lncRNA_genes.classifier.final.1.txt Candidate_lncRNA_genes.classifier.final.txt 
## 统计数目和平均距离
perl temp10.pl
# 1-1	6317	6317	12911
# 1-10-1	288	288	12328
# 1-10-2	15	50	9727
# 10-1	8	8	9752
# 10-2	1	11	632

# 统计平均距离及integenic LncRNA数目
perl temp12.pl
# <=1000	968
# >1000 & <=5000	2094
# >5000 & <=10000 1144
# >10000	2423
perl temp15.pl >raw_distance_to_gene.txt

# filter according to reads number
perl temp13.pl >Candidate_lncRNA_genes.classifier.final.filter_reads.txt

## the LncRNA gene number in each scaffold
perl temp6.pl >LncRNA_gene_location_scaffold.nb.txt
```

## Use IGV
Load gtf file: File --> load from file (select the gtf file);   
Load genome: genome --> load Genome from file (select the fasta file)    

## check the GO functions
```bash
# kangjingliang@kangjingliangdeMacBook-Pro 一  9 12 07:59:47 ~/Documents/2022/Apoly_lncRNA
less HE_lincRNAs_enrichment_nonfilter.txt|perl -alne 'next if ! /^\[OVER\]/;@F=split /\t/;@F=split /\t/;print "$F[2]\t$F[-3]" if $F[3] eq "BIOLOGICAL_PROCESS"'|perl -alne '@F=split /\t/;$hash{$F[0]}=$F[1];END{foreach $key (sort {$hash{$b}<=>$hash{$a}} keys %hash){print "$key\t$hash{$key}"}}'|less
extract_gene_functions -i HE_lincRNAs_nofilter_enrichment.txt -a Apoly_gene_annotation.txt --gene_column 1 --func_column 3 --functions HE_lincRNAs_over_function.txt --output HE_lincRNAs_over_function_gene

less Final_enrichment.txt|perl -alne 'next if ! /^\[OVER\]/;@F=split /\t/;@F=split /\t/;print "$F[2]\t$F[3]\t$F[-2]"'|perl -alne '@F=split /\t/;$name=$F[0]."\t".$F[1];$hash{$name}=$F[2];END{foreach $key (sort {$hash{$b}<=>$hash{$a}} keys %hash){print "$key\t$hash{$key}"}}'|less

extract_gene_functions -i HE_lincRNAs_nofilter_enrichment.txt -a Apoly_gene_annotation.txt --gene_column 1 --func_column 3 --functions ion_transport_funcion.txt --output ion_transport_funcion_gene

cat turquoise_lincRNAs_turquoise_neighbouring_genes.txt wgcna_lincRNAs_DE_neighbouring_genes.txt DE_lincRNAs_DE_neighbouring_genes.txt >Final_lincRNAs_neighbouring_gene_CO2.txt
```

## lincRNAs information
```bash
# Kang@fishlab3 Mon Sep 12 08:09:14 ~/Desktop/LncRNA/classifier
cat Candidate_lncRNA_genes.classifier.final.txt Candidate_lncRNA.no_interaction.location.txt >Candidate_lncRNA_genes.classifier.final.1.txt
mv Candidate_lncRNA_genes.classifier.final.1.txt Candidate_lncRNA_genes.classifier.final.txt

# Kang@fishlab3 Mon Sep 12 08:16:33 ~/Desktop/LncRNA/classifier
perl temp17.pl >Candidate_lncRNA_genes.classifier.final.filter_reads.txt

# kangjingliang@kangjingliangdeMacBook-Pro 四  6 09 11:16:38 ~/Documents/2021/lncRNA/LncRNA_write/2021-11-18
less Candidate_lncRNA_genes.classifier.final.filter_reads.txt # list for final lncRNA
# orientation
# Kang@fishlab3 Fri Sep 09 11:14:22 ~/Desktop/LncRNA/classifier
perl temp16.pl|perl -alne '$hash{$F[-1]}++;END{foreach my $a(keys %hash){print "$a\t$hash{$a}"}}'
#antisense	354
#strand_unknow	1057
#sense	243

# Kang@fishlab3 Mon Sep 12 08:23:54 ~/Desktop/LncRNA/classifier
perl temp16.pl|perl -alne '$aa=$F[-2]."\t".$F[-1];$hash{$aa}++;END{foreach my $a(keys %hash){print "$a\t$hash{$a}"}}'
# antisense	convergent;divergent	2
# antisense	divergent;convergent	2
# strand_unknow	strand(s)	1057
# antisense	divergent	146
# antisense	convergent	204
# sense	same_strand	243

# summary the exon number of lincRNAs
# Kang@fishlab3 Mon Sep 12 16:52:20 ~/Desktop/LncRNA/classifier
perl temp18.pl|perl -alne '$hash{$F[1]}++;END{foreach my $key (sort {$hash{$b}<=>$hash{$a}}keys %hash){print "$key\t$hash{$key}"}}'
# 1	2143
# 2	892
# 3	155
# 4	51
# 5	22
# 6	7
# 8	1
# 7	1
perl temp18.pl|sort -k2n >lincRNAs_exon.txt
less lincRNAs_exon.txt|perl -alne '$hash{$F[1]}++;END{foreach my $a(sort keys %hash){print "$a\t$hash{$a}"}}' >lincRNAs_exon_nb.txt
perl temp25.pl >Coding_exon.txt
perl temp26.pl >coding_lincRNAs_exon.txt

# the length of lincRNAs; if the lincRNA contain many transcripts, select the longest one
# Kang@fishlab3 Tue Sep 13 10:08:57 ~/Desktop/LncRNA/classifier
perl temp19.pl|sort -k3n >lincRNAs_length.txt
perl temp20.pl >lincRNAs_len_nb.txt
perl temp23.pl|sort -k3n >coding_length.txt
perl temp24.pl >coding_lincRNAs_length.txt

# count GC content
# Kang@fishlab3 Tue Sep 13 14:12:04 ~/Desktop/LncRNA/classifier
perl temp21.pl|sort -k3n >lincRNAs_gc.txt
perl temp22.pl|sort -k3n >coding_gc.txt
cat lincRNAs_gc.txt coding_gc.txt >gc_content.txt

# the number of coding genes and lincRNAs along genome
# kangjingliang@kangjingliangdeMacBook-Pro 二  9 13 21:28:26 ~/Documents/2021/lncRNA/LncRNA_write/circos
perl temp7.pl >coding_lincRNAs_nb.txt
```
## Expression of lincRNAs and coding genes
TPM normalization   
```R
countdata<-read.table("featureCounts_all.txt",skip = 1,sep="\t",header = T,row.names = 1)
metadata <- countdata[,1:5]
countdata <- countdata[,6:ncol(countdata)]
prefix<-"featureCounts_all"
kb <- metadata$Length / 1000
rpk <- countdata / kb
tpm <- t(t(rpk)/colSums(rpk) * 1000000)
write.csv(tpm,paste0(prefix,"_tpm.csv"))
q()
```
filtering low expressed transcripts    
```bash
# kangjingliang@kangjingliangdeMacBook-Pro 三  9 14 10:03:38 ~/Documents/2022/Apoly_lncRNA
less featureCounts_all_tpm.csv|perl -alne 'if (/bam/i){s/\.sorted\.bam//ig;print}else{@a=split /\,/;my $j;for($i=1;$i<@a;$i++){$j++ if $a[$i]>=1};print if $j>=204}' >featureCounts_all_tpm_filter.csv
less featureCounts_all_tpm_filter.csv|perl -alne 's/\"//g;print"Transcript\tAverage_TPM" if /dark/i; next if /dark/i;@a=split /\,/;my ($total, $ave);for($i=1;$i<@a;$i++){$total+=$a[$i]};$ave=$total/226;$ave=sprintf("%.0f",$ave);print "$a[0]\t$ave"' >featureCounts_all_tpm_filter_average.txt
# kangjingliang@kangjingliangdeMacBook-Pro 三  9 14 10:33:50 ~/Documents/2022/Apoly_lncRNA/revison
perl temp1.pl >coding_lincRNAs_average_tpm.txt

# check the expression correlation between lincRNAs and coding genes
# Kang@fishlab3 Wed Sep 14 22:49:20 ~/Desktop/LncRNA/classifier
nohup perl temp27.pl >lincRNAs_coding_spearman.txt &
#[1] 12368
# too long time
# only compute the lincRNAs and coding genes in the WGCNA module
# Kang@fishlab3 Thu Sep 15 09:31:27 ~/Desktop/LncRNA/classifier
mkdir spearman
cp ../run_testR.pl ./; cp ../featureCounts_all_tpm_filter.csv ./
# Kang@fishlab3 Thu Sep 15 10:44:28 ~/Desktop/LncRNA/classifier/spearman
perl temp27.pl Turqouise_lincRNAs.txt Turqouise_coding.txt >lincRNAs_coding_spearman.txt &
# [1] 14797
# divide Turqouise_coding.txt into 20 small files
split -l 313 Turqouise_coding.txt Turqouise_coding_split
# Kang@fishlab3 Thu Sep 15 11:34:53 ~/Desktop/LncRNA/classifier/spearman
nohup perl temp1.pl >run_spearman.process 2>&1 &
# [1] 22359
# only use PNG inds
# Kang@fishlab3 Thu Sep 15 12:17:34 ~/Desktop/LncRNA/classifier/spearman
perl temp2.pl >featureCounts_all_tpm_filter_apoly.csv
nohup perl temp1.pl >run_spearman.process 2>&1 &
# [1] 28865

# Check the correlation between DEGs and DElincRNAs
# SNORLAX
# (base) kang1234@celia-PowerEdge-T640 Thu Sep 15 16:06:58 ~/lncRNA_revison
split -l 172 DE_codings.txt DE_codings_
nohup perl temp1.pl >run_spearman.process 2>&1 &
# [1] 5857
# the number of coding genes that are significantly related with lincRNAs
# (base) kang1234@celia-PowerEdge-T640 Fri Sep 16 06:45:54 ~/lncRNA_revison
for i in *_spearman.txt;do less ${i}|perl -alne 'print if $F[-2]<=0.01 && abs($F[-1])>=0.9';done|perl -alne '$hash{$F[1]}.=$F[0]."\t";END{foreach my $key(sort keys %hash){print "$key\t$hash{$key}"}}'|wc -l
# 148

# Kang@fishlab3 Fri Sep 16 06:48:58 ~/Desktop/LncRNA/classifier/spearman
for i in *_spearman.txt;do less ${i}|perl -alne 'print if $F[-2]<=0.01 && abs($F[-1])>=0.9';done|perl -alne '$hash{$F[1]}.=$F[0]."\t";END{foreach my $key(sort keys %hash){print "$key\t$hash{$key}"}}'|wc -l
# 219

# The left two modules 
# 1. skyblue
# Kang@fishlab3 Fri Sep 16 06:54:39 ~/Desktop/LncRNA/classifier/spearman
mkdir skyblue
# Kang@fishlab3 Fri Sep 16 06:58:20 ~/Desktop/LncRNA/classifier/spearman/skyblue
cp ../featureCounts_all_tpm_filter_apoly.csv ../temp27.pl ../temp1.pl ../test.R ./
vi skyblue_lincRNAs.txt; vi skyblue_codings.txt
split -l 14 skyblue_codings.txt skyblue_codings_
# Kang@fishlab3 Fri Sep 16 07:04:39 ~/Desktop/LncRNA/classifier/spearman/skyblue
nohup perl temp1.pl >run_spearman.process 2>&1 &
# [1] 25903
# 2. darkgreen
# (base) kang1234@celia-PowerEdge-T640 Fri Sep 16 07:05:50 ~/lncRNA_revison
mkdir darkgreen
# (base) kang1234@celia-PowerEdge-T640 Fri Sep 16 07:06:31 ~/lncRNA_revison/darkgreen
cp ../featureCounts_all_tpm_filter_apoly.csv ../temp27.pl ../temp1.pl ../test.R ./
vi darkgreen_codings.txt; vi darkgreen_lincRNAs.txt
split -l 24 darkgreen_codings.txt darkgreen_codings_
# (base) kang1234@celia-PowerEdge-T640 Fri Sep 16 07:10:32 ~/lncRNA_revison/darkgreen
nohup perl temp1.pl >run_spearman.process 2>&1 &
# [1] 27408

# p<=0.01; rho>=0.9
# 1. DE lincRNAs and coding genes
# (base) kang1234@celia-PowerEdge-T640 Fri Sep 16 10:18:28 ~/lncRNA_revison
perl temp2.pl DE_all_spearman.txt >DE_coding_lincRNA_correlation.txt # 129
perl temp3.pl DE_all_spearman.txt >DE_lincRNA_coding_correlation.txt # 14 lincRNAs in total; a DE lincRNA might be correlated with maximum 53 DEGs
perl temp4.pl DE_coding_lincRNA_correlation.txt >DE_coding_lincRNA_correlation_neighbour.txt
# 2. WGCNA: darkgreen
# (base) kang1234@celia-PowerEdge-T640 Fri Sep 16 10:44:03 ~/lncRNA_revison/darkgreen
cp temp2.pl temp3.pl Apoly_gene_annotation.txt darkgreen/
cat *_spearman.txt >darkgreen_all_spearman.txt
perl temp2.pl darkgreen_all_spearman.txt # 0 lincRNAs


# p<=0.01; rho>=0.9
# 3. WGCNA: Turqouise
# Kang@fishlab3 Fri Sep 16 10:48:08 ~/Desktop/LncRNA/classifier/spearman
scp kang1234@147.8.76.177:~/lncRNA_revison/darkgreen/temp2.pl temp3.pl
scp kang1234@147.8.76.177:~/lncRNA_revison/darkgreen/temp3.pl temp4.pl
scp kang1234@147.8.76.177:~/lncRNA_revison/darkgreen/Apoly_gene_annotation.txt
cat *_spearman.txt >Turqouise_all_spearman.txt
perl temp3.pl Turqouise_all_spearman.txt >Turqouise_coding_lincRNA_correlation.txt # 129
perl temp4.pl Turqouise_all_spearman.txt >Turqouise_lincRNA_coding_correlation.txt # 14 lincRNAs in total; a lincRNA in Turqouise might be correlated with maximum 87 DEGs
# The coding gene and lincRNA is pair or not 
perl temp5.pl Turqouise_coding_lincRNA_correlation.txt >Turqouise_coding_lincRNA_correlation_neighbour.txt

# 4. WGCNA: skyblue
# Kang@fishlab3 Fri Sep 16 10:59:51 ~/Desktop/LncRNA/classifier/spearman
cp temp3.pl temp4.pl Apoly_gene_annotation.txt skyblue/
# Kang@fishlab3 Thu Sep 22 11:56:17 ~/Desktop/LncRNA/classifier/spearman/skyblue
cat *_spearman.txt > skyblue_all_spearman.txt
perl temp3.pl skyblue_all_spearman.txt # 0 lincRNAs
```

## Enrichment
```bash
# kangjingliang@kangjingliangdeMacBook-Pro 三  9 21 09:57:25 ~/Documents/2022/Apoly_lncRNA/Enrichment
less Turqouise_coding_lincRNA_correlation_neighbour.txt|cut -f 1|perl -alne '$a=$_."-mRNA-1";print $a' >Turqouise_coding_lincRNA_correlation.gene.txt
# kangjingliang@kangjingliangdeMacBook-Pro 三  9 21 09:57:54 ~/Documents/2022/Apoly_lncRNA/Enrichment
less DE_coding_lincRNA_correlation_neighbour.txt|cut -f 1|perl -alne '$a=$_."-mRNA-1";print $a' >DE_coding_lincRNA_correlation.gene.txt
# combine the two results
perl temp1.pl >combine_correlation.txt
less combine_correlation.txt|perl -alne 'next if /^Coding_gene/;$a=$F[0]."-mRNA-1";print $a' >combine_correlation_genes.txt
```
