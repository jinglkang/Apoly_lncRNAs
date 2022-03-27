# Apoly lncRNA
## Putative lncRNA detected by CPC, CPAT, FEELnc
in Kang@fishlab3 Sat Mar 26 16:57:49 ~/Desktop/LncRNA   
**CPC_lncRNA_gene.txt**   
**CPAT_lncRNA_gene.txt**   
**FEELnc_lncRNA_gene.txt**   
***
## The final lincRNAs detected by which tools
The final lincRNAs: Kang@fishlab3 Sat Mar 26 16:57:49 ~/Desktop/LncRNA/classifier   
Candidate_lncRNA_genes.classifier.final.filter_reads.txt   
```bash
# Kang@fishlab3 Sat Mar 26 17:35:28 ~/Desktop/LncRNA
perl lincRNA_in_tools.pl|perl -alne '$hash{$F[-1]}++;END{foreach $key (keys %hash){print "$key\t$hash{$key}"}}'
# CPC;-;FEELnc	101
# CPC;CPAT;FEELnc	405
# CPC;CPAT;-	2174
# -;CPAT;FEELnc	592
```
***
## lincRNA detection by the three tools
### CPC
```bash
# Kang@fishlab3 Sat Mar 26 17:49:18 ~/Desktop/LncRNA
mv temp10.pl make_new_gtf.pl
perl make_new_gtf.pl CPC_lncRNA_gene.txt >CPC_lncRNA.gtf
mkdir CPC_classifier
FEELnc_classifier.pl --lncrna=CPC_lncRNA.gtf --mrna=apoly_primary_geneannotation_v1.gtf --log=CPC_classifier/CPC_lncRNA.log > CPC_classifier/CPC_lncRNA.classifier.txt
cd CPC_classifier
less CPC_lncRNA.classifier.txt|perl -alne 'my $info;for (my $i = 1; $i < @F; $i++){$info.=$F[$i]."\t"};$info=~s/\s+$//;print $info if /^isBest/ || /^1/'|perl -alne '@a=split /\t/;print if /intergenic/' >CPC_lincRNA.txt # 6538
```
### CPAT
```bash
# Kang@fishlab3 Sat Mar 26 17:49:18 ~/Desktop/LncRNA
perl make_new_gtf.pl CPAT_lncRNA_gene.txt >CPAT_lncRNA.gtf
mkdir CPAT_classifier
FEELnc_classifier.pl --lncrna=CPAT_lncRNA.gtf --mrna=apoly_primary_geneannotation_v1.gtf --log=CPAT_classifier/CPAT_lncRNA.log > CPAT_classifier/CPAT_lncRNA.classifier.txt
cd CPAT_classifier
less CPAT_lncRNA.classifier.txt|perl -alne 'my $info;for (my $i = 1; $i < @F; $i++){$info.=$F[$i]."\t"};$info=~s/\s+$//;print $info if /^isBest/ || /^1/'|perl -alne '@a=split /\t/;print if /intergenic/' >CPAT_lincRNA.txt # 14213
```
### FEELnc
```bash
# Kang@fishlab3 Sat Mar 26 17:49:18 ~/Desktop/LncRNA
perl make_new_gtf.pl FEELnc_lncRNA_gene.txt >FEELnc_lncRNA.gtf
mkdir FEELnc_classifier
FEELnc_classifier.pl --lncrna=FEELnc_lncRNA.gtf --mrna=apoly_primary_geneannotation_v1.gtf --log=FEELnc_classifier/FEELnc_lncRNA.log > FEELnc_classifier/FEELnc_lncRNA.classifier.txt
cd FEELnc_classifier
less FEELnc_lncRNA.classifier.txt|perl -alne 'my $info;for (my $i = 1; $i < @F; $i++){$info.=$F[$i]."\t"};$info=~s/\s+$//;print $info if /^isBest/ || /^1/'|perl -alne '@a=split /\t/;print if /intergenic/' >FEELnc_lincRNA.txt # 4318
```
