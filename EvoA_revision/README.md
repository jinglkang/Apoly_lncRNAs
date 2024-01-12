# Evolutionary applications: lncRNAs revision
```bash
# RNAcentral: The non-coding RNA sequence database (https://rnacentral.org/)
# rnacentral_species_specific_ids_lncRNAs.fasta
# index all rnacentral sequences
makeblastdb -in rnacentral_species_specific_ids.fasta -dbtype nucl -parse_seqids -out RNAcenDB

# only select lncRNAs sequences
# Kang@fishlab3 Wed Dec 27 21:47:19 /media/HDD
nohup perl temp1.pl > rnacentral_species_specific_ids_lncRNAs.fasta &

# index the lncRNAs nuc sequences
# Kang@fishlab3 Thu Dec 28 09:18:29 /media/HDD
makeblastdb -in rnacentral_species_specific_ids_lncRNAs.fasta -dbtype nucl -parse_seqids -out lncRNAsDB

# Kang@fishlab3 Wed Dec 27 11:42:44 ~/Desktop/LncRNA
vi lincRNAs_id.txt
# merged_all_ind.fa
# Kang@fishlab3 Wed Dec 27 14:34:19 ~/Desktop/LncRNA
perl temp13.pl > Final_lincRNAs.fas
scp Final_lincRNAs.fas jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/lncRNAs

# transfer the data to HPC
# jlkang@hpc2021 Wed Dec 27 18:59:44 /lustre1/g/sbs_schunter/Kang
mkdir lncRNAs
# Kang@fishlab3 Wed Dec 27 18:58:06 /media/HDD
scp rnacentral_species_specific_ids.fasta.gz jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/lncRNAs
# jlkang@hpc2021 Wed Dec 27 19:03:46 /lustre1/g/sbs_schunter/Kang/lncRNAs
module load blast-plus/2.13.0
# blastn -outfmt 6 -query query.fasta -out query.fa.bla -db Cleaner_wrasse -evalue 1e-10 -num_threads 30
# split as four files for blastn
split -l 1636 Final_lincRNAs.fas Final_lincRNAs_split
cat Final_lincRNAs_splita*.bla > Final_lincRNAs_Total.bla

# 558 of 3272 lincRNAs were blasted to lncRNAs
# jlkang@hpc2021 Thu Dec 28 11:21:30 /lustre1/g/sbs_schunter/Kang/lncRNAs
less Final_lincRNAs_bestblast.txt|perl -alne '@a=split /\t/;print $a[-1]'|perl -alne 's/\s+lncRNA\s+\(.*\)//;$hash{$_}++;END{foreach my $key (sort keys %hash){print "$key\t$hash{$key}"}}' > lncRNAs_blast_spe.txt

# kangjingliang@kangjingliangdeMacBook-Pro 四 12 28 11:26:40 ~/Documents/2023/lncRNA/EA_revision
scp jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/lncRNAs/Final_lincRNAs_bestblast.txt ./
scp jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/lncRNAs/lncRNAs_blast_spe.txt ./

# the lncRNAs of Apoly can be mapped to 34 fish species
# extract the lncRNAs of the fish species to detect the orthologous genes with Apoly
# jlkang@hpc2021 Thu Dec 28 22:49:51 /lustre1/g/sbs_schunter/Kang/lncRNAs
mkdir Orthlogous
cp rnacentral_species_specific_ids_lncRNAs.fasta Orthlogous/
# jlkang@hpc2021 Thu Dec 28 23:37:32 /lustre1/g/sbs_schunter/Kang/lncRNAs/Orthlogous
perl temp1.pl
# jlkang@hpc2021 Thu Dec 28 23:42:07 /lustre1/g/sbs_schunter/Kang/lncRNAs/Orthlogous
less ../Final_lincRNAs.fas|perl -alne 'if (/\>/){$nb++;my $id=Apoly.$nb;s/>//;print ">$id\t$_"}else{print}' > Apoly.fas
# Kang@fishlab3 Fri Dec 29 00:09:17 ~/Desktop/LncRNA
scp -r jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/lncRNAs/Orthlogous ./
# orf dectection
# Kang@fishlab3 Fri Dec 29 00:17:47 ~/Desktop/LncRNA/Orthlogous
ll *.fas|perl -alne '(my $spe)=$F[-1]=~/(.*)\.fas/;my $dir=$spe."_orf";my $cmd="TransDecoder.LongOrfs -t $F[-1] -O $dir";print $cmd' > 1.sh
sh 1.sh
# extract the longest_orfs.pep in the result folder for orthologous detection
# Kang@fishlab3 Fri Dec 29 00:27:48 ~/Desktop/LncRNA/Orthlogous
mkdir Total_peps
less Mapped_spe.txt|perl -alne 'my $dir=$F[-1]."_orf";my $nm=$F[-1].".fas";my $cmd="cp $dir/longest_orfs.pep Total_peps/$nm";system($cmd)'
cp Apoly_orf/longest_orfs.pep Total_peps/Apoly.fas
nohup orthofinder -f Total_peps -a 30 >orthofinder.process 2>&1 &
# [1] 20151
# Kang@fishlab3 Fri Dec 29 01:00:38 ~/Desktop/LncRNA/Orthlogous/Total_peps/OrthoFinder/Results_Dec29/Orthogroups
less Orthogroups.GeneCount.tsv|perl -alne 'next if /^Orthogroup/;print if $F[2]>0'|wc -l # 164 orthgroups include Apoly
less Orthogroups.GeneCount.tsv|perl -alne 'next if /^Orthogroup/;print if $F[2]>0'|perl -alne 'my $nb;for ($i=1;$i<@F-1;$i++){$nb++ if $F[$i]>=1};print $nb if $nb>=2'|perl -alne '$hash{$_}++;END{foreach my $key (sort keys %hash){my $nb=$key-1;print "$nb\t$hash{$key}"}}' > Orthogroup_shared.txt
less Orthogroup_shared.txt|perl -alne '$nb+=$F[1] if $F[0]==1;END{print $nb}' # 36 orthogroups were shared by Apoly and one another fish species
less Orthogroup_shared.txt|perl -alne '$nb+=$F[1] if $F[0]>=2 && $F[0]<=10;END{print $nb}' # 85 orthogroups were shared by Apoly and 2-10 another fish species
less Orthogroup_shared.txt|perl -alne '$nb+=$F[1] if $F[0]>10;END{print $nb}' # 43 orthogroups were shared by Apoly and more than 10 another fish species

# these orthologous might be not the lincRNAs
# Kang@fishlab3 Tue Jan 02 10:39:21 ~/Desktop/LncRNA
mkdir Orthlogous_classifier
# Kang@fishlab3 Tue Jan 02 11:06:14 ~/Desktop/LncRNA/Orthlogous_classifier
less Ame-1.0.2.110.gtf.gz|grep -i 'lncrna' > Ame_lncRNAs.gtf # Ame lncRNAs gtf
less Ame-1.0.2.110.gtf | grep -i 'protein_coding' > Ame_codingRNAs.gtf # Ame coding RNAs gtf
FEELnc_classifier.pl --lncrna=Ame_lncRNAs.gtf --mrna=Ame_codingRNAs.gtf > Ame.classifier.txt # lncRNAs classify
less Ame.classifier.txt|perl -alne 'if (/^isBest/){print}elsif($F[0]==1 && /intergenic/){print}'|perl -alne 'next if /^isBest/;print $F[1]'|sort -u|wc -l # 2289 lincRNAs
for i in *.gz;do echo -e ${i};less ${i}|grep -i 'lncrna'|wc -l;done
# Danio_rerio.GRCz11.110.gtf.gz, Eptatretus_burgeri.Eburgeri_3.2.110.gtf.gz, Lepisosteus_oculatus.LepOcu1.110.gtf.gz: no lncRNAs annotation, are lincRNAs annotation

less Cgo-3.1.110.gtf|grep -i 'lncrna' > Cgo_lncRNAs.gtf
less Cgo-3.1.110.gtf|grep -i 'protein_coding' > Cgo_codingRNAs.gtf
FEELnc_classifier.pl --lncrna=Cgo_lncRNAs.gtf --mrna=Cgo_codingRNAs.gtf > Cgo.classifier.txt # lncRNAs classify

less Ate-1.2.110.gtf|grep -i 'lncrna' > Ate_lncRNAs.gtf
less Ate-1.2.110.gtf|grep -i 'protein_coding' > Ate_codingRNAs.gtf
FEELnc_classifier.pl --lncrna=Ate_lncRNAs.gtf --mrna=Ate_codingRNAs.gtf > Ate.classifier.txt # lncRNAs classify

# Kang@fishlab3 Tue Jan 02 14:34:27 ~/Desktop/LncRNA/Orthlogous_classifier
perl temp1.pl > classify.sh

# it seems that we can not use ORFs to detect the orthologous lincRNAs
# Kang@fishlab3 Wed Jan 03 02:10:17 ~/Desktop/LncRNA/Orthlogous_classifier
cat *classifier.txt|perl -alne 'next if /^isBest/;print $F[1] if $F[0]==1'|sort -u > All_lincRNAs.txt # all lincRNAs across the 34 fish species
# Kang@fishlab3 Wed Jan 03 02:12:50 ~/Desktop/LncRNA/Orthlogous
cp ~/Desktop/LncRNA/Orthlogous_classifier/All_lincRNAs.txt ./
# Kang@fishlab3 Wed Jan 03 02:20:38 ~/Desktop/LncRNA/Orthlogous
mkdir lincRNAs_ortho
perl temp2.pl
# index for all these fish species
# makeblastdb -in Ame.fas -dbtype nucl -parse_seqids -out Ame
# Kang@fishlab3 Wed Jan 03 02:51:14 ~/Desktop/LncRNA/Orthlogous/lincRNAs_ortho
mkdir blastn_fromApoly; mkdir blastn_toApoly;
# reciprocal blast: Apoly blast to each species; each species blast to Apoly
# blastn -outfmt 6 -query Apoly.fas -out Apoly_Ame.bla -db Ame -evalue 1e-3 -num_threads 30
# Kang@fishlab3 Wed Jan 03 03:08:25 ~/Desktop/LncRNA/Orthlogous/lincRNAs_ortho
nohup perl temp2.pl &
# [1] 28588
# reciprocal best
# Kang@fishlab3 Wed Jan 03 03:15:25 ~/Desktop/LncRNA/Orthlogous/lincRNAs_ortho
mkdir blastn_rp_best
# Kang@fishlab3 Wed Jan 03 03:50:31 ~/Desktop/LncRNA/Orthlogous/lincRNAs_ortho/blastn_fromApoly
perl temp1.pl|perl -alne '$hash1{$F[0]}++;$hash2{$F[0]}.=$F[1].";";END{foreach my $key (sort keys %hash2){$hash2{$key}=~s/\;$//;print "$key\t$hash1{$key}\t$hash2{$key}"}}' > Final_lincRNAs_Ortho.txt
# 228 orthologous lincRNAs between Apoly and at least one another species
# Kang@fishlab3 Wed Jan 03 09:51:28 ~/Desktop/LncRNA/Orthlogous/lincRNAs_ortho
cp blastn_fromApoly/Final_lincRNAs_Ortho.txt ./
# Kang@fishlab3 Wed Jan 03 10:06:56 ~/Desktop/LncRNA/Orthlogous/lincRNAs_ortho
perl temp3.pl > Final_lincRNAs_Ortho_list.txt

# download protein sequences to check whether these ortholgous lincRNAs have same neighbour coding genes
# extract the longest transcript to represent the gene
# concatenate pep sequences
# Kang@fishlab3 Wed Jan 03 14:48:42 ~/Desktop/LncRNA/Orthlogous_classifier
mkdir Conca
# Kang@fishlab3 Wed Jan 03 15:01:39 ~/Desktop/LncRNA/Orthlogous_classifier/Conca
perl temp1.pl > all_34spe_pep_longest.all.fa
perl temp3.pl >Conca/all_34spe_pep.all.fa
cat 
# extract all neighbouring coding genes of the orthologous lincRNAs
# Kang@fishlab3 Wed Jan 03 15:17:37 ~/Desktop/LncRNA/Orthlogous_classifier
perl temp2.pl > Conca/lincRNAs_neighbourCoding.txt
less 1.txt |perl -alne 's/\r//g;my @a=split;print "$F[0]\t$F[-1]"' > Conca/Apoly_lincRNAs_neighbourCoding.txt
cat all_34spe_pep_longest.all.fa Apoly_pep_longest.all.fa > Final_pep_longest.all.fa
cat lincRNAs_neighbourCoding.txt Apoly_lincRNAs_neighbourCoding.txt > Final_lincRNAs_neighbourCoding.txt
# Annotate the pep sequences of the ortholgous lincRNAs have same neighbour coding genes
perl temp3.pl > Final_lincRNAs_Ortho_neighbour_coding_pep.fa
# Annotate based on uniprot
# diamond blastp -q makeblastdb_input.2.fa -e 1e-5 --sensitive -d ./blastdb --out blast_output.txt
# (base) kang1234@celia-PowerEdge-T640 Wed Jan 03 16:32:59 ~/swiss-prot
diamond blastp -q Final_lincRNAs_Ortho_neighbour_coding_pep.fa -e 1e-5 --sensitive -d uniprot-filtered-reviewed_yes.fasta --out Final_lincRNAs_Ortho_neighbour_coding_pep_blast_output.txt
# (base) kang1234@celia-PowerEdge-T640 Wed Jan 03 16:39:28 ~/swiss-prot
# (base) kang1234@celia-PowerEdge-T640 Wed Jan 03 16:40:15 ~/swiss-prot
perl extract_best_blastp.pl Final_lincRNAs_Ortho_neighbour_coding_pep_blast_output.txt|sort -u > Final_lincRNAs_Ortho_neighbour_coding_pep_blast_output_best_ano.txt
# Kang@fishlab3 Wed Jan 03 16:36:32 ~/Desktop/LncRNA/Orthlogous_classifier/Conca
scp kang1234@147.8.76.177:~/swiss-prot/Final_lincRNAs_Ortho_neighbour_coding_pep_blast_output_best_ano.txt ./
# Kang@fishlab3 Wed Jan 03 17:42:15 ~/Desktop/LncRNA/Orthlogous_classifier/Conca
perl temp4.pl > syntenic_results.txt

# about the annotation of coding genes in Apoly
```
