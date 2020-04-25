
###########################################################################################################
######## Mapping of the raw reads on Wuhan-hu-1 reference genome (GenBank NC_045512.2 ) ###################
###########################################################################################################

# SNP calling pipeline, to be done on a HPC cluster
mkdir 02.results_bwa
mkdir 023.depths
mkdir 07.SNP_filtering_temp_files

for files in ./01.sra_fastq/*.fastq
do
fileout=$(echo $files | perl -F'\/' -ane '$F[-1] =~ s/fastq/sam/; print $F[-1]')
echo $fileout

# mapping
bwa mem -M -t 6 refseq_sars.fasta $files > ./02.results_bwa/${fileout}

# sorting bam file (needed in order to get MarkDuplicates to work)
samtools sort $files -o ./02.results_bwa/${fileout}sorted.sam

# removing of PCR duplicates in order not to study a illegit high sequencing depth
java -jar /media/bacterio/DATA2/DATA/picard2.22/picard.jar MarkDuplicates \
      I=./02.results_bwa/${fileout}sorted.sam \
      O=./02.results_bwa/${fileout}rmdup.bam \
      VALIDATION_STRINGENCY=LENIENT \
      REMOVE_DUPLICATES=true \
      M=./02.results_bwa/${fileout}.delme.txt

rm -f ./02.results_bwa/${fileout}.delme.txt

samtools index ./02.results_bwa/${fileout}rmdup.bam

# For all sample, I determine the depth of sequencing for all the nucleotide positions of the ref 
samtools depth -aa ./02.results_bwa/${fileout}rmdup.bam > ./023.depths/${fileout}

# i will now identify and not study further the SRA projects that led to a sequencing 
# depth of less than 5 on more than 1% of the length of the reference genome
env fileout=$fileout perl -F'\t' -ane 'BEGIN{use List::Util qw(sum);
$name = $ENV{fileout}}
 $F[2] =~ s/\r?\n//;
 push @values, $F[2];
if($F[2] > 5){ $a++ }
END{
if(@values == 0){$average = 0 }else{$average = sum(@values) / @values;}  # calcul de la moyenne
$pcent_sup_5 = $a * 100 / 29903;
 print "$name\t$average\t$pcent_sup_5\n" } ' ./023.depths/${fileout} >> 023.depth.all_indivs

# keeping track of the positions of the reference genome not covered for at least one individual 
perl -F'\t' -ane '$F[-1] =~ s/\r?\n//; if($F[-1] == 0 ){ print $F[1] . "\n" } ' ./023.depths/${fileout} >> 05.1.low_depth_positions
sort -n 05.1.low_depth_positions | uniq > ./07.SNP_filtering_temp_files/05.1.low_depth_positions_uniq

done

rm -f 05.1.low_depth_positions


mkdir 02.results_bwa_selected
# proper copy of the interesting files in 02.results_bwa_selected
perl -F'\t' -ane 'BEGIN{use File::Copy qw(copy);}; if($F[-1] > 99){
$old_file = "./02.results_bwa/" . $F[0] . "rmdup.bam";
$new_file = "./02.results_bwa_selected/" . $F[0] . "rmdup.bam";
copy  $old_file, $new_file}else{ $count_deleted++ };END{print STDERR $count_deleted . " deleted strains\n"}' 023.depth.all_indivs

# it's possible to keep only the mapped reads of each bam file to save space if needed
# for files in ./02.results_bwa_selected/*; do
# fileout=$(echo $files | perl -F'\/' -ane 'print $F[-1]')
# samtools view -b -F 4 $files > ./02.results_bwa_selected_mapped/${fileout} ; done


###########################################################################################################
####################  Retrieving the biosample metadata for the selected SRA projects ###################
###########################################################################################################

# i restrict the SRA metadata table to the kept projects :
perl -F'\t' -ane 'if($. ==1 ){ print }elsif( -f "./02.results_bwa_selected/" . $F[0] . ".rmdup.bam"){print}' 00.sra_results.tab > 00.sra_and_biosample_results_selected_indivs.tab

# Extract metadata from biosample # production of a list of the biosample metadata to retrieve & delete  the header line
cut -f20 -d$'\t' 00.sra_and_biosample_results_selected_indivs.tab | sed '1d' > 00.biosamples

# method to retrieve large numbers of biosample metadata even though NCBI is not happy about such large querries.
mkdir 00.temp_biosamples
split -l 100 00.biosamples ./00.temp_biosamples/00.biosamples_

for block in `ls ./00.temp_biosamples/00.biosamples_*`; 
  do 
  tr "\n" "," < $block | sed 's/,$//' | sed 's/^/\nefetch -db biosample -mode xml -id /'; 
done > 00.get_biosamples_info.sh

bash 00.get_biosamples_info.sh > 00.biosample.xml
perl -ne 's/</\n</g; print' 00.biosample.xml > 00.biosample-reformat.xml
rm -rf 00.temp_biosamples ; rm -f 00.biosample.xml

perl 00.biosample_xml_to_tab.pl > 00.biosample_results.tab

# Addition of the biosample information to the SRA metadata table
perl -F'\t' -ane 'BEGIN{open (GB, "<", "00.biosample_results.tab");while ($line = <GB>) {
$line =~ s/\r?\n//;
@fields = split /\t/,$line;
$h{$fields[1]} = $fields[2] . "\t" . $fields[3] . "\t" . $fields[4] }
close GB};
s/\r?\n//;
$F[-1] =~ s/\r?\n//;
if($. == 1){print $_ . "\tcollection_date\tgeo_loc_name\tisolation_source\n"}
else{
print $_ . "\t" . $h{$F[19]} . "\n" } ' 00.sra_results.tab > 00.sra_and_biosample_results.tab

###########################################################################################################
######################################  SNP Calling and filtering  ########################################
###########################################################################################################

# to be done on a HPC cluster
mkdir 03.results_vcf
rm -rf 04.1.filtrated_vcf

for files in ./02.results_bwa_selected/*
do
fileout=$(echo $files | perl -F'\/' -ane '$F[-1] =~ s/rmdup.bam//; print $F[-1]')
freebayes --min-alternate-fraction 0.2 --ploidy 1 --no-indels --haplotype-length 1 --no-complex --no-mnps --min-alternate-count 2 -f ../00.ref/refseq_sars.fasta --bam ${files} > ./03.results_vcf/${fileout}vcf
perl ./scripts/filter_vcf.pl --vcf ./03.results_vcf/${fileout}vcf --output_folder 04.1.filtrated_vcf --sliding 0 --min_alt_allele_freq 0.65 --min_depth 1 --min_quality 20
done

# keeping all positions outputed in ALL vcf (will be used to output the "variant only" wgs fasta file). These positions include positions that WILL be filtered OUT.
perl -F'\t' -ane 'if(! /^#/){ print $F[1] . "\n" }' ./04.1.filtrated_vcf/*.vcf > 05.2.variable_positions
perl -ne 's/\r?\n//; $h{$_} = "ok";END{ foreach $key (sort { $a <=> $b } keys %h){print $key . "\n"} }' 05.2.variable_positions > ./07.SNP_filtering_temp_files/05.2.variable_positions_uniq
rm -f 05.2.variable_positions

###########################################################################################################
###################  Creation of aligned fasta SNP Calling and filtering  #################################
###########################################################################################################

for files in ./04.1.filtrated_vcf/*.vcf
do
perl ./scripts/annotated_vcf2WGS_perl.pl --vcf ${files} \
--output_folder_fasta_all_sites 07.filtrated_WGS_all_sites \
--output_folder_fasta_only_variant 07.filtrated_WGS_only_variant \
--fasta_reference ./01.ref/refseq_sars.fasta \
--low_depth_positions ./07.SNP_filtering_temp_files/05.1.low_depth_positions_uniq \
--SNP_positions_to_extract ./07.SNP_filtering_temp_files/05.2.variable_positions_uniq
done

# i will create the two multifasta with N on positions
cat ./07.filtrated_WGS_all_sites/* > 07.filtrated_WGS_all_sites.fasta
cat ./07.filtrated_WGS_only_variant/* > 07.filtrated_WGS_only_variant.fasta








