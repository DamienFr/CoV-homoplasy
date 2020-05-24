
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

perl ./scripts/annotated_vcf2WGS_perl.pl --vcf 04.1.filtrated_vcf \
--output_fasta_all_sites 07.filtrated_WGS_all_sites.fasta \
--output_fasta_only_variant 07.filtrated_WGS_only_variant.fasta \
--fasta_reference ./01.ref/refseq_sars.fasta \
--SNP_positions_to_exclude ./07.SNP_filtering_temp_files/05.1.low_depth_positions_uniq \
--SNP_positions_to_extract ./07.SNP_filtering_temp_files/05.2.variable_positions_uniq \
--SNP_stats_folder 08.SNP_filtering_stats \
--max_N 100








