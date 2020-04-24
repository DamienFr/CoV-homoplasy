#!/usr/bin/perl
use Getopt::Long;

=head1 DESCRIPTION

	This script produces a fasta file from a VCF and a reference genome from which the SNPs have been called.
	ONLY works with annotated VCF files originating from vcf_filter.pl. This is because SNPs have already been annotated based on their quality (depth, phred quality, coverage, allele frequency and presence in SNP clusters).
	As vcf_filter.pl, this script ONLY works with biallelic SNPS and 1 nucleotide long REf and ALT alleles.

=head1 USAGE

	perl annotated_vcf2WGS.pl --vcf XXX.vcf --output_folder_fasta_all_sites all_sites --output_folder_fasta_only_variant only_variant --fasta_reference ref.fasta --low_depth_positions low_depth.txt --SNP_positions_to_extract variable_positions.txt

=head1 OPTIONS

	--vcf=s" => \$files,
	--output_folder_fasta_all_sites	Output folder for complete fasta sequence (eg. the length of the reference genome). Created if non existing. (REQUIRED)
	--output_folder_fasta_only_variant	Output folder for "variant-only" fasta sequence (eg. the length of the nb of SNPs in SNP_positions_to_extract file). Created if non existing. (REQUIRED)
	--fasta_reference	Reference genome on whoch reads have been mapped to produce the VCFs. Multifasta files will make the script die (REQUIRED)
	--low_depth_positions	File containing the list of positions not covered by at least one individual. These positions will NOT be considered (REQUIRED)
	--SNP_positions_to_extract	File containing the list of variable SNP positions. Only those will be considered.  Usefull to not take into account sites detecetd as recombinant by ClonalFrameML for ex. (REQUIRED)
	-h or --help:	This Documentation

=head1 AUTHOR

	Damien Richard, 2020

=cut

GetOptions(
	"vcf=s" => \$files,
	"output_folder_fasta_all_sites=s" => \$folder_out,
	"output_folder_fasta_only_variant=s" => \$folder_out_variant_only,
	"fasta_reference=s" => \$ref,	
	"low_depth_positions=s" => \$low_depth,
	"SNP_positions_to_extract=s" => \$var_positions,
	"h|help" => \$hl
);

die `pod2text $0` unless $files && $folder_out && $folder_out_variant_only && $ref && $low_depth && $var_positions ;
die `pod2text $0` if $hl;

$folder_out =~ s/^\.|\///g;
if (! -d $folder_out){ mkdir $folder_out }

$folder_out_variant_only =~ s/^\.|\///g;
if (! -d $folder_out_variant_only){ mkdir $folder_out_variant_only }


$fileout = (split /\//, $files)[-1] ;
$fileout =~ s/.vcf//g;


####################################################################


open (VAR, "<", $var_positions) || die "Cannot open $var_positions: $!";
while ($line2 = <VAR>) {
$line2 =~ s/\r?\n//;
$pos_var{$line2} = "1";
}
close VAR;

open (LOW, "<", $low_depth) || die "Cannot open $low_depth: $!";
while ($line3 = <LOW>) {
$line3 =~ s/\r?\n//;
$pos_low{$line3} = "1";
}
close LOW;

open (GB, "<", $ref) || die "Cannot open $ref: $!";
while ($line = <GB>) {
if($line =~ /^>/){ $a++; 
if($a > 1 ){die "Script killed : cannot work on multi-molecule fasta file\n";}}else{
$line =~ s/\s|\r?\n//g; $seq .= $line;
}}
close GB;

# building of the shortref, already including "N" at the positions of low coverage for at least one individual (not modificable behaviour at the time)
# at the same time building the table that keeps track of the nucleotide position of the first, second etc ... artificial bases of the fake shortref: $position_table
foreach $key (sort { $a <=> $b } keys %pos_var){
$logg++;
if(!exists($pos_low{$key})){
$short_ref.= substr $seq, ($key-1), 1;
}else{
$short_ref.="N"
};
$k++; $position_table{$key}=$k;
}

open (VC, "<", $files) || die "Cannot open $files: $!";
while ($line4 = <VC>) {

$line_number++;
# checking that the file has been annotated using vcf_filter.pl
if($line_number == 1){ 
if( $line4 !~ /.*\/annotated_using_vcf_filter.pl$/ ){
die "This programm ONLY works with annotated VCF files originating from vcf_filter.pl.\n This is because SNPs have already been annotated based on their quality (depth, phred quality, coverage, allele frequency and presence in SNP clusters.\n\n";
}
 }

@F = split /\t/, $line4 ;

if(!/^#/ && exists($pos_var{$F[1]})){ # skip header # la condition que la position soit indiquee comme variable dans %pos_var est pas obligatoire pour le premier passage du programme (avant clonalframe, mais l est pour le deuxieme
# car c est cela qui permet de virer les positions recombinantes
$h{$F[0]} = "pouet"; if(scalar %h > 1){die "Script killed : Vcf issue: multiple molecules are present in the vcf\n"}

if($F[6] eq ""){
substr($seq,$F[1]-1,1) = $F[4];
substr($short_ref,$position_table{$F[1]}-1,1) = $F[4];
 }elsif($F[6] eq "Heteroplastic;"){
 substr($seq,$F[1]-1,1) = "N";# "H";
substr($short_ref,$position_table{$F[1]}-1,1) = "N";# "H";
  }elsif($F[6] eq "Low_depth;"){
   substr($seq,$F[1]-1,1) = "N";# "V";
substr($short_ref,$position_table{$F[1]}-1,1) = "N";# "V";
  }elsif($F[6] eq "Bad_mapping;"){
     substr($seq,$F[1]-1,1) = "N";# "M";
substr($short_ref,$position_table{$F[1]}-1,1) = "N";# "M";
    }elsif($F[6] eq "SnpCluster;"){
         substr($seq,$F[1]-1,1) = "N";# "S";
substr($short_ref,$position_table{$F[1]}-1,1) = "N";# "S";
      }elsif(exists($hpos_low{$F[1]})){
          substr($seq,$F[1]-1,1) = "N"; # il est deja egal a N ... voir ou je le laisse
substr($short_ref,$position_table{$F[1]}-1,1) = "N";     
       }else{
substr($seq,$F[1]-1,1) = "N";# "B";
substr($short_ref,$position_table{$F[1]}-1,1) = "N";}# "B";}
};

}
close VC;

open (OUT, ">", "./" . $folder_out_variant_only . "/" . $fileout . "_simplified.fasta");
print OUT ">" . $fileout . "\n" . $short_ref . "\n";
close OUT;


open (OUTALL, ">", "./" . $folder_out . "/" . $fileout . ".fasta");
print OUTALL ">" . $fileout . "\n" . $seq . "\n";
close OUTALL;




























