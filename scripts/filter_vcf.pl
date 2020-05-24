#!/usr/bin/perl
use Getopt::Long;

=head1 DESCRIPTION

	This script annotates SNPs in a VCF file based on their quality (depth, phred quality, coverage, allele frequency and presence in SNP clusters).
	This script ONLY works with biallelic SNPS and 1 nucleotide long REF and ALT alleles.

=head1 USAGE

	perl filter_vcf.pl --vcf XXX.vcf --output_folder filtrated_vcf [--sliding 10] [--min_alt_allele_freq 0.8] [--min_depth 5] [--min_quality 20]

=head1 OPTIONS

	--vcf=s" => \$files,
	--output_folder_fasta_all_sites	Output folder for complete fasta sequence (eg. the length of the reference genome). Created if non existing. (REQUIRED)
	--sliding	Degfault : 10
	--min_alt_allele_freq	Default : 0.8
	--min_depth	Min total sequencing depth of the position in order to be considered for SNP Default : 5
	--min_quality	Min quality of the position in order to be considered for SNP Default : 20
	-h or --help:	This Documentation

=head1 AUTHOR

	Damien Richard, 2020

=cut

$sliding=10; 
$af=0.8 ;
$mincov=5;
$qual=20;

GetOptions(
	"vcf=s" => \$files,
	"output_folder=s" => \$folder_out,
	"sliding=s" => \$sliding,
	"min_alt_allele_freq=s" => \$af,	
	"min_depth=s" => \$mincov,
	"min_quality=s" => \$qual,
	"h|help" => \$hl
);

die `pod2text $0` unless $files && $folder_out ;
die `pod2text $0` if $hl;

$folder_out =~ s/^\.|\///g;
if (! -d $folder_out){ mkdir $folder_out }

$fileout = (split /\//, $files)[-1] ;
print STDERR "[INFO] Input file name $fileout \n";
$fileout =~ s/.vcf//g;

open (OUT, ">", "./" . $folder_out . "/" . $fileout . ".vcf");

open (VC, "<", $files) || die "Cannot open $files: $!";
while ($line4 = <VC>) {

	$line_number++;

	if($line4 =~ /^#/){
		if($line_number == 1){
			$line4 =~ s/\r?\n//; print OUT $line4 . "/annotated_using_vcf_filter.pl\n";
		}else{
			print OUT $line4 ;
		}
		if($line4 =~ /^#[^#]/){
			@F = split /\t/, $line4 ; $F[-1] =~ s/\r?\n//;
			if(length($F[-1]) == 0){
			print STDERR "[INFO] No Read group (optional) found. VCF produced by GATK should have read groups.";
			}else{
			print STDERR "[INFO] Read group (optional) is \"" . $F[-1] . "\"\n";
			}
		} # #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SRR11059940.#
	}else{

		@F = split /\t/, $line4 ;
	
		if(length($F[3]) == 1 && length($F[4]) == 1){ # if to only take into acount bi-allelic SNP (no multiallelic snp and no indels)
			@nb_fields_tmp = split /\:/ , $F[8] ; $nb_fields = scalar(@nb_fields_tmp);
			$nb_info_field{$nb_fields} = "1"; # i keep track of the number of info field so that i can print it out on STDERR and warn user if multiple values are found
			if( $nb_fields == 8 && (split /\:/,$F[8])[1] eq "DP" && (split /\:/,$F[8])[5] eq "AO" ){
			# GT:DP:AD:RO:QR:AO:QA:GL
				$dp = (split /:/,$F[9])[1]; # profondeur totale
				$oa = (split /:/,$F[9])[5]; # profondeur ALT
			}elsif( $nb_fields == 7 && (split /\:/,$F[8])[1] eq "DP" && (split /\:/,$F[8])[4] eq "AO" ){
			# GT:DP:RO:QR:AO:QA:GL
				$dp = (split /:/,$F[9])[1]; # profondeur totale
				$oa = (split /:/,$F[9])[4]; # profondeur ALT
			}elsif( $nb_fields == 5 && (split /\:/,$F[8])[1] eq "AD" && (split /\:/,$F[8])[2] eq "DP" ){
			#GT:AD:DP:GQ:PL	1/1:0,9:9:26:323,26,0
				$dp = (split /:/,$F[9])[2]; # profondeur totale
				$oa = (split /,/,(split /:/,$F[9])[1])[1] ; # profondeur ALT 
				print STDERR $F[1] . "\t" . $oa . "\n";
			}
			else{die "cannot read properly the info field of sample " . $fileout . ", will die now\n" }
		
			$F[6] = "";
			if($dp < $mincov){ $F[6] .= "Low_depth;" } 
			if( $F[5] <= $qual){ $F[6] .= "Bad_mapping;" }
			if(($oa / $dp) <= $af ){ $F[6] .= "Heteroplastic;" }
		
			if( $F[6] ne ""){
			$SNP{$F[1]} = [@F];
			}else{
		#if($last_pass && ($F[1] - $last_pass) < $sliding ){$F[6] .= "SnpCluster;"; ${$SNP{$last_pass}}[6] .= "SnpCluster;" }
		$SNP{$F[1]} = [@F];
		$last_pass = $F[1];
		}
		} # end of # if to only take into acount bi-allelic SNP 
} # fin du premier else

} # end of reading vcf file

foreach $nb_field (sort { $a <=> $b } keys %nb_info_field){
	print STDERR "[INFO] VCF file has $nb_field info fields\n";
}

$nb_of_different_info_fields = keys %nb_info_field;
if($nb_of_different_info_fields > 1){  print STDERR "WARNING\nYou should consider revising your VCF, as the number of info fields appears variable along the file\n\n"}


$count_deleted = 0;
foreach $position (sort { $a <=> $b } keys %SNP){
	print OUT join("\t",@{$SNP{$position}});
	if( ${$SNP{$position}}[6] eq "" ){$count_normal++}else{$count_normal++; $count_deleted++; }
}
print STDERR "[INFO] " . $count_deleted . " out of " . $count_normal . " SNPs were annotated as low quality\n\n" ;

close OUT;



