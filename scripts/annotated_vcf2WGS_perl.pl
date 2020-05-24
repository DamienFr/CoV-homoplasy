#!/usr/bin/perl
use Getopt::Long;
use Cwd qw(cwd);

=head1 DESCRIPTION

	This script produces a fasta file from a VCF and a reference genome from which the SNPs have been called.
	ONLY works with annotated VCF files originating from vcf_filter.pl. This is because SNPs have already been annotated based on their quality (depth, phred quality, coverage, allele frequency and presence in SNP clusters).
	As vcf_filter.pl, this script ONLY works with biallelic SNPS and 1 nucleotide long REf and ALT alleles.

=head1 USAGE

	perl annotated_vcf2WGS.pl --vcf XXX.vcf --output_fasta_all_sites all_sites --output_fasta_only_variant only_variant --fasta_reference ref.fasta --SNP_positions_to_exclude low_depth.txt --SNP_positions_to_extract variable_positions.txt --outgroup strain01 --max_N 2 --SNP_stats_folder snp_metrics

=head1 OPTIONS

	--vcf=s" => Folder containing vcf files to process (all .vcf files will be processed) (REQUIRED)
	--output_fasta_all_sites	Output for complete fasta sequence (eg. the length of the reference genome) (REQUIRED)
	--output_fasta_only_variant	Output for "variant-only" fasta sequence (eg. the length of the nb of SNPs in SNP_positions_to_extract file). (REQUIRED)
	--fasta_reference	Reference genome on whoch reads have been mapped to produce the VCFs. Multifasta files will make the script die (REQUIRED)
	--SNP_positions_to_exclude	File containing the list of positions to exclude. Usefull to specify low depth positions and/or positions not covered and/or recombinant positions. These positions will be turned to "N" if and only if they are part of the variable positions.
	--SNP_positions_to_extract	File containing the list of variable SNP positions. Only those will be considered.  Usefull to not take into account sites detecetd as recombinant by ClonalFrameML for ex. (REQUIRED)
	--outgroup	Outgroup. That's it. It's the outgroup. Not more not less. (OPTIONAL)
	--max_N	Maximum number of N per nucleotide position across all individuals. N in the sequence of the outgroup DO NOT count. Positions with more N will be deleted for the "only variant" fasta file and converted to "N" for all individuals for the "all sites" fasta file. Surviving variable positions will be outputed. 
	--SNP_stats_folder	Folder to which multiple files providing SNP number metrics will be written. Created if non-existant (REQUIRED)
	-h or --help:	This Documentation

=head1 AUTHOR

	Damien Richard, 2020

=cut

GetOptions(
	"vcf=s" => \$vcf_folder,
	"output_fasta_all_sites=s" => \$out_all_sites,
	"output_fasta_only_variant=s" => \$out_variant_only,
	"fasta_reference=s" => \$ref,	
	"SNP_positions_to_exclude=s" => \$pos_to_exclude,
	"SNP_positions_to_extract=s" => \$var_positions,
	"SNP_stats_folder=s" => \$stats_folder,
	"outgroup=s" => \$outgroup,
	"max_N=s" => \$max_N,
	"h|help" => \$hl,
    "debug" => \$debug
);

die `pod2text $0` unless $vcf_folder && $out_all_sites && $out_variant_only && $ref && $pos_to_exclude && $var_positions && $stats_folder && defined($max_N) ;
die `pod2text $0` if $hl;

####################################################################

open (VAR, "<", $var_positions) || die "Cannot open $var_positions: $!";
while ($line2 = <VAR>) {
$line2 =~ s/\r?\n//;
$pos_var{$line2} = "1";
}
close VAR;

open (EXC, "<", $pos_to_exclude) || die "Cannot open $pos_to_exclude: $!";
while ($line3 = <EXC>) {
$line3 =~ s/\r?\n//;
$h_pos_exclude{$line3} = "1";
$h_N_all_sites{$line3} = $max_N+1; # To ensure that positions to exclude are outputed as "N" and do not appear in variable position files i need them to be excluded based on their number of N (eg. all indiv)
}
close EXC;

open (GB, "<", $ref) || die "Cannot open $ref: $!";
while ($line = <GB>) {
if($line =~ /^>/){ $a++; 
if($a > 1 ){die "Script killed : cannot work on multi-molecule fasta file\n";}}else{
$line =~ s/\s|\r?\n//g; $seq_orig .= $line;
}}
close GB;


# building of the shortref, already including "N" at the positions of low coverage for at least one individual (not modificable behaviour at the time)
# at the same time building the table that keeps track of the nucleotide position of the first, second etc ... artificial bases of the fake shortref: $position_table
foreach $key (sort { $a <=> $b } keys %pos_var){
if(!exists($h_pos_exclude{$key})){
$short_ref_orig .= substr $seq_orig, ($key-1), 1;
}else{
$short_ref_orig.="N" # a voir si je laisse ca
};
$k++; $position_table{$key}=$k;
}

$unfiltrated_total_SNP = keys %pos_var;

# version without the possibility of using positions to exclude
# which maybe should be implemented because in order to exclude positions we just have to exclude them from varpos ... do we ? 
# foreach $key (sort { $a <=> $b } keys %pos_var){
# $short_ref_orig .= substr $seq_orig, ($key-1), 1;
# $k++; $position_table{$key}=$k;
# }

$dir = cwd;
$vcf_folder =~ s/^\.\/|^\/|\/$//g;
$stats_folder =~ s/^\.\/|^\/|\/$//g;

mkdir($stats_folder);

opendir $dh, $vcf_folder;
print STDERR "Opening folder " . $vcf_folder . " from current directory $dir \n";
@vcf_files = readdir $dh;

foreach $files (@vcf_files)
{
    # skip . and ..
    if ($files eq '.' or $files eq '..') {
        next;
    }

if($debug){ print STDERR "Working on $files\n"; }

$short_ref = $short_ref_orig;
$seq = $seq_orig;
$line_number=0;

open ($vc, "<", $dir . "/" . $vcf_folder . "/" . $files) || die "Cannot open $files: $!";
$files =~ s/\.vcf// ;
while ($line4 = <$vc>) {
$line_number++;
# checking that the file has been annotated using vcf_filter.pl
if($line_number == 1){ 
if( $line4 !~ /.*annotated_using_vcf_filter.pl$/ ){
die "$files This programm ONLY works with annotated VCF files originating from vcf_filter.pl.\n This is because SNPs have already been annotated based on their quality (depth, phred quality, coverage, allele frequency and presence in SNP clusters.\n\n";
}
}

@F = split /\t/, $line4 ;
;
if($line4 !~ /^#/ && exists($pos_var{$F[1]}) && $line4 !~ /^$/)
{ # skip header # la condition que la position soit indiquee comme variable dans %pos_var est pas obligatoire pour le premier passage du programme (avant clonalframe, mais l est pour le deuxieme
# car c est cela qui permet de virer les positions recombinantes
$h_nb_molecules{$F[0]} = "pouet"; if(scalar %h_nb_molecules > 1){die "Script killed : Vcf issue: multiple molecules are present in the vcf or different references were used to generate the vcf of the input folder\n"}

if($files eq $outgroup){ $h_var_pos_only_outgroup{$F[1]} = "ok"}else{ $h_var_pos_no_outgroup{$F[1]} = "ok" }

if($F[6] eq ""){
substr($seq,$F[1]-1,1) = $F[4];
substr($short_ref,$position_table{$F[1]}-1,1) = $F[4];
 }elsif($F[6] eq "Heteroplastic;"){
 substr($seq,$F[1]-1,1) = "N"; # "H"
 substr($short_ref,$position_table{$F[1]}-1,1) = "N"; # "H"
 if($files ne $outgroup ){ $h_N_all_sites{$F[1]}++ ; $h_N_variant_only{$position_table{$F[1]}}++ ; }else{$h_pos_outgroup_N{$F[1]}++ }
  }elsif($F[6] eq "Low_depth;"){
   substr($seq,$F[1]-1,1) = "N"; # "V"
   substr($short_ref,$position_table{$F[1]}-1,1) = "N"; # "V"
   if($files ne $outgroup ){ $h_N_all_sites{$F[1]}++ ; $h_N_variant_only{$position_table{$F[1]}}++ }else{$h_pos_outgroup_N{$F[1]}++ }
  }elsif($F[6] eq "Bad_mapping;"){
     substr($seq,$F[1]-1,1) = "N"; # "M"
     substr($short_ref,$position_table{$F[1]}-1,1) = "N";  # "M"
     if($files ne $outgroup ){ $h_N_all_sites{$F[1]}++ ; $h_N_variant_only{$position_table{$F[1]}}++ }else{$h_pos_outgroup_N{$F[1]}++ }
    }elsif($F[6] eq "SnpCluster;"){
         substr($seq,$F[1]-1,1) = "N"; # "S"
         substr($short_ref,$position_table{$F[1]}-1,1) = "N"; # "S"
         if($files ne $outgroup ){ $h_N_all_sites{$F[1]}++ ; $h_N_variant_only{$position_table{$F[1]}}++ }else{$h_pos_outgroup_N{$F[1]}++ }
      }elsif(exists($h_pos_exclude{$F[1]})){
          substr($seq,$F[1]-1,1) = "N"; 
          substr($short_ref,$position_table{$F[1]}-1,1) = "N";
          if($files ne $outgroup ){ $h_N_all_sites{$F[1]}++ ; $h_N_variant_only{$position_table{$F[1]}}++ }else{ $h_pos_outgroup_N{$F[1]}++ }
       }else{
           substr($seq,$F[1]-1,1) = "N"; # "B";
           substr($short_ref,$position_table{$F[1]}-1,1) = "N"; # "B"
           if($files ne $outgroup ){ $h_N_all_sites{$F[1]}++ ; $h_N_variant_only{$position_table{$F[1]}}++ }else{$h_pos_outgroup_N{$F[1]}++ }
           }
}

# end of individual loop

} # end of files line reading

$h_result_short{$files} = $short_ref ;
$h_result_long{$files} = $seq ;


close $vc;
} # end of foreach file

closedir $dh;

# i delete from the hashes the positions that do NOT have too many "N"s. so that later i delete positions that remain in the hashes
for $k (sort { $a <=> $b } keys %h_N_variant_only){
	if($h_N_variant_only{$k} <= $max_N ){
		delete $h_N_variant_only{$k};
	}
}

for $k (sort { $a <=> $b } keys %h_N_all_sites){
if($debug){ print STDERR "$h_N_all_sites{$k} $k \n" }
	if($h_N_all_sites{$k} <= $max_N ){
		delete $h_N_all_sites{$k};
	}else{
		delete($pos_var{$k}) ; delete($h_var_pos_no_outgroup{$k}); # i delete the record of the position as being variable (both in total count and in total without outgroup counts) because it will be replaced by N for all individual
		delete($h_var_pos_only_outgroup{$k}); # including the outgroup so i delete also the record for the outgroup
	}
}

# for all variable positions where the outgroup is N and other strains do not have a SNP, i will
# * delete the position for "only variant" file. To do that i add the position in h_N_variant_only hash
# * do nothing for the "all sites" file, because outgroup will be N and others will have same genotype, that's okay.
foreach $pos_outgroup_N (keys %h_pos_outgroup_N){
if( !exists($h_var_pos_no_outgroup{$pos_outgroup_N}) ){
delete $h_var_pos_only_outgroup{$pos_outgroup_N};
delete $pos_var{$pos_outgroup_N};
$h_N_variant_only{$position_table{$pos_outgroup_N}} = "1";

}}



open (OUT, ">", "./" . $stats_folder . "/outgroup_variable_positions.txt");
print OUT $_ . "\n" for sort { $a <=> $b } keys %h_var_pos_only_outgroup;
close OUT;

$v_var_pos_only_outgroup = keys %h_var_pos_only_outgroup;

open (NOUT, ">",  "./" . $stats_folder . "/no_outgroup_variable_positions.txt");
print NOUT $_ . "\n" and delete $h_var_pos_only_outgroup{$_} for sort { $a <=> $b } keys %h_var_pos_no_outgroup;
close NOUT;

$v_var_pos_no_outgroup = keys %h_var_pos_no_outgroup;

open (TOT, ">",  "./" . $stats_folder . "/total_variable_positions.txt");
print TOT $_ . "\n" for sort { $a <=> $b } keys %pos_var;
close TOT;

$v_pos_var = keys %pos_var;

open (OUTP, ">",  "./" . $stats_folder . "/outgroup_private_variable_positions.txt");
print OUTP $_ . "\n"  for sort { $a <=> $b } keys %h_var_pos_only_outgroup;
close OUTP;

$v_var_pos_private_outgroup = keys %h_var_pos_only_outgroup;


print STDERR "\nUnfiltrated :\nTotal SNP (with outgroup ): $unfiltrated_total_SNP \n\nFiltrated :\nTotal filtrated SNP : $v_pos_var\nSNP without outgroup : $v_var_pos_no_outgroup\nSNP of the outgroup : $v_var_pos_only_outgroup\nSNP private to the outgroup : $v_var_pos_private_outgroup\n";
open (GLO, ">",  "./" . $stats_folder . "/global_stats.txt");
print GLO "\nUnfiltrated :\nTotal SNP (with outgroup ): $unfiltrated_total_SNP \n\nFiltrated :\nTotal filtrated SNP : $v_pos_var\nSNP without outgroup : $v_var_pos_no_outgroup\nSNP of the outgroup : $v_var_pos_only_outgroup\nSNP private to the outgroup : $v_var_pos_private_outgroup\n";
close GLO;

open (ALL, ">", $out_all_sites);
open (VAR, ">", $out_variant_only);

if($debug){ $tmppp = keys %position_table;
print STDERR "Before deletion of filtered SNPs, position_table has $tmppp entries\n";} # 3643

%position_to_del_shrt_aln = %h_N_variant_only ;

# keys of position_table are real genome positions from the vcf
foreach $input_position ( keys %position_table){ 
if( !exists($pos_var{$input_position}) ){  $position_to_del_shrt_aln{$position_table{$input_position}} = "to_delete" ; delete($position_table{$input_position}) }
}

# new :

if($debug){ $tmppp = keys %position_to_del_shrt_aln;
print STDERR "Deletion table has $tmppp entries\n";} # 20
# end new

if($debug){$tmppp = keys %position_table;
print STDERR "After deletion of filtered SNPs, position_table has $tmppp entries\n";} # 3147

foreach $key ( keys %h_result_short){
# new :
	$decalage = 1;
	foreach $position (sort { $a <=> $b } keys %position_to_del_shrt_aln){
		substr( $h_result_short{$key} , ($position - $decalage), 1 ) =~ s/.// ; $decalage++; 
	}


# old :
# substr( $h_result_short{$key} , ($_ - 1), 1 ) =~ s/./N/ for keys %h_N_variant_only;
# we can't delete positions now because as soon as you delete a position, the other ones get wrong (-1 each time)
# substr( $h_result_short{$key} , ($_ - 1), 1 ) =~ s/.// for keys %h_N_variant_only;

print VAR ">" . $key . "\n" . $h_result_short{$key} . "\n";

substr( $h_result_long{$key} , ($_ - 1), 1 ) =~ s/./N/ for keys %h_N_all_sites;
print ALL ">" . $key . "\n" . $h_result_long{$key} . "\n";
}

close VAR;
close ALL;


