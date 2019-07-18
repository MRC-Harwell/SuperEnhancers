#!/usr/bin/perl
# EPU comparison analysis / TADs comparison analysis
#Author: Sid

use Getopt::Long;
use Cwd ;

my ($path,$help,$out,$in) = "" ;
GetOptions(
    'help'   => \$help,
    'o=s' => \$out ,
    'i=s' => \$in ,
    'p=s' => \$path ,

) or die "\n**********  Incorrect usage!  ***********\nrun with -help option to see the useage\n ";

sub useage { die(qq/
	USAGE : perl <script> <arguments>
	ARGUMENTS :
                    REQUIRED
                    -i -> input file
                    -o -> name of the output files\/directory
                    OPTIONAL
                    -p -> OUTPUT path where the files need to be saved [default = cwd ]
                    -help -> prints this help message
               \n/);
}

if($help) { &useage ;}
#if() { print "\n MISSING ARGUMENTS : Give all the required options\n" ; &useage ;}
#unless( ) { print "\n \n" ; &useage ;}
if( !$path || $path =~ /\./ ) { $cwd = getcwd ; $path = "$cwd" ; }


# epu file
my $epu = "/NGS/users/Sid/ENHANCERS/SuperEnh_27ac/Great_targets/Strict_targets/Epus/Data/epu_cord.txt" ;
# tad file
#my $epu = "/NGS/users/Sid/ENHANCERS/SuperEnh_27ac/Great_targets/Strict_targets/Epus/Data_tads/cortex_tads.txt" ;
#my $epu = "/NGS/users/Sid/ENHANCERS/SuperEnh_27ac/Great_targets/Strict_targets/Epus/Data_tads/mESC_tads.txt" ;

# enhancer file for extracting enhancer ids
my $enh = "/NGS/users/Sid/ENHANCERS/SuperEnh_27ac/all.bed" ;
my %enhIds ;
open (enh, $enh)or die "Cannot open file $enh \n";
foreach(<enh>) {
	$_ =~ s/\n|\r|\s+$//g ;
	my($chr,$start,$end,$id) = split("\t",$_) ;
	$enhIds{"$chr-$start-$end"} = $id ;
}
close enh ;


my %geneCordinates ;
open (go,"<gene.cord.txt")or die "Cannot open file gene.cord.txt \n";
foreach(<go>) {
	$_ =~ s/\n|\r|\s+$//g ;
	my($chr,$start,$end,$name) = split("\t",$_) ;
	$geneCordinates{$name} = "$chr:$start-$end" ;
}
close go ;



my %genEpu ; my %EpuData ; 
open (o1,"<gene.epu.intersect.txt")or die "Cannot open file gene.epu.intersect.txt\n";
foreach(<o1>) {
	$_ =~ s/\n|\r|\s+$//g ;
	my($chr,$start,$end,$name,$ec,$es,$ee,$eid,$tmp) = split("\t",$_) ;
	$genEpu{$name} = $eid ;
	$EpuData{$eid} = "$ec:$es-$ee" ;
}
close o1 ;

my %enhEpu ; my %enhData ;
open (o2,"<enh.epu.intersect.txt")or die "Cannot open file enh.epu.intersect.txt\n";
foreach(<o2>) {
	$_ =~ s/\n|\r|\s+$//g ;
	my($chr,$start,$end,$name,$ec,$es,$ee,$eid) = split("\t",$_) ;
	$enhEpu{$name} = $eid ;
	$enhData{$eid} = "$ec:$es-$ee" ;
}
close o2 ;


open (in, $in)or die "Cannot open file $in \n";
open(out,">RESULTS.txt") ;
foreach(<in>) {
  $_ =~ s/\n|\r|\s+$//g ;
  my($chr,$start,$end,$name,$dist) = split("\t",$_) ;
  $gCordinates = $geneCordinates{$name} ;

  my $id ;
	if(exists $enhIds{"$chr-$start-$end"}){
		$id = $enhIds{"$chr-$start-$end"} ;
	}
	else { $id = "missing_id" ;}

  my $enh_epu ; my $enh_epu_cord ;
	if(exists $enhEpu{$name}){
		$enh_epu = $enhEpu{$name} ;
		$enh_epu_cord = $enhData{$enh_epu} ;
	}
	else { $enh_epu = "No_overlap_enh" ; $enh_epu_cord = "No_overlap_enh" ;}

	my $gene_epu ; my $gene_epu_cord ;
	if(exists $genEpu{$name}){
		$gene_epu = $genEpu{$name} ;
		$gene_epu_cord = $EpuData{$gene_epu} ;
	}
	else { $gene_epu = "No_overlap_gene" ; $gene_epu_cord = "No_overlap_gene" ;}

  my $match ;
	if($enh_epu eq $gene_epu){
		$match = "yes" ;
	}
  elsif($gene_epu eq "No_overlap_gene" && $enh_epu eq "No_overlap_enh"){
		$match = "No_overlap_both" ;
	}
	elsif($gene_epu eq "No_overlap_gene"){
		$match = "No_overlap_gene" ;
	}
  elsif($enh_epu eq "No_overlap_enh"){
		$match = "No_overlap_enh" ;
	}
	else{ $match = "no" ;}

	print out "$chr\t$start\t$end\t$id\t$enh_epu\t$enh_epu_cord\t$name\t$gene_epu\t$gCordinates\t$gene_epu_cord\t$match\n" ;

}
close in ; close out ;

