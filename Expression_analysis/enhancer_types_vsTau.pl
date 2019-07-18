#!/usr/bin/perl
# EXPRESSION ANALYSIS TAU
#Author: Sid

use Getopt::Long;
use Cwd ;
use List::Util qw( min max sum);
#use Statistics::R;

my ($path,$help,$out,$reg,$exp) = "" ;
GetOptions(
    'help'   => \$help,
    'o=s' => \$out ,
    'exp=s' => \$exp ,
    'reg=s' => \$reg ,
    'p=s' => \$path ,

) or die "\n**********  Incorrect usage!  ***********\nrun with -help option to see the useage\n ";

sub useage { die(qq/
	USAGE : perl <script> <arguments>
	ARGUMENTS :
                    REQUIRED

                    OPTIONAL
                    -p -> OUTPUT path where the files need to be saved [default = cwd ]
                    -help -> prints this help message
               \n/);
}

if($help) { &useage ;}
#if() { print "\n MISSING ARGUMENTS : Give all the required options\n" ; &useage ;}
#unless( ) { print "\n \n" ; &useage ;}
if( !$path || $path =~ /\./ ) { $cwd = getcwd ; $path = "$cwd" ; }


#TAU results file
my $tau_file = "/NGS/users/Sid/ENHANCERS/SuperEnh_27ac/Great_targets/Strict_targets/Tau_analysis/RESULTS.txt" ;

# all GREAT genes
my $allGenes = "/NGS/users/Sid/ENHANCERS/SuperEnh_27ac/Great_targets/Strict_targets/CellTypeSpecific_analysis/VsTau/all.genes.cellTypeSpec.txt" ;


open(out,">all.genes.tau.txt") ;
open(nfound,">not_found.txt") ;

open(in,$allGenes) or die "Cannot open $allGenes" ;
foreach(<in>){
	$_ =~ s/\n|\r|\s+$//g ;
	my($gene,$cellTypeNum,$group) = split("\t",$_,3) ;
  my $lines = `grep -w $gene $tau_file | grep $group | wc -l` ; chomp($lines) ;
  if($lines == 0) {
    print nfound "$_\n" ;
  }
  else{
    my $grep = `grep -w $gene $tau_file | grep $group` ;
    my @grep = split("\n",$grep) ;
    foreach my $l(@grep){
      my($g,$tauf,$name,$tis) = split("\t",$l) ;
      print out "$gene\t$cellTypeNum\t$group\t$tauf\n" ;
    }
  }
}
close in ; close out ; close nfound ;
