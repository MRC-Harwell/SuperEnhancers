#!/usr/bin/perl
#Calculate genomic coverage of chromatin states in 22 tissues
#Author: Sid

use Getopt::Long;
use Cwd ;

my ($path,$help,$out,$dir) = "" ;
GetOptions(
    'help'   => \$help,
    'o=s' => \$out ,
    'dir=s' => \$dir ,
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

my @tissues = (
	"BAT",
	"BmarrowDm",
	"Bmarrow",
	"Cerebellum",
	"CH12",
	"Cortex",
	"Esb4",
	"Es-E14",
	"Heart",
	"Kidney",
	"Limb",
	"Liver",
	"Lung",
	"MEF",
	"MEL",
	"OlfactoryBulb",
	"Placenta",
	"SmallIntestine",
	"Spleen",
	"Testis",
	"Thymus",
	"Wbrain"
) ;


my $dir = "/NGS/users/Sid/ChromHmm/LearnModel/State_6/Dense_modified" ;

foreach my $tis(@tissues){
  open(in, "<$dir/$tis\_6_dense.mod.bed") or die "Cannot open file $dir/$tis\_6_dense.mod.bed" ;
  my ($st1,$st2,$st3,$st4,$st5,$st6) = 0 ;
  foreach my $line(<in>){
    if ($line =~ /^track/) {next ;}
    my($chr,$start,$end,$state,$rest) = split("\t",$line,5) ;
    my $diff = $end - $start ;
    if($state eq 1) { $st1 = $diff + $st1 ;}
    elsif($state eq 2) {$st2 = $diff + $st2 ;}
    elsif($state eq 3) {$st3 = $diff + $st3 ;}
    elsif($state eq 4) {$st4 = $diff + $st4 ;}
    elsif($state eq 5) {$st5 = $diff + $st5 ;}
    elsif($state eq 6) {$st6 = $diff + $st6 ;}
  }
  close in ;

  my $total_bp = 2654911517 ;
  my $st1_per = ($st1/$total_bp)*100 ;
  my $st2_per = ($st2/$total_bp)*100 ;
  my $st3_per = ($st3/$total_bp)*100 ;
  my $st4_per = ($st4/$total_bp)*100 ;
  my $st5_per = ($st5/$total_bp)*100 ;
  my $st6_per = ($st6/$total_bp)*100 ;
 
  print "$tis\t$st1_per\t$st2_per\t$st3_per\t$st4_per\t$st5_per\t$st6_per\n" ;

}
