#!/usr/bin/perl
# calculating the tau (core calculation)
#Author: Sid

use Getopt::Long;
use Cwd ;
use List::Util qw( min max );

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
                    -i -> input matrix file from part 1 script
                    OPTIONAL
                    Change the ensemble file if you want the annotation as well
                    -p -> OUTPUT path where the files need to be saved [default = cwd ]
                    -help -> prints this help message
               \n/);
}

if($help) { &useage ;}
#if() { print "\n MISSING ARGUMENTS : Give all the required options\n" ; &useage ;}
#unless( ) { print "\n \n" ; &useage ;}
if( !$path || $path =~ /\./ ) { $cwd = getcwd ; $path = "$cwd" ; }

########## If the first column of the matrix needs to be annotated ###############
# ensembl file
my $anno = "/NGS/users/Sid/ChromHmm/Posterior/State_6/POSTERIOR/Tau_expression/Ens+refGene_mm9/annotation.txt" ;
open(anno,$anno) or die "cannot open $anno" ;
my %annotation ;
foreach(<anno>) {
  $_ =~ s/\n|\r|\s+$//g ;
  if($_ =~ /^Ensembl/) { next ;}
  #my($id,$name,$type,$status) = split("\t",$_) ;
  my($id,$name,$length) = split("\t",$_) ;
  $annotation{$id} = $name ;
}
close anno ;


open(in,$in) or die "cannot open $in" ;
my @header;
foreach(<in>) {
    $_ =~ s/\n|\r|\s+$//g ;
    if($_ =~ /^\t/) {
        my @tissues = split("\t",$_) ;
        splice(@tissues,0,1) ;
        foreach my $t(@tissues) {
            push(@header,$t) ;
        }
        next ;
    }
    else {
        my ($cord,$rest) = split("\t",$_,2) ;
        my @elements = split("\t",$rest) ;
        my $max = max @elements ;
        my $max_1 = $max - 1 ;

        # get the indices of max values
        my @indices ; my $count = 0 ; my $sum ;
        foreach my $el(@elements) {
            $count++ ;
            #if ($el == $max || $el == $max_1) { # max and max-1
            if ($el == $max) {                   # max
                if($el != 0) {
                    push(@indices,$count - 1) ;
                }
            }

            # tau calculation
            my $norm = 1 - ($el/$max) ;
            $sum = $norm + $sum ;
        }

        my $total_tissues = scalar @elements ;
        my $tau = $sum/($total_tissues-1) ;

        #getting tissues with max bin
        my @m_tis ;
        foreach my $i(@indices) {
            my $m = $header[$i] ;
            if($m eq "Es.E14") { $m = "Es-E14" ;}
            push(@m_tis,$m) ;
        }
        my $m_tis_names = join(",",@m_tis) ;
        #print "$cord\t$tau\t$m_tis_names\n" ;
        ## if annotate the cord
        if(exists $annotation{$cord}){
          my $value = $annotation{$cord} ;
          print "$cord\t$value\t$tau\t$m_tis_names\n" ;
        }
        else { print "$cord\tNOT_FOUND\n" ;}
    }
}
close in ;
