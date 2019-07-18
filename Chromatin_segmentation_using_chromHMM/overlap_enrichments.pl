#!/usr/bin/perl
#Calculate average coverage of chromatin states in 22 tissues
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
                    -dir -> input file
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


system "mkdir -m a=rwx $dir/OverlapEnrichment" ;
system "cat $dir/*overlap.txt | grep -v state | grep -v Base > $dir/OverlapEnrichment/all_overlap.txt" ;
system "Rscript /NGS/users/Sid/SCRIPTS/Chromatin_segmentation/pivot_table.R all_overlap.txt $dir/OverlapEnrichment/" ;


open(in, "<$dir/OverlapEnrichment/all_overlap_pivot.txt") or die "Cannot open file $dir/OverlapEnrichment/all_overlap_pivot.txt" ;
open(out,">$dir/OverlapEnrichment/all_overlap_table.txt") ;

foreach(<in>) {
    $_ =~ s/\n|\r|\"|\s+$//g ;
    if($_ =~ /^X1/) { next ;}
    else {
        my($tmp,$type,$state,$value) = split("\t",$_) ;
        my $orig_state = $state ;
        $state =~ s/X/E/ ;
        $state =~ s/\..+// ;
        print out "$type\t$state\t$value\n" ;
    }
}
close in ; close out ;

my @categories = ("Genome","CpGIsland","RefSeqExon","RefSeqGene","RefSeqTES","RefSeqTSS","RefSeqTSS2kb") ;
foreach(@categories) {
    system "grep -w $_ $dir/OverlapEnrichment/all_overlap_table.txt > $dir/OverlapEnrichment/$_.table.txt" ;
    system "Rscript /NGS/users/Sid/SCRIPTS/Chromatin_segmentation/overlap_enrichment_boxPlot.R $_.table.txt $dir/OverlapEnrichment" ;
}










     
