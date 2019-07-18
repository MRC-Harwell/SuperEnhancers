#!/usr/bin/perl
# tau results post processing
#Author: Sid 

use Getopt::Long;
use Cwd ;

my ($path,$help,$out,$in,$matrix) = "" ;
GetOptions(
    'help'   => \$help,
    'o=s' => \$out ,
    'i=s' => \$in ,
    'm=s' => \$matrix ,
    'p=s' => \$path ,
    
) or die "\n**********  Incorrect usage!  ***********\nrun with -help option to see the useage\n "; 

sub useage { die(qq/
	USAGE : perl <script> <arguments>
	ARGUMENTS : 
                    REQUIRED
                    -i -> input tau results file
                    -m -> imput matrix file
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
    "Es.E14",
    "Heart",
    "Kidney",
    "Limb",
    "Liver",
    "Lung",
    "MEF",
    "MEL" ,
    "OlfactoryBulb",
    "Placenta",
    "SmallIntestine",
    "Spleen" ,
    "Testis" ,
    "Thymus" ,
    "Wbrain",
) ;

my %matrix ;
open(mat,$matrix) or die "cannot open $matrix" ;
foreach my $line(<mat>) {
    $line =~ s/\n|\r|\s+$//g ;
    if($line =~ /^Coordinates/) {
        next ;
    }
    my($cord,$values) = split("\t",$line,2) ;
    $matrix{$cord} = $values ;
}
close mat ;


# highly tissue-specific
foreach my $tis(@tissues) {
    print "extracting for $tis\n" ;
    open(in,$in) or die "cannot open $in" ;
    open(out,">$tis.matrix.txt") ;
    foreach my $line(<in>) {
        $line =~ s/\n|\r|\s+$//g ;
        my($cord,$tau,$tissue) = split("\t",$line) ;
        my @array = split(",",$tissue) ;
        if($tau > 0.85 && $tissue eq "$tis") {
        	if( exists $matrix{$cord}) {
            	print out "$cord\t$matrix{$cord}\n";
            }
            else { print "Error: Could not find $cord in the matrix file $matrix\n";}
        }
    }
    close in ; close out ;
    system "perl /NGS/users/Sid/SCRIPTS/Chromatin_segmentation/matrixToBed.pl -i $tis.matrix.txt > $tis.bed" ;
}


# intermediate
open(in,$in) or die "cannot open $in" ;
open(out,">intermediate.matrix.txt") ;
foreach my $line(<in>) {
    $line =~ s/\n|\r|\s+$//g ;
    my($cord,$tau,$tissue) = split("\t",$line) ;
    my @array = split(",",$tissue) ;
    if($tau < 0.85 && $tau > 0.15) {
    	if( exists $matrix{$cord}) {
        	print out "$cord\t$matrix{$cord}\n";
        }
        else { print "Error: Could not find $cord in the matrix file $matrix\n";}
    }
}
close in ; close out ;


# low
open(in,$in) or die "cannot open $in" ;
open(out,">low.matrix.txt") ;
foreach my $line(<in>) {
    $line =~ s/\n|\r|\s+$//g ;
    my($cord,$tau,$tissue) = split("\t",$line) ;
    my @array = split(",",$tissue) ;
    if($tau <= 0.15) {
    	if( exists $matrix{$cord}) {
        	print out "$cord\t$matrix{$cord}\n";
        }
        else { print "Error: Could not find $cord in the matrix file $matrix\n";}
    }
}
close in ; close out ;



foreach my $tis(@tissues) {
    print "extracting for $tis\n" ;
    open(in,$in) or die "cannot open $in" ;
    open(out,">$tis.matrix.all.txt") ;
    foreach my $line(<in>) {
        $line =~ s/\n|\r|\s+$//g ;
        my($cord,$tau,$tissue) = split("\t",$line) ;
        my @array = split(",",$tissue) ;
        if($tau > 0.85) {
            foreach my $t(@array) {
                if($t eq "$tis") {
                    if( exists $matrix{$cord}) {
                        print out "$cord\t$matrix{$cord}\n";
                    }
                    else { print "Error: Could not find $cord in the matrix file $matrix\n";}
                }
            }
        }
    }
    close in ; close out ;
    #system "perl /NGS/users/Sid/SCRIPTS/Chromatin_segmentation/matrixToBed.pl -i $tis.matrix.txt > $tis.bed" ;
    last ;
}

