#!/usr/bin/perl
# modifying gene targets association from GREAT into a gene list, also removing replicates
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
                    -i -> input assoc file from GREAT
                    OPTIONAL
                    -p -> OUTPUT path where the files need to be saved [default = cwd ]
                    -help -> prints this help message
               \n/);
}

if($help) { &useage ;}
#if() { print "\n MISSING ARGUMENTS : Give all the required options\n" ; &useage ;}
#unless( ) { print "\n \n" ; &useage ;}
if( !$path || $path =~ /\./ ) { $cwd = getcwd ; $path = "$cwd" ; }


# Extracting all, no filtering on distance
my %unique ;
open(in,"<$in") or die "Cannot open $in\n" ;
foreach(<in>) {
    $_ =~ s/\r|\n|\s+$//g ;
    if($_ =~ /^#/) { next ; }
    my($id,$genes) = split("\t",$_) ;
    my @targets = split(",",$genes) ;
    foreach my $g(@targets) {
        $g =~ s/\s+\(.*\d+\)//g ;
        $g =~ s/^\s+// ;
        #print "$g\n" ;
        if( !exists $unique{$g}) {
            print "$g\n" ;
            $unique{$g} = $id ;
        }

    }

}
close in ;


####### For filtering with a distance threshold (+-50kb)
my %unique ;
open(in,"<$in") or die "Cannot open $in\n" ;
foreach(<in>) {
    $_ =~ s/\r|\n|\s+$//g ;
    if($_ =~ /^#/) { next ; }
    my($id,$genes) = split("\t",$_) ;
    my @targets = split(",",$genes) ;
    foreach my $g(@targets) {
        $g =~ s/^\s+// ;
        my($name,$dist) = split(/\s+/,$g) ;
        $dist =~ s/\(|\)|\+|-//g ;
        if ($dist <= 50000) {
            if( !exists $unique{$name}) {
                print "$name\n" ;
                $unique{$name} = $dist ;
            }
        }
    }
}
close in ;


# ####### For Plotting the DISTRIBUTION

open(in,"<$in") or die "Cannot open $in\n" ;
my $cluster = $in ;
$cluster =~ s/.enhancer-gene\.txt//g ;
open(out,">Distribution/$cluster.txt") ;
foreach(<in>) {
    $_ =~ s/\r|\n|\s+$//g ;
    if($_ =~ /^#/) { next ; }
    my($id,$genes) = split("\t",$_) ;
    my @targets = split(",",$genes) ;
    foreach my $g(@targets) {
        $g =~ s/^\s+// ;
        my($name,$dist) = split(/\s+/,$g) ;
        $dist =~ s/\(|\)//g ;
        $dist = $dist/1000 ;
        print out "$name\t$dist\t$cluster\tEnhancer\n" ;
    }
}
close in ; close out ;


