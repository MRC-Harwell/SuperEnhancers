#!/usr/bin/perl
# Matrix to BED file
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
                    -i -> input matrix file
                    OPTIONAL
                    -p -> OUTPUT path where the files need to be saved [default = cwd ]
                    -help -> prints this help message
               \n/);
}

if($help) { &useage ;}
#if() { print "\n MISSING ARGUMENTS : Give all the required options\n" ; &useage ;}
#unless( ) { print "\n \n" ; &useage ;}
if( !$path || $path =~ /\./ ) { $cwd = getcwd ; $path = "$cwd" ; }


#open(out,">$path/cluster_$i.bed");
open (in, $in)or die "Cannot open file $in \n";	
my $count = 0 ;
foreach(<in>){
	$_ =~ s/\r|\n|\s+$|\"//g ;
	unless($_ =~ /^Coordinates/) {
		$count++ ;
		my @array = split("\t",$_) ;
		my $cordinate = $array[0] ;			
		my($chr,$start,$end) = split("-",$cordinate) ;				
		print "$chr\t$start\t$end\tenhancer_$count\n" ;							
	}
}
close in ;





















