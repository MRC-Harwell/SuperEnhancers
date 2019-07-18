#!/usr/bin/perl
# run ROSE algorithm
#Author: Sid

use Getopt::Long;
use Cwd ;

my ($path,$help,$out,$bed,$bam,$control) = "" ;
GetOptions(
    'help'   => \$help,
    'o=s' => \$out ,
    'i=s' => \$bed ,
    'r=s' => \$bam ,
    'c=s' => \$control ,
    'p=s' => \$path ,

) or die "\n**********  Incorrect usage!  ***********\nrun with -help option to see the useage\n ";
if( !$path || $path =~ /\./ ) { $cwd = getcwd ; $path = "$cwd" ; }


# Running ROSE
system "python /NGS/Software/rose/ROSE_main.py -i $bed -g MM9 -o $path -r $bam -c $control -t 2000" ;

#For making points plot of rose
#R --no-save /NGS/users/Sid/ENHANCERS/Rose/Bmarrow/Pol2/ Bmarrow_12KB_STITCHED_TSS_DISTAL_ENHANCER_REGION_MAP.txt Bmarrow NONE < /NGS/Software/rose/ROSE_callSuper.R


# For Metagene profile matrix ####
system "python /NGS/Software/rose/ROSE_bamToGFF.py -i $bed -b $bam -o $path/$out -m 300 -r" ;
