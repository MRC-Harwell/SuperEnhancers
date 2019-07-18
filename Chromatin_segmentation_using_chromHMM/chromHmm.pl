#!/usr/bin/perl
# runs chromatin segmentation algorithm
#Author: Sid

use Getopt::Long;
use Cwd ;

my ($path,$help,$out,$fileTable,$num) = "" ;
GetOptions(
    'help'   => \$help,
    'o=s' => \$out ,
    'fileTable=s' => \$fileTable ,
    'p=s' => \$path ,
    'n=i' => \$num ,

) or die "\n**********  Incorrect usage!  ***********\nrun with -help option to see the useage\n ";

sub useage { die(qq/
	USAGE : perl <script> <arguments>
	ARGUMENTS :
                    REQUIRED
                    -fileTable -> input text file with file names [ cell\\tmark\\tbedFile\\tcontrol_bedFile(optional)\\n]
                    -o -> name of the output files [ example : give cell\/tissue name ]
                    -n <int> -> number of states
                    OPTIONAL
                    -p -> OUTPUT path where the files need to be saved [default = cwd ]
                    -help -> prints this help message
               \n/);
}

if($help) { &useage ;}
#if( !$fileTable || !$num || !$out) { print "\n MISSING ARGUMENTS : Give all the required options\n" ; &useage ;}
#unless( ) { print "\n \n" ; &useage ;}
if( !$path || $path =~ /\./ ) { $cwd = getcwd ; $path = "$cwd" ; }


$inputDir = "/NGS/working_projects/ENCODE/LicrHistoneMouse/BroadPeak/" ;
system "mkdir -m a=rwx $path/$out.binaryFiles " ;
print "converting input bed files to binarized output .....\n" ;

system "java -mx16g -jar /NGS/Software/ChromHMM/ChromHMM.jar BinarizeBed -c /NGS/working_projects/ENCODE/LicrHistoneMouse/BAM/Pooled -o outputcontrolSignal -t outputsignal /NGS/Software/ChromHMM/CHROMSIZES/mm9.txt $inputDir $fileTable $path/$out.binaryFiles" ;


################ CHROMHMM LEARNMODEL REQUIRES X11 FORWARDING ####################
for(my $i=2;$i<=8;$i++){
  print "STATE_$i.........\n" ;
  system "mkdir -m a=rwx /$path/LearnModel/State_$i" ;
  system "java -mx16g -jar /NGS/Software/ChromHMM/ChromHMM.jar LearnModel -color 153,153,0 -l /NGS/Software/ChromHMM/CHROMSIZES/mm9.txt -p 10 -printposterior $path/$out.binaryFiles $path/LearnModel/State_$i $i mm9" ;
}
