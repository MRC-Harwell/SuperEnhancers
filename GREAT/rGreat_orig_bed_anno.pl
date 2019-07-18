#!/usr/bin/perl
# annotate originl BED with rGREAT results
#Author: Sid

use Getopt::Long;
use Cwd ;

my ($path,$help,$out,$in,$bed) = "" ;
GetOptions(
    'help'   => \$help,
    'o=s' => \$out ,
    'i=s' => \$in ,
	'bed=s' => \$bed ,
    'p=s' => \$path ,

) or die "\n**********  Incorrect usage!  ***********\nrun with -help option to see the useage\n ";

sub useage { die(qq/
	USAGE : perl <script> <arguments>
	ARGUMENTS :
                    REQUIRED
                    -i -> input .anno file from rGREAT script
					-bed -> original BED
                    OPTIONAL
                    -p -> OUTPUT path where the files need to be saved [default = cwd ]
                    -help -> prints this help message
               \n/);
}


open(in,"<$in") or die "Cannot open $in\n" ;
my %data ;
foreach(<in>) {
  $_ =~ s/\r|\n|\s+$//g ;
  if($_ =~ /^chr\t/) { next ; }
  my($chr,$start,$end,$gene,$dist) = split("\t",$_) ;
  $data{"$chr-$start-$end"} = "$gene\t$dist" ;
}
close in ;

open(bed,"<$bed") or die "Cannot open $bed\n" ;
foreach(<bed>){
    $_ =~ s/\r|\n|\s+$//g ;
    if($_ =~ /^chr\t/) { next ; }
    my($chr,$start,$end,$rest) = split("\t",$_,4) ;
	my $value ;
  	if(exists $data{"$chr-$start-$end"}){
    	$value = $data{"$chr-$start-$end"} ;
	}
	else{$value = "NA\tNA" ;}
    print "$_\t$value\n" ;
}
close bed ;







