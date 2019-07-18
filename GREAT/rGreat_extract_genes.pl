#!/usr/bin/perl
# annotate ROSE results back to the constituent enhacer elements
#Author: Sid

use Getopt::Long;
use Cwd ;
use List::Util qw( min max sum);

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
                    -i -> input results file from rGREAT
                    OPTIONAL
                    -p -> OUTPUT path where the files need to be saved [default = cwd ]
                    -help -> prints this help message
               \n/);
}

open(in,"<$in") or die "Cannot open $in\n" ;
my %data ; my %genes ; my %geneDist ;
foreach(<in>) {
  $_ =~ s/\r|\n|\s+$//g ;
  if($_ =~ /^chr\t/) { next ; }
  my($chr,$start,$end,$gene,$dist) = split("\t",$_) ;
  if(exists $data{"$chr-$start-$end"}){
    my $previous = $data{"$chr-$start-$end"} ;
    my $new = $dist ;
    $previous =~ s/\+|-//g ;
    $new =~ s/\+|-//g ;

    if($new < $previous){

      $data{"$chr-$start-$end"} = $dist ;
      $genes{"$chr-$start-$end"} = $gene ;
      $geneDist{"$chr-$start-$end"} = "$gene\t$dist" ;
    }
  }
  else{
    $data{"$chr-$start-$end"} = $dist ;
    $genes{"$chr-$start-$end"} = $gene ;
    $geneDist{"$chr-$start-$end"} = "$gene\t$dist" ;
  }
}
close in ;

my $out = $in ;
$out =~ s/\.rGreat.txt//g ;

my %unique ;
open(o1,">$out.genes.txt") ;
my @names = values %genes ;
foreach(@names){
  if(!exists $unique{$_}){
    $unique{$_} = undef ;
    print o1 "$_\n" ;
  }
}
close o1 ;

open(o2,">$out.dist.txt") ;
my @d = values %geneDist ;
foreach(@d){
  print o2 "$_\n" ;
}
close o2 ;


open(o3,">$out.anno.txt") ;
my @a = keys %geneDist ;
foreach(@a){
  my $key = $_ ;
  my $gene = $geneDist{$key} ;
  my($chr,$start,$end) = split("-",$key) ;
  print o3 "$chr\t$start\t$end\t$gene\n" ;
}
close o3 ;
