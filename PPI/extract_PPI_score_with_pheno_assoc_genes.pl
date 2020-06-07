#!/usr/bin/perl
# extracting protein interaction score with phenotype associated genes from STRING database data (locally downloaded)
# Author: Sid

use Getopt::Long;
use Cwd ;
use List::Util qw( min max sum);
use threads;
use List::MoreUtils qw(part) ;
use Data::Dumper ;

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
                    -i -> input gene list file
			-o -> output file name

                    OPTIONAL
                    -p -> OUTPUT path where the files need to be saved [default = cwd ]
                    -help -> prints this help message
               \n/);
}

if($help) { &useage ;}
#if() { print "\n MISSING ARGUMENTS : Give all the required options\n" ; &useage ;}
#unless( ) { print "\n \n" ; &useage ;}
if( !$path || $path =~ /\./ ) { $cwd = getcwd ; $path = "$cwd" ; }

my $mp_dir = "Mgi_14June_17/Dissection" ;
my @mps = ("MP:0005388") ;

# gene to protein ID STRING mapping file
my $mapping = "/NGS/users/Sid/StringDB_interactions/10090.protein.aliases.v10.5.txt";

# STRING interactions data
my $interactions = "/NGS/users/Sid/StringDB_interactions/10090.protein.links.detailed.v10.5.txt";

my %mapping ;
open(mapp,$mapping) or die "Cannot open $mapping" ;
foreach my $l(<mapp>){
  $l =~ s/\n|\r|\s+$//g ;
  if($l =~ /^##/){ next ;}
  my($id,$gene,$rest) = split("\t",$l,3) ;
  $mapping{$gene} = $id ;
}
close mapp ;

my %mp_genes ;
foreach my $m(@mps){
  open(mgi,"<$mp_dir/$m.genes.txt") or die "Cannot open $mp_dir/$m.genes.txt" ;
  my $total = 0 ; my $count = 0 ;
  foreach my $g(<mgi>){
    $total++ ;
    $g =~ s/\n|\r|\s+$//g ;
    if(exists $mapping{$g}){
      my $protein_id = $mapping{$g} ;
      $mp_genes{$protein_id} = undef ;
    }
    else { $count++ ;}
  }
  close mgi ;
  my $failed ;
  if($count == 0){
    $failed = 0 ;
  }
  else { $failed = ($count/$total)*100 ;}
  print "\n$failed % of genes names in $m failed to map with STRING protein IDs\n" ;
}


my @input ;
open(in,$in) or die "Cannot open $in" ;
foreach my $gene(<in>){
  $gene =~ s/\n|\r|\s+$//g ;
  push(@input,"$gene") ;
}
close in ;


# dividing my input data into small parts to run on multiple threads
my $partitions = 5 ;
my $i = 0 ;
my @data_part = part {$partitions * $i++ /@input} @input ; # this outputs an array of array references


my $j = 0 ; my @threads ;
foreach(@data_part){
	my @part = @{$_} ;
	$j++ ;
	$j = threads->create(\&analysis,\@part);
	push(@threads,$j) ;
}

open(out,">$out") ;
foreach(@threads){
   my @res = @{$_->join()};
   foreach my $r(@res) {
	   print out "$r\n" ;
   }
}
close out ;


sub analysis {
  my $id = threads->tid();
  my @array = @{$_[0]} ;
  my @RESULTS ;
  foreach my $gene(@array){
     my $final; my @score_val ;
    if(exists $mapping{$gene}){
      my $protein_id = $mapping{$gene} ;
      open(data,$interactions) or die "Cannot open $interactions" ;
      foreach my $line(<data>){
        if($line =~ /^protein/) { next ;}
        my ($pr1,$pr2,$neigh,$fusion,$coo,$coex,$exp,$db,$text,$score) = split(" ",$line) ;
        if ($score >= 900){
          if($pr1 eq $protein_id){
            if(exists $mp_genes{$pr2}){
              push(@score_val,$score) ;
            }
          }
          elsif($pr2 eq $protein_id){
            if(exists $mp_genes{$pr1}){
              push(@score_val,$score) ;
            }
          }
        }
      }
      close data ;
  		if (scalar @score_val == 0) {
  			$final = 0 ;
  		}
      	else {$final = sum (@score_val) ;}
    }
    else{
      $final = 0 ;
    }
    push(@RESULTS,"$gene\t$final") ;
	#print "$gene\t$final\n" ;
  }
  return \@RESULTS ;
}
