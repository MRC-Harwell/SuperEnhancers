#!/usr/bin/perl
# making of the matrix of super enhancers and typical enh for cell type specificity (distinct enhancer types)
# this script was also used to make enhacer matrix for Random forest
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

                    OPTIONAL
                    -p -> OUTPUT path where the files need to be saved [default = cwd ]
                    -help -> prints this help message
               \n/);
}

if($help) { &useage ;}
#if() { print "\n MISSING ARGUMENTS : Give all the required options\n" ; &useage ;}
#unless( ) { print "\n \n" ; &useage ;}
if( !$path || $path =~ /\./ ) { $cwd = getcwd ; $path = "$cwd" ; }

# Gene list - all super enhancers genes/ typical enh genes
my $gene_file = "/NGS/users/Sid/ENHANCERS/SuperEnh_27ac/Great_targets/Strict_targets/CellTypeSpecific_analysis/all.te.uniq.genes.txt" ;

#my $gene_file = " /NGS/users/Sid/ChromHmm/Posterior/State_6/POSTERIOR/Tau_enhancers/TauResults/Highly_specific/SuperEnh/AllGenes_super/new_pc_genes_list.txt" ;

# super enhancers /typical enhancers
my $dir= "/NGS/users/Sid/ENHANCERS/SuperEnh_27ac/Great_targets/Strict_targets" ;
#my $dir = "/NGS/users/Sid/ENHANCERS/SuperEnh_27ac/Great_targets/Strict_targets/Ppi_scores" ;

my @tissues = (
	"BAT",
	"BmarrowDm",
	"Bmarrow",
	"Cerebellum",
	"CH12",
	"Cortex",
	"Esb4",
	"Es-E14",
	"Heart",
	"Kidney",
	"Limb",
	"Liver",
	"Lung",
	"MEF",
	"MEL",
	"OlfactoryBulb",
	"Placenta",
	"SmallIntestine",
	"Spleen",
	"Testis",
	"Thymus",
	"Wbrain"
) ;


open(in,"<$gene_file") or die "Cannot open $gene_file" ;

foreach my $gene(<in>){
	$gene =~ s/\n|\r|\s+$//g ;
	print "$gene\t" ;
	my $count = 0 ;
	foreach my $tis(@tissues) {
		#print "$tis ...\n" ;
		$count++ ;
		my $result ;
		#my $lines = `grep -w $gene $dir/$tis.ppiScore.txt | wc -l` ; chomp($lines) ;
    my $lines = `grep -w $gene $dir/$tis.te.uniq.genes.txt | wc -l` ; chomp($lines) ;
    if($lines == 0) {
			$result = 0 ;
		}
		else {
			$result = 1 ;
      # my $grep = `grep -w $gene $dir/$tis.ppiScore.txt` ;
      # my @grep = split("\n",$grep) ;
	    # foreach my $l(@grep){
      #   $l =~ s/\n|\r|\s+$//g ;
      #   my($g,$sc) = split("\t",$l) ;
      #   $result = $sc ;
      #}
		}
		if($count == 22){
			print "$result" ;
		}
		else{ print "$result\t" ; }
	}
	print "\n" ;
}
close in ;
