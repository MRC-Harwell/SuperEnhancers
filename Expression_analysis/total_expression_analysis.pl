#!/usr/bin/perl
# EXPRESSION ANALYSIS
#Author: Sid

use Getopt::Long;
use Cwd ;
use List::Util qw( min max sum);
#use Statistics::R;

my ($path,$help,$out,$reg,$exp) = "" ;
GetOptions(
    'help'   => \$help,
    'o=s' => \$out ,
    'exp=s' => \$exp ,
    'reg=s' => \$reg ,
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

my %tissueReplace = (
	"1" => "BAT",
	"2" => "BmarrowDm",
	"3" =>"Bmarrow",
	"4" =>"Cerebellum",
	"5" =>"CH12",
	"6" =>"Cortex",
	"7" =>"Esb4",
	"8" =>"Es-E14",
	"9" =>"Heart",
	"10" =>"Kidney",
	"11" =>"Limb",
	"12" =>"Liver",
	"13" =>"Lung",
	"14" =>"MEF",
	"15" =>"MEL",
	"16" =>"OlfactoryBulb",
	"17" =>"Placenta",
	"18" =>"SmallIntestine",
	"19" =>"Spleen",
	"20" =>"Testis",
	"21" =>"Thymus",
	"22" =>"Wbrain"
) ;

#expression matrix protein coding
my $exp_file = "/NGS/users/Sid/ENHANCERS/SuperEnh_27ac/AllGenes_exp/mat.txt" ;
# super enhancers and typical enhancers
my $dir_1 = "/NGS/users/Sid/ENHANCERS/SuperEnh_27ac/Great_targets/Strict_targets" ;
# weak enhancers
my $dir_2 = "/NGS/users/Sid/ChromHmm/Posterior/State_6/POSTERIOR/Tau_weak_enhancer/Results/Highly_specific/Great_targets/Strict_targets" ;
# all GREAT genes
my $allGenes = "/NGS/users/Sid/RNA_seq/gtfTsre/Great_genes_mm9.txt" ;

my %expression ;
open(ex,$exp_file) or die "Cannot open $exp_file" ;
foreach(<ex>){
	$_ =~ s/\n|\r|\s+$//g ;
	my($gene,$rest) = split("\t",$_,2) ;
	$expression{$gene} = $rest ;
}
close ex ;

open(out,">RESULTS.txt") ;
open(nfound,">not_found.txt") ;
foreach my $tis(@tissues){
	print "$tis\n";

	my %geneCheck ; # to make sure no genes repeat in each group
	# super enhancers
	open(se,"$dir_1/$tis.se.genes.txt") or die "Cannot open $dir_1/$tis.se.genes.txt" ;
	foreach(<se>){
		$_ =~ s/\n|\r|\s+$//g ;
		my($gene,$tmp) = split("\t",$_,2) ;
		$geneCheck{$gene} = "super" ;
		if(exists $expression{$gene}){
			my $value = $expression{$gene} ;
			my @elements = split("\t",$value) ;
			my $tisCount = 0 ;
			foreach my $el(@elements) {
				$tisCount++ ;
				my $tissue = $tissueReplace{$tisCount} ;
				if($tis eq $tissue){
					print out "Super\t$el\t$gene\t$tissue\n" ;
				}
			}
		}
		else { print nfound "$gene\t$tis\tSuper\n" ;}
	}
	close se ;

	# typical enhancers
	open(te,"$dir_1/$tis.te.uniq.genes.txt") or die "Cannot open $dir_1/$tis.te.uniq.genes.txt" ;
	foreach(<te>){
		$_ =~ s/\n|\r|\s+$//g ;
		my($gene,$tmp) = split("\t",$_,2) ;
		$geneCheck{$gene} = "typical" ;
		if(exists $expression{$gene}){
			my $value = $expression{$gene} ;
			my @elements = split("\t",$value) ;
			my $tisCount = 0 ;
			foreach my $el(@elements) {
				$tisCount++ ;
				my $tissue = $tissueReplace{$tisCount} ;
				if($tis eq $tissue){
					print out "Typical\t$el\t$gene\t$tissue\n" ;
				}
			}
		}
		else { print nfound "$gene\t$tis\tTypical\n" ;}
	}
	close te ;

	# weak enhancers
	open(we,"$dir_2/$tis.genes.txt") or die "Cannot open $dir_2/$tis.genes.txt" ;
	foreach(<we>){
		$_ =~ s/\n|\r|\s+$//g ;
		my($gene,$tmp) = split("\t",$_,2) ;
		if (exists $geneCheck{$gene}){ next ;}
		else {
			$geneCheck{$gene} = "weak" ;
			if(exists $expression{$gene}){
				my $value = $expression{$gene} ;
				my @elements = split("\t",$value) ;
				my $tisCount = 0 ;
				foreach my $el(@elements) {
					$tisCount++ ;
					my $tissue = $tissueReplace{$tisCount} ;
					if($tis eq $tissue){
						print out "Weak\t$el\t$gene\t$tissue\n" ;
					}
				}
			}
			else { print nfound "$gene\t$tis\tWeak\n" ;}
		}
	}
	close we ;

	# Absent
	open(all,$allGenes) or die "Cannot open $allGenes" ;
	foreach(<all>){
		$_ =~ s/\n|\r|\s+$//g ;
		my($id,$chr,$pos,$strand,$gene) = split("\t",$_) ;
		if (exists $geneCheck{$gene}){ next ;}
		else {
			$geneCheck{$gene} = "absent" ;
			if(exists $expression{$gene}){
				my $value = $expression{$gene} ;
				my @elements = split("\t",$value) ;
				my $tisCount = 0 ;
				foreach my $el(@elements) {
					$tisCount++ ;
					my $tissue = $tissueReplace{$tisCount} ;
					if($tis eq $tissue){
						print out "Absent\t$el\t$gene\t$tissue\n" ;
					}
				}
			}
			else { print nfound "$gene\t$tis\tAbsent\n" ;}
		}
	}
	close all ;

}
close out ; close nfound ;
