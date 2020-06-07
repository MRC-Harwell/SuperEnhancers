#!/usr/bin/perl
# to dissect chromatin segmentation bed files into 200 bps windows , replacing selected annotation with "1" and others with "0"
# doing this to use for clustering promoter and enhancer regions in different tissues
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
               
                    OPTIONAL
                    -p -> OUTPUT path where the files need to be saved [default = cwd ]
                    -help -> prints this help message
               \n/);
}

if($help) { &useage ;}
#if() { print "\n MISSING ARGUMENTS : Give all the required options\n" ; &useage ;}
#unless( ) { print "\n \n" ; &useage ;}
if( !$path || $path =~ /\./ ) { $cwd = getcwd ; $path = "$cwd" ; }


my $annotation_1 = "E2" ;    #annotation to keep
my $dir = "/NGS/users/Sid/ChromHmm/LearnModel/State_6" ;
my $max_states = 6 ;
my @files = (
	"BAT_$max_states\_segments.bed",
	"Bmarrow_$max_states\_segments.bed",
	"BmarrowDm_$max_states\_segments.bed",
	"Cerebellum_$max_states\_segments.bed",
	"CH12_$max_states\_segments.bed",
	"Cortex_$max_states\_segments.bed",
	"Esb4_$max_states\_segments.bed",
	"Es-E14_$max_states\_segments.bed",
	"Heart_$max_states\_segments.bed",
	"Kidney_$max_states\_segments.bed",
	"Limb_$max_states\_segments.bed",
	"Liver_$max_states\_segments.bed",
	"Lung_$max_states\_segments.bed",
	"MEF_$max_states\_segments.bed",
	"MEL_$max_states\_segments.bed",
	"OlfactoryBulb_$max_states\_segments.bed",
	"Placenta_$max_states\_segments.bed",
	"SmallIntestine_$max_states\_segments.bed",
	"Spleen_$max_states\_segments.bed",
	"Testis_$max_states\_segments.bed",
	"Thymus_$max_states\_segments.bed",
	"Wbrain_$max_states\_segments.bed"
) ;


system "mkdir -m a=rwx $dir/Clustering_matrix" ;

foreach my $seg(@files){
	print "$seg\n" ; 
	my $tissue = $seg ;
	$tissue =~ s/_6_segments.bed// ;
	open (in,"<$dir/$seg")or die "Cannot open file $dir/$seg \n";
	open(out,">$dir/Clustering_matrix/$tissue.$annotation_1.txt") ;
	print out "Coordinates\t$tissue\n" ;
	foreach(<in>) {
		my $value ;
		$_ =~ s/\n|\r|\s+$//g ;
		my($chr,$start,$end,$anno) = split("\t",$_) ;

		#if($anno eq $annotation_1 || $anno eq $annotation_2) { $value = 1 ; }
		if($anno eq $annotation_1) { $value = 1 ; }
		else { $value = 0 ; } 

		my $diff = $end - $start ;
		if ($diff != 200) {
			my $i ;
			for($i = $start ; $i < $end ; $i = $i + 200) {
				my $e = $i + 200 ;
				print out "$chr-$i-$e\t$value\n" ;
			}
		}
		else {print out "$chr-$start-$end\t$value\n" ;}
	}
	close in; close out ;
}





		
    
