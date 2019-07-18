#!/usr/bin/perl
#Annotates the emmision states from ChromHmm [BED file ] 
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
	USAGE : perl <script> <arguments> <BED files>
	ARGUMENTS : 
                    REQUIRED
                    LOOK usage
                    
                    OPTIONAL
                    -p -> OUTPUT path where the files need to be saved
                    -help -> prints this help message
               \n/);
}

if($help) { &useage ;}
#if(!$in) { print "\n MISSING ARGUMENTS : Give all the required options\n" ; &useage ;}
#unless( ) { print "\n \n" ; &useage ;}
if( !$path || $path =~ /\./ ) { $cwd = getcwd ; $path = "$cwd" ; }


foreach my $bed(@ARGV) {
	print "\nprocessing file $bed\n" ;
	my $name = $bed ;
	$name =~ s/\.bed// ;

	open (in, $bed)or die "Cannot open file $bed \n";
	open(out, ">$path/$name.mod.bed");
	foreach(<in>){		
		my ($emision,$colour) ;
		if($_ =~ /^track.*/) { print out "$_" ; next ; }
		else {
			$_ =~ s/\n|\r|\s+$//g ;
			my($chr,$start,$end,$state,$tmp1,$tmp2,$tmp3,$tmp4,$rgb) = split("\t",$_) ;
			if ($state eq 1) {$emision = "promoter" ; $colour = "255,102,102";  }
			elsif ($state eq 2) {$emision = "active_promoter" ; $colour = "255,0,0" ; }
			elsif ($state eq 3) {$emision = "active_regulatory" ; $colour = "176,224,230" ; }
			elsif ($state eq 4) {$emision = "strong_enhancer" ; $colour = "255,165,0" ;}
			elsif ($state eq 5) {$emision = "heterochromatin" ; $colour = "255,255,255" ; }
			elsif ($state eq 6) {$emision = "weak_enhancer" ; $colour = "255,255,153" ; }			
			else { print "ERROR : $_\n" ; }
			print out "$chr\t$start\t$end\t$state\t$tmp1\t$tmp2\t$tmp3\t$tmp4\t$colour\n" ;			
		}
	} 

	close in ; close out ;
	print"output BED wriiten to $name.mod.bed\n" ;

	print "indexing $bed for Igv .....\n";
	system "qsub -cwd -j y -b yes -N $name -P NGS -o $path/GridLogs/$prefix -pe big 1 java -Xmx1200M -jar /NGS/Software/IGVTools/igvtools.jar index $path/$name.mod.bed" ;
	print"-----------------------------------------------------------\n" ;
}


          
