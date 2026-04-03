#!/usr/bin/perl

$f='';
$option='';
$fileout='';

my($f,$fileout,$option,$k)=@ARGV;

if($f eq '' || $option eq '' || $fileout eq '' || $k ne ''){
	  print "usage: stripsdf <file.sdf> <Outputfile.sdf> <Remained DATA BLOCK NAMES in quotes >\n";
	  exit;
};
$option=" ".$option." ";
$ok=1;
$nb=0;
open(OUT,">$fileout");
open(IN,"<$f");
while(<IN>){
	if($_=~/^>/){
		$tmp=$_;
		$tmp=~s/<//g;
		$tmp=~s/>//g;
		@get=split(' ',$tmp);
		#print "$get[0]\n";
		if($option=~/ $get[0] /){
			$ok=1;
		}
		else{
			$ok=0;
		};
	};
	$ok=1 if($_=~/^\$\$\$\$/);
 	print OUT $_ if($ok);
 	$nb++ if($_=~/^\$\$\$\$/);
};
close(IN);
close(OUT);
#print "$nb molecules\n";

