#!/usr/bin/perl

	local($file,$motif,$file2)=@ARGV;
	
       	if($file eq "" || $motif eq "" || $file2 eq ""){
                  die "usage: adddatablock.pl <file.sdf> <datablock name to add> <file with one data per line>\n";
        };
	
	@data="";
	$i=1;
	open(IN,"<$file2");
	while (<IN>){
		if($_!~/^#/){
			@get=split(' ',$_);
			$data[$i]=$get[1];
			$i++;
		};
	};
	close(IN);
	$i=$i-1;
	print "$i data in $file2\n";

	open(OUT1,">adddatablock.sdf");
	$flagnew=1;
	@ligne='';
	$i=0;
	$j=0;
	$name="";
	open(IN,"<$file");
	while (<IN>){
		$conv2=$_;
		@get=split(' ',$_);
		if($flagnew){
			$flagnew=0;
			$flagname=0;
			@ligne='';
			$i=0;
			$j++;
			$name="";
		};

		$ligne[$i]=$_;
		$i++;
		if($conv2 =~/(\$\$\$\$)/){
			$flagnew=1;
			foreach $k (0..@ligne-2){
				print OUT1 $ligne[$k];
			};
			if($data[$j] ne ""){
				print OUT1 "> <$motif>\n";
				print OUT1 "$data[$j]\n";
				print OUT1 "\n";
			};
			print OUT1 "$conv2";
		};
	};
	close(IN);
  	close(OUT1);
	print "see adddatablock.sdf\n";

