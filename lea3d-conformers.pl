#!/usr/bin/perl

	$file='';
	$option='';
	$nbconf=1;

	my($file,$option,$nbconf)=@ARGV;
	if($file eq '' || $option eq ''){
		print "usage: \n lea3d-conformers.pl <file.sdf> <rdkit or rdkitmin> <nb_max_conformers (default value = 1); rdkitmin with only 1 conformer >\n";
        	exit;
	};

	$leaexe=$0;
        #windows requires 2 steps
        $leaexe=~s/lea3d-conformers\.pl$//;
        $leaexe=~s/\/$//;
        #print "perl scripts in $leaexe\n";

	$nbconf=1 if($option eq "rdkitmin");
	$nomfile3d=$file;
	$nomfile3d=~s/.sdf//;
	$nomfile3d=$nomfile3d.'3D'.".sdf";

	if($option eq "rdkit" || $option eq "rdkitmin"){
		unlink "rdkit.sdf" if(-e "rdkit.sdf");
		open(OUT,">rdkit.sdf");
		close(OUT);
		$nb=&nbsdf($file);
		foreach $i (1..$nb){
			&searchsdfi($file,$i,"tmp.sdf");
			if(-e "tmp.sdf"){
				if($option eq "rdkit"){
					chop($rdkit=`python $leaexe/rdkit-confs.py tmp.sdf $nbconf tmp3D.sdf`);
				}
				else{
					chop($rdkit=`python $leaexe/rdkit-min-iter-long.py tmp.sdf tmp3D.sdf`);
				};
				#print "$rdkit\n";
				if(-e "tmp3D.sdf"){
					open(OUT,"<tmp3D.sdf");
					open(IN,">>rdkit.sdf");
					while(<OUT>){
						print IN $_;
					};
					close(OUT);
					close(IN);
				}
				else{
					print "$file no $i - rdkit failed\n";
				};
				unlink "tmp3D.sdf";
				unlink "tmp.sdf";
			}
			else{
				print "$file no $i not found ?\n";
			};
		};
		if(-e "rdkit.sdf" && !-z "rdkit.sdf"){
			rename "rdkit.sdf", "$nomfile3d";
			print "$nomfile3d Done\n";
		}
		else{
			print "$nomfile3d failed\n$rdkit\n";
		};
	};

###########################################################################################
#SUBROUTINES

sub nbsdf{
        local($fsdf)=@_;

        $nbligne=0;
        $nbd=0;
        open(IN,"<$fsdf");
        while(<IN>){
                $nbd++ if($_=~/\$\$\$\$/);
                $nbligne++;
        };
        close(IN);
        if($nbligne >=5 && $nbd==0){#case of .mol
                $nbd=1;
        };
        #print "$nbd\n";
        $nbd;
};

###################################################################################
###################################################################################

sub searchsdfi{
        local($fsdf,$motif,$fsdfo)=@_;

        $flag=1;
        $flagname=0;
        $flagnew=1;
        $fit=0;
        @ligne='';
        $i=0;
        open(IN,"<$fsdf");
        open(OUT,">$fsdfo");
        while (<IN>){
                $conv2=$_;
                if ($flagnew){
                        $fit++;
                        $flagnew=0;
                        $flagname=0;
                        @ligne='';
                        $i=0;
                };
                $ligne[$i]=$_;
                $i++;
                if ($fit == $motif){
                        $flagname=1;
                };
                if ($conv2 =~/(\$\$\$\$)/){
                        $flagnew=1;
                        if ($flagname){
                                foreach $i (0..@ligne-1){
                                        print OUT $ligne[$i];
                                };
                                $flagname=0;
                                last;
                        };
                };
        };
        close(IN);
        close(OUT);
};

###################################################################################
###################################################################################

