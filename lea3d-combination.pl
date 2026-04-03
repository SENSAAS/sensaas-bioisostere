#!/usr/bin/perl

   local($init,$type,$full,$isolated)= @ARGV;

   $leaexe=$0;
   #windows requires 2 steps
   $leaexe=~s/lea3d-combination\.pl$//;
   $leaexe=~s/\/$//;
   if($leaexe eq ""){
   	$leaexe=".";
   };
   #print "perl scripts in $leaexe\n";
   
   	#set 1 if want to keep original sdf datablock
	$ifdatablock=1;

   	#set the maximum combinations allowed (time cost)
   	$nbmaxcombin=1000;

	#windows conda search user_generic_name file until $nbfgtmax
        $nbfgtmax=20;

	if($init eq '' || $type eq ''){
		die "usage: <>.pl <file.sdf> <type: fgts or parts or user_generic_name (eg: part_)> 2 0 \n"; 
	};
	
	$nbinput=&nbsdf($init);
	if($nbinput > 1 && ($type eq "fgts" || $type eq "parts")){
		#die "$init must contain only one structure\n";
		#print "$init contains $nbinput structures. Combinations will be generated for each\n";
	}
	elsif($nbinput > 1){
		$nbinput=1; #$init file is not used in case of user_generic_name 
	};

	if($full==2){
		$isolated=0;
	};
	$filesdf="combination.sdf";

#Loop over $nbinput
foreach $moli (1..$nbinput){
	&searchsdfi($init,$moli,"init-tmp.sdf");
	$init2="init-tmp.sdf";
	
	foreach $rmi (0..$nbmaxcombin){
                $namerm="combin_".$rmi.".sdf";
                if(-e $namerm){
                        unlink $namerm;
                };
        };

### retrieve datablock info
	if($ifdatablock==1){
		$datablock="";
		$gosdf=0;
		open(SDF,"<$init2");
		while(<SDF>){
			$conv2=$_;
			$gosdf=0 if($conv2 =~/(\$\$\$\$)/);
			$gosdf=1 if ($conv2 =~/^>/);
			if($gosdf){
				$datablock=$datablock."$_";
			};
		};
		close(SDF);
	};

### Fragmente le SDF

	$list6="";
	#windows conda:
        $nblist6=0;
        foreach $ni6 (1..$nbfgtmax){
                $namelist6=$type.$ni6.".sdf";
                if(-e "$namelist6" && !-z "$namelist6"){
                        $nblist6++;
                        $list6=$list6." $namelist6";
                }
                else{
                        last;
                };
        };

	$list=$list6;
	@fgt=split(' ',$list);
	#print "mol$moli fgts= @fgt\n";

### GENERE LES COMBINAISONS

	#fill with single fragments
	@combin_num="";
	@fgtorigin="";
	$reference="";
	foreach $i (0..@fgt-1){
		$combin_num[$i]=" $i ";
		$reference=$reference." $i ";
		$fgtorigin[$i]=&readorigin($fgt[$i]);
		#print "\tfgt indice no $i => $combin_num[$i] with -$fgtorigin[$i]-\n";
	};
	$k=@fgt;
	$nbfgtmax=$k;
	#print "$nbfgtmax fragments\n";

#if($nbfgtmax < $nbmaxfgt){
	$countcombin=0;
	$fin=@fgt-1;
	foreach $i (0..$fin){
		$combini=$i;
		foreach $m (0..$fin){
		$combinj=$m;

		@pilei="";
		@pilej="";
		$pi=0;
		$depart=1;
		$depart=0 if($nbfgtmax==1);
		while(($pilei[0] ne "" || $depart==1) && $countcombin <= $nbmaxcombin){
			if($combinj <= $fin){
				$match=&matchorigin($fgtorigin[$combini],$fgtorigin[$combinj]);
				#print "test $combini ($fgtorigin[$combini]) $combinj ($fgtorigin[$combinj]) match=$match\n";
				if($combin_num[$combini]!~/ $combinj / && $match==1){
					#print "combinj=$combinj ($fgtorigin[$combini]) ; combin_num[$combini] ($fgtorigin[$combinj]) =$combin_num[$combini];\n";
					unlink "combin.sdf" if(-e "combin.sdf");
					chop($out3 =`perl $leaexe/lea3d-LINK_2MOL.pl $fgt[$combini] 0 $fgt[$combinj] 0 1`);
					if(-e "combin.sdf" && !-z "combin.sdf"){
						rename "combin.sdf", "combin_$k.sdf";
						$fgt[$k]="combin_$k.sdf";
						$combin_num[$k]="$combin_num[$combini] $combin_num[$combinj]";
						$fgtorigin[$k]=$fgtorigin[$combini]." ".$fgtorigin[$combinj];
						@count=split(' ',$combin_num[$k]);
						$nbcount=@count;
						last if($full==2 && $nbcount==$nbfgtmax);#stop once rebuild
						$pilej[$pi]=$combinj+1;
						$pilei[$pi]=$combini;
						#print "success add pile i=@pilei j=@pilej\n";
						$pi++;
						$depart=0;
						$combini=$k;
						$combinj=0;
						$k++;
						$countcombin++;
					}
					else{
						#print "pilej j=@pilej\n";
						$combinj++;
						#print "1- fail to link change combinj to $combinj\n";
					};
				}
				else{
					#print "pilej j=@pilej\n";
					$combinj++;
					#print "2- already tested change combinj to $combinj\n";
				};
			}
			elsif($pilei[0] ne ""){
				#print "avant pile i=@pilei j=@pilej\n";
				$combini=$pilei[$pi-1];
				$pilei[$pi-1]="";
				$combinj=$pilej[$pi-1];
				$pilej[$pi-1]="";
				$pi=$pi-1;
				#print "apres pile i=@pilei j=@pilej\n";
				#print "test combin_num no combini $combini et fgt combinj $combinj\n";
			}
			else{#$pilei[0] eq ""
				$depart=0;
			};
		};
		};#$m
	};	
	if($countcombin > $nbmaxcombin ){
		#print "Warning: stops after $countcombin tested combinations (maximum = $nbmaxcombin)\n";
	};

### ENLEVE LES DOUBLONS

	$countcombinunique=0;
	@combin_unique="";
	foreach $q (0..@combin_num-1){
		@range=split(' ',$combin_num[$q]);
		$nbrange=@range;
		#print "R @range\n";
		@range2=sort @range;
		#print "R2 @range2\n";
		$combin_tmp=join(' ',@range2);
		$vu=1;
		foreach $p (0..@combin_num-1){
			if($combin_unique[$p] ne ''){
				#$vu=0 if($combin_unique[$p] eq $combin_tmp || $combin_tmp eq $reference);
				$vu=0 if($combin_unique[$p] eq $combin_tmp);
			};
		};
		if($vu){
			if($nbrange == $nbfgtmax && $full){#in case of full reconstruction only
				#if($nbrange == $nbmaxcombin && $full){
				$combin_unique[$q]=$combin_tmp;
			}
			#elsif($nbrange < $nbmaxcombin){
			elsif($nbrange < $nbfgtmax){# in case of combinatory
				$combin_unique[$q]=$combin_tmp;
			};
			$countcombinunique++;	
		};
	};
	#print "$countcombinunique unique combinations (out of $countcombin (includes single fragments)) in the table combin_unique\n";

### AFFICHAGE SOLUTIONS

	$nbfilenew=0;
	foreach $i (0..@fgt-1){
		if($combin_unique[$i] ne '' ){
 	               	@nbfgt=split(' ',$combin_unique[$i]);
        	        $nbfgts=@nbfgt;
			if($full==2 && $nbfgts==$nbfgtmax){
				#print "write the complete structure only\n";
				$nbfilenew++;
				$lignesdf="";
				$lignesdf=$lignesdf."> <COMBINATION>\nFGT=$fgt[$i] NB_FGTS=$nbfgts\n";
				&printsdf($fgt[$i],$lignesdf,$filesdf);
			}
			elsif($full<2){
				if(($isolated==0 && $nbfgts > 1) || ($isolated && $nbfgts >0)){
					$nbfilenew++;
					$lignesdf="";
					$listname="";
					@listnametab=split(" ",$combin_unique[$i]);
					foreach $op (0..@listnametab-1){
						$listname=$listname." ".$fgt[$listnametab[$op]];
					};
					$lignesdf=$lignesdf."> <COMBINATION>\nFGT=$fgt[$i] NB_FGTS= $nbfgts LIST=$listname\n";
					&printsdf($fgt[$i],$lignesdf,$filesdf);
					#print "combination $fgt[$i] \t$nbfgts fragments\n";
				}
				else{
					#to remove file see below
					$combin_unique[$i]="";
				};	
			};
		};
	};

	if($nbfilenew==0){
		open(IN,">>toosmall.sdf");
		open(OP,"<$init2");
		while(<OP>){
			print IN $_;
		};
		close(OP);
		close(IN);
		#print "No combinations: does the input molecule too small? ($nbfgtmax fragments (full building option=$full) See toosmall.sdf\n";
	};	

#####################################################################################
### Nettoyage du repertoire ####

	unlink "init-tmp.sdf";

	$m=0;
	foreach $q (0..@combin_num-1){
		if($full<2){
			if($combin_unique[$q] eq ""){
				unlink "combin_$q.sdf";
			}
			else{
				$m++;
				#print "$m combin_$q.sdf $combin_unique[$q]\n";
				unlink "combin_$q.sdf";
			};
		}
		else{
			unlink "combin_$q.sdf";
		};
	};

	foreach $rmi (0..$nbmaxcombin){
                $namerm="combin_".$rmi.".sdf";
                if(-e $namerm){
                        unlink $namerm;
                };
        };
	unlink "control_fgt";

};#Loop over $nbinput

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

sub printsdf {

	local($fileopen,$phrase,$fileout)=@_;

	$flagnew=1;
	open(SDF,"<$fileopen");
	open(SDT,">>$fileout");
	while(<SDF>){
		$conv2=$_;
		if($flagnew){
			$flagnew=0;
			@ligne='';
			$si=0;
			$gosdf=1;
		};
		#comment following line to write info about points and origin
		#$gosdf=0 if ($conv2 =~/^>/);
		if($gosdf && $conv2 !~/(\$\$\$\$)/){
			$ligne[$si]=$_;
			$si++;
		};
		if ($conv2 =~/(\$\$\$\$)/){
			$flagnew=1;
			foreach $si (0..@ligne-1){
				print SDT $ligne[$si];
			};
		};
	};
	print SDT "$datablock" if($ifdatablock==1 && $datablock ne "");
	print SDT "$phrase\n";
	print SDT "\$\$\$\$\n";
	close(SDT);
	close(SDF);
};

##############
sub readorigin {
	
	local($fileopen)=@_;

        $flagnew=1;
	$no=0;
        open(IN,"<$fileopen");
        while(<IN>){
                if($flagnew){
                        $flagname=0;
                        $flagname2=0;
                        $no++;
                        $ig=1;
                        $istratom=0;
                        $points="";
                        $origins="";
                        $compt=0;
                        $flagnew=0;
                };
                @getstr = split(' ',$_);
                $compt++;
                if (($compt > 4) && ($ig <= $istratom)){
                        $atom=$getstr[3];
                        $ig++;
                };
                if ($compt == 4){
                        $istratom=$getstr[0];
                        $istrbond=$getstr[1];

                        @coller=split(' *',$istratom);
                        if(@coller>3 && @coller==6){
                                $istratom=$coller[0].$coller[1].$coller[2];
                                $istrbond=$coller[3].$coller[4].$coller[5];
                        }
                        elsif(@coller>3 && @coller==5){
                                if($_=~/^\s/){
                                        $istratom=$coller[0].$coller[1];
                                        $istrbond=$coller[2].$coller[3].$coller[4];
                                }
                                else{
                                        $istratom=$coller[0].$coller[1].$coller[2];
                                        $istrbond=$coller[3].$coller[4];
                                };
                        };
                };
                if ($_=~/\$\$\$\$/){
                        $flagnew=1;
                        $points=~s/-/ /g;
                        $points=" ".$points." ";
                        $origins=~s/-/ /g;
                        $origins=" ".$origins." ";
		};
                if($flagname && $getstr[0] ne ''){
                        $points=$getstr[0];
                        $flagname=0;
                };
                if($flagname2 && $getstr[0] ne ''){
                        $origins=$getstr[0];
                        $flagname2=0;
                };
                if($_=~/^>/ && $_=~/POINTS/){
                        $flagname=1;
                };
                if($_=~/^>/ && $_=~/ORIGIN/){
                        $flagname2=1;
                };
        };
        close(IN);
	
	$origins;	
};

########################################

sub matchorigin {

	local($a,$b)=@_;
	
	$meme=0;
	@gata=split(' ',$a);
	foreach $d (0..@gata-1){
		@gat=split('_',$gata[$d]);
		$inverse=$gat[1]."_".$gat[0];
		if($b=~/ $gata[$d] / || $b=~/ $inverse /){
			$meme=1;
		};
	};
	$meme;
};

########################################
