#!/usr/bin/perl

	local($filesdf,$optionsplit)= @ARGV;

	$leaexe=$0;
	#windows requires 2 steps
	$leaexe=~s/lea3d-MAKE_FGTS_AGGREG\.pl$//;
	$leaexe=~s/\/$//;
	if($leaexe eq ""){
                $leaexe=".";
        };
	#print "perl scripts in $leaexe\n";

	if($filesdf eq ''){
		die "file ?\nusage: <>.pl <file.sdf with only one molecule> <option split acyclic (default) or not '0'> creates aggregated fgt if small substituents on a ring and split acyclic fgt into 2 smaller fgts; output file make_fgts_aggreg.sdf \n";
	};

	# used in sensaas-bioisostere $optionsplit=0
        if($optionsplit eq ""){
                $optionsplit=1;
        }
        else{
                $optionsplit=0;
        };

	#nbheavymax defines which size of substituent to aggregate 
	#if 1 (CH3, Cl, Br, H, I, F, OH) ; if 2 (CH3CH3, OCH3, ...)
	$nbheavymax=2;	

	unlink "make_fgts.sdf" if(-e "make_fgts.sdf");
	unlink "ring.sdf" if(-e "ring.sdf");
	unlink "fused_rings.sdf" if(-e "fused_rings.sdf");
	unlink "special.sdf" if(-e "special.sdf");
	unlink "substituent.sdf" if(-e "substituent.sdf");
       	unlink "acyclic.sdf" if(-e "acyclic.sdf");
	unlink "linker.sdf" if(-e "linker.sdf");

	if($optionsplit==0){
		chop($log=` perl $leaexe/lea3d-MAKE_FGTS.pl $filesdf `);
	}
	else{
		#create fgts then split large acyclic fgts
		&make_fgts_split($filesdf);
		if(-e "make_fgts_split.sdf" && !-z "make_fgts_split.sdf"){
			rename "make_fgts_split.sdf", "make_fgts.sdf";
		};
	
	};
	#clean
	unlink "make_fgts_split.sdf" if(-e "make_fgts_split.sdf");
	unlink "acyclic-parts.sdf" if(-e "acyclic-parts.sdf");
	#read ring.sdf fused_rings.sdf special.sdf linker.sdf substituent.sdf
	unlink "ring.sdf" if(-e "ring.sdf");
	unlink "fused_rings.sdf" if(-e "fused_rings.sdf");
	unlink "special.sdf" if(-e "special.sdf");
	unlink "linker.sdf" if(-e "linker.sdf");
	unlink "substituent.sdf" if(-e "substituent.sdf");
	unlink "acyclic.sdf" if(-e "acyclic.sdf");

	$nb2=&nbsdf("make_fgts.sdf");	
	#print "$nb2 fragments in make_fgts.sdf \n";

if($nb2 > 1){

	#read make_fgts.sdf and combin substituents to rings if nb heavy atom <= $nbheavymax
	@name="";#name of <fgt>
	@point="";#atom connected to X
	@nbpoints="";#nb X
	@heavy="";
	@origin="";#origin datablock of atom connected to X
	$flagnew=1;
	$no=0;
	open(IN,"<make_fgts.sdf");
	while(<IN>){
		if($flagnew){
			$flagname=0;
			$flagname2=0;
			$no++;
			$ig=1;
			$istratom=0;
			$nbheavy=0;
			$points="";
			$origins="";
			$compt=0;
			$flagnew=0;
			$line="";
			$fgt="make_fgts_".$no.".sdf";
		};
		$line=$line.$_;
		@getstr = split(' ',$_);
		$compt++;
		if (($compt > 4) && ($ig <= $istratom)){
			$atom=$getstr[3];
			$nbheavy++ if($atom ne 'H' && $atom ne 'X' && $atom ne 'Cl' && $atom ne 'Br' && $atom ne 'I' && $atom ne 'F');
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
			@get=split('-',$points);
			$nbp=@get;
			$points=~s/-/ /g;
			$points=" ".$points." ";
			$origins=~s/-/ /g;
			$origins=" ".$origins." ";
			#@get1=split('-',$origins);
			$name[$no]=$fgt;
			$nbpoints[$no]=$nbp;
			$point[$no]=$points;
			$origin[$no]=$origins;
			@heavy[$no]=$nbheavy;
			open(OUT,">$fgt");
			print OUT "$line";
			close(OUT);
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

	foreach $i (1..$nb2){
		if($heavy[$i] <= $nbheavymax && $nbpoints[$i]==1){
			@get=split(' ',$origin[$i]);
			@get2=split('_',$get[0]);
			$inverse=$get2[1]."_".$get2[0];
			foreach $j (1..$nb2){
				if($j!=$i && ($origin[$j]=~/ $get[0] / || $origin[$j]=~/ $inverse /)){
					chop($log=` perl $leaexe/lea3d-LINK_2MOL.pl $name[$j] 0 $name[$i] 0 1`);
					if(-e "combin.sdf" && !-z "combin.sdf"){
						rename "combin.sdf", $name[$j];
						unlink "$name[$i]";
						$name[$i]="";
					};
				};
			};
		};
	};
	open(OUT,">make_fgts_aggreg.sdf");
	foreach $i (1..$nb2){
		if($name[$i] ne ""){
			open(IN,"<$name[$i]");
			while(<IN>){
				print OUT $_;
			};
			close(IN);
			unlink $name[$i];
		};
	};
	close(OUT);

};

###################################################################################
###################################################################################

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

#usage: sub () <file.sdf> creates fgts and additional fragments when linkers and substituents are large (acyclic.sdf=substituent.sdf+linker.sdf)
sub make_fgts_split{

	local($fsdf)=@ARGV;

	#rule 1: if can cut bond between C and atom !=H
	#rule 2: $nbbondsmin defines which size of acyclic to split 
	$nbbondsmin=4;	

	chop($log=` perl $leaexe/lea3d-MAKE_FGTS.pl $fsdf `); 

	#make_fgts.sdf = ring.sdf fused_rings.sdf special.sdf linker.sdf substituent.sdf
	#if molecule without ring then acyclic.sdf=molecule
	unlink "acyclic-parts.sdf" if(-e "acyclic-parts.sdf");
	unlink "make_fgts_split.sdf" if(-e "make_fgts_split.sdf");

	open(OUT,">make_fgts_split.sdf");
	if(-e "ring.sdf" && !-z "ring.sdf"){
		open(IN,"<ring.sdf");
		while(<IN>){
			print OUT $_;
		};
		close(IN);
	};
	if(-e "fused_rings.sdf" && !-z "fused_rings.sdf"){
		open(IN,"<fused_rings.sdf");
		while(<IN>){
			print OUT $_;
		};
		close(IN);
	};
	if(-e "special.sdf" && !-z "special.sdf"){
		open(IN,"<special.sdf");
		while(<IN>){
			print OUT $_;
		};
		close(IN);
	};
	close(OUT);
	#Keep only acyclic.sdf
	unlink "substituent.sdf" if(-e "substituent.sdf");
	unlink "linker.sdf" if(-e "linker.sdf");
	unlink "ring.sdf" if(-e "ring.sdf");
	unlink "fused_rings.sdf" if(-e "fused_rings.sdf");
	unlink "special.sdf" if(-e "special.sdf");

	if(!-z "acyclic.sdf" && -e "acyclic.sdf"){
		#check size of each and split into small parts before cat with make_fgts_split.sdf
		#read each acyclic
		$flagnew=1;
		$moli=0;
		open(IN,"<acyclic.sdf");
		while(<IN>){
			if($flagnew){
				$moli++;
				$compt=0;
                        	$flagnew=0;
				$ig=1;
				$jg=0;
				@strx='';
				@stry='';
				@strz='';
				@atom='';
	                        @coval='';
        	                @fonc='';
                	        @ifonc='';
                        	@covfonc='';
                        	@bond='';
                        	@listb='';
                        	@typeb='';
                        	$blanc=' ';
                        	@lignebond='';
				$cut=0;
				$bondcut="";
				@sdf='';
				$nbheavy=0;
			};
			@getstr = split(' ',$_);
			$sdf[$compt]=$_;
			$compt++;
			if (($compt > 4) && ($ig <= $istratom)){
				$strx[$ig]=$getstr[0];
				$stry[$ig]=$getstr[1];
				$strz[$ig]=$getstr[2];
				$atom[$ig]=$getstr[3];
				$nbheavy++ if($atom[$ig] ne 'H' && $atom[$ig] ne 'X' && $atom[$ig] ne 'Cl' && $atom[$ig] ne 'Br' && $atom[$ig] ne 'I' && $atom[$ig] ne 'F');
				$ig++;
			};
			if (($compt > 4) && ($ig > $istratom) && ($jg <=$istrbond)){
				if ($jg == 0){
                                	$jg++;
                        	}
                        	else{
                                	@coller=split(' *',$getstr[0]);
                                	@coller2=split(' *',$getstr[1]);
                                	if(@coller==6 && $getstr[1] ne ""){
                                        	$getstr[0]=$coller[0].$coller[1].$coller[2];
                                        	$getstr[2]=$getstr[1];
                                        	$getstr[1]=$coller[3].$coller[4].$coller[5];
                                	}
                                	elsif(@coller==6 && $getstr[1] eq ""){
                                        	$getstr[0]=$coller[0].$coller[1];
                                        	$getstr[1]=$coller[2].$coller[3].$coller[4];
                                        	$getstr[2]=$coller[5];
                                	}
                                	elsif(@coller==5){
                                        	if($_=~/^\s/){
                                                	$getstr[0]=$coller[0].$coller[1];
                                                	$getstr[2]=$getstr[1];
                                                	$getstr[1]=$coller[2].$coller[3].$coller[4];
                                        	}
                                        	else{
                                                	$getstr[0]=$coller[0].$coller[1].$coller[2];
                                                	$getstr[2]=$getstr[1];
                                                	$getstr[1]=$coller[3].$coller[4];
                                        	};
                                	}
                                	elsif(@coller==4){
                                        	if($_=~/^\s/){
                                                	$getstr[0]=$coller[0];
                                                	$getstr[2]=$getstr[1];
                                                	$getstr[1]=$coller[1].$coller[2].$coller[3];
                                        	}
                                        	else{
                                                	$getstr[0]=$coller[0].$coller[1].$coller[2];
                                                	$getstr[2]=$getstr[1];
                                                	$getstr[1]=$coller[3];
                                        	};
                               	 	}
                                	elsif(@coller2==4){
                                        	$getstr[1]=$coller2[0].$coller2[1].$coller2[2];
                                        	$getstr[2]=$coller2[3];
                                	}
                                	elsif(@coller==7){
                                        	$getstr[0]=$coller[0].$coller[1].$coller[2];
                                        	$getstr[1]=$coller[3].$coller[4].$coller[5];
                                        	$getstr[2]=$coller[6];
                                	};

					if($getstr[2]==1 && ($atom[$getstr[1]] eq 'C' || $atom[$getstr[0]] eq 'C') && $atom[$getstr[1]] ne 'H' && $atom[$getstr[0]] ne 'H' && $atom[$getstr[1]] ne 'X' && $atom[$getstr[0]] ne 'X'){
						$cut++;
						$compt2=$compt-1;
						$bondcut=$bondcut." ".$compt2;#to break if not printed from @line
						#print "bond $jg\n";
					};

					$bond[$getstr[0]]=$bond[$getstr[0]].$blanc.$getstr[1].$blanc.$getstr[2];
                                	$listb[$getstr[0]]=$listb[$getstr[0]].$blanc.$getstr[1];
                                	$typeb[$getstr[0]]=$typeb[$getstr[0]].$blanc.$getstr[2];

                                	$bond[$getstr[1]]=$bond[$getstr[1]].$blanc.$getstr[0].$blanc.$getstr[2];
                                	$listb[$getstr[1]]=$listb[$getstr[1]].$blanc.$getstr[0];
                                	$typeb[$getstr[1]]=$typeb[$getstr[1]].$blanc.$getstr[2];

                                	$fonc[$getstr[0]]=$fonc[$getstr[0]].$blanc.$getstr[2].'-'.$atom[$getstr[1]].$blanc;
                                	$ifonc[$getstr[0]]=$ifonc[$getstr[0]].$blanc.$getstr[1].$blanc;
                                	$covfonc[$getstr[0]]=$covfonc[$getstr[0]].$blanc.$getstr[2];
                                	$coval[$getstr[0]]=$coval[$getstr[0]]+$getstr[2];

                                	$fonc[$getstr[1]]=$fonc[$getstr[1]].$blanc.$getstr[2].'-'.$atom[$getstr[0]].$blanc;
                                	$ifonc[$getstr[1]]=$ifonc[$getstr[1]].$blanc.$getstr[0].$blanc;
                                	$covfonc[$getstr[1]]=$covfonc[$getstr[1]].$blanc.$getstr[2];
                                	$coval[$getstr[1]]=$coval[$getstr[1]]+$getstr[2];
                                	$lignebond[$jg]=$_;
                                	$jg++;
                        	};
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
				if($cut >= $nbbondsmin){
					#print "acyclic.sdf no $moli has $cut single bonds that can be splitted\n";
					&splitacyclic;
				}
				else{#print in make_fgts_split.sdf
					open(OUT,">>make_fgts_split.sdf");
					foreach $i (0..@sdf-1){
						print OUT "$sdf[$i]";
					};
					close(OUT);
				};
			};
		};
		close(IN);

	};
	if(-e "acyclic-parts.sdf" && !-z "acyclic-parts.sdf"){
		open(OUT,">>make_fgts_split.sdf");
		open(IN,"<acyclic-parts.sdf");
		while(<IN>){
			print OUT $_;
		};
		close(IN);
		close(OUT);
	};
};

###################################################################################
###################################################################################
	
sub splitacyclic{
	
	@listcut=split(' ',$bondcut);
	#$half=int($istratom/2);
	#print "half=$half\n";
	@sizefgt1="";
	@sizefgt2="";
	@suite1="";
	@tabifx="";
	@tabpoint="";
	$diffsize=1000;
	$diffsizei=0;
	foreach $j (0..@listcut-1){
		@ij=split(' ',$sdf[$listcut[$j]]);
		#print "cut between @ij ?\n";
		#try each bond
		$one=0;
		$suite=" 1 ";
		$continue=1;
		$bi=1;
		$pos=0;
		$ifx="";
		$point="";
		while($continue){
			@get=split(' ',$ifonc[$bi]);
	                foreach $k (0..@get-1){
				#exclude bond $j
				if(($bi==$ij[0] && $get[$k]==$ij[1]) || ($bi==$ij[1] && $get[$k]==$ij[0])){
					$ifx=$get[$k];
					$point=$bi;
				}
				else{
        	                	$suite=$suite."$get[$k] " if ($suite!~/ $get[$k] /);
				};
                	};
                	@get2=split(' ',$suite);
                	$pos++;
                	if($pos == @get2){
                        	 $continue=0;
                	}
                	else{
                        	$bi=$get2[$pos];
                	};
        	};
		@get2=split(' ',$suite);
		$heavy1=0;
		foreach $k (0..@get2-1){
			$heavy1++ if($atom[$get2[$k]] ne 'H' && $atom[$get2[$k]] ne 'X' && $atom[$get2[$k]] ne 'F' && $atom[$get2[$k]] ne 'Cl' && $atom[$get2[$k]] ne 'I' && $atom[$get2[$k]] ne 'Br');
		};
		$sizefgt1[$j]=$heavy1;
		$sizefgt2[$j]=$nbheavy-$heavy1;
		$diff=sqrt(($sizefgt1[$j]-$sizefgt2[$j])*($sizefgt1[$j]-$sizefgt2[$j]));
		#print "$j $sizefgt1[$j] $sizefgt2[$j] diff=$diff\n";
		if($diffsize > $diff){
			$diffsize=$diff;
			$diffsizei=$j;
		};
		$suite1[$j]=$suite;
		$tabifx[$j]=$ifx;
		$tabpoint[$j]=$point;
	};
	#print "Best cut at no $diffsizei size1=$sizefgt1[$diffsizei] size2=$sizefgt2[$diffsizei]\n";
	@get2=split(' ',$suite1[$diffsizei]);
	$nbacc=@get2;
	$acc=$suite1[$diffsizei];
	@ij=split(' ',$sdf[$listcut[$diffsizei]]);
	$originij=$ij[0]."_".$ij[1];
        #check ORIGIN tag - it must be different from what is already in make_fgts.sdf for other fgts
	$vu=0;
	open(OUP,"<make_fgts.sdf");
	while(<OUP>){
		if($_=~/$originij/){
			$vu=1;
		};
	};
	close(OUP);
	if($vu==1){
        	#change
                $originij=($ij[0]+100)."_".($ij[1]+100);
                #print "create new origin tag here $originij\n";
	}
	else{
		$originijrev=$ij[1]."_".$ij[0];
		$vu1=0;
		open(OUP,"<make_fgts.sdf");
		while(<OUP>){
			if($_=~/$originijrev/){
				$vu1=1;
			};
		};
		close(OUP);
		if($vu1==1){
			#change
			$originij=($ij[0]+100)."_".($ij[1]+100);
			#print "create new origin tag here $originij\n";
		};
	};
	&printpart($suite1[$diffsizei],$tabifx[$diffsizei],$tabpoint[$diffsizei]);

	#print second part
	$continue=1;
	while($continue){
		$suite=" ";
		$bi='';
		foreach $k (1..$istratom){
			$bi=$k if($acc!~/ $k /);
		};
		$suite=$suite."$bi ";
		$cont=1;
                $pos=0;
		$ifx="";
		$point="";
                while($cont){
                	@get=split(' ',$ifonc[$bi]);
   	             	foreach $k (0..@get-1){
				if(($bi==$ij[0] && $get[$k]==$ij[1]) || ($bi==$ij[1] && $get[$k]==$ij[0])){
					$ifx=$get[$k];
					$point=$bi;
				}
				else{
                        		$suite=$suite."$get[$k] " if ($suite!~/ $get[$k] /);
				};
                        };
                        @get2=split(' ',$suite);
                        $pos++;
                        if($pos == @get2){
                        	$cont=0;
                        }
                        else{
                        	$bi=$get2[$pos];
                        };
		};
		$acc=$acc.' '.$suite;#allow to print more salts until @getd >= $istratom
		@getd=split(' ',$acc);
		$continue=0 if(@getd >= $istratom);
		$continue=0 if($bi eq '');
		&printpart($suite,$ifx,$point);
	};
};

###################################################################################
###################################################################################

sub printpart{
	local($atomes,$ifx,$point)=@_;

	$atomes=' '.$atomes.' ';
	@getp=split(' ',$atomes);
	$longa=@getp;
	$longb=0;
	$gh=1;
	@ligneb='';
	foreach $f (($istratom+5-1)..@sdf-1){ # 5-1 to begin index 0 see $compt++
		@getb=split(' ',$sdf[$f]);
		@coller=split(' *',$getb[0]);
		@coller2=split(' *',$getb[1]);
		if(@coller==6 && $getb[1] ne ""){
			$getb[0]=$coller[0].$coller[1].$coller[2];
			$getb[2]=$getb[1];
			$getb[1]=$coller[3].$coller[4].$coller[5];
		}
		elsif(@coller==6 && $getb[1] eq ""){
			$getb[0]=$coller[0].$coller[1];
			$getb[1]=$coller[2].$coller[3].$coller[4];
                        $getb[2]=$coller[5];
		}
		elsif(@coller==5){
	                if($sdf[$f]=~/^\s/){
        	                $getb[0]=$coller[0].$coller[1];
                                $getb[2]=$getb[1];
                                $getb[1]=$coller[2].$coller[3].$coller[4];
                	}
                	else{
                		$getb[0]=$coller[0].$coller[1].$coller[2];
                                $getb[2]=$getb[1];
                                $getb[1]=$coller[3].$coller[4];
                      	};
            	}
		elsif(@coller==4){
			if($sdf[$f]=~/^\s/){
	                        $getb[0]=$coller[0];
                                $getb[2]=$getb[1];
                                $getb[1]=$coller[1].$coller[2].$coller[3];
        		}
                	else{
                		$getb[0]=$coller[0].$coller[1].$coller[2];
                                $getb[2]=$getb[1];
                                $getb[1]=$coller[3];
           		};
          	}
                elsif(@coller2==4){
                		$getb[1]=$coller2[0].$coller2[1].$coller2[2];
                                $getb[2]=$coller2[3];
              	}
		elsif(@coller==7){
			$getb[0]=$coller[0].$coller[1].$coller[2];
                        $getb[1]=$coller[3].$coller[4].$coller[5];
                 	$getb[2]=$coller[6];
		};
		$gh=0 if($sdf[$f]=~/^>/ || $getb[0] eq '' || $sdf[$f]=~/^M/ || $sdf[$f]=~/^\$\$\$\$/);
		if($gh){
			if($atomes=~/ $getb[0] / || $atomes=~/ $getb[1] /){
				$a1='';
				$a2='';
				foreach $p (0..@getp-1){
					$a1=($p+1) if($getb[0] eq $getp[$p]);
					$a2=($p+1) if($getb[1] eq $getp[$p]);
					$a1=$longa+1 if($getb[0] eq $ifx);
					$a2=$longa+1 if($getb[1] eq $ifx);
				};
				$z="  0  0  0  0";
				$ligneb[$longb]=sprintf"%3s%3s%3s  0  0  0  0",$a1,$a2,$getb[2];
				$longb++;
			};
		};
	};
	open(OUT,">>acyclic-parts.sdf");
	foreach $p (0..2){
		print OUT "$sdf[$p]";
	};
	$longa=$longa+1;#add -X ($ifx)
	printf OUT "%3s%3s  0  0  0  0  0  0  0  0999 V2000\n",$longa,$longb;
	$pi=0;
	@newnumber='';
	foreach $f (0..@getp-1){
		$pi++;
		print OUT "$sdf[$getp[$f]+3]";
		$newnumber[$getp[$f]]=$pi;
		$origini=$pi if($point==$getp[$f]);
	};
	#change connected atom $ifx in X
	@get6=split(' ',$sdf[$ifx+3]);
	$get6[3]='X';
	$atomx=sprintf"%10s%10s%10s%1s%1s%1s%3s%3s%3s\n",$get6[0],$get6[1],$get6[2],$blanc,$get6[3],$blanc,$get6[4],$get6[5],$get6[6];
	print OUT "$atomx";
	foreach $f (0..@ligneb-1){
		print OUT "$ligneb[$f]\n";
	};
	print OUT "M  END\n";
#> <POINTS>
#2-3-3-5
#
#> <ORIGIN>
#9_10-6_7-6_8-4_3
	$ecritfin=0;
	$flagpoint=0;
	$flagorigin=0;
	$flagvu=0;
	@skiporigin="";
	foreach $mp (0..@sdf-1){
		$ecritfin=1 if($sdf[$mp]=~/^>/);
		if($flagpoint){
			@get6=split(' ',$sdf[$mp]);
			#replace old number by new
			@get7=split('-',$get6[0]);
			$tmppoint='';
			foreach $k (0..@get7-1){
				if($newnumber[$get7[$k]]ne ''){#is empty if in the other acyclic part: to skip
					if($tmppoint eq ''){
						$tmppoint=$newnumber[$get7[$k]];
					}
					else{
						$tmppoint=$tmppoint.'-'.$newnumber[$get7[$k]];
					};
				}
				else{#do not print ORIGIN linked to this -X (is in the other part)
					$skiporigin[$k]=1;
				};
			};
			$tmppoint=$tmppoint."-".$origini;
			$tmppoint=~s/ //g;
			$tmppoint=~s/--/-/g;
			$tmppoint=~s/^-//;
			$tmppoint=~s/-$//;
			print OUT "$tmppoint\n";
			$flagpoint=0;
		}
		elsif($flagorigin){
			@get6=split(' ',$sdf[$mp]);
			#keep only if in this part
			@get7=split('-',$get6[0]);
			$tmporigin='';
			foreach $k (0..@get7-1){
				if($skiporigin[$k]!=1){#print if not ==1
					if($tmporigin eq ''){
						$tmporigin=$get7[$k];
					}
					else{
						$tmporigin=$tmporigin.'-'.$get7[$k];
					};
				};
			};
			$tmporigin=$tmporigin.'-'.$originij;
			$tmporigin=~s/ //g;
			$tmporigin=~s/--/-/g;
			$tmporigin=~s/^-//;
			$tmporigin=~s/-$//;
			print OUT "$tmporigin\n";
			$flagorigin=0;
		}
		elsif($sdf[$mp]=~/\$\$\$\$/ && $flagvu==0){
			$origini=~s/ //g;
			$origini=~s/--/-/g;
			$origini=~s/^-//;
			$origini=~s/-$//;
			$originij=~s/ //g;
			$originij=~s/--/-/g;
			$originij=~s/^-//;
			$originij=~s/-$//;
			print OUT "> <POINTS>\n";
                        print OUT "$origini\n";
                        print OUT "\n";
                        print OUT "> <ORIGIN>\n";
                        print OUT "$originij\n";
                        print OUT "\n";
                        print OUT "\$\$\$\$\n";
		}
		elsif($ecritfin){
			print OUT "$sdf[$mp]";
		};
		if($sdf[$mp]=~/^>/ && $sdf[$mp]=~/<POINTS>/){
			$flagpoint=1;
			$flagvu=1;
		};
		if($sdf[$mp]=~/^>/ && $sdf[$mp]=~/<ORIGIN>/){
			$flagorigin=1;
			$flagvu=1;
		};
	};
	print OUT "\$\$\$\$\n" if($ecritfin==0);
	close(OUT);
};

######################################################################
######################################################################
