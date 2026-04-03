#!/usr/bin/perl

	my($filesdf,$filelistatom)=@ARGV;

	$leaexe=$0;
	#windows requires 2 steps
	$leaexe=~s/lea3d-MAKE_FGTS_LIB\.pl$//;
	$leaexe=~s/\/$//;
	#print "perl scripts in $leaexe\n";
	
	if($filesdf eq ''){
                print "usage: \n lea3d-MAKE_FGTS_LIB.pl <file sdf> <file with list atom> \n";
                exit;
        };

	#IMPORTANT
	# $filesdf is in 3D
	# be sure that the previous step either select the 3D coordinates file or generate a conformer if it comes from the sketcher

	@lib="";
        $libnofgt="";
        $li=0;

	if($filelistatom eq ""){
		#use makefgtagg with the option or not to split acyclics in smaller fgt
		#suggest list at next step
		# option '1' in order to split acyclics:
		chop($log=`perl $leaexe/lea3d-MAKE_FGTS_AGGREG.pl $filesdf 1 `);#modify as .pl

		#output file make_fgts_aggreg.sdf
                unlink "make_fgts.sdf" if(-e "make_fgts.sdf");

		#suggest other libraries with combination of 2
		#read make_fgts_aggreg.sdf 
		#if 1 -X then nothing
		#if 2 or more -X then combinations

	}
	else{
		#read $filelistatom
		#makefgtagg without splitting acyclics (because the datablock LISTATOM cannot be correctly splitted - info is lost)
		#match and try to recombin if required
	
		@listatom="";
		open(IN,"<$filelistatom");
		while(<IN>){
			if($_!=/^#/){
				@get=split(' ',$_);
				foreach $i (0..@get-1){
					$listatom[$i]=$get[$i] if($get[$i]=~/\d/);
				};
			};
		};
		close(IN);
		#print "list to match: @listatom\n";

		# option '0' for not to split acyclics:
		chop($log=`perl $leaexe/lea3d-MAKE_FGTS_AGGREG.pl $filesdf 0 `);#modify iy as .pl

		#output file make_fgts_aggreg.sdf
		unlink "make_fgts.sdf" if(-e "make_fgts.sdf");


		####### Retrieve information on fragments
		
		chop($log=`python $leaexe/listsdf.py make_fgts_aggreg.sdf ORIGIN `);
		unlink "datablock.txt";
		@origin="";
		$o=1;
		@get=split('\n',$log);
		foreach $i (0..@get-1){
			@get2=split(' ',$get[$i]);
			if($get2[0] ne ""){
				$origin[$o]=$get2[0];
				#print "fgt $o - $origin[$o]\n";
				$o++;
			};
		};
		$listcombin="";
		$listsmallpart="";
		@lista="";
		$o=1;
		chop($log=`python $leaexe/listsdf.py make_fgts_aggreg.sdf LISTATOM `);
		unlink "datablock.txt";
		@get=split('\n',$log);
		@match="";
		foreach $i (0..@get-1){
			@get2=split(' ',$get[$i]);
			$long=@get2;
			if($get2[0] ne ""){
				$m=join(' ',@get2);
				$m=" ".$m." ";
				$lista[$o]=$m;
				#print "fgt $o - $lista[$o] / $origin[$o]\n";

				#check if match @listatom
				$p=0;
				foreach $k (0..@listatom-1){
					if($m=~/ $listatom[$k] /){
						$p++;
						#print "same atom in different fgt?\n" if($match[$k] ne "");
						$match[$k]=$o;
					};
				};
				$ratio=$p/$long;
				if($p==$long){
					#print "fgt $o full match listatom\n";
					$listcombin=$listcombin." ".$o;
				}
				elsif($ratio > 0.5){
					#print "fgt $o has $p atoms in $o but not all (ie $long; ratio = $ratio) is selected\n";
					$listcombin=$listcombin." ".$o;
				}
				elsif($ratio<=0.5){
					#print "fgt $o has $p atoms in $o less than 50% but used if no fgt in $listcombin \n";
					$listsmallpart=$o;
				};
				$o++;
			};
		};
		#print "@match\n";
		if($listcombin eq ""){
			$listcombin=$listsmallpart;
		};
		@get2=split(' ',$listcombin);
		$nb=@get2;
		#print "$nb fragments match the list of atoms\n";	

		#######################################

		$newlib=0;
		if($nb==1){ 
			$lib[$li]="make_fgts_aggreg.sdf";
                        $libnofgt[$li]=$get2[0];
                        $li++;

		}
		elsif($nb>1){
			#print "$nb fgt full match to combin - fgts: $listcombin\n";
			$nb2=&nbsdf("make_fgts_aggreg.sdf");
			#print "$nb2 fragments in make_fgts_aggreg.sdf\n";
			chop($log=`python $leaexe/splitsdf.py make_fgts_aggreg.sdf libfgt `);

			#link 2 first fgts
			$name1="libfgt_".$get2[0].".sdf";
			$name2="libfgt_".$get2[1].".sdf";
			if($nb>=3){
				$name3="libfgt_".$get2[2].".sdf";
			};
			#first link fgt 1 to 2
			chop($log =`perl $leaexe/lea3d-LINK_2MOL.pl $name1 0 $name2 0 1 `);
			if(-e "combin.sdf" && !-z "combin.sdf"){
				$newlib=1;
				rename "combin.sdf", $name1;
				unlink "$name2";
				#if 3 fgts
				if($nb>=3){
					chop($log =`perl $leaexe/lea3d-LINK_2MOL.pl $name1 0 $name3 0 1 `);
					if(-e "combin.sdf" && !-z "combin.sdf"){
						rename "combin.sdf", $name1;
						unlink "$name3";
					};
				};
			}
			elsif($nb>=3){# name1 is not directly connected to name2: try name3
				chop($log =`perl $leaexe/lea3d-LINK_2MOL.pl $name1 0 $name3 0 1 `);
				if(-e "combin.sdf" && !-z "combin.sdf"){
					$newlib=1;
					rename "combin.sdf", $name1;
					unlink "$name3";
					#try to link fgt 2
					chop($log =`perl $leaexe/lea3d-LINK_2MOL.pl $name1 0 $name2 0 1 `);
					if(-e "combin.sdf" && !-z "combin.sdf"){
						rename "combin.sdf", $name1;
						unlink "$name2";
					};
				};
			}
			else{
				print "It not possible to create a combination with fragments $get2[0] and $get2[1]\n";
			};
			if($newlib){
				#create the new lib file
                                unlink "atmp.sdf" if(-e "atmp.sdf");
                                $nb3=0;
                                foreach $i (1..$nb2){
                                	$namei="libfgt_".$i.".sdf";
                                        if(-e "$namei" && !-z "$namei"){
                                        	$nb3++;
                                                open(OUT,">>atmp.sdf");
                                                open(IN,"<$namei");
                                                while(<IN>){
                                                	print OUT $_;
                                                };
                                                close(IN);
                                                close(OUT);
                                        };
                                        unlink "$namei" if(-e "$namei");
                                };
                                rename "atmp.sdf", "lib1.sdf";
				print "lib1.sdf contains $nb3 fgts\n";
                                $lib[$li]="lib1.sdf";
                                $libnofgt[$li]=$get2[0];
                                $li++;
			}
			else{
				#use the first identified fragment;
				$lib[$li]="make_fgts_aggreg.sdf";
		                $libnofgt[$li]=$get2[0];
                		$li++;
			};
		};

	};

        #print to suggest list at next step
        # +  suggest to select the entire molecule $filesdf?
	# 2- make_fgts_aggreg.sdf
	# 3- lib1.sdf
	# 4- ...
	
	print "list of atoms to match: @listatom in $filesdf (fragmented in make_fgts_aggreg.sdf)\n";
	print "$nb fragments match the list of atoms\n";
        if($nb>=3){
 	print "The number of fragments that match the list of atoms is $nb (fragments number @get2 in make_fgts_aggreg.sdf). However, only the 3 first will be connected if it is possible\n";
        };

	print "Results:\n";
	if($lib[0] ne ""){
		#print "\nInput list of atoms matches:\n";
		#print "- $filesdf (here the entire molecule will be used to find fragments having a similar shape\n";
		#print "- make_fgts_aggreg.sdf and you will select the fragment you want to replace at the next step\n";
		foreach $i (0..@lib-1){
			if(-e "$lib[$i]" && !-z "$lib[$i]" && $libnofgt[$i] ne ""){
				print "$lib[$i] fgt $libnofgt[$i] \n";
			};
		};
	};
	#print "\nUse make_fgts_aggreg.sdf and you will select the fragment you want to replace at the next step\n";
	print "make_fgts_aggreg.sdf \n";

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

