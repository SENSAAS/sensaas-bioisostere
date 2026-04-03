#!/usr/bin/perl

	my($f,$sdf,$ref,$nbmatchx)=@ARGV;
	if($f eq '' || $sdf eq '' || $ref eq ''){
        	die "usage: .pl <matrix-sensaas.txt> <catsensaas.sdf> <reference.sdf file for identifying X positions> <(option) number of required corresponding X (default = count the number of X in reference.sdf)>\n<>.pl matrix-sensaas.txt catsensaas.sdf reference.sdf ";
	};

	$leaexe=$0;
	#windows requires 2 steps
	$leaexe=~s/lea3d-ordered-catsensaas-filtered\.pl$//;
	$leaexe=~s/\/$//;
	#print "perl scripts in $leaexe\n";
	
	$distmax=1.0;
	
	$nbmatchx=&nbax($ref);
	if($nbmatchx > 1){#permissive in those cases
		#$distmax=3.0;
		$distmax=3.5;
	};

	print "inputs are $f and $sdf and matching $nbmatchx X dummy atoms with reference $ref (maximum distance = $distmax)\n";

	$j=0;
        @data="";
        open(IN,"<$f");
        while(<IN>){
                @get=split(' ',$_);
		if(@get==1){
			#if one value per line
			$data[$j]=$get[0];
			$j++;
		}
		else{
			#if several value on a line
			foreach $i (0..@get-1){
				$data[$j]=$get[$i];
				$j++;
			};
		};
	};
	close(IN);

####	#POINTS and ORIGIN from reference
	$point="";
	$origin="";
	$flagname=0;
	$flagname2=0;
	$flagnew=1;
	open(IN,"<$ref");
	while(<IN>){
		if($flagnew){
			$ig=1;
			$jg=0;
                        @strx='';
                        @stry='';
                        @strz='';
			@atomxlinked="";
                        @atom='';
			$compt=0;
			$flagnew=0;
			$istratom=0;
			$istrbond=0;
		};
		@getstr = split(' ',$_);
		$compt++;

		if (($compt > 4) && ($ig <= $istratom)){
			$strx[$ig]=sprintf "%4.2f",$getstr[0];
                        $stry[$ig]=sprintf "%4.2f",$getstr[1];
                        $strz[$ig]=sprintf "%4.2f",$getstr[2];
			$atom[$ig]=$getstr[3];
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
				if($atom[$getstr[0]] eq 'X'){
					$atomxlinked[$getstr[0]]=$getstr[1];
				}
				elsif($atom[$getstr[1]] eq 'X'){
					$atomxlinked[$getstr[1]]=$getstr[0];
				};
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
		};
		if($flagname && $getstr[0] ne ''){
			$point=$getstr[0];
			$flagname=0;
		};
		if($flagname2 && $getstr[0] ne ''){
			$origin=$getstr[0];
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
	#print "$point and $origin\n";
	@points="";
	@origins="";
	@get=split('-',$point);
	@get1=split('-',$origin);
	foreach $o (0..@get-1){
		$points[$get[$o]]=$get[$o];
		$origins[$get[$o]]=$get1[$o];
	};

####    #ORDERING
	$l=1; 
	$lmax=100; #extract max best solutions

	unlink "ordered-filtered-catsensaas.sdf" if(-e "ordered-filtered-catsensaas.sdf");
	unlink "ordered-filtered-matrix-sensaas.txt" if(-e "ordered-filtered-matrix-sensaas.txt");
	open(OUT,">ordered-filtered-matrix-sensaas.txt");
	foreach $i (0..@data-1){
		$min=-1;
		$mini=-1;
		foreach $j (0..@data-1){
			if($data[$j] ne ""){
				if($data[$j] > $min){
					$min=$data[$j];
					$mini=$j;
				};
			};
		};
		$k=$mini+1;

		&searchsdfi($sdf,$k,"atmp.sdf");

####	        #POINTS from sol	
		@lignesdf="";
		$pointsol="";
        	$originsol="";
        	$flagname=0;
        	$flagname2=0;
        	$flagnew=1;
        	open(IN,"<atmp.sdf");
		while(<IN>){
	                if($flagnew){
        	                $ig=1;
                	        $jg=0;
                        	@strxsol='';
                        	@strysol='';
                        	@strzsol='';
                        	@atomxlinkedsol="";
                        	@atomsol='';
                        	$compt=0;
                        	$flagnew=0;
                        	$istratom=0;
                        	$istrbond=0;
                	};
                	@getstr = split(' ',$_);
                	$compt++;
			$lignesdf[$compt]=$_;

			if (($compt > 4) && ($ig <= $istratom)){
				$strxsol[$ig]=sprintf "%4.2f",$getstr[0];
        	                $strysol[$ig]=sprintf "%4.2f",$getstr[1];
                	        $strzsol[$ig]=sprintf "%4.2f",$getstr[2];
                        	$atomsol[$ig]=$getstr[3];
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
                                if($atomsol[$getstr[0]] eq 'X'){
                                        $atomxlinkedsol[$getstr[0]]=$getstr[1];
                                }
                                elsif($atomsol[$getstr[1]] eq 'X'){
                                        $atomxlinkedsol[$getstr[1]]=$getstr[0];
                                };
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
        	        };
	                if($flagname && $getstr[0] ne ''){
				#to ease the replacement
				$lignesdf[$compt]="pointstoreplace\n";
        	        	$pointsol=$getstr[0];
 	                	$flagname=0;
                	};
                	if($flagname2 && $getstr[0] ne ''){
				#to ease the replacement
				$lignesdf[$compt]="origintoreplace\n";
                        	$originsol=$getstr[0];
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
		@pointsols="";
		@originsols="";
		@get=split('-',$pointsol);
		@get1=split('-',$originsol);
		foreach $o (0..@get-1){
			$pointsols[$get[$o]]=$get[$o];
			$originsols[$get[$o]]=$get1[$o];
		};

		#check if X dummy atoms are close
		$nbx=0;
		
		$lignepoint="";
		$ligneorigin="";
		$already="";	
		$alreadyo="";
		foreach $bi (1..@strx-1){
			if($atom[$bi] eq 'X'){
				$smallestdist=1000;
				$smallestdisti=0;
				foreach $bis (1..@strxsol-1){
					if($atomsol[$bis] eq 'X'){
						$dist=sqrt(($strxsol[$bis]-$strx[$bi])*($strxsol[$bis]-$strx[$bi]) + ($strysol[$bis]-$stry[$bi])*($strysol[$bis]-$stry[$bi]) + ($strzsol[$bis]-$strz[$bi])*($strzsol[$bis]-$strz[$bi]));
						if($dist<$smallestdist && $already!~/ $bis / && $alreadyo!~/ $bi /){
							#print "$bi $bis $dist\n";
							$smallestdist=$dist;
							$smallestdisti=$bis;
						};
					};
				};
				if($smallestdist<1000 && $smallestdisti>0){
					if($smallestdist < $distmax && $already!~/ $smallestdisti / && $alreadyo!~/ $bi /){
						#print "$bi $smallestdisti $smallestdist\n";
						$nbx++;
						if($lignepoint eq ""){
							$ligneorigin=$origins[$atomxlinked[$bi]];
							$lignepoint=$atomxlinkedsol[$smallestdisti];
						}
						else{
							$ligneorigin=$ligneorigin."-".$origins[$atomxlinked[$bi]];
							$lignepoint=$lignepoint."-".$atomxlinkedsol[$smallestdisti];
						};
						$already=$already." $smallestdisti ";
						$alreadyo=$alreadyo." $bi ";
						#last;
					};
				};
			};
		};

		if($nbx >= $nbmatchx){ 
		#if($nbx >= $nbmatchx && $l <= $lmax){
			print OUT "$k $data[$mini]\n";

			#print "sol $k:\n$lignepoint\n$ligneorigin\n";
			#if POINTS and ORIGIN exist in atmp.sdf
			open(OUP,">atmpnew.sdf");
			$vupoint=0;
			$vuorigin=0;
			foreach $bi (0..@lignesdf-1){
				$vupoint=1 if($lignesdf[$bi]=~/POINTS/);
				$vuorigin=1 if($lignesdf[$bi]=~/ORIGIN/);
				if($lignesdf[$bi]=~/pointstoreplace/){
					$lignesdf[$bi]=~s/pointstoreplace/$lignepoint/;
				};
				if($lignesdf[$bi]=~/origintoreplace/){
					$lignesdf[$bi]=~s/origintoreplace/$ligneorigin/;
				};
				if($lignesdf[$bi]=~/\$\$\$\$/ && $vupoint==0){
					print OUP "> <POINTS>\n";
					print OUP "$lignepoint\n";
					print OUP "\n";
				};
				if($lignesdf[$bi]=~/\$\$\$\$/ && $vuorigin==0){
					print OUP "> <ORIGIN>\n";
					print OUP "$ligneorigin\n";
					print OUP "\n";
				};
				print OUP "$lignesdf[$bi]";
			};
			close(OUP);
	
			open(IN,">>ordered-filtered-catsensaas.sdf");
			open(OP,"<atmpnew.sdf");
                        while(<OP>){
                            print IN $_;
                        };
                	close(OP);
                	close(IN);
			
			$l++;
		};

		$data[$mini]="";
	};
	unlink "atmp.sdf" if(-e "atmp.sdf");
	unlink "atmpnew.sdf" if(-e "atmpnew.sdf");
	close(OUT);
	$l=$l-1;
	print "see ordered-filtered-matrix-sensaas.txt and ordered-filtered-catsensaas.sdf with $l solutions\n";

###########################################################################################
#SUBROUTINES

sub nbax{
       local($fsdf)=@_;

	$nbxatom=0;
	$flagnew=1;
        open(MOL,"<$fsdf");
        while(<MOL>){
		if($flagnew){
			$compt=0;
                        $ig=1;
                        $jg=0;
			$flagnew=0;
                };
                @getstr = split(' ',$_);
                $compt++;
                if (($compt > 4) && ($ig <= $istratom)){
			$nbxatom++ if($getstr[3] eq 'X');
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
                        print "$nbxatom\n";
                };
        };
        close(MOL);
	#print "$nbxatom\n";
	$nbxatom;
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
        open(INI,"<$fsdf");
        open(OU,">$fsdfo");
        while (<INI>){
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
                                        print OU $ligne[$i];
                                };
                                $flagname=0;
                                last;
                        };
                };
        };
        close(INI);
        close(OU);
};

###################################################################################
###################################################################################


