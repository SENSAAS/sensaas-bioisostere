#!/usr/bin/perl

#list of libraries in db
#drug, drug2, drug3, pdb, pdb2, pdb3 and test

        my($filesdf,$filefgts,$nofgt,$db)=@ARGV;

        if($filesdf eq '' || $filefgts eq '' || $nofgt eq '' || $db eq ''){
                print "usage: \n sensaas-bioisostere.pl <file sdf of the entire molecule> <file sdf containing fragments> <no of the fgt to replace in filefgts> <library of fragment to screen >\n";
                exit;
        };

	if(!-e "$db" || -z "$db"){
		print "Library of fragments $db not found\n";
		exit;
	};

	$database=$db;
	$db="lib";
	print "database $db address $database \n";
	
	$leaexe=$0;
	#windows requires 2 steps
	$leaexe=~s/sensaas-bioisostere\.pl$//;
	$leaexe=~s/\/$//;
	print "perl scripts in $leaexe\n";
	
	print "SENSAAS-Bioisostere:\nTarget $filesdf \n";

	$nb2=&nbsdf($filefgts);

	#clean
	foreach $i (1..$nb2){
		$fgtna="libfgt_".$i.".sdf";
		unlink "$fgtna" if(-e "$fgtna");
	};

	chop($log=`python $leaexe/splitsdf.py $filefgts libfgt `);
	$referencex="libfgt_".$nofgt.".sdf";
	open(IN,">referencex.sdf");
        open(OP,"<$referencex");
                while(<OP>){
                        print IN $_;
                };
        close(OP);
	close(IN);


	#how many subtitutions
	chop($nbx=`python $leaexe/nbx.py $referencex `);
	print "$nb2 fragments in $filefgts - referencex.sdf = $referencex contains the fragment to replace with $nbx substituants\n";

	#reduce computing time
	unlink "filterx.sdf" if (-e "filterx.sdf");
	chop($log=`python $leaexe/filterx.py $database $nbx `);
	#output filterx.sdf
	$nb5=&nbsdf("filterx.sdf");
	print "Library: $nb5 fragments in $database having >= $nbx substitution positions\n";
	
	#execute similarity calculations
	print "Screening using metasensaas (option -s mean -l 2) \n";
	chop($log=`python $leaexe/meta-sensaas.py referencex.sdf filterx.sdf –s mean -l 2 `);
	open(OYO,">meta-sensaas_fgt.out");
	print OYO "$log\n";
	close(OYO);
	
	#rank and check substitution are appropriately positioned 
	chop($log=`perl $leaexe/lea3d-ordered-catsensaas-filtered.pl matrix-sensaas.txt catsensaas.sdf referencex.sdf $nbx `);
	$m1="fragments-".$db.".sdf";
        $m2="fragments-".$db."-score.txt";
        $m3="ordered-filtered-".$m1;
        $m4="ordered-filtered-".$m2;
	rename "matrix-sensaas.txt", "$m2";
        rename "catsensaas.sdf", "$m1";
        rename "ordered-filtered-catsensaas.sdf", "$m3";
        rename "ordered-filtered-matrix-sensaas.txt", "$m4";

	#Build the entire bioisoster molecule
	chop($log=`python $leaexe/build-bioisostere.py $m3 libfgt_ $referencex `);
	# X dummy atoms must be replaced by H
	chop($log=`python $leaexe/replacex.py build-bioisostere.sdf `);
	rename "replacex.sdf", "build-bioisostereH.sdf";

	
	#Minimize 3D structures
	chop($log=`perl $leaexe/lea3d-conformers.pl build-bioisostereH.sdf rdkitmin 1 `);
	rename "build-bioisostereH3D.sdf", "build-bioisostereH-min.sdf";
	$nb3=&nbsdf("build-bioisostereH-min.sdf");
	print "$nb3 built bioisosteres\n";

	#Realign solutions on the Target
	chop($log=`python $leaexe/meta-sensaas.py $filesdf build-bioisostereH-min.sdf -l 2 `);
	open(OYO,">meta-sensaas_bioisostere.out");
        print OYO "$log\n";
        close(OYO);

	#rank
	chop($log=`python $leaexe/ordered-catsensaas.py matrix-sensaas.txt catsensaas.sdf `);
	
	$m5="bioisosteres-".$db.".sdf";
        $m6="bioisosteres-".$db."-score.txt";
        $m7="ordered-".$m5;
        $m8="ordered-".$m6;
        rename "matrix-sensaas.txt", "$m6";
        rename "catsensaas.sdf", "$m5";
        rename "ordered-catsensaas.sdf", "$m7";
        rename "ordered-scores.txt", "$m8";

	#remove datablocks
	chop($log=`perl $leaexe/stripsdf.pl $m7 stripsdf.sdf \"-\" `);

        #add scores in datablock
        chop($log=`perl $leaexe/adddatablocksensaas.pl stripsdf.sdf SENSAAS-SCORE $m8 `);
        rename "adddatablock.sdf", "$m7";
        unlink "stripsdf.sdf";

	#clean
        foreach $i (1..$nb2){
                $fgtna="libfgt_".$i.".sdf";
                unlink "$fgtna" if(-e "$fgtna");
        };
	unlink "referencex.sdf" if(-e "referencex.sdf");
	unlink "filterx.sdf" if(-e "filterx.sdf");
	unlink "combination.sdf" if(-e "combination.sdf");
	unlink "bestsensaas.sdf" if(-e "bestsensaas.sdf");
	unlink "copyfgt.sdf" if(-e "copyfgt.sdf");

	$nb4=&nbsdf($m7);
	 print "Visualize the $nb4 aligned bioisosteres in $m7 with PyMol (scores in $m8) - Selected fragments from $db are in $m3 (scores in $m4) \n";

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

