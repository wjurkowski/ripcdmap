#!/usr/bin/perl -w
use strict;
use warnings;

if ($#ARGV < 0 ) {die "Run type available \
-r reactions\
-s species\
-m modules\n";
}


if (($ARGV[0] eq "-r") and ($#ARGV != 2)) {die "Program used with parameters [-r] [CellDesigner xml] [choose reaction type: -ge,-reg, -met, -all]\n";}
elsif (($ARGV[0] eq "-s") and ($#ARGV != 1)) {die "Program used with parameters [-s] [CellDesigner xml]\n";}
elsif (($ARGV[0] eq "-m") and ($#ARGV != 1)) {die "Program used with parameters [-m] [CellDesigner xml]\n";}
elsif ($ARGV[0] ne "-s" and $ARGV[0] ne "-r" and $ARGV[0] ne "-m") {die "Wrong run type. Availabe options are:\
-r reactions\
-s species\
-m modules\n";}

my @map=open_file($ARGV[1]);
my $runtype=$ARGV[0];
my (%SAiR, %SiR);

if($runtype eq "-r" or $runtype eq "-m") {
my (@reactions,@modifications,@modifiers,@substrates,@products,%species_names,%species_localization,@np,@ns,@nm,@pubmed);
my $rn=0;
my $rtswitch=$ARGV[2];
my $outt=$ARGV[3];
my $lOReact=0;
my $lOSpec=0;

for my $k (0..$#map){#parse the map
	my $lin=$map[$k];
#species description from listOfSpecies
	if($lin =~ /.*<listOfSpecies>/){$lOSpec = 1;}
        if($lin =~ /.*<\/listOfSpecies>/){$lOSpec = 0;}
	if($lOSpec == 1 and $lin =~ /.+species metaid.+id.{2}(.+).{2}name.{2}(.+).{2}compartment.{2}(\w+)".+$/){
		$species_names{$1}=$2;
		$species_localization{$1}=$3;
	}
#reactions
	if($lin =~ /.*<listOfReactions>/){$lOReact = 1;}
        if($lin =~ /.*<\/listOfReactions>/){$lOReact = 0;}
	if($lOReact == 1 and $lin =~ /.*<reaction.+ id="(\w+)".+reversible="(\w+)".$/){
		$rn++;
		my $revers=$2;
		$reactions[$rn][0] = $1;#rid
		$reactions[$rn][1] = $revers;
		$ns[$rn]=0;
		$np[$rn]=0;
		$nm[$rn]=0;
		$pubmed[$rn]="Reference unknown";
	}
	elsif($lOReact == 1 and $lin =~ /.*<reaction.+ id="(\w+)".+$/){
                $rn++;
                my $revers="true";
                $reactions[$rn][0] = $1;#rid
                $reactions[$rn][1] = $revers;
                $ns[$rn]=0;
                $np[$rn]=0;
                $nm[$rn]=0;
                $pubmed[$rn]="Reference unknown";
        }
	if($lOReact == 1 and $lin =~ /.*<celldesigner.reactionType.{1}(\w+).{2}celldesigner.+$/) {
		my $rtype = $1;
		$reactions[$rn][2] = $rtype;
	}
	if($lOReact == 1 and $lin =~ /.*<celldesigner.modification type.{2}(\w+).{2}modifiers.{2}(\w+).{2}aliases.{2}(\w+)".+$/) {
		$nm[$rn]++;
		$modifications[$rn][$nm[$rn]] = $1;#modification
		$modifiers[$rn][$nm[$rn]] = $2;#modifier
		my $salias = $3;
		push(@{ $SAiR{$salias} }, $reactions[$rn][0]);
		$SiR{$2}++;
	}
	if($lOReact == 1 and $lin =~ /.*<celldesigner:baseReactant alias.{2}(\w+).{2}species.{2}(\w+)".+$/) {
		$ns[$rn]++;
		my $substrate = $2;
		$substrates[$rn][$ns[$rn]] = $substrate;
		my $salias = $1;
		push(@{ $SAiR{$salias} }, $reactions[$rn][0]);
		#print "$salias: @{ $SiR{$salias} }\n";
		$SiR{$substrate}++;
	}  
	if($lOReact == 1 and $lin =~ /.*<celldesigner:reactantLink alias.{2}(\w+).{2}reactant.{2}(\w+)".+$/) {
		$ns[$rn]++;
                my $substrate = $2;
                $substrates[$rn][$ns[$rn]] = $substrate;
		my $salias = $1;
		push(@{ $SAiR{$salias} }, $reactions[$rn][0]);
		$SiR{$substrate}++;
	}
	if($lOReact == 1 and $lin =~ /.*<celldesigner:baseProduct alias.{2}(\w+).{2}species.{2}(\w+)".+$/) {
		$np[$rn]++;
		my $product = $2;	
		$products[$rn][$np[$rn]] = $product;
		my $salias = $1;
		push(@{ $SAiR{$salias} }, $reactions[$rn][0]);
		$SiR{$product}++;
	}
	if($lOReact == 1 and $lin =~ /.*<celldesigner:productLink alias.{2}(\w+).{2}product.{2}(\w+)".+$/) {
		$np[$rn]++;
		my $product = $2;	
		$products[$rn][$np[$rn]] = $product;
		my $salias = $1;
		push(@{ $SAiR{$salias} }, $reactions[$rn][0]);
		$SiR{$product}++;
	}
	if($lOReact == 1 and $lin =~ /.+rdf.resource.+pubmed.{1}(\w+).{3}/){
		$pubmed[$rn]=$1;
	}
}#end of map parsing
#merge multiple species states? 

if($runtype eq "-r"){#analyse reactions
#do not print in -m mode
my $reakcje=substr($ARGV[1],0,rindex($ARGV[1],".")).".reactions.txt";
open(OUT, ">$reakcje");

for my $i (1..$rn){#iterate reactions
  if ($rtswitch eq "-ge"){#regulation of gene expression
#print "$reactions[$i][2] \n";
   if($reactions[$i][2] eq "TRANSCRIPTION" or $reactions[$i][2] eq "TRANSLATION") {
	for my $m (1..$np[$i]){
	  if($ns[$i] > 1 or $nm[$i] > 0){
		for my $l (1..$ns[$i]){
		    printf OUT "$reactions[$i][0]\t$reactions[$i][1]\t$reactions[$i][2]\t$substrates[$i][$l] -> $products[$i][$m].ps\t$pubmed[$i]\t$species_names{$substrates[$i][$l]}-$substrates[$i][$l] -> $species_names{$products[$i][$m]}-$products[$i][$m].ps\n";
		    #$graph{$reactions[$i][0]}="$substrates[$i][$l] -> $products[$i][$m].ps";
		}
		for my $l (1..$nm[$i]){
 		 my $tab0=$modifications[$i][$l];
 		 my $tab1=$modifiers[$i][$l];
		 if($tab0 eq "CATALYSIS" or $tab0 eq "TRIGGER"){
		    printf OUT "$reactions[$i][0]\t$reactions[$i][1]\t$reactions[$i][2]\t$tab1 -> $products[$i][$m].ps\t$pubmed[$i]\t$species_names{$tab1}-$tab1 -> $species_names{$products[$i][$m]}-$products[$i][$m].ps\n";
		    #$graph{$reactions[$i][0]}="$tab1 -> $products[$i][$m].ps";
		 }
		 elsif($tab0 eq "INHIBITION"){
		    printf OUT "$reactions[$i][0]\t$reactions[$i][1]\t$reactions[$i][2]\t$tab1 -| $products[$i][$m].ps\t$pubmed[$i]\t$species_names{$tab1}-$tab1 -| $species_names{$products[$i][$m]}-$products[$i][$m].ps\n";
		    #$graph{$reactions[$i][0]}="$tab1 -| $products[$i][$m].ps";
		 }
		}
		printf OUT "$reactions[$i][0]\t$reactions[$i][1]\t$reactions[$i][2]\t$products[$i][$m].ps -> $products[$i][$m]\t$pubmed[$i]\t$species_names{$products[$i][$m]}-$products[$i][$m].ps -> $species_names{$products[$i][$m]}-$products[$i][$m]\n";
		#$graph{$reactions[$i][0]}="$products[$i][$m].ps -> $products[$i][$m]";
	  }
	  else{
                printf OUT "$reactions[$i][0]\t$reactions[$i][1]\t$reactions[$i][2]\t$substrates[$i][1] -> $products[$i][$m]\t$pubmed[$i]\t$species_names{$substrates[$i][1]}-$substrates[$i][1] -> $species_names{$products[$i][$m]}-$products[$i][$m]\n";
		#$graph{$reactions[$i][0]}="$substrates[$i][1] -> $products[$i][$m]";
          }
	}
   }
  }

  if ($rtswitch eq "-reg" or $rtswitch eq "-all"){#all types of regulatory interactions
   if($reactions[$i][2] eq "TRANSCRIPTION" or $reactions[$i][2] eq "TRANSLATION"  or $reactions[$i][2] eq "POSITIVE_INFLUENCE" or $reactions[$i][2] eq "UNKNOWN_POSITIVE_INFLUENCE" or $reactions[$i][2] eq "TRIGGER" or $reactions[$i][2] eq "UNKNOWN_REDUCED_TRIGGER") {
	for my $m (1..$np[$i]){
	  if($ns[$i] > 1 or $nm[$i] > 0){
		for my $l (1..$ns[$i]){
		   printf OUT "$reactions[$i][0]\t$reactions[$i][1]\t$reactions[$i][2]\t$substrates[$i][$l] -> $products[$i][$m].ps\t$pubmed[$i]\t$species_names{$substrates[$i][$l]}-$substrates[$i][$l] -> $species_names{$products[$i][$m]}-$substrates[$i][$l].ps\n";
		   #$graph{$reactions[$i][0]}="$substrates[$i][$l] -> $products[$i][$m].ps";
		}
          	for my $l (1..$nm[$i]){
		 my $tab0=$modifications[$i][$l];
 		 my $tab1=$modifiers[$i][$l];
                 if($tab0 eq "CATALYSIS" or $tab0 eq "TRIGGER" or $tab0 eq "UNKNOWN_CATALYSIS"){
		   printf OUT "$reactions[$i][0]\t$reactions[$i][1]\t$reactions[$i][2]\t$tab1 -> $products[$i][$m].ps\t$pubmed[$i]\t$species_names{$tab1}-$tab1 -> $species_names{$products[$i][$m]}-$products[$i][$m].ps\n";
		   #$graph{$reactions[$i][0]}="$tab1 -> $products[$i][$m].ps";
                 }
                 elsif($tab0 eq "INHIBITION" or $tab0 eq "UNKNOWN_INHIBITION"){
		    printf OUT "$reactions[$i][0]\t$reactions[$i][1]\t$reactions[$i][2]\t$tab1 -| $products[$i][$m].ps\t$pubmed[$i]\t$species_names{$tab1}-$tab1 -| $species_names{$products[$i][$m]}-$products[$i][$m].ps\n";
		    #$graph{$reactions[$i][0]}="$tab1 -| $products[$i][$m].ps";
		 }
          	}
		printf OUT "$reactions[$i][0]\t$reactions[$i][1]\t$reactions[$i][2]\t$products[$i][$m].ps -> $products[$i][$m]\t$pubmed[$i]\t$species_names{$products[$i][$m]}-$products[$i][$m].ps -> $species_names{$products[$i][$m]}-$products[$i][$m]\n";
		#$graph{$reactions[$i][0]}="$products[$i][$m].ps -> $products[$i][$m]";
          }
	  else{
		printf OUT "$reactions[$i][0]\t$reactions[$i][1]\t$reactions[$i][2]\t$substrates[$i][1] -> $products[$i][$m]\t$pubmed[$i]\t$species_names{$substrates[$i][1]}-$substrates[$i][1] -> $species_names{$products[$i][$m]}-$products[$i][$m]\n";
		#$graph{$reactions[$i][0]}="$substrates[$i][1] -> $products[$i][$m]";
	  }	
	}
   }
   if($reactions[$i][2] eq "NEGATIVE_INFLUENCE" or $reactions[$i][2] eq "UNKNOWN_NEGATIVE_INFLUENCE") {
   	for my $m (1..$np[$i]){
	  if($ns[$i]>1 or $nm[$i]>0){
		for my $l (1..$ns[$i]){
		   printf OUT "$reactions[$i][0]\t$reactions[$i][1]\t$reactions[$i][2]\t$substrates[$i][$l] -| $products[$i][$m].ps\t$pubmed[$i]\t$species_names{$substrates[$i][$l]}-$substrates[$i][$l] -| $species_names{$products[$i][$m]}-$products[$i][$m].ps\n";
		   #$graph{$reactions[$i][0]}="$substrates[$i][$l] -| $products[$i][$m].ps";
                }
                for my $l (1..$nm[$i]){
		  my $tab0=$modifications[$i][$l];
 		  my $tab1=$modifiers[$i][$l];
		  if($tab0 eq "CATALYSIS" or $tab0 eq "TRIGGER" or $tab0 eq "UNKNOWN_CATALYSIS"){
		    printf OUT "$reactions[$i][0]\t$reactions[$i][1]\t$reactions[$i][2]\t$tab1 -| $products[$i][$m].ps\t$pubmed[$i]\t$species_names{$tab1}-$tab1 -| $species_names{$products[$i][$m]}-$products[$i][$m].ps\n";
		    #$graph{$reactions[$i][0]}="$tab1 -| $products[$i][$m].ps";
		  }
		  elsif($tab0 eq "INHIBITION" or $tab0 eq "UNKNOWN_INHIBITION"){
		    printf OUT "$reactions[$i][0]\t$reactions[$i][1]\t$reactions[$i][2]\t$tab1 -> $products[$i][$m].ps\t$pubmed[$i]\t$species_names{$tab1}-$tab1 -> $species_names{$products[$i][$m]}-$products[$i][$m].ps\n";
		    #$graph{$reactions[$i][0]}="$tab1 -> $products[$i][$m].ps";
		  }
                }
		printf OUT "$reactions[$i][0]\t$reactions[$i][1]\t$reactions[$i][2]\t$products[$i][$m].ps -> $products[$i][$m]\t$pubmed[$i]\t$species_names{$products[$i][$m]}-$products[$i][$m].ps -> $species_names{$products[$i][$m]}-$products[$i][$m]\n";
		#$graph{$reactions[$i][0]}="$products[$i][$m].ps -> $products[$i][$m]";
	  }
	  else{
                printf OUT "$reactions[$i][0]\t$reactions[$i][1]\t$reactions[$i][2]\t$substrates[$i][1] -| $products[$i][$m]\t$pubmed[$i]\t$species_names{$substrates[$i][1]}-$substrates[$i][1] -| $species_names{$products[$i][$m]}-$products[$i][$m]\n";
		#$graph{$reactions[$i][0]}="$substrates[$i][1] -| $products[$i][$m]";
          }
	}
   }
  }
  if ($rtswitch eq "-met" or $rtswitch eq "-all"){#metabolic network and other physical transitions
   if($reactions[$i][2] eq "STATE_TRANSITION" or $reactions[$i][2] eq "KNOWN_TRANSITION_OMITTED" or $reactions[$i][2] eq "TRANSPORT" or $reactions[$i][2] eq "DISSOCIATION" or $reactions[$i][2] eq "TRUNCATION" or $reactions[$i][2] eq "HETERODIMER_ASSOCIATION" ) {
	for my $m (1..$np[$i]){
print "dupa\n";
	  my $prod=$products[$i][$m];
	  if (not exists $species_names{$prod}){$species_names{$prod}="no name";}
	  if($ns[$i]>1 or $nm[$i]>0){
        	for my $l (1..$ns[$i]){
		   my $subs=$substrates[$i][$l]; 
		   if (not exists $species_names{$subs}){$species_names{$subs}="no name";}
		   printf OUT "$reactions[$i][0]\t$reactions[$i][1]\t$reactions[$i][2]\t$subs -> $prod.ps\t$pubmed[$i]\t$species_names{$subs}-$subs -> $species_names{$prod}-$prod.ps\n";
		   #$graph{$reactions[$i][0]}="$substrates[$i][$l] -> $products[$i][$m].ps"; 
                }
                for my $l (1..$nm[$i]){
		  my $tab0=$modifications[$i][$l];
 		  my $tab1=$modifiers[$i][$l];
		  if (not exists $species_names{$tab1}){$species_names{$tab1}="no name";}
		  if($tab0 eq "CATALYSIS" or $tab0 eq "TRIGGER" or $tab0 eq "UNKNOWN_CATALYSIS"){
		    my $subs=$tab1;
		    printf OUT "$reactions[$i][0]\t$reactions[$i][1]\t$reactions[$i][2]\t$subs -> $prod.ps\t$pubmed[$i]\t$species_names{$subs}-$subs -> $species_names{$prod}-$prod.ps\n";
		    #$graph{$reactions[$i][0]}="$tab1 -> $products[$i][$m].ps"; 
		  }
		  elsif($tab0 eq "INHIBITION" or $tab0 eq "UNKNOWN_INHIBITION"){
		    my $subs=$tab1;
		    printf OUT "$reactions[$i][0]\t$reactions[$i][1]\t$reactions[$i][2]\t$subs -> $prod.ps\t$pubmed[$i]\t$species_names{$subs}-$subs -> $species_names{$prod}-$prod.ps\n";
		    #$graph{$reactions[$i][0]}="$tab1 -| $products[$i][$m].ps"; 
		  }
                }
		printf OUT "$reactions[$i][0]\t$reactions[$i][1]\t$reactions[$i][2]\t$prod.ps -> $prod\t$pubmed[$i]\t$species_names{$prod}-$prod.ps -> $species_names{$prod}-$prod\n";
		#$graph{$reactions[$i][0]}="$products[$i][$m].ps -> $products[$i][$m]";
	  }
	  else{
		my $subs=$substrates[$i][1];
		if (not exists $species_names{$substrates[$i][1]}){$species_names{$substrates[$i][1]}="no name";}
        	printf OUT "$reactions[$i][0]\t$reactions[$i][1]\t$reactions[$i][2]\t$subs -> $prod\t$pubmed[$i]\t$species_names{$subs}-$subs -> $species_names{$prod}-$prod\n";
		#$graph{$reactions[$i][0]}="$substrates[$i][1] -> $products[$i][$m]";
          }
	}
	if($reactions[$i][1] eq "true"){#resolve reversibility of reactions
	for my $l (1..$ns[$i]){
	  my $prod=$substrates[$i][$l];
	  if (not exists $species_names{$prod}){$species_names{$prod}="no name";}
	  if($np[$i]>1 or $nm[$i]>0){
		for my $m (1..$np[$i]){
		   my $subs=$products[$i][$m];
		   if (not exists $species_names{$subs}){$species_names{$subs}="no name";}
		   printf OUT "$reactions[$i][0]\t$reactions[$i][1]\t$reactions[$i][2]\t$subs -> $prod.ps\t$pubmed[$i]\t$species_names{$subs}-$subs -> $species_names{$prod}-$prod.ps\n";
		   #$graph{$reactions[$i][0]}="$substrates[$i][$l] -> $products[$i][$m].ps"; 
                }
                for my $l (1..$nm[$i]){
		  my $tab0=$modifications[$i][$l];
 		  my $tab1=$modifiers[$i][$l];
		  if (not exists $species_names{$tab1}){$species_names{$tab1}="no name";}
		  if($tab0 eq "CATALYSIS" or $tab0 eq "TRIGGER" or $tab0 eq "UNKNOWN_CATALYSIS"){
		    my $subs=$tab1;
		    printf OUT "$reactions[$i][0]\t$reactions[$i][1]\t$reactions[$i][2]\t$subs -> $prod.ps\t$pubmed[$i]\t$species_names{$subs}-$subs -> $species_names{$prod}-$prod.ps\n";
		    #$graph{$reactions[$i][0]}="$tab1 -> $products[$i][$m].ps"; 
		  }
		  elsif($tab0 eq "INHIBITION" or $tab0 eq "UNKNOWN_INHIBITION"){
		    my $subs=$tab1;
		    printf OUT "$reactions[$i][0]\t$reactions[$i][1]\t$reactions[$i][2]\t$subs -> $prod.ps\t$pubmed[$i]\t$species_names{$subs}-$subs -> $species_names{$prod}-$prod.ps\n";
		    #$graph{$reactions[$i][0]}="$tab1 -| $products[$i][$m].ps"; 
		  }
                }
		printf OUT "$reactions[$i][0]\t$reactions[$i][1]\t$reactions[$i][2]\t$prod.ps -> $prod\t$pubmed[$i]\t$species_names{$prod}-$prod.ps -> $species_names{$prod}-$prod\n";
		#$graph{$reactions[$i][0]}="$products[$i][$m].ps -> $products[$i][$m]";
	  }
	  else{
		my $subs=$products[$i][1];
		if (not exists $species_names{$substrates[$i][1]}){$species_names{$substrates[$i][1]}="no name";}
        	printf OUT "$reactions[$i][0]\t$reactions[$i][1]\t$reactions[$i][2]\t$subs -> $prod\t$pubmed[$i]\t$species_names{$subs}-$subs -> $species_names{$prod}-$prod\n";
		#$graph{$reactions[$i][0]}="$substrates[$i][1] -> $products[$i][$m]";
          }
	}  
	}
   }
  }
}#end of reaction iteration
}#end of printing switch
}#end of reaction parsing block
if($runtype eq "-s" or $runtype eq "-m") {#parse species 
 my (@species,%species_name,%species_localization,%species_position,%species_type,@fcompx1,@fcompx2,@fcompy1,@fcompy2,@fcomp_name);
 my (@incspecies,%incspecies_name,%incspecies_withinof,$spec,%compx,%compy,%spx,%spy,%spinfm,$specal,%specalias);
 my (%compartment_name,%compartment_outside);
 my $lOInclSpec=0;
 my $lOCSpecAl=0;
 my $lOSpecAl=0;
 my $lOLayers=0;
 my $lOComp=0;
 my $lOSpec=0;
 my $ns=0;
 my $nl=0;
 my $nis=0;

  for my $k (0..$#map){
        	my $lin=$map[$k];
#from celldesigner:listOfIncludedSpecies
		if($lin =~ /.*<celldesigner:listOfIncludedSpecies>/){$lOInclSpec = 1;}
                if($lin =~ /.*<\/celldesigner:listOfIncludedSpecies>/){$lOInclSpec = 0;}
 	#species id
		if($lOInclSpec == 1 and $lin =~ /.*<celldesigner:species.+id.{2}(.+).{2}name.{2}(.+)">/){
			$nis++;
			$incspecies[$nis]=$1;
			$incspecies_name{$1}=$2;
		}
	#complexSpecies
		if($lOInclSpec == 1 and $lin =~ /.*<celldesigner:complexSpecies>(.+)<.celldesigner:complexSpecies>/){
			$incspecies_withinof{$incspecies[$nis]}=$1;#hash of uniq species id being part of a complex $1
		}
	#proteinReference
	#listOfModifications
		#To be done - VERY IMPORTANT!
#celldesigner:listOfCompartmentAliases
	#id, coordinates and color
#celldesigner:listOfComplexSpeciesAliases (species that are complexes)
		if($lin =~ /.*<celldesigner:listOfComplexSpeciesAliases>/){$lOCSpecAl = 1;}
                if($lin =~ /.*<\/celldesigner:listOfComplexSpeciesAliases>/){$lOCSpecAl = 0;}
	#complexSpeciesAlias - mapping to species id
                if($lOCSpecAl == 1 and $lin =~ /.*<celldesigner:complexSpeciesAlias.+id.{2}(\w+).+species.{2}(\w+)".+/){
			$spec=$1;
			$specalias{$1}=$2;
		}
	#coordinates
		if($lOCSpecAl == 1 and $lin =~ /.*<celldesigner:bounds x.{2}(.+).{2}y.{2}(.+).{2}w.{2}(.+).{2}h.{2}(.+).{3}$/){
			# x and y of upper-left corner of the box
			$spx{$spec} = ($1+$1+$3)/2;
                        $spy{$spec} = ($2+$2+$4)/2;
		}
#celldesigner:listOfSpeciesAliases> (species and includes species)
		if($lin =~ /.*<celldesigner:listOfSpeciesAliases>/){$lOSpecAl = 1;}
                if($lin =~ /.*<\/celldesigner:listOfSpeciesAliases>/){$lOSpecAl = 0;}
	#speciesAlias mapping to species id
		if($lOSpecAl == 1 and $lin =~ /.*<celldesigner:speciesAlias.+id.{2}(\w+).+species.{2}(\w+)".+/){
			$specal=$1;
			$specalias{$1}=$2;#list of species ID for given alias
		}
	#coordinates
		if($lOSpecAl == 1 and $lin =~ /.*<celldesigner:bounds x.{2}(.+).{2}y.{2}(.+).{2}w.{2}(.+).{2}h.{2}(.+).{3}$/){
			# x and y of upper-left corner of the box
			$spx{$specal} = ($1+$1+$3)/2;
			$spy{$specal} = ($2+$2+$4)/2;
		}
	#activity
#celldesigner:listOfGroups>
#celldesigner:listOfProteins
	#protein id/reference
#celldesigner:listOfGenes
#celldesigner:listOfRNAs
#celldesigner:listOfLayers
		if($lin =~ /.*<celldesigner:listOfLayers>/){$lOLayers = 1;}
		if($lin =~ /.*<\/celldesigner:listOfLayers>/){$lOLayers = 0;}
#<celldesigner:layerSpeciesAlias contains all info
		if($lOLayers == 1 and $lin =~ /.*<celldesigner:layerSpeciesAlias.+/){$nl++;}#iterate number of functional module
		if($lOLayers == 1 and $lin =~ /.*<celldesigner:bounds x.{2}(.+).{2}y.{2}(.+).{2}w.{2}(.+).{2}h.{2}(.+).{3}$/){
			# x and y of upper-left corner of the box
			$fcompx1[$nl] = $1;
			$fcompy1[$nl] = $2;
			$fcompx2[$nl] = $fcompx1[$nl]+$3;
			$fcompy2[$nl] = $fcompy1[$nl]+$4;
		}
		if($lOLayers == 1 and not $lin =~ /.+celldesigner.+/){
			$fcomp_name[$nl] = $lin;
		}
#celldesigner:listOfBlockDiagrams
#from listOfCompartments:
		if($lin =~ /.*<listOfCompartments>/){$lOComp = 1;}
                if($lin =~ /.*<\/listOfCompartments>/){$lOComp = 0;}
        	if($lOComp == 1){
		  $compartment_name{"default"}="extracellular";
		  if($lin =~ /.*compartment metaid.+id.{2}(.+).{2}name.{2}(.+).{2}size.+outside.{2}(\w+)".+$/){
			$compartment_name{$1}=$2;#name of compartment
			$compartment_outside{$1}=$3;
		  }
		}
#from listOfSpecies
		if($lin =~ /.*<listOfSpecies>/){$lOSpec = 1;}
		if($lin =~ /.*<\/listOfSpecies>/){$lOSpec = 0;}
        	if($lOSpec == 1 and $lin =~ /.+species metaid.+id.{2}(.+).{2}name.{2}(.+).{2}compartment.{2}(\w+)".+$/){
			$species_name{$1}=$2;
			$species_localization{$1}=$3;
			$ns++;
			$species[$ns]=$1;
			$spinfm{$1}="none";
		}
		if($lOSpec == 1 and $lin =~ /.+positionToCompartment.{1}(.+).{2}celldesigner.+$/){
			$species_position{$species[$ns]}=$1;
		}
		 if($lOSpec == 1 and $lin =~ /.+class.{1}(.+).{2}celldesigner.+$/){
                        $species_type{$species[$ns]}=$1;
                }
  }
	chomp @fcomp_name;# remove trailing characters 
  if($runtype eq "-s") {#print species
	my $species=substr($ARGV[1],0,rindex($ARGV[1],".")).".species.txt";
	my $complexes=substr($ARGV[1],0,rindex($ARGV[1],".")).".complexes.txt";
	open(OUT, ">$species");
	open(OUT2, ">$complexes");
	for my $i (1..$ns){
	  unless (exists $incspecies_withinof{$species[$i]}) {
		$incspecies_withinof{$species[$i]} = "not_in_complex";
	  }
	  if($species_localization{$species[$i]} ne "default"){
		printf OUT "$species[$i]\t$species_name{$species[$i]}\t$species_type{$species[$i]}\t$species_localization{$species[$i]}\t$compartment_name{$species_localization{$species[$i]}}\t$species_position{$species[$i]}\t$compartment_outside{$species_localization{$species[$i]}}\n";
	  }
	  elsif($species_localization{$species[$i]} eq "default"){
		printf OUT "$species[$i]\t$species_name{$species[$i]}\t$species_type{$species[$i]}\t$species_localization{$species[$i]}\tnone\t$species_position{$species[$i]}\tnone\n";
	  }
	}

	foreach my $key (keys %incspecies_withinof){
	  unless($incspecies_withinof{$key} eq "not_in_complex"){
		my $name="noname";
		if(exists $species_name{$incspecies_withinof{$key}}){$name=$species_name{$incspecies_withinof{$key}};}
                if(exists $incspecies_name{$incspecies_withinof{$key}}){$name=$incspecies_name{$incspecies_withinof{$key}};}
		printf OUT2 "$key\t$incspecies_name{$key}\t$incspecies_withinof{$key}\t$name\n";
	  }
	}
	close(OUT);
	close(OUT2);
  }

  if($runtype eq "-m") {#parse functional compartments
#two types of modules printing: 
#	1) prints all entities (by species and complex alias) to facilitate search of additional connections within or between modules (include localization) OUT3
#	2) prints only entities with species or complex species ID that are included in reactions (include localization) OUT5

          my $modules=substr($ARGV[1],0,rindex($ARGV[1],".")).".all_in_fm.txt";
	  my $modules2=substr($ARGV[1],0,rindex($ARGV[1],".")).".fm_reactions.txt";
	  my $modules3=substr($ARGV[1],0,rindex($ARGV[1],".")).".fm.txt";
          open(OUT3, ">$modules");
	  open(OUT4, ">$modules2");
	  open(OUT5, ">$modules3");
	  my(@react, %fmodule, %fmodulen,%fmofsp,%fmofspn);
	  foreach my $w (keys %specalias){#iterate species aliases including complex species
		if(exists $incspecies_withinof{$specalias{$w}}){
			$fmodulen{$w}=0;
			$fmodule{$w}="in complex";#TODO - resolve membership of inlcuded species. 
		}
		else{
			my $mn=$species_localization{$specalias{$w}};
			$fmodule{$w}=$compartment_name{$mn};
			if($mn eq "default"){$fmodulen{$w}=$nl+1;}
			else{
				$mn=~s/\D//;
				$fmodulen{$w}=$nl+1+$mn;
			}
		}
		#print "$species_localization{$specalias{$w}}\n";
		my $sx = $spx{$w};
		my $sy = $spy{$w};
#print "dupsko $w $specalias{$w} $sx $sy\n";
		for my $l(1..$nl){
		  my $mx1 = $fcompx1[$l];
		  my $mx2 = $fcompx2[$l];
		  my $my1 = $fcompy1[$l];
		  my $my2 = $fcompy2[$l];
			  if($sx > $mx1 and $sx < $mx2){
				if($sy > $my1 and $sy < $my2){
					 $fmodule{$w}=$fcomp_name[$l];
					 $fmodulen{$w}=$l;
				}
			  }	
		}
		if(exists $fmofsp{$specalias{$w}}){
			my $old=$fmofsp{$specalias{$w}};
			my $oldn=$fmofspn{$specalias{$w}};
			my $new=$old."-".$fmodule{$w};
			my $newn=$oldn."-".$fmodulen{$w};
			$fmofsp{$specalias{$w}}=$new;
			$fmofspn{$specalias{$w}}=$newn;
	  	}
	  	else{
			$fmofsp{$specalias{$w}}=$fmodule{$w};
			$fmofspn{$specalias{$w}}=$fmodulen{$w};
	  	}
	  	my $name = "noname";
	  	if(exists $species_name{$specalias{$w}}){$name=$species_name{$specalias{$w}};}
	  	if(exists $incspecies_name{$specalias{$w}}){$name=$incspecies_name{$specalias{$w}};}
	  	printf OUT3 "$specalias{$w}\t$w\t$name\t$fmodule{$w}\t$fmodulen{$w}\n";
	  	for my $i ( 0 .. $#{$SAiR{$w}} ) {
#my @k = split('\+',$w);
#print "$w $k[0] $i $SAiR{$w}[$i] $#{$SAiR{$w}}\n";
	      		print OUT4 "$w\t$SAiR{$w}[$i]\t$fmodule{$w}\t$fmodulen{$w}\n";
	  	}
	  }
	  close(OUT3);
	  close(OUT4);
	  foreach my $w (keys %fmofsp){#iterate species (including complexes)
	    if(exists $SiR{$w}){
		my $name = "noname";
		if(exists $species_name{$w}){$name=$species_name{$w};}
		printf OUT5 "$w\t$name\t$fmofsp{$w}\t$fmofspn{$w}\n";
	    }
	  }
  }

}#end of species parsing


sub open_file{
        my ($file_name)=@_;
        open(INP1, "< $file_name") or die "Can not open an input file: $!";
        my @file1=<INP1>;
        close (INP1);
        chomp @file1;
        return @file1;
}
