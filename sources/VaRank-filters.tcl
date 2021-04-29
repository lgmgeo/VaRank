############################################################################################################
# VaRank 1.5.2                                                                                             #
#                                                                                                          #
# VaRank: a simple and powerful tool for ranking genetic variants                                          #
#                                                                                                          #
# Copyright (C) 2016-2021 Veronique Geoffroy (veronique.geoffroy@inserm.fr)                                #
#                         Jean Muller (jeanmuller@unistra.fr)                                              #
#                                                                                                          #
# Please cite the following article:                                                                       #
# Geoffroy V.*, Pizot C.*, Redin C., Piton A., Vasli N., Stoetzel C., Blavier A., Laporte J. and Muller J. #
# VaRank: a simple and powerful tool for ranking genetic variants.                                         #
# PeerJ. 2015 (10.7717/peerj.796)                                                                          #
#                                                                                                          #
# This is part of VaRank source code.                                                                      #
#                                                                                                          #
# This program is free software; you can redistribute it and/or                                            #
# modify it under the terms of the GNU General Public License                                              #
# as published by the Free Software Foundation; either version 3                                           #
# of the License, or (at your option) any later version.                                                   #
#                                                                                                          #
# This program is distributed in the hope that it will be useful,                                          #
# but WITHOUT ANY WARRANTY; without even the implied warranty of                                           #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                            #
# GNU General Public License for more details.                                                             #
#                                                                                                          #
# You should have received a copy of the GNU General Public License                                        #
# along with this program; If not, see <http://www.gnu.org/licenses/>.                                     #
############################################################################################################
 
##
## Filtering of all the existing ranking files (not already filtered).
##
## OUTPUTS: 	g_VaRank(vcfDir)/"family"_"patient"_filteredVariants.rankingByVar.tsv (1 par patient)
##

proc executeFilters {} {

    global g_VaRank
    global g_vcfINFOS
    global g_ANNOTATION
    global env 

    puts "...writing all filtered ranking files ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"

    ## Remove all old filtered files
    ################################
    foreach filteredFile [glob -nocomplain $g_VaRank(vcfDir)/*filteredVariants*.tsv] {file delete -force $filteredFile}
    
    ## List all frequency colum names
    #################################
    set L_freqColumns {}
    foreach colName $g_VaRank(L_outputColHeader) {
	if {[regexp "AF$|Freq|FREQ" $colName] && ![regexp "alleleFrequency|rsMAF" $colName]} {
	    lappend L_freqColumns $colName
	}
    }

    ## Filter every line with the user's options.
    ## Create the filtered files.
    #############################################################
    foreach rankFile [glob -nocomplain $g_VaRank(vcfDir)/*Variants.rankingBy*.tsv] {
	regsub "all(Variants)" $rankFile "filtered\\1" outputfile
	set OutputText ""
	foreach L [LinesFromFile $rankFile] {
	    if {[regexp "^## Barcode:" $L] || $L == ""} {append OutputText $L; continue}
	    if {[regexp "^##"          $L] || $L == ""} {append OutputText "\n$L"; continue}
	    set Ls [split $L "\t"]

	    if {[regexp "^variantID" $L]} {
		append OutputText "\n$L"

		set i_clinic          [lsearch -exact $Ls "^rsClinicalSignificance"];
		set i_clinic2         [lsearch -exact $Ls "clinVarClinSignifs"];      
		set i_valid           [lsearch -exact $Ls "^clinVarReviewStatus"];  
		#set i_valid           [lsearch -regexp $Ls "^rsValidations"];          

		set i_cov             [lsearch -exact $Ls "^totalReadDepth"];  
		set i_read            [lsearch -exact $Ls "^varReadDepth"];    
		set i_percent         [lsearch -exact $Ls "^varReadPercent"]; 
		
		foreach colName $g_VaRank(L_outputColHeader) {
		    set iCol    "i_$colName"
		    set $iCol [lsearch -exact $Ls "$colName"]
		}

		continue        
	    }
	    set cov     [lindex $Ls $i_cov] 
	    set read    [lindex $Ls $i_read]  
	    set percent [lindex $Ls $i_percent]

	    # Filtering on read coverage
	    if {$cov != "" && $cov != "NA" && $cov != "." && $cov != "0"} { 
		if {$cov < $g_VaRank(depthFilter)} {continue}

		if {$read != "" && $read != "NA" && $read != "."} {
		    if {$read    < $g_VaRank(readFilter)}        {continue}
		    if {$percent < $g_VaRank(readPercentFilter)} {continue}
		}
	    }

	    # Filtering on frequency data
	    if {$g_VaRank(freqFilter) != ""} {
		set TooFrequent 0

		# rsMAF
		if {[info exists rsMAFAllele] && $i_rsMAF ne -1} {			
		    if {$rsMAFAllele eq $alt} {		   
			set freq [lindex $Ls $i_rsMAF] 
			if {[set g_VaRank(metrics)]=="fr"} {
			    #puts "Filters replacing $freq"
			    regsub {\,} $freq "." freq
			    #puts "Filters: $freq"
			}
			if {$freq != "NA" && $freq != "" && $freq>=$g_VaRank(freqFilter)} {
			    continue
			}
		    }
		}

		# gnomadFilter 
		if {[info exists i_gnomadFilter]} {
		    set gnomadFilter [lindex $Ls $i_gnomadFilter]
		} else {
		    set gnomadFilter "PASS"; # If I haven't this information, I can't filter on it
		}

		# filtering with all freq
		foreach colName $L_freqColumns {		  
		    # Don't use gnomad frequencies if not "PASS"
		    if {[regexp "^gnomad" $colName] && $gnomadFilter ne "PASS"} {continue}
		    # Don't use "EVS/ESP" frequencies (none column with these frequencies)
		    set iCol "i_$colName"	    	    
		    if {[set $iCol] eq -1} {
			continue
		    } else {
			set freq [lindex $Ls [set $iCol]] 
		    }	    
		    if {[set g_VaRank(metrics)]=="fr"} {
			#puts "Filters replacing $freq"
			regsub {\,} $freq "." freq
			#puts "Filters: $freq"
		    }

		    if {$freq != "NA" && $freq != "" && $freq>=$g_VaRank(freqFilter)} {set TooFrequent 1;break}
		}
		if {$TooFrequent} {continue}
	    }

	    #Filtering on sequencing data
	    #
	    #Here we need to improve filter and add also frequency
	    if {[info exists env(ALAMUT)]} {
		if {$g_VaRank(rsFilter) == "removeNonPathoRS"} {
		    set    valid  [lindex $Ls $i_valid]  
		    set    clinic [lindex $Ls $i_clinic] 
		    append clinic " [lindex $Ls $i_clinic2]"
		    
		    # Vero, 2014/11/07  ----  "Clinical significance"
		    # Alamut developpers communication:   
		    # "Depuis le 9 septembre 2014, notre source de données dbSNP n’est plus le NCBI mais Ensembl.
		    # En raison de problèmes de qualité rencontrés dans les fichiers XML fournis par le NCBI, nous avons décidé d’utiliser l’API d’Ensembl comme source de données dbSNP."
		    # La signification clinique fournie par Ensembl correspond à la signification clinique donnée par la base de données ClinVar (NCBI).
		    #
		    # http://www.ensembl.org/info/genome/variation/data_description.html#clin_significance
		    # http://www.ncbi.nlm.nih.gov/clinvar/docs/clinsig/
		    #
		    # [association] 	         association	
		    # [benign] 	                 benign	        
		    # [confers sensitivity] 	 confers sensitivity
		    # [drug response] 	 	 drug response
		    # [likely benign] 	 	 likely benign	
		    # [likely pathogenic]  		 likely pathogenic
		    # [not provided] 	         not provided	
		    # [other] 	       		 other
		    # [pathogenic] 	                 pathogenic
		    # [protective] 	        	 protective
		    # [risk factor] 	         risk factor
		    # [uncertain significance] uncertain significance
		    
		    # 2015/08/24
		    # WARNING : We can have both benign and pathogenic annotation in the "valid" field:
		    # 22_51064039_G_C benign,pathogenic -- Cluster/Frequency/Submitter/DoubleHit/HapMap/1000G
		    # 7_151478406_C_T benign,likelypathogenic -- Cluster/Frequency/1000G
		    # 12_121437382_A_G benign,pathogenic -- Cluster/Frequency/HapMap/1000G
		    # 3_8775661_C_T notprovided,benign,pathogenic -- Cluster/Frequency/Submitter/DoubleHit/HapMap/1000G
		    
		    # 2016/02/29
		    # Values for rsClinicalSignificance
		    # not_provided, pathogenic, likely_benign, benign, likely_pathogenic, uncertain_significance, risk_factor, association, protective, other, drug_response, confers_sensitivity
		    # Values for clinVarClinSignifs
		    # Uncertain significance, Likely benign, Likely pathogenic, Benign, Pathogenic, association, not provided, drug response, other, protective, risk factor, conflicting data from submitters, confers sensitivity, Conflicting interpretations of pathogenicity, Affects
		    #
		    # We explicitely remove the conflicting ones and highligh whenever there is pathogenic or likely pathogenic annotations
		    
		    # Excluding validated benign variants (Keeping the potential pathogenic ones)
		    #
		    #regsub -all -nocase "conflicting interpretations of pathogenicity" $clinicChanged "" clinicChanged
		    #if {[regexp "benign|^NA$" $clinic] && ![regexp "pathogenic" $clinic]} {}

		    if {[regexp -nocase "benign" $clinic] && ![regexp -nocase "pathogenic|risk" $clinic]} {
			#Keeping only the ones with 2 validation or more
			#if {[regsub -all "/|;" $valid "" toto] >= 1 || [regexp -nocase "^yes$" $valid]} {continue}
			if {$valid > 1} {continue}
		    }
		}
	    }

	    append OutputText "\n$L"
	}
	if {[regexp "rankingByVar" $outputfile]} {
	    ReplaceTextInFile $OutputText $outputfile
	    continue
	}

	## Ranking by gene have to be done once again after filtering
	regexp "fam\[0-9\]+_(.+)_allVariants.rankingBy" $rankFile match patient
	set OutputText2 ""
	set lGenes ""
	set lVariantsRankingByGene ""
	set lID ""
	foreach L [split $OutputText "\n"] {
	    if {$L == ""} {continue}
	    if {[regexp "^## Barcode:" $L]} {append OutputText2     $L; continue}
	    if {[regexp "^##"          $L]} {append OutputText2 "\n$L"; continue}
            set Ls [split $L "\t"]

	    if {[regexp "^variantID" $L]} {
		append OutputText2 "\n$L"
		set i_gene  [lsearch -exact $Ls "gene"];           if {$i_gene  == -1} {puts "gene: column number not found - Exit"; exit}
		set i_score [lsearch -exact $Ls "varankVarScore"]; if {$i_score == -1} {puts "varankVarScore: column number not found - Exit"; exit}
		continue
	    }
	    set id   [lindex $Ls 0]     
	    set gene [lindex $Ls $i_gene] 
	    if {$gene == "NA"} {continue}

	    if {[lsearch -exact $lGenes "$gene"] == -1} {
		lappend lGenes $gene
	    }
	    set l($id) $L
		
	    set score [lindex $Ls $i_score] 
	    if {$score eq ""} {continue} 
 
	    lappend lVariants($gene) "$id $score"
	}
	
	foreach gene $lGenes {
	    set liste {}
	    foreach duoIDScore $lVariants($gene) {
		set id [lindex $duoIDScore 0]
		lappend liste $id
	    }
	    set maxScore1 [lindex [lindex $lVariants($gene) 0] 1]
	    if {[llength $lVariants($gene)] == 1} {
		set BestDuoScore "[expr {$maxScore1*2}]"
	    } else {
		set bestID [lindex [lindex $lVariants($gene) 0] 0]
		if {[regexp "$patient:(\[^: \]+):" $g_vcfINFOS($bestID) match homhtz] && [regexp -nocase "hom" $homhtz]} {
		    set BestDuoScore "[expr {$maxScore1*2}]"
		} else {
		    set maxScore2    [lindex [lindex $lVariants($gene) 1] 1]
		    set BestDuoScore "[expr {$maxScore1+$maxScore2}]"
		}
	    }
	    lappend lVariantsRankingByGene "{$liste} $BestDuoScore"
	}
	if {[info exists lVariantsRankingByGene]} {
	    set lVariantsRankingByGene [lsort -command DescendingSortOnElement1 $lVariantsRankingByGene]
	}    
	foreach val $lVariantsRankingByGene {
	    set lID [lindex $val 0]
	    foreach id $lID {
		append OutputText2 "\n$l($id)"
	    }
	}

	ReplaceTextInFile "$OutputText2" $outputfile
	
	## Very important 
	if {[info exists l]} {array unset l "*"}
	if {[info exists lVariants]} {array unset lVariants "*"}
    }
    return
}


