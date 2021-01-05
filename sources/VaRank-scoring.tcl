############################################################################################################
# VaRank 1.4.4                                                                                             #
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


## Return the bonus value to add to the VaRank score
proc addBonus {siftPred siftMed pph2 phastcons} {

    global g_VaRank

    if {[info exists g_VaRank(DEBUG)]} {puts "Add Bonus values: $siftPred $siftMed $pph2 $phastcons"}
    
    set bonus 0
    
    if {[info exists g_VaRank(snpeffDir)]} {
	if {[regexp -nocase "deleterious" $siftPred]} {
	    incr bonus $g_VaRank(B_SIFT)
	}
	if {[regexp "probably damaging|deleterious" $pph2]} {
	    incr bonus $g_VaRank(B_PPH2)
	}	
        if {[regexp "probablydamaging|deleterious" $pph2]} {
            incr bonus $g_VaRank(B_PPH2)
        }
    } else {
	if {[string is double $siftMed]} {
	    if {[regexp "Deleterious" $siftPred] && $siftMed >= 2.75 && $siftMed <= 3.5} {
		incr bonus $g_VaRank(B_SIFT)
	    }
	}
	if {[regexp "deleterious" $pph2]} {
	    incr bonus $g_VaRank(B_PPH2)
	}
    }
    if {[string is double $phastcons] && $phastcons > 0.95} {incr bonus $g_VaRank(B_phastCons)}
    
    if {[info exists g_VaRank(DEBUG)]} {puts "Bonus $bonus"}
    
    return $bonus
}

proc scoreAllTheID {} {

    global g_VaRank
    global g_ANNOTATION
    global g_PPH2
    global g_lScore
    global g_deltaSSF
    global g_deltaMES
    global g_deltaNNS
    global env

    # WARNING: If alamut annotation files of different versions have been merged (with different number of annotation columns), 
    # it implies a bug here. Correction: remove all alamut annotation files and rerun VaRank.
    

    ## Loading PPH2 data (neutral or deleterious or unknown) in the pph2($access,$pos,$aa1,$aa2) variable
    #####################################################################################################
    set pph2HumVarOutput "$g_VaRank(vcfDir)/PPH2/PPH2humVar_all.txt"   
    if {[file exists $pph2HumVarOutput]} {	
	puts "...loading PPH2 annotation ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"	
	set PPH2Header 0
	foreach L [LinesFromFile $pph2HumVarOutput] {
	    if {[regexp "^#o_acc" $L]} {
		set L [split $L "\t"]
		set i_pos      [lsearch -regexp $L "o_pos"];      if {$i_pos      == -1} {puts "o_pos: column number not found. Something must be wrong with PPH2 run - Exit"; exit}
		set i_aa1      [lsearch -regexp $L "o_aa1"];      if {$i_aa1      == -1} {puts "o_aa1: column number not found. Something must be wrong with PPH2 run - Exit"; exit}
		set i_aa2      [lsearch -regexp $L "o_aa2"];      if {$i_aa2      == -1} {puts "o_aa2: column number not found. Something must be wrong with PPH2 run - Exit"; exit}
		set i_pph2     [lsearch -regexp $L "pph2_class"]; if {$i_pph2     == -1} {puts "pph2class: column number not found. Something must be wrong with PPH2 run - Exit"; exit}
		set PPH2Header 1		
		continue
	    }
	    set L [split $L "\t"]
	    
	    if {$PPH2Header!="1"} {puts "PPH2 output header not found. Something must be wrong with PPH2 run - Exit"; exit}
	    
	    regsub -all " " [lindex $L 0] "" access
	    ## We remove the version number at the end of the ID (eg: refseqP ID = NP_001185763.1)
	    regsub "\\..*$" $access "" access
	    regsub -all " " [lindex $L $i_pos]  "" pos
	    regsub -all " " [lindex $L $i_aa1]  "" aa1
	    regsub -all " " [lindex $L $i_aa2]  "" aa2
	    regsub -all " " [lindex $L $i_pph2] "" actpph2

	    if {![info exists pph2($access,$pos,$aa1,$aa2)] || $pph2($access,$pos,$aa1,$aa2) == "unknown" || $actpph2 == "deleterious"} {
		set pph2($access,$pos,$aa1,$aa2) $actpph2
	    }
	}
    }
	
    ## Define the "g_lScore" global variable with the genetic variants analysed 
    ###########################################################################
    puts "...scoring each genetic variant ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
    #regsub -all "^\\\\?{|\\\\?}$" [split $g_ANNOTATION(#id) "\t"] "" L
    set L [split $g_ANNOTATION(#id) "\t"]
    
    # for transfert of pph2 annotation in the  "g_PPH2" global variable
    set i_refseq    [lsearch -regexp $L "^protein$"]; # This annotation doesn't exists with SNPEFF annotation.
    set i_uniprot   [lsearch -regexp $L "^Uniprot$"]; if {$i_uniprot == -1} {puts "column number not found for uniprot - Exit"; exit}
    set i_pos       [lsearch -regexp $L "^posAA$"];   if {$i_pos     == -1} {puts "column number not found for posAA - Exit"; exit}
    set i_AA        [lsearch -regexp $L "^wtAA_1$"];  if {$i_AA      == -1} {puts "column number not found for wtAA_1 - Exit"; exit}
    set i_var       [lsearch -regexp $L "^varAA_1$"]; if {$i_var     == -1} {puts "column number not found for varAA_1 - Exit"; exit}

    # for scoring
    set i_gene      [lsearch -exact $L "gene"];           if {$i_gene   == -1} {puts "column number not found for gene - Exit"; exit}
    set i_trans     [lsearch -exact $L "transcript"];     
    set i_length    [lsearch -exact $L "transLen"];      
    set i_clinic    [lsearch -exact $L "rsClinicalSignificance"];
    set i_clinic2   [lsearch -exact $L "clinVarClinSignifs"];
    set i_effect    [lsearch -exact $L "codingEffect"];   if {$i_effect == -1} {puts "column number not found for codingEffect - Exit"; exit}
    set i_varType   [lsearch -regexp $L "^varType$"];     if {$i_varType  == -1} {puts "column number not found for varType - Exit"; exit}

    set i_distSS    [lsearch -exact $L "distNearestSS"];  
    set i_nearSS    [lsearch -exact $L "nearestSSType"]; 
    set i_lse       [lsearch -exact $L "localSpliceEffect"];
    set i_varlocation    [lsearch -exact $L "varLocation"];    if {$i_varlocation == -1} {puts "column number not found for varLocation - Exit"; exit}
    set i_wtSSF     [lsearch -exact $L "wtSSFScore"];    
    set i_wtMES     [lsearch -exact $L "wtMaxEntScore"]; 
    set i_wtNNS     [lsearch -exact $L "wtNNSScore"];     
    set i_varSSF    [lsearch -exact $L "varSSFScore"];   
    set i_varMES    [lsearch -exact $L "varMaxEntScore"]; 
    set i_varNNS    [lsearch -exact $L "varNNSScore"];   

    set i_siftPred  [lsearch -exact $L "SIFTprediction"]; if {$i_siftPred  == -1} {puts "column number not found for SIFTprediction - Exit"; exit}
    set i_siftMed   [lsearch -exact $L "SIFTmedian"];     
    set i_phastcons [lsearch -exact $L "phastCons"];      if {$i_phastcons == -1} {puts "column number not found for phastCons - Exit"; exit}

    #set i_phylop [lsearch -exact $L "phyloP"]; if {$i_phylop == -1} {puts "column number not found for phyloP - Exit"; exit}

    
    foreach ID [array name g_ANNOTATION] {
	if {$ID == "#id"} {continue}
	set score      0
	set bestScore  0
	set bestLength 0
	set allGene   {}
	set bestGene ".."

	if {[info exists g_VaRank(DEBUG)]} {puts "Variant: $ID"}

	# For a single variant, multiple isoforms can give different score. Each is evaluated.
	set L_pathogenic {}
	foreach L $g_ANNOTATION($ID) {
	    set L [split $L "\t"]

	    # Reset score between isoforms
	    set score 0

	    set trans     ""
	    set length    ""
	    set gene      ""
	    set effect    ""
	    set varType   ""
	    set clinic    ""
	    set siftPred  ""
	    set siftMed   ""
	    set thePPH2   "NA"
	    set phastcons ""


	    if {[info exists g_VaRank(DEBUG)]} {puts "Data: $L"}
	    
	    set varType   [lindex $L $i_varType]
	    set trans  [lindex $L $i_trans]  

	    set length [lindex $L $i_length]   
	    set gene   [lindex $L $i_gene]     
	    if {![regexp -nocase "$gene" $allGene]} {lappend allGene $gene}
	    
	    set effect [lindex $L $i_effect]   
	    set clinic [lindex $L $i_clinic]   
	    append clinic " [lindex $L $i_clinic2]"

	    set siftPred [lindex $L $i_siftPred] 
	    set siftMed  [lindex $L $i_siftMed]  
	    set phastcons [lindex $L $i_phastcons]
	    set varlocation [lindex $L $i_varlocation]


	    # General scoring scheme for big categories

	    if {[regexp "stop gain|stop_gained" $effect]} {
		set score $g_VaRank(S_StopGain)
		incr score [addBonus "NA" "NA" "NA" $phastcons]
	    } elseif {[regexp -nocase "frameshift" $effect]} {
		set score $g_VaRank(S_Fs)
	    } elseif {[regexp "startloss|start_lost|start loss" $effect]} {
		set score $g_VaRank(S_StartLoss)
		incr score [addBonus "NA" "NA" "NA" $phastcons]
	    } elseif {[regexp "stoploss|stop_lost|stop loss" $effect]} {
		set score $g_VaRank(S_StopLoss)
		incr score [addBonus "NA" "NA" "NA" $phastcons]
	    } elseif {[regexp "missense|rare_amino_acid" $effect] && [regexp -nocase "substitution|SNP" $varType]} {
		if {[file exists $pph2HumVarOutput]} {	    		    
		    # Report of some "$ID $codingEffect $varType" values:
		    #####################################################	    
		    # 1_909238_G_C missense substitution
		    # 1_909242_A_G synonymous substitution
		    # 1_911595_A_G NA substitution
		    # 1_12888523_T_G start loss substitution
		    # 10_27687225_A_G stop loss substitution
		    # 11_20949960_C_T missense substitution
		    # 8_142446923_C_T missense substitution
		    # 1_914333_C_G missense substitution
		    #
		    # Some CNV from Pindel have a missense coding effect, but are not substitution (SNV). For example:
		    # id	                          chrom  pos      gene	    protein	        Uniprot	 varType     codingEffect
		    # 19_4511674_CACAGCA_TACGGTG   19     4511674  PLIN4     NP_001073869.1	Q96Q06	 delins	     missense   
		    #
		    # Selection of the lines with: ("codingEffect" = missense) && ("varType" = substitution)  => ("varLocation" = exon) 
		    set IDuniprot [lindex $L $i_uniprot]
		    set Pos       [lindex $L $i_pos]
		    set AA1       [lindex $L $i_AA]
		    set Var       [lindex $L $i_var]
		    if {$i_refseq ne "-1"} {
			set IDrefseqp_V [lindex $L $i_refseq]
			regsub "\\..*$" $IDrefseqp_V "" IDrefseqp
		    } else {
			set IDrefseqp_V "NA"
			set IDrefseqp   "NA"
		    }
		    regsub "\\..*$" $IDrefseqp_V "" IDrefseqp		    
		    if {[info exists pph2($IDuniprot,$Pos,$AA1,$Var)]} {
			set thePPH2 $pph2($IDuniprot,$Pos,$AA1,$Var)
		    } elseif {[info exists pph2($IDrefseqp_V,$Pos,$AA1,$Var)]} {
			set thePPH2 $pph2($IDrefseqp_V,$Pos,$AA1,$Var)
		    } elseif {[info exists pph2($IDrefseqp,$Pos,$AA1,$Var)]} {
			set thePPH2 $pph2($IDrefseqp,$Pos,$AA1,$Var)
		    } 
		}		
		set score $g_VaRank(S_Missense)
		# Add SIFT & PPH2 bonus for missense only.
		incr score [addBonus $siftPred $siftMed $thePPH2 $phastcons]
	    } elseif {[regexp "in-frame|inframe" $effect]} {
		set score $g_VaRank(S_Inframe)
	    } elseif {[regexp "synonymous|start_retained|stop_retained" $effect]} {
		set score $g_VaRank(S_Synonymous)
		incr score [addBonus "NA" "NA" "NA" $phastcons]
	    } elseif {[regexp "intron|exon" $varlocation]} {
		set score $g_VaRank(S_ExonIntron)
	    } elseif {[regexp "UTR" $varlocation]} {
		set score $g_VaRank(S_UTR)
	    }

	    if {[info exists g_VaRank(DEBUG)]} {puts "Score $score"}



	    
	    # Scoring for splice effect, only with ALAMUT annotation
	    ########################################################
	    if {[info exists env(ALAMUT)]} {
		
		set distSS [lindex $L $i_distSS]
		set nearSS [lindex $L $i_nearSS] 
		set lse    [lindex $L $i_lse]
		set wtSSF  [lindex $L $i_wtSSF] 
		set wtMES  [lindex $L $i_wtMES] 
		set wtNNS  [lindex $L $i_wtNNS] 
		set varSSF [lindex $L $i_varSSF]
		set varMES [lindex $L $i_varMES] 
		set varNNS [lindex $L $i_varNNS] 
		
		set n   0
		set tot 0
		
		
		if {[info exists g_VaRank(DEBUG)]} {puts "Splice Scores: $wtSSF, $wtMES, $wtNNS, $varSSF, $varMES, $varNNS"}
		
		
		if {$wtSSF != "" && $wtSSF != "NA" && $wtSSF != "0"} {
		    if {$varSSF=="" || $varSSF=="NA"} {set varSSF 0}
		    set  g_deltaSSF($ID,$trans) [format "%.1f" [expr {($varSSF-$wtSSF)/$wtSSF*100.0}]]
		    if {$g_deltaSSF($ID,$trans) < $g_VaRank(SSFcutoff)} {incr n}
		    incr tot
		} else {set g_deltaSSF($ID,$trans) "NA"}
		
		if {[info exists g_VaRank(DEBUG)]} {puts "g_deltaSSF $trans [set g_deltaSSF($ID,$trans)]"}
		
		if {$wtMES != "" && $wtMES != "NA" && $wtMES != "0"} {
		    if {$varMES == "" || $varMES=="NA"} {set varMES 0}
		    
		    set  g_deltaMES($ID,$trans) [format "%.1f" [expr {($varMES-$wtMES)/$wtMES*100.0}]]
		    if {$g_deltaMES($ID,$trans) < $g_VaRank(MEScutoff)} {incr n}
		    incr tot
		} else {set g_deltaMES($ID,$trans) "NA"}
		
		if {[info exists g_VaRank(DEBUG)]} {puts "g_deltaMES $trans [set g_deltaMES($ID,$trans)]"}

		if {$wtNNS != "" && $wtNNS != "NA"  && $wtNNS != "0"} {
		    if {$varNNS == "" || $varNNS=="NA"} {set varNNS 0}
		    
		    set  g_deltaNNS($ID,$trans) [format "%.1f" [expr {($varNNS-$wtNNS)/$wtNNS*100.0}]]
		    if {$g_deltaNNS($ID,$trans) < $g_VaRank(NNScutoff)} {incr n}
		    incr tot
		} else {set g_deltaNNS($ID,$trans) "NA"}
		
		if {[info exists g_VaRank(DEBUG)]} {puts "g_deltaNNS $trans [set g_deltaNNS($ID,$trans)]"}
		
		set splicingEffect 0
		if {$tot != 0 && [expr {$n*1.0/$tot}] >= 0.5} {set splicingEffect 1}
		
		if {[info exists g_VaRank(DEBUG)]} {puts "splicingEffect $splicingEffect"}

		if {$splicingEffect} {
		    if {$score <= $g_VaRank(S_EssentialSplice)} {
			# Essential splice site
			if {$nearSS == "5'" && [regexp "^1$|^2$" $distSS]} {
			    set  score $g_VaRank(S_EssentialSplice)
			    incr score [addBonus "NA" "NA" "NA" $phastcons]
			} elseif {$nearSS == "3'" && [regexp "^-1$|^-2$" $distSS]} {
			    set  score $g_VaRank(S_EssentialSplice)
			    incr score [addBonus "NA" "NA" "NA" $phastcons]
			}
		    } 
		    if {$score <= $g_VaRank(S_CloseSplice)} {
			# Intron-Exon boundary
			if {$nearSS == "5'" && $distSS >= -3 && $distSS <= 6} {
			    set  score $g_VaRank(S_CloseSplice)
			    incr score [addBonus "NA" "NA" "NA" $phastcons]
			} elseif {$nearSS == "3'" && $distSS >= -12 && $distSS <= 2} {
			    set  score $g_VaRank(S_CloseSplice)
			    incr score [addBonus "NA" "NA" "NA" $phastcons]
			}
		    } 
		    if {$score <= $g_VaRank(S_DeepSplice)} {
			# Deep intron-exon boundary
			if {$varlocation == "intron"} {
			    set score $g_VaRank(S_DeepSplice)
			    incr score [addBonus "NA" "NA" "NA" $phastcons]
			}
		    }
		}
		## Use of "LocalSpliceEffect" field from Alamut (even if $splicingEffect == 0):
		## VaRank score = 40 : New Donor Site, New Acceptor Site, Cryptic Donor Strongly Activated, Cryptic Acceptor Strongly Activated,
		## VaRank score = 35 : Cryptic Donor Weakly Activated, Cryptic Acceptor Weakly Activated
		if {[regexp "New|Strongly" $lse] && $score<$g_VaRank(S_LSEstrong)} {set score $g_VaRank(S_LSEstrong); incr score [addBonus "NA" "NA" "NA" $phastcons]}
		if {[regexp "Weakly" $lse] && $score<$g_VaRank(S_LSEweak)} {set score $g_VaRank(S_LSEweak); incr score [addBonus "NA" "NA" "NA" $phastcons]}
	    }

	    if {[info exists g_VaRank(DEBUG)]} {puts "Score $ID $trans: $score / Best $bestScore"}
	    
	    # Within the different isoforms we keep the one for which the effect is maximum
	    # If similar score between isoforms we keep the longest one

	    # puts "$ID $gene $trans - $score - $bestScore"

	    if {$score > $bestScore} {

		#puts "Better"

		set bestGene   $gene
		set bestL      $L
		set bestScore  $score
		set bestLength $length

		set g_PPH2($ID) $thePPH2
	
	    } elseif {$score == $bestScore} {
		#puts "equal"

		if {$length > $bestLength || $length eq ""} {
		    #puts "Longer"
		    set bestGene   $gene
		    set bestL      $L
		    set bestLength $length

		    set g_PPH2($ID) $thePPH2
		   
		}
	    }

	    # 2016/02/29
	    #
	    # Values for rsClinicalSignificance:
	    # ----------------------------------
	    # not_provided, pathogenic, likely_benign, benign, likely_pathogenic, uncertain_significance, risk_factor, association, protective, other, drug_response, confers_sensitivity
	    # 
	    # Values for clinVarClinSignifs:
	    # ------------------------------
	    #Uncertain significance, Likely benign, Likely pathogenic, Benign, Pathogenic, association, not provided, drug response, other, protective, risk factor, conflicting data from submitters, confers sensitivity, Conflicting interpretations of pathogenicity, Affects
	    
	    # To notice: With SnpEff annotations, clinic = ""

	    # We explicitely remove the conflicting ones and highlight whenever there is pathogenic or likely pathogenic annotations
	    set clinicChanged ""
	    regsub -all "non-pathogenic" $clinic "" clinicChanged
	    regsub -all -nocase "conflicting interpretations of pathogenicity" $clinicChanged "" clinicChanged

	    if {[regexp -nocase "pathogenic" $clinicChanged]} {
		lappend L_pathogenic "$gene $score $length"
		set LineOfThePathogenicGene($gene,$length) $L
	        set thePPH2OfThePathogenicGene($gene,$length) $thePPH2
	    } 


	}

	if {$L_pathogenic ne {}} {
	    set L_pathogenic [lsort -command DescendingSortOnElement2 $L_pathogenic]
	    set L_pathogenic [lsort -command DescendingSortOnElement1 $L_pathogenic]
	    set top [lindex $L_pathogenic 0]
	    set bestGene [lindex $top 0]
	    set bestScore $g_VaRank(S_Known)
	    set bestLength [lindex $top 2]
            set g_PPH2($ID) $thePPH2OfThePathogenicGene($bestGene,$bestLength)
	    array unset thePPH2OfThePathogenicGene "*"
	    set bestL $LineOfThePathogenicGene($bestGene,$bestLength)
	    array unset LineOfThePathogenicGene "*"
	} 

	if {$bestGene != $allGene} {
	    regsub "$bestGene" $allGene "" allGene
	    set allGene "$bestGene/[join $allGene "/"]"
	}

	#puts "$allGene"


	set g_ANNOTATION($ID) [join [lreplace $bestL $i_gene $i_gene "$allGene"] "\t"]
	
	if {[info exists g_VaRank(DEBUG)]} {puts "Best $ID lreplace $bestL $i_gene $i_gene $allGene"}
	
	lappend g_lScore "$ID $bestScore"

    }

    puts "...classifying each genetic variant ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"

    if {![info exists g_lScore] || $g_lScore == {}} {
	puts "No score could be retrieved something must be wrong with variation identifiers (g_lScore) - Exit";exit
	if {[info exists g_VaRank(DEBUG)]} {
	    puts "List of uploaded variation from alamut: [array name g_ANNOTATION]"
	    puts "g_lScore $g_lScore"
	}
    }
    
    set g_lScore [lsort -command DescendingSortOnElement1 $g_lScore]

    return
}



