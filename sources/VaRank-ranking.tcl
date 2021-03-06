############################################################################################################
# VaRank 2.0                                                                                               #
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


proc ExternalAnnotations args {

    #To add external annotation to a gene in particular
    #
    # Format is tab separated values, 1st line is a header, 1st "genes" column is the gene name, rest is free
    # Typical use would be a gene file containing specific annotations such as tranmission mode, disease, expression...
    # WARNING: shouldn't contain "{" and "}". Replaced here by "(" and ")"
    
    # args :
    #   - $F
    #   - "L_Files"    -return-> Annotation files list
    #   - $F,"L_ID"    -return-> List of all genes from $F
    #   - $F,Header    -return-> Header from $F (without the "genes" column)
    #   - $F,$ID       -return-> Annotation for $ID (= gene) from $F ([lrange $L 1 end])
    
    global g_VaRank
    global g_ExtAnnotation

    set What [join $args ","]
    if {[info exists g_ExtAnnotation($What)]} {return [set g_ExtAnnotation($What)]}

    if {![info exists g_VaRank(extann)] || $g_VaRank(extann) eq ""} {return}
    
    set g_ExtAnnotation(L_Files) $g_VaRank(extann)

#    if {[file exists [lindex $args 0]]} {
#	set L_Files_Anno [lindex $args 0]
#    } else {
#	set L_Files_Anno [set g_VaRank(extann)]
#	set What [join [concat [list $L_Files_Anno] $args] ","]
 #   }
    
    #puts "--args $args"
    #puts "--$L_Files_Anno"
    #puts "--$What"

#    if {[info exists g_ExtAnnotation($What)]} {return [set g_ExtAnnotation($What)]}


    foreach File_Anno [set g_VaRank(extann)] {
	if {[info exists g_ExtAnnotation($File_Anno,Loaded)]} {continue}
	if {[regexp "README$|results.txt$" $File_Anno]} {continue}

	puts "...Loading [file tail $File_Anno] ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
	set g_ExtAnnotation($File_Anno,Loaded) 1
		
	if {![info exists g_ExtAnnotation($File_Anno,L_ID)]} {
	    set g_ExtAnnotation($File_Anno,L_ID) {}
	}
	
	set NbHeaderColumns 0
	set First 1
	if {[regexp ".gz$" $File_Anno]} {
	    set lLines [LinesFromGZFile $File_Anno]
	} else {
	    set lLines [LinesFromFile $File_Anno]
	}
	foreach L $lLines {
	    ## Bug during WritingByGene if "{ or }" are presents.
	    ## To be done here, inside each element of the list (not on "$lLines").
	    regsub -all "{" $L "(" L
	    regsub -all "}" $L ")" L
	    
	    if {![regexp -nocase {[a-z0-9]+} $L] || [regexp -nocase {^[\#]} $L]} {continue}
	    
	    if {$First} {
		set First 0
		if {![regexp -nocase "^( +)?gene(s)?$" [string trim [lindex [split $L "\t"] 0]]]} {
		    puts "\tReading $File_Anno:\n\tWARNING: Header ($L) does not contain a \"genes\" named column. Not used."
		    set g_ExtAnnotation($File_Anno,Header) ""
		    # => removed this file from the list 
		    set i [lsearch -exact $g_VaRank(extann) "$File_Anno"]
		    set g_VaRank(extann) [lreplace $g_VaRank(extann) $i $i]
		    break
		} else {
		    # Registering the header Without keeping the "genes" column name 
		    regsub "^( +)?gene(s)?\t" $L "" g_ExtAnnotation($File_Anno,Header)
		    set NbHeaderColumns [llength [split $L "\t"]]
		}
		continue
	    } 
	    
	    set Ls [split $L "\t"]
	    set ID [string trim [lindex $Ls 0]]
	    if {[llength $Ls] ne $NbHeaderColumns} {puts "WARNING: $ID, format (nb columns) is not good ([llength $L] vs $NbHeaderColumns columns)"}
	    
	    if {![info exists g_ExtAnnotation($File_Anno,$ID)]} {
		lappend g_ExtAnnotation($File_Anno,L_ID) $ID
	    } else {
		puts "WARNING: Reading Annotation file $File_Anno and $ID is seen multiple times."
	    }
	    set g_ExtAnnotation($File_Anno,$ID) [join [lrange $Ls 1 end] "\t"]
	}
	
	puts "\t...[llength [set g_ExtAnnotation($File_Anno,L_ID)]] gene identifiers and [llength [split [set g_ExtAnnotation($File_Anno,Header)] "\t"]] annotations columns ([join [split [set g_ExtAnnotation($File_Anno,Header)] "\t"] ", "])."
	
	if {[info exists g_ExtAnnotation($What)]} {return [set g_ExtAnnotation($What)]}
    }
    return ""
}


# DEFINITION OF:
################
# - g_perso($patient) "$zygosity\t$totalReadDepth\t$varReadDepth\t$varReadPercent\t$QUALphred"
# - g_famBarcode($fam)

# RETURN OF:
############
# "$barcode $homCount $hetCount $alleleCount $sampleCount $avgVariantDepth $sdVariantDepth $countVariantDepth $avgTotalDepth $sdTotalDepth $countTotalDepth"

# Barcode for each variation observed 
#####################################
# 0 homozygous/no information for the variation
# 1 heterozygous for the variation
# 2 homozygous for the variation

# Compute the statistics for coverage for all patients for one SNV
##################################################################
proc findBarcodesAndStatFor {ID} {
    global g_lPatientsOf
    global g_famBarcode
    global g_allPatients
    global g_perso
    global g_Statistics
    global g_VaRank

    # totalReadDepth = DP information from the VCF
    # varReadDepth = NR information from the VCF

    catch {array unset g_famBarcode "*"}
    catch {array unset g_perso "*"}
    

    # Definition of $g_perso($patient) and $zygous($patient)
    ########################################################
    set L_TotalDepth {}
    set L_VariantDepth   {}    
    foreach fam [array names g_lPatientsOf] {
	foreach patient $g_lPatientsOf($fam) {
	    set test [db_vcfData eval {SELECT homhet, dp, nr, qual FROM vcfSampleData AS vs \
					   INNER JOIN sampleName AS sn ON vs.sampleName_id = sn.sampleName_id \
					   WHERE ID = $ID AND sampleName = $patient}]
	    if {$test ne ""} {
		set zygosity       [lindex $test 0]
		set totalReadDepth [lindex $test 1]
		set varReadDepth   [lindex $test 2]
		set QUALphred      [lindex $test 3]

		##DEBUG
		#puts "$ID $zygosity\t$totalReadDepth\t$varReadDepth"

		#Depth of coverage
		if {$totalReadDepth ne "NA" && [regexp {^[0-9]} $totalReadDepth]} {lappend L_TotalDepth $totalReadDepth}
		if {$varReadDepth ne "NA" && [regexp {^[0-9]} $varReadDepth]}  {lappend L_VariantDepth   $varReadDepth}

		#Ratio
		if {[regexp "\[^0-9\]" $varReadDepth] || [regexp "\[^0-9\]" $totalReadDepth] || $totalReadDepth eq "0"} {
		    set varReadPercent "NA"
		} else {
		    set varReadPercent [format "%.0f" [expr {$varReadDepth*100.0/$totalReadDepth}]]
		}
		set g_perso($patient) "$zygosity\t$totalReadDepth\t$varReadDepth\t$varReadPercent\t$QUALphred"
		set zygous($patient) "$zygosity"
	    }
	}
    }

    # Definition of $avgVariantDepth $sdVariantDepth $countVariantDepth $avgTotalDepth $sdTotalDepth $countTotalDepth
    #################################################################################################################
    set avgTotalDepth   "-1"
    set sdTotalDepth     "-1"
    set countTotalDepth [llength $L_TotalDepth]
    if {$countTotalDepth > 1 && $L_TotalDepth ne 0 && $L_TotalDepth ne {}} {
	set MVSD [BasicStatistics $L_TotalDepth]
	set avgTotalDepth [format "%.0f" [lindex $MVSD 0]]
	set sdTotalDepth   [format "%.0f" [lindex $MVSD 2]]
    } else {
	if {$L_TotalDepth ne {}} {
	    set avgTotalDepth [format "%.0f" $L_TotalDepth]
	} 
	set sdTotalDepth   "0"
    }

    set avgVariantDepth   "-1"
    set sdVariantDepth     "-1"
    set countVariantDepth [llength $L_VariantDepth]
    if {$countVariantDepth>1 && $L_VariantDepth ne 0 && $L_VariantDepth ne {}} {
	set MVSD [BasicStatistics $L_VariantDepth]
	set avgVariantDepth [format "%.0f" [lindex $MVSD 0]]
	set sdVariantDepth   [format "%.0f" [lindex $MVSD 2]]
    } else {
	if {$L_VariantDepth ne {}} {
	    set avgVariantDepth [format "%.0f" $L_VariantDepth]
	} else {
	    set avgVariantDepth "-1"
	}
	set sdVariantDepth   "0"
    }
    #puts "$ID $avgVariantDepth $sdVariantDepth $countVariantDepth $avgTotalDepth $sdTotalDepth $countTotalDepth"


    # Definition of $barcode, $g_famBarcode($fam), $homCount, $hetCount
    #################################################################################################################
    set nbTotPatient [llength $g_allPatients]
    set barcode ""
    set nbPatientHetAveclID 0
    set nbPatientHomAveclID 0
    foreach fam [array names g_lPatientsOf] {
	set g_famBarcode($fam) ""
	foreach patient $g_lPatientsOf($fam) {
	    if {[info exists zygous($patient)]} {
		incr nbPatientAveclID
		if {$zygous($patient) eq "hom" || $zygous($patient) eq "hom?"} {
		    incr nbPatientHomAveclID
		    append g_famBarcode($fam) "2"
		    append barcode "2"
		} elseif {$zygous($patient) eq "het" || $zygous($patient) eq "het?"} {
		    incr nbPatientHetAveclID
		    append g_famBarcode($fam) "1"
		    append barcode "1"
		} else {
		    #incr nbPatientHetAveclID
		    append g_famBarcode($fam) "0"
		    append barcode "0"
		    #append g_famBarcode($fam) "*"
		    #append barcode "*"
		}
	    } else {
		append g_famBarcode($fam) "0"
		append barcode "0"
	    }
	}
    }
    set homCount "$nbPatientHomAveclID"
    set hetCount "$nbPatientHetAveclID"


    # Definition of $alleleCount and $sampleCount
    #################################################################################################################
    # alleleCount is now the allele count (the homozygous counts as 2)
    set alleleCount [expr {(2*$homCount)+$hetCount}]
    set sampleCount "$nbTotPatient"


    # Global statistics
    #################################################################################################################
    if {$alleleCount eq "0"} {
	incr g_Statistics(all,Null) 
    } elseif {$homCount eq "0"} {
	incr g_Statistics(all,Het) 
    } elseif {$hetCount eq "0"} {
    	incr g_Statistics(all,Hom) 
    } else {
	incr g_Statistics(all,Both) 
    }

    #puts "$ID $barcode $homCount $hetCount $alleleCount $sampleCount"
    #puts "$ID $avgVariantDepth $sdVariantDepth $countVariantDepth $avgTotalDepth $sdTotalDepth $countTotalDepth"
    return [list $barcode $homCount $hetCount $alleleCount $sampleCount $avgVariantDepth $sdVariantDepth $countVariantDepth $avgTotalDepth $sdTotalDepth $countTotalDepth]
}


##
## Ranking by variant for all the variants, creation of 1 output file by patient.
## No filter applied on these output files.
##
## OUTPUTS: g_VaRank(vcfDir)/"family"_"patient"_"genomeBuild"_allVariants.rankingByVar.tsv (1 by patient)
##
proc writeAllVariantsRankingByVar {} {

    global g_vcfINFOS_Supp
    global g_VaRank
    global g_ANNOTATION
    global g_PPH2
    global g_lScore
    global g_allPatients
    global g_lPatientsOf
    global g_famBarcode
    global g_deltaMES
    global g_deltaSSF
    global g_deltaNNS
    global g_perso
    global g_Statistics
    global env

    
    ## Delete output files if they already exist: 
    #############################################
    set genomeBuild ""
    if {[info exist g_VaRank(alamutHumanDB)]} {
	set genomeBuild $g_VaRank(alamutHumanDB)
    } elseif {[info exist g_VaRank(snpeffHumanDB)]} {
	set genomeBuild $g_VaRank(snpeffHumanDB)
    }
    foreach fam [array names g_lPatientsOf] {
	foreach patient $g_lPatientsOf($fam) {
	    if {$g_VaRank(SamOut) ne "all" && [lsearch -exact -nocase $g_VaRank(SamOut) $patient] eq -1} {continue}
	    set outputfile "$g_VaRank(vcfDir)/[set fam]_[set patient]_${genomeBuild}_allVariants.rankingByVar.tsv"
	    if {[file exists $outputfile]} {file delete -force $outputfile}
	}
    }
   

    ## Search column numbers in "g_ANNOTATION(#id)" to parse in g_ANNOTATION
    ########################################################################         
    # output column names calculated by VaRank:
    set L_colNamesFromVarank {variantID chr start end ref alt zygosity totalReadDepth varReadDepth varReadPercent QUALphred deltaMaxEntScorePercent deltaSSFscorePercent deltaNNSscorePercent localSpliceAnnotation PPH2pred varankVarScore annotationAnalysis avgTotalDepth sdTotalDepth countTotalDepth avgVariantDepth sdVariantDepth countVariantDepth familyBarcode barcode homCount hetCount alleleCount sampleCount alleleFrequency samVa}
    # search for indices in [split $g_ANNOTATION(#id) "\t"]:                               
    set L_colNamesFromAnnotation {}
    set L_colNamesToRemove {}
    foreach colName $g_VaRank(L_outputColHeader) {
	# no need of indices for columns calculated by VaRank
	if {[lsearch -exact $L_colNamesFromVarank $colName] ne -1} {continue}
	# SEARCH for indices in g_ANNOTATION(#id),
	# puts WARNINGS if annotation column is missing
	set iCol  "i_$colName"
	set $iCol [lsearch -regexp [split $g_ANNOTATION(#id) "\t"] "^$colName"]
	if {[set $iCol] eq -1} {
	    puts "\twarning: column name \"$colName\" not found in annotations. Removed."
	    lappend L_colNamesToRemove $colName
	    continue
	}
	lappend L_colNamesFromAnnotation $colName
    }
    # Remove the column names "L_colNamesToRemove" from $g_VaRank(L_outputColHeader)
    set new_L_outputColHeader {}
    foreach cntr $L_colNamesToRemove {
	set i [lsearch -exact $g_VaRank(L_outputColHeader) $cntr]
	set g_VaRank(L_outputColHeader) [lreplace $g_VaRank(L_outputColHeader) $i $i]
    }

    
    ## Statistics
    #############
    if {[info exists g_VaRank(snpeffDir)]} {
	set L_codingEffect $g_VaRank(codingEffect)
	set L_varLocation  "unknown"
    } else {
	# Alamut coding effect useful for statistics
	set L_codingEffect [list synonymous missense "stop gain" in-frame frameshift "start loss" "stop loss"]
	set L_varLocation  [list intron upstream "5'UTR" "3'UTR" downstream "splice site"]
    }

    
    ## Define the 3 headlines for each output ($outputLine($patient)):
    ###################################################################
    foreach fam [array names g_lPatientsOf] {
	foreach patient $g_lPatientsOf($fam) {

	  
	    if {$g_VaRank(SamOut) ne "all" && [lsearch -exact -nocase $g_VaRank(SamOut) $patient] eq -1} {continue}		    
	    set outputfile "$g_VaRank(vcfDir)/[set fam]_[set patient]_${genomeBuild}_allVariants.rankingByVar.tsv"
	    file delete -force $outputfile	
	

	    # Continue if 'patient' is not in the output samples list:
	    if {$g_VaRank(SamOut) ne "all" && [lsearch -exact -nocase $g_VaRank(SamOut) $patient] eq -1} {continue}			     
	    
	    # Compulsory header columns
	    set    outputLine($patient) "## Barcode: $g_allPatients"
	    append outputLine($patient) "\n## FamilyBarcode: $g_lPatientsOf($fam)"
	    append outputLine($patient) "\n[join $g_VaRank(L_outputColHeader) "\t"]"
	    
	    # Optional header columns
	    if {[info exists g_vcfINFOS_Supp(Header)] && [set g_vcfINFOS_Supp(Header)] ne {}} {
		set rajoutVCFinfos "[join [set g_vcfINFOS_Supp(Header)] "\t"]"
		append outputLine($patient) "\t$rajoutVCFinfos"
	    } 
	    if {[info exists g_VaRank(extann)] && $g_VaRank(extann) ne ""} {
		ExternalAnnotations
		append outputLine($patient) "\tgenes"
		foreach F [ExternalAnnotations L_Files] {
		    #puts $F 
		    #puts "Header >[ExternalAnnotations $F Header]<"
		    append outputLine($patient) "\t[ExternalAnnotations $F Header]"
		}
	    }
	    # Exomiser header columns
	    append outputLine($patient) "\tExomiser_gene_pheno_score\tHuman_pheno_evidence\tMouse_pheno_evidence\tFish_pheno_evidence"
	    
	    # End of header line
	    lappend L_outputLine($patient) $outputLine($patient)
	    set outputLine($patient) ""
	}
    }


    ## Downloading genetic variants annotated by snpeff or alamut.
    ## AlamutAnalysis = "yes"
    ## g_lScore is sorted in descending order of scores (with no more redundancy).
    ##############################################################################
    puts "...writing output files: all variants, ranking by var ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
    puts "\t...from annotated data ([llength $g_lScore] variants) ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"

    set outputSize 0

    foreach duoIDscore $g_lScore {
	set ID [lindex $duoIDscore 0]
	set varankVarScore [lindex $duoIDscore 1]
	set L [split $g_ANNOTATION($ID) "\t"]
	#puts "g_ANNOTATION($ID): $g_ANNOTATION($ID)"
	#puts "$ID - $L\n"

	set chromStartRefAltRs [db_vcfData eval {SELECT chrom, pos, ref, alt, rsID, rsValid FROM vcfSVdata WHERE ID = $ID}]
	if {$chromStartRefAltRs eq ""} {continue}
	
	# => Variation is computed and at least 1 patient has it...	
	set variantID $ID
	set chr   [lindex $chromStartRefAltRs 0]
	set start [lindex $chromStartRefAltRs 1]
	set ref   [lindex $chromStartRefAltRs 2]
	set alt   [lindex $chromStartRefAltRs 3]
	set end   [expr {$start+[string length $ref]-1}]
	#DEBUG
	#puts "startrefalt:$ID - $start - $ref - $alt"

	# Attributing the value for each annotation category 
	foreach colName $L_colNamesFromAnnotation {
	    set iCol "i_$colName"	    	    
	    if {[set $iCol] eq -1} {
		set $colName NA
	    } else {
		set $colName [lindex $L [set $iCol]]   ;# varankVarScore defined here
	    }	    
	    #DEBUG
	    #puts "iCol $iCol - [set $iCol] - [set $colName]"
	}

	# If rsID is given by VCF, we can keep the rsID and the rsValidated informations
	# (Else, "rsID and rsValidations" are already defined above with "$L_colNamesFromAnnotation")
	if {$g_VaRank(rsFromVCF) eq "yes"} {
	    set rsIdVCF [lindex $chromStartRefAltRs 4]
	    # Rs from VCF is GOOD
	    if {![isNotAnRS $rsIdVCF]} {
		set rsId  $rsIdVCF
		set rsValidations [lindex $chromStartRefAltRs 5]
	    } else {
		# rsID from vcf is not a good rsID testing from Alamut
		set rsId "NA"; set rsValidations "NA"
	    }
	} 
	
	# Definition of the following values:
	# barcode homCount hetCount alleleCount sampleCount avgVariantDepth sdVariantDepth countVariantDepth avgTotalDepth sdTotalDepth countTotalDepth
	# alleleFrequency
        set infos ""
        set infos [findBarcodesAndStatFor $ID]
	lassign $infos barcode homCount hetCount alleleCount sampleCount avgVariantDepth sdVariantDepth countVariantDepth avgTotalDepth sdTotalDepth countTotalDepth

	set alleleFrequency [format "%.4f" [expr {$alleleCount*1.0/2/$sampleCount}]]

	# Definition of SamVa
	if {[regexp -nocase "samVa" $g_VaRank(L_outputColHeader)]} {
	    set nb    0
	    set samVa ""
	    foreach b [split $barcode ""] s $g_allPatients {
		if {$b ne "0"} {append samVa "${b}_$s "; incr nb}
		if {$nb eq 10} {append samVa "..."; break}
	    }
	}

	# Converting empty "count" or "gnomadReadDepth" into -1 allow to filter for below 0 including missing values (For rsMAFCount, gnomadHomCount_*, ...)
	# Converting empty "frequencies" into -1 allow to filter for below 0.01 including missing values (For rsMAF, 1000g_AF, exacAltFreq_sas, allFrequency...)
	foreach colName $g_VaRank(L_outputColHeader) {
	    if {[regexp "Count|gnomadReadDepth" $colName]} {
		if {[set $colName] eq "NA"} {set $colName "-1";continue}
	    }

	    if {![regexp "AF$|Freq|FREQ" $colName]} {continue}
	    if {[set $colName] eq "NA"} {set $colName "-1";continue}	    
	    # Change also metrics to "." to ","
	    if {[set g_VaRank(metrics)] eq "fr"} {regsub {\.} [set $colName] "," $colName}
	}

	# Splice effect annotation (only with ALAMUT annotation)
	if {[info exists env(ALAMUT)]} {
	    set localSpliceAnnotation "NA"
	    if {$nearestSSType eq "5'" && [regexp "^1$|^2$" $distNearestSS]} {
		set localSpliceAnnotation "essential splice donor"
	    } elseif {$nearestSSType eq "3'" && [regexp "^-1$|^-2$" $distNearestSS]} {
		set localSpliceAnnotation "essential splice acceptor"
	    } elseif {$nearestSSType eq "5'" && $distNearestSS >= -3  && $distNearestSS <= 6} {
		set localSpliceAnnotation "close splice donor"
	    } elseif {$nearestSSType eq "3'" && $distNearestSS >= -12 && $distNearestSS <= 2} {
		set localSpliceAnnotation "close splice acceptor"
	    }
	    if {[catch {set deltaMaxEntScorePercent "$g_deltaMES($ID,$transcript)"} Message]} {set deltaMaxEntScorePercent "NA"}
	    if {[catch {set deltaSSFscorePercent "$g_deltaSSF($ID,$transcript)"   } Message]} {set deltaSSFscorePercent "NA"}
	    if {[catch {set deltaNNSscorePercent "$g_deltaNNS($ID,$transcript)"   } Message]} {set deltaNNSscorePercent "NA"}
	    #DEBUG
	    #puts "$ID $transcript : $g_deltaMES($ID,$transcript)\t$varMaxEntScore\t$g_deltaSSF($ID,$transcript)"
	}

	# PPH2
	set PPH2pred $g_PPH2($ID)
	
	# Collecting data for the global statistics
	if {![info exists g_Statistics(All)]} {
	    foreach e [concat $L_codingEffect $L_varLocation unknown] {
		set g_Statistics(All,$e)   0
	    }
	    set g_Statistics(All) 0
	} 
	incr g_Statistics(All) 

	# Counting for coding effect
	if {$codingEffect ne "NA" && $codingEffect ne ""} {
	    incr g_Statistics(All,$codingEffect)  
	} elseif {$varLocation ne "NA" && $varLocation ne ""} {
	    incr g_Statistics(All,$varLocation)  
	} else {
	    incr g_Statistics(All,unknown)
	} 
	
	set annotationAnalysis "yes"
	set barcode "'$barcode'"

	
        # Searching for annotation specific of each patient #################################################
	# Then append all annotations (non specific + specific) to the "outputLine(patient)" variable #######
	
	# For the above given ID, tmp definition of vcfINFOS_Supp_Tmp($patient)
	if {[set g_VaRank(vcfInfo)]=="yes"} {
	    catch {unset vcfINFOS_Supp_Tmp}
	    foreach fam [array names g_lPatientsOf] {
		foreach patient $g_lPatientsOf($fam) {
		    set theInfo [db_vcfData eval {SELECT INFO FROM vcfSampleData AS vs \
						   INNER JOIN sampleName AS sn ON vs.sampleName_id = sn.sampleName_id \
						   WHERE ID = $ID AND sampleName = $patient}]
		    if {$theInfo eq ""} {continue}
		    set theInfo [lindex $theInfo 0]; # SQLite return a "1 element" list (that we don't want)
		    if {$g_VaRank(vcfFields)=="all"} {
			set vcfINFOS_Supp_Tmp($patient) [split $theInfo ";"]
		    } else {
			foreach duoInfos_VCF [split $theInfo ";"] {
			    set duo_Header_Infos_VCF [split $duoInfos_VCF "="]
			    set Header [lindex $duo_Header_Infos_VCF 0]
			    if {[lsearch -exact $g_VaRank(vcfFields) $Header]!=-1} {
				lappend vcfINFOS_Supp_Tmp($patient) $duoInfos_VCF
			    }
			}
		    }
		}
	    }
	}

	foreach fam [array names g_lPatientsOf] {
	    foreach patient $g_lPatientsOf($fam) {

		if {$g_VaRank(SamOut) ne "all" && [lsearch -exact -nocase $g_VaRank(SamOut) $patient] eq -1} {continue}
		
		set theInfo [db_vcfData eval {SELECT INFO FROM vcfSampleData AS vs \
						  INNER JOIN sampleName AS sn ON vs.sampleName_id = sn.sampleName_id \
						  WHERE ID = $ID AND sampleName = $patient}]
		if {$theInfo eq ""} {continue}

		set familyBarcode "'$g_famBarcode($fam)'"

		set zygosity               [lindex $g_perso($patient) 0]
		set totalReadDepth         [lindex $g_perso($patient) 1]
		set varReadDepth           [lindex $g_perso($patient) 2]
		set varReadPercent         [lindex $g_perso($patient) 3]
		set QUALphred              [lindex $g_perso($patient) 4]

		# Collecting data for the statistics per patient	   
		if {![info exists g_Statistics($patient)]} {
		    foreach e [concat $L_codingEffect $L_varLocation unknown] {
			set g_Statistics($patient,$e)      0
			set g_Statistics($patient,$e,Hom)  0
			set g_Statistics($patient,$e,Het)  0
			set g_Statistics($patient,$e,Null) 0
		    }
		    set g_Statistics($patient) 0
		} 
		incr g_Statistics($patient) 

		set HomHet [lindex [split $g_perso($patient) "\t"] 0]

		# Counting for coding effect
		if {$codingEffect ne "NA" && $codingEffect ne ""} {
		    incr g_Statistics($patient,$codingEffect)  
		    if {[regexp "hom" $HomHet]} {
			incr g_Statistics($patient,$codingEffect,Hom)  
		    } elseif {[regexp "het" $HomHet]} {
			incr g_Statistics($patient,$codingEffect,Het)  
		    } else {
			incr g_Statistics($patient,$codingEffect,Null)
		    }
		} elseif {$varLocation ne "NA" && $varLocation ne ""} {
		    incr g_Statistics($patient,$varLocation)  
		    if {[regexp "hom" $HomHet]} {
			incr g_Statistics($patient,$varLocation,Hom)  
		    } elseif {[regexp "het" $HomHet]} {
			incr g_Statistics($patient,$varLocation,Het)  
		    } else {
			incr g_Statistics($patient,$varLocation,Null)
		    }
		} else {
		    incr g_Statistics($patient,unknown)
		    if {[regexp "hom" $HomHet]} {
			incr g_Statistics($patient,unknown,Hom)  
		    } elseif {[regexp "het" $HomHet]} {
			incr g_Statistics($patient,unknown,Het)  
		    } else {
			incr g_Statistics($patient,unknown,Null)
		    }
		}
		

		# append values for compulsory columns (the ones from the configfile)
		set L_rankText {} 
		foreach colName $g_VaRank(L_outputColHeader) {
		    lappend L_rankText [set $colName]
		}
		append outputLine($patient) [join $L_rankText "\t"] ;# add varankVarScore here

		# Append values for optional columns (the ones from the VCF files):
		###################################################################
		# -> from the INFO column of the VCF:
		set l_Infos_VCF  {}
		set l_headers_ID {}
		set l_data_ID    {}		
		set rajoutVCFinfos {}
		if {[set g_VaRank(vcfInfo)] eq "yes"} {
		    set l_Infos_VCF ""
		    if {[info exists vcfINFOS_Supp_Tmp($patient)]} {set l_Infos_VCF [set vcfINFOS_Supp_Tmp($patient)]}
		    foreach infos_VCF $l_Infos_VCF {
			set l_Header_Infos_VCF ""
			set Header ""
			set Data   ""			
			if {[regexp "=" $infos_VCF]} {
			    set l_Header_Infos_VCF [split $infos_VCF "="]			    
			    set Header [lindex $l_Header_Infos_VCF 0]
			    set Data   [lindex $l_Header_Infos_VCF 1]
			} else {
			    set l_Header_Infos_VCF $infos_VCF			    
			    set Header $l_Header_Infos_VCF
			    set Data   $l_Header_Infos_VCF
			}			
			lappend l_headers_ID $Header
			lappend l_data_ID    $Data
		    }
		    foreach NewHeaders_VCF [set g_vcfINFOS_Supp(Header)] {
			set  i_header [lsearch -exact  $l_headers_ID $NewHeaders_VCF]
			if {$i_header eq -1} {lappend rajoutVCFinfos "NA"} else {lappend rajoutVCFinfos [lindex $l_data_ID $i_header]}
		    }
		    set rajoutVCFinfos [join $rajoutVCFinfos "\t"]
		}		
		if {$rajoutVCFinfos ne {}} {
		    append outputLine($patient) "\t$rajoutVCFinfos"
		} 

		# Adding external user annotations at the end of each output files line
		#######################################################################
		if {[info exists g_VaRank(extann)] && $g_VaRank(extann) ne ""} {
		    ## Adding external annotations only on the first gene (decided the 19/07/2018 with Jean)		    
		    set g [lindex [split $gene "/"] 0] 
		    append outputLine($patient) "\t$g"
		    foreach F [ExternalAnnotations L_Files] {
			if {[ExternalAnnotations $F $g] ne ""} {
			    append outputLine($patient) "\t[ExternalAnnotations $F $g]"
			} else {			    
			    set NbHeader [llength [split [ExternalAnnotations $F Header] "\t"]]
			    append outputLine($patient) "\t[join [lrepeat $NbHeader ""] "\t"]"
			}
		    }
		}

		# Adding exomiser annotation
		############################
		## Adding exomiser annotations only on the first gene
		set g [lindex [split $gene "/"] 0] 
		append outputLine($patient) "\t[ExomiserAnnotation "$patient" "$g" "all"]" 
		
		# End of annotation line
		########################
		lappend L_outputLine($patient) $outputLine($patient)
		incr outputSize [string length $outputLine($patient)]
		set outputLine($patient) ""

		# Intermediate writing to save RAM
		##################################
		if {$outputSize > 8000000} {
		    set outputSize 0
		    foreach famHere [array names g_lPatientsOf] {
			foreach patientHere $g_lPatientsOf($famHere) {			    
			    if {$g_VaRank(SamOut) ne "all" && [lsearch -exact -nocase $g_VaRank(SamOut) $patientHere] eq -1} {continue}		    
			    set outputfileHere "$g_VaRank(vcfDir)/[set famHere]_[set patientHere]_${genomeBuild}_allVariants.rankingByVar.tsv"
			    if {[info exists L_outputLine($patientHere)]} {
				WriteTextInFile [join "$L_outputLine($patientHere)" "\n"] $outputfileHere
				unset L_outputLine($patientHere)
			    }
			}
		    }
		}
	    }
	}
    }

    ## Downloading genetic variants not analysed by alamut.
    ## AlamutAnalysis = "no"
    #######################################################

    # Calculate the number of unannotated variants
    set n_notannotated 0
    set L_IDs [db_vcfData eval {SELECT ID FROM vcfSVdata}]

    foreach ID $L_IDs {
	if {[info exists g_ANNOTATION($ID)]} {continue}
        # => Variation is in VCF but not annotated
	incr n_notannotated
    }
    puts "\t...from data not annotated ($n_notannotated variants) ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"

    set outputSize 0

    # Needed if there are only variants not analysed by Alamut
    set i_gene [lsearch -exact [split $g_ANNOTATION(#id) "\t"] "gene"]
    foreach ID $L_IDs {
	if {[info exists g_ANNOTATION($ID)]} {continue}
        # => Variation is in VCF but not annotated

	set chromStartRefAltRs [db_vcfData eval {SELECT chrom, pos, ref, alt, rsID, rsValid FROM vcfSVdata WHERE ID = $ID}]
	if {$chromStartRefAltRs eq ""} {continue}
	
	# => Variation is computed and at least 1 patient has it...	
	set variantID $ID
	set chr   [lindex $chromStartRefAltRs 0]
	set start [lindex $chromStartRefAltRs 1]
	set ref   [lindex $chromStartRefAltRs 2]
	set alt   [lindex $chromStartRefAltRs 3]
	set end   [expr {$start+[string length $ref]-1}]

	#DEBUG
	#puts "chrstartrefalt:$ID - $chr - $start - $ref - $alt"
	
	# Attributing the "NA" or "-1" value for each annotation category 
	if {![info exists attribution1]} {  ;# Values are either empty or equal to "NA"/"-1". No need to redo this values attribution several times.
	    set attribution1 "done"
	    set deltaMaxEntScorePercent "NA"
	    set deltaSSFscorePercent "NA"
	    set deltaNNSscorePercent "NA"
	    foreach colName $L_colNamesFromAnnotation {
		if {[regexp "AF$|Freq|Count|gnomadReadDepth" $colName]} {
		    set $colName "-1"
		} else {
		    set $colName "NA"; # varankVarScore defined here with "NA". Need to be set to "0" below.
		}	    
	    }
	    set varankVarScore 0
	}
	
	## Attribution of a gene name for upstream deletion (from VCF 4.2, ALT='*')
	if {[regexp "_\\*$" $ID]} {
	    set gene    [geneForStarID $ID $i_gene]
	    set varType "upstream deletion"
	} else {
	    set gene    "NA"
	    set varType "NA"
	}

	# If rsID is given by VCF, we can keep the rsID and the rsValidated informations
	if {[set g_VaRank(rsFromVCF)] eq "yes"} {
	    set rsIdVCF   [lindex $chromStartRefAltRs 4]
	    #rs from VCF is GOOD
	    if {![isNotAnRS $rsIdVCF]} {
		set rsId $rsIdVCF
		set rsValidations [lindex $chromStartRefAltRs 5]
	    } else {
		set rsId "NA"; set rsValidations "NA"
	    }
	}

	# Definition of the following values:
	# barcode homCount hetCount alleleCount sampleCount avgVariantDepth sdVariantDepth countVariantDepth avgTotalDepth sdTotalDepth countTotalDepth
	# alleleFrequency
	set infos   ""
	set infos   [findBarcodesAndStatFor $ID]
	lassign $infos barcode homCount hetCount alleleCount sampleCount avgVariantDepth sdVariantDepth countVariantDepth avgTotalDepth sdTotalDepth countTotalDepth

	set alleleFrequency [format "%.4f" [expr {$alleleCount*1.0/2/$sampleCount}]]
	# Change metrics from "." to ","
	if {[set g_VaRank(metrics)] eq "fr"} {regsub {\.} $alleleFrequency "," alleleFrequency}

	# Definition of SamVa
	if {[regexp -nocase "samVa" $g_VaRank(L_outputColHeader)]} {
	    set nb 0
	    set samVa ""
	    foreach b [split $barcode ""] s $g_allPatients {
		if {$b ne "0"} {append samVa "${b}_$s "; incr nb}
		if {$nb eq 10} {append samVa "..."; break}
	    }
	}

	# Collecting data for the global statistics
	if {![info exists g_Statistics(All)]} {
	    foreach e [concat $L_codingEffect $L_varLocation unknown] {
		set g_Statistics(All,$e)   0
	    }
	    set g_Statistics(All) 0
	} 
	incr g_Statistics(All) 
	incr g_Statistics(All,unknown)

	set annotationAnalysis "no"
	set barcode "'$barcode'"


        # Searching for annotation specific of each patient #################################################
	# Then append all annotations (non specific + specific) to the "outputLine(patient)" variable #######
	
	# For the above given ID, tmp definition of vcfINFOS_Supp_Tmp($patient)
	if {[set g_VaRank(vcfInfo)]=="yes"} {
	    catch {unset vcfINFOS_Supp_Tmp}
	    foreach fam [array names g_lPatientsOf] {
		foreach patient $g_lPatientsOf($fam) {
		    set theInfo [db_vcfData eval {SELECT INFO FROM vcfSampleData AS vs \
						   INNER JOIN sampleName AS sn ON vs.sampleName_id = sn.sampleName_id \
						   WHERE ID = $ID AND sampleName = $patient}]
		    if {$theInfo eq ""} {continue}
		    set theInfo [lindex $theInfo 0]; # SQLite return a "1 element" list (that we don't want)
		    if {$g_VaRank(vcfFields)=="all"} {
			set vcfINFOS_Supp_Tmp($patient) [split $theInfo ";"]
		    } else {
			foreach duoInfos_VCF [split $theInfo ";"] {
			    set duo_Header_Infos_VCF [split $duoInfos_VCF "="]
			    set Header [lindex $duo_Header_Infos_VCF 0]
			    if {[lsearch -exact $g_VaRank(vcfFields) $Header]!=-1} {
				lappend vcfINFOS_Supp_Tmp($patient) $duoInfos_VCF
			    }
			}
		    }
		}		             
	    }
	}

	catch {unset attributionForThisID}
	foreach fam [array names g_lPatientsOf] {
	    foreach patient $g_lPatientsOf($fam) {
		
		if {$g_VaRank(SamOut) ne "all" && [lsearch -exact -nocase $g_VaRank(SamOut) $patient] eq -1} {continue}		

		set theInfo [db_vcfData eval {SELECT INFO FROM vcfSampleData AS vs \
						  INNER JOIN sampleName AS sn ON vs.sampleName_id = sn.sampleName_id \
						  WHERE ID = $ID AND sampleName = $patient}]
		if {$theInfo eq ""} {continue}

		set familyBarcode "'$g_famBarcode($fam)'"
		
		set zygosity               [lindex $g_perso($patient) 0]
		set totalReadDepth         [lindex $g_perso($patient) 1]
		set varReadDepth           [lindex $g_perso($patient) 2]
		set varReadPercent         [lindex $g_perso($patient) 3]
		set QUALphred              [lindex $g_perso($patient) 4]

		# Collecting data for the statistics per patient
		if {![info exists g_Statistics($patient)]} {
		    foreach e [concat $L_codingEffect $L_varLocation unknown] {
			set g_Statistics($patient,$e)      0
			set g_Statistics($patient,$e,Hom)  0
			set g_Statistics($patient,$e,Het)  0
			set g_Statistics($patient,$e,Null) 0
		    }
		    set g_Statistics($patient) 0
		} 
		incr g_Statistics($patient) 
		
		set HomHet [lindex [split $g_perso($patient) "\t"] 0]

		# $codingEffect eq "NA" && $varLocation eq "NA"
		incr g_Statistics($patient,unknown)
		if {[regexp "hom" $HomHet]} {
		    incr g_Statistics($patient,unknown,Hom)
		} elseif {[regexp "het" $HomHet]} {
		    incr g_Statistics($patient,unknown,Het)
		} else {
		    incr g_Statistics($patient,unknown,Null)
		}
		
		# append values for compulsory columns (the ones from the configfile)
		set L_rankText {}
		foreach colName $g_VaRank(L_outputColHeader) {
		    lappend L_rankText [set $colName]
		}
		append outputLine($patient) [join $L_rankText "\t"]
		#puts "!! [join $L_rankText "\t"]"		
	
		# Append values for optional columns (the ones from the VCF files):
		###################################################################
		if {![info exists attributionForThisID]} {
		    set attributionForThisID "done"
		    
		    set rajoutVCFinfos {}		

		    # -> from the INFO column of the VCF:
		    set l_Infos_VCF  {}
		    set l_headers_ID {}
		    set l_data_ID    {}		
		    if {[set g_VaRank(vcfInfo)] eq "yes"} {
			set l_Infos_VCF ""
		        if {[info exists vcfINFOS_Supp_Tmp($patient)]} {set l_Infos_VCF [set vcfINFOS_Supp_Tmp($patient)]}
			foreach infos_VCF $l_Infos_VCF {
			    set l_Header_Infos_VCF ""
			    set Header ""
			    set Data   ""			
			    if {[regexp "=" $infos_VCF]} {
				set l_Header_Infos_VCF [split $infos_VCF "="]			    
				set Header [lindex $l_Header_Infos_VCF 0]
				set Data   [lindex $l_Header_Infos_VCF 1]
			    } else {
				set l_Header_Infos_VCF $infos_VCF			    
				set Header $l_Header_Infos_VCF
				set Data   $l_Header_Infos_VCF
			    }			
			    lappend l_headers_ID $Header
			    lappend l_data_ID    $Data
			}		    
			foreach NewHeaders_VCF [set g_vcfINFOS_Supp(Header)] {
			    set  i_header [lsearch -exact  $l_headers_ID $NewHeaders_VCF]
			    if {$i_header eq -1} {lappend rajoutVCFinfos "NA"} else {lappend rajoutVCFinfos [lindex $l_data_ID $i_header]}
			}
			set rajoutVCFinfos [join $rajoutVCFinfos "\t"]
		    }
		}
		if {$rajoutVCFinfos ne {}} {
		    append outputLine($patient) "\t$rajoutVCFinfos"
		} 
		    
		
		## Adding external user annotations at the end of the output files lines
		########################################################################
		if {![info exists attribution2]} {  ;# Values are always the same. No need to redo this values attribution several times.
		    set attribution2 "done"
		    
		    set ExtAnnPat ""
		    if {[info exists g_VaRank(extann)] && $g_VaRank(extann) ne ""} {
			# We have: $gene eq "NA"		    
			append ExtAnnPat "\t"
			foreach F [ExternalAnnotations L_Files] {
			    #puts $F 
			    set NbHeader [llength [split [ExternalAnnotations $F Header] "\t"]]
			    append ExtAnnPat "\t[join [lrepeat $NbHeader ""] "\t"]"
			    #puts "Header >[ExternalAnnotations $F Header]< ---> $NbHeader"
			    #puts "\t[join [lrepeat $NbHeader NA] "\t"]"
			}		    
			
			# Adding exomiser annotation
			############################
			append ExtAnnPat "\t-1.0\tNA\tNA\tNA"
		    }
		}
		append outputLine($patient) "$ExtAnnPat"
		
		# End of annotation line
		########################
		lappend L_outputLine($patient) $outputLine($patient)
		incr outputSize [string length $outputLine($patient)]
		set outputLine($patient) ""
		
		# Intermediate writing to save RAM
		##################################
		if {$outputSize > 8000000} {
		    set outputSize 0
		    foreach famHere [array names g_lPatientsOf] {
			foreach patientHere $g_lPatientsOf($famHere) {			    
			    if {$g_VaRank(SamOut) ne "all" && [lsearch -exact -nocase $g_VaRank(SamOut) $patientHere] eq -1} {continue}		    
			    set outputfileHere "$g_VaRank(vcfDir)/[set famHere]_[set patientHere]_${genomeBuild}_allVariants.rankingByVar.tsv"
			    if {[info exists L_outputLine($patientHere)]} {
				WriteTextInFile [join "$L_outputLine($patientHere)" "\n"] $outputfileHere
				unset L_outputLine($patientHere)
			    }
			}
		    }
		}
            }
	}
    }	
    
    ## Writing "*_allVariants.rankingByVar" output files
    ####################################################
    foreach fam [array names g_lPatientsOf] {
	foreach patient $g_lPatientsOf($fam) {

	    if {$g_VaRank(SamOut) ne "all" && [lsearch -exact -nocase $g_VaRank(SamOut) $patient] eq -1} {continue}

	    set outputfile "$g_VaRank(vcfDir)/[set fam]_[set patient]_${genomeBuild}_allVariants.rankingByVar.tsv"
	    if {[info exists L_outputLine($patient)]} {
		WriteTextInFile [join "$L_outputLine($patient)" "\n"] $outputfile
		unset L_outputLine($patient)
	    }
	}
    }

    ## NOT BY DEFAULT!!
    ## Replacing sequences with length > $g_VaRank(maxSizeOfSeq) bp.
    ## Replace long sequences into the ID and from the ref and the alt:
    ## ID = chr_pos_'refLength'bp_'altLength'bp
    ## ref = 'refLength'bp
    ## alt = 'altLength'bp
    ##
    ## Remark: If sequence variations are too long, they can't be open in excel => bug. 
    ## Excel: maximum of character by cell = 32 767 char
    if {[info exists g_VaRank(maxSizeOfSeq)]} {
	puts "\t...updating VariantID nomenclature for indels > $g_VaRank(maxSizeOfSeq) bp ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
	
	foreach fam [array names g_lPatientsOf] {
	    foreach patient $g_lPatientsOf($fam) {
		
		if {$g_VaRank(SamOut) ne "all" && [lsearch -exact -nocase $g_VaRank(SamOut) $patient] eq -1} {continue}

		set outputfile "$g_VaRank(vcfDir)/[set fam]_[set patient]_${genomeBuild}_allVariants.rankingByVar.tsv"

		set newText ""
		foreach L [LinesFromFile $outputfile] {
		    set Ls [split $L "\t"]
		    if {[regexp "^#" $L]} {append newText "$L\n"; continue}
		    if {[regexp "^variantID" $L]} {
			set i_variantID [lsearch -exact $Ls "variantID"]; if {$i_variantID eq -1} {puts "variantID column not found for VaRank output - Exit"; exit}
			set i_ref       [lsearch -exact $Ls "ref"];       if {$i_ref       eq -1} {puts "ref column not found for VaRank output - Exit"; exit}
			set i_alt       [lsearch -exact $Ls "alt"];       if {$i_alt       eq -1} {puts "alt column not found for VaRank output - Exit"; exit}
			append newText "$L\n"
			continue
		    }
		    set variantID [lindex $Ls $i_variantID]
		    set ref       [lindex $Ls $i_ref] 
		    set alt       [lindex $Ls $i_alt]
		    set lengthRef [string length $ref] 
		    set lengthAlt [string length $alt]
		    
		    if {$lengthRef > $g_VaRank(maxSizeOfSeq)} {
			set ref "${lengthRef}bp"
			set Ls [lreplace $Ls $i_ref $i_ref $ref]
			
			## chrom = 1-22, X,Y, M, MT
			if {![regexp "(\[a-zA-Z0-9\]+)_(\[0-9\]+)_" $variantID match chr pos]} {
			    puts "WARNING: no chromosome matching into $variantID"
			    append newText "[join $Ls "\t"]\n"
			    continue
			}
			if {$lengthAlt > $g_VaRank(maxSizeOfSeq)} {			 
			    set alt "${lengthAlt}bp"
			    set Ls [lreplace $Ls $i_alt $i_alt $alt]
			    set newVariantID "${chr}_${pos}_${lengthRef}bp_${lengthAlt}bp"
			    
			} else {
			    set newVariantID "${chr}_${pos}_${lengthRef}bp_$alt"
			}
			set i 0
			while {[info exists correspondance($newVariantID)] && $correspondance($newVariantID) ne $variantID} {
			    incr i
			    set newVariantID "$newVariantID.$i"			    
			}
			set correspondance($newVariantID) $variantID
			set Ls [lreplace $Ls $i_variantID $i_variantID $newVariantID]
			## For the writing of output files: all variants, ranking by gene:
		    } elseif {$lengthAlt > $g_VaRank(maxSizeOfSeq)} {
			set alt "${lengthAlt}bp"
			set Ls [lreplace $Ls $i_alt $i_alt $alt]
			
			if {![regexp "(\[a-zA-Z0-9\]+)_(\[0-9\]+)_" $variantID match chr pos]} {
			    puts "WARNING: no chromosome matching into $variantID"
			    append newText "[join $Ls "\t"]\n"
			    continue
			}
			set newVariantID "${chr}_${pos}_${lengthRef}bp_$alt"
			set i 0
			while {[info exists correspondance($newVariantID)] && $correspondance($newVariantID) ne $variantID} {
			    incr i
			    set newVariantID "$newVariantID.$i"	
			}
			set correspondance($newVariantID) $variantID
			set Ls [lreplace $Ls $i_variantID $i_variantID $newVariantID]
			## For the writing of output files: all variants, ranking by gene:
		    }
		    append newText "[join $Ls "\t"]\n"
		}
		ReplaceTextInFile $newText $outputfile
	    }
	}
    }

    return
}


##
## Ranking by gene for all variants, creation of 1 output file by patient.
## No filter applied on these output files.
##
## OUTPUTS: g_VaRank(vcfDir)/"family"_"patient"_allVariants.rankingByGene.tsv (1 by patient)
##
proc writeAllVariantsRankingByGene {} {

    global g_VaRank
    global g_lPatientsOf
    global g_allPatients

    puts "...writing output files: all variants, ranking by gene ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"

    set genomeBuild ""
    if {[info exist g_VaRank(alamutHumanDB)]} {
        set genomeBuild $g_VaRank(alamutHumanDB)
    } elseif {[info exist g_VaRank(snpeffHumanDB)]} {
        set genomeBuild $g_VaRank(snpeffHumanDB)
    }
    foreach fam [array names g_lPatientsOf] {
	foreach patient $g_lPatientsOf($fam) {
	    #puts "\t...$patient ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
	    
	    ## Delete output file if it already exists
	    set outputfile "$g_VaRank(vcfDir)/[set fam]_[set patient]_${genomeBuild}_allVariants.rankingByGene.tsv"
	    if {[file exists $outputfile]} {file delete -force $outputfile} 

	    ## Check if this output patient file is required in $g_VaRank(SamOut)
	    if {$g_VaRank(SamOut) ne "all" && [lsearch -exact -nocase $g_VaRank(SamOut) $patient] eq -1} {continue}

	    ## Searches for all htz and hom SNV, for each patient and for each gene:
	    ## --> creation of lVariants($gene) variables
	    ## --> creation of lGenes and lID variables
	    catch {unset lVariants}
	    set lGenes {}
	    set lID    {}
	    set rankingByVarFile "$g_VaRank(vcfDir)/[set fam]_[set patient]_${genomeBuild}_allVariants.rankingByVar.tsv"
	    if {![file exists $rankingByVarFile]} {
		puts "WARNING: $rankingByVarFile doesn't exist."
		continue
	    }
	    #puts "\t\t...searching for all variants ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
	    catch {unset Tab}
	    foreach L [LinesFromFile $rankingByVarFile] {
		if {$L eq ""} {continue}
		if {[regexp "^## Barcode: (.*)" $L match allExomes]} {set textAllEx "$L"; continue}
		if {[regexp "^##" $L]} {continue}
		if {[regexp "^variantID" $L]} {
		    set HeaderText "$L"
		    set L [split $L "\t"]
		    set i_gene     [lsearch -regexp $L "gene$"];            if {$i_gene  eq -1} {puts "gene: column number not found - Exit"; exit}
		    set i_score    [lsearch -regexp $L "^varankVarScore$"]; if {$i_score eq -1} {puts "varankVarScore: column number not found - Exit"; exit}
		    set i_zygosity [lsearch -regexp $L "^zygosity$"];       if {$i_zygosity eq -1} {puts "zygosity: column number not found - Exit"; exit}
		    continue
		}
		set L [split $L "\t"]
		set id   [lindex $L 0]  
		set gene [lindex $L $i_gene]
		if {$gene eq "NA"} {continue}
		set zygosity [lindex $L $i_zygosity]
		## In case of value like: gene = "MC1R/TUBB3", each gene have to be treated separately.
		## (Else, bug to assembly variants from gene = "MC1R/TUBB3" with variants from gene = "MC1R")
		## => but, by this way, we introduce redondancy in "*allVariants.rankingByGene.tsv" files.
		foreach gene [split $gene "/"] {
		    if {![info exists Tab($gene)]} {set Tab($gene) 1;lappend lGenes $gene}
		    if {![info exists Tab($id)]}   {set Tab($id)   1;lappend lID    $id}

		    #if {[lsearch -exact $lID "$id"]      eq -1} {lappend lID $id}
		    #if {[lsearch -exact $lGenes "$gene"] eq -1} {lappend lGenes $gene}
		    set score [lindex $L $i_score] 
		    
		    lappend lVariants($gene) "$id $score $zygosity"
		}
	    }

	    ## Calculate the score for each gene:
	    ## - If the most deleterious variant is hom ("scoreMDHom"), gene score = "scoreMDhom" x 2
	    ## - If the most deleterious variant is het ("scoreMDHet"): 
	    ##		- if there is another variant in the same gene, gene score = "scoreMDHet" + "following score"
	    ##		- if there isn't another variant in the same gene, gene score = "scoreMDHet" x 2
	    ##
	    ## --> creation of the variable lVariantsRankingByGene
	    catch {unset lVariantsRankingByGene}
	    #puts "\t\t...scoring of each gene  ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
	    foreach gene $lGenes {
		if {![info exists lVariants($gene)]} {continue}
		set lVariants($gene) [lsort -command DescendingSortOnElement1 $lVariants($gene)]
		set liste {}
		foreach duoIDScore $lVariants($gene) {
		    set id [lindex $duoIDScore 0]
		    lappend liste $id
		}
		set maxScore1 [lindex [lindex $lVariants($gene) 0] 1]
		if {[llength $lVariants($gene)] eq 1} {
		    set BestDuoScore "[expr {$maxScore1*2}]"
		} else {
		    set bestID [lindex [lindex $lVariants($gene) 0] 0]
		    set zygosityOfTheBesTID [lindex [lindex $lVariants($gene) 0] 2]
		    if {[regexp -nocase "hom" $zygosityOfTheBesTID]} {
			set BestDuoScore "[expr {$maxScore1*2}]"
		    } else {
			set maxScore2 [lindex [lindex $lVariants($gene) 1] 1]
			set BestDuoScore "[expr {$maxScore1+$maxScore2}]"
		    }
		}
		lappend lVariantsRankingByGene "{$liste} $BestDuoScore"
	    }
	    if {[info exists lVariantsRankingByGene]} {
		set lVariantsRankingByGene [lsort -command DescendingSortOnElement1 $lVariantsRankingByGene]
	    } else {
		continue
	    }
	    
	    ## Downloading each ranking file line.
	    ## set rankingByVarFile "$g_VaRank(vcfDir)/[set fam]_[set patient]_${genomeBuild}_allVariants.rankingByVar.tsv"
	    foreach L [LinesFromFile $rankingByVarFile] {
		set id [lindex $L 0]
		if {[info exists Tab($id)]} {set ligne($patient,$id) $L}
	    }
	    
	    ## Writing of the patient output
	    set outputfile "$g_VaRank(vcfDir)/[set fam]_[set patient]_${genomeBuild}_allVariants.rankingByGene.tsv"
	    ReplaceTextInFile "$textAllEx\n## FamilyBarcode: $g_lPatientsOf($fam)\n$HeaderText" $outputfile
	    #puts "\t\t...writing \"$outputfile\" output files ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"	    
	    set i 0
	    set L_Lines {}
	    foreach infoForOneGene $lVariantsRankingByGene {
		foreach id [lindex $infoForOneGene 0] {
		    lappend L_Lines $ligne($patient,$id)
		    incr i
		    if {$i > 10000} {
			WriteTextInFile [join $L_Lines "\n"] $outputfile
			set i 0
			set L_Lines {}
		    }
		}
	    }
	    if {$L_Lines ne ""} {
		WriteTextInFile [join $L_Lines "\n"] $outputfile
	    }
	    
	    # Memory safe!!!
	    array unset ligne "*"
	}
    }

    return
}

