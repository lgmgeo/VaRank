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

proc AssignID {chrom pos ref alt} {

    global g_VaRank

    #VCF describing variation
    #
    #SNP VCF record
    #20     3 .         C      T    .   PASS  DP=100
    #
    #Insertion VCF record
    #20     3 .         C      CTAG    .   PASS  DP=100
    #This is a insertion since the reference base C is being replaced by C [the reference base] plus three insertion bases TAG. Again there are only two alleles so I have the two following segregating haplotypes:
    
    #Deletion VCF record
    #20     2 .         TCG      T    .   PASS  DP=100
    #This is a deletion of two reference bases since the reference allele TCG is being replaced by just the T [the reference base]. Again there are only two alleles so I have the two following segregating haplotypes:
    
    #Mixed VCF record for a microsatellite
    #20     2 .         TCGCG      TCG,TCGCGCG    .   PASS  DP=100
    #This is a mixed type record containing a 2 base insertion and a 2 base deletion. There are three segregating alleles so I have the three following haplotypes:
    #Ref: a t c g c g - - a // C is the reference base
    #Ref: a t c g - - - - a // following the C base is a deletion  of 2 bases
    #Ref: a t c g c g c g a // following the C base is a insertion of 2 bases
    
    #test sur la taille si alt > ref = insertion
    #test sur la taille si alt < ref = deletion


    ## Save for each VaRank ID the VCF coresponding information
    set saveIDfile "$g_VaRank(vcfDir)/VCF_Coordinates_Conversion.tsv"
    set nativepos $pos
    set nativeref $ref
    set nativealt $alt
   
    if {![regexp "\\-" $ref] && ![regexp "\\-" $alt]} {
	set ref_length [string length $ref]
	set alt_length [string length $alt]
	
	if {$ref_length > $alt_length} {
	    #Here we have a deletion
	    if {[string range $ref 0 [expr {$alt_length - 1}]] eq $alt} {
		set ref [string range $ref $alt_length end]
		set alt "-"
		incr pos $alt_length
	    } else {
		if {[info exists g_VaRank(DEBUG)]} {puts "WARNING: could not solve this deletion chr$chrom at position $pos $ref>$alt, continue as it is."}
	    } 				    
	    #if {[catch {regsub "^[set alt]" $ref "" ref} Message]} {}
	} elseif {$ref_length < $alt_length} {
	    #Here we have a insertion
	    if {[string range $alt 0 [expr {$ref_length - 1}]] eq $ref} {
		set alt [string range $alt $ref_length end]
		set ref "-"
		incr pos [expr {$ref_length-1}]
	    } else {
		if {[info exists g_VaRank(DEBUG)]} {puts "WARNING: could not solve this insertion chr$chrom at position $pos $ref>$alt, continue as it is."}
	    } 				    
	    
	    #if {[catch {regsub "^[set ref]" $alt "" alt} Message]} {puts $Message;puts "Skipping $chrom $pos $ref, too big";continue}
	}
    }

    set simples [SimplifyVariation $pos $ref $alt]
    set pos [lindex $simples 0]
    set ref [lindex $simples 1]
    set alt [lindex $simples 2]

    set ID "${chrom}_${pos}_${ref}_${alt}"
    WriteTextInFile "$ID\t$chrom\t$nativepos\t$nativeref\t$nativealt" $saveIDfile
    
    return $ID
}

proc SimplifyVariation {pos ref alt} {

    # SNV like...: pos 12345, ATGG > ATTG (length-ref = length-alt)
    # ==> Transformed with: pos 12347, G>T (removes the two identical ends)

    # INS like... : pos 12345, CA > CATT
    # ==> Transformed with: pos 12346, A > ATT

    # DEL like... : pos 12345, CATT > CA
    # ==> Transformed with: pos 12346, ATT > A

    # No simplification for something like: TGCTGCTGCCAC > CTGCTGC
    #                                          -------

    # No simplification if ref or alt > 400 NT


    if {[string length $ref] > 400 || [string length $alt] > 400} {
	# Special case not treated:
	return "$pos $ref $alt"
    }

    ## VCF V4.2 specification: The '*' allele is reserved to indicate that the allele is missing due to a upstream deletion
    if {$alt == "*"} {
	return "$pos $ref $alt"
    }

    set newalt "$alt"
    set newref "$ref"
    set newpos $pos

    if {[string length $ref] == [string length $alt]} {
        ### Case of a SNV with two identical ends to remove:
	set test 0
	foreach NTref [split $ref ""] NTalt [split $alt ""] {
	    if {$test} {
		append addref $NTref
		append addalt $NTalt
		if {$NTref != $NTalt} {
		    append newref $addref
		    append newalt $addalt
		    set addref ""
		    set addalt ""
		}
		continue
	    }
	    if {$NTref == $NTalt} {
		incr newpos
	    } else {
		set test 1
		set newref $NTref
		set newalt $NTalt
		set addref ""
		set addalt ""
	    }
	}
    } elseif {[string length $ref] < [string length $alt]} {
        ### Case of an INS:
	if {[regexp "^${ref}(.*)$" $alt match ins]} {
	    set commun [string index $ref end]
	    set newref $commun
	    set newalt "${commun}${ins}"
	    set newpos [expr {$pos+[string length $ref]-1}]
	}
    } else {
        ### Case of a DEL:
	if {[regexp "^${alt}(.*)$" $ref match del]} {
	    set commun [string index $alt end]
	    set newalt $commun
	    set newref ${commun}${del}
	    set newpos [expr {$pos+[string length $alt]-1}]
	}
    }

    return "$newpos $newref $newalt"

    # puts "42456670 C > CTCTT ==> [SimplifyVariation 42456670 C CTCTT] // 42456670 C > CTCTT" 
    # puts "12345 A > G ==> [SimplifyVariation 12345 A G] // 12345 A > G"
    # puts "12345 ATGG > ATTG ==> [SimplifyVariation 12345 ATGG ATTG] // 12347 G > T"
    # puts "12345 CA > TTCA ==> [SimplifyVariation 12345 CA TTCA] // 12345 CA > TTCA"
    # puts "12345 CA > CATT ==> [SimplifyVariation 12345 CA CATT] // 12346 A > ATT"
    # puts "12345 CATT > CA ==> [SimplifyVariation 12345 CATT CA] // 12346 ATT > A"
    # puts "12 GCCT > GC ==> [SimplifyVariation 12 GCCT GC] // 13 CCT > C"

}



## Parsing of VCF input file(s).
## Output:
## 	2 global variables:
##	- g_allPatients = "patient1 patient2 ..." 
##	- g_vcfINFOS(ID) = "chrom pos ref alt rsID rsValidation patient1:homhet:dp:nr:qual patient2:homhet:dp:nr:qual ..."
##	  If ID is absent from patient1, so "patient1:homhet:dp:nr:qual" is not put in memory.
##      - $g_VaRank(vcfDir)/VCF_Coordinates_Conversion.tsv
##
##  Empty data (./.) and wild type variation are filtered out for every patient individually
##
proc parseVCFfiles {} {

    global g_VaRank
    global g_vcfINFOS
    global g_vcfINFOS_Supp
    global g_allPatients
    global g_lPatientsOf

    set g_allPatients {}
    set patientsDir   $g_VaRank(vcfDir)
    
    set nbPatients_Total   0
    set nbVariations_Total 0
    set nbFiles 0

    set L_NewHeaders_VCF {}
    set L_AllHeaders_VCF {}
    set g_vcfINFOS(L_IDs) {}
		    
    set saveIDfile "$g_VaRank(vcfDir)/VCF_Coordinates_Conversion.tsv"
    file delete -force $saveIDfile
    WriteTextInFile "VariantID\tCHROM\tPOS\tREF\tALT" $saveIDfile

    foreach vcfFile [glob -nocomplain $patientsDir/*.vcf $patientsDir/*.vcf.gz] {

	incr nbFiles

	set DP_once 1
	set NR_once 1

	set l_patients {}

	puts "...parsing the VCF file ($vcfFile) ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
	
	set FirstTime 0

	set NbDPSaved    0
	set NbDPNotFound 0
	set NbDPEmpty    0
	set NbADSaved    0 
	set NbADNotFound 0
	set NbADEmpty    0

	if {[regexp ".vcf.gz$" $vcfFile]} {
	    set lLines [LinesFromGZFile $vcfFile]
	} else {
	    set lLines [LinesFromFile $vcfFile]
	}
	foreach L $lLines {
	    if {[regexp "^##" $L]} {continue}
	    if {[regexp "^#CHROM" $L]} {
		if {![regexp "^#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" $L]} {
	            puts "######################################################################"
        	    puts "Bad header line syntax of the VCF:"
		    puts "$L"
		    puts "The header line doesn't name the 8 fixed, mandatory columns."
		    puts "Should begin with \"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\""
        	    puts "######################################################################"
        	    exit

		}
		if {$FirstTime} {puts "WARNING: [file tail $vcfFile] seems to contain multiple headers";continue}

		set L [split $L "\t"]
		set i_chr    [lsearch -exact $L "#CHROM" ]; if {$i_chr    == -1} {puts "Bad header line syntax into $vcfFile. \n#CHROM column not found - Exit"; exit}
		set i_pos    [lsearch -exact $L "POS"];     if {$i_pos    == -1} {puts "Bad header line syntax into $vcfFile. \nPOS column not found - Exit"; exit}
		set i_ref    [lsearch -exact $L "REF"];     if {$i_ref    == -1} {puts "Bad header line syntax into $vcfFile. \nREF column not found - Exit"; exit}
		set i_alt    [lsearch -exact $L "ALT"];     if {$i_alt    == -1} {puts "Bad header line syntax into $vcfFile. \nALT column not found - Exit"; exit}
		set i_id     [lsearch -exact $L "ID"];      if {$i_id     == -1} {puts "Bad header line syntax into $vcfFile. \nID column not found - Exit"; exit}
		set i_qual   [lsearch -exact $L "QUAL"];    if {$i_qual   == -1} {puts "Bad header line syntax into $vcfFile. \nQUAL column not found - Exit"; exit}
		set i_valid  [lsearch -exact $L "VALID"]
		set i_filter [lsearch -exact $L "FILTER"];  if {$i_filter == -1} {puts "Bad header line syntax into $vcfFile. \nFILTER column not found - Exit"; exit}
		set i_info   [lsearch -exact $L "INFO"];    if {$i_info   == -1} {puts "Bad header line syntax into $vcfFile. \nINFO column not found - Exit"; exit}
		set i_format [lsearch -exact $L "FORMAT"];  if {$i_format == -1} {puts "Bad header line syntax into $vcfFile. \nFORMAT column not found - Exit"; exit}

		#We need to have a global list for multiple files and multiple patients per file

		set nbPatients_File   0
		set nbVariations_File 0

		for {set i [expr {$i_format+1}]} {$i < [llength $L]} {incr i} {
		    set patient [lindex $L $i]
		    incr nbPatients_File 
		    #puts $patient
		    
		    if {[lsearch -exact $g_allPatients $patient] != -1} {
			puts "\tWARNING: $patient seems to be present in different VCF files"
			lappend l_patients    $patient
		    } else {
			lappend g_allPatients $patient
			lappend l_patients    $patient
		    }

		    set i_$patient $i
		}
		set FirstTime 1

		continue
	    }

	    ## Checking that the VCF file contain an header line beginning with "#CHROM"
	    if {! $FirstTime} {
		puts "WARNING: [file tail $vcfFile] not contain a good header line syntax (not beginning with \"#CHROM\"). Exit."
		exit
	    }

	    set L [split $L "\t"]

	    ## A ameliorer!! Vero 2015/11/20
	    ## If there is only normal hom in the line (GT = 0|0 or 0/0), so there isn't mutation to analyse by VaRank
	    if {![regexp "1/0|0/1|2/0|0/2|1/1|2/2|2/1|1/2" $L] && ![regexp "1\\\|0|0\\\|1|2\\\|0|0\\\|2|1\\\|1|2\\\|2|2\\\|1|1\\\|2" $L]} {continue}

	    set chrom [lindex $L $i_chr] 
	    regsub -all "^chr" $chrom "" chrom
	    
	    set pos [lindex $L $i_pos] 
	    set ref [lindex $L $i_ref] 
	    set alt [lindex $L $i_alt] 
	    set rs  [lindex $L $i_id]  

	    if {[set g_VaRank(rsFromVCF)]=="yes"} {
		if {[isNotAnRS $rs]} {
		    set rs "NA"; set valid "NA"
		} else {
		    if {$i_valid == -1} {
			set valid "NA"
		    } else {
			set valid [lindex $L $i_valid] 
		    }
		}
	    } else {
		set rs "NA"; set valid "NA"
	    }

	    set filter [lindex $L $i_filter]
	    set info   [lindex $L $i_info] 
	    regsub -all {\"} $info "" info

	    set format [lindex $L $i_format] 
	    set qual   [lindex $L $i_qual]

	    # ALT column
	    ############
	    # ex: ALT = "A"
	    # ex: ALT = "A,AAC,AACACAC,AACACACACACACAC,*"

	    #Analysing the FORMAT COLUMN
	    ############################
	    set l_format [split $format ":"]

	    # GT: Genotype 
	    ##############
	    # 0 is for reference allele.
	    # 1 is for the first mutated allele of the ALT column
	    # 2 is for the second mutated allele of the ALT column
	    # 3 is for the third mutated allele of the ALT column
	    # ...etc
	    # GT value: "0|0 = hmz reference", "1|0 or 0|1 or 2|0 or 0|2 = htz", "1|1 or 2|2 = mutated hmz", "2|1 or 1|2 = double htz".
	    # In case with GT = 1|2 or 2|1: represents for VaRank 2 different ID (2 htz). Ok, will be see if deleterious in the compound htz.
	    #
	    set  j_gt [lsearch -exact $l_format "GT"]
	    if {$j_gt == -1} {
		puts "FORMAT:$format";puts $L
		puts "\tERROR: GT (genotype) absent at least once from the FORMAT column - Exit"; exit}

	    # DP: Read Depth
	    ################
	    # The field describes the total depth of reads that passed the caller's internal quality control metrics 
	    #
	    set  j_dp [lsearch -exact $l_format "DP"]
	    if {$j_dp == -1} {
		if {$DP_once} {
		    puts "\tWARNING: DP (total read depth) absent at least once from the FORMAT column - continue"
		    if {[info exists g_VaRank(DEBUG)]} {puts "FORMAT:$format";puts $L}
		    set DP_once 0
		}
	    }

	    # NR: Read Number
	    #################
	    # Alternatively, AD can be used this one for getting the number of read for the variant
	    # Alternatively, AC can be used this one for getting the number of read for the variant
	    # Alternatively, AO can be used this one for getting the number of read for the variant
	    # Alternatively, DV can be used this one for getting the number of read for the variant (samtools mpileup, 2016, used by Quentin)
	    #                DV = "Number of high-quality non-reference bases"
	    set L_j_nr {}
	    set j_nr [lsearch -exact $l_format "NR"]
	    if {$j_nr ne -1} {
		lappend L_j_nr $j_nr
	    }
	    set  j_nr [lsearch -exact $l_format "AD"]
	    if {$j_nr ne -1} {
		lappend L_j_nr $j_nr
	    }
	    set j_nr  [lsearch -exact $l_format "AC"]
	    if {$j_nr ne -1} {
		lappend L_j_nr $j_nr
	    }
	    #USE AO for allele read depth;AO from Life technologies (PGM)
	    set  j_nr [lsearch -exact $l_format "AO"]
	    if {$j_nr ne -1} {
		lappend L_j_nr $j_nr
	    }
	    set j_nr  [lsearch -exact $l_format "DV"]
	    if {$j_nr ne -1} {
		lappend L_j_nr $j_nr
	    }
	    if {$L_j_nr eq {}} {
		if {$NR_once} {
		    puts "\tWARNING: NR, AD, AC, AO and DV (read depth) absent at least once from the FORMAT column - continue"
		    set NR_once 0
		    #puts "FORMAT:$format"
		}
	    }
	

	    # Analysing the INFO COLUMN
	    # Example of value for the INFO column:
	    # AC=2,1;AF=0.045,0.023;AN=44;BaseQRankSum=-7.360e-01;ClippingRankSum=-3.800e-01;DP=174;FS=12.506;InbreedingCoeff=-0.0226;MQ=60.72;MQRankSum=-3.400e-01;POSITIVE_TRAIN_SITE;QD=16.81;ReadPosRankSum=-7.360e-01;SOR=1.188;VQSLOD=1.16;culprit=FS
	    # We need to extract the INFO column and additional information from there.
	    # Storing is done here but extraction is done while ranking to ensure the big picture
	    #
	    # To integrate the info column from the VCF we need to generate a list of cumulated headers from all patients and all input files
	    #
	    set l_duoInfos_VCF ""
	    set   DPINFO    ""
	    set i_DPINFO    ""
	    
	    if {$info!="" && $info!="." && $info!="NA"} {
		set l_duoInfos_VCF [split $info ";"]
		
		#No DP information from the FORMAT column so we try to rescue somehow
		if {$j_dp == -1} {
		    set  i_DPINFO [lsearch -regexp $l_duoInfos_VCF "^DP="]
		    if {$i_DPINFO!="-1"} {
			set DPINFO [split [lindex $l_duoInfos_VCF $i_DPINFO] "="]
			set DPINFO [lindex $DPINFO 1]
		    }
		}

		#save the different headers from INFO column (L_NewHeaders_VCF)
		if {[set g_VaRank(vcfInfo)]=="yes"} {
		    foreach duoInfos_VCF $l_duoInfos_VCF {
			#Info is usually Header=Data
			#Some field could be only data and no header, in that case we will get data, and set data as header too (set as header=header)
			set duo_Header_Infos_VCF [split $duoInfos_VCF "="]
			set Header [lindex $duo_Header_Infos_VCF 0]
			
			if {[llength $duo_Header_Infos_VCF]>2} {puts ">>>>>>>>>$ID $duo_Header_Infos_VCF ----- $Header"}
			
			if {[lsearch -exact $L_AllHeaders_VCF $Header]==-1} {
			    lappend L_AllHeaders_VCF $Header
			    if {$g_VaRank(vcfFields)=="all" || [lsearch -exact $g_VaRank(vcfFields) $Header]!=-1} {
				lappend L_NewHeaders_VCF $Header
			    }
			}    
		    }
		}
	    }
 

	    #Variant/allele position, in case of multiple allele to get the correct data
	    #The counting includes the fact that position at index 0 is the reference position
	    set k 0 ; #ALT=G,C : GT 0/0 is for WT, GT 1/1 for hmz G, GT 2/2 for hmz C.
	    set naltn [llength [split $alt ","]]
	    foreach altn [split $alt ","] {
		incr k
		if {$altn == ""} {continue}
		## Structural variants (ALT=<ID=type,Description=description>) from the VCF input file are not analysed by VaRank
		if {[regexp "<.+>" $altn]} {continue}
		set ID   [AssignID $chrom $pos $ref $altn]
		set sID  [split $ID "_"]
		set pos_tmp  [lindex $sID 1]
		set ref_tmp  [lindex $sID 2]
		set altn_tmp [lindex $sID 3]

		if {![info exists g_vcfINFOS($ID)]} {
		    set g_vcfINFOS($ID) "$chrom $pos_tmp $ref_tmp $altn_tmp $rs $valid" ; # To add only if there is at least 1 altn ok.
		    lappend g_vcfINFOS(L_IDs) $ID
		}

		set firstTime 1 ; # some variant are present in VCF but absent from samples (GT always like 0/0)

		# Analysing the 'sample' COLUMNS
		foreach patient $l_patients {

		    if {![info exists i_$patient]} {continue}; # In case of several VCF files in input 
		    set value [lindex $L [set i_$patient]]
 
		    #Empty data are represented by ./. and filtered out
		    #if {$value == "./."} {continue} ; # 2015/11/01 seems not to be necessary...

		    #List of values
		    set lvalue [split $value ":"]
		    #Genotype for the patient
		    set gt [lindex $lvalue $j_gt]

		    # Wild type variation (actually not variation) are filtered out
		    #
		    # If a call cannot be made for a sample at a given locus, 
		    # '.' is specified for each missing allele in the GT field (example: GT='./.')
		    # WARNING: VCF with at least 10 ALT? $k =10? If yes, taken into account in the following regexp.
		    if {![regexp "^$k/" $gt] && ![regexp "/$k$" $gt] && ![regexp "^$k\\\|" $gt] && ![regexp "\\\|$k$" $gt]} {continue}

		    # some variant are present in VCF but absent from samples (GT always like 0/0)
		    if {$firstTime} {
			set firstTime 0
			incr nbVariations_File
		    }
		    
		    ## 2015/11/01: Keep the './.' info???


		    #We need to extract the INFO column and additional information from there and this for each patient individually
		    #Storing is done here but extraction is done while ranking to ensure the big picture
		    #Analysing the INFO COLUMN
		    #
		    if {[set g_VaRank(vcfInfo)]=="yes"} {
			if {![info exists g_vcfINFOS_Supp($ID,$patient)]} {
			    if {$g_VaRank(vcfFields)=="all"} {
				set g_vcfINFOS_Supp($ID,$patient) [split $info ";"]
			    } else {
				foreach duoInfos_VCF [split $info ";"] {
				    set duo_Header_Infos_VCF [split $duoInfos_VCF "="]
				    set Header [lindex $duo_Header_Infos_VCF 0]
				    if {[lsearch -exact $g_VaRank(vcfFields) $Header]!=-1} {
					lappend g_vcfINFOS_Supp($ID,$patient) $duoInfos_VCF
				    }
				}
			    }
			}
		    }
		    

		    if {$j_dp != -1} {
			set  dp [lindex $lvalue $j_dp]
			# 2018-07-27: debuggage for the old VCF from the diag (DP present in the format column but not in the patient column)
			if {$dp eq ""} {set dp "0"}
			if {$dp=="0"} {incr NbDPEmpty} else {incr NbDPSaved}
		    } elseif {$DPINFO != "" && [llength $l_patients]<=1} {
			set  dp $DPINFO
			if {$dp=="0"} {incr NbDPEmpty} else {incr NbDPSaved}
		    } else {
			set  dp "NA"
			incr NbDPNotFound
		    }

		    #In case of multiple variants (REF=A and ALT=G,T => AD=9,20,5) this should be retrieved by the indexed k
		    ## WARNING: in some old VCF, AD is not given for the REF (REF=A and ALT=G,T => AD=20,5). 
		    if {$L_j_nr ne {}} {
			foreach j_nr $L_j_nr {
			    set nr [lindex $lvalue $j_nr] 
			    set nr [split $nr ","]
			    if {[llength $nr] eq [expr {$naltn+1}]} {
				set nr [lindex $nr $k]
			    } elseif {[llength $nr] eq $naltn} {
				#set nr [lindex [split $nr ","] [expr {$k-1}]]			
				set nr [lindex $nr [expr {$k-1}]]			
			    } else {
				set nr "NA"
				#puts "FORMAT ERROR:\n$naltn ALT given for [llength $nr] allelic depth.$L\nPlease check the read count: [lindex $lvalue $j_nr]\nExit"; exit 
			    }
			    if {$nr eq "0" || $nr eq "."} {continue}
			    break
			}
			if {$nr eq "0" || $nr eq "."} {
			    incr NbADEmpty
			} else {
			    incr NbADSaved
			}
		    } else {
			set nr "NA"
			incr NbADNotFound
		    }

		    #Getting the genotype from the variants.
		    #
		    #Multiple ways to get the information, from the GT field, from the info column and various filed in there or recalculate from depth and read depth
		    #
		    set StatutHomHet "NA"
		    
		    #puts "$gt-$dp-$nr-$k"
		    #puts "[llength $nr]==[expr {$naltn+1}]"

		    #WARNING: We can imagine to have 10 different ALT. e.g. GT=11/11
		    if {[regexp "^$k\(/|\\\|\)$k$" $gt]} {
			set n 2
		    } elseif {[regexp "^$k\(/|\\\|\)" $gt]} {
			set n 1
		    } elseif {[regexp "\(/|\\\|\)$k$" $gt]} {
			set n 1
		    } else {
			set n 0
		    }

		    #puts "$gt - $n"
		    if {[set g_VaRank(Homstatus)]=="yes"} {
			#Recompute the ratio for homozygosity/heterozygosity detection
			#
			#Ratio for homozygosity
			set Ratio    0
			set RatioHom [set g_VaRank(Homcutoff)]
			if {($nr !="." && $nr != "NA" && $nr != "0") && ($dp !="." && $dp != "NA" && $dp != "0")} {
			    if {[info exists g_VaRank(DEBUG)]} {puts "NR:$nr, DP:$dp"}

			    set Ratio [format "%.0f" [expr {$nr*100.0/$dp}]]

			    if {[info exists g_VaRank(DEBUG)]} {puts "$nr $dp gives a ratio $Ratio"}

			    if {$Ratio>=$RatioHom} {set StatutHomHet "hom"} else {set StatutHomHet "het"}
			} else {
			    #If no ratio can be computed (eg, BGI files with no coverage), we keep the old way
			    if {$n == 1} {set StatutHomHet "het"} elseif {$n == 2} {set StatutHomHet "hom"}
			}
		    } else {
			#Get the homozygosity/heterozygosity detection from the files
			if {$n == 1} {
			    set StatutHomHet "het"
			} elseif {$n == 2} {
			    set StatutHomHet "hom"
			} else {
			    #In case the ratio is not given properly (bug observed in some tools), recompute the ratio for homozygosity/heterozygosity detection
			    #Many programs have trouble in given GT calls for multiple allele in vcf
			    #
			    set Ratio    0
			    set RatioHom [set g_VaRank(Homcutoff)]
			    
			    if {($nr !="." && $nr != "NA" && $nr != "0") && ($dp !="." && $dp != "NA" && $dp != "0")} {

				if {[info exists g_VaRank(DEBUG)]} {puts "NR:$nr, DP:$dp"}
				
				set Ratio [format "%.0f" [expr {$nr*100.0/$dp}]]
				
				if {[info exists g_VaRank(DEBUG)]} {puts "$nr $dp gives a ratio $Ratio"}
				
				if {$Ratio>=$RatioHom} {set StatutHomHet "hom?"} else {set StatutHomHet "het?"}
			    }
			}
		    }

		    #Now that we have defined that NA is by default for zygosity status, one has to selectively not include variants in other patients when no coverage is available
		    # Vero 2015/11/03 - BUG!!!
		    #if {$StatutHomHet == "NA" && ($nr =="." || $nr == "NA" || $nr == "0") && ($dp =="." || $dp == "NA" || $dp == "0")} {continue}
		    if {$StatutHomHet == "NA"} {continue}
		    #if {$nr =="." || $nr == "NA" || $nr == "0"} {continue}
		    if {$nr =="." || $nr == "NA"} {continue}

		    #puts "$patient:$StatutHomHet:$dp:$nr:$qual"

		    lappend g_vcfINFOS($ID) "$patient:$StatutHomHet:$dp:$nr:$qual"
		    set     g_vcfINFOS($ID,$patient) 1
		}
	    }
	}
	incr   nbPatients_Total $nbPatients_File
	incr nbVariations_Total $nbVariations_File

	puts "\tFile loaded: $nbVariations_File non redundant variation(s) in $nbPatients_File sample(s) ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])."
	puts "\t             Total Read Depth:  $NbDPSaved variant(s) are different to 0, $NbDPEmpty variant(s) are equal to 0, $NbDPNotFound are empty"
	puts "\t             Allele Read Depth: $NbADSaved variant(s) are different to 0, $NbADEmpty variant(s) are equal to 0, $NbADNotFound are empty"
    }

    if {[set g_VaRank(vcfInfo)]=="yes"} {
	set g_vcfINFOS_Supp(Header) [lsort -dictionary $L_NewHeaders_VCF]
    }

    puts "...VCF file(s) loaded: $nbFiles file(s) for $nbVariations_Total variation(s) in $nbPatients_Total sample(s) ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
    
    #foreach L [array get g_vcfINFOS] {puts $L}

    ## Check if all the patient given by family in the config file are presents in the VCF file.
    if {[info exists g_lPatientsOf]} {
	foreach fam [array names g_lPatientsOf] {

	    set lBadPatients  {}
	    set lPatients_Tmp {}

	    foreach patient $g_lPatientsOf($fam) {
		if {[lsearch -exact $g_allPatients $patient] == -1} {lappend lBadPatients $patient} else {lappend lPatients_Tmp $patient}
	    }
	    
	    if {$lBadPatients != {}} {
		puts "\t...patient names \"[join $lBadPatients " "]\" (given in config file for family study) absent from the VCF file. Updating $fam patient list ($lPatients_Tmp)."
		#puts "\t...Families not taken into account in this study"
			
		set g_lPatientsOf($fam) $lPatients_Tmp
	    }
	}
    }

    ## Patients without an attributed family in the config file have here a new family.
    ## If no family have been defined in the config file, we attribute a family for each patient of the VCF files.
    set lPasTrouve {}
    if {[info exists g_lPatientsOf]} {
	foreach patient $g_allPatients {
	    set pastrouve 1

	    foreach fam [array names g_lPatientsOf] {
		if {[lsearch -exact $g_lPatientsOf($fam) $patient] != -1} {set pastrouve 0}
	    }
	    if {$pastrouve} {lappend lPasTrouve $patient}
	}
    } else {
	foreach patient $g_allPatients {
	    lappend lPasTrouve $patient
	}
    }

    set lPasTrouve [lsort -dictionary $lPasTrouve]

    set i 1
    foreach patient $lPasTrouve {
	while {[info exists g_lPatientsOf(fam$i)]} {incr i}
	set g_lPatientsOf(fam$i) $patient
    }

    ## We re-sort the list of patients in g_allPatients to match those in g_lPatientsOf
    ## -> Used to write header files of ranking (barcodes lines)
    unset g_allPatients
    set   g_allPatients {}
    foreach fam [array names g_lPatientsOf] {
	foreach patient $g_lPatientsOf($fam) {
	    lappend g_allPatients $patient
	}
    }

    ## If SamOut option is set, check that the list of good samples is not empty. Else, use all samples.
    if {$g_VaRank(SamOut)!="all"} {
	set l {}    
	foreach patient $g_VaRank(SamOut) {
	    if {[lsearch -exact -nocase $g_allPatients $patient]!=-1} {
		lappend l $patient
	    }
	}
	if {[llength $l]==0} {
	    puts "...WARNING:"
	    puts "\t-SamOut option is badly defined (\"$g_VaRank(SamOut)\")."
	    puts "\tSamples should be choosen into \"$g_allPatients\""
	    puts "\tOption not used"
	    set g_VaRank(SamOut) "all" 
	} else {
	    set g_VaRank(SamOut) "$l"
	}
    }

    if {[set g_VaRank(vcfInfo)]=="yes"} {
	if {$L_AllHeaders_VCF != $L_NewHeaders_VCF} {
	    puts "...List of cumulated headers (from INFO column), available from all patients and all VCF input files:"
	    puts "\t$L_AllHeaders_VCF"
	    
	    puts "...List of headers (from INFO column) reported in output files:"
	    puts "\t$L_NewHeaders_VCF"
	}
    }

    return
}

