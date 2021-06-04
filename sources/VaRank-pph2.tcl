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

proc changeModeOfPPH2lockFiles {} {
    global g_VaRank

    ## In case of use by different users (else, pph2 will crash...)
    foreach f [glob -nocomplain $g_VaRank(pph2Dir)/scratch/lock/*] {
	catch {file attributes $f -permissions 0777} Message
    }

    return
}

proc searchProteicSequence {IDprot database} {
    
    ## Cleaning:
    foreach f [glob -nocomplain wgetz\?-id*] {
	file delete -force $f
    }
    file delete -force Wget.log

    ## Searching for a RefSeqP sequence on a regularly updated site.
    set command "/usr/bin/wget -nv -a Wget.log -N \"http://srs.ebi.ac.uk/srsbin/cgi-bin/wgetz?-id+_bIV1er7R7+-e+\\\[$database:'$IDprot'\\\]+-qnum+2+-enum+1\""

    if {[catch {set Seq [eval exec $command]} Message]} {
	set Seq ""
    } else {
	set test 1
	foreach f [glob -nocomplain wgetz\?-id*] {
	    set Seq ""
	    foreach L [LinesFromFile $f] {
		if {[regexp "^ >>refseqp" $L]} {set test 0; continue}
		if {!$test && [regexp "^<" $L]} {set test 1; break}
		if {$test} {continue}
		regsub -all " " $L "" L
		regsub -all "\[0-9\]" $L "" L
		append Seq "[string toupper $L]"
	    }
	}
    }
    
    ## Cleaning:
    foreach f [glob -nocomplain wgetz\?-id**] {
	file delete -force $f
    }
    file delete -force Wget.log
    
    # Return the sequence
    return $Seq
}



proc searchBadAA1position {IDprot Seq patientsDir} {

    ## Search into the created PPH2input file the SNV corresponding to the ID given in fasta.
    ## Check if SNV positions correspond to AA1 provided.
    ## Returns 1 if the positions are not all OK, 0 otherwise.
    foreach L [LinesFromFile $patientsDir/PPH2/PPH2input_all.txt] {
	if {![regexp "$IDprot" $L]} {continue}
	set Pos [expr {[lindex $L 1]-1}]
	set AA1 [lindex $L 2]
	if {[string index $Seq $Pos] != $AA1} {
	    #puts "\t$L"
	    #puts "\t***** Colle pas: [string index $Seq $Pos] != $AA1"
	    return 1
	} else {
	    #puts "\tOk for [lindex $L 0] - [string index $Seq $Pos] eq $AA1"
	}
    }
    
    return 0
}

proc FastaSequence_HumanUniProt args {

    #Store the uniprot identifiers and sequences

    global g_VaRank
    global SeqFasta

    regsub "sources" $g_VaRank(sourcesDir) "pph2Databases" pph2Databases
    if {$pph2Databases eq "" || ![file exists $pph2Databases] } {return ""}
    if {$g_VaRank(uniprot) eq ""} {return ""}

    set InputFile "$g_VaRank(uniprot)"
    if {![file exists $InputFile]} {puts "Inputfile UNIPROT $g_VaRank(uniprot)";return ""}

    set What [join [concat $InputFile $args] ","]

    if {[info exists SeqFasta($What)]} {return [set SeqFasta($What)]}
    if {![info exists SeqFasta($InputFile,Loaded)]} {
	set SeqFasta($InputFile,Loaded) 1

	set ID  ""
	set Seq ""
	set SeqFasta($InputFile,L_ID) {}

	if {[regexp ".gz$" $InputFile]} {
	    set F [open "|gzip -cd $InputFile"] 
	} else {
	    set F [open "$InputFile"]
	}
	while {[gets $F Line]>=0} {
	    if {[string first ">" $Line] eq "0"} {
		if {$Seq ne ""} {
		    #if {[info exists SeqFasta($InputFile,$ID)]} {puts "$ID already seen"}
		    
		    regsub -all " "  $Seq "" Seq
		    #regsub -all {\-} $Seq "" Seq
		    
		    set SeqFasta($InputFile,$ID) $Seq
		    #puts "$ID [set SeqFasta($InputFile,$ID)]"
		    set ID  ""
		    set Seq ""
		}

		#>tr|A0JP02|A0JP02_HUMAN PLEKHA5 protein OS=Homo sapiens GN=PLEKHA5 PE=2 SV=1
		#>sp|A0M8Q6|LAC7_HUMAN Ig lambda-7 chain C region OS=Homo sapiens GN=IGLC7 PE=1 SV=2

		set  ID [string trim [lindex [split $Line "|"] 1]]
		if {$ID eq ""} {puts "WARNING $Line"}
		lappend SeqFasta($InputFile,L_ID) $ID
	    } else {
		append Seq [string trim $Line]
	    }
	}
	close $F
	
	if {$Seq ne ""} {
	    #if {[info exists SeqFasta($InputFile,$ID)]} {puts "$ID already seen"}
	    
	    regsub -all " "  $Seq "" Seq
	    #regsub -all {\-} $Seq "" Seq

	    set SeqFasta($InputFile,$ID) $Seq
	    set ID  ""
	    set Seq ""
	}
	
	set SeqFasta($InputFile,Loaded) 1

	if {[info exists SeqFasta($What)]} {return [set SeqFasta($What)]}
    } else {
	return ""
    }
}

proc FastaSequence_HumanRefSeq args {

    # Store the refseq identifiers and sequences

    global g_VaRank
    global SeqFasta

    regsub "sources" $g_VaRank(sourcesDir) "pph2Databases" pph2Databases
    if {$pph2Databases eq "" || ![file exists $pph2Databases] } {return ""}
    if {$g_VaRank(refseq) eq ""} {return ""}
    
    set InputFile "$g_VaRank(refseq)"
    if {![file exists $InputFile]} {return ""}

    set What [join [concat $InputFile $args] ","]

    if {  [info exists SeqFasta($What)]} {return [set SeqFasta($What)]}
    if {! [info exists SeqFasta($InputFile,Loaded)]} {
	set SeqFasta($InputFile,Loaded) 1

	set ID  ""
	set Seq ""
	set SeqFasta($InputFile,L_ID) {}

	if {[regexp ".gz$" $InputFile]} {
	    set F [open "|gzip -cd $InputFile"] 
	} else {
	    set F [open "$InputFile"]
	}
	while {[gets $F Line]>=0} {
	    if {[string first ">" $Line] eq "0"} {
		if {$Seq ne ""} {
		    #if {[info exists SeqFasta($InputFile,$ID)]} {puts "$ID already seen"}

		    regsub -all " "  $Seq "" Seq
		    #regsub -all {\-} $Seq "" Seq

		    set SeqFasta($InputFile,$ID) $Seq
		    #puts "$ID [set SeqFasta($InputFile,$ID)]"
		    set ID  ""
		    set Seq ""
		}

		# Old files:
		#>gi|53292629|ref|NP_001005405.1| keratin-associated protein 5-11 [Homo sapiens]
		#>gi|52317162|ref|NP_001004713.1| olfactory receptor 1I1 [Homo sapiens]

		# New files (2018/04/17):
		#>NP_001078836.1 zinc finger protein 655 isoform b [Homo sapiens]
		#>NP_001340251.1 immunoglobulin superfamily member 11 isoform f precursor [Homo sapiens]
		if {[regexp "^>gi" $Line]} {
		    set ID [string trim [lindex [split $Line "|"] 3]]
		} else {
		    regexp "^>(\[^ \]+) " $Line match ID
		}

		if {$ID eq ""} {puts "WARNING $Line"}
		lappend SeqFasta($InputFile,L_ID) $ID
	    } else {
		append Seq [string trim $Line]
	    }
	}
	close $F
	
	if {$Seq!=""} {
	    #if {[info exists SeqFasta($InputFile,$ID)]} {puts "$ID already seen"}
	    
	    regsub -all " "  $Seq "" Seq
	    #regsub -all {\-} $Seq "" Seq

	    set SeqFasta($InputFile,$ID) $Seq
	    set ID  ""
	    set Seq ""
	}
	
	set SeqFasta($InputFile,Loaded) 1

	if {[info exists SeqFasta($What)]} {return [set SeqFasta($What)]}
    } else {
	return ""
    }
}


##
## - Create the PPH2 input file from g_ANNOTATION
##   ->  Format of the input file created (corresponding to V2.2.2 of PPH2): ID-UniProt/RefSeqP position AA1 AA2
## - Create a unique fasta file (RefSeqPsequences.fasta) for the RefSeqP ID of all the patients.
##
## OUTPUT:
##	- $patientsDir/PPH2/PPH2input_all.txt
## 	- $patientsDir/PPH2/RefSeqPsequences.fasta
##
## Use of the wgetz command at: "http://srs.ebi.ac.uk/srsbin/cgi-bin/wgetz?" wich seems to be regularly updataed. Site is retired now since 19/12/2013
##				"http://bioinfo.ceinge.unina.it/srs7131bin/cgi-bin/wgetz?" which seems not to be regularly updated.
##
## TO KNOW : In the PPH2 input file and in the "RefSeqPsequences.fasta" file:
##            - The RefSeqP ID found with ebi haven't the version number at the end (eg: NP_001185763)
##            - The RefSeqP ID found with bioinfo.ceinge.unina.it have the version number at the end (eg: NP_001185763.1)
##
##	      The SNV without associated sequence are written into the PPH2 input.
##	      => This means that if you do not manually add it in "RefSeqPsequences.fasta", PPH2 will crashed on this SNV.
proc createPPH2Input {} {

    global g_VaRank
    global g_ANNOTATION

    #set lRefSeqID ""
    
    set lRefSeqID(L_IDs) {}
    
    ## pph2Dir must be given in the config file or in the command line.
    ###################################################################
    if {$g_VaRank(pph2Dir) eq ""} {return}

    set patientsDir $g_VaRank(vcfDir)
    if {![file exists "$patientsDir/PPH2"]} {file mkdir "$patientsDir/PPH2"}

    set PPH2InputFile "$patientsDir/PPH2/PPH2input_all.txt"
    puts "...creation of the PPH2 input file ($PPH2InputFile) ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"

    ## If the RefSeqPsequences.fasta file already exists:
    ## - Load all the sequences of this file (without headline, without "\n") into $Seq($IDprot).
    ## - Load the IDprot of these sequences into $InfosRefSeqPsequences

    #set InfosRefSeqPsequences ""

    set fastaFile "$patientsDir/PPH2/RefSeqPsequences.fasta"
    if {[file exists $fastaFile]} {
	foreach L [LinesFromFile $fastaFile] {
	    if {$L eq ""} {continue}
	    if {[regexp "^>(.*)" $L match IDprot]} {
		#append InfosRefSeqPsequences "$IDprot "

		set InfosRefSeqPsequences($IDprot) 1
		
		#Not sure we need to store the short version see with Vero
		set IDprotShort ""
		set IDprotShort $IDprot
		
		regsub {\.[0-9]+$} $IDprotShort "" IDprotShort

		set InfosRefSeqPsequences($IDprotShort) 1
		
		continue
	    }
	    append Seq($IDprot) "$L"
	}
    }

    set PPH2InputFile "$patientsDir/PPH2/PPH2input_all.txt"
    file delete -force $PPH2InputFile
    set Info ""
    puts "\t...Extracting identifiers and sequences from annotation data, Uniprot and RefSeqP fasta files."
    set L [split $g_ANNOTATION(#id) "\t"]

    set i_refseq  [lsearch -regexp $L "^protein$"]; # This annotation doesn't exists with SNPEFF annotation.
    set i_uniprot [lsearch -regexp $L "^Uniprot$"]; if {$i_uniprot eq -1} {puts "column number not found for Uniprot - Exit"; exit}
    set i_pos     [lsearch -regexp $L "^posAA$"];   if {$i_pos     eq -1} {puts "column number not found for posAA - Exit"; exit}
    set i_AA      [lsearch -regexp $L "^wtAA_1$"];  if {$i_AA      eq -1} {puts "column number not found for wtAA_1 - Exit"; exit}
    set i_var     [lsearch -regexp $L "^varAA_1$"]; if {$i_var     eq -1} {puts "column number not found for varAA_1 - Exit"; exit}
    set i_varType [lsearch -regexp $L "^varType$"]; if {$i_varType  eq -1} {puts "column number not found for varType - Exit"; exit}
    set i_effect  [lsearch -regexp $L "^codingEffect$"]; if {$i_effect  eq -1} {puts "column number not found for codingEffect - Exit"; exit}

    set nbTestedUni   0
    set nbNotFoundUni 0
    set nbFoundUni    0
    set L_NotFoundUni {}

    set nbTestedRef   0
    set nbNotFoundRef 0
    set nbFoundRef    0
    set L_NotFoundRef {}

    #puts "From Uniprot [join [FastaSequence_HumanUniProt L_ID] ","]"
    #puts "From Refseq [join [FastaSequence_HumanRefSeq L_ID] ","]"
    
    foreach ID [array name g_ANNOTATION] {

	if {$ID eq "#id"} {continue}
	foreach L $g_ANNOTATION($ID) {
	    set L [split $L "\t"] ; # do not combined the last 2 lines!

	    # WARNING: If alamut annotation files of different versions have been merged (with different number of annotation columns), 
	    # it implies a bug here. Correction: remove all alamut annotation files and rerun VaRank.

	    # Selection of the lines with: "varType" = substitution; "varLocation" = exon; "codingEffect" = missense
	    set codingEffect [lindex $L $i_effect]
	    set varType [lindex $L $i_varType]

	    # Some CNV from Pindel have a missense coding effect, but are not substitution (SNV). For example:
	    #id	                      chrom  pos      gene	protein	        Uniprot	 varType     codingEffect
	    #19_4511674_CACAGCA_TACGGTG   19     4511674  PLIN4     NP_001073869.1	Q96Q06	 delins	     missense
	    if {![regexp "missense|rare_amino_acid" $codingEffect] || ![regexp -nocase "substitution|SNP" $varType]} {continue}

	    # Selection of the UniProt ID
	    set SNVid     [lindex $L 0]
	    set IDuniprot [lindex $L $i_uniprot]
	    set Pos       [expr {[lindex $L $i_pos]-1}]
	    set AA1       [lindex $L $i_AA]

	    if {$IDuniprot != "" && $IDuniprot != "NA"} {
		# Check the position and the AA1 on the sequence of this UniProt ID.
		
		incr nbTestedUni

		if {![info exists Seq($IDuniprot)]} {
		    set Sequni ""
		    set Sequni [FastaSequence_HumanUniProt $IDuniprot]
		    set  Seq($IDuniprot) $Sequni
		    
		    if {$Seq($IDuniprot) eq ""} {
			lappend L_NotFoundUni $IDuniprot
			
			#EBI SRS webserver is retired since 19/12/2013
			#FIND ANOTHER WAY????

			#set Seq($IDuniprot) [searchProteicSequence $IDuniprot UNIPROT]
		    } 
		} 
		
		if {$Seq($IDuniprot) eq ""} {
		    set IDuniprot ""
		    incr nbNotFoundUni
		} else {
		    if {[string index $Seq($IDuniprot) $Pos] != $AA1} {
			set IDuniprot ""
			incr nbNotFoundUni
		    } else {
			incr nbFoundUni
		    }
		}
	    }


	    ## - No UniProt ID or...
	    ## - UniProt ID not found or...
	    ## - UniProt ID doesn't match for AA1...
	    if {$IDuniprot eq "" || $IDuniprot eq "NA"} {
		
		# Selection of the RefSeqP ID (with the number version at the end)
		if {$i_refseq ne "-1"} {
		    set IDrefseqp_V [lindex $L $i_refseq]
		    regsub "\\..*$" $IDrefseqp_V "" IDrefseqp
		} else {
		    set IDrefseqp_V "NA"
		    set IDrefseqp   "NA"
		}

		if {$IDrefseqp eq "" || $IDrefseqp eq "NA"} {
		    if {[info exists g_VaRank(DEBUG)]} {		    
			puts "\t...ERROR for $SNVid: No UniProt ID or RefSeqP ID !! SNV NOT WRITTEN INTO THE INPUT FILE :"
			puts "\t$L"
		    }
		    continue
		} else {
		    incr nbTestedRef

		    # Get the sequence if possible
		    # Check the position and the AA1 on the sequence of this UniProt ID.

		    set NotFound 1
		    foreach IDrefseqTmp [list $IDrefseqp_V $IDrefseqp] {

			if {![info exists Seq($IDrefseqTmp)]} {
			    set Seqref ""
			    set Seqref [FastaSequence_HumanRefSeq $IDrefseqTmp]
			    set Seq($IDrefseqTmp) $Seqref
			}
			# Check the position and the AA1 on the sequence of this RefSeqP ID
			if {$Seq($IDrefseqTmp) != ""} {
			    if {[string index $Seq($IDrefseqTmp) $Pos] != $AA1} {
				set Seq($IDrefseqTmp) ""
			    } else {
				set NotFound 0
				incr nbFoundRef
				break
			    }
			} 
		    }
		    set IDrefseqp $IDrefseqTmp
		    
		    if {$NotFound} {
			incr    nbNotFoundRef
			lappend L_NotFoundRef $IDrefseqp
		    }
		    
		    #if {![regexp $IDrefseqp $lRefSeqID] && ![regexp $IDrefseqp $InfosRefSeqPsequences]} {lappend lRefSeqID $IDrefseqp}

		    if {![info exists lRefSeqID($IDrefseqp)] && ![info exists InfosRefSeqPsequences($IDrefseqp)]} {set lRefSeqID($IDrefseqp) 1;lappend lRefSeqID(L_IDs) $IDrefseqp}
		}
		
		set l "$IDrefseqp\t[lindex $L $i_pos]\t$AA1\t[lindex $L $i_var]"
	    } else {
		set l "$IDuniprot\t[lindex $L $i_pos]\t$AA1\t[lindex $L $i_var]"
	    }
	    # Elimination of the redondance
	    if {[regexp $l $Info]} {continue}
	    append Info "$l\n"
	}
    }

    if {[info exists g_VaRank(DEBUG)]} {
	puts "\t...UniProt: total tested $nbTestedUni for $nbFoundUni found and $nbNotFoundUni not found ([join $L_NotFoundUni ","])"
	puts "\t...RefSeqP: total tested $nbTestedRef for $nbFoundRef found and $nbNotFoundRef not found ([join $L_NotFoundRef ","])"
    } else {
	puts "\t...UniProt: total tested $nbTestedUni for $nbFoundUni found and $nbNotFoundUni not found"
	puts "\t...RefSeqP: total tested $nbTestedRef for $nbFoundRef found and $nbNotFoundRef not found"
    }

    if {$Info eq ""} {
	puts "\t...no missense found, skipping PPH2 step."
	return 0
    } else {
	puts "\t...writing the PPH2 input file"
	WriteTextInFile $Info $PPH2InputFile
    }

    puts "\t...creation of the sequences file for PPH2 ($fastaFile)"

    # We work in $patientsDir. Creation/Complementation of the "RefSeqPsequences.fasta" file.
    foreach IDrefseqp [set lRefSeqID(L_IDs)] {
	#foreach IDrefseqp $lRefSeqID {}

	if {[info exists InfosRefSeqPsequences($IDrefseqp)]} {continue}

	#if {[regexp $IDrefseqp $InfosRefSeqPsequences]} {continue}
	
	## Load the sequence if it is not.
	# (<=> Possible if VaRank was running with an already existing RefSeqPsequences.fasta,
	#       and if the sequence is not present in the "RefSeqPsequences.fasta" file)
	if {![info exists Seq($IDrefseqp)]} {
	    set Seq($IDrefseqp) [searchProteicSequence $IDrefseqp REFSEQP]
	}

	# Search if the sequence is OK or not.
	if {$Seq($IDrefseqp) eq "" || [searchBadAA1position $IDrefseqp $Seq($IDrefseqp) $patientsDir]} {
	    ## The RefSeqP sequences not found have their SNV written in the PPH2 input files.
	    ## Sequences should be added manually (if found on the net) in RefSeqPsequences.fasta file
	    ## Sequences should match whith the AA1 position.
	    ## Else, remove the SNV line in PPH2 inputfile.
	    # puts "\tRefSeqP ID ($IDrefseqp) not found or doesn't match for all AA1. SEQUENCE NOT WRITTEN INTO RefSeqPsequences.fasta!!"
	    continue
	}

	## Sequences loaded are not in the good format. No "\n".
	set fastaOk ">$IDrefseqp"
	set l [string length $Seq($IDrefseqp)]
	for {set i 0} {$i < $l} {incr i 60} {
	    append fastaOk "\n[string range $Seq($IDrefseqp) $i [expr {$i+59}]]"
	}

	WriteTextInFile $fastaOk $fastaFile
    }

    return
}

## Return 0 if all the variations of the input file have been analysed or if the input file is empty.
## Else return 1 and create a inputFile.tmp with the variations not yet analysed.
proc pph2Step1IsNotCompleted {inputFile featFile errorFile} {

    if {[file size $inputFile] eq 0} {
	puts "Empty input file for PPH2: $inputFile"
	return 0
    }

    set lIDinput [LinesFromFile $inputFile]

    file delete -force $inputFile.tmp

    if {![file exists $featFile] && ![file exists $errorFile]} {
	file copy -force $inputFile $inputFile.tmp
	return 1
    }
    
    ## Load in lIDana the ID already analysed by PPH2.
    set lIDana {}
    if {[file exists $featFile]} {
	foreach L [LinesFromFile $featFile] {
	    if {$L eq ""} {continue}
	    #lappend lIDana [join [lrange $L 0 3] " "]
	    set ID [join [lrange $L 0 3] " "]
	    if {![info exists TabIDana($ID)]} {set TabIDana($ID) 1;lappend lIDana $ID} else {continue}
	}
    }
    if {[file exists $errorFile]} {
	foreach L [LinesFromFile $errorFile] {
	    if {$L eq ""} {continue}
	    if {![regexp ":" $L]} {
		if {[llength $L] eq 4} {
		    #lappend lIDana [join [lrange $L 0 3] " "]

		    set ID [join [lrange $L 0 3] " "]
		    if {![info exists TabIDana($ID)]} {set TabIDana($ID) 1;lappend lIDana $ID} else {continue}
		}
	    }
	}
    }

    if {$lIDana eq ""} {
	file copy -force $inputFile $inputFile.tmp
	return 1
    }

    set test 0
    foreach L [LinesFromFile $inputFile] {
	if {$L eq ""} {continue}
	set ID [join [lrange $L 0 3] " "]
	if {![info exists TabIDana($ID)]} {
	    #if {[lsearch -exact $lIDana "$ID"] eq -1} {}
	    set test 1
	    WriteTextInFile $L $inputFile.tmp
	}
    }

    return $test
}

## Return 0 if all the variations of the input file have been analysed or if the input file is empty.
## Else return 1.
proc pph2Step2IsNotCompleted {inputFile HumVarOutput errorFile} {

    if {![file exists $HumVarOutput]} {
	return 1
    }

    set lIDinput [LinesFromFile $inputFile]

    ## Load in lIDana the ID already analysed by PPH2.
    set lIDana {}
    foreach L [LinesFromFile $HumVarOutput] {
	#lappend lIDana [join [lrange $L 0 3] " "]
	set ID [join [lrange $L 0 3] " "]
	if {![info exists TabIDana($ID)]} {set TabIDana($ID) 1;lappend lIDana $ID} else {continue}
    }

    if {[file exists $errorFile]} {
	foreach L [LinesFromFile $errorFile] {
	    if {$L eq ""} {continue}
	    if {![regexp ":" $L]} {
		if {[llength $L] eq 4} {
		    #lappend lIDana [join [lrange $L 0 3] " "]
		    set ID [join [lrange $L 0 3] " "]
		    if {![info exists TabIDana($ID)]} {set TabIDana($ID) 1;lappend lIDana $ID} else {continue}
		}
	    } 
	}
    }

    set test 0
    foreach L [LinesFromFile $inputFile] {
	if {$L eq ""} {continue}
	set ID [join [lrange $L 0 3] " "]
	if {![info exists TabIDana($ID)]} {
	    #if {[lsearch -exact $lIDana "$ID"] eq -1} {}
	    return 1
	}
    }
    
    return 0
}

## CONTEXT:
## With a PPH2 input file containing several SNV, if PPH2 crashes on a SNV line so PPH2 crashes.
## <=> PPH2 doesn't run on all the SNV lines.
##
## USE OF THIS PROC :
## runPPH2-1by1 runs PPH2 independently for each SNV line of the PPH2 input file
## So, if PPH2 crashes, PPH2 will continue to run with the following SNV line.
##
## OUTPUT :
## - $patientsDir/PPH2/PPH2features_all.txt
## - $patientsDir/PPH2/PPH2errors_all.txt (containing all the SNV lines where PPH2 has crashed)
##
proc runPPH2-1by1 {} {

    global g_VaRank

    set patientsDir $g_VaRank(vcfDir)
    
    ## pph2Dir must be given in the config file or in the command line
    if {$g_VaRank(pph2Dir) eq ""} {
	puts "...PPH environment variable not specified, not running local installation of PPH2"
	return
    }

    set PPH2inputFile "$patientsDir/PPH2/PPH2input_all.txt"
    if {![file exists $PPH2inputFile]} {puts "$PPH2inputFile doesn't exist. Exit."; exit}

    puts "...running PPH2 ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
    set inputUnite        "$patientsDir/PPH2/unite.PPH2input"
    set pph2FeaturesUnite "$patientsDir/PPH2/unite.PPH2features"
    set pph2FeaturesAll   "$patientsDir/PPH2/PPH2features_all.txt"
    set pph2ErrorAll      "$patientsDir/PPH2/PPH2errors_all.txt"

    foreach file "$inputUnite $pph2FeaturesUnite" {file delete -force $file}
    
    ################################
    # First step of PPH2 (long step)
    ################################
    if {![pph2Step1IsNotCompleted $PPH2inputFile $pph2FeaturesAll $pph2ErrorAll]} {
	puts "\t...pph2 step 1 already done and completed, continue"
    } else {
	if {![file exists $pph2FeaturesAll]} {
	    puts "\t...running pph2 step 1 ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
	} else {
	    puts "\t...pph2 step 1 already done but is not completed. Running again ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
	}
	if {[file exists $patientsDir/PPH2/RefSeqPsequences.fasta]} {        
	    set PPH2command "$g_VaRank(pph2Dir)/bin/run_pph.pl -s $patientsDir/PPH2/RefSeqPsequences.fasta $inputUnite > $pph2FeaturesUnite"
	} else {
	    set PPH2command "$g_VaRank(pph2Dir)/bin/run_pph.pl $inputUnite > $pph2FeaturesUnite"
	}
	## The first SNV of the input is analysed to keep the headline of the PPH2 output: "#o_acc  o_acc  o_pos  o_aa1   o_aa2   snp_id ..."
	foreach L [LinesFromFile $PPH2inputFile.tmp] {
	    if {$L eq ""} {continue}
	    ReplaceTextInFile $L $inputUnite

	    catch {eval exec $PPH2command} Message

	    if {[regexp -nocase "error" $Message]} {
		if {![file exists $pph2FeaturesAll]} {
		    WriteTextInFile [lindex [LinesFromFile $pph2FeaturesUnite] 0] $pph2FeaturesAll
		}
		WriteTextInFile [ContentFromFile $inputUnite] $pph2ErrorAll
		set error [lindex [LinesFromFile $pph2FeaturesUnite] end]
		if {![regexp "#o_snp_id|#o_acc" $error]} {
		    WriteTextInFile $error $pph2ErrorAll
		}
		WriteTextInFile "$Message\n" $pph2ErrorAll
	    } else {
		if {![file exists $pph2FeaturesAll]} {
		    WriteTextInFile [ContentFromFile $pph2FeaturesUnite] $pph2FeaturesAll
		} else {
		    WriteTextInFile [lindex [LinesFromFile $pph2FeaturesUnite] end] $pph2FeaturesAll
		}
	    }	
	    break		
	}
	
	# Treatment of the following SNV
	foreach L [lrange [LinesFromFile $PPH2inputFile.tmp] 1 end] {
	    if {$L eq ""} {continue}
	    ReplaceTextInFile $L $inputUnite
	    catch {eval exec $PPH2command} Message
	    if {[regexp -nocase "error" $Message]} {
		WriteTextInFile [ContentFromFile $inputUnite] $pph2ErrorAll
		set error [lindex [LinesFromFile $pph2FeaturesUnite] end]
		if {![regexp "#o_snp_id|#o_acc" $error]} {
		    WriteTextInFile $error $pph2ErrorAll
		}
		WriteTextInFile "$Message\n" $pph2ErrorAll
	    } else {
		WriteTextInFile [lindex [LinesFromFile $pph2FeaturesUnite] end] $pph2FeaturesAll
	    }
	}

	# Cleaning
	foreach file "$inputUnite $pph2FeaturesUnite" {file delete -force $file}
    } 

    ##############################
    # Second step of PPH2 (fast) :
    ##############################
    set pph2HumVarOutput "$patientsDir/PPH2/PPH2humVar_all.txt"
    if {![pph2Step2IsNotCompleted $PPH2inputFile $pph2HumVarOutput $pph2ErrorAll]} {
	puts "\t...pph2 step 2 already done and completed, continue "
    } else {
	if {![file exists $pph2HumVarOutput]} {
	    puts "\t...running pph2 step 2 ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
	} else {
	    puts "\t...pph2 step 2 already done but is not completed. Running again ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
	}
	set PPH2command "$g_VaRank(pph2Dir)/bin/run_weka.pl -l $g_VaRank(pph2Dir)/models/HumVar.UniRef100.NBd.f11.model $pph2FeaturesAll >& $pph2HumVarOutput"
	if {[catch {eval exec $PPH2command} Message]} {
	    puts "$PPH2command"
	    puts "ERROR : $Message"
	}
    }

    return
}

