############################################################################################################
# VaRank 1.5                                                                                               #
#                                                                                                          #
# VaRank: a simple and powerful tool for ranking genetic variants                                          #
#                                                                                                          #
# Copyright (C) 2016-2021 Veronique Geoffroy (veronique.geoffroy@inserm.fr)                                # 
#                         Jean Muller (jeanmuller@unistra.fr)                                              # 
#                                                                                                          #
# Please cite the following article:                                                                       #
# Geoffroy V.*, Pizot C.*, Redin C., Piton A., Vasli N., Stoetzel C., Blavier A., Laporte J. and Muller J. #
# VaRank: a simple and powerful tool for ranking genetic variants.                                         #
# PeerJ. 2015 (10.7717/peerj.796)			       		                                   #
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

proc checkIfCorrupted {} {

    global g_VaRank

    ## Checking if pre-existing "AlamutAnnotations_all.txt" file is corrupted (not the same length on all lines)
    set VCFdir $g_VaRank(vcfDir)
    set annFile "$VCFdir/Alamut/AlamutAnnotations_all.txt"

    if {![file exists $annFile]} {return}

    puts "...checking pre-existing AlamutAnnotations_all.txt file ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
    set lengthL 0
    set f [open $annFile]
    while {![eof $f]} {
	set L [gets $f]
	if {$L == ""} {continue}
	# calcul of the header length
	if {[regexp "^#id" $L]} {
	    set lengthL [llength [split $L "\t"]]
	    continue
	}
	if {[regexp "^#" $L]} {continue}
	if {[llength [split $L "\t"]] != $lengthL} {
	    if {$lengthL == 0} {
		puts "\tHeader is not present on the first lines"
	    } else {
		puts "AlamutAnnotations_all.txt file seems to be corrupted:"
		puts "\tHeader length = $lengthL"
		puts "\tAnnotations line not with the same length ([llength [split $L "\t"]]):\n\t$L"
	    }
	    puts "Exit"
	    exit
	}
    }
    close $f
    return
}

proc checkAlamut {} {

    global g_VaRank
    global g_donnsplice

    puts "\t...checking Alamut Batch ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
    
    # If Alamut annotation is already done for all variants, no need to check Alamut
    if {$g_VaRank(skipAlamutChecks) eq "yes"} {
	puts "\t\t...skip (no need)"
	if {![info exists g_donnsplice]} {set g_donnsplice ""}
	return
    }

    if {![file exists $g_VaRank(vcfDir)/Alamut]} {file mkdir $g_VaRank(vcfDir)/Alamut}

    ## Creating an input test file for Alamut Batch
    set AlamutInputFile [file join $g_VaRank(vcfDir) Alamut test.input]
    set annFile         [file join $g_VaRank(vcfDir) Alamut test.ann]
    set unannFile       [file join $g_VaRank(vcfDir) Alamut test.unann]
    set outputFile      [file join $g_VaRank(vcfDir) Alamut test.output]

    ReplaceTextInFile "id\t11\t66293652\tT\tG" $AlamutInputFile

    ## Checking if the licence is accepted
    set test 0
    foreach f [glob -nocomplain $g_VaRank(alamutDir)/alamut-*.ini] {
	if {[file exists $f]} {
	    if {[regexp "LicenseAccepted=true" [ContentFromFile $f]]} {
		set test 1
	    }
	}
    }
    if {! $test} {
	puts "#############################################################################################################"
	puts "\t...Alamut Batch License agreement is apparently not accepted ($f)."
	puts "\t   Please, run Alamut Batch on a single example and accept the terms of licence before running VaRank."
	puts "\t   You can run the two following command lines:"
	puts "\t   cd $g_VaRank(vcfDir)/Alamut"
	puts "\t   $g_VaRank(alamutDir)/alamut-batch --assbly $g_VaRank(alamutHumanDB) --in test.input --ann test.ann --unann test.unann --alltrans --outputannonly --ssIntronicRange 2 --outputDataVersions"
	puts "\t   EXIT of VaRank!"
	puts "#############################################################################################################"
	exit
    }

    ## checking if the --donnsplice option is available (version dependant)
    catch {eval exec $g_VaRank(alamutDir)/alamut-batch} Message
    if {[regexp "donnsplice" $Message]} {
	set g_donnsplice "--donnsplice"
    } else {
	set g_donnsplice ""
    }

    ## Building the Alamut Batch command
    set alamutCmd "$g_VaRank(alamutDir)/alamut-batch --outputDataVersions --assbly $g_VaRank(alamutHumanDB) --in $AlamutInputFile --ann $annFile --unann $unannFile --outputannonly --ssIntronicRange 2"

    if {$g_VaRank(AlamutAlltrans)=="yes"} {append alamutCmd " --alltrans"}
    if {$g_VaRank(AlamutProcesses)}       {append alamutCmd " --processes $g_VaRank(AlamutProcesses)"}


    if {([info exists g_VaRank(proxyUser)] && $g_VaRank(proxyUser) != "") && ([info exists g_VaRank(proxyServer)] && $g_VaRank(proxyServer) != "")} {
	append alamutCmd " --proxyuser $g_VaRank(proxyUser) --proxypasswd $g_VaRank(proxyPasswd) --proxyserver $g_VaRank(proxyServer) --proxyport $g_VaRank(proxyPort)"
    }
    append alamutCmd " >>& $outputFile"

    if {[catch {eval exec $alamutCmd} Message]} {	
	puts "\nAlamut Batch unexpected stop\n$Message" 

	set What "Invalid parameter: --processes"
	if {[regexp $What $Message] || [regexp $What [ContentFromFile $outputFile]]} {
	    puts "################################################################################################"
	    puts " ...Alamut Batch doesn't accept \"--processes\" option (seems not to be a standalone licence)."
	    puts "    Please rerun VaRank without \"-AlamutProcesses $g_VaRank(AlamutProcesses)\" option."
	    puts "    EXIT of VaRank!"
	    puts "################################################################################################"
	    exit
	} 
	
	set What "Sorry. Access denied: key .* in use by"
	if {[regexp $What $Message] || [regexp $What [ContentFromFile $outputFile]]} {
	    puts "##############################################################"
	    puts " ...Alamut Batch access denied because key in use."
	    puts "    Please rerun VaRank when Alamut Batch will be available."
	    puts "    EXIT of VaRank!"
	    puts "##############################################################"
	    exit
	} 
	
	set What "expired license"
	if {[regexp $What $Message] || [regexp $What [ContentFromFile $outputFile]]} {	    
	    puts "##################################################################"
	    puts " ...Alamut Batch expired license!"
	    puts "    Please make sure your Alamut Batch license is still active."
	    puts "    EXIT of VaRank!"
	    puts "##################################################################"
	    exit
	} 	
	
    } else {
	#Checking Alamut Batch version compatibility both upper and lower authorized versions
	#UpperVersion
	set A 1
	set B 9
	set C 0
	
	#LowerVersion
	set A_l 1
	set B_l 6
	set C_l 0
	
	if {[regexp "alamut-.*? version (\[0-9\]+)\\.(\[0-9\]+)\\.(\[0-9\]+)" [ContentFromFile $outputFile] match vA vB vC]} {
	    set test 0
	    if {$vA > $A || $vA < $A_l} {
		set test 1
	    } elseif {$vB > $B || $vB < $B_l} {
		set test 1
	    } elseif {$vC > $C || $vC < $C_l} {
		set test 1
	    }
	    if {$test} {
		puts "#################################################################################"
		puts " ...WARNING: VaRank $g_VaRank(Version) is compatible with Alamut Batch $A_l.$B_l.$C_l up to $A.$B.$C"
		puts "    You are using Alamut Batch $vA.$vB.$vC"
		puts "#################################################################################"
	    }
	}
    }
    puts "\t   checks ok"
    puts "\t...running annotations"
}
  

## Use of the g_vcfINFOS global variable to create the Alamut input file (one file for all patients)
proc createAlamutInputFile {} {

    global g_vcfINFOS
    global g_VaRank

    set VCFdir $g_VaRank(vcfDir)
    set AlamutInputFile "$VCFdir/Alamut/AlamutInputFile_all.txt"

    if {![file exists "$VCFdir/Alamut"]} {file mkdir "$VCFdir/Alamut"}

    puts "...creation of the alamut input file ($AlamutInputFile) ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
    file delete -force $AlamutInputFile

    set nbSNV 0
    set AlamutText ""

    foreach ID [set g_vcfINFOS(L_IDs)] {
	# If ALT=*, then variant is not wrote into the input file for annotation.
	if {[regexp "_\\*$" $ID]} {continue}

	append AlamutText "$ID\t[join [lrange $g_vcfINFOS($ID) 0 3] "\t"]\n"		
	incr nbSNV
    }	
    if {[info exists g_VaRank(DEBUG)]} {puts "\t...$nbSNV variations created"}

    regsub "\n$" $AlamutText "" AlamutText
    WriteTextInFile $AlamutText $AlamutInputFile

    return
}


## Return 0 if all the input file variations have been analysed or if the input file is empty.
## Else, return 1 and create an TmpAlamutInputFile with the variations not yet analysed.
proc AlamutIsNotCompleted {inputFile annFile unannFile outputFile} {
    
    regsub "InputFile_all.txt" $inputFile "CannotLoadGene.txt" cannotloadgeneFile

    if {[file size $inputFile] == 0} {
	puts "\t...Empty input file for Alamut: $inputFile"
	return 0
    }

    regsub "InputFile_all.txt" $inputFile "InputFile_remainingID.tmp" TmpAlamutInputFile
    file delete -force $TmpAlamutInputFile

    ## Load into lIDana the ID already analysed by Alamut.
    set lIDana {}
    if {[file exists $annFile]} {
	## annFile can be a big data file
	set f [open $annFile]
	while {![eof $f]} {
	    set L [gets $f]
	    if {[regexp "^#" $L]} {continue}
	    set ID [lindex $L 0]
	    if {![info exists TabIDana($ID)]} {set TabIDana($ID) 1;lappend lIDana $ID} else {continue}
	}
	close $f	
    }

    if {[file exists $unannFile]} {
	# Line example:
	# 1:887801/1_887801_A_G/A/G       Transcript NM_015658.3: Genome/Transcript discrepancy: Alternate genomic nucleotide (G) same as transcript nucleotide (Assembly: GRCh37)

	foreach L [LinesFromFile $unannFile] {
	    set ID [lindex [split [lindex $L 0] "/"] 1]
	    if {![info exists TabIDana($ID)]} {set TabIDana($ID) 1;lappend lIDana $ID} else {continue}
	}
    } 
    if {$lIDana == ""} {
	file copy -force $inputFile $TmpAlamutInputFile
	return 1
    }

    ## Load into TabIDinvalid the ID that need to be updated by Alamut developer:
    ## Genes defined with "Invalid data while loading gene" or "Cannot load gene" must be requested by mail to the Alamut developers.
    ## These ID don't work with Alamut, we do not put them in the inputfile.
    if {[file exists $outputFile]} {
	foreach L [LinesFromFile $outputFile] {
	    #puts $L
	    if {[regexp "Cannot load gene (.*)" $L match gene]} {
		regsub ":" $L "" L
		#Room for improvment, use lsearch -index 0 or add the identifiers from lIDinvalid into a Tab like for lIDana
		#JEAN
		#lappend lIDinvalid {[lindex $L 0] $gene}
		set ID [lindex $L 0]
		if {![info exists TabIDinvalid($ID)]} {set TabIDinvalid($ID) $gene} else {continue}
	    }
	}
    }

    set test 0
    ## Some ID (With Knome for example) are too long to be treated by a lsearch:
    ## "couldn't compile regular expression pattern: out of memory"
    ## They are not analyzed and are written in $unannFile.strange
    file delete -force $unannFile.strange

    ## TmpAlamutInputFile creation:
    file delete -force $cannotloadgeneFile
    foreach L [LinesFromFile $inputFile] {
	if {$L == ""} {continue}

	set ID [lindex $L 0]
	
	if {[info exists TabIDinvalid($ID)]} {
	    WriteTextInFile "Cannot load gene [set TabIDinvalid($ID)]!!! $ID not analysed" $cannotloadgeneFile
	    #WriteTextInFile "Cannot load gene [lindex [lindex $lIDinvalid $i] end]!!! $ID not analysed" $cannotloadgeneFile
	    continue
	}
	if {![info exists TabIDana($ID)]} {
	    set test 1
	    WriteTextInFile $L $TmpAlamutInputFile
	}
    }
    
    if {[file exists $cannotloadgeneFile]} {
	WriteTextInFile "\nPlease, ask for an update for these genes to the Alamut Batch developers." $cannotloadgeneFile
    }

    return $test
}


#################################################################################################
# OUTPUT:
#    - $VCFdir/Alamut/AlamutAnnotations_all.txt
#    - $VCFdir/Alamut/AlamutUnannotated_all.txt
#    - $VCFdir/Alamut/AlamutOutput_all.txt
#
# RETURN :
#    - return "1" if Alamut has been run for all the variations
#    - else "exit" (if alamut has been run 10 times and is still not completed)
#     	 	   (or if alamut license expired)
#################################################################################################
proc runAlamut {} {

    global g_VaRank
    global g_donnsplice

    set VCFdir $g_VaRank(vcfDir)
    set AlamutInputFile "$VCFdir/Alamut/AlamutInputFile_all.txt"
    regsub "AlamutInputFile"   $AlamutInputFile "AlamutAnnotations"  annFile
    regsub "AlamutInputFile"   $AlamutInputFile "AlamutUnannotated"  unannFile
    regsub "AlamutInputFile"   $AlamutInputFile "AlamutOutput"       outputFile
    regsub "InputFile_all.txt" $AlamutInputFile "CannotLoadGene.txt" cannotloadgeneFile

    file delete -force $outputFile
    file delete -force $cannotloadgeneFile

    ## Running Alamut for the first time if the output files don't exist.
    #
    puts "...running Alamut Batch ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"

    if {[catch {checkAlamut} Message]} {
        puts "VaRank seems to have a problem while checking Alamut Batch. Exit"
        puts "######################################################################"
        puts "$Message"
        puts "######################################################################"
        exit
    }
    if {![file exists $annFile] && ![file exists $unannFile]} {
	set DateDuJour [clock seconds]
	ReplaceTextInFile "Alamut Batch started: [clock format $DateDuJour -format "%B %d %Y - %H:%M"]" $outputFile
	set doitagain 1
	while {$doitagain} {

	    #Building the Alamut Batch command
	    set alamutCmd "$g_VaRank(alamutDir)/alamut-batch --outputDataVersions --assbly $g_VaRank(alamutHumanDB) $g_donnsplice --in $AlamutInputFile --ann $annFile --unann $unannFile --outputannonly --ssIntronicRange 2"
	    if {$g_VaRank(AlamutAlltrans)=="yes"} {append alamutCmd " --alltrans"}
	    if {$g_VaRank(AlamutProcesses)}       {append alamutCmd " --processes $g_VaRank(AlamutProcesses)"}
	    if {([info exists g_VaRank(proxyUser)] && $g_VaRank(proxyUser) != "") && ([info exists g_VaRank(proxyServer)] && $g_VaRank(proxyServer) != "")} {
		append alamutCmd " --proxyuser $g_VaRank(proxyUser) --proxypasswd $g_VaRank(proxyPasswd) --proxyserver $g_VaRank(proxyServer) --proxyport $g_VaRank(proxyPort)"
	    }
	    append alamutCmd " >>& $outputFile"


	    if {[catch {eval exec $alamutCmd} Message]} {
		## Alamut has implemented a token system so that we can run only x instance(s) at once.
		## Else we get the error message "Sorry. ... Access denied."

		WriteTextInFile "\nAlamut unexpected stop\n$Message" $outputFile

		set ErrorLines ""
		if {[file exists $outputFile]} {
		    set ErrorLines [join [LinesFromFile $outputFile] " "]
		}
		#file delete -force $outputFile
		

		#JEAN NECESSARY now that we have
		set What "Sorry. Access denied: key .* in use by"
		if {[regexp $What $Message]||[regexp $What $ErrorLines]} {
		    puts "\t...Alamut Batch access denied because key in use, waiting."
		    after 600000; # on attend 10 min qu'une instance se libère
		} else {
		    set doitagain 0
		}

	    } else {
		set doitagain 0
	    }
	}
	set DateDuJour [clock seconds]
	WriteTextInFile "Alamut Batch finished: [clock format $DateDuJour -format "%B %d %Y - %H:%M"]\n\n" $outputFile
    } elseif {![AlamutIsNotCompleted $AlamutInputFile $annFile $unannFile $outputFile]} {
	puts "\t...Alamut Batch has already analysed all variations, continue ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
	return 1
    }
    
    if {![info exists g_VaRank(AlamutHeader)]} {
	#BUG HEADER PROBLEM one way to solve it

	if {[file exists $annFile]} {
	    set Header [FirstLineFromFile $annFile]
	    if {![regexp {^#} $Header]} {puts "\t...WARNING Trying to save the Alamut Batch header but already missing"} else {set g_VaRank(AlamutHeader) $Header}
	}
    }
    
    ## Rerun Alamut until all variations have been analysed (max executed 10 times)
    #
    set i 1
    set i_Limit 11

    regsub "InputFile_all.txt" $AlamutInputFile "InputFile_remainingID.tmp" TmpAlamutInputFile

    puts "\t...checking if Alamut Batch has already run on all ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
    set doitagain    1
    while {$i < $i_Limit && [AlamutIsNotCompleted $AlamutInputFile $annFile $unannFile $outputFile] && $doitagain} {
	puts "\t...Alamut Batch didn't finish to analyse all variations. New run is launched (Try $i) ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"

	set DateDuJour [clock seconds]
	WriteTextInFile "Alamut Batch started: [clock format $DateDuJour -format "%B %d %Y - %H:%M"]" $outputFile

	#set j 0
	#set j_Limit 3
	#while {$doitagain && $j<$j_Limit} {}
	while {$doitagain} {

	    #Building the Alamut Batch command
	    set alamutCmd "$g_VaRank(alamutDir)/alamut-batch --outputDataVersions --assbly $g_VaRank(alamutHumanDB) $g_donnsplice --in $TmpAlamutInputFile --ann $annFile.tmp --unann $unannFile.tmp --outputannonly --ssIntronicRange 2"

	    if {$g_VaRank(AlamutAlltrans)=="yes"} {append alamutCmd " --alltrans"}
	    if {$g_VaRank(AlamutProcesses)}       {append alamutCmd " --processes $g_VaRank(AlamutProcesses)"}
	    
	    if {([info exists g_VaRank(proxyUser)] && $g_VaRank(proxyUser) != "") && ([info exists g_VaRank(proxyServer)] && $g_VaRank(proxyServer) != "")} {
		append alamutCmd " --proxyuser $g_VaRank(proxyUser) --proxypasswd $g_VaRank(proxyPasswd) --proxyserver $g_VaRank(proxyServer) --proxyport $g_VaRank(proxyPort)"
	    }
	    append alamutCmd " >>& $outputFile"
	    WriteTextInFile "\n$alamutCmd\n" $outputFile
	    if {[catch {eval exec $alamutCmd} Message]} {
		#set Message "Sorry. Access denied: key 1321321 in use by"
		WriteTextInFile "\nAlamut Batch unexpected stop\n$Message" $outputFile

		set ErrorLines ""
		if {[file exists $outputFile]} {
		    set ErrorLines [join [LinesFromFile $outputFile] " "]
		}
		
		set What "expired license"
		if {[regexp $What $Message]||[regexp $What $ErrorLines]} {
		    #Sorry. Access denied: expired license

		    puts "################################################################"
		    puts " ...Alamut Batch expired license!"
		    puts "    Please make sure your Alamut Batch license is still active."
		    puts "    EXIT of VaRank!"
		    puts "################################################################"
		    exit
		} 

		set What "child killed: segmentation violation"
		if {[regexp $What $Message]||[regexp $What $ErrorLines]} {
		    #Some ID caused "segmentation violation" in the alamut program only on debian. Not successfully debbugged by alamut.

		    puts "###############################################"
		    puts " ...Alamut Batch segmentation violation"
		    puts "###############################################"
		    
		    # Keep annotations already done 
		    if {[file exists $annFile.tmp] && [file size $annFile.tmp] != 0} {
			# several header lines since adding "--outputDataVersions" in the command line. 
			set f [open $annFile.tmp]
			set nb 0
			set infosAnn {}
			while {![eof $f]} {
			    set L [gets $f]
			    if {[regexp "^#" $L]} {continue}
			    incr nb
			    if {$nb > 100000} {
				WriteTextInFile [join $infosAnn "\n"] $annFile
				set nb 0
				set infosAnn {}
			    }
			    lappend infosAnn $L
			}
			close $f	
			if {$infosAnn != ""} {
			    WriteTextInFile [join $infosAnn "\n"] $annFile
			}
		    }
		    file delete -force $annFile.tmp

		    if {[file exists $unannFile.tmp] && [file size $unannFile.tmp] != 0} {
			WriteTextInFile [join [LinesFromFile $unannFile.tmp] "\n"] $unannFile
		    }
		    file delete -force $unannFile.tmp

		    ## Add the corruptedID in the unannFile
		    set ID ""
		    set corruptedID ""
		    set lastIDAnalysed ""
		    foreach L [LinesFromFile $outputFile] {
			set ID [lindex $L 0]
			if {[regexp "^\[0-9XY\]_\[0-9\]+_" $ID]} {set lastIDAnalysed $ID}
		    }
		    regsub ":" $lastIDAnalysed "" lastIDAnalysed
		    if {$lastIDAnalysed != ""} {
			set test 0
			foreach L [LinesFromFile $TmpAlamutInputFile] {
			    if {$test} {set corruptedID [lindex $L 0]; break}
			    if {[regexp $lastIDAnalysed $L]} {set test 1}			
			}
		    } else {
			# This is the first ID of TmpAlamutInputFile that is corrupted
			set corruptedID [lindex [lindex [LinesFromFile $TmpAlamutInputFile] 0] 0]
		    }
		    WriteTextInFile "$corruptedID child killed: segmentation violation" $unannFile
		    puts "$corruptedID : child killed: segmentation violation"
		    puts "Must be requested by mail to the Alamut developers."
		    		    
		    ## Should be deleted to permit the regexp into "$ErrorLines"
		    file delete -force $outputFile

		    # Waiting 10 minutes before reruning and expecting a free token
		    after 600000

		    AlamutIsNotCompleted $AlamutInputFile $annFile $unannFile $outputFile
		    continue
		} 

		set What "Sorry. Access denied: key .* in use by"
		if {[regexp $What $Message]||[regexp $What $ErrorLines]} {
		    set DateDuJour [clock seconds]
		    puts "\t\t...Alamut Batch licence in use, waiting some time before retry [clock format $DateDuJour -format "%B %d %Y - %H:%M"]"
		    WriteTextInFile "... Licence in use waiting some time before retry [clock format $DateDuJour -format "%B %d %Y - %H:%M"]" $outputFile

		    #Waiting 10 minutes before reruning and expecting a free token
		    after 600000

		    if {[file exists $annFile.tmp] && [file size $annFile.tmp] != 0} {
			set  infosAnn [lrange [LinesFromFile $annFile.tmp] 1 end]
			if {$infosAnn != ""} {
			    WriteTextInFile [join $infosAnn "\n"] $annFile
			}
		    }
		    if {[file exists $unannFile.tmp] && [file size $unannFile.tmp] != 0} {
			WriteTextInFile [join [LinesFromFile $unannFile.tmp] "\n"] $unannFile
		    }
		    
		    file delete -force $annFile.tmp
		    file delete -force $unannFile.tmp

		    ## Should be deleted to permit the regexp into "$ErrorLines"
		    #WriteTextInFile [join [LinesFromFile $outputFile] "\n"] $outputFile.$i
		    file delete -force $outputFile
		    
		    #We keep unlimited because it can be very long before
		    #incr i
		} 

		# Alamut error message when the connection is lost:
		#     ***BATCH FINISHED WITH ERRORS***
		#     Couldn't cleanly close connection
		#     Alamut Batch unexpected stop
		#     child process exited abnormally
		set What "***BATCH FINISHED WITH ERRORS***"
		if {[regexp $What $Message]||[regexp $What $ErrorLines]} {
		    set DateDuJour [clock seconds]
		    puts "\t\t...Alamut Batch finished with errors, waiting some time before retry [clock format $DateDuJour -format "%B %d %Y - %H:%M"]"
		    WriteTextInFile "... Alamut Batch finished with errors, waiting some time before retry [clock format $DateDuJour -format "%B %d %Y - %H:%M"]" $outputFile

		    # Waiting 10 minutes before reruning and expecting a free token
		    after 600000

		    if {[file exists $annFile.tmp] && [file size $annFile.tmp] != 0} {
			set  infosAnn [lrange [LinesFromFile $annFile.tmp] 1 end]
			if {$infosAnn != ""} {
			    WriteTextInFile [join $infosAnn "\n"] $annFile
			}
		    }
		    if {[file exists $unannFile.tmp] && [file size $unannFile.tmp] != 0} {
			WriteTextInFile [join [LinesFromFile $unannFile.tmp] "\n"] $unannFile
		    }
		    
		    file delete -force $annFile.tmp
		    file delete -force $unannFile.tmp

		    ## Should be deleted to permit the regexp into "$ErrorLines"
		    #WriteTextInFile [join [LinesFromFile $outputFile] "\n"] $outputFile.$i
		    file delete -force $outputFile
		    
		    # We keep unlimited because it can be very long before to retrieve connection
		    # incr i
		} else {
		    set doitagain 0
		}
		
	    } else {
		WriteTextInFile "\n$Message" $outputFile
		set doitagain 0
	    }
	}
	set DateDuJour [clock seconds]

	incr i
	if {[file exists $annFile.tmp] && [file size $annFile.tmp] != 0} {
	    set f [open $annFile.tmp]
	    set nb 0
	    set texteL {}
    	    while {![eof $f]} {
                set L [gets $f]
		if {[regexp "^#" $L]} {continue}
		incr nb
		if {$nb > 100000} {
			WriteTextInFile [join $texteL "\n"] $annFile
			set nb 0
			set texteL {}
		}
		lappend texteL $L
	    }
	    close $f
	    WriteTextInFile [join $texteL "\n"] $annFile
	}

	if {[file exists $unannFile.tmp] && [file size $unannFile.tmp] != 0} {
	    WriteTextInFile [join [LinesFromFile $unannFile.tmp] "\n"] $unannFile
	}
	file delete -force $annFile.tmp
	file delete -force $unannFile.tmp
    }

    if {$i > $i_Limit || [AlamutIsNotCompleted $AlamutInputFile $annFile $unannFile $outputFile]} {
	# Alamut did not run on all variations
	if {[file exists $cannotloadgeneFile]} {
	    puts "[ContentFromFile $cannotloadgeneFile]\n"
	}

	puts "###########################################################"
	puts "   WARNING: Alamut Batch is not finished after $i run!"
	puts "   EXIT of VaRank!"
	puts "###########################################################"
	exit
    }

    puts "\t...Alamut Batch finished: [clock format $DateDuJour -format "%B %d %Y - %H:%M"]"
    WriteTextInFile "Alamut Batch finished: [clock format $DateDuJour -format "%B %d %Y - %H:%M"]\n\n" $outputFile
    
    ## Alamut has been run for all the variations
    if {[file exists $cannotloadgeneFile]} {
	puts "[ContentFromFile $cannotloadgeneFile]\n"
    }
    
    return 1
} 


proc parseAlamutFile {} {

    global g_ANNOTATION
    global g_VaRank
    global g_vcfINFOS

    set VCFdir $g_VaRank(vcfDir)
    set annFile "$VCFdir/Alamut/AlamutAnnotations_all.txt"

    puts "...parsing Alamut Batch results ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"

    set f [open $annFile]
    set g_vcfINFOS(#id) ""; # not to have a bug! Unset just below after the while.

    # Search for the alamut header line
    while {![eof $f]} {
        set L [gets $f]
	if {[regexp "^#id" $L]} {set headerAlamutLine [split $L "\t"]; break} 
    }
    close $f

    # Search for colName from alamut to keep in memory (no need to keep useless annotation in the RAM!)
    set L_alamutColNameToKeep {}
    foreach colName $g_VaRank(L_outputColHeader) {
	if {[lsearch -exact "$headerAlamutLine" $colName] ne -1} {
	    lappend L_alamutColNameToKeep "$colName"
	}
    }
    set g_ANNOTATION(#id) "[join $L_alamutColNameToKeep "\t"]"

    # parsing
    set f [open $annFile]
    while {![eof $f]} {
        set L [gets $f]
	if {[regexp "^#" $L]} {continue}

	# Load only variants from the patients analyzed
	set ID [lindex $L 0]
	if {![info exists g_vcfINFOS($ID)]} {continue}

	# Replacing the empty columns by NA
	# Note regsub is faster
	regsub -all "\t\t" $L "\tNA\t" L
	regsub -all "\t\t" $L "\tNA\t" L; # Same command used twice because of the possible succession of 3 "\t" ("\t\t\t") 
	regsub "\t$" $L "\tNA" L

	set Ls [split $L "\t"]
	eval lassign \$Ls $headerAlamutLine

	## ID corresponds to 1 genetic variant.
	## There may be multiple transcripts per genetic variant (strand +, strand -, different splicing isoforms).
	## This need to be handled specifically later to ensure coherence of the specific data for each isoform

	set newL "" 
	foreach col $L_alamutColNameToKeep {
	    lappend newL [set $col]
	}
	lappend g_ANNOTATION($ID) "[join $newL "\t"]"
    }
    close $f
    unset g_vcfINFOS(#id)

    # Analyzing the alamut header
    if {![info exists g_ANNOTATION(#id)]} {
        if {[info exists g_VaRank(AlamutHeader)] && [set g_VaRank(AlamutHeader)]!=""} {
            lappend g_ANNOTATION(#id) "$L"
            puts "\t...WARNING Alamut Batch output file seems to miss the header. Rescue done!"
        } else {
            puts "\t...WARNING Alamut Batch output file seems to miss the header. Exit of VaRank!"
            exit
        }
    }

    return 
}



