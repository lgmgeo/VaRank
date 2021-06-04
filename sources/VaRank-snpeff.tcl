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

proc createSnpEffInputFile {} {

    global g_VaRank

    ## Keep SnpEff input files already existing
    ## No deletion

    ## Removing snpEff input files not already annotated
    foreach inputFile [glob -nocomplain $g_VaRank(vcfDir)/SnpEff/InputDir/*.vcf] {
	regsub "InputDir" $inputFile "OutputDir" annFile
	regsub ".vcf" $annFile ".varType.dbsnp.dbnsfp.phastCons.vcf" annFile
	if {![file exists $annFile]} {file delete -force $inputFile}
    }

    ## Creation of the SnpEff input files (containing non redundant variants)
    #
    puts "...creation of the SnpEff input files ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
    file mkdir $g_VaRank(vcfDir)/SnpEff
    file mkdir $g_VaRank(vcfDir)/SnpEff/Input

    ## If some SnpEff outputFiles already exists, their ID should not be reported in the input file (because already annotated)
    ## Load into lIDana the ID already analysed by SnpEff
    set lIDana {}
    foreach annFile [glob -nocomplain $g_VaRank(vcfDir)/SnpEff/Output/*.varType.dbsnp.dbnsfp.phastCons.vcf] {
	foreach L [LinesFromFile $annFile] {
	    if {[regexp "^#" $L]} {continue}
	    lappend lIDana "[lindex $L 0]_[lindex $L 1]_[lindex $L 3]_[lindex $L 4]"		
	}
    }	

    set i 1
    # Warning: Check the incrementation ($i) for inputfile names (should not be the same as already existed output file)
    while {[file exists $g_VaRank(vcfDir)/SnpEff/Output/allinone.$i.eff.varType.dbsnp.dbnsfp.phastCons.vcf]} {
	incr i
    }
    set nbL 0
    set snpeffInputFile $g_VaRank(vcfDir)/SnpEff/Input/allinone.$i.vcf
    ReplaceTextInFile "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttutu" $snpeffInputFile

    set lVariantLines {}
    set allIDfromVCFs {}
    foreach vcfFile [glob -nocomplain $g_VaRank(vcfDir)/*.vcf $g_VaRank(vcfDir)/*.vcf.gz] {
	puts "\t...parsing $vcfFile"

        if {[regexp ".vcf.gz$" $vcfFile]} {
            set lLines [LinesFromGZFile $vcfFile]
        } else {
            set lLines [LinesFromFile $vcfFile]
        }
        foreach L $lLines {
            if {[regexp "^##" $L]} {continue}
            if {[regexp "^#CHROM" $L]} {
                set Ls [split $L "\t"]
                set i_chr [lsearch -exact $Ls "#CHROM" ]; if {$i_chr != 0} {puts "Bad header line syntax: $L\nExit"; exit}
                set i_pos [lsearch -exact $Ls "POS"];     if {$i_pos != 1} {puts "Bad header line syntax: $L\nExit"; exit}
                set i_id  [lsearch -exact $Ls "ID"];      if {$i_id  != 2} {puts "Bad header line syntax: $L\nExit"; exit}
                set i_ref [lsearch -exact $Ls "REF"];     if {$i_ref != 3} {puts "Bad header line syntax: $L\nExit"; exit}
                set i_alt [lsearch -exact $Ls "ALT"];     if {$i_alt != 4} {puts "Bad header line syntax: $L\nExit"; exit}
             
		continue
            }

            if {$L == ""} {continue}

            ## If there is only wild type hom in the line (GT = 0|0 or 0/0), so there isn't mutation to analyse by VaRank
            if {![regexp "1/0|0/1|2/0|0/2|1/1|2/2|2/1|1/2" $L] && ![regexp "1\\\|0|0\\\|1|2\\\|0|0\\\|2|1\\\|1|2\\\|2|2\\\|1|1\\\|2" $L]} {continue}

            set Ls [split $L "\t"]
            set chrom [lindex $Ls 0] 
            regsub -all "chr" $chrom "" chrom

            set pos [lindex $Ls 1]
            set ref [lindex $Ls 3]
            set alt [lindex $Ls 4] 

	    foreach altn [split $alt ","] {
		if {![regexp "^\[ACGTN*\]+$" $altn]} { 
		    puts "############################################################################"
		    puts "\t   WARNING: illegal character in VCF's ALT fields ('$altn')"
		    puts "\t   Please, check your $vcfFile file."
		    puts "\t   EXIT of VaRank!"
		    puts "############################################################################"
		    puts ""
		    exit
		}
		
		# If ALT=*, then variant is not wrote into the input file for annotation.
		if {$altn == "*"} {continue}

		## ID already annotated
		if {[lsearch -exact $lIDana "${chrom}_${pos}_${ref}_${altn}"] != -1} {continue}
           	set ID [AssignID $chrom $pos $ref $altn]
	    	if {[lsearch -exact $allIDfromVCFs $ID] == -1} {
		    lappend allIDfromVCFs $ID
		    ## Remove genotype like 1/2, 0/3, 0|2 (2 doesn't exist anymore)
		    ## All genotypes are set to 0/1
		    set Ln [concat $chrom [lrange $Ls 1 3] $altn "30\tPASS\tAC=1;DP=30\tGT:AD:GQ:DP\t0/1:20,10:40:30"]
		    set Ln [join $Ln "\t"]
		    lappend lVariantLines $Ln
		    incr nbL
		    if {$nbL > 40000} {
			set lVariantLines [lsort -command AscendingSortOnElement1 $lVariantLines]
			set lVariantLines [lsort -command AscendingSortOnElement0 $lVariantLines]
			WriteTextInFile [join $lVariantLines "\n"] $snpeffInputFile
			incr i
			while {[file exists $g_VaRank(vcfDir)/SnpEff/Output/allinone.$i.eff.varType.dbsnp.dbnsfp.phastCons.vcf]} {
			    incr i
			}
			set nbL 0
			set snpeffInputFile $g_VaRank(vcfDir)/SnpEff/Input/allinone.$i.vcf
			ReplaceTextInFile "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttutu" $snpeffInputFile			
			set lVariantLines {}
		    }
	    	} 
	    }
	}
    }
    set lVariantLines [lsort -command AscendingSortOnElement1 $lVariantLines]
    set lVariantLines [lsort -command AscendingSortOnElement0 $lVariantLines]

    if {[llength $lVariantLines] != 0} {
	WriteTextInFile [join $lVariantLines "\n"] $snpeffInputFile
    } else {
	file delete -force $snpeffInputFile
    }
}

proc ErrorWithinSnpEffFile {vcfOutputFile} {

    set test 1

    regexp ".(eff)(.varType)?(.dbsnp)?(.dbnsfp)?(.phastCons)?(.msigDb)?.vcf$" $vcfOutputFile match eff varType dbsnp dbnsfp phastCons msigDb
    
    set lTags {"##INFO=<ID=ANN"}
    if {$varType != ""} {lappend lTags "##INFO=<ID=VARTYPE"}
    if {$dbsnp != ""} {lappend lTags "##INFO=<ID=CAF"}
    if {$dbnsfp != ""} {lappend lTags "##INFO=<ID=dbNSFP"}
    if {$phastCons != ""} {lappend lTags "##INFO=<ID=PhastCons"}
    if {$msigDb != ""} {lappend lTags "##INFO=<ID=MSigDb"}
    
    # vcfOutputFile should exists
    if {![file exists $vcfOutputFile]} {return $test}

    # logfile should exists, even if empty
    if {![file exists $vcfOutputFile.log]} {return $test}

    # Check java logfile
    foreach L [LinesFromFile $vcfOutputFile.log] {
	if {[regexp "^Exception in thread|java.lang." $L]} {return $test}
    }

    # Check the vcf output file
    set lReads {}
    foreach L [LinesFromFile $vcfOutputFile] {
	if {![regexp "^#" $L]} {break}
	foreach tag $lTags {
	    if {[regexp "^$tag" $L]} {
		lappend lReads "$tag"
		break
	    }
	}
    }
    set lReads [lsort -unique $lReads]
    if {[llength $lReads] == [llength $lTags]} {set test 0}

    return $test
}

proc runOneSnpEffCommand {vcfOutputFile effCommand} {

    ReplaceTextInFile "SnpEff started: [clock format [clock seconds] -format "%B %d %Y - %H:%M"]\n" $vcfOutputFile.log
    WriteTextInFile "$effCommand\n" $vcfOutputFile.log

    if {[catch {eval exec $effCommand} Message]} {
	WriteTextInFile "$Message" $vcfOutputFile.log
	WriteTextInFile "SnpEff finished: [clock format [clock seconds] -format "%B %d %Y - %H:%M"]\n" $vcfOutputFile.log
	if {[ErrorWithinSnpEffFile $vcfOutputFile]} {
	    puts "############################################################################"
	    puts "\t   SnpEff error!"
	    puts "\t   Please, check the $vcfOutputFile.log file."
	    puts "\t   EXIT of VaRank!"
	    puts "############################################################################"
	    exit
	}
    }
}

proc checkSnpEff {} {

    global g_VaRank

    puts "...checking SnpEff ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"    

    # VCF input file for testing
    set vcfFile "$g_VaRank(SnpEffTestsDir)/example.vcf"

    if {![file exists $vcfFile]} {
	ReplaceTextInFile "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1" $vcfFile
	WriteTextInFile "1\t69511\t.\tA\tG\t136\tPASS\tAC=2;DP=11;FREQ_HTZ=1.14;FREQ_HOM=34.66;DB;MG;EVS\tGT:AD:GQ:DP\t1/1:0,11:30:11" $vcfFile
	WriteTextInFile "1\t762273\t.\tG\tA\t36\tPASS\tAC=1;DP=28;FREQ_HTZ=11.93;FREQ_HOM=26.7;DB;MG\tGT:AD:GQ:DP\t0/1:19,9:36:28" $vcfFile
	WriteTextInFile "4\t1795039\t.\tG\tA\t35\tPASS\tAC=1;DP=28;FREQ_HTZ=11.93;FREQ_HOM=26.7;DB;MG\tGT:AD:GQ:DP\t0/1:19,9:36:28\t0/0:19,9:36:19" $vcfFile
	WriteTextInFile "19\t13010643\t.\tG\tT\t33\tPASS\tAC=1;DP=28;FREQ_HTZ=11.93;FREQ_HOM=26.7;DB;MG\tGT:AD:GQ:DP\t0/1:19,9:36:28\t0/1:17,11:28:23" $vcfFile
    }

    # checking
    set effCommand "$g_VaRank(javaSnpEffCommand)/snpEff.jar"
    if {[catch {eval exec $effCommand} Message]} {
	if {[regexp "SnpEff version SnpEff (\[0-9\]+)\\.(\[0-9\]+)" $Message match version subversion]} {
	    if {$version < 4 || "$version.$subversion" == "4.0"} {
		puts "############################################################################"
		puts "\t   VaRank is compatible with SnpEff from version 4.1 "
		puts "\t   You are using version $version.$subversion ($g_VaRank(snpeffDir)/snpEff.jar)"
		puts "\t   Please install a new version of SnpEff"
		puts "\t   EXIT of VaRank!"
		puts "############################################################################"
		exit
	    }
	} else {
	    puts "############################################################################"
	    puts "\t   SnpEff error!"
	    puts "\t   $Message"
	    puts "\t   Please, check your SnpEff installation"
	    puts "\t   EXIT of VaRank!"
	    puts "############################################################################"
	    exit
	}
    } 
    puts "\t...using SnpEff version $version.$subversion"

    # If SnpEff checking has already been done previously, continue
    regsub ".vcf$" $vcfFile ".eff.varType.dbsnp.dbnsfp.phastCons.vcf" vcfOutputFile
    if {[file exists $vcfOutputFile] && ([file mtime $vcfOutputFile] > [file mtime $g_VaRank(snpeffDir)/snpEff.jar])} {
	return
    }

    # Remove outputfile if they already exist
    foreach outputFile [glob -nocomplain $g_VaRank(SnpEffTestsDir)/example.eff.*] {
	file delete -force $outputFile
    }

    puts "\t...checking step1: snpEff.jar eff annotation ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
    regsub "vcf" $vcfFile "eff.vcf" vcfOutputFile
    set effCommand "$g_VaRank(javaSnpEffCommand)/snpEff.jar eff -c $g_VaRank(snpeffDir)/snpEff.config -v $g_VaRank(snpeffHumanDB) $vcfFile > $vcfOutputFile"
    runOneSnpEffCommand $vcfOutputFile $effCommand 
    set vcfFile $vcfOutputFile
    

    puts "\t...checking step2: SnpSift.jar varType annotation ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
    regsub "vcf" $vcfFile "varType.vcf" vcfOutputFile
    set effCommand "$g_VaRank(javaSnpEffCommand)/SnpSift.jar varType $vcfFile > $vcfOutputFile"
    runOneSnpEffCommand $vcfOutputFile $effCommand 
    set vcfFile $vcfOutputFile

    puts "\t...checking step3: SnpSift.jar dbSNP annotation ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
    regsub "vcf" $vcfFile "dbsnp.vcf" vcfOutputFile
    set effCommand "$g_VaRank(javaSnpEffCommand)/SnpSift.jar annotate $g_VaRank(dbSNP) $vcfFile > $vcfOutputFile"
    runOneSnpEffCommand $vcfOutputFile $effCommand 
    set vcfFile $vcfOutputFile

    puts "\t...checking step4: SnpSift.jar dbnsfp annotation ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
    regsub "vcf" $vcfFile "dbnsfp.vcf" vcfOutputFile
    set effCommand "$g_VaRank(javaSnpEffCommand)/SnpSift.jar dbnsfp -db $g_VaRank(dbNSFP) $vcfFile > $vcfOutputFile"
    runOneSnpEffCommand $vcfOutputFile $effCommand 
    set vcfFile $vcfOutputFile

    puts "\t...checking step5: SnpSift.jar phastCons annotation ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
    regsub "vcf" $vcfFile "phastCons.vcf" vcfOutputFile
    set effCommand "$g_VaRank(javaSnpEffCommand)/SnpSift.jar phastCons $g_VaRank(phastConsDB) $vcfFile > $vcfOutputFile"
    runOneSnpEffCommand $vcfOutputFile $effCommand 
    set vcfFile $vcfOutputFile

    #    puts "\t...checking step6: SnpSift.jar geneSets annotation ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
    #    regsub "vcf" $vcfFile "msigDb.vcf" vcfOutputFile
    #    set effCommand "$g_VaRank(javaSnpEffCommand)/SnpSift.jar geneSets -v $g_VaRank(msigdb) $vcfFile > $vcfOutputFile"
    #    runOneSnpEffCommand $vcfOutputFile $effCommand 
    #    set vcfFile $vcfOutputFile

    file delete -force snpEff_genes.txt 
    file delete -force snpEff_summary.html 

    return
}

#################################################################################################
# OUTPUTS:
#    - $g_VaRank(vcfDir)/SnpEff/'vcfFile'.eff.varType.dbsnp.dbnsfp.phastCons.vcf
#    (- $g_VaRank(vcfDir)/SnpEff/'vcfFile'.eff.varType.dbsnp.dbnsfp.phastCons.msigDb.vcf)
#
# RETURN :
#    - return "1" if SnpEff has been run for all the variations
#    - else "exit" 
#################################################################################################
proc runSnpEff {} {

    global g_VaRank

    ## Running SnpEff 
    #
    puts "...running SnpEff ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
    file mkdir $g_VaRank(vcfDir)/SnpEff/Output
    foreach vcfFile [glob -nocomplain $g_VaRank(vcfDir)/SnpEff/Input/*.vcf $g_VaRank(vcfDir)/SnpEff/Input/*.vcf.gz] {
	puts "\t...Annotation of $vcfFile"
	
	## Checking if all annotations are already done 
	#regsub "vcf(.gz)?" $vcfFile "eff.varType.dbsnp.dbnsfp.phastCons.msigDb.vcf" vcfOutputFile
	regsub "vcf(.gz)?$" $vcfFile "eff.varType.dbsnp.dbnsfp.phastCons.vcf" vcfOutputFile	
	regsub "/Input/" $vcfOutputFile "/Output/" vcfOutputFile

	if {[ErrorWithinSnpEffFile $vcfOutputFile]} {

	    puts "\tstep1: snpEff.jar eff annotation ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
	    regsub "vcf(.gz)?$" $vcfFile "eff.vcf" vcfOutputFile
	    regsub "/Input/" $vcfOutputFile "/Output/" vcfOutputFile
	    set effCommand "$g_VaRank(javaSnpEffCommand)/snpEff.jar eff -c $g_VaRank(snpeffDir)/snpEff.config -v $g_VaRank(snpeffHumanDB) $vcfFile > $vcfOutputFile"
	    runOneSnpEffCommand $vcfOutputFile $effCommand 
	    set vcfFile $vcfOutputFile

	    puts "\tstep2: SnpSift.jar varType annotation ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
	    regsub "vcf$" $vcfFile "varType.vcf" vcfOutputFile
	    set effCommand "$g_VaRank(javaSnpEffCommand)/SnpSift.jar varType $vcfFile > $vcfOutputFile"
	    runOneSnpEffCommand $vcfOutputFile $effCommand 
	    file delete -force $vcfFile 
	    set vcfFile $vcfOutputFile

	    puts "\tstep3: SnpSift.jar dbSNP annotation ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
	    regsub "vcf$" $vcfFile "dbsnp.vcf" vcfOutputFile
	    set effCommand "$g_VaRank(javaSnpEffCommand)/SnpSift.jar annotate $g_VaRank(dbSNP) $vcfFile > $vcfOutputFile"
	    runOneSnpEffCommand $vcfOutputFile $effCommand 
	    file delete -force $vcfFile 
	    set vcfFile $vcfOutputFile

	    puts "\tstep4: SnpSift.jar dbnsfp annotation ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
	    regsub "vcf$" $vcfFile "dbnsfp.vcf" vcfOutputFile
	    set effCommand "$g_VaRank(javaSnpEffCommand)/SnpSift.jar dbnsfp -db $g_VaRank(dbNSFP) $vcfFile > $vcfOutputFile"
	    runOneSnpEffCommand $vcfOutputFile $effCommand 
	    file delete -force $vcfFile 
	    set vcfFile $vcfOutputFile

	    puts "\tstep5: SnpSift.jar phastCons annotation ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
	    regsub "vcf$" $vcfFile "phastCons.vcf" vcfOutputFile
	    set effCommand "$g_VaRank(javaSnpEffCommand)/SnpSift.jar phastCons $g_VaRank(phastConsDB) $vcfFile > $vcfOutputFile"
	    runOneSnpEffCommand $vcfOutputFile $effCommand 
	    file delete -force $vcfFile 
	    set vcfFile $vcfOutputFile

	    #	    puts "\tstep6: SnpSift.jar geneSets annotation ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
	    #	    regsub "vcf$" $vcfFile "msigDb.vcf" vcfOutputFile
	    #	    set effCommand "$g_VaRank(javaSnpEffCommand)/SnpSift.jar geneSets -v $g_VaRank(msigdb) $vcfFile > $vcfOutputFile"
	    #	    runOneSnpEffCommand $vcfOutputFile $effCommand 
	    #	    file delete -force $vcfFile 
	    #	    set vcfFile $vcfOutputFile

	    file delete -force snpEff_genes.txt
	    file delete -force snpEff_summary.html

	} else {
	    puts "\t   already done."
	}
    }

    return 1
}




   
# Colum names existing with Alamut, but set to "NA" with SnpEff
###############################################################
# protein (refseq ID)
# transcript
# transLen
# rsClinicalSignificance
# clinVarClinSignifs
# SIFTmedian
# distNearestSS
# nearestSSType
# localSpliceEffect
# wtSSFScore
# wtMaxEntScore
# wtNNSScore
# varSSFScore
# varMaxEntScore
# varNNSScore                       

                                                
proc parseSnpEffFile {} {

    global g_ANNOTATION
    global g_VaRank

    ## Determine the header of g_ANNOTATION() (= g_ANNOTATION(#id))
    ###############################################################
    # Here, we report column names from "SnpEff Annotation" that we will use and report in VaRank outputs.
    # 
    # Column order is not important.
    set header {gene rsID Uniprot codingEffect varLocation exon intron varType Annotation_Impact Gene_ID Feature_Type Feature_ID Transcript_BioType cNomen pNomen wtAA_1 posAA varAA_1 cDNA.pos cDNA.length CDS.pos CDS.length AA.pos AA.length Distance dbNSFP_1000Gp1_AF dbNSFP_1000Gp1_AFR_AF dbNSFP_1000Gp1_AMR_AF dbNSFP_1000Gp1_ASN_AF dbNSFP_1000Gp1_EUR_AF dbNSFP_CADD_phred dbNSFP_ExAC_AC dbNSFP_ExAC_AF dbNSFP_ExAC_AFR_AC dbNSFP_ExAC_AFR_AF dbNSFP_ExAC_AMR_AC dbNSFP_ExAC_AMR_AF dbNSFP_ExAC_Adj_AC dbNSFP_ExAC_Adj_AF dbNSFP_ExAC_EAS_AC dbNSFP_ExAC_EAS_AF dbNSFP_ExAC_FIN_AC dbNSFP_ExAC_FIN_AF dbNSFP_ExAC_NFE_AC dbNSFP_ExAC_NFE_AF dbNSFP_ExAC_SAS_AC dbNSFP_ExAC_SAS_AF dbNSFP_FATHMM_pred dbNSFP_GERP___NR dbNSFP_GERP___RS dbNSFP_LRT_pred dbNSFP_MetaSVM_pred dbNSFP_MutationAssessor_pred dbNSFP_MutationTaster_pred dbNSFP_PROVEAN_pred dbNSFP_Polyphen2_HDIV_pred dbNSFP_Polyphen2_HVAR_pred dbNSFP_SIFT_pred phastCons SIFTprediction PPH2class LOF NMD}

    set g_ANNOTATION(#id) [join $header "\t"]


    puts "...parsing SnpEff annotations ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
    #foreach annFile [glob -nocomplain $g_VaRank(vcfDir)/SnpEff/*.eff.varType.dbsnp.dbnsfp.phastCons.msigDb.vcf] 
    foreach annFile [glob -nocomplain $g_VaRank(vcfDir)/SnpEff/Output/*.eff.varType.dbsnp.dbnsfp.phastCons.vcf] {
	puts "\t$annFile"
	foreach L [LinesFromFile $annFile] {
	    if {[regexp "^##" $L]} {continue}
	    if {[regexp "^#CHROM" $L]} {
		set Ls [split $L "\t"]
		set i_chr    [lsearch -exact $Ls "#CHROM"];  if {$i_chr    == -1} {puts "Bad header line syntax into $vcfFile. \n#CHROM column not found - Exit"; exit}
		set i_pos    [lsearch -exact $Ls "POS"];     if {$i_pos    == -1} {puts "Bad header line syntax into $vcfFile. \nPOS column not found - Exit"; exit}
		set i_ref    [lsearch -exact $Ls "REF"];     if {$i_ref    == -1} {puts "Bad header line syntax into $vcfFile. \nREF column not found - Exit"; exit}
		set i_alt    [lsearch -exact $Ls "ALT"];     if {$i_alt    == -1} {puts "Bad header line syntax into $vcfFile. \nALT column not found - Exit"; exit}
		set i_id     [lsearch -exact $Ls "ID"];      if {$i_id     == -1} {puts "Bad header line syntax into $vcfFile. \nID column not found - Exit"; exit}
		set i_info   [lsearch -exact $Ls "INFO"];    if {$i_info   == -1} {puts "Bad header line syntax into $vcfFile. \nINFO column not found - Exit"; exit}
		continue
	    }

	    if {$L == ""} {continue}

	    ## If there is only normal hom in the line (GT = 0|0 or 0/0), so there isn't mutation to analyse by VaRank
	    if {![regexp "1/0|0/1|2/0|0/2|1/1|2/2|2/1|1/2" $L] && ![regexp "1\\\|0|0\\\|1|2\\\|0|0\\\|2|1\\\|1|2\\\|2|2\\\|1|1\\\|2" $L]} {continue}

	    ## Variables initialisation:
	    ############################
	    #  - variables from SnpEff (not "ANN") that are not in the header:
	    set dbNSFP_Uniprot_acc         "NA"
	    set VARTYPE                    "NA"
	    set PhastCons                  "NA"
	    # - variables from the header
	    foreach col $header {
		set $col "NA"
	    }

	    ## Parsing 
	    ##########
	    set Ls [split $L "\t"]
	    set chrom [lindex $Ls $i_chr] 
	    regsub -all "^chr" $chrom "" chrom

	    set pos [lindex $Ls $i_pos] 
	    set ref [lindex $Ls $i_ref] 
	    set alt [lindex $Ls $i_alt] 
	    set rsID [lindex $Ls $i_id] 
	    if {[isNotAnRS $rsID]} {set rsID "NA"} 


	    ## ALT=T --> ANN=T|...
	    ## ALT=G,T --> ANN=G|... , T|...
	    ## ALT=C,G (cancer) --> ANN=G-C|... (format = ALT-REFERENCE)
	    ## ALT=G (into MNP)--> ANN=G-chr1:123456_A>T|...

	    ## Multiple alt (G,T) have been split before SnpEff annotation: 1 alt per line in annotated VCF
	    set INFO [lindex $Ls $i_info] 
	    regsub -all {\"} $INFO "" INFO
	    set L_INFO [split $INFO ";"]
	    
	    ## Attribution of the value for each $colName
	    foreach inf $L_INFO {
		if {[regexp "(.*)=(.*)" $inf match colName value]} {
		    # Treatment for frequences and phastCons values: keep only the first value
		    if {[regexp "AF$|FREQ|phastCons" $colName]} {
			if {[regexp "_ESP6500_" $colName]} {continue}
			set $colName $value
			foreach el [split $value ","] {
			    if {$el == "" || $el == "."} {
				set $colName "NA"
			    } else {
				set $colName $el
				break
			    }
			}
		    } else {
			set $colName $value
		    }
		}
	    }
	    set varType        "$VARTYPE"
	    set PPH2class      "$dbNSFP_Polyphen2_HVAR_pred"
	    set SIFTprediction "$dbNSFP_SIFT_pred"
	    set Uniprot        "$dbNSFP_Uniprot_acc"
	    set phastCons      "$PhastCons"
	    ### In case of several prediction (T,D,D), keep the deleterious effect
            if {[regexp "T" $SIFTprediction]} {
		set SIFTprediction "Tolerated"
		if {[regexp "D" $SIFTprediction]} {set SIFTprediction "Deleterious"}		
	    } elseif {[regexp "D" $SIFTprediction]} {
		set SIFTprediction "Deleterious"
	    } else {
		set SIFTprediction "NA"
	    }

	    ## If $PPH environment variable is not defined
	    if {$g_VaRank(pph2Dir) == ""} {
		## D=probably damaging, P=possibly damaging, B=benign
		if {[regexp "^D" $PPH2class]} {set PPH2class "probably damaging"}
		if {[regexp "^P" $PPH2class]} {set PPH2class "possibly damaging"}
		if {[regexp "^B" $PPH2class]} {set PPH2class "neutral"}
	    } else {
		## $PPH environment variable is defined: don't use PPH2class
		set PPH2class "NA"
	    }

	    ## Listing of all genes at this position (e.g. gene="MROHS/FLJ43860")
	    set L_Ann [split $ANN ","]
	    set lGenes {}
	    foreach ann $L_Ann {
		set ann [split $ann "|"]
		lappend lGenes [lindex $ann 3]
	    }
	    set firstGN [lindex $lGenes 0]
	    set lGenes  [lsort -unique $lGenes]
	    # Keep the first annotated gene (the “most deleterious”) in first position
	    set gene "$firstGN"
	    foreach g $lGenes {
		if {$g == ""} {continue}
		if {$g != $firstGN} {lappend gene $g}
	    }
	    set gene [join $gene "/"]

	    ## ALT=TAA --> ANN=TAA|...
	    ## ALT=G,T --> ANN=G|... , T|...
	    ## ALT=C,G (cancer) --> ANN=G-C|... (format = ALT-REFERENCE)
	    ## ALT=G (into MNP)--> ANN=G-chr1:123456_A>T|...	    	

	    # Keep annotations only of the first annotated variant (the “most deleterious”) (only if $codingEffect is not equal to "sequence_feature")
	    foreach ann $L_Ann {
		# cf "http://snpeff.sourceforge.net/SnpEff_manual.html#ann"
		if {[regexp "(\[^|\]*)\\|(\[^|\]*)\\|(\[^|\]*)\\|(\[^|\]*)\\|(\[^|\]*)\\|(\[^|\]*)\\|(\[^|\]*)\\|(\[^|\]*)\\|(\[^|\]*)\\|(\[^|\]*)\\|(\[^|\]*)\\|(\[^|/\]*)/?(\[^|/\]*)\\|(\[^|/\]*)/?(\[^|/\]*)\\|(\[^|/\]*)/?(\[^|/\]*)\\|(\[^|\]*)\\|(\[^|,/\]*)/?(\[^|,/\]*)/?(\[^|/,;\]*)" $ann match Allele Annotation Annotation_Impact Gene_Name Gene_ID Feature_Type Feature_ID Transcript_BioType Rank HGVS.c HGVS.p cDNA.pos cDNA.length CDS.pos CDS.length AA.pos AA.length Distance ERRORS WARNINGS INFOsup]} {
		    #if {![eval regexp "$expression" $INFO match $lVariables]} {} ; #doesn't work

		    foreach colName {Allele Annotation Annotation_Impact Gene_Name Gene_ID Feature_Type Feature_ID Transcript_BioType Rank HGVS.c HGVS.p cDNA.pos cDNA.length CDS.pos CDS.length AA.pos AA.length Distance ERRORS WARNINGS INFOsup} {
			if {[set $colName] == ""} {set $colName "NA"}
		    }

		    if {[info exists g_VaRank(DEBUG)]} {
			puts "INFO: $INFO"
			puts "Allele: $Allele"
			puts "Annotation: $Annotation"
			puts "Annotation_Impact: $Annotation_Impact"
			puts "Gene_Name: $Gene_Name"
			puts "Gene_ID: $Gene_ID"
			puts "Feature_Type: $Feature_Type"
			puts "Feature_ID: $Feature_ID"
			puts "Transcript_BioType: $Transcript_BioType"
			puts "Rank: $Rank"
			puts "HGVS.c: ${HGVS.c}"
			puts "HGVS.p: ${HGVS.p}"
			puts "cDNA.pos: ${cDNA.pos}"
			puts "cDNA.length: ${cDNA.length}"
			puts "CDS.pos: ${CDS.pos}"
			puts "CDS.length: ${CDS.length}"
			puts "AA.pos: ${AA.pos}"
			puts "AA.length: ${AA.length}"
			puts "Distance: $Distance"
			puts "ERRORS: $ERRORS"
			puts "WARNINGS: $WARNINGS"
			puts "INFOsup: $INFOsup"
		    }
		    
		    set codingEffect "$Annotation"
		    set cNomen "$Feature_ID:${HGVS.c}"
		    set pNomen "${HGVS.p}"
		    set wtAA_1 "NA"; set posAA "NA"; set varAA_1 "NA"
		    ## p.His309fs*  p.Phe508del p.His7dup p.Trp33_Lys35delinsArg p.Met1ext-5
		    if {[regexp "p\\.(\[A-Za-z\]{3})(\[0-9\]+)(\[A-Za-z*\]+)" ${HGVS.p} match wtAA_1 posAA varAA_1]} {
			if {[regexp "fs" $varAA_1]} {set varAA_1 "*"}
			set wtAA_1  [correspondanceAA3inAA1 $wtAA_1]
			set varAA_1 [correspondanceAA3inAA1 $varAA_1]
		    }


		    ## Effect into Sequence Ontology (codingEffect)
		    ## Exon: exon_variant exon_loss_variant frameshift_variant missense_variant initiator_codon_variant stop_retained_variant rare_amino_acid_variant stop_lost start_lost stop_gained synonymous_variant start_retained stop_retained_variant 
		    ## Intron: intron_variant conserved_intron_variant splice_acceptor_variant splice_donnor_variant
		    ## Exon/Intron: splice_region_variant

		    set varLocation "NA"
		    set exon "NA"
		    set intron "NA"
		    if {[regexp "exon_variant|exon_loss_variant|frameshift_variant|missense_variant|initiator_codon_variant|stop_retained_variant|rare_amino_acid_variant|stop_lost|start_lost|stop_gained|synonymous_variant|start_retained|stop_retained_variant" $codingEffect]} {
			set varLocation "exon"
			set exon '$Rank'
		    } elseif {[regexp "intron_variant|conserved_intron_variant|splice_acceptor_variant|splice_donnor_variant" $codingEffect]} {
			set varLocation "intron"
			set intron '$Rank'
		    } elseif {[regexp "5_prime" $codingEffect]} {
			set varLocation "5'UTR"
		    } elseif {[regexp "3_prime" $codingEffect]} {
			set varLocation "3'UTR"
		    } elseif {[regexp "downstream" $codingEffect]} {
			set varLocation "downstream"
		    } elseif {[regexp "upstream" $codingEffect]} {
			set varLocation "upstream"
		    }
		    
		    ## codingEffect == frameshift_variant&missense_variant
		    ## keep only the first annotation (important for statistics)
		    set codingEffect [lindex [split $codingEffect "&"] 0]
		    lappend g_VaRank(codingEffect) $codingEffect

		    ## Determine/simplify the ID
		    set ID [AssignID $chrom $pos $ref $alt]		    


		    ## Definition of g_ANNOTATION($ID)
		    ##################################
		    if {[info exists g_ANNOTATION($ID)]} {unset g_ANNOTATION($ID)}
		    set newline {}
		    foreach colName $header {
			lappend newline [set $colName]
		    }
		    ## Don't put a "set" in place of the "lappend"! (else bug during the scoring)
		    lappend g_ANNOTATION($ID) "[join $newline "\t"]"	 
		    #puts "g_ANNOTATION($ID): $g_ANNOTATION($ID)"


		    if {$codingEffect == "sequence_feature"} {continue}
		    break
		}
	    }
	}
    }

    return 
}

proc correspondanceAA3inAA1 {AA3} {
    set c(Gly) G
    set c(Ala) A
    set c(Val) V 
    set c(Leu) L
    set c(Ile) I
    set c(Met) M
    set c(Phe) F
    set c(Trp) W
    set c(Pro) P
    set c(Ser) S
    set c(Thr) T
    set c(Cys) C
    set c(Tyr) Y
    set c(Asn) N
    set c(Gln) Q
    set c(Asp) D
    set c(Glu) E
    set c(Lys) K
    set c(Arg) R
    set c(His) H
    set c(*) *
    
    if {[info exists c($AA3)]} {return $c($AA3)} else {return "NA"}
}
