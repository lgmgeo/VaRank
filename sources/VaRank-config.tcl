############################################################################################################
# VaRank 1.5.1                                                                                             #
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

## Loading the different options in the following order:
## - Default options
## - Config file options (if file exists)
## - Options given in arguments
#
## Please note: the case in the name of options is not important. Ex: "vcfdir" = "vcfDir"
##
 
proc VaRank_Version {} {

    #VaRank Version 

    global g_VaRank

    set g_VaRank(Version) "1.4.3"

    return $g_VaRank(Version)
}

proc configureVaRank {argv} {

    global g_VaRank
    global g_lPatientsOf
    global env
    global L_hpo


    puts "...downloading the configuration data ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"

    ## We must have a vcf input file
    ################################
    set  k [lsearch -regexp -nocase $argv "vcfDir"]
    if {$k ne -1} {
	set g_VaRank(vcfDir) [lindex $argv [expr {$k+1}]]
    } else {
	puts "VaRank needs in argument a study directory containing your vcf input file (-vcfDir ...) - Exit."
	exit
    }
    set  liste_vcfFile [glob -nocomplain $g_VaRank(vcfDir)/*.vcf $g_VaRank(vcfDir)/*.vcf.gz]
    if {$liste_vcfFile eq ""} {
	puts "VaRank needs a VCF file in input. No VCF file present in $g_VaRank(vcfDir) - Exit."; exit
    }

    ## Load default options
    #######################
    puts "\t...configuration data by default"
    regsub "sources" $g_VaRank(sourcesDir) "pph2Databases" pph2Dir

    set g_VaRank(AlamutAlltrans) "yes"
    set g_VaRank(alamutHumanDB) "GRCh37"
    set g_VaRank(AlamutProcesses) 0
    
    set g_VaRank(B_phastCons) 5 
    set g_VaRank(B_SIFT)      5
    set g_VaRank(B_PPH2)      5

    set g_VaRank(dbNSFP)        ""
    set g_VaRank(dbSNP)         ""
    set g_VaRank(depthFilter)   10
    set g_VaRank(extann)        ""
    set g_VaRank(freqFilter)  0.01
    set g_VaRank(Homstatus)   "no"
    set g_VaRank(Homcutoff)     80
    set g_VaRank(javaPath)  "java"

    set g_VaRank(MEScutoff)          -15
    set g_VaRank(metrics)           "us"
    set g_VaRank(msigdb)              ""
    set g_VaRank(NNScutoff)          -10
    set g_VaRank(phastConsCutoff)   0.95
    set g_VaRank(phastConsDB)         ""
    set g_VaRank(proxyPasswd) "password"
    set g_VaRank(proxyPort)       "8080"
    set g_VaRank(proxyServer)         ""
    set g_VaRank(proxyUser)           ""
    set g_VaRank(readFilter)          10
    set g_VaRank(readPercentFilter)   15
    set g_VaRank(refseq)     "$pph2Dir/human.protein.faa.gz"
    set g_VaRank(rsFilter)   removeNonPathoRS
    set g_VaRank(rsFromVCF)          "no"

    set g_VaRank(S_Known)          110
    set g_VaRank(S_CloseSplice)     70
    set g_VaRank(S_DeepSplice)      25
    set g_VaRank(S_EssentialSplice) 90
    set g_VaRank(S_ExonIntron)       2
    set g_VaRank(S_Fs)             100
    set g_VaRank(S_Inframe)         30
    set g_VaRank(S_LSEstrong)       40
    set g_VaRank(S_LSEweak)         35
    set g_VaRank(S_Missense)        50
    set g_VaRank(S_StartLoss)       80
    set g_VaRank(S_StopGain)       100
    set g_VaRank(S_StopLoss)        80
    set g_VaRank(S_Synonymous)      10
    set g_VaRank(S_UTR)              1

    set g_VaRank(skipAlamutChecks) "no"

    set g_VaRank(snpeffHumanDB)     "" 
    set g_VaRank(SnpEffTestsDir) "Tests/SnpEffTests/SnpEff/checkSnpEff"
    set g_VaRank(SSFcutoff)         -5

    set g_VaRank(SamOut)         "all"
    set g_VaRank(SamVa)           "no"

    set g_VaRank(uniprot) "$pph2Dir/HUMAN.fasta.gz"

    set g_VaRank(vcfFields)  {all}
    set g_VaRank(vcfInfo)    "no"


    set g_VaRank(L_outputColHeader) "variantID gene rsID chr start end ref alt zygosity totalReadDepth varReadDepth varReadPercent QUALphred Uniprot codingEffect varLocation exon intron varType Annotation_Impact Gene_ID Feature_Type Feature_ID Transcript_BioType cNomen pNomen wtAA_1 posAA varAA_1 cDNA.pos cDNA.length CDS.pos CDS.length AA.pos AA.length Distance LOF NMD dbNSFP_1000Gp1_AF dbNSFP_1000Gp1_AFR_AF dbNSFP_1000Gp1_AMR_AF dbNSFP_1000Gp1_ASN_AF dbNSFP_1000Gp1_EUR_AF dbNSFP_CADD_phred dbNSFP_ExAC_AC dbNSFP_ExAC_AF dbNSFP_ExAC_AFR_AC dbNSFP_ExAC_AFR_AF dbNSFP_ExAC_AMR_AC dbNSFP_ExAC_AMR_AF dbNSFP_ExAC_Adj_AC dbNSFP_ExAC_Adj_AF dbNSFP_ExAC_EAS_AC dbNSFP_ExAC_EAS_AF dbNSFP_ExAC_FIN_AC dbNSFP_ExAC_FIN_AF dbNSFP_ExAC_NFE_AC dbNSFP_ExAC_NFE_AF dbNSFP_ExAC_SAS_AC dbNSFP_ExAC_SAS_AF dbNSFP_FATHMM_pred dbNSFP_GERP___NR dbNSFP_GERP___RS dbNSFP_LRT_pred dbNSFP_MetaSVM_pred dbNSFP_MutationAssessor_pred dbNSFP_MutationTaster_pred dbNSFP_PROVEAN_pred dbNSFP_Polyphen2_HDIV_pred dbNSFP_Polyphen2_HVAR_pred dbNSFP_SIFT_pred phastCons SIFTprediction PPH2class varankVarScore annotationAnalysis avgTotalDepth sdTotalDepth countTotalDepth avgVariantDepth sdVariantDepth countVariantDepth familyBarcode barcode homCount hetCount alleleCount sampleCount alleleFrequency samVa"
    if {[info exists env(ALAMUT)]} {
	set g_VaRank(L_outputColHeader) "variantID gene geneDesc omimId transcript strand transLen cdsLen chr start end ref alt protein Uniprot posAA wtAA_1 wtCodon varAA_1 zygosity totalReadDepth varReadDepth varReadPercent QUALphred varType codingEffect varLocation assembly exon intron gDNAstart gDNAend gNomen cDNAstart cDNAend cNomen pNomen clinVarIds clinVarOrigins clinVarMethods clinVarClinSignifs clinVarReviewStatus clinVarPhenotypes rsId rsValidations rsValidationNumber rsAncestralAllele rsHeterozygosity rsClinicalSignificance rsMAF rsMAFAllele rsMAFCount 1000g_AF 1000g_AFR_AF 1000g_SAS_AF 1000g_EAS_AF 1000g_EUR_AF 1000g_AMR_AF gnomadAltFreq_all gnomadAltFreq_afr gnomadAltFreq_amr gnomadAltFreq_asj gnomadAltFreq_eas gnomadAltFreq_sas gnomadAltFreq_nfe gnomadAltFreq_fin gnomadAltFreq_oth gnomadAltFreq_popmax gnomadHomCount_all gnomadHomCount_afr gnomadHomCount_amr gnomadHomCount_asj gnomadHomCount_eas gnomadHomCount_sas gnomadHomCount_nfe gnomadHomCount_fin gnomadHomCount_oth gnomadHetCount_all gnomadHetCount_afr gnomadHetCount_amr gnomadHetCount_asj gnomadHetCount_eas gnomadHetCount_sas gnomadHetCount_nfe gnomadHetCount_fin gnomadHetCount_oth gnomadHemCount_all gnomadHemCount_afr gnomadHemCount_amr gnomadHemCount_asj gnomadHemCount_eas gnomadHemCount_sas gnomadHemCount_nfe gnomadHemCount_fin gnomadHemCount_oth gnomadFilter gnomadReadDepth gnomadOrigin deltaMaxEntScorePercent wtMaxEntScore varMaxEntScore deltaSSFscorePercent wtSSFScore varSSFScore deltaNNSscorePercent wtNNSScore varNNSScore nearestSSChange distNearestSS nearestSSType localSpliceEffect localSpliceAnnotation localSS_pos localSS_wtMaxEntScore localSS_varMaxEntScore localSS_wtNNSScore localSS_varNNSScore localSS_wtSSFScore localSS_varSSFScore branchPointPos branchPointChange proteinDomain1 proteinDomain2 proteinDomain3 proteinDomain4 SIFTprediction SIFTweight SIFTmedian PPH2pred phyloP phastCons granthamDist AGVGDclass AGVGDgv AGVGDgd varankVarScore annotationAnalysis avgTotalDepth sdTotalDepth countTotalDepth avgVariantDepth sdVariantDepth countVariantDepth familyBarcode barcode homCount hetCount alleleCount sampleCount alleleFrequency samVa"
    } 
    
    ## Load config file options and output column names
    ###################################################
    set lOptionsOk "AlamutAlltrans alamutHumanDB AlamutProcesses B_phastCons B_PPH2 B_SIFT dbNSFP dbSNP depthFilter extann freqFilter Homcutoff Homstatus javaPath MEScutoff metrics msigdb NNScutoff phastConsCutoff phastConsDB proxyPasswd proxyPort proxyServer proxyUser readFilter readPercentFilter refseq rsFilter rsFromVCF S_CloseSplice S_DeepSplice S_EssentialSplice S_ExonIntron S_Fs S_Inframe S_Known S_LSEstrong S_LSEweak S_Missense S_StartLoss S_StopGain S_StopLoss S_Synonymous  S_UTR skipAlamutChecks SamOut SamVa snpeffHumanDB SnpEffTestsDir SSFcutoff uniprot vcfDir vcfFields vcfInfo"
    set L_outputColHeaderBis ""
    set configFile "$g_VaRank(vcfDir)/configfile"
    if {[file exists $configFile]} {
	puts "\t...configuration data from $configFile"
	foreach L [LinesFromFile $configFile] {
	    if {[regexp "^#" $L]} {continue}
	    if {$L eq "" || [regexp "^\[ \t\]+$" $L]} {continue}
	    
	    # extracting the families
	    if {[regexp "^(fam\[0-9\]+) *: *(.+)" $L match family lpat]} {
		foreach val [split [string trim $lpat] " "] {
		    lappend g_lPatientsOf($family) $val
		}
		continue
	    }
	    # extracting the phenotype (HPO)
	    if {[regexp "^sample = +(\[^ \t\]+)" $L match sample]} {
		continue
	    }
	    if {[regexp "^HP:" $L]} {
		regsub -all " |;" $L "," L
		set L_hpo($sample) $L
		continue
	    }
	    # extracting the options
	    if {[regexp "^-" $L]} {
		regsub -all "^-|:" $L "" L
		set optionName  [lindex $L 0]
		set optionValue [lindex $L 1]
		set k [lsearch -exact -nocase $lOptionsOk $optionName]
		if {$k ne -1} {
		    set optionName [lindex $lOptionsOk $k]
		    set g_VaRank($optionName) $optionValue
		} else {
		    puts "############################################################################"
		    puts "\"$optionName\" option not known."
		    puts "For more information on the arguments, please use the -help option"
		    puts "############################################################################"
		    exit
		}
		continue
	    }	
	    # extracting the output column header
	    lappend L_outputColHeaderBis "[lindex $L 0]"
	}
	if {$L_outputColHeaderBis ne {}} {set g_VaRank(L_outputColHeader) "$L_outputColHeaderBis"}
    }

    ## Load options given in arguments
    ##################################
    puts "\t...configuration data given in arguments"
    regsub -all "^-" $argv "" argv
    set i 0
    set j 1
    while {$j < [llength $argv]} {
	set optionName [lindex $argv $i]
	regsub -all "^-|:" $optionName "" optionName

	set optionValue [lindex $argv $j]
	set  k [lsearch -exact -nocase $lOptionsOk $optionName]
	if {$k ne -1} {
	    set optionName [lindex $lOptionsOk $k]
	    set g_VaRank($optionName) $optionValue
	} else {
	    puts "\"$optionName\" option not known."
	    puts "For more information on the arguments, please use the -help option"
	    exit
	}

	incr i 2
	incr j 2
    }

    puts "\t...checking configuration data"
    ##########################
    ##########################

    ## It must be existing directories
    foreach option "SnpEffTestsDir" {
	regsub "sources" $g_VaRank(sourcesDir) $g_VaRank($option) g_VaRank($option)
	if {![file isdirectory $g_VaRank($option)]} {file mkdir $g_VaRank($option)}
    }


    ## It must be an integer value for the -readFilter -depthFilter -S_... -B_... options.
    foreach option "readFilter depthFilter S_Known S_StopGain S_Fs S_EssentialSplice S_StartLoss S_StopLoss S_Missense S_CloseSplice S_LSEstrong S_LSEweak S_Inframe S_Synonymous S_DeepSplice S_ExonIntron S_UTR B_phastCons B_SIFT B_PPH2 AlamutProcesses" {
	if {![regexp "^\[0-9\]+$" $g_VaRank($option)]} {
	    puts "############################################################################"
	    puts "Bad option value : -$option = $g_VaRank($option)"	
	    puts "Should be an integer. Exit"
	    puts "############################################################################"
	    exit
	}
    }

    ## It must be an existing files    
    foreach option "uniprot refseq" {
	if {![file exists "$g_VaRank($option)"] && $g_VaRank($option) ne ""} {
	    puts "############################################################################"
	    puts "Bad option value for $option option, file does not exists ($g_VaRank($option))."
	    puts "############################################################################"
	    set g_VaRank($option) ""
	}
    }
    if {$g_VaRank(uniprot) eq "" && [file exists "$pph2Dir/HUMAN.fasta.gz"]} {set g_VaRank(uniprot) "$pph2Dir/HUMAN.fasta.gz"}
    if {$g_VaRank(refseq) eq "" && [file exists "$pph2Dir/human.protein.faa.gz"]} {set g_VaRank(refseq) "$pph2Dir/human.protein.faa.gz"}
    
    ## It must be external annotation files (one or several, tsv) for the "-extann" option
    ## Check that extann files exist and are tsv.
    foreach extFile $g_VaRank(extann) {
	if {![file exists $extFile] && $extFile ne ""} {
	    set g_VaRank(extann) ""
	    puts "############################################################################"
	    puts "Bad option value for extann, file does not exists $extFile."
	    puts "############################################################################"
	    break
	}
    }
    ## Automaticly add the files contained in the /ExtAnn directory
    regsub "sources" $g_VaRank(sourcesDir) "ExtAnn" extAnnDir
    foreach annotFile [glob -nocomplain $extAnnDir/*] {
	if {[regexp "^README$|^results.txt$" [file tail $annotFile]]} {continue}
	lappend g_VaRank(extann) $annotFile
    }
    ## Sort -unique on all these files
    set g_VaRank(extann) [UniqueSortingOfFiles $g_VaRank(extann)]

    ## It must be "removeNonPathoRS" or "none" for the -rsFilter option.
    if {![regexp -nocase "removeNonPathoRS|none" $g_VaRank(rsFilter)]} {
	puts "############################################################################"
	puts "Bad option value: -rsFilter = $g_VaRank(rsFilter)"
	puts "Should be \"removeNonPathoRS\" or \"none\""
	puts "############################################################################"
	exit
    }
    
    ## It must be a value comprised between 0 and 100 for the -Homcutoff
    foreach option "Homcutoff readPercentFilter" {
	if {$g_VaRank($option) > 100 || $g_VaRank($option) < 0} {
	    puts "Bad option value: -$option = $g_VaRank($option)"
	    puts "Should be an integer comprised into the range values: \[0,100\]"
	    exit
	}
    }

    ## It must be a value comprised between -100 and 0 for the -MEScutoff -SSFcutoff -NNScutoff options.
    foreach option "MEScutoff SSFcutoff NNScutoff" {
	if {$g_VaRank($option) < -100 || $g_VaRank($option) > 0} {
	    puts "############################################################################"
	    puts "Bad option value: -$option = $g_VaRank($option)"
	    puts "Should be an integer comprised into the range values: \[-100,0\]"
	    puts "############################################################################"
	    exit
	}
    }

    ## It must be a value comprised between 0 and 1 for the -phastConsCutoff option.
    if {$g_VaRank(phastConsCutoff) < 0 ||$g_VaRank(phastConsCutoff) > 1} {
	puts "############################################################################"
	puts "Bad option value: -phastConsCutoff = $g_VaRank(phastConsCutoff)"
	puts "Should be a value comprised into the range values: \[0,1\]"
	puts "############################################################################"
	exit
    }

    ## It must be a value comprised between 0 and 1 for the freqFilter option.
    if {$g_VaRank(freqFilter) < 0 || $g_VaRank(freqFilter) > 1} {
	puts "############################################################################"
	puts "Bad option value: -freqFilter = $g_VaRank(freqFilter)"
	puts "Should be a value comprised into the range values: \[0,1\]"
	puts "############################################################################"
	exit
    }

    ## It must be "fr" or "us" for the metrics option.
    if {![regexp -nocase "^(fr)|(us)$" $g_VaRank(metrics)]} {
	puts "############################################################################"
	puts "Bad option value: -metrics = $g_VaRank(metrics)"
	puts "Should be \"fr\" or \"us\""
	puts "############################################################################"
	exit
    }

    if {[info exists g_VaRank(snpeffDir)]} {

	## javaPath file be correct
 	set g_VaRank(javaSnpEffCommand) "$g_VaRank(javaPath) -Xmx4g -jar $g_VaRank(snpeffDir)"

        ## It must be "NCBI36|GRCh37|GRCh38" for the -alamutHumanDB option.
        if {![regexp "NCBI36|GRCh37|GRCh38" $g_VaRank(alamutHumanDB)]} {
	    puts "############################################################################"
            puts "Bad option value: -alamutHumanDB = $g_VaRank(alamutHumanDB)"
            puts "Value should be something like NCBI36, GRCh37 or GRCh38"
	    puts "############################################################################"
            exit
        }


	if {[info exists g_VaRank(snpeffDir)]} {
	    ## It must be something like "GRCh37.75|hg19" for the -snpeffHumanDB option.
	    if {![regexp -nocase "^GRCh|^hg" $g_VaRank(snpeffHumanDB)]} {
	        puts "############################################################################"
		puts "Bad option value: -snpeffHumanDB = $g_VaRank(snpeffHumanDB)"
		puts "Value should be something like \"GRCh37.75\" or \"hg19\""
	        puts "############################################################################"
		exit
	    }
	    ## It must be something like "/.../.../dbNSFP2.4.txt.gz" for the -dbNSFP option.
	    if {![regexp -nocase "dbNSFP.*.txt.gz$" $g_VaRank(dbNSFP)]} {
	        puts "############################################################################"
		puts "Bad option value: -dbNSFP = $g_VaRank(dbNSFP)"
		puts "Value should be a file named like \"*dbNSFP*.txt.gz\""
	        puts "############################################################################"
		exit
	    }
	    if {![file exists $g_VaRank(dbNSFP)]} {
	        puts "############################################################################"
		puts "Bad option value: -dbNSFP = $g_VaRank(dbNSFP)"
		puts "$g_VaRank(dbNSFP) file doesn't exist"
	        puts "############################################################################"
		exit
	    }
	    ## It must be something like "db/dbSNP.2015-01-09_00-All.vcf" for the -dbSNP option.
	    if {![regexp -nocase ".*.vcf(.gz)?$" $g_VaRank(dbSNP)]} {
	        puts "############################################################################"
		puts "Bad option value: -dbSNP = $g_VaRank(dbSNP)"
		puts "Should be a file.vcf"
	        puts "############################################################################"
		exit
	    }
	    if {![file exists $g_VaRank(dbSNP)]} {
	        puts "############################################################################"
		puts "Bad option value: -dbSNP = $g_VaRank(dbSNP)"
		puts "$g_VaRank(dbSNP) file doesn't exist"
	        puts "############################################################################"
		exit
	    }
	    ## It must be a directory for the -phastConsDB option.
	    if {![file isdirectory $g_VaRank(phastConsDB)]} {
	        puts "############################################################################"
		puts "Bad option value: -phastConsDB = $g_VaRank(phastConsDB)"
		puts "$g_VaRank(phastConsDB) is not a directory"
	        puts "############################################################################"
		exit
	    }
	    ## It must be something like "msigdb.v4.0.symbols.gmt" for the -msigdb option.
#	    if {![regexp -nocase "msigdb.*.symbols.gmt$" $g_VaRank(msigdb)]} {
#	        puts "############################################################################"
#		puts "Bad option value: -msigdb = $g_VaRank(msigdb)"
#		puts "Value should be something like \"$g_VaRank(snpeffDir)/msigdb.v4.0.symbols.gmt\""
#	        puts "############################################################################"
#		exit
#	    }
#	    if {![file exists $g_VaRank(msigdb)]} {
#	        puts "############################################################################"
#		puts "Bad option value: -msigdb = $g_VaRank(msigdb)"
#		puts "$g_VaRank(msigdb) file doesn't exist"
#	        puts "############################################################################"
#		exit
#	    }
	}
    }

    puts "\t******************************************"
    puts "\tVaRank has been run with these arguments :"
    puts "\t******************************************"
    set lKey [array names g_VaRank]
    foreach key [lsort $lKey] {
	if {[regexp "L_outputColHeader|skipAlamutChecks" $key]} {continue}
	puts "\t-[format "%-20s %s" $key $g_VaRank($key)]"
    }
    puts "\t******************************************\n"
    

    return	
}
