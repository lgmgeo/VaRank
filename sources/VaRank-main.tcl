#!/usr/bin/env tclsh

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

global g_VaRank
global g_allPatients
global g_lPatientsOf
global g_vcfINFOS
global g_ANNOTATION
global g_PPH2
global g_lScore
global g_deltaSSF
global g_deltaMES
global g_deltaNNS
global g_Statistics
global env

## Checking for environment variables needed (if all are defined).
if {![info exists env(VARANK)]} {
    puts "\"VARANK\" environment variable not specified. Please defined it before running VaRank. Exit"; exit
}
if {![info exists env(SNPEFF)]} {
    if {![info exists env(ALAMUT)]} {
	puts "No annotation engine defined:"
        puts "Please defined one of the \"SNPEFF\" or \"ALAMUT\" environment variables before running VaRank."
	puts "Exit"; exit
    }
}


set g_VaRank(sourcesDir) "$env(VARANK)/sources"


if {[info exists env(PPH)] && $env(PPH)!=""} {
    set g_VaRank(pph2Dir) "$env(PPH)"
} else {set g_VaRank(pph2Dir) ""}

puts "Tcl/Tk version: [info tclversion]"
source $g_VaRank(sourcesDir)/VaRank-alamut.tcl
source $g_VaRank(sourcesDir)/VaRank-config.tcl
source $g_VaRank(sourcesDir)/VaRank-exomiser.tcl
source $g_VaRank(sourcesDir)/VaRank-filters.tcl
source $g_VaRank(sourcesDir)/VaRank-general.tcl
source $g_VaRank(sourcesDir)/VaRank-help.tcl
source $g_VaRank(sourcesDir)/VaRank-pph2.tcl
source $g_VaRank(sourcesDir)/VaRank-ranking.tcl
source $g_VaRank(sourcesDir)/VaRank-scoring.tcl
source $g_VaRank(sourcesDir)/VaRank-snpeff.tcl
source $g_VaRank(sourcesDir)/VaRank-stats.tcl
source $g_VaRank(sourcesDir)/VaRank-vcf.tcl
source $g_VaRank(sourcesDir)/VaRank-upstreamDeletion.tcl

puts "VaRank [VaRank_Version]"
puts "VaRank is a program for Ranking genetic Variation from NGS data"
puts ""
puts "Copyright (C) 2016-2021 GEOFFROY Veronique and MULLER Jean"
puts ""
puts "Please feel free to contact us for any suggestions or bug reports"
puts "email: veronique.geoffroy@inserm.fr; jeanmuller@unistra.fr"
puts ""

if {[info exists env(ALAMUT)]} {
    set g_VaRank(alamutDir) "$env(ALAMUT)"
    if {[info exists env(SNPEFF)]} {
	puts "...WARNING:"
	puts "\tBoth ALAMUT and SNPEFF environment variables exist."
	puts "\tRunning VaRank with ALAMUT.\n"
    }
} else {
    set g_VaRank(snpeffDir) "$env(SNPEFF)"
}

#To activate the debug mode and see
#set g_VaRank(DEBUG) 1

## Initialisation of the statistics
if {![info exists g_Statistics]} {
    set g_Statistics(all,Null) 0
    set g_Statistics(all,Het)  0
    set g_Statistics(all,Hom)  0
    set g_Statistics(all,Both) 0
}

## No argument given:
if {$argv == ""} {
    puts "Arguments are missing see help below\n"
    showHelp; exit
}

## Needing help?
if {[regexp -nocase "help" $argv]} {showHelp; exit}

## Downloading configuration:
if {[catch {configureVaRank $argv} Message]} {
    puts "VaRank seems to have a problem during configuration step. Exit"
    puts "######################################################################"
    puts "$Message"
    puts "######################################################################"
    exit
}
if {[catch {ExternalAnnotations} Message]} {
    puts "VaRank seems to have a problem while downloading external annotation. Exit"
    puts "######################################################################"
    puts "$Message"
    puts "######################################################################"
    exit
}

## Downloading VCF files data:
if {[catch {parseVCFfiles} Message]} {
    puts "VaRank seems to have a problem while parsing VCF files. Exit"
    puts "######################################################################"
    puts "$Message"
    puts "######################################################################"
    exit
}

if {[info exists g_VaRank(snpeffDir)]} {
    ## Creation of the SnpEff input files (containing non redundant variants)
    if {[catch {createSnpEffInputFile} Message]} {
	puts "VaRank seems to have a problem while creating the SnpEff input files. Exit"
	puts "######################################################################"
	puts "$Message"
	puts "######################################################################"
	exit
    }
    
    ## Checking SnpEff program:
    if {[catch {checkSnpEff} Message]} {
	puts "VaRank seems to have a problem while checking SnpEff. Exit"
	puts "######################################################################"
	puts "$Message"
	puts "######################################################################"
	exit
    } 
    
    ## Running SnpEff:
    if {[catch {runSnpEff} Message]} {
	puts "VaRank seems to have a problem while running SnpEff. Exit"
	puts "######################################################################"
	puts "$Message"
	puts "######################################################################"
	exit
    } 
    
    ## Downloading SnpEff data:
    if {[catch {parseSnpEffFile} Message]} {
	puts "VaRank seems to have a problem while parsing the SnpEff annotation file. Exit"
	puts "######################################################################"
	puts "$Message"
	puts "######################################################################"
	exit
    }   
} else {
    ## Creation of the alamut input file (1 for all patients)
    if {[catch {createAlamutInputFile} Message]} {
	puts "VaRank seems to have a problem while creating the Alamut Batch input file. Exit"
	puts "######################################################################"
	puts "$Message"
	puts "######################################################################"
	exit
    }
    ## Checking if pre-existing "AlamutAnnotations_all.txt" file is corrupted 
    if {[catch {checkIfCorrupted} Message]} {
	puts "VaRank seems to have a problem while running checks of the AlamutAnnotations_all.txt file. Exit"
	puts "######################################################################"
	puts "$Message"
	puts "######################################################################"
	exit
    }     
    ## Running alamut:
    if {[catch {runAlamut} Message]} {
	puts "VaRank seems to have a problem while running Alamut Batch. Exit"
	puts "######################################################################"
	puts "$Message"
	puts "######################################################################"
	exit
    } 
    ## Downloading alamut data:
    if {[catch {parseAlamutFile} Message]} {
	puts "VaRank seems to have a problem while parsing the Alamut Batch annotation file. Exit"
	puts "######################################################################"
	puts "$Message"
	puts "######################################################################"
	exit
    } 
}

## Creation of the PolyPhen-2 input file:
## If no missense to be analyzed we skip this step
if {[createPPH2Input]!=0} {
    ## Runnning PPH2:
    if {[catch {runPPH2-1by1} Message]} {
	changeModeOfPPH2lockFiles
	puts "VaRank seems to have a problem while running PPH2. Exit"
	puts "######################################################################"
	puts "$Message"
	puts "######################################################################"
	exit
    } 
    changeModeOfPPH2lockFiles
} 

## Preparation of the phenotype-driven analysis (Exomiser)
set L_allGenes [searchForAllGenesContainingVariants]
if {$g_VaRank(hpo) ne "" && $L_allGenes ne ""} {
    checkExomiserInstallation "$L_allGenes"
    runExomiser "$L_allGenes" "$g_VaRank(hpo)"
}

## Scoring genetic variants:
if {[catch {scoreAllTheID} Message]} {
    puts "VaRank seems to have a problem while scoring the variations. Exit"
    puts "######################################################################"
    puts "$Message"
    puts "######################################################################"
    exit
}

## Writing output files
if {[catch {writeAllVariantsRankingByVar} Message]} {
    puts "VaRank seems to have a problem while writing output files (by var). Exit"
    puts "######################################################################"
    puts "$Message"
    puts "######################################################################"
    exit
}
if {[catch {writeAllVariantsRankingByGene} Message]} {
    puts "VaRank seems to have a problem while writing output files (by gene). Exit"
    puts "######################################################################"
    puts "$Message"
    puts "######################################################################"
    exit
}

## Filtering the output files:
if {[catch {executeFilters} Message]} {
    puts "VaRank seems to have a problem while filtering the output files. Exit"
    puts "######################################################################"
    puts "$Message"
    puts "######################################################################"
    exit
}

# Statistics
if {[catch {writeAllStatistics} Message]} {
    puts "VaRank seems to have a problem while writing statistics. Exit"
    puts "######################################################################"
    puts "$Message"
    puts "######################################################################"
    exit
}

if {[info exists g_Statistics(All)]} {
    if {$g_VaRank(SamOut)=="all"} {
	set tot "[llength $g_allPatients]"
    } else {
	set tot "[llength $g_VaRank(SamOut)]"
    }
    puts "...VaRank Global Statistics on $tot samples: [set g_Statistics(All)] variation(s):"
    puts "\t\t[set g_Statistics(all,Both)] found with homozygous and heterozygous status"
    puts "\t\t[set g_Statistics(all,Het)] found with only the heterozygous status"
    puts "\t\t[set g_Statistics(all,Hom)] found with only the homozygous status"
    if {$g_Statistics(all,Null)!=0} {
	puts "\tWARNING: [set g_Statistics(all,Null)] found in none sample (GT=0/0 in all samples or no depth of coverage)."
	puts "\tWARNING: These variations are not included in the annotated output files."
    }
}
puts "...VaRank is done with the analysis ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"


