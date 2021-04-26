#!/usr/bin/env tclsh

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
global L_hpo

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

puts "Tcl version: [info tclversion]"
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
puts "...running on [exec hostname]"

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
configureVaRank $argv
ExternalAnnotations

## Downloading VCF files data:
parseVCFfiles

if {[info exists g_VaRank(snpeffDir)]} {
    ## Creation of the SnpEff input files (containing non redundant variants)
    createSnpEffInputFile
    checkSnpEff
    runSnpEff
    parseSnpEffFile
} else {
    createAlamutInputFile
    checkIfCorrupted
    runAlamut
    parseAlamutFile
}

## Creation of the PolyPhen-2 input file:
## If no missense to be analyzed we skip this step
if {[createPPH2Input]!=0} {
    ## Runnning PPH2:
    runPPH2-1by1
    changeModeOfPPH2lockFiles
} 

## Phenotype-driven analysis (Exomiser)
if {[info exists L_hpo]} {
    set L_3infos [checkExomiserInstallation] ;# {$port $applicationPropertiesTmpFile $idService}
    foreach sample [array names L_hpo] {
	set L_allGenes [searchForAllGenesContainingVariants $sample]
	if {$L_hpo($sample) ne "" && $L_allGenes ne ""} {
	    runExomiser "$sample" "$L_allGenes" "$L_hpo($sample)" "$L_3infos"
	}
    }
    # End the REST service
    set idService [lindex $L_3infos 2]
    if {[catch {exec kill -9 $idService} Message]} {
	puts "End the REST service:"
	puts $Message
    }
}

## Scoring genetic variants:
scoreAllTheID


## Writing output files
writeAllVariantsRankingByVar
writeAllVariantsRankingByGene

## Filtering the output files:
executeFilters

# Statistics
writeAllStatistics

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


