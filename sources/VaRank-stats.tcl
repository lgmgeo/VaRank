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

proc writeAllStatistics {} {

    global g_VaRank
    global g_lPatientsOf
    global g_Statistics
    global g_allPatients
    
    puts "...writing patients and global statistics ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"

    set globalOutputfile "$g_VaRank(vcfDir)/SNV_global_statistics.tsv"

    if {[info exists g_VaRank(snpeffDir)]} {
	set L_codingEffect $g_VaRank(codingEffect)
	set L_varLocation "unknown"
    } else {
	#Alamut coding effect useful for statistics
	set L_codingEffect [list synonymous missense "stop gain" in-frame frameshift "start loss" "stop loss"]
	set L_varLocation  [list intron upstream "5'UTR" "3'UTR" downstream "splice site"]
    }
 
    set lEffects [lsort [concat $L_codingEffect $L_varLocation]]
    lappend lEffects "unknown"

    # need to format for the output
    set maxlength 0
    foreach effect $lEffects {	
	set len [string length $effect]
	if {$len > $maxlength} {set maxlength $len}
    }
    incr maxlength

    foreach effect $lEffects {set L_$effect {}}
    set L_Total {}
    foreach fam [array names g_lPatientsOf] {
	foreach patient $g_lPatientsOf($fam) {	 
   
	    if {$g_VaRank(SamOut)!="all" && [lsearch -exact -nocase $g_VaRank(SamOut) $patient]==-1} {continue}

	    if {![info exists g_Statistics($patient)]} {puts "\t...WARNING: $patient statistics could not be written because VaRank outputdata already existed."; continue}
	    set Line {}
	    ## Header of the patientOutputfile
	    lappend Line [join [list [format "%-${maxlength}s" What] Total Homozygous Heterozygous Null] "\t"]
	    foreach effect $lEffects {
		## Important if SnpEff annotation is used
		if {![info exists g_Statistics($patient,$effect)]} {set g_Statistics($patient,$effect) 0}
		if {![info exists g_Statistics($patient,$effect,Hom)]} {set g_Statistics($patient,$effect,Hom) 0}
		if {![info exists g_Statistics($patient,$effect,Het)]} {set g_Statistics($patient,$effect,Het) 0}
		if {![info exists g_Statistics($patient,$effect,Null)]} {set g_Statistics($patient,$effect,Null) 0}
		lappend L_$effect  [set g_Statistics($patient,$effect)]
		lappend Line [join [list [format "%-${maxlength}s" $effect] [set g_Statistics($patient,$effect)] [set g_Statistics($patient,$effect,Hom)] [set g_Statistics($patient,$effect,Het)] [set g_Statistics($patient,$effect,Null)]] "\t"]
	    }
	    lappend Line "\nTotal:  [set g_Statistics($patient)] variations."
	    
	    lappend L_Total  [set g_Statistics($patient)]

	    set patientOutputfile "$g_VaRank(vcfDir)/[set fam]_[set patient]_statistics.tsv"
	    if {[file exists $patientOutputfile]} {file delete -force $patientOutputfile}
	    WriteTextInFile [join $Line "\n"] $patientOutputfile
	}
    }	

    set     Line {}
    if {$g_VaRank(SamOut)=="all"} {
	lappend Line "Global statistics (non redundant) on [llength $g_allPatients] samples:\n"
    } else {
	lappend Line "Global statistics (non redundant) on [llength $g_VaRank(SamOut)] samples:\n"
    }
    lappend Line [join [list [format "%-${maxlength}s" What] "Total" Mean SD] "\t"]
    foreach effect $lEffects {	
	if {![info exists g_Statistics(All,$effect)]} {continue}

	set Mean NA
	set SD   NA

	if {[llength [set L_$effect]>1] && [set L_$effect]!=0 && [set L_$effect]!={}} {

	    set MVSD [BasicStatistics [set L_$effect]]
	    set Mean [format "%.0f" [lindex $MVSD 0]]
	    set SD   [format "%.0f" [lindex $MVSD 2]]
	}
	lappend Line [join [list [format "%-${maxlength}s" $effect] [set g_Statistics(All,$effect)] $Mean $SD] "\t"]
    }

    set Mean NA
    if {[llength [set L_Total]]>1} {
	set MVSD [BasicStatistics [set L_Total]]
	set Mean [format "%.0f" [lindex $MVSD 0]]
    }
    if {[info exists g_Statistics(All)]} {
	lappend Line "\n[join [list [format "%-${maxlength}s" Total] [set g_Statistics(All)] $Mean] "\t"]"
    }
    

    if {[llength $Line]>2} {
	if {[file exists $globalOutputfile]} {file delete -force $globalOutputfile}
	WriteTextInFile [join $Line "\n"] $globalOutputfile
    } else {
	puts "\t...WARNING: $patient statistics could not be written because VaRank outputdata already existed."
    }
    return 
}

proc BasicStatistics {Liste} {

    #Allow the calculation of Basic Statistics
    #Mean, Variance, Standart deviation
    
    #SV  Sum of values
    #SC  Sum of Squares
    #MV  Mean of Values
    #MC  Mean of Squares
    #Var Variance

    set SC  "0.0"
    set SV  "0.0"
    set Var "0.0"
    set SD  "0.0"

    set NV [llength $Liste]
    foreach V $Liste {
	set SV [expr {$SV+$V}]
	set SC [expr {($SC+pow($V,2))*1.0}]
    }    
    
    #puts "SV [join $SV ","]"
    #puts "NV [join $NV ","]"

    set MV  [expr {$SV/$NV}]
    set MC  [expr {$SC/$NV}]
    set Var [expr {$MC-pow($MV,2)}]

    if {$Var>0.0} {
	set SD  [expr {sqrt($Var)}]
    } else {
	set SD 0.0
    }
    #puts "Nb: $NV - SV: $SV - SC: $SC - MV: $MV - MC: $MC - Var: $Var - SD: $SD"

    return [list $MV $Var $SD]
}
