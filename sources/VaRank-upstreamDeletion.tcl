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

proc geneForStarID {starID i_gene} {

    global g_VaRank
    global g_ANNOTATION
    # In case of ID with ALT="*" (stands for a upstream deletion), retrieve gene name annotation

    set gene "NA"

    if {![regexp "(.+)_(.+)_(.+)_\\*" $starID match chrom pos ref]} {return "$gene"}

    set saveIDfile "$g_VaRank(vcfDir)/VCF_Coordinates_Conversion.tsv"
    if {![file exists $saveIDfile]} {return "$gene"}

    set consensusID ""
    foreach L [LinesFromFile $saveIDfile] {
	if {[regexp "$chrom\t$pos\t$ref\t\[^*\]" $L]} {
	    set consensusID [lindex $L 0]
	    break
	}
    }
    if {[info exists g_ANNOTATION($consensusID)]} {
	set gene [lindex $g_ANNOTATION($consensusID) $i_gene]
    }

    return "$gene"
}

