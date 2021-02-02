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

proc DescendingSortOnElement1 {X Y {N 1}} {
    return [expr {[lindex $X $N]<[lindex $Y $N]}]
}
proc DescendingSortOnElement2 {X Y {N 2}} {
    return [expr {[lindex $X $N]<[lindex $Y $N]}]
}
proc AscendingSortOnElement0 {X Y {N 0}} {
    return [expr {[lindex $X $N]>[lindex $Y $N]}]
}
proc AscendingSortOnElement1 {X Y {N 1}} {
    return [expr {[lindex $X $N]>[lindex $Y $N]}]
}

##############################################################################
#                          WORKING WITH FILES
##############################################################################
proc FirstLineFromFile {{File ""}} {

    if {[regexp ".gz$" $File]} {
	set F [open "|gzip -cd $File"] 
    } else {
	set F [open "$File"]
    }

    set L_Lines {}
    set First 1

    while {[gets $F Line]>=0} {
	if {$First} {set First 0;break}
    }
    close $F

    return $Line
}

proc ContentFromFile {{File ""}} {
    if {[string equal $File ""]} {return ""}
    set f     [open $File r]
    set Text [read -nonewline $f]
    close $f
    return $Text
}

proc LinesFromFile {{File ""}} {
    return [split [ContentFromFile $File] "\n"]
}

proc ContentFromGZFile {{File ""}} {
    if {[string equal $File ""]} {return ""}
    set f     [open "|gzip -cd $File" r]
    set Texte [read -nonewline $f]
    close $f
    return $Texte
}

proc LinesFromGZFile {{File ""}} {
    return [split [ContentFromGZFile $File] "\n"]
}

proc ReplaceTextInFile {NewText fichier} {
    set    fifi [open $fichier w]
    puts  $fifi $NewText
    close $fifi 
    return 1
}


proc WriteTextInFile {text fichier} {
    set    fifi [open $fichier a]
    puts  $fifi $text
    close $fifi 
    return 1
}

proc UniqueSortingOfFiles {L_Files} {
    # Return a list without the same file listed several times.
    # <=> Search for the absolute path and name of each file.

    set new_L_Files {}
    foreach f $L_Files {
	# While $f is a link:
	while {![catch {set f [file readlink $f]} Message]} {
	    continue
	}
	lappend new_L_Files "[file normalize $f]" 
    }

    set new_L_Files [lsort -unique $new_L_Files]

    return $new_L_Files
}



##############################################################################
#                          WORKING WITH rsID
##############################################################################

proc isNotAnRS {rs} {
        if {[regexp "^rs\[0-9\]+$" $rs]} {return 0} else {return 1}
}

