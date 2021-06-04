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


proc createTheSQLiteDB {} {
    
    global g_VaRank

    # Creation of the "db_vcfData" SQLite database to store the VCF data
    puts "...creation of the tmp SQLite database to load the VCF data"
    regsub "sources" $g_VaRank(sourcesDir) "bash" bashDir
    set dbFileTail "VaRank_vcfData_tmp_[clock format [clock seconds] -format "%Y%m%d-%H%M%S"].db"
    set sqLiteDBfile "$g_VaRank(vcfDir)/$dbFileTail"
    puts "\t$sqLiteDBfile"
    
    set diskNFS ""
    catch {set diskNFS [exec bash $bashDir/isNFSdisk.bash $g_VaRank(vcfDir)]}
    if {[regexp -nocase "nfs" $diskNFS]} {
	if {![catch {set dfkh [exec df -kh [file dirname $sqLiteDBfile] | tail -1]} Message]} {
	    puts "\t...on the filesystem: [lindex $dfkh 0]"
	}
	puts "\t   WARNING: Using NFS volumes ($diskNFS) to write the SQLite database can lead to an increase of the running time."
    }
    file delete -force "$sqLiteDBfile"
    file delete -force "${sqLiteDBfile}-journal"
    sqlite3 db_vcfData "$sqLiteDBfile"

    db_vcfData eval {PRAGMA foreign_keys = ON;}
    
    db_vcfData eval {
	CREATE TABLE vcfSVdata (
				chrom    TEXT,
				pos      INT,
				ref      TEXT,
				alt      TEXT,
				rsID     TEXT,
				rsValid  TEXT,
				ID       TEXT UNIQUE
				)
    }    
    db_vcfData eval {
	CREATE TABLE sampleName (
				 sampleName_id INTEGER NOT NULL PRIMARY KEY,
				 sampleName    TEXT NOT NULL
				 )
    }
    db_vcfData eval {
	CREATE TABLE vcfSampleData (
				ID            TEXT,
				homhet        TEXT,
				dp            INTEGER,
				nr            INTEGER,
				qual          TEXT,
				INFO          TEXT,
				sampleName_id INTEGER,
				FOREIGN KEY(sampleName_id) REFERENCES sampleName(sampleName_id),
				FOREIGN KEY(ID) REFERENCES vcfSVdata(ID)
				)
    }

    return $sqLiteDBfile
}


proc closeTheSQLiteDB {sqLiteDBfile} {
    
    # Close and delete the SQLite db
    puts "...deletion of the tmp SQLite database"
    puts "\t($sqLiteDBfile)"
    db_vcfData close    
    file delete -force "$sqLiteDBfile"

    return
}
