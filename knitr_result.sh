#    knitr_result.sh - part of RIEMS - Reliable Information Extraction from Metagenomic Sequence datasets
#    Copyright (C) 2009-2016  Ariane Belka, Maria Jenckel, Matthias Scheuch, Dirk Hoeper
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

#!/bin/bash

. ${installdir}/Config.txt    


md5complete=`md5sum ${arbeitsvz}/resultprotocol-complete.txt | cut -d " " -f1`                  # calculate md5sum for resultprotocol-complete.txt and get first column (2nd is filename)
md5compact=`md5sum ${arbeitsvz}/resultprotocol-compact.txt | cut -d " " -f1`                    # calculate md5sum for resultprotocol-compact.txt and get first column (2nd is filename)

if [[ $resultprotocol == "y" ]]                                                                 # if the creation of resultprotocol was selected (default), then ...
    then
        restreads=`cut -f8 ${arbeitsvz}/asignedAcc.txt | grep -c "unclassified"`                # get number of unclassified reads (in asignedAcc.txt as unclassified)
        lowqual=`cut -f8 ${arbeitsvz}/asignedAcc.txt | grep -c "low quality"`                   # get number of low quality reads (in asignedAcc.txt as low quality)
        totalreads=`wc -l < ${arbeitsvz}/asignedAcc.txt`                                        # get number of all reads (number of lines in asignedAcc.txt)
        
        export md5complete                                                                      # export infos neede for subsequent R script
        export md5compact                                                                       # |
        export restreads                                                                        # |
        export lowqual                                                                          # |
        export totalreads                                                                       # |
        export projektname                                                                      # V
        
        echo "Construct: *.tex" >> zeiten.log                                                   # info in zeiten.log
        date >> zeiten.log                                                                      # date to zeiten.log
        ${R}/Rscript -e "library(knitr); knit('${scriptT}')" 12 &>>${arbeitsvz}/console.log     # start R script by using knitR (resultprotocol.Rnw)
        date >> zeiten.log                                                                      # date to zeiten.log
        
        tex2pdf=${arbeitsvz}/resultprotocol.tex                                                 # add filename to variable
        echo "Construct: *.pdf" >> zeiten.log                                                   # info to zeiten.log
        date >> zeiten.log                                                                      # date to zeiten.log
        pdflatex ${tex2pdf} 12 &>>${arbeitsvz}/console.log                                      # write tex-file to pdf
        pdflatex ${tex2pdf} 12 &>>${arbeitsvz}/console.log                                      # write tex-file to pdf (again as errors can occur in first attempt)
        date >> zeiten.log                                                                      # date to zeiten.log
        
        mv resultprotocol.pdf ${projektname}-resultprotocol.pdf                                 # rename resultprotocol.pdf by adding projectname
        mv resultprotocol.tex ${projektname}-resultprotocol.tex                                 # rename resultprotocol.tex by adding projectname
        mv resultprotocol.aux ${projektname}-resultprotocol.aux                                 # rename resultprotocol.aux by adding projectname
        mv resultprotocol.log ${projektname}-resultprotocol.log                                 # rename resultprotocol.log by adding projectname
        
        echo -ne "Done"                                                                         # user info
            
fi

