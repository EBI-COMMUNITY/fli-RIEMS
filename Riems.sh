#    RIEMS.sh - part of RIEMS - Reliable Information Extraction from Metagenomic Sequence datasets
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
    
. /hps/nobackup/nucleotide/rahman/development/fli-RIEMS/Config.txt                                                                                            # <- Path of config file
. /hps/nobackup/nucleotide/rahman/development/fli-RIEMS/functions.sh                                                                                          # and path to functions
. /hps/nobackup/nucleotide/rahman/development/fli-RIEMS/TaxidDetermination.sh                                                                                 # and path to more functions
    
echo -e "\n
***************** RIEMS 4.0 - Reliable Information Extraction of Metagenomic Sequence datasets ******************    
*                                                                                                               *
*                 Copyright (C) 2009-2016  Ariane Belka, Maria Jenckel, Matthias Scheuch, Dirk Hoeper           *
*                 This program comes with ABSOLUTELY NO WARRANTY                                                *
*                 This is free software, and you are welcome to redistribute it                                 *
*                 under certain conditions.                                                                     *
*                                                                                                               *
*****************************************************************************************************************\n    
RIEMS 4.0 Workflow started at `date`\n\n"                                                                       # user info

resultprotocol=y
furtherAnalysis=n
viriems=n
unset projektname
unset arbeitsvz
while getopts f:p:o:j:i:t:x:r:v: opt 2>/dev/null                                                                # input options
    do                                                                                                          
        case ${opt} in                                                                                          # get input parameters
            o)                              
                outputordner=$OPTARG;;                                                                          # enter working directory (absolute path)
            p)                              
                projektname=$OPTARG;;                                                                           # if using multiple sequence files, a name for the output folder must be specified
            j)                                                                                                             
                unspecified_input+=("$OPTARG");;                                                                # path and name to input date if not format is specified 
            f)                              
                furtherAnalysis=$OPTARG;;                                                                       # if -c y, then copy input file to output folder
            i)          
                illumina_input+=("$OPTARG");;                                                                   # specify illumina input files
            t)          
                ionTorrent_input+=("$OPTARG");;                                                                 # specify ionTorrent input files
            x)          
                tax=(${OPTARG//,/ });;                                                                          # input of taxid for prior screening (comma-seperated)
            r)
                resultprotocol=$OPTARG;;
            v)  
                viriems=$OPTARG;;
            *)                                                                                                  # at unvalid input                                                                                                                        
                cat ${installdir}/input_info.txt                                                                # show info and 
                exit;;                                                                                          # exit program                                                                                     
        esac                            
    done                        
    
if [ $# -lt 2 ]                         
    then                                                                                                        # if false input in arguments than show info and exit program 
        cat ${installdir}/input_info.txt                    
        exit
fi

for i in ${unspecified_input[@]}; do unspecified_input1+=(`readlink -m $i`); done
for i in ${illumina_input[@]}; do illumina_input1+=(`readlink -m $i`); done
for i in ${ionTorrent_input[@]}; do ionTorrent_input1+=(`readlink -m $i`); done

input=(${unspecified_input1[@]} ${illumina_input1[@]} ${ionTorrent_input1[@]})                                  # summarize all inputdata in one array

if [[ -z $input ]]                                                                                              # if no input was given, then ...
    then
        echo -ne "You have no input file specified.\nPlease add an input file by using parameters -i/-t/-j and start the analysis again.\n\n"  # user info
        return 0
fi        

if [[ -z $projektname ]]                                                                                        # if projectname is empty, then ...
    then
        projektname=`echo $input[0] | sed 's/.*\///g' | sed 's/\.[^\.]*$//' | sed 's/_/-/g'`                    # use the first input as projectname and delete everything before backslash, all behind last point (file-format) and substitut "_" by "-"
    else
        projektname=`echo $projektname | sed 's/_/-/g'`                                                         # else tkae projectname and substitut "_" by "-" ("_" makes problems in final latex skript)
fi
arbeitsvz=${outputordner}/${projektname}                                                                        # define working directory with outputfolder and project


if [ -d ${arbeitsvz} ]      # If the working directory already exists, delete it
    then
        rm -r ${arbeitsvz}
fi


if [ -d ${arbeitsvz} ]
    then                                                                                                        # if working directory already exists
        echo -e "\n\nThere is already a Folder named '${arbeitsvz}' in your current output directory."          # show info to user and exit program
    #    echo -e "Please delete or rename this Folder and start the analysis again.\n" ; exit
    else
        mkdir -p ${arbeitsvz}/Ausgangsdateien                                                                   # creat working directory and directory "Ausgangsdateien"
        mkdir ${arbeitsvz}/MultiBlast                                                                           # create Blast directory                                                                          
        touch ${arbeitsvz}/AnBlast.txt                                                                          # create empty Anblast.txt file
        touch ${arbeitsvz}/TrimmedReads.fastq                                                                   # create empty fastq-file
        touch ${arbeitsvz}/ReadStatus.txt                                                                       # create empty ReadStatus-file
        touch ${arbeitsvz}/TrimStatus.txt                                                                       # create empty TrimStatus-file
        touch ${arbeitsvz}/asignedAcc.txt                                                                       # create empty file for assigned Reads
        touch ${arbeitsvz}/report.txt
        for i in ${input[*]}                                                                                    # for each element in tax, do ...
           do  
               echo $i | sed 's/.*\///g' | sed 's/\.[^\.]*$//' | sed 's/_/-/g' >> ${arbeitsvz}/input.txt        # write all input-file names in one file to use later in R
           done 
        mkdir -p ${arbeitsvz}/DB                                                                                # create a directory for all reference sequences
fi 
 
touch ${arbeitsvz}/input_report.txt                                                                             # create empty input_report.txt (needed for resultprotocol)
if [[ -n ${unspecified_input} ]]                                                                                # if unspecified input is specified, then...
    then
        echo "unspecified input" >> ${arbeitsvz}/input_report.txt                                               # write "unspecified input" to input_report.txt
fi

if [[ -n ${illumina_input} ]]                                                                                   # if illumina input is specified, then...
    then
        echo "illumina input" >> ${arbeitsvz}/input_report.txt                                                  # write "illumina input" to input_report.txt
fi

if [[ -n ${ionTorrent_input} ]]                                                                                 # if ionTorrent input is specified, then...
    then
        echo "ionTorrent input" >> ${arbeitsvz}/input_report.txt                                                # write "ionTorrent input" to input_report.txt
fi
    
if [[ ${copy} == y ]]                                                                                           # copy was set to yes
    then    
        cp ${input} ${arbeitsvz}/.                                                                              # copy input files
        ln -s ${input[*]} ${arbeitsvz}/Ausgangsdateien/.                                                        # and create a link
    else                                                        
        ln -s ${input[*]} ${arbeitsvz}/Ausgangsdateien/.                                                        # otherwise just create a link
fi  
   
mkdir -p ${refseqdir}                                                                                           # refseqdir from Config.txt, if it does not exist

if [[ -n ${tax} ]]                                                                                              # if taxids were given, then ...
    then    
        touch ${arbeitsvz}/tmp.txt                                                                              # create a temporary file
        for i in ${tax[*]}                                                                                      # for each element in tax, do ...
            do  
                echo $i >> ${arbeitsvz}/tmp.txt                                                                 # write given taxids into temporaty file
            done    
        sort ${arbeitsvz}/tmp.txt | uniq > ${arbeitsvz}/uniq-tids.txt                                           # sort taxids and make them uniq and write to new file
        rm ${arbeitsvz}/tmp.txt                                                                                 # remove temporary file
fi  

if ! [[ -s ${installdir}/gi_exclude.txt ]]                                                                      # check if there is already exists a list of gis to exclude from blast analysis, if not then...
    then
        . ${installdir}/TaxID_exclude.sh                                                                        # start TaxID_exclude.sh to create that list
fi
echo -ne "For progress please check ${arbeitsvz}/progress.log.\n\n"                                             # user info

. ${installdir}/metagenomanalyse_neu.sh 1>${arbeitsvz}/progress.log 2>${arbeitsvz}/error.log                    # start metagenomic workflow and save errors in error.log and terminal output in progress.log
