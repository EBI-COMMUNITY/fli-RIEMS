#!/bin/bash
# Schnellere Grep-Funktion zur Zuordnung der viralen-Taxids aus den gesplitteten .dmp-Dateien über Case  
function GiExtraction()                                                                                       
{
  case "${line}" in                                                                                       # line: jeweilige Viridae-Gi
    1*) case "${line}" in                                                                                 # beginnt diese mit der 1, schau nach der
          1)    echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-1.dmp`\n" >> tmp.dmp ;;
          10*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-10.dmp`\n" >> tmp.dmp ;;          # 2. Ziffer und geh dann in die 
          11*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-11.dmp`\n" >> tmp.dmp ;;          # entsprechende Datei und grep sie dort
          12*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-12.dmp`\n" >> tmp.dmp ;;
          13*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-13.dmp`\n" >> tmp.dmp ;;
          14*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-14.dmp`\n" >> tmp.dmp ;;
          15*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-15.dmp`\n" >> tmp.dmp ;;
          16*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-16.dmp`\n" >> tmp.dmp ;;
          17*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-17.dmp`\n" >> tmp.dmp ;;
          18*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-18.dmp`\n" >> tmp.dmp ;;
          19*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-19.dmp`\n" >> tmp.dmp ;;
        esac ;;
    2*) case "${line}" in
          2)    echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-2.dmp`\n" >> tmp.dmp ;;
          20*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-20.dmp`\n" >> tmp.dmp ;;           
          21*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-21.dmp`\n" >> tmp.dmp ;;          
          22*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-22.dmp`\n" >> tmp.dmp ;;
          23*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-23.dmp`\n" >> tmp.dmp ;;
          24*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-24.dmp`\n" >> tmp.dmp ;;
          25*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-25.dmp`\n" >> tmp.dmp ;;
          26*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-26.dmp`\n" >> tmp.dmp ;;
          27*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-27.dmp`\n" >> tmp.dmp ;;
          28*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-28.dmp`\n" >> tmp.dmp ;;
          29*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-29.dmp`\n" >> tmp.dmp ;;
        esac ;;
    3*) case "${line}" in
          3)    echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-3.dmp`\n" >> tmp.dmp ;;
          30*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-30.dmp`\n" >> tmp.dmp ;;           
          31*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-31.dmp`\n" >> tmp.dmp ;;          
          32*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-32.dmp`\n" >> tmp.dmp ;;
          33*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-33.dmp`\n" >> tmp.dmp ;;
          34*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-34.dmp`\n" >> tmp.dmp ;;
          35*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-35.dmp`\n" >> tmp.dmp ;;
          36*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-36.dmp`\n" >> tmp.dmp ;;
          37*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-37.dmp`\n" >> tmp.dmp ;;
          38*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-38.dmp`\n" >> tmp.dmp ;;
          39*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-39.dmp`\n" >> tmp.dmp ;;
        esac ;;
    4*) case "${line}" in
          4)    echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-4.dmp`\n" >> tmp.dmp ;;
          40*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-40.dmp`\n" >> tmp.dmp ;;           
          41*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-41.dmp`\n" >> tmp.dmp ;;          
          42*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-42.dmp`\n" >> tmp.dmp ;;
          43*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-43.dmp`\n" >> tmp.dmp ;;
          44*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-44.dmp`\n" >> tmp.dmp ;;
          45*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-45.dmp`\n" >> tmp.dmp ;;
          46*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-46.dmp`\n" >> tmp.dmp ;;
          47*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-47.dmp`\n" >> tmp.dmp ;;
          48*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-48.dmp`\n" >> tmp.dmp ;;
          49*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-49.dmp`\n" >> tmp.dmp ;;
        esac ;;
    5*) case "${line}" in
          5)    echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-5.dmp`\n" >> tmp.dmp ;;
          50*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-50.dmp`\n" >> tmp.dmp ;;           
          51*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-51.dmp`\n" >> tmp.dmp ;;          
          52*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-52.dmp`\n" >> tmp.dmp ;;
          53*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-53.dmp`\n" >> tmp.dmp ;;
          54*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-54.dmp`\n" >> tmp.dmp ;;
          55*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-55.dmp`\n" >> tmp.dmp ;;
          56*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-56.dmp`\n" >> tmp.dmp ;;
          57*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-57.dmp`\n" >> tmp.dmp ;;
          58*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-58.dmp`\n" >> tmp.dmp ;;
          59*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-59.dmp`\n" >> tmp.dmp ;;
        esac ;;
    6*) case "${line}" in
          6)    echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-6.dmp`\n" >> tmp.dmp ;;
          60*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-60.dmp`\n" >> tmp.dmp ;;           
          61*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-61.dmp`\n" >> tmp.dmp ;;          
          62*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-62.dmp`\n" >> tmp.dmp ;;
          63*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-63.dmp`\n" >> tmp.dmp ;;
          64*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-64.dmp`\n" >> tmp.dmp ;;
          65*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-65.dmp`\n" >> tmp.dmp ;;
          66*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-66.dmp`\n" >> tmp.dmp ;;
          67*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-67.dmp`\n" >> tmp.dmp ;;
          68*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-68.dmp`\n" >> tmp.dmp ;;
          69*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-69.dmp`\n" >> tmp.dmp ;;
        esac ;;
    7*) case "${line}" in
          7)    echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-7.dmp`\n" >> tmp.dmp ;;
          70*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-70.dmp`\n" >> tmp.dmp ;;           
          71*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-71.dmp`\n" >> tmp.dmp ;;          
          72*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-72.dmp`\n" >> tmp.dmp ;;
          73*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-73.dmp`\n" >> tmp.dmp ;;
          74*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-74.dmp`\n" >> tmp.dmp ;;
          75*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-75.dmp`\n" >> tmp.dmp ;;
          76*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-76.dmp`\n" >> tmp.dmp ;;
          77*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-77.dmp`\n" >> tmp.dmp ;;
          78*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-78.dmp`\n" >> tmp.dmp ;;
          79*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-79.dmp`\n" >> tmp.dmp ;;
        esac ;;
    8*) case "${line}" in
          8)    echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-8.dmp`\n" >> tmp.dmp ;;
          80*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-80.dmp`\n" >> tmp.dmp ;;           
          81*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-81.dmp`\n" >> tmp.dmp ;;          
          82*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-82.dmp`\n" >> tmp.dmp ;;
          83*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-83.dmp`\n" >> tmp.dmp ;;
          84*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-84.dmp`\n" >> tmp.dmp ;;
          85*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-85.dmp`\n" >> tmp.dmp ;;
          86*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-86.dmp`\n" >> tmp.dmp ;;
          87*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-87.dmp`\n" >> tmp.dmp ;;
          88*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-88.dmp`\n" >> tmp.dmp ;;
          89*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-89.dmp`\n" >> tmp.dmp ;;
        esac ;;
    9*) case "${line}" in
          9)    echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-9.dmp`\n" >> tmp.dmp ;;
          90*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-90.dmp`\n" >> tmp.dmp ;;           
          91*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-91.dmp`\n" >> tmp.dmp ;;          
          92*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-92.dmp`\n" >> tmp.dmp ;;
          93*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-93.dmp`\n" >> tmp.dmp ;;
          94*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-94.dmp`\n" >> tmp.dmp ;;
          95*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-95.dmp`\n" >> tmp.dmp ;;
          96*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-96.dmp`\n" >> tmp.dmp ;;
          97*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-97.dmp`\n" >> tmp.dmp ;;
          98*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-98.dmp`\n" >> tmp.dmp ;;
          99*)  echo -e "`grep "^\<$line\>" taxid_acc_nucl/taxid_acc_nucl-99.dmp`\n" >> tmp.dmp ;;
        esac ;;
  esac
}

