#!/bin/bash
# get the daily AVISO merged dataset to produce nc file

USER=*              #This is the FTP user that has access to the server.
PASS=*           #This is the password for the FTP user.
HOST1=my.cmems-du.eu    #Host for DT
HOST2=nrt.cmems-du.eu   #Host for NRT

n=$(cat date.txt | wc -l)
yyyymmdd1=$(sed -n "1 p" date.txt)
yyyymmdd2=$(sed -n "$n p" date.txt)

mkdir AVISO
cd AVISO
mkdir phy_l4
cd phy_l4
mkdir allsat
cd allsat

while read -r line
do
    
    YEAR=`echo $line | cut -b1-4`
    MONTH=`echo $line | cut -b5-6`
  
    hostdir1=Core/SEALEVEL_MED_PHY_L4_REP_OBSERVATIONS_008_051/dataset-duacs-rep-medsea-merged-allsat-phy-l4/${YEAR}/${MONTH}
    hostdir2=Core/SEALEVEL_MED_PHY_L4_NRT_OBSERVATIONS_008_050/dataset-duacs-nrt-medsea-merged-allsat-phy-l4/${YEAR}/${MONTH}
    
    # remove old file
    rm dt_med_allsat_phy_l4_${line}*

    # get the daily DT file from CMEMS
     wget -r -nd ftp://$USER:$PASS@${HOST1}/${hostdir1}/dt_med_allsat_phy_l4_${line}_*.nc

    if [ ! -f "dt_med_allsat_phy_l4_${line}*" ]; then

        # remove old file
        rm nrt_med_allsat_phy_l4_${line}*

        # get the daily NRT file from CMEMS
        wget -r -nd ftp://$USER:$PASS@${HOST2}/${hostdir2}/nrt_med_allsat_phy_l4_${line}_*.nc
                
    fi

done < ../../../date.txt

/opt/cdo-1.6.3/bin/cdo -f nc copy dt_med_allsat_phy_l4_${YEAR}* ../dt_med_allsat_phy_l4_${yyyymmdd1}_${yyyymmdd2}.nc
/opt/cdo-1.6.3/bin/cdo -f nc copy nrt_med_allsat_phy_l4_${YEAR}* ../nrt_med_allsat_phy_l4_${yyyymmdd1}_${yyyymmdd2}.nc

cd ../../..
