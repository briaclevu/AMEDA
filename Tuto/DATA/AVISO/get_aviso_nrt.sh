#!/bin/bash
# get the daily AVISO merged dataset to produce nc file

USER=*              #This is the FTP user that has access to the server.
PASS=*           #This is the password for the FTP user.
HOST2=nrt.cmems-du.eu   #Host for NRT

n=$(cat date.txt | wc -l)
yyyymmdd1=$(sed -n "1 p" date.txt)
yyyymmdd2=$(sed -n "$n p" date.txt)

mkdir phy_l4
cd phy_l4
mkdir allsat
cd allsat

while read -r line
do
    
    YEAR=`echo $line | cut -b1-4`
    MONTH=`echo $line | cut -b5-6`
  
    hostdir2=Core/SEALEVEL_EUR_PHY_L4_NRT_OBSERVATIONS_008_060/dataset-duacs-nrt-europe-merged-allsat-phy-l4/${YEAR}/${MONTH}
    
    # remove old file
    rm nrt_europe_allsat_phy_l4_${line}*

    # get the daily NRT file from CMEMS
    wget -r -nd ftp://$USER:$PASS@${HOST2}/${hostdir2}/nrt_europe_allsat_phy_l4_${line}_*.nc
                
done < ../../date.txt

cdo -f nc copy nrt_europe_allsat_phy_l4_${YEAR}* ../nrt_europe_allsat_phy_l4_${yyyymmdd1}_${yyyymmdd2}.nc

cd ../..
