#!/bin/bash

export USER="$(id -u -n)"
export LOGNAME=${USER}
export HOME=/sphenix/user/${LOGNAME}/analysis/JS-Jet/JetValidation/macro
export MYINSTALL=/sphenix/u/${USER}/install

source /opt/sphenix/core/bin/sphenix_setup.sh -n
source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL
 
# print the environment - needed for debugging
#printenv

#echo $dataFileList
#simFileList=/sphenix/user/bkimelman/analysis/CaloEmbedAnalysis/macro/fileLists/simLists/simDST_$(printf "%04d" $1).list
#geomList=/sphenix/user/bkimelman/simGlobalLists/simGlobalDST_$(printf "%04d" $1).list
#truthList=/sphenix/user/bkimelman/simTruthLists/simTruthDST_$(printf "%04d" $1).list
#dataFileList=/sphenix/user/bkimelman/Run23745_ana399_DSTs/Run23745_ana399_DST_$(printf "%04d" $1).list
#dataFileList=$HOME/dst_calo_list.list
#dataFileList=$HOME/lists/dst_calo_list_49219_$(printf "%04d" $1).list

#dataFileList=$HOME/lists/dst_calo_list_49219_$(printf "%04d" $1).list
#truthList=$HOME/lists/dst_calo_fitting_49219_$(printf "%04d" $1).list

#dataFileList="/sphenix/lustre01/sphnxpro/production/physics/run2pp/caloy2calib/ana450_2024p009/run_00049200_00049300/DST_CALO_run2pp_ana446_2024p007-00049216-$(printf "%05d" $1).root"
#truthList="/sphenix/lustre01/sphnxpro/production/physics/run2pp/caloy2fitting/ana446_2024p007/run_00049200_00049300/DST_CALOFITTING_run2pp_ana446_2024p007-00049219-$(printf "%05d" $1).root"

jetcaloList="/sphenix/lustre01/sphnxpro/production/run2pp/physics/ana462_2024p010_v001/DST_JETCALO/run_00049200_00049300/dst/DST_JETCALO_run2pp_ana462_2024p010_v001-00049219-$(printf "%05d" $1).root"
jetList="/sphenix/lustre01/sphnxpro/production/run2pp/physics/ana462_2024p010_v001/DST_JETCALO/run_00049200_00049300/dst/DST_JET_run2pp_ana462_2024p010_v001-00049219-$(printf "%05d" $1).root"
#geomList="/sphenix/lustre01/sphnxpro/production/physics/run2pp/caloy2fitting/ana446_2024p007/run_00049200_00049300/DST_CALOFITTING_run2pp_ana446_2024p007-00049219-$(printf "%05d" $1).root"
caloList="/sphenix/lustre01/sphnxpro/production/run2pp/physics/ana462_2024p010_v001/DST_CALO/run_00049200_00049300/dst/DST_CALO_run2pp_ana446_2024p007-00049219-$(printf "%05d" $1).root"



#dataFileList=$HOME/dst_calo_run2pp-00042183.list
#dataFileList="/sphenix/lustre01/sphnxpro/physics/slurp/caloy2calib/ana437_2024p007/run_00053300_00053400/DST_CALO_run2pp_ana437_2024p007-00053376-$(printf "%05d" $1).root"
#dataFileList="/sphenix/lustre01/sphnxpro/physics/slurp/caloy2calib/ana437_2024p007/run_00053300_00053400/DST_CALO_run2pp_ana437_2024p007-00053376-$(printf "%05d" $1).root"

#echo simFileList: $simFileList
#echo geomList: $geomList
echo jetcaloList: $jetcaloList
echo jetList: $jetList
echo caloList: $caloList

 # use this for sim
#root.exe -q -b Fun4All_EMJetVal.C\(\"$truthList\",\"$simFileList\",\"$geomList\",\"/sphenix/u/jamesj3j3/analysis/JS-Jet/JetValidation/macro/condorTest/EMTree_$(printf "%04d" $1).root\"\)

 
 # use this for data

#root.exe -q -b Fun4All_JetVal.C\(\"\",\"$dataFileList\",\"\",\"/sphenix/user/$USER/analysis/JS-Jet/JetValidation/macro/Run23745_ana399_DST_$(printf "%04d" $1).root\"\)

#root.exe -q -b Fun4All_JetVal.C\(\"\",\"$dataFileList\",\"\",\"/sphenix/user/$USER/analysis/JS-Jet/JetValidation/macro/out/JetVal_run2pp_ana437_2024p007_00049270_$(printf "%04d" $1).root\"\)

root.exe -q -b Fun4All_JetVal.C\(\"$jetcaloList\",\"$jetList\",\"$caloList\",\"/sphenix/user/$USER/analysis/JS-Jet/JetValidation/macro/out/JetVal_run2pp_ana462_2024p007_00049219_$(printf "%04d" $1).root\"\)


#root.exe -q -b Fun4All_JetVal.C\(\"$jetcaloList\",\"$dataFileList\",\"$geomList\",\"$rawList\",\"/sphenix/user/$USER/analysis/JS-Jet/JetValidation/macro/out/JetVal_run2pp_ana446_2024p007_00049219_$(printf "%04d" $1).root\"\)

#root.exe -q -b Fun4All_JetVal.C\(\"$truthlist\",\"\",\"\",\"/sphenix/user/$USER/analysis/JS-Jet/JetValidation/macro/out/JetVal_run2pp_ana437_2024p007_00049219_$(printf "%04d" $1).root\"\)
#root.exe -q -b Fun4All_JetVal.C\(\"\",\"$dataFileList\",\"\",\"/sphenix/user/$USER/analysis/JS-Jet/JetValidation/macro/out/JetVal_run2pp_ana437_2024p007_00042183_$(printf "%04d" $1).root\"\)

echo all done
