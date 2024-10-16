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
#globalList=/sphenix/user/bkimelman/simGlobalLists/simGlobalDST_$(printf "%04d" $1).list
#truthList=/sphenix/user/bkimelman/simTruthLists/simTruthDST_$(printf "%04d" $1).list
#dataFileList=/sphenix/user/bkimelman/Run23745_ana399_DSTs/Run23745_ana399_DST_$(printf "%04d" $1).list
dataFileList=$HOME/dst_calo_list.list
#dataFileList="/sphenix/lustre01/sphnxpro/physics/slurp/caloy2calib/ana437_2024p007/run_00053300_00053400/DST_CALO_run2pp_ana437_2024p007-00053376-$(printf "%05d" $1).root"
#dataFileList="/sphenix/lustre01/sphnxpro/physics/slurp/caloy2calib/ana437_2024p007/run_00053300_00053400/DST_CALO_run2pp_ana437_2024p007-00053376-"$(printf "%05d" $1).root

#echo simFileList: $simFileList
#echo globalList: $globalList
#echo truthList: $truthList
echo dataFileList: $dataFileList

 # use this for sim
#root.exe -q -b Fun4All_EMJetVal.C\(\"$truthList\",\"$simFileList\",\"$globalList\",\"/sphenix/u/jamesj3j3/analysis/JS-Jet/JetValidation/macro/condorTest/EMTree_$(
#printf "%04d" $1).root\"\)

 
 # use this for data

#root.exe -q -b Fun4All_JetVal.C\(\"\",\"$dataFileList\",\"\",\"/sphenix/user/$USER/analysis/JS-Jet/JetValidation/macro/Run23745_ana399_DST_$(printf "%04d" $1).root\"\)
root.exe -q -b Fun4All_JetVal.C\(\"\",\"$dataFileList\",\"\",\"/sphenix/user/$USER/analysis/JS-Jet/JetValidation/macro/out/JetVal_run2pp_ana437_2024p007_00053049_$(printf "%04d" $1).root\"\)

echo all done
