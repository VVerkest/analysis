#!/bin/bash                                                                                                                                                                                                                                                                   
source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.410
# source /opt/sphenix/core/bin/setup_local.sh /sphenix/user/ChengWei/sPH_dNdeta/build_folder/build_HR_ntuple/install

data_type=${1}
control_file=${2}
core_i=${3}
output_sub_folder=${4}
used_zvtx_folder=${5}
root -l -b -q ${control_file}/evtXY_mother.C\(${core_i},\"${output_sub_folder}\",\"${used_zvtx_folder}\"\)

# core_i=${1}

# for i in {0..19};
# do
#     # root -l -b -q /sphenix/user/ChengWei/INTT/INTT_dNdeta_repo/DST_MC/data20869_xy_geo1.C\(${i},${core_i}\)
#     root -l -b -q /sphenix/user/ChengWei/INTT/INTT_dNdeta_repo/DST_MC/control_files/MC_HIJING398/Geo_study/test_xy_398_000_geo.C\(${i},${core_i}\)
#     sleep 3
# done

# i=${1}

# root.exe -l -b -q /sphenix/user/ChengWei/sPH_dNdeta/dNdEta_Run2023/macros/Fun4All_G4_sPHENIX.C\(0,0,\"/sphenix/user/ChengWei/sPH_dNdeta/dNdEta_INTT_MC_388_000/INTTRecoClusters_test_400_$i.root\",\"/sphenix/user/ChengWei/sPH_dNdeta/dNdEta_Run2023/macros/list/dNdEta_INTT/list_folder_388_000/dst_INTTdNdEta_400_$i.list\",0\)
