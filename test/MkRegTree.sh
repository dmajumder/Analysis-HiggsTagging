#!/bin/bash

eos=/afs/cern.ch/project/eos/installation/cms/bin/eos.select
batchdir=${PWD}
workdir=/afs/cern.ch/work/d/devdatta/CMSREL/Analysis/CMSSW_8_0_20/src/Analysis/HiggsTagging/test
cd ${workdir}
eval `scramv1 runtime -sh`
cd ${batchdir}
echo $PWD
cp ${workdir}/MkRegTree.py ${workdir}/inputFiles.py ${batchdir}
cd ${batchdir}
python MkRegTree.py 
${eos} cp TrainingTree_BGToHH4b_narrow_B2GAnaFW_v80x_v2p2.root /eos/cms/store/user/devdatta/B2GEDMNTuples/CMSSW_8_0_X_V2_AK8Reg/TrainingTree_BGToHH4b_narrow_B2GAnaFW_v80x_v2p2.root
rm -f TrainingTree_BGToHH4b_narrow_B2GAnaFW_v80x_v2p2.root MkRegTree.py inputFiles.py
