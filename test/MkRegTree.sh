#!/bin/bash

batchdir=${PWD}
workdir=/afs/cern.ch/work/d/devdatta/CMSREL/Analysis/CMSSW_8_0_20/src/Analysis/HiggsTagging/test
cd ${workdir}
eval `scramv1 runtime -sh`
export PATH=/afs/cern.ch/user/d/devdatta/.local/bin:/afs/cern.ch/cms/lumi/brilconda-1.0.1/bin:${PATH}
export PYTHONPATH=/afs/cern.ch/user/d/devdatta/bin/rootpy:${PYTHONPATH}
echo "PYTHONPATH", $PYTHONPATH
cd ${batchdir}
echo $PWD
python ${workdir}/MkRegTree.py
eos cp TrainingTree_BGToHH4b_narrow_B2GAnaFW_v80x_v2p2.root /eos/cms/store/user/devdatta/B2GEDMNTuples/CMSSW_8_0_X_V2_AK8Reg/TrainingTree_BGToHH4b_narrow_B2GAnaFW_v80x_v2p2.root
