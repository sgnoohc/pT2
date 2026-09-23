# Source CMS default setup
source /cvmfs/cms.cern.ch/cmsset_default.sh

# Check SCRAM_ARCH
if [ "$SCRAM_ARCH" != "el9_amd64_gcc12" ]; then
    echo "SCRAM_ARCH not known, please manually set CMSSW"
    exit 1
fi

# Set up CMSSW environment
cd /cvmfs/cms.cern.ch/el9_amd64_gcc12/cms/cmssw/CMSSW_15_1_0/src/ || exit 1
cmsenv
cd -
