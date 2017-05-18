import os
from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'cascadev2ppb2013InvMass'
config.General.workArea = 'cascadev2ppb2013InvMass'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
#config.JobType.pluginName = 'PrivateMC'
config.JobType.psetName = os.path.expandvars('$CMSSW_BASE/src/XiAnalyzer/XiAnalyzer/test/xianalysis_cfg.py')

config.section_("Data")
config.Data.inputDataset = '/PAMinBiasUPC/davidlw-PA2013_FlowCorr_PromptReco_MB_Gplus_Rereco_ReTracking_v18-25c9a89be536a41c8ccb3c75e9fd9358/USER'
#config.Data.inputDataset = '/MinBias_TuneZ2_7TeV-pythia6/davidlw-Skim_v15-1709136df11f28c8f8a8a944d51c46a6/USER'
#config.Data.userInputFiles = list(open('HMMC90.txt'))
config.Data.inputDBS = 'phys03'
#config.Data.primaryDataset = 'MinBias_TuneZ2star_7TeV_pythia6'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
config.Data.outLFNDirBase = '/store/group/phys_heavyions/btran/'
config.Data.publication = False
config.Data.publishDBS = 'phys03'
config.Data.outputDatasetTag = 'pPb_Cascade_Rereco_MB'

config.section_("Site")
config.Site.storageSite = 'T2_CH_CERN'
config.Site.whitelist = ['T2_US_MIT']

