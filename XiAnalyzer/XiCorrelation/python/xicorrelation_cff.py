import FWCore.ParameterSet.Config as cms

#from RiceHig.V0Analysis.xicorrelation_cfi import *
# For interactive
from XiAnalyzer.XiCorrelation.xicorrelation_cfi import *



xiCorrelator = xiCorrelation.clone();
hltHM = hltHM.clone();
