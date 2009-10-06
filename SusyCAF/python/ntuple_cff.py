import FWCore.ParameterSet.Config as cms

from SUSYBSMAnalysis.SusyCAF.theBigNtuple_cfi import *

#Schedule
ntuplize_step = cms.Path( theBigNtuple )
schedule = cms.Schedule( ntuplize_step )
