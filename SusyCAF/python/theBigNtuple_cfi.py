import FWCore.ParameterSet.Config as cms

from SUSYBSMAnalysis.SusyCAF.SusyCAF_Event_cfi import *
from SUSYBSMAnalysis.SusyCAF.SusyCAF_MET_cfi import *
from SUSYBSMAnalysis.SusyCAF.SusyCAF_Electron_cfi import *

susyTree = cms.EDAnalyzer("SusyTree",
                             outputCommands = cms.untracked.vstring(
    'drop *',
    'keep *_susycafevent_*_*',
    'keep *_susycafmet_*_*',
    'keep *_susycafelectron_*_*',
    ))

theBigNtuple = cms.Sequence( (susycafevent +
                              susycafmet +
                              susycafelectron ) *
                              susyTree
                             )
