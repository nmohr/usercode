#
#  SUSY-PAT configuration fragment
#
#  PAT configuration for the SUSY group - 37X series
#  More information here:
#  https://twiki.cern.ch/twiki/bin/view/CMS/SusyPatLayer1DefV9


import FWCore.ParameterSet.Config as cms

def addDefaultSUSYPAT(process, mcInfo=True, HLTMenu='HLT', JetMetCorrections='Spring10', mcVersion='' ,theJetNames = ['IC5Calo','AK5JPT'],doValidation=False, extMatch=True):
    loadPF2PAT(process,mcInfo,JetMetCorrections,'PF')
    if not mcInfo:
	removeMCDependence(process)
    loadMCVersion(process,mcVersion,mcInfo)
    loadPAT(process,JetMetCorrections)
    addJetMET(process,theJetNames,mcVersion)
    loadPATTriggers(process,HLTMenu)

    #-- Counter for the number of processed events --------------------------------
    process.eventCountProducer = cms.EDProducer("EventCountProducer")

    # Full path
    process.susyPatDefaultSequence = cms.Sequence( process.eventCountProducer 
                                                   * process.patDefaultSequence * process.patPF2PATSequencePF
                                                   #* process.patTrigger * process.patTriggerEvent
                                                    )
    if mcInfo and extMatch:
        extensiveMatching(process)
        process.susyPatDefaultSequence.replace(process.patDefaultSequence, process.extensiveMatching+process.patDefaultSequence)

    if mcVersion == '35x' and 'JPT' in ''.join(theJetNames): 
    	process.susyPatDefaultSequence.replace(process.eventCountProducer, process.eventCountProducer * process.recoJPTJets)
    if doValidation:
        loadSusyValidation(process)
        process.susyPatDefaultSequence.replace(process.patPF2PATSequencePF, process.patPF2PATSequencePF * process.ak5CaloJetsL2L3 * process.metJESCorAK5CaloJet  * process.RecoSusyValidation * process.PatSusyValidation*process.MEtoEDMConverter)

def extensiveMatching(process):
    process.load("SimGeneral.TrackingAnalysis.trackingParticlesNoSimHits_cfi")    # On RECO
    process.load("SimMuon.MCTruth.MuonAssociatorByHits_cfi")  # On RECO
    process.mergedTruth = cms.EDProducer("GenPlusSimParticleProducer",
        src           = cms.InputTag("g4SimHits"), # use "famosSimHits" for FAMOS
        setStatus     = cms.int32(5),             # set status = 8 for GEANT GPs
        filter        = cms.vstring("pt > 3.0"),  # just for testing (optional)
        genParticles   = cms.InputTag("genParticles") # original genParticle list
    )
    process.load("MuonAnalysis.MuonAssociators.muonClassificationByHits_cfi")

    from MuonAnalysis.MuonAssociators.muonClassificationByHits_cfi import addUserData as addClassByHits
    addClassByHits(process.patMuons,labels=['classByHitsGlb'],extraInfo=True)
    addClassByHits(process.patMuonsPF,labels=['classByHitsGlb'],extraInfo=True)
    
    #process.load("ElectronAnalysis.ElectronAssociators.electronClassificationByHits_cfi")
    #from ElectronAnalysis.ElectronAssociators.electronClassificationByHits_cfi import addUserData as addEClassByHits
    #addEClassByHits(process.patElectrons,labels=['classByHitsGsf'],extraInfo=True)
    #addEClassByHits(process.patElectronsPF,labels=['classByHitsGsf'],extraInfo=True)

    process.extensiveMatching = cms.Sequence(process.mergedTruth+process.muonClassificationByHits) #+process.electronClassificationByHits)


def loadMCVersion(process, mcVersion, mcInfo):
    #-- To be able to run on 35X input samples ---------------------------------------
    from PhysicsTools.PatAlgos.tools.cmsswVersionTools import run36xOn35xInput
    if not mcVersion:
	return
    elif mcVersion == '35x': 
	run36xOn35xInput(process)
	if mcInfo:
		run36xOnReRecoMC(process)
    	#-- Jet plus tracks are in RECO in 36X, but not in 35X-----------------------
	process.load("RecoJets.Configuration.RecoJPTJets_cff")
    else: raise ValueError, "Unknown MC version: %s" % (mcVersion)


def loadPAT(process,JetMetCorrections):
    #-- Changes for electron and photon ID ----------------------------------------
    # Turn off photon-electron cleaning (i.e., flag only)
    process.cleanPatPhotons.checkOverlaps.electrons.requireNoOverlaps = False

    # Remove embedding of superClusters, will keep entire superCluster collection
    process.patElectrons.embedSuperCluster = False
    process.patPhotons.embedSuperCluster   = False
    
    #-- Tuning of Monte Carlo matching --------------------------------------------
    # Also match with leptons of opposite charge
    process.electronMatch.checkCharge = False
    process.electronMatch.maxDeltaR   = cms.double(0.2)
    process.electronMatch.maxDPtRel   = cms.double(999999.)
    process.muonMatch.checkCharge     = False
    process.muonMatch.maxDeltaR       = cms.double(0.2)
    process.muonMatch.maxDPtRel       = cms.double(999999.)
    process.tauMatch.checkCharge      = False
    process.tauMatch.maxDeltaR        = cms.double(0.3)
    process.patJetPartonMatch.maxDeltaR  = cms.double(0.25)
    process.patJetPartonMatch.maxDPtRel  = cms.double(999999.)
    process.patJetGenJetMatch.maxDeltaR  = cms.double(0.25)
    process.patJetGenJetMatch.maxDPtRel  = cms.double(999999.)

    #-- Jet corrections -----------------------------------------------------------
    process.patJetCorrFactors.corrSample = JetMetCorrections 

def loadPF2PAT(process,mcInfo,JetMetCorrections,postfix):
    #-- PAT standard config -------------------------------------------------------
    process.load("PhysicsTools.PatAlgos.patSequences_cff")
    #-- Jet corrections -----------------------------------------------------------
    process.patJetCorrFactors.corrSample = JetMetCorrections 
    #-- PF2PAT config -------------------------------------------------------------
    from PhysicsTools.PatAlgos.tools.pfTools import usePF2PAT
    usePF2PAT(process,runPF2PAT=True, jetAlgo='AK5',runOnMC=mcInfo,postfix=postfix)
    process.patJetsPF.embedGenJetMatch = False
    process.patJetsPF.embedPFCandidates = False
    #-- Relax isolation -----------------------------------------------------------
    process.pfNonIsolatedElectronsPF = process.pfIsolatedElectronsPF.clone(combinedIsolationCut = 3.)
    process.pfNonIsolatedMuonsPF = process.pfIsolatedMuonsPF.clone(combinedIsolationCut = 3.)
    process.pfElectronSequencePF.replace(process.pfIsolatedElectronsPF,process.pfNonIsolatedElectronsPF+process.pfIsolatedElectronsPF)
    process.pfMuonSequencePF.replace(process.pfIsolatedMuonsPF,process.pfNonIsolatedMuonsPF+process.pfIsolatedMuonsPF)
    process.pfIsolatedMuonsPF.combinedIsolationCut = 0.25
    process.pfIsolatedElectronsPF.combinedIsolationCut = 0.3
    process.patElectronsPF.pfElectronSource = "pfNonIsolatedElectronsPF"
    process.patMuonsPF.pfMuonSource = "pfNonIsolatedMuonsPF"
    process.electronMatchPF.checkCharge = False
    process.muonMatchPF.checkCharge = False

def loadPATTriggers(process,HLTMenu):
    #-- Trigger matching ----------------------------------------------------------
    from PhysicsTools.PatAlgos.tools.trigTools import switchOnTrigger
    switchOnTrigger( process )
    #process.patTriggerSequence.remove( process.patTriggerMatcher )
    #process.patTriggerEvent.patTriggerMatches  = []
    # If we have to rename the default trigger menu
    process.patTrigger.processName = HLTMenu
    process.patTriggerEvent.processName = HLTMenu

def addSUSYJetCollection(process,jets = 'IC5Calo',mcVersion='',doJTA=False,doType1MET=False,doJetID=True,jetIdLabel=None):
    if mcVersion == '35x':
        from PhysicsTools.PatAlgos.tools.cmsswVersionTools import addJetCollection35X as addJetCollection
    else:
	from PhysicsTools.PatAlgos.tools.jetTools import addJetCollection
    algorithm = jets[0:3]
    type = jets[3:len(jets)]
    jetCorrLabel = (algorithm,type)
    if 'IC' in algorithm: collection = algorithm.replace('IC','iterativeCone')
    elif 'SC' in algorithm: collection = algorithm.replace('SC','sisCone')
    elif 'AK' in algorithm: collection = algorithm.replace('AK','ak')
    elif 'KT' in algorithm: collection = algorithm.replace('KT','kt')
    else: raise ValueError, "Unknown jet algorithm: %s" % (jets)
    jetIdLabel = algorithm.lower()
    if type == 'Calo':
	jetCollection = '%(collection)sCaloJets' % locals()
        doJTA = True
	if not 'AK7' in algorithm:
		doType1MET = True
    elif type == 'PF':
	jetCollection = '%(collection)sPFJets' % locals()
        doJTA = True
	doJetID = False
    elif type == 'JPT':
        if 'IC' in algorithm: collectionJPT = algorithm.replace('IC','Icone')
        elif 'SC' in algorithm: collectionJPT = algorithm.replace('SC','Siscone')
        elif 'AK' in algorithm: collectionJPT = algorithm.replace('AK','AntiKt')
        else: raise ValueError, "Unknown jet algorithm: %s" % (jets)
        jetCollection = 'JetPlusTrackZSPCorJet%(collectionJPT)s' % locals()
    elif type == 'Track':
	jetCollection = '%(collection)sTrackJets' % locals()
    	jetCorrLabel = None
	doJetID = False
    else: raise ValueError, "Unknown jet type: %s" % (jets)

    addJetCollection(process, cms.InputTag(jetCollection),
                     algorithm, type,
                     doJTA            = doJTA,
                     doBTagging       = True,
                     jetCorrLabel     = jetCorrLabel,
                     doType1MET       = doType1MET,
                     doL1Cleaning     = True,
                     doL1Counters     = True,
                     doJetID          = doJetID,
		     jetIdLabel       = jetIdLabel,
                     genJetCollection = cms.InputTag('%(collection)sGenJets' % locals())
                     )

def addJetMET(process,theJetNames,mcVersion):
    #-- Extra Jet/MET collections -------------------------------------------------
    # Add a few jet collections...
    for jetName in theJetNames:
    	addSUSYJetCollection(process,jetName,mcVersion)
    
    #-- Tune contents of jet collections  -----------------------------------------
    theJetNames.append('')
    for jetName in theJetNames:
        module = getattr(process,'patJets'+jetName)
        module.addTagInfos = False    # Remove tag infos
        module.embedGenJetMatch = False # Only keep reference, since we anyway keep the genJet collections
        #module.embedCaloTowers = True # To drop calo towers
    theJetNames.pop()
    
    # Add tcMET
    from PhysicsTools.PatAlgos.tools.metTools import addTcMET #, addPfMET
    addTcMET(process,'TC')
    #addPfMET(process,'PF') #is in PF2PAT

    # Rename default jet collection for uniformity
    process.cleanPatJetsAK5Calo = process.cleanPatJets
    process.patMETsAK5Calo      = process.patMETs
    process.patMHTsAK5Calo      = process.patMHTs

    # Modify subsequent modules
    process.patHemispheres.patJets = process.cleanPatJetsAK5Calo.label()
    process.countPatJets.src       = process.cleanPatJetsAK5Calo.label()
    
    # Add MHT (inserted until officially suported)
    from PhysicsTools.PatAlgos.producersLayer1.mhtProducer_cff import makePatMHTs, patMHTs
    process.countPatCandidates.replace(process.countPatJets, process.countPatJets + process.makePatMHTs)
    process.patMHTs.jetTag      = 'patJets'
    process.patMHTs.electronTag = 'patElectrons'
    process.patMHTs.muonTag     = 'patMuons'
    process.patMHTs.tauTag      = 'patTaus'
    process.patMHTs.photonTag   = 'patPhotons'

    # Modify counters' input
    process.patCandidateSummary.candidates.remove(cms.InputTag('patMETs'))
    process.patCandidateSummary.candidates.append(cms.InputTag('patMETsAK5Calo'))
    process.patCandidateSummary.candidates.append(cms.InputTag('patMHTsAK5Calo'))
    process.cleanPatCandidateSummary.candidates.remove(cms.InputTag('cleanPatJets'))
    process.cleanPatCandidateSummary.candidates.append(cms.InputTag('cleanPatJetsAK5Calo'))
    # Add new jet collections to counters (MET done automatically)
    for jets in theJetNames: 
        process.patCandidateSummary.candidates.append(cms.InputTag('patJets'+jets))
        process.selectedPatCandidateSummary.candidates.append(cms.InputTag('selectedPatJets'+jets))
        process.cleanPatCandidateSummary.candidates.append(cms.InputTag('cleanPatJets'+jets))
	

def removeMCDependence( process ):
    #-- Remove MC dependence ------------------------------------------------------
    from PhysicsTools.PatAlgos.tools.coreTools import removeMCMatching
    removeMCMatching(process, ['All'])

def loadSusyValidation(process):
    process.load("JetMETCorrections.Configuration.JetCorrectionProducers_cff")
    process.load("DQM.Physics.susyValidation_cfi")
    process.load("DQMServices.Components.MEtoEDMConverter_cfi")
    process.load("DQMServices.Core.DQM_cfg")
    process.load("DQMServices.Components.DQMEnvironment_cfi")
    process.DQMStore = cms.Service("DQMStore")
    process.DQMStore.collateHistograms = cms.untracked.bool(True)
    process.options = cms.untracked.PSet(
 	fileMode = cms.untracked.string('NOMERGE')
    )

def getSUSY_pattuple_outputCommands( process ):
	from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent, patExtraAodEventContent, patTriggerEventContent, patTriggerStandAloneEventContent, patEventContentTriggerMatch
	keepList = []
    	susyAddEventContent = [ # PAT Objects
	# Keep PF2PAT output
        'keep *_patTaus_*_*',         
        'keep *_selectedPatMuonsPF_*_*',         
        'keep *_selectedPatElectronsPF_*_*',         
        'keep *_selectedPatTausPF_*_*',         
        'keep *_selectedPatJetsPF_*_*',
	#L1 trigger info         
	'keep L1GlobalTriggerObjectMapRecord_*_*_*',
        'keep L1GlobalTriggerReadoutRecord_*_*_*',
        # Generator information
        'keep recoGenJets_*GenJets*_*_*',
        'keep recoGenMETs_*_*_*',
	#Number of processed events
        'keep edmMergeableCounter_eventCountProducer_*_*',
	'keep recoRecoChargedRefCandidates_trackRefsForJets_*_*',
	'keep recoTrackJets_ak5TrackJets_*_*',
	'keep *_electronMergedSeeds_*_*',
	'keep *_Conversions_*_*',
	'keep recoPFCandidates_particleFlow_*_*',
        'keep recoSuperClusters_corrected*_*_*',
	'keep recoSuperClusters_pfElectronTranslator_*_*',
        'keep *_gsfElectronCores_*_*',    #Keep electron core
        'keep *_photonCore_*_*',        #Keep electron core
        'keep recoConversions_conversions_*_*',
        'keep recoTracks_*onversions_*_*',
        'keep HcalNoiseSummary_*_*_*', #Keep the one in RECO
	'keep *BeamHaloSummary_*_*_*',
	#DQM
	'keep *_MEtoEDMConverter_*_PAT',
	#'drop *_towerMaker_*_*'
        ] 
	keepList.extend(patEventContent)
	keepList.extend(patExtraAodEventContent)
	keepList.extend(patTriggerEventContent)
	keepList.extend(patEventContentTriggerMatch)
	keepList.extend(susyAddEventContent)
	return keepList


def run36xOnReRecoMC( process ):
    """
    ------------------------------------------------------------------
    running GenJets for ak5 and ak7

    process : process
    genJets : which gen jets to run
    ------------------------------------------------------------------    
    """
    print "*********************************************************************"
    print "NOTE TO USER: when running on 31X samples re-recoed in 3.5.6         "
    print "              with this CMSSW version of PAT                         "
    print "              it is required to re-run the GenJet production for     "
    print "              anti-kT since that is not part of the re-reco          "
    print "*********************************************************************"
    process.load("RecoJets.Configuration.GenJetParticles_cff")
    process.load("RecoJets.JetProducers.ak5GenJets_cfi")
    process.ak7GenJets = process.ak5GenJets.clone( rParam = 0.7 )
    process.makePatJets.replace( process.patJetCharge, process.genParticlesForJets+process.ak5GenJets+process.ak7GenJets+process.patJetCharge)
    #-- Remove changes for GenJets ------------------------------------------------
    process.genParticlesForJets.ignoreParticleIDs = cms.vuint32(1000022, 1000012, 1000014, 1000016, 2000012,
        2000014, 2000016, 1000039, 5100039, 4000012,
        4000014, 4000016, 9900012, 9900014, 9900016,
        39)
    process.genParticlesForJets.excludeResonances = True


