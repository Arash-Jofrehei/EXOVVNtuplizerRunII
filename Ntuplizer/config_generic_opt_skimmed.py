#cmsRun  config_generic_opt_skimmed.py  RunPeriod="Fall17" # for MC
#cmsRun  config_generic_opt_skimmed.py  RunPeriod="Run2017B" # for Data from 2017
#cmsRun  config_generic_opt_skimmed.py  RunPeriod="Run2018B" # for Data from 2018
#cmsRun  config_generic_opt_skimmed.py  RunPeriod="Autumn18" # for MC from 2018
#cmsRun  config_generic_opt_skimmed.py  RunPeriod="Summer16" # for MC from 2016


###### Process initialization ##########

import FWCore.ParameterSet.Config as cms

process = cms.Process("Ntuple")

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.Geometry.GeometryRecoDB_cff')



# Add PbPb centrality
process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.load('RecoHI.HiCentralityAlgos.HiCentrality_cfi')
process.hiCentrality.produceHFhits = False
process.hiCentrality.produceHFtowers = False
process.hiCentrality.produceEcalhits = False
process.hiCentrality.produceZDChits = True
process.hiCentrality.produceETmidRapidity = False
process.hiCentrality.producePixelhits = False
process.hiCentrality.produceTracks = False
process.hiCentrality.producePixelTracks = False
process.hiCentrality.reUseCentrality = True
process.hiCentrality.srcZDChits = cms.InputTag("QWzdcreco")
process.hiCentrality.srcReUse = cms.InputTag("hiCentrality","","reRECO")
process.centralityBin.Centrality = cms.InputTag("hiCentrality")
process.centralityBin.centralityVariable = cms.string("HFtowers")
process.centralityBin.nonDefaultGlauberModel = cms.string("")
process.cent_seq = cms.Sequence(process.hiCentrality * process.centralityBin)




process.TFileService = cms.Service("TFileService",
                                    fileName = cms.string('flatTuple.root')
#                                    fileName = cms.string('flatTuple_mc.root')
                                   )

#from EXOVVNtuplizerRunII.Ntuplizer.ntuplizerOptions_data_cfi import config
from EXOVVNtuplizerRunII.Ntuplizer.ntuplizerOptions_generic_cfi import config

# change from its original value
#config["DZCUT"] = 0.25
#config["FSIGCUT"] = 3
#config["VPROBCUT"] = 0.1
#config["DNNCUT"] = 0.2

				   
####### Config parser ##########

import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing ('analysis')

options.register( 'RunPeriod',
                  '',
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.string,          # string, int, or float
                  "RunNumber (Default Run2017B)")

options.register( 'runUpToEarlyF',
                  'false',
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list                                                                                                                                 
                  VarParsing.VarParsing.varType.bool,          # string, int, or float                                                                                                                        
                  "false")# https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideAboutPythonConfigFile

####



#options.maxEvents = 2000
options.maxEvents = 10000

#2015 data file
options.inputFiles = [
  'root://cms-xrd-global.cern.ch//store/hidata/HIRun2018A/HIForward/AOD/04Apr2019-v1/00000/19CFC2EF-7F8F-194C-92D9-01E0F59802AD.root' # 2018 data
#  '/store/hidata/HIRun2018A/HIForward/AOD/04Apr2019-v1/00000/00DB29C1-E3F8-6546-BB32-120BFD2C7FEA.root' # 2018 data
#  'file:/work/ajofrehe/gtau/AOD/CMSSW_10_3_2/src/EXOVVNtuplizerRunII/Ntuplizer/HIN-HINPbPbWinter16DR-00315.root'
#  '/store/user/ajofrehe/gtau_v01/ggTauTau_TuneCUEP_5p02TeV_MG5_aMCatNLO_pythia8-pLHE/ggTauTau_TuneCUEP_5p02TeV_MG5_aMCatNLO_pythia8-AODSIM-v2/210316_105241/0000/HIN-HINPbPbWinter16DR-00315_1.root'
#  '/store/user/ajofrehe/gtau_v01/ggTauTau_TuneCUEP_5p02TeV_MG5_aMCatNLO_pythia8-pLHE/ggTauTau_TuneCUEP_5p02TeV_MG5_aMCatNLO_pythia8-AODSIM-v2/210316_105241/0000/HIN-HINPbPbWinter16DR-00315_10.root' # 2015 Signal
#  '/store/himc/HINPbPbWinter16DR/ggBBbar_4f_TuneCUETP8M1_5p02TeV_MG5_aMCatNLO_pythia8/AODSIM/NoPU_75X_mcRun2_HeavyIon_v14_ext1-v2/280000/00359902-EF3E-EB11-BA40-FA163ED5170D.root' #2015 BBbar
#  '/store/himc/HINPbPbWinter16DR/ggCCbar_TuneCUETP8M1_5p02TeV_MG5_aMCatNLO_pythia8/AODSIM/NoPU_75X_mcRun2_HeavyIon_v14_ext1-v3/40000/00FA85F6-4673-EB11-8D21-842B2B6890DE.root' #2015 CCbar
#  '/store/himc/HINPbPbWinter16DR/ggCCbar_TuneCUETP8M1_5p02TeV_MG5_aMCatNLO_pythia8/AODSIM/NoPU_75X_mcRun2_HeavyIon_v14_ext1-v3/40000/F6AD22F8-7173-EB11-89BE-BC97E17B3080.root' #2015 CCbar
#  '/store/user/gkrintir/ggTauTau_TuneCUEP_5p02TeV_MG5_aMCatNLO_pythia8-pLHE/ggTauTau_TuneCUEP_5p02TeV_MG5_aMCatNLO_pythia8-AODSIM-v2/210319_031327/0000/HIN-HINPbPbWinter16DR-00315_10.root'
#  '/store/himc/HINPbPbAutumn18DR/ggTauTau_TuneCP5_5p02TeV_SuperChic_pythia8/AODSIM/NoPUlowPtPhotonReg_LbyL_103X_upgrade2018_realistic_HI_LowPtPhotonReg_v2-v2/10000/05D69AFE-D2D0-6647-803A-E2F4FA2D87F1.root' #SuperChic Signal 2018
#  '/store/hidata/HIRun2015/HIForward/AOD/02May2016-v1/50000/5C634C6F-1019-E611-BCBF-D4AE528FF351.root'
#  '/store/himc/HINPbPbWinter16DR/ggBBbar_4f_TuneCUETP8M1_5p02TeV_MG5_aMCatNLO_pythia8/AODSIM/NoPU_75X_mcRun2_HeavyIon_v14_ext1-v2/280000/00359902-EF3E-EB11-BA40-FA163ED5170D.root'
#  '/store/group/phys_heavyions/rchudasa/ggtautu/GammaGammaTauTau_PbPb_2015_MC/ggTauTau_TuneCP5_5p02TeV_amcatnlo_pythia8_pLHE/ggTauTau_TuneCP5_5p02TeV_amcatnlo_pythia8_RECO/201018_084927/0000/ggtautau_reco_1.root'
#  '/store/himc/HINPbPbAutumn18DR/ggTauTau_TuneCP5_5p02TeV_SuperChic_pythia8/AODSIM/NoPUlowPtPhotonReg_LbyL_103X_upgrade2018_realistic_HI_LowPtPhotonReg_v2-v2/10000/05D69AFE-D2D0-6647-803A-E2F4FA2D87F1.root'
#  '/store/group/phys_heavyions/rchudasa/ggtautu/GammaGammaTauTau_PbPb_2015_MC/ggTauTau_TuneCP5_5p02TeV_amcatnlo_pythia8_pLHE/ggTauTau_TuneCP5_5p02TeV_amcatnlo_pythia8_RECO/201018_084927/0000/ggtautau_reco_10.root'
#  '/store/himc/HINPbPbWinter16DR/ggBBbar_4f_TuneCUETP8M1_5p02TeV_MG5_aMCatNLO_pythia8/AODSIM/NoPU_75X_mcRun2_HeavyIon_v14-v1/20000/2A31F974-711E-EB11-92FB-20040FEABE68.root',
#  'file:/scratch/ytakahas/2A31F974-711E-EB11-92FB-20040FEABE68.root'
#  '/store/hidata/HIRun2015/HIForward/AOD/02May2016-v1/00000/0038049F-0D25-E611-B57F-F01FAFD691F4.root'
]

options.parseArguments()

process.options  = cms.untracked.PSet( 
                     wantSummary = cms.untracked.bool(True),
                     SkipEvent = cms.untracked.vstring('ProductNotFound'),
                     allowUnscheduled = cms.untracked.bool(True),
                     )

#process.options.numberOfThreads=cms.untracked.uint32(2)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

# run, lumi, event ID
#                              eventsToProcess = cms.untracked.VEventRange('1:94:182460'),

if config["RUNONMC"]:
  process.source = cms.Source("PoolSource",
                              fileNames = cms.untracked.vstring(options.inputFiles),
#                              eventsToProcess = cms.untracked.VEventRange('1:417:334860'),
                              duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),                              
                              ) 
else:                    
  process.source = cms.Source("PoolSource",
                              fileNames = cms.untracked.vstring(options.inputFiles),
#                              skipEvents=cms.untracked.uint32(23000)
                              ) 

print " process source filenames %s" %(process.source) 
######## Sequence settings ##########

hltFiltersProcessName = 'RECO'
#import pdb; pdb.set_trace()
if config["RUNONMC"] or config["JSONFILE"].find('reMiniAOD') != -1:
  hltFiltersProcessName = 'PAT'

# ####### Logger ##########
process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('Ntuple')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit = cms.untracked.int32(1)
)

process.MessageLogger.cerr.FwkReport.reportEvery = 10000

####### Define conditions ##########
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag import GlobalTag

GT = ''


jetcorr_levels=[]
jetcorr_levels_groomed=[]
if config["RUNONMC"]:
  jetcorr_levels = cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute'])
  jetcorr_levels_groomed = cms.vstring(['L2Relative', 'L3Absolute']) # NO L1 corretion for groomed jets
else:
  jetcorr_levels = cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual'])
  jetcorr_levels_groomed = cms.vstring(['L2Relative', 'L3Absolute', 'L2L3Residual'])

   
######### read JSON file for data ##########					                                                             
if not(config["RUNONMC"]) and config["USEJSON"]:

  import FWCore.PythonUtilities.LumiList as LumiList
  import FWCore.ParameterSet.Types as CfgTypes
  process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
  myLumis = LumiList.LumiList(filename = config["JSONFILE"]).getCMSSWString().split(',')
  process.source.lumisToProcess.extend(myLumis) 

  

####### Redo Jet clustering sequence ##########
betapar = cms.double(0.0)
fatjet_ptmin = 100.0

from RecoJets.Configuration.RecoPFJets_cff import *
from RecoJets.JetProducers.AnomalousCellParameters_cfi import *
from RecoJets.JetProducers.PFJetParameters_cfi import *

from PhysicsTools.PatAlgos.tools.helpers import *
pattask = getPatAlgosToolsTask(process)
                                                                                                          
process.chs = cms.EDFilter("CandPtrSelector",
  src = cms.InputTag('packedPFCandidates'),
  cut = cms.string('fromPV')
)

process.ak4PFJetsCHS = ak4PFJetsCHS.clone( src = 'chs' )
process.ak4PFJetsCHS.doAreaFastjet = True


from PhysicsTools.SelectorUtils.tools.vid_id_tools import *

#dataFormat=DataFormat.MiniAOD
#switchOnVIDElectronIdProducer(process,dataFormat,task=pattask)

#process.egmGsfElectronIDSequence = cms.Sequence(process.egmGsfElectronIDs)
     
#my_id_modules = [
#                 'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Fall17_94X_V1_cff',
#                 'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronHLTPreselecition_Summer16_V1_cff',
#                 'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff',
#                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V1_cff',
#                 ]
           
#add them to the VID producer
#for idmod in my_id_modules:
#    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection,task=pattask)

####### Ntuplizer initialization ##########
jetsAK4 = "slimmedJets"


METS = "slimmedMETs"
METS_EGclean = "slimmedMETsEGClean"
METS_MEGclean = "slimmedMETsMuEGClean"
METS_uncorr = "slimmedMETsUncorrected"

if config["USENOHF"]: METS = "slimmedMETsNoHF"  

##___________________ MET significance and covariance matrix ______________________##


##___________________ Jets ______________________##

  

######## JEC ########
jecLevelsAK8chs = []
jecLevelsAK8Groomedchs = []
jecLevelsAK4chs = []
jecLevelsAK4 = []
jecLevelsAK8Puppi = []
jecLevelsForMET = []

print "1. options.RunPeriod ", options.RunPeriod
if options.RunPeriod=="" : options.RunPeriod=options.inputFiles[0]

if  config["RUNONMC"] :
  JECprefix = ""
  if ("Fall17" in options.RunPeriod):
    JECprefix = "Fall17_17Nov2017_V32"
    GT ='94X_mcRun2_asymptotic_v3'
  elif ("Summer16" in options.RunPeriod):
    JECprefix = "Summer16_07Aug2017_V11"
    GT = '94X_mc2017_realistic_v17'
  elif (("Autumn18" in options.RunPeriod) or ("Fall18" in options.RunPeriod)):
    JECprefix = "Autumn18_V8"
    GT = '102X_upgrade2018_realistic_v18'
  else:
    JECprefix = "Autumn18_V8"
    GT = '102X_upgrade2018_realistic_v18'

  #jecAK8chsUncFile = "JEC/%s_MC_Uncertainty_AK8PFchs.txt"%(JECprefix)
  jecAK4chsUncFile = "JEC/%s_MC_Uncertainty_AK4PFchs.txt"%(JECprefix)
 



else : #Data
   JECprefix = ""
   JEC_runDependent_suffix= ""

   if ("2017" in options.RunPeriod):
     if ("Run2017B" in  options.RunPeriod): JEC_runDependent_suffix= "B"
     elif ("Run2017C" in  options.RunPeriod): JEC_runDependent_suffix= "C"
     elif ("Run2017D" in  options.RunPeriod): JEC_runDependent_suffix= "D"
     elif ("Run2017E" in  options.RunPeriod): JEC_runDependent_suffix= "E"
     elif ("Run2017F" in  options.RunPeriod): JEC_runDependent_suffix= "F"
     
     JECprefix = "Fall17_17Nov2017"+JEC_runDependent_suffix+"_V32"
     GT = '94X_dataRun2_v11'
     
     
   elif ("2016" in options.RunPeriod):
     if ("Run2016D" in  options.RunPeriod or "Run2016B" in  options.RunPeriod  or "Run2016C" in  options.RunPeriod  ): JEC_runDependent_suffix= "ABC"
     elif ("Run2016E" in  options.RunPeriod): JEC_runDependent_suffix= "EF"
     elif ("Run2016G" in  options.RunPeriod): JEC_runDependent_suffix= "GH"
     elif ("Run2016F" in  options.RunPeriod and  not options.runUpToEarlyF): JEC_runDependent_suffix= "GH"
     elif ("Run2016F" in  options.RunPeriod and   options.runUpToEarlyF): JEC_runDependent_suffix= "EF"


     JECprefix = "Summer16_07Aug2017"+JEC_runDependent_suffix+"_V11"
     GT ='94X_dataRun2_v10'

   elif ("2018" in options.RunPeriod):
     if ("Run2018A" in  options.RunPeriod ): 
       JEC_runDependent_suffix= "A"
       GT="102X_dataRun2_Sep2018ABC_v2" 
     elif ("Run2018B" in  options.RunPeriod): 
       JEC_runDependent_suffix= "B"
       GT="102X_dataRun2_Sep2018ABC_v2"
     elif ("Run2018C" in  options.RunPeriod): 
       JEC_runDependent_suffix= "C"
       GT="102X_dataRun2_Sep2018ABC_v2"
     elif ("Run2018D" in  options.RunPeriod): 
       JEC_runDependent_suffix= "D"
       GT = '102X_dataRun2_Prompt_v16' 

     JECprefix = "Autumn18_Run"+JEC_runDependent_suffix+"_V8"
    
   #jecAK8chsUncFile = "JEC/%s_DATA_Uncertainty_AK8PFchs.txt"%(JECprefix)
   jecAK4chsUncFile = "JEC/%s_DATA_Uncertainty_AK4PFchs.txt"%(JECprefix)
 
#   GT = '106X_dataRun2_v27' 
   print "jec JEC_runDependent_suffix %s ,  prefix %s " %(JEC_runDependent_suffix,JECprefix)

print "jec prefix ", JECprefix

print "doing corrections  to met on the fly %s" ,config["CORRMETONTHEFLY"]

print "*************************************** GLOBAL TAG *************************************************" 
print GT
print "****************************************************************************************************" 
process.GlobalTag = GlobalTag(process.GlobalTag, GT)




if config["CORRMETONTHEFLY"]:  
   if config["RUNONMC"]:
     jecLevelsForMET = [				       
     	 'JEC/%s_MC_L1FastJet_AK4PFchs.txt'%(JECprefix),
     	 'JEC/%s_MC_L2Relative_AK4PFchs.txt'%(JECprefix),
     	 'JEC/%s_MC_L3Absolute_AK4PFchs.txt'%(JECprefix)
       ]
   else:       					       
     jecLevelsForMET = [
     	 'JEC/%s_DATA_L1FastJet_AK4PFchs.txt'%(JECprefix),
     	 'JEC/%s_DATA_L2Relative_AK4PFchs.txt'%(JECprefix),
     	 'JEC/%s_DATA_L3Absolute_AK4PFchs.txt'%(JECprefix),
         'JEC/%s_DATA_L2L3Residual_AK4PFchs.txt'%(JECprefix)
       ]	
      			    

                                                                       
################## Ntuplizer ###################
process.ntuplizer = cms.EDAnalyzer("Ntuplizer",
    runOnMC	      = cms.bool(config["RUNONMC"]),
#    useHammer	      = cms.bool(config["USEHAMMER"]),
    doGenParticles    = cms.bool(config["DOGENPARTICLES"]),
    doGenEvent	      = cms.bool(config["DOGENEVENT"]),
    doPileUp	      = cms.bool(config["DOPILEUP"]),
    doBsTauTau	      = cms.bool(config["DOBSTAUTAU"]),
    dotwoPiBsTauTau	      = cms.bool(config["DOtwoPiBSTAUTAU"]),
    isTruth           = cms.bool(config["ISTRUTH"]),
    doVertices	      = cms.bool(config["DOVERTICES"]),
    doMissingEt       = cms.bool(config["DOMISSINGET"]),
    doGenHist         = cms.bool(config["DOGENHIST"]),
    verbose           = cms.bool(config["VERBOSE"]),
    dzcut             = cms.double(config['DZCUT']),
    fsigcut           = cms.double(config['FSIGCUT']),
    vprobcut          = cms.double(config['VPROBCUT']),
    tau_charge        = cms.uint32(config['TAU_CHARGE']),
    vertices = cms.InputTag("offlinePrimaryVertices"),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    taus = cms.InputTag("slimmedTaus"),
    muons = cms.InputTag("muons"),
    electrons = cms.InputTag("slimmedElectrons"),
    ebRecHits = cms.InputTag("reducedEgamma","reducedEBRecHits"),
    centralitySrc = cms.InputTag("hiCentrality"),

#    eleHEEPId51Map = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV51"),
#    eleHEEPIdMap = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV60"),
#    eleVetoIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-veto"),
#    eleLooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
#    eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
#    eleTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight"),

#    eleVetoIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V1-veto"),
#    eleLooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V1-loose"),
#    eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V1-medium"),
#    eleTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V1-tight"),

#    eleHLTIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronHLTPreselection-Summer16-V1"), 
#    eleHEEPIdMap = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV70"),
                                   
#    eleMVAMediumIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-noIso-V1-wp90"),
#    eleMVATightIdMap  = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-noIso-V1-wp80"),
#    mvaValuesMap     = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17NoIsoV1Values"),
#    mvaCategoriesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Categories"),
    dupCluster          = cms.InputTag("particleFlowEGammaGSFixed:dupECALClusters"),
    hitsNotReplaced     = cms.InputTag("ecalMultiAndGSGlobalRecHitEB:hitsNotReplaced"),
    mets = cms.InputTag(METS),
    mets_EGclean = cms.InputTag(METS_EGclean),
    mets_MEGclean = cms.InputTag(METS_MEGclean),
    mets_uncorr = cms.InputTag(METS_uncorr),
    mets_puppi = cms.InputTag("slimmedMETsPuppi"),
    mets_mva = cms.InputTag("MVAMET","MVAMET"),
    corrMetPx = cms.string("+0.1166 + 0.0200*Nvtx"),
    corrMetPy = cms.string("+0.2764 - 0.1280*Nvtx"),
    jecAK4forMetCorr = cms.vstring( jecLevelsForMET ),
    jetsForMetCorr = cms.InputTag(jetsAK4),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    genparticles = cms.InputTag("prunedGenParticles"),
    gentaus = cms.InputTag("tauGenJets"),
    PUInfo = cms.InputTag("addPileupInfo"),
    genEventInfo = cms.InputTag("generator"),
    externallheProducer = cms.InputTag("externalLHEProducer"),
    HLT = cms.InputTag("TriggerResults","","HLT"),
    triggerobjects = cms.InputTag("slimmedPatTrigger"),
    triggerprescales = cms.InputTag("patTrigger"),
    noiseFilter = cms.InputTag('TriggerResults','', hltFiltersProcessName),
    jecpath = cms.string(''),
   
    
    ## Noise Filters ###################################
    # defined here: https://github.com/cms-sw/cmssw/blob/CMSSW_7_4_X/PhysicsTools/PatAlgos/python/slimming/metFilterPaths_cff.py
    noiseFilterSelection_HBHENoiseFilter = cms.string('Flag_HBHENoiseFilter'),   # both data and MC for 2018,
    noiseFilterSelection_HBHENoiseFilterLoose = cms.InputTag("HBHENoiseFilterResultProducer", "HBHENoiseFilterResultRun2Loose"),
    noiseFilterSelection_HBHENoiseFilterTight = cms.InputTag("HBHENoiseFilterResultProducer", "HBHENoiseFilterResultRun2Tight"),
    noiseFilterSelection_HBHENoiseIsoFilter = cms.InputTag("HBHENoiseFilterResultProducer", "HBHEIsoNoiseFilterResult"),    # both data and MC for 2018,  
    noiseFilterSelection_ecalBadCalibReducedMINIAODFilter = cms.InputTag("ecalBadCalibReducedMINIAODFilter"),  # both data and MC for 2018,
    noiseFilterSelection_CSCTightHaloFilter = cms.string('Flag_CSCTightHaloFilter'),
    noiseFilterSelection_CSCTightHalo2015Filter = cms.string('Flag_CSCTightHalo2015Filter'),
    noiseFilterSelection_hcalLaserEventFilter = cms.string('Flag_hcalLaserEventFilter'),
    noiseFilterSelection_EcalDeadCellTriggerPrimitiveFilter = cms.string('Flag_EcalDeadCellTriggerPrimitiveFilter'),  # both data and MC for 2018,
    noiseFilterSelection_goodVertices = cms.string('Flag_goodVertices'),  # both data and MC for 2018,
    noiseFilterSelection_trackingFailureFilter = cms.string('Flag_trackingFailureFilter'),
    noiseFilterSelection_eeBadScFilter = cms.string('Flag_eeBadScFilter'),
    noiseFilterSelection_ecalLaserCorrFilter = cms.string('Flag_ecalLaserCorrFilter'),
    noiseFilterSelection_trkPOGFilters = cms.string('Flag_trkPOGFilters'),
    
    #New for ICHEP 2016
    noiseFilterSelection_CSCTightHaloTrkMuUnvetoFilter = cms.string('Flag_CSCTightHaloTrkMuUnvetoFilter'),
    noiseFilterSelection_globalTightHalo2016Filter = cms.string('Flag_globalTightHalo2016Filter'),
    noiseFilterSelection_globalSuperTightHalo2016Filter = cms.string('Flag_globalSuperTightHalo2016Filter'), # both data and MC for 2018,  
    noiseFilterSelection_HcalStripHaloFilter = cms.string('Flag_HcalStripHaloFilter'),
    noiseFilterSelection_chargedHadronTrackResolutionFilter = cms.string('Flag_chargedHadronTrackResolutionFilter'),
    noiseFilterSelection_muonBadTrackFilter = cms.string('Flag_muonBadTrackFilter'),
    
    #New for Moriond
    noiseFilterSelection_badMuonsFilter = cms.string('Flag_BadPFMuonFilter'),    #('Flag_badMuons'),  # both data and MC for 2018, 
    noiseFilterSelection_duplicateMuonsFilter = cms.string('Flag_duplicateMuons'),
    noiseFilterSelection_nobadMuonsFilter = cms.string('Flag_nobadMuons'),

    # and the sub-filters
    noiseFilterSelection_trkPOG_manystripclus53X = cms.string('Flag_trkPOG_manystripclus53X'),
    noiseFilterSelection_trkPOG_toomanystripclus53X = cms.string('Flag_trkPOG_toomanystripclus53X'),
    noiseFilterSelection_trkPOG_logErrorTooManyClusters = cms.string('Flag_trkPOG_logErrorTooManyClusters'),
    # summary
    noiseFilterSelection_metFilters = cms.string('Flag_METFilters'),

    packedpfcandidates = cms.InputTag('particleFlow'),
    SecondaryVertices = cms.InputTag('slimmedSecondaryVertices'),
#    losttrack = cms.InputTag('lostTracks')
)

process.load('RecoMET.METFilters.ecalBadCalibFilter_cfi')

baddetEcallist = cms.vuint32(
    [872439604,872422825,872420274,872423218,
     872423215,872416066,872435036,872439336,
     872420273,872436907,872420147,872439731,
     872436657,872420397,872439732,872439339,
     872439603,872422436,872439861,872437051,
     872437052,872420649,872422436,872421950,
     872437185,872422564,872421566,872421695,
     872421955,872421567,872437184,872421951,
     872421694,872437056,872437057,872437313])


process.ecalBadCalibReducedMINIAODFilter = cms.EDFilter(
    "EcalBadCalibFilter",
    EcalRecHitSource = cms.InputTag("reducedEgamma:reducedEERecHits"),
    ecalMinEt        = cms.double(50.),
    baddetEcal    = baddetEcallist, 
    taggingMode = cms.bool(True),
    debug = cms.bool(True)#False
    )




process.TransientTrackBuilderESProducer = cms.ESProducer("TransientTrackBuilderESProducer",
    ComponentName = cms.string('TransientTrackBuilder')
)


####### Final path ##########
process.p = cms.Path()

#process.p += process.ecalBadCalibReducedMINIAODFilter

if config["RUNONMC"]:
  process.load("PhysicsTools.JetMCAlgos.TauGenJets_cfi")
  process.tauGenJets.GenParticles = cms.InputTag('genParticles')
  process.p += process.tauGenJets

  process.ntuplizer.vertices = cms.InputTag("hiSelectedVertex")
  process.ntuplizer.packedpfcandidates = cms.InputTag('particleFlowTmp')
  process.ntuplizer.genparticles = cms.InputTag("genParticles")

process.p += process.ntuplizer
process.p.associate(pattask)

print pattask

#  LocalWords:  tauIdMVAIsoDBoldDMwLT
