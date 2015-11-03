import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST2")

process.load("FWCore.MessageService.MessageLogger_cfi")
#process.load("Configuration.StandardSequences.Services_cff")
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('EventFilter.SiPixelRawToDigi.SiPixelRawToDigi_cfi')

#process.load('Geometry.TrackerGeometryBuilder.customTrackerParametersRun2.py')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(
            #'/store/relval/CMSSW_7_4_10_patch2/SingleMuon/RAW/74X_dataRun2_HLT_frozen_v1_resub_RelVal_singleMuon2015C-v1/00000/1072FBB8-4D53-E511-A496-00261894397E.root'
            # '/store/user/abaty/2010VirginRaw/HICorePhysics/VirginRAW_2010_HICorePhysics_SKIM_Cent_0_5/150615_203732/0000/VirginRAW_2010_HICorePhysics_SKIM_Cent_0_5_1.root'
            'file:/afs/cern.ch/work/k/katatar/public/PixelReadoutDQM/0CEFE112-8E63-E511-93F7-0025905A60D0.root'  
             '/store/user/abaty/2010VirginRaw/HICorePhysics/VirginRAW_2010_HICorePhysics_SKIM_Cent_0_5/150615_203732/0000/VirginRAW_2010_HICorePhysics_SKIM_Cent_0_5_1.root'
            #'/store/relval/CMSSW_7_4_10_patch2/SingleMuon/RAW/74X_dataRun2_HLT_frozen_v1_resub_RelVal_singleMuon2015C-v1/00000/1072FBB8-4D53-E511-A496-00261894397E.root'
            #'/store/relval/CMSSW_7_5_3_patch1/DoubleEG/RAW/75X_dataRun2_HLT_withOfflineCustomisation_v0_newCond_RelVal_doubEG2015B-v1/00000/0CEFE112-8E63-E511-93F7-0025905A60D0.root'

                 )
                            )

process.TFileService = cms.Service("TFileService",
      fileName = cms.string("histoPixelAnalyzer.root"),
      closeFileFast = cms.untracked.bool(True)
  )

#added--------------------------------
import Geometry.HcalEventSetup.hcalTopologyIdeal_cfi
import Geometry.HcalEventSetup.hcalTopologyConstants_cfi as hcalTopologyConstants_cfi
hcalTopologyIdeal = Geometry.HcalEventSetup.hcalTopologyIdeal_cfi.hcalTopologyIdeal.clone()
hcalTopologyIdeal.hcalTopologyConstants = cms.PSet(hcalTopologyConstants_cfi.hcalTopologyConstants)



#added-------------------------------


process.GlobalTag = cms.ESSource( "PoolDBESSource",    
#     globaltag = cms.string( "74X_dataRun2_HLT_frozen_v1" ),
    globaltag = cms.string( "75X_mcRun2_asymptotic_v6" ),
    #globaltag = cms.string( "75X_mcRun2_asymptotic_v6" ),
    #globaltag = cms.string( "75X_dataRun2_HLT_withOfflineCustomisation_v0" ),



    RefreshEachRun = cms.untracked.bool( True ),
    RefreshOpenIOVs = cms.untracked.bool( False ),
    toGet = cms.VPSet(
    ),
    DBParameters = cms.PSet(
      authenticationPath = cms.untracked.string( "." ),
      connectionRetrialTimeOut = cms.untracked.int32( 60 ),
      idleConnectionCleanupPeriod = cms.untracked.int32( 10 ),
      messageLevel = cms.untracked.int32( 0 ),
      enablePoolAutomaticCleanUp = cms.untracked.bool( False ),
      enableConnectionSharing = cms.untracked.bool( True ),
      enableReadOnlySessionOnUpdateConnection = cms.untracked.bool( False ),
      connectionTimeOut = cms.untracked.int32( 0 ),
      connectionRetrialPeriod = cms.untracked.int32( 10 )
    ),
    RefreshAlways = cms.untracked.bool( False ),
    connect = cms.string( "frontier://(proxyurl=http://localhost:3128)(serverurl=http://localhost:800\
0/FrontierOnProd)(serverurl=http://localhost:8000/FrontierOnProd)(retrieve-ziplevel=0)/CMS_CONDITIONS\
" ),
    ReconnectEachRun = cms.untracked.bool( True ),
    BlobStreamerName = cms.untracked.string( "TBufferBlobStreamingService" ),
    DumpStat = cms.untracked.bool( False )
)


if 'GlobalTag' in process.__dict__:
    from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag as customiseGlobalTag
    process.GlobalTag = customiseGlobalTag(process.GlobalTag, globaltag = '75X_mcRun2_asymptotic_v6')
#    process.GlobalTag = customiseGlobalTag(process.GlobalTag, globaltag = '74X_dataRun2_HLT_frozen_v1')
    process.GlobalTag.connect   = 'frontier://FrontierProd/CMS_CONDITIONS'
    process.GlobalTag.pfnPrefix = cms.untracked.string('frontier://FrontierProd/')
    for pset in process.GlobalTag.toGet.value():
        pset.connect = pset.connect.value().replace('frontier://FrontierProd/', 'frontier://FrontierProd/')
    # fix for multi-run processing                                                                    
    process.GlobalTag.RefreshEachRun = cms.untracked.bool( False )
    process.GlobalTag.ReconnectEachRun = cms.untracked.bool( False )



process.TrackerDigiGeometryESModule = cms.ESProducer( "TrackerDigiGeometryESModule",
  appendToDataLabel = cms.string( "" ),
  fromDDD = cms.bool( False ),
  #trackerGeometryConstants = cms.PSet(
  #  ROCS_X = cms.int32( 0 ),
  #  ROCS_Y = cms.int32( 0 ),
  #  upgradeGeometry = cms.bool( False ),
  #  BIG_PIX_PER_ROC_Y = cms.int32( 2 ),
  #  BIG_PIX_PER_ROC_X = cms.int32( 1 ),
  #  ROWS_PER_ROC = cms.int32( 80 ),
  #  COLS_PER_ROC = cms.int32( 52 )
  #),
  #removed for 7_5_3_patch1 running
  applyAlignment = cms.bool( True ),
  alignmentsLabel = cms.string( "" )
)
process.TrackerGeometricDetESModule = cms.ESProducer( "TrackerGeometricDetESModule",
  appendToDataLabel = cms.string( "" ),
  fromDDD = cms.bool( False )
)

process.siPixelQualityESProducer = cms.ESProducer( "SiPixelQualityESProducer",
  ListOfRecordToMerge = cms.VPSet(
    cms.PSet(  record = cms.string( "SiPixelQualityFromDbRcd" ),
      tag = cms.string( "" )
    ),
    cms.PSet(  record = cms.string( "SiPixelDetVOffRcd" ),
      tag = cms.string( "" )
    )
  )
)



#HF Reco----------------

process.load("RecoLocalCalo.HcalRecProducers.HcalHitReconstructor_hf_cfi")

process.hcalReco = cms.Path(process.hcalDigis*process.hfreco)
process.raw2digi_step = cms.Path(process.siPixelDigis)

process.test = cms.EDAnalyzer('SiPixelAnalyzer',
                              src = cms.InputTag("siPixelDigis"),
                              srcHFhits = cms.InputTag("hfreco"),
                              outputFile = cms.untracked.string("quickTestPixelAna.root")
                              )

process.p = cms.Path(process.test)

process.schedule = cms.Schedule(process.hcalReco,process.raw2digi_step)#,process.p)

