import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST2")

process.load("FWCore.MessageService.MessageLogger_cfi")
#process.load("Configuration.StandardSequences.Services_cff")
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('EventFilter.SiPixelRawToDigi.SiPixelRawToDigi_cfi')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(50) )

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(
            'file:/mnt/hadoop/cms/store/user/cmcginn/hltJet/HI/20150520_502MB_RAWSkim/rawSkim_L1Cent05_144_1_IQr.root'
                )
                            )

process.TFileService = cms.Service("TFileService",
      fileName = cms.string("histoPixelAnalyzer.root"),
      closeFileFast = cms.untracked.bool(True)
  )


process.GlobalTag = cms.ESSource( "PoolDBESSource",
    globaltag = cms.string( "GR_H_V39" ),
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
    process.GlobalTag = customiseGlobalTag(process.GlobalTag, globaltag = 'MCHI2_74_V3')
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
  trackerGeometryConstants = cms.PSet(
    ROCS_X = cms.int32( 0 ),
    ROCS_Y = cms.int32( 0 ),
    upgradeGeometry = cms.bool( False ),
    BIG_PIX_PER_ROC_Y = cms.int32( 2 ),
    BIG_PIX_PER_ROC_X = cms.int32( 1 ),
    ROWS_PER_ROC = cms.int32( 80 ),
    COLS_PER_ROC = cms.int32( 52 )
  ),
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




process.raw2digi_step = cms.Path(process.siPixelDigis)

process.test = cms.EDAnalyzer('SiPixelAnalyzer',
                              src = cms.InputTag("siPixelDigis"),
                              outputFile = cms.untracked.string("quickTestPixelAna.root")
                              )

process.p = cms.Path(process.test)

process.schedule = cms.Schedule(process.raw2digi_step,process.p)

