# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: reco -s RAW2DIGI --filein file:/afs/cern.ch/work/k/katatar/public/PixelReadoutDQM/0CEFE112-8E63-E511-93F7-0025905A60D0.root --conditions 75X_dataRun2_HLT_withOfflineCustomisation_v0 --no_exec --data -n 4
import FWCore.ParameterSet.Config as cms

process = cms.Process('RAW2DIGI')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
#process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.RawToDigi_Repacked_cff')  # for /HIHighPt/HIRun2011-v1/RAW
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
                                      #'file:/afs/cern.ch/work/k/katatar/public/PixelReadoutDQM/0CEFE112-8E63-E511-93F7-0025905A60D0.root'
#                                       'file:/data/abaty/VirginRaw_CentralitySkims/VirginRAW_2010_HICorePhysics_SKIM_Cent_0_5_1.root',
#                                       'file:/data/abaty/VirginRaw_CentralitySkims/VirginRAW_2010_HICorePhysics_SKIM_Cent_0_5_10.root',
#                                       'file:/data/abaty/VirginRaw_CentralitySkims/VirginRAW_2010_HICorePhysics_SKIM_Cent_25_50_102.root',
#                                       'file:/data/abaty/VirginRaw_CentralitySkims/VirginRAW_2010_HICorePhysics_SKIM_Cent_50_100_104.root'
#                                         '/store/hidata/HIRun2010/HIAllPhysics/RAW/v1/000/152/698A4EE0EAD-D8FB-DF11-A668-003048F1C832.root'
#                                         '/store/hidata/HIRun2010/HIAllPhysics/RAW/v1/000/152/698/00EF189B-CAFB-DF11-B3E7-003048F024FE.root',
#                                         '/store/hidata/HIRun2010/HIAllPhysics/RAW/v1/000/152/698/961F4067-D9FB-DF11-8211-001D09F23A3E.root'
#                                         '/store/hidata/HIRun2011/HIHighPt/RAW/v1/000/181/531/4E33D225-0A0D-E111-8C91-0025901D5DB8.root'

# some files from dataset = /HIHighPt/HIRun2011-v1/RAW
'file:/mnt/hadoop/cms/store/hidata/HIRun2011/HIHighPt/RAW/v1/000/181/611/0856BEB4-380E-E111-9495-00215AEDFCCC.root',
'file:/mnt/hadoop/cms/store/hidata/HIRun2011/HIHighPt/RAW/v1/000/181/611/0A834C89-150E-E111-9D66-BCAEC518FF44.root',
'file:/mnt/hadoop/cms/store/hidata/HIRun2011/HIHighPt/RAW/v1/000/181/611/0C7D9789-280E-E111-9ABE-BCAEC518FF7A.root',
'file:/mnt/hadoop/cms/store/hidata/HIRun2011/HIHighPt/RAW/v1/000/181/611/2CB7DDF9-300E-E111-87E2-003048F1BF66.root',
'file:/mnt/hadoop/cms/store/hidata/HIRun2011/HIHighPt/RAW/v1/000/181/611/30A94F53-240E-E111-8F06-BCAEC5329726.root',
'file:/mnt/hadoop/cms/store/hidata/HIRun2011/HIHighPt/RAW/v1/000/181/611/3253ECE5-1B0E-E111-ACF2-BCAEC532971C.root',
'file:/mnt/hadoop/cms/store/hidata/HIRun2011/HIHighPt/RAW/v1/000/181/611/34B60206-120E-E111-B1AF-BCAEC532970D.root',
'file:/mnt/hadoop/cms/store/hidata/HIRun2011/HIHighPt/RAW/v1/000/181/611/384FD322-140E-E111-A0AB-BCAEC532971F.root',
'file:/mnt/hadoop/cms/store/hidata/HIRun2011/HIHighPt/RAW/v1/000/181/611/38B93B17-200E-E111-96E2-0025901D626C.root',
'file:/mnt/hadoop/cms/store/hidata/HIRun2011/HIHighPt/RAW/v1/000/181/611/483F0EAD-170E-E111-838B-003048F24A04.root',
'file:/mnt/hadoop/cms/store/hidata/HIRun2011/HIHighPt/RAW/v1/000/181/611/48FB70BC-250E-E111-A2E7-BCAEC5329701.root',
'file:/mnt/hadoop/cms/store/hidata/HIRun2011/HIHighPt/RAW/v1/000/181/611/4CD0ADBE-120E-E111-A998-BCAEC53296FB.root',
'file:/mnt/hadoop/cms/store/hidata/HIRun2011/HIHighPt/RAW/v1/000/181/611/5486A159-110E-E111-9E89-BCAEC518FF80.root',
'file:/mnt/hadoop/cms/store/hidata/HIRun2011/HIHighPt/RAW/v1/000/181/611/54B04F5C-180E-E111-A6F4-BCAEC5329717.root',
'file:/mnt/hadoop/cms/store/hidata/HIRun2011/HIHighPt/RAW/v1/000/181/611/5C4454F4-160E-E111-841F-003048F1C424.root',
'file:/mnt/hadoop/cms/store/hidata/HIRun2011/HIHighPt/RAW/v1/000/181/611/5CF98311-330E-E111-99BD-BCAEC5364CED.root',
'file:/mnt/hadoop/cms/store/hidata/HIRun2011/HIHighPt/RAW/v1/000/181/611/60C6E5CD-200E-E111-A79F-BCAEC5364CED.root',
'file:/mnt/hadoop/cms/store/hidata/HIRun2011/HIHighPt/RAW/v1/000/181/611/6488D178-2D0E-E111-804D-BCAEC5364CFB.root',
'file:/mnt/hadoop/cms/store/hidata/HIRun2011/HIHighPt/RAW/v1/000/181/611/6EBB2EA1-3D0E-E111-95B7-003048F118E0.root',
'file:/mnt/hadoop/cms/store/hidata/HIRun2011/HIHighPt/RAW/v1/000/181/611/70449657-2B0E-E111-B19A-002481E0E56C.root'

                                      ),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('reco nevts:1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string(''),
        filterName = cms.untracked.string('')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    fileName = cms.untracked.string('reco_RAW2DIGI.root'),
    outputCommands = process.RECOSIMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run1_data', '')

# Path and EndPath definitions
process.RawToDigi_custom = cms.Sequence(#process.csctfDigis
                         #+process.dttfDigis
                         #+process.gctDigis
                          #+process.gtDigis
                        #+process.gtEvmDigis
                         process.siPixelDigis
                         #+process.siStripDigis
                         +process.hcalDigis
)

#process.raw2digi_step = cms.Path(process.RawToDigi)
#process.L1Reco_step = cms.Path(process.L1Reco)
process.raw2digi_step = cms.Path(process.RawToDigi_custom)

process.load("RecoLocalCalo.HcalRecProducers.HcalHitReconstructor_hf_cfi")
process.hfrecoPath = cms.Path(process.hfreco)

process.test = cms.EDAnalyzer('SiPixelAnalyzer',
                              src = cms.InputTag("siPixelDigis"),
                              srcHFhits = cms.InputTag("hfreco"),
                              outputFile = cms.untracked.string("quickTestPixelAna.root")
                              )

process.p = cms.Path(process.test)

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step, process.hfrecoPath, process.p)
#process.schedule = cms.Schedule(process.raw2digi_step, process.L1Reco_step)#,process.reconstruction_step,process.endjob_step,process.RECOSIMoutput_step)
