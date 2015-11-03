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
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(2)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/k/katatar/public/PixelReadoutDQM/0CEFE112-8E63-E511-93F7-0025905A60D0.root'),
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
process.GlobalTag = GlobalTag(process.GlobalTag, '75X_dataRun2_HLT_withOfflineCustomisation_v0', '')

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
