import FWCore.ParameterSet.Config as cms

PixelAnalyzer = cms.EDAnalyzer("SiPixelAnalyzer",
                               src = cms.InputTag("siPixelDigis"),
                               outputFile = cms.untracked.string("testhist.root")     
)

