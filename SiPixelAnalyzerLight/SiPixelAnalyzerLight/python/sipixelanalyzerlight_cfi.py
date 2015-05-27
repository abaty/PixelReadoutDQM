import FWCore.ParameterSet.Config as cms

PixelAnalyzerLight = cms.EDAnalyzer("SiPixelAnalyzerLight",
    src = cms.InputTag("siPixelDigis"),
	SelectBx = cms.uint32(0), #0 = no selection, 1= select good bx, 2 = select bad bx
	CollisionBx = cms.vuint32(2001,
   2201,
   2401,
   2601,
   10911,
   11111,
   11311,
   11511,
   19821,
   20021,
   20221,
   20421,
   28731,
   28931,
   29131,
   29331)
)

