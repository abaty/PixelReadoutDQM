// -*- C++ -*-
//
// Package:    SiPixelAnalyzerLight
// Class:      SiPixelAnalyzerLight
//  
/**\class SiPixelAnalyzerLight SiPixelAnalyzerLight.cc Validation/SiPixelAnalyzerLight/src/SiPixelAnalyzerLight.cc

 Description: <one line class summary>
  
 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Ivan Amos Cali
//         Created:  Wed Jun  4 18:45:43 CEST 2008
// $Id$
//
//
#include <memory>
#include <string>
#include <iostream>

// system include files
#include <memory>

// user include files 
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
//#include "DataFormats/Common/interface/Handle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"

//edit here for 740
#include "DataFormats/Common/interface/EDProductfwd.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerNumberingBuilder/interface/GeometricDet.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetType.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
//#include "Geometry/TrackerTopology/interface/RectangularPixelTopology.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// For ROOT
#include "TROOT.h"
#include "TF1.h"
#include "TH2F.h"
#include "TH1F.h"


using namespace edm;
//
// class decleration
//

class SiPixelAnalyzerLight : public edm::EDAnalyzer {
   public:
      explicit SiPixelAnalyzerLight(const edm::ParameterSet&);
      ~SiPixelAnalyzerLight();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
  


  // ----------member data ---------------------------
  std::vector<uint32_t> CollBunches_;
  edm::InputTag src_;
  const TrackerGeometry* trGeo_;
  
  uint32_t SelectBx_;   //0 = no selection, 1= select good bx, 2 = select bad bx
   
  TH1F *DcolumnHits_;
  TH1F *DcolumnHits1_;
  TH1F *DcolumnHits2_;
  TH1F *DcolumnHits3_;  
  TH1F *BarrelSlinkHitsDistrib_;
  TH1F *BarrelSlinkHitsDistrib1_;
  TH1F *BarrelSlinkHitsDistrib2_;
  TH1F *BarrelSlinkHitsDistrib3_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
SiPixelAnalyzerLight::SiPixelAnalyzerLight(const edm::ParameterSet& iConfig):
CollBunches_(iConfig.getParameter<std::vector<uint32_t> >("CollisionBx"))

{
   //now do what ever initialization is needed 
   
   src_ =  iConfig.getParameter<edm::InputTag>( "src" );
   SelectBx_ = iConfig.getParameter<uint32_t>( "SelectBx" );
   
   /*
   CollBunches_.clear();
   CollBunches_.push_back(2001);
   CollBunches_.push_back(2201);
   CollBunches_.push_back(2401);
   CollBunches_.push_back(2601);
   CollBunches_.push_back(10911);
   CollBunches_.push_back(11111);
   CollBunches_.push_back(11311);
   CollBunches_.push_back(11511);
   CollBunches_.push_back(19821);
   CollBunches_.push_back(20021);
   CollBunches_.push_back(20221);
   CollBunches_.push_back(20421);
   CollBunches_.push_back(28731);
   CollBunches_.push_back(28931);
   CollBunches_.push_back(29131);
   CollBunches_.push_back(29331);
  */

    edm::Service<TFileService> fs;
 
 
     //General histos---------------------------------------------------------
	 
     BarrelSlinkHitsDistrib_ = fs->make<TH1F>("BarrelSlinkHitsDistrib", "Slink Hits -Barrel-", 1001, -0.5, 1000.5);
     BarrelSlinkHitsDistrib_->SetYTitle("Entries");
     BarrelSlinkHitsDistrib_->SetXTitle("Hits per RO Link");
     //BarrelSlinkHitsDistrib_->SetDirectory(oFile_->GetDirectory(0));
	 BarrelSlinkHitsDistrib1_ = fs->make<TH1F>("BarrelSlinkHitsDistrib1", "Slink Hits -Barrel L1-", 1001, -0.5, 1000.5);
     BarrelSlinkHitsDistrib1_->SetYTitle("Entries");
     BarrelSlinkHitsDistrib1_->SetXTitle("Hits per RO Link");
     //BarrelSlinkHitsDistrib1_->SetDirectory(oFile_->GetDirectory(0));
	 BarrelSlinkHitsDistrib2_ = fs->make<TH1F>("BarrelSlinkHitsDistrib2", "Slink Hits -Barrel L2-", 1001, -0.5, 1000.5);
     BarrelSlinkHitsDistrib2_->SetYTitle("Entries");
     BarrelSlinkHitsDistrib2_->SetXTitle("Hits per RO Link");
     //BarrelSlinkHitsDistrib2_->SetDirectory(oFile_->GetDirectory(0));
	 BarrelSlinkHitsDistrib3_ = fs->make<TH1F>("BarrelSlinkHitsDistrib3", "Slink Hits -Barrel L3-", 1001, -0.5, 1000.5);
	 BarrelSlinkHitsDistrib3_->SetYTitle("Entries");
     BarrelSlinkHitsDistrib3_->SetXTitle("Hits per RO Link");
     //BarrelSlinkHitsDistrib3_->SetDirectory(oFile_->GetDirectory(0));
	 
	 
     DcolumnHits_ = fs->make<TH1F>("DcolumnHits", "Hits per Dcolumn", 161, -0.5, 160.5);
     DcolumnHits_ ->SetYTitle("Entries");
     DcolumnHits_->SetXTitle("Hits per Dcolumn");
     //DcolumnHits_->SetDirectory(oFile_->GetDirectory(0));
 
     DcolumnHits1_ = fs->make<TH1F>("DcolumnHits1", "Hits per Dcolumn L1", 161, -0.5, 160.5);
     DcolumnHits1_ ->SetYTitle("Entries");
     DcolumnHits1_->SetXTitle("Hits per Dcolumn");
     //DcolumnHits1_->SetDirectory(oFile_->GetDirectory(0));
	 
	 DcolumnHits2_ = fs->make<TH1F>("DcolumnHits2", "Hits per Dcolumn L2", 161, -0.5, 160.5);
     DcolumnHits2_ ->SetYTitle("Entries");
     DcolumnHits2_->SetXTitle("Hits per Dcolumn");
     //DcolumnHits2_->SetDirectory(oFile_->GetDirectory(0));
	 
	 DcolumnHits3_ = fs->make<TH1F>("DcolumnHits3", "Hits per Dcolumn L3", 161, -0.5, 160.5);
     DcolumnHits3_ ->SetYTitle("Entries");
     DcolumnHits3_->SetXTitle("Hits per Dcolumn");
     //DcolumnHits3_->SetDirectory(oFile_->GetDirectory(0));
   
   
    
}


SiPixelAnalyzerLight::~SiPixelAnalyzerLight()
{
 
   
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
SiPixelAnalyzerLight::analyze(const edm::Event& e, const edm::EventSetup& iSetup)
{
  //int run       = e.id().run();
  //int event     = e.id().event();

  //int lumiBlock = e.luminosityBlock();
  //  uint32_t bx        = e.bunchCrossing();
  //int orbit     = e.orbitNumber();

  //COMMENT OUT TILL SELECT HERE  
  //  bool GoodBx = false; //SelectBx_ : 0 = no selection, 1= select good bx, 2 = select bad bx
  /*
  for(size_t itBXl =0; itBXl <  CollBunches_.size(); ++itBXl){
	if(CollBunches_[itBXl] == bx){
		GoodBx = true;
		itBXl = CollBunches_.size();
	}
  }
  */
  
  //  if((SelectBx_==1&&!GoodBx)||(SelectBx_==2&&GoodBx)) return;
    
 
 // using namespace sipixelobjects;
   edm::Handle<edm::DetSetVector<PixelDigi> > pixelDigis;
   e.getByLabel(src_, pixelDigis);

   
   edm::ESHandle<TrackerGeometry> tracker;
   iSetup.get<TrackerDigiGeometryRecord>().get( tracker );    
   trGeo_ = tracker.product();
  

   //looping over FED and Links
    edm::DetSetVector<PixelDigi>::const_iterator DSViter;
    
    
     
    //looping over the DIgis
   //-------------------------------------------------------------------------   
    //-------------------------------------------------------------------------
    for( DSViter = pixelDigis->begin() ; DSViter != pixelDigis->end(); DSViter++) {
          uint32_t id = DSViter->id;
          DetId  detId(id);
                    
          edm::DetSet<PixelDigi>::const_iterator iter;          
          
 
          //detector Geometry-- the coordinate are of the module center
          const PixelGeomDetUnit* PixelModuleGeom = dynamic_cast<const PixelGeomDetUnit*> (trGeo_->idToDet(id));   //detector geometry -> it returns the center of the module
         // double detZ = PixelModuleGeom->surface().position().z();        //module z      
         // double detR = PixelModuleGeom->surface().position().perp();        //module R                                
         // double detEta = PixelModuleGeom->surface().position().eta();    //module eta                           
         // double detPhi = PixelModuleGeom->surface().position().phi();    //module phi 
          uint32_t ncolumns = PixelModuleGeom->specificTopology().ncolumns(); //n of columns
          uint32_t nrows = PixelModuleGeom->specificTopology().nrows();       //n of rows
          
         // uint32_t nROCs = ncolumns * nrows / (80*52);                   
          uint32_t nHModule = 1;
          if (nrows >80) nHModule = 2;
          


	  //Selct the Barrel
          //---------------------------------------------------------------------------------------------------------
         if(detId.subdetId() ==PixelSubdetector::PixelBarrel ) {                //selcting barrel modules
             PXBDetId  bdetid(id);
             
             uint32_t layer  = bdetid.layer();   // Layer:1,2,3.
           //  uint32_t ladder = bdetid.ladder();  // Ladeer: 1-20, 32, 44. 
           //  uint32_t ring = bdetid.module();  // Z-index: 1-8. (ring)
       
           
             uint32_t BarrelModuleHits=0, HModHits[2];
             HModHits[0] = 0;
             HModHits[1] = 0;     
             
             
             std::vector<uint32_t>  rocColumnsHits, rocDColumnsHits;
			 rocColumnsHits.clear();
			 rocColumnsHits.insert(rocColumnsHits.begin(), nHModule * ncolumns, 0);
			 rocDColumnsHits.clear();
             rocDColumnsHits.insert(rocDColumnsHits.begin(),(nHModule *  ncolumns)/2, 0);
                                   
            
            //std::cout << "Looping over Digis" << std::endl;             
             // loop over module Digis
             //------------------------------------------------------------------
             for ( iter = DSViter->data.begin() ; iter != DSViter->data.end(); iter++ ){  
                ++BarrelModuleHits;
                uint32_t HModule = 0;
                if ((*iter).row() >= 80) HModule =1;
                ++rocColumnsHits[HModule * ncolumns + (*iter).column()];
		        ++rocDColumnsHits[(HModule * ncolumns +(*iter).column())/2];
                ++HModHits[HModule];  
                
             }

			// std::cout << "Writing Dcol" << std::endl;              
             for(size_t dcolumns = 0; dcolumns < rocDColumnsHits.size();  ++dcolumns){
			   if(layer ==1) DcolumnHits1_->Fill(rocDColumnsHits[dcolumns]);
			   if(layer ==2) DcolumnHits2_->Fill(rocDColumnsHits[dcolumns]);
			   if(layer ==3) DcolumnHits3_->Fill(rocDColumnsHits[dcolumns]);
			   DcolumnHits_->Fill(rocDColumnsHits[dcolumns]);
			 }

             //std::cout << "Writing link" << std::endl;             
             if(layer ==1){
			   BarrelSlinkHitsDistrib1_->Fill(HModHits[0]);
			   BarrelSlinkHitsDistrib_->Fill(HModHits[0]);
               if(nHModule > 1){
					BarrelSlinkHitsDistrib1_->Fill(HModHits[1]);           
					BarrelSlinkHitsDistrib_->Fill(HModHits[1]);
				}
			}else if(layer ==2){
	           BarrelSlinkHitsDistrib2_->Fill(HModHits[0]);
			   BarrelSlinkHitsDistrib_->Fill(HModHits[0]);
               if(nHModule > 1){
					BarrelSlinkHitsDistrib2_->Fill(HModHits[1]);           
					BarrelSlinkHitsDistrib_->Fill(HModHits[1]);
				}
	         
	         }else if (layer ==3){
			   BarrelSlinkHitsDistrib3_->Fill(BarrelModuleHits);
	           BarrelSlinkHitsDistrib_->Fill(BarrelModuleHits);
             }

            


             //std::cout << "End" << std::endl;
             
             
         }

 

}
}

// ------------ method called once each job just before starting event loop  ------------
void 
SiPixelAnalyzerLight::beginJob(const edm::EventSetup&)
{
      
    
    
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SiPixelAnalyzerLight::endJob() {
 //  oFile_->Write();
 //  oFile_->Close();
}


//define this as a plug-in
DEFINE_FWK_MODULE(SiPixelAnalyzerLight);
