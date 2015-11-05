// -*- C++ -*-
//
// Package:    SiPixelAnalyzer
// Class:      SiPixelAnalyzer
//  
/**\class SiPixelAnalyzer SiPixelAnalyzer.cc Validation/SiPixelAnalyzer/src/SiPixelAnalyzer.cc

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


// system include files
#include <memory>
#include <iostream>
#include <ctime>        // std::clock()
#include <algorithm>    // std::find()

// user include files 
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerNumberingBuilder/interface/GeometricDet.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetType.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/RectangularPixelTopology.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"

#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"

#include "CondFormats/SiPixelObjects/interface/SiPixelFedCablingMap.h"
#include "CondFormats/SiPixelObjects/interface/SiPixelFedCablingTree.h"
#include "CondFormats/DataRecord/interface/SiPixelFedCablingMapRcd.h"
#include "CondFormats/SiPixelObjects/interface/PixelFEDCabling.h"
#include "CondFormats/SiPixelObjects/interface/SiPixelFrameConverter.h"
#include "CondFormats/SiPixelObjects/interface/PixelFEDLink.h"
#include "CondFormats/SiPixelObjects/interface/PixelROC.h"
#include "CondFormats/SiPixelObjects/interface/ElectronicIndex.h"
#include "CondFormats/SiPixelObjects/interface/DetectorIndex.h"

//ROOT inclusion
#include "TROOT.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TH1F.h"
#include "TH2F.h"

using namespace std;
using namespace edm;
//
// class decleration
//

class SiPixelAnalyzer : public edm::EDAnalyzer {
   public:
      explicit SiPixelAnalyzer(const edm::ParameterSet&);
      ~SiPixelAnalyzer();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
  


  // ----------member data ---------------------------
  std::string outputFile_;
  edm::InputTag src_;

  edm::EDGetTokenT<HFRecHitCollection> srcHFhits_;

  unsigned int BarrelColumnsOffset_;
  unsigned int gnHModuleBarrel_;
  unsigned int EndcapColumnsOffset_;
  unsigned int gnHModuleEndcap_;
  const TrackerGeometry* trGeo_;
  TFile* oFile_;
  TNtuple* BarrelModuleNt_;
  TNtuple* BarrelDigisNt_;
  TNtuple* BarrelColumnsNt_;
  TNtuple* BarrelDColumnsNt_;
  

  TNtuple* EndcapModuleNt_;
  TNtuple* EndcapDigisNt_;
  TNtuple* EndcapColumnsNt_;
  TNtuple* EndcapDColumnsNt_;  
    
  TNtuple* FEDNt_;
  TNtuple* LinksNt_;
  TNtuple* ROCOccNt_;
  
  TH1F* DcolumnHits_; 
  TH1D* Occupancy_;
  TH2D* OccupancyZ_;
  TH1F * BarrelSlinkHitsDistrib_;

  TH2D* LinkOcc_[6];
  TH1D* LinkByLinkOcc_;  
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
SiPixelAnalyzer::SiPixelAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed 
   outputFile_ = iConfig.getUntrackedParameter<string>("outputFile", "pixeldigihisto.root");
   src_ =  iConfig.getParameter<edm::InputTag>( "src" );
   
   srcHFhits_ = consumes<HFRecHitCollection>(iConfig.getParameter<edm::InputTag>("srcHFhits"));

    
     oFile_ = new TFile((const char*)outputFile_.c_str(), "RECREATE");
     //FEDs studies---------------------------------------------------------------------
     FEDNt_ = new TNtuple("FED", "FED","fedid:links:roc",100000);
     LinksNt_= new TNtuple("Links", "Links","fedid:linkn:nHits",100000);
     ROCOccNt_= new TNtuple("ROCOcc", "ROCOcc","idH:idL:ROCHits:SubDet:gHModule:ROCn:ModuleHits:nROCs",100000);
     
     
     //Barrel Ntuples----------------------------------------------------------------------------
     BarrelModuleNt_ = new TNtuple("BarrelPixelDistrib", "BarrelPixelDistri","idH:idL:layer:ladder:ring:z:R:phi:eta:Mhits:HMod0Hits:HMod1Hits:ncolumns:nrows",100000);
     BarrelDigisNt_ = new TNtuple("BarrelDigis", "BarrelPixelDigisDist","idH:idL:layer:ladder:ring:adc:column:row",100000);
     BarrelColumnsNt_ = new TNtuple("BarrelColumnOcc", "BarrelColumnOcc","idH:idL:layer:ladder:ring:gHmodule:column:colHits:gcolumn",100000);
     BarrelDColumnsNt_ = new TNtuple("BarrelDColumnOcc", "BarrelDColumnOcc","idH:idL:layer:ladder:ring:gHmodule:Dcolumn:DcolHits:gDcolumn",100000);
   
     //histos---------------------------------------------------------
     BarrelSlinkHitsDistrib_ = new TH1F("BarrelSlinkHitsDistrib", "Slink Hits -Barrel-", 801, 0., 800.);
     BarrelSlinkHitsDistrib_->GetYaxis()->SetTitle("Entries");
     BarrelSlinkHitsDistrib_->GetXaxis()->SetTitle("Hits per RO Link");
     BarrelSlinkHitsDistrib_->SetDirectory(oFile_->GetDirectory(0));
 
     for(int i = 0; i<6; i++)
     {    
       if(i<3) LinkOcc_[i] = new TH2D(Form("FEDOccupancy_B_L%d",i+1),Form("FED Occupancy (Barrel Layer%d);HF Energy Sum;Average FED Occupancy",i+1),40,0,300000,100,0,3);
       else    LinkOcc_[i] = new TH2D(Form("FEDOccupancy_EC_D%d",i-2),Form("FED Occupancy (Endcap Disk%d);HF Energy Sum;Average FED Occupancy",i-2),40,0,300000,100,0,3);
       LinkOcc_[i]->SetDirectory(oFile_->GetDirectory(0));
     }
     LinkByLinkOcc_ = new TH1D("LinkByLinkOcc",";36*detid+linkid;hits",1500,0,1500);
     LinkByLinkOcc_->SetDirectory(oFile_->GetDirectory(0));
 
    //Endcap Ntuples----------------------------------------------------------------------------
    EndcapModuleNt_ = new TNtuple("EndcapPixelDistrib", "EndcapPixelDistri","idH:idL:side:disk:blade:panel:module:z:R:phi:eta:HMod0Hits:HMod1Hits:ncolumns:nrows",100000);
    EndcapDigisNt_ = new TNtuple("EndcapDigis", "EndcapPixelDigisDist","idH:idL:side:disk:blade:panel:module:adc:column:row",100000);
    EndcapColumnsNt_ = new TNtuple("EndcapColumnOcc", "EndcapColumnOcc","idH:idL:side:disk:blade:panel:gHmodule:column:colHits:gcolumn",100000);
    EndcapDColumnsNt_ = new TNtuple("EndcapDColumnOcc", "EndcapDColumnOcc","idH:idL:side:disk:blade:panel:gHmodule:Dcolumn:DcolHits:gDcolumn",100000);

     //General histos---------------------------------------------------------
     DcolumnHits_ = new TH1F("DcolumnHits", "Hits per Dcolumn", 161, 0, 160);
     DcolumnHits_ ->GetYaxis()->SetTitle("Entries");
     DcolumnHits_->GetXaxis()->SetTitle("Hits per Dcolumn");
     DcolumnHits_->SetDirectory(oFile_->GetDirectory(0));
 
     Occupancy_ = new TH1D("Occupancy", "Occupancy", 3001, 0., 3.);
     Occupancy_->GetYaxis()->SetTitle("Entries");
     Occupancy_->GetXaxis()->SetTitle("Occupancy [%]");
     Occupancy_->SetDirectory(oFile_->GetDirectory(0));

     OccupancyZ_ = new TH2D("OccupancyZ", "Occupancy", 1020, -51., 51., 3001, 0., 3.);
     OccupancyZ_->GetYaxis()->SetTitle("Occupancy [%]");
     OccupancyZ_->GetXaxis()->SetTitle("Z");
     OccupancyZ_->SetDirectory(oFile_->GetDirectory(0));

}


SiPixelAnalyzer::~SiPixelAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
SiPixelAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    std::clock_t    start_analyze, end_analyze;
    std::clock_t    start_analyze1, start_analyze2;
    start_analyze = std::clock();
    std::cout << "started, CLOCK_REALTIME         : " << (double)CLOCK_REALTIME << std::endl;

   //centrality first
   double HFRecHitSum = 0;
   Handle<HFRecHitCollection> hits;
   iEvent.getByToken(srcHFhits_,hits);
   for( size_t ihit = 0; ihit<hits->size(); ++ ihit){
     const HFRecHit & rechit = (*hits)[ ihit ];
     HFRecHitSum += rechit.energy();
   }
   std::cout << "hits->size() = " <<hits->size() << std::endl;
   std::cout << "HFRecHitSum  = " <<HFRecHitSum << std::endl;

   using namespace sipixelobjects;

   gnHModuleBarrel_=0;
   BarrelColumnsOffset_ =0;
   gnHModuleEndcap_=0;
   EndcapColumnsOffset_ =0;

   start_analyze1 = std::clock();

   edm::Handle<edm::DetSetVector<PixelDigi> > pixelDigis;
   iEvent.getByLabel(src_, pixelDigis);

   edm::ESHandle<SiPixelFedCablingMap> map;
   iSetup.get<SiPixelFedCablingMapRcd>().get( map );
//   const SiPixelFedCablingMap * theCablingMap = map.product();
   std::cout <<"Map number " << map->version() << std::endl;
   
   edm::ESHandle<TrackerGeometry> tracker;
   iSetup.get<TrackerDigiGeometryRecord>().get( tracker );    
   trGeo_ = tracker.product();

   //looping over FED and Links
   //-------------------------------------------------------------------------
   //EDIT HERE FOR CONVERSION ERROR UNIQUE POINTER TO SIPIXELFEDCBLINGTREE
   std::unique_ptr<SiPixelFedCablingTree> CablingTree = map->cablingTree();
   const std::vector<const PixelFEDCabling *>  cabling = CablingTree->fedList();
   typedef std::vector<const PixelFEDCabling *>::const_iterator FEDiter;
   edm::DetSetVector<PixelDigi>::const_iterator DSViter;   

   start_analyze2 = std::clock();

   int numFEDs = cabling.size();
   std::cout << "numFEDs = " << numFEDs << std::endl;
   int numPixelDigis_global = pixelDigis->size();

   std::vector<sipixelobjects::PixelROC> pixelROCs;
   std::vector<unsigned int> rawIDsROC;
   std::vector<unsigned int> fedIDs;
   for (std::vector<const PixelFEDCabling *>::const_iterator ifed=cabling.begin();
    ifed != cabling.end(); ifed++) {
     unsigned int fed = (**ifed).id();
     unsigned int numLink = (**ifed).numberOfLinks();
     for (unsigned int link=1; link <= numLink; link++) {
       const PixelFEDLink * pLink = (**ifed).link(link);
       if (pLink==0) continue;
       unsigned int linkId = pLink->id();
       if (linkId != 0 && linkId!= link)
           std::cout << "PROBLEM WITH LINK NUMBER!!!!" << std::endl;
       unsigned int numberROC = pLink->numberOfROCs();
       for (unsigned int roc=1; roc <= numberROC; roc++) {
         const PixelROC * pROC = pLink->roc(roc);
         if (pROC==0) continue;
         if (pROC->idInLink() != roc)
             std::cout << "PROBLEM WITH ROC NUMBER!!!!" << std::endl;
//         Key key = {fed, link, roc};
//         theMap[key] = (*pROC);
         pixelROCs.push_back(*pROC);
         rawIDsROC.push_back(pROC->rawId());
         fedIDs.push_back(fed);
       }
     }
   }

   int numRawIDsROC = rawIDsROC.size();
   std::cout<< "numRawIDsROC = " << numRawIDsROC <<std::endl;
   uint32_t nhits[numFEDs][6] = { {0} };
   uint32_t totalPix[numFEDs][6] = { {0} };

   for( DSViter = pixelDigis->begin() ; DSViter != pixelDigis->end(); DSViter++) {
       unsigned int id = DSViter->id;
       unsigned int fedID = 999999;

       std::vector<unsigned int>::iterator it = std::find(rawIDsROC.begin(),rawIDsROC.end(), id);
       if(it == rawIDsROC.end()) continue;
       auto pos = it - rawIDsROC.begin();
       fedID = fedIDs.at(pos);

       edm::DetSet<PixelDigi>::const_iterator  begin = DSViter->data.begin();
       edm::DetSet<PixelDigi>::const_iterator  end   = DSViter->data.end();
       edm::DetSet<PixelDigi>::const_iterator iter;

       DetId  detId(id);
       int detLayer = -1;
       if(detId.subdetId() ==PixelSubdetector::PixelBarrel ) {                //selcting barrel modules
           PXBDetId  bdetid(id);
       detLayer  = bdetid.layer()-1;   // Layer:1,2,3. -1
       }
       else{
           PXFDetId  fdetid(id);
           detLayer  = fdetid.disk()+2; //1, 2, 3 +2 for endcap
       }

       const PixelGeomDetUnit* PixelModuleGeom = dynamic_cast<const PixelGeomDetUnit*> (trGeo_->idToDet(id));   //detector geometry -> it returns the center of the module
       uint32_t ncolumns = PixelModuleGeom->specificTopology().ncolumns(); //n of columns
       uint32_t nrows = PixelModuleGeom->specificTopology().nrows();       //n of rows

       int nhits_PixelDigi[6] = {0};

       totalPix[fedID][detLayer] += ncolumns*nrows;
       for ( iter = begin ; iter != end; iter++ ){  //llop over digi
           GlobalPixel global = {(*iter).row(), (*iter).column()};

           for(std::vector<sipixelobjects::PixelROC>::const_iterator it_pixelROCs = pixelROCs.begin(); it_pixelROCs!=pixelROCs.end(); ++it_pixelROCs){
               PixelROC roc = *it_pixelROCs;
               LocalPixel local = roc.toLocal(global);

               if(local.valid())
               {
                   nhits[fedID][detLayer] ++;
                   nhits_PixelDigi[detLayer]++;
                   break;
               }
           }
       }
   }

   for(int j=0; j<numFEDs; ++j)
   {
       std::cout << "Fed id: " << j << std::endl;
       for(int i = 0; i<6; i++){
           if(totalPix[j][i]!=0){
               LinkOcc_[i]->Fill(HFRecHitSum,100.0*nhits[j][i]/((double)totalPix[j][i]));
               std::cout << "Layer: " << i << " Occupancy: " << 100.0*nhits[j][i]/((double)totalPix[j][i]) << std::endl;
           }
       }
   }
   std::cout << "numPixelDigis_global       = " << numPixelDigis_global << std::endl;

   end_analyze = std::clock();
   std::cout.precision(6);      // get back to default precision
   std::cout << "analyzer finished in             : " << (end_analyze - start_analyze) / (double)CLOCKS_PER_SEC << " seconds" << std::endl;
   std::cout << "analyzer1 finished in            : " << (end_analyze - start_analyze1) / (double)CLOCKS_PER_SEC << " seconds" << std::endl;
   std::cout << "analyzer2 finished in            : " << (end_analyze - start_analyze2) / (double)CLOCKS_PER_SEC << " seconds" << std::endl;
   std::cout << "finished, CLOCK_REALTIME         : " << (double)CLOCK_REALTIME << std::endl;
     
     
    //looping over the DIgis
   //-------------------------------------------------------------------------   
    //-------------------------------------------------------------------------

    /* 
    for( DSViter = pixelDigis->begin() ; DSViter != pixelDigis->end(); DSViter++) {
          uint32_t id = DSViter->id;
          float idL = id & 0xFFFF;
          float idH = (id & 0xFFFF0000) >> 16;
          DetId  detId(id);
          //uint32_t SubDet = detId.subdetId();
          
          edm::DetSet<PixelDigi>::const_iterator  begin = DSViter->data.begin();
          edm::DetSet<PixelDigi>::const_iterator  end   = DSViter->data.end();
          edm::DetSet<PixelDigi>::const_iterator iter;          
          
 
          //detector Geometry-- the coordinate are of the module center
          const PixelGeomDetUnit* PixelModuleGeom = dynamic_castconst PixelGeomDetUnit*> (trGeo_->idToDet(id));   //detector geometry -> it returns the center of the module
          //double detZ = PixelModuleGeom->surface().position().z();        //module z      
          //double detR = PixelModuleGeom->surface().position().perp();        //module R                                
          //double detEta = PixelModuleGeom->surface().position().eta();    //module eta                           
          //double detPhi = PixelModuleGeom->surface().position().phi();    //module phi 
          uint32_t ncolumns = PixelModuleGeom->specificTopology().ncolumns(); //n of columns
          uint32_t nrows = PixelModuleGeom->specificTopology().nrows();       //n of rows
          
          //uint32_t nROCs = ncolumns * nrows / (80*52);                   
          uint32_t nHModule = 1;
          if (nrows >80) nHModule = 2;
          


	  //Selct the Barrel
          //---------------------------------------------------------------------------------------------------------
         if(detId.subdetId() ==PixelSubdetector::PixelBarrel ) {                //selcting barrel modules
             PXBDetId  bdetid(id);
             
             uint32_t layer  = bdetid.layer();   // Layer:1,2,3.
             uint32_t ladder = bdetid.ladder();  // Ladeer: 1-20, 32, 44. 
             uint32_t ring = bdetid.module();  // Z-index: 1-8. (ring)
       
           
             unsigned int BarrelModuleHits=0, HModHits[2];
             HModHits[0] = 0;
             HModHits[1] = 0;     
             
             
             uint32_t* rocColumnsHits = new unsigned int [nHModule * ncolumns];
             uint32_t* rocDColumnsHits = new unsigned int [nHModule *  ncolumns/2];
           
             
             
             for(uint32_t j =0; j < nHModule; ++j){
                for (uint32_t i = 0; i< ncolumns; ++i){
		         rocColumnsHits[j* ncolumns +i] = 0;
		         rocDColumnsHits[(j* ncolumns + i)/2] =0;
               }
             }
            
  
             // loop over module Digis
             //------------------------------------------------------------------
             for ( iter = begin ; iter != end; iter++ ){  
                ++BarrelModuleHits;
       
                uint32_t HModule = 0;
                if ((*iter).row() > 80) HModule =1;
                ++rocColumnsHits[HModule * ncolumns + (*iter).column()];
		        ++rocDColumnsHits[(HModule * ncolumns +(*iter).column())/2];
                ++HModHits[HModule];  
   
          //      BarrelDigisNt_->Fill(idH, idL, layer, ladder, ring,(*iter).adc(), (*iter).column(), (*iter).row());
                 
             }

            uint32_t ROCn=0, ROCHits =0;
            for(uint32_t j=0; j < nHModule; ++j){
               for(uint32_t i =0; i < ncolumns; ++i){
                        ROCHits += rocColumnsHits[j* ncolumns + i];
                        BarrelColumnsNt_->Fill(idH, idL, layer, ladder, ring,gnHModuleBarrel_+j,j*ncolumns + i, rocColumnsHits[j* ncolumns + i], BarrelColumnsOffset_ + i + j*ncolumns);
                        if((i+1)%52 ==0){ // end of the chip 
                            ROCOccNt_->Fill(idH, idL, ROCHits,SubDet,gnHModuleBarrel_+j,ROCn , BarrelModuleHits, nROCs );
                            ++ROCn;
                            ROCHits=0;
                        }
               }
                        
	           for(uint32_t i =0; i < ncolumns/2; ++i)  BarrelDColumnsNt_->Fill(idH, idL,layer, ladder, ring,gnHModuleBarrel_+j,j*ncolumns/2 + i, rocDColumnsHits[j* ncolumns/2 +i], BarrelColumnsOffset_/2 + i + j*ncolumns/2);
               for(uint32_t i =0; i < ncolumns/2; ++i) DcolumnHits_->Fill(rocDColumnsHits[j* ncolumns/2 +i]);
	         }

             //BarrelModuleNt_->Fill(idH, idL,layer,ladder, ring,detZ,detR, detPhi, detEta, BarrelModuleHits, HModHits[0], HModHits[1], ncolumns, nrows);
             if(layer <3){
	           //BarrelSlinkHitsDistrib_->Fill(HModHits[0]);
                   if(layer==1) LinkOcc_[0]->Fill(HFRecHitSum,(double)(BarrelModuleHits / (ncolumns * nrows) * 100));
                   else if(layer==2) LinkOcc_[1]->Fill(HFRecHitSum,(double)(BarrelModuleHits / (ncolumns * nrows) * 100));
               //if(nHModule > 1) BarrelSlinkHitsDistrib_->Fill(HModHits[1]);           
	         }else if (layer ==3){
	           //BarrelSlinkHitsDistrib_->Fill(BarrelModuleHits);
                   LinkOcc_[2]->Fill(HFRecHitSum,(double)(BarrelModuleHits / (ncolumns * nrows) * 100));
             }

            // Occupancy_->Fill((double)(BarrelModuleHits / (ncolumns * nrows) * 100));
            // Occupancy_->Fill((double)detZ, (double)(BarrelModuleHits / (ncolumns * nrows) * 100));


            // BarrelColumnsOffset_ += ncolumns * nHModule;  //ncolumns/Hmodule * n HModules
            // gnHModuleBarrel_+=nHModule;
             
             delete [] rocColumnsHits;
             delete [] rocDColumnsHits;
             
         }

         //Selct the Endcap
         //---------------------------------------------------------------------------------------------------------
         if(detId.subdetId()==PixelSubdetector::PixelEndcap ){ 
           PXFDetId  fdetid(id);
           uint32_t side  = fdetid.side(); //size =1 for -z, 2 for +z
           uint32_t disk  = fdetid.disk(); //1, 2, 3
           uint32_t blade = fdetid.blade(); //1-23
           uint32_t panel = fdetid.panel();//panel =1,2
           uint32_t mod   = fdetid.module();
           
           uint32_t EndcapModuleHits=0, HModHits[2];
            HModHits[0] = 0;
            HModHits[1] = 0;     
             
             
             uint32_t* rocColumnsHits = new uint32_t [nHModule * ncolumns];
             uint32_t* rocDColumnsHits = new uint32_t [nHModule *  ncolumns/2];
             
             for(uint32_t j =0; j < nHModule; ++j){
	       for (uint32_t i = 0; i< ncolumns; ++i){
		 rocColumnsHits[j* ncolumns +i] = 0;
		 rocDColumnsHits[(j* ncolumns + i)/2] =0;
               }
             }
      
           // loop over module Digis
           //------------------------------------------------------------------
           for ( iter = begin ; iter != end; iter++ ){  
	          ++EndcapModuleHits;
                 
              uint32_t HModule = 0;
              if ((*iter).row() > 80) HModule =1;
	          ++rocColumnsHits[HModule * ncolumns + (*iter).column()];
		      ++rocDColumnsHits[(HModule * ncolumns +(*iter).column())/2];
              ++HModHits[HModule]; 
            //  EndcapDigisNt_->Fill(idH, idL,side, disk, blade,panel,mod,(*iter).adc(), (*iter).column(), (*iter).row());
           }
           
          // EndcapModuleNt_->Fill(idH, idL, side,disk, blade,panel, mod, detZ, detR,detPhi, detEta, EndcapModuleHits,HModHits[0], HModHits[1], ncolumns, nrows);
           //EndcapModuleNt_->Fill(idH, idL, side,disk, blade,panel, mod, detZ, detR,detPhi, detEta, HModHits[0], HModHits[1], ncolumns, nrows);
           if(disk==1) LinkOcc_[3]->Fill(HFRecHitSum,(double)(EndcapModuleHits / (ncolumns * nrows) * 100));
           else if(disk==2) LinkOcc_[4]->Fill(HFRecHitSum,(double)(EndcapModuleHits / (ncolumns * nrows) * 100));
           else if(disk==3) LinkOcc_[5]->Fill(HFRecHitSum,(double)(EndcapModuleHits / (ncolumns * nrows) * 100));

           
           Occupancy_->Fill((double)(EndcapModuleHits / (ncolumns * nrows) * 100));
           OccupancyZ_->Fill((double)detZ,(double)( EndcapModuleHits / (ncolumns * nrows) * 100));

           uint32_t ROCn=0, ROCHits =0;
           for(uint32_t j=0; j < nHModule; ++j){
	           for(uint32_t i =0; i < ncolumns; ++i){
                    ROCHits += rocColumnsHits[j* ncolumns + i];   
                    EndcapColumnsNt_->Fill(idH, idL,side,disk, blade, panel,gnHModuleEndcap_+j,j*ncolumns + i, rocColumnsHits[j* ncolumns + i], EndcapColumnsOffset_ + i + j*ncolumns);
                    if((i+1)%52 ==0){ // end of the chip 
                            ROCOccNt_->Fill(idH, idL,ROCHits,SubDet,gnHModuleEndcap_+j,ROCn, EndcapModuleHits, nROCs);
                            ++ROCn;
                            ROCHits=0;
                    }
               }
               for(uint32_t i =0; i < ncolumns/2; ++i)  EndcapDColumnsNt_->Fill(idH, idL,side, disk, blade, panel, gnHModuleEndcap_+j,j*ncolumns/2 + i, rocDColumnsHits[j* ncolumns/2 +i], EndcapColumnsOffset_/2 + i + j*ncolumns/2);
               for(uint32_t i =0; i < ncolumns/2; ++i) DcolumnHits_->Fill(rocDColumnsHits[j* ncolumns/2 +i]);
           } 


             EndcapColumnsOffset_ += ncolumns * nHModule;  //ncolumns/Hmodule * n HModules
             gnHModuleEndcap_+=nHModule;
             
             delete [] rocColumnsHits;
             delete [] rocDColumnsHits;
         }
 

    }*/
   
}

// ------------ method called once each job just before starting event loop  ------------
void 
SiPixelAnalyzer::beginJob(const edm::EventSetup&)
{
    
     oFile_ = new TFile((const char*)outputFile_.c_str(), "RECREATE");
     //FEDs studies---------------------------------------------------------------------
     FEDNt_ = new TNtuple("FED", "FED","fedid:links:roc",100000);
     LinksNt_= new TNtuple("Links", "Links","fedid:linkn:nHits",100000);
     ROCOccNt_= new TNtuple("ROCOcc", "ROCOcc","idH:idL:ROCHits:SubDet:gHModule:ROCn:ModuleHits:nROCs",100000);
     
     
     //Barrel Ntuples----------------------------------------------------------------------------
     BarrelModuleNt_ = new TNtuple("BarrelPixelDistrib", "BarrelPixelDistri","idH:idL:layer:ladder:ring:z:R:phi:eta:Mhits:HMod0Hits:HMod1Hits:ncolumns:nrows",100000);
     BarrelDigisNt_ = new TNtuple("BarrelDigis", "BarrelPixelDigisDist","idH:idL:layer:ladder:ring:adc:column:row",100000);
     BarrelColumnsNt_ = new TNtuple("BarrelColumnOcc", "BarrelColumnOcc","idH:idL:layer:ladder:ring:gHmodule:column:colHits:gcolumn",100000);
     BarrelDColumnsNt_ = new TNtuple("BarrelDColumnOcc", "BarrelDColumnOcc","idH:idL:layer:ladder:ring:gHmodule:Dcolumn:DcolHits:gDcolumn",100000);
   
     //histos---------------------------------------------------------
     BarrelSlinkHitsDistrib_ = new TH1F("BarrelSlinkHitsDistrib", "Slink Hits -Barrel-", 801, 0., 800.);
     BarrelSlinkHitsDistrib_->GetYaxis()->SetTitle("Entries");
     BarrelSlinkHitsDistrib_->GetXaxis()->SetTitle("Hits per RO Link");
     BarrelSlinkHitsDistrib_->SetDirectory(oFile_->GetDirectory(0));
 
    //Endcap Ntuples----------------------------------------------------------------------------
    EndcapModuleNt_ = new TNtuple("EndcapPixelDistrib", "EndcapPixelDistri","idH:idL:side:disk:blade:panel:module:z:R:phi:eta:HMod0Hits:HMod1Hits:ncolumns:nrows",100000);
    EndcapDigisNt_ = new TNtuple("EndcapDigis", "EndcapPixelDigisDist","idH:idL:side:disk:blade:panel:module:adc:column:row",100000);
    EndcapColumnsNt_ = new TNtuple("EndcapColumnOcc", "EndcapColumnOcc","idH:idL:side:disk:blade:panel:gHmodule:column:colHits:gcolumn",100000);
    EndcapDColumnsNt_ = new TNtuple("EndcapDColumnOcc", "EndcapDColumnOcc","idH:idL:side:disk:blade:panel:gHmodule:Dcolumn:DcolHits:gDcolumn",100000);

     //General histos---------------------------------------------------------
     DcolumnHits_ = new TH1F("DcolumnHits", "Hits per Dcolumn", 161, 0, 160);
     DcolumnHits_ ->GetYaxis()->SetTitle("Entries");
     DcolumnHits_->GetXaxis()->SetTitle("Hits per Dcolumn");
     DcolumnHits_->SetDirectory(oFile_->GetDirectory(0));
 
     Occupancy_ = new TH1D("Occupancy", "Occupancy", 3001, 0., 3.);
     Occupancy_->GetYaxis()->SetTitle("Entries");
     Occupancy_->GetXaxis()->SetTitle("Occupancy [%]");
     Occupancy_->SetDirectory(oFile_->GetDirectory(0));

     OccupancyZ_ = new TH2D("OccupancyZ", "Occupancy", 1020, -51., 51., 3001, 0., 3.);
     OccupancyZ_->GetYaxis()->SetTitle("Occupancy [%]");
     OccupancyZ_->GetXaxis()->SetTitle("Z");
     OccupancyZ_->SetDirectory(oFile_->GetDirectory(0));
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SiPixelAnalyzer::endJob() {
   oFile_->Write();
   oFile_->Close();
}


//define this as a plug-in
DEFINE_FWK_MODULE(SiPixelAnalyzer);
