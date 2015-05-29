#include "TNtuple.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TH1F.h"

#include <string>
#include <iostream>

void makePixelHist(const std::string inName, const std::string outName)
{
  TFile* inFile_p = new TFile(inName.c_str(), "READ");
  TNtuple* links_p = (TNtuple*)inFile_p->Get("Links");

  Float_t fedid; 
  Float_t linkn; 
  Float_t nHits; 

  links_p->SetBranchAddress("fedid", &fedid);
  links_p->SetBranchAddress("linkn", &linkn);
  links_p->SetBranchAddress("nHits", &nHits);

  const Int_t maxFedId = 39;
  const Int_t maxLink = 36;
  TH1F* nHits_h[maxFedId][maxLink];

  for(Int_t iter = 0; iter < maxFedId; iter++){
    for(Int_t iter2 = 0; iter2 < maxLink; iter2++){
      nHits_h[iter][iter2] = new TH1F(Form("nHits_fed%d_link%d_h", iter+1, iter2+1), Form(";nHits_{fed%d,link%d};Events", iter+1, iter2+1), 100, 0, 400);
      nHits_h[iter][iter2]->GetXaxis()->CenterTitle();
    }
  }

  const Long64_t nEntries = links_p->GetEntries();

  for(Long64_t entry = 0; entry < nEntries; entry++){
    links_p->GetEntry(entry);

    if(fedid < 1 || fedid > maxFedId){
      std::cout << "ERROR: fedid out of bounds" << std::endl;
      break;
    }

    if(linkn < 1 || linkn > maxLink){
      std::cout << "ERROR: linkn out of bounds" << std::endl;
      break;
    }

    nHits_h[(Int_t(fedid))-1][(Int_t(linkn))-1]->Fill(nHits);
  }

  TFile* outFile_p = new TFile(outName.c_str(), "RECREATE");
  for(Int_t iter = 0; iter < maxFedId; iter++){
    TDirectory* tempDir_p = outFile_p->mkdir(Form("fed%d", iter+1));
    tempDir_p->cd();
    for(Int_t iter2 = 0; iter2 < maxLink; iter2++){
      nHits_h[iter][iter2]->Write("", TObject::kOverwrite);
    }
  }
  outFile_p->Close();
  delete outFile_p;

  for(Int_t iter = 0; iter < maxFedId; iter++){
    for(Int_t iter2 = 0; iter2 < maxLink; iter2++){
      delete nHits_h[iter][iter2];
    }
  }

  inFile_p->Close();
  delete inFile_p;

  return;
}
