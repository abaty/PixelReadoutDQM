#include "TNtuple.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TH1F.h"

#include <string>
#include <vector>
#include <iostream>
#include <fstream>

int makePixelHist(const std::string inName)
{
  std::string buffer;
  std::vector<std::string> listOfFiles;
  Int_t nLines = 0;
  std::ifstream inList(inName.data());

  if(!inList.is_open()){
    std::cout << "Error opening file. Exiting." <<std::endl;
    return 1;
  }
  else{
    while(true){
      inList >> buffer;
      if(inList.eof()) break;
      listOfFiles.push_back(buffer);
      nLines++;
    }
  }

  std::string outName = inName;
  const std::string cutString = ".txt";
  const std::string addString = "_HIST.root";

  std::size_t strIndex = outName.find(cutString);
  if(!(strIndex == std::string::npos)) outName.replace(strIndex, cutString.length(), addString);

  TFile* outFile_p = new TFile(outName.c_str(), "RECREATE");

  const Int_t maxFedId = 39;
  const Int_t maxLink = 36;
  TH1F* nHitsPerFed_h[maxFedId];
  TH1F* nHitsPerLink_h[maxFedId][maxLink];

  for(Int_t iter = 0; iter < maxFedId; iter++){
    nHitsPerFed_h[iter] = new TH1F(Form("nHitsPerFed_fed%d_h", iter+1), Form(";nHits_{fed%d};Events", iter+1), 20, 0.001, 7999.999);
    nHitsPerFed_h[iter]->GetXaxis()->CenterTitle();

    for(Int_t iter2 = 0; iter2 < maxLink; iter2++){
      nHitsPerLink_h[iter][iter2] = new TH1F(Form("nHitsPerLink_fed%d_link%d_h", iter+1, iter2+1), Form(";nHits_{fed%d,link%d};Events", iter+1, iter2+1), 20, 0.001, 399.999);
      nHitsPerLink_h[iter][iter2]->GetXaxis()->CenterTitle();
    }
  }


  for(Int_t fileIter = 0; fileIter < nLines; fileIter++){
    std::cout << "Running file " << listOfFiles[fileIter] << std::endl;
    TFile* inFile_p = new TFile(listOfFiles[fileIter].c_str(), "READ");
    TNtuple* links_p = (TNtuple*)inFile_p->Get("Links");
  
    Float_t fedid; 
    Float_t linkn; 
    Float_t nHits; 
    
    links_p->SetBranchAddress("fedid", &fedid);
    links_p->SetBranchAddress("linkn", &linkn);
    links_p->SetBranchAddress("nHits", &nHits);
    

    const Long64_t nEntries = links_p->GetEntries();

    Float_t fedCounter = 0;

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
      
      nHitsPerLink_h[(Int_t(fedid))-1][(Int_t(linkn))-1]->Fill(nHits);
      if(linkn != 36) fedCounter += nHits;
      else{
	nHitsPerFed_h[(Int_t(fedid))-1]->Fill(fedCounter);
	fedCounter = 0;
      }
    }

    inFile_p->Close();
    delete inFile_p;
  }    

  for(Int_t iter = 0; iter < maxFedId; iter++){
    TDirectory* tempDir_p = outFile_p->mkdir(Form("fed%d", iter+1));
    tempDir_p->cd();
    nHitsPerFed_h[iter]->Write("", TObject::kOverwrite);
    for(Int_t iter2 = 0; iter2 < maxLink; iter2++){
      nHitsPerLink_h[iter][iter2]->Write("", TObject::kOverwrite);
    }
  }
  
  for(Int_t iter = 0; iter < maxFedId; iter++){
    delete nHitsPerFed_h[iter];
    for(Int_t iter2 = 0; iter2 < maxLink; iter2++){
      delete nHitsPerLink_h[iter][iter2];
    }
  }
  
  outFile_p->Close();
  delete outFile_p;

  return 0;
}


int main(int argc, char *argv[])
{
  if(argc != 2){
    std::cout << "Usage: makePixelHist <inName>" << std::endl;
    std::cout << "argNum: " << argc << std::endl;
    for(Int_t iter = 0; iter < argc; iter++){
      std::cout << "arg " << iter << ": " << argv[iter] << std::endl;
    }
  }

  int rStatus = -1;
  rStatus = makePixelHist(argv[1]);
  return rStatus;
}
