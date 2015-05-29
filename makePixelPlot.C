#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TString.h"
#include "TDatime.h"
#include "TMath.h"
#include "TLatex.h"
#include "TStyle.h"

#include <string>
#include <vector>
#include <iostream>
#include <fstream>

void claverCanvasSaving(TCanvas* c, TString s,TString format="gif"){
  TDatime* date = new TDatime();
  c->SaveAs(Form("%s_%d.%s",s.Data(),date->GetDate(), format.Data()));
  return;
}

int makePixelPlot(const std::string inName)
{
  gStyle->SetOptStat(0);
  
  std::string outName = inName;
  const std::string cutString = "_HIST.root";
  const std::string addString = "_PLOT.root";

  std::size_t strIndex = outName.find(cutString);
  if(!(strIndex == std::string::npos)) outName.replace(strIndex, cutString.length(), addString);


  TFile* inFile_p = new TFile(inName.c_str(), "READ");

  const Int_t maxFedId = 39;
  const Int_t maxLink = 36;
  TCanvas* fedCanv_p[maxFedId];
  TCanvas* meanFedCanv_p;
  TCanvas* linkCanv_p[maxFedId];
  TCanvas* meanLinkCanv_p[maxFedId];

  TH1F* nHitsPerFed_h[maxFedId];
  TH1F* meanHitsPerFed_h;
  TH1F* nHitsPerLink_h[maxFedId][maxLink];
  TH1F* meanHitsPerLink_h[maxFedId];

  Float_t maxPerFed[maxFedId];
  Float_t minPerFed[maxFedId];

  for(Int_t iter = 0; iter < maxFedId; iter++){
    maxPerFed[iter] = -1;
    minPerFed[iter] = 1000000;
  }

  meanFedCanv_p = new TCanvas(Form("meanHitsPerFed_c"), Form("meanHitsPerFed_c"), 700, 700);
  meanHitsPerFed_h = (TH1F*)inFile_p->Get(Form("meanHitsPerFed_h"));
  meanHitsPerFed_h->DrawCopy("E1");

  for(Int_t iter = 0; iter < maxFedId; iter++){
    fedCanv_p[iter] = new TCanvas(Form("nHitsPerFed_fed%d_c", iter+1), Form("nHitsPerFed_fed%d_c", iter+1), 700, 700);
    nHitsPerFed_h[iter] = (TH1F*)inFile_p->Get(Form("fed%d/nHitsPerFed_fed%d_h", iter+1, iter+1));
    nHitsPerFed_h[iter]->DrawCopy("E1");

    meanLinkCanv_p[iter] = new TCanvas(Form("meanHitsPerLink_fed%d_c", iter+1), Form("meanHitsPerLink_fed%d_c", iter+1), 700, 700);
    meanHitsPerLink_h[iter] = (TH1F*)inFile_p->Get(Form("fed%d/meanHitsPerLink_fed%d_h", iter+1, iter+1));
    meanHitsPerLink_h[iter]->DrawCopy("E1");


    linkCanv_p[iter] = new TCanvas(Form("nHitsPerLink_fed%d_c", iter+1), Form("nHitsPerLink_fed%d_c", iter+1), 700, 700);
    linkCanv_p[iter]->Divide(9, 4, 0.0, 0.0);

    for(Int_t iter2 = 0; iter2 < maxLink; iter2++){
      nHitsPerLink_h[iter][iter2] = (TH1F*)inFile_p->Get(Form("fed%d/nHitsPerLink_fed%d_link%d_h", iter+1, iter+1, iter2+1));

      Float_t tempMax = nHitsPerLink_h[iter][iter2]->GetMaximum();
      Float_t tempMin = nHitsPerLink_h[iter][iter2]->GetMinimum();
      if(tempMax + TMath::Sqrt(tempMax) > maxPerFed[iter]) maxPerFed[iter] = tempMax + TMath::Sqrt(tempMax);
      if(tempMin - TMath::Sqrt(tempMin) < minPerFed[iter]) minPerFed[iter] = tempMin - TMath::Sqrt(tempMin);
						  
      linkCanv_p[iter]->cd(iter2+1);
      nHitsPerLink_h[iter][iter2]->DrawCopy("E1");
    }
    if(minPerFed[iter] < 0) minPerFed[iter] = 0;
  }

  for(Int_t iter = 0; iter < maxFedId; iter++){
    TLatex* temp = new TLatex();
    temp->SetNDC();
    temp->SetTextFont(43);
    temp->SetTextSizePixels(14);

    for(Int_t iter2 = 0; iter2 < maxLink; iter2++){
      linkCanv_p[iter]->cd(iter2+1);
      nHitsPerLink_h[iter][iter2]->SetMaximum(maxPerFed[iter]);
      nHitsPerLink_h[iter][iter2]->SetMinimum(minPerFed[iter]);
      nHitsPerLink_h[iter][iter2]->DrawCopy("E1");

      if(iter2 == 0) temp->DrawLatex(.15, .8, Form("fedId == %d", iter+1)); 
      temp->DrawLatex(.15, .9, Form("linkn == %d", iter2+1));
    }
  } 


  TFile* outFile_p = new TFile(outName.c_str(), "RECREATE");
  meanFedCanv_p->Write("", TObject::kOverwrite);
  claverCanvasSaving(meanFedCanv_p, Form("pdfDir/meanHitsPerFed"), "pdf");

  for(Int_t iter = 0; iter < maxFedId; iter++){
    fedCanv_p[iter]->Write("", TObject::kOverwrite);
    claverCanvasSaving(fedCanv_p[iter], Form("pdfDir/nHitsPerFed_fed%d", iter+1), "pdf");

    meanLinkCanv_p[iter]->Write("", TObject::kOverwrite);
    claverCanvasSaving(meanLinkCanv_p[iter], Form("pdfDir/meanHitsPerLink_fed%d", iter+1), "pdf");

    linkCanv_p[iter]->Write("", TObject::kOverwrite);
    claverCanvasSaving(linkCanv_p[iter], Form("pdfDir/nHitsPerLink_fed%d", iter+1), "pdf");
  }
  outFile_p->Close();
  delete outFile_p;

  delete meanFedCanv_p;
  for(Int_t iter = 0; iter < maxFedId; iter++){
    delete fedCanv_p[iter];
    delete meanLinkCanv_p[iter];
    delete linkCanv_p[iter];
  }

  inFile_p->Close();
  delete inFile_p;

  return 0;
}

int main(int argc, char *argv[])
{
  if(argc != 2){
    std::cout << "Usage: makePixelPlot <inName>" << std::endl;
    std::cout << "argNum: " << argc << std::endl;
    for(Int_t iter = 0; iter < argc; iter++){
      std::cout << "arg " << iter << ": " << argv[iter] << std::endl;
    }
  }

  int rStatus = -1;
  rStatus = makePixelPlot(argv[1]);
  return rStatus;
}
