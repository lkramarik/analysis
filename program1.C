#include "TFile.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TH2F.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TROOT.h"
#include <iostream>
using namespace std;

void Fitting(TH1F* h) {
    TF1 *fit = new TF1("fit0","gaus",0,100);
    //     TF1 *fit = new TF1("fit0","[0]*exp(-0.5*((x-[1])/[2])**2)",0,100);
    fit -> SetParameter(0,5);
    fit -> SetParameter(1,50);
    h -> Fit(fit,"L");
    mean = fit -> GetParameter(1);
    meanError = fit -> GetParError(1);
    cout<<mean<<endl;
    cout<<meanError<<endl;
}


void program() {
    TFile* fin = new TFile("results.root" ,"r");

    
    Float_t flag, mass_hs;
    
    TH1D *hUS = (TH1D*)fin -> Get("h_inv_mass_US"); 
    TH1D *hLS = (TH1D*)fin -> Get("h_inv_mass_LS"); 
    hUS -> Rebin(10);
    hLS -> Rebin(10);

    
//     Long64_t number = tripleTree -> GetEntries(); //pocet vstupov v NTuple ("zozname" trojic)
//     for (Int_t i = 0; i < number; i++) { //cyklus, ktory postupne prejde vsetky vstupy v zozname
//         tripleTree -> GetEntry(i);   //vyberem jeden vstup, ktory ma cislo "i"
//         if ((flag == 0)||(flag == 1)) {
//             h1 -> Fill(mass_hs); 
//         }
//         if ((flag == 2)||(flag == 3)||(flag == 4)||(flag == 5)) {
//             h2 -> Fill(mass_hs);   
//             
//         }
//         
//         
//         
//     }
// //     h1 -> Add(h2, -1);
// 
//     Double_t h1Integral = h1 -> Integral(1, h1 -> FindBin(1.8), "") + h1 -> Integral(h1 -> FindBin(1.95), h1 -> GetNbinsX(), "");
//     Double_t h2Integral = h2 -> Integral(1, h2 -> FindBin(1.8), "") + h2 -> Integral(h1 -> FindBin(1.95), h2 -> GetNbinsX(), "");
//     
//     h2 -> Scale(h1Integral/h2Integral);
//         h1 -> Add(h2, -1);

    hUS -> Draw();
    hLS -> SetLineColor(45);
    hLS -> Draw("same");
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
}