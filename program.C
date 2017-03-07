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

void Fitting(TH1D* h) {
    TF1 *fit = new TF1("fit0","gaus",0,100);
    //     TF1 *fit = new TF1("fit0","[0]*exp(-0.5*((x-[1])/[2])**2)",0,100);
    fit -> SetParameter(0,5);
    fit -> SetParameter(1,1.8);
    h -> Fit(fit,"L");
    mean = fit -> GetParameter(1);
    meanError = fit -> GetParError(1);
    cout<<mean<<endl;
    cout<<meanError<<endl;
}


void program() {
//     TFile* fin = new TFile("results_vpd1.root" ,"r");
//     TFile* fin = new TFile("results_bht1.root" ,"r");
    
    TFile* fin = new TFile("results_all.root" ,"r");
    

    
    Float_t flag, mass_hs;
    
    TH1D *hUS = (TH1D*)fin -> Get("h_inv_mass_US"); 
    TH1D *hLS = (TH1D*)fin -> Get("h_inv_mass_LS"); 
    hUS -> Rebin(3);
    hLS -> Rebin(3);

    TTree *t3 = (TTree*)f->Get("t_output");
    t_output -> Scan();
    
    
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
    
    Double_t hUSIntegral = hUS -> Integral(hUS -> FindBin(1.65), hUS -> FindBin(1.8), "");
    Double_t hLSIntegral = hLS -> Integral(hLS -> FindBin(1.65), hLS -> FindBin(1.8), "");
    
    Double_t hUSIntegral = hUSIntegral + (hUS -> Integral(hUS -> FindBin(1.9), hUS -> FindBin(2), ""));
    Double_t hLSIntegral = hLSIntegral + (hLS -> Integral(hLS -> FindBin(1.9), hLS -> FindBin(2), ""));
        

    
//     Double_t hUSIntegral = hUS -> Integral(hUS -> FindBin(1.5), hUS -> FindBin(1.7), "") + hUS -> Integral(hUS -> FindBin(1.99), hUS -> GetNbinsX(), "");
//     Double_t hLSIntegral = hLS -> Integral(hLS -> FindBin(1.5), hLS -> FindBin(1.7), "") + hLS -> Integral(hLS -> FindBin(1.99), hLS -> GetNbinsX(), "");
//     
    hLS -> Scale(hUSIntegral/hLSIntegral);
    hUS -> Add(hLS, -1);
  
    TF1 *fit = new TF1("fit0","gaus",1.8,1.9);
        TF1 *fit = new TF1("fit0","[0]*exp(-0.5*((x-[1])/[2])**2)",1.75,1.95);
    fit -> SetParameter(0,5000);
    fit -> SetParameter(1,1.86);
    fit -> SetParameter(2,0.06);
    
    hUS -> Fit(fit,"R");
/*
    cout<<fit -> GetParameter(1)<<endl;
    cout<<fit -> GetParError(1)<<endl;
    */
    
    hUS -> Draw();
//     hLS -> SetLineColor(45);
//     hLS -> Draw("same");
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
}