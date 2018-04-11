#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TMath.h"

//-----------------------------------------------------------------------------------------------------------------------------

//from Matt's Skimming Branch of LHCb Analysis

//# Data
//itype=11   sqrts=7   year=2011  name=Data2011         fname=/r01/lhcb/mkenzie/Bc2Dmunu/Data2011.root  tname=AnalysisTree
//itype=12   sqrts=8   year=2012  name=Data2012         fname=/r01/lhcb/mkenzie/Bc2Dmunu/Data2012.root  tname=AnalysisTree
//itype=15   sqrts=13  year=2015  name=Data2015         fname=/r01/lhcb/mkenzie/Bc2Dmunu/Data2015.root  tname=AnalysisTree
//itype=16   sqrts=13  year=2016  name=Data2016         fname=/r01/lhcb/mkenzie/Bc2Dmunu/Data2016.root  tname=AnalysisTree
//itype=17   sqrts=13  year=2017  name=Data2017         fname=/r01/lhcb/mkenzie/Bc2Dmunu/Data2017.root  tname=AnalysisTree
//
//# Norm Data
//#itype=111   sqrts=7  year=2011  name=DataNorm2011     fname=/r01/lhcb/mkenzie/Bc2Dmunu/DataNorm2011.root  tname=AnalysisTree
//#itype=112   sqrts=8  year=2012  name=DataNorm2012     fname=/r01/lhcb/mkenzie/Bc2Dmunu/DataNorm2012.root  tname=AnalysisTree
//#itype=115   sqrts=13 year=2015  name=DataNorm2015     fname=/r01/lhcb/mkenzie/Bc2Dmunu/DataNorm2015.root  tname=AnalysisTree
//#itype=116   sqrts=13 year=2016  name=DataNorm2016     fname=/r01/lhcb/mkenzie/Bc2Dmunu/DataNorm2016.root  tname=AnalysisTree
//#itype=117   sqrts=13 year=2017  name=DataNorm2017     fname=/r01/lhcb/mkenzie/Bc2Dmunu/DataNorm2017.root  tname=AnalysisTree
//
//
//# MC
//# SIG
//itype=-20 sqrts=8  year=2012  name=MCBc2DMuNu         fname=/r01/lhcb/mkenzie/Bc2Dmunu/MCBc2DMuNu.root tname=AnalysisTree
//# BKG
//itype=-29 sqrts=8 year=2012 name=MCBu2DMuNu fname=/r01/lhcb/mkenzie/Bc2Dmunu/MCBu2DMuNu.root tname=AnalysisTree
//-----------------------------------------------------------------------------------------------------------------------------


//load, read, hists
void BDTscan()
{
    //from macro, use 'new' to keep TTree and TFile alive after function
    TFile *f = new TFile("/var/pcfst/r01/lhcb/mkenzie/Bc2Dmunu/AnalysisOut.root");
    TTree *AnalysisTree = (TTree*)f->Get("AnalysisTree");

    //declare branches of interest
    Double_t Bmcorr, Bmcorrerr, bdt;
    Int_t Mucharge, Kcharge;
   
    AnalysisTree->SetBranchStatus("*", 0);
    AnalysisTree->SetBranchStatus("bu_rejection_bdtoutput",1);
    AnalysisTree->SetBranchStatus("B_plus_MCORR",1);
    AnalysisTree->SetBranchStatus("B_plus_MCORRERR",1);
    AnalysisTree->SetBranchStatus("Mu_plus_ID",1);
    AnalysisTree->SetBranchStatus("K_minus_ID",1);
    AnalysisTree->SetBranchStatus("itype",1);
    AnalysisTree->SetBranchStatus("B_plus_MCORR",1);


    AnalysisTree->SetBranchAddress("bu_rejection_bdtoutput",&bdt);
    AnalysisTree->SetBranchAddress("B_plus_MCORR",&Bmcorr);
    AnalysisTree->SetBranchAddress("B_plus_MCORRERR",&Bmcorrerr);
    AnalysisTree->SetBranchAddress("Mu_plus_ID",&Mucharge);
    AnalysisTree->SetBranchAddress("K_minus_ID",&Kcharge);

    
    Double_t nentries = AnalysisTree->GetEntries();
    Double_t totSig, partialSig, EffSig, BDT_cutval, lowMassVal, highMassVal, maxBDT, minBDT, significance;
    
    //sigmas
    significance=5.0;
    
    //signal MC
    totSig =  AnalysisTree->GetEntries("itype==-20");

    //signal mass window
    lowMassVal=4500.0;
    highMassVal=10000.0;
    
    //BDT extrema - note that this assumer BDT range [-1, +1]
    maxBDT=+1.;
    minBDT=-1;

    //sloppy 
    Int_t target_nentries;
    target_nentries = 100;

    Double_t increment;
    increment = (maxBDT-minBDT)/target_nentries;
    

    Double_t BDTvals[target_nentries];
    Double_t FoMvals[target_nentries];
    Double_t SigEffs[target_nentries];
    
    Double_t bkg, punzi;
    
    for (Int_t i=0; i<target_nentries; i++){
        
        BDT_cutval = minBDT+(i*increment);
        BDTvals[i] = BDT_cutval;
    

        partialSig = AnalysisTree->GetEntries(Form("itype==-20 && bu_rejection_bdtoutput>(%f)", BDT_cutval));
        
        //note that particle ID convention means oppositely charged hadron lepton required posive
        //branch_ID product
        bkg =  AnalysisTree->GetEntries(Form("(itype==11 || itype==12) && B_plus_MCORR > %f && B_plus_MCORR < %f && bu_rejection_bdtoutput > %f && ((K_minus_ID*Mu_plus_ID)>0)", lowMassVal, highMassVal, BDT_cutval ) ); 
        
        
        EffSig = partialSig/totSig;
        punzi = EffSig/(significance/2 + sqrt(bkg));
        
        FoMvals[i] = punzi;
        SigEffs[i] = EffSig;
    }

    TCanvas *c1 = new TCanvas("c1","Evaluation of Punzi FoM vs BDT");
    //c1->SetCanvasSize(600, 400);
    //c1->SetWindowSize(200, 200);
    //c1->Divide(2,1);
    //c1->cd(1);
    
    TGraph *gr1 = new TGraph (target_nentries, BDTvals, FoMvals);
    gr1->Draw("ALP*");
    gr1->SetMarkerStyle(21);
    gr1->SetMarkerColor(38);
    gr1->SetTitle("Punzi FoM vs BDT Output");
    gr1->GetYaxis()->SetTitle("#frac{#epsilon_{s} }{ #frac{a}{2} + #sqrt{N_{B}} }");
    gr1->GetXaxis()->SetTitle("BDT Output [Arbitrary Units]");
    c1->SaveAs("FoM.pdf");

    TCanvas *c2 = new TCanvas("c2","Signal Effiency vs BDT");
    TGraph *gr2 = new TGraph (target_nentries, BDTvals, SigEffs);
    gr2->Draw("ALP*");
    gr2->SetMarkerStyle(21);
    gr2->SetMarkerColor(38);
    gr2->SetTitle("#epsilon_{s} vs BDT Output");
    gr2->GetYaxis()->SetTitle("#epsilon_{s}");
    gr2->GetXaxis()->SetTitle("BDT Output [Arbitrary Units]");
    c2->SaveAs("Efficiency.C");



//https://root.cern.ch/doc/master/classTColor.html
}
