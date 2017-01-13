#include "TFile.h"
#include "TDirectoryFile.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TDatime.h"
#include "TKey.h"
#include "TLine.h"
#include "TAxis.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TMath.h"
#include "TStyle.h"
#include "TF1.h"

#include "TTree.h"
#include "TDirectory.h"
#include "TLorentzVector.h"
#include <math.h>

#include <fstream>

#include <string>
#include <iostream>
#include <vector>

#define	centbin 4
#define trkptbin 8

char saythis[500];

TString cent[centbin+1] = {"0","10","30","50","100"};
TString trkpt[trkptbin+1] = {"0p7","1","2","3","4","8","12","16","20"};

Double_t int_bin_bounds[9] = {0.7,1.,2.,3.,4.,8.,12.,16.,20.};
int integral_bins = 9;

void JFF_unc_old(){

	TFile *file_q_rg = TFile::Open("/home/dhanush/Documents/Kurt_histos/PbPb_5TeV_MC_PythiaHydjet_MixHistos_RecoGenReduced_QuarkJet_sube0_Merged_noMix_newJetTrackCorrections_fineBin.root");
	TFile *file_g_rg = TFile::Open("/home/dhanush/Documents/Kurt_histos/PbPb_5TeV_MC_PythiaHydjet_MixHistos_RecoGenReduced_GluonJet_sube0_Merged_noMix_newJetTrackCorrections_fineBin.root");
    TFile *file_rg = TFile::Open("/home/dhanush/Documents/Kurt_histos/PbPb_5TeV_MC_PythiaHydjet_MixHistos_RecoGenReduced_sube0_Merged_refpt_newJetTrackCorrections_fineBin.root"); 

    TH2D *q_RecoJet_GenTrack_hJetTrackSignalBackground[4][10];
    TH2D *g_RecoJet_GenTrack_hJetTrackSignalBackground[4][10];
    TH2D *RecoJet_GenTrack_hJetTrackSignalBackground[4][10];

    TH1D *q_RecoJet_GenTrack_hJetTrackSignalBackground_eta[4][10];
    TH1D *g_RecoJet_GenTrack_hJetTrackSignalBackground_eta[4][10];
    TH1D *RecoJet_GenTrack_hJetTrackSignalBackground_eta[4][10];

    TH1D *RecoJet_GenTrack_hJetTrackME_eta[4][10];
    TH1D *q_RecoJet_GenTrack_hJetTrackME_eta[4][10];
    TH1D *g_RecoJet_GenTrack_hJetTrackME_eta[4][10]; 

    TFile *file_q_gg = TFile::Open("/home/dhanush/Documents/Kurt_histos/PbPb_5TeV_MC_PythiaHydjet_MixHistos_GenGenReduced_QuarkJet_sube0_Merged_noMix_newJetTrackCorrections_fineBin.root");
	TFile *file_g_gg = TFile::Open("/home/dhanush/Documents/Kurt_histos/PbPb_5TeV_MC_PythiaHydjet_MixHistos_GenGenReduced_GluonJet_sube0_Merged_noMix_newJetTrackCorrections_fineBin.root");
    TFile *file_gg = TFile::Open("/home/dhanush/Documents/Kurt_histos/PbPb_5TeV_MC_PythiaHydjet_MixHistos_GenGenReduced_sube0_Merged_refpt_newJetTrackCorrections_fineBin.root");

    //TFile f_yield("/home/dhanush/Documents/JetTrack2016/subtractedyield_integral.root", "RECREATE"); 

    TH2D *q_GenJet_GenTrack_hJetTrackSignalBackground[4][10];
    TH2D *g_GenJet_GenTrack_hJetTrackSignalBackground[4][10];
    TH2D *GenJet_GenTrack_hJetTrackSignalBackground[4][10];

    TH1D *q_GenJet_GenTrack_hJetTrackSignalBackground_eta[4][10];
    TH1D *g_GenJet_GenTrack_hJetTrackSignalBackground_eta[4][10];
    TH1D *GenJet_GenTrack_hJetTrackSignalBackground_eta[4][10];
 
    TH1D *recogen_gengen[4][10];
    TH1D *recogen_gengen_q[4][10];
    TH1D *recogen_gengen_g[4][10];

    TH1D *hintegral_eta[4]; 
    TH1D *hintegral_eta_q[4];
    TH1D *hintegral_eta_g[4];

    Double_t integral_all[4][10];
    Double_t integral_q[4][10];
    Double_t integral_g[4][10];

    TH1F *reco_jets[4];
    TH1F *q_reco_jets[4];
    TH1F *g_reco_jets[4];
    TH1F *gen_jets[4];
    TH1F *q_gen_jets[4];
    TH1F *g_gen_jets[4];

    for(int ibin=0;ibin<4;ibin++){

      Double_t integral_reco = 0;
      Double_t integral_reco_q = 0;
      Double_t integral_reco_g = 0;
      Double_t integral_gen = 0;
      Double_t integral_gen_q = 0;
      Double_t integral_gen_g = 0;  
      
      sprintf(saythis,"hintegral_eta_%d",ibin);
      hintegral_eta[ibin] = new TH1D(saythis,"",integral_bins-1,int_bin_bounds);
      hintegral_eta[ibin]->Sumw2();

      sprintf(saythis,"hintegral_eta_q_%d",ibin);
      hintegral_eta_q[ibin] = new TH1D(saythis,"",integral_bins-1,int_bin_bounds);
      hintegral_eta_q[ibin]->Sumw2();  

      sprintf(saythis,"hintegral_eta_g_%d",ibin);
      hintegral_eta_g[ibin] = new TH1D(saythis,"",integral_bins-1,int_bin_bounds);
      hintegral_eta_g[ibin]->Sumw2();

      sprintf(saythis,"reco_jets_%d",ibin);
      reco_jets[ibin] = (TH1F*)file_rg->Get((TString)("RecoJet_GenTrack_all_jets_corrpTCent"+cent[ibin]+"_Cent"+cent[ibin+1]+"_Pt100_Pt300"))->Clone(saythis);
      //reco_jets[ibin]->Scale(1./reco_jets[ibin]->GetXaxis()->GetBinWidth(1));
      integral_reco = reco_jets[ibin]->Integral();

      sprintf(saythis,"q_reco_jets_%d",ibin);
      q_reco_jets[ibin] = (TH1F*)file_q_rg->Get((TString)("RecoJet_GenTrack_all_jets_corrpTCent"+cent[ibin]+"_Cent"+cent[ibin+1]+"_Pt100_Pt300"))->Clone(saythis);
      //q_reco_jets[ibin]->Scale(1./q_reco_jets[ibin]->GetXaxis()->GetBinWidth(1));
      integral_reco_q = q_reco_jets[ibin]->Integral(); 

      sprintf(saythis,"g_reco_jets_%d",ibin);
      g_reco_jets[ibin] = (TH1F*)file_g_rg->Get((TString)("RecoJet_GenTrack_all_jets_corrpTCent"+cent[ibin]+"_Cent"+cent[ibin+1]+"_Pt100_Pt300"))->Clone(saythis);
      //g_reco_jets[ibin]->Scale(1./g_reco_jets[ibin]->GetXaxis()->GetBinWidth(1));
      integral_reco_g = g_reco_jets[ibin]->Integral(); 

      sprintf(saythis,"gen_jets_%d",ibin);  
      gen_jets[ibin] = (TH1F*)file_gg->Get((TString)("GenJet_GenTrack_all_jets_corrpTCent"+cent[ibin]+"_Cent"+cent[ibin+1]+"_Pt100_Pt300"))->Clone(saythis);
      //gen_jets[ibin]->Scale(1./gen_jets[ibin]->GetXaxis()->GetBinWidth(1));
      integral_gen = gen_jets[ibin]->Integral(); 

      sprintf(saythis,"q_gen_jets_%d",ibin);  
      q_gen_jets[ibin] = (TH1F*)file_q_gg->Get((TString)("GenJet_GenTrack_all_jets_corrpTCent"+cent[ibin]+"_Cent"+cent[ibin+1]+"_Pt100_Pt300"))->Clone(saythis);
      //q_gen_jets[ibin]->Scale(1./q_gen_jets[ibin]->GetXaxis()->GetBinWidth(1)); 
      integral_gen_q = q_gen_jets[ibin]->Integral();

      sprintf(saythis,"g_gen_jets_%d",ibin);  
      g_gen_jets[ibin] = (TH1F*)file_g_gg->Get((TString)("GenJet_GenTrack_all_jets_corrpTCent"+cent[ibin]+"_Cent"+cent[ibin+1]+"_Pt100_Pt300"))->Clone(saythis);
      //g_gen_jets[ibin]->Scale(1./g_gen_jets[ibin]->GetXaxis()->GetBinWidth(1));  
      integral_gen_q = q_gen_jets[ibin]->Integral();  

      for(int ibin3=0;ibin3<8;ibin3++){

        int y1 = 0;
        int y2 = 0;

        // reco jet gen tracks

        sprintf(saythis,"RecoJet_GenTrack_hJetTrackSignalBackground_cent%d_trkpt%d",ibin,ibin3);
        RecoJet_GenTrack_hJetTrackSignalBackground[ibin][ibin3] = (TH2D*)file_rg->Get((TString)("RecoJet_GenTrack_hJetTrackSignalBackgroundCent"+cent[ibin]+"_Cent"+cent[ibin+1]+"_Pt100_Pt300_TrkPt"+trkpt[ibin3]+"_TrkPt"+trkpt[ibin3+1]))->Clone(saythis);
        
        y1 = RecoJet_GenTrack_hJetTrackSignalBackground[ibin][ibin3]->GetYaxis()->FindBin(TMath::Pi()/2);
        y2 = RecoJet_GenTrack_hJetTrackSignalBackground[ibin][ibin3]->GetYaxis()->FindBin(3*TMath::Pi()/2);
/*
        sprintf(saythis,"RecoJet_GenTrack_hJetTrackME_eta_cent%d_trkpt%d",ibin,ibin3);
        RecoJet_GenTrack_hJetTrackME_eta[ibin][ibin3] = RecoJet_GenTrack_hJetTrackSignalBackground[ibin][ibin3] -> ProjectionX("",y1,y2);
        RecoJet_GenTrack_hJetTrackME_eta[ibin][ibin3]->Scale(1./RecoJet_GenTrack_hJetTrackME_eta[ibin][ibin3]->GetMaximum());
*/
        sprintf(saythis,"RecoJet_GenTrack_hJetTrackSignalBackground_eta_cent%d_trkpt%d",ibin,ibin3);
        RecoJet_GenTrack_hJetTrackSignalBackground_eta[ibin][ibin3] = RecoJet_GenTrack_hJetTrackSignalBackground[ibin][ibin3] -> ProjectionX();   
        //RecoJet_GenTrack_hJetTrackSignalBackground_eta[ibin][ibin3]->Sumw2();
        RecoJet_GenTrack_hJetTrackSignalBackground_eta[ibin][ibin3]->Scale(1./RecoJet_GenTrack_hJetTrackSignalBackground[ibin][ibin3]->GetXaxis()->GetBinWidth(1));
        //RecoJet_GenTrack_hJetTrackSignalBackground_eta[ibin][ibin3]->Scale(1./RecoJet_GenTrack_hJetTrackSignalBackground[ibin][ibin3]->GetXaxis()->GetBinWidth(1));
        
        RecoJet_GenTrack_hJetTrackSignalBackground_eta[ibin][ibin3]->Scale(1./integral_reco);   
        //RecoJet_GenTrack_hJetTrackSignalBackground_eta[ibin][ibin3]->Write();
        RecoJet_GenTrack_hJetTrackSignalBackground_eta[ibin][ibin3]->SetLineColor(kBlack);

        sprintf(saythis,"q_RecoJet_GenTrack_hJetTrackSignalBackground_cent%d_trkpt%d",ibin,ibin3);
        q_RecoJet_GenTrack_hJetTrackSignalBackground[ibin][ibin3] = (TH2D*)file_q_rg->Get((TString)("RecoJet_GenTrack_hJetTrackSignalBackgroundCent"+cent[ibin]+"_Cent"+cent[ibin+1]+"_Pt100_Pt300_TrkPt"+trkpt[ibin3]+"_TrkPt"+trkpt[ibin3+1]))->Clone(saythis);
/*
        sprintf(saythis,"q_RecoJet_GenTrack_hJetTrackME_eta_cent%d_trkpt%d",ibin,ibin3);
        q_RecoJet_GenTrack_hJetTrackME_eta[ibin][ibin3] = q_RecoJet_GenTrack_hJetTrackSignalBackground[ibin][ibin3] -> ProjectionX("",y1,y2);
        q_RecoJet_GenTrack_hJetTrackME_eta[ibin][ibin3]->Scale(1./q_RecoJet_GenTrack_hJetTrackME_eta[ibin][ibin3]->GetMaximum());
*/
        sprintf(saythis,"q_RecoJet_GenTrack_hJetTrackSignalBackground_eta_cent%d_trkpt%d",ibin,ibin3);
        q_RecoJet_GenTrack_hJetTrackSignalBackground_eta[ibin][ibin3] = q_RecoJet_GenTrack_hJetTrackSignalBackground[ibin][ibin3] -> ProjectionX();  
        //q_RecoJet_GenTrack_hJetTrackSignalBackground_eta[ibin][ibin3]->Sumw2();
        q_RecoJet_GenTrack_hJetTrackSignalBackground_eta[ibin][ibin3]->Scale(1./q_RecoJet_GenTrack_hJetTrackSignalBackground[ibin][ibin3]->GetXaxis()->GetBinWidth(1));  
        
        q_RecoJet_GenTrack_hJetTrackSignalBackground_eta[ibin][ibin3]->Scale(1./integral_reco_q);
        //q_RecoJet_GenTrack_hJetTrackSignalBackground_eta[ibin][ibin3]->Write();
        q_RecoJet_GenTrack_hJetTrackSignalBackground_eta[ibin][ibin3]->SetLineColor(kBlue);
        q_RecoJet_GenTrack_hJetTrackSignalBackground_eta[ibin][ibin3]->GetXaxis()->SetTitle("dEta");
        q_RecoJet_GenTrack_hJetTrackSignalBackground_eta[ibin][ibin3]->GetXaxis()->SetTitleSize(0.05); 

        sprintf(saythis,"g_RecoJet_GenTrack_hJetTrackSignalBackground_cent%d_trkpt%d",ibin,ibin3);
        g_RecoJet_GenTrack_hJetTrackSignalBackground[ibin][ibin3] = (TH2D*)file_g_rg->Get((TString)("RecoJet_GenTrack_hJetTrackSignalBackgroundCent"+cent[ibin]+"_Cent"+cent[ibin+1]+"_Pt100_Pt300_TrkPt"+trkpt[ibin3]+"_TrkPt"+trkpt[ibin3+1]))->Clone(saythis);
/*
        sprintf(saythis,"g_RecoJet_GenTrack_hJetTrackME_eta_cent%d_trkpt%d",ibin,ibin3);
        g_RecoJet_GenTrack_hJetTrackME_eta[ibin][ibin3] = g_RecoJet_GenTrack_hJetTrackSignalBackground[ibin][ibin3] -> ProjectionX("",y1,y2);
        g_RecoJet_GenTrack_hJetTrackME_eta[ibin][ibin3]->Scale(1./g_RecoJet_GenTrack_hJetTrackME_eta[ibin][ibin3]->GetMaximum()); 
*/
        sprintf(saythis,"g_RecoJet_GenTrack_hJetTrackSignalBackground_eta_cent%d_trkpt%d",ibin,ibin3); 
        g_RecoJet_GenTrack_hJetTrackSignalBackground_eta[ibin][ibin3] = g_RecoJet_GenTrack_hJetTrackSignalBackground[ibin][ibin3] -> ProjectionX();
        //g_RecoJet_GenTrack_hJetTrackSignalBackground_eta[ibin][ibin3]->Sumw2();
        g_RecoJet_GenTrack_hJetTrackSignalBackground_eta[ibin][ibin3]->Scale(1./g_RecoJet_GenTrack_hJetTrackSignalBackground[ibin][ibin3]->GetXaxis()->GetBinWidth(1)); 
        
        g_RecoJet_GenTrack_hJetTrackSignalBackground_eta[ibin][ibin3]->Scale(1./integral_reco_g);  
        //g_RecoJet_GenTrack_hJetTrackSignalBackground_eta[ibin][ibin3]->Write();
        g_RecoJet_GenTrack_hJetTrackSignalBackground_eta[ibin][ibin3]->SetLineColor(kRed);  

        // gen jet gen tracks

        sprintf(saythis,"GenJet_GenTrack_hJetTrackSignalBackground_cent%d_trkpt%d",ibin,ibin3);
        GenJet_GenTrack_hJetTrackSignalBackground[ibin][ibin3] = (TH2D*)file_gg->Get((TString)("GenJet_GenTrack_hJetTrackSignalBackgroundCent"+cent[ibin]+"_Cent"+cent[ibin+1]+"_Pt100_Pt300_TrkPt"+trkpt[ibin3]+"_TrkPt"+trkpt[ibin3+1]))->Clone(saythis);

        sprintf(saythis,"GenJet_GenTrack_hJetTrackSignalBackground_eta_cent%d_trkpt%d",ibin,ibin3);
        GenJet_GenTrack_hJetTrackSignalBackground_eta[ibin][ibin3] = GenJet_GenTrack_hJetTrackSignalBackground[ibin][ibin3] -> ProjectionX(); 
        //GenJet_GenTrack_hJetTrackSignalBackground_eta[ibin][ibin3]->Sumw2();
        GenJet_GenTrack_hJetTrackSignalBackground_eta[ibin][ibin3]->Scale(1./GenJet_GenTrack_hJetTrackSignalBackground[ibin][ibin3]->GetXaxis()->GetBinWidth(1));
        
        GenJet_GenTrack_hJetTrackSignalBackground_eta[ibin][ibin3]->Scale(1./integral_gen); 
        //GenJet_GenTrack_hJetTrackSignalBackground_eta[ibin][ibin3]->Write();
        GenJet_GenTrack_hJetTrackSignalBackground_eta[ibin][ibin3]->SetLineColor(kBlack);
        GenJet_GenTrack_hJetTrackSignalBackground_eta[ibin][ibin3]->GetXaxis()->SetTitle("dEta");
        GenJet_GenTrack_hJetTrackSignalBackground_eta[ibin][ibin3]->GetXaxis()->SetTitleSize(0.05);

        sprintf(saythis,"q_GenJet_GenTrack_hJetTrackSignalBackground_cent%d_trkpt%d",ibin,ibin3);
        q_GenJet_GenTrack_hJetTrackSignalBackground[ibin][ibin3] = (TH2D*)file_q_gg->Get((TString)("GenJet_GenTrack_hJetTrackSignalBackgroundCent"+cent[ibin]+"_Cent"+cent[ibin+1]+"_Pt100_Pt300_TrkPt"+trkpt[ibin3]+"_TrkPt"+trkpt[ibin3+1]))->Clone(saythis);

        sprintf(saythis,"q_GenJet_GenTrack_hJetTrackSignalBackground_eta_cent%d_trkpt%d",ibin,ibin3);
        q_GenJet_GenTrack_hJetTrackSignalBackground_eta[ibin][ibin3] = q_GenJet_GenTrack_hJetTrackSignalBackground[ibin][ibin3] -> ProjectionX(); 
        //q_GenJet_GenTrack_hJetTrackSignalBackground_eta[ibin][ibin3]->Sumw2();
        q_GenJet_GenTrack_hJetTrackSignalBackground_eta[ibin][ibin3]->Scale(1./q_GenJet_GenTrack_hJetTrackSignalBackground[ibin][ibin3]->GetXaxis()->GetBinWidth(1)); 
        
        q_GenJet_GenTrack_hJetTrackSignalBackground_eta[ibin][ibin3]->Scale(1./integral_gen_q); 
        //q_GenJet_GenTrack_hJetTrackSignalBackground_eta[ibin][ibin3]->Write(); 
        q_GenJet_GenTrack_hJetTrackSignalBackground_eta[ibin][ibin3]->SetLineColor(kBlue);
 
        sprintf(saythis,"g_GenJet_GenTrack_hJetTrackSignalBackground_cent%d_trkpt%d",ibin,ibin3);
        g_GenJet_GenTrack_hJetTrackSignalBackground[ibin][ibin3] = (TH2D*)file_g_gg->Get((TString)("GenJet_GenTrack_hJetTrackSignalBackgroundCent"+cent[ibin]+"_Cent"+cent[ibin+1]+"_Pt100_Pt300_TrkPt"+trkpt[ibin3]+"_TrkPt"+trkpt[ibin3+1]))->Clone(saythis);

        sprintf(saythis,"g_GenJet_GenTrack_hJetTrackSignalBackground_eta_cent%d_trkpt%d",ibin,ibin3); 
        g_GenJet_GenTrack_hJetTrackSignalBackground_eta[ibin][ibin3] = g_GenJet_GenTrack_hJetTrackSignalBackground[ibin][ibin3] -> ProjectionX(); 
        //g_GenJet_GenTrack_hJetTrackSignalBackground_eta[ibin][ibin3]->Sumw2();
        g_GenJet_GenTrack_hJetTrackSignalBackground_eta[ibin][ibin3]->Scale(1./g_GenJet_GenTrack_hJetTrackSignalBackground[ibin][ibin3]->GetXaxis()->GetBinWidth(1));
        
        g_GenJet_GenTrack_hJetTrackSignalBackground_eta[ibin][ibin3]->Scale(1./integral_gen_g);  
        //g_GenJet_GenTrack_hJetTrackSignalBackground_eta[ibin][ibin3]->Write();
        g_GenJet_GenTrack_hJetTrackSignalBackground_eta[ibin][ibin3]->SetLineColor(kRed);

        //recogen - gengen

        sprintf(saythis,"recogen_gengen_cent%d_trkpt%d",ibin,ibin3);
        recogen_gengen[ibin][ibin3] = (TH1D*)RecoJet_GenTrack_hJetTrackSignalBackground_eta[ibin][ibin3]->Clone(saythis);
        recogen_gengen[ibin][ibin3] -> Add(GenJet_GenTrack_hJetTrackSignalBackground_eta[ibin][ibin3],-1); 
        recogen_gengen[ibin][ibin3] -> SetLineColor(kBlack);
        recogen_gengen[ibin][ibin3] -> Rebin(5);
        recogen_gengen[ibin][ibin3] -> GetYaxis()->SetRangeUser(-10,5);
        recogen_gengen[ibin][ibin3] -> GetYaxis()->SetLabelSize(0.07);
        recogen_gengen[ibin][ibin3] -> GetXaxis()->SetTitle("dEta");
        recogen_gengen[ibin][ibin3] -> GetXaxis()->SetTitleSize(0.05);

        sprintf(saythis,"recogen_gengen_q_cent%d_trkpt%d",ibin,ibin3);
        recogen_gengen_q[ibin][ibin3] = (TH1D*)q_RecoJet_GenTrack_hJetTrackSignalBackground_eta[ibin][ibin3]->Clone(saythis);
        recogen_gengen_q[ibin][ibin3] -> Add(GenJet_GenTrack_hJetTrackSignalBackground_eta[ibin][ibin3],-1);
        recogen_gengen_q[ibin][ibin3] -> Rebin(5);
        recogen_gengen_q[ibin][ibin3] -> SetLineColor(kBlue);
 
        sprintf(saythis,"recogen_gengen_g_cent%d_trkpt%d",ibin,ibin3);
        recogen_gengen_g[ibin][ibin3] = (TH1D*)g_RecoJet_GenTrack_hJetTrackSignalBackground_eta[ibin][ibin3]->Clone(saythis);
        recogen_gengen_g[ibin][ibin3] -> Add(GenJet_GenTrack_hJetTrackSignalBackground_eta[ibin][ibin3],-1);  
        recogen_gengen_g[ibin][ibin3] -> Rebin(5);
        recogen_gengen_g[ibin][ibin3] -> SetLineColor(kRed);

        integral_all[ibin][ibin3] = recogen_gengen[ibin][ibin3]->Integral();
        //cout<<"integral of all "<<ibin<<"  "<<ibin3<<"  "<<integral_all[ibin][ibin3]<<endl;
        Double_t binwidth = int_bin_bounds[ibin3+1] - int_bin_bounds[ibin3];
        hintegral_eta[ibin]->SetBinContent(ibin3+1,integral_all[ibin][ibin3]/binwidth);
        hintegral_eta[ibin]->SetMarkerStyle(20);
        hintegral_eta[ibin]->SetMarkerColor(kBlack);
        //cout<<hintegral_eta[ibin]->GetXaxis()->GetBinWidth(ibin3+1)<<endl;  
        
        integral_q[ibin][ibin3] = recogen_gengen_q[ibin][ibin3]->Integral();
        //cout<<"integral of q "<<ibin<<"  "<<ibin3<<"  "<<integral_q[ibin][ibin3]<<endl;
        hintegral_eta_q[ibin]->SetBinContent(ibin3+1,integral_q[ibin][ibin3]/binwidth);  
        hintegral_eta_q[ibin]->SetMarkerStyle(20); 
        hintegral_eta_q[ibin]->SetMarkerColor(kBlue);

        integral_g[ibin][ibin3] = recogen_gengen_g[ibin][ibin3]->Integral();
        //cout<<"integral of g "<<ibin<<"  "<<ibin3<<"  "<<integral_g[ibin][ibin3]<<endl; 
        hintegral_eta_g[ibin]->SetBinContent(ibin3+1,integral_g[ibin][ibin3]/binwidth);
        //hintegral_eta_g[ibin]->GetYaxis()->SetRangeUser();
        hintegral_eta_g[ibin]->SetMarkerStyle(20);
        hintegral_eta_g[ibin]->SetMarkerColor(kRed); 

      } 

      hintegral_eta_g[ibin]->GetYaxis()->SetRangeUser(hintegral_eta_q[ibin]->GetBinContent(hintegral_eta_q[ibin]->GetMinimumBin())*1.2,hintegral_eta_g[ibin]->GetBinContent(hintegral_eta_g[ibin]->GetMaximumBin())*5);
    }     
    //cout<<"integral is"<<RecoJet_GenTrack_hJetTrackSignalBackground_eta[2][3]->Integral()<<endl;
    //q_RecoJet_GenTrack_hJetTrackSignalBackground_eta[1][5]->Draw();

    //f_yield.Close();

/// subtracted yield integral canvas

    TCanvas *c1 = new TCanvas("c1","",1200,300);
    c1->Divide(4,1);

    c1->cd(1)->SetPad(0.01,0.01,0.25,0.99);

    hintegral_eta_g[3]->SetLineColor(kRed);
    hintegral_eta_g[3]->Draw("e1 P same");

    hintegral_eta_q[3]->SetLineColor(kBlue);
    hintegral_eta_q[3]->Draw("e1 P same");
 
    hintegral_eta[3]->SetLineColor(kBlack);
    hintegral_eta[3]->SetLineWidth(2);
    hintegral_eta[3]->Draw("e1 P same");
    
    c1->cd(2)->SetPad(0.25,0.01,0.5,0.99);

    hintegral_eta_g[2]->SetLineColor(kRed);
    hintegral_eta_g[2]->Draw("e1 P same");

    hintegral_eta_q[2]->SetLineColor(kBlue);
    hintegral_eta_q[2]->Draw("e1 P same");
    
    hintegral_eta[2]->SetLineColor(kBlack);
    hintegral_eta[2]->SetLineWidth(2);
    hintegral_eta[2]->Draw("e1 P same");

    c1->cd(3)->SetPad(0.5,0.01,0.75,0.99);
    
    hintegral_eta_g[1]->SetLineColor(kRed);
    hintegral_eta_g[1]->Draw("e1 P same");

    hintegral_eta_q[1]->SetLineColor(kBlue);
    hintegral_eta_q[1]->Draw("e1 P same");
    
    hintegral_eta[1]->SetLineColor(kBlack);
    hintegral_eta[1]->SetLineWidth(2);
    hintegral_eta[1]->Draw("e1 P same");

    c1->cd(4)->SetPad(0.75,0.01,0.99,0.99);

    hintegral_eta_g[0]->SetLineColor(kRed);
    hintegral_eta_g[0]->Draw("e1 P same");

    hintegral_eta_q[0]->SetLineColor(kBlue);
    hintegral_eta_q[0]->Draw("e1 P same");
    
    hintegral_eta[0]->SetLineColor(kBlack);
    hintegral_eta[0]->SetLineWidth(2);
    hintegral_eta[0]->Draw("e1 P same");

//////subtracted deta yields

    TCanvas *c2 = new TCanvas("c2","",1200,1200);
    c2->Divide(4,4);
    gStyle->SetOptStat(0); 

    c2->cd(1)->SetPad(0.01,0.75,0.25,0.99);
    recogen_gengen[3][4]->Draw("same");
    recogen_gengen_q[3][4]->Draw("same");
    recogen_gengen_g[3][4]->Draw("same");
    c2->cd(2)->SetPad(0.25,0.75,0.5,0.99);
    recogen_gengen[2][4]->Draw("same");
    recogen_gengen_q[2][4]->Draw("same");
    recogen_gengen_g[2][4]->Draw("same");
    c2->cd(3)->SetPad(0.5,0.75,0.75,0.99);
    recogen_gengen[1][4]->Draw("same");
    recogen_gengen_q[1][4]->Draw("same");
    recogen_gengen_g[1][4]->Draw("same");
    c2->cd(4)->SetPad(0.75,0.75,0.99,0.99);
    recogen_gengen[0][4]->Draw("same");
    recogen_gengen_q[0][4]->Draw("same");
    recogen_gengen_g[0][4]->Draw("same");
    c2->cd(5)->SetPad(0.01,0.5,0.25,0.75);
    recogen_gengen[3][5]->Draw("same");
    recogen_gengen_q[3][5]->Draw("same");
    recogen_gengen_g[3][5]->Draw("same");
    c2->cd(6)->SetPad(0.25,0.5,0.5,0.75);
    recogen_gengen[2][5]->Draw("same");
    recogen_gengen_q[2][5]->Draw("same");
    recogen_gengen_g[2][5]->Draw("same");
    c2->cd(7)->SetPad(0.5,0.5,0.75,0.75);
    recogen_gengen[1][5]->Draw("same");
    recogen_gengen_q[1][5]->Draw("same");
    recogen_gengen_g[1][5]->Draw("same");
    c2->cd(8)->SetPad(0.75,0.5,0.99,0.75);
    recogen_gengen[0][5]->Draw("same");
    recogen_gengen_q[0][5]->Draw("same");
    recogen_gengen_g[0][5]->Draw("same");
    c2->cd(9)->SetPad(0.01,0.25,0.25,0.5);
    recogen_gengen[3][6]->Draw("same");
    recogen_gengen_q[3][6]->Draw("same");
    recogen_gengen_g[3][6]->Draw("same");
    c2->cd(10)->SetPad(0.25,0.25,0.5,0.5);
    recogen_gengen[2][6]->Draw("same");
    recogen_gengen_q[2][6]->Draw("same");
    recogen_gengen_g[2][6]->Draw("same");
    c2->cd(11)->SetPad(0.5,0.25,0.75,0.5);
    recogen_gengen[1][6]->Draw("same");
    recogen_gengen_q[1][6]->Draw("same");
    recogen_gengen_g[1][6]->Draw("same");
    c2->cd(12)->SetPad(0.75,0.25,0.99,0.5);
    recogen_gengen[0][6]->Draw("same");
    recogen_gengen_q[0][6]->Draw("same");
    recogen_gengen_g[0][6]->Draw("same");
    c2->cd(13)->SetPad(0.01,0.01,0.25,0.25);
    recogen_gengen[3][7]->Draw("same");
    recogen_gengen_q[3][7]->Draw("same");
    recogen_gengen_g[3][7]->Draw("same");
    c2->cd(14)->SetPad(0.25,0.01,0.5,0.25);
    recogen_gengen[2][7]->Draw("same");
    recogen_gengen_q[2][7]->Draw("same");
    recogen_gengen_g[2][7]->Draw("same");   
    c2->cd(15)->SetPad(0.5,0.01,0.75,0.25);
    recogen_gengen[1][7]->Draw("same");
    recogen_gengen_q[1][7]->Draw("same");
    recogen_gengen_g[1][7]->Draw("same");
    c2->cd(16)->SetPad(0.75,0.01,0.99,0.25);
    recogen_gengen[0][7]->Draw("same");
    recogen_gengen_q[0][7]->Draw("same");
    recogen_gengen_g[0][7]->Draw("same");

    TCanvas *c3 = new TCanvas("c3","",1200,1200);
    c3->Divide(4,4);
    gStyle->SetOptStat(0); 

    c3->cd(1)->SetPad(0.01,0.75,0.25,0.99);
    recogen_gengen[3][0]->Draw("same");
    recogen_gengen_q[3][0]->Draw("same");
    recogen_gengen_g[3][0]->Draw("same");
    c3->cd(2)->SetPad(0.25,0.75,0.5,0.99);
    recogen_gengen[2][0]->Draw("same");
    recogen_gengen_q[2][0]->Draw("same");
    recogen_gengen_g[2][0]->Draw("same");
    c3->cd(3)->SetPad(0.5,0.75,0.75,0.99);
    recogen_gengen[1][0]->Draw("same");
    recogen_gengen_q[1][0]->Draw("same");
    recogen_gengen_g[1][0]->Draw("same");
    c3->cd(4)->SetPad(0.75,0.75,0.99,0.99);
    recogen_gengen[0][0]->Draw("same");
    recogen_gengen_q[0][0]->Draw("same");
    recogen_gengen_g[0][0]->Draw("same");
    c3->cd(5)->SetPad(0.01,0.5,0.25,0.75);
    recogen_gengen[3][1]->Draw("same");
    recogen_gengen_q[3][1]->Draw("same");
    recogen_gengen_g[3][1]->Draw("same");
    c3->cd(6)->SetPad(0.25,0.5,0.5,0.75);
    recogen_gengen[2][1]->Draw("same");
    recogen_gengen_q[2][1]->Draw("same");
    recogen_gengen_g[2][1]->Draw("same");
    c3->cd(7)->SetPad(0.5,0.5,0.75,0.75);
    recogen_gengen[1][1]->Draw("same");
    recogen_gengen_q[1][1]->Draw("same");
    recogen_gengen_g[1][1]->Draw("same");
    c3->cd(8)->SetPad(0.75,0.5,0.99,0.75);
    recogen_gengen[0][1]->Draw("same");
    recogen_gengen_q[0][1]->Draw("same");
    recogen_gengen_g[0][1]->Draw("same");
    c3->cd(9)->SetPad(0.01,0.25,0.25,0.5);
    recogen_gengen[3][2]->Draw("same");
    recogen_gengen_q[3][2]->Draw("same");
    recogen_gengen_g[3][2]->Draw("same");
    c3->cd(10)->SetPad(0.25,0.25,0.5,0.5);
    recogen_gengen[2][2]->Draw("same");
    recogen_gengen_q[2][2]->Draw("same");
    recogen_gengen_g[2][2]->Draw("same");
    c3->cd(11)->SetPad(0.5,0.25,0.75,0.5);
    recogen_gengen[1][2]->Draw("same");
    recogen_gengen_q[1][2]->Draw("same");
    recogen_gengen_g[1][2]->Draw("same");
    c3->cd(12)->SetPad(0.75,0.25,0.99,0.5);
    recogen_gengen[0][2]->Draw("same");
    recogen_gengen_q[0][2]->Draw("same");
    recogen_gengen_g[0][2]->Draw("same");
    c3->cd(13)->SetPad(0.01,0.01,0.25,0.25);
    recogen_gengen[3][3]->Draw("same");
    recogen_gengen_q[3][3]->Draw("same");
    recogen_gengen_g[3][3]->Draw("same");
    c3->cd(14)->SetPad(0.25,0.01,0.5,0.25);
    recogen_gengen[2][3]->Draw("same");
    recogen_gengen_q[2][3]->Draw("same");
    recogen_gengen_g[2][3]->Draw("same");   
    c3->cd(15)->SetPad(0.5,0.01,0.75,0.25);
    recogen_gengen[1][3]->Draw("same");
    recogen_gengen_q[1][3]->Draw("same");
    recogen_gengen_g[1][3]->Draw("same");
    c3->cd(16)->SetPad(0.75,0.01,0.99,0.25);
    recogen_gengen[0][3]->Draw("same");
    recogen_gengen_q[0][3]->Draw("same");
    recogen_gengen_g[0][3]->Draw("same");

/*
    c2->cd(1)->SetPad(0.01,0.75,0.25,0.99);
    //g_RecoJet_GenTrack_hJetTrackSignalBackground_eta[3][0]->Draw("e1 same");
    g_RecoJet_GenTrack_hJetTrackSignalBackground_eta[3][0]->Draw("e1 same");
    q_RecoJet_GenTrack_hJetTrackSignalBackground_eta[3][0]->Draw("e1 same");
    RecoJet_GenTrack_hJetTrackSignalBackground_eta[3][0]->Draw("e1 same");
    
    c2->cd(2)->SetPad(0.25,0.75,0.5,0.99);
    //g_RecoJet_GenTrack_hJetTrackSignalBackground_eta[2][0]->Draw("e1 same");
    g_RecoJet_GenTrack_hJetTrackSignalBackground_eta[2][0]->Draw("e1 same");
    q_RecoJet_GenTrack_hJetTrackSignalBackground_eta[2][0]->Draw("e1 same");
    RecoJet_GenTrack_hJetTrackSignalBackground_eta[2][0]->Draw("e1 same");
    
    c2->cd(3)->SetPad(0.5,0.75,0.75,0.99);
    //g_RecoJet_GenTrack_hJetTrackSignalBackground_eta[1][0]->Draw("e1 same");
    g_RecoJet_GenTrack_hJetTrackSignalBackground_eta[1][0]->Draw("e1 same");
    q_RecoJet_GenTrack_hJetTrackSignalBackground_eta[1][0]->Draw("e1 same");
    RecoJet_GenTrack_hJetTrackSignalBackground_eta[1][0]->Draw("e1 same");    
    
    c2->cd(4)->SetPad(0.75,0.75,0.99,0.99);
    //g_RecoJet_GenTrack_hJetTrackSignalBackground_eta[0][0]->Draw("e1 same");
    g_RecoJet_GenTrack_hJetTrackSignalBackground_eta[0][0]->Draw("e1 same");
    q_RecoJet_GenTrack_hJetTrackSignalBackground_eta[0][0]->Draw("e1 same");
    RecoJet_GenTrack_hJetTrackSignalBackground_eta[0][0]->Draw("e1 same");
    c2->cd(5)->SetPad(0.01,0.5,0.25,0.75);
    //g_RecoJet_GenTrack_hJetTrackSignalBackground_eta[3][1]->Draw("e1 same");
    g_RecoJet_GenTrack_hJetTrackSignalBackground_eta[3][1]->Draw("e1 same");
    q_RecoJet_GenTrack_hJetTrackSignalBackground_eta[3][1]->Draw("e1 same");
    RecoJet_GenTrack_hJetTrackSignalBackground_eta[3][1]->Draw("e1 same");
    
    c2->cd(6)->SetPad(0.25,0.5,0.5,0.75);
    //g_RecoJet_GenTrack_hJetTrackSignalBackground_eta[2][1]->Draw("e1 same");
    g_RecoJet_GenTrack_hJetTrackSignalBackground_eta[2][1]->Draw("e1 same");
    q_RecoJet_GenTrack_hJetTrackSignalBackground_eta[2][1]->Draw("e1 same");
    RecoJet_GenTrack_hJetTrackSignalBackground_eta[2][1]->Draw("e1 same");
    
    c2->cd(7)->SetPad(0.5,0.5,0.75,0.75);
    //g_RecoJet_GenTrack_hJetTrackSignalBackground_eta[1][1]->Draw("e1 same");
    g_RecoJet_GenTrack_hJetTrackSignalBackground_eta[1][1]->Draw("e1 same");
    q_RecoJet_GenTrack_hJetTrackSignalBackground_eta[1][1]->Draw("e1 same");
    RecoJet_GenTrack_hJetTrackSignalBackground_eta[1][1]->Draw("e1 same");
    
    c2->cd(8)->SetPad(0.75,0.5,0.99,0.75);
    //g_RecoJet_GenTrack_hJetTrackSignalBackground_eta[0][1]->Draw("e1 same");
    g_RecoJet_GenTrack_hJetTrackSignalBackground_eta[0][1]->Draw("e1 same");
    q_RecoJet_GenTrack_hJetTrackSignalBackground_eta[0][1]->Draw("e1 same");
    RecoJet_GenTrack_hJetTrackSignalBackground_eta[0][1]->Draw("e1 same");
    
    c2->cd(9)->SetPad(0.01,0.25,0.25,0.5);
    //g_RecoJet_GenTrack_hJetTrackSignalBackground_eta[3][2]->Draw("e1 same");
    g_RecoJet_GenTrack_hJetTrackSignalBackground_eta[3][2]->Draw("e1 same");
    q_RecoJet_GenTrack_hJetTrackSignalBackground_eta[3][2]->Draw("e1 same");
    RecoJet_GenTrack_hJetTrackSignalBackground_eta[3][2]->Draw("e1 same");
    
    c2->cd(10)->SetPad(0.25,0.25,0.5,0.5);
    //g_RecoJet_GenTrack_hJetTrackSignalBackground_eta[2][2]->Draw("e1 same");
    g_RecoJet_GenTrack_hJetTrackSignalBackground_eta[2][2]->Draw("e1 same");
    q_RecoJet_GenTrack_hJetTrackSignalBackground_eta[2][2]->Draw("e1 same");
    RecoJet_GenTrack_hJetTrackSignalBackground_eta[2][2]->Draw("e1 same");
    
    c2->cd(11)->SetPad(0.5,0.25,0.75,0.5);
    //g_RecoJet_GenTrack_hJetTrackSignalBackground_eta[1][2]->Draw("e1 same");
    g_RecoJet_GenTrack_hJetTrackSignalBackground_eta[1][2]->Draw("e1 same");
    q_RecoJet_GenTrack_hJetTrackSignalBackground_eta[1][2]->Draw("e1 same");
    RecoJet_GenTrack_hJetTrackSignalBackground_eta[1][2]->Draw("e1 same");
    
    c2->cd(12)->SetPad(0.75,0.25,0.99,0.5);
    //g_RecoJet_GenTrack_hJetTrackSignalBackground_eta[0][2]->Draw("e1 same");
    g_RecoJet_GenTrack_hJetTrackSignalBackground_eta[0][2]->Draw("e1 same");
    q_RecoJet_GenTrack_hJetTrackSignalBackground_eta[0][2]->Draw("e1 same");
    RecoJet_GenTrack_hJetTrackSignalBackground_eta[0][2]->Draw("e1 same");
    
    c2->cd(13)->SetPad(0.01,0.01,0.25,0.25);
    //g_RecoJet_GenTrack_hJetTrackSignalBackground_eta[3][3]->Draw("e1 same");
    g_RecoJet_GenTrack_hJetTrackSignalBackground_eta[3][3]->Draw("e1 same");
    q_RecoJet_GenTrack_hJetTrackSignalBackground_eta[3][3]->Draw("e1 same");
    RecoJet_GenTrack_hJetTrackSignalBackground_eta[3][3]->Draw("e1 same");
    
    c2->cd(14)->SetPad(0.25,0.01,0.5,0.25);
    //g_RecoJet_GenTrack_hJetTrackSignalBackground_eta[2][3]->Draw("e1 same");
    g_RecoJet_GenTrack_hJetTrackSignalBackground_eta[2][3]->Draw("e1 same");    
    q_RecoJet_GenTrack_hJetTrackSignalBackground_eta[2][3]->Draw("e1 same");    
    RecoJet_GenTrack_hJetTrackSignalBackground_eta[2][3]->Draw("e1 same");
        
    c2->cd(15)->SetPad(0.5,0.01,0.75,0.25);
    //g_RecoJet_GenTrack_hJetTrackSignalBackground_eta[1][3]->Draw("e1 same");
    g_RecoJet_GenTrack_hJetTrackSignalBackground_eta[1][3]->Draw("e1 same");
    q_RecoJet_GenTrack_hJetTrackSignalBackground_eta[1][3]->Draw("e1 same");
    RecoJet_GenTrack_hJetTrackSignalBackground_eta[1][3]->Draw("e1 same");
    
    c2->cd(16)->SetPad(0.75,0.01,0.99,0.25);
    //g_RecoJet_GenTrack_hJetTrackSignalBackground_eta[0][3]->Draw("e1 same");
    g_RecoJet_GenTrack_hJetTrackSignalBackground_eta[0][3]->Draw("e1 same");
    q_RecoJet_GenTrack_hJetTrackSignalBackground_eta[0][3]->Draw("e1 same");    
    RecoJet_GenTrack_hJetTrackSignalBackground_eta[0][3]->Draw("e1 same");
*/     
/*
    TCanvas *c3 = new TCanvas("c3","",1600,1600);
    c3->Divide(4,4);
    c3->cd(1)->SetPad(0.01,0.75,0.25,0.99);
    c3->cd(2)->SetPad(0.25,0.75,0.5,0.99);
    c3->cd(3)->SetPad(0.5,0.75,0.75,0.99);
    c3->cd(4)->SetPad(0.75,0.75,0.99,0.99);
    c3->cd(5)->SetPad(0.01,0.5,0.25,0.75);
    c3->cd(6)->SetPad(0.25,0.5,0.5,0.75);
    c3->cd(7)->SetPad(0.5,0.5,0.75,0.75);
    c3->cd(8)->SetPad(0.75,0.5,0.99,0.75);
    c3->cd(9)->SetPad(0.01,0.25,0.25,0.5);
    c3->cd(10)->SetPad(0.25,0.25,0.5,0.5);
    c3->cd(11)->SetPad(0.5,0.25,0.75,0.5);
    c3->cd(12)->SetPad(0.75,0.25,0.99,0.5);
    c3->cd(13)->SetPad(0.01,0.01,0.25,0.25);
    c3->cd(14)->SetPad(0.25,0.01,0.5,0.25);
    c3->cd(15)->SetPad(0.5,0.01,0.75,0.25);
    c3->cd(16)->SetPad(0.75,0.01,0.99,0.25);
    TCanvas *c4 = new TCanvas("c4","",1600,1600);
    c4->Divide(4,4);
    c4->cd(1)->SetPad(0.01,0.75,0.25,0.99);
    c4->cd(2)->SetPad(0.25,0.75,0.5,0.99);
    c4->cd(3)->SetPad(0.5,0.75,0.75,0.99);
    c4->cd(4)->SetPad(0.75,0.75,0.99,0.99);
    c4->cd(5)->SetPad(0.01,0.5,0.25,0.75);
    c4->cd(6)->SetPad(0.25,0.5,0.5,0.75);
    c4->cd(7)->SetPad(0.5,0.5,0.75,0.75);
    c4->cd(8)->SetPad(0.75,0.5,0.99,0.75);
    c4->cd(9)->SetPad(0.01,0.25,0.25,0.5);
    c4->cd(10)->SetPad(0.25,0.25,0.5,0.5);
    c4->cd(11)->SetPad(0.5,0.25,0.75,0.5);
    c4->cd(12)->SetPad(0.75,0.25,0.99,0.5);
    c4->cd(13)->SetPad(0.01,0.01,0.25,0.25);
    c4->cd(14)->SetPad(0.25,0.01,0.5,0.25);
    c4->cd(15)->SetPad(0.5,0.01,0.75,0.25);
    c4->cd(16)->SetPad(0.75,0.01,0.99,0.25);
*/
}