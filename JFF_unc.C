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

using namespace std;

#define	centbin 4
#define trkptbin 8

char saythis[500];

TString cent[centbin+1] = {"0","10","30","50","100"};
TString trkpt[trkptbin+1] = {"0p7","1","2","3","4","8","12","16","20"};

Double_t int_bin_bounds[9] = {0.7,1.,2.,3.,4.,8.,12.,16.,20.};
int integral_bins = 9;

void JFF_unc(){

	TFile *file_q_rg = TFile::Open("/home/dhanush/Documents/Kurt_histos/PbPb_5TeV_MC_PythiaHydjet_MixHistos_RecoGenReduced_QuarkJet_sube0_Merged_noMix_newJetTrackCorrections_fineBin.root");
	TFile *file_g_rg = TFile::Open("/home/dhanush/Documents/Kurt_histos/PbPb_5TeV_MC_PythiaHydjet_MixHistos_RecoGenReduced_GluonJet_sube0_Merged_noMix_newJetTrackCorrections_fineBin.root");
    TFile *file_rg = TFile::Open("/home/dhanush/Documents/Kurt_histos/PbPb_5TeV_MC_PythiaHydjet_MixHistos_RecoGenReduced_sube0_Merged_refpt_newJetTrackCorrections_fineBin.root"); 

    TFile *file_q_gg = TFile::Open("/home/dhanush/Documents/Kurt_histos/PbPb_5TeV_MC_PythiaHydjet_MixHistos_GenGenReduced_QuarkJet_sube0_Merged_noMix_newJetTrackCorrections_fineBin.root");
	TFile *file_g_gg = TFile::Open("/home/dhanush/Documents/Kurt_histos/PbPb_5TeV_MC_PythiaHydjet_MixHistos_GenGenReduced_GluonJet_sube0_Merged_noMix_newJetTrackCorrections_fineBin.root");
    TFile *file_gg = TFile::Open("/home/dhanush/Documents/Kurt_histos/PbPb_5TeV_MC_PythiaHydjet_MixHistos_GenGenReduced_sube0_Merged_refpt_newJetTrackCorrections_fineBin.root");

    //TFile f_yield("/home/dhanush/Documents/JetTrack2016/subtractedyield_integral.root", "RECREATE"); 

    TH2D *q_RecoJet_GenTrack_hJetTrackSignalBackground[centbin][trkptbin];
    TH2D *g_RecoJet_GenTrack_hJetTrackSignalBackground[centbin][trkptbin];
    TH2D *RecoJet_GenTrack_hJetTrackSignalBackground[centbin][trkptbin];

    TH2D *q_RecoJet_GenTrack_hJetTrackME[centbin][trkptbin];
    TH2D *g_RecoJet_GenTrack_hJetTrackME[centbin][trkptbin];
    TH2D *RecoJet_GenTrack_hJetTrackME[centbin][trkptbin];

    TH1D *q_RecoJet_GenTrack_hJetTrackSignalBackground_eta[centbin][trkptbin];
    TH1D *g_RecoJet_GenTrack_hJetTrackSignalBackground_eta[centbin][trkptbin];
    TH1D *RecoJet_GenTrack_hJetTrackSignalBackground_eta[centbin][trkptbin];

    TH1D *RecoJet_GenTrack_hJetTrackME_eta[centbin][trkptbin];
    TH1D *q_RecoJet_GenTrack_hJetTrackME_eta[centbin][trkptbin];
    TH1D *g_RecoJet_GenTrack_hJetTrackME_eta[centbin][trkptbin]; 

    TH2D *q_GenJet_GenTrack_hJetTrackSignalBackground[centbin][trkptbin];
    TH2D *g_GenJet_GenTrack_hJetTrackSignalBackground[centbin][trkptbin];
    TH2D *GenJet_GenTrack_hJetTrackSignalBackground[centbin][trkptbin];

    TH2D *GenJet_GenTrack_hJetTrackME[centbin][trkptbin];
    TH1D *GenJet_GenTrack_hJetTrackME_eta[centbin][trkptbin];

    TH1D *q_GenJet_GenTrack_hJetTrackSignalBackground_eta[centbin][trkptbin];
    TH1D *g_GenJet_GenTrack_hJetTrackSignalBackground_eta[centbin][trkptbin];
    TH1D *GenJet_GenTrack_hJetTrackSignalBackground_eta[centbin][trkptbin];
 
    TH1D *recogen_gengen[centbin][trkptbin];
    TH1D *recogen_gengen_q[centbin][trkptbin];
    TH1D *recogen_gengen_g[centbin][trkptbin];

    TH1D *hintegral_eta[centbin]; 
    TH1D *hintegral_eta_q[centbin];
    TH1D *hintegral_eta_g[centbin];

    Double_t integral_all[centbin][trkptbin];
    Double_t integral_q[centbin][trkptbin];
    Double_t integral_g[centbin][trkptbin];

    TH1F *reco_jets[centbin];
    TH1F *q_reco_jets[centbin];
    TH1F *g_reco_jets[centbin];
    TH1F *gen_jets[centbin];
    TH1F *q_gen_jets[centbin];
    TH1F *g_gen_jets[centbin];

    //TLegend *legend[4][10];

    for(int ibin=0;ibin<centbin;ibin++){

      Double_t integral_reco = 0;
      Double_t integral_reco_q = 0;
      Double_t integral_reco_g = 0;
      Double_t integral_gen = 0;
      
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
      integral_reco = reco_jets[ibin]->Integral();
      //cout<<integral_reco<<endl;

      sprintf(saythis,"q_reco_jets_%d",ibin);
      q_reco_jets[ibin] = (TH1F*)file_q_rg->Get((TString)("RecoJet_GenTrack_all_jets_corrpTCent"+cent[ibin]+"_Cent"+cent[ibin+1]+"_Pt100_Pt300"))->Clone(saythis);
      integral_reco_q = q_reco_jets[ibin]->Integral(); 

      sprintf(saythis,"g_reco_jets_%d",ibin);
      g_reco_jets[ibin] = (TH1F*)file_g_rg->Get((TString)("RecoJet_GenTrack_all_jets_corrpTCent"+cent[ibin]+"_Cent"+cent[ibin+1]+"_Pt100_Pt300"))->Clone(saythis);
      integral_reco_g = g_reco_jets[ibin]->Integral(); 

      sprintf(saythis,"gen_jets_%d",ibin);  
      gen_jets[ibin] = (TH1F*)file_gg->Get((TString)("GenJet_GenTrack_all_jets_corrpTCent"+cent[ibin]+"_Cent"+cent[ibin+1]+"_Pt100_Pt300"))->Clone(saythis);
      integral_gen = gen_jets[ibin]->Integral();   
      //cout<<integral_gen<<endl;

      for(int ibin3=0;ibin3<trkptbin;ibin3++){

        // gen jets gen tracks 

        sprintf(saythis,"GenJet_GenTrack_hJetTrackSignalBackground_cent%d_trkpt%d",ibin,ibin3);
        GenJet_GenTrack_hJetTrackSignalBackground[ibin][ibin3] = (TH2D*)file_gg->Get((TString)("GenJet_GenTrack_hJetTrackSignalBackgroundCent"+cent[ibin]+"_Cent"+cent[ibin+1]+"_Pt100_Pt300_TrkPt"+trkpt[ibin3]+"_TrkPt"+trkpt[ibin3+1]))->Clone(saythis); 
        //GenJet_GenTrack_hJetTrackSignalBackground[ibin][ibin3] -> Scale(1./GenJet_GenTrack_hJetTrackSignalBackground[ibin][ibin3]->GetXaxis()->GetBinWidth(2));
        GenJet_GenTrack_hJetTrackSignalBackground[ibin][ibin3] -> Scale(1./gen_jets[ibin]->Integral());
        
        int y1 = GenJet_GenTrack_hJetTrackSignalBackground[ibin][ibin3]->GetYaxis()->FindBin(TMath::Pi()/2);
        int y2 = GenJet_GenTrack_hJetTrackSignalBackground[ibin][ibin3]->GetYaxis()->FindBin(3*TMath::Pi()/2); 

        sprintf(saythis,"GenJet_GenTrack_hJetTrackME_eta_cent%d_trkpt%d",ibin,ibin3);
        GenJet_GenTrack_hJetTrackME_eta[ibin][ibin3] = GenJet_GenTrack_hJetTrackSignalBackground[ibin][ibin3] -> ProjectionX(saythis,y1,y2);
        GenJet_GenTrack_hJetTrackME_eta[ibin][ibin3]->Scale(1./GenJet_GenTrack_hJetTrackME_eta[ibin][ibin3]->GetMaximum());
        
        // gen jet gen track ME

        sprintf(saythis,"GenJet_GenTrack_hJetTrackME_cent%d_trkpt%d",ibin,ibin3);
        GenJet_GenTrack_hJetTrackME[ibin][ibin3] = new TH2D(saythis,"",500,-5.,5.,200,-TMath::Pi()/2,3*TMath::Pi()/2);

        for (int ixbin = 1; ixbin < 501; ixbin++){
            Double_t bincontent = 0.;
            Double_t binerror = 0.; 
            bincontent = GenJet_GenTrack_hJetTrackME_eta[ibin][ibin3]->GetBinContent(ixbin);     
            binerror = GenJet_GenTrack_hJetTrackME_eta[ibin][ibin3]->GetBinError(ixbin);       
            for (int iybin = 1; iybin < 201; iybin++){
              GenJet_GenTrack_hJetTrackME[ibin][ibin3]->SetBinContent(ixbin,iybin,bincontent);
              GenJet_GenTrack_hJetTrackME[ibin][ibin3]->SetBinError(ixbin,iybin,binerror); 
            }
        }

        GenJet_GenTrack_hJetTrackSignalBackground[ibin][ibin3]->Divide(GenJet_GenTrack_hJetTrackME[ibin][ibin3]);     

        // reco jet gen tracks

        sprintf(saythis,"RecoJet_GenTrack_hJetTrackSignalBackground_cent%d_trkpt%d",ibin,ibin3);
        RecoJet_GenTrack_hJetTrackSignalBackground[ibin][ibin3] = (TH2D*)file_rg->Get((TString)("RecoJet_GenTrack_hJetTrackSignalBackgroundCent"+cent[ibin]+"_Cent"+cent[ibin+1]+"_Pt100_Pt300_TrkPt"+trkpt[ibin3]+"_TrkPt"+trkpt[ibin3+1]))->Clone(saythis);
        //RecoJet_GenTrack_hJetTrackSignalBackground[ibin][ibin3] -> Scale(1./RecoJet_GenTrack_hJetTrackSignalBackground[ibin][ibin3]->GetXaxis()->GetBinWidth(2));
        RecoJet_GenTrack_hJetTrackSignalBackground[ibin][ibin3] -> Scale(1./reco_jets[ibin]->Integral());

        sprintf(saythis,"RecoJet_GenTrack_hJetTrackME_eta_cent%d_trkpt%d",ibin,ibin3);
        RecoJet_GenTrack_hJetTrackME_eta[ibin][ibin3] = RecoJet_GenTrack_hJetTrackSignalBackground[ibin][ibin3] -> ProjectionX(saythis,y1,y2);
        RecoJet_GenTrack_hJetTrackME_eta[ibin][ibin3]->Scale(1./RecoJet_GenTrack_hJetTrackME_eta[ibin][ibin3]->GetMaximum());

        // reco jet gen track ME

        sprintf(saythis,"RecoJet_GenTrack_hJetTrackME_cent%d_trkpt%d",ibin,ibin3);
        RecoJet_GenTrack_hJetTrackME[ibin][ibin3] = new TH2D(saythis,"",500,-5.,5.,200,-TMath::Pi()/2,3*TMath::Pi()/2);

        for (int ixbin = 1; ixbin < 501; ixbin++){
            Double_t bincontent = 0.;
            Double_t binerror = 0.;
            bincontent = RecoJet_GenTrack_hJetTrackME_eta[ibin][ibin3]->GetBinContent(ixbin);            
            binerror = RecoJet_GenTrack_hJetTrackME_eta[ibin][ibin3]->GetBinError(ixbin);
            for (int iybin = 1; iybin < 201; iybin++){
              RecoJet_GenTrack_hJetTrackME[ibin][ibin3]->SetBinContent(ixbin,iybin,bincontent);
              RecoJet_GenTrack_hJetTrackME[ibin][ibin3]->SetBinError(ixbin,iybin,binerror); 
            }
        }

        RecoJet_GenTrack_hJetTrackSignalBackground[ibin][ibin3] -> Divide(RecoJet_GenTrack_hJetTrackME[ibin][ibin3]);  
        RecoJet_GenTrack_hJetTrackSignalBackground[ibin][ibin3] -> Add(GenJet_GenTrack_hJetTrackSignalBackground[ibin][ibin3],-1); 

        int phi_lowBin = RecoJet_GenTrack_hJetTrackSignalBackground[ibin][ibin3]->GetYaxis()->FindBin(-1);
        int phi_highBin = RecoJet_GenTrack_hJetTrackSignalBackground[ibin][ibin3]->GetYaxis()->FindBin(1);

        // reco q jet gen track

        sprintf(saythis,"q_RecoJet_GenTrack_hJetTrackSignalBackground_cent%d_trkpt%d",ibin,ibin3);
        q_RecoJet_GenTrack_hJetTrackSignalBackground[ibin][ibin3] = (TH2D*)file_q_rg->Get((TString)("RecoJet_GenTrack_hJetTrackSignalBackgroundCent"+cent[ibin]+"_Cent"+cent[ibin+1]+"_Pt100_Pt300_TrkPt"+trkpt[ibin3]+"_TrkPt"+trkpt[ibin3+1]))->Clone(saythis);
        //q_RecoJet_GenTrack_hJetTrackSignalBackground[ibin][ibin3] -> Scale(1./q_RecoJet_GenTrack_hJetTrackSignalBackground[ibin][ibin3]->GetXaxis()->GetBinWidth(2));
        q_RecoJet_GenTrack_hJetTrackSignalBackground[ibin][ibin3] -> Scale(1./q_reco_jets[ibin]->Integral()); 
        
        sprintf(saythis,"q_RecoJet_GenTrack_hJetTrackME_eta_cent%d_trkpt%d",ibin,ibin3);
        q_RecoJet_GenTrack_hJetTrackME_eta[ibin][ibin3] = q_RecoJet_GenTrack_hJetTrackSignalBackground[ibin][ibin3] -> ProjectionX(saythis,y1,y2);
        q_RecoJet_GenTrack_hJetTrackME_eta[ibin][ibin3]->Scale(1./q_RecoJet_GenTrack_hJetTrackME_eta[ibin][ibin3]->GetMaximum()); 

        // reco q jet gen track ME

        sprintf(saythis,"q_RecoJet_GenTrack_hJetTrackME_cent%d_trkpt%d",ibin,ibin3);
        q_RecoJet_GenTrack_hJetTrackME[ibin][ibin3] = new TH2D(saythis,"",500,-5.,5.,200,-TMath::Pi()/2,3*TMath::Pi()/2);

        for (int ixbin = 1; ixbin < 501; ixbin++){
            Double_t bincontent = 0.;
            Double_t binerror = 0.;
            bincontent = q_RecoJet_GenTrack_hJetTrackME_eta[ibin][ibin3]->GetBinContent(ixbin);
            binerror = q_RecoJet_GenTrack_hJetTrackME_eta[ibin][ibin3]->GetBinError(ixbin);            
            for (int iybin = 1; iybin < 201; iybin++){
              q_RecoJet_GenTrack_hJetTrackME[ibin][ibin3]->SetBinContent(ixbin,iybin,bincontent); 
              q_RecoJet_GenTrack_hJetTrackME[ibin][ibin3]->SetBinError(ixbin,iybin,binerror);
            }
        }
 
        q_RecoJet_GenTrack_hJetTrackSignalBackground[ibin][ibin3] -> Divide(q_RecoJet_GenTrack_hJetTrackME[ibin][ibin3]); 
        q_RecoJet_GenTrack_hJetTrackSignalBackground[ibin][ibin3] -> Add(GenJet_GenTrack_hJetTrackSignalBackground[ibin][ibin3],-1);
        
        // reco g jet gen track 

        sprintf(saythis,"g_RecoJet_GenTrack_hJetTrackSignalBackground_cent%d_trkpt%d",ibin,ibin3);
        g_RecoJet_GenTrack_hJetTrackSignalBackground[ibin][ibin3] = (TH2D*)file_g_rg->Get((TString)("RecoJet_GenTrack_hJetTrackSignalBackgroundCent"+cent[ibin]+"_Cent"+cent[ibin+1]+"_Pt100_Pt300_TrkPt"+trkpt[ibin3]+"_TrkPt"+trkpt[ibin3+1]))->Clone(saythis);
        //g_RecoJet_GenTrack_hJetTrackSignalBackground[ibin][ibin3] -> Scale(1./g_RecoJet_GenTrack_hJetTrackSignalBackground[ibin][ibin3]->GetXaxis()->GetBinWidth(2));
        g_RecoJet_GenTrack_hJetTrackSignalBackground[ibin][ibin3] -> Scale(1./g_reco_jets[ibin]->Integral());
        
        sprintf(saythis,"g_RecoJet_GenTrack_hJetTrackME_eta_cent%d_trkpt%d",ibin,ibin3);
        g_RecoJet_GenTrack_hJetTrackME_eta[ibin][ibin3] = g_RecoJet_GenTrack_hJetTrackSignalBackground[ibin][ibin3] -> ProjectionX(saythis,y1,y2);
        g_RecoJet_GenTrack_hJetTrackME_eta[ibin][ibin3]->Scale(1./g_RecoJet_GenTrack_hJetTrackME_eta[ibin][ibin3]->GetMaximum());

        // reco g jet gen track ME

        sprintf(saythis,"g_RecoJet_GenTrack_hJetTrackME_cent%d_trkpt%d",ibin,ibin3);
        g_RecoJet_GenTrack_hJetTrackME[ibin][ibin3] = new TH2D(saythis,"",500,-5.,5.,200,-TMath::Pi()/2,3*TMath::Pi()/2);

        for (int ixbin = 1; ixbin < 501; ixbin++){
            Double_t bincontent = 0.;
            Double_t binerror = 0.;
            bincontent = g_RecoJet_GenTrack_hJetTrackME_eta[ibin][ibin3]->GetBinContent(ixbin);
            binerror = g_RecoJet_GenTrack_hJetTrackME_eta[ibin][ibin3]->GetBinError(ixbin);            
            for (int iybin = 1; iybin < 201; iybin++){
              g_RecoJet_GenTrack_hJetTrackME[ibin][ibin3]->SetBinContent(ixbin,iybin,bincontent);
              g_RecoJet_GenTrack_hJetTrackME[ibin][ibin3]->SetBinError(ixbin,iybin,binerror); 
            }
        }

        g_RecoJet_GenTrack_hJetTrackSignalBackground[ibin][ibin3] -> Divide(g_RecoJet_GenTrack_hJetTrackME[ibin][ibin3]);
        g_RecoJet_GenTrack_hJetTrackSignalBackground[ibin][ibin3] -> Add(GenJet_GenTrack_hJetTrackSignalBackground[ibin][ibin3],-1);

        //cout<<g_RecoJet_GenTrack_hJetTrackSignalBackground[ibin][ibin3]->GetYaxis()->GetBinWidth(1)<<endl;

        sprintf(saythis,"recogen_gengen_cent%d_trkpt%d",ibin,ibin3);
        recogen_gengen[ibin][ibin3] = RecoJet_GenTrack_hJetTrackSignalBackground[ibin][ibin3] -> ProjectionX(saythis,phi_lowBin, phi_highBin,"e");
        recogen_gengen[ibin][ibin3] -> Scale(1./RecoJet_GenTrack_hJetTrackSignalBackground[ibin][ibin3]->GetXaxis()->GetBinWidth(2));
        recogen_gengen[ibin][ibin3] -> Rebin(5);
        recogen_gengen[ibin][ibin3] -> Scale(1./5.);
        recogen_gengen[ibin][ibin3] -> SetLineColor(kBlack);

        sprintf(saythis,"recogen_gengen_q_cent%d_trkpt%d",ibin,ibin3);
        recogen_gengen_q[ibin][ibin3] = q_RecoJet_GenTrack_hJetTrackSignalBackground[ibin][ibin3] -> ProjectionX(saythis,phi_lowBin, phi_highBin,"e");  
        recogen_gengen_q[ibin][ibin3] -> Scale(1./q_RecoJet_GenTrack_hJetTrackSignalBackground[ibin][ibin3]->GetXaxis()->GetBinWidth(2));
        recogen_gengen_q[ibin][ibin3] -> Rebin(5);
        recogen_gengen_q[ibin][ibin3] -> Scale(1./5.);
        recogen_gengen_q[ibin][ibin3] -> SetLineColor(kBlue);
        recogen_gengen_q[ibin][ibin3] -> GetXaxis()->SetTitle("dEta");
        recogen_gengen_q[ibin][ibin3] -> GetXaxis()->SetTitleSize(0.05);  

        sprintf(saythis,"recogen_gengen_g_cent%d_trkpt%d",ibin,ibin3); 
        recogen_gengen_g[ibin][ibin3] = g_RecoJet_GenTrack_hJetTrackSignalBackground[ibin][ibin3] -> ProjectionX(saythis,phi_lowBin, phi_highBin,"e");
        recogen_gengen_g[ibin][ibin3] -> Scale(1./g_RecoJet_GenTrack_hJetTrackSignalBackground[ibin][ibin3]->GetXaxis()->GetBinWidth(2));
        recogen_gengen_g[ibin][ibin3] -> Rebin(5);
        recogen_gengen_g[ibin][ibin3] -> Scale(1./5.);
        recogen_gengen_g[ibin][ibin3] -> SetLineColor(kRed); 
        recogen_gengen_g[ibin][ibin3] -> GetYaxis()->SetLabelSize(0.07);
        recogen_gengen_g[ibin][ibin3] -> GetXaxis()->SetTitle("dEta");
        recogen_gengen_g[ibin][ibin3] -> GetXaxis()->SetTitleSize(0.05); 
        recogen_gengen_g[ibin][ibin3] -> GetYaxis() -> SetRangeUser((recogen_gengen_q[ibin][ibin3]->GetBinContent(50))*2,(recogen_gengen_g[ibin][ibin3]->GetBinContent(50))*3);  

        ///integral of deta yields
        
        int eta_lowBin = recogen_gengen[ibin][ibin3]->GetXaxis()->FindBin(-1);
        int eta_highBin = recogen_gengen[ibin][ibin3]->GetXaxis()->FindBin(1); 

        //cout<<eta_highBin<<"  "<<eta_lowBin<<endl;

        integral_all[ibin][ibin3] = recogen_gengen[ibin][ibin3]->Integral(eta_lowBin,eta_highBin,"e");
        //cout<<"integral of all "<<ibin<<"  "<<ibin3<<"  "<<integral_all[ibin][ibin3]<<endl;
        Double_t binwidth = int_bin_bounds[ibin3+1] - int_bin_bounds[ibin3];
        hintegral_eta[ibin]->SetBinContent(ibin3+1,integral_all[ibin][ibin3]/binwidth);
        hintegral_eta[ibin]->SetMarkerStyle(20);
        hintegral_eta[ibin]->SetMarkerColor(kBlack);
        //cout<<hintegral_eta[ibin]->GetXaxis()->GetBinWidth(ibin3+1)<<endl;  
        
        integral_q[ibin][ibin3] = recogen_gengen_q[ibin][ibin3]->Integral(eta_lowBin,eta_highBin,"e");
        //cout<<"integral of q "<<ibin<<"  "<<ibin3<<"  "<<integral_q[ibin][ibin3]<<endl;
        hintegral_eta_q[ibin]->SetBinContent(ibin3+1,integral_q[ibin][ibin3]/binwidth);  
        hintegral_eta_q[ibin]->SetMarkerStyle(20); 
        hintegral_eta_q[ibin]->SetMarkerColor(kBlue);

        integral_g[ibin][ibin3] = recogen_gengen_g[ibin][ibin3]->Integral(eta_lowBin,eta_highBin,"e");
        //cout<<"integral of g "<<ibin<<"  "<<ibin3<<"  "<<integral_g[ibin][ibin3]<<endl; 
        hintegral_eta_g[ibin]->SetBinContent(ibin3+1,integral_g[ibin][ibin3]/binwidth);
        //hintegral_eta_g[ibin]->GetYaxis()->SetRangeUser();
        hintegral_eta_g[ibin]->SetMarkerStyle(20);
        hintegral_eta_g[ibin]->SetMarkerColor(kRed); 
/*
        legend[ibin][ibin3]->AddEntry((TObject*)0,(TString)(cent[ibin]+"-"+cent[ibin+1]), "");
        legend[ibin][ibin3]->AddEntry((TObject*)0,(TString)(trkpt[ibin3]+"-"+trkpt[ibin3+1]), "");
*/
      } 

      hintegral_eta_g[ibin]->GetYaxis()->SetRangeUser(hintegral_eta_q[ibin]->GetBinContent(hintegral_eta_q[ibin]->GetMinimumBin())*2,hintegral_eta_g[ibin]->GetBinContent(hintegral_eta_g[ibin]->GetMaximumBin())*3);
      //hintegral_eta_g[ibin]->GetYaxis()->SetRangeUser(-6,4);
    }     
    //cout<<"integral is"<<RecoJet_GenTrack_hJetTrackSignalBackground_eta[2][3]->Integral()<<endl;
  
    TCanvas *c_temp = new TCanvas("c_temp","",400,400);
    c_temp->cd(1);
    GenJet_GenTrack_hJetTrackSignalBackground[1][1]->Draw("SURF2");

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
    recogen_gengen_g[3][4]->Draw("same");
    recogen_gengen_q[3][4]->Draw("same");

    TLegend *leg2 = new TLegend(0.2,0.5,0.5,0.9);
    leg2->AddEntry((TObject*)0, "RecoGen - GenGen", "");
    leg2->Draw("same"); 

    recogen_gengen[3][4]->Draw("same");
    c2->cd(2)->SetPad(0.25,0.75,0.5,0.99);
    recogen_gengen_g[2][4]->Draw("same");
    recogen_gengen_q[2][4]->Draw("same");
    recogen_gengen[2][4]->Draw("same");
    
    TLegend *leg1 = new TLegend(0.2,0.5,0.5,0.9);
    leg1->AddEntry(recogen_gengen[0][5],"all jets","lp");
    leg1->AddEntry(recogen_gengen_q[0][5],"q jets","lp");
    leg1->AddEntry(recogen_gengen_g[0][5],"g jets","lp");
    leg1->Draw("same");

    c2->cd(3)->SetPad(0.5,0.75,0.75,0.99);
    recogen_gengen_g[1][4]->Draw("same");
    recogen_gengen_q[1][4]->Draw("same");
    recogen_gengen[1][4]->Draw("same");
    c2->cd(4)->SetPad(0.75,0.75,0.99,0.99);
    recogen_gengen_g[0][4]->Draw("same");
    recogen_gengen_q[0][4]->Draw("same");
    recogen_gengen[0][4]->Draw("same");
    c2->cd(5)->SetPad(0.01,0.5,0.25,0.75);
    recogen_gengen_g[3][5]->Draw("same");
    recogen_gengen_q[3][5]->Draw("same");
    recogen_gengen[3][5]->Draw("same");
    c2->cd(6)->SetPad(0.25,0.5,0.5,0.75);
    recogen_gengen_g[2][5]->Draw("same");
    recogen_gengen_q[2][5]->Draw("same");
    recogen_gengen[2][5]->Draw("same");
    c2->cd(7)->SetPad(0.5,0.5,0.75,0.75);
    recogen_gengen_g[1][5]->Draw("same");
    recogen_gengen_q[1][5]->Draw("same");
    recogen_gengen[1][5]->Draw("same");
    c2->cd(8)->SetPad(0.75,0.5,0.99,0.75);
    recogen_gengen_g[0][5]->Draw("same");
    recogen_gengen_q[0][5]->Draw("same");
    recogen_gengen[0][5]->Draw("same");
    c2->cd(9)->SetPad(0.01,0.25,0.25,0.5);
    recogen_gengen_g[3][6]->Draw("same");
    recogen_gengen_q[3][6]->Draw("same");
    recogen_gengen[3][6]->Draw("same");
    c2->cd(10)->SetPad(0.25,0.25,0.5,0.5);
    recogen_gengen_g[2][6]->Draw("same");
    recogen_gengen_q[2][6]->Draw("same");
    recogen_gengen[2][6]->Draw("same");
    c2->cd(11)->SetPad(0.5,0.25,0.75,0.5);
    recogen_gengen_g[1][6]->Draw("same");
    recogen_gengen_q[1][6]->Draw("same");
    recogen_gengen[1][6]->Draw("same");
    c2->cd(12)->SetPad(0.75,0.25,0.99,0.5);
    recogen_gengen_g[0][6]->Draw("same");
    recogen_gengen_q[0][6]->Draw("same");
    recogen_gengen[0][6]->Draw("same");
    c2->cd(13)->SetPad(0.01,0.01,0.25,0.25);
    recogen_gengen_g[3][7]->Draw("same");
    recogen_gengen_q[3][7]->Draw("same");
    recogen_gengen[3][7]->Draw("same");
    c2->cd(14)->SetPad(0.25,0.01,0.5,0.25);
    recogen_gengen_g[2][7]->Draw("same");
    recogen_gengen_q[2][7]->Draw("same");
    recogen_gengen[2][7]->Draw("same");   
    c2->cd(15)->SetPad(0.5,0.01,0.75,0.25);
    recogen_gengen_g[1][7]->Draw("same");
    recogen_gengen_q[1][7]->Draw("same");
    recogen_gengen[1][7]->Draw("same");
    c2->cd(16)->SetPad(0.75,0.01,0.99,0.25);
    recogen_gengen_g[0][7]->Draw("same");
    recogen_gengen_q[0][7]->Draw("same");
    recogen_gengen[0][7]->Draw("same");

    TCanvas *c3 = new TCanvas("c3","",1200,1200);
    c3->Divide(4,4);
    gStyle->SetOptStat(0); 

    c3->cd(1)->SetPad(0.01,0.75,0.25,0.99);
    recogen_gengen_g[3][0]->Draw("same");
    recogen_gengen_q[3][0]->Draw("same");
    recogen_gengen[3][0]->Draw("same");

    leg2->Draw("same"); 

    c3->cd(2)->SetPad(0.25,0.75,0.5,0.99);
    recogen_gengen_g[2][0]->Draw("same");
    recogen_gengen_q[2][0]->Draw("same");

    leg1->Draw("same");

    recogen_gengen[2][0]->Draw("same");
    c3->cd(3)->SetPad(0.5,0.75,0.75,0.99);
    recogen_gengen_g[1][0]->Draw("same");
    recogen_gengen_q[1][0]->Draw("same");
    recogen_gengen[1][0]->Draw("same");
    c3->cd(4)->SetPad(0.75,0.75,0.99,0.99);
    recogen_gengen_g[0][0]->Draw("same");
    recogen_gengen_q[0][0]->Draw("same");
    recogen_gengen[0][0]->Draw("same");
    c3->cd(5)->SetPad(0.01,0.5,0.25,0.75);
    recogen_gengen_g[3][1]->Draw("same");
    recogen_gengen_q[3][1]->Draw("same");
    recogen_gengen[3][1]->Draw("same");
    c3->cd(6)->SetPad(0.25,0.5,0.5,0.75);
    recogen_gengen_g[2][1]->Draw("same");
    recogen_gengen_q[2][1]->Draw("same");
    recogen_gengen[2][1]->Draw("same");
    c3->cd(7)->SetPad(0.5,0.5,0.75,0.75);
    recogen_gengen_g[1][1]->Draw("same");
    recogen_gengen_q[1][1]->Draw("same");
    recogen_gengen[1][1]->Draw("same");
    c3->cd(8)->SetPad(0.75,0.5,0.99,0.75);
    recogen_gengen_g[0][1]->Draw("same");
    recogen_gengen_q[0][1]->Draw("same");
    recogen_gengen[0][1]->Draw("same");
    c3->cd(9)->SetPad(0.01,0.25,0.25,0.5);
    recogen_gengen_g[3][2]->Draw("same");
    recogen_gengen_q[3][2]->Draw("same");
    recogen_gengen[3][2]->Draw("same");
    c3->cd(10)->SetPad(0.25,0.25,0.5,0.5);
    recogen_gengen_g[2][2]->Draw("same");
    recogen_gengen_q[2][2]->Draw("same");
    recogen_gengen[2][2]->Draw("same");
    c3->cd(11)->SetPad(0.5,0.25,0.75,0.5);
    recogen_gengen_g[1][2]->Draw("same");
    recogen_gengen_q[1][2]->Draw("same");
    recogen_gengen[1][2]->Draw("same");
    c3->cd(12)->SetPad(0.75,0.25,0.99,0.5);
    recogen_gengen_g[0][2]->Draw("same");
    recogen_gengen_q[0][2]->Draw("same");
    recogen_gengen[0][2]->Draw("same");
    c3->cd(13)->SetPad(0.01,0.01,0.25,0.25);
    recogen_gengen_g[3][3]->Draw("same");
    recogen_gengen_q[3][3]->Draw("same");
    recogen_gengen[3][3]->Draw("same");
    c3->cd(14)->SetPad(0.25,0.01,0.5,0.25);
    recogen_gengen_g[2][3]->Draw("same");
    recogen_gengen_q[2][3]->Draw("same");
    recogen_gengen[2][3]->Draw("same");   
    c3->cd(15)->SetPad(0.5,0.01,0.75,0.25);
    recogen_gengen_g[1][3]->Draw("same");
    recogen_gengen_q[1][3]->Draw("same");
    recogen_gengen[1][3]->Draw("same");
    c3->cd(16)->SetPad(0.75,0.01,0.99,0.25);
    recogen_gengen_g[0][3]->Draw("same");
    recogen_gengen_q[0][3]->Draw("same");
    recogen_gengen[0][3]->Draw("same");

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
