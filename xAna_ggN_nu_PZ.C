#include <fstream>
#include <vector>
#include <iostream>
#include <TLine.h>
#include <TH1D.h>
#include <TRandom.h>
#include <TLorentzVector.h>
#include <math.h>
#include <algorithm> 
#include "Utilities.h"
#include <TProfile.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <string>

//#include "Utilities.h"
//#include "puweicalc.h"
#include "untuplizer.h"
//#include "New_untuplizer/untuplizer.h"
//#include "ElectronSelections.h"
//#include "MuonSelections.h"
//#include "PhotonSelections.h"
//#include "sigmaEff.C"
#include "TMVA_regression_nu.h"

void xAna_ggN_nu_PZ
(const char* outpath1 = "/data3/PZLai/BJetRegression/mini_tree/re_Leading_PZ.root",const char* outpath2 = "/data3/PZLai/BJetRegression/mini_tree/re_Trailing_PZ.root")
{
  // open tree(s) with events to be processed
  // print out which type of variables should be associated with tree branches
  //data.Print();
  //int ShowCout = 0;
  //int remove999 = 1;
	
  TreeReader data("/data1/pwang/MC/V08_00_11_01/job_spring16_TT_powheg_ext3.root","ggNtuplizer/EventTree");
  //Prepare Tree
	
  Int_t jet1nPVs;
  Float_t jet1Pt, jet1Corr, jet1Eta, jet1Mt, jet1LeadTrackPt, jet1LepPtRel,
    jet1LepPt, jet1LepdR, jet1neHEF, jet1neEmEF, jet1vtxPt, jet1vtxMass,
    jet1vtx3dL, jet1vtxNtrk, jet1vtx3deL, jet1PFMET, jet1METDPhi;
  Float_t jet1resolution;
    
  TTree* jet1 = new TTree("jet1","mini_tree");
  //18 training var
  jet1->Branch("jet1Pt",           &jet1Pt,          "jet1Pt/F");
  jet1->Branch("jet1Corr",         &jet1Corr,        "jet1Corr/F");
  jet1->Branch("jet1nPVs",         &jet1nPVs,        "jet1nPVs/I");
  jet1->Branch("jet1Eta",          &jet1Eta,         "jet1Eta/F");
  jet1->Branch("jet1Mt",           &jet1Mt,          "jet1Mt/F");
  jet1->Branch("jet1LeadTrackPt",  &jet1LeadTrackPt, "jet1LeadTrackPt/F");
  jet1->Branch("jet1LepPtRel",     &jet1LepPtRel,    "jet1LepPtRel/F");
  jet1->Branch("jet1LepPt",        &jet1LepPt,       "jet1LepPt/F");
  jet1->Branch("jet1LepdR",        &jet1LepdR,       "jet1LepdR/F");
  jet1->Branch("jet1neHEF",        &jet1neHEF,       "jet1neHEF/F");
  jet1->Branch("jet1neEmEF",       &jet1neEmEF,      "jet1neEmEF/F");
  jet1->Branch("jet1vtxPt",        &jet1vtxPt,       "jet1vtxPt/F");
  jet1->Branch("jet1vtxMass",      &jet1vtxMass,     "jet1vtxMass/F");
  jet1->Branch("jet1vtx3dL",       &jet1vtx3dL,      "jet1vtx3dL/F");
  jet1->Branch("jet1vtxNtrk",      &jet1vtxNtrk,     "jet1vtxNtrk/F");
  jet1->Branch("jet1vtx3deL",      &jet1vtx3deL,     "jet1vtx3deL/F");
  jet1->Branch("jet1PFMET",        &jet1PFMET,       "jet1PFMET/F");
  jet1->Branch("jet1METDPhi",      &jet1METDPhi,     "jet1METDPhi/F");
  jet1->Branch("jet1resolution",   &jet1resolution,      "jet1resolution");

  Int_t jet2nPVs;
  Float_t jet2Pt, jet2Corr, jet2Eta, jet2Mt, jet2LeadTrackPt, jet2LepPtRel,
    jet2LepPt, jet2LepdR, jet2neHEF, jet2neEmEF, jet2vtxPt, jet2vtxMass,
    jet2vtx3dL, jet2vtxNtrk, jet2vtx3deL, jet2PFMET, jet2METDPhi;
  Float_t jet2resolution;
    
  TTree* jet2 = new TTree("jet2","mini_tree");
  //18 training var
  jet2->Branch("jet2Pt",           &jet2Pt,          "jet2Pt/F");
  jet2->Branch("jet2Corr",         &jet2Corr,        "jet2Corr/F");
  jet2->Branch("jet2nPVs",         &jet2nPVs,        "jet2nPVs/I");
  jet2->Branch("jet2Eta",          &jet2Eta,         "jet2Eta/F");
  jet2->Branch("jet2Mt",           &jet2Mt,          "jet2Mt/F");
  jet2->Branch("jet2LeadTrackPt",  &jet2LeadTrackPt, "jet2LeadTrackPt/F");
  jet2->Branch("jet2LepPtRel",     &jet2LepPtRel,    "jet2LepPtRel/F");
  jet2->Branch("jet2LepPt",        &jet2LepPt,       "jet2LepPt/F");
  jet2->Branch("jet2LepdR",        &jet2LepdR,       "jet2LepdR/F");
  jet2->Branch("jet2neHEF",        &jet2neHEF,       "jet2neHEF/F");
  jet2->Branch("jet2neEmEF",       &jet2neEmEF,      "jet2neEmEF/F");
  jet2->Branch("jet2vtxPt",        &jet2vtxPt,       "jet2vtxPt/F");
  jet2->Branch("jet2vtxMass",      &jet2vtxMass,     "jet2vtxMass/F");
  jet2->Branch("jet2vtx3dL",       &jet2vtx3dL,      "jet2vtx3dL/F");
  jet2->Branch("jet2vtxNtrk",      &jet2vtxNtrk,     "jet2vtxNtrk/F");
  jet2->Branch("jet2vtx3deL",      &jet2vtx3deL,     "jet2vtx3deL/F");
  jet2->Branch("jet2PFMET",        &jet2PFMET,       "jet2PFMET/F");
  jet2->Branch("jet2METDPhi",      &jet2METDPhi,     "jet2METDPhi/F");
  jet2->Branch("jet2resolution",   &jet2resolution,      "jet2resolution");
    

    
  //Event loop
  //for (Long64_t ev = data.GetEntriesFast(); ev > 0; --ev) {
  for (Long64_t ev = 0; ev < data.GetEntriesFast()/20; ev++) {
    //for (Long64_t ev = 0; ev < 1000000; ev++) {
    if (ev % 50000 == 0)
      fprintf(stderr, "Processing event %lli of %lli\n", ev+1 , data.GetEntriesFast());
        
    data.GetEntry(ev);
    //start process
    // b_quark
    Int_t		nMC			=	data.GetInt("nMC");
    Int_t*        mcPID		=	data.GetPtrInt("mcPID");
    Float_t*	mcPt		=	data.GetPtrFloat("mcPt");
    Float_t*	mcE			=	data.GetPtrFloat("mcE");
    Float_t*	mcEta		=	data.GetPtrFloat("mcEta");
    Float_t*	mcPhi		=	data.GetPtrFloat("mcPhi");
    Float_t*	mcMass		=	data.GetPtrFloat("mcMass");
    Int_t*        mcMomPID	=	data.GetPtrInt("mcMomPID");
    Float_t*	mcMomMass	=	data.GetPtrFloat("mcMomMass");
    // Genjet
    Float_t*	jetGenEn	=	data.GetPtrFloat("jetGenEn");
    Float_t*	jetGenPt	=	data.GetPtrFloat("jetGenPt");
    Float_t*	jetGenEta	=	data.GetPtrFloat("jetGenEta");
    Float_t*	jetGenPhi	=	data.GetPtrFloat("jetGenPhi");
    // Recojet
    Int_t		nJet		=	data.GetInt("nJet");
    Float_t*	jetEn		=	data.GetPtrFloat("jetEn");
    Float_t*	jetPt		=	data.GetPtrFloat("jetPt");
    Float_t*	jetEta		=	data.GetPtrFloat("jetEta");
    Float_t*	jetPhi		=	data.GetPtrFloat("jetPhi");
    Int_t*      jetPartonID	=	data.GetPtrInt("jetPartonID");
    // JetGenJet
    Float_t*	jetGenJetEn		=	data.GetPtrFloat("jetGenJetEn");
    Float_t*	jetGenJetPt		=	data.GetPtrFloat("jetGenJetPt");
    Float_t*	jetGenJetEta	=	data.GetPtrFloat("jetGenJetEta");
    Float_t*	jetGenJetPhi	=	data.GetPtrFloat("jetGenJetPhi");
    Int_t*      jetGenPartonID	=	data.GetPtrInt("jetGenPartonID");
    Int_t*      jetGenPartonMomID	=	data.GetPtrInt("jetGenPartonMomID");
        
    // 18 variables
    Int_t		nPVs      			  =	data.GetInt("nVtx");
    //Float_t* jetPt                    = data.GetPtrFloat("jetPt");
    Float_t* jetRawPt                 = data.GetPtrFloat("jetRawPt");
    Float_t  rho                      = data.GetFloat("rho");
    //Float_t* jetEta                   = data.GetPtrFloat("jetEta");
    //Float_t* jetPhi                   = data.GetPtrFloat("jetPhi");
    Float_t* jetMt                    = data.GetPtrFloat("jetMt");
    Float_t* jetLeadTrackPt           = data.GetPtrFloat("jetLeadTrackPt");
    Float_t* jetLepTrackPt            = data.GetPtrFloat("jetLepTrackPt");
    Float_t* jetLepTrackEta           = data.GetPtrFloat("jetLepTrackEta");
    Float_t* jetLepTrackPhi           = data.GetPtrFloat("jetLepTrackPhi");
    Float_t* jetNHF                   = data.GetPtrFloat("jetNHF");
    Float_t* jetNEF                   = data.GetPtrFloat("jetNEF");
    Int_t*   jetNCH                   = data.GetPtrInt("jetNCH");
    Float_t* jetVtxPt                 = data.GetPtrFloat("jetVtxPt");
    Float_t* jetVtxMass               = data.GetPtrFloat("jetVtxMass");
    Float_t* jetVtx3DVal              = data.GetPtrFloat("jetVtx3DVal");
    Float_t* jetVtxNtrks              = data.GetPtrFloat("jetVtxNtrks");
    Float_t* jetVtx3DSig              = data.GetPtrFloat("jetVtx3DSig");
    Float_t  pfMET                    = data.GetFloat("pfMET");
    Float_t  pfMETPhi                 = data.GetFloat("pfMETPhi");
        
    //jet_Ptcor_jetGenJet_18	= TMVA_18_jetGenJet_nu(data, i) * jetPt;
     
      
      vector<int> mc1List;
      mc1List.clear();
      vector<int> nuList;
      nuList.clear();
      vector<int> nuList_usd[2];
      nuList_usd[0].clear();
      nuList_usd[1].clear();
      for (int i=0; i < nMC; i++){
          if( abs(mcPID[i])== 12 || abs(mcPID[i])== 14 || abs(mcPID[i])== 16) nuList.push_back(i);
          if (abs(mcPID[i])!=5) continue;//b quark
          if (abs(mcMomPID[i])!=6) continue;//t quark
          
          mc1List.push_back(i);
      }
      
        //sort
        if (mc1List.size()!=0 ){
            for (size_t n=0; n<=mc1List.size() - 1 ;n++){
                int temp;
                for (size_t j = 0; j < mc1List.size() - 1 - n; j++){
                    if (mcPt[mc1List[j]] < mcPt[mc1List[j+1]]) {
                        temp = mc1List[j];
                        mc1List[j] = mc1List[j+1];
                        mc1List[j+1] = temp;
                    }
                }
            }
        }
        
      float myResolution1 = -1;
      float myResolution2 = -1;
      vector<int> jetList;
      jetList.clear();
      int matchNuCout=0;
      for (int i=0; i < nJet ; i++){
          if( fabs(jetGenPartonID[i]) != 5) continue;
          if( fabs(jetGenPartonMomID[i]) != 6) continue;
          if( jetPt[i] < 0 || jetGenJetPt[i] < 0) continue;
          if( ev == 66450047) continue;
          if(fabs(jetEta[i])>2.4) continue;
          if (jetGenJetEn[i]<0) matchNuCout=1;
          jetList.push_back(i);
      }//nJetloop
      //=====================================================================================================================================
      
      if( jetList.size() < 2) continue;
      if (mc1List.size() < 2) continue;
      
      
      vector<int> acc_jet;
      acc_jet.clear();
      vector<int> acc_b1;
      acc_b1.clear();
      for(size_t i=0;i<jetList.size();i++){
          //pair 1 matching
          int flag_p1 = 0;
          for(size_t j=0;j<mc1List.size();j++){
              if (deltaR(mcEta[mc1List[j]],mcPhi[mc1List[j]],jetEta[jetList[i]],jetPhi[jetList[i]])<0.35){
                  int flag=0;
                  for(size_t k=0;k<acc_b1.size();k++){
                      if(acc_b1[k]==mc1List[j])flag=1;
                  }
                  if (flag==1)continue;
                  acc_jet.push_back(jetList[i]);
                  acc_b1.push_back(mc1List[j]);
                  flag_p1=1;
                  break;
              }
          }
      }
      
      if( acc_jet.size() < 2) continue;
      if (mc1List.size()<2)continue;
      
      float jet1_Ptcor_jetGenJet_18, jet2_Ptcor_jetGenJet_18;
      jet1_Ptcor_jetGenJet_18	=  TMVA_18_jetGenJet_Leading_nu(data, acc_jet[0]) * jetPt[acc_jet[0]];
      jet2_Ptcor_jetGenJet_18	=  TMVA_18_jetGenJet_Trailing_nu(data, acc_jet[1]) * jetPt[acc_jet[1]];
      if( jet1_Ptcor_jetGenJet_18 < 20.) continue;
      if( jet2_Ptcor_jetGenJet_18 < 20.) continue;

      
      //cout << "Hello1" << endl;
      TLorentzVector thisGenjet[2];
      TLorentzVector thisRecojet[2];
      
      thisGenjet[0].SetPtEtaPhiE(jetGenJetPt[acc_jet[0]], jetGenJetEta[acc_jet[0]], jetGenJetPhi[acc_jet[0]], jetGenJetEn[acc_jet[0]]);
      thisGenjet[1].SetPtEtaPhiE(jetGenJetPt[acc_jet[1]], jetGenJetEta[acc_jet[1]], jetGenJetPhi[acc_jet[1]], jetGenJetEn[acc_jet[1]]);
      thisRecojet[0].SetPtEtaPhiE(jet1_Ptcor_jetGenJet_18,jetEta[acc_jet[0]],jetPhi[acc_jet[0]],jetEn[acc_jet[0]]);
      thisRecojet[1].SetPtEtaPhiE(jet2_Ptcor_jetGenJet_18,jetEta[acc_jet[1]],jetPhi[acc_jet[1]],jetEn[acc_jet[1]]);
      //cout << "Hello2" << endl;
      //adding nutrino
      //cout << "Hello2" << endl;
      if (nuList.size()!=0 && acc_jet.size()>1){
          for (size_t j=0;j<nuList.size();j++){
              int matched_nuJet = 999;
              vector<int> pass_nuJet;
              pass_nuJet.clear();
              if(acc_jet.size()>1){
                  if (deltaR(jetGenJetEta[acc_jet[0]],jetGenJetPhi[acc_jet[0]],mcEta[nuList[j]],mcPhi[nuList[j]])<0.4) pass_nuJet.push_back(acc_jet[0]);
                  if (deltaR(jetGenJetEta[acc_jet[1]],jetGenJetPhi[acc_jet[1]],mcEta[nuList[j]],mcPhi[nuList[j]])<0.4) pass_nuJet.push_back(acc_jet[1]);
              }
              if (pass_nuJet.size() > 0){
                  //Find the closer Jets!
                  float DR1 = deltaR(jetGenJetEta[pass_nuJet[0]],jetGenJetPhi[pass_nuJet[0]],mcEta[nuList[j]],mcPhi[nuList[j]]);
                  matched_nuJet = pass_nuJet[0];
                  for (size_t n=1;n<pass_nuJet.size();n++){
                      float DR2 = deltaR(jetGenJetEta[pass_nuJet[n]],jetGenJetPhi[pass_nuJet[n]],mcEta[nuList[j]],mcPhi[nuList[j]]);
                      if (DR2<DR1) {
                          matched_nuJet = pass_nuJet[n];
                          DR1 = DR2;
                      }
                  }
                  if(acc_jet.size()>1) if (matched_nuJet == acc_jet[0] )nuList_usd[0].push_back(nuList[j]);
                  if(acc_jet.size()>1) if (matched_nuJet == acc_jet[1] )nuList_usd[1].push_back(nuList[j]);
              }
              pass_nuJet.clear();
          }
      }
      if (matchNuCout==1)for (size_t n=0;n<nuList_usd[0].size();n++) cout<<nuList_usd[0][n]<<" "<<deltaR(jetGenJetEta[acc_jet[0]],jetGenJetPhi[acc_jet[0]],mcEta[nuList_usd[0][n]],mcPhi[nuList_usd[0][n]])<<endl;
      if (matchNuCout==1)for (size_t n=0;n<nuList_usd[1].size();n++) cout<<nuList_usd[1][n]<<" "<<deltaR(jetGenJetEta[acc_jet[1]],jetGenJetPhi[acc_jet[1]],mcEta[nuList_usd[1][n]],mcPhi[nuList_usd[1][n]])<<endl;
      //Adding neutrinos
      if(acc_jet.size()>1){
          for (size_t n=0;n<nuList_usd[0].size();n++){
              TLorentzVector nu;
              nu.SetPtEtaPhiE(mcPt[nuList_usd[0][n]],mcEta[nuList_usd[0][n]],mcPhi[nuList_usd[0][n]],mcE[nuList_usd[0][n]]);
              thisGenjet[0] += nu;
          }
          for (size_t n=0;n<nuList_usd[1].size();n++){
              TLorentzVector nu;
              nu.SetPtEtaPhiE(mcPt[nuList_usd[1][n]],mcEta[nuList_usd[1][n]],mcPhi[nuList_usd[1][n]],mcE[nuList_usd[1][n]]);
              thisGenjet[1] += nu;
          }
      }

//mycode
      
      
    //=================================================================================================================
        
        
/*
      //jet energy regression
    vector<int> acc_mc;
    acc_mc.clear();
    vector<int> nuList;
    nuList.clear();
    vector<int> nuList_usd[2];
    nuList_usd[0].clear();
    nuList_usd[1].clear();
    for (int i=0; i < nMC; i++){
      if( abs(mcPID[i])== 12 || abs(mcPID[i])== 14 || abs(mcPID[i])== 16) nuList.push_back(i);
      if (abs(mcPID[i])!=5) continue;//b quark
      if (abs(mcMomPID[i])!=6) continue;//t quark
                
      acc_mc.push_back(i);
    }
    //sort
    if (acc_mc.size()!=0 ){
        for (size_t n=0; n<=acc_mc.size() - 1 ;n++){
            int temp;
            for (size_t j = 0; j < acc_mc.size() - 1 - n; j++){
                if (mcPt[acc_mc[j]] < mcPt[acc_mc[j+1]]) {
                    temp = acc_mc[j];
                    acc_mc[j] = acc_mc[j+1];
                    acc_mc[j+1] = temp;
                }
            }
        }
    }
        
    //===========================================================================================================================
    //b jet selection ver 2.0 9/5
    vector<int> Bjet1;
    vector<int> Bjet2;
    Bjet1.clear();
    Bjet2.clear();
    //loop all pairs
    for (int n1=0;n1<nJet;n1++){
        for(int n2=n1+1;n2<nJet;n2++){
            Bjet1.push_back(n1);
            Bjet2.push_back(n2);
        }
    }
    //define leading & trailing
    for (size_t i=0;i<Bjet1.size();i++){
        if(jetPt[Bjet1[i]]<jetPt[Bjet2[i]]){
            int tmp;
            tmp = Bjet1[i];
            Bjet1[i] = Bjet2[i];
            Bjet2[i] = tmp;
        }
    }
    vector<float> bjRegfactor[2];
    bjRegfactor[0].clear();
    bjRegfactor[1].clear();
        for (size_t n=0;n<Bjet1.size();n++){
        bjRegfactor[0].push_back(TMVA_18_jetGenJet_Leading_nu(data, Bjet1[n]));
        
        }
        
        for (size_t n=0;n<Bjet2.size();n++){
        bjRegfactor[1].push_back(TMVA_18_jetGenJet_Trailing_nu(data, Bjet2[n]));
            
        }
        
    vector<int> mcList; mcList.clear();
    vector<int> acc_jet; acc_jet.clear();
    vector<float> acc_jet_reg; acc_jet_reg.clear();
        for (size_t n=0;n<Bjet1.size();n++){
            //general cut
            int cut_flag=0; int i;
            for (int a=0;a<2;a++){
                if (a==0) i=Bjet1[n];
                if (a==1) i=Bjet2[n];
                
                if( fabs(jetGenPartonID[i]) != 5)       cut_flag = 1;//add
                if( fabs(jetGenPartonMomID[i]) != 6)    cut_flag = 1;//add
                if(jetPt[i]*bjRegfactor[a][n]<20)		cut_flag = 1;
                if(fabs(jetEta[i])>2.4)					cut_flag = 1;
                if(jetGenPt[i]<0 || jetGenJetPt[i]<0)	cut_flag = 1;
            }
            if (cut_flag==1) continue;
            cut_flag = 0;
        
            // b matching
            int mc[2];
            for (int a=0;a<2;a++){
                if (a==0) i=Bjet1[n];
                if (a==1) i=Bjet2[n];
                vector<float> mc_dR;
                mc_dR.clear();
                vector<int> mc_dR_index;
                mc_dR_index.clear();
                for(size_t m=0;m<acc_mc.size();m++){
                    if (deltaR(mcEta[acc_mc[m]],mcPhi[acc_mc[m]],jetEta[i],jetPhi[i])>0.35)continue;
                    mc_dR_index.push_back(acc_mc[m]);
                    mc_dR.push_back(deltaR(mcEta[acc_mc[m]],mcPhi[acc_mc[m]],jetEta[i],jetPhi[i]));
                }
                if(mc_dR.size()>0){
                    //dR sorting
                    for (size_t m=0; m<mc_dR.size() - 1 ;m++){
                        float temp;
                        for (size_t j = 0; j < mc_dR.size() - 1 - m; j++){
                            if (mc_dR[j] > mc_dR[j+1]) {
                                temp = mc_dR[j];
                                mc_dR[j] = mc_dR[j+1];
                                mc_dR[j+1] = temp;
                                temp = mc_dR_index[j];
                                mc_dR_index[j] = mc_dR_index[j+1];
                                mc_dR_index[j+1] = temp;
                            }
                        }
                    }
                    mc[a]= mc_dR_index[0]; //select minimum dR
                }else{cut_flag=1;}
            }
            
            if (cut_flag==1) continue;
            mcList.push_back(mc[0]);
            mcList.push_back(mc[1]);
            acc_jet.push_back(Bjet1[n]);
            acc_jet.push_back(Bjet2[n]);
            acc_jet_reg.push_back(bjRegfactor[0][n]);
            acc_jet_reg.push_back(bjRegfactor[1][n]);
            break;
        }
  
    if(acc_jet.size() < 2) continue;
        
    TLorentzVector thisGenjet[2];
    thisGenjet[0].SetPtEtaPhiE(jetGenJetPt[acc_jet[0]], jetGenJetEta[acc_jet[0]], jetGenJetPhi[acc_jet[0]], jetGenJetEn[acc_jet[0]]);
    thisGenjet[1].SetPtEtaPhiE(jetGenJetPt[acc_jet[1]], jetGenJetEta[acc_jet[1]], jetGenJetPhi[acc_jet[1]], jetGenJetEn[acc_jet[1]]);
    TLorentzVector thisRecojet[2];
    thisRecojet[0].SetPtEtaPhiE(jetPt[acc_jet[0]],jetEta[acc_jet[0]],jetPhi[acc_jet[0]],jetEn[acc_jet[0]]);
    thisRecojet[1].SetPtEtaPhiE(jetPt[acc_jet[1]],jetEta[acc_jet[1]],jetPhi[acc_jet[1]],jetEn[acc_jet[1]]);

    if (nuList.size()!=0 && acc_jet.size()>1){
        
      for (size_t j=0;j<nuList.size();j++){
	int matched_nuJet = 999;
	vector<int> pass_nuJet;
	pass_nuJet.clear();
	if(acc_jet.size()>1){
	  if (deltaR(jetGenJetEta[acc_jet[0]],jetGenJetPhi[acc_jet[0]],mcEta[nuList[j]],mcPhi[nuList[j]])<0.4) pass_nuJet.push_back(acc_jet[0]);
	  if (deltaR(jetGenJetEta[acc_jet[1]],jetGenJetPhi[acc_jet[1]],mcEta[nuList[j]],mcPhi[nuList[j]])<0.4) pass_nuJet.push_back(acc_jet[1]);
	}
	if (pass_nuJet.size() > 0){
	  //Find the closer Jets!
	  float DR1 = deltaR(jetGenJetEta[pass_nuJet[0]],jetGenJetPhi[pass_nuJet[0]],mcEta[nuList[j]],mcPhi[nuList[j]]);
	  matched_nuJet = pass_nuJet[0];
	  for (size_t n=1;n<pass_nuJet.size();n++){
	    float DR2 = deltaR(jetGenJetEta[pass_nuJet[n]],jetGenJetPhi[pass_nuJet[n]],mcEta[nuList[j]],mcPhi[nuList[j]]);
	    if (DR2<DR1) {
	      matched_nuJet = pass_nuJet[n];
	      DR1 = DR2;
	    }
	  }
	  if(acc_jet.size()>1) if (matched_nuJet == acc_jet[0] )nuList_usd[0].push_back(nuList[j]);
	  if(acc_jet.size()>1) if (matched_nuJet == acc_jet[1] )nuList_usd[1].push_back(nuList[j]);
	}
	pass_nuJet.clear();
      }
    }
    //Adding neutrinos
    if(acc_jet.size()>1){
      for (size_t n=0;n<nuList_usd[0].size();n++){
	TLorentzVector nu;
	nu.SetPtEtaPhiE(mcPt[nuList_usd[0][n]],mcEta[nuList_usd[0][n]],mcPhi[nuList_usd[0][n]],mcE[nuList_usd[0][n]]);
	thisGenjet[0] += nu;
      }
      for (size_t n=0;n<nuList_usd[1].size();n++){
	TLorentzVector nu;
	nu.SetPtEtaPhiE(mcPt[nuList_usd[1][n]],mcEta[nuList_usd[1][n]],mcPhi[nuList_usd[1][n]],mcE[nuList_usd[1][n]]);
	thisGenjet[1] += nu;
      }
    }
*/
        
      
    /*
      float myResolution1 = -1;
      float myResolution2 = -1;
      TLorentzVector thisGenjet[2];
      TLorentzVector thisRecojet[2];
      vector<int> acc_jet;
      acc_jet.clear();
      int matchNuCout=0;
      for (int i=0; i < nJet ; i++){
      //selct bjet
      if( ev == 66450047) continue;
      if( jetGenPartonID[i] == -99) continue;
      if( fabs(jetGenPartonMomID[i]) != 6) continue;
      if( fabs(jetGenPartonID[i]) != 5) continue;
      if( jetGenJetPt[i] < 0 || jetPt[i] < 0) continue;
      if (jetGenJetEn[i]<0) matchNuCout=1;
            
      acc_jet.push_back(i);

      }//nJetloop
            
      if( acc_jet.size() < 2) continue;
      float jet1_Ptcor_jetGenJet_18, jet2_Ptcor_jetGenJet_18;
      jet1_Ptcor_jetGenJet_18	=  TMVA_18_jetGenJet_Leading_nu(data, acc_jet[0]) * jetPt[acc_jet[0]];
      if( jet1_Ptcor_jetGenJet_18 < 20.) continue;
      jet2_Ptcor_jetGenJet_18	=  TMVA_18_jetGenJet_Trailing_nu(data, acc_jet[1]) * jetPt[acc_jet[1]];
      if( jet2_Ptcor_jetGenJet_18 < 20.) continue;
        
      //cout << "hello2" << endl;
      thisGenjet[0].SetPtEtaPhiE(jetGenJetPt[acc_jet[0]], jetGenJetEta[acc_jet[0]], jetGenJetPhi[acc_jet[0]], jetGenJetEn[acc_jet[0]]);
      thisGenjet[1].SetPtEtaPhiE(jetGenJetPt[acc_jet[1]], jetGenJetEta[acc_jet[1]], jetGenJetPhi[acc_jet[1]], jetGenJetEn[acc_jet[1]]);
      thisRecojet[0].SetPtEtaPhiE(jet1_Ptcor_jetGenJet_18,jetEta[acc_jet[0]],jetPhi[acc_jet[0]],jetEn[acc_jet[0]]);
      thisRecojet[1].SetPtEtaPhiE(jet2_Ptcor_jetGenJet_18,jetEta[acc_jet[1]],jetPhi[acc_jet[1]],jetEn[acc_jet[1]]);
      //cout << "hello3" << endl;
      //adding nutrino
      vector<int> nuList;
      nuList.clear();
      vector<int> nuList_usd[2];
      nuList.clear();
      nuList_usd[0].clear();
      nuList_usd[1].clear();
      for (int nupar=0; nupar < nMC; nupar++){
      if( abs(mcPID[nupar])== 12 || abs(mcPID[nupar])== 14 || abs(mcPID[nupar])== 16){
      nuList.push_back(nupar);
      }
      }
        
      if (nuList.size()!=0 && acc_jet.size()>1){
      for (size_t j=0;j<nuList.size();j++){
      int matched_nuJet = 999;
      vector<int> pass_nuJet;
      pass_nuJet.clear();
      if(acc_jet.size()>1){
      if (deltaR(jetGenJetEta[acc_jet[0]],jetGenJetPhi[acc_jet[0]],mcEta[nuList[j]],mcPhi[nuList[j]])<0.4) pass_nuJet.push_back(acc_jet[0]);
      if (deltaR(jetGenJetEta[acc_jet[1]],jetGenJetPhi[acc_jet[1]],mcEta[nuList[j]],mcPhi[nuList[j]])<0.4) pass_nuJet.push_back(acc_jet[1]);
      }
                
      if (pass_nuJet.size() > 0){
      //Find the closer Jets!
      float DR1 = deltaR(jetGenJetEta[pass_nuJet[0]],jetGenJetPhi[pass_nuJet[0]],mcEta[nuList[j]],mcPhi[nuList[j]]);
      matched_nuJet = pass_nuJet[0];
      for (size_t n=1;n<pass_nuJet.size();n++){
      float DR2 = deltaR(jetGenJetEta[pass_nuJet[n]],jetGenJetPhi[pass_nuJet[n]],mcEta[nuList[j]],mcPhi[nuList[j]]);
      if (DR2<DR1) {
      matched_nuJet = pass_nuJet[n];
      DR1 = DR2;
      }
      }
      if(acc_jet.size()>1) if (matched_nuJet == acc_jet[0] )nuList_usd[0].push_back(nuList[j]);
      if(acc_jet.size()>1) if (matched_nuJet == acc_jet[1] )nuList_usd[1].push_back(nuList[j]);
      }
      }
      }
      //cout << "hello4" << endl;
      if (matchNuCout==1)for (size_t n=0;n<nuList_usd[0].size();n++) cout<<nuList_usd[0][n]<<" "<<deltaR(jetGenJetEta[acc_jet[0]],jetGenJetPhi[acc_jet[0]],mcEta[nuList_usd[0][n]],mcPhi[nuList_usd[0][n]])<<endl;
      if (matchNuCout==1)for (size_t n=0;n<nuList_usd[1].size();n++) cout<<nuList_usd[1][n]<<" "<<deltaR(jetGenJetEta[acc_jet[1]],jetGenJetPhi[acc_jet[1]],mcEta[nuList_usd[1][n]],mcPhi[nuList_usd[1][n]])<<endl;
      //Adding neutrinos
      if(acc_jet.size()>1){
      for (size_t n=0;n<nuList_usd[0].size();n++){
      TLorentzVector nu;
      nu.SetPtEtaPhiE(mcPt[nuList_usd[0][n]],mcEta[nuList_usd[0][n]],mcPhi[nuList_usd[0][n]],mcE[nuList_usd[0][n]]);
      thisGenjet[0] += nu;
      }
      for (size_t n=0;n<nuList_usd[1].size();n++){
      TLorentzVector nu;
      nu.SetPtEtaPhiE(mcPt[nuList_usd[1][n]],mcEta[nuList_usd[1][n]],mcPhi[nuList_usd[1][n]],mcE[nuList_usd[1][n]]);
      thisGenjet[1] += nu;
      }
      }*/
      //
      
      
    jet1resolution =  (thisRecojet[0].Pt() - thisGenjet[0].Pt()) / thisGenjet[0].Pt();
    jet2resolution =  (thisRecojet[1].Pt() - thisGenjet[1].Pt()) / thisGenjet[1].Pt();
        
    if (acc_jet.size()<2)continue;
    jet1Pt           = thisRecojet[0].Pt();//Jet pT
    jet1Corr         = thisRecojet[0].Pt()/jetRawPt[acc_jet[0]];//JEC
    jet1nPVs           = nPVs;
    jet1Eta          = jetEta[acc_jet[0]];//Jet η
    jet1Mt           = jetMt[acc_jet[0]];//Jet transverse mass
    jet1LeadTrackPt  = (jetLeadTrackPt[acc_jet[0]] < 0) ? 0. : jetLeadTrackPt[acc_jet[0]];//pTLeadTrk, transverse momentum of the leading track in the jet
    jet1LepPtRel  = (jetLepTrackPt[acc_jet[0]] < 0) ? 0. : jetLepTrackPt[acc_jet[0]]/jetPt[acc_jet[0]];//Soft Lepton pTRel, relative transverse momentum of soft lepton candidate in the jet;
    jet1LepPt     = (jetLepTrackPt[acc_jet[0]] < 0) ? 0. : jetLepTrackPt[acc_jet[0]];//Soft Lepton pT, transverse momentum of soft lepton candidate in the jet;
    jet1LepdR = (jetLepTrackPt[acc_jet[0]] < 0) ? 0. : deltaR(jetEta[acc_jet[0]], jetPhi[acc_jet[0]], jetLepTrackEta[acc_jet[0]], jetLepTrackPhi[acc_jet[0]]);//Soft Lepton dR, relative distance in the η-phi space of soft lepton candidate in the jet and the jet;
    jet1neHEF        = jetNHF[acc_jet[0]];//Neutral hadron energy fraction
    jet1neEmEF       = jetNEF[acc_jet[0]];//Photon energy fraction
    jet1vtxPt        = jetVtxPt[acc_jet[0]];//SecVtxPt, pT of the jet secondary vertex
    jet1vtxMass      = jetVtxMass[acc_jet[0]];//SecVtxM, Mass of the jet secondary vertex
    jet1vtx3dL       = jetVtx3DVal[acc_jet[0]];//SecVtxdL, the 3-d flight length of the jet secondary vertex
    jet1vtxNtrk      = (Int_t) jetVtxNtrks[acc_jet[0]];//SecVtxNtrk, track multiplicity of the reconstructed secondary vertex
    jet1vtx3deL      = jetVtx3DSig[acc_jet[0]];//SecVtxdeL, Error on the 3-d flight length of the jet secondary vertex
    jet1PFMET        = pfMET;
    jet1METDPhi      = deltaPhi(pfMETPhi,jetPhi[acc_jet[0]]);
    jet1->Fill();
        
    jet2Pt           = thisRecojet[1].Pt();//Jet pT
    jet2Corr         = thisRecojet[1].Pt()/jetRawPt[acc_jet[1]];//JEC
    jet2nPVs           = nPVs;
    jet2Eta          = jetEta[acc_jet[1]];//Jet η
    jet2Mt           = jetMt[acc_jet[1]];//Jet transverse mass
    jet2LeadTrackPt  = (jetLeadTrackPt[acc_jet[1]] < 0) ? 0. : jetLeadTrackPt[acc_jet[1]];//pTLeadTrk, transverse momentum of the leading track in the jet
    jet2LepPtRel  = (jetLepTrackPt[acc_jet[1]] < 0) ? 0. : jetLepTrackPt[acc_jet[1]]/jetPt[acc_jet[1]];//Soft Lepton pTRel, relative transverse momentum of soft lepton candidate in the jet;
    jet2LepPt     = (jetLepTrackPt[acc_jet[1]] < 0) ? 0. : jetLepTrackPt[acc_jet[1]];//Soft Lepton pT, transverse momentum of soft lepton candidate in the jet;
    jet2LepdR = (jetLepTrackPt[acc_jet[1]] < 0) ? 0. : deltaR(jetEta[acc_jet[1]], jetPhi[acc_jet[1]], jetLepTrackEta[acc_jet[1]], jetLepTrackPhi[acc_jet[1]]);//Soft Lepton dR, relative distance in the η-phi space of soft lepton candidate in the jet and the jet;
    jet2neHEF        = jetNHF[acc_jet[1]];//Neutral hadron energy fraction
    jet2neEmEF       = jetNEF[acc_jet[1]];//Photon energy fraction
    jet2vtxPt        = jetVtxPt[acc_jet[1]];//SecVtxPt, pT of the jet secondary vertex
    jet2vtxMass      = jetVtxMass[acc_jet[1]];//SecVtxM, Mass of the jet secondary vertex
    jet2vtx3dL       = jetVtx3DVal[acc_jet[1]];//SecVtxdL, the 3-d flight length of the jet secondary vertex
    jet2vtxNtrk      = (Int_t) jetVtxNtrks[acc_jet[1]];//SecVtxNtrk, track multiplicity of the reconstructed secondary vertex
    jet2vtx3deL      = jetVtx3DSig[acc_jet[1]];//SecVtxdeL, Error on the 3-d flight length of the jet secondary vertex
    jet2PFMET        = pfMET;
    jet2METDPhi      = deltaPhi(pfMETPhi,jetPhi[acc_jet[1]]);
    jet2->Fill();
  
  }//eventloop
    
  TFile* fo1 = TFile::Open(outpath1, "RECREATE");
  jet1->Write("training");
  fo1->Write();
    
  delete jet1;
  delete fo1;
    
  TFile* fo2 = TFile::Open(outpath2, "RECREATE");
  jet2->Write("training");
  fo2->Write();

  fprintf(stderr, "Processed all events\n");
  //jet->Write("", TObject::kOverwrite);
    
  delete jet2;
  delete fo2;
    
}//main

