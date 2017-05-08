/*C/C++ Includes{{{*/
#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <map>
#include <cmath>
/*}}}*/

/*ROOT Includes{{{*/
#include <TSystem.h>
#include <TString.h>
#include <TStyle.h>
#include <Riostream.h>
#include "TObjString.h"
#include <TNamed.h>
#include <TPRegexp.h>
#include <TObjArray.h>
#include <TChain.h>
#include <TMath.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TROOT.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TDatime.h>
#include <TError.h>
#include <TVirtualFitter.h>
#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TSQLRow.h>
#include <TCut.h>
#include <TMultiGraph.h>
#include <TCutG.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TPaletteAxis.h>
#include <TRandom3.h>
//#include <TMatrix.h>
/*}}}*/

#include "SIDIS_Acceptance.h"
const double MeV2GeV=0.001;
using namespace std;

const double Deg2Rad = TMath::DegToRad();
const double Rad2Deg = TMath::RadToDeg();
/*Detector Resolutions{{{*/
const double Sigma_Dp_E = 0.02; //2%, energy resolution for electron from GEM tracking
const double Sigma_Theta_E = 0.6/1000.; //0.6mrad, Angular resolution for electron, determined by GEM tracking
const double Sigma_Phi_E =  5.0/1000.; //5mrad, Angular resolution for electron, determined by GEM tracking
const double Sigma_Fermi =  200.0/1000.; //GeV,The average Fermi motion for a nucleon in a nucleus is about 200MeV

const double eMass = 0.511/1000;//electron mass 
const double piMass = 0.13957018;
const double pMass = 0.938272;
const double nMass = 0.939565;
const double he3Mass = 2.*pMass+nMass;

Int_t CheckLaws(TLorentzVector* P_t,TLorentzVector* P_E0, TLorentzVector* P_e, TLorentzVector* P_g, TLorentzVector* P_pro);
Double_t GetMM( TLorentzVector* P_t,TLorentzVector* P_E0, TLorentzVector* P_e, TLorentzVector* P_g);
Double_t GetMP( TLorentzVector* P_t,TLorentzVector* P_E0, TLorentzVector* P_e, TLorentzVector* P_g);
TLorentzVector* GetMV(TLorentzVector* P_t,TLorentzVector* P_E0, TLorentzVector* P_e, TLorentzVector* P_pim);

/*}}}*/

int main(){
    //int main(){
    gStyle->SetOptStat(0);
    //const double Lumi = 1.0e36; // cm-2*s-1, for He3 nuclear not for nucleons
    //const double nBcm2 = 1e-33;
    //const double PI = 3.1415926;
    const TString Target = "He3";
    const TString particle="pim";
    double EBeam = 11.0;
    double time = 48 * 24 * 3600.0;
    
    Int_t CASE = 0; cout<<"--- What Configuration (1->NONE, 2->EL+MS, 3->EL+MS+FSI, 4->Fermi, 5->EL+MS-Fermi)? "; cin>>CASE;/*{{{*/

    Int_t target_direction=1; 
    cout<<"--- Target Direction (1->Up, 2->Down): "; cin>>target_direction;
    TString targetname="";
    if(target_direction==1) targetname="up";
    if(target_direction==2) targetname="down";/*}}}*/

    //CLEO Acceptance
    SIDIS_Acceptance *accpt = new SIDIS_Acceptance();

    TString filename, new_filename;/*{{{*/
    TChain *t1 = new TChain("t1");
    
    //const int fileNO = 1000;
    int fileNO = 999;
    //for(int i=0;i<fileNO;i++){
    for(int i=1;i<fileNO;i++){
        if(CASE==1){
            //if(target_direction==1) filename = Form("./RootFiles_Up_Simple/DEMP_Ee_11_Events_1000000_File_%d.root", i);
            if(target_direction==1) filename = Form("./RootFiles_Up_Simple/DEMP_Ee_11_Events_1000000_File_%d.root", i+1000);
            if(target_direction==2) filename = Form("./RootFiles_Down_Simple/DEMP_Ee_11_Events_1000000_File_%d.root", i);
        }
        
        if(CASE==2){
            if(target_direction==1) filename = Form("./RootFiles_Up_EL/DEMP_Ee_11_Events_1000000_File_%d_Fermi__Eloss__MS_.root", i);
            if(target_direction==2) filename = Form("./RootFiles_Down_EL/DEMP_Ee_11_Events_1000000_File_%d_Fermi__Eloss__MS_.root", i);
        }

        if(CASE==3){
            if(target_direction==1) filename = Form("./RootFiles_Up_EL_FSI/DEMP_Ee_11_Events_1000000_File_%d_Fermi__Eloss__FSI__MS_.root", i);
            if(target_direction==2) filename = Form("./RootFiles_Down_EL_FSI/DEMP_Ee_11_Events_1000000_File_%d_Fermi__Eloss__FSI__MS_.root", i);
        }

        if(CASE==4){
            if(target_direction==1) filename = Form("./RootFiles_Up_Fermi/DEMP_Ee_11_Events_1000000_File_%d_Fermi_.root", i);
            if(target_direction==2) filename = Form("./RootFiles_Down_Fermi/DEMP_Ee_11_Events_1000000_File_%d_Fermi_.root", i);
        }

        if(CASE==5){
            if(target_direction==1) filename = Form("./RootFiles_Up_EL_noFermi/DEMP_Ee_11_Events_1000000_File_%d_Eloss__MS_.root", i);
            if(target_direction==2) filename = Form("./RootFiles_Down_EL_noFermi/DEMP_Ee_11_Events_1000000_File_%d_Eloss__MS_.root", i);
        }
       
        cout<<"--- Reading in file = "<<filename.Data()<<"\r";
        t1->Add(filename.Data());
    }
    int N_entries = t1->GetEntries();
    cout<<"--- Total Events = "<< N_entries<<endl;/*}}}*/
   
    /*Define Variables, Branches{{{*/
    Int_t NRecorded, NGenerated; 
    Double_t W,Qsq, t, t_Para, Epsilon, x, y, z;/*{{{*/
    Double_t t_cor, Qsq_cor, W_cor, x_cor, y_cor, z_cor;
    Double_t tgt_theta, tgt_phi, tgt_ene, tgt_mom, tgt_px, tgt_py, tgt_pz;
    Double_t beam_theta, beam_phi, beam_ene, beam_mom, beam_px, beam_py, beam_pz;
    Double_t beam_theta_cor, beam_phi_cor, beam_ene_cor, beam_mom_cor, beam_px_cor, beam_py_cor, beam_pz_cor;
    Double_t ele_theta, ele_phi, ele_ene, ele_mom, ele_px, ele_py, ele_pz;
    Double_t ele_theta_cor, ele_phi_cor, ele_ene_cor, ele_mom_cor, ele_px_cor, ele_py_cor, ele_pz_cor;
    Double_t ele_theta_res, ele_phi_res, ele_ene_res, ele_mom_res, ele_px_res, ele_py_res, ele_pz_res;
    Double_t pim_theta, pim_phi, pim_ene, pim_mom, pim_px, pim_py, pim_pz;
    Double_t pim_theta_cor, pim_phi_cor, pim_ene_cor, pim_mom_cor, pim_px_cor, pim_py_cor, pim_pz_cor;
    Double_t pim_theta_res, pim_phi_res, pim_ene_res, pim_mom_res, pim_px_res, pim_py_res, pim_pz_res;
    Double_t pro_theta, pro_phi, pro_ene, pro_mom, pro_px, pro_py, pro_pz;
    Double_t pro_theta_cor, pro_phi_cor, pro_ene_cor, pro_mom_cor, pro_px_cor, pro_py_cor, pro_pz_cor;
    Double_t pro_theta_res, pro_phi_res, pro_ene_res, pro_mom_res, pro_px_res, pro_py_res, pro_pz_res;
    Double_t EventWeight, WilliamsWeight, DedrickWeight, CatchenWeight;
    Double_t Theta_Pion_Photon, Phi, PhiS, Phi_cor, PhiS_cor, Vertex_X, Vertex_Y, Vertex_Z;
    Double_t Asym_PhiS, Asym_PhiMinusPhiS, Asym_2PhiMinusPhiS, Asym_3PhiMinusPhiS,Asym_PhiPlusPhiS,Asym_2PhiPlusPhiS; 
    Double_t Sigma_PhiS, Sigma_PhiMinusPhiS, Sigma_2PhiMinusPhiS, Sigma_3PhiMinusPhiS, Sigma_PhiPlusPhiS, Sigma_2PhiPlusPhiS;
    Double_t Sigma_Lab, Sigma_UU, Sigma_UT, SSAsym, SineAsym;
    Double_t Sig_L, Sig_T, Sig_LT, Sig_TT;
    Double_t Flux_Factor_RF, Flux_Factor_Col, Jacobian_CM, Jacobian_CM_RF, Jacobian_CM_Col, A_Factor, Photon_Factor, Photon_Theta_cor, Photon_Theta;
    /*}}}*/

    /*Branches{{{*/
    //Target Quantities/*{{{*/
    t1->SetBranchAddress("Target_Theta_Col",                          &tgt_theta     );
    t1->SetBranchAddress("Target_Phi_Col",                            &tgt_phi       );
    t1->SetBranchAddress("Target_Energy_Col_GeV",                     &tgt_ene       );
    t1->SetBranchAddress("Target_Mom_Col_GeV",                        &tgt_mom       );
    t1->SetBranchAddress("Target_MomZ_Col_GeV",                       &tgt_pz        );
    t1->SetBranchAddress("Target_MomX_Col_GeV",                       &tgt_px        );
    t1->SetBranchAddress("Target_MomY_Col_GeV",                       &tgt_py        );/*}}}*/

    //Beam Electron Quantities/*{{{*/
    t1->SetBranchAddress("Electron_Theta_Col",                        &beam_theta      );
    t1->SetBranchAddress("Electron_Phi_Col",                          &beam_phi        );
    t1->SetBranchAddress("Electron_Energy_Col_GeV",                   &beam_ene        );
    t1->SetBranchAddress("Electron_Mom_Col_GeV",                      &beam_mom        );
    t1->SetBranchAddress("Electron_MomX_Col_GeV",                     &beam_px         );
    t1->SetBranchAddress("Electron_MomY_Col_GeV",                     &beam_py         );
    t1->SetBranchAddress("Electron_MomZ_Col_GeV",                     &beam_pz         );
    t1->SetBranchAddress("Electron_Corrected_Theta_Col",              &beam_theta_cor );
    t1->SetBranchAddress("Electron_Corrected_Phi_Col",                &beam_phi_cor   );
    t1->SetBranchAddress("Electron_Corrected_Energy_Col_GeV",         &beam_ene_cor   );
    t1->SetBranchAddress("Electron_Corrected_Mom_Col_GeV",            &beam_mom_cor   );
    t1->SetBranchAddress("Electron_Corrected_MomX_Col_GeV",           &beam_px_cor    );
    t1->SetBranchAddress("Electron_Corrected_MomY_Col_GeV",           &beam_py_cor    );
    t1->SetBranchAddress("Electron_Corrected_MomZ_Col_GeV",           &beam_pz_cor    );/*}}}*/

    //Scattered Electron Quantities/*{{{*/
    t1->SetBranchAddress("ScatElec_Theta_Col",                        &ele_theta     );
    t1->SetBranchAddress("ScatElec_Phi_Col",                          &ele_phi       );
    t1->SetBranchAddress("ScatElec_Energy_Col_GeV",                   &ele_ene       );
    t1->SetBranchAddress("ScatElec_Mom_Col_GeV",                      &ele_mom       );
    t1->SetBranchAddress("ScatElec_MomX_Col_GeV",                     &ele_px        );
    t1->SetBranchAddress("ScatElec_MomY_Col_GeV",                     &ele_py        );
    t1->SetBranchAddress("ScatElec_MomZ_Col_GeV",                     &ele_pz        );
    t1->SetBranchAddress("ScatElec_Corrected_Theta_Col",              &ele_theta_cor);
    t1->SetBranchAddress("ScatElec_Corrected_Phi_Col",                &ele_phi_cor  );
    t1->SetBranchAddress("ScatElec_Corrected_Energy_Col_GeV",         &ele_ene_cor  );
    t1->SetBranchAddress("ScatElec_Corrected_Mom_Col_GeV",            &ele_mom_cor  );
    t1->SetBranchAddress("ScatElec_Corrected_MomX_Col_GeV",           &ele_px_cor   );
    t1->SetBranchAddress("ScatElec_Corrected_MomY_Col_GeV",           &ele_py_cor   );
    t1->SetBranchAddress("ScatElec_Corrected_MomZ_Col_GeV",           &ele_pz_cor   );/*}}}*/

    //Pion- Quantities/*{{{*/
    t1->SetBranchAddress("Pion_Theta_Col",                            &pim_theta     );
    t1->SetBranchAddress("Pion_Phi_Col",                              &pim_phi       );
    t1->SetBranchAddress("Pion_Energy_Col_GeV",                       &pim_ene       );
    t1->SetBranchAddress("Pion_Mom_Col_GeV",                          &pim_mom       );
    t1->SetBranchAddress("Pion_MomX_Col_GeV",                         &pim_px        );
    t1->SetBranchAddress("Pion_MomY_Col_GeV",                         &pim_py        );
    t1->SetBranchAddress("Pion_MomZ_Col_GeV",                         &pim_pz        );
    t1->SetBranchAddress("Pion_Corrected_Theta_Col",                  &pim_theta_cor);
    t1->SetBranchAddress("Pion_Corrected_Phi_Col",                    &pim_phi_cor  );
    t1->SetBranchAddress("Pion_Corrected_Energy_Col_GeV",             &pim_ene_cor  );
    t1->SetBranchAddress("Pion_Corrected_Mom_Col_GeV",                &pim_mom_cor  );
    t1->SetBranchAddress("Pion_Corrected_MomX_Col_GeV",               &pim_px_cor   );
    t1->SetBranchAddress("Pion_Corrected_MomY_Col_GeV",               &pim_py_cor   );
    t1->SetBranchAddress("Pion_Corrected_MomZ_Col_GeV",               &pim_pz_cor   );/*}}}*/

    //Recoil Proton Quantities/*{{{*/
    t1->SetBranchAddress("RecoilProton_Theta_Col",                    &pro_theta          );
    t1->SetBranchAddress("RecoilProton_Phi_Col",                      &pro_phi            );
    t1->SetBranchAddress("RecoilProton_Energy_Col_GeV",               &pro_ene            );
    t1->SetBranchAddress("RecoilProton_Mom_Col_GeV",                  &pro_mom            );
    t1->SetBranchAddress("RecoilProton_MomX_Col_GeV",                 &pro_px             );
    t1->SetBranchAddress("RecoilProton_MomY_Col_GeV",                 &pro_py             );
    t1->SetBranchAddress("RecoilProton_MomZ_Col_GeV",                 &pro_pz             );
    t1->SetBranchAddress("RecoilProton_Corrected_Theta_Col",          &pro_theta_cor     );
    t1->SetBranchAddress("RecoilProton_Corrected_Phi_Col",            &pro_phi_cor       );
    t1->SetBranchAddress("RecoilProton_Corrected_Energy_Col_GeV",     &pro_ene_cor       );
    t1->SetBranchAddress("RecoilProton_Corrected_Mom_Col_GeV",        &pro_mom_cor       );
    t1->SetBranchAddress("RecoilProton_Corrected_MomX_Col_GeV",       &pro_px_cor        );
    t1->SetBranchAddress("RecoilProton_Corrected_MomY_Col_GeV",       &pro_py_cor        );
    t1->SetBranchAddress("RecoilProton_Corrected_MomZ_Col_GeV",       &pro_pz_cor        );/*}}}*/

    t1->SetBranchAddress("Phi",                                       &Phi                 );/*{{{*/
    t1->SetBranchAddress("PhiS",                                      &PhiS                );
    t1->SetBranchAddress("Phi_Corrected",                             &Phi_cor            );
    t1->SetBranchAddress("PhiS_Corrected",                            &PhiS_cor           );
    t1->SetBranchAddress("Photon_Theta_Col",                          &Photon_Theta        );
    t1->SetBranchAddress("Photon_Corrected_Theta_Col",                &Photon_Theta_cor        );
    //t1->SetBranchAddress("Factor_Col",                                &Photon_Factor       );
    /*}}}*/

    //Different XS Quantities/
    //In the lab frame/*{{{*/
    t1->SetBranchAddress("ZASigma_Lab",                               &Sigma_Lab         );//Pleaes explain here
    t1->SetBranchAddress("ZASigma_UU_Col",                            &Sigma_UU          );//Pleaes explain here
    t1->SetBranchAddress("RorySigma_UT_Col",                          &Sigma_UT        );//Pleaes explain here

    //XSs corresponding to different Asymmetries
    //dSigma_k_Sin(...) = kFactor * Sin(...) * fAsym_k * dSigma_UU, in lab frame
    //t1->SetBranchAddress("Sig_PhiS_Col",                                 &Sigma_PhiS           );
    //t1->SetBranchAddress("Sig_Phi_Minus_PhiS_Col",                       &Sigma_PhiMinusPhiS );
    //t1->SetBranchAddress("Sig_2Phi_Minus_PhiS_Col",                      &Sigma_2PhiMinusPhiS);
    //t1->SetBranchAddress("Sig_3Phi_Minus_PhiS_Col",                      &Sigma_3PhiMinusPhiS);
    //t1->SetBranchAddress("Sig_Phi_Plus_PhiS_Col",                        &Sigma_PhiPlusPhiS  );
    //t1->SetBranchAddress("Sig_2Phi_Plus_PhiS_Col",                       &Sigma_2PhiPlusPhiS );
    t1->SetBranchAddress("Term_Phi_S_Col",                                 &Sigma_PhiS           );
    t1->SetBranchAddress("Term_PhiMinusPhi_S_Col",                       &Sigma_PhiMinusPhiS );
    t1->SetBranchAddress("Term_2PhiMinusPhi_S_Col",                      &Sigma_2PhiMinusPhiS);
    t1->SetBranchAddress("Term_3PhiMinusPhi_S_Col",                      &Sigma_3PhiMinusPhiS);
    t1->SetBranchAddress("Term_PhiPlusPhi_S_Col",                        &Sigma_PhiPlusPhiS  );
    t1->SetBranchAddress("Term_2PhiPlusPhi_S_Col",                       &Sigma_2PhiPlusPhiS );

    //6 Asymmetries, "_Col" means in the lab frame, w/o that means the rest frame
    t1->SetBranchAddress("SSAsym",                                       &SSAsym           );
    t1->SetBranchAddress("SineAsym",                                     &SineAsym           );
    t1->SetBranchAddress("AsymPhi_S_Col",                                &Asym_PhiS           );
    t1->SetBranchAddress("AsymPhiMinusPhi_S_Col",                        &Asym_PhiMinusPhiS   );
    t1->SetBranchAddress("Asym2PhiMinusPhi_S_Col",                       &Asym_2PhiMinusPhiS  );/*}}}*/
    t1->SetBranchAddress("Asym3PhiMinusPhi_S_Col",                       &Asym_3PhiMinusPhiS  );
    t1->SetBranchAddress("AsymPhiMinusPhi_S_Col",                        &Asym_PhiMinusPhiS  );
    t1->SetBranchAddress("AsymPhiPlusPhi_S_Col",                         &Asym_PhiPlusPhiS    );
    t1->SetBranchAddress("Asym2PhiPlusPhi_S_Col",                        &Asym_2PhiPlusPhiS    );

    //LT XS in the rest frame/*{{{*/
    t1->SetBranchAddress("ZASig_T",                                   &Sig_T             );//Pleaes explain here
    t1->SetBranchAddress("ZASig_L",                                   &Sig_L             );//Pleaes explain here
    t1->SetBranchAddress("ZASig_LT",                                  &Sig_LT            );//Pleaes explain here
    t1->SetBranchAddress("ZASig_TT",                                  &Sig_TT            );/*}}}*/

    //Weights/*{{{*/
    t1->SetBranchAddress("EventWeight",                               &EventWeight   );//if chaining more root files, weight=weight/N_file
    t1->SetBranchAddress("WilliamsWeight",                            &WilliamsWeight);//FSI weight
    t1->SetBranchAddress("DedrickWeight",                             &DedrickWeight );//FSI weight
    t1->SetBranchAddress("CatchenWeight",                             &CatchenWeight );//FSI weight  /*}}}*/

    //Other Quantities/*{{{*/
    t1->SetBranchAddress("NRecorded",                                 &NRecorded            ); 
    t1->SetBranchAddress("NGenerated",                                &NGenerated           );

    t1->SetBranchAddress("W_GeV",                                     &W                    );
    t1->SetBranchAddress("Qsq_GeV",                                   &Qsq                  );
    t1->SetBranchAddress("T_Para_GeV",                                &t_Para               );
    t1->SetBranchAddress("T_GeV",                                     &t                    );
    t1->SetBranchAddress("Epsilon",                                   &Epsilon              );
    t1->SetBranchAddress("x",                                         &x                    );
    t1->SetBranchAddress("y",                                         &y                    );
    t1->SetBranchAddress("z",                                         &z                    );

    t1->SetBranchAddress("T_Corrected_GeV",                           &t_cor               );
    t1->SetBranchAddress("Qsq_Corrected_GeV",                         &Qsq_cor             );
    t1->SetBranchAddress("W_Corrected_GeV",                           &W_cor               );
    t1->SetBranchAddress("x_Corrected",                               &x_cor               );
    t1->SetBranchAddress("y_Corrected",                               &y_cor               );
    t1->SetBranchAddress("z_Corrected",                               &z_cor               );

    t1->SetBranchAddress("Vertex_X",                                  &Vertex_X             );
    t1->SetBranchAddress("Vertex_Y",                                  &Vertex_Y             );
    t1->SetBranchAddress("Vertex_Z",                                  &Vertex_Z             );
    t1->SetBranchAddress("Theta_Pion_Photon_Col",                     &Theta_Pion_Photon    );

    t1->SetBranchAddress("A",                                         &A_Factor              );
    t1->SetBranchAddress("Flux_Factor_RF",                            &Flux_Factor_RF           );
    t1->SetBranchAddress("Flux_Factor_Col",                           &Flux_Factor_Col       );
    t1->SetBranchAddress("Jacobian_CM",                               &Jacobian_CM           );
    t1->SetBranchAddress("Jacobian_CM_RF",                            &Jacobian_CM_RF        );
    t1->SetBranchAddress("Jacobian_CM_Col",                           &Jacobian_CM_Col        );
    /*}}}*/
    /*}}}*/
    /*}}}*/

    if(CASE==1) new_filename = Form("./rootfiles/DEMP_%s_simple.root", targetname.Data());
    if(CASE==2) new_filename = Form("./rootfiles/DEMP_%s_mult.root", targetname.Data());
    if(CASE==3) new_filename = Form("./rootfiles/DEMP_%s_mult_fsi.root", targetname.Data());
    if(CASE==4) new_filename = Form("./rootfiles/DEMP_%s_fermi.root", targetname.Data());
    if(CASE==5) new_filename = Form("./rootfiles/DEMP_%s_mult_nofermi.root", targetname.Data());
    cout<<"---  Saving in file = "<<new_filename.Data()<<endl;
    TFile* f2=new TFile(new_filename.Data(),"recreate"); 
   
    /*Define new Tree and new Branch{{{*/
    double ele_acc_f, ele_acc_l, pim_acc_f,pim_acc_l, pro_acc_f,pro_acc_l,total_acc, total_acc_cor, total_acc_res;
    double tp,tp_cor, Lumi_PSF, MM, MM_cor, MM_res,MP, MP_cor, MP_res,dilute;
    double weight, weight_uu, weight_ut, weight_3m1, weight_2m1, weight_1m1, weight_0p1, weight_1p1, weight_2p1;
    //double ele_mom_res, ele_theta_res, ele_phi_res;
    //double pim_mom_res, pim_theta_res, pim_phi_res;
    //double pro_mom_res, pro_theta_res, pro_phi_res;

    //TTree *t2 = t1->CloneTree(0);
    TTree *t2 = new TTree("T","a new tree");

    t2->Branch("NRecorded",       &NRecorded,     "NRecorded/I");/*{{{*/
    t2->Branch("NGenerated",       &NGenerated,     "NGenerated/I");

    t2->Branch("Epsilon", &Epsilon, "Epsilon/D" );
    t2->Branch("Qsq", &Qsq ,"Qsq/D");
    t2->Branch("t", &t ,"t/D");
    t2->Branch("tp", &tp ,"tp/D");
    t2->Branch("t_Para", &t_Para ,"t_Para/D");
    t2->Branch("W", &W ,"W/D");
    t2->Branch("x", &x ,"x/D");
    t2->Branch("y", &y ,"y/D");
    t2->Branch("z", &z ,"z/D");

    t2->Branch("Qsq_cor", &Qsq_cor ,"Qsq_cor/D");
    t2->Branch("t_cor", &t_cor ,"t_cor/D");
    t2->Branch("tp_cor", &tp_cor ,"tp_cor/D");
    t2->Branch("W_cor", &W_cor ,"W_cor/D");
    t2->Branch("x_cor", &x_cor ,"x_cor/D");
    t2->Branch("y_cor", &y_cor ,"y_cor/D");
    t2->Branch("z_cor", &z_cor ,"z_cor/D"); 

    t2->Branch("Vertex_X",   &Vertex_X   ,"Vertex_X/D");
    t2->Branch("Vertex_Y",   &Vertex_Y   ,"Vertex_Y/D");
    t2->Branch("Vertex_Z",   &Vertex_Z   ,"Vertex_Z/D");
    t2->Branch("Theta_Pion_Photon",   &Theta_Pion_Photon   ,"Theta_Pion_Photon/D");

    t2->Branch("A_Factor",                                  &A_Factor          ,"data/D"    );
    t2->Branch("Flux_Factor_RF",                            &Flux_Factor_RF       ,"data/D"    );
    t2->Branch("Flux_Factor_Col",                           &Flux_Factor_Col   ,"data/D"    );
    t2->Branch("Jacobian_CM",                               &Jacobian_CM       ,"data/D"    );
    t2->Branch("Jacobian_CM_RF",                            &Jacobian_CM_RF    ,"data/D"    );
    t2->Branch("Jacobian_CM_Col",                           &Jacobian_CM_Col   ,"data/D"    );
    /*}}}*/

    t2->Branch("Phi",       &Phi,       "Phi/D");/*{{{*/
    t2->Branch("PhiS",      &PhiS,      "PhiS/D");
    t2->Branch("Phi_cor",  &Phi_cor,  "Phi_cor/D");
    t2->Branch("PhiS_cor", &PhiS_cor, "PhiS_cor/D");
    t2->Branch("Photon_Theta",      &Photon_Theta,      "Photon_Theta/D");
    t2->Branch("Photon_Theta_cor",      &Photon_Theta_cor,      "Photon_Theta_cor/D");
    t2->Branch("Photon_Factor",      &Photon_Factor,      "Photon_Factor/D");
    /*}}}*/

    t2->Branch("Sigma_Lab",     &Sigma_Lab, "data/D");                              /*{{{*/
    t2->Branch("Sigma_UU",      &Sigma_UU,   "data/D");                              
    t2->Branch("Sigma_UT",      &Sigma_UT,   "data/D");                              
    t2->Branch("Sig_T",         &Sig_T, "data/D");                                  
    t2->Branch("Sig_L",         &Sig_L, "data/D");                                  
    t2->Branch("Sig_LT",        &Sig_LT, "data/D");                                 
    t2->Branch("Sig_TT",        &Sig_TT, "data/D");                                      /*}}}*/

    //Six Asymmetries and polarized XS/*{{{*/
    t2->Branch("SSAsym",              &SSAsym, "data/D");
    t2->Branch("SineAsym",            &SineAsym, "data/D");
    t2->Branch("Asym_PhiS",           &Asym_PhiS, "data/D");
    t2->Branch("Asym_PhiPlusPhiS",    &Asym_PhiPlusPhiS, "data/D");
    t2->Branch("Asym_2PhiPlusPhiS",   &Asym_2PhiPlusPhiS, "data/D");
    t2->Branch("Asym_PhiMinusPhiS",   &Asym_PhiMinusPhiS, "data/D");
    t2->Branch("Asym_2PhiMinusPhiS",  &Asym_2PhiMinusPhiS, "data/D");
    t2->Branch("Asym_3PhiMinusPhiS",  &Asym_3PhiMinusPhiS, "data/D");

    t2->Branch("Sigma_PhiS",           &Sigma_PhiS, "data/D");
    t2->Branch("Sigma_PhiPlusPhiS",    &Sigma_PhiPlusPhiS, "data/D");
    t2->Branch("Sigma_2PhiPlusPhiS",   &Sigma_2PhiPlusPhiS, "data/D");
    t2->Branch("Sigma_PhiMinusPhiS",   &Sigma_PhiMinusPhiS, "data/D");
    t2->Branch("Sigma_2PhiMinusPhiS",  &Sigma_2PhiMinusPhiS, "data/D");
    t2->Branch("Sigma_3PhiMinusPhiS",  &Sigma_3PhiMinusPhiS, "data/D");/*}}}*/

    t2->Branch("EventWeight",         &EventWeight, "data/D");                              /*{{{*/
    t2->Branch("WilliamsWeight",      &WilliamsWeight, "data/D");                              
    t2->Branch("DedrickWeight",       &DedrickWeight, "data/D");                              
    t2->Branch("CatchenWeight",       &CatchenWeight, "data/D");                              /*}}}*/

    t2->Branch("tgt_px", &tgt_px, "tgt_px/D");/*{{{*/
    t2->Branch("tgt_py", &tgt_py, "tgt_py/D");
    t2->Branch("tgt_pz", &tgt_pz, "tgt_pz/D");
    t2->Branch("tgt_theta", &tgt_theta, "tgt_theta/D");
    t2->Branch("tgt_phi", &tgt_phi, "tgt_phi/D");
    t2->Branch("tgt_ene", &tgt_ene, "tgt_ene/D");
    t2->Branch("tgt_mom", &tgt_mom, "tgt_mom/D");/*}}}*/

    t2->Branch("beam_px", &beam_px, "beam_px/D");/*{{{*/
    t2->Branch("beam_py", &beam_py, "beam_py/D");
    t2->Branch("beam_pz", &beam_pz, "beam_pz/D");
    t2->Branch("beam_theta", &beam_theta, "beam_theta/D");
    t2->Branch("beam_phi", &beam_phi, "beam_phi/D");
    t2->Branch("beam_ene", &beam_ene, "beam_ene/D");
    t2->Branch("beam_mom", &beam_mom, "beam_mom/D");

    t2->Branch("beam_px_cor", &beam_px_cor, "beam_px_cor/D");
    t2->Branch("beam_py_cor", &beam_py_cor, "beam_py_cor/D");
    t2->Branch("beam_pz_cor", &beam_pz_cor, "beam_pz_cor/D");
    t2->Branch("beam_theta_cor", &beam_theta_cor, "beam_theta_cor/D");
    t2->Branch("beam_phi_cor", &beam_phi_cor, "beam_phi_cor/D");
    t2->Branch("beam_ene_cor", &beam_ene_cor, "beam_ene_cor/D");
    t2->Branch("beam_mom_cor", &beam_mom_cor, "beam_mom_cor/D");/*}}}*/

    t2->Branch("ele_px",    &ele_px    ,"ele_px/D");/*{{{*/
    t2->Branch("ele_py",    &ele_py    ,"ele_py/D");
    t2->Branch("ele_pz",    &ele_pz    ,"ele_pz/D");
    t2->Branch("ele_mom",   &ele_mom   ,"ele_mom/D");
    t2->Branch("ele_ene",   &ele_ene   ,"ele_ene/D");
    t2->Branch("ele_theta", &ele_theta ,"ele_theta/D");
    t2->Branch("ele_phi",   &ele_phi   ,"ele_phi/D");

    t2->Branch("ele_px_cor",    &ele_px_cor    ,"ele_px_cor/D");
    t2->Branch("ele_py_cor",    &ele_py_cor    ,"ele_py_cor/D");
    t2->Branch("ele_pz_cor",    &ele_pz_cor    ,"ele_pz_cor/D");
    t2->Branch("ele_mom_cor",   &ele_mom_cor   ,"ele_mom_cor/D");
    t2->Branch("ele_ene_cor",   &ele_ene_cor   ,"ele_ene_cor/D");
    t2->Branch("ele_theta_cor", &ele_theta_cor ,"ele_theta_cor/D");
    t2->Branch("ele_phi_cor",   &ele_phi_cor   ,"ele_phi_cor/D");
   
    t2->Branch("ele_ene_res",   &ele_ene_res   ,"ele_ene_res/D");
    t2->Branch("ele_px_res",    &ele_px_res    ,"ele_px_res/D");
    t2->Branch("ele_py_res",    &ele_py_res    ,"ele_py_res/D");
    t2->Branch("ele_pz_res",    &ele_pz_res    ,"ele_pz_res/D");
    t2->Branch("ele_mom_res",   &ele_mom_res   ,"ele_mom_res/D");
    t2->Branch("ele_theta_res", &ele_theta_res ,"ele_theta_res/D");
    t2->Branch("ele_phi_res",   &ele_phi_res   ,"ele_phi_res/D");/*}}}*/
    
    t2->Branch("pim_ene",   &pim_ene   ,"pim_ene/D");/*{{{*/
    t2->Branch("pim_px",    &pim_px    ,"pim_px/D");
    t2->Branch("pim_py",    &pim_py    ,"pim_py/D");
    t2->Branch("pim_pz",    &pim_pz    ,"pim_pz/D");
    t2->Branch("pim_mom",   &pim_mom   ,"pim_mom/D");
    t2->Branch("pim_theta", &pim_theta ,"pim_theta/D");
    t2->Branch("pim_phi",   &pim_phi   ,"pim_phi/D");

    t2->Branch("pim_ene_cor",   &pim_ene_cor   ,"pim_ene_cor/D");
    t2->Branch("pim_px_cor",    &pim_px_cor    ,"pim_px_cor/D");
    t2->Branch("pim_py_cor",    &pim_py_cor    ,"pim_py_cor/D");
    t2->Branch("pim_pz_cor",    &pim_pz_cor    ,"pim_pz_cor/D");
    t2->Branch("pim_mom_cor",   &pim_mom_cor   ,"pim_mom_cor/D");
    t2->Branch("pim_theta_cor", &pim_theta_cor ,"pim_theta_cor/D");
    t2->Branch("pim_phi_cor",   &pim_phi_cor   ,"pim_phi_cor/D");
 
    t2->Branch("pim_ene_res",   &pim_ene_res   ,"pim_ene_res/D");
    t2->Branch("pim_px_res",    &pim_px_res    ,"pim_px_res/D");
    t2->Branch("pim_py_res",    &pim_py_res    ,"pim_py_res/D");
    t2->Branch("pim_pz_res",    &pim_pz_res    ,"pim_pz_res/D");
    t2->Branch("pim_mom_res",   &pim_mom_res   ,"pim_mom_res/D");
    t2->Branch("pim_theta_res", &pim_theta_res ,"pim_theta_res/D");
    t2->Branch("pim_phi_res",   &pim_phi_res   ,"pim_phi_res/D");
    /*}}}*/

    t2->Branch("pro_ene",   &pro_ene   ,"pro_ene/D"); /*{{{*/
    t2->Branch("pro_px",    &pro_px    ,"pro_px/D");
    t2->Branch("pro_py",    &pro_py    ,"pro_py/D");
    t2->Branch("pro_pz",    &pro_pz    ,"pro_pz/D");
    t2->Branch("pro_mom",   &pro_mom   ,"pro_mom/D");
    t2->Branch("pro_theta", &pro_theta ,"pro_theta/D");
    t2->Branch("pro_phi",   &pro_phi   ,"pro_phi/D");

    t2->Branch("pro_ene_cor",   &pro_ene_cor   ,"pro_ene_cor/D");
    t2->Branch("pro_px_cor",    &pro_px_cor    ,"pro_px_cor/D");
    t2->Branch("pro_py_cor",    &pro_py_cor    ,"pro_py_cor/D");
    t2->Branch("pro_pz_cor",    &pro_pz_cor    ,"pro_pz_cor/D");
    t2->Branch("pro_mom_cor",   &pro_mom_cor   ,"pro_mom_cor/D");
    t2->Branch("pro_theta_cor", &pro_theta_cor ,"pro_theta_cor/D");
    t2->Branch("pro_phi_cor",   &pro_phi_cor   ,"pro_phi_cor/D");

    t2->Branch("pro_ene_res",   &pro_ene_res   ,"pro_ene_res/D");
    t2->Branch("pro_px_res",    &pro_px_res    ,"pro_px_res/D");
    t2->Branch("pro_py_res",    &pro_py_res    ,"pro_py_res/D");
    t2->Branch("pro_pz_res",    &pro_pz_res    ,"pro_pz_res/D");
    t2->Branch("pro_mom_res",   &pro_mom_res   ,"pro_mom_res/D");
    t2->Branch("pro_theta_res", &pro_theta_res ,"pro_theta_res/D");
    t2->Branch("pro_phi_res",   &pro_phi_res   ,"pro_phi_res/D");
    /*}}}*/

    //Add SoLID acceptance/*{{{*/
    t2->Branch("ele_acc_f",     &ele_acc_f,     "ele_acc_f/D");
    t2->Branch("ele_acc_l",     &ele_acc_l,     "ele_acc_l/D");
    t2->Branch("pim_acc_f",     &pim_acc_f,     "pim_acc_f/D");
    t2->Branch("pim_acc_l",     &pim_acc_l,     "pim_acc_l/D");
    t2->Branch("pro_acc_f",     &pro_acc_f,     "pro_acc_f/D");
    t2->Branch("pro_acc_l",     &pro_acc_l,     "pro_acc_l/D");/*}}}*/

    //Add other quantities/*{{{*/
    t2->Branch("weight",        &weight,        "weight/D");
    t2->Branch("weight_uu",        &weight_uu,        "weight_uu/D"); //weight for unpolarized XS
    t2->Branch("weight_ut",        &weight_ut,        "weight_ut/D"); //weight for polarized XS
    t2->Branch("weight_3m1",        &weight_3m1,        "weight_3m1/D"); //weight for Sin(3Phi-PhiS) module
    t2->Branch("weight_2m1",        &weight_2m1,        "weight_2m1/D"); //weight for Sin(2Phi-PhiS) module
    t2->Branch("weight_1m1",        &weight_1m1,        "weight_1m1/D"); //weight for Sin(Phi-PhiS) module
    t2->Branch("weight_0p1",        &weight_0p1,        "weight_0p1/D"); //weight for Sin(PhiS) module
    t2->Branch("weight_1p1",        &weight_1p1,        "weight_1p1/D"); //weight for Sin(Phi+PhiS) module
    t2->Branch("weight_2p1",        &weight_2p1,        "weight_2p1/D"); //weight for Sin(2Phi+PhiS) module
    t2->Branch("dilute",        &dilute,        "dilute/D");
    //t2->Branch("PSF",           &PSF,           "PSF/D");
    t2->Branch("MM",            &MM,            "MM/D");
    t2->Branch("MM_cor",        &MM_cor,        "MM_cor/D");
    t2->Branch("MM_res",        &MM_res,        "MM_res/D");
    t2->Branch("MP",            &MP,            "MP/D");
    t2->Branch("MP_cor",        &MP_cor,        "MP_cor/D");
    t2->Branch("MP_res",        &MP_res,        "MP_res/D");
    /*}}}*/
    
    t2->Branch("Lumi_PSF",            &Lumi_PSF,            "Lumi_PSF/D");
    t2->Branch("time",            &time,            "time/D");
    t2->Branch("fileNO",          &fileNO,          "fileNO/I");
    t2->Branch("total_acc",       &total_acc,       "total_acc/D");
    t2->Branch("total_acc_cor",  &total_acc_cor,  "total_acc_cor/D");
    t2->Branch("total_acc_res",   &total_acc_res,   "total_acc_res/D");
    //t2->Branch("",            &,            "/D");
    /*}}}*/

    /*Vectors & Histograms{{{*/
    new_filename.ReplaceAll("DEMP_", "histo_");
    TFile* histo =new TFile(new_filename.Data(),"recreate"); 
    TLorentzVector *P_E0 = new TLorentzVector();//incoming electron
    TLorentzVector *P_e = new TLorentzVector();//scattered electron with Eloss
    TLorentzVector *P_pim = new TLorentzVector();//photon
    TLorentzVector *P_pro = new TLorentzVector();//proton or neutron
    TLorentzVector *P_t = new TLorentzVector();//target, either proton or neutron
    TLorentzVector *MV = new TLorentzVector();//scattered electron with Eloss

    TH1F *hMM_dvmp = new TH1F("hMM_dvmp",Form("DVMP Missing Mass at %3.1f GeV",EBeam), 200,0.0, 2.5);/*{{{*/
    hMM_dvmp->SetXTitle("Hadron Missing Mass (GeV)");
    hMM_dvmp->GetXaxis()->CenterTitle(1);
    hMM_dvmp->SetYTitle("Rate (Hz)");
    hMM_dvmp->GetYaxis()->CenterTitle(1);
    TH1F *hMM_dvmp_cor = new TH1F("hMM_dvmp_cor",Form("DVMP Missing Mass at %3.1f GeV",EBeam), 200,0.0, 2.5);
    hMM_dvmp_cor->SetXTitle("Hadron Missing Mass (GeV)");
    hMM_dvmp_cor->GetXaxis()->CenterTitle(1);
    hMM_dvmp_cor->SetYTitle("Rate (Hz)");
    hMM_dvmp_cor->GetYaxis()->CenterTitle(1);
    TH1F *hMM_dvmp_res = new TH1F("hMM_dvmp_res",Form("DVMP Missing Mass at %3.1f GeV",EBeam), 200,0.0, 2.5);
    hMM_dvmp_res->SetXTitle("Hadron Missing Mass (GeV)");
    hMM_dvmp_res->GetXaxis()->CenterTitle(1);
    hMM_dvmp_res->SetYTitle("Rate (Hz)");
    hMM_dvmp_res->GetYaxis()->CenterTitle(1);
    /*}}}*/

    TH1F *hMM_dvmp_cut = new TH1F("hMM_dvmp_cut",Form("DVMP Missing Mass at %3.1f GeV",EBeam), 200,0.0, 2.5);/*{{{*/
    hMM_dvmp_cut->SetXTitle("Hadron Missing Mass (GeV)");
    hMM_dvmp_cut->GetXaxis()->CenterTitle(1);
    hMM_dvmp_cut->SetYTitle("Rate (Hz)");
    hMM_dvmp_cut->GetYaxis()->CenterTitle(1);
    TH1F *hMM_dvmp_cut_res = new TH1F("hMM_dvmp_cut_res",Form("DVMP Missing Mass at %3.1f GeV",EBeam), 200,0.0, 2.5);
    hMM_dvmp_cut_res->SetXTitle("Hadron Missing Mass (GeV)");
    hMM_dvmp_cut_res->GetXaxis()->CenterTitle(1);
    hMM_dvmp_cut_res->SetYTitle("Rate (Hz)");
    hMM_dvmp_cut_res->GetYaxis()->CenterTitle(1);
    TH1F *hMM_dvmp_cut_cor = new TH1F("hMM_dvmp_cut_cor",Form("DVMP Missing Mass at %3.1f GeV",EBeam), 200,0.0, 2.5);
    hMM_dvmp_cut_cor->SetXTitle("Hadron Missing Mass (GeV)");
    hMM_dvmp_cut_cor->GetXaxis()->CenterTitle(1);
    hMM_dvmp_cut_cor->SetYTitle("Rate (Hz)");
    hMM_dvmp_cut_cor->GetYaxis()->CenterTitle(1);
    /*}}}*/

    TH1F *hMP_dvmp = new TH1F("hMP_dvmp",Form("DVMP Missing Momentum at %3.1f GeV",EBeam), 200,0.0, 2.5);/*{{{*/
    hMP_dvmp->SetXTitle("Hadron Missing Momentum (GeV/c)");
    hMP_dvmp->GetXaxis()->CenterTitle(1);
    hMP_dvmp->SetYTitle("Rate (Hz)");
    hMP_dvmp->GetYaxis()->CenterTitle(1);
    TH1F *hMP_dvmp_res = new TH1F("hMP_dvmp_res",Form("DVMP Missing Momentum at %3.1f GeV",EBeam), 200,0.0, 2.5);
    hMP_dvmp_res->SetXTitle("Hadron Missing Momentum (GeV/c)");
    hMP_dvmp_res->GetXaxis()->CenterTitle(1);
    hMP_dvmp_res->SetYTitle("Rate (Hz)");
    hMP_dvmp_res->GetYaxis()->CenterTitle(1);
    TH1F *hMP_dvmp_cor = new TH1F("hMP_dvmp_cor",Form("DVMP Missing Momentum at %3.1f GeV",EBeam), 200,0.0, 2.5);
    hMP_dvmp_cor->SetXTitle("Hadron Missing Momentum (GeV/c)");
    hMP_dvmp_cor->GetXaxis()->CenterTitle(1);
    hMP_dvmp_cor->SetYTitle("Rate (Hz)");
    hMP_dvmp_cor->GetYaxis()->CenterTitle(1);
    /*}}}*/

    TH1F *hMT_dvmp = new TH1F("hMT_dvmp",Form("DVMP Missing Angle at %3.1f GeV",EBeam), 200,0.0, 65.);/*{{{*/
    hMT_dvmp->SetXTitle("Hadron Missing Angle (Degree)");
    hMT_dvmp->GetXaxis()->CenterTitle(1);
    hMT_dvmp->SetYTitle("Rate (Hz)");
    hMT_dvmp->GetYaxis()->CenterTitle(1);
    TH1F *hMT_dvmp_res = new TH1F("hMT_dvmp_res",Form("DVMP Missing Angle at %3.1f GeV",EBeam), 200,0.0, 65.);
    hMT_dvmp_res->SetXTitle("Hadron Missing Angle (Degree)");
    hMT_dvmp_res->GetXaxis()->CenterTitle(1);
    hMT_dvmp_res->SetYTitle("Rate (Hz)");
    hMT_dvmp_res->GetYaxis()->CenterTitle(1);
    TH1F *hMT_dvmp_cor = new TH1F("hMT_dvmp_cor",Form("DVMP Missing Angle at %3.1f GeV",EBeam), 200,0.0, 65.);
    hMT_dvmp_cor->SetXTitle("Hadron Missing Angle (Degree)");
    hMT_dvmp_cor->GetXaxis()->CenterTitle(1);
    hMT_dvmp_cor->SetYTitle("Rate (Hz)");
    hMT_dvmp_cor->GetYaxis()->CenterTitle(1);
    /*}}}*/
    /*}}}*/

    /*Fill in DVMP Missing Mass{{{*/
    double total_rate_dvmp=0.0; //Make sure the value is right
    double total_rate_dvmp_4GeV=0.0; //Make sure the value is right
    for(Long64_t i=0;i<N_entries;i++){
        t1->GetEntry(i);
        //if(ele_theta>7.5&&ele_theta<24.5&&ele_mom>1.&&ele_mom<11
        //&&pim_theta>7.5&&pim_theta<24.5&&pim_mom>1.&&pim_mom<11
        //&& W>=2.&&Qsq>1.0&&SigmaPara>1e-33){//any additional cuts should be added in here
        Lumi_PSF = EventWeight / Sigma_Lab;
        dilute = Epsilon*Sig_L/(Epsilon*Sig_L + Sig_T);
        tp = t-t_Para;
        tp_cor = t_cor-t_Para;
        //weight = EventWeight * total_acc_res / fileNO * time;//XS*Lumi*PSF*Acc, give the real rates

        ///////////////////////////////////////////////////////////////////////////
        //Consider the detector resolution here*{{{*/
        //make sure to use the corrected quantitie for multiple scattering and eloss
        //Electron/*{{{*/
        ele_mom_res = gRandom->Gaus(ele_mom_cor, Sigma_Dp_E*ele_mom_cor);//GeV, for electron, E ~= P
        ele_theta_res = gRandom->Gaus(ele_theta_cor*Deg2Rad, Sigma_Theta_E);//rad
        ele_phi_res = gRandom->Gaus(ele_phi_cor*Deg2Rad, Sigma_Phi_E);//rad

        ele_px_res = ele_mom_res * sin(ele_theta_res)*cos(ele_phi_res); 
        ele_py_res = ele_mom_res * sin(ele_theta_res)*sin(ele_phi_res); 
        ele_pz_res = ele_mom_res * cos(ele_theta_res);
        ele_ene_res = sqrt(ele_mom_res*ele_mom_res + eMass*eMass);
        ele_phi_res *=Rad2Deg; ele_theta_res *= Rad2Deg;
        /*}}}*/

        //Pion/*{{{*/
        pim_mom_res = gRandom->Gaus(pim_mom_cor, Sigma_Dp_E*pim_mom_cor);//GeV, for electron, E ~= P
        pim_theta_res = gRandom->Gaus(pim_theta_cor*Deg2Rad, Sigma_Theta_E);//rad
        pim_phi_res = gRandom->Gaus(pim_phi_cor*Deg2Rad, Sigma_Phi_E);//rad

        pim_px_res = pim_mom_res * sin(pim_theta_res)*cos(pim_phi_res); 
        pim_py_res = pim_mom_res * sin(pim_theta_res)*sin(pim_phi_res); 
        pim_pz_res = pim_mom_res * cos(pim_theta_res);
        pim_ene_res = sqrt(pim_mom_res*pim_mom_res + piMass*piMass);
        pim_phi_res *=Rad2Deg; pim_theta_res *= Rad2Deg;
        /*}}}*/

        //Proton/*{{{*/
        pro_mom_res = gRandom->Gaus(pro_mom_cor, Sigma_Dp_E*pro_mom_cor);//GeV, for electron, E ~= P
        pro_theta_res = gRandom->Gaus(pro_theta_cor*Deg2Rad, Sigma_Theta_E);//rad
        pro_phi_res = gRandom->Gaus(pro_phi_cor*Deg2Rad, Sigma_Phi_E);//rad

        pro_px_res = pro_mom_res * sin(pro_theta_res)*cos(pro_phi_res); 
        pro_py_res = pro_mom_res * sin(pro_theta_res)*sin(pro_phi_res); 
        pro_pz_res = pro_mom_res * cos(pro_theta_res);
        pro_ene_res = sqrt(pro_mom_res*pro_mom_res + piMass*piMass);
        pro_phi_res *=Rad2Deg; pro_theta_res *= Rad2Deg;
        /*}}}*/
        /*}}}*/

        /*Get acceptance of e and pi-{{{*/
        //Make sure to use the corrected quantities for multipile scattering and eloss effects
        //Do not use the smeared quantities since we are about whether particles are in the accepntace or not, but not how good we measure
        /*Elec Acc {{{*/
        ele_acc_f = accpt->GetAcc("e-","forward", ele_mom, ele_theta);
        ele_acc_l = accpt->GetAcc("e-","large", ele_mom, ele_theta);
        if(ele_mom<1.0||ele_theta>14.8||ele_theta<8.0)//GeV, CLEO
            ele_acc_f=0.0;//Farward-Angle EC Cut at 1 GeV
        if(ele_mom<3.5||ele_theta<16.0||ele_theta>24)//GeV,CLEO
            ele_acc_l=0.0; //Larger-Angle EC Cut at 3 GeV
        if(ele_acc_f>1.) 
            ele_acc_f=1.0; 
        if(ele_acc_l>1.) 
            ele_acc_l=1.0; 
        total_acc = (ele_acc_l+ele_acc_f);

        ele_acc_f = accpt->GetAcc("e-","forward", ele_mom_cor, ele_theta_cor);
        ele_acc_l = accpt->GetAcc("e-","large", ele_mom_cor, ele_theta_cor);
        if(ele_mom_cor<1.0||ele_theta_cor>14.8||ele_theta_cor<8.0)//GeV, CLEO
            ele_acc_f=0.0;//Farward-Angle EC Cut at 1 GeV
        if(ele_mom_cor<3.5||ele_theta_cor<16.0||ele_theta_cor>24)//GeV,CLEO
            ele_acc_l=0.0; //Larger-Angle EC Cut at 3 GeV
        if(ele_acc_f>1.) 
            ele_acc_f=1.0; 
        if(ele_acc_l>1.) 
            ele_acc_l=1.0; 
        total_acc_cor = (ele_acc_l+ele_acc_f);

        ele_acc_f = accpt->GetAcc("e-","forward", ele_mom_res, ele_theta_res);
        ele_acc_l = accpt->GetAcc("e-","large", ele_mom_res, ele_theta_res);
        if(ele_mom_res<1.0||ele_theta_res>14.8||ele_theta_res<8.0)//GeV, CLEO
            ele_acc_f=0.0;//Farward-Angle EC Cut at 1 GeV
        if(ele_mom_res<3.5||ele_theta_res<16.0||ele_theta_res>24)//GeV,CLEO
            ele_acc_l=0.0; //Larger-Angle EC Cut at 3 GeV
        if(ele_acc_f>1.) 
            ele_acc_f=1.0; 
        if(ele_acc_l>1.) 
            ele_acc_l=1.0; 
        total_acc_res = (ele_acc_l+ele_acc_f);
        /*}}}*/

        /*Pion Acc{{{*/
        pim_acc_f = accpt->GetAcc("pi-","forward", pim_mom, pim_theta);
        pim_acc_l = accpt->GetAcc("pi-","large", pim_mom, pim_theta);
        if(pim_theta>14.8||pim_theta<8.0||pim_mom<0.||pim_mom>11.)//GeV, CLEO
            pim_acc_f=0.0;
        if(pim_theta<16.0||pim_theta>24.0||pim_mom<0.||pim_mom>11.)//GeV, CLEO
            pim_acc_l=0.0; 
        if(pim_acc_f>1.) 
            pim_acc_f=1.0; 
        if(pim_acc_l>1.) 
            pim_acc_l=1.0; 
        total_acc *= (pim_acc_f);

        pim_acc_f = accpt->GetAcc("pi-","forward", pim_mom_cor, pim_theta_cor);
        pim_acc_l = accpt->GetAcc("pi-","large", pim_mom_cor, pim_theta_cor);
        if(pim_theta_cor>14.8||pim_theta_cor<8.0||pim_mom_cor<0.||pim_mom_cor>11.)//GeV, CLEO
            pim_acc_f=0.0;
        if(pim_theta_cor<16.0||pim_theta_cor>24.0||pim_mom_cor<0.||pim_mom_cor>11.)//GeV, CLEO
            pim_acc_l=0.0; 
        if(pim_acc_f>1.) 
            pim_acc_f=1.0; 
        if(pim_acc_l>1.) 
            pim_acc_l=1.0; 
        total_acc_cor *= (pim_acc_f);

        pim_acc_f = accpt->GetAcc("pi-","forward", pim_mom_res, pim_theta_res);
        pim_acc_l = accpt->GetAcc("pi-","large", pim_mom_res, pim_theta_res);
        if(pim_theta_res>14.8||pim_theta<8.0||pim_mom_res<0.||pim_mom_res>11.)//GeV, CLEO
            pim_acc_f=0.0;
        if(pim_theta_res<16.0||pim_theta_res>24.0||pim_mom_res<0.||pim_mom_res>11.)//GeV, CLEO
            pim_acc_l=0.0; 
        if(pim_acc_f>1.) 
            pim_acc_f=1.0; 
        if(pim_acc_l>1.) 
            pim_acc_l=1.0; 
        total_acc_res *= (pim_acc_f);

        /*}}}*/

        /*Proton Acc {{{*/
        //The momentum cut is applied by EC while for proton, we don't reply on EC to tell the acceptance.
        //What I assume here is that we can detecto all energy range of protons,
        //unlike electrons which need to be separated from pions
        //   pro_acc_f = accpt->GetAcc("p","forward", pro_mom, pro_theta);
        //   pro_acc_l = accpt->GetAcc("p","large", pro_mom, pro_theta);
        //   pro_acc_f = accpt->GetThetaAcc("p","forward", pro_theta);
        //   pro_acc_l = accpt->GetThetaAcc("p","large", pro_theta);

        pro_acc_f = 1.0;
        pro_acc_l = 1.0;
        if(pro_theta>14.8||pro_theta<8.0||pro_mom<0.||pro_mom>11.)//GeV, CLEO
            pro_acc_f=0.0;
        if(pro_theta<16.0||pro_theta>24.0||pro_mom<0.||pro_mom>11.)//GeV, CLEO
            pro_acc_l=0.0; 
        if(pro_acc_f>1.) 
            pro_acc_f=1.0; 
        if(pro_acc_l>1.) 
            pro_acc_l=1.0; 
        total_acc *= (pro_acc_l+pro_acc_f);

        pro_acc_f = 1.0;
        pro_acc_l = 1.0;
        if(pro_theta_cor>14.8||pro_theta_cor<8.0||pro_mom_cor<0.||pro_mom_cor>11.)//GeV, CLEO
            pro_acc_f=0.0;
        if(pro_theta_cor<16.0||pro_theta_cor>24.0||pro_mom_cor<0.||pro_mom_cor>11.)//GeV, CLEO
            pro_acc_l=0.0; 
        if(pro_acc_f>1.) 
            pro_acc_f=1.0; 
        if(pro_acc_l>1.) 
            pro_acc_l=1.0; 
        total_acc_cor *= (pro_acc_l+pro_acc_f);

        pro_acc_f = 1.0;
        pro_acc_l = 1.0;
        if(pro_theta_res>14.8||pro_theta_res<8.0 ||pro_mom_res<0.||pro_mom_res>11.)//GeV, CLEO
            pro_acc_f=0.0;
        if(pro_theta_res<16.0||pro_theta_res>24.0||pro_mom_res<0.||pro_mom_res>11.)//GeV, CLEO
            pro_acc_l=0.0; 
        if(pro_acc_f>1.) 
            pro_acc_f=1.0; 
        if(pro_acc_l>1.) 
            pro_acc_l=1.0; 
        total_acc_res *= (pro_acc_l+pro_acc_f);
        /*}}}*/

        /*}}}*/

        /*define weights{{{*/
        //Didn't save in the Root file so have to reconstruct it back
        /*double Photon_Factor_temp = Sigma_PhiS / sin(PhiS*Deg2Rad) / (Asym_PhiS*Sigma_UU) / 0.865;*/
        /*double Photon_Theta_temp = asin( sqrt( (1-1./Photon_Factor)/pow(sin(PhiS*Deg2Rad),2) ) )*Rad2Deg;*/
        if(targetname=="up")
            Photon_Factor =  1.0 ;// sqrt(1.0 - pow(sin(Photon_Theta*Deg2Rad) ,2) * pow(PhiS*Deg2Rad ,2));
        if(targetname=="down")
            Photon_Factor = -1.0 ;// sqrt(1.0 - pow(sin(Photon_Theta*Deg2Rad) ,2) * pow(PhiS*Deg2Rad ,2));

        //Remove the polarization values from the polarized XS, assuming 100% polarization
        //Note that the angular module sin(m*Phi_cor + l*PhiS_cor) are using the corrected azimuthal angles
        Sigma_3PhiMinusPhiS = Photon_Factor * Sigma_UU * Asym_3PhiMinusPhiS * sin(3.*Phi_cor*Deg2Rad - PhiS_cor*Deg2Rad);
        Sigma_2PhiMinusPhiS = Photon_Factor * Sigma_UU * Asym_2PhiMinusPhiS * sin(2.*Phi_cor*Deg2Rad - PhiS_cor*Deg2Rad);
        Sigma_PhiMinusPhiS  = Photon_Factor * Sigma_UU * Asym_PhiMinusPhiS  * sin(1.*Phi_cor*Deg2Rad - PhiS_cor*Deg2Rad);
        Sigma_PhiS          = Photon_Factor * Sigma_UU * Asym_PhiS          * sin(                     PhiS_cor*Deg2Rad);
        Sigma_PhiPlusPhiS   = Photon_Factor * Sigma_UU * Asym_PhiPlusPhiS   * sin(1.*Phi_cor*Deg2Rad + PhiS_cor*Deg2Rad);
        Sigma_2PhiPlusPhiS  = Photon_Factor * Sigma_UU * Asym_2PhiPlusPhiS  * sin(2.*Phi_cor*Deg2Rad + PhiS_cor*Deg2Rad);

        Sigma_UT = Sigma_3PhiMinusPhiS 
            + Sigma_2PhiMinusPhiS
            + Sigma_PhiMinusPhiS
            + Sigma_PhiS
            + Sigma_PhiPlusPhiS
            + Sigma_2PhiPlusPhiS;
    
        //correction from the generator: a minus sign here, note that when using He3, we also need to have its 60% polarization, eg. 0.6*0.865*Sigma_UT
        Sigma_Lab = Sigma_UU - 0.6 * 0.865 * Sigma_UT; 
        
        weight     = Lumi_PSF * total_acc_res / fileNO * time * Sigma_Lab;//XS*Lumi*PSF*Acc, give the real rates
        weight_uu  =  Lumi_PSF * total_acc_res / fileNO * time * Sigma_UU;
        weight_ut  =  Lumi_PSF * total_acc_res / fileNO * time * Sigma_UT;
        weight_3m1 =  Lumi_PSF * total_acc_res / fileNO * time * Sigma_3PhiMinusPhiS;
        weight_2m1 =  Lumi_PSF * total_acc_res / fileNO * time * Sigma_2PhiMinusPhiS;
        weight_1m1 =  Lumi_PSF * total_acc_res / fileNO * time * Sigma_PhiMinusPhiS;
        weight_0p1 =  Lumi_PSF * total_acc_res / fileNO * time * Sigma_PhiS;
        weight_1p1 =  Lumi_PSF * total_acc_res / fileNO * time * Sigma_PhiPlusPhiS;
        weight_2p1 =  Lumi_PSF * total_acc_res / fileNO * time * Sigma_2PhiPlusPhiS;
        /*}}}*/

        if(i<10)
            cout<<Form("-- Lumi_PSF =%f, Photon_Factor=%f, Sigma_UU=%f, Sigma_UT=%f, Sigma_Lab=%f", Lumi_PSF, Photon_Factor, Sigma_UU, Sigma_UT, Sigma_Lab)<<endl;

        total_rate_dvmp += weight/time;
        /*Q2>4 Cut: Fill Trees and Missing Mass{{{*/
        if(Qsq>4){
            total_rate_dvmp_4GeV += weight/time;
        
            /*Missing w/o resolutions{{{*/
            P_E0->SetPxPyPzE(beam_px,beam_py,beam_pz, beam_ene);
            P_t->SetPxPyPzE(tgt_px,tgt_py,tgt_pz,tgt_ene);
            
            P_e->SetPxPyPzE(ele_px,ele_py,ele_pz,ele_ene);	
            P_pim->SetPxPyPzE(pim_px,pim_py,pim_pz,pim_ene);	
            P_pro->SetPxPyPzE(pro_px,pro_py,pro_pz,pro_ene);	
            int err = CheckLaws(P_t, P_E0, P_e, P_pim, P_pro);//Check whether momentum and energy conserve first
            if (err < 1e-33){
                cerr<<"---- Momentum and Energy Conservation Laws are broken!! Something is wrong!!!"<<endl;
            }
            //After checking everying is fine, now put back to the realistic case where we don't know the individual target motion
            P_t->SetPxPyPzE(0., 0., 0., nMass);
            MM = GetMM(P_t, P_E0, P_e, P_pim);	
            MP = GetMP(P_t, P_E0, P_e, P_pim);	
            MV = GetMV(P_t, P_E0, P_e, P_pim);	

            hMP_dvmp->Fill(MP, weight/time*total_acc);
            hMT_dvmp->Fill(MV->Theta()*Rad2Deg, weight/time*total_acc);
            hMM_dvmp->Fill(MM, weight/time*total_acc);
            if(MP>0.&&MP<1.2&&MV->Theta()*Rad2Deg>8 && MV->Theta()*Rad2Deg<24.0)
                hMM_dvmp_cut->Fill(MM, weight/time*total_acc);
            /*}}}*/
 
            /*Missing w/o resolutions+nuclear effect{{{*/
            P_E0->SetPxPyPzE(beam_px_cor,beam_py_cor,beam_pz_cor, beam_ene_cor);
            //P_t->SetPxPyPzE(tgt_px,tgt_py,tgt_pz,tgt_ene);
            P_t->SetPxPyPzE(0., 0., 0., nMass);
            
            P_e->SetPxPyPzE(ele_px_cor,ele_py_cor,ele_pz_cor,ele_ene_cor);	
            P_pim->SetPxPyPzE(pim_px_cor,pim_py_cor,pim_pz_cor,pim_ene_cor);	
            P_pro->SetPxPyPzE(pro_px_cor,pro_py_cor,pro_pz_cor,pro_ene_cor);	

            MM_cor = GetMM(P_t, P_E0, P_e, P_pim);	
            MP_cor = GetMP(P_t, P_E0, P_e, P_pim);	
            MV = GetMV(P_t, P_E0, P_e, P_pim);	

            hMP_dvmp_cor->Fill(MP_cor, weight/time*total_acc_cor);
            hMT_dvmp_cor->Fill(MV->Theta()*Rad2Deg, weight/time*total_acc_cor);
            hMM_dvmp_cor->Fill(MM_cor, weight/time*total_acc_cor);
            if(MP>0.&&MP<1.2&&MV->Theta()*Rad2Deg>8 && MV->Theta()*Rad2Deg<24.0)
                hMM_dvmp_cut_cor->Fill(MM_cor, weight/time*total_acc_cor);
            /*}}}*/
            
            /*Missing w/ resolutions{{{*/
            P_E0->SetPxPyPzE(0.,0.,beam_ene_cor, beam_ene_cor);
            //P_t->SetPxPyPzE(tgt_px,tgt_py,tgt_pz,tgt_ene);
            P_t->SetPxPyPzE(0., 0., 0., nMass);
            
            P_e->SetPxPyPzE(ele_px_res,ele_py_res,ele_pz_res,ele_ene_res);	
            P_pim->SetPxPyPzE(pim_px_res,pim_py_res,pim_pz_res,pim_ene_res);	
            P_pro->SetPxPyPzE(pro_px_res,pro_py_res,pro_pz_res,pro_ene_res);	

            MM_res = GetMM(P_t, P_E0, P_e, P_pim);	
            MP_res = GetMP(P_t, P_E0, P_e, P_pim);	
            MV = GetMV(P_t, P_E0, P_e, P_pim);	

            hMP_dvmp_res->Fill(MP_res, weight/time*total_acc_res);
            hMT_dvmp_res->Fill(MV->Theta()*Rad2Deg, weight/time*total_acc_res);
            hMM_dvmp_res->Fill(MM_res, weight/time*total_acc_res);
            if(MP_res>0.&&MP_res<1.2&&MV->Theta()*Rad2Deg>8 && MV->Theta()*Rad2Deg<24.0)
                hMM_dvmp_cut_res->Fill(MM_res, weight/time*total_acc_res);
            /*}}}*/
 
            t2->Fill(); 
            if(!(i%1000))
                cout<<"--- Processed events = "<<i<<"\r";
        }/*Q2>4 cut}}}*/

    }// events loop ends here
    /*}}}*/

    cout<<Form("--- Total DVMP Rate= %f Hz (%f Hz at Q2>4)",total_rate_dvmp,total_rate_dvmp_4GeV)<<endl;
    TString outf_name = new_filename.ReplaceAll(".root",".txt");
    ofstream outf(outf_name.Data());
    outf<<Form("--- Total DVMP Rate= %f Hz (%f Hz at Q2>4)",total_rate_dvmp,total_rate_dvmp_4GeV)<<endl;

    histo->cd();
    hMM_dvmp->Write(); hMM_dvmp_cor->Write(); hMM_dvmp_res->Write(); 
    hMM_dvmp_cut->Write(); hMM_dvmp_cut_cor->Write(); hMM_dvmp_cut_res->Write(); 
    hMP_dvmp->Write(); hMP_dvmp_cor->Write(); hMP_dvmp_res->Write(); 
    hMT_dvmp->Write(); hMT_dvmp_cor->Write(); hMT_dvmp_res->Write(); 

    f2->cd();
    hMM_dvmp->Write(); hMM_dvmp_cor->Write(); hMM_dvmp_res->Write(); 
    hMM_dvmp_cut->Write(); hMM_dvmp_cut_cor->Write(); hMM_dvmp_cut_res->Write(); 
    hMP_dvmp->Write(); hMP_dvmp_cor->Write(); hMP_dvmp_res->Write(); 
    hMT_dvmp->Write(); hMT_dvmp_cor->Write(); hMT_dvmp_res->Write(); 
    t2->Write();
    
    f2->Close();
    histo->Close();

    return 0;
}

/*Double_t GetMM(TLorentzVector* P_E0, TLorentzVector* P_e, TLorentzVector* P_pim){{{*/
Double_t GetMM(TLorentzVector* P_t,TLorentzVector* P_E0, TLorentzVector* P_e, TLorentzVector* P_pim){
    //Calculate the Missing Mass
    double E_miss = (P_t->E() + P_E0->E()) - (P_e->E()+P_pim->E());
    double px_miss = (P_t->Px() + P_E0->Px()) - (P_e->Px()+P_pim->Px()); 
    double py_miss = (P_t->Py() + P_E0->Py()) - (P_e->Py()+P_pim->Py()); 
    double pz_miss = (P_t->Pz() + P_E0->Pz()) - (P_e->Pz()+P_pim->Pz()); 
    double kMM = sqrt( E_miss*E_miss - (pow(px_miss,2)+pow(py_miss,2)+pow(pz_miss,2)) );
    return kMM;
}
/*}}}*/

/*Double_t GetMP(TLorentzVector* P_E0, TLorentzVector* P_e, TLorentzVector* P_pim){{{*/
Double_t GetMP(TLorentzVector* P_t,TLorentzVector* P_E0, TLorentzVector* P_e, TLorentzVector* P_pim){
    //Calculate the missing momentum
    double px_miss = (P_t->Px() + P_E0->Px()) - (P_e->Px()+P_pim->Px()); 
    double py_miss = (P_t->Py() + P_E0->Py()) - (P_e->Py()+P_pim->Py()); 
    double pz_miss = (P_t->Pz() + P_E0->Pz()) - (P_e->Pz()+P_pim->Pz()); 
    double kMP = sqrt( pow(px_miss,2)+pow(py_miss,2)+pow(pz_miss,2) );
    return kMP;
}
/*}}}*/

/*TLorentzVector* GetMV(TLorentzVector* P_E0, TLorentzVector* P_e, TLorentzVector* P_pim){{{*/
TLorentzVector* GetMV(TLorentzVector* P_t,TLorentzVector* P_E0, TLorentzVector* P_e, TLorentzVector* P_pim){
    //Calculate the missing momentum
    double E_miss = (P_t->E() + P_E0->E()) - (P_e->E()+P_pim->E());
    double px_miss = (P_t->Px() + P_E0->Px()) - (P_e->Px()+P_pim->Px()); 
    double py_miss = (P_t->Py() + P_E0->Py()) - (P_e->Py()+P_pim->Py()); 
    double pz_miss = (P_t->Pz() + P_E0->Pz()) - (P_e->Pz()+P_pim->Pz()); 

    TLorentzVector *vec= new TLorentzVector;
    vec->SetPxPyPzE(px_miss, py_miss,pz_miss,E_miss);

    return vec;
}
/*}}}*/

/*Double_t CheckLaws(TLorentzVector* P_t,TLorentzVector* P_E0, TLorentzVector* P_e, TLorentzVector* P_pim, TLorentzVector* P_pro){{{*/
Int_t CheckLaws(TLorentzVector* P_t,TLorentzVector* P_E0, TLorentzVector* P_e, TLorentzVector* P_pim, TLorentzVector* P_pro){
    double energy_check = (P_t->E() + P_E0->E()) - (P_e->E()+P_pim->E()+P_pro->E());
    double px_check =(P_t->Px() + P_E0->Px()) - (P_e->Px()+P_pim->Px()+P_pro->Px()); 
    double py_check =(P_t->Py() + P_E0->Py()) - (P_e->Py()+P_pim->Py()+P_pro->Py()); 
    double pz_check =(P_t->Pz() + P_E0->Pz()) - (P_e->Pz()+P_pim->Pz()+P_pro->Pz()); 

    Int_t err = -1;
    if(fabs(energy_check)<0.01 && fabs(px_check)<0.01 && fabs(py_check)<0.01 && fabs(pz_check)<0.01)
        err = 1;
    else{
        cerr<<Form("*** dE = %f,  dPx = %f, dPy = %f, dPz = %f", energy_check, px_check, py_check, pz_check)<<endl;
        err = -1;
    }
    return err;

}
/*}}}*/
