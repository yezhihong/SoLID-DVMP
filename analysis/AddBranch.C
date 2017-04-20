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
    const Int_t time = 48 * 24 * 3600;

    //CLEO Acceptance
    SIDIS_Acceptance *accpt = new SIDIS_Acceptance();

    TString filename, new_filename;
    TChain *t1 = new TChain("t1");
    
    //const int fileNO = 1000;
    const int fileNO = 999;
    //for(int i=0;i<fileNO;i++){
    for(int i=1;i<fileNO;i++){
        filename = Form("./RootFiles_Up_Simple/DEMP_Ee_11_Events_100000_File_%d.root", i);
        //filename = Form("./RootFiles_Down_Simple/DEMP_Ee_11_Events_100000_File_%d.root", i);
        cout<<"--- Reading in file = "<<filename.Data()<<"\r";
        t1->Add(filename.Data());
    }
    int N_entries = t1->GetEntries();
    cout<<"--- Total Events = "<< N_entries<<endl;
   
    /*Define Variables, Branches{{{*/
    Int_t NRecorded, NGenerated; 
    Double_t W,Qsq, t, t_Para, Epsilon, x, y, z;/*{{{*/
    Double_t t_corr, Qsq_corr, W_corr, x_corr, y_corr, z_corr;
    Double_t tgt_theta, tgt_phi, tgt_ene, tgt_mom, tgt_px, tgt_py, tgt_pz;
    Double_t beam_theta, beam_phi, beam_ene, beam_mom, beam_px, beam_py, beam_pz;
    Double_t beam_corr_theta, beam_corr_phi, beam_corr_ene, beam_corr_mom, beam_corr_px, beam_corr_py, beam_corr_pz;
    Double_t ele_theta, ele_phi, ele_ene, ele_mom, ele_px, ele_py, ele_pz;
    Double_t ele_corr_theta, ele_corr_phi, ele_corr_ene, ele_corr_mom, ele_corr_px, ele_corr_py, ele_corr_pz;
    Double_t pim_theta, pim_phi, pim_ene, pim_mom, pim_px, pim_py, pim_pz;
    Double_t pim_corr_theta, pim_corr_phi, pim_corr_ene, pim_corr_mom, pim_corr_px, pim_corr_py, pim_corr_pz;
    Double_t pro_theta, pro_phi, pro_ene, pro_mom, pro_px, pro_py, pro_pz;
    Double_t pro_corr_theta, pro_corr_phi, pro_corr_ene, pro_corr_mom, pro_corr_px, pro_corr_py, pro_corr_pz;
    Double_t EventWeight, WilliamsWeight, DedrickWeight, CatchenWeight;
    Double_t Theta_Pion_Photon, Phi, PhiS, Phi_corr, PhiS_corr, Vertex_X, Vertex_Y, Vertex_Z;
    Double_t Asym_PhiS, Asym_PhiMinusPhiS, Asym_2PhiMinusPhiS, Asym_3PhiMinusPhiS,Asym_PhiPlusPhiS,Asym_2PhiPlusPhiS; 
    Double_t Sigma_PhiS, Sigma_PhiMinusPhiS, Sigma_2PhiMinusPhiS, Sigma_3PhiMinusPhiS, Sigma_PhiPlusPhiS, Sigma_2PhiPlusPhiS;
    Double_t Sigma_Lab, Sigma_UU, Sigma_UT, SSAsym, SineAsym;
    Double_t Sig_L, Sig_T, Sig_LT, Sig_TT;
    Double_t Flux_Factor_RF, Flux_Factor_Col, Jacobian_CM, Jacobian_CM_RF, Jacobian_CM_Col, A_Factor;
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
    t1->SetBranchAddress("Electron_Corrected_Theta_Col",              &beam_corr_theta );
    t1->SetBranchAddress("Electron_Corrected_Phi_Col",                &beam_corr_phi   );
    t1->SetBranchAddress("Electron_Corrected_Energy_Col_GeV",         &beam_corr_ene   );
    t1->SetBranchAddress("Electron_Corrected_Mom_Col_GeV",            &beam_corr_mom   );
    t1->SetBranchAddress("Electron_Corrected_MomX_Col_GeV",           &beam_corr_px    );
    t1->SetBranchAddress("Electron_Corrected_MomY_Col_GeV",           &beam_corr_py    );
    t1->SetBranchAddress("Electron_Corrected_MomZ_Col_GeV",           &beam_corr_pz    );/*}}}*/

    //Scattered Electron Quantities/*{{{*/
    t1->SetBranchAddress("ScatElec_Theta_Col",                        &ele_theta     );
    t1->SetBranchAddress("ScatElec_Phi_Col",                          &ele_phi       );
    t1->SetBranchAddress("ScatElec_Energy_Col_GeV",                   &ele_ene       );
    t1->SetBranchAddress("ScatElec_Mom_Col_GeV",                      &ele_mom       );
    t1->SetBranchAddress("ScatElec_MomX_Col_GeV",                     &ele_px        );
    t1->SetBranchAddress("ScatElec_MomY_Col_GeV",                     &ele_py        );
    t1->SetBranchAddress("ScatElec_MomZ_Col_GeV",                     &ele_pz        );
    t1->SetBranchAddress("ScatElec_Corrected_Theta_Col",              &ele_corr_theta);
    t1->SetBranchAddress("ScatElec_Corrected_Phi_Col",                &ele_corr_phi  );
    t1->SetBranchAddress("ScatElec_Corrected_Energy_Col_GeV",         &ele_corr_ene  );
    t1->SetBranchAddress("ScatElec_Corrected_Mom_Col_GeV",            &ele_corr_mom  );
    t1->SetBranchAddress("ScatElec_Corrected_MomX_Col_GeV",           &ele_corr_px   );
    t1->SetBranchAddress("ScatElec_Corrected_MomY_Col_GeV",           &ele_corr_py   );
    t1->SetBranchAddress("ScatElec_Corrected_MomZ_Col_GeV",           &ele_corr_pz   );/*}}}*/

    //Pion- Quantities/*{{{*/
    t1->SetBranchAddress("Pion_Theta_Col",                            &pim_theta     );
    t1->SetBranchAddress("Pion_Phi_Col",                              &pim_phi       );
    t1->SetBranchAddress("Pion_Energy_Col_GeV",                       &pim_ene       );
    t1->SetBranchAddress("Pion_Mom_Col_GeV",                          &pim_mom       );
    t1->SetBranchAddress("Pion_MomX_Col_GeV",                         &pim_px        );
    t1->SetBranchAddress("Pion_MomY_Col_GeV",                         &pim_py        );
    t1->SetBranchAddress("Pion_MomZ_Col_GeV",                         &pim_pz        );
    t1->SetBranchAddress("Pion_Corrected_Theta_Col",                  &pim_corr_theta);
    t1->SetBranchAddress("Pion_Corrected_Phi_Col",                    &pim_corr_phi  );
    t1->SetBranchAddress("Pion_Corrected_Energy_Col_GeV",             &pim_corr_ene  );
    t1->SetBranchAddress("Pion_Corrected_Mom_Col_GeV",                &pim_corr_mom  );
    t1->SetBranchAddress("Pion_Corrected_MomX_Col_GeV",               &pim_corr_px   );
    t1->SetBranchAddress("Pion_Corrected_MomY_Col_GeV",               &pim_corr_py   );
    t1->SetBranchAddress("Pion_Corrected_MomZ_Col_GeV",               &pim_corr_pz   );/*}}}*/

    //Recoil Proton Quantities/*{{{*/
    t1->SetBranchAddress("RecoilProton_Theta_Col",                    &pro_theta          );
    t1->SetBranchAddress("RecoilProton_Phi_Col",                      &pro_phi            );
    t1->SetBranchAddress("RecoilProton_Energy_Col_GeV",               &pro_ene            );
    t1->SetBranchAddress("RecoilProton_Mom_Col_GeV",                  &pro_mom            );
    t1->SetBranchAddress("RecoilProton_MomX_Col_GeV",                 &pro_px             );
    t1->SetBranchAddress("RecoilProton_MomY_Col_GeV",                 &pro_py             );
    t1->SetBranchAddress("RecoilProton_MomZ_Col_GeV",                 &pro_pz             );
    t1->SetBranchAddress("RecoilProton_Corrected_Theta_Col",          &pro_corr_theta     );
    t1->SetBranchAddress("RecoilProton_Corrected_Phi_Col",            &pro_corr_phi       );
    t1->SetBranchAddress("RecoilProton_Corrected_Energy_Col_GeV",     &pro_corr_ene       );
    t1->SetBranchAddress("RecoilProton_Corrected_Mom_Col_GeV",        &pro_corr_mom       );
    t1->SetBranchAddress("RecoilProton_Corrected_MomX_Col_GeV",       &pro_corr_px        );
    t1->SetBranchAddress("RecoilProton_Corrected_MomY_Col_GeV",       &pro_corr_py        );
    t1->SetBranchAddress("RecoilProton_Corrected_MomZ_Col_GeV",       &pro_corr_pz        );/*}}}*/

    t1->SetBranchAddress("Phi",                                       &Phi                 );/*{{{*/
    t1->SetBranchAddress("PhiS",                                      &PhiS                );
    t1->SetBranchAddress("Phi_Corrected",                             &Phi_corr            );
    t1->SetBranchAddress("PhiS_Corrected",                            &PhiS_corr           );/*}}}*/

    //Different XS Quantities/
    //In the lab frame/*{{{*/
    t1->SetBranchAddress("ZASigma_Lab",                               &Sigma_Lab         );//Pleaes explain here
    t1->SetBranchAddress("ZASigma_UU_Col",                            &Sigma_UU          );//Pleaes explain here
    t1->SetBranchAddress("RorySigma_UT_Col",                          &Sigma_UT        );//Pleaes explain here

    //XSs corresponding to different Asymmetries
    //dSigma_k_Sin(...) = kFactor * Sin(...) * fAsym_k * dSigma_UU, in lab frame
    t1->SetBranchAddress("Sig_PhiS_Col",                                 &Sigma_PhiS           );
    t1->SetBranchAddress("Sig_Phi_Minus_PhiS_Col",                       &Sigma_PhiMinusPhiS );
    t1->SetBranchAddress("Sig_2Phi_Minus_PhiS_Col",                      &Sigma_2PhiMinusPhiS);
    t1->SetBranchAddress("Sig_3Phi_Minus_PhiS_Col",                      &Sigma_3PhiMinusPhiS);
    t1->SetBranchAddress("Sig_Phi_Plus_PhiS_Col",                        &Sigma_PhiPlusPhiS  );
    t1->SetBranchAddress("Sig_2Phi_Plus_PhiS_Col",                       &Sigma_2PhiPlusPhiS );

    //6 Asymmetries, "_Col" means in the lab frame, w/o that means the rest frame
    t1->SetBranchAddress("SSAsym",                                       &SSAsym           );
    t1->SetBranchAddress("SineAsym",                                     &SineAsym           );
    t1->SetBranchAddress("AsymPhi_S_Col",                                &Asym_PhiS           );
    t1->SetBranchAddress("AsymPhiMinusPhi_S_Col",                        &Asym_PhiMinusPhiS   );
    t1->SetBranchAddress("AsymPhiPlusPhi_S_Col",                         &Asym_PhiPlusPhiS    );
    t1->SetBranchAddress("Asym3PhiMinusPhi_S_Col",                       &Asym_3PhiMinusPhiS  );
    t1->SetBranchAddress("AsymPhiMinusPhi_S_Col",                        &Asym_2PhiMinusPhiS  );
    t1->SetBranchAddress("Asym2PhiMinusPhi_S_Col",                       &Asym_2PhiMinusPhiS  );/*}}}*/

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

    t1->SetBranchAddress("T_Corrected_GeV",                           &t_corr               );
    t1->SetBranchAddress("Qsq_Corrected_GeV",                         &Qsq_corr             );
    t1->SetBranchAddress("W_Corrected_GeV",                           &W_corr               );
    t1->SetBranchAddress("x_Corrected",                               &x_corr               );
    t1->SetBranchAddress("y_Corrected",                               &y_corr               );
    t1->SetBranchAddress("z_Corrected",                               &z_corr               );

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

    new_filename = "./rootfiles/DEMP_Ee_11_4_up_simple.root";
    //new_filename = "./rootfiles/DEMP_Ee_11_4_down_simple.root";
    cout<<"---  Saving in file = "<<new_filename.Data()<<endl;
    TFile* f2=new TFile(new_filename.Data(),"recreate"); 
   
    /*Define new Tree and new Branch{{{*/
    double ele_acc_f, ele_acc_l, pim_acc_f,pim_acc_l, pro_acc_f,pro_acc_l;
    double MM, MM_res,dilute,weight;
    double ele_mom_res, ele_theta_res, ele_phi_res;
    double pim_mom_res, pim_theta_res, pim_phi_res;
    double pro_mom_res, pro_theta_res, pro_phi_res;

    //TTree *t2 = t1->CloneTree(0);
    TTree *t2 = new TTree("T","a new tree");

    t2->Branch("NRecorded",       &NRecorded,     "NRecorded/I");/*{{{*/
    t2->Branch("NGenerated",       &NGenerated,     "NGenerated/I");

    t2->Branch("Epsilon", &Epsilon, "Epsilon/D" );
    t2->Branch("Qsq", &Qsq ,"Qsq/D");
    t2->Branch("t", &t ,"t/D");
    t2->Branch("W", &W ,"W/D");
    t2->Branch("x", &x ,"x/D");
    t2->Branch("y", &y ,"y/D");
    t2->Branch("z", &z ,"z/D");

    t2->Branch("Qsq_corr", &Qsq_corr ,"Qsq_corr/D");
    t2->Branch("t_corr", &t_corr ,"t_corr/D");
    t2->Branch("W_corr", &W_corr ,"W_corr/D");
    t2->Branch("x_corr", &x_corr ,"x_corr/D");
    t2->Branch("y_corr", &y_corr ,"y_corr/D");
    t2->Branch("z_corr", &z_corr ,"z_corr/D"); 

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
    t2->Branch("Phi_corr",  &Phi_corr,  "Phi_corr/D");
    t2->Branch("PhiS_corr", &PhiS_corr, "PhiS_corr/D");/*}}}*/

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
    t2->Branch("tgt_ene", &tgt_mom, "tgt_mom/D");/*}}}*/

    t2->Branch("beam_px", &beam_px, "beam_px/D");/*{{{*/
    t2->Branch("beam_py", &beam_py, "beam_py/D");
    t2->Branch("beam_pz", &beam_pz, "beam_pz/D");
    t2->Branch("beam_theta", &beam_theta, "beam_theta/D");
    t2->Branch("beam_phi", &beam_phi, "beam_phi/D");
    t2->Branch("beam_ene", &beam_ene, "beam_ene/D");
    t2->Branch("beam_ene", &beam_mom, "beam_mom/D");

    t2->Branch("beam_corr_px", &beam_corr_px, "beam_corr_px/D");
    t2->Branch("beam_corr_py", &beam_corr_py, "beam_corr_py/D");
    t2->Branch("beam_corr_pz", &beam_corr_pz, "beam_corr_pz/D");
    t2->Branch("beam_corr_theta", &beam_corr_theta, "beam_corr_theta/D");
    t2->Branch("beam_corr_phi", &beam_corr_phi, "beam_corr_phi/D");
    t2->Branch("beam_corr_ene", &beam_corr_ene, "beam_corr_ene/D");
    t2->Branch("beam_corr_ene", &beam_corr_mom, "beam_corr_mom/D");/*}}}*/

    t2->Branch("pim_ene",   &pim_ene   ,"pim_ene/D");/*{{{*/
    t2->Branch("pim_px",    &pim_px    ,"pim_px/D");
    t2->Branch("pim_py",    &pim_py    ,"pim_py/D");
    t2->Branch("pim_pz",    &pim_pz    ,"pim_pz/D");
    t2->Branch("pim_mom",   &pim_mom   ,"pim_mom/D");
    t2->Branch("pim_theta", &pim_theta ,"pim_theta/D");
    t2->Branch("pim_phi",   &pim_phi   ,"pim_phi/D");

    t2->Branch("pim_corr_ene",   &pim_corr_ene   ,"pim_corr_ene/D");
    t2->Branch("pim_corr_px",    &pim_corr_px    ,"pim_corr_px/D");
    t2->Branch("pim_corr_py",    &pim_corr_py    ,"pim_corr_py/D");
    t2->Branch("pim_corr_pz",    &pim_corr_pz    ,"pim_corr_pz/D");
    t2->Branch("pim_corr_mom",   &pim_corr_mom   ,"pim_corr_mom/D");
    t2->Branch("pim_corr_theta", &pim_corr_theta ,"pim_corr_theta/D");
    t2->Branch("pim_corr_phi",   &pim_corr_phi   ,"pim_corr_phi/D");/*}}}*/

    t2->Branch("ele_ene",   &ele_ene   ,"ele_ene/D");/*{{{*/
    t2->Branch("ele_px",    &ele_px    ,"ele_px/D");
    t2->Branch("ele_py",    &ele_py    ,"ele_py/D");
    t2->Branch("ele_pz",    &ele_pz    ,"ele_pz/D");
    t2->Branch("ele_mom",   &ele_mom   ,"ele_mom/D");
    t2->Branch("ele_theta", &ele_theta ,"ele_theta/D");
    t2->Branch("ele_phi",   &ele_phi   ,"ele_phi/D");

    t2->Branch("ele_corr_ene",   &ele_corr_ene   ,"ele_corr_ene/D");
    t2->Branch("ele_corr_px",    &ele_corr_px    ,"ele_corr_px/D");
    t2->Branch("ele_corr_py",    &ele_corr_py    ,"ele_corr_py/D");
    t2->Branch("ele_corr_pz",    &ele_corr_pz    ,"ele_corr_pz/D");
    t2->Branch("ele_corr_mom",   &ele_corr_mom   ,"ele_corr_mom/D");
    t2->Branch("ele_corr_theta", &ele_corr_theta ,"ele_corr_theta/D");
    t2->Branch("ele_corr_phi",   &ele_corr_phi   ,"ele_corr_phi/D");/*}}}*/

    t2->Branch("pro_ene",   &pro_ene   ,"pro_ene/D"); /*{{{*/
    t2->Branch("pro_px",    &pro_px    ,"pro_px/D");
    t2->Branch("pro_py",    &pro_py    ,"pro_py/D");
    t2->Branch("pro_pz",    &pro_pz    ,"pro_pz/D");
    t2->Branch("pro_mom",   &pro_mom   ,"pro_mom/D");
    t2->Branch("pro_theta", &pro_theta ,"pro_theta/D");
    t2->Branch("pro_phi",   &pro_phi   ,"pro_phi/D");

    t2->Branch("pro_corr_ene",   &pro_corr_ene   ,"pro_corr_ene/D");
    t2->Branch("pro_corr_px",    &pro_corr_px    ,"pro_corr_px/D");
    t2->Branch("pro_corr_py",    &pro_corr_py    ,"pro_corr_py/D");
    t2->Branch("pro_corr_pz",    &pro_corr_pz    ,"pro_corr_pz/D");
    t2->Branch("pro_corr_mom",   &pro_corr_mom   ,"pro_corr_mom/D");
    t2->Branch("pro_corr_theta", &pro_corr_theta ,"pro_corr_theta/D");
    t2->Branch("pro_corr_phi",   &pro_corr_phi   ,"pro_corr_phi/D");/*}}}*/

    //Add SoLID acceptance/*{{{*/
    t2->Branch("ele_acc_f",     &ele_acc_f,     "ele_acc_f/D");
    t2->Branch("ele_acc_l",     &ele_acc_l,     "ele_acc_l/D");
    t2->Branch("pim_acc_f",     &pim_acc_f,     "pim_acc_f/D");
    t2->Branch("pim_acc_l",     &pim_acc_l,     "pim_acc_l/D");
    t2->Branch("pro_acc_f",     &pro_acc_f,     "pro_acc_f/D");
    t2->Branch("pro_acc_l",     &pro_acc_l,     "pro_acc_l/D");/*}}}*/

    //Add detector resolutions/*{{{*/
    t2->Branch("ele_mom_res",   &ele_mom_res,   "ele_mom_res/D");
    t2->Branch("ele_theta_res", &ele_theta_res, "ele_theta_res/D");
    t2->Branch("ele_phi_res",   &ele_phi_res,   "ele_phi_res/D");
    t2->Branch("pim_mom_res",   &pim_mom_res,   "pim_mom_res/D");
    t2->Branch("pim_theta_res", &pim_theta_res, "pim_theta_res/D");
    t2->Branch("pim_phi_res",   &pim_phi_res,   "pim_phi_res/D");
    t2->Branch("pro_mom_res",   &pro_mom_res,   "pro_mom_res/D");
    t2->Branch("pro_theta_res", &pro_theta_res, "pro_theta_res/D");
    t2->Branch("pro_phi_res",   &pro_phi_res,   "pro_phi_res/D");/*}}}*/

    //Add other quantities/*{{{*/
    t2->Branch("weight",        &weight,        "weight/D");
    t2->Branch("dilute",        &dilute,        "dilute/D");
    //t2->Branch("PSF",           &PSF,           "PSF/D");
    t2->Branch("MM",            &MM,            "MM/D");
    t2->Branch("MM_res", &MM_res, "MM_res/D");/*}}}*/
    /*}}}*/

    /*Vectors & Histograms{{{*/
    TLorentzVector *P_E0 = new TLorentzVector();//incoming electron
    TLorentzVector *P_e = new TLorentzVector();//scattered electron with Eloss
    TLorentzVector *P_pim = new TLorentzVector();//photon

    TLorentzVector *P_pro = new TLorentzVector();//proton or neutron
    TLorentzVector *P_t = new TLorentzVector();//target, either proton or neutron
    TLorentzVector *P_e_res = new TLorentzVector();//scattered electron with resolution
    TLorentzVector *P_pim_res = new TLorentzVector();//photon with resolution
    TLorentzVector *P_pro_res = new TLorentzVector();//photon with resolution

    TH1F *hMM = new TH1F("hMM",Form("DVMP Missing Mass at %3.1f GeV",EBeam), 400,0.0, 6.0);
    hMM->SetXTitle("Hadron Missing Mass (GeV)");
    hMM->GetXaxis()->CenterTitle(1);
    hMM->SetYTitle("Rate (Hz)");
    hMM->GetYaxis()->CenterTitle(1);
    TH1F *hMM_res = new TH1F("hMM_res",Form("DVMP Missing Mass at %3.1f GeV",EBeam), 400,0.0, 6.0);
    hMM_res->SetXTitle("Hadron Missing Mass (GeV)");
    hMM_res->GetXaxis()->CenterTitle(1);
    hMM_res->SetYTitle("Rate (Hz)");
    hMM_res->GetYaxis()->CenterTitle(1);
    /*}}}*/

    /*Fill in DVMP Missing Mass{{{*/
    double total_rate_dvmp=0.0; //Make sure the value is right
    for(Long64_t i=0;i<N_entries;i++){
        t1->GetEntry(i);
        if(Qsq>4){
            //if(ele_theta>7.5&&ele_theta<24.5&&ele_mom>1.&&ele_mom<11
            //&&pim_theta>7.5&&pim_theta<24.5&&pim_mom>1.&&pim_mom<11
            //&& W>=2.&&Qsq>1.0&&SigmaPara>1e-33){//any additional cuts should be added in here

            ///////////////////////////////////////////////////////////////////////////
            //Consider the detector resolution here*{{{*/
            //make sure to use the corrected quantitie for multiple scattering and eloss
            //Electron/*{{{*/
            ele_mom_res = gRandom->Gaus(ele_corr_mom, Sigma_Dp_E*ele_corr_mom);//GeV, for electron, E ~= P
            ele_theta_res = gRandom->Gaus(ele_corr_theta*Deg2Rad, Sigma_Theta_E);//rad
            ele_phi_res = gRandom->Gaus(ele_corr_phi*Deg2Rad, Sigma_Phi_E);//rad

            double ePx_res = ele_mom_res * sin(ele_theta_res)*cos(ele_phi_res); 
            double ePy_res = ele_mom_res * sin(ele_theta_res)*sin(ele_phi_res); 
            double ePz_res = ele_mom_res * cos(ele_theta_res);
            double eE_res = sqrt(ele_mom_res*ele_mom_res + eMass*eMass);
            ele_phi_res *=Rad2Deg; ele_theta_res *= Rad2Deg;
            P_e_res->SetPxPyPzE(ePx_res, ePy_res, ePz_res, eE_res);	/*}}}*/

            //Pion/*{{{*/
            pim_mom_res = gRandom->Gaus(pim_corr_mom, Sigma_Dp_E*pim_corr_mom);//GeV, for electron, E ~= P
            pim_theta_res = gRandom->Gaus(pim_corr_theta*Deg2Rad, Sigma_Theta_E);//rad
            pim_phi_res = gRandom->Gaus(pim_corr_phi*Deg2Rad, Sigma_Phi_E);//rad

            double pim_Px_res = pim_mom_res * sin(pim_theta_res)*cos(pim_phi_res); 
            double pim_Py_res = pim_mom_res * sin(pim_theta_res)*sin(pim_phi_res); 
            double pim_Pz_res = pim_mom_res * cos(pim_theta_res);
            double pimE_res = sqrt(pim_mom_res*pim_mom_res + piMass*piMass);
            pim_phi_res *=Rad2Deg; pim_theta_res *= Rad2Deg;
            P_pim_res->SetPxPyPzE(pim_Px_res, pim_Py_res, pim_Pz_res, pimE_res);	/*}}}*/

            //Proton/*{{{*/
            pro_mom_res = gRandom->Gaus(pro_corr_mom, Sigma_Dp_E*pro_corr_mom);//GeV, for electron, E ~= P
            pro_theta_res = gRandom->Gaus(pro_corr_theta*Deg2Rad, Sigma_Theta_E);//rad
            pro_phi_res = gRandom->Gaus(pro_corr_phi*Deg2Rad, Sigma_Phi_E);//rad

            double pro_Px_res = pro_mom_res * sin(pro_theta_res)*cos(pro_phi_res); 
            double pro_Py_res = pro_mom_res * sin(pro_theta_res)*sin(pro_phi_res); 
            double pro_Pz_res = pro_mom_res * cos(pro_theta_res);
            double proE_res = sqrt(pro_mom_res*pro_mom_res + piMass*piMass);
            pro_phi_res *=Rad2Deg; pro_theta_res *= Rad2Deg;
            P_pro_res->SetPxPyPzE(pro_Px_res, pro_Py_res, pro_Pz_res, proE_res);	/*}}}*/
            /*}}}*/

            /*Get acceptance of e and pi-{{{*/
            //Make sure to use the corrected quantities for multipile scattering and eloss effects
            //Do not use the smeared quantities since we are about whether particles are in the accepntace or not, but not how good we measure
            /*Elec Acc {{{*/
            ele_acc_f = accpt->GetAcc("e-","forward", ele_corr_mom, ele_corr_theta);
            ele_acc_l = accpt->GetAcc("e-","large", ele_corr_mom, ele_corr_theta);
            if(ele_corr_mom<1.0||ele_corr_theta>14.8||ele_corr_theta<8.0)//GeV, CLEO
                ele_acc_f=0.0;//Farward-Angle EC Cut at 1 GeV
            if(ele_corr_mom<3.5||ele_corr_theta<16.0||ele_corr_theta>24)//GeV,CLEO
                ele_acc_l=0.0; //Larger-Angle EC Cut at 3 GeV
            if(ele_acc_f>1.) 
                ele_acc_f=1.0; 
            if(ele_acc_l>1.) 
                ele_acc_l=1.0; 

            //ele_acc_f = accpt->GetAcc("e-","forward", ele_mom, ele_theta);
            //ele_acc_l = accpt->GetAcc("e-","large", ele_mom, ele_theta);
            //if(ele_mom<1.0||ele_theta>14.8||ele_theta<8.0)//GeV, CLEO
                //ele_acc_f=0.0;//Farward-Angle EC Cut at 1 GeV
            //if(ele_mom<3.5||ele_theta<16.0||ele_theta>24)//GeV,CLEO
                //ele_acc_l=0.0; //Larger-Angle EC Cut at 3 GeV
            //if(ele_acc_f>1.) 
                //ele_acc_f=1.0; 
            //if(ele_acc_l>1.) 
                //ele_acc_l=1.0; /*}}}*/
            
            /*Pion Acc{{{*/
            pim_acc_f = accpt->GetAcc("pi-","forward", pim_corr_mom, pim_corr_theta);
            pim_acc_l = accpt->GetAcc("pi-","large", pim_corr_mom, pim_corr_theta);
            if(pim_corr_theta>14.8||pim_corr_theta<8.0||pim_corr_mom<0.||pim_corr_mom>11.)//GeV, CLEO
            pim_acc_f=0.0;
            if(pim_corr_theta<16.0||pim_corr_theta>24.0||pim_corr_mom<0.||pim_corr_mom>11.)//GeV, CLEO
            pim_acc_l=0.0; 
            if(pim_acc_f>1.) 
            pim_acc_f=1.0; 
            if(pim_acc_l>1.) 
            pim_acc_l=1.0; 

            //pim_acc_f = accpt->GetAcc("pi-","forward", pim_mom, pim_theta);
            //pim_acc_l = accpt->GetAcc("pi-","large", pim_mom, pim_theta);
            //if(pim_theta>14.8||pim_theta<8.0||pim_mom<0.||pim_mom>11.)//GeV, CLEO
                //pim_acc_f=0.0;
            //if(pim_theta<16.0||pim_theta>24.0||pim_mom<0.||pim_mom>11.)//GeV, CLEO
                //pim_acc_l=0.0; 
            //if(pim_acc_f>1.) 
                //pim_acc_f=1.0; 
            //if(pim_acc_l>1.) 
                //pim_acc_l=1.0; 
            /*}}}*/

            /*Proton Acc {{{*/
            //The momentum cut is applied by EC while for proton, we don't reply on EC to tell the acceptance.
            //What I assume here is that we can detecto all energy range of protons,
            //unlike electrons which need to be separated from pions
            //     pro_acc_f = accpt->GetAcc("p","forward", pro_mom, pro_theta);
            //     pro_acc_l = accpt->GetAcc("p","large", pro_mom, pro_theta);
            //   pro_acc_f = accpt->GetThetaAcc("p","forward", pro_theta);
            //   pro_acc_l = accpt->GetThetaAcc("p","large", pro_theta);

            pro_acc_f = 1.0;
            pro_acc_l = 1.0;
            if(pro_corr_theta>14.8||pro_corr_theta<8.0||pro_corr_mom<0.||pro_corr_mom>11.)//GeV, CLEO
                pro_acc_f=0.0;
            if(pro_corr_theta<16.0||pro_corr_theta>24.0||pro_corr_mom<0.||pro_corr_mom>11.)//GeV, CLEO
                pro_acc_l=0.0; 
            if(pro_acc_f>1.) 
                pro_acc_f=1.0; 
            if(pro_acc_l>1.) 
                pro_acc_l=1.0; 

            //pro_acc_f = 1.0;
            //pro_acc_l = 1.0;
            //if(pro_theta>14.8||pro_theta<8.0||pro_mom<0.||pro_mom>11.)//GeV, CLEO
                //pro_acc_f=0.0;
            //if(pro_theta<16.0||pro_theta>24.0||pro_mom<0.||pro_mom>11.)//GeV, CLEO
                //pro_acc_l=0.0; 
            //if(pro_acc_f>1.) 
                //pro_acc_f=1.0; 
            //if(pro_acc_l>1.) 
                //pro_acc_l=1.0; 
            /*}}}*/

            double ele_acceptance=(ele_acc_f+ele_acc_l);
            //double pim_acceptance=(pim_acc_l+pim_acc_f);
            double pim_acceptance=pim_acc_f;
            double pro_acceptance=pro_acc_f+pro_acc_l;

            double total_acceptance=ele_acceptance*pim_acceptance*pro_acceptance;
            /*}}}*/

            dilute = Epsilon*Sig_L/(Epsilon*Sig_L + Sig_T);
            weight = EventWeight * total_acceptance/fileNO * time;//XS*Lumi*PSF*Acc, give the real rates
            total_rate_dvmp += weight/time;

            /*Missing Mass w/o resolutions{{{*/
            P_E0->SetPxPyPzE(0.,0.,EBeam, EBeam);
            P_t->SetPxPyPzE(0.,0.,0., nMass);
            P_e->SetPxPyPzE(ele_px,ele_py,ele_pz,ele_ene);	
            P_pim->SetPxPyPzE(pim_px,pim_py,pim_pz,pim_ene);	
            P_pro->SetPxPyPzE(pro_px,pro_py,pro_pz,pro_ene);	

            int err = CheckLaws(P_t, P_E0, P_e, P_pim, P_pro);//Check whether momentum and energy conserve first
            if (err < 1e-33){
            cerr<<"---- Momentum and Energy Conservation Laws are broken!! Something is wrong!!!"<<endl;
            }

            MM = GetMM(P_t, P_E0, P_e, P_pim);	
            if(ele_theta>7.5&&ele_theta<24.5&&ele_mom>1.&&ele_mom<11
                    &&pim_theta>7.5&&pim_theta<24.5&&pim_mom>1.&&pim_mom<11
                    && W>=2.&&Qsq>1.0&&Sigma_Lab>1e-33)//any additional cuts should be added in here
                hMM->Fill(MM, weight*total_acceptance);
            /*}}}*/

            /*Missing Mass w/ detector resolutions{{{*/
            /////////////////////////////////////////
            MM_res = GetMM(P_t, P_E0, P_e_res, P_pim_res);	
            if(ele_theta>7.5&&ele_theta<24.5&&ele_mom>1.&&ele_mom<11
                    &&pim_theta>7.5&&pim_theta<24.5&&pim_mom>1.&&pim_mom<11
                    && W>=2.&&Qsq>1.0&&Sigma_Lab>1e-33)//any additional cuts should be added in here
                hMM_res->Fill(MM_res, weight*total_acceptance);	
            //////////////////////////////////////////*}}}*/
            t2->Fill(); 
            if(!(i%1000))
                cout<<"--- Processed events = "<<i<<"\r";
        }

    }// events loop ends here
    cout<<"--- Total DVMP Rate="<<total_rate_dvmp<<endl;
    /*}}}*/

    f2->cd();
    hMM->Write(); hMM_res->Write(); t2->Write();
    f2->Close();

    return 0;
}

/*Double_t GetMM(TLorentzVector* P_E0, TLorentzVector* P_e, TLorentzVector* P_pim){{{*/
Double_t GetMM(TLorentzVector* P_t,TLorentzVector* P_E0, TLorentzVector* P_e, TLorentzVector* P_pim){
    Double_t kMM = 0.0;

    double E_miss = (P_t->E() + P_E0->E()) - (P_e->E()+P_pim->E());
    double px_miss = (P_t->Px() + P_E0->Px()) - (P_e->Px()+P_pim->Px()); 
    double py_miss = (P_t->Py() + P_E0->Py()) - (P_e->Py()+P_pim->Py()); 
    double pz_miss = (P_t->Pz() + P_E0->Pz()) - (P_e->Pz()+P_pim->Pz()); 

    kMM = sqrt( E_miss*E_miss - (pow(px_miss,2)+pow(py_miss,2)+pow(pz_miss,2)) );

    return kMM;
}
/*}}}*/

/*Double_t CheckLaws(TLorentzVector* P_t,TLorentzVector* P_E0, TLorentzVector* P_e, TLorentzVector* P_pim, TLorentzVector* P_h){{{*/
Int_t CheckLaws(TLorentzVector* P_t,TLorentzVector* P_E0, TLorentzVector* P_e, TLorentzVector* P_pim, TLorentzVector* P_h){
    Int_t err = -1;

    double energy_check = (P_t->E() + P_E0->E()) - (P_e->E()+P_pim->E()+P_h->E());
    double px_check =(P_t->Px() + P_E0->Px()) - (P_e->Px()+P_pim->Px()+P_h->Px()); 
    double py_check =(P_t->Py() + P_E0->Py()) - (P_e->Py()+P_pim->Py()+P_h->Py()); 
    double pz_check =(P_t->Pz() + P_E0->Pz()) - (P_e->Pz()+P_pim->Pz()+P_h->Pz()); 

    if(fabs(energy_check)<0.01 && fabs(px_check)<0.01 && fabs(py_check)<0.01 && fabs(pz_check)<0.01)
        err = 1;
    else{

        //	cerr<<Form("*** dE = %f,  dPx = %f, dPy = %f, dPz = %f", energy_check, px_check, py_check, pz_check)<<endl;

        err = -1;
    }
    return err;

}
/*}}}*/
