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

using namespace std;
const double Deg2Rad = TMath::DegToRad();
const double Rad2Deg = TMath::RadToDeg();
const double nMass = 0.939565;
const double EBeam = 11.0;
Int_t CheckLaws(TLorentzVector* P_t,TLorentzVector* P_E0, TLorentzVector* P_e, TLorentzVector* P_g, TLorentzVector* P_pro);
Double_t GetMM( TLorentzVector* P_t,TLorentzVector* P_E0, TLorentzVector* P_e, TLorentzVector* P_g);
Double_t GetMP( TLorentzVector* P_t,TLorentzVector* P_E0, TLorentzVector* P_e, TLorentzVector* P_g);
TLorentzVector* GetMV(TLorentzVector* P_t,TLorentzVector* P_E0, TLorentzVector* P_e, TLorentzVector* P_pim);

const double MP_CUT = 1.0;
int main(){
    //int main(){
    gStyle->SetOptStat(0);

    TString type_name = "";/*{{{*/
    Int_t type=0;
    cout<<"--- type = (1->simple, 2->mult, 3->mult_fsi, 4->fermi, 5->mult-fermi)"; cin>>type;
    if(type==1) type_name ="simple";
    if(type==2) type_name ="mult";
    if(type==3) type_name ="mult_fsi";
    if(type==4) type_name ="fermi";
    if(type==5) type_name ="mult_nofermi";

    TString pol_name="";
    int pol=0;
    cout<<"--- Pol = (1->up, 2->down) "; cin >> pol;
    if(pol==1) pol_name="up";
    if(pol==2) pol_name="down";/*}}}*/

    TString prefix = "";
    TString finalfile = Form("../rootfiles/DEMP_%s_%s.root", pol_name.Data(), type_name.Data());
    cout<<"-- Reading in file  "<<finalfile.Data()<<endl;
    TFile *file = new TFile(finalfile.Data(),"r");
    TTree *t0 = (TTree*) file->Get("T");
    int N_entries = t0->GetEntries();
    
    TString histoname = Form("histo_%s_%s_cut1.root", pol_name.Data(), type_name.Data());
    
    Double_t W,Qsq, t,tp, t_Para, Epsilon, x, y, z;/*{{{*/
    Double_t t_cor, tp_cor, Qsq_cor, W_cor, x_cor, y_cor, z_cor;
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
    Double_t Flux_Factor_RF, Flux_Factor_Col, Jacobian_CM, Jacobian_CM_RF, Jacobian_CM_Col, A_Factor, time, Photon_Factor, Photon_Theta;
    Int_t NRecorded, NGenerated, fileNO; 

    double ele_acc_f, ele_acc_l, pim_acc_f,pim_acc_l, pro_acc_f,pro_acc_l,total_acc, total_acc_cor, total_acc_res;
    double MM, MM_res, MM_cor, MP, MP_res, MP_cor,Lumi_PSF, dilute,weight, weight_uu, weight_ut, weight_3m1, weight_2m1, weight_1m1, weight_0p1, weight_1p1, weight_2p1;
    /*}}}*/
    
    /*Define old Tree and old Branch{{{*/
    t0->SetBranchAddress("NRecorded",       &NRecorded);/*{{{*/
    t0->SetBranchAddress("NGenerated",       &NGenerated);

    t0->SetBranchAddress("Epsilon", &Epsilon);
    t0->SetBranchAddress("Qsq", &Qsq );
    t0->SetBranchAddress("t", &t );
    t0->SetBranchAddress("tp", &tp );
    t0->SetBranchAddress("t_Para", &t_Para );
    t0->SetBranchAddress("W", &W );
    t0->SetBranchAddress("x", &x );
    t0->SetBranchAddress("y", &y );
    t0->SetBranchAddress("z", &z );

    t0->SetBranchAddress("Qsq_cor", &Qsq_cor );
    t0->SetBranchAddress("t_cor", &t_cor);
    t0->SetBranchAddress("tp_cor", &tp_cor);
    t0->SetBranchAddress("W_cor", &W_cor);
    t0->SetBranchAddress("x_cor", &x_cor);
    t0->SetBranchAddress("y_cor", &y_cor);
    t0->SetBranchAddress("z_cor", &z_cor); 

    t0->SetBranchAddress("Vertex_X",   &Vertex_X );
    t0->SetBranchAddress("Vertex_Y",   &Vertex_Y );
    t0->SetBranchAddress("Vertex_Z",   &Vertex_Z );
    t0->SetBranchAddress("Theta_Pion_Photon",   &Theta_Pion_Photon );
    t0->SetBranchAddress("Photon_Theta",   &Photon_Theta );
    t0->SetBranchAddress("Photon_Factor",  &Photon_Factor);

    t0->SetBranchAddress("A_Factor",                                  &A_Factor        );
    t0->SetBranchAddress("Flux_Factor_RF",                            &Flux_Factor_RF  );
    t0->SetBranchAddress("Flux_Factor_Col",                           &Flux_Factor_Col );
    t0->SetBranchAddress("Jacobian_CM",                               &Jacobian_CM     );
    t0->SetBranchAddress("Jacobian_CM_RF",                            &Jacobian_CM_RF  );
    t0->SetBranchAddress("Jacobian_CM_Col",                           &Jacobian_CM_Col );
    /*}}}*/

    t0->SetBranchAddress("Phi",       &Phi     );/*{{{*/
    t0->SetBranchAddress("PhiS",      &PhiS    );
    t0->SetBranchAddress("Phi_cor",  &Phi_cor);
    t0->SetBranchAddress("PhiS_cor", &PhiS_cor);/*}}}*/

    t0->SetBranchAddress("Sigma_Lab",     &Sigma_Lab);                              /*{{{*/
    t0->SetBranchAddress("Sigma_UU",      &Sigma_UU);                              
    t0->SetBranchAddress("Sigma_UT",      &Sigma_UT);                              
    t0->SetBranchAddress("Sig_T",         &Sig_T);
    t0->SetBranchAddress("Sig_L",         &Sig_L);
    t0->SetBranchAddress("Sig_LT",        &Sig_LT);
    t0->SetBranchAddress("Sig_TT",        &Sig_TT);                                      /*}}}*/

    //Six Asymmetries and polarized XS/*{{{*/
    t0->SetBranchAddress("SSAsym",              &SSAsym  );
    t0->SetBranchAddress("SineAsym",            &SineAsym  );
    t0->SetBranchAddress("Asym_PhiS",           &Asym_PhiS  );
    t0->SetBranchAddress("Asym_PhiPlusPhiS",    &Asym_PhiPlusPhiS  );
    t0->SetBranchAddress("Asym_2PhiPlusPhiS",   &Asym_2PhiPlusPhiS  );
    t0->SetBranchAddress("Asym_PhiMinusPhiS",   &Asym_PhiMinusPhiS  );
    t0->SetBranchAddress("Asym_2PhiMinusPhiS",  &Asym_2PhiMinusPhiS  );
    t0->SetBranchAddress("Asym_3PhiMinusPhiS",  &Asym_3PhiMinusPhiS  );

    t0->SetBranchAddress("Sigma_PhiS",           &Sigma_PhiS  );
    t0->SetBranchAddress("Sigma_PhiPlusPhiS",    &Sigma_PhiPlusPhiS  );
    t0->SetBranchAddress("Sigma_2PhiPlusPhiS",   &Sigma_2PhiPlusPhiS  );
    t0->SetBranchAddress("Sigma_PhiMinusPhiS",   &Sigma_PhiMinusPhiS  );
    t0->SetBranchAddress("Sigma_2PhiMinusPhiS",  &Sigma_2PhiMinusPhiS  );
    t0->SetBranchAddress("Sigma_3PhiMinusPhiS",  &Sigma_3PhiMinusPhiS  );/*}}}*/

    t0->SetBranchAddress("EventWeight",         &EventWeight  );                              /*{{{*/
    t0->SetBranchAddress("WilliamsWeight",      &WilliamsWeight  );                              
    t0->SetBranchAddress("DedrickWeight",       &DedrickWeight  );                              
    t0->SetBranchAddress("CatchenWeight",       &CatchenWeight  );                              /*}}}*/

    t0->SetBranchAddress("tgt_px", &tgt_px);/*{{{*/
    t0->SetBranchAddress("tgt_py", &tgt_py);
    t0->SetBranchAddress("tgt_pz", &tgt_pz);
    t0->SetBranchAddress("tgt_theta", &tgt_theta);
    t0->SetBranchAddress("tgt_phi", &tgt_phi);
    t0->SetBranchAddress("tgt_ene", &tgt_ene);
    t0->SetBranchAddress("tgt_mom", &tgt_mom);/*}}}*/

    t0->SetBranchAddress("beam_px", &beam_px);/*{{{*/
    t0->SetBranchAddress("beam_py", &beam_py);
    t0->SetBranchAddress("beam_pz", &beam_pz);
    t0->SetBranchAddress("beam_theta", &beam_theta);
    t0->SetBranchAddress("beam_phi", &beam_phi);
    t0->SetBranchAddress("beam_ene", &beam_ene);
    t0->SetBranchAddress("beam_mom", &beam_mom);

    t0->SetBranchAddress("beam_px_cor", &beam_px_cor);
    t0->SetBranchAddress("beam_py_cor", &beam_py_cor);
    t0->SetBranchAddress("beam_pz_cor", &beam_pz_cor);
    t0->SetBranchAddress("beam_theta_cor", &beam_theta_cor);
    t0->SetBranchAddress("beam_phi_cor", &beam_phi_cor);
    t0->SetBranchAddress("beam_ene_cor", &beam_ene_cor);
    t0->SetBranchAddress("beam_mom_cor", &beam_mom_cor);
    /*}}}*/

    t0->SetBranchAddress("pim_ene",   &pim_ene   );/*{{{*/
    t0->SetBranchAddress("pim_px",    &pim_px    );
    t0->SetBranchAddress("pim_py",    &pim_py    );
    t0->SetBranchAddress("pim_pz",    &pim_pz    );
    t0->SetBranchAddress("pim_mom",   &pim_mom   );
    t0->SetBranchAddress("pim_theta", &pim_theta );
    t0->SetBranchAddress("pim_phi",   &pim_phi   );

    t0->SetBranchAddress("pim_ene_cor",   &pim_ene_cor  );
    t0->SetBranchAddress("pim_px_cor",    &pim_px_cor   );
    t0->SetBranchAddress("pim_py_cor",    &pim_py_cor   );
    t0->SetBranchAddress("pim_pz_cor",    &pim_pz_cor   );
    t0->SetBranchAddress("pim_mom_cor",   &pim_mom_cor  );
    t0->SetBranchAddress("pim_theta_cor", &pim_theta_cor);
    t0->SetBranchAddress("pim_phi_cor",   &pim_phi_cor  ); 
    
    t0->SetBranchAddress("pim_ene_res",   &pim_ene_res  );
    t0->SetBranchAddress("pim_px_res",    &pim_px_res   );
    t0->SetBranchAddress("pim_py_res",    &pim_py_res   );
    t0->SetBranchAddress("pim_pz_res",    &pim_pz_res   );
    t0->SetBranchAddress("pim_mom_res",   &pim_mom_res  );
    t0->SetBranchAddress("pim_theta_res", &pim_theta_res);
    t0->SetBranchAddress("pim_phi_res",   &pim_phi_res  );
    /*}}}*/

    t0->SetBranchAddress("ele_ene",   &ele_ene );/*{{{*/
    t0->SetBranchAddress("ele_px",    &ele_px  );
    t0->SetBranchAddress("ele_py",    &ele_py  );
    t0->SetBranchAddress("ele_pz",    &ele_pz  );
    t0->SetBranchAddress("ele_mom",   &ele_mom );
    t0->SetBranchAddress("ele_theta", &ele_theta);
    t0->SetBranchAddress("ele_phi",   &ele_phi );

    t0->SetBranchAddress("ele_ene_cor",   &ele_ene_cor  );
    t0->SetBranchAddress("ele_px_cor",    &ele_px_cor   );
    t0->SetBranchAddress("ele_py_cor",    &ele_py_cor   );
    t0->SetBranchAddress("ele_pz_cor",    &ele_pz_cor   );
    t0->SetBranchAddress("ele_mom_cor",   &ele_mom_cor  );
    t0->SetBranchAddress("ele_theta_cor", &ele_theta_cor);
    t0->SetBranchAddress("ele_phi_cor",   &ele_phi_cor  );
   
    t0->SetBranchAddress("ele_ene_res",   &ele_ene_res  );
    t0->SetBranchAddress("ele_px_res",    &ele_px_res   );
    t0->SetBranchAddress("ele_py_res",    &ele_py_res   );
    t0->SetBranchAddress("ele_pz_res",    &ele_pz_res   );
    t0->SetBranchAddress("ele_mom_res",   &ele_mom_res  );
    t0->SetBranchAddress("ele_theta_res", &ele_theta_res);
    t0->SetBranchAddress("ele_phi_res",   &ele_phi_res  );
    /*}}}*/

    t0->SetBranchAddress("pro_ene",   &pro_ene  ); /*{{{*/
    t0->SetBranchAddress("pro_px",    &pro_px   );
    t0->SetBranchAddress("pro_py",    &pro_py   );
    t0->SetBranchAddress("pro_pz",    &pro_pz   );
    t0->SetBranchAddress("pro_mom",   &pro_mom  );
    t0->SetBranchAddress("pro_theta", &pro_theta);
    t0->SetBranchAddress("pro_phi",   &pro_phi  );

    t0->SetBranchAddress("pro_ene_cor",   &pro_ene_cor   );
    t0->SetBranchAddress("pro_px_cor",    &pro_px_cor    );
    t0->SetBranchAddress("pro_py_cor",    &pro_py_cor    );
    t0->SetBranchAddress("pro_pz_cor",    &pro_pz_cor    );
    t0->SetBranchAddress("pro_mom_cor",   &pro_mom_cor   );
    t0->SetBranchAddress("pro_theta_cor", &pro_theta_cor );
    t0->SetBranchAddress("pro_phi_cor",   &pro_phi_cor   );
 
    t0->SetBranchAddress("pro_ene_res",   &pro_ene_res   );
    t0->SetBranchAddress("pro_px_res",    &pro_px_res    );
    t0->SetBranchAddress("pro_py_res",    &pro_py_res    );
    t0->SetBranchAddress("pro_pz_res",    &pro_pz_res    );
    t0->SetBranchAddress("pro_mom_res",   &pro_mom_res   );
    t0->SetBranchAddress("pro_theta_res", &pro_theta_res );
    t0->SetBranchAddress("pro_phi_res",   &pro_phi_res   );/*}}}*/

    //Add SoLID acceptance/*{{{*/
    t0->SetBranchAddress("ele_acc_f",     &ele_acc_f);
    t0->SetBranchAddress("ele_acc_l",     &ele_acc_l);
    t0->SetBranchAddress("pim_acc_f",     &pim_acc_f);
    t0->SetBranchAddress("pim_acc_l",     &pim_acc_l);
    t0->SetBranchAddress("pro_acc_f",     &pro_acc_f);
    t0->SetBranchAddress("pro_acc_l",     &pro_acc_l);/*}}}*/

    //Add detector resolutions/*{{{*/
    t0->SetBranchAddress("ele_mom_res",   &ele_mom_res  );
    t0->SetBranchAddress("ele_theta_res", &ele_theta_res);
    t0->SetBranchAddress("ele_phi_res",   &ele_phi_res  );
    t0->SetBranchAddress("pim_mom_res",   &pim_mom_res  );
    t0->SetBranchAddress("pim_theta_res", &pim_theta_res);
    t0->SetBranchAddress("pim_phi_res",   &pim_phi_res  );
    t0->SetBranchAddress("pro_mom_res",   &pro_mom_res  );
    t0->SetBranchAddress("pro_theta_res", &pro_theta_res);
    t0->SetBranchAddress("pro_phi_res",   &pro_phi_res  );/*}}}*/

    //Add other quantities/*{{{*/
    t0->SetBranchAddress("weight",            &weight        );
    t0->SetBranchAddress("weight_uu",         &weight_uu     ); //weight for unpolarized XS
    t0->SetBranchAddress("weight_ut",         &weight_ut     ); //weight for polarized XS
    t0->SetBranchAddress("weight_3m1",        &weight_3m1    ); //weight for Sin(3Phi-PhiS) module
    t0->SetBranchAddress("weight_2m1",        &weight_2m1    ); //weight for Sin(2Phi-PhiS) module
    t0->SetBranchAddress("weight_1m1",        &weight_1m1    ); //weight for Sin(Phi-PhiS) module
    t0->SetBranchAddress("weight_0p1",        &weight_0p1    ); //weight for Sin(PhiS) module
    t0->SetBranchAddress("weight_1p1",        &weight_1p1    ); //weight for Sin(Phi+PhiS) module
    t0->SetBranchAddress("weight_2p1",        &weight_2p1    ); //weight for Sin(2Phi+PhiS) module
    t0->SetBranchAddress("dilute",            &dilute        );
    t0->SetBranchAddress("MM",     &MP);
    t0->SetBranchAddress("MM_res", &MM_res);
    t0->SetBranchAddress("MM_cor", &MM_cor);
    t0->SetBranchAddress("MP",     &MP);
    t0->SetBranchAddress("MP_cor", &MP_cor);
    t0->SetBranchAddress("MP_res", &MP_res);
    t0->SetBranchAddress("Lumi_PSF", &Lumi_PSF);
 
    t0->SetBranchAddress("time",            &time);
    t0->SetBranchAddress("fileNO",          &fileNO);
    t0->SetBranchAddress("total_acc",       &total_acc);
    t0->SetBranchAddress("total_acc_cor",  &total_acc_cor);
    t0->SetBranchAddress("total_acc_res",   &total_acc_res);
    /*}}}*/
/*}}}*/
    
    /*Vectors & Histograms{{{*/
    TFile* histo =new TFile(histoname.Data(),"recreate"); 
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
        t0->GetEntry(i);
        
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
            if(MP>0.&&MP<MP_CUT&&MV->Theta()*Rad2Deg>8 && MV->Theta()*Rad2Deg<24.0)
                hMM_dvmp_cut->Fill(MM, weight/time*total_acc);
            /*}}}*/

            /*Missing w/o resolutions + w/ nuclear effect{{{*/
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
            if(MP>0.&&MP<MP_CUT&&MV->Theta()*Rad2Deg>8 && MV->Theta()*Rad2Deg<24.0)
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
            if(MP_res>0.&&MP_res<MP_CUT&&MV->Theta()*Rad2Deg>8 && MV->Theta()*Rad2Deg<24.0)
                hMM_dvmp_cut_res->Fill(MM_res, weight/time*total_acc_res);
            /*}}}*/

            if(!(i%1000))
                cout<<"--- Processed events = "<<i<<"\r";
        }/*Q2>4 cut}}}*/

    }// events loop ends here
    /*}}}*/

    TString outf_name = Form("rate_%s_%s.txt", pol_name.Data(), type_name.Data());
    cout<<Form("--- Total DVMP Rate= %f Hz (%f Hz at Q2>4)",total_rate_dvmp,total_rate_dvmp_4GeV)<<endl;
    ofstream outf(outf_name.Data());
    outf<<Form("--- Total DVMP Rate= %f Hz (%f Hz at Q2>4)",total_rate_dvmp,total_rate_dvmp_4GeV)<<endl;

    histo->cd();
    hMM_dvmp->Write(); hMM_dvmp_cor->Write(); hMM_dvmp_res->Write(); 
    hMM_dvmp_cut->Write(); hMM_dvmp_cut_cor->Write(); hMM_dvmp_cut_res->Write(); 
    hMP_dvmp->Write(); hMP_dvmp_cor->Write(); hMP_dvmp_res->Write(); 
    hMT_dvmp->Write(); hMT_dvmp_cor->Write(); hMT_dvmp_res->Write(); 

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
