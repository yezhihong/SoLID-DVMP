/*C/C++ Includes{{{*/
#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <map>
#include <cmath>
//#include "Rtypes.h"
//#include "math.h"

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
#include <TH1D.h>
#include <TH2D.h>
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
#include <TRandom3.h>
#include <TApplication.h>
#include <Rtypes.h>
#include <TTree.h>
//#include <TMatrix.h>
/*}}}*/

using namespace std;
const Double_t PI = 3.1415926;
const Double_t Deg2Rad = PI/180.0;
Int_t makebin();

/*Int_t main{{{*/
Int_t main(){
    Int_t err = -1000;
    err = makebin();
    if(err<0){
        cerr<<"***** Return error code = "<<err <<endl;
    }
    return err;
}
/*}}}*/

Int_t makebin(){
    gStyle->SetOptStat(0);

    TString pol_name="";
    int pol=0;
    cout<<"--- Pol = (1->up, 2->down) "; cin >> pol;
    if(pol==1) pol_name="up";
    if(pol==2) pol_name="down";

    TString prefix = "";
    TString finalfile = Form("./rootfiles/DEMP_Ee_11_4_%s_simple.root", pol_name.Data());
    TFile *file = new TFile(finalfile.Data(),"r");
    TTree *t0 = (TTree*) file->Get("T");
    int N_entries = t0->GetEntries();
   
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
    Double_t Flux_Factor_RF, Flux_Factor_Col, Jacobian_CM, Jacobian_CM_RF, Jacobian_CM_Col, A_Factor, time;
    Int_t NRecorded, NGenerated, fileNO; 

    double ele_acc_f, ele_acc_l, pim_acc_f,pim_acc_l, pro_acc_f,pro_acc_l,total_acceptance;
    double MM, MM_res,dilute,weight, weight_uu, weight_ut, weight_3m1, weight_2m1, weight_1m1, weight_0p1, weight_1p1, weight_2p1;
    double ele_mom_res, ele_theta_res, ele_phi_res;
    double pim_mom_res, pim_theta_res, pim_phi_res;
    double pro_mom_res, pro_theta_res, pro_phi_res;
    Int_t Q2BIN;
    /*}}}*/
    
    /*Define old Tree and old Branch{{{*/
    t0->SetBranchAddress("NRecorded",       &NRecorded);/*{{{*/
    t0->SetBranchAddress("NGenerated",       &NGenerated);

    t0->SetBranchAddress("Epsilon", &Epsilon);
    t0->SetBranchAddress("Qsq", &Qsq );
    t0->SetBranchAddress("t", &t );
    t0->SetBranchAddress("W", &W );
    t0->SetBranchAddress("x", &x );
    t0->SetBranchAddress("y", &y );
    t0->SetBranchAddress("z", &z );

    t0->SetBranchAddress("Qsq_corr", &Qsq_corr );
    t0->SetBranchAddress("t_corr", &t_corr);
    t0->SetBranchAddress("W_corr", &W_corr);
    t0->SetBranchAddress("x_corr", &x_corr);
    t0->SetBranchAddress("y_corr", &y_corr);
    t0->SetBranchAddress("z_corr", &z_corr); 

    t0->SetBranchAddress("Vertex_X",   &Vertex_X );
    t0->SetBranchAddress("Vertex_Y",   &Vertex_Y );
    t0->SetBranchAddress("Vertex_Z",   &Vertex_Z );
    t0->SetBranchAddress("Theta_Pion_Photon",   &Theta_Pion_Photon );

    t0->SetBranchAddress("A_Factor",                                  &A_Factor        );
    t0->SetBranchAddress("Flux_Factor_RF",                            &Flux_Factor_RF  );
    t0->SetBranchAddress("Flux_Factor_Col",                           &Flux_Factor_Col );
    t0->SetBranchAddress("Jacobian_CM",                               &Jacobian_CM     );
    t0->SetBranchAddress("Jacobian_CM_RF",                            &Jacobian_CM_RF  );
    t0->SetBranchAddress("Jacobian_CM_Col",                           &Jacobian_CM_Col );
    /*}}}*/

    t0->SetBranchAddress("Phi",       &Phi     );/*{{{*/
    t0->SetBranchAddress("PhiS",      &PhiS    );
    t0->SetBranchAddress("Phi_corr",  &Phi_corr);
    t0->SetBranchAddress("PhiS_corr", &PhiS_corr);/*}}}*/

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
    t0->SetBranchAddress("tgt_ene", &tgt_mom);/*}}}*/

    t0->SetBranchAddress("beam_px", &beam_px);/*{{{*/
    t0->SetBranchAddress("beam_py", &beam_py);
    t0->SetBranchAddress("beam_pz", &beam_pz);
    t0->SetBranchAddress("beam_theta", &beam_theta);
    t0->SetBranchAddress("beam_phi", &beam_phi);
    t0->SetBranchAddress("beam_ene", &beam_ene);
    t0->SetBranchAddress("beam_ene", &beam_mom);

    t0->SetBranchAddress("beam_corr_px", &beam_corr_px);
    t0->SetBranchAddress("beam_corr_py", &beam_corr_py);
    t0->SetBranchAddress("beam_corr_pz", &beam_corr_pz);
    t0->SetBranchAddress("beam_corr_theta", &beam_corr_theta);
    t0->SetBranchAddress("beam_corr_phi", &beam_corr_phi);
    t0->SetBranchAddress("beam_corr_ene", &beam_corr_ene);
    t0->SetBranchAddress("beam_corr_ene", &beam_corr_mom);/*}}}*/

    t0->SetBranchAddress("pim_ene",   &pim_ene   );/*{{{*/
    t0->SetBranchAddress("pim_px",    &pim_px    );
    t0->SetBranchAddress("pim_py",    &pim_py    );
    t0->SetBranchAddress("pim_pz",    &pim_pz    );
    t0->SetBranchAddress("pim_mom",   &pim_mom   );
    t0->SetBranchAddress("pim_theta", &pim_theta );
    t0->SetBranchAddress("pim_phi",   &pim_phi   );

    t0->SetBranchAddress("pim_corr_ene",   &pim_corr_ene  );
    t0->SetBranchAddress("pim_corr_px",    &pim_corr_px   );
    t0->SetBranchAddress("pim_corr_py",    &pim_corr_py   );
    t0->SetBranchAddress("pim_corr_pz",    &pim_corr_pz   );
    t0->SetBranchAddress("pim_corr_mom",   &pim_corr_mom  );
    t0->SetBranchAddress("pim_corr_theta", &pim_corr_theta);
    t0->SetBranchAddress("pim_corr_phi",   &pim_corr_phi  );/*}}}*/

    t0->SetBranchAddress("ele_ene",   &ele_ene );/*{{{*/
    t0->SetBranchAddress("ele_px",    &ele_px  );
    t0->SetBranchAddress("ele_py",    &ele_py  );
    t0->SetBranchAddress("ele_pz",    &ele_pz  );
    t0->SetBranchAddress("ele_mom",   &ele_mom );
    t0->SetBranchAddress("ele_theta", &ele_theta);
    t0->SetBranchAddress("ele_phi",   &ele_phi );

    t0->SetBranchAddress("ele_corr_ene",   &ele_corr_ene  );
    t0->SetBranchAddress("ele_corr_px",    &ele_corr_px   );
    t0->SetBranchAddress("ele_corr_py",    &ele_corr_py   );
    t0->SetBranchAddress("ele_corr_pz",    &ele_corr_pz   );
    t0->SetBranchAddress("ele_corr_mom",   &ele_corr_mom  );
    t0->SetBranchAddress("ele_corr_theta", &ele_corr_theta);
    t0->SetBranchAddress("ele_corr_phi",   &ele_corr_phi  );/*}}}*/

    t0->SetBranchAddress("pro_ene",   &pro_ene  ); /*{{{*/
    t0->SetBranchAddress("pro_px",    &pro_px   );
    t0->SetBranchAddress("pro_py",    &pro_py   );
    t0->SetBranchAddress("pro_pz",    &pro_pz   );
    t0->SetBranchAddress("pro_mom",   &pro_mom  );
    t0->SetBranchAddress("pro_theta", &pro_theta);
    t0->SetBranchAddress("pro_phi",   &pro_phi  );

    t0->SetBranchAddress("pro_corr_ene",   &pro_corr_ene   );
    t0->SetBranchAddress("pro_corr_px",    &pro_corr_px    );
    t0->SetBranchAddress("pro_corr_py",    &pro_corr_py    );
    t0->SetBranchAddress("pro_corr_pz",    &pro_corr_pz    );
    t0->SetBranchAddress("pro_corr_mom",   &pro_corr_mom   );
    t0->SetBranchAddress("pro_corr_theta", &pro_corr_theta );
    t0->SetBranchAddress("pro_corr_phi",   &pro_corr_phi   );/*}}}*/

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
    t0->SetBranchAddress("MM",            &MM);
    t0->SetBranchAddress("MM_res", &MM_res);
    
    //t0->SetBranchAddress("time",            &time);
    //t0->SetBranchAddress("fileNO",            &fileNO);
    //t0->SetBranchAddress("total_acceptance",            &total_acceptance);[>}}}<]
    /*}}}*/

    /*Define #1 new Tree and new Branch{{{*/
    TFile *f1 = new TFile(Form("./rootfiles/dvmp_%s_t1_newt.root", pol_name.Data()), "recreate");
    TTree *t1 = new TTree("T","a new tree");

    t1->Branch("NRecorded",       &NRecorded,     "NRecorded/I");/*{{{*/
    t1->Branch("NGenerated",       &NGenerated,     "NGenerated/I");

    t1->Branch("Epsilon", &Epsilon, "Epsilon/D" );
    t1->Branch("Qsq", &Qsq ,"Qsq/D");
    t1->Branch("t", &t ,"t/D");
    t1->Branch("W", &W ,"W/D");
    t1->Branch("x", &x ,"x/D");
    t1->Branch("y", &y ,"y/D");
    t1->Branch("z", &z ,"z/D");

    t1->Branch("Qsq_corr", &Qsq_corr ,"Qsq_corr/D");
    t1->Branch("t_corr", &t_corr ,"t_corr/D");
    t1->Branch("W_corr", &W_corr ,"W_corr/D");
    t1->Branch("x_corr", &x_corr ,"x_corr/D");
    t1->Branch("y_corr", &y_corr ,"y_corr/D");
    t1->Branch("z_corr", &z_corr ,"z_corr/D"); 

    t1->Branch("Vertex_X",   &Vertex_X   ,"Vertex_X/D");
    t1->Branch("Vertex_Y",   &Vertex_Y   ,"Vertex_Y/D");
    t1->Branch("Vertex_Z",   &Vertex_Z   ,"Vertex_Z/D");
    t1->Branch("Theta_Pion_Photon",   &Theta_Pion_Photon   ,"Theta_Pion_Photon/D");

    t1->Branch("A_Factor",                                  &A_Factor          ,"data/D"    );
    t1->Branch("Flux_Factor_RF",                            &Flux_Factor_RF       ,"data/D"    );
    t1->Branch("Flux_Factor_Col",                           &Flux_Factor_Col   ,"data/D"    );
    t1->Branch("Jacobian_CM",                               &Jacobian_CM       ,"data/D"    );
    t1->Branch("Jacobian_CM_RF",                            &Jacobian_CM_RF    ,"data/D"    );
    t1->Branch("Jacobian_CM_Col",                           &Jacobian_CM_Col   ,"data/D"    );
    /*}}}*/

    t1->Branch("Phi",       &Phi,       "Phi/D");/*{{{*/
    t1->Branch("PhiS",      &PhiS,      "PhiS/D");
    t1->Branch("Phi_corr",  &Phi_corr,  "Phi_corr/D");
    t1->Branch("PhiS_corr", &PhiS_corr, "PhiS_corr/D");/*}}}*/

    t1->Branch("Sigma_Lab",     &Sigma_Lab, "data/D");                              /*{{{*/
    t1->Branch("Sigma_UU",      &Sigma_UU,   "data/D");                              
    t1->Branch("Sigma_UT",      &Sigma_UT,   "data/D");                              
    t1->Branch("Sig_T",         &Sig_T, "data/D");                                  
    t1->Branch("Sig_L",         &Sig_L, "data/D");                                  
    t1->Branch("Sig_LT",        &Sig_LT, "data/D");                                 
    t1->Branch("Sig_TT",        &Sig_TT, "data/D");                                      /*}}}*/

    //Six Asymmetries and polarized XS/*{{{*/
    t1->Branch("SSAsym",              &SSAsym, "data/D");
    t1->Branch("SineAsym",            &SineAsym, "data/D");
    t1->Branch("Asym_PhiS",           &Asym_PhiS, "data/D");
    t1->Branch("Asym_PhiPlusPhiS",    &Asym_PhiPlusPhiS, "data/D");
    t1->Branch("Asym_2PhiPlusPhiS",   &Asym_2PhiPlusPhiS, "data/D");
    t1->Branch("Asym_PhiMinusPhiS",   &Asym_PhiMinusPhiS, "data/D");
    t1->Branch("Asym_2PhiMinusPhiS",  &Asym_2PhiMinusPhiS, "data/D");
    t1->Branch("Asym_3PhiMinusPhiS",  &Asym_3PhiMinusPhiS, "data/D");

    t1->Branch("Sigma_PhiS",           &Sigma_PhiS, "data/D");
    t1->Branch("Sigma_PhiPlusPhiS",    &Sigma_PhiPlusPhiS, "data/D");
    t1->Branch("Sigma_2PhiPlusPhiS",   &Sigma_2PhiPlusPhiS, "data/D");
    t1->Branch("Sigma_PhiMinusPhiS",   &Sigma_PhiMinusPhiS, "data/D");
    t1->Branch("Sigma_2PhiMinusPhiS",  &Sigma_2PhiMinusPhiS, "data/D");
    t1->Branch("Sigma_3PhiMinusPhiS",  &Sigma_3PhiMinusPhiS, "data/D");/*}}}*/

    t1->Branch("EventWeight",         &EventWeight, "data/D");                              /*{{{*/
    t1->Branch("WilliamsWeight",      &WilliamsWeight, "data/D");                              
    t1->Branch("DedrickWeight",       &DedrickWeight, "data/D");                              
    t1->Branch("CatchenWeight",       &CatchenWeight, "data/D");                              /*}}}*/

    t1->Branch("tgt_px", &tgt_px, "tgt_px/D");/*{{{*/
    t1->Branch("tgt_py", &tgt_py, "tgt_py/D");
    t1->Branch("tgt_pz", &tgt_pz, "tgt_pz/D");
    t1->Branch("tgt_theta", &tgt_theta, "tgt_theta/D");
    t1->Branch("tgt_phi", &tgt_phi, "tgt_phi/D");
    t1->Branch("tgt_ene", &tgt_ene, "tgt_ene/D");
    t1->Branch("tgt_ene", &tgt_mom, "tgt_mom/D");/*}}}*/

    t1->Branch("beam_px", &beam_px, "beam_px/D");/*{{{*/
    t1->Branch("beam_py", &beam_py, "beam_py/D");
    t1->Branch("beam_pz", &beam_pz, "beam_pz/D");
    t1->Branch("beam_theta", &beam_theta, "beam_theta/D");
    t1->Branch("beam_phi", &beam_phi, "beam_phi/D");
    t1->Branch("beam_ene", &beam_ene, "beam_ene/D");
    t1->Branch("beam_ene", &beam_mom, "beam_mom/D");

    t1->Branch("beam_corr_px", &beam_corr_px, "beam_corr_px/D");
    t1->Branch("beam_corr_py", &beam_corr_py, "beam_corr_py/D");
    t1->Branch("beam_corr_pz", &beam_corr_pz, "beam_corr_pz/D");
    t1->Branch("beam_corr_theta", &beam_corr_theta, "beam_corr_theta/D");
    t1->Branch("beam_corr_phi", &beam_corr_phi, "beam_corr_phi/D");
    t1->Branch("beam_corr_ene", &beam_corr_ene, "beam_corr_ene/D");
    t1->Branch("beam_corr_ene", &beam_corr_mom, "beam_corr_mom/D");/*}}}*/

    t1->Branch("pim_ene",   &pim_ene   ,"pim_ene/D");/*{{{*/
    t1->Branch("pim_px",    &pim_px    ,"pim_px/D");
    t1->Branch("pim_py",    &pim_py    ,"pim_py/D");
    t1->Branch("pim_pz",    &pim_pz    ,"pim_pz/D");
    t1->Branch("pim_mom",   &pim_mom   ,"pim_mom/D");
    t1->Branch("pim_theta", &pim_theta ,"pim_theta/D");
    t1->Branch("pim_phi",   &pim_phi   ,"pim_phi/D");

    t1->Branch("pim_corr_ene",   &pim_corr_ene   ,"pim_corr_ene/D");
    t1->Branch("pim_corr_px",    &pim_corr_px    ,"pim_corr_px/D");
    t1->Branch("pim_corr_py",    &pim_corr_py    ,"pim_corr_py/D");
    t1->Branch("pim_corr_pz",    &pim_corr_pz    ,"pim_corr_pz/D");
    t1->Branch("pim_corr_mom",   &pim_corr_mom   ,"pim_corr_mom/D");
    t1->Branch("pim_corr_theta", &pim_corr_theta ,"pim_corr_theta/D");
    t1->Branch("pim_corr_phi",   &pim_corr_phi   ,"pim_corr_phi/D");/*}}}*/

    t1->Branch("ele_ene",   &ele_ene   ,"ele_ene/D");/*{{{*/
    t1->Branch("ele_px",    &ele_px    ,"ele_px/D");
    t1->Branch("ele_py",    &ele_py    ,"ele_py/D");
    t1->Branch("ele_pz",    &ele_pz    ,"ele_pz/D");
    t1->Branch("ele_mom",   &ele_mom   ,"ele_mom/D");
    t1->Branch("ele_theta", &ele_theta ,"ele_theta/D");
    t1->Branch("ele_phi",   &ele_phi   ,"ele_phi/D");

    t1->Branch("ele_corr_ene",   &ele_corr_ene   ,"ele_corr_ene/D");
    t1->Branch("ele_corr_px",    &ele_corr_px    ,"ele_corr_px/D");
    t1->Branch("ele_corr_py",    &ele_corr_py    ,"ele_corr_py/D");
    t1->Branch("ele_corr_pz",    &ele_corr_pz    ,"ele_corr_pz/D");
    t1->Branch("ele_corr_mom",   &ele_corr_mom   ,"ele_corr_mom/D");
    t1->Branch("ele_corr_theta", &ele_corr_theta ,"ele_corr_theta/D");
    t1->Branch("ele_corr_phi",   &ele_corr_phi   ,"ele_corr_phi/D");/*}}}*/

    t1->Branch("pro_ene",   &pro_ene   ,"pro_ene/D"); /*{{{*/
    t1->Branch("pro_px",    &pro_px    ,"pro_px/D");
    t1->Branch("pro_py",    &pro_py    ,"pro_py/D");
    t1->Branch("pro_pz",    &pro_pz    ,"pro_pz/D");
    t1->Branch("pro_mom",   &pro_mom   ,"pro_mom/D");
    t1->Branch("pro_theta", &pro_theta ,"pro_theta/D");
    t1->Branch("pro_phi",   &pro_phi   ,"pro_phi/D");

    t1->Branch("pro_corr_ene",   &pro_corr_ene   ,"pro_corr_ene/D");
    t1->Branch("pro_corr_px",    &pro_corr_px    ,"pro_corr_px/D");
    t1->Branch("pro_corr_py",    &pro_corr_py    ,"pro_corr_py/D");
    t1->Branch("pro_corr_pz",    &pro_corr_pz    ,"pro_corr_pz/D");
    t1->Branch("pro_corr_mom",   &pro_corr_mom   ,"pro_corr_mom/D");
    t1->Branch("pro_corr_theta", &pro_corr_theta ,"pro_corr_theta/D");
    t1->Branch("pro_corr_phi",   &pro_corr_phi   ,"pro_corr_phi/D");/*}}}*/

    //Add SoLID acceptance/*{{{*/
    t1->Branch("ele_acc_f",     &ele_acc_f,     "ele_acc_f/D");
    t1->Branch("ele_acc_l",     &ele_acc_l,     "ele_acc_l/D");
    t1->Branch("pim_acc_f",     &pim_acc_f,     "pim_acc_f/D");
    t1->Branch("pim_acc_l",     &pim_acc_l,     "pim_acc_l/D");
    t1->Branch("pro_acc_f",     &pro_acc_f,     "pro_acc_f/D");
    t1->Branch("pro_acc_l",     &pro_acc_l,     "pro_acc_l/D");/*}}}*/

    //Add detector resolutions/*{{{*/
    t1->Branch("ele_mom_res",   &ele_mom_res,   "ele_mom_res/D");
    t1->Branch("ele_theta_res", &ele_theta_res, "ele_theta_res/D");
    t1->Branch("ele_phi_res",   &ele_phi_res,   "ele_phi_res/D");
    t1->Branch("pim_mom_res",   &pim_mom_res,   "pim_mom_res/D");
    t1->Branch("pim_theta_res", &pim_theta_res, "pim_theta_res/D");
    t1->Branch("pim_phi_res",   &pim_phi_res,   "pim_phi_res/D");
    t1->Branch("pro_mom_res",   &pro_mom_res,   "pro_mom_res/D");
    t1->Branch("pro_theta_res", &pro_theta_res, "pro_theta_res/D");
    t1->Branch("pro_phi_res",   &pro_phi_res,   "pro_phi_res/D");/*}}}*/

    //Add other quantities/*{{{*/
    t1->Branch("weight",        &weight,        "weight/D");
    t1->Branch("weight_uu",        &weight_uu,        "weight_uu/D"); //weight for unpolarized XS
    t1->Branch("weight_ut",        &weight_ut,        "weight_ut/D"); //weight for polarized XS
    t1->Branch("weight_3m1",        &weight_3m1,        "weight_3m1/D"); //weight for Sin(3Phi-PhiS) module
    t1->Branch("weight_2m1",        &weight_2m1,        "weight_2m1/D"); //weight for Sin(2Phi-PhiS) module
    t1->Branch("weight_1m1",        &weight_1m1,        "weight_1m1/D"); //weight for Sin(Phi-PhiS) module
    t1->Branch("weight_0p1",        &weight_0p1,        "weight_0p1/D"); //weight for Sin(PhiS) module
    t1->Branch("weight_1p1",        &weight_1p1,        "weight_1p1/D"); //weight for Sin(Phi+PhiS) module
    t1->Branch("weight_2p1",        &weight_2p1,        "weight_2p1/D"); //weight for Sin(2Phi+PhiS) module
    t1->Branch("dilute",        &dilute,        "dilute/D");
    //t1->Branch("PSF",           &PSF,           "PSF/D");
    t1->Branch("MM",            &MM,            "MM/D");
    t1->Branch("MM_res", &MM_res, "MM_res/D");/*}}}*/
    
    t1->Branch("time",            &time,            "time/D");
    //t1->Branch("fileNO",            &fileNO,            "fileNO/I");
    t1->Branch("total_acceptance",            &total_acceptance,            "total_acceptance/D");
    t1->Branch("Q2BIN", &Q2BIN, "Q2BIN/I");
    //t1->Branch("",            &,            "/D");
    /*}}}*/

    /*Define #2 new Tree and new Branch{{{*/
    TFile *f2 = new TFile(Form("./rootfiles/dvmp_%s_t2_newt.root", pol_name.Data()), "recreate");
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
    t2->Branch("MM_res", &MM_res, "MM_res/D");/*}}}*/
    
    t2->Branch("time",            &time,            "time/D");
    //t2->Branch("fileNO",            &fileNO,            "fileNO/I");
    t2->Branch("total_acceptance",            &total_acceptance,            "total_acceptance/D");
    t2->Branch("Q2BIN", &Q2BIN, "Q2BIN/I");
    //t2->Branch("",            &,            "/D");
    /*}}}*/
    
    /*Define #3 new Tree and new Branch{{{*/
    TFile *f3 = new TFile(Form("./rootfiles/dvmp_%s_t3_newt.root", pol_name.Data()), "recreate");
    TTree *t3 = new TTree("T","a new tree");

    t3->Branch("NRecorded",       &NRecorded,     "NRecorded/I");/*{{{*/
    t3->Branch("NGenerated",       &NGenerated,     "NGenerated/I");

    t3->Branch("Epsilon", &Epsilon, "Epsilon/D" );
    t3->Branch("Qsq", &Qsq ,"Qsq/D");
    t3->Branch("t", &t ,"t/D");
    t3->Branch("W", &W ,"W/D");
    t3->Branch("x", &x ,"x/D");
    t3->Branch("y", &y ,"y/D");
    t3->Branch("z", &z ,"z/D");

    t3->Branch("Qsq_corr", &Qsq_corr ,"Qsq_corr/D");
    t3->Branch("t_corr", &t_corr ,"t_corr/D");
    t3->Branch("W_corr", &W_corr ,"W_corr/D");
    t3->Branch("x_corr", &x_corr ,"x_corr/D");
    t3->Branch("y_corr", &y_corr ,"y_corr/D");
    t3->Branch("z_corr", &z_corr ,"z_corr/D"); 

    t3->Branch("Vertex_X",   &Vertex_X   ,"Vertex_X/D");
    t3->Branch("Vertex_Y",   &Vertex_Y   ,"Vertex_Y/D");
    t3->Branch("Vertex_Z",   &Vertex_Z   ,"Vertex_Z/D");
    t3->Branch("Theta_Pion_Photon",   &Theta_Pion_Photon   ,"Theta_Pion_Photon/D");

    t3->Branch("A_Factor",                                  &A_Factor          ,"data/D"    );
    t3->Branch("Flux_Factor_RF",                            &Flux_Factor_RF       ,"data/D"    );
    t3->Branch("Flux_Factor_Col",                           &Flux_Factor_Col   ,"data/D"    );
    t3->Branch("Jacobian_CM",                               &Jacobian_CM       ,"data/D"    );
    t3->Branch("Jacobian_CM_RF",                            &Jacobian_CM_RF    ,"data/D"    );
    t3->Branch("Jacobian_CM_Col",                           &Jacobian_CM_Col   ,"data/D"    );
    /*}}}*/

    t3->Branch("Phi",       &Phi,       "Phi/D");/*{{{*/
    t3->Branch("PhiS",      &PhiS,      "PhiS/D");
    t3->Branch("Phi_corr",  &Phi_corr,  "Phi_corr/D");
    t3->Branch("PhiS_corr", &PhiS_corr, "PhiS_corr/D");/*}}}*/

    t3->Branch("Sigma_Lab",     &Sigma_Lab, "data/D");                              /*{{{*/
    t3->Branch("Sigma_UU",      &Sigma_UU,   "data/D");                              
    t3->Branch("Sigma_UT",      &Sigma_UT,   "data/D");                              
    t3->Branch("Sig_T",         &Sig_T, "data/D");                                  
    t3->Branch("Sig_L",         &Sig_L, "data/D");                                  
    t3->Branch("Sig_LT",        &Sig_LT, "data/D");                                 
    t3->Branch("Sig_TT",        &Sig_TT, "data/D");                                      /*}}}*/

    //Six Asymmetries and polarized XS/*{{{*/
    t3->Branch("SSAsym",              &SSAsym, "data/D");
    t3->Branch("SineAsym",            &SineAsym, "data/D");
    t3->Branch("Asym_PhiS",           &Asym_PhiS, "data/D");
    t3->Branch("Asym_PhiPlusPhiS",    &Asym_PhiPlusPhiS, "data/D");
    t3->Branch("Asym_2PhiPlusPhiS",   &Asym_2PhiPlusPhiS, "data/D");
    t3->Branch("Asym_PhiMinusPhiS",   &Asym_PhiMinusPhiS, "data/D");
    t3->Branch("Asym_2PhiMinusPhiS",  &Asym_2PhiMinusPhiS, "data/D");
    t3->Branch("Asym_3PhiMinusPhiS",  &Asym_3PhiMinusPhiS, "data/D");

    t3->Branch("Sigma_PhiS",           &Sigma_PhiS, "data/D");
    t3->Branch("Sigma_PhiPlusPhiS",    &Sigma_PhiPlusPhiS, "data/D");
    t3->Branch("Sigma_2PhiPlusPhiS",   &Sigma_2PhiPlusPhiS, "data/D");
    t3->Branch("Sigma_PhiMinusPhiS",   &Sigma_PhiMinusPhiS, "data/D");
    t3->Branch("Sigma_2PhiMinusPhiS",  &Sigma_2PhiMinusPhiS, "data/D");
    t3->Branch("Sigma_3PhiMinusPhiS",  &Sigma_3PhiMinusPhiS, "data/D");/*}}}*/

    t3->Branch("EventWeight",         &EventWeight, "data/D");                              /*{{{*/
    t3->Branch("WilliamsWeight",      &WilliamsWeight, "data/D");                              
    t3->Branch("DedrickWeight",       &DedrickWeight, "data/D");                              
    t3->Branch("CatchenWeight",       &CatchenWeight, "data/D");                              /*}}}*/

    t3->Branch("tgt_px", &tgt_px, "tgt_px/D");/*{{{*/
    t3->Branch("tgt_py", &tgt_py, "tgt_py/D");
    t3->Branch("tgt_pz", &tgt_pz, "tgt_pz/D");
    t3->Branch("tgt_theta", &tgt_theta, "tgt_theta/D");
    t3->Branch("tgt_phi", &tgt_phi, "tgt_phi/D");
    t3->Branch("tgt_ene", &tgt_ene, "tgt_ene/D");
    t3->Branch("tgt_ene", &tgt_mom, "tgt_mom/D");/*}}}*/

    t3->Branch("beam_px", &beam_px, "beam_px/D");/*{{{*/
    t3->Branch("beam_py", &beam_py, "beam_py/D");
    t3->Branch("beam_pz", &beam_pz, "beam_pz/D");
    t3->Branch("beam_theta", &beam_theta, "beam_theta/D");
    t3->Branch("beam_phi", &beam_phi, "beam_phi/D");
    t3->Branch("beam_ene", &beam_ene, "beam_ene/D");
    t3->Branch("beam_ene", &beam_mom, "beam_mom/D");

    t3->Branch("beam_corr_px", &beam_corr_px, "beam_corr_px/D");
    t3->Branch("beam_corr_py", &beam_corr_py, "beam_corr_py/D");
    t3->Branch("beam_corr_pz", &beam_corr_pz, "beam_corr_pz/D");
    t3->Branch("beam_corr_theta", &beam_corr_theta, "beam_corr_theta/D");
    t3->Branch("beam_corr_phi", &beam_corr_phi, "beam_corr_phi/D");
    t3->Branch("beam_corr_ene", &beam_corr_ene, "beam_corr_ene/D");
    t3->Branch("beam_corr_ene", &beam_corr_mom, "beam_corr_mom/D");/*}}}*/

    t3->Branch("pim_ene",   &pim_ene   ,"pim_ene/D");/*{{{*/
    t3->Branch("pim_px",    &pim_px    ,"pim_px/D");
    t3->Branch("pim_py",    &pim_py    ,"pim_py/D");
    t3->Branch("pim_pz",    &pim_pz    ,"pim_pz/D");
    t3->Branch("pim_mom",   &pim_mom   ,"pim_mom/D");
    t3->Branch("pim_theta", &pim_theta ,"pim_theta/D");
    t3->Branch("pim_phi",   &pim_phi   ,"pim_phi/D");

    t3->Branch("pim_corr_ene",   &pim_corr_ene   ,"pim_corr_ene/D");
    t3->Branch("pim_corr_px",    &pim_corr_px    ,"pim_corr_px/D");
    t3->Branch("pim_corr_py",    &pim_corr_py    ,"pim_corr_py/D");
    t3->Branch("pim_corr_pz",    &pim_corr_pz    ,"pim_corr_pz/D");
    t3->Branch("pim_corr_mom",   &pim_corr_mom   ,"pim_corr_mom/D");
    t3->Branch("pim_corr_theta", &pim_corr_theta ,"pim_corr_theta/D");
    t3->Branch("pim_corr_phi",   &pim_corr_phi   ,"pim_corr_phi/D");/*}}}*/

    t3->Branch("ele_ene",   &ele_ene   ,"ele_ene/D");/*{{{*/
    t3->Branch("ele_px",    &ele_px    ,"ele_px/D");
    t3->Branch("ele_py",    &ele_py    ,"ele_py/D");
    t3->Branch("ele_pz",    &ele_pz    ,"ele_pz/D");
    t3->Branch("ele_mom",   &ele_mom   ,"ele_mom/D");
    t3->Branch("ele_theta", &ele_theta ,"ele_theta/D");
    t3->Branch("ele_phi",   &ele_phi   ,"ele_phi/D");

    t3->Branch("ele_corr_ene",   &ele_corr_ene   ,"ele_corr_ene/D");
    t3->Branch("ele_corr_px",    &ele_corr_px    ,"ele_corr_px/D");
    t3->Branch("ele_corr_py",    &ele_corr_py    ,"ele_corr_py/D");
    t3->Branch("ele_corr_pz",    &ele_corr_pz    ,"ele_corr_pz/D");
    t3->Branch("ele_corr_mom",   &ele_corr_mom   ,"ele_corr_mom/D");
    t3->Branch("ele_corr_theta", &ele_corr_theta ,"ele_corr_theta/D");
    t3->Branch("ele_corr_phi",   &ele_corr_phi   ,"ele_corr_phi/D");/*}}}*/

    t3->Branch("pro_ene",   &pro_ene   ,"pro_ene/D"); /*{{{*/
    t3->Branch("pro_px",    &pro_px    ,"pro_px/D");
    t3->Branch("pro_py",    &pro_py    ,"pro_py/D");
    t3->Branch("pro_pz",    &pro_pz    ,"pro_pz/D");
    t3->Branch("pro_mom",   &pro_mom   ,"pro_mom/D");
    t3->Branch("pro_theta", &pro_theta ,"pro_theta/D");
    t3->Branch("pro_phi",   &pro_phi   ,"pro_phi/D");

    t3->Branch("pro_corr_ene",   &pro_corr_ene   ,"pro_corr_ene/D");
    t3->Branch("pro_corr_px",    &pro_corr_px    ,"pro_corr_px/D");
    t3->Branch("pro_corr_py",    &pro_corr_py    ,"pro_corr_py/D");
    t3->Branch("pro_corr_pz",    &pro_corr_pz    ,"pro_corr_pz/D");
    t3->Branch("pro_corr_mom",   &pro_corr_mom   ,"pro_corr_mom/D");
    t3->Branch("pro_corr_theta", &pro_corr_theta ,"pro_corr_theta/D");
    t3->Branch("pro_corr_phi",   &pro_corr_phi   ,"pro_corr_phi/D");/*}}}*/

    //Add SoLID acceptance/*{{{*/
    t3->Branch("ele_acc_f",     &ele_acc_f,     "ele_acc_f/D");
    t3->Branch("ele_acc_l",     &ele_acc_l,     "ele_acc_l/D");
    t3->Branch("pim_acc_f",     &pim_acc_f,     "pim_acc_f/D");
    t3->Branch("pim_acc_l",     &pim_acc_l,     "pim_acc_l/D");
    t3->Branch("pro_acc_f",     &pro_acc_f,     "pro_acc_f/D");
    t3->Branch("pro_acc_l",     &pro_acc_l,     "pro_acc_l/D");/*}}}*/

    //Add detector resolutions/*{{{*/
    t3->Branch("ele_mom_res",   &ele_mom_res,   "ele_mom_res/D");
    t3->Branch("ele_theta_res", &ele_theta_res, "ele_theta_res/D");
    t3->Branch("ele_phi_res",   &ele_phi_res,   "ele_phi_res/D");
    t3->Branch("pim_mom_res",   &pim_mom_res,   "pim_mom_res/D");
    t3->Branch("pim_theta_res", &pim_theta_res, "pim_theta_res/D");
    t3->Branch("pim_phi_res",   &pim_phi_res,   "pim_phi_res/D");
    t3->Branch("pro_mom_res",   &pro_mom_res,   "pro_mom_res/D");
    t3->Branch("pro_theta_res", &pro_theta_res, "pro_theta_res/D");
    t3->Branch("pro_phi_res",   &pro_phi_res,   "pro_phi_res/D");/*}}}*/

    //Add other quantities/*{{{*/
    t3->Branch("weight",        &weight,        "weight/D");
    t3->Branch("weight_uu",        &weight_uu,        "weight_uu/D"); //weight for unpolarized XS
    t3->Branch("weight_ut",        &weight_ut,        "weight_ut/D"); //weight for polarized XS
    t3->Branch("weight_3m1",        &weight_3m1,        "weight_3m1/D"); //weight for Sin(3Phi-PhiS) module
    t3->Branch("weight_2m1",        &weight_2m1,        "weight_2m1/D"); //weight for Sin(2Phi-PhiS) module
    t3->Branch("weight_1m1",        &weight_1m1,        "weight_1m1/D"); //weight for Sin(Phi-PhiS) module
    t3->Branch("weight_0p1",        &weight_0p1,        "weight_0p1/D"); //weight for Sin(PhiS) module
    t3->Branch("weight_1p1",        &weight_1p1,        "weight_1p1/D"); //weight for Sin(Phi+PhiS) module
    t3->Branch("weight_2p1",        &weight_2p1,        "weight_2p1/D"); //weight for Sin(2Phi+PhiS) module
    t3->Branch("dilute",        &dilute,        "dilute/D");
    //t3->Branch("PSF",           &PSF,           "PSF/D");
    t3->Branch("MM",            &MM,            "MM/D");
    t3->Branch("MM_res", &MM_res, "MM_res/D");/*}}}*/
    
    t3->Branch("time",            &time,            "time/D");
    //t3->Branch("fileNO",            &fileNO,            "fileNO/I");
    t3->Branch("total_acceptance",            &total_acceptance,            "total_acceptance/D");
    t3->Branch("Q2BIN", &Q2BIN, "Q2BIN/I");
    //t3->Branch("",            &,            "/D");
    /*}}}*/

    /*Define #4 new Tree and new Branch{{{*/
    TFile *f4 = new TFile(Form("./rootfiles/dvmp_%s_t4_newt.root", pol_name.Data()), "recreate");
    TTree *t4 = new TTree("T","a new tree");

    t4->Branch("NRecorded",       &NRecorded,     "NRecorded/I");/*{{{*/
    t4->Branch("NGenerated",       &NGenerated,     "NGenerated/I");

    t4->Branch("Epsilon", &Epsilon, "Epsilon/D" );
    t4->Branch("Qsq", &Qsq ,"Qsq/D");
    t4->Branch("t", &t ,"t/D");
    t4->Branch("W", &W ,"W/D");
    t4->Branch("x", &x ,"x/D");
    t4->Branch("y", &y ,"y/D");
    t4->Branch("z", &z ,"z/D");

    t4->Branch("Qsq_corr", &Qsq_corr ,"Qsq_corr/D");
    t4->Branch("t_corr", &t_corr ,"t_corr/D");
    t4->Branch("W_corr", &W_corr ,"W_corr/D");
    t4->Branch("x_corr", &x_corr ,"x_corr/D");
    t4->Branch("y_corr", &y_corr ,"y_corr/D");
    t4->Branch("z_corr", &z_corr ,"z_corr/D"); 

    t4->Branch("Vertex_X",   &Vertex_X   ,"Vertex_X/D");
    t4->Branch("Vertex_Y",   &Vertex_Y   ,"Vertex_Y/D");
    t4->Branch("Vertex_Z",   &Vertex_Z   ,"Vertex_Z/D");
    t4->Branch("Theta_Pion_Photon",   &Theta_Pion_Photon   ,"Theta_Pion_Photon/D");

    t4->Branch("A_Factor",                                  &A_Factor          ,"data/D"    );
    t4->Branch("Flux_Factor_RF",                            &Flux_Factor_RF       ,"data/D"    );
    t4->Branch("Flux_Factor_Col",                           &Flux_Factor_Col   ,"data/D"    );
    t4->Branch("Jacobian_CM",                               &Jacobian_CM       ,"data/D"    );
    t4->Branch("Jacobian_CM_RF",                            &Jacobian_CM_RF    ,"data/D"    );
    t4->Branch("Jacobian_CM_Col",                           &Jacobian_CM_Col   ,"data/D"    );
    /*}}}*/

    t4->Branch("Phi",       &Phi,       "Phi/D");/*{{{*/
    t4->Branch("PhiS",      &PhiS,      "PhiS/D");
    t4->Branch("Phi_corr",  &Phi_corr,  "Phi_corr/D");
    t4->Branch("PhiS_corr", &PhiS_corr, "PhiS_corr/D");/*}}}*/

    t4->Branch("Sigma_Lab",     &Sigma_Lab, "data/D");                              /*{{{*/
    t4->Branch("Sigma_UU",      &Sigma_UU,   "data/D");                              
    t4->Branch("Sigma_UT",      &Sigma_UT,   "data/D");                              
    t4->Branch("Sig_T",         &Sig_T, "data/D");                                  
    t4->Branch("Sig_L",         &Sig_L, "data/D");                                  
    t4->Branch("Sig_LT",        &Sig_LT, "data/D");                                 
    t4->Branch("Sig_TT",        &Sig_TT, "data/D");                                      /*}}}*/

    //Six Asymmetries and polarized XS/*{{{*/
    t4->Branch("SSAsym",              &SSAsym, "data/D");
    t4->Branch("SineAsym",            &SineAsym, "data/D");
    t4->Branch("Asym_PhiS",           &Asym_PhiS, "data/D");
    t4->Branch("Asym_PhiPlusPhiS",    &Asym_PhiPlusPhiS, "data/D");
    t4->Branch("Asym_2PhiPlusPhiS",   &Asym_2PhiPlusPhiS, "data/D");
    t4->Branch("Asym_PhiMinusPhiS",   &Asym_PhiMinusPhiS, "data/D");
    t4->Branch("Asym_2PhiMinusPhiS",  &Asym_2PhiMinusPhiS, "data/D");
    t4->Branch("Asym_3PhiMinusPhiS",  &Asym_3PhiMinusPhiS, "data/D");

    t4->Branch("Sigma_PhiS",           &Sigma_PhiS, "data/D");
    t4->Branch("Sigma_PhiPlusPhiS",    &Sigma_PhiPlusPhiS, "data/D");
    t4->Branch("Sigma_2PhiPlusPhiS",   &Sigma_2PhiPlusPhiS, "data/D");
    t4->Branch("Sigma_PhiMinusPhiS",   &Sigma_PhiMinusPhiS, "data/D");
    t4->Branch("Sigma_2PhiMinusPhiS",  &Sigma_2PhiMinusPhiS, "data/D");
    t4->Branch("Sigma_3PhiMinusPhiS",  &Sigma_3PhiMinusPhiS, "data/D");/*}}}*/

    t4->Branch("EventWeight",         &EventWeight, "data/D");                              /*{{{*/
    t4->Branch("WilliamsWeight",      &WilliamsWeight, "data/D");                              
    t4->Branch("DedrickWeight",       &DedrickWeight, "data/D");                              
    t4->Branch("CatchenWeight",       &CatchenWeight, "data/D");                              /*}}}*/

    t4->Branch("tgt_px", &tgt_px, "tgt_px/D");/*{{{*/
    t4->Branch("tgt_py", &tgt_py, "tgt_py/D");
    t4->Branch("tgt_pz", &tgt_pz, "tgt_pz/D");
    t4->Branch("tgt_theta", &tgt_theta, "tgt_theta/D");
    t4->Branch("tgt_phi", &tgt_phi, "tgt_phi/D");
    t4->Branch("tgt_ene", &tgt_ene, "tgt_ene/D");
    t4->Branch("tgt_ene", &tgt_mom, "tgt_mom/D");/*}}}*/

    t4->Branch("beam_px", &beam_px, "beam_px/D");/*{{{*/
    t4->Branch("beam_py", &beam_py, "beam_py/D");
    t4->Branch("beam_pz", &beam_pz, "beam_pz/D");
    t4->Branch("beam_theta", &beam_theta, "beam_theta/D");
    t4->Branch("beam_phi", &beam_phi, "beam_phi/D");
    t4->Branch("beam_ene", &beam_ene, "beam_ene/D");
    t4->Branch("beam_ene", &beam_mom, "beam_mom/D");

    t4->Branch("beam_corr_px", &beam_corr_px, "beam_corr_px/D");
    t4->Branch("beam_corr_py", &beam_corr_py, "beam_corr_py/D");
    t4->Branch("beam_corr_pz", &beam_corr_pz, "beam_corr_pz/D");
    t4->Branch("beam_corr_theta", &beam_corr_theta, "beam_corr_theta/D");
    t4->Branch("beam_corr_phi", &beam_corr_phi, "beam_corr_phi/D");
    t4->Branch("beam_corr_ene", &beam_corr_ene, "beam_corr_ene/D");
    t4->Branch("beam_corr_ene", &beam_corr_mom, "beam_corr_mom/D");/*}}}*/

    t4->Branch("pim_ene",   &pim_ene   ,"pim_ene/D");/*{{{*/
    t4->Branch("pim_px",    &pim_px    ,"pim_px/D");
    t4->Branch("pim_py",    &pim_py    ,"pim_py/D");
    t4->Branch("pim_pz",    &pim_pz    ,"pim_pz/D");
    t4->Branch("pim_mom",   &pim_mom   ,"pim_mom/D");
    t4->Branch("pim_theta", &pim_theta ,"pim_theta/D");
    t4->Branch("pim_phi",   &pim_phi   ,"pim_phi/D");

    t4->Branch("pim_corr_ene",   &pim_corr_ene   ,"pim_corr_ene/D");
    t4->Branch("pim_corr_px",    &pim_corr_px    ,"pim_corr_px/D");
    t4->Branch("pim_corr_py",    &pim_corr_py    ,"pim_corr_py/D");
    t4->Branch("pim_corr_pz",    &pim_corr_pz    ,"pim_corr_pz/D");
    t4->Branch("pim_corr_mom",   &pim_corr_mom   ,"pim_corr_mom/D");
    t4->Branch("pim_corr_theta", &pim_corr_theta ,"pim_corr_theta/D");
    t4->Branch("pim_corr_phi",   &pim_corr_phi   ,"pim_corr_phi/D");/*}}}*/

    t4->Branch("ele_ene",   &ele_ene   ,"ele_ene/D");/*{{{*/
    t4->Branch("ele_px",    &ele_px    ,"ele_px/D");
    t4->Branch("ele_py",    &ele_py    ,"ele_py/D");
    t4->Branch("ele_pz",    &ele_pz    ,"ele_pz/D");
    t4->Branch("ele_mom",   &ele_mom   ,"ele_mom/D");
    t4->Branch("ele_theta", &ele_theta ,"ele_theta/D");
    t4->Branch("ele_phi",   &ele_phi   ,"ele_phi/D");

    t4->Branch("ele_corr_ene",   &ele_corr_ene   ,"ele_corr_ene/D");
    t4->Branch("ele_corr_px",    &ele_corr_px    ,"ele_corr_px/D");
    t4->Branch("ele_corr_py",    &ele_corr_py    ,"ele_corr_py/D");
    t4->Branch("ele_corr_pz",    &ele_corr_pz    ,"ele_corr_pz/D");
    t4->Branch("ele_corr_mom",   &ele_corr_mom   ,"ele_corr_mom/D");
    t4->Branch("ele_corr_theta", &ele_corr_theta ,"ele_corr_theta/D");
    t4->Branch("ele_corr_phi",   &ele_corr_phi   ,"ele_corr_phi/D");/*}}}*/

    t4->Branch("pro_ene",   &pro_ene   ,"pro_ene/D"); /*{{{*/
    t4->Branch("pro_px",    &pro_px    ,"pro_px/D");
    t4->Branch("pro_py",    &pro_py    ,"pro_py/D");
    t4->Branch("pro_pz",    &pro_pz    ,"pro_pz/D");
    t4->Branch("pro_mom",   &pro_mom   ,"pro_mom/D");
    t4->Branch("pro_theta", &pro_theta ,"pro_theta/D");
    t4->Branch("pro_phi",   &pro_phi   ,"pro_phi/D");

    t4->Branch("pro_corr_ene",   &pro_corr_ene   ,"pro_corr_ene/D");
    t4->Branch("pro_corr_px",    &pro_corr_px    ,"pro_corr_px/D");
    t4->Branch("pro_corr_py",    &pro_corr_py    ,"pro_corr_py/D");
    t4->Branch("pro_corr_pz",    &pro_corr_pz    ,"pro_corr_pz/D");
    t4->Branch("pro_corr_mom",   &pro_corr_mom   ,"pro_corr_mom/D");
    t4->Branch("pro_corr_theta", &pro_corr_theta ,"pro_corr_theta/D");
    t4->Branch("pro_corr_phi",   &pro_corr_phi   ,"pro_corr_phi/D");/*}}}*/

    //Add SoLID acceptance/*{{{*/
    t4->Branch("ele_acc_f",     &ele_acc_f,     "ele_acc_f/D");
    t4->Branch("ele_acc_l",     &ele_acc_l,     "ele_acc_l/D");
    t4->Branch("pim_acc_f",     &pim_acc_f,     "pim_acc_f/D");
    t4->Branch("pim_acc_l",     &pim_acc_l,     "pim_acc_l/D");
    t4->Branch("pro_acc_f",     &pro_acc_f,     "pro_acc_f/D");
    t4->Branch("pro_acc_l",     &pro_acc_l,     "pro_acc_l/D");/*}}}*/

    //Add detector resolutions/*{{{*/
    t4->Branch("ele_mom_res",   &ele_mom_res,   "ele_mom_res/D");
    t4->Branch("ele_theta_res", &ele_theta_res, "ele_theta_res/D");
    t4->Branch("ele_phi_res",   &ele_phi_res,   "ele_phi_res/D");
    t4->Branch("pim_mom_res",   &pim_mom_res,   "pim_mom_res/D");
    t4->Branch("pim_theta_res", &pim_theta_res, "pim_theta_res/D");
    t4->Branch("pim_phi_res",   &pim_phi_res,   "pim_phi_res/D");
    t4->Branch("pro_mom_res",   &pro_mom_res,   "pro_mom_res/D");
    t4->Branch("pro_theta_res", &pro_theta_res, "pro_theta_res/D");
    t4->Branch("pro_phi_res",   &pro_phi_res,   "pro_phi_res/D");/*}}}*/

    //Add other quantities/*{{{*/
    t4->Branch("weight",        &weight,        "weight/D");
    t4->Branch("weight_uu",        &weight_uu,        "weight_uu/D"); //weight for unpolarized XS
    t4->Branch("weight_ut",        &weight_ut,        "weight_ut/D"); //weight for polarized XS
    t4->Branch("weight_3m1",        &weight_3m1,        "weight_3m1/D"); //weight for Sin(3Phi-PhiS) module
    t4->Branch("weight_2m1",        &weight_2m1,        "weight_2m1/D"); //weight for Sin(2Phi-PhiS) module
    t4->Branch("weight_1m1",        &weight_1m1,        "weight_1m1/D"); //weight for Sin(Phi-PhiS) module
    t4->Branch("weight_0p1",        &weight_0p1,        "weight_0p1/D"); //weight for Sin(PhiS) module
    t4->Branch("weight_1p1",        &weight_1p1,        "weight_1p1/D"); //weight for Sin(Phi+PhiS) module
    t4->Branch("weight_2p1",        &weight_2p1,        "weight_2p1/D"); //weight for Sin(2Phi+PhiS) module
    t4->Branch("dilute",        &dilute,        "dilute/D");
    //t4->Branch("PSF",           &PSF,           "PSF/D");
    t4->Branch("MM",            &MM,            "MM/D");
    t4->Branch("MM_res", &MM_res, "MM_res/D");/*}}}*/
    
    t4->Branch("time",            &time,            "time/D");
    //t4->Branch("fileNO",            &fileNO,            "fileNO/I");
    t4->Branch("total_acceptance",            &total_acceptance,            "total_acceptance/D");
    t4->Branch("Q2BIN", &Q2BIN, "Q2BIN/I");
    //t4->Branch("",            &,            "/D");
    /*}}}*/
    
    /*Define #5 new Tree and new Branch{{{*/
    TFile *f5 = new TFile(Form("./rootfiles/dvmp_%s_t5_newt.root", pol_name.Data()), "recreate");
    TTree *t5 = new TTree("T","a new tree");

    t5->Branch("NRecorded",       &NRecorded,     "NRecorded/I");/*{{{*/
    t5->Branch("NGenerated",       &NGenerated,     "NGenerated/I");

    t5->Branch("Epsilon", &Epsilon, "Epsilon/D" );
    t5->Branch("Qsq", &Qsq ,"Qsq/D");
    t5->Branch("t", &t ,"t/D");
    t5->Branch("W", &W ,"W/D");
    t5->Branch("x", &x ,"x/D");
    t5->Branch("y", &y ,"y/D");
    t5->Branch("z", &z ,"z/D");

    t5->Branch("Qsq_corr", &Qsq_corr ,"Qsq_corr/D");
    t5->Branch("t_corr", &t_corr ,"t_corr/D");
    t5->Branch("W_corr", &W_corr ,"W_corr/D");
    t5->Branch("x_corr", &x_corr ,"x_corr/D");
    t5->Branch("y_corr", &y_corr ,"y_corr/D");
    t5->Branch("z_corr", &z_corr ,"z_corr/D"); 

    t5->Branch("Vertex_X",   &Vertex_X   ,"Vertex_X/D");
    t5->Branch("Vertex_Y",   &Vertex_Y   ,"Vertex_Y/D");
    t5->Branch("Vertex_Z",   &Vertex_Z   ,"Vertex_Z/D");
    t5->Branch("Theta_Pion_Photon",   &Theta_Pion_Photon   ,"Theta_Pion_Photon/D");

    t5->Branch("A_Factor",                                  &A_Factor          ,"data/D"    );
    t5->Branch("Flux_Factor_RF",                            &Flux_Factor_RF       ,"data/D"    );
    t5->Branch("Flux_Factor_Col",                           &Flux_Factor_Col   ,"data/D"    );
    t5->Branch("Jacobian_CM",                               &Jacobian_CM       ,"data/D"    );
    t5->Branch("Jacobian_CM_RF",                            &Jacobian_CM_RF    ,"data/D"    );
    t5->Branch("Jacobian_CM_Col",                           &Jacobian_CM_Col   ,"data/D"    );
    /*}}}*/

    t5->Branch("Phi",       &Phi,       "Phi/D");/*{{{*/
    t5->Branch("PhiS",      &PhiS,      "PhiS/D");
    t5->Branch("Phi_corr",  &Phi_corr,  "Phi_corr/D");
    t5->Branch("PhiS_corr", &PhiS_corr, "PhiS_corr/D");/*}}}*/

    t5->Branch("Sigma_Lab",     &Sigma_Lab, "data/D");                              /*{{{*/
    t5->Branch("Sigma_UU",      &Sigma_UU,   "data/D");                              
    t5->Branch("Sigma_UT",      &Sigma_UT,   "data/D");                              
    t5->Branch("Sig_T",         &Sig_T, "data/D");                                  
    t5->Branch("Sig_L",         &Sig_L, "data/D");                                  
    t5->Branch("Sig_LT",        &Sig_LT, "data/D");                                 
    t5->Branch("Sig_TT",        &Sig_TT, "data/D");                                      /*}}}*/

    //Six Asymmetries and polarized XS/*{{{*/
    t5->Branch("SSAsym",              &SSAsym, "data/D");
    t5->Branch("SineAsym",            &SineAsym, "data/D");
    t5->Branch("Asym_PhiS",           &Asym_PhiS, "data/D");
    t5->Branch("Asym_PhiPlusPhiS",    &Asym_PhiPlusPhiS, "data/D");
    t5->Branch("Asym_2PhiPlusPhiS",   &Asym_2PhiPlusPhiS, "data/D");
    t5->Branch("Asym_PhiMinusPhiS",   &Asym_PhiMinusPhiS, "data/D");
    t5->Branch("Asym_2PhiMinusPhiS",  &Asym_2PhiMinusPhiS, "data/D");
    t5->Branch("Asym_3PhiMinusPhiS",  &Asym_3PhiMinusPhiS, "data/D");

    t5->Branch("Sigma_PhiS",           &Sigma_PhiS, "data/D");
    t5->Branch("Sigma_PhiPlusPhiS",    &Sigma_PhiPlusPhiS, "data/D");
    t5->Branch("Sigma_2PhiPlusPhiS",   &Sigma_2PhiPlusPhiS, "data/D");
    t5->Branch("Sigma_PhiMinusPhiS",   &Sigma_PhiMinusPhiS, "data/D");
    t5->Branch("Sigma_2PhiMinusPhiS",  &Sigma_2PhiMinusPhiS, "data/D");
    t5->Branch("Sigma_3PhiMinusPhiS",  &Sigma_3PhiMinusPhiS, "data/D");/*}}}*/

    t5->Branch("EventWeight",         &EventWeight, "data/D");                              /*{{{*/
    t5->Branch("WilliamsWeight",      &WilliamsWeight, "data/D");                              
    t5->Branch("DedrickWeight",       &DedrickWeight, "data/D");                              
    t5->Branch("CatchenWeight",       &CatchenWeight, "data/D");                              /*}}}*/

    t5->Branch("tgt_px", &tgt_px, "tgt_px/D");/*{{{*/
    t5->Branch("tgt_py", &tgt_py, "tgt_py/D");
    t5->Branch("tgt_pz", &tgt_pz, "tgt_pz/D");
    t5->Branch("tgt_theta", &tgt_theta, "tgt_theta/D");
    t5->Branch("tgt_phi", &tgt_phi, "tgt_phi/D");
    t5->Branch("tgt_ene", &tgt_ene, "tgt_ene/D");
    t5->Branch("tgt_ene", &tgt_mom, "tgt_mom/D");/*}}}*/

    t5->Branch("beam_px", &beam_px, "beam_px/D");/*{{{*/
    t5->Branch("beam_py", &beam_py, "beam_py/D");
    t5->Branch("beam_pz", &beam_pz, "beam_pz/D");
    t5->Branch("beam_theta", &beam_theta, "beam_theta/D");
    t5->Branch("beam_phi", &beam_phi, "beam_phi/D");
    t5->Branch("beam_ene", &beam_ene, "beam_ene/D");
    t5->Branch("beam_ene", &beam_mom, "beam_mom/D");

    t5->Branch("beam_corr_px", &beam_corr_px, "beam_corr_px/D");
    t5->Branch("beam_corr_py", &beam_corr_py, "beam_corr_py/D");
    t5->Branch("beam_corr_pz", &beam_corr_pz, "beam_corr_pz/D");
    t5->Branch("beam_corr_theta", &beam_corr_theta, "beam_corr_theta/D");
    t5->Branch("beam_corr_phi", &beam_corr_phi, "beam_corr_phi/D");
    t5->Branch("beam_corr_ene", &beam_corr_ene, "beam_corr_ene/D");
    t5->Branch("beam_corr_ene", &beam_corr_mom, "beam_corr_mom/D");/*}}}*/

    t5->Branch("pim_ene",   &pim_ene   ,"pim_ene/D");/*{{{*/
    t5->Branch("pim_px",    &pim_px    ,"pim_px/D");
    t5->Branch("pim_py",    &pim_py    ,"pim_py/D");
    t5->Branch("pim_pz",    &pim_pz    ,"pim_pz/D");
    t5->Branch("pim_mom",   &pim_mom   ,"pim_mom/D");
    t5->Branch("pim_theta", &pim_theta ,"pim_theta/D");
    t5->Branch("pim_phi",   &pim_phi   ,"pim_phi/D");

    t5->Branch("pim_corr_ene",   &pim_corr_ene   ,"pim_corr_ene/D");
    t5->Branch("pim_corr_px",    &pim_corr_px    ,"pim_corr_px/D");
    t5->Branch("pim_corr_py",    &pim_corr_py    ,"pim_corr_py/D");
    t5->Branch("pim_corr_pz",    &pim_corr_pz    ,"pim_corr_pz/D");
    t5->Branch("pim_corr_mom",   &pim_corr_mom   ,"pim_corr_mom/D");
    t5->Branch("pim_corr_theta", &pim_corr_theta ,"pim_corr_theta/D");
    t5->Branch("pim_corr_phi",   &pim_corr_phi   ,"pim_corr_phi/D");/*}}}*/

    t5->Branch("ele_ene",   &ele_ene   ,"ele_ene/D");/*{{{*/
    t5->Branch("ele_px",    &ele_px    ,"ele_px/D");
    t5->Branch("ele_py",    &ele_py    ,"ele_py/D");
    t5->Branch("ele_pz",    &ele_pz    ,"ele_pz/D");
    t5->Branch("ele_mom",   &ele_mom   ,"ele_mom/D");
    t5->Branch("ele_theta", &ele_theta ,"ele_theta/D");
    t5->Branch("ele_phi",   &ele_phi   ,"ele_phi/D");

    t5->Branch("ele_corr_ene",   &ele_corr_ene   ,"ele_corr_ene/D");
    t5->Branch("ele_corr_px",    &ele_corr_px    ,"ele_corr_px/D");
    t5->Branch("ele_corr_py",    &ele_corr_py    ,"ele_corr_py/D");
    t5->Branch("ele_corr_pz",    &ele_corr_pz    ,"ele_corr_pz/D");
    t5->Branch("ele_corr_mom",   &ele_corr_mom   ,"ele_corr_mom/D");
    t5->Branch("ele_corr_theta", &ele_corr_theta ,"ele_corr_theta/D");
    t5->Branch("ele_corr_phi",   &ele_corr_phi   ,"ele_corr_phi/D");/*}}}*/

    t5->Branch("pro_ene",   &pro_ene   ,"pro_ene/D"); /*{{{*/
    t5->Branch("pro_px",    &pro_px    ,"pro_px/D");
    t5->Branch("pro_py",    &pro_py    ,"pro_py/D");
    t5->Branch("pro_pz",    &pro_pz    ,"pro_pz/D");
    t5->Branch("pro_mom",   &pro_mom   ,"pro_mom/D");
    t5->Branch("pro_theta", &pro_theta ,"pro_theta/D");
    t5->Branch("pro_phi",   &pro_phi   ,"pro_phi/D");

    t5->Branch("pro_corr_ene",   &pro_corr_ene   ,"pro_corr_ene/D");
    t5->Branch("pro_corr_px",    &pro_corr_px    ,"pro_corr_px/D");
    t5->Branch("pro_corr_py",    &pro_corr_py    ,"pro_corr_py/D");
    t5->Branch("pro_corr_pz",    &pro_corr_pz    ,"pro_corr_pz/D");
    t5->Branch("pro_corr_mom",   &pro_corr_mom   ,"pro_corr_mom/D");
    t5->Branch("pro_corr_theta", &pro_corr_theta ,"pro_corr_theta/D");
    t5->Branch("pro_corr_phi",   &pro_corr_phi   ,"pro_corr_phi/D");/*}}}*/

    //Add SoLID acceptance/*{{{*/
    t5->Branch("ele_acc_f",     &ele_acc_f,     "ele_acc_f/D");
    t5->Branch("ele_acc_l",     &ele_acc_l,     "ele_acc_l/D");
    t5->Branch("pim_acc_f",     &pim_acc_f,     "pim_acc_f/D");
    t5->Branch("pim_acc_l",     &pim_acc_l,     "pim_acc_l/D");
    t5->Branch("pro_acc_f",     &pro_acc_f,     "pro_acc_f/D");
    t5->Branch("pro_acc_l",     &pro_acc_l,     "pro_acc_l/D");/*}}}*/

    //Add detector resolutions/*{{{*/
    t5->Branch("ele_mom_res",   &ele_mom_res,   "ele_mom_res/D");
    t5->Branch("ele_theta_res", &ele_theta_res, "ele_theta_res/D");
    t5->Branch("ele_phi_res",   &ele_phi_res,   "ele_phi_res/D");
    t5->Branch("pim_mom_res",   &pim_mom_res,   "pim_mom_res/D");
    t5->Branch("pim_theta_res", &pim_theta_res, "pim_theta_res/D");
    t5->Branch("pim_phi_res",   &pim_phi_res,   "pim_phi_res/D");
    t5->Branch("pro_mom_res",   &pro_mom_res,   "pro_mom_res/D");
    t5->Branch("pro_theta_res", &pro_theta_res, "pro_theta_res/D");
    t5->Branch("pro_phi_res",   &pro_phi_res,   "pro_phi_res/D");/*}}}*/

    //Add other quantities/*{{{*/
    t5->Branch("weight",        &weight,        "weight/D");
    t5->Branch("weight_uu",        &weight_uu,        "weight_uu/D"); //weight for unpolarized XS
    t5->Branch("weight_ut",        &weight_ut,        "weight_ut/D"); //weight for polarized XS
    t5->Branch("weight_3m1",        &weight_3m1,        "weight_3m1/D"); //weight for Sin(3Phi-PhiS) module
    t5->Branch("weight_2m1",        &weight_2m1,        "weight_2m1/D"); //weight for Sin(2Phi-PhiS) module
    t5->Branch("weight_1m1",        &weight_1m1,        "weight_1m1/D"); //weight for Sin(Phi-PhiS) module
    t5->Branch("weight_0p1",        &weight_0p1,        "weight_0p1/D"); //weight for Sin(PhiS) module
    t5->Branch("weight_1p1",        &weight_1p1,        "weight_1p1/D"); //weight for Sin(Phi+PhiS) module
    t5->Branch("weight_2p1",        &weight_2p1,        "weight_2p1/D"); //weight for Sin(2Phi+PhiS) module
    t5->Branch("dilute",        &dilute,        "dilute/D");
    //t5->Branch("PSF",           &PSF,           "PSF/D");
    t5->Branch("MM",            &MM,            "MM/D");
    t5->Branch("MM_res", &MM_res, "MM_res/D");/*}}}*/
    
    t5->Branch("time",            &time,            "time/D");
    //t5->Branch("fileNO",            &fileNO,            "fileNO/I");
    t5->Branch("total_acceptance",            &total_acceptance,            "total_acceptance/D");
    t5->Branch("Q2BIN", &Q2BIN, "Q2BIN/I");
    //t5->Branch("",            &,            "/D");
    /*}}}*/

    /*Define #6 new Tree and new Branch{{{*/
    TFile *f6 = new TFile(Form("./rootfiles/dvmp_%s_t6_newt.root", pol_name.Data()), "recreate");
    TTree *t6 = new TTree("T","a new tree");

    t6->Branch("NRecorded",       &NRecorded,     "NRecorded/I");/*{{{*/
    t6->Branch("NGenerated",       &NGenerated,     "NGenerated/I");

    t6->Branch("Epsilon", &Epsilon, "Epsilon/D" );
    t6->Branch("Qsq", &Qsq ,"Qsq/D");
    t6->Branch("t", &t ,"t/D");
    t6->Branch("W", &W ,"W/D");
    t6->Branch("x", &x ,"x/D");
    t6->Branch("y", &y ,"y/D");
    t6->Branch("z", &z ,"z/D");

    t6->Branch("Qsq_corr", &Qsq_corr ,"Qsq_corr/D");
    t6->Branch("t_corr", &t_corr ,"t_corr/D");
    t6->Branch("W_corr", &W_corr ,"W_corr/D");
    t6->Branch("x_corr", &x_corr ,"x_corr/D");
    t6->Branch("y_corr", &y_corr ,"y_corr/D");
    t6->Branch("z_corr", &z_corr ,"z_corr/D"); 

    t6->Branch("Vertex_X",   &Vertex_X   ,"Vertex_X/D");
    t6->Branch("Vertex_Y",   &Vertex_Y   ,"Vertex_Y/D");
    t6->Branch("Vertex_Z",   &Vertex_Z   ,"Vertex_Z/D");
    t6->Branch("Theta_Pion_Photon",   &Theta_Pion_Photon   ,"Theta_Pion_Photon/D");

    t6->Branch("A_Factor",                                  &A_Factor          ,"data/D"    );
    t6->Branch("Flux_Factor_RF",                            &Flux_Factor_RF       ,"data/D"    );
    t6->Branch("Flux_Factor_Col",                           &Flux_Factor_Col   ,"data/D"    );
    t6->Branch("Jacobian_CM",                               &Jacobian_CM       ,"data/D"    );
    t6->Branch("Jacobian_CM_RF",                            &Jacobian_CM_RF    ,"data/D"    );
    t6->Branch("Jacobian_CM_Col",                           &Jacobian_CM_Col   ,"data/D"    );
    /*}}}*/

    t6->Branch("Phi",       &Phi,       "Phi/D");/*{{{*/
    t6->Branch("PhiS",      &PhiS,      "PhiS/D");
    t6->Branch("Phi_corr",  &Phi_corr,  "Phi_corr/D");
    t6->Branch("PhiS_corr", &PhiS_corr, "PhiS_corr/D");/*}}}*/

    t6->Branch("Sigma_Lab",     &Sigma_Lab, "data/D");                              /*{{{*/
    t6->Branch("Sigma_UU",      &Sigma_UU,   "data/D");                              
    t6->Branch("Sigma_UT",      &Sigma_UT,   "data/D");                              
    t6->Branch("Sig_T",         &Sig_T, "data/D");                                  
    t6->Branch("Sig_L",         &Sig_L, "data/D");                                  
    t6->Branch("Sig_LT",        &Sig_LT, "data/D");                                 
    t6->Branch("Sig_TT",        &Sig_TT, "data/D");                                      /*}}}*/

    //Six Asymmetries and polarized XS/*{{{*/
    t6->Branch("SSAsym",              &SSAsym, "data/D");
    t6->Branch("SineAsym",            &SineAsym, "data/D");
    t6->Branch("Asym_PhiS",           &Asym_PhiS, "data/D");
    t6->Branch("Asym_PhiPlusPhiS",    &Asym_PhiPlusPhiS, "data/D");
    t6->Branch("Asym_2PhiPlusPhiS",   &Asym_2PhiPlusPhiS, "data/D");
    t6->Branch("Asym_PhiMinusPhiS",   &Asym_PhiMinusPhiS, "data/D");
    t6->Branch("Asym_2PhiMinusPhiS",  &Asym_2PhiMinusPhiS, "data/D");
    t6->Branch("Asym_3PhiMinusPhiS",  &Asym_3PhiMinusPhiS, "data/D");

    t6->Branch("Sigma_PhiS",           &Sigma_PhiS, "data/D");
    t6->Branch("Sigma_PhiPlusPhiS",    &Sigma_PhiPlusPhiS, "data/D");
    t6->Branch("Sigma_2PhiPlusPhiS",   &Sigma_2PhiPlusPhiS, "data/D");
    t6->Branch("Sigma_PhiMinusPhiS",   &Sigma_PhiMinusPhiS, "data/D");
    t6->Branch("Sigma_2PhiMinusPhiS",  &Sigma_2PhiMinusPhiS, "data/D");
    t6->Branch("Sigma_3PhiMinusPhiS",  &Sigma_3PhiMinusPhiS, "data/D");/*}}}*/

    t6->Branch("EventWeight",         &EventWeight, "data/D");                              /*{{{*/
    t6->Branch("WilliamsWeight",      &WilliamsWeight, "data/D");                              
    t6->Branch("DedrickWeight",       &DedrickWeight, "data/D");                              
    t6->Branch("CatchenWeight",       &CatchenWeight, "data/D");                              /*}}}*/

    t6->Branch("tgt_px", &tgt_px, "tgt_px/D");/*{{{*/
    t6->Branch("tgt_py", &tgt_py, "tgt_py/D");
    t6->Branch("tgt_pz", &tgt_pz, "tgt_pz/D");
    t6->Branch("tgt_theta", &tgt_theta, "tgt_theta/D");
    t6->Branch("tgt_phi", &tgt_phi, "tgt_phi/D");
    t6->Branch("tgt_ene", &tgt_ene, "tgt_ene/D");
    t6->Branch("tgt_ene", &tgt_mom, "tgt_mom/D");/*}}}*/

    t6->Branch("beam_px", &beam_px, "beam_px/D");/*{{{*/
    t6->Branch("beam_py", &beam_py, "beam_py/D");
    t6->Branch("beam_pz", &beam_pz, "beam_pz/D");
    t6->Branch("beam_theta", &beam_theta, "beam_theta/D");
    t6->Branch("beam_phi", &beam_phi, "beam_phi/D");
    t6->Branch("beam_ene", &beam_ene, "beam_ene/D");
    t6->Branch("beam_ene", &beam_mom, "beam_mom/D");

    t6->Branch("beam_corr_px", &beam_corr_px, "beam_corr_px/D");
    t6->Branch("beam_corr_py", &beam_corr_py, "beam_corr_py/D");
    t6->Branch("beam_corr_pz", &beam_corr_pz, "beam_corr_pz/D");
    t6->Branch("beam_corr_theta", &beam_corr_theta, "beam_corr_theta/D");
    t6->Branch("beam_corr_phi", &beam_corr_phi, "beam_corr_phi/D");
    t6->Branch("beam_corr_ene", &beam_corr_ene, "beam_corr_ene/D");
    t6->Branch("beam_corr_ene", &beam_corr_mom, "beam_corr_mom/D");/*}}}*/

    t6->Branch("pim_ene",   &pim_ene   ,"pim_ene/D");/*{{{*/
    t6->Branch("pim_px",    &pim_px    ,"pim_px/D");
    t6->Branch("pim_py",    &pim_py    ,"pim_py/D");
    t6->Branch("pim_pz",    &pim_pz    ,"pim_pz/D");
    t6->Branch("pim_mom",   &pim_mom   ,"pim_mom/D");
    t6->Branch("pim_theta", &pim_theta ,"pim_theta/D");
    t6->Branch("pim_phi",   &pim_phi   ,"pim_phi/D");

    t6->Branch("pim_corr_ene",   &pim_corr_ene   ,"pim_corr_ene/D");
    t6->Branch("pim_corr_px",    &pim_corr_px    ,"pim_corr_px/D");
    t6->Branch("pim_corr_py",    &pim_corr_py    ,"pim_corr_py/D");
    t6->Branch("pim_corr_pz",    &pim_corr_pz    ,"pim_corr_pz/D");
    t6->Branch("pim_corr_mom",   &pim_corr_mom   ,"pim_corr_mom/D");
    t6->Branch("pim_corr_theta", &pim_corr_theta ,"pim_corr_theta/D");
    t6->Branch("pim_corr_phi",   &pim_corr_phi   ,"pim_corr_phi/D");/*}}}*/

    t6->Branch("ele_ene",   &ele_ene   ,"ele_ene/D");/*{{{*/
    t6->Branch("ele_px",    &ele_px    ,"ele_px/D");
    t6->Branch("ele_py",    &ele_py    ,"ele_py/D");
    t6->Branch("ele_pz",    &ele_pz    ,"ele_pz/D");
    t6->Branch("ele_mom",   &ele_mom   ,"ele_mom/D");
    t6->Branch("ele_theta", &ele_theta ,"ele_theta/D");
    t6->Branch("ele_phi",   &ele_phi   ,"ele_phi/D");

    t6->Branch("ele_corr_ene",   &ele_corr_ene   ,"ele_corr_ene/D");
    t6->Branch("ele_corr_px",    &ele_corr_px    ,"ele_corr_px/D");
    t6->Branch("ele_corr_py",    &ele_corr_py    ,"ele_corr_py/D");
    t6->Branch("ele_corr_pz",    &ele_corr_pz    ,"ele_corr_pz/D");
    t6->Branch("ele_corr_mom",   &ele_corr_mom   ,"ele_corr_mom/D");
    t6->Branch("ele_corr_theta", &ele_corr_theta ,"ele_corr_theta/D");
    t6->Branch("ele_corr_phi",   &ele_corr_phi   ,"ele_corr_phi/D");/*}}}*/

    t6->Branch("pro_ene",   &pro_ene   ,"pro_ene/D"); /*{{{*/
    t6->Branch("pro_px",    &pro_px    ,"pro_px/D");
    t6->Branch("pro_py",    &pro_py    ,"pro_py/D");
    t6->Branch("pro_pz",    &pro_pz    ,"pro_pz/D");
    t6->Branch("pro_mom",   &pro_mom   ,"pro_mom/D");
    t6->Branch("pro_theta", &pro_theta ,"pro_theta/D");
    t6->Branch("pro_phi",   &pro_phi   ,"pro_phi/D");

    t6->Branch("pro_corr_ene",   &pro_corr_ene   ,"pro_corr_ene/D");
    t6->Branch("pro_corr_px",    &pro_corr_px    ,"pro_corr_px/D");
    t6->Branch("pro_corr_py",    &pro_corr_py    ,"pro_corr_py/D");
    t6->Branch("pro_corr_pz",    &pro_corr_pz    ,"pro_corr_pz/D");
    t6->Branch("pro_corr_mom",   &pro_corr_mom   ,"pro_corr_mom/D");
    t6->Branch("pro_corr_theta", &pro_corr_theta ,"pro_corr_theta/D");
    t6->Branch("pro_corr_phi",   &pro_corr_phi   ,"pro_corr_phi/D");/*}}}*/

    //Add SoLID acceptance/*{{{*/
    t6->Branch("ele_acc_f",     &ele_acc_f,     "ele_acc_f/D");
    t6->Branch("ele_acc_l",     &ele_acc_l,     "ele_acc_l/D");
    t6->Branch("pim_acc_f",     &pim_acc_f,     "pim_acc_f/D");
    t6->Branch("pim_acc_l",     &pim_acc_l,     "pim_acc_l/D");
    t6->Branch("pro_acc_f",     &pro_acc_f,     "pro_acc_f/D");
    t6->Branch("pro_acc_l",     &pro_acc_l,     "pro_acc_l/D");/*}}}*/

    //Add detector resolutions/*{{{*/
    t6->Branch("ele_mom_res",   &ele_mom_res,   "ele_mom_res/D");
    t6->Branch("ele_theta_res", &ele_theta_res, "ele_theta_res/D");
    t6->Branch("ele_phi_res",   &ele_phi_res,   "ele_phi_res/D");
    t6->Branch("pim_mom_res",   &pim_mom_res,   "pim_mom_res/D");
    t6->Branch("pim_theta_res", &pim_theta_res, "pim_theta_res/D");
    t6->Branch("pim_phi_res",   &pim_phi_res,   "pim_phi_res/D");
    t6->Branch("pro_mom_res",   &pro_mom_res,   "pro_mom_res/D");
    t6->Branch("pro_theta_res", &pro_theta_res, "pro_theta_res/D");
    t6->Branch("pro_phi_res",   &pro_phi_res,   "pro_phi_res/D");/*}}}*/

    //Add other quantities/*{{{*/
    t6->Branch("weight",        &weight,        "weight/D");
    t6->Branch("weight_uu",        &weight_uu,        "weight_uu/D"); //weight for unpolarized XS
    t6->Branch("weight_ut",        &weight_ut,        "weight_ut/D"); //weight for polarized XS
    t6->Branch("weight_3m1",        &weight_3m1,        "weight_3m1/D"); //weight for Sin(3Phi-PhiS) module
    t6->Branch("weight_2m1",        &weight_2m1,        "weight_2m1/D"); //weight for Sin(2Phi-PhiS) module
    t6->Branch("weight_1m1",        &weight_1m1,        "weight_1m1/D"); //weight for Sin(Phi-PhiS) module
    t6->Branch("weight_0p1",        &weight_0p1,        "weight_0p1/D"); //weight for Sin(PhiS) module
    t6->Branch("weight_1p1",        &weight_1p1,        "weight_1p1/D"); //weight for Sin(Phi+PhiS) module
    t6->Branch("weight_2p1",        &weight_2p1,        "weight_2p1/D"); //weight for Sin(2Phi+PhiS) module
    t6->Branch("dilute",        &dilute,        "dilute/D");
    //t6->Branch("PSF",           &PSF,           "PSF/D");
    t6->Branch("MM",            &MM,            "MM/D");
    t6->Branch("MM_res", &MM_res, "MM_res/D");/*}}}*/
    
    t6->Branch("time",            &time,            "time/D");
    //t6->Branch("fileNO",            &fileNO,            "fileNO/I");
    t6->Branch("total_acceptance",            &total_acceptance,            "total_acceptance/D");
    t6->Branch("Q2BIN", &Q2BIN, "Q2BIN/I");
    //t6->Branch("",            &,            "/D");
    /*}}}*/
    
    /*Define #7 new Tree and new Branch{{{*/
    TFile *f7 = new TFile(Form("./rootfiles/dvmp_%s_t7_newt.root", pol_name.Data()), "recreate");
    TTree *t7 = new TTree("T","a new tree");

    t7->Branch("NRecorded",       &NRecorded,     "NRecorded/I");/*{{{*/
    t7->Branch("NGenerated",       &NGenerated,     "NGenerated/I");

    t7->Branch("Epsilon", &Epsilon, "Epsilon/D" );
    t7->Branch("Qsq", &Qsq ,"Qsq/D");
    t7->Branch("t", &t ,"t/D");
    t7->Branch("W", &W ,"W/D");
    t7->Branch("x", &x ,"x/D");
    t7->Branch("y", &y ,"y/D");
    t7->Branch("z", &z ,"z/D");

    t7->Branch("Qsq_corr", &Qsq_corr ,"Qsq_corr/D");
    t7->Branch("t_corr", &t_corr ,"t_corr/D");
    t7->Branch("W_corr", &W_corr ,"W_corr/D");
    t7->Branch("x_corr", &x_corr ,"x_corr/D");
    t7->Branch("y_corr", &y_corr ,"y_corr/D");
    t7->Branch("z_corr", &z_corr ,"z_corr/D"); 

    t7->Branch("Vertex_X",   &Vertex_X   ,"Vertex_X/D");
    t7->Branch("Vertex_Y",   &Vertex_Y   ,"Vertex_Y/D");
    t7->Branch("Vertex_Z",   &Vertex_Z   ,"Vertex_Z/D");
    t7->Branch("Theta_Pion_Photon",   &Theta_Pion_Photon   ,"Theta_Pion_Photon/D");

    t7->Branch("A_Factor",                                  &A_Factor          ,"data/D"    );
    t7->Branch("Flux_Factor_RF",                            &Flux_Factor_RF       ,"data/D"    );
    t7->Branch("Flux_Factor_Col",                           &Flux_Factor_Col   ,"data/D"    );
    t7->Branch("Jacobian_CM",                               &Jacobian_CM       ,"data/D"    );
    t7->Branch("Jacobian_CM_RF",                            &Jacobian_CM_RF    ,"data/D"    );
    t7->Branch("Jacobian_CM_Col",                           &Jacobian_CM_Col   ,"data/D"    );
    /*}}}*/

    t7->Branch("Phi",       &Phi,       "Phi/D");/*{{{*/
    t7->Branch("PhiS",      &PhiS,      "PhiS/D");
    t7->Branch("Phi_corr",  &Phi_corr,  "Phi_corr/D");
    t7->Branch("PhiS_corr", &PhiS_corr, "PhiS_corr/D");/*}}}*/

    t7->Branch("Sigma_Lab",     &Sigma_Lab, "data/D");                              /*{{{*/
    t7->Branch("Sigma_UU",      &Sigma_UU,   "data/D");                              
    t7->Branch("Sigma_UT",      &Sigma_UT,   "data/D");                              
    t7->Branch("Sig_T",         &Sig_T, "data/D");                                  
    t7->Branch("Sig_L",         &Sig_L, "data/D");                                  
    t7->Branch("Sig_LT",        &Sig_LT, "data/D");                                 
    t7->Branch("Sig_TT",        &Sig_TT, "data/D");                                      /*}}}*/

    //Six Asymmetries and polarized XS/*{{{*/
    t7->Branch("SSAsym",              &SSAsym, "data/D");
    t7->Branch("SineAsym",            &SineAsym, "data/D");
    t7->Branch("Asym_PhiS",           &Asym_PhiS, "data/D");
    t7->Branch("Asym_PhiPlusPhiS",    &Asym_PhiPlusPhiS, "data/D");
    t7->Branch("Asym_2PhiPlusPhiS",   &Asym_2PhiPlusPhiS, "data/D");
    t7->Branch("Asym_PhiMinusPhiS",   &Asym_PhiMinusPhiS, "data/D");
    t7->Branch("Asym_2PhiMinusPhiS",  &Asym_2PhiMinusPhiS, "data/D");
    t7->Branch("Asym_3PhiMinusPhiS",  &Asym_3PhiMinusPhiS, "data/D");

    t7->Branch("Sigma_PhiS",           &Sigma_PhiS, "data/D");
    t7->Branch("Sigma_PhiPlusPhiS",    &Sigma_PhiPlusPhiS, "data/D");
    t7->Branch("Sigma_2PhiPlusPhiS",   &Sigma_2PhiPlusPhiS, "data/D");
    t7->Branch("Sigma_PhiMinusPhiS",   &Sigma_PhiMinusPhiS, "data/D");
    t7->Branch("Sigma_2PhiMinusPhiS",  &Sigma_2PhiMinusPhiS, "data/D");
    t7->Branch("Sigma_3PhiMinusPhiS",  &Sigma_3PhiMinusPhiS, "data/D");/*}}}*/

    t7->Branch("EventWeight",         &EventWeight, "data/D");                              /*{{{*/
    t7->Branch("WilliamsWeight",      &WilliamsWeight, "data/D");                              
    t7->Branch("DedrickWeight",       &DedrickWeight, "data/D");                              
    t7->Branch("CatchenWeight",       &CatchenWeight, "data/D");                              /*}}}*/

    t7->Branch("tgt_px", &tgt_px, "tgt_px/D");/*{{{*/
    t7->Branch("tgt_py", &tgt_py, "tgt_py/D");
    t7->Branch("tgt_pz", &tgt_pz, "tgt_pz/D");
    t7->Branch("tgt_theta", &tgt_theta, "tgt_theta/D");
    t7->Branch("tgt_phi", &tgt_phi, "tgt_phi/D");
    t7->Branch("tgt_ene", &tgt_ene, "tgt_ene/D");
    t7->Branch("tgt_ene", &tgt_mom, "tgt_mom/D");/*}}}*/

    t7->Branch("beam_px", &beam_px, "beam_px/D");/*{{{*/
    t7->Branch("beam_py", &beam_py, "beam_py/D");
    t7->Branch("beam_pz", &beam_pz, "beam_pz/D");
    t7->Branch("beam_theta", &beam_theta, "beam_theta/D");
    t7->Branch("beam_phi", &beam_phi, "beam_phi/D");
    t7->Branch("beam_ene", &beam_ene, "beam_ene/D");
    t7->Branch("beam_ene", &beam_mom, "beam_mom/D");

    t7->Branch("beam_corr_px", &beam_corr_px, "beam_corr_px/D");
    t7->Branch("beam_corr_py", &beam_corr_py, "beam_corr_py/D");
    t7->Branch("beam_corr_pz", &beam_corr_pz, "beam_corr_pz/D");
    t7->Branch("beam_corr_theta", &beam_corr_theta, "beam_corr_theta/D");
    t7->Branch("beam_corr_phi", &beam_corr_phi, "beam_corr_phi/D");
    t7->Branch("beam_corr_ene", &beam_corr_ene, "beam_corr_ene/D");
    t7->Branch("beam_corr_ene", &beam_corr_mom, "beam_corr_mom/D");/*}}}*/

    t7->Branch("pim_ene",   &pim_ene   ,"pim_ene/D");/*{{{*/
    t7->Branch("pim_px",    &pim_px    ,"pim_px/D");
    t7->Branch("pim_py",    &pim_py    ,"pim_py/D");
    t7->Branch("pim_pz",    &pim_pz    ,"pim_pz/D");
    t7->Branch("pim_mom",   &pim_mom   ,"pim_mom/D");
    t7->Branch("pim_theta", &pim_theta ,"pim_theta/D");
    t7->Branch("pim_phi",   &pim_phi   ,"pim_phi/D");

    t7->Branch("pim_corr_ene",   &pim_corr_ene   ,"pim_corr_ene/D");
    t7->Branch("pim_corr_px",    &pim_corr_px    ,"pim_corr_px/D");
    t7->Branch("pim_corr_py",    &pim_corr_py    ,"pim_corr_py/D");
    t7->Branch("pim_corr_pz",    &pim_corr_pz    ,"pim_corr_pz/D");
    t7->Branch("pim_corr_mom",   &pim_corr_mom   ,"pim_corr_mom/D");
    t7->Branch("pim_corr_theta", &pim_corr_theta ,"pim_corr_theta/D");
    t7->Branch("pim_corr_phi",   &pim_corr_phi   ,"pim_corr_phi/D");/*}}}*/

    t7->Branch("ele_ene",   &ele_ene   ,"ele_ene/D");/*{{{*/
    t7->Branch("ele_px",    &ele_px    ,"ele_px/D");
    t7->Branch("ele_py",    &ele_py    ,"ele_py/D");
    t7->Branch("ele_pz",    &ele_pz    ,"ele_pz/D");
    t7->Branch("ele_mom",   &ele_mom   ,"ele_mom/D");
    t7->Branch("ele_theta", &ele_theta ,"ele_theta/D");
    t7->Branch("ele_phi",   &ele_phi   ,"ele_phi/D");

    t7->Branch("ele_corr_ene",   &ele_corr_ene   ,"ele_corr_ene/D");
    t7->Branch("ele_corr_px",    &ele_corr_px    ,"ele_corr_px/D");
    t7->Branch("ele_corr_py",    &ele_corr_py    ,"ele_corr_py/D");
    t7->Branch("ele_corr_pz",    &ele_corr_pz    ,"ele_corr_pz/D");
    t7->Branch("ele_corr_mom",   &ele_corr_mom   ,"ele_corr_mom/D");
    t7->Branch("ele_corr_theta", &ele_corr_theta ,"ele_corr_theta/D");
    t7->Branch("ele_corr_phi",   &ele_corr_phi   ,"ele_corr_phi/D");/*}}}*/

    t7->Branch("pro_ene",   &pro_ene   ,"pro_ene/D"); /*{{{*/
    t7->Branch("pro_px",    &pro_px    ,"pro_px/D");
    t7->Branch("pro_py",    &pro_py    ,"pro_py/D");
    t7->Branch("pro_pz",    &pro_pz    ,"pro_pz/D");
    t7->Branch("pro_mom",   &pro_mom   ,"pro_mom/D");
    t7->Branch("pro_theta", &pro_theta ,"pro_theta/D");
    t7->Branch("pro_phi",   &pro_phi   ,"pro_phi/D");

    t7->Branch("pro_corr_ene",   &pro_corr_ene   ,"pro_corr_ene/D");
    t7->Branch("pro_corr_px",    &pro_corr_px    ,"pro_corr_px/D");
    t7->Branch("pro_corr_py",    &pro_corr_py    ,"pro_corr_py/D");
    t7->Branch("pro_corr_pz",    &pro_corr_pz    ,"pro_corr_pz/D");
    t7->Branch("pro_corr_mom",   &pro_corr_mom   ,"pro_corr_mom/D");
    t7->Branch("pro_corr_theta", &pro_corr_theta ,"pro_corr_theta/D");
    t7->Branch("pro_corr_phi",   &pro_corr_phi   ,"pro_corr_phi/D");/*}}}*/

    //Add SoLID acceptance/*{{{*/
    t7->Branch("ele_acc_f",     &ele_acc_f,     "ele_acc_f/D");
    t7->Branch("ele_acc_l",     &ele_acc_l,     "ele_acc_l/D");
    t7->Branch("pim_acc_f",     &pim_acc_f,     "pim_acc_f/D");
    t7->Branch("pim_acc_l",     &pim_acc_l,     "pim_acc_l/D");
    t7->Branch("pro_acc_f",     &pro_acc_f,     "pro_acc_f/D");
    t7->Branch("pro_acc_l",     &pro_acc_l,     "pro_acc_l/D");/*}}}*/

    //Add detector resolutions/*{{{*/
    t7->Branch("ele_mom_res",   &ele_mom_res,   "ele_mom_res/D");
    t7->Branch("ele_theta_res", &ele_theta_res, "ele_theta_res/D");
    t7->Branch("ele_phi_res",   &ele_phi_res,   "ele_phi_res/D");
    t7->Branch("pim_mom_res",   &pim_mom_res,   "pim_mom_res/D");
    t7->Branch("pim_theta_res", &pim_theta_res, "pim_theta_res/D");
    t7->Branch("pim_phi_res",   &pim_phi_res,   "pim_phi_res/D");
    t7->Branch("pro_mom_res",   &pro_mom_res,   "pro_mom_res/D");
    t7->Branch("pro_theta_res", &pro_theta_res, "pro_theta_res/D");
    t7->Branch("pro_phi_res",   &pro_phi_res,   "pro_phi_res/D");/*}}}*/

    //Add other quantities/*{{{*/
    t7->Branch("weight",        &weight,        "weight/D");
    t7->Branch("weight_uu",        &weight_uu,        "weight_uu/D"); //weight for unpolarized XS
    t7->Branch("weight_ut",        &weight_ut,        "weight_ut/D"); //weight for polarized XS
    t7->Branch("weight_3m1",        &weight_3m1,        "weight_3m1/D"); //weight for Sin(3Phi-PhiS) module
    t7->Branch("weight_2m1",        &weight_2m1,        "weight_2m1/D"); //weight for Sin(2Phi-PhiS) module
    t7->Branch("weight_1m1",        &weight_1m1,        "weight_1m1/D"); //weight for Sin(Phi-PhiS) module
    t7->Branch("weight_0p1",        &weight_0p1,        "weight_0p1/D"); //weight for Sin(PhiS) module
    t7->Branch("weight_1p1",        &weight_1p1,        "weight_1p1/D"); //weight for Sin(Phi+PhiS) module
    t7->Branch("weight_2p1",        &weight_2p1,        "weight_2p1/D"); //weight for Sin(2Phi+PhiS) module
    t7->Branch("dilute",        &dilute,        "dilute/D");
    //t7->Branch("PSF",           &PSF,           "PSF/D");
    t7->Branch("MM",            &MM,            "MM/D");
    t7->Branch("MM_res", &MM_res, "MM_res/D");/*}}}*/
    
    t7->Branch("time",            &time,            "time/D");
    //t7->Branch("fileNO",            &fileNO,            "fileNO/I");
    t7->Branch("total_acceptance",            &total_acceptance,            "total_acceptance/D");
    t7->Branch("Q2BIN", &Q2BIN, "Q2BIN/I");
    //t7->Branch("",            &,            "/D");
    /*}}}*/
    
    /*Define #8 new Tree and new Branch{{{*/
    TFile *f8 = new TFile(Form("./rootfiles/dvmp_%s_t8_newt.root", pol_name.Data()), "recreate");
    TTree *t8 = new TTree("T","a new tree");

    t8->Branch("NRecorded",       &NRecorded,     "NRecorded/I");/*{{{*/
    t8->Branch("NGenerated",       &NGenerated,     "NGenerated/I");

    t8->Branch("Epsilon", &Epsilon, "Epsilon/D" );
    t8->Branch("Qsq", &Qsq ,"Qsq/D");
    t8->Branch("t", &t ,"t/D");
    t8->Branch("W", &W ,"W/D");
    t8->Branch("x", &x ,"x/D");
    t8->Branch("y", &y ,"y/D");
    t8->Branch("z", &z ,"z/D");

    t8->Branch("Qsq_corr", &Qsq_corr ,"Qsq_corr/D");
    t8->Branch("t_corr", &t_corr ,"t_corr/D");
    t8->Branch("W_corr", &W_corr ,"W_corr/D");
    t8->Branch("x_corr", &x_corr ,"x_corr/D");
    t8->Branch("y_corr", &y_corr ,"y_corr/D");
    t8->Branch("z_corr", &z_corr ,"z_corr/D"); 

    t8->Branch("Vertex_X",   &Vertex_X   ,"Vertex_X/D");
    t8->Branch("Vertex_Y",   &Vertex_Y   ,"Vertex_Y/D");
    t8->Branch("Vertex_Z",   &Vertex_Z   ,"Vertex_Z/D");
    t8->Branch("Theta_Pion_Photon",   &Theta_Pion_Photon   ,"Theta_Pion_Photon/D");

    t8->Branch("A_Factor",                                  &A_Factor          ,"data/D"    );
    t8->Branch("Flux_Factor_RF",                            &Flux_Factor_RF       ,"data/D"    );
    t8->Branch("Flux_Factor_Col",                           &Flux_Factor_Col   ,"data/D"    );
    t8->Branch("Jacobian_CM",                               &Jacobian_CM       ,"data/D"    );
    t8->Branch("Jacobian_CM_RF",                            &Jacobian_CM_RF    ,"data/D"    );
    t8->Branch("Jacobian_CM_Col",                           &Jacobian_CM_Col   ,"data/D"    );
    /*}}}*/

    t8->Branch("Phi",       &Phi,       "Phi/D");/*{{{*/
    t8->Branch("PhiS",      &PhiS,      "PhiS/D");
    t8->Branch("Phi_corr",  &Phi_corr,  "Phi_corr/D");
    t8->Branch("PhiS_corr", &PhiS_corr, "PhiS_corr/D");/*}}}*/

    t8->Branch("Sigma_Lab",     &Sigma_Lab, "data/D");                              /*{{{*/
    t8->Branch("Sigma_UU",      &Sigma_UU,   "data/D");                              
    t8->Branch("Sigma_UT",      &Sigma_UT,   "data/D");                              
    t8->Branch("Sig_T",         &Sig_T, "data/D");                                  
    t8->Branch("Sig_L",         &Sig_L, "data/D");                                  
    t8->Branch("Sig_LT",        &Sig_LT, "data/D");                                 
    t8->Branch("Sig_TT",        &Sig_TT, "data/D");                                      /*}}}*/

    //Six Asymmetries and polarized XS/*{{{*/
    t8->Branch("SSAsym",              &SSAsym, "data/D");
    t8->Branch("SineAsym",            &SineAsym, "data/D");
    t8->Branch("Asym_PhiS",           &Asym_PhiS, "data/D");
    t8->Branch("Asym_PhiPlusPhiS",    &Asym_PhiPlusPhiS, "data/D");
    t8->Branch("Asym_2PhiPlusPhiS",   &Asym_2PhiPlusPhiS, "data/D");
    t8->Branch("Asym_PhiMinusPhiS",   &Asym_PhiMinusPhiS, "data/D");
    t8->Branch("Asym_2PhiMinusPhiS",  &Asym_2PhiMinusPhiS, "data/D");
    t8->Branch("Asym_3PhiMinusPhiS",  &Asym_3PhiMinusPhiS, "data/D");

    t8->Branch("Sigma_PhiS",           &Sigma_PhiS, "data/D");
    t8->Branch("Sigma_PhiPlusPhiS",    &Sigma_PhiPlusPhiS, "data/D");
    t8->Branch("Sigma_2PhiPlusPhiS",   &Sigma_2PhiPlusPhiS, "data/D");
    t8->Branch("Sigma_PhiMinusPhiS",   &Sigma_PhiMinusPhiS, "data/D");
    t8->Branch("Sigma_2PhiMinusPhiS",  &Sigma_2PhiMinusPhiS, "data/D");
    t8->Branch("Sigma_3PhiMinusPhiS",  &Sigma_3PhiMinusPhiS, "data/D");/*}}}*/

    t8->Branch("EventWeight",         &EventWeight, "data/D");                              /*{{{*/
    t8->Branch("WilliamsWeight",      &WilliamsWeight, "data/D");                              
    t8->Branch("DedrickWeight",       &DedrickWeight, "data/D");                              
    t8->Branch("CatchenWeight",       &CatchenWeight, "data/D");                              /*}}}*/

    t8->Branch("tgt_px", &tgt_px, "tgt_px/D");/*{{{*/
    t8->Branch("tgt_py", &tgt_py, "tgt_py/D");
    t8->Branch("tgt_pz", &tgt_pz, "tgt_pz/D");
    t8->Branch("tgt_theta", &tgt_theta, "tgt_theta/D");
    t8->Branch("tgt_phi", &tgt_phi, "tgt_phi/D");
    t8->Branch("tgt_ene", &tgt_ene, "tgt_ene/D");
    t8->Branch("tgt_ene", &tgt_mom, "tgt_mom/D");/*}}}*/

    t8->Branch("beam_px", &beam_px, "beam_px/D");/*{{{*/
    t8->Branch("beam_py", &beam_py, "beam_py/D");
    t8->Branch("beam_pz", &beam_pz, "beam_pz/D");
    t8->Branch("beam_theta", &beam_theta, "beam_theta/D");
    t8->Branch("beam_phi", &beam_phi, "beam_phi/D");
    t8->Branch("beam_ene", &beam_ene, "beam_ene/D");
    t8->Branch("beam_ene", &beam_mom, "beam_mom/D");

    t8->Branch("beam_corr_px", &beam_corr_px, "beam_corr_px/D");
    t8->Branch("beam_corr_py", &beam_corr_py, "beam_corr_py/D");
    t8->Branch("beam_corr_pz", &beam_corr_pz, "beam_corr_pz/D");
    t8->Branch("beam_corr_theta", &beam_corr_theta, "beam_corr_theta/D");
    t8->Branch("beam_corr_phi", &beam_corr_phi, "beam_corr_phi/D");
    t8->Branch("beam_corr_ene", &beam_corr_ene, "beam_corr_ene/D");
    t8->Branch("beam_corr_ene", &beam_corr_mom, "beam_corr_mom/D");/*}}}*/

    t8->Branch("pim_ene",   &pim_ene   ,"pim_ene/D");/*{{{*/
    t8->Branch("pim_px",    &pim_px    ,"pim_px/D");
    t8->Branch("pim_py",    &pim_py    ,"pim_py/D");
    t8->Branch("pim_pz",    &pim_pz    ,"pim_pz/D");
    t8->Branch("pim_mom",   &pim_mom   ,"pim_mom/D");
    t8->Branch("pim_theta", &pim_theta ,"pim_theta/D");
    t8->Branch("pim_phi",   &pim_phi   ,"pim_phi/D");

    t8->Branch("pim_corr_ene",   &pim_corr_ene   ,"pim_corr_ene/D");
    t8->Branch("pim_corr_px",    &pim_corr_px    ,"pim_corr_px/D");
    t8->Branch("pim_corr_py",    &pim_corr_py    ,"pim_corr_py/D");
    t8->Branch("pim_corr_pz",    &pim_corr_pz    ,"pim_corr_pz/D");
    t8->Branch("pim_corr_mom",   &pim_corr_mom   ,"pim_corr_mom/D");
    t8->Branch("pim_corr_theta", &pim_corr_theta ,"pim_corr_theta/D");
    t8->Branch("pim_corr_phi",   &pim_corr_phi   ,"pim_corr_phi/D");/*}}}*/

    t8->Branch("ele_ene",   &ele_ene   ,"ele_ene/D");/*{{{*/
    t8->Branch("ele_px",    &ele_px    ,"ele_px/D");
    t8->Branch("ele_py",    &ele_py    ,"ele_py/D");
    t8->Branch("ele_pz",    &ele_pz    ,"ele_pz/D");
    t8->Branch("ele_mom",   &ele_mom   ,"ele_mom/D");
    t8->Branch("ele_theta", &ele_theta ,"ele_theta/D");
    t8->Branch("ele_phi",   &ele_phi   ,"ele_phi/D");

    t8->Branch("ele_corr_ene",   &ele_corr_ene   ,"ele_corr_ene/D");
    t8->Branch("ele_corr_px",    &ele_corr_px    ,"ele_corr_px/D");
    t8->Branch("ele_corr_py",    &ele_corr_py    ,"ele_corr_py/D");
    t8->Branch("ele_corr_pz",    &ele_corr_pz    ,"ele_corr_pz/D");
    t8->Branch("ele_corr_mom",   &ele_corr_mom   ,"ele_corr_mom/D");
    t8->Branch("ele_corr_theta", &ele_corr_theta ,"ele_corr_theta/D");
    t8->Branch("ele_corr_phi",   &ele_corr_phi   ,"ele_corr_phi/D");/*}}}*/

    t8->Branch("pro_ene",   &pro_ene   ,"pro_ene/D"); /*{{{*/
    t8->Branch("pro_px",    &pro_px    ,"pro_px/D");
    t8->Branch("pro_py",    &pro_py    ,"pro_py/D");
    t8->Branch("pro_pz",    &pro_pz    ,"pro_pz/D");
    t8->Branch("pro_mom",   &pro_mom   ,"pro_mom/D");
    t8->Branch("pro_theta", &pro_theta ,"pro_theta/D");
    t8->Branch("pro_phi",   &pro_phi   ,"pro_phi/D");

    t8->Branch("pro_corr_ene",   &pro_corr_ene   ,"pro_corr_ene/D");
    t8->Branch("pro_corr_px",    &pro_corr_px    ,"pro_corr_px/D");
    t8->Branch("pro_corr_py",    &pro_corr_py    ,"pro_corr_py/D");
    t8->Branch("pro_corr_pz",    &pro_corr_pz    ,"pro_corr_pz/D");
    t8->Branch("pro_corr_mom",   &pro_corr_mom   ,"pro_corr_mom/D");
    t8->Branch("pro_corr_theta", &pro_corr_theta ,"pro_corr_theta/D");
    t8->Branch("pro_corr_phi",   &pro_corr_phi   ,"pro_corr_phi/D");/*}}}*/

    //Add SoLID acceptance/*{{{*/
    t8->Branch("ele_acc_f",     &ele_acc_f,     "ele_acc_f/D");
    t8->Branch("ele_acc_l",     &ele_acc_l,     "ele_acc_l/D");
    t8->Branch("pim_acc_f",     &pim_acc_f,     "pim_acc_f/D");
    t8->Branch("pim_acc_l",     &pim_acc_l,     "pim_acc_l/D");
    t8->Branch("pro_acc_f",     &pro_acc_f,     "pro_acc_f/D");
    t8->Branch("pro_acc_l",     &pro_acc_l,     "pro_acc_l/D");/*}}}*/

    //Add detector resolutions/*{{{*/
    t8->Branch("ele_mom_res",   &ele_mom_res,   "ele_mom_res/D");
    t8->Branch("ele_theta_res", &ele_theta_res, "ele_theta_res/D");
    t8->Branch("ele_phi_res",   &ele_phi_res,   "ele_phi_res/D");
    t8->Branch("pim_mom_res",   &pim_mom_res,   "pim_mom_res/D");
    t8->Branch("pim_theta_res", &pim_theta_res, "pim_theta_res/D");
    t8->Branch("pim_phi_res",   &pim_phi_res,   "pim_phi_res/D");
    t8->Branch("pro_mom_res",   &pro_mom_res,   "pro_mom_res/D");
    t8->Branch("pro_theta_res", &pro_theta_res, "pro_theta_res/D");
    t8->Branch("pro_phi_res",   &pro_phi_res,   "pro_phi_res/D");/*}}}*/

    //Add other quantities/*{{{*/
    t8->Branch("weight",        &weight,        "weight/D");
    t8->Branch("weight_uu",        &weight_uu,        "weight_uu/D"); //weight for unpolarized XS
    t8->Branch("weight_ut",        &weight_ut,        "weight_ut/D"); //weight for polarized XS
    t8->Branch("weight_3m1",        &weight_3m1,        "weight_3m1/D"); //weight for Sin(3Phi-PhiS) module
    t8->Branch("weight_2m1",        &weight_2m1,        "weight_2m1/D"); //weight for Sin(2Phi-PhiS) module
    t8->Branch("weight_1m1",        &weight_1m1,        "weight_1m1/D"); //weight for Sin(Phi-PhiS) module
    t8->Branch("weight_0p1",        &weight_0p1,        "weight_0p1/D"); //weight for Sin(PhiS) module
    t8->Branch("weight_1p1",        &weight_1p1,        "weight_1p1/D"); //weight for Sin(Phi+PhiS) module
    t8->Branch("weight_2p1",        &weight_2p1,        "weight_2p1/D"); //weight for Sin(2Phi+PhiS) module
    t8->Branch("dilute",        &dilute,        "dilute/D");
    //t8->Branch("PSF",           &PSF,           "PSF/D");
    t8->Branch("MM",            &MM,            "MM/D");
    t8->Branch("MM_res", &MM_res, "MM_res/D");/*}}}*/
    
    t8->Branch("time",            &time,            "time/D");
    //t8->Branch("fileNO",            &fileNO,            "fileNO/I");
    t8->Branch("total_acceptance",            &total_acceptance,            "total_acceptance/D");
    t8->Branch("Q2BIN", &Q2BIN, "Q2BIN/I");
    //t8->Branch("",            &,            "/D");
    /*}}}*/
    
    //const Int_t tbin = 7;
    //const Double_t t_cut[8] = {0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.75, 1.2};
    //const Int_t tbin = 8;
    const Double_t t_cut[9] = {0.0, 0.30, 0.40, 0.50, 0.60, 0.75, 0.95, 1.2, 2.0};
 
    const double_t par[4] = {3.293, 6.989, -3.174, 0.4937}; //for Q2 cut

    for(Long64_t i=0;i<N_entries;i++){
        t0->GetEntry(i);
        if(Epsilon>0.55&&Epsilon<0.75&&W>2&&Qsq>4){

            double Q2Cut = par[0]+par[1]*t + par[2]*t*t + par[3]*t*t*t;
            if(Qsq<=Q2Cut) Q2BIN = 1;
            if(Qsq>Q2Cut) Q2BIN = 2;

            time = 48 * 24 * 60 * 60;
            total_acceptance = (ele_acc_f+ele_acc_l);
            total_acceptance *=(pim_acc_f+pim_acc_l);
            total_acceptance *=(pro_acc_f+pro_acc_l);

            weight *=total_acceptance * time;
            weight_uu *=total_acceptance * time;
            weight_ut *=total_acceptance * time;

            if( t>t_cut[0] && t<t_cut[1])  t1->Fill();
            if( t>t_cut[1] && t<t_cut[2])  t2->Fill();
            if( t>t_cut[2] && t<t_cut[3])  t3->Fill();
            if( t>t_cut[3] && t<t_cut[4])  t4->Fill();
            if( t>t_cut[4] && t<t_cut[5])  t5->Fill();
            if( t>t_cut[5] && t<t_cut[6])  t6->Fill();
            if( t>t_cut[6] && t<t_cut[7])  t7->Fill();
            if( t>t_cut[7] && t<t_cut[8])  t8->Fill();
        }
    }


    f1->cd(); t1->Write(); f1->Close();
    f2->cd(); t2->Write(); f2->Close();
    f3->cd(); t3->Write(); f3->Close();
    f4->cd(); t4->Write(); f4->Close();
    f5->cd(); t5->Write(); f5->Close();
    f6->cd(); t6->Write(); f6->Close();
    f7->cd(); t7->Write(); f7->Close();
    f8->cd(); t8->Write(); f8->Close();
    file->Close();

    return 0;
}
