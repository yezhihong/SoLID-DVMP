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
Double_t max(Double_t a, Double_t b);
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
    TString target = "n";
    TString energy = "11";
    Int_t Q2BIN = 0;
    /*Initial Factors{{{*/
    if(target!="p" && target!="n"){
        cerr<<"I don't know this target flag!"<<endl;
        return -1;
    }
    if(energy!="11" && energy!="8" && energy!="11p8"){
        cerr<<"I don't know this energy flag!"<<endl;
        return -2;
    }

    const Double_t dilute_factor = 0.9;
    const Double_t target_polarization = 0.6; //60% polarization
    //const Double_t beam_polarization = 1.0; //60% polarization
    const Double_t det_eff = 0.85; //85% detector efficiency for electrons and hadrons
    const Double_t target_factor = 0.865; //neutron polarization 86.5%
    //	const Double_t Asys = 0.006; 
    //const Double_t Nsys = 2e+6; 
    /*}}}*/

    TString prefix = "";
    TString finalfile = Form("../rootfiles/new_DEMP_Ee_11_4.root");
    TFile *file = new TFile(finalfile.Data(),"r");
    TTree *t0 = (TTree*) file->Get("T");
   
    /*Define new Tree and new Branch{{{*/
    Double_t Epsilon, Qsq, T, W, x, y, z, vertexz;
    Double_t ele_ene,ele_px,ele_py,ele_pz,ele_mom; 
    Double_t ele_phi,ele_theta;
    Double_t pim_ene,pim_px,pim_py,pim_pz,pim_mom; 
    Double_t pim_phi,pim_theta;
    Double_t pro_ene,pro_px,pro_py,pro_pz,pro_mom; 
    Double_t pro_phi,pro_theta;
    Double_t ele_acc_f, ele_acc_l, pim_acc_f,pim_acc_l, pro_acc_f,pro_acc_l;
    Double_t t,MM, MM_res,dilute,weight;
    Double_t ele_mom_res, ele_theta_res, ele_phi_res;
    Double_t pim_mom_res, pim_theta_res, pim_phi_res;
    Double_t tgt_px,tgt_py,tgt_pz;//to simulate the Fermi motion
    Double_t phi_h, phi_S, ZASig_T, ZASig_L, ZASig_LT, ZASig_TT, SSAsym, SineAsym, ZASigma_Lab, EventWeight, PSF;
    Int_t N_Total;

    /*Old Tree{{{*/
    t0->SetBranchAddress("Epsilon", &Epsilon);
    t0->SetBranchAddress("Qsq", &Qsq);
    t0->SetBranchAddress("t", &t);
    t0->SetBranchAddress("W", &W);
    t0->SetBranchAddress("x", &x);
    t0->SetBranchAddress("y", &y);
    t0->SetBranchAddress("z", &z);
    t0->SetBranchAddress("phi_h", &phi_h );
    t0->SetBranchAddress("phi_S", &phi_S );

    t0->SetBranchAddress("ZASigma_Lab", &ZASigma_Lab    );                              
    t0->SetBranchAddress("EventWeight" ,&EventWeight    );                              
    t0->SetBranchAddress("ZASig_T",     &ZASig_T        );                                  
    t0->SetBranchAddress("ZASig_L",     &ZASig_L        );                                  
    t0->SetBranchAddress("ZASig_LT",     &ZASig_LT      );                                 
    t0->SetBranchAddress("ZASig_TT",     &ZASig_TT      );                                 
    t0->SetBranchAddress("SSAsym",     &SSAsym         );                                   
    t0->SetBranchAddress("SineAsym",     &SineAsym       );                                 
    
    t0->SetBranchAddress("vertexz",   &vertexz   );

    t0->SetBranchAddress("pim_ene",   &pim_ene   );
    t0->SetBranchAddress("pim_px",    &pim_px    );
    t0->SetBranchAddress("pim_py",    &pim_py    );
    t0->SetBranchAddress("pim_pz",    &pim_pz    );
    t0->SetBranchAddress("pim_mom",   &pim_mom   );
    t0->SetBranchAddress("pim_theta", &pim_theta );
    t0->SetBranchAddress("pim_phi",   &pim_phi   );

    t0->SetBranchAddress("ele_ene",   &ele_ene   );
    t0->SetBranchAddress("ele_px",    &ele_px    );
    t0->SetBranchAddress("ele_py",    &ele_py    );
    t0->SetBranchAddress("ele_pz",    &ele_pz    );
    t0->SetBranchAddress("ele_mom",   &ele_mom   );
    t0->SetBranchAddress("ele_theta", &ele_theta );
    t0->SetBranchAddress("ele_phi",   &ele_phi   );

    t0->SetBranchAddress("pro_ene",   &pro_ene   );
    t0->SetBranchAddress("pro_px",    &pro_px    );
    t0->SetBranchAddress("pro_py",    &pro_py    );
    t0->SetBranchAddress("pro_pz",    &pro_pz    );
    t0->SetBranchAddress("pro_mom",   &pro_mom   );
    t0->SetBranchAddress("pro_theta", &pro_theta );
    t0->SetBranchAddress("pro_phi",   &pro_phi   );

    t0->SetBranchAddress("ele_acc_f",     &ele_acc_f     );
    t0->SetBranchAddress("ele_acc_l",     &ele_acc_l     );
    t0->SetBranchAddress("pim_acc_f",     &pim_acc_f     );
    t0->SetBranchAddress("pim_acc_l",     &pim_acc_l     );
    t0->SetBranchAddress("pro_acc_f",     &pro_acc_f     );
    t0->SetBranchAddress("pro_acc_l",     &pro_acc_l     );
    t0->SetBranchAddress("ele_mom_res",   &ele_mom_res   );
    t0->SetBranchAddress("ele_theta_res", &ele_theta_res );
    t0->SetBranchAddress("ele_phi_res",   &ele_phi_res   );
    t0->SetBranchAddress("pim_mom_res",   &pim_mom_res   );
    t0->SetBranchAddress("pim_theta_res", &pim_theta_res );
    t0->SetBranchAddress("pim_phi_res",   &pim_phi_res   );
    t0->SetBranchAddress("weight",        &weight        );
    t0->SetBranchAddress("dilute",        &dilute        );
    t0->SetBranchAddress("PSF",           &PSF           );
    t0->SetBranchAddress("N_Total",       &N_Total       );
    t0->SetBranchAddress("MM",            &MM            );
    t0->SetBranchAddress("MM_res", &MM_res);
    t0->SetBranchAddress("tgt_px", &tgt_px);
    t0->SetBranchAddress("tgt_py", &tgt_py);
    t0->SetBranchAddress("tgt_pz", &tgt_pz);
    Long64_t N_entries = t0->GetEntries();
    /*}}}*/

    TFile *f1 = new TFile("../rootfiles/dvmp_t1.root", "recreate");/*{{{*/
    TTree *t1 = new TTree("T","a new tree");

    t1->Branch("Epsilon", &Epsilon, "Epsilon/D" );
    t1->Branch("Qsq", &Qsq ,"Qsq/D");
    t1->Branch("t", &t ,"t/D");
    t1->Branch("W", &W ,"W/D");
    t1->Branch("x", &x ,"x/D");
    t1->Branch("y", &y ,"y/D");
    t1->Branch("z", &z ,"z/D");
    t1->Branch("phi_h", &phi_h ,"phi_h/D");
    t1->Branch("phi_S", &phi_S ,"phi_S/D");

    t1->Branch("ZASigma_Lab",     &ZASigma_Lab, "data/D");                              
    t1->Branch("EventWeight",     &EventWeight, "data/D");                              
    t1->Branch("ZASig_T",         &ZASig_T, "data/D");                                  
    t1->Branch("ZASig_L",         &ZASig_L, "data/D");                                  
    t1->Branch("ZASig_LT",        &ZASig_LT, "data/D");                                 
    t1->Branch("ZASig_TT",        &ZASig_TT, "data/D");                                 
    t1->Branch("SSAsym",          &SSAsym, "data/D");                                   
    t1->Branch("SineAsym",        &SineAsym, "data/D");                                 
    
    t1->Branch("vertexz",   &vertexz   ,"vertexz/D");

    t1->Branch("pim_ene",   &pim_ene   ,"pim_ene/D");
    t1->Branch("pim_px",    &pim_px    ,"pim_px/D");
    t1->Branch("pim_py",    &pim_py    ,"pim_py/D");
    t1->Branch("pim_pz",    &pim_pz    ,"pim_pz/D");
    t1->Branch("pim_mom",   &pim_mom   ,"pim_mom/D");
    t1->Branch("pim_theta", &pim_theta ,"pim_theta/D");
    t1->Branch("pim_phi",   &pim_phi   ,"pim_phi/D");

    t1->Branch("ele_ene",   &ele_ene   ,"ele_ene/D");
    t1->Branch("ele_px",    &ele_px    ,"ele_px/D");
    t1->Branch("ele_py",    &ele_py    ,"ele_py/D");
    t1->Branch("ele_pz",    &ele_pz    ,"ele_pz/D");
    t1->Branch("ele_mom",   &ele_mom   ,"ele_mom/D");
    t1->Branch("ele_theta", &ele_theta ,"ele_theta/D");
    t1->Branch("ele_phi",   &ele_phi   ,"ele_phi/D");

    t1->Branch("pro_ene",   &pro_ene   ,"pro_ene/D");
    t1->Branch("pro_px",    &pro_px    ,"pro_px/D");
    t1->Branch("pro_py",    &pro_py    ,"pro_py/D");
    t1->Branch("pro_pz",    &pro_pz    ,"pro_pz/D");
    t1->Branch("pro_mom",   &pro_mom   ,"pro_mom/D");
    t1->Branch("pro_theta", &pro_theta ,"pro_theta/D");
    t1->Branch("pro_phi",   &pro_phi   ,"pro_phi/D");

    t1->Branch("ele_acc_f",     &ele_acc_f,     "ele_acc_f/D");
    t1->Branch("ele_acc_l",     &ele_acc_l,     "ele_acc_l/D");
    t1->Branch("pim_acc_f",     &pim_acc_f,     "pim_acc_f/D");
    t1->Branch("pim_acc_l",     &pim_acc_l,     "pim_acc_l/D");
    t1->Branch("pro_acc_f",     &pro_acc_f,     "pro_acc_f/D");
    t1->Branch("pro_acc_l",     &pro_acc_l,     "pro_acc_l/D");
    t1->Branch("ele_mom_res",   &ele_mom_res,   "ele_mom_res/D");
    t1->Branch("ele_theta_res", &ele_theta_res, "ele_theta_res/D");
    t1->Branch("ele_phi_res",   &ele_phi_res,   "ele_phi_res/D");
    t1->Branch("pim_mom_res",   &pim_mom_res,   "pim_mom_res/D");
    t1->Branch("pim_theta_res", &pim_theta_res, "pim_theta_res/D");
    t1->Branch("pim_phi_res",   &pim_phi_res,   "pim_phi_res/D");
    t1->Branch("weight",        &weight,        "weight/D");
    t1->Branch("dilute",        &dilute,        "dilute/D");
    t1->Branch("PSF",           &PSF,           "PSF/D");
    t1->Branch("N_Total",       &N_Total,       "N_Total/I");
    t1->Branch("MM",            &MM,            "MM/D");
    t1->Branch("MM_res", &MM_res, "MM_res/D");
    t1->Branch("tgt_px", &tgt_px, "tgt_px/D");
    t1->Branch("tgt_py", &tgt_py, "tgt_py/D");
    t1->Branch("tgt_pz", &tgt_pz, "tgt_pz/D");/*}}}*/
    
    TFile *f2 = new TFile("../rootfiles/dvmp_t2.root", "recreate");/*{{{*/
    TTree *t2 = new TTree("T","a new tree");
    t2->Branch("Epsilon", &Epsilon, "Epsilon/D" );
    t2->Branch("Qsq", &Qsq ,"Qsq/D");
    t2->Branch("t", &t ,"t/D");
    t2->Branch("W", &W ,"W/D");
    t2->Branch("x", &x ,"x/D");
    t2->Branch("y", &y ,"y/D");
    t2->Branch("z", &z ,"z/D");
    t2->Branch("phi_h", &phi_h ,"phi_h/D");
    t2->Branch("phi_S", &phi_S ,"phi_S/D");

    t2->Branch("ZASigma_Lab",     &ZASigma_Lab, "data/D");                              
    t2->Branch("EventWeight",     &EventWeight, "data/D");                              
    t2->Branch("ZASig_T",         &ZASig_T, "data/D");                                  
    t2->Branch("ZASig_L",         &ZASig_L, "data/D");                                  
    t2->Branch("ZASig_LT",        &ZASig_LT, "data/D");                                 
    t2->Branch("ZASig_TT",        &ZASig_TT, "data/D");                                 
    t2->Branch("SSAsym",          &SSAsym, "data/D");                                   
    t2->Branch("SineAsym",        &SineAsym, "data/D");                                 
    
    t2->Branch("vertexz",   &vertexz   ,"vertexz/D");

    t2->Branch("pim_ene",   &pim_ene   ,"pim_ene/D");
    t2->Branch("pim_px",    &pim_px    ,"pim_px/D");
    t2->Branch("pim_py",    &pim_py    ,"pim_py/D");
    t2->Branch("pim_pz",    &pim_pz    ,"pim_pz/D");
    t2->Branch("pim_mom",   &pim_mom   ,"pim_mom/D");
    t2->Branch("pim_theta", &pim_theta ,"pim_theta/D");
    t2->Branch("pim_phi",   &pim_phi   ,"pim_phi/D");

    t2->Branch("ele_ene",   &ele_ene   ,"ele_ene/D");
    t2->Branch("ele_px",    &ele_px    ,"ele_px/D");
    t2->Branch("ele_py",    &ele_py    ,"ele_py/D");
    t2->Branch("ele_pz",    &ele_pz    ,"ele_pz/D");
    t2->Branch("ele_mom",   &ele_mom   ,"ele_mom/D");
    t2->Branch("ele_theta", &ele_theta ,"ele_theta/D");
    t2->Branch("ele_phi",   &ele_phi   ,"ele_phi/D");

    t2->Branch("pro_ene",   &pro_ene   ,"pro_ene/D");
    t2->Branch("pro_px",    &pro_px    ,"pro_px/D");
    t2->Branch("pro_py",    &pro_py    ,"pro_py/D");
    t2->Branch("pro_pz",    &pro_pz    ,"pro_pz/D");
    t2->Branch("pro_mom",   &pro_mom   ,"pro_mom/D");
    t2->Branch("pro_theta", &pro_theta ,"pro_theta/D");
    t2->Branch("pro_phi",   &pro_phi   ,"pro_phi/D");

    t2->Branch("ele_acc_f",     &ele_acc_f,     "ele_acc_f/D");
    t2->Branch("ele_acc_l",     &ele_acc_l,     "ele_acc_l/D");
    t2->Branch("pim_acc_f",     &pim_acc_f,     "pim_acc_f/D");
    t2->Branch("pim_acc_l",     &pim_acc_l,     "pim_acc_l/D");
    t2->Branch("pro_acc_f",     &pro_acc_f,     "pro_acc_f/D");
    t2->Branch("pro_acc_l",     &pro_acc_l,     "pro_acc_l/D");
    t2->Branch("ele_mom_res",   &ele_mom_res,   "ele_mom_res/D");
    t2->Branch("ele_theta_res", &ele_theta_res, "ele_theta_res/D");
    t2->Branch("ele_phi_res",   &ele_phi_res,   "ele_phi_res/D");
    t2->Branch("pim_mom_res",   &pim_mom_res,   "pim_mom_res/D");
    t2->Branch("pim_theta_res", &pim_theta_res, "pim_theta_res/D");
    t2->Branch("pim_phi_res",   &pim_phi_res,   "pim_phi_res/D");
    t2->Branch("weight",        &weight,        "weight/D");
    t2->Branch("dilute",        &dilute,        "dilute/D");
    t2->Branch("PSF",           &PSF,           "PSF/D");
    t2->Branch("N_Total",       &N_Total,       "N_Total/I");
    t2->Branch("MM",            &MM,            "MM/D");
    t2->Branch("MM_res", &MM_res, "MM_res/D");
    t2->Branch("tgt_px", &tgt_px, "tgt_px/D");
    t2->Branch("tgt_py", &tgt_py, "tgt_py/D");
    t2->Branch("tgt_pz", &tgt_pz, "tgt_pz/D");/*}}}*/
    
    TFile *f3 = new TFile("../rootfiles/dvmp_t3.root", "recreate");/*{{{*/
    TTree *t3 = new TTree("T","a new tree");
    t3->Branch("Epsilon", &Epsilon, "Epsilon/D" );
    t3->Branch("Qsq", &Qsq ,"Qsq/D");
    t3->Branch("t", &t ,"t/D");
    t3->Branch("W", &W ,"W/D");
    t3->Branch("x", &x ,"x/D");
    t3->Branch("y", &y ,"y/D");
    t3->Branch("z", &z ,"z/D");
    t3->Branch("phi_h", &phi_h ,"phi_h/D");
    t3->Branch("phi_S", &phi_S ,"phi_S/D");

    t3->Branch("ZASigma_Lab",     &ZASigma_Lab, "data/D");                              
    t3->Branch("EventWeight",     &EventWeight, "data/D");                              
    t3->Branch("ZASig_T",         &ZASig_T, "data/D");                                  
    t3->Branch("ZASig_L",         &ZASig_L, "data/D");                                  
    t3->Branch("ZASig_LT",        &ZASig_LT, "data/D");                                 
    t3->Branch("ZASig_TT",        &ZASig_TT, "data/D");                                 
    t3->Branch("SSAsym",          &SSAsym, "data/D");                                   
    t3->Branch("SineAsym",        &SineAsym, "data/D");                                 
    
    t3->Branch("vertexz",   &vertexz   ,"vertexz/D");

    t3->Branch("pim_ene",   &pim_ene   ,"pim_ene/D");
    t3->Branch("pim_px",    &pim_px    ,"pim_px/D");
    t3->Branch("pim_py",    &pim_py    ,"pim_py/D");
    t3->Branch("pim_pz",    &pim_pz    ,"pim_pz/D");
    t3->Branch("pim_mom",   &pim_mom   ,"pim_mom/D");
    t3->Branch("pim_theta", &pim_theta ,"pim_theta/D");
    t3->Branch("pim_phi",   &pim_phi   ,"pim_phi/D");

    t3->Branch("ele_ene",   &ele_ene   ,"ele_ene/D");
    t3->Branch("ele_px",    &ele_px    ,"ele_px/D");
    t3->Branch("ele_py",    &ele_py    ,"ele_py/D");
    t3->Branch("ele_pz",    &ele_pz    ,"ele_pz/D");
    t3->Branch("ele_mom",   &ele_mom   ,"ele_mom/D");
    t3->Branch("ele_theta", &ele_theta ,"ele_theta/D");
    t3->Branch("ele_phi",   &ele_phi   ,"ele_phi/D");

    t3->Branch("pro_ene",   &pro_ene   ,"pro_ene/D");
    t3->Branch("pro_px",    &pro_px    ,"pro_px/D");
    t3->Branch("pro_py",    &pro_py    ,"pro_py/D");
    t3->Branch("pro_pz",    &pro_pz    ,"pro_pz/D");
    t3->Branch("pro_mom",   &pro_mom   ,"pro_mom/D");
    t3->Branch("pro_theta", &pro_theta ,"pro_theta/D");
    t3->Branch("pro_phi",   &pro_phi   ,"pro_phi/D");

    t3->Branch("ele_acc_f",     &ele_acc_f,     "ele_acc_f/D");
    t3->Branch("ele_acc_l",     &ele_acc_l,     "ele_acc_l/D");
    t3->Branch("pim_acc_f",     &pim_acc_f,     "pim_acc_f/D");
    t3->Branch("pim_acc_l",     &pim_acc_l,     "pim_acc_l/D");
    t3->Branch("pro_acc_f",     &pro_acc_f,     "pro_acc_f/D");
    t3->Branch("pro_acc_l",     &pro_acc_l,     "pro_acc_l/D");
    t3->Branch("ele_mom_res",   &ele_mom_res,   "ele_mom_res/D");
    t3->Branch("ele_theta_res", &ele_theta_res, "ele_theta_res/D");
    t3->Branch("ele_phi_res",   &ele_phi_res,   "ele_phi_res/D");
    t3->Branch("pim_mom_res",   &pim_mom_res,   "pim_mom_res/D");
    t3->Branch("pim_theta_res", &pim_theta_res, "pim_theta_res/D");
    t3->Branch("pim_phi_res",   &pim_phi_res,   "pim_phi_res/D");
    t3->Branch("weight",        &weight,        "weight/D");
    t3->Branch("dilute",        &dilute,        "dilute/D");
    t3->Branch("PSF",           &PSF,           "PSF/D");
    t3->Branch("N_Total",       &N_Total,       "N_Total/I");
    t3->Branch("MM",            &MM,            "MM/D");
    t3->Branch("MM_res", &MM_res, "MM_res/D");
    t3->Branch("tgt_px", &tgt_px, "tgt_px/D");
    t3->Branch("tgt_py", &tgt_py, "tgt_py/D");
    t3->Branch("tgt_pz", &tgt_pz, "tgt_pz/D");/*}}}*/

    TFile *f4 = new TFile("../rootfiles/dvmp_t4.root", "recreate");/*{{{*/
    TTree *t4 = new TTree("T","a new tree");
    t4->Branch("Epsilon", &Epsilon, "Epsilon/D" );
    t4->Branch("Qsq", &Qsq ,"Qsq/D");
    t4->Branch("t", &t ,"t/D");
    t4->Branch("W", &W ,"W/D");
    t4->Branch("x", &x ,"x/D");
    t4->Branch("y", &y ,"y/D");
    t4->Branch("z", &z ,"z/D");
    t4->Branch("phi_h", &phi_h ,"phi_h/D");
    t4->Branch("phi_S", &phi_S ,"phi_S/D");

    t4->Branch("ZASigma_Lab",     &ZASigma_Lab, "data/D");                              
    t4->Branch("EventWeight",     &EventWeight, "data/D");                              
    t4->Branch("ZASig_T",         &ZASig_T, "data/D");                                  
    t4->Branch("ZASig_L",         &ZASig_L, "data/D");                                  
    t4->Branch("ZASig_LT",        &ZASig_LT, "data/D");                                 
    t4->Branch("ZASig_TT",        &ZASig_TT, "data/D");                                 
    t4->Branch("SSAsym",          &SSAsym, "data/D");                                   
    t4->Branch("SineAsym",        &SineAsym, "data/D");                                 
    
    t4->Branch("vertexz",   &vertexz   ,"vertexz/D");

    t4->Branch("pim_ene",   &pim_ene   ,"pim_ene/D");
    t4->Branch("pim_px",    &pim_px    ,"pim_px/D");
    t4->Branch("pim_py",    &pim_py    ,"pim_py/D");
    t4->Branch("pim_pz",    &pim_pz    ,"pim_pz/D");
    t4->Branch("pim_mom",   &pim_mom   ,"pim_mom/D");
    t4->Branch("pim_theta", &pim_theta ,"pim_theta/D");
    t4->Branch("pim_phi",   &pim_phi   ,"pim_phi/D");

    t4->Branch("ele_ene",   &ele_ene   ,"ele_ene/D");
    t4->Branch("ele_px",    &ele_px    ,"ele_px/D");
    t4->Branch("ele_py",    &ele_py    ,"ele_py/D");
    t4->Branch("ele_pz",    &ele_pz    ,"ele_pz/D");
    t4->Branch("ele_mom",   &ele_mom   ,"ele_mom/D");
    t4->Branch("ele_theta", &ele_theta ,"ele_theta/D");
    t4->Branch("ele_phi",   &ele_phi   ,"ele_phi/D");

    t4->Branch("pro_ene",   &pro_ene   ,"pro_ene/D");
    t4->Branch("pro_px",    &pro_px    ,"pro_px/D");
    t4->Branch("pro_py",    &pro_py    ,"pro_py/D");
    t4->Branch("pro_pz",    &pro_pz    ,"pro_pz/D");
    t4->Branch("pro_mom",   &pro_mom   ,"pro_mom/D");
    t4->Branch("pro_theta", &pro_theta ,"pro_theta/D");
    t4->Branch("pro_phi",   &pro_phi   ,"pro_phi/D");

    t4->Branch("ele_acc_f",     &ele_acc_f,     "ele_acc_f/D");
    t4->Branch("ele_acc_l",     &ele_acc_l,     "ele_acc_l/D");
    t4->Branch("pim_acc_f",     &pim_acc_f,     "pim_acc_f/D");
    t4->Branch("pim_acc_l",     &pim_acc_l,     "pim_acc_l/D");
    t4->Branch("pro_acc_f",     &pro_acc_f,     "pro_acc_f/D");
    t4->Branch("pro_acc_l",     &pro_acc_l,     "pro_acc_l/D");
    t4->Branch("ele_mom_res",   &ele_mom_res,   "ele_mom_res/D");
    t4->Branch("ele_theta_res", &ele_theta_res, "ele_theta_res/D");
    t4->Branch("ele_phi_res",   &ele_phi_res,   "ele_phi_res/D");
    t4->Branch("pim_mom_res",   &pim_mom_res,   "pim_mom_res/D");
    t4->Branch("pim_theta_res", &pim_theta_res, "pim_theta_res/D");
    t4->Branch("pim_phi_res",   &pim_phi_res,   "pim_phi_res/D");
    t4->Branch("weight",        &weight,        "weight/D");
    t4->Branch("dilute",        &dilute,        "dilute/D");
    t4->Branch("PSF",           &PSF,           "PSF/D");
    t4->Branch("N_Total",       &N_Total,       "N_Total/I");
    t4->Branch("MM",            &MM,            "MM/D");
    t4->Branch("MM_res", &MM_res, "MM_res/D");
    t4->Branch("tgt_px", &tgt_px, "tgt_px/D");
    t4->Branch("tgt_py", &tgt_py, "tgt_py/D");
    t4->Branch("tgt_pz", &tgt_pz, "tgt_pz/D");/*}}}*/
    
    TFile *f5 = new TFile("../rootfiles/dvmp_t5.root", "recreate");/*{{{*/
    TTree *t5 = new TTree("T","a new tree");
    t5->Branch("Epsilon", &Epsilon, "Epsilon/D" );
    t5->Branch("Qsq", &Qsq ,"Qsq/D");
    t5->Branch("t", &t ,"t/D");
    t5->Branch("W", &W ,"W/D");
    t5->Branch("x", &x ,"x/D");
    t5->Branch("y", &y ,"y/D");
    t5->Branch("z", &z ,"z/D");
    t5->Branch("phi_h", &phi_h ,"phi_h/D");
    t5->Branch("phi_S", &phi_S ,"phi_S/D");

    t5->Branch("ZASigma_Lab",     &ZASigma_Lab, "data/D");                              
    t5->Branch("EventWeight",     &EventWeight, "data/D");                              
    t5->Branch("ZASig_T",         &ZASig_T, "data/D");                                  
    t5->Branch("ZASig_L",         &ZASig_L, "data/D");                                  
    t5->Branch("ZASig_LT",        &ZASig_LT, "data/D");                                 
    t5->Branch("ZASig_TT",        &ZASig_TT, "data/D");                                 
    t5->Branch("SSAsym",          &SSAsym, "data/D");                                   
    t5->Branch("SineAsym",        &SineAsym, "data/D");                                 
    
    t5->Branch("vertexz",   &vertexz   ,"vertexz/D");

    t5->Branch("pim_ene",   &pim_ene   ,"pim_ene/D");
    t5->Branch("pim_px",    &pim_px    ,"pim_px/D");
    t5->Branch("pim_py",    &pim_py    ,"pim_py/D");
    t5->Branch("pim_pz",    &pim_pz    ,"pim_pz/D");
    t5->Branch("pim_mom",   &pim_mom   ,"pim_mom/D");
    t5->Branch("pim_theta", &pim_theta ,"pim_theta/D");
    t5->Branch("pim_phi",   &pim_phi   ,"pim_phi/D");

    t5->Branch("ele_ene",   &ele_ene   ,"ele_ene/D");
    t5->Branch("ele_px",    &ele_px    ,"ele_px/D");
    t5->Branch("ele_py",    &ele_py    ,"ele_py/D");
    t5->Branch("ele_pz",    &ele_pz    ,"ele_pz/D");
    t5->Branch("ele_mom",   &ele_mom   ,"ele_mom/D");
    t5->Branch("ele_theta", &ele_theta ,"ele_theta/D");
    t5->Branch("ele_phi",   &ele_phi   ,"ele_phi/D");

    t5->Branch("pro_ene",   &pro_ene   ,"pro_ene/D");
    t5->Branch("pro_px",    &pro_px    ,"pro_px/D");
    t5->Branch("pro_py",    &pro_py    ,"pro_py/D");
    t5->Branch("pro_pz",    &pro_pz    ,"pro_pz/D");
    t5->Branch("pro_mom",   &pro_mom   ,"pro_mom/D");
    t5->Branch("pro_theta", &pro_theta ,"pro_theta/D");
    t5->Branch("pro_phi",   &pro_phi   ,"pro_phi/D");

    t5->Branch("ele_acc_f",     &ele_acc_f,     "ele_acc_f/D");
    t5->Branch("ele_acc_l",     &ele_acc_l,     "ele_acc_l/D");
    t5->Branch("pim_acc_f",     &pim_acc_f,     "pim_acc_f/D");
    t5->Branch("pim_acc_l",     &pim_acc_l,     "pim_acc_l/D");
    t5->Branch("pro_acc_f",     &pro_acc_f,     "pro_acc_f/D");
    t5->Branch("pro_acc_l",     &pro_acc_l,     "pro_acc_l/D");
    t5->Branch("ele_mom_res",   &ele_mom_res,   "ele_mom_res/D");
    t5->Branch("ele_theta_res", &ele_theta_res, "ele_theta_res/D");
    t5->Branch("ele_phi_res",   &ele_phi_res,   "ele_phi_res/D");
    t5->Branch("pim_mom_res",   &pim_mom_res,   "pim_mom_res/D");
    t5->Branch("pim_theta_res", &pim_theta_res, "pim_theta_res/D");
    t5->Branch("pim_phi_res",   &pim_phi_res,   "pim_phi_res/D");
    t5->Branch("weight",        &weight,        "weight/D");
    t5->Branch("dilute",        &dilute,        "dilute/D");
    t5->Branch("PSF",           &PSF,           "PSF/D");
    t5->Branch("N_Total",       &N_Total,       "N_Total/I");
    t5->Branch("MM",            &MM,            "MM/D");
    t5->Branch("MM_res", &MM_res, "MM_res/D");
    t5->Branch("tgt_px", &tgt_px, "tgt_px/D");
    t5->Branch("tgt_py", &tgt_py, "tgt_py/D");
    t5->Branch("tgt_pz", &tgt_pz, "tgt_pz/D");/*}}}*/
   
    TFile *f6 = new TFile("../rootfiles/dvmp_t6.root", "recreate");/*{{{*/
    TTree *t6 = new TTree("T","a new tree");
    t6->Branch("Epsilon", &Epsilon, "Epsilon/D" );
    t6->Branch("Qsq", &Qsq ,"Qsq/D");
    t6->Branch("t", &t ,"t/D");
    t6->Branch("W", &W ,"W/D");
    t6->Branch("x", &x ,"x/D");
    t6->Branch("y", &y ,"y/D");
    t6->Branch("z", &z ,"z/D");
    t6->Branch("phi_h", &phi_h ,"phi_h/D");
    t6->Branch("phi_S", &phi_S ,"phi_S/D");

    t6->Branch("ZASigma_Lab",     &ZASigma_Lab, "data/D");                              
    t6->Branch("EventWeight",     &EventWeight, "data/D");                              
    t6->Branch("ZASig_T",         &ZASig_T, "data/D");                                  
    t6->Branch("ZASig_L",         &ZASig_L, "data/D");                                  
    t6->Branch("ZASig_LT",        &ZASig_LT, "data/D");                                 
    t6->Branch("ZASig_TT",        &ZASig_TT, "data/D");                                 
    t6->Branch("SSAsym",          &SSAsym, "data/D");                                   
    t6->Branch("SineAsym",        &SineAsym, "data/D");                                 
    
    t6->Branch("vertexz",   &vertexz   ,"vertexz/D");

    t6->Branch("pim_ene",   &pim_ene   ,"pim_ene/D");
    t6->Branch("pim_px",    &pim_px    ,"pim_px/D");
    t6->Branch("pim_py",    &pim_py    ,"pim_py/D");
    t6->Branch("pim_pz",    &pim_pz    ,"pim_pz/D");
    t6->Branch("pim_mom",   &pim_mom   ,"pim_mom/D");
    t6->Branch("pim_theta", &pim_theta ,"pim_theta/D");
    t6->Branch("pim_phi",   &pim_phi   ,"pim_phi/D");

    t6->Branch("ele_ene",   &ele_ene   ,"ele_ene/D");
    t6->Branch("ele_px",    &ele_px    ,"ele_px/D");
    t6->Branch("ele_py",    &ele_py    ,"ele_py/D");
    t6->Branch("ele_pz",    &ele_pz    ,"ele_pz/D");
    t6->Branch("ele_mom",   &ele_mom   ,"ele_mom/D");
    t6->Branch("ele_theta", &ele_theta ,"ele_theta/D");
    t6->Branch("ele_phi",   &ele_phi   ,"ele_phi/D");

    t6->Branch("pro_ene",   &pro_ene   ,"pro_ene/D");
    t6->Branch("pro_px",    &pro_px    ,"pro_px/D");
    t6->Branch("pro_py",    &pro_py    ,"pro_py/D");
    t6->Branch("pro_pz",    &pro_pz    ,"pro_pz/D");
    t6->Branch("pro_mom",   &pro_mom   ,"pro_mom/D");
    t6->Branch("pro_theta", &pro_theta ,"pro_theta/D");
    t6->Branch("pro_phi",   &pro_phi   ,"pro_phi/D");

    t6->Branch("ele_acc_f",     &ele_acc_f,     "ele_acc_f/D");
    t6->Branch("ele_acc_l",     &ele_acc_l,     "ele_acc_l/D");
    t6->Branch("pim_acc_f",     &pim_acc_f,     "pim_acc_f/D");
    t6->Branch("pim_acc_l",     &pim_acc_l,     "pim_acc_l/D");
    t6->Branch("pro_acc_f",     &pro_acc_f,     "pro_acc_f/D");
    t6->Branch("pro_acc_l",     &pro_acc_l,     "pro_acc_l/D");
    t6->Branch("ele_mom_res",   &ele_mom_res,   "ele_mom_res/D");
    t6->Branch("ele_theta_res", &ele_theta_res, "ele_theta_res/D");
    t6->Branch("ele_phi_res",   &ele_phi_res,   "ele_phi_res/D");
    t6->Branch("pim_mom_res",   &pim_mom_res,   "pim_mom_res/D");
    t6->Branch("pim_theta_res", &pim_theta_res, "pim_theta_res/D");
    t6->Branch("pim_phi_res",   &pim_phi_res,   "pim_phi_res/D");
    t6->Branch("weight",        &weight,        "weight/D");
    t6->Branch("dilute",        &dilute,        "dilute/D");
    t6->Branch("PSF",           &PSF,           "PSF/D");
    t6->Branch("N_Total",       &N_Total,       "N_Total/I");
    t6->Branch("MM",            &MM,            "MM/D");
    t6->Branch("MM_res", &MM_res, "MM_res/D");
    t6->Branch("tgt_px", &tgt_px, "tgt_px/D");
    t6->Branch("tgt_py", &tgt_py, "tgt_py/D");
    t6->Branch("tgt_pz", &tgt_pz, "tgt_pz/D");/*}}}*/

    TFile *f7 = new TFile("../rootfiles/dvmp_t7.root", "recreate");/*{{{*/
    TTree *t7 = new TTree("T","a new tree");
    t7->Branch("Epsilon", &Epsilon, "Epsilon/D" );
    t7->Branch("Qsq", &Qsq ,"Qsq/D");
    t7->Branch("t", &t ,"t/D");
    t7->Branch("W", &W ,"W/D");
    t7->Branch("x", &x ,"x/D");
    t7->Branch("y", &y ,"y/D");
    t7->Branch("z", &z ,"z/D");
    t7->Branch("phi_h", &phi_h ,"phi_h/D");
    t7->Branch("phi_S", &phi_S ,"phi_S/D");

    t7->Branch("ZASigma_Lab",     &ZASigma_Lab, "data/D");                              
    t7->Branch("EventWeight",     &EventWeight, "data/D");                              
    t7->Branch("ZASig_T",         &ZASig_T, "data/D");                                  
    t7->Branch("ZASig_L",         &ZASig_L, "data/D");                                  
    t7->Branch("ZASig_LT",        &ZASig_LT, "data/D");                                 
    t7->Branch("ZASig_TT",        &ZASig_TT, "data/D");                                 
    t7->Branch("SSAsym",          &SSAsym, "data/D");                                   
    t7->Branch("SineAsym",        &SineAsym, "data/D");                                 
    
    t7->Branch("vertexz",   &vertexz   ,"vertexz/D");

    t7->Branch("pim_ene",   &pim_ene   ,"pim_ene/D");
    t7->Branch("pim_px",    &pim_px    ,"pim_px/D");
    t7->Branch("pim_py",    &pim_py    ,"pim_py/D");
    t7->Branch("pim_pz",    &pim_pz    ,"pim_pz/D");
    t7->Branch("pim_mom",   &pim_mom   ,"pim_mom/D");
    t7->Branch("pim_theta", &pim_theta ,"pim_theta/D");
    t7->Branch("pim_phi",   &pim_phi   ,"pim_phi/D");

    t7->Branch("ele_ene",   &ele_ene   ,"ele_ene/D");
    t7->Branch("ele_px",    &ele_px    ,"ele_px/D");
    t7->Branch("ele_py",    &ele_py    ,"ele_py/D");
    t7->Branch("ele_pz",    &ele_pz    ,"ele_pz/D");
    t7->Branch("ele_mom",   &ele_mom   ,"ele_mom/D");
    t7->Branch("ele_theta", &ele_theta ,"ele_theta/D");
    t7->Branch("ele_phi",   &ele_phi   ,"ele_phi/D");

    t7->Branch("pro_ene",   &pro_ene   ,"pro_ene/D");
    t7->Branch("pro_px",    &pro_px    ,"pro_px/D");
    t7->Branch("pro_py",    &pro_py    ,"pro_py/D");
    t7->Branch("pro_pz",    &pro_pz    ,"pro_pz/D");
    t7->Branch("pro_mom",   &pro_mom   ,"pro_mom/D");
    t7->Branch("pro_theta", &pro_theta ,"pro_theta/D");
    t7->Branch("pro_phi",   &pro_phi   ,"pro_phi/D");

    t7->Branch("ele_acc_f",     &ele_acc_f,     "ele_acc_f/D");
    t7->Branch("ele_acc_l",     &ele_acc_l,     "ele_acc_l/D");
    t7->Branch("pim_acc_f",     &pim_acc_f,     "pim_acc_f/D");
    t7->Branch("pim_acc_l",     &pim_acc_l,     "pim_acc_l/D");
    t7->Branch("pro_acc_f",     &pro_acc_f,     "pro_acc_f/D");
    t7->Branch("pro_acc_l",     &pro_acc_l,     "pro_acc_l/D");
    t7->Branch("ele_mom_res",   &ele_mom_res,   "ele_mom_res/D");
    t7->Branch("ele_theta_res", &ele_theta_res, "ele_theta_res/D");
    t7->Branch("ele_phi_res",   &ele_phi_res,   "ele_phi_res/D");
    t7->Branch("pim_mom_res",   &pim_mom_res,   "pim_mom_res/D");
    t7->Branch("pim_theta_res", &pim_theta_res, "pim_theta_res/D");
    t7->Branch("pim_phi_res",   &pim_phi_res,   "pim_phi_res/D");
    t7->Branch("weight",        &weight,        "weight/D");
    t7->Branch("dilute",        &dilute,        "dilute/D");
    t7->Branch("PSF",           &PSF,           "PSF/D");
    t7->Branch("N_Total",       &N_Total,       "N_Total/I");
    t7->Branch("MM",            &MM,            "MM/D");
    t7->Branch("MM_res", &MM_res, "MM_res/D");
    t7->Branch("tgt_px", &tgt_px, "tgt_px/D");
    t7->Branch("tgt_py", &tgt_py, "tgt_py/D");
    t7->Branch("tgt_pz", &tgt_pz, "tgt_pz/D");/*}}}*/
    /*}}}*/
    
    const Int_t tbin = 7;
    const Double_t t_cut[8] = {0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.75, 1.2};
    const Int_t time = 48 * 24 * 3600;
    
    const Double_t Norm_Fact = time * (pow(target_factor * dilute_factor,2) * det_eff);
  
    for(Long64_t i=0;i<N_entries;i++){
        t0->GetEntry(i);
        //the weight now is the combination of weight*acceptance*lumi
        weight = EventWeight*(ele_acc_f+ele_acc_l)*(pim_acc_f+pim_acc_l)*(pro_acc_f+pro_acc_l)*Norm_Fact;
        if(Epsilon>0.55&&Epsilon<0.75&&W>2&&Qsq>4){

            if( t>t_cut[0] && t<t_cut[1])  t1->Fill();
            if( t>t_cut[1] && t<t_cut[2])  t2->Fill();
            if( t>t_cut[2] && t<t_cut[3])  t3->Fill();
            if( t>t_cut[3] && t<t_cut[4])  t4->Fill();
            if( t>t_cut[4] && t<t_cut[5])  t5->Fill();
            if( t>t_cut[5] && t<t_cut[6])  t6->Fill();
            if( t>t_cut[6] && t<t_cut[7])  t7->Fill();
        }
    }


    f1->cd(); t1->Write(); f1->Close();
    f2->cd(); t2->Write(); f2->Close();
    f3->cd(); t3->Write(); f3->Close();
    f4->cd(); t4->Write(); f4->Close();
    f5->cd(); t5->Write(); f5->Close();
    f6->cd(); t6->Write(); f6->Close();
    f7->cd(); t7->Write(); f7->Close();
    file->Close();

    return 0;
}

Double_t max(Double_t a, Double_t b){/*{{{*/
    if(a>b)
        return a;
    else
        return b;
}/*}}}*/
