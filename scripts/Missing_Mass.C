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
const double Sigma_EC_E = 0.02; //2%, energy resolution for electron from GEM tracking
const double Sigma_EC_G = 0.05; //10%, energy resolution for photon
const double Sigma_Theta_E = 0.6/1000.; //0.6mrad, Angular resolution for electron, determined by GEM tracking
const double Sigma_Phi_E =  5.0/1000.; //5mrad, Angular resolution for electron, determined by GEM tracking
const double Sigma_X_G = 1.0; //cm, Position resolution for photon, determined by EC cluster reconstruction 
const double Sigma_Y_G = 1.0; //cm, Position resolution for photon, determined by EC cluster reconstruction
const double Sigma_VZ = 0.5; //cm, VertexZ resolution for photon, determined by electron GEM tracking
const double Length = 790.0; //cm, FAEC to target distance

const double eMass = 0.511/1000;//electron mass 
const double piMass = 0.13957018;
const double pMass = 0.938272;
const double nMass = 0.939565;

Int_t CheckLaws(TLorentzVector* P_t,TLorentzVector* P_E0, TLorentzVector* P_e, TLorentzVector* P_g, TLorentzVector* P_h);
Double_t GetMM( TLorentzVector* P_t,TLorentzVector* P_E0, TLorentzVector* P_e, TLorentzVector* P_g);
/*}}}*/

void Missing_Mass(){
//int main(){
	gStyle->SetOptStat(0);
	const double Lumi = 1.0e36; // cm-2*s-1, for He3 nuclear not for nucleons
	const double KHz = 1e-3;
	const double nBcm2 = 1e-33;
	const double PI = 3.1415926;
    const double DEG2RAD = PI/180.0;
	const TString Target = "He3";
    const TString particle="pim";
    const double hMass = nMass; //GeV, Proton Mass

    //CLEO Acceptance
    SIDIS_Acceptance *accpt = new SIDIS_Acceptance();
       
    //const Double_t EBeam = 11.0;
    //const double Ratio = 20.7598; //11GeV
    //TFile* f1=new TFile("./RootFiles/SoLID_Excl_Targ_N_Ee_11_Events_20e6.root","r"); 
	
    const Double_t EBeam = 8.80;
    const double Ratio = 25.1826; //8.8GeV
    TFile* f1=new TFile("./RootFiles/SoLID_Excl_Targ_N_Ee_8P8_Events_20e6.root","r"); 
    
    int energy_flag =int(EBeam);//GeV
	/*Define DVMP Rootfile and variables{{{*/
    TTree* t1 = (TTree*) f1->Get("t1");

    Double_t Epsilon, Qsq, T, W, x, y, z, ZASig_L, ZASig_T, ZASigmaPara,vertexz;
    Double_t Neutron_Energy_Col,Neutron_MomX_Col,Neutron_MomY_Col,Neutron_MomZ_Col,Neutron_Mom_Col; 
    Double_t Neutron_Phi_Col,Neutron_Theta_Col;
    Double_t ScatElec_Energy_Col,ScatElec_MomX_Col,ScatElec_MomY_Col,ScatElec_MomZ_Col,ScatElec_Mom_Col; 
    Double_t ScatElec_Phi_Col,ScatElec_Theta_Col;
    Double_t Pion_Energy_Col,Pion_MomX_Col,Pion_MomY_Col,Pion_MomZ_Col,Pion_Mom_Col; 
    Double_t Pion_Phi_Col,Pion_Theta_Col;
    Double_t Proton_Energy_Col,Proton_MomX_Col,Proton_MomY_Col,Proton_MomZ_Col,Proton_Mom_Col; 
    Double_t Proton_Phi_Col,Proton_Theta_Col;

    t1->SetBranchAddress("Epsilon", &Epsilon );
    t1->SetBranchAddress("Qsq", &Qsq );
    t1->SetBranchAddress("T", &T );
    t1->SetBranchAddress("W", &W );
    t1->SetBranchAddress("x", &x );
    t1->SetBranchAddress("y", &y );
    t1->SetBranchAddress("z", &z );
    t1->SetBranchAddress("ZASig_L", &ZASig_L );
    t1->SetBranchAddress("ZASig_T", &ZASig_T );
    t1->SetBranchAddress("ZASigmaPara", &ZASigmaPara );

    t1->SetBranchAddress("Vertex_Z", &vertexz );
    t1->SetBranchAddress("Pion_Energy_Col", &Pion_Energy_Col );
    t1->SetBranchAddress("Pion_MomX_Col", &Pion_MomX_Col );
    t1->SetBranchAddress("Pion_MomY_Col", &Pion_MomY_Col );
    t1->SetBranchAddress("Pion_MomZ_Col", &Pion_MomZ_Col );
    t1->SetBranchAddress("Pion_Mom_Col", &Pion_Mom_Col );
    t1->SetBranchAddress("Pion_Theta_Col", &Pion_Theta_Col );
    t1->SetBranchAddress("Pion_Phi_Col", &Pion_Phi_Col);

    t1->SetBranchAddress("ScatElec_Energy_Col", &ScatElec_Energy_Col );
    t1->SetBranchAddress("ScatElec_MomX_Col", &ScatElec_MomX_Col );
    t1->SetBranchAddress("ScatElec_MomY_Col", &ScatElec_MomY_Col );
    t1->SetBranchAddress("ScatElec_MomZ_Col", &ScatElec_MomZ_Col );
    t1->SetBranchAddress("ScatElec_Mom_Col", &ScatElec_Mom_Col );
    t1->SetBranchAddress("ScatElec_Theta_Col", &ScatElec_Theta_Col );
    t1->SetBranchAddress("ScatElec_Phi_Col", &ScatElec_Phi_Col );
	
    t1->SetBranchAddress("Neutron_Energy_Col", &Proton_Energy_Col );
    t1->SetBranchAddress("Neutron_MomX_Col", &Proton_MomX_Col );
    t1->SetBranchAddress("Neutron_MomY_Col", &Proton_MomY_Col );
    t1->SetBranchAddress("Neutron_MomZ_Col", &Proton_MomZ_Col );
    t1->SetBranchAddress("Neutron_Mom_Col", &Proton_Mom_Col );
    t1->SetBranchAddress("Neutron_Theta_Col", &Proton_Theta_Col );
    t1->SetBranchAddress("Neutron_Phi_Col", &Proton_Phi_Col );

   double PSF =(EBeam*0.9-EBeam*0.1) * (28.0-5.0)*DEG2RAD * 2*PI
               *(EBeam*0.9-EBeam*0.1) * (60.0-0.0)*DEG2RAD * 2*PI;
	Long64_t N_entries=t1->GetEntries();
    Long64_t N_Total = N_entries * Ratio; 
    /*}}}*/
    
    /*Vectors & Histograms{{{*/
	TLorentzVector *P_E0 = new TLorentzVector();//incoming electron
	TLorentzVector *P_e = new TLorentzVector();//scattered electron with Eloss
	TLorentzVector *P_pim = new TLorentzVector();//photon
	TLorentzVector *P_h = new TLorentzVector();//proton or neutron
	TLorentzVector *P_t = new TLorentzVector();//target, either proton or neutron
	TLorentzVector *P_e_res = new TLorentzVector();//scattered electron with resolution
	TLorentzVector *P_pim_res = new TLorentzVector();//photon with resolution

    TH1F *hMM_sidis = new TH1F("hMM_sidis",Form("SIDIS Missing Mass of %s at %3.1f GeV", particle.Data(),EBeam), 200,0.0, 2.8);
	hMM_sidis->SetXTitle("Hadron Missing Mass (GeV)");
	hMM_sidis->GetXaxis()->CenterTitle(1);
	hMM_sidis->SetYTitle("Rate (Hz)");
	hMM_sidis->GetYaxis()->CenterTitle(1);
	TH1F *hMM_sidis_res = new TH1F("hMM_sidis_res",Form("SIDIS Missing Mass of %s at %3.1f GeV", particle.Data(),EBeam), 200,0.0, 2.8);
	hMM_sidis_res->SetXTitle("Hadron Missing Mass (GeV)");
	hMM_sidis_res->GetXaxis()->CenterTitle(1);
	hMM_sidis_res->SetYTitle("Rate (Hz)");
	hMM_sidis_res->GetYaxis()->CenterTitle(1);
    
	TH1F *hMM_dvmp = new TH1F("hMM_dvmp",Form("DVMP Missing Mass at %3.1f GeV",EBeam), 200,0.0, 2.8);
	hMM_dvmp->SetXTitle("Hadron Missing Mass (GeV)");
	hMM_dvmp->GetXaxis()->CenterTitle(1);
	hMM_dvmp->SetYTitle("Rate (Hz)");
	hMM_dvmp->GetYaxis()->CenterTitle(1);
	TH1F *hMM_dvmp_res = new TH1F("hMM_dvmp_res",Form("DVMP Missing Mass at %3.1f GeV",EBeam), 200,0.0, 2.8);
	hMM_dvmp_res->SetXTitle("Hadron Missing Mass (GeV)");
	hMM_dvmp_res->GetXaxis()->CenterTitle(1);
	hMM_dvmp_res->SetYTitle("Rate (Hz)");
	hMM_dvmp_res->GetYaxis()->CenterTitle(1);
    /*}}}*/

    /*Fill in DVMP Missing Mass{{{*/
    double total_rate_dvmp=0.0; //Make sure the value is right
    for(Long64_t i=0;i<N_entries;i++){
		t1->GetEntry(i);
        ScatElec_Mom_Col*=MeV2GeV;//into GeV
        ScatElec_MomX_Col*=MeV2GeV;//into GeV
        ScatElec_MomY_Col*=MeV2GeV;//into GeV
        ScatElec_MomZ_Col*=MeV2GeV;//into GeV
        ScatElec_Energy_Col*=MeV2GeV;//into GeV
        Proton_Mom_Col*=MeV2GeV;//into GeV
        Proton_MomX_Col*=MeV2GeV;//into GeV
        Proton_MomY_Col*=MeV2GeV;//into GeV
        Proton_MomZ_Col*=MeV2GeV;//into GeV
        Proton_Energy_Col*=MeV2GeV;//into GeV
        Pion_Mom_Col*=MeV2GeV;//into GeV
        Pion_MomX_Col*=MeV2GeV;//into GeV
        Pion_MomY_Col*=MeV2GeV;//into GeV
        Pion_MomZ_Col*=MeV2GeV;//into GeV
        Pion_Energy_Col*=MeV2GeV;//into GeV
        Qsq *= MeV2GeV * MeV2GeV;
        T *= MeV2GeV * MeV2GeV;
        W *= MeV2GeV;

        if(ScatElec_Theta_Col>7.5&&ScatElec_Theta_Col<24.5&&ScatElec_Mom_Col>1.&&ScatElec_Mom_Col<11
          &&Pion_Theta_Col>7.5&&Pion_Theta_Col<24.5&&Pion_Mom_Col>1.&&Pion_Mom_Col<11
          && W>=2.&&Qsq>1.0&&ZASigmaPara>1e-33){//any additional cuts should be added in here
            /*Get acceptance of e and pi-{{{*/
            double ele_forward_acceptance = accpt->GetAcc("e-","forward", ScatElec_Mom_Col, ScatElec_Theta_Col);
            double ele_large_acceptance = accpt->GetAcc("e-","large", ScatElec_Mom_Col, ScatElec_Theta_Col);
			if(ScatElec_Mom_Col<1.0||ScatElec_Theta_Col>14.5||ScatElec_Theta_Col<8.0)//GeV, CLEO
				ele_forward_acceptance=0.0;//Farward-Angle EC Cut at 1 GeV
			if(ScatElec_Mom_Col<3.5||ScatElec_Theta_Col<16.0||ScatElec_Theta_Col>24)//GeV,CLEO
				ele_large_acceptance=0.0; //Larger-Angle EC Cut at 3 GeV
			if(ele_forward_acceptance>1.) 
				ele_forward_acceptance=1.0; 
			if(ele_large_acceptance>1.) 
				ele_large_acceptance=1.0; 

            double pim_forward_acceptance = accpt->GetAcc("pi-","forward", Pion_Mom_Col, Pion_Theta_Col);
            double pim_large_acceptance = accpt->GetAcc("pi-","large", Pion_Mom_Col, Pion_Theta_Col);
            if(Pion_Theta_Col>14.8||Pion_Theta_Col<8.0||Pion_Mom_Col<1.||Pion_Mom_Col>11.)//GeV, CLEO
                pim_forward_acceptance=0.0;
			if(Pion_Theta_Col<16.0||Pion_Theta_Col>24.0||Pion_Mom_Col<3.5||Pion_Mom_Col>11.)//GeV, CLEO
				pim_large_acceptance=0.0; 
			if(pim_forward_acceptance>1.) 
				pim_forward_acceptance=1.0; 
			if(pim_large_acceptance>1.) 
				pim_large_acceptance=1.0; 
            
            double event_weight=ZASigmaPara*PSF/N_Total*nBcm2*Lumi ;   //in Hz

			double ele_acceptance=(ele_forward_acceptance+ele_large_acceptance);
			//double pim_acceptance=(pim_large_acceptance+pim_forward_acceptance);
			double pim_acceptance=pim_forward_acceptance;
			double forward_acceptance=ele_forward_acceptance*pim_acceptance;
			double large_acceptance=ele_large_acceptance*pim_acceptance;
			double total_acceptance=ele_acceptance*pim_acceptance;
            /*}}}*/
            total_rate_dvmp += event_weight * total_acceptance;

            /*Missing w/o resolutions{{{*/
            P_E0->SetPxPyPzE(0.,0.,EBeam, EBeam);
            P_t->SetPxPyPzE(0.,0.,0., hMass);
            P_e->SetPxPyPzE(ScatElec_MomX_Col,ScatElec_MomY_Col,ScatElec_MomZ_Col,ScatElec_Energy_Col);	
            P_pim->SetPxPyPzE(Pion_MomX_Col,Pion_MomY_Col,Pion_MomZ_Col,Pion_Energy_Col);	
            P_h->SetPxPyPzE(Proton_MomX_Col,Proton_MomY_Col,Proton_MomZ_Col,Proton_Energy_Col);	

            int err = CheckLaws(P_t, P_E0, P_e, P_pim, P_h);//Check whether momentum and energy conserve first
            if (err < 1e-33){
                cerr<<"---- Momentum and Energy Conservation Laws are broken!! Something is wrong!!!"<<endl;
            }

            double MM = GetMM(P_t, P_E0, P_e, P_pim);	
            hMM_dvmp->Fill(MM, event_weight*total_acceptance);
            /*}}}*/

        ///////////////////////////////////////////////////////////////////////////
        //Now consider the detector resolution here
        //////////////////////////////////////////////////////////////////////////*{{{*/
        //
        //Electron/*{{{*/
        double eP_res = gRandom->Gaus(ScatElec_Mom_Col, Sigma_EC_E*ScatElec_Mom_Col);//GeV, for electron, E ~= P
        double eTheta_res = gRandom->Gaus(ScatElec_Theta_Col*Deg2Rad, Sigma_Theta_E);//rad
        double ePhi_res = gRandom->Gaus(ScatElec_Phi_Col*Deg2Rad, Sigma_Phi_E);//rad

        double ePx_res = eP_res * sin(eTheta_res)*cos(ePhi_res); 
        double ePy_res = eP_res * sin(eTheta_res)*sin(ePhi_res); 
        double ePz_res = eP_res * cos(eTheta_res);
        double eE_res = sqrt(eP_res*eP_res + eMass*eMass);

        P_e_res->SetPxPyPzE(ePx_res, ePy_res, ePz_res, eE_res);	/*}}}*/
        
        //Pion/*{{{*/
        double pimP_res = gRandom->Gaus(Pion_Mom_Col, Sigma_EC_E*Pion_Mom_Col);//GeV, for electron, E ~= P
        double pimTheta_res = gRandom->Gaus(Pion_Theta_Col*Deg2Rad, Sigma_Theta_E);//rad
        double pimPhi_res = gRandom->Gaus(Pion_Phi_Col*Deg2Rad, Sigma_Phi_E);//rad

        double pimPx_res = pimP_res * sin(pimTheta_res)*cos(pimPhi_res); 
        double pimPy_res = pimP_res * sin(pimTheta_res)*sin(pimPhi_res); 
        double pimPz_res = pimP_res * cos(pimTheta_res);
        double pimE_res = sqrt(pimP_res*pimP_res + piMass*piMass);

        P_pim_res->SetPxPyPzE(pimPx_res, pimPy_res, pimPz_res, pimE_res);	/*}}}*/
        

        /////////////////////////////////////////
        double MM_res = GetMM(P_t, P_E0, P_e_res, P_pim_res);	
        hMM_dvmp_res->Fill(MM_res, event_weight*total_acceptance);	
        //////////////////////////////////////////*}}}*/
        
        }

	}// events loop ends here
    cout<<"--- Total DVMP Rate="<<total_rate_dvmp<<endl;
    /*}}}*/

	/*Define SIDIS root file{{{*/
    TString prefix = "./sidis_rootfiles/sidis_3he_";
    TString posfix,new_filename;
    TChain *t2 = new TChain("T","T");
    for (Int_t i=1; i<=2;i++){
        if(energy_flag==11){
            posfix.Form("_11_0_%d_0.root",i);
            new_filename = prefix + particle + posfix;
            cerr<<Form(" @@@ Adding Root File: %s", new_filename.Data())<<endl;
            t2->AddFile(new_filename);
        }
        if(energy_flag==8){
            posfix.Form("_8_0_%d_0.root",i);
            new_filename = prefix + particle + posfix;
            cerr<<Form("     Adding Root File: %s", new_filename.Data())<<endl;
            t2->AddFile(new_filename);
        }
    }
	cerr<<Form("   Got total number of events = %d", (int)(t2->GetEntries()))<<endl;

	Double_t theta_gen, phi_gen,mom_gen;
	//Double_t Q2,W,Wp,x,y,z,pt,nu,s;
	Double_t Q2,Wp,pt,nu,s;
	Double_t theta_q,phi_q;
	Double_t theta_s,phi_s,phi_h;
	Double_t jacoF,dxs_hp,dxs_hm;
	Double_t mom_ele,mom_had, weight_hp, weight_hm;
	Double_t theta_ele,theta_had;
	Double_t phi_ele,phi_had;
	Double_t mom_pro,energy_pro;
	Double_t mom_ini_ele,energy_ini_ele;
	Double_t dilute[2];
	Int_t nsim;

	t2->SetBranchAddress("Q2",&Q2);
	t2->SetBranchAddress("W",&W);
	t2->SetBranchAddress("Wp",&Wp);
	t2->SetBranchAddress("x",&x);
	t2->SetBranchAddress("y",&y);
	t2->SetBranchAddress("z",&z);
	t2->SetBranchAddress("nu",&nu);
	t2->SetBranchAddress("s",&s);
	t2->SetBranchAddress("pt",&pt);
	t2->SetBranchAddress("theta_q",&theta_q);
	t2->SetBranchAddress("theta_s",&theta_s);
	t2->SetBranchAddress("phi_h",&phi_h);
	t2->SetBranchAddress("phi_s",&phi_s);
	t2->SetBranchAddress("jacoF",&jacoF);
	t2->SetBranchAddress("dxs_hm",&dxs_hm);
	t2->SetBranchAddress("dxs_hp",&dxs_hp);
	t2->SetBranchAddress("mom_ele",&mom_ele);
	t2->SetBranchAddress("mom_had",&mom_had);
	t2->SetBranchAddress("theta_ele",&theta_ele);
	t2->SetBranchAddress("theta_had",&theta_had);
	t2->SetBranchAddress("phi_ele",&phi_ele);
	t2->SetBranchAddress("phi_had",&phi_had);
	t2->SetBranchAddress("nsim",&nsim);
	t2->SetBranchAddress("dilute_p",&dilute[0]);
	t2->SetBranchAddress("dilute_m",&dilute[1]);

	double N_sidis=t2->GetEntries();
	t2->GetEntry(N_entries-1);          //get nsim for this rootfile
	double N_simulate=(double)(nsim);
	cout<<"N_simulate: "<<N_simulate<<endl;
	double electron_phase_space=(cos(7/180.*3.1415926) - cos(30/180.*3.1415926))*2*3.14159265*(EBeam-0.5);   // theta: 7~30 degree,  2pi phi coverage, 0.5~11 GeV Momentum coverage 	
	double hadron_phase_space=(cos(7/180.*3.1415926) - cos(30/180.*3.1415926))*2*3.14159265*(6-0.5);  //theta, 7~30 degree,  2pi phi coverage, 0.5~6 GeV Momentum coverage
	double Phase_space=electron_phase_space*hadron_phase_space;           //electron*hadron phase space eg, for electron: delta_cos_theta*delta_phi*delta_energy
	cout<<"Phase_space: "<<electron_phase_space<<"	"<<hadron_phase_space<<"	"<<Phase_space<<endl;
	    
    /*}}}*/

    /*Fill in SIDIS Missing Mass{{{*/
    double total_rate_sidis=0.0; //Make sure the value is right
    for(Long64_t i=0;i<N_sidis;i++){
		t2->GetEntry(i);
        theta_ele *= Rad2Deg;
        phi_ele *= Rad2Deg;
        theta_had *= Rad2Deg;
        phi_had *= Rad2Deg;

        //if(theta_ele>7.5&&theta_ele<24.5&&mom_ele>1.&&mom_ele<11
                //&&theta_had>7.5&&theta_had<24.5&&mom_had>1.&&mom_had<11
          //&&Q2>1.0&&W>2.0){//any additional cuts should be added in here
		if(Q2>=1&& W>=2.){ //zhihong's cut
            /*Get acceptance of e and pi-{{{*/
            double ele_forward_acceptance = accpt->GetAcc("e-","forward", mom_ele, theta_ele);
            double ele_large_acceptance = accpt->GetAcc("e-","large", mom_ele, theta_ele);
			if(mom_ele<1.0||theta_ele>14.5||theta_ele<8.0)//GeV, CLEO
				ele_forward_acceptance=0.0;//Farward-Angle EC Cut at 1 GeV
			if(mom_ele<3.5||theta_ele<16.0||theta_ele>24)//GeV,CLEO
				ele_large_acceptance=0.0; //Larger-Angle EC Cut at 3 GeV
			if(ele_forward_acceptance>1.) 
				ele_forward_acceptance=1.0; 
			if(ele_large_acceptance>1.) 
				ele_large_acceptance=1.0; 

            double pim_forward_acceptance = accpt->GetAcc("pi-","forward", mom_had, theta_had);
            double pim_large_acceptance = accpt->GetAcc("pi-","large", mom_had, theta_had);
            if(theta_had>14.8||theta_had<8.0||mom_had<1.||mom_had>11.)//GeV, CLEO
                pim_forward_acceptance=0.0;
			if(theta_had<16.0||theta_had>24.0||mom_had<3.5||mom_had>11.)//GeV, CLEO
				pim_large_acceptance=0.0; 
			if(pim_forward_acceptance>1.) 
				pim_forward_acceptance=1.0; 
			if(pim_large_acceptance>1.) 
				pim_large_acceptance=1.0; 
            
			weight_hm=dxs_hm*Phase_space/N_simulate;
            double event_weight=weight_hm*nBcm2*Lumi ;   //in Hz

			double ele_acceptance=(ele_forward_acceptance+ele_large_acceptance);
			//double pim_acceptance=(pim_large_acceptance+pim_forward_acceptance);
			double pim_acceptance=pim_forward_acceptance;
			double forward_acceptance=ele_forward_acceptance*pim_acceptance;
			double large_acceptance=ele_large_acceptance*pim_acceptance;
			double total_acceptance=ele_acceptance*pim_acceptance;

            total_rate_sidis += event_weight * total_acceptance;
            /*}}}*/

            /*Missing w/o resolutions{{{*/
            P_E0->SetPxPyPzE(0.,0.,EBeam, EBeam);
            P_t->SetPxPyPzE(0.,0.,0., hMass);
            double mom_ele_x = mom_ele * sin(theta_ele*Deg2Rad) * cos(phi_ele*Deg2Rad);
            double mom_ele_y = mom_ele * sin(theta_ele*Deg2Rad) * sin(phi_ele*Deg2Rad);
            double mom_ele_z = mom_ele * cos(theta_ele*Deg2Rad);
            double energy_ele = sqrt(mom_ele*mom_ele + eMass*eMass);
            P_e->SetPxPyPzE(mom_ele_x,mom_ele_y,mom_ele_z,energy_ele);	

            double mom_had_x = mom_had * sin(theta_had*Deg2Rad) * cos(phi_had*Deg2Rad);
            double mom_had_y = mom_had * sin(theta_had*Deg2Rad) * sin(phi_had*Deg2Rad);
            double mom_had_z = mom_had * cos(theta_had*Deg2Rad);
            double energy_had = sqrt(mom_had*mom_had + piMass*piMass);
            P_pim->SetPxPyPzE(mom_had_x,mom_had_y,mom_had_z,energy_had);	

            double MM = GetMM(P_t, P_E0, P_e, P_pim);	
            hMM_sidis->Fill(MM, event_weight*total_acceptance);
            /*}}}*/

        ///////////////////////////////////////////////////////////////////////////
        //Now consider the detector resolution here
        //////////////////////////////////////////////////////////////////////////*{{{*/
        //
        //Electron/*{{{*/
        double eP_res = gRandom->Gaus(mom_ele, Sigma_EC_E*mom_ele);//GeV, for electron, E ~= P
        double eTheta_res = gRandom->Gaus(theta_ele*Deg2Rad, Sigma_Theta_E);//rad
        double ePhi_res = gRandom->Gaus(phi_ele*Deg2Rad, Sigma_Phi_E);//rad

        double ePx_res = eP_res * sin(eTheta_res)*cos(ePhi_res); 
        double ePy_res = eP_res * sin(eTheta_res)*sin(ePhi_res); 
        double ePz_res = eP_res * cos(eTheta_res);
        double eE_res = sqrt(eP_res*eP_res + eMass*eMass);

        P_e_res->SetPxPyPzE(ePx_res, ePy_res, ePz_res, eE_res);	/*}}}*/
        
        /*Pion{{{*/
        double pimP_res = gRandom->Gaus(mom_had, Sigma_EC_E*mom_had);//GeV, for hadctron, E ~= P
        double pimTheta_res = gRandom->Gaus(theta_had*Deg2Rad, Sigma_Theta_E);//rad
        double pimPhi_res = gRandom->Gaus(phi_had*Deg2Rad, Sigma_Phi_E);//rad

        double pimPx_res = pimP_res * sin(pimTheta_res)*cos(pimPhi_res); 
        double pimPy_res = pimP_res * sin(pimTheta_res)*sin(pimPhi_res); 
        double pimPz_res = pimP_res * cos(pimTheta_res);
        double pimE_res = sqrt(pimP_res*pimP_res + piMass*piMass);

        P_pim_res->SetPxPyPzE(pimPx_res, pimPy_res, pimPz_res, pimE_res);	/*}}}*/


        /////////////////////////////////////////
        double MM_res = GetMM(P_t, P_E0, P_e_res, P_pim_res);	
        hMM_sidis_res->Fill(MM_res, event_weight*total_acceptance);	
        //////////////////////////////////////////*}}}*/
        
        }

	}// events loop ends here
    cout<<"--- Total SIDIS Rate="<<total_rate_sidis<<endl;
    /*}}}*/

    TCanvas *c1 =new TCanvas("c1","c1", 800,600);
    hMM_sidis_res->SetLineColor(4);
    hMM_dvmp_res->Draw();
    hMM_sidis_res->SetLineColor(2);
    hMM_sidis_res->Draw("same");

    TFile *outf = new TFile(Form("MM_n_%d.root",energy_flag), "recreate");
    hMM_sidis->Write();
    hMM_sidis_res->Write();
    hMM_dvmp->Write();
    hMM_dvmp_res->Write();
    outf->Close();
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

		cerr<<Form("*** dE = %f,  dPx = %f, dPy = %f, dPz = %f", energy_check, px_check, py_check, pz_check)<<endl;

		err = -1;
	}
	return err;

}
/*}}}*/
