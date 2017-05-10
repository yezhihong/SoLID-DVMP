///////////////2//////////////////////
//     DVMP Asymmetries Fitting     //
//     TMinuit Minimization         //
//    ---  Zhihong Ye 04/26/2017    //
//////////////////////////////////////
/*Headers{{{*/

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>
#include <iomanip>
#include <string>
#include <cstdio>

#include "math.h"
#include "stdio.h"
#include <stdlib.h>

#include "TGraph.h"
#include "TGraphErrors.h"
#include "TPaveStats.h"
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2F.h>
#include <TAxis.h>
#include <TLine.h>
#include <TLatex.h>
#include <TLegend.h>
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TObject.h"
#include "TSystem.h"
#include <TVirtualFitter.h>
#include <TMinuit.h>
//#include <TFitterMinuit.h>
#include <TApplication.h>
#include <vector>
using namespace std;

/*}}}*/

//Set to one but keep in mind that the weight_ut should set pol=1 as well.
static const Double_t POL = 0.6*0.865 * 0.9;//He3-Pol, Neutron-Effective-Pol, and Dilution
static const Double_t pi = 3.1415926;
static const Double_t Deg2Rad = pi/180.;
static const Double_t Rad2Deg = 180./pi;

TString type_name = "";
TString bin_name = "";
TString fit_pars = "fit5";
ofstream outf;
Double_t asym_1m1_avg,asym_2m1_avg,asym_3m1_avg,asym_0p1_avg,asym_1p1_avg;
Double_t asym_1m1_fit,asym_2m1_fit,asym_3m1_fit,asym_0p1_fit,asym_1p1_fit;
Double_t asym_1m1_err,asym_2m1_err,asym_3m1_err,asym_0p1_err,asym_1p1_err;
Double_t t_avg, tp_avg, Q2_avg, xb_avg,W_avg, dilute_avg,Asym,Astat;
void LoadData( Int_t IT, Int_t IQ);

/*Main{{{*/
Int_t main()
{ 
    gStyle->SetOptFit(1);  
    gStyle->SetOptStat(0);

    Int_t iType = 0;
    cout<<"--- Which file ? (1->simple, 2->mult, 3->mult_fsi, 4->fermi, 5->mult_nofermi)  "; cin >> iType;
    if(iType==1) type_name = "simple"; 
    if(iType==2) type_name = "mult"; 
    if(iType==3) type_name = "mult_fsi"; 
    if(iType==4) type_name = "fermi"; 
    if(iType==5) type_name = "mult_nofermi"; 

    Int_t bin_type = 0;
    cout<<"--- Which Bining? (1->t, 2->tp, 3->log(tp) )"; cin>> bin_type;
    if(bin_type==1) bin_name ="t";
    if(bin_type==2) bin_name ="tp";
    if(bin_type==3) bin_name ="logtp";

    int BINS = 0;
    if(bin_type==1) {BINS= 7;}
    if(bin_type==2) {BINS= 9;}
    if(bin_type==3) {BINS= 10;}

    ifstream inputf; 
    inputf.open(Form("../asym_extr/results/%s_dvmp_par_%s_%s.dat",bin_name.Data(), type_name.Data(), fit_pars.Data()));
    TString com;
    double temp;
    int t_bin, Q2_bin;
    inputf>>com>>com>>com>>com>>com
          >>com>>com>>com>>com>>com
          >>com>>com>>com>>com>>com
          >>com>>com>>com>>com>>com
          >>com>>com>>com>>com>>com>>com>>com;
    for(int i=1; i<=BINS;i++){
        for(int j=0; j<1;j++){

            cout<<Form("--- working on IT = %d, IQ = %d", i,j)<<endl;
            inputf >> t_bin >> Q2_bin
                >> asym_1m1_avg>> asym_1m1_fit>> asym_1m1_err
                >> asym_0p1_avg>> asym_0p1_fit>> asym_0p1_err
                >> asym_2m1_avg>> asym_2m1_fit>> asym_2m1_err
                >> asym_3m1_avg>> asym_3m1_fit>> asym_3m1_err
                >> asym_1p1_avg>> asym_1p1_fit>> asym_1p1_err
                >> temp >> temp >> Asym >> Astat
                >> t_avg >> tp_avg >> xb_avg >> Q2_avg >> W_avg >> dilute_avg;

            if(i!=t_bin || j!=Q2_bin){
               cout<<Form("Something wrong?!  T=%d/%d,  Q=%d/%d", i, t_bin, j, Q2_bin)<<endl;
            
            }

            outf.open(Form("./database/BIN_%s_dvmp_par_%s_%s_%d_%d.dat",bin_name.Data(), type_name.Data(), fit_pars.Data(),i,j));
            outf<<Form("%4s %4s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s",
                    "#t", "#Q2", 
                    "A1M1_AVG", "A_1M1", "dA_1M1", 
                    "A0P1_AVG", "A_0P1", "dA_0P1", 
                    "A2M1_AVG", "A_2M1", "dA_2M1", 
                    "A3M1_AVG", "A_3M1", "dA_3M1", 
                    "A1P1_AVG", "A_1P1", "dA_1P1",
                    "Asym", "Astat",
                    "t", "tp","xb", "Q2", "W", "Dilute"
                    )
                <<endl;

            outf<<Form("%4d %4d %10.3e %10.3e %10.4e %10.3e %10.3e %10.4e %10.3e %10.3e %10.4e %10.3e %10.3e %10.4e %10.3e %10.3e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e",
                    i, j, 
                    asym_1m1_avg, asym_1m1_fit, asym_1m1_err,
                    asym_0p1_avg, asym_0p1_fit, asym_0p1_err,
                    asym_2m1_avg, asym_2m1_fit, asym_2m1_err,
                    asym_3m1_avg, asym_3m1_fit, asym_3m1_err,
                    asym_1p1_avg, asym_1p1_fit, asym_1p1_err,
                    Asym, Astat,
                    t_avg, tp_avg, xb_avg, Q2_avg, W_avg, dilute_avg 
                    )
                <<endl;

            outf <<Form("%6s %6s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s",
                    "Phi", "PhiS",
                    "N_U_avg", "N_D_avg", "Asym_avg","dAsym_avg",
                    "N_U_fit", "N_D_fit", "Asym_fit","dAsym_fit",
                    "Asym_cal","dAsym_cal",
                    "GOOD_fit", "GOOD_cal"
                    )<<endl;

            /////////////////////////////
            //Load Data from MC Rootfiles to bin on phi and phiS
            /////////////////////////////
            LoadData(i, j);
            outf.close();
        }
    }
    return 0;
}
/*}}}*/

/*LoadData{{{*/
void LoadData( Int_t IT, Int_t IQ){
    /*Define{{{*/
    Double_t phi_h, phi_s, MP_res;
    Double_t asym_1m1, asym_2m1, asym_3m1, asym_0p1, asym_1p1, sigma_uu;
    Double_t weight, weight_uu, weight_ut;
    Double_t weight_1m1,weight_2m1,weight_3m1,weight_0p1,weight_1p1;
    Int_t Q2BIN;

    Double_t weight_fit, weight_ut_fit;
    Double_t weight_1m1_fit,weight_2m1_fit,weight_3m1_fit,weight_0p1_fit,weight_1p1_fit;

    const double PhiH_BIN[13] = {0,30,60,90,120,150,180,210,240,270,300,330,360};
    const double PhiS_BIN[13] = {0,30,60,90,120,150,180,210,240,270,300,330,360};
    double Ncnt_U_avg[12][12], Ncnt_D_avg[12][12], Ncnt_U_fit[12][12], Ncnt_D_fit[12][12];
    double Asym_avg[12][12], Asym_fit[12][12], dAsym_avg[12][12], dAsym_fit[12][12];
    double Asym_cal[12][12], dAsym_cal[12][12];
    for(int i=0;i<12;i++){
        for(int j=0;j<12;j++){
            Ncnt_U_avg[i][j]=0.0;   Ncnt_D_avg[i][j]=0.0;
            Ncnt_U_fit[i][j]=0.0;   Ncnt_D_fit[i][j]=0.0;
            Asym_avg[i][j] = 0.0;   dAsym_avg[i][j] = 0.0;
            Asym_fit[i][j] = 0.0;   dAsym_fit[i][j] = 0.0;
            Asym_cal[i][j] = 0.0;   dAsym_cal[i][j] = 0.0;
        }
    }/*}}}*/

   /*Target Polarization is up{{{*/
   TFile *fu = new TFile(Form("../rootfiles/%s_dvmp_up_t%d_%s.root",bin_name.Data(), IT, type_name.Data()));
   TTree *Tu = (TTree*) gDirectory->Get("T");
   Tu->SetBranchAddress("Phi_cor", &phi_h);/*{{{*/
   Tu->SetBranchAddress("PhiS_cor", &phi_s);//should use the corrected quantity which accounts for fermi-motion etc.
   Tu->SetBranchAddress("Q2BIN", &Q2BIN);

   //Tu->SetBranchAddress("t_cor", &t);
      Tu->SetBranchAddress("MP_res", &MP_res);
   Tu->SetBranchAddress("weight", &weight);
   Tu->SetBranchAddress("weight_uu", &weight_uu);
   Tu->SetBranchAddress("weight_ut", &weight_ut);
   Tu->SetBranchAddress("weight_1m1", &weight_1m1);
   Tu->SetBranchAddress("weight_2m1", &weight_2m1);
   Tu->SetBranchAddress("weight_3m1", &weight_3m1);
   Tu->SetBranchAddress("weight_0p1", &weight_0p1);
   Tu->SetBranchAddress("weight_1p1", &weight_1p1);
  
   Tu->SetBranchAddress("Asym_PhiMinusPhiS",  &asym_1m1);
   Tu->SetBranchAddress("Asym_2PhiMinusPhiS", &asym_2m1);
   Tu->SetBranchAddress("Asym_3PhiMinusPhiS", &asym_3m1);
   Tu->SetBranchAddress("Asym_PhiS",          &asym_0p1);
   Tu->SetBranchAddress("Asym_PhiPlusPhiS",   &asym_1p1);
   Tu->SetBranchAddress("Sigma_UU", &sigma_uu);/*}}}*/
   
   Int_t Nu = Tu->GetEntries();
   for(int i=0;i<Nu; i++){
     Tu->GetEntry(i);
     if(Q2BIN!=IQ && IQ!=0) continue; //choose the right Q2 bin

     //last chance to apply whatever cuts here
     if(MP_res>1.2) continue;
   
     if(isnan(weight_uu) || isinf(weight_uu)) weight_uu=0.0;
     if(weight_uu<1e-12 && weight_uu>1e12) weight_uu=0.0;
  
     phi_h *= Deg2Rad;
     phi_s *= Deg2Rad;
    
     /*From MC DATA{{{*/
     weight_1m1 = weight_uu * (asym_1m1* sin(1.*phi_h - phi_s) );
     weight_0p1 = weight_uu * (asym_0p1* sin(0.*phi_h + phi_s) );
     weight_2m1 = weight_uu * (asym_2m1* sin(2.*phi_h - phi_s) );
     weight_3m1 = weight_uu * (asym_3m1* sin(3.*phi_h - phi_s) );
     weight_1p1 = weight_uu * (asym_1p1* sin(1.*phi_h + phi_s) );
     
     weight_ut = weight_1m1
               + weight_0p1;
     if(fit_pars=="fit5"){
         weight_ut += weight_2m1;
         weight_ut += weight_3m1;
         weight_ut += weight_1p1;
     }

     weight_ut *= POL;
     weight = weight_uu - weight_ut;/*}}}*/
     
     /*From Fit{{{*/
     weight_1m1_fit = weight_uu * (asym_1m1_fit* sin(1.*phi_h - phi_s) );
     weight_0p1_fit = weight_uu * (asym_0p1_fit* sin(0.*phi_h + phi_s) );
     weight_2m1_fit = weight_uu * (asym_2m1_fit* sin(2.*phi_h - phi_s) );
     weight_3m1_fit = weight_uu * (asym_3m1_fit* sin(3.*phi_h - phi_s) );
     weight_1p1_fit = weight_uu * (asym_1p1_fit* sin(1.*phi_h + phi_s) );
     
     weight_ut_fit = weight_1m1_fit
                   + weight_0p1_fit;
     if(fit_pars=="fit5"){
         weight_ut_fit += weight_2m1_fit;
         weight_ut_fit += weight_3m1_fit;
         weight_ut_fit += weight_1p1_fit;
     }

     weight_ut_fit *= POL;
     weight_fit = weight_uu - weight_ut_fit;/*}}}*/

     for(int i=0;i<12;i++){
         for(int j=0;j<12;j++){
             if( ( phi_h>PhiH_BIN[i]*Deg2Rad && phi_h<=PhiH_BIN[i+1]*Deg2Rad) && ( phi_s>PhiS_BIN[j] *Deg2Rad&& phi_s<=PhiS_BIN[j+1]*Deg2Rad)){
                 Ncnt_U_avg[i][j] += weight;
                 Ncnt_U_fit[i][j] += weight_fit;
             }
         }
     }

   }
   fu->Close();
   /*}}}*/

   /*Target Polarization is down{{{*/
   TFile *fd = new TFile(Form("../rootfiles/%s_dvmp_down_t%d_%s.root",bin_name.Data(), IT, type_name.Data()));
   TTree *Td = (TTree*) gDirectory->Get("T");
   Td->SetBranchAddress("Phi_cor", &phi_h);/*{{{*/
   Td->SetBranchAddress("PhiS_cor", &phi_s);//should use the corrected quantity which accounts for fermi-motion etc.
   Td->SetBranchAddress("Q2BIN", &Q2BIN);

   Td->SetBranchAddress("MP_res", &MP_res);
   Td->SetBranchAddress("weight", &weight);
   Td->SetBranchAddress("weight_uu", &weight_uu);
   Td->SetBranchAddress("weight_ut", &weight_ut);
   Td->SetBranchAddress("weight_1m1", &weight_1m1);
   Td->SetBranchAddress("weight_2m1", &weight_2m1);
   Td->SetBranchAddress("weight_3m1", &weight_3m1);
   Td->SetBranchAddress("weight_0p1", &weight_0p1);
   Td->SetBranchAddress("weight_1p1", &weight_1p1);
   Td->SetBranchAddress("Asym_PhiMinusPhiS",  &asym_1m1);
   Td->SetBranchAddress("Asym_2PhiMinusPhiS", &asym_2m1);
   Td->SetBranchAddress("Asym_3PhiMinusPhiS", &asym_3m1);
   Td->SetBranchAddress("Asym_PhiS",          &asym_0p1);
   Td->SetBranchAddress("Asym_PhiPlusPhiS",   &asym_1p1);
   Td->SetBranchAddress("Sigma_UU", &sigma_uu);/*}}}*/
   
   Int_t Nd = Td->GetEntries();
   for(int i=0;i<Nd; i++){
     Td->GetEntry(i);
     if(Q2BIN!=IQ && IQ!=0) continue; //choose the right Q2 bin
     //last chance to apply whatever cuts here
     if(MP_res>1.2) continue;
    
     //Note: In the generator Ahmed switch the sign of the polarization. 
     //In this fit, I fix the absolute polarization values, and rotate the phi_S
     //phi_s += 180.0; 
     phi_h *= Deg2Rad; phi_s *= Deg2Rad;
    
     if(isnan(weight_uu) || isinf(weight_uu)) weight_uu=0.0;
     if(weight_uu<1e-12 && weight_uu>1e12) weight_uu=0.0;
     
     /*From MC DATA{{{*/
     weight_1m1 = weight_uu * (asym_1m1* sin(1.*phi_h - phi_s) );
     weight_0p1 = weight_uu * (asym_0p1* sin(0.*phi_h + phi_s) );
     weight_2m1 = weight_uu * (asym_2m1* sin(2.*phi_h - phi_s) );
     weight_3m1 = weight_uu * (asym_3m1* sin(3.*phi_h - phi_s) );
     weight_1p1 = weight_uu * (asym_1p1* sin(1.*phi_h + phi_s) );
     
     weight_ut = weight_1m1
         + weight_0p1;
     if(fit_pars=="fit5"){
         weight_ut += weight_2m1;
         weight_ut += weight_3m1;
         weight_ut += weight_1p1;
     }

     weight_ut *= -1.*POL;
     weight = weight_uu - weight_ut; //follow the HERMES thesis to put a minus sign here/*}}}*/

     /*From Fit{{{*/
     weight_1m1_fit = weight_uu * (asym_1m1_fit* sin(1.*phi_h - phi_s) );
     weight_0p1_fit = weight_uu * (asym_0p1_fit* sin(0.*phi_h + phi_s) );
     weight_2m1_fit = weight_uu * (asym_2m1_fit* sin(2.*phi_h - phi_s) );
     weight_3m1_fit = weight_uu * (asym_3m1_fit* sin(3.*phi_h - phi_s) );
     weight_1p1_fit = weight_uu * (asym_1p1_fit* sin(1.*phi_h + phi_s) );
     
     weight_ut_fit = weight_1m1_fit
                   + weight_0p1_fit;
     if(fit_pars=="fit5"){
         weight_ut_fit += weight_2m1_fit;
         weight_ut_fit += weight_3m1_fit;
         weight_ut_fit += weight_1p1_fit;
     }

     weight_ut_fit *= POL;
     weight_fit = weight_uu - weight_ut_fit;/*}}}*/

     for(int i=0;i<12;i++){
         for(int j=0;j<12;j++){
             if( ( phi_h>PhiH_BIN[i]*Deg2Rad && phi_h<=PhiH_BIN[i+1]*Deg2Rad  ) && ( phi_s>PhiS_BIN[j] *Deg2Rad&& phi_s<=PhiS_BIN[j+1]*Deg2Rad)  ){
                 Ncnt_D_avg[i][j] += weight;
                 Ncnt_D_fit[i][j] += weight_fit;
             }
         }
     }  
  }
   fd->Close();
   /*}}}*/

   double dNU = 0., dND=0.0;
   double phi_h_temp = 0.0, phi_s_temp=0.0;
   double GOOD_asym_fit[12][12],GOOD_asym_cal[12][12];

   TH2F * h_avg = new TH2F("h_avg","", 12,0,360, 12,0,360);
   TH2F * h_fit = new TH2F("h_fit","", 12,0,360, 12,0,360);
   TH2F * h_cal = new TH2F("h_cal","", 12,0,360, 12,0,360);
   TH2F * g_fit = new TH2F("g_fit","", 144,0,144, 144,-1,1);
   TH2F * g_cal = new TH2F("g_cal","", 144,0,144, 144,-1,1);

   for(int i=0;i<12;i++){
         for(int j=0;j<12;j++){
             //Asymmetries and error from the MC data/*{{{*/
            Asym_avg[i][j] = (Ncnt_U_avg[i][j]-Ncnt_D_avg[i][j]) / (Ncnt_U_avg[i][j]+Ncnt_D_avg[i][j]);
            dNU = sqrt( Ncnt_U_avg[i][j]);
            if (dNU<1e-33) dNU = 1;
            dND = sqrt( Ncnt_D_avg[i][j]);
            if (dND<1e-33) dND = 1;
            dAsym_avg[i][j] = sqrt( pow((2*dND*Ncnt_U_avg[i][j]), 2)+pow((2*dNU*Ncnt_D_avg[i][j]),2))/(Ncnt_U_avg[i][j]+Ncnt_D_avg[i][j]) ;/*}}}*/

             //Asymmetries and error from the fitted values /*{{{*/
            dNU = sqrt( Ncnt_U_fit[i][j]);
            if (dNU<1e-33) dNU = 1;
            dND = sqrt( Ncnt_D_fit[i][j]);
            if (dND<1e-33) dND = 1;
            Asym_fit[i][j] = (Ncnt_U_fit[i][j]-Ncnt_D_fit[i][j]) / (Ncnt_U_fit[i][j]+Ncnt_D_fit[i][j]);
            dAsym_fit[i][j] = sqrt( pow((2*dND*Ncnt_U_fit[i][j]), 2)+pow((2*dNU*Ncnt_D_fit[i][j]),2))/(Ncnt_U_fit[i][j]+Ncnt_D_fit[i][j]) ;/*}}}*/
             
            //Asymmetries and error from page#84 of the HERMES thesis /*{{{*/
            phi_h_temp = 0.5*(PhiH_BIN[i]*Deg2Rad+PhiH_BIN[i+1]*Deg2Rad);
            phi_s_temp = 0.5*(PhiS_BIN[j]*Deg2Rad+PhiS_BIN[j+1]*Deg2Rad);

            Asym_cal[i][j]=0.0;
            Asym_cal[i][j]+= (asym_1m1_fit* sin(1.*phi_h_temp - phi_s_temp) );
            Asym_cal[i][j]+= (asym_0p1_fit* sin(0.*phi_h_temp + phi_s_temp) );
            Asym_cal[i][j]+= (asym_2m1_fit* sin(2.*phi_h_temp - phi_s_temp) );
            Asym_cal[i][j]+= (asym_3m1_fit* sin(3.*phi_h_temp - phi_s_temp) );
            Asym_cal[i][j]+= (asym_1p1_fit* sin(1.*phi_h_temp + phi_s_temp) );

            dAsym_cal[i][j]=0.0;
            dAsym_cal[i][j]+= pow(asym_1m1_err* sin(1.*phi_h_temp - phi_s_temp),2);
            dAsym_cal[i][j]+= pow(asym_0p1_err* sin(0.*phi_h_temp + phi_s_temp),2);
            dAsym_cal[i][j]+= pow(asym_2m1_err* sin(2.*phi_h_temp - phi_s_temp),2);
            dAsym_cal[i][j]+= pow(asym_3m1_err* sin(3.*phi_h_temp - phi_s_temp),2);
            dAsym_cal[i][j]+= pow(asym_1p1_err* sin(1.*phi_h_temp + phi_s_temp),2);
            dAsym_cal[i][j]=sqrt(dAsym_cal[i][j]);/*}}}*/

            //Goodness of the fit from page#84 of the HERMES thesis 
            GOOD_asym_fit[i][j] = (Asym_avg[i][j] - Asym_fit[i][j])/sqrt( pow(dAsym_avg[i][j] ,2)+pow(dAsym_fit[i][j] ,2));
            GOOD_asym_cal[i][j] = (Asym_avg[i][j] - Asym_cal[i][j])/sqrt( pow(dAsym_avg[i][j] ,2)+pow(dAsym_cal[i][j] ,2));

            outf <<Form("%6.3f %6.3f %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e",
                    phi_h_temp*Rad2Deg, phi_s_temp*Rad2Deg,
                    Ncnt_U_avg[i][j], Ncnt_D_avg[i][j], Asym_avg[i][j], dAsym_avg[i][j],
                    Ncnt_U_fit[i][j], Ncnt_D_fit[i][j], Asym_fit[i][j], dAsym_fit[i][j],
                    Asym_cal[i][j], dAsym_cal[i][j],
                    GOOD_asym_fit[i][j], GOOD_asym_cal[i][j]
                    )<<endl;

            h_avg->Fill(phi_h_temp*Rad2Deg, phi_s_temp*Rad2Deg, Asym_avg[i][j] );
            h_fit->Fill(phi_h_temp*Rad2Deg, phi_s_temp*Rad2Deg, Asym_fit[i][j] );
            h_cal->Fill(phi_h_temp*Rad2Deg, phi_s_temp*Rad2Deg, Asym_cal[i][j] );
            g_fit->Fill(i*j, GOOD_asym_fit[i][j] );
            g_cal->Fill(i*j, GOOD_asym_cal[i][j] );
         }
   }

   TFile *histo = new TFile(Form("./database/histo_%s_dvmp_par_%s_%s_%d_%d.root",bin_name.Data(), type_name.Data(), fit_pars.Data(),IT,IQ),"recreate");
   histo->cd();
   h_avg->Write(); h_fit->Write(); h_cal->Write();
   g_fit->Write(); g_cal->Write();
   histo->Close();

   delete h_avg;
   delete h_fit;
   delete h_cal;
   delete g_cal;
   delete g_fit;
}
/*}}}*/
