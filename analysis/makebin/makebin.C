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

    Int_t target_direction=1; 
    cout<<"--- Target Direction (1->Up, 2->Down): "; cin>>target_direction;
    TString targetname="";
    if(target_direction==1) targetname="up";
    if(target_direction==2) targetname="down";

    const Double_t dilute_factor = 0.9;
    const Double_t target_polarization = 0.6; //60% polarization
    //const Double_t beam_polarization = 1.0; //60% polarization
    const Double_t det_eff = 0.85; //85% detector efficiency for electrons and hadrons
    const Double_t target_factor = 0.865; //neutron polarization 86.5%
    //	const Double_t Asys = 0.006; 
    //const Double_t Nsys = 2e+6; 
    /*}}}*/

       
    /*Binning{{{*/
    TString histoname;
    Double_t N_raw= 0.0,N_out = 0.0, Asym= 0.0, Astat=0.0;
    Double_t t_min = 0.0, t_max = 0.0;
    
    TCanvas *c1 = new TCanvas("c1","c1",900,900);
    TCanvas *c3 = new TCanvas("c3","c3",1000,300);

    TString filename = "";
    TString prefix = Form("./database/");
    filename = Form("%s_asym_%s_1.dat",target.Data(),  targetname.Data());
    TString new_filename = prefix + filename;
    ofstream outf1(new_filename);

    filename = Form("%s_asym_%s_2.dat",target.Data(),  targetname.Data());
    new_filename = prefix + filename;
    ofstream outf2(new_filename);

    const Int_t time = 48 * 24 * 3600;
    //const Double_t Norm_Fact = (pow(target_factor * dilute_factor,2) * det_eff);
    const Double_t det_eff_total = pow(sqrt(det_eff), 3);
    const Int_t tbin = 8;
    const Double_t t_cut[9] = {0.0, 0.30, 0.40, 0.50, 0.60, 0.75, 0.95, 1.2, 2.0};

    const TString cut0=Form("(weight_uu*(1.0-0.1*sin((Phi-PhiS)*3.1415926/180.) *%2.1f*%2.1f)*%3.2f)", target_polarization, dilute_factor, det_eff_total);//weight= (Sig_UU+Sig_UT)*Lumi*PSF * Acc, 
    const TString cut1=Form("(weight_uu*(1.0+0.1*sin(PhiS*3.1415926/180.)*%2.1f*%2.1f)*%3.2f)", target_polarization, dilute_factor, det_eff_total);//weight= (Sig_UU+Sig_UT)*Lumi*PSF * Acc, 
    TString cutXS,cutAsym;
    const Int_t PhiBin = 12;
 
    histoname = Form("./database/histo_%s_all.root",  targetname.Data());
    TFile *outroot = new TFile(histoname.Data(),"recreate");

    for (Int_t i=0;i<tbin;i++){
        TString finalfile = Form("../rootfiles/dvmp_%s_t%d_newt.root",  targetname.Data(), i+1);
        TFile *file = new TFile(finalfile.Data(),"r");
        TTree *t0 = (TTree*) file->Get("T");
        t_min = t_cut[i];
        t_max = t_cut[i+1];

        for (Int_t j=1;j<=2;j++){
            TString cutKin = Form("(Epsilon>0.55&&Epsilon<0.75&&W>2&&Qsq>4&&t>%3.2f&&t<%3.2f&&(Q2BIN==%d))", t_min,t_max, j);
            cutXS = cut0 + "*" + cutKin;
            cutAsym= cut1 + "*" + cutKin;
            cerr<<endl<<"--- CutXS   = "<<cutXS.Data()<<endl;
            cerr<<"--- CutAsym = "<<cutAsym.Data()<<endl;
            TCut cut = (TCut) (cutXS);

            /*Histograms{{{*/
            //Q2
            TH1D *h1Q2  =new TH1D(Form("h1Q2"),"h1Q2",500,0.,10.);
            //x 
            TH1D *h1x = new TH1D(Form("h1x"),"h1x",500,0.1,0.7);
            //W 
            TH1D *h1W = new TH1D(Form("h1W"),"h1W",500,0.,10.);
            //Epsilon 
            TH1D *h1Ep = new TH1D(Form("h1Ep"),"h1Ep",500,0.,1.);
            //t 
            TH1D *h1t = new TH1D(Form("h1t"),"h1t",500, 0.0, 2.2);
            //R=Sigma_L/Sigma_T
            TH1D *h1R = new TH1D(Form("h1R"),"h1R",500,0.,10.);
            //F, dilution from L/T separation
            TH1D *h1F = new TH1D(Form("h1F"),"h1F",500,0.,1.);
            //TSA 
            TH1D *h1TSA = new TH1D(Form("h1TSA"),"h1TSA",500,-1.0,0.2);
            //Sigma
            TH1D *h1XS = new TH1D(Form("h1XS"),"h1XS (log10)",500,-12.0,2.);

            //Weight
            TH1D *h1Wt = new TH1D(Form("h1Wt"),"h1Wt",500,0.0,0.01);


            TH1D *h1PhiS = new TH1D(Form("h1PhiS"),"#phi_{S}",PhiBin,0.0,360.);
            TH1D *h1PhiH = new TH1D(Form("h1PhiH"),"#phi",PhiBin,0.0,360.);
            TH2D *h2Phi = new TH2D(Form("h2Phi"),"#phi_{S}:#phi",PhiBin,0.0,360.,PhiBin,0,360);
            TH2D *h2PhiAsym = new TH2D(Form("h2PhiAsym"),"#phi_{S}:#phi",PhiBin,0.0,360.,PhiBin,0,360);
            /*}}}*/

            /*t Binning{{{*/
            c1->Clear();c1->Divide(3,3);
            c1->cd(1); h1Q2->SetLineColor(1);t0->Draw("Qsq>>h1Q2",cut);
            c1->cd(2); h1x->SetLineColor(1); t0->Draw("x>>h1x",cut);
            c1->cd(3); h1W->SetLineColor(1); t0->Draw("W>>h1W",cut);
            c1->cd(4); h1t->SetLineColor(1); t0->Draw("t>>h1t",cut);
            c1->cd(5); h1TSA->SetLineColor(1);t0->Draw("-dilute*SSAsym>>h1TSA",cut);
            c1->cd(6); t0->Draw("log10(Sigma_Lab)>>h1XS","");
            c1->cd(7); t0->Draw("PhiS>>h1PhiS",cut,"");
            c1->cd(8); t0->Draw("Phi>>h1PhiH",cut,"");
            c1->cd(9); t0->Draw("PhiS:Phi>>h2Phi",cut,"colz");
            c1->Print(Form("./figure/plot_%s_%d_%d.png",  targetname.Data(),i+1,j));
            c1->Print(Form("./figure/plot_%s_%d_%d.pdf",  targetname.Data(),i+1,j));

            c3->cd(); h2Phi->Draw("LEGO2"); c3->Print(Form("./figure/PhiS_h_%s_%d_%d.png",targetname.Data(), i+1,j));

            c1->cd(5); h1Ep->SetLineColor(1); t0->Draw("Epsilon>>h1Ep",cut);
            c1->cd(6); h1R->SetLineColor(i); t0->Draw("Sig_L/Sig_T>>h1R",cut);
            c1->cd(7); h1F->SetLineColor(i); t0->Draw("dilute>>h1F",cut);
            c1->cd(8); h1Wt->SetLineColor(1); t0->Draw("weight>>h1Ep",cut);
            c1->cd(9); t0->Draw("PhiS:Phi>>h2PhiAsym",cutAsym,"colz");

            N_raw = h1x->GetSum()/det_eff_total * time;
            //N_raw = h1x->GetSum()/Norm_Fact * time;
            N_out = h1x->GetSum();
            Asym = h1TSA->GetMean();
            if(N_out>10)
                //Astat = 1./sqrt(N_out);
                Astat = 1./target_polarization/dilute_factor/target_factor * sqrt(1-pow(target_polarization*dilute_factor*target_factor*Asym, 2))/sqrt(N_out);
            else
                Astat = -1.0;

            if(j==1)
                outf1<<Form("%4d %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %12.4e %12.4e %12.4e %12.4e %12.4e",/*{{{*/
                        i,
                        h1Q2->GetMean(),
                        h1x->GetMean(),
                        h1W->GetMean(),
                        h1t->GetMean(),
                        h1Ep->GetMean(),
                        h1R->GetMean(),
                        h1F->GetMean(),
                        h1XS->GetMean(),
                        Asym,
                        Astat,
                        N_out,
                        N_raw)
                    <<endl;/*}}}*/
            if(j==2)
                outf2<<Form("%4d %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %12.4e %12.4e %12.4e %12.4e %12.4e",/*{{{*/
                        i,
                        h1Q2->GetMean(),
                        h1x->GetMean(),
                        h1W->GetMean(),
                        h1t->GetMean(),
                        h1Ep->GetMean(),
                        h1R->GetMean(),
                        h1F->GetMean(),
                        h1XS->GetMean(),
                        Asym,
                        Astat,
                        N_out,
                        N_raw)
                    <<endl;/*}}}*/

            /*Get phi_S and phi_H bin{{{*/
            Double_t Ncnt[PhiBin][PhiBin], Nstat[PhiBin][PhiBin], PhiS[PhiBin][PhiBin], PhiH[PhiBin][PhiBin];
            Double_t PhiBinSize = (360-0.0)/PhiBin;
            for (int l=0;l<PhiBin;l++){
                for (int k=0;k<PhiBin;k++){
                    Ncnt[l][k] = h2PhiAsym->GetBinContent(l+1,k+1);//[phi_h][phi_S]
                    Nstat[l][k] = h2Phi->GetBinContent(l+1,k+1);
                    PhiS[l][k] = 0.5*PhiBinSize + k*PhiBinSize;
                    PhiH[l][k] = 0.5*PhiBinSize + l*PhiBinSize;

                    if(j==1) outf1<<Form("%d  %d  %10.4f  %10.4f  %12.4e  %12.4e", l,k, PhiS[l][k], PhiH[l][k], Ncnt[l][k], Nstat[l][k])<<endl;
                    if(j==2) outf2<<Form("%d  %d  %10.4f  %10.4f  %12.4e  %12.4e", l,k, PhiS[l][k], PhiH[l][k], Ncnt[l][k], Nstat[l][k])<<endl;
                }
            }
            /*}}}*/
            /*}}}*/
            cout<<Form("----- Bin#%d: x=%5.4f, Q2=%5.3f, t=%5.4f, N=%f/%d ",i, h1x->GetMean(), h1Q2->GetMean(), h1t->GetMean(), double (N_out/time), int(N_raw))<<endl;

            h1Q2->SetName(Form("h1Q2_t%d_Q%d", i+1, j));
            h1x->SetName(Form("h1x_t%d_Q%d", i+1, j));
            h1t->SetName(Form("h1t_t%d_Q%d", i+1, j));
            h1W->SetName(Form("h1W_t%d_Q%d", i+1, j));
            h1Ep->SetName(Form("h1Ep_t%d_Q%d", i+1, j));
            h1R->SetName(Form("h1R_t%d_Q%d", i+1, j));
            h1F->SetName(Form("h1F_t%d_Q%d", i+1, j));
            h1TSA->SetName(Form("h1TSA_t%d_Q%d", i+1, j));
            h1XS->SetName(Form("h1XS_t%d_Q%d", i+1, j));
            h1PhiH->SetName(Form("h1PhiH_t%d_Q%d", i+1, j));
            h1PhiS->SetName(Form("h1PhiS_t%d_Q%d", i+1, j));
            h2Phi->SetName(Form("h2Phi_t%d_Q%d", i+1, j));
            h2PhiAsym->SetName(Form("h2PhiAsym_t%d_Q%d", i+1, j));
            outroot->cd(); 
            h1Q2->Write();  h1x->Write(); h1W->Write(); h1Ep->Write(); 
            h1t->Write();   h1R->Write(); h1F->Write(); h1TSA->Write(); h1XS->Write(); 
            h2Phi->Write(); h2PhiAsym->Write(); h1PhiS->Write(); h1PhiH->Write();
        }
        file->Close();
    }
    /*}}}*/
    outf1.close();
    outf2.close();
    outroot->Close();

    return 0;
}

Double_t max(Double_t a, Double_t b){/*{{{*/
    if(a>b)
        return a;
    else
        return b;
}/*}}}*/
