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

#include "RooRealVar.h"/*{{{*/
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooFormulaVar.h"
#include "RooGenericPdf.h"
#include "RooPolynomial.h"
#include "RooChi2Var.h"
#include "RooMinimizer.h"
#include "RooPlot.h"
#include "RooFitResult.h"
using namespace RooFit ;
using namespace std;

/*}}}*/

//void asym_extr()
int main()
{
   Int_t I = 0; cout<<"--- Which -t bin? "; cin >> I;
   TFile *f1 = new TFile(Form("../rootfiles/dvmp_t%d.root", I));
   TTree *T = (TTree*) gDirectory->Get("T");

   Double_t tphis, tphih, tssasym, tdilute, tweight;
   T->SetBranchAddress("phi_S", &tphis);
   T->SetBranchAddress("phi_h", &tphih);
   T->SetBranchAddress("SSAsym", &tssasym);
   T->SetBranchAddress("dilute", &tdilute);
   T->SetBranchAddress("weight", &tweight);
   Int_t Nevt = T->GetEntries();

   // create observables and read data from ROOT file
   RooRealVar phi_S("phi_S","phi_S", 0.0, 360) ;
   RooRealVar phi_h("phi_h","phi_h", 0.0, 360) ;
   //RooRealVar SSAsym("SSAsym","SSAsym", 0.0, 0.6) ;
   //RooRealVar dilute("dilute","dilute", 0.,1.) ;
   //RooRealVar weight("weight","weight",0., 1000.);
   RooRealVar beta1("beta1","beta1",0., 1000.);
   RooRealVar beta2("beta2","beta2",0., 1000.);
   RooRealVar w("w","w",0., 1000.);
  
   //read data from the tree
   RooDataSet data("data", "data",RooArgSet(phi_S, phi_h, w), w.GetName());
   for(int i=0;i<Nevt; i++){
     T->GetEntry(i);
     //double tw = tweight*(1.0-tdilute*tssasym*sin(3.0*tphis*3.1415926/180+tphih*3.1415926/180));
     double tw = tweight*(1.0-0.3*sin(3.0*tphis*3.1415926/180+tphih*3.1415926/180));
     beta1= sin(3.0*tphis*3.1415926/180+tphih*3.1415926/180);
     beta2= sin(1.0*tphis*3.1415926/180-tphih*3.1415926/180);
     phi_S = tphis;
     phi_h = tphih;
     w = tw;
     data.add( RooArgSet(phi_S, phi_h, w), tw,0);
     
     //cout<<Form("--- phi_s=%f, phi_h=%f,  w=%f/%f", tphis, tphih, tweight, tw)<<endl;
   }

   w.Print();

   // Dataset d is now a dataset with two observable (x,w) 
   data.Print() ;
   //data->setWeightVar(w);
   //RooDataSet wdata(data->GetName(), data->etTitle(), data, data->get(), 0, "w");


   // U n b i n n e d   M L   f i t   t o   w e i g h t e d   d a t a 
   // ---------------------------------------------------------------
   // Construction asymmetry p.d.f. for fitting
   RooRealVar a0("a0","a0",1, 0.,2.) ;
   RooRealVar a1("a1","a1",-0.3,-1,0.) ;
   RooGenericPdf AsymPDF("AsymPDF", " a0*(1.+a1*sin(3*phi_S*3.1415926/180.+phi_h*3.1415926/180.))", RooArgSet(phi_S, phi_h, a0, a1));

   // Fit asymmetry PDF to weighted data
   // NOTE: A plain Maximum likelihood fit to weighted data does in general 
   //       NOT result in correct error estimates, unless individual
   //       event weights represent Poisson statistics themselves.
   //       
   // Fit with 'wrong' errors
   //RooFitResult* r_ml_wgt = AsymPDF.fitTo(data,Save()) ;
   
   // A first order correction to estimated parameter errors in an 
   // (unbinned) ML fit can be obtained by calculating the
   // covariance matrix as
   //
   //    V' = V C-1 V
   //
   // where V is the covariance matrix calculated from a fit
   // to -logL = - sum [ w_i log f(x_i) ] and C is the covariance
   // matrix calculated from -logL' = -sum [ w_i^2 log f(x_i) ] 
   // (i.e. the weights are applied squared)
   //
   // A fit in this mode can be performed as follows:
   RooFitResult* r_ml_wgt_corr = AsymPDF.fitTo(data,Save(),SumW2Error(kTRUE)) ;

   // P l o t   w e i g h e d   d a t a   a n d   f i t   r e s u l t 
   // ---------------------------------------------------------------
   // Construct plot frame
   RooPlot* frame_S = phi_S.frame(Title("Unbinned ML fit, binned chi^2 fit to weighted data")) ;
   RooPlot* frame_h = phi_h.frame(Title("Unbinned ML fit, binned chi^2 fit to weighted data")) ;

   // Plot data using sum-of-weights-squared error rather than Poisson errors
   data.plotOn(frame_S,DataError(RooAbsData::SumW2)) ;
   // Overlay result of the fit to weighted data
   AsymPDF.plotOn(frame_S) ;
   data.plotOn(frame_S,DataError(RooAbsData::SumW2)) ;

   // Plot data using sum-of-weights-squared error rather than Poisson errors
   data.plotOn(frame_h,DataError(RooAbsData::SumW2)) ;
   // Overlay result of the fit to weighted data
   AsymPDF.plotOn(frame_h) ;
   data.plotOn(frame_h,DataError(RooAbsData::SumW2)) ;

   TCanvas* c1 = new TCanvas("c1","c1", 1200,400);
   c1->Divide(2);
   c1->cd(1); frame_S->Draw();
   c1->cd(2); frame_h->Draw();

   c1->Print(Form("fit_%d.png",I));

   //// ML Fit of pdf to equivalent unweighted dataset
   //// -----------------------------------------------------------------------------------------
   //// Construct a pdf with the same shape as p0 after weighting
   //RooGenericPdf AsymPDF_copy("AsymPDF", " a0*(1.+a1*sin(3*phi_S*3.1415926/180.+phi_h*3.1415926/180.))", RooArgSet(phi_S, phi_h, a0, a1));

   //// Sample a dataset with the same number of events as data
   //RooDataSet* data2 = AsymPD_copy.generate(RooArgSet(phi_S, phi_h),1000) ;

   //// Sample a dataset with the same number of weights as data
   //RooDataSet* data3 = AsymPD_copy.generate(RooArgSet(phi_S, phi_h),5000) ;

   //// Fit the 2nd order polynomial to both unweighted datasets and save the results for comparison
   //RooFitResult* r_ml_unw10 = AsymPDF.fitTo(*data2,Save()) ;
   //RooFitResult* r_ml_unw43 = AsymPDF.fitTo(*data3,Save()) ;

f1->Close(); 
}
