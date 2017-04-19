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

void plot(){
    int energy=11;
    TString filename=""; 
    filename="noFR_E11_0_0.root";

    TString pro_cut0 = "W>2.&&(pro_theta>24.&&pro_theta<70.)";
   // TString pro_cut0 = "W>2.&&((pro_theta>8.&&pro_theta<14.8) || (pro_theta>24.&&pro_theta<70.))";
    TString pro_cut = "W>2.&&((pro_theta>8.&&pro_theta<14.8) || (pro_theta>16.&&pro_theta<24.))";


    TFile* f1 = new TFile(filename.Data(),"r");
    //TFile* f1 = new TFile(filename.Data(),"update");
    TTree* T = (TTree*) f1->Get("t2");

    gStyle->SetOptStat(0);

    TCanvas *c1 = new TCanvas("c1","c1", 600,900);
    c1->Divide(1,3);

    TCanvas *c2 = new TCanvas("c2","c2", 600,900);
    c2->Divide(1,3);

    TH2F* k1 = new TH2F("k1","Electron Acceptance (e_acc*pi_acc)", 100, 5.0, 28.0, 300, 0.5, 8.5);/*{{{*/
    TH2F* k2 = new TH2F("k2","Pion Acceptance (e_acc*pi_acc)", 100, 5.0, 28.0, 300, 0.5, 8.5);
    TH2F* k3 = new TH2F("k3","Proton Acceptance (e_acc*pi_acc)", 100, 5.0, 55.0, 300, 0., 1.6);
    k1->SetXTitle("#theta_{e} (degree)");
    k1->SetYTitle("P_{e} (GeV/c)");
    k1->GetXaxis()->CenterTitle(1);
    k1->GetYaxis()->CenterTitle(1);
    k2->SetXTitle("#theta_{#pi^{-}} (degree)");
    k2->SetYTitle("P_{#pi^{-}} (GeV/c)");
    k2->GetXaxis()->CenterTitle(1);
    k2->GetYaxis()->CenterTitle(1);
    k3->SetXTitle("#theta_{p} (degree)");
    k3->SetYTitle("P_{p} (GeV/c)");
    k3->GetXaxis()->CenterTitle(1);
    k3->GetYaxis()->CenterTitle(1);
     
    k2->GetXaxis()->SetTitleOffset(0.7);
    k2->GetXaxis()->SetTitleSize(0.08);
    k2->GetYaxis()->SetTitleOffset(0.5);
    k2->GetYaxis()->SetTitleSize(0.08);
 
    k1->GetXaxis()->SetTitleOffset(0.7);
    k1->GetXaxis()->SetTitleSize(0.08);
    k1->GetYaxis()->SetTitleOffset(0.5);
    k1->GetYaxis()->SetTitleSize(0.08);

    k3->GetXaxis()->SetTitleOffset(0.7);
    k3->GetXaxis()->SetTitleSize(0.08);
    k3->GetYaxis()->SetTitleOffset(0.5);
    k3->GetYaxis()->SetTitleSize(0.08);
/*}}}*/

    k1->GetZaxis()->SetRangeUser(0.000, 7000.01);
    k2->GetZaxis()->SetRangeUser(0.000, 7000.01);
    k3->GetZaxis()->SetRangeUser(0.000, 7000.01);

    c2->cd(1);
    T->Draw("ele_mom:ele_theta>>k1",Form("time*weight*(ele_acc_f+ele_acc_l)*(pim_acc_f)*(Q2>4.0&&W>2)*(%s)", pro_cut0.Data()),"colz");
    gPad->SetLogz(1);
    c2->cd(2);
    T->Draw("pim_mom:pim_theta>>k2",Form("time*weight*(ele_acc_f+ele_acc_l)*(pim_acc_f)*(Q2>4.0&&W>2)*(%s)", pro_cut0.Data()),"colz");
    gPad->SetLogz(1);
    c2->cd(3);
    T->Draw("pro_mom:pro_theta>>k3",Form("time*weight*(ele_acc_f+ele_acc_l)*(pim_acc_f)*(Q2>4.0&&W>2)*(%s)", pro_cut0.Data()),"colz");
    gPad->SetLogz(1);
    c2->Modified();
    c2->Update();

    c2->Print(Form("E%d_acc_epi_Q2gt4_noP.png", energy));
    c2->Print(Form("E%d_acc_epi_Q2gt4_noP.pdf", energy));

    TCanvas *c3 = new TCanvas("c3","c3", 600,900);
    c3->Divide(1,3);

    TH2F* j1 = new TH2F("j1","Electron Acceptance (e_acc*pi_acc*p_acc)", 300, 5.0, 28.0, 300, 0.5, 8.5);/*{{{*/
    TH2F* j2 = new TH2F("j2","Pion Acceptance (e_acc*pi_acc*p_acc)", 300, 5.0, 28.0, 300, 0.5, 8.5);
    TH2F* j3 = new TH2F("j3","Proton Acceptance (e_acc*pi_acc*p_acc)", 300, 5.0, 28.0, 300, 0., 1.6);
    j1->SetXTitle("#theta_{e} (degree)");
    j1->SetYTitle("P_{e} (GeV/c)");
    j1->GetXaxis()->CenterTitle(1);
    j1->GetYaxis()->CenterTitle(1);
    j2->SetXTitle("#theta_{#pi^{-}} (degree)");
    j2->SetYTitle("P_{#pi^{-}} (GeV/c)");
    j2->GetXaxis()->CenterTitle(1);
    j2->GetYaxis()->CenterTitle(1);
    j3->SetXTitle("#theta_{p} (degree)");
    j3->SetYTitle("P_{p} (GeV/c)");
    j3->GetXaxis()->CenterTitle(1);
    j3->GetYaxis()->CenterTitle(1);
 
    j2->GetXaxis()->SetTitleOffset(0.7);
    j2->GetXaxis()->SetTitleSize(0.08);
    j2->GetYaxis()->SetTitleOffset(0.5);
    j2->GetYaxis()->SetTitleSize(0.08);
 
    j1->GetXaxis()->SetTitleOffset(0.7);
    j1->GetXaxis()->SetTitleSize(0.08);
    j1->GetYaxis()->SetTitleOffset(0.5);
    j1->GetYaxis()->SetTitleSize(0.08);

    j3->GetXaxis()->SetTitleOffset(0.7);
    j3->GetXaxis()->SetTitleSize(0.08);
    j3->GetYaxis()->SetTitleOffset(0.5);
    j3->GetYaxis()->SetTitleSize(0.08);
/*}}}*/

    j1->GetZaxis()->SetRangeUser(0.000, 800.01);
    j2->GetZaxis()->SetRangeUser(0.000, 800.01);
    j3->GetZaxis()->SetRangeUser(0.000, 800.01);

    c3->cd(1);
    T->Draw("ele_mom:ele_theta>>j1",Form("time*weight*(ele_acc_f+ele_acc_l)*(pim_acc_f)*(%s)*(Q2>4.0&&W>2)", pro_cut.Data()),"colz");
    gPad->SetLogz(1);
    c3->cd(2);
    T->Draw("pim_mom:pim_theta>>j2",Form("time*weight*(ele_acc_f+ele_acc_l)*(pim_acc_f)*(%s)*(Q2>4.0&&W>2)", pro_cut.Data()),"colz");
    gPad->SetLogz(1);
    c3->cd(3);
    T->Draw("pro_mom:pro_theta>>j3",Form("time*weight*(ele_acc_f+ele_acc_l)*(pim_acc_f)*(%s)*(Q2>4.0&&W>2)", pro_cut.Data()),"colz");
    gPad->SetLogz(1);
    c3->Modified();
    c3->Update();

    c3->Print(Form("E%d_acc_epip_Q2gt4.png", energy));
    c3->Print(Form("E%d_acc_epip_Q2gt4.pdf", energy));

    cout<<Form("%f   %f   %f", k1->GetMaximum(), k2->GetMaximum(), k3->GetMaximum() )<<endl;
    cout<<Form("%f   %f   %f", j1->GetMaximum(), j2->GetMaximum(), j3->GetMaximum() )<<endl;

    //f1->cd(); 
    //h1->Write(); h2->Write(); h3->Write();
    //k1->Write(); k2->Write(); k3->Write();
    //j1->Write(); j2->Write(); j3->Write();
    //f1->Close();
}
