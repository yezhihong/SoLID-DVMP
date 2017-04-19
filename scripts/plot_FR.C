{
  gStyle->SetOptStat(0);
  TFile *f1 = new TFile("noFR_E11_0_0.root");
  TTree* T = (TTree*) f1->Get("t2");

  TCanvas *c1 = new TCanvas("c1","c1", 800,600);

  TH2F *h1 = new TH2F("h1","", 90, 0.0, 1.2, 90, 3.5, 9.0);
  h1->SetXTitle("-t (GeV^{2})");
  h1->SetYTitle("Q^{2} (GeV^{2})");
  h1->GetXaxis()->CenterTitle(1);
  h1->GetYaxis()->CenterTitle(1);

  h1->GetZaxis()->SetRangeUser(0, 6500);

  T->Draw("Q2:t>>h1","weight*time*(ele_acc_f+ele_acc_l)*pim_acc_f*(W>2.0&&Epsilon>0.55&&Epsilon<0.75&&(pro_mom_cor<1.0&&(pro_theta_cor>8.&&pro_theta_cor<14.8)||(pro_theta_cor>16&pro_theta_cor<24)))","colz");

 gPad->SetLogz(1);

}
