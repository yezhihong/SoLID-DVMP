{
    gStyle->SetOptStat(0);

  TCanvas *c1 = new TCanvas("c1","c1", 800,800);
  c1->Divide(2,2);

  TFile *f1 =new TFile("DEMP_Ee_11_4_up_simple.root");
 
  TH2D *h1u = new TH2D("h1u","Target Up: t vs #phi_{h}",12,0.,360,50, 0.0,1.5);
  h1u->SetXTitle("#phi_{h} (Degree)");
  h1u->SetYTitle("-t (GeV^{2})");
  h1u->GetXaxis()->CenterTitle(1);
  h1u->GetYaxis()->CenterTitle(1);

  TH2D *h2u = new TH2D("h2u","Target Up: t vs #phi_{S}",12,0.,360,50, 0.0,1.5);
  h2u->SetXTitle("#phi_{S} (Degree)");
  h2u->SetYTitle("-t (GeV^{2})");
  h2u->GetXaxis()->CenterTitle(1);
  h2u->GetYaxis()->CenterTitle(1);

  c1->cd(1); T->Draw("t:Phi>>h1u", "weight*(Epsilon>0.55&&Epsilon<0.75&&W>2&&Qsq>4)", "colz");
  c1->cd(2); T->Draw("t:PhiS>>h2u", "weight*(Epsilon>0.55&&Epsilon<0.75&&W>2&&Qsq>4)", "colz");

  TFile *f2 =new TFile("DEMP_Ee_11_4_down_simple.root");

  TH2D *h1d = new TH2D("h1d","Target Down: t vs #phi_{h}",12,0.,360,50, 0.0,1.5);
  h1d->SetXTitle("#phi_{h} (Degree)");
  h1d->SetYTitle("-t (GeV^{2})");
  h1d->GetXaxis()->CenterTitle(1);
  h1d->GetYaxis()->CenterTitle(1);

  TH2D *h2d = new TH2D("h2d","Target Down: t vs #phi_{S}",12,0.,360,50, 0.0,1.5);
  h2d->SetXTitle("#phi_{S} (Degree)");
  h2d->SetYTitle("-t (GeV^{2})");
  h2d->GetXaxis()->CenterTitle(1);
  h2d->GetYaxis()->CenterTitle(1);


  c1->cd(3); T->Draw("t:Phi>>h1d", "weight*(Epsilon>0.55&&Epsilon<0.75&&W>2&&Qsq>4)", "colz");
  c1->cd(4); T->Draw("t:PhiS>>h2d", "weight*(Epsilon>0.55&&Epsilon<0.75&&W>2&&Qsq>4)", "colz");
  c1->Print("phi_t.pdf");
  c1->Print("phi_t.png");

  TFile* histo = new TFile("histo.root","recreate");
  histo->cd();
  h1u->Write(); h2u->Write(); h1d->Write(); h2d->Write();
  histo->Close();

}
