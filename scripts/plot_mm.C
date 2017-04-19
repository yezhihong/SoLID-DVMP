//void plot_mm(){
{

    gStyle->SetOptStat(0);

    //TFile *f1 = new TFile("MM_p_11_Q2gt4_fermi_11m.root");
    TFile *f1 = new TFile("MM_p_11_Q2gt4_fermi_mplt1.0.root");
    //TFile *f2 = new TFile("MM_3he_11_Q2gt4_fermi_11m.root");
    TFile *f2 = new TFile("MM_3he_11_Q2gt4_fermi_mplt1.0.root");

    //TFile *f1 = new TFile("MM_p_11_Q2gt4_prd_fermi.root");
    //TFile *f2 = new TFile("MM_3he_11_Q2gt4_prd_fermi.root");
 
    TH1F *hMM_n_dvmp = (TH1F*) f1->Get("hMM_dvmp");
    TH1F *hMM_n_dvmp_res = (TH1F*) f1->Get("hMM_dvmp_res");
    TH1F *hMM_n_dvmp_cut = (TH1F*) f1->Get("hMM_dvmp_cut");
    TH1F *hMM_n_dvmp_cut_res = (TH1F*) f1->Get("hMM_dvmp_cut_res");

    TH1F *hMP_n_dvmp = (TH1F*) f1->Get("hMP_dvmp");
    TH1F *hMP_n_dvmp_res = (TH1F*) f1->Get("hMP_dvmp_res");
    TH1F *hMT_n_dvmp = (TH1F*) f1->Get("hMT_dvmp");
    TH1F *hMT_n_dvmp_res = (TH1F*) f1->Get("hMT_dvmp_res");
 
    double scale = 0.23/hMP_n_dvmp->GetSum();
   
    TH1F *hMM_n_sidis = (TH1F*) f1->Get("hMM_sidis");
    TH1F *hMM_n_sidis_res = (TH1F*) f1->Get("hMM_sidis_res");
    TH1F *hMM_n_sidis_cut = (TH1F*) f1->Get("hMM_sidis_cut");
    TH1F *hMM_n_sidis_cut_res = (TH1F*) f1->Get("hMM_sidis_cut_res");
 
    TH1F *hMP_n_sidis = (TH1F*) f1->Get("hMP_sidis");
    TH1F *hMP_n_sidis_res = (TH1F*) f1->Get("hMP_sidis_res");
    TH1F *hMT_n_sidis = (TH1F*) f1->Get("hMT_sidis");
    TH1F *hMT_n_sidis_res = (TH1F*) f1->Get("hMT_sidis_res");

    TH1F *hMM_p_dvmp = (TH1F*) f2->Get("hMM_dvmp");
    TH1F *hMM_p_dvmp_res = (TH1F*) f2->Get("hMM_dvmp_res");
    TH1F *hMM_p_dvmp_cut = (TH1F*) f2->Get("hMM_dvmp_cut");
    TH1F *hMM_p_dvmp_cut_res = (TH1F*) f2->Get("hMM_dvmp_cut_res");
   
    TH1F *hMP_p_dvmp = (TH1F*) f2->Get("hMP_dvmp");
    TH1F *hMP_p_dvmp_res = (TH1F*) f2->Get("hMP_dvmp_res");
    TH1F *hMT_p_dvmp = (TH1F*) f2->Get("hMT_dvmp");
    TH1F *hMT_p_dvmp_res = (TH1F*) f2->Get("hMT_dvmp_res");

    TH1F *hMM_p_sidis = (TH1F*) f2->Get("hMM_sidis");
    TH1F *hMM_p_sidis_res = (TH1F*) f2->Get("hMM_sidis_res");
    TH1F *hMM_p_sidis_cut = (TH1F*) f2->Get("hMM_sidis_cut");
    TH1F *hMM_p_sidis_cut_res = (TH1F*) f2->Get("hMM_sidis_cut_res");

    TH1F *hMP_p_sidis = (TH1F*) f2->Get("hMP_sidis");
    TH1F *hMP_p_sidis_res = (TH1F*) f2->Get("hMP_sidis_res");
    TH1F *hMT_p_sidis = (TH1F*) f2->Get("hMT_sidis");
    TH1F *hMT_p_sidis_res = (TH1F*) f2->Get("hMT_sidis_res");

    TH1F *hMP_3he_sidis = (TH1F*) hMP_n_sidis->Clone();
    hMP_3he_sidis->Add(hMP_p_sidis, 2.0);//two protons in he3
    TH1F *hMP_3he_sidis_res = (TH1F*) hMP_n_sidis_res->Clone();
    hMP_3he_sidis_res->Add(hMP_p_sidis_res, 2.0);//two protons in he3

    TH1F *hMM_3he_sidis = (TH1F*) hMM_n_sidis->Clone();
    hMM_3he_sidis->Add(hMM_p_sidis, 2.0);//two protons in he3
    TH1F *hMM_3he_sidis_res = (TH1F*) hMM_n_sidis_res->Clone();
    hMM_3he_sidis_res->Add(hMM_p_sidis_res, 2.0);//two protons in he3

    TH1F *hMM_3he_sidis_cut = (TH1F*) hMM_n_sidis_cut->Clone();
    hMM_3he_sidis_cut->Add(hMM_p_sidis_cut, 2.0);//two protons in he3
    TH1F *hMM_3he_sidis_cut_res = (TH1F*) hMM_n_sidis_cut_res->Clone();
    hMM_3he_sidis_cut_res->Add(hMM_p_sidis_cut_res, 2.0);//two protons in he3

    cerr<<"--- Residual total SIDIS rate: "<<hMM_3he_sidis_cut_res->GetSum()<<endl;


    hMP_n_dvmp->Scale(scale);
    hMP_n_dvmp_res->Scale(scale);
    hMM_n_dvmp->Scale(scale);
    hMM_n_dvmp_res->Scale(scale);
    hMM_n_dvmp_cut->Scale(scale);
    hMM_n_dvmp_cut_res->Scale(scale);

    TCanvas *c1 = new TCanvas("c1","c1", 800,600);
    hMP_3he_sidis->SetLineColor(3); hMP_3he_sidis->Draw();
    hMP_3he_sidis_res->SetLineColor(2); hMP_3he_sidis_res->Draw("same");
    hMP_n_dvmp->SetLineColor(6); hMP_n_dvmp->Draw("same");
    hMP_n_dvmp_res->SetLineColor(4); hMP_n_dvmp_res->Draw("same");

    TCanvas *c2 = new TCanvas("c2","c2",1200,600);
    c2->Divide(2,1);
    c2->cd(1);

    hMM_3he_sidis->SetLineColor(3); hMM_3he_sidis->Draw();
    hMM_3he_sidis_res->SetLineColor(2); hMM_3he_sidis_res->Draw("same");
    hMM_n_dvmp->SetLineColor(6); hMM_n_dvmp->Draw("same");
    hMM_n_dvmp_res->SetLineColor(4); hMM_n_dvmp_res->Draw("same");
    
    c2->cd(2);
    hMM_n_dvmp_cut_res->SetLineColor(4); hMM_n_dvmp_cut_res->Draw();
    hMM_n_dvmp_cut->SetLineColor(6); hMM_n_dvmp_cut->Draw("same");
    
   // hMM_n_dvmp_cut->SetLineColor(6); hMM_n_dvmp_cut->Draw();
   // hMM_n_dvmp_cut_res->SetLineColor(4); hMM_n_dvmp_cut_res->Draw("same");
    hMM_3he_sidis_cut_res->SetLineColor(2); hMM_3he_sidis_cut_res->Draw("same");
    hMM_3he_sidis_cut->SetLineColor(3); hMM_3he_sidis_cut->Draw("same");
    

}
