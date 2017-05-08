void plot_mm(){
    gStyle->SetOptStat(0);

    TFile *f1 = new TFile("SIDIS_MM_p_11_Q2gt4.root");
    TFile *f2 = new TFile("SIDIS_MM_n_11_Q2gt4.root");
    TFile *f3 = new TFile("./histo_up_mult.root");

    TH1F *hMM_n_dvmp = (TH1F*) f3->Get("hMM_dvmp");
    TH1F *hMM_n_dvmp_cut = (TH1F*) f3->Get("hMM_dvmp_cut");
    TH1F *hMM_n_dvmp_cor = (TH1F*) f3->Get("hMM_dvmp_cor");
    TH1F *hMM_n_dvmp_cut_cor = (TH1F*) f3->Get("hMM_dvmp_cut_cor");
    TH1F *hMM_n_dvmp_res = (TH1F*) f3->Get("hMM_dvmp_res");
    TH1F *hMM_n_dvmp_cut_res = (TH1F*) f3->Get("hMM_dvmp_cut_res");

    TH1F *hMP_n_dvmp = (TH1F*) f3->Get("hMP_dvmp");
    TH1F *hMT_n_dvmp = (TH1F*) f3->Get("hMT_dvmp");
    TH1F *hMP_n_dvmp_cor = (TH1F*) f3->Get("hMP_dvmp_cor");
    TH1F *hMT_n_dvmp_cor = (TH1F*) f3->Get("hMT_dvmp_cor");
    TH1F *hMP_n_dvmp_res = (TH1F*) f3->Get("hMP_dvmp_res");
    TH1F *hMT_n_dvmp_res = (TH1F*) f3->Get("hMT_dvmp_res");
    
    TH1F *hMM_n_sidis = (TH1F*) f1->Get("hMM_sidis");
    TH1F *hMM_n_sidis_res = (TH1F*) f1->Get("hMM_sidis_res");
    TH1F *hMM_n_sidis_cut = (TH1F*) f1->Get("hMM_sidis_cut");
    TH1F *hMM_n_sidis_cut_res = (TH1F*) f1->Get("hMM_sidis_cut_res");

    TH1F *hMP_n_sidis = (TH1F*) f1->Get("hMP_sidis");
    TH1F *hMP_n_sidis_res = (TH1F*) f1->Get("hMP_sidis_res");
    TH1F *hMT_n_sidis = (TH1F*) f1->Get("hMT_sidis");
    TH1F *hMT_n_sidis_res = (TH1F*) f1->Get("hMT_sidis_res");

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

    TCanvas *c1 = new TCanvas("c1","c1", 800,600);
    hMP_3he_sidis_res->SetLineColor(2); hMP_3he_sidis_res->Draw();
    hMP_n_dvmp_res->SetLineColor(1); hMP_n_dvmp_res->Draw("same");
    hMP_n_dvmp_cor->SetLineStyle(9); hMP_n_dvmp_cor->SetLineColor(4); hMP_n_dvmp_cor->Draw("same");

    TCanvas *c2 = new TCanvas("c2","c2",800,600);
    hMM_3he_sidis_res->SetLineColor(2); hMM_3he_sidis_res->Draw();
    hMM_n_dvmp_res->SetLineColor(1); hMM_n_dvmp_res->Draw("same");
    hMM_n_dvmp_cor->SetLineStyle(3); hMM_n_dvmp_cor->SetLineColor(4); hMM_n_dvmp_cor->Draw("same");

    TCanvas *c3 = new TCanvas("c3","c3",800,600);
    hMM_n_dvmp_cut_res->SetLineColor(1); hMM_n_dvmp_cut_res->Draw("same");
    hMM_n_dvmp_cut_cor->SetLineStyle(9); hMM_n_dvmp_cut_cor->SetLineColor(4); hMM_n_dvmp_cut_cor->Draw("same");
    hMM_3he_sidis_cut_res->SetLineColor(2); hMM_3he_sidis_cut_res->Draw("same");
}
