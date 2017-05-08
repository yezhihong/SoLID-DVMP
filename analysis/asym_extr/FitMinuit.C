///////////////2//////////////////////
//     DVMP Asymmetries Fitting     //
//     TMinuit Minimization         //
//    ---  Zhihong Ye 04/26/2017    //
//////////////////////////////////////
#include "FitMinuit.h"

//Set to one but keep in mind that the weight_ut should set pol=1 as well.
static const Double_t POL = 0.6*0.865 * 0.9;//He3-Pol, Neutron-Effective-Pol, and Dilution

TString type_name = "";
TString bin_name = "";
TString fit_pars = "fit5";
ofstream outf;
Double_t asym_1m1_avg,asym_2m1_avg,asym_3m1_avg,asym_0p1_avg,asym_1p1_avg,Asym,Astat;
Double_t t_avg, tp_avg, Q2_avg, xb_avg,W_avg, dilute_avg;

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

    Int_t IT = 0, IQ = 0;

    Int_t bin_type = 0;
    cout<<"--- Which Bining? (1->t, 2->tp)"; cin>> bin_type;
    if(bin_type==1) bin_name ="t";
    if(bin_type==2) bin_name ="tp";

    int BINS = 0;
    if(bin_type==1){ cout<<"--- Which t-bin? (1--8, 0 for all)  "; cin >> IT;  BINS= 8;}
    if(bin_type==2){ cout<<"--- Which tp-bin? (1--11, 0 for all)  "; cin >> IT;  BINS= 11;}

    if(IT==0){
        outf.open(Form("%s_dvmp_par_%s_%s.dat",bin_name.Data(), type_name.Data(), fit_pars.Data()));
        outf<<Form("%4s %4s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s",
                "#t", "#Q2", 
                "A1M1_AVG", "A_1M1", "dA_1M1", 
                "A0P1_AVG", "A_0P1", "dA_0P1", 
                "A2M1_AVG", "A_2M1", "dA_2M1", 
                "A3M1_AVG", "A_3M1", "dA_3M1", 
                "A1P1_AVG", "A_1P1", "dA_1P1",
                "Norm", "dNorm", "Asym", "Astat",
                "t", "tp","xb", "Q2", "W", "Dilute"
                )
            <<endl;


        for(int i=1; i<=BINS;i++){
            //for(int j=0; j<=2;j++){
            for(int j=0; j<1;j++){
                IT = i;
                IQ = j;

                asym_1m1_avg=0.0;
                asym_2m1_avg=0.0;
                asym_3m1_avg=0.0;
                asym_0p1_avg=0.0;
                asym_1p1_avg=0.0;
                Asym=0.0;
                Astat=0.0;
                t_avg=0.0; 
                tp_avg=0.0; 
                Q2_avg=0.0; 
                xb_avg=0.0;
                W_avg=0.0; 
                dilute_avg=0.0;

                /////////////////////////////
                //Load Data from MC Rootfiles
                /////////////////////////////
                LoadData(IT, IQ);

                /////////////////////////////
                //TMinuit Minimization
                /////////////////////////////
                double fyPar[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
                DoMinuit(fyPar, IT, IQ, 1);
            }
        }
        outf.close();
    }else{
        cout<<"--- Which Q2-bin? (1--2, 0 for no Q2 binning)  "; cin >> IQ;
        LoadData(IT, IQ);

        /////////////////////////////
        //TMinuit Minimization
        /////////////////////////////
        double fyPar[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        DoMinuit(fyPar, IT, IQ, 0);
    }

        return 0;
    }
/*}}}*/

/*DoMunuit{{{*/
void DoMinuit(double *par, Int_t IT, Int_t IQ, Int_t SAVE)
{
	cout << "      *************************************************" << endl; 
	cout << "      *          Minimization Fitting for             *" << endl;
	cout << "      *              DVMP Asymmetries                 *" << endl;
	cout << "      *             Z. Ye 04/26/2017                  *" << endl;
	cout << "      *************************************************" << endl; 
	cout << endl;
	gSystem->Load("libMinuit");
	static const Int_t iParaNum=6;

	TMinuit *pMinuit = new TMinuit(iParaNum);  //initialize TMinuit with a maximum of iParaNum params
	pMinuit->SetFCN(myUML);

	Double_t arglist[10];
	Int_t ierflg = 0;

	/// set print level
	/*  SET PRIntout  <level>
		Sets the print level, determining how much output will be
		produced. Allowed values and their meanings are displayed
		after a SHOw PRInt command, and are currently <level>=:
		[-1]  no output except from SHOW commands
		[0]  minimum output
		[1]  default value, normal output
		[2]  additional output giving intermediate results.
		[3]  maximum output, showing progress of minimizations. */
	arglist[0] = 1;
	pMinuit->mnexcm("SET PRIntout",arglist,1,ierflg);

	/*
	   SET NOWarnings
	   Supresses Minuit warning messages.
	   SET WARnings
	   Instructs Minuit to output warning messages when suspicious
	   conditions arise which may indicate unreliable results.
	   This is the default.
	   */
	arglist[0] = 1;
	pMinuit->mnexcm("SET NOWarnings",arglist,0,ierflg);

	//Set Name of Parameters
	TString namePar[iParaNum]={"norm", "a1m1","a0p1","a2m1","a3m1","a1p1"};
	
        //Set Initial Values of Parameters
        Double_t inPar[iParaNum]={par[0],par[1],par[2],par[3],par[4],par[5]};   

	//Set Stepsize of Parameters
	Double_t step[iParaNum]={1.0e-5,1.0e-5,1.0e-5,1.0e-5,1.0e-5,1.0e-5};

	//Set Min of Parameters, value==0 for No-Bound
	Double_t minVal[iParaNum]={0,0,0,0,0,0};

	//Set Max of Parameters, value==0 for No-Bound
	Double_t maxVal[iParaNum]={0,0,0,0,0,0};

	//Initializing Parameters
	for(int ii=0;ii<iParaNum;ii++){
            pMinuit->DefineParameter(ii,namePar[ii], inPar[ii], step[ii], minVal[ii], maxVal[ii]);
	}

	//Fix parameters, starting from 1 not 0
	arglist[0] = 1;   pMinuit->mnexcm("FIX ", arglist ,1,ierflg);
	//arglist[0] = 2;   pMinuit->mnexcm("FIX ", arglist ,1,ierflg);
	//arglist[0] = 3;   pMinuit->mnexcm("FIX ", arglist ,1,ierflg);
    if(fit_pars=="fit2"){
        arglist[0] = 4;   pMinuit->mnexcm("FIX ", arglist ,1,ierflg);
        arglist[0] = 5;   pMinuit->mnexcm("FIX ", arglist ,1,ierflg);
        arglist[0] = 6;   pMinuit->mnexcm("FIX ", arglist ,1,ierflg);
    }

	arglist[0] = 1;
	pMinuit->mnexcm("SET ERR", arglist,1,ierflg);

	/*
	//SET LIMits  0->[parno]  1->[lolim]  2->[uplim]
	//	  arglist[0] = 0;	   //this means set all parameters with the same limit
	//	  arglist[0] = 1;	   //this means set 1st parameters with the specified limit
	//	  arglist[0] = 2;	   //this means set 2nd parameters with the sepecified limit
	arglist[0] = 0;
	arglist[1] = 0.;
	arglist[2] = 0.;
	pMinuit->mnexcm("SET LIMits",arglist,3,ierflg);
	*/

	//Set Iteration Times and Use Migrad to start minimization
	arglist[0] = 500;
	//arglist[1] = 0.001;
	pMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);

	//Get Results
	for(int ii = 0; ii < iParaNum; ii++)
	{
		pMinuit->GetParameter(ii,outPar[ii],err[ii]);    
		par[ii] = outPar[ii];
	}
	//Put the results into a file
        cout<<endl<<"-----------------------------------------------------------------------"<<endl;
        cout<<"-----------------------------------------------------------------------"<<endl;

        cout<<Form("Real:         A1=%6.3f,  A2=%6.3f,  A3=%6.3f,  A4=%6.3f   A5=%6.3f, Asym=%10.4e, Astat=%10.4e",
                asym_1m1_avg,asym_0p1_avg,asym_2m1_avg,asym_3m1_avg,asym_1p1_avg, Asym, Astat)<<endl;

        cout<<Form("Asymmetries:  A1=%6.3f,  A2=%6.3f,  A3=%6.3f,  A4=%6.3f,  A5=%6.3f, Norm=%6.3f", 
                par[1], par[2], par[3], par[4], par[5],par[0])<<endl;
        cout<<Form("     Errors: dA1=%6.3f, dA2=%6.3f, dA3=%6.3f, dA4=%6.3f, dA5=%6.3f, dNorm=%6.3f", 
        err[1], err[2], err[3], err[4], err[5], err[0])<<endl;
        cout<<"-----------------------------------------------------------------------"<<endl;
        cout<<"-----------------------------------------------------------------------"<<endl<<endl;
        
        if(SAVE==1){
            outf<<Form("%4d %4d %10.3e %10.3e %10.4e %10.3e %10.3e %10.4e %10.3e %10.3e %10.4e %10.3e %10.3e %10.4e %10.3e %10.3e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e",
                    IT, IQ, 
                    asym_1m1_avg, par[1], err[1],
                    asym_0p1_avg, par[2], err[2],
                    asym_2m1_avg, par[3], err[3],
                    asym_3m1_avg, par[4], err[4],
                    asym_1p1_avg, par[5], err[5],
                    par[0],err[0], Asym, Astat,
                    t_avg, tp_avg, xb_avg, Q2_avg, W_avg, dilute_avg 
                    )
                <<endl;
        }
}
/*}}}*/

/*myUML{{{*/
void myUML(Int_t& npar, Double_t* deriv, Double_t& L, Double_t *par, Int_t flag){
    //Set Parameters
    unsigned int aNu = vPhiH_U.size();
    double lnL_U = 0.0;
    for ( unsigned int ii=0; ii< aNu; ii++ ){
        lnL_U += vWeight_U[ii] * log(func(vFactor_U[ii], vPhiH_U[ii], vPhiS_U[ii], par)) ;
    }

    unsigned int aNd = vPhiH_D.size();
    double lnL_D = 0.0;
    for ( unsigned int ii=0; ii< aNd; ii++ ){
        lnL_D += vWeight_D[ii] * log(func(vFactor_D[ii], vPhiH_D[ii], vPhiS_D[ii], par)) ;
    }

    L = -1.0 * (lnL_U + lnL_D);
}
/*}}}*/

/*func{{{*/
Double_t func(Double_t factor, Double_t phiH, Double_t phiS, Double_t *par){

    double norm = par[0];
    double asym_1m1 = par[1];
    double asym_0p1 = par[2];
    double asym_2m1 = par[3];
    double asym_3m1 = par[4];
    double asym_1p1 = par[5];

    double f = 1.0+norm - POL * factor * (
              asym_1m1*sin(1.*phiH - phiS) 
            + asym_0p1*sin(0.*phiH + phiS) 
            + asym_2m1*sin(2.*phiH - phiS) 
            + asym_3m1*sin(3.*phiH - phiS) 
            + asym_1p1*sin(1.*phiH + phiS) 
            );

    return f;
}
/*}}}*/

/*LoadData{{{*/
void LoadData( Int_t IT=0, Int_t IQ=0){

   Double_t phi_h, phi_s, Photon_Factor;
   Double_t total_acc, total_acc_cor, total_acc_res, MP_res;
   Double_t weight, weight_uu, weight_ut;
   Double_t weight_1m1,weight_2m1,weight_3m1,weight_0p1,weight_1p1;
   Double_t asym_1m1, asym_2m1, asym_3m1, asym_0p1, asym_1p1, sigma_uu;
   Double_t t, tp, t_Para, Q2, xb, W, dilute;
   Int_t Q2BIN;
  
   /*Target Polarization is up{{{*/
   TFile *fu = new TFile(Form("../rootfiles/%s_dvmp_up_t%d_%s.root",bin_name.Data(), IT, type_name.Data()));
   TTree *Tu = (TTree*) gDirectory->Get("T");
   Tu->SetBranchAddress("Phi_cor", &phi_h);/*{{{*/
   Tu->SetBranchAddress("PhiS_cor", &phi_s);//should use the corrected quantity which accounts for fermi-motion etc.
   Tu->SetBranchAddress("Photon_Factor", &Photon_Factor);
   Tu->SetBranchAddress("Q2BIN", &Q2BIN);

   Tu->SetBranchAddress("t_cor", &t);
   Tu->SetBranchAddress("t_Para", &t_Para);
   Tu->SetBranchAddress("W_cor", &W);
   Tu->SetBranchAddress("x_cor", &xb);
   Tu->SetBranchAddress("Qsq_cor", &Q2);
   Tu->SetBranchAddress("dilute", &dilute);

   Tu->SetBranchAddress("total_acc", &total_acc);
   Tu->SetBranchAddress("total_acc_cor", &total_acc_cor);
   Tu->SetBranchAddress("total_acc_res", &total_acc_res);
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
   Double_t weight_sum_up = 0;
   for(int i=0;i<Nu; i++){
     Tu->GetEntry(i);
     if(Q2BIN!=IQ && IQ!=0) continue; //choose the right Q2 bin

     //last chance to apply whatever cuts here
     if(MP_res>1.2) continue;

     tp = t - t_Para;
     phi_h *= Deg2Rad;
     phi_s *= Deg2Rad;
     Photon_Factor = 1.0;
     
     if(isnan(weight_uu) || isinf(weight_uu)) weight_uu=0.0;
     if(weight_uu<1e-12 && weight_uu>1e12) weight_uu=0.0;

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

     weight_ut *= POL*Photon_Factor;
     weight = weight_uu - weight_ut;

     vPhiH_U.push_back(phi_h);
     vPhiS_U.push_back(phi_s);
     vFactor_U.push_back(Photon_Factor);
     vWeight_U.push_back(weight);

     asym_1m1_avg += asym_1m1 * weight_uu;
     asym_0p1_avg += asym_0p1 * weight_uu;
     asym_2m1_avg += asym_2m1 * weight_uu;
     asym_3m1_avg += asym_3m1 * weight_uu;
     asym_1p1_avg += asym_1p1 * weight_uu;

     t_avg += t * weight_uu;
     tp_avg += tp * weight_uu;
     W_avg += W * weight_uu;
     xb_avg += xb * weight_uu;
     Q2_avg += Q2 * weight_uu;
     dilute_avg += dilute * weight_uu;

     weight_sum_up += weight;
 }
   fu->Close();
   /*}}}*/

   /*Target Polarization is down{{{*/
   TFile *fd = new TFile(Form("../rootfiles/%s_dvmp_down_t%d_%s.root",bin_name.Data(), IT, type_name.Data()));
   TTree *Td = (TTree*) gDirectory->Get("T");
   Td->SetBranchAddress("Phi_cor", &phi_h);/*{{{*/
   Td->SetBranchAddress("PhiS_cor", &phi_s);//should use the corrected quantity which accounts for fermi-motion etc.
   Td->SetBranchAddress("Photon_Factor", &Photon_Factor);
   Td->SetBranchAddress("Q2BIN", &Q2BIN);

   Td->SetBranchAddress("t_cor", &t);
   Td->SetBranchAddress("t_Para", &t_Para);
   Td->SetBranchAddress("W_cor", &W);
   Td->SetBranchAddress("x_cor", &xb);
   Td->SetBranchAddress("Qsq_cor", &Q2);
   Td->SetBranchAddress("dilute", &dilute);
   
   Td->SetBranchAddress("total_acc", &total_acc);
   Td->SetBranchAddress("total_acc_cor", &total_acc_cor);
   Td->SetBranchAddress("total_acc_res", &total_acc_res);
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
   Double_t weight_sum_down = 0;
   for(int i=0;i<Nd; i++){
     Td->GetEntry(i);
     if(Q2BIN!=IQ && IQ!=0) continue; //choose the right Q2 bin
     //last chance to apply whatever cuts here
     if(MP_res>1.2) continue;
     tp = t - t_Para;
    
     //Note: In the generator Ahmed switch the sign of the polarization. 
     //In this fit, I fix the absolute polarization values, and rotate the phi_S
     phi_s += 180.0;
     phi_h *= Deg2Rad;
     phi_s *= Deg2Rad;
     Photon_Factor = 1.0;
    
     if(isnan(weight_uu) || isinf(weight_uu)) weight_uu=0.0;
     if(weight_uu<1e-12 && weight_uu>1e12) weight_uu=0.0;
 
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

     weight_ut *= POL*Photon_Factor;
     weight = weight_uu - weight_ut; //follow the HERMES thesis to put a minus sign here

     vPhiH_D.push_back(phi_h);
     vPhiS_D.push_back(phi_s);
     vFactor_D.push_back(Photon_Factor);
     vWeight_D.push_back(weight);
   
     asym_1m1_avg += asym_1m1 * weight_uu;
     asym_0p1_avg += asym_0p1 * weight_uu;
     asym_2m1_avg += asym_2m1 * weight_uu;
     asym_3m1_avg += asym_3m1 * weight_uu;
     asym_1p1_avg += asym_1p1 * weight_uu;

     t_avg += t * weight_uu;
     tp_avg += tp * weight_uu;
     W_avg += W * weight_uu;
     xb_avg += xb * weight_uu;
     Q2_avg += Q2 * weight_uu;
     dilute_avg += dilute * weight_uu;

     weight_sum_down += weight;
}
   fd->Close();
   /*}}}*/

   Double_t weight_sum = weight_sum_up;
   weight_sum +=weight_sum_down;
   
   asym_1m1_avg /= weight_sum;
   asym_0p1_avg /= weight_sum;
   asym_2m1_avg /= weight_sum;
   asym_3m1_avg /= weight_sum;
   asym_1p1_avg /= weight_sum;

   t_avg /= weight_sum;
   tp_avg /= weight_sum;
   W_avg /= weight_sum;
   xb_avg /= weight_sum;
   Q2_avg /= weight_sum;
   dilute_avg /= weight_sum;
   
   Asym = (weight_sum_up - weight_sum_down)/weight_sum * 1./POL;
   Astat = 1./sqrt(weight_sum*0.5);

   cout<<Form("D:  A_1m1=%10.4e A_0p1=%10.4e   A_2m1=%10.4e  A_3m1=%10.4e  A_1p1=%10.4e,, Asym=%10.4e, Astat=%10.4e",
           asym_1m1_avg,asym_0p1_avg,asym_2m1_avg,asym_3m1_avg,asym_1p1_avg, Asym, Astat)<<endl;

}
/*}}}*/
