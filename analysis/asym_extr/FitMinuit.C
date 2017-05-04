///////////////2//////////////////////
//     DVMP Asymmetries Fitting     //
//     TMinuit Minimization         //
//    ---  Zhihong Ye 04/26/2017    //
//////////////////////////////////////
#include "FitMinuit.h"

//Set to one but keep in mind that the weight_ut should set pol=1 as well.
static const Double_t POL = 0.6*0.865 * 0.9;//He3-Pol, Neutron-Effective-Pol, and Dilution

TString type_name = "simple";
//TString type_name = "mult";
ofstream outf;
Double_t asym_1m1_avg,asym_2m1_avg,asym_3m1_avg,asym_0p1_avg,asym_1p1_avg,Astat;

/*Main{{{*/
Int_t main()
{ 
    gStyle->SetOptFit(1);  
    gStyle->SetOptStat(0);
    
    Int_t iType = 0;
    cout<<"--- Which file ? (1->simple, 2->mult, 3->mult_fsi)  "; cin >> iType;
    if(iType==1) type_name = "simple"; 
    if(iType==2) type_name = "mult"; 
    if(iType==3) type_name = "mult_fsi"; 

    Int_t IT = 0, IQ = 0;
    cout<<"--- Which t-bin? (1--8, 0 for all)  "; cin >> IT;
    if(IT>0){
        cout<<"--- Which Q2-bin? (1--2)  "; cin >> IQ;
        LoadData(IT, IQ);

        /////////////////////////////
        //TMinuit Minimization
        /////////////////////////////
        double fyPar[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
        DoMinuit(fyPar, IT, IQ, 0);
    }
    else{
        outf.open(Form("dvmp_par_%s_real.dat", type_name.Data()));
        outf<<Form("%4s %4s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s ",
                "#t", "#Q2", 
                "A_AVG", "A_1M1", "dA_1M1", 
                "A_AVG", "A_2M1", "dA_2M1", 
                "A_AVG", "A_3M1", "dA_3M1", 
                "A_AVG", "A_0P1", "dA_0P1", 
                "A_AVG", "A_1P1", "dA_1P1")
            <<endl;


        for(int i=1; i<=7;i++){
            //for(int j=0; j<=2;j++){
            for(int j=0; j<1;j++){
                IT = i;
                IQ = j;

                /////////////////////////////
                //Load Data from MC Rootfiles
                /////////////////////////////
                LoadData(IT, IQ);

                /////////////////////////////
                //TMinuit Minimization
                /////////////////////////////
                double fyPar[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
                DoMinuit(fyPar, IT, IQ, 1);
            }
        }

        outf.close();
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
	static const Int_t iParaNum=5;

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
	TString namePar[iParaNum]={"a1m1","a2m1","a3m1","a0p1","a1p1"};
	
        //Set Initial Values of Parameters
        Double_t inPar[iParaNum]={par[0],par[1],par[2],par[3],par[4]};//,par[5]};   

	//Set Stepsize of Parameters
	Double_t step[iParaNum]={1.0e-5,1.0e-5,1.0e-5,1.0e-5,1.0e-5};

	//Set Min of Parameters, value==0 for No-Bound
	Double_t minVal[iParaNum]={0,0,0,0,0};

	//Set Max of Parameters, value==0 for No-Bound
	Double_t maxVal[iParaNum]={0,0,0,0,0};

	//Initializing Parameters
	for(int ii=0;ii<iParaNum;ii++){
            pMinuit->DefineParameter(ii,namePar[ii], inPar[ii], step[ii], minVal[ii], maxVal[ii]);
	}

	//Fix parameters, starting from 1 not 0
	//arglist[0] = 2;   pMinuit->mnexcm("FIX ", arglist ,1,ierflg);
	//arglist[0] = 3;   pMinuit->mnexcm("FIX ", arglist ,1,ierflg);
	//arglist[0] = 4;   pMinuit->mnexcm("FIX ", arglist ,1,ierflg);
	//arglist[0] = 5;   pMinuit->mnexcm("FIX ", arglist ,1,ierflg);

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

        cout<<Form("Real:         A1=%6.3f,  A2=%6.3f,  A3=%6.3f,  A4=%6.3f   A5=%6.3f, Astat=%10.4e",
                asym_1m1_avg,asym_2m1_avg,asym_3m1_avg,asym_0p1_avg,asym_1p1_avg, Astat)<<endl;

        cout<<Form("Asymmetries:  A1=%6.3f,  A2=%6.3f,  A3=%6.3f,  A4=%6.3f,  A5=%6.3f", 
                par[0], par[1], par[2], par[3], par[4])<<endl;
        cout<<Form("     Errors: dA1=%6.3f, dA2=%6.3f, dA3=%6.3f, dA4=%6.3f, dA5=%6.3f", 
        err[0], err[1], err[2], err[3], err[4])<<endl;
        cout<<"-----------------------------------------------------------------------"<<endl;
        cout<<"-----------------------------------------------------------------------"<<endl<<endl;
        
        if(SAVE==1){
            outf<<Form("%4d %4d %10.3e %10.3e %10.4e %10.3e %10.3e %10.4e %10.3e %10.3e %10.4e %10.3e %10.3e %10.4e %10.3e %10.3e %10.4e %10.4e",
                    IT, IQ, 
                    asym_1m1_avg, par[0], err[0],
                    asym_2m1_avg, par[1], err[1],
                    asym_3m1_avg, par[2], err[2],
                    asym_0p1_avg, par[3], err[3],
                    asym_1p1_avg, par[4], err[4],
                    Astat)
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

    double asym_1m1 = par[0];
    double asym_2m1 = par[1];
    double asym_3m1 = par[2];
    double asym_0p1 = par[3];
    double asym_1p1 = par[4];

    double f = 1.0 - POL * factor * (
              asym_1m1*sin(1.*phiH - phiS) 
            + asym_2m1*sin(2.*phiH - phiS) 
            + asym_3m1*sin(3.*phiH - phiS) 
            + asym_0p1*sin(0.*phiH + phiS) 
            + asym_1p1*sin(1.*phiH + phiS) 
            );

    return f;
}
/*}}}*/

/*LoadData{{{*/
void LoadData( Int_t IT=0, Int_t IQ=0){

   Double_t phi_h, phi_s, Photon_Factor;
   Double_t weight, weight_uu, weight_ut;
   Double_t weight_1m1,weight_2m1,weight_3m1,weight_0p1,weight_1p1;
   Double_t asym_1m1, asym_2m1, asym_3m1, asym_0p1, asym_1p1, sigma_uu;
   Int_t Q2BIN;
  
  //Target Polarization is up/*{{{*/
   TFile *fu = new TFile(Form("../rootfiles/dvmp_up_t%d_%s.root", IT, type_name.Data()));
   TTree *Tu = (TTree*) gDirectory->Get("T");
   
   Tu->SetBranchAddress("Phi", &phi_h);
   Tu->SetBranchAddress("PhiS_corr", &phi_s);//should use the corrected quantity which accounts for fermi-motion etc.
   Tu->SetBranchAddress("Photon_Factor", &Photon_Factor);
   Tu->SetBranchAddress("Q2BIN", &Q2BIN);
   
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
   Tu->SetBranchAddress("Sigma_UU", &sigma_uu);
   
   Int_t Nu = Tu->GetEntries();
   int weight_sum = 0.0;
   for(int i=0;i<Nu; i++){
     Tu->GetEntry(i);
     if(Q2BIN!=IQ && IQ!=0) continue; //choose the right Q2 bin
    
     phi_h *= Deg2Rad;
     phi_s *= Deg2Rad;
     Photon_Factor = 1.0;

     weight_1m1 = weight_uu * (asym_1m1* sin(1.*phi_h - phi_s) );
     weight_2m1 = weight_uu * (asym_2m1* sin(2.*phi_h - phi_s) );
     weight_3m1 = weight_uu * (asym_3m1* sin(3.*phi_h - phi_s) );
     weight_0p1 = weight_uu * (asym_0p1* sin(0.*phi_h + phi_s) );
     weight_1p1 = weight_uu * (asym_1p1* sin(1.*phi_h + phi_s) );
     
     weight_ut = weight_1m1
               + weight_2m1
               + weight_3m1
               + weight_0p1
               + weight_1p1;

     weight_ut *= POL*Photon_Factor;
     weight = weight_uu - weight_ut;

     vPhiH_U.push_back(phi_h);
     vPhiS_U.push_back(phi_s);
     vFactor_U.push_back(Photon_Factor);
     vWeight_U.push_back(weight);
     vWeight_UU_U.push_back(weight_uu);
     vWeight_UT_U.push_back(weight_ut);
     vWeight_1M1_U.push_back(weight_1m1);
     vWeight_2M1_U.push_back(weight_2m1);
     vWeight_3M1_U.push_back(weight_3m1);
     vWeight_0P1_U.push_back(weight_0p1);
     vWeight_1P1_U.push_back(weight_1p1);
  
     asym_1m1_avg += asym_1m1 * weight_uu;
     asym_2m1_avg += asym_2m1 * weight_uu;
     asym_3m1_avg += asym_3m1 * weight_uu;
     asym_0p1_avg += asym_0p1 * weight_uu;
     asym_1p1_avg += asym_1p1 * weight_uu;
     weight_sum += weight_uu;
 }
   fu->Close();
   asym_1m1_avg /= weight_sum;
   asym_2m1_avg /= weight_sum;
   asym_3m1_avg /= weight_sum;
   asym_0p1_avg /= weight_sum;
   asym_1p1_avg /= weight_sum;
   Astat = 1./sqrt(weight_sum);

   cout<<Form("U:  A_1m1=%10.4e  A_2m1=%10.4e  A_3m1=%10.4e  A_0p1=%10.4e  A_1p1=%10.4e, Astat=%10.4e",
           asym_1m1_avg,asym_2m1_avg,asym_3m1_avg,asym_0p1_avg,asym_1p1_avg, Astat)<<endl;

   //outf<<Form("U:  A_1m1=%10.4e  A_2m1=%10.4e  A_3m1=%10.4e  A_0p1=%10.4e  A_1p1=%10.4e, Astat=%10.4e",
           //asym_1m1_avg,asym_2m1_avg,asym_3m1_avg,asym_0p1_avg,asym_1p1_avg, 1./sqrt(weight_sum))<<endl;

   /*}}}*/

  //Target Polarization is down/*{{{*/
   TFile *fd = new TFile(Form("../rootfiles/dvmp_down_t%d_%s.root", IT, type_name.Data()));
   TTree *Td = (TTree*) gDirectory->Get("T");
   
   Td->SetBranchAddress("Phi", &phi_h);
   Td->SetBranchAddress("PhiS_corr", &phi_s);//should use the corrected quantity which accounts for fermi-motion etc.
   Td->SetBranchAddress("Photon_Factor", &Photon_Factor);
   Td->SetBranchAddress("Q2BIN", &Q2BIN);
   
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
   Td->SetBranchAddress("Sigma_UU", &sigma_uu);

   Int_t Nd = Td->GetEntries();
   weight_sum = 0.0;
   for(int i=0;i<Nd; i++){
     Td->GetEntry(i);
     if(Q2BIN!=IQ && IQ!=0) continue; //choose the right Q2 bin
     //Note: In the generator Ahmed switch the sign of the polarization. 
     //In this fit, I fix the absolute polarization values, and rotate the phi_S
     phi_s += 180.0;
     phi_h *= Deg2Rad;
     phi_s *= Deg2Rad;
     Photon_Factor = 1.0;
 
     weight_1m1 = weight_uu * (asym_1m1* sin(1.*phi_h - phi_s) );
     weight_2m1 = weight_uu * (asym_2m1* sin(2.*phi_h - phi_s) );
     weight_3m1 = weight_uu * (asym_3m1* sin(3.*phi_h - phi_s) );
     weight_0p1 = weight_uu * (asym_0p1* sin(0.*phi_h + phi_s) );
     weight_1p1 = weight_uu * (asym_1p1* sin(1.*phi_h + phi_s) );
     
     weight_ut = weight_1m1
               + weight_2m1
               + weight_3m1
               + weight_0p1
               + weight_1p1;

     weight_ut *= POL*Photon_Factor;
     weight = weight_uu - weight_ut; //follow the HERMES thesis to put a minus sign here

     vPhiH_D.push_back(phi_h);
     vPhiS_D.push_back(phi_s);
     vFactor_D.push_back(Photon_Factor);
     vWeight_D.push_back(weight);
     vWeight_UU_D.push_back(weight_uu);
     vWeight_UT_D.push_back(weight_ut);
     vWeight_1M1_D.push_back(weight_1m1);
     vWeight_2M1_D.push_back(weight_2m1);
     vWeight_3M1_D.push_back(weight_3m1);
     vWeight_0P1_D.push_back(weight_0p1);
     vWeight_1P1_D.push_back(weight_1p1);
   
     asym_1m1_avg += asym_1m1 * weight_uu;
     asym_2m1_avg += asym_2m1 * weight_uu;
     asym_3m1_avg += asym_3m1 * weight_uu;
     asym_0p1_avg += asym_0p1 * weight_uu;
     asym_1p1_avg += asym_1p1 * weight_uu;
     weight_sum += weight_uu;
}
   fd->Close();
   asym_1m1_avg /= weight_sum;
   asym_2m1_avg /= weight_sum;
   asym_3m1_avg /= weight_sum;
   asym_0p1_avg /= weight_sum;
   asym_1p1_avg /= weight_sum;
   Astat = 1./sqrt(weight_sum);

   cout<<Form("D:  A_1m1=%10.4e  A_2m1=%10.4e  A_3m1=%10.4e  A_0p1=%10.4e  A_1p1=%10.4e, Astat=%10.4e",
           asym_1m1_avg,asym_2m1_avg,asym_3m1_avg,asym_0p1_avg,asym_1p1_avg, Astat)<<endl;

   //outf<<Form("D:  A_1m1=%10.4e  A_2m1=%10.4e  A_3m1=%10.4e  A_0p1=%10.4e  A_1p1=%10.4e, Astat=%10.4e",
           //asym_1m1_avg,asym_2m1_avg,asym_3m1_avg,asym_0p1_avg,asym_1p1_avg, 1./sqrt(weight_sum))<<endl;

   /*}}}*/

}
/*}}}*/
