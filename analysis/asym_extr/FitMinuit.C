///////////////2//////////////////////
//     DVMP Asymmetries Fitting     //
//     TMinuit Minimization         //
//    ---  Zhihong Ye 04/26/2017    //
//////////////////////////////////////
#include "FitMinuit.h"
//Fake Asymmetries for testing only
const double A1 = -0.3;//1m1
const double A2 = -0.2;//2m1
const double A3 = -0.1;//3m1
const double A4 = 0.1;//0p1
const double A5 = 0.15;//1p1

/*Main{{{*/
Int_t main()
{ 
	gStyle->SetOptFit(1);  
	gStyle->SetOptStat(0);

        Int_t I = 0;
        cout<<"--- Which t-bin? (1--7)"; cin >> I;

        //Load Data from MC Rootfiles
	LoadData(I);

	/////////////////////////////
	//TMinuit Minimization
	/////////////////////////////
	double fyPar[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
        DoMinuit(fyPar);
	return 0;
}
/*}}}*/

/*DoMunuit{{{*/
void DoMinuit(double *par)
{
	cout << "      *************************************************" << endl; 
	cout << "      *          Minimization Fitting for             *" << endl;
	cout << "      *              DVMP Asymmetries                 *" << endl;
	cout << "      *             Z. Ye 09/20/2010                  *" << endl;
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
	TString namePar[iParaNum]={"a1m1","a2m1","a3m1","a0p1","a1p1"};//, "a2p1"};
	
        //Set Initial Values of Parameters
        Double_t inPar[iParaNum]={par[0],par[1],par[2],par[3],par[4]};//,par[5]};   

	//Set Stepsize of Parameters
	Double_t step[iParaNum]={ 0.00000001,0.00000001,0.00001,0.00000001,0.00000001};//,0.00000001};

	//Set Min of Parameters, value==0 for No-Bound
	// Double_t minVal[iParaNum]={ 0.0001,0.0001,100.0,0.0001,0.0001};
	Double_t minVal[iParaNum]={ 0.0,0.0,0.0,0.0,0.0};

	//Set Max of Parameters, value==0 for No-Bound
	Double_t maxVal[iParaNum]={ 0.0,0.0,0.0,0.0,0.0};//,0.0};

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
	arglist[1] = 0.01;
	pMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);

	//Get Results
	for(int ii = 0; ii < iParaNum; ii++)
	{
		pMinuit->GetParameter(ii,outPar[ii],err[ii]);    
		par[ii] = outPar[ii];
	}
	//Put the results into a file
	cerr<<Form("Asymmetries:  A1=%6.3f,  A2=%6.3f,  A3=%6.3f,  A4=%6.3f,  A5=%6.3f", 
                par[0], par[1], par[2], par[3], par[4])<<endl;
	cerr<<Form("     Errors: dA1=%6.3f, dA2=%6.4f, dA3=%6.4f, dA4=%6.4f, dA5=%6.4f", 
                err[0], err[1], err[2], err[3], err[4])<<endl;

}
/*}}}*/

/*myUML{{{*/
void myUML(Int_t& npar, Double_t* deriv, Double_t& L, Double_t *par, Int_t flag){
    //Set Parameters
    int aNu = vPhiH_U.size();
    double lnL_U = 0.0;
    for ( unsigned int ii=0; ii< aNu; ii++ ){
        lnL_U += vWeight_UT_U[ii] * log(func(vPhiH_U[ii], vPhiS_U[ii], par)) ;
    }

    int aNd = vPhiH_D.size();
    double lnL_D = 0.0;
    for ( unsigned int ii=0; ii< aNd; ii++ ){
        lnL_D += vWeight_UT_D[ii] * log(func(vPhiH_D[ii], vPhiS_D[ii], par)) ;
    }

    L = -1.0 * (lnL_U + lnL_D);
}
/*}}}*/

/*func{{{*/
Double_t func(Double_t phiH, Double_t phiS, Double_t *par){

    double asym_1m1 = par[0];
    double asym_2m1 = par[1];
    double asym_3m1 = par[2];
    double asym_0p1 = par[3];
    double asym_1p1 = par[4];

    double f = 1.0 + POL * (
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
void LoadData( Int_t I=0){

   Double_t phi_h, phi_s;
   Double_t weight, weight_uu, weight_ut;
   Double_t weight_1m1,weight_2m1,weight_3m1,weight_0p1,weight_1p1;
  
  //Target Polarization is up/*{{{*/
   TFile *fu = new TFile(Form("../rootfiles/dvmp_up_t%d.root", I));
   TTree *Tu = (TTree*) gDirectory->Get("T");
   
   Tu->SetBranchAddress("Phi", &phi_h);
   Tu->SetBranchAddress("PhiS", &phi_s);
   Tu->SetBranchAddress("weight", &weight);
   Tu->SetBranchAddress("weight_uu", &weight_uu);
   Tu->SetBranchAddress("weight_ut", &weight_ut);
   Tu->SetBranchAddress("weight_1m1", &weight_1m1);
   Tu->SetBranchAddress("weight_2m1", &weight_2m1);
   Tu->SetBranchAddress("weight_3m1", &weight_3m1);
   Tu->SetBranchAddress("weight_0p1", &weight_0p1);
   Tu->SetBranchAddress("weight_1p1", &weight_1p1);

   Int_t Nu = Tu->GetEntries();
   for(int i=0;i<Nu; i++){
     Tu->GetEntry(i);
     phi_h *= Deg2Rad;
     phi_s *= Deg2Rad;

     //fake weights based on fake asymmetries for testing
     weight_1m1 = A1 * sin(1.0*phi_h - phi_s );
     weight_2m1 = A2 * sin(2.0*phi_h - phi_s );
     weight_3m1 = A3 * sin(3.0*phi_h - phi_s );
     weight_0p1 = A4 * sin(0.0*phi_h + phi_s );
     weight_1p1 = A5 * sin(1.0*phi_h + phi_s );
     weight_ut = weight_uu * (1.0+weight_1m1+weight_2m1+weight_3m1+weight_0p1+weight_1p1);

     vPhiH_U.push_back(phi_h);
     vPhiS_U.push_back(phi_s);
     vWeight_UU_U.push_back(weight_uu);
     vWeight_UT_U.push_back(weight_ut);
     vWeight_1M1_U.push_back(weight_1m1);
     vWeight_2M1_U.push_back(weight_2m1);
     vWeight_3M1_U.push_back(weight_3m1);
     vWeight_0P1_U.push_back(weight_0p1);
     vWeight_1P1_U.push_back(weight_1p1);
   }
   fu->Close();/*}}}*/

  //Target Polarization is down/*{{{*/
   TFile *fd = new TFile(Form("../rootfiles/dvmp_down_t%d.root", I));
   TTree *Td = (TTree*) gDirectory->Get("T");
   
   Td->SetBranchAddress("Phi", &phi_h);
   Td->SetBranchAddress("PhiS", &phi_s);
   Td->SetBranchAddress("weight", &weight);
   Td->SetBranchAddress("weight_uu", &weight_uu);
   Td->SetBranchAddress("weight_ut", &weight_ut);
   Td->SetBranchAddress("weight_1m1", &weight_1m1);
   Td->SetBranchAddress("weight_2m1", &weight_2m1);
   Td->SetBranchAddress("weight_3m1", &weight_3m1);
   Td->SetBranchAddress("weight_0p1", &weight_0p1);
   Td->SetBranchAddress("weight_1p1", &weight_1p1);

   Int_t Nd = Td->GetEntries();
   for(int i=0;i<Nd; i++){
     Td->GetEntry(i);
     phi_s += 180.0;
     phi_h *= Deg2Rad;
     phi_s *= Deg2Rad;
      
     //fake weights based on fake asymmetries for testing
     weight_1m1 = A1 * sin(1.0*phi_h - phi_s );
     weight_2m1 = A2 * sin(2.0*phi_h - phi_s );
     weight_3m1 = A3 * sin(3.0*phi_h - phi_s );
     weight_0p1 = A4 * sin(0.0*phi_h + phi_s );
     weight_1p1 = A5 * sin(1.0*phi_h + phi_s );
     weight_ut = weight_uu * (1.0+weight_1m1+weight_2m1+weight_3m1+weight_0p1+weight_1p1);
     
     vPhiH_D.push_back(phi_h);
     vPhiS_D.push_back(phi_s);
     vWeight_UU_D.push_back(weight_uu);
     vWeight_UT_D.push_back(weight_ut);
     vWeight_1M1_D.push_back(weight_1m1);
     vWeight_2M1_D.push_back(weight_2m1);
     vWeight_3M1_D.push_back(weight_3m1);
     vWeight_0P1_D.push_back(weight_0p1);
     vWeight_1P1_D.push_back(weight_1p1);
   }
   fd->Close();/*}}}*/

}
/*}}}*/
