/*Headers{{{*/

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>
#include <iomanip>
#include <string>
#include <cstdio>

#include "math.h"
#include "stdio.h"
#include <stdlib.h>

#include "TGraph.h"
#include "TGraphErrors.h"
#include "TPaveStats.h"
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TF1.h>
#include <TH1.h>
#include <TAxis.h>
#include <TLine.h>
#include <TLatex.h>
#include <TLegend.h>
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TObject.h"
#include "TSystem.h"
#include <TVirtualFitter.h>
#include <TMinuit.h>
//#include <TFitterMinuit.h>
#include <TApplication.h>
#include <vector>
using namespace std;

/*}}}*/

void myUML(Int_t& npar, Double_t* deriv, Double_t& f, Double_t *par, Int_t flag);
Double_t func(Double_t, Double_t , Double_t *par);
void LoadData(Int_t, Int_t);
void DoMinuit(double *par, Int_t, Int_t, Int_t);

static const Double_t POL = 1.0;
static const Double_t pi = 3.1415926;
static const Double_t Deg2Rad = pi/180.;
static const Double_t Rad2Deg = 180./pi;
Double_t outPar[5], err[5];

vector <double> vPhiH_U;
vector <double> vPhiS_U;
vector <double> vWeight_UU_U;
vector <double> vWeight_UT_U;
vector <double> vWeight_1M1_U;
vector <double> vWeight_2M1_U;
vector <double> vWeight_3M1_U;
vector <double> vWeight_0P1_U;
vector <double> vWeight_1P1_U;
vector <double> vWeight_2P1_U;

vector <double> vPhiH_D;
vector <double> vPhiS_D;
vector <double> vWeight_UU_D;
vector <double> vWeight_UT_D;
vector <double> vWeight_1M1_D;
vector <double> vWeight_2M1_D;
vector <double> vWeight_3M1_D;
vector <double> vWeight_0P1_D;
vector <double> vWeight_1P1_D;
vector <double> vWeight_2P1_D;
