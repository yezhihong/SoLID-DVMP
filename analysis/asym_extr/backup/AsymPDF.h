/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
  * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/

#ifndef ASYMPDF
#define ASYMPDF

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
 
class AsymPDF : public RooAbsPdf {
public:
  AsymPDF() {} ; 
  AsymPDF(const char *name, const char *title,
	      RooAbsReal& _phi_S,
	      RooAbsReal& _ phi_h,
	      RooAbsReal& _ a,
	      RooAbsReal& _ b);
  AsymPDF(const AsymPDF& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new AsymPDF(*this,newname); }
  inline virtual ~AsymPDF() { }

protected:

  RooRealProxy phi_S ;
  RooRealProxy  phi_h ;
  RooRealProxy  a ;
  RooRealProxy  b ;
  
  Double_t evaluate() const ;

private:

  ClassDef(AsymPDF,1) // Your description goes here...
};
 
#endif
