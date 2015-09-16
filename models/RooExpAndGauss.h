/*****************************************************************************
 * Project: Erasmus                                                           *
 * Package: SomeMassModels                                                    *
 *    File: $Id: RooExpAndGauss.h,v 1.13 2007/07/12 20:30:49 diegoms Exp $
 * Authors:                                                                  *
 *           Diego Martinez Santos          diego.martinez.santos@cern.ch
 * Purpose:
 * Description of the B->JPsiX background
 *                                                                           *
 *****************************************************************************/


 // -- CLASS DESCRIPTION [PDF] --
 // This function describes approximately the bkg from B->JpsiX in channels
// like B+-->JpsiK+.


#ifndef ROOEXPANDGAUSS
#define ROOEXPANDGAUSS

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"

class RooExpAndGauss : public RooAbsPdf {
public:
  RooExpAndGauss() { } ;
  RooExpAndGauss(const char *name, const char *title,
	      RooAbsReal& _x,
	      RooAbsReal& _sh_mean,
	      RooAbsReal& _sh_sigma,
	      RooAbsReal& _sh_trans);
  RooExpAndGauss(const RooExpAndGauss& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooExpAndGauss(*this,newname); }
  inline virtual ~RooExpAndGauss() { }

protected:

  RooRealProxy x ;
  RooRealProxy sh_mean ;
  RooRealProxy sh_sigma ;
  RooRealProxy sh_trans ;

  Double_t evaluate() const ;

private:

  ClassDef(RooExpAndGauss,1) // Your description goes here...
};

#endif
