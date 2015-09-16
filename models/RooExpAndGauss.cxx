
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

 #include <iostream>

 #include "RooExpAndGauss.h"
 #include "RooAbsReal.h"

 ClassImp(RooExpAndGauss)

 RooExpAndGauss::RooExpAndGauss(const char *name, const char *title,
                        RooAbsReal& _x,
                        RooAbsReal& _sh_mean,
                        RooAbsReal& _sh_sigma,
                        RooAbsReal& _sh_trans) :
   RooAbsPdf(name,title),
   x("x","x",this,_x),
   sh_mean("sh_mean","sh_mean",this,_sh_mean),
   sh_sigma("sh_sigma","sh_sigma",this,_sh_sigma),
   sh_trans("sh_trans","sh_trans",this,_sh_trans)
 {
 }


 RooExpAndGauss::RooExpAndGauss(const RooExpAndGauss& other, const char* name) :
   RooAbsPdf(other,name),
   x("x",this,other.x),
   sh_mean("sh_mean",this,other.sh_mean),
   sh_sigma("sh_sigma",this,other.sh_sigma),
   sh_trans("sh_trans",this,other.sh_trans)
 {
 }



 Double_t RooExpAndGauss::evaluate() const
 {
   // ENTER EXPRESSION IN TERMS OF VARIABLE ARGUMENTS HERE
   double beta = (sh_mean - sh_trans)/(sh_sigma*sh_sigma);
   double c = exp(-0.5*pow((sh_trans-sh_mean)/sh_sigma,2))*exp(-beta*sh_trans);
//   if (x <= sh_trans){
//    return c*exp(beta*x);
//  };
//
//  return exp(-0.5*pow((x-sh_mean)/sh_sigma,2)); //think a 1/2 was missed before BsMuMuPy v4.4

  return double(x <= sh_trans) * c*exp(beta*x) + double(x > sh_trans) * exp(-0.5*pow((x-sh_mean)/sh_sigma,2));

 }



