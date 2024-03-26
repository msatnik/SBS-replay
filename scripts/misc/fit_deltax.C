#include "TH1F.h"
#include "TF1.h"
#include "TAxis.h"
#include "TStyle.h"

void fit_deltax(TH1F *hdx, double Np, double Nn, double mup=-0.8, double mun=0.0, double sigmap=0.1, double sigman=0.1, int bgflag=0, double fit_xmin=-2, double fit_xmax=1.0 ){
  //bgflag = 0: pol4 for the background shape
  //bgflag = 1: gaus for the background shape:

  gStyle->SetOptFit();
  
  double xmin = hdx->GetXaxis()->GetXmin();
  double xmax = hdx->GetXaxis()->GetXmin();

  TF1 *fitfunc_pol4 = new TF1("fitfunc_pol4", "[0]*exp(-0.5*pow((x-[1])/[2],2))+[3]*exp(-0.5*pow((x-[4])/[5],2))+[6]+[7]*x+[8]*x*x+[9]*x*x*x+[10]*pow(x,4)",-10,10);
  TF1 *fitfunc_gaus = new TF1("fitfunc_gaus", "[0]*exp(-0.5*pow((x-[1])/[2],2))+[3]*exp(-0.5*pow((x-[4])/[5],2))+[6]*exp(-0.5*pow((x-[7])/[8],2))",-10,10);
  
  fitfunc_pol4->SetNpx(1000);
  
  double startpar[11] = { Np, mup, sigmap, Nn, mun, sigman, 0, 0, 0, 0, 0 };

  fitfunc_pol4->SetParameters(startpar);

  if( Nn <= 0.1 ){
    fitfunc_pol4->FixParameter(3,0.0);
    fitfunc_pol4->FixParameter(4,0.0);
    fitfunc_pol4->FixParameter(5,1.0);
  }

  double startpar_gaus[9] = {Np, mup, sigmap, Nn, mun, sigman, 0.1*Np, 0.38*mun+0.62*mup, 1.0};

  fitfunc_gaus->SetParameters(startpar_gaus);

  if( bgflag == 0 ){
    hdx->Fit(fitfunc_pol4,"S","",fit_xmin,fit_xmax);
  } else {
    hdx->Fit(fitfunc_gaus, "S","",fit_xmin,fit_xmax);
  }
    

  TF1 *ngaus = new TF1( "ngaus", "[0]*exp(-0.5*pow((x-[1])/[2],2))",-10,10);
  TF1 *pgaus = new TF1( "pgaus", "[0]*exp(-0.5*pow((x-[1])/[2],2))",-10,10);

  for( int i=0; i<=2; i++ ){
    if( bgflag == 0 ){
      ngaus->SetParameter( i, fitfunc_pol4->GetParameter(i+3) );
      pgaus->SetParameter( i, fitfunc_pol4->GetParameter(i) );
    } else {
      ngaus->SetParameter( i, fitfunc_gaus->GetParameter(i+3) );
      pgaus->SetParameter( i, fitfunc_gaus->GetParameter(i) );
    }
  } 

  double Nproton, Nneutron;
  
  double xminp = pgaus->GetParameter(1)-10.0*pgaus->GetParameter(2);
  double xmaxp = pgaus->GetParameter(1)+10.0*pgaus->GetParameter(2);

  double xminn = ngaus->GetParameter(1)-10.0*ngaus->GetParameter(2);
  double xmaxn = ngaus->GetParameter(1)+10.0*ngaus->GetParameter(2);
  
  
  Nproton = pgaus->Integral(xminp,xmaxp) / hdx->GetBinWidth(10);
  Nneutron = ngaus->Integral(xminn, xmaxn) / hdx->GetBinWidth(10);
  
  cout << "Number of neutrons = " << Nneutron << " +/- " << sqrt(Nneutron) << endl;
  cout << "Number of protons = " << Nproton << " +/- " << sqrt(Nproton) << endl;

  double ratio = Nneutron/Nproton;
  double dratio = fabs(ratio) * sqrt( 1.0/Nneutron + 1.0/Nproton );

  cout << "n/p ratio = " << ratio << " +/- " << dratio << endl;

  double FFratio = sqrt(ratio) * 2.793/1.913;
  double dFFratio = dratio / (2.0*FFratio) * 2.793/1.913;

  cout << "GMn/GMp (neglect GEp, GEn) = " << FFratio << " +/- " << dFFratio << endl;
  
  //y = x^(1/2), dy = 1/(2y) dx

}
