#include "TH1D.h"
#include "TF1.h"
#include "TChain.h"
#include "TTree.h"
#include "TCut.h"
#include "TCanvas.h"

TH1D *hW2elastic;

void FitGaus_FWHM( TH1D *htest, double thresh=0.5 ){
  int binmax = htest->GetMaximumBin();
  int binlow = binmax, binhigh = binmax;

  double max = htest->GetBinContent(binmax);

  while( htest->GetBinContent(binlow) >= thresh*max && binlow > 1 ){binlow--;}
  while( htest->GetBinContent(binhigh) >= thresh*max && binhigh < htest->GetNbinsX() ){ binhigh++; }

  double xlow = htest->GetBinLowEdge(binlow);
  double xhigh = htest->GetBinLowEdge(binhigh+1);

  htest->Fit("gaus","q0S","",xlow, xhigh);
}

double fitfunc( double *x, double *par ){
  double W2 = x[0];

  double Norm = par[0];
  // Let's use (up to) 6th order polynomial for the background:
  double bg = 0.0;
  for( int i = 0; i<7; i++ ){
    bg += par[i+1]*pow(W2,i);
  }

  double signal = Norm * hW2elastic->Interpolate(W2);

  return signal + bg; 
  
}

//For the dx efficiency method, fit a Gaussian 
double fitfunc_dx( double *x, double *par ){
  double dx = x[0];

  //Fit 4th-order polynomial to sidebands to get background:

  //These should ordinarily be fixed:
  double dxpeak = par[0];
  double dxsigma = par[1];

  if( fabs( dx - dxpeak ) < 3.5*dxsigma ) TF1::RejectPoint();
  
  double poly = 0.0;
  for( int i=0; i<5; i++ ){
    poly += par[i+2]*pow(dx,i);
  }

  return poly;

}

void Fit_W2_inclusive( TH1D *hW2all, TH1D *hW2elastic_temp, int order_bg=4, double W2min=-1.0, double W2max=2.0 ){
  TF1 *FitFunc = new TF1( "FitFunc", fitfunc, -1, 3, 8 );

  FitFunc->SetNpx(1000);
  
  if( order_bg > 6 ) order_bg = 6;

  hW2elastic = new TH1D( *hW2elastic_temp );
  
  if( hW2all != nullptr && hW2elastic != nullptr ){
    double startpar[] = { 1.0, 0, 0, 0, 0, 0, 0, 0 };

    FitFunc->SetParameters(startpar);

   
    
    int ipar=0;
    while( ipar <= order_bg ){
      ipar++;
    }

    while( ipar <= 6 ){
      FitFunc->FixParameter( ipar+1, 0.0 );
      ipar++;
    }

    hW2all->Fit( FitFunc, "", "", W2min, W2max );

    hW2all->SetLineColor(1);
    hW2all->SetLineWidth(2);
    
    hW2all->DrawCopy();
    //hW2elastic->Scale( FitFunc->GetParameter(0) );

    TH1D *hW2elastic_copy = new TH1D(*hW2elastic); 
    //    hW2elastic->Draw("SAME");

    hW2elastic_copy->SetName("hW2elastic_copy");

    hW2elastic_copy->Scale( FitFunc->GetParameter(0) );
    
    TH1D *hW2_bg = new TH1D( *hW2all );

    hW2_bg->SetName("hW2_bg");

    hW2_bg->Add( hW2all, hW2elastic_copy, 1., -1. );

    hW2_bg->SetMarkerStyle(21);
    hW2_bg->SetMarkerColor(4);
    hW2_bg->SetLineColor(4);
    
    hW2_bg->DrawCopy("SAME");

    hW2elastic_copy->SetMarkerStyle(20);
    hW2elastic_copy->SetMarkerColor(2);
    hW2elastic_copy->SetLineColor(2);
    
    hW2elastic_copy->DrawCopy("same");

    *hW2elastic_temp = *hW2elastic_copy;
    
  }


}

void GetElasticCounts( const char *rootfilename, double W2min=0.4, double W2max = 1.4, TCut cut="", TCut fiducialcut="", double W2min_fit=0, double W2max_fit=1.4, double dymin=-0.75, double dymax=0.3 ){
  TChain *C = new TChain("Tout");

  C->Add(rootfilename);

  TH1D *hW2_nocut = new TH1D("hW2_nocut", ";W^{2} (GeV^{2});", 100, -1.0,3.0 );
  TH1D *hW2_HCALcut = new TH1D("hW2_HCALcut", ";W^{2} (GeV^{2});", 100, -1.0, 3.0 );
  TH1D *hW2_anticut = new TH1D("hW2_anticut", ";W^{2} (GeV^{2});", 100, -1.0, 3.0 );

  TH1D *hdx_Wcut = new TH1D("hdx_Wcut", ";#Deltax (m);", 300,-2,1);

  
  C->Project("hW2_nocut", "W2", fiducialcut);
  C->Project("hW2_HCALcut", "W2", cut&&fiducialcut);
  C->Project("hW2_anticut", "W2", !cut&&fiducialcut);
 
  TString cutstring; 

  cutstring.Form("%g<W2&&W2<%g",W2min,W2max);
  //  cutstring.Form("%g<W2&&W2<%g&&%g<deltay&&deltay<%g",W2min,W2max, dymin, dymax);

  TCut Wcut = cutstring.Data();

  Wcut.Print();

  //Add loose deltay cut:
  C->Project("hdx_Wcut", "deltax",Wcut&&fiducialcut); 

  TCanvas *c1 = new TCanvas("c1","c1",1200,750);
  c1->Divide(3,1,.001,.001);

  c1->cd(1);
  
  Fit_W2_inclusive( hW2_nocut, hW2_HCALcut, 6, W2min_fit, W2max_fit );

  double Counts,dCounts;

  Counts = hW2_HCALcut->IntegralAndError( hW2_HCALcut->FindBin(W2min), hW2_HCALcut->FindBin(W2max), dCounts );

  cout << "Elastic yield (total with fiducial cut) = " << Counts << " +/- " << dCounts << endl;

  double CountsAcut,dCountsAcut;

  c1->cd(2);
  
  Fit_W2_inclusive( hW2_anticut, hW2_HCALcut, 6, W2min_fit, W2max_fit );

  CountsAcut = hW2_HCALcut->IntegralAndError( hW2_HCALcut->FindBin(W2min), hW2_HCALcut->FindBin(W2max), dCountsAcut );

  cout << "Elastic yield (passed fiducial cut, failed HCAL cut) = " << CountsAcut << " +/- " << dCountsAcut << endl;

  //Total inclusive counts  

  double efficiency = 1.-CountsAcut/Counts;
  double defficiency = sqrt(efficiency*(1.-efficiency)/Counts);
  
  cout << "HCAL effective proton efficiency = " << efficiency << " +/- " << defficiency << endl;

  c1->cd(3);
  hdx_Wcut->Draw();

  FitGaus_FWHM( hdx_Wcut, 0.5 );
  
  double dxmean = ( (TF1*) hdx_Wcut->GetListOfFunctions()->FindObject("gaus") )->GetParameter("Mean");
  double dxsigma = ( (TF1*) hdx_Wcut->GetListOfFunctions()->FindObject("gaus") )->GetParameter("Sigma");

  TF1 *func_dxsideband = new TF1("func_dxsideband", fitfunc_dx, -2, 1, 7 );
  func_dxsideband->FixParameter(0,dxmean);
  func_dxsideband->FixParameter(1,dxsigma);

  for( int ipar=2; ipar<7; ipar++ ){
    func_dxsideband->SetParameter(ipar,0.0);
  }
  
  hdx_Wcut->Fit(func_dxsideband,"S","",-2,1);

  double bg, dbg;

  bg = func_dxsideband->Integral(-2,1) / hdx_Wcut->GetBinWidth(1);
  dbg = func_dxsideband->IntegralError(-2,1) / hdx_Wcut->GetBinWidth(1);

  double tot,dtot;
  //tot = hdx_Wcut->IntegralAndError(1,hdx_Wcut->GetNbinsX(),dtot);
  double elastic = hdx_Wcut->IntegralAndError(1,hdx_Wcut->GetNbinsX(),dtot) - bg;

  //VERY rough error = 
  double delastic = sqrt(dtot*dtot+dbg*dbg);

  double eff_dxside = elastic/Counts;
  double deff_dxside = sqrt(eff_dxside*(1.-eff_dxside)/Counts);
  
  cout << "Elastic counts from dx sideband fit = " << elastic << " +/- " << delastic << ", background = " 
       << bg << " +/- " << dbg << ", naive efficiency = " << elastic/Counts 
       << " +/- " << deff_dxside << endl;
}
