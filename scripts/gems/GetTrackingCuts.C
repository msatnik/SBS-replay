#include "TH1D.h"
#include "TH2D.h"
#include "TChain.h"
#include "TTree.h"
#include "TCut.h"
#include "TEventList.h"
#include "TFile.h"
#include "TROOT.h"
#include "TF1.h"
#include <iostream>
#include <fstream>
//#include <ROOT/RDataFrame.hxx>

//using namespace ROOT;

vector<double> threshsamp_module;
vector<double> threshstrip_module;
vector<double> threshclust_module;
vector<TH1D*> hsamp;
vector<TH1D*> hstrip;
vector<TH1D*> hclust;
vector<TH1D*> htmaxU;
vector<TH1D*> htmaxV;


void DrawHist(TH1D *hist, int flag = 0, double value = 0){

  //flag = 0 for calorimeter plots
  //flag = 1 for ADC plots

  hist->Draw();

  TPaveText *pt;
  if(flag == 0){
    pt = new TPaveText(0.12,0.75,0.35,0.88,"ndc");
    pt->AddText(Form("Center = %g",hist->GetMean()));
    pt->AddText(Form("Width = %g",hist->GetRMS() * 4.5));
  }
  else if(flag == 1){
    pt = new TPaveText(0.55,0.75,0.75,0.88,"ndc");
    pt->AddText(Form("Threshold = %g",value));
  }
  pt->SetFillColor(0);
  pt->Draw("same");


}

void DrawThreshPlots(TCanvas *c,vector<TH1D*> hsamp,vector<TH1D*> hstrip,vector<TH1D*> hclust,int imod){

  c->cd(1);
  hsamp[imod]->Draw();
  
  c->cd(2);
  hstrip[imod]->Draw();

  c->cd(3);
  hclust[imod]->Draw();
}

void FitGaus_FWHM( TH1D *htest, double thresh=0.5 ){
  int binmax = htest->GetMaximumBin();
  int binlow = binmax, binhigh = binmax;

  double max = htest->GetBinContent(binmax);

  while( htest->GetBinContent(binlow) >= thresh*max && binlow > 1 ){binlow--;}
  while( htest->GetBinContent(binhigh) >= thresh*max && binhigh < htest->GetNbinsX() ){ binhigh++; }

  double xlow = htest->GetBinCenter(binlow);
  double xhigh = htest->GetBinCenter(binhigh);

  htest->Fit("gaus","qS","",xlow, xhigh);
}

void GetTrackingCuts( const char *rootfilename, const char *outfilename="GENtrackingcuts.root", const int nmodules=2 ){

  ROOT::EnableImplicitMT(10);
  //We want this macro to do several things:
  // 1. Set tracking constraints
  // 2. Set thresholds

  //TFile *fin = new TFile(rootfilename,"READ");

  TChain *T = new TChain("T");
  T->Add(rootfilename);

 
  
  //ROOT::RDataFrame F(*C);


  TCut global_cut =  "bb.tr.n>=1&&abs(bb.tr.vz[0])<0.27&&bb.gem.track.chi2ndf[0]<10&&bb.gem.track.nhits[0]>3&&bb.ps.e>0.2&&abs(bb.etot_over_p-1.)<0.25";
  
  //TCut global_cut = "bb.tr.n==1&&abs(bb.tr.vz)<0.08&&bb.gem.track.chi2ndf<10&&bb.gem.track.nhits>3&&bb.ps.e>0.25&&abs(bb.etot_over_p-1.)<0.25";

  TEventList *elist = new TEventList("elist");

  T->Draw(">>elist",global_cut);

  T->SetEventList(elist);

  cout << "number of events passing global cut = " << elist->GetN() << endl;
  
  TFile *fout = new TFile( outfilename, "RECREATE" );
  
  TH1D *hdxfcp = new TH1D("hdxfcp", "Front Constraint X;x_{track}-x_{fcp} (m);", 300, -0.75, 0.75 );
  TH1D *hdyfcp = new TH1D("hdyfcp", "Front Constraint Y;y_{track}-y_{fcp} (m);", 300, -0.30, 0.30 );
  TH1D *hdxbcp = new TH1D("hdxbcp", "Back Constraint X;x_{track}+x'_{track}z_{bcp}-x_{bcp} (m);", 300, -0.30, 0.50 );
  TH1D *hdybcp = new TH1D("hdybcp", "Back Constrant Y;y_{track}+y'_{track}z_{bcp}-y_{bcp} (m);", 300, -0.30, 0.30 );

  TH2D *hdxdyfcp = new TH2D("hdxdyfcp", ";y_{track}-y_{fcp} (m);x_{track}-x_{fcp} (m)", 150, -0.15, 0.15, 150, -0.15, 0.15 );
  TH2D *hdxdybcp = new TH2D("hdxdybcp", ";y_{track}+y'_{track}z_{bcp}-y_{bcp} (m);x_{track}+x'_{track}z_{bcp}-x_{bcp} (m)", 150, -0.15, 0.15, 150, -0.15, 0.15 );

  TH1D *hdthcp = new TH1D("hdthcp", "Track #theta Constraint; x'_{track}-x'_{constraint};", 300, -0.1, 0.1 );
  TH1D *hdphcp = new TH1D("hdphcp", "Track #phi Constraint; y'_{track}-y'_{constraint};", 300, -0.1, 0.1 );
  
  T->Project( "hdxbcp", "bb.tr.x+bb.tr.th*bb.z_bcp-bb.x_bcp" );
  T->Project( "hdybcp", "bb.tr.y+bb.tr.ph*bb.z_bcp-bb.y_bcp" );
  T->Project( "hdxdybcp", "bb.tr.x+bb.tr.th*bb.z_bcp-bb.x_bcp:bb.tr.y+bb.tr.ph*bb.z_bcp-bb.y_bcp" );
  T->Project( "hdxfcp", "bb.tr.x+bb.tr.th*bb.z_fcp-bb.x_fcp" );
  T->Project( "hdyfcp", "bb.tr.y+bb.tr.ph*bb.z_fcp-bb.y_fcp" );
  T->Project( "hdxdyfcp", "bb.tr.x+bb.tr.th*bb.z_fcp-bb.x_fcp:bb.tr.y+bb.tr.ph*bb.z_fcp-bb.y_fcp" );
  T->Project( "hdthcp", "bb.tr.th-(bb.x_bcp-bb.x_fcp)/(bb.z_bcp-bb.z_fcp)" );
  T->Project( "hdphcp", "bb.tr.ph-(bb.y_bcp-bb.y_fcp)/(bb.z_bcp-bb.z_fcp)" );

  TString fname_db = outfilename;
  fname_db.ReplaceAll(".root",".dat");

  ofstream dbfile(fname_db.Data());


  fname_db.ReplaceAll(".dat", "_timecuts.dat");
  ofstream dbfile_tcut(fname_db.Data() );
  
  dbfile << "bb.frontconstraint_x0 = " << hdxfcp->GetMean() << endl;
  dbfile << "bb.frontconstraint_y0 = " << hdyfcp->GetMean() << endl;
  dbfile << "bb.backconstraint_x0 = " << hdxbcp->GetMean() << endl;
  dbfile << "bb.backconstraint_y0 = " << hdybcp->GetMean() << endl << endl;

  dbfile << "bb.frontconstraintwidth_x = " << hdxfcp->GetRMS() * 4.5 << endl;
  dbfile << "bb.frontconstraintwidth_y = " << hdyfcp->GetRMS() * 4.5 << endl;
  dbfile << "bb.backconstraintwidth_x = " << hdxbcp->GetRMS() * 4.5 << endl;
  dbfile << "bb.backconstraintwidth_y = " << hdybcp->GetRMS() * 4.5 << endl << endl;

  dbfile << "bb.gem.constraintwidth_theta = " << hdthcp->GetRMS() * 4.5 << endl;
  dbfile << "bb.gem.constraintwidth_phi = " << hdphcp->GetRMS() * 4.5 << endl;

  

  threshsamp_module.resize(nmodules);
  threshstrip_module.resize(nmodules);
  threshclust_module.resize(nmodules);
  hsamp.resize(nmodules);
  hstrip.resize(nmodules);
  hclust.resize(nmodules);
  htmaxU.resize(nmodules);
  htmaxV.resize(nmodules);

  
  dbfile << endl;

  
  
  
  //clustercut += ;
  
  //get individual sample width for each module and check it against track-based sample threshold; apply 5sigma or track-based, whichever is higher:
  for( int imod=0; imod<nmodules; imod++ ){
    TCut clustercut = (TString::Format("bb.gem.hit.nstripu>1&&bb.gem.hit.nstripv>1&&bb.gem.hit.module==%d",imod)).Data();
    
    TString histname;
    histname.Form("hADCpedsubU_allstrips_bb_gem_m%d", imod );

    TFile *fin = T->GetFile();
    
    TH1D *htemp;
    fin->GetObject( histname.Data(), htemp );


    int binmax, binlow, binhigh;
    double max;
    double sigmaU = 20.0, sigmaV = 20.0;
    
    if( htemp ){
    
      binmax = htemp->GetMaximumBin();
      max = htemp->GetBinContent(binmax);
      
      binlow=binmax;
      binhigh=binmax;
      
      while( htemp->GetBinContent(binlow) >= 0.4*max && binlow > 1 ){binlow--;}
      while( htemp->GetBinContent(binhigh) >= 0.4*max && binhigh < htemp->GetNbinsX() ){binhigh++; }
      
      
      
      if( htemp->GetEntries() >= 100 ){
	htemp->Fit("gaus", "q", "", htemp->GetBinCenter(binlow), htemp->GetBinCenter(binhigh) );
	
	sigmaU = ( ( (TF1*) htemp->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Sigma");
      }
    }
    
    histname.Form("hADCpedsubV_allstrips_bb_gem_m%d", imod );
      
    TH1D *htempV;
    fin->GetObject( histname.Data(), htempV );
    
    if( htempV ){
      
      binmax = htempV->GetMaximumBin();

      max = htempV->GetBinContent(binmax);
     
      binlow=binmax;
      binhigh=binmax;
      
      while( htempV->GetBinContent(binlow) >= 0.4*max && binlow > 1 ){binlow--;}
      while( htempV->GetBinContent(binhigh) >= 0.4*max && binhigh < htempV->GetNbinsX() ){binhigh++; }
      
      if( htempV->GetEntries() >= 100 ){
	htempV->Fit("gaus", "q", "", htempV->GetBinCenter(binlow), htempV->GetBinCenter(binhigh) );
	
	sigmaV = ( (TF1*) (htempV->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Sigma");
      }
    }
    
    double sigma = 0.5*(sigmaU+sigmaV);

    threshsamp_module[imod] =  5.0*sigma;
      
    TString hsampname,hstripname,hclustname;
    hsampname.Form("hADCmaxsamp_m%d",imod);
    hstripname.Form("hADCmaxstrip_m%d",imod);
    hclustname.Form("hADCclust_m%d",imod);
    
    //fout->cd();
    
    TH1D *hADCmaxsamp = new TH1D(hsampname.Data(),Form("Module %i Sample ADC;max ADC sample;",imod), 1500, 0, 3000 );
    TH1D *hADCmaxstrip = new TH1D(hstripname.Data(), Form("Module %i Strip Sum ADC;strip ADC sum;",imod), 1500, 0, 15000 );
    TH1D *hADCclust = new TH1D(hclustname.Data(), Form("Module %i Cluster Sum ADC;cluster ADC sum;",imod), 1500, 0, 30000 );
    
    T->Project( hsampname.Data(), "0.5*(bb.gem.hit.ADCmaxsampU+bb.gem.hit.ADCmaxsampV)", clustercut );

    if( hADCmaxsamp->GetEntries() < 300 ) continue;
    
    double fraction = 0.003;

    int ibin=1;

    ibin=1;
    while( hADCmaxsamp->Integral(1,ibin) < fraction * hADCmaxsamp->GetEntries() ){ ibin++; }
    double threshsamp_avg = hADCmaxsamp->GetBinCenter( ibin );

    threshsamp_module[imod] = std::max( 3.5*sigma, std::min( threshsamp_avg, 5.0*sigma ) );
    
    TString samplecut;
    samplecut.Form("0.5*(bb.gem.hit.ADCmaxsampU+bb.gem.hit.ADCmaxsampV)>%g",threshsamp_module[imod]);

    //cout << "\"" << samplecut << "\"" << endl;
    
    TCut modulecut = clustercut;
    modulecut += samplecut;

    //cout << "\"" << modulecut << "\"" << endl;
    
					
    T->Project( hstripname.Data(), "0.5*(bb.gem.hit.ADCmaxstripU+bb.gem.hit.ADCmaxstripV)", modulecut );
    T->Project( hclustname.Data(), "bb.gem.hit.ADCavg", modulecut );

    //hADCmaxsamp->Print();
    //hADCmaxstrip->Print();
    //hADCclust->Print();
    
    hsamp[imod] = hADCmaxsamp;
    hstrip[imod] = hADCmaxstrip;
    hclust[imod] = hADCclust;


    fraction = 0.002;
    
    ibin=1;
    while( hADCmaxstrip->Integral(1,ibin) < fraction * hADCmaxstrip->GetEntries() ){ ibin++; }
    double threshstrip_avg = hADCmaxstrip->GetBinCenter( ibin );
    
    ibin=1;
    while( hADCclust->Integral(1,ibin) < fraction * hADCclust->GetEntries() ){ ibin++; }
    
    double threshclust_avg = hADCclust->GetBinCenter( ibin );
    
    threshstrip_module[imod] = threshstrip_avg;
    threshclust_module[imod] = threshclust_avg;

    TString dbline;
    dbfile << dbline.Format( "bb.gem.m%d.threshold_sample = %g", imod, threshsamp_module[imod] ) << endl;
    dbfile << dbline.Format( "bb.gem.m%d.threshold_stripsum = %g", imod, threshstrip_module[imod] ) << endl;
    dbfile << dbline.Format( "bb.gem.m%d.threshold_clustersum = %g", imod, threshclust_module[imod] ) << endl << endl;

    TString hname;

    TCut hitqualitycut = (hname.Format("bb.gem.hit.module==%d&&bb.gem.track.nhits[bb.gem.hit.trackindex]>3&&bb.gem.track.chi2ndf[bb.gem.hit.trackindex]<10.&&abs(bb.tr.vz[bb.gem.hit.trackindex])<0.1",imod)).Data();
    
    //Let's do the following now:
    TH1D *hdeltat = new TH1D( hname.Format( "hdeltat_mod%d", imod ), "Standard hit time U-V ;#Deltat (ns);", 200,-50,50);
    TH1D *hdeltat_deconv = new TH1D( hname.Format( "hdeltat_deconv_mod%d", imod ), "Deconvoluted hit time U-V; #Delta t (ns);", 200, -100, 100 );
    TH1D *hdeltat_fit = new TH1D( hname.Format( "hdeltat_fit_mod%d", imod ), "Fit hit time U - V; #Delta t (ns);", 200, -50., 50. );

    TH1D *htavg = new TH1D( hname.Format("htavg_mod%d", imod ), "Standard hit time (U+V)/2; t_{avg} (ns);", 150, 0.0, 150.0 );
    TH1D *htdeconv = new TH1D( hname.Format("htdeconv_mod%d", imod ), "Deconvoluted hit time (U+V)/2; t_{avg} (ns);", 150, -100, 200 );
    TH1D *htfit = new TH1D( hname.Format("htfit_mod%d", imod ), "Fit hit time (U+V)/2; t_{avg} (ns);",150,-100, 200 );

    TH1D *hADCratio = new TH1D( hname.Format( "hADCratio_mod%d", imod ), "ADC ratio; ADCV/ADCU;", 200, 0, 2 );
    TH1D *hADCratio_deconv = new TH1D( hname.Format( "hADCratio_deconv_mod%d", imod ), "ADC ratio deconvoluted; ADCV/ADCU;", 200, 0, 2 );

    TH1D *hADCasym = new TH1D( hname.Format( "hADCasym_mod%d", imod ), "ADC asymmetry;", 200, -1., 1. );
    TH1D *hADCasym_deconv = new TH1D( hname.Format( "hADCasym_deconv_mod%d", imod ), "Deconv ADC asymmetry;", 200, -1., 1. );

    T->Project( hADCasym->GetName(), "bb.gem.hit.ADCasym", hitqualitycut );
    T->Project( hADCasym_deconv->GetName(), "bb.gem.hit.ADCasym_deconv", hitqualitycut );

    FitGaus_FWHM( hADCasym, 0.3 );

    double sigmaASYM = ( (TF1*) (hADCasym->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Sigma");

    
    
    T->Project( hdeltat->GetName(), "bb.gem.hit.deltat", hitqualitycut );
    T->Project( hdeltat_deconv->GetName(), "bb.gem.hit.deltat_deconv", hitqualitycut );
    T->Project( hdeltat_fit->GetName(), "bb.gem.hit.deltat_fit", hitqualitycut );

    hdeltat->Fit( "gaus", "qS", "", hdeltat->GetMean()-hdeltat->GetRMS(), hdeltat->GetMean() + hdeltat->GetRMS() );
    hdeltat_deconv->Fit( "gaus", "qS", "", hdeltat_deconv->GetMean()-hdeltat_deconv->GetRMS(), hdeltat_deconv->GetMean() + hdeltat_deconv->GetRMS() );
    hdeltat_fit->Fit( "gaus", "qS", "", hdeltat_fit->GetMean()-hdeltat_fit->GetRMS(), hdeltat_fit->GetMean() + hdeltat_fit->GetRMS() );

    FitGaus_FWHM( hdeltat, 0.4 );
    FitGaus_FWHM( hdeltat_deconv, 0.4 );
    FitGaus_FWHM( hdeltat_fit, 0.4 );
    
    T->Project( htavg->GetName(), "0.5*(bb.gem.hit.Utime+bb.gem.hit.Vtime)", hitqualitycut );
    T->Project( htdeconv->GetName(), "0.5*(bb.gem.hit.UtimeDeconv+bb.gem.hit.VtimeDeconv)", hitqualitycut );
    T->Project( htfit->GetName(), "0.5*(bb.gem.hit.UtimeFit+bb.gem.hit.VtimeFit)", hitqualitycut );

    

    FitGaus_FWHM( htavg, 0.25 );
    FitGaus_FWHM( htdeconv, 0.25 );
    FitGaus_FWHM( htfit, 0.25 );
    
    T->Project( hADCratio->GetName(), "bb.gem.hit.ADCV/bb.gem.hit.ADCU", hitqualitycut );
    T->Project( hADCratio_deconv->GetName(), "bb.gem.hit.ADCV_deconv/bb.gem.hit.ADCU_deconv", hitqualitycut );

    FitGaus_FWHM( hADCratio, 0.4 );
    FitGaus_FWHM( hADCratio_deconv, 0.4 );
    
    TH1D *hmaxstrip_tavg = new TH1D(hname.Format("hmaxstrip_tavg_mod%d",imod), "Max. strip time (U+V)/2; t_{avg} (ns);", 150,0,150);
    TH1D *hmaxstrip_tavg_deconv = new TH1D(hname.Format("hmaxstrip_tavg_deconv_mod%d", imod ), "Max. strip Deconv time (U+V)/2; t_{avg} (ns);", 150, -100, 200 );
    TH1D *hmaxstrip_tavg_fit = new TH1D(hname.Format("hmaxstrip_tavg_fit_mod%d",imod), "Max. strip Fit time (U+V)/2; t_{avg} (ns);", 150, -100, 200 );

    TH1D *hmaxstrip_dt = new TH1D(hname.Format("hmaxstrip_dt_mod%d",imod), "Max. strip time (U-V); #Delta t (ns);", 100, -50, 50 );
    TH1D *hmaxstrip_dtdeconv = new TH1D(hname.Format("hmaxstrip_dt_deconv_mod%d", imod ), "Max. strip deconv. time (U-V); #Delta t (ns);", 150,-75,75);
    TH1D *hmaxstrip_dtfit = new TH1D(hname.Format("hmaxstrip_dt_fit_mod%d", imod ), "Max strip fit time (U-V); #Delta t (ns);", 100, -50, 50 );
    
    T->Project( hmaxstrip_tavg->GetName(), "0.5*(bb.gem.hit.UtimeMaxStrip+bb.gem.hit.VtimeMaxStrip)", hitqualitycut );
    T->Project( hmaxstrip_tavg_deconv->GetName(), "0.5*(bb.gem.hit.UtimeMaxStripDeconv+bb.gem.hit.VtimeMaxStripDeconv)", hitqualitycut );
    T->Project( hmaxstrip_tavg_fit->GetName(), "0.5*(bb.gem.hit.UtimeMaxStripFit+bb.gem.hit.VtimeMaxStripFit)", hitqualitycut );

    T->Project( hmaxstrip_dt->GetName(), "bb.gem.hit.UtimeMaxStrip-bb.gem.hit.VtimeMaxStrip", hitqualitycut );
    T->Project( hmaxstrip_dtdeconv->GetName(), "bb.gem.hit.UtimeMaxStripDeconv-bb.gem.hit.VtimeMaxStripDeconv", hitqualitycut );
    T->Project( hmaxstrip_dtfit->GetName(), "bb.gem.hit.UtimeMaxStripFit-bb.gem.hit.VtimeMaxStripFit", hitqualitycut );
    
    FitGaus_FWHM( hmaxstrip_tavg, 0.25 );
    FitGaus_FWHM( hmaxstrip_tavg_deconv, 0.25 );
    FitGaus_FWHM( hmaxstrip_tavg_fit, 0.5 );

    FitGaus_FWHM( hmaxstrip_dt, 0.25 );
    FitGaus_FWHM( hmaxstrip_dtdeconv, 0.25 );
    FitGaus_FWHM( hmaxstrip_dtfit, 0.25 );

    TH1D *hmaxstrip_tU = new TH1D(hname.Format( "hmaxstrip_tU_mod%d", imod ), "Max. U strip time;t (ns);", 150,0,150 );
    TH1D *hmaxstrip_tV = new TH1D(hname.Format( "hmaxstrip_tV_mod%d", imod ), "Max. V strip time;t (ns);", 150,0,150 );

    TH1D *hmaxstrip_tU_deconv = new TH1D(hname.Format("hmaxstrip_tU_deconv_mod%d", imod), "Max. U strip time; t_{deconv} (ns);", 150, -100, 200 );
    TH1D *hmaxstrip_tV_deconv = new TH1D(hname.Format("hmaxstrip_tV_deconv_mod%d", imod), "Max. V strip time; t_{deconv} (ns);", 150, -100, 200 );

    TH1D *hmaxstrip_tU_fit = new TH1D(hname.Format("hmaxstrip_tU_fit_mod%d", imod), "Max. U strip time; t_{fit} (ns);", 150, -100, 200 );
    TH1D *hmaxstrip_tV_fit = new TH1D(hname.Format("hmaxstrip_tV_fit_mod%d", imod), "Max. V strip time; t_{fit} (ns);", 150, -100, 200 );

    T->Project( hmaxstrip_tU->GetName(), "bb.gem.hit.UtimeMaxStrip", hitqualitycut );
    T->Project( hmaxstrip_tV->GetName(), "bb.gem.hit.VtimeMaxStrip", hitqualitycut );
    
    T->Project( hmaxstrip_tU_deconv->GetName(), "bb.gem.hit.UtimeMaxStripDeconv", hitqualitycut );
    T->Project( hmaxstrip_tV_deconv->GetName(), "bb.gem.hit.VtimeMaxStripDeconv", hitqualitycut );
    T->Project( hmaxstrip_tU_fit->GetName(), "bb.gem.hit.UtimeMaxStripFit", hitqualitycut );
    T->Project( hmaxstrip_tV_fit->GetName(), "bb.gem.hit.VtimeMaxStripFit", hitqualitycut );

    FitGaus_FWHM( hmaxstrip_tU, 0.25 );
    FitGaus_FWHM( hmaxstrip_tV, 0.25 );
    FitGaus_FWHM( hmaxstrip_tU_deconv, 0.25 );
    FitGaus_FWHM( hmaxstrip_tV_deconv, 0.25 );
    FitGaus_FWHM( hmaxstrip_tU_fit, 0.25 );
    FitGaus_FWHM( hmaxstrip_tV_fit, 0.25 );
    
    double tmeanU = ( (TF1*) (hmaxstrip_tU->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Mean");
    double tsigmaU = ( (TF1*) (hmaxstrip_tU->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Sigma");

    double tmeanV = ( (TF1*) (hmaxstrip_tV->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Mean");
    double tsigmaV = ( (TF1*) (hmaxstrip_tV->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Sigma");
    //TString dbline;
    //dbfile_tcut << endl;
    dbfile_tcut << dbline.Format( "bb.gem.m%d.maxstrip_t0 = %9.4g %9.4g", imod, tmeanU, tmeanV ) << endl;
    dbfile_tcut << dbline.Format( "bb.gem.m%d.maxstrip_tsigma = %9.4g %9.4g", imod, tsigmaU, tsigmaV ) << endl;
    dbfile_tcut << dbline.Format( "bb.gem.m%d.maxstrip_tcut = %9.4g %9.4g", imod, 4.5, 4.5 ) << endl << endl;

    double tmeanU_deconv = ( (TF1*) (hmaxstrip_tU_deconv->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Mean");
    double tsigmaU_deconv = ( (TF1*) (hmaxstrip_tU_deconv->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Sigma");

    double tmeanV_deconv = ( (TF1*) (hmaxstrip_tV_deconv->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Mean");
    double tsigmaV_deconv = ( (TF1*) (hmaxstrip_tV_deconv->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Sigma");

    //    dbfile_tcut << endl;
    dbfile_tcut << dbline.Format( "bb.gem.m%d.maxstrip_t0_deconv = %9.4g %9.4g", imod, tmeanU_deconv, tmeanV_deconv ) << endl;
    dbfile_tcut << dbline.Format( "bb.gem.m%d.maxstrip_tsigma_deconv = %9.4g %9.4g", imod, tsigmaU_deconv, tsigmaV_deconv ) << endl;
    dbfile_tcut << dbline.Format( "bb.gem.m%d.maxstrip_tcut_deconv = %9.4g %9.4g", imod, 4.5, 4.5 ) << endl << endl;

    double tmeanU_fit = ( (TF1*) (hmaxstrip_tU_fit->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Mean");
    double tsigmaU_fit = ( (TF1*) (hmaxstrip_tU_fit->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Sigma");

    double tmeanV_fit = ( (TF1*) (hmaxstrip_tV_fit->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Mean");
    double tsigmaV_fit = ( (TF1*) (hmaxstrip_tV_fit->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Sigma");

    //dbfile_tcut << endl;
    dbfile_tcut << dbline.Format( "bb.gem.m%d.maxstrip_t0_fit = %9.4g %9.4g", imod, tmeanU_fit, tmeanV_fit ) << endl;
    dbfile_tcut << dbline.Format( "bb.gem.m%d.maxstrip_tsigma_fit = %9.4g %9.4g", imod, tsigmaU_fit, tsigmaV_fit ) << endl;
    dbfile_tcut << dbline.Format( "bb.gem.m%d.maxstrip_tcut_fit = %9.4g %9.4g", imod, 4.5, 4.5 ) << endl << endl;


    TH1D *htUhit = new TH1D(hname.Format( "htUhit_mod%d", imod ), "U hit time; t (ns);", 150, 0, 150 );
    TH1D *htVhit = new TH1D(hname.Format( "htVhit_mod%d", imod ), "V hit time; t (ns);", 150, 0, 150 );
    TH1D *htUhitDeconv = new TH1D(hname.Format( "htUhitDeconv_mod%d", imod ), "Deconv U hit time; t (ns);", 150, -100, 200 );
    TH1D *htVhitDeconv = new TH1D(hname.Format( "htVhitDeconv_mod%d", imod ), "Deconv V hit time; t (ns);", 150, -100, 200 );
    TH1D *htUhitFit = new TH1D(hname.Format( "htUhitFit_mod%d", imod ), "Fit U hit time; t (ns);", 150, -100, 200 );
    TH1D *htVhitFit = new TH1D(hname.Format( "htVhitFit_mod%d", imod ), "Fit V hit time; t (ns);", 150, -100, 200 );

    //Unfortunately we also apparently need mean cluster times separately for U and V strips as well: 
    T->Project( htUhit->GetName(), "bb.gem.hit.Utime", hitqualitycut );
    T->Project( htVhit->GetName(), "bb.gem.hit.Vtime", hitqualitycut );
    T->Project( htUhitDeconv->GetName(), "bb.gem.hit.UtimeDeconv", hitqualitycut );
    T->Project( htVhitDeconv->GetName(), "bb.gem.hit.VtimeDeconv", hitqualitycut );
    T->Project( htUhitFit->GetName(), "bb.gem.hit.UtimeFit", hitqualitycut );
    T->Project( htVhitFit->GetName(), "bb.gem.hit.VtimeFit", hitqualitycut );

    FitGaus_FWHM( htUhit, 0.25 );
    FitGaus_FWHM( htVhit, 0.25 );
    FitGaus_FWHM( htUhitDeconv, 0.25 );
    FitGaus_FWHM( htVhitDeconv, 0.25 );
    FitGaus_FWHM( htUhitFit, 0.25 );
    FitGaus_FWHM( htVhitFit, 0.25 );

    tmeanU = ( (TF1*) (htUhit->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Mean");
    tmeanV = ( (TF1*) (htVhit->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Mean");
    tsigmaU = ( (TF1*) (htUhit->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Sigma");
    tsigmaV = ( (TF1*) (htVhit->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Sigma");

    dbfile_tcut << dbline.Format( "bb.gem.m%d.HitTimeMean = %10.5g %10.5g", imod, tmeanU, tmeanV ) << endl;
    dbfile_tcut << dbline.Format( "bb.gem.m%d.HitTimeSigma = %10.5g %10.5g", imod, tsigmaU, tsigmaV ) << endl;

    tmeanU = ( (TF1*) (htUhitDeconv->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Mean");
    tmeanV = ( (TF1*) (htVhitDeconv->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Mean");
    tsigmaU = ( (TF1*) (htUhitDeconv->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Sigma");
    tsigmaV = ( (TF1*) (htVhitDeconv->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Sigma");

    dbfile_tcut << dbline.Format( "bb.gem.m%d.HitTimeMeanDeconv = %10.5g %10.5g", imod, tmeanU, tmeanV ) << endl;
    dbfile_tcut << dbline.Format( "bb.gem.m%d.HitTimeSigmaDeconv = %10.5g %10.5g", imod, tsigmaU, tsigmaV ) << endl;
    
    tmeanU = ( (TF1*) (htUhitFit->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Mean");
    tmeanV = ( (TF1*) (htVhitFit->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Mean");
    tsigmaU = ( (TF1*) (htUhitFit->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Sigma");
    tsigmaV = ( (TF1*) (htVhitFit->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Sigma");

    dbfile_tcut << dbline.Format( "bb.gem.m%d.HitTimeMeanFit = %10.5g %10.5g", imod, tmeanU, tmeanV ) << endl;
    dbfile_tcut << dbline.Format( "bb.gem.m%d.HitTimeSigmaFit = %10.5g %10.5g", imod, tsigmaU, tsigmaV ) << endl << endl;

    double dtsigma = ( (TF1*) ( hdeltat->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Sigma");
    double dtsigma_deconv = ( (TF1*) ( hdeltat_deconv->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Sigma");
    double dtsigma_fit = ( (TF1*) ( hdeltat_fit->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Sigma");

    dbfile_tcut << dbline.Format("bb.gem.m%d.deltat_sigma = %10.5g", imod, dtsigma ) << endl;
    dbfile_tcut << dbline.Format("bb.gem.m%d.deltat_sigma_deconv = %10.5g", imod, dtsigma_deconv ) << endl;
    dbfile_tcut << dbline.Format("bb.gem.m%d.deltat_sigma_fit = %10.5g", imod, dtsigma_fit ) << endl << endl;

    dbfile_tcut << dbline.Format("bb.gem.m%d.deltat_cut = %10.5g", imod, dtsigma*3.5 ) << endl;
    dbfile_tcut << dbline.Format("bb.gem.m%d.deltat_cut_deconv = %10.5g", imod, dtsigma_deconv*3.5 ) << endl;
    dbfile_tcut << dbline.Format("bb.gem.m%d.deltat_cut_fit = %10.5g", imod, dtsigma_fit*3.5 ) << endl << endl;

    double ADCratio_sigma = ( (TF1*) (hADCratio->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Sigma");
    double ADCratio_sigma_deconv = ( (TF1*) (hADCratio_deconv->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Sigma");

    dbfile_tcut << dbline.Format("bb.gem.m%d.ADCratio_sigma = %10.5g", imod, ADCratio_sigma ) << endl;
    dbfile_tcut << dbline.Format( "bb.gem.m%d.ADCasym_sigma = %10.5g", imod, sigmaASYM ) << endl << endl;
    dbfile_tcut << dbline.Format("bb.gem.m%d.ADCratio_sigma_deconv = %10.5g", imod, ADCratio_sigma_deconv ) << endl << endl;


    htmaxU = hmaxstrip_tU
    
    
  }
  
  fout->Write();
  

  ///// Here we pring out all the pdfs for the calorimeter and thresholds /////
  TString outfilepdf = outfilename;
  outfilepdf.ReplaceAll(".root",".pdf");

  TCanvas *c1 = new TCanvas("c1","",1600,1200);
  c1->Divide(2,3);

  c1->cd(1);
  DrawHist(hdxfcp);

  c1->cd(2);
  DrawHist(hdyfcp);

  c1->cd(3);
  DrawHist(hdxbcp);
  
  c1->cd(4);
  DrawHist(hdybcp);

  c1->cd(5);
  DrawHist(hdthcp);
  
  c1->cd(6);
  DrawHist(hdphcp);

  c1->Print(outfilepdf + "(");
  
  for(int imod=0; imod < nmodules; imod++){
    TCanvas *c= new TCanvas("c","",1200,1200);
    c->Divide(1,3);
    
    c->cd(1);
    DrawHist(hsamp[imod],1,threshsamp_module[imod]);
    
    c->cd(2);
    DrawHist(hstrip[imod],1,threshstrip_module[imod]);

    c->cd(3);
    DrawHist(hclust[imod],1,threshclust_module[imod]);

    if(imod == nmodules - 1) c->Print(outfilepdf + ")");
    else c->Print(outfilepdf);
  }
  
  ///// Here we print out all the plots for the timing cuts ///////
  TString outfilepdf_time = outfilename;
  outfilepdf_time.ReplaceAll(".root","_timecuts.pdf");

  TCanvas *c2 = new TCanvas("c2","",1600,1200);
  

}
