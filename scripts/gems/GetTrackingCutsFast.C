#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TCut.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TTreeFormula.h"
#include <iostream>
#include <fstream>
#include "TString.h"
#include "TClonesArray.h"


void DrawHist(TH1D *hist, vector<double> value, int flag = 0){

  //flag = 0 for calorimeter plots
  //flag = 1 for ADC plots
  //flag = 2 for time plots

  hist->Draw();

  TPaveText *pt;
  if(flag == 0){
    pt = new TPaveText(0.12,0.75,0.35,0.88,"ndc");
    pt->AddText(Form("Center = %g",value[0]));
    pt->AddText(Form("Width = %g",value[1]));
  }
  else if(flag == 1){
    pt = new TPaveText(0.55,0.75,0.75,0.88,"ndc");
    pt->AddText(Form("Threshold = %g",value[0]));
  }
  else if(flag == 2){
    pt = new TPaveText(0.12,0.75,0.35,0.88,"ndc");
    if(value.size() == 2){
      pt->AddText(Form("Mean = %g",value[0]));
      pt->AddText(Form("#sigma = %g",value[1]));
    }
    if(value.size() == 1){
      pt->AddText(Form("Cut = %g",value[0]));
    }
 }
  pt->SetFillColor(0);
  pt->Draw("same");


}

void FitGaus_FWHM( TH1D *htest, double thresh=0.5 ){
  int binmax = htest->GetMaximumBin();
  int binlow = binmax, binhigh = binmax;

  double max = htest->GetBinContent(binmax);

  while( htest->GetBinContent(binlow) >= thresh*max && binlow > 1 ){binlow--;}
  while( htest->GetBinContent(binhigh) >= thresh*max && binhigh < htest->GetNbinsX() ){ binhigh++; }

  double xlow = htest->GetBinCenter(binlow);
  double xhigh = htest->GetBinCenter(binhigh);

  htest->Fit("gaus","q0S","",xlow, xhigh);
}

void GetTrackingCutsFast( const char *configfilename, const char *outfilename="GENtrackingcuts.root", int nmodules=8 ){

  ifstream infile(configfilename);
  
  TChain *C = new TChain("T");

  TString currentline;
  while( infile >> currentline && !currentline.BeginsWith("endlist") ){
    if( !currentline.BeginsWith("#") ){
      C->Add(currentline.Data());
    }
  }

  TCut globalcut = "";
  
  while( infile >> currentline && !currentline.BeginsWith("endcut") ){
    if( !currentline.BeginsWith("#") ){
      globalcut += currentline.Data();
    }
  }

  //Initialize branches:
  
  //what branches do we need? Tracking stuff:
  //int ntracks=0;
  const int MAXNTRACKS=100;
  //int nhits=0;
  const int MAXNHITS = 1000;
  const int MAXNCLUST = 1;


  C->SetBranchStatus("*",0);
  
  double EPS, ESH;
  C->SetBranchStatus("bb.etot_over_p",1);
  C->SetBranchStatus("bb.ps.e",1);
  C->SetBranchStatus("bb.sh.e",1);
  C->SetBranchAddress("bb.ps.e",&EPS);
  C->SetBranchAddress("bb.sh.e",&ESH);

  double ntracks;
  double tracknhits[MAXNTRACKS];
  double trackngoodhits[MAXNTRACKS];
  double trackchi2ndf[MAXNTRACKS];
  double trackchi2ndf_hitquality[MAXNTRACKS];
  double trackt0[MAXNTRACKS];
  double xfp[MAXNTRACKS], yfp[MAXNTRACKS], thfp[MAXNTRACKS], phfp[MAXNTRACKS];
  double xtar[MAXNTRACKS], ytar[MAXNTRACKS], thtar[MAXNTRACKS], phtar[MAXNTRACKS];
  double rxfp[MAXNTRACKS], ryfp[MAXNTRACKS], rthfp[MAXNTRACKS], rphfp[MAXNTRACKS];
  double dxfp[MAXNTRACKS], dyfp[MAXNTRACKS], dthfp[MAXNTRACKS], dphfp[MAXNTRACKS];
  double p[MAXNTRACKS],px[MAXNTRACKS],py[MAXNTRACKS],pz[MAXNTRACKS],vz[MAXNTRACKS]; 
  double xfcp[MAXNCLUST],yfcp[MAXNCLUST],zfcp[MAXNCLUST],xbcp[MAXNCLUST],ybcp[MAXNCLUST],zbcp[MAXNCLUST];

  C->SetBranchStatus("bb.gem.track.nhits",1);
  C->SetBranchStatus("bb.gem.track.ngoodhits",1);
  C->SetBranchStatus("bb.gem.track.chi2ndf",1);
  C->SetBranchStatus("bb.gem.track.chi2ndf_hitquality",1);
  C->SetBranchStatus("bb.gem.track.t0",1);
  C->SetBranchStatus("bb.tr.x",1);
  C->SetBranchStatus("bb.tr.y",1);
  C->SetBranchStatus("bb.tr.th",1);
  C->SetBranchStatus("bb.tr.ph",1);
  C->SetBranchStatus("bb.tr.tg_x",1);
  C->SetBranchStatus("bb.tr.tg_y",1);
  C->SetBranchStatus("bb.tr.tg_th",1);
  C->SetBranchStatus("bb.tr.tg_ph",1);
  C->SetBranchStatus("bb.tr.r_x",1);
  C->SetBranchStatus("bb.tr.r_y",1);
  C->SetBranchStatus("bb.tr.r_th",1);
  C->SetBranchStatus("bb.tr.r_ph",1);
  C->SetBranchStatus("bb.tr.d_x",1);
  C->SetBranchStatus("bb.tr.d_y",1);
  C->SetBranchStatus("bb.tr.d_th",1);
  C->SetBranchStatus("bb.tr.d_ph",1);
  C->SetBranchStatus("bb.tr.n",1);
  C->SetBranchStatus("bb.tr.p",1);
  C->SetBranchStatus("bb.tr.px",1);
  C->SetBranchStatus("bb.tr.py",1);
  C->SetBranchStatus("bb.tr.pz",1);
  C->SetBranchStatus("bb.tr.vz",1);
  C->SetBranchStatus("bb.x_fcp",1);
  C->SetBranchStatus("bb.y_fcp",1);
  C->SetBranchStatus("bb.z_fcp",1);
  C->SetBranchStatus("bb.x_bcp",1);
  C->SetBranchStatus("bb.y_bcp",1);
  C->SetBranchStatus("bb.z_bcp",1);
  
  C->SetBranchAddress("bb.gem.track.nhits",tracknhits);
  C->SetBranchAddress("bb.gem.track.ngoodhits",trackngoodhits);
  C->SetBranchAddress("bb.gem.track.chi2ndf",trackchi2ndf);
  C->SetBranchAddress("bb.gem.track.chi2ndf_hitquality",trackchi2ndf_hitquality);
  C->SetBranchAddress("bb.gem.track.t0",trackt0);
  C->SetBranchAddress("bb.tr.x",xfp);
  C->SetBranchAddress("bb.tr.y",yfp);
  C->SetBranchAddress("bb.tr.th",thfp);
  C->SetBranchAddress("bb.tr.ph",phfp);
  C->SetBranchAddress("bb.tr.tg_x",xtar);
  C->SetBranchAddress("bb.tr.tg_y",ytar);
  C->SetBranchAddress("bb.tr.tg_th",thtar);
  C->SetBranchAddress("bb.tr.tg_ph",phtar);
  C->SetBranchAddress("bb.tr.r_x",rxfp);
  C->SetBranchAddress("bb.tr.r_y",ryfp);
  C->SetBranchAddress("bb.tr.r_th",rthfp);
  C->SetBranchAddress("bb.tr.r_ph",rphfp);
  C->SetBranchAddress("bb.tr.d_x",dxfp);
  C->SetBranchAddress("bb.tr.d_y",dyfp);
  C->SetBranchAddress("bb.tr.d_th",dthfp);
  C->SetBranchAddress("bb.tr.d_ph",dphfp);
  C->SetBranchAddress("bb.tr.n",&ntracks);
  C->SetBranchAddress("bb.tr.p",p);
  C->SetBranchAddress("bb.tr.px",px);
  C->SetBranchAddress("bb.tr.py",py);
  C->SetBranchAddress("bb.tr.pz",pz);
  C->SetBranchAddress("bb.tr.vz",vz);
  C->SetBranchAddress("bb.x_fcp",xfcp);
  C->SetBranchAddress("bb.y_fcp",yfcp);
  C->SetBranchAddress("bb.z_fcp",zfcp);
  C->SetBranchAddress("bb.x_bcp",xbcp);
  C->SetBranchAddress("bb.y_bcp",ybcp);
  C->SetBranchAddress("bb.z_bcp",zbcp);

  double ngoodhits;
  //Track hit variables:
  double module[MAXNHITS];
  double nstripu[MAXNHITS], nstripv[MAXNHITS];
  double trackindex[MAXNHITS];

  C->SetBranchStatus("bb.gem.hit.trackindex",1);
  C->SetBranchStatus("bb.gem.hit.ngoodhits",1);
  C->SetBranchStatus("bb.gem.hit.module",1);
  C->SetBranchStatus("bb.gem.hit.nstripu",1);
  C->SetBranchStatus("bb.gem.hit.nstripv",1);

  C->SetBranchAddress("bb.gem.hit.ngoodhits",&ngoodhits);
  C->SetBranchAddress("bb.gem.hit.trackindex",trackindex);
  C->SetBranchAddress("bb.gem.hit.module",module);
  C->SetBranchAddress("bb.gem.hit.nstripu",nstripu);
  C->SetBranchAddress("bb.gem.hit.nstripv",nstripv);
  
  double ADCmaxsampU[MAXNHITS], ADCmaxsampV[MAXNHITS], ADCmaxstripU[MAXNHITS], ADCmaxstripV[MAXNHITS], ADCU[MAXNHITS], ADCV[MAXNHITS], ADCavg[MAXNHITS], ADCasym[MAXNHITS], ADCasymDeconv[MAXNHITS],DeconvADCU[MAXNHITS], DeconvADCV[MAXNHITS];

  C->SetBranchStatus("bb.gem.hit.ADCmaxsampU",1);
  C->SetBranchStatus("bb.gem.hit.ADCmaxsampV",1);
  C->SetBranchStatus("bb.gem.hit.ADCmaxstripU",1);
  C->SetBranchStatus("bb.gem.hit.ADCmaxstripV",1);
  C->SetBranchStatus("bb.gem.hit.ADCU",1);
  C->SetBranchStatus("bb.gem.hit.ADCV",1);
  C->SetBranchStatus("bb.gem.hit.ADCavg",1);
  C->SetBranchStatus("bb.gem.hit.ADCasym",1);
  C->SetBranchStatus("bb.gem.hit.ADCU_deconv",1);
  C->SetBranchStatus("bb.gem.hit.ADCV_deconv",1);
  C->SetBranchStatus("bb.gem.hit.ADCasym_deconv",1);

  C->SetBranchAddress("bb.gem.hit.ADCmaxsampU",ADCmaxsampU);
  C->SetBranchAddress("bb.gem.hit.ADCmaxsampV",ADCmaxsampV);
  C->SetBranchAddress("bb.gem.hit.ADCmaxstripU",ADCmaxstripU);
  C->SetBranchAddress("bb.gem.hit.ADCmaxstripV",ADCmaxstripV);
  C->SetBranchAddress("bb.gem.hit.ADCU",ADCU);
  C->SetBranchAddress("bb.gem.hit.ADCV",ADCV);
  C->SetBranchAddress("bb.gem.hit.ADCavg",ADCavg);
  C->SetBranchAddress("bb.gem.hit.ADCasym",ADCasym);
  C->SetBranchAddress("bb.gem.hit.ADCU_deconv",DeconvADCU);
  C->SetBranchAddress("bb.gem.hit.ADCV_deconv",DeconvADCV);
  C->SetBranchAddress("bb.gem.hit.ADCasym_deconv",ADCasymDeconv);

  double UtimeMaxStrip[MAXNHITS],VtimeMaxStrip[MAXNHITS];
  double UtimeMaxStripDeconv[MAXNHITS],VtimeMaxStripDeconv[MAXNHITS];
  double UtimeMaxStripFit[MAXNHITS],VtimeMaxStripFit[MAXNHITS];

  double Utime[MAXNHITS],Vtime[MAXNHITS];
  double UtimeDeconv[MAXNHITS],VtimeDeconv[MAXNHITS];
  double UtimeFit[MAXNHITS],VtimeFit[MAXNHITS];

  double HitTavg_corr[MAXNHITS];

  C->SetBranchStatus("bb.gem.hit.Tavg_corr",1);
  C->SetBranchAddress("bb.gem.hit.Tavg_corr",HitTavg_corr);
  
  C->SetBranchStatus("bb.gem.hit.UtimeMaxStrip",1);
  C->SetBranchStatus("bb.gem.hit.VtimeMaxStrip",1);
  C->SetBranchStatus("bb.gem.hit.UtimeMaxStripDeconv",1);
  C->SetBranchStatus("bb.gem.hit.VtimeMaxStripDeconv",1);
  C->SetBranchStatus("bb.gem.hit.UtimeMaxStripFit",1);
  C->SetBranchStatus("bb.gem.hit.VtimeMaxStripFit",1);

  C->SetBranchAddress("bb.gem.hit.UtimeMaxStrip",UtimeMaxStrip);
  C->SetBranchAddress("bb.gem.hit.VtimeMaxStrip",VtimeMaxStrip);
  C->SetBranchAddress("bb.gem.hit.UtimeMaxStripDeconv",UtimeMaxStripDeconv);
  C->SetBranchAddress("bb.gem.hit.VtimeMaxStripDeconv",VtimeMaxStripDeconv);
  C->SetBranchAddress("bb.gem.hit.UtimeMaxStripFit",UtimeMaxStripFit);
  C->SetBranchAddress("bb.gem.hit.VtimeMaxStripFit",VtimeMaxStripFit);

  C->SetBranchStatus("bb.gem.hit.Utime",1);
  C->SetBranchStatus("bb.gem.hit.Vtime",1);
  C->SetBranchStatus("bb.gem.hit.UtimeDeconv",1);
  C->SetBranchStatus("bb.gem.hit.VtimeDeconv",1);
  C->SetBranchStatus("bb.gem.hit.UtimeFit",1);
  C->SetBranchStatus("bb.gem.hit.VtimeFit",1);

  C->SetBranchAddress("bb.gem.hit.Utime",Utime);
  C->SetBranchAddress("bb.gem.hit.Vtime",Vtime);
  C->SetBranchAddress("bb.gem.hit.UtimeDeconv",UtimeDeconv);
  C->SetBranchAddress("bb.gem.hit.VtimeDeconv",VtimeDeconv);
  C->SetBranchAddress("bb.gem.hit.UtimeFit",UtimeFit);
  C->SetBranchAddress("bb.gem.hit.VtimeFit",VtimeFit);
  
  double deltat[MAXNHITS], deltatFit[MAXNHITS], deltatDeconv[MAXNHITS];

  C->SetBranchStatus("bb.gem.hit.deltat",1);
  C->SetBranchStatus("bb.gem.hit.deltat_deconv",1);
  C->SetBranchStatus("bb.gem.hit.deltat_fit",1);

  C->SetBranchAddress("bb.gem.hit.deltat",deltat);
  C->SetBranchAddress("bb.gem.hit.deltat_deconv",deltatDeconv);
  C->SetBranchAddress("bb.gem.hit.deltat_fit",deltatFit);

  TFile *fout = new TFile(outfilename,"RECREATE");
  
  //Next: initialize histograms and output file:
  TH1D *hdxfcp = new TH1D("hdxfcp", "Front Constraint X;x_{track}-x_{fcp} (m);", 500, -0.15, 0.15 );
  TH1D *hdyfcp = new TH1D("hdyfcp", "Front Constraint Y;y_{track}-y_{fcp} (m);", 500, -0.15, 0.15 );
  TH1D *hdxbcp = new TH1D("hdxbcp", "Back Constraint X;x_{track}+x'_{track}z_{bcp}-x_{bcp} (m);", 500, -0.15, 0.15 );
  TH1D *hdybcp = new TH1D("hdybcp", "Back Constraint Y;y_{track}+y'_{track}z_{bcp}-y_{bcp} (m);", 500, -0.15, 0.15 );

  TH2D *hdxdyfcp = new TH2D("hdxdyfcp", ";y_{track}-y_{fcp} (m);x_{track}-x_{fcp} (m)", 150, -0.15, 0.15, 150, -0.15, 0.15 );
  TH2D *hdxdybcp = new TH2D("hdxdybcp", ";y_{track}+y'_{track}z_{bcp}-y_{bcp} (m);x_{track}+x'_{track}z_{bcp}-x_{bcp} (m)", 150, -0.15, 0.15, 150, -0.15, 0.15 );

  TH1D *hdthcp = new TH1D("hdthcp", "Track #theta Constraint; x'_{track}-x'_{constraint};", 300, -0.1, 0.1 );
  TH1D *hdphcp = new TH1D("hdphcp", "Track #phi Constraint; y'_{track}-y'_{constraint};", 300, -0.1, 0.1 );

  //Forget TClonesArrays, use TH2D for module-dependent stuff:
  TH2D *hADCmaxsamp_vs_module = new TH2D("hADCmaxsamp_vs_module",";module;max strip max ADC sample", nmodules, -0.5, nmodules-0.5, 1500,0,3000);
  TH2D *hADCmaxstrip_vs_module = new TH2D("hADCmaxstrip_vs_module",";module;max strip ADC sum", nmodules, -0.5, nmodules-0.5, 1500,0,15000);

  TH2D *hADCclust_vs_module = new TH2D("hADCclust_vs_module",";module;cluster sum", nmodules, -0.5, nmodules-0.5, 1500,0,30000);


  TH2D *hdeltat_mod = new TH2D("hdeltat_mod", "Standard hit times ;module;#Delta t (ns)", nmodules, -0.5, nmodules-0.5, 200,-30,30);
  TH2D *hdeltat_deconv_mod = new TH2D("hdeltat_deconv_mod", "Deconvoluted hit times ;module;#Delta t (ns)", nmodules, -0.5, nmodules-0.5, 200,-80,80);
  TH2D *hdeltat_fit_mod = new TH2D("hdeltat_fit_mod", "Fitted hit times ;module;#Delta t (ns)", nmodules, -0.5, nmodules-0.5, 200,-30,30);

  TH2D *hADCasym_mod = new TH2D("hADCasym_mod", ";module; ADC asymmetry", nmodules,-0.5,nmodules-0.5, 200,-0.8,0.8);
  TH2D *hADCasym_deconv_mod = new TH2D("hADCasym_deconv_mod", ";module; ADC asymmetry", nmodules,-0.5,nmodules-0.5, 200,-1.0,1.0);

  TH2D *hADCratio_mod = new TH2D("hADCratio_mod", ";module; ADC V/U ratio", nmodules,-0.5,nmodules-0.5, 200,0.0,2.0);
  TH2D *hADCratio_deconv_mod = new TH2D("hADCratio_deconv_mod", ";module; ADC V/U ratio", nmodules,-0.5,nmodules-0.5, 200,0.0,2.0);

  TH2D *hmaxstrip_tU_mod = new TH2D("hmaxstrip_tU_mod", "Max U strip ;module; time (ns)", nmodules, -0.5, nmodules-0.5, 150,0,150); 
  TH2D *hmaxstrip_tV_mod = new TH2D("hmaxstrip_tV_mod", "Max V strip ;module; time (ns)", nmodules, -0.5, nmodules-0.5, 150,0,150); 

  TH2D *hmaxstrip_tU_deconv_mod = new TH2D("hmaxstrip_tU_deconv_mod", "Max U strip ;module; deconv time (ns)", nmodules, -0.5, nmodules-0.5, 150,-100,200); 
  TH2D *hmaxstrip_tV_deconv_mod = new TH2D("hmaxstrip_tV_deconv_mod", "Max V strip ;module; deconv time (ns)", nmodules, -0.5, nmodules-0.5, 150,-100,200); 

  TH2D *hmaxstrip_tU_fit_mod = new TH2D("hmaxstrip_tU_fit_mod", "Max U strip ;module; fit time (ns)", nmodules, -0.5, nmodules-0.5, 150,-100,200); 
  TH2D *hmaxstrip_tV_fit_mod = new TH2D("hmaxstrip_tV_fit_mod", "Max V strip ;module; fit time (ns)", nmodules, -0.5, nmodules-0.5, 150,-100,200); 

  TH2D *htU_mod = new TH2D("htU_mod", "U clusters;module; time (ns)", nmodules,-0.5,nmodules-0.5, 150,0,150);
  TH2D *htV_mod = new TH2D("htV_mod", "V clusters;module; time (ns)", nmodules,-0.5,nmodules-0.5, 150,0,150);

  TH2D *htU_deconv_mod = new TH2D("htU_deconv_mod", "U clusters;module; time deconv (ns)", nmodules,-0.5,nmodules-0.5, 150,-100,200);
  TH2D *htV_deconv_mod = new TH2D("htV_deconv_mod", "V clusters;module; time deconv (ns)", nmodules,-0.5,nmodules-0.5, 150,-100,200);
  
  TH2D *htU_fit_mod = new TH2D("htU_fit_mod", "U clusters;module; fit time (ns)", nmodules,-0.5,nmodules-0.5, 150,-100,200);
  TH2D *htV_fit_mod = new TH2D("htV_fit_mod", "V clusters;module; fit time (ns)", nmodules,-0.5,nmodules-0.5, 150,-100,200);

  TH2D *hdt_corr_mod = new TH2D("hdt_corr_mod", "Hit average corrected time minus track t0;T_{corr}-t_{0}^{track} (ns);", nmodules,-0.5,nmodules-0.5, 100,-50,50);
  
  TString fname_db = outfilename;
  fname_db.ReplaceAll(".root",".dat");
  
  ofstream dbfile(fname_db.Data());

  fname_db.ReplaceAll(".dat", "_timecuts.dat");
  ofstream dbfile_tcuts(fname_db.Data() );

  long nevent=0;

  TTreeFormula *GlobalCut = new TTreeFormula( "GlobalCut", globalcut, C );

  int treenum=0, currenttreenum=0;

  while( C->GetEntry( nevent++ ) ){
    currenttreenum = C->GetTreeNumber();
    if( nevent == 1 || currenttreenum != treenum ){
      treenum = currenttreenum;
      GlobalCut->UpdateFormulaLeaves();
    }

    if( nevent % 100000 == 0 ) cout << nevent << endl;

    bool passedcut = GlobalCut->EvalInstance(0) != 0;

    // cout << "passed cut, ntracks, ngoodhits = " << passedcut << ", "
    // 	 << ntracks << ", " << ngoodhits << endl;


    
    if( passedcut && int(ntracks) >= 1){

      hdxfcp->Fill( xfp[0] - xfcp[0] );
      hdyfcp->Fill( yfp[0] - yfcp[0] );
      hdxbcp->Fill( xfp[0]+thfp[0]*(zbcp[0]-zfcp[0]) - xbcp[0] );
      hdybcp->Fill( yfp[0]+phfp[0]*(zbcp[0]-zfcp[0]) - ybcp[0] );

      hdxdyfcp->Fill( yfp[0]-yfcp[0], xfp[0]-xfcp[0] );
      hdxdybcp->Fill( yfp[0]+phfp[0]*(zbcp[0]-zfcp[0]) - ybcp[0],
		      xfp[0]+thfp[0]*(zbcp[0]-zfcp[0]) - xbcp[0] );
      
      double thcp = (xbcp[0]-xfcp[0])/(zbcp[0]-zfcp[0]);
      double phcp = (ybcp[0]-yfcp[0])/(zbcp[0]-zfcp[0]);

      hdthcp->Fill( thfp[0]-thcp );
      hdphcp->Fill( phfp[0]-phcp );
      
      for( int ihit=0; ihit<int(ngoodhits); ihit++ ){
	if( int(trackindex[ihit]) == 0 && nstripu[ihit]>1&&nstripv[ihit]>1 ){
	  
	  
	  hADCmaxsamp_vs_module->Fill( module[ihit], 0.5*(ADCmaxsampU[ihit]+ADCmaxsampV[ihit]) );
	  hADCmaxstrip_vs_module->Fill( module[ihit], 0.5*(ADCmaxstripU[ihit]+ADCmaxstripV[ihit]) );

	  hADCclust_vs_module->Fill( module[ihit], ADCavg[ihit] );

	  if( ADCavg[ihit] >= 1500.0 ){
	    hADCasym_mod->Fill( module[ihit], ADCasym[ihit] );
	    hADCasym_deconv_mod->Fill( module[ihit], ADCasymDeconv[ihit] );
	    hADCratio_mod->Fill( module[ihit], ADCV[ihit]/ADCU[ihit] );
	    hADCratio_deconv_mod->Fill( module[ihit], DeconvADCV[ihit]/DeconvADCU[ihit] );

	    hdeltat_mod->Fill( module[ihit], deltat[ihit] );
	    hdeltat_deconv_mod->Fill( module[ihit], deltatDeconv[ihit] );
	    hdeltat_fit_mod->Fill( module[ihit], deltatFit[ihit] );
	  }

	  if( ADCavg[ihit] >= 750.0 ){
	    hmaxstrip_tU_mod->Fill( module[ihit], UtimeMaxStrip[ihit] );
	    hmaxstrip_tV_mod->Fill( module[ihit], VtimeMaxStrip[ihit] );
	    
	    hmaxstrip_tU_deconv_mod->Fill( module[ihit], UtimeMaxStripDeconv[ihit] );
	    hmaxstrip_tV_deconv_mod->Fill( module[ihit], VtimeMaxStripDeconv[ihit] );
	    
	    hmaxstrip_tU_fit_mod->Fill( module[ihit], UtimeMaxStripFit[ihit] );
	    hmaxstrip_tV_fit_mod->Fill( module[ihit], VtimeMaxStripFit[ihit] );
	    
	    htU_mod->Fill( module[ihit], Utime[ihit] );
	    htV_mod->Fill( module[ihit], Vtime[ihit] );
	    
	    htU_deconv_mod->Fill( module[ihit], UtimeDeconv[ihit] );
	    htV_deconv_mod->Fill( module[ihit], VtimeDeconv[ihit] );
	    
	    htU_fit_mod->Fill( module[ihit], UtimeFit[ihit] );
	    htV_fit_mod->Fill( module[ihit], VtimeFit[ihit] );

	    hdt_corr_mod->Fill( module[ihit], HitTavg_corr[ihit]-trackt0[0] );
	  }
	}
      }
    }
  }

  //// Here we create a pdf which will be filled in the following loop ////

  TString outfilepdf = outfilename;
  outfilepdf.ReplaceAll(".root",".pdf");

  TString outfilepdf_time = outfilename;
  outfilepdf_time.ReplaceAll(".root","_timecuts.pdf");

  TCanvas *c1 = new TCanvas("c1","",1600,1200);
  c1->Divide(2,3);
  c1->cd(1);
  DrawHist(hdxfcp,{hdxfcp->GetMean(),hdxfcp->GetRMS()*4.5});
  c1->cd(2);
  DrawHist(hdyfcp,{hdyfcp->GetMean(),hdyfcp->GetRMS()*4.5});
  c1->cd(3);
  DrawHist(hdxbcp,{hdxbcp->GetMean(),hdxbcp->GetRMS()*4.5});
  c1->cd(4);
  DrawHist(hdybcp,{hdybcp->GetMean(),hdybcp->GetRMS()*4.5});
  c1->cd(5);
  DrawHist(hdthcp,{hdthcp->GetMean(),hdthcp->GetRMS()*4.5});
  c1->cd(6);
  DrawHist(hdphcp,{hdphcp->GetMean(),hdphcp->GetRMS()*4.5});

  c1->Print(outfilepdf + "(");  //Open the pdf and make the first page

  ///////////////////////////////////////////////////////////////////
  double nsigma = 4.5;
  
  double maxstrip_t0[nmodules][2];
  double maxstrip_tsigma[nmodules][2];
  //double maxstrip_tcut[nmodules][2];

  double maxstrip_t0_deconv[nmodules][2];
  double maxstrip_tsigma_deconv[nmodules][2];
  //double maxstrip_tcut_deconv[nmodules][2];

  double maxstrip_t0_fit[nmodules][2];
  double maxstrip_tsigma_fit[nmodules][2];

  double tmean[nmodules][6];
  double tsigma[nmodules][6];
 
  //double maxstrip_tcut_fit[nmodules][2];
  
  //Let's start with timing cuts:
  for( int imod=0; imod<nmodules; imod++ ){

    TH1D *hADCtemp = hADCmaxsamp_vs_module->ProjectionY("hADCtemp",imod+1,imod+1);
    ///Create canvas which will hold the ADC plots //////
    TCanvas *c= new TCanvas("c","",1200,1200);
    c->Divide(1,3);

    if( hADCtemp->GetEntries() >= 300 ){

      int binlow=1;

      while( hADCtemp->Integral(1,binlow) < 0.001*hADCtemp->GetEntries() ){binlow++;}

      vector<double> threshsamp;
      threshsamp.push_back(hADCtemp->GetBinLowEdge(binlow));
      
      TString dbline;
      dbfile << dbline.Format("bb.gem.m%d.threshold_sample = %12.5g", imod, threshsamp[0] ) << endl;

      ///// Add this histogram to the pdf output
      TH1D *htemp = hADCtemp;
      htemp->SetName(Form("hsample_m%i",imod));
      htemp->SetTitle(Form("Module %i Sample ADC",imod));
      c->cd(1);
      DrawHist(htemp,threshsamp,1);
    }

    hADCtemp = hADCmaxstrip_vs_module->ProjectionY("hADCtemp",imod+1,imod+1);

    if( hADCtemp->GetEntries() >= 300 ){
      int binlow=1;
      
      while( hADCtemp->Integral(1,binlow) < 0.001*hADCtemp->GetEntries() ){binlow++;}
      
      vector<double> threshstrip;
      threshstrip.push_back(hADCtemp->GetBinLowEdge(binlow));
      
      TString dbline;
      dbfile << dbline.Format("bb.gem.m%d.threshold_stripsum = %12.5g", imod, threshstrip[0] ) << endl;

      ///// Add this histogram to the pdf output
      TH1D *htemp = hADCtemp;
      htemp->SetName(Form("hstrip_m%i",imod));
      htemp->SetTitle(Form("Module %i Strip Sum ADC",imod));
      c->cd(2);
      DrawHist(htemp,threshstrip,1);
    }

    hADCtemp = hADCclust_vs_module->ProjectionY("hADCtemp",imod+1,imod+1);

    if( hADCtemp->GetEntries() >= 300 ){
      int binlow=1;
      
      while( hADCtemp->Integral(1,binlow) < 0.001*hADCtemp->GetEntries() ){binlow++;}
      
      vector<double> thresh;
      thresh.push_back(hADCtemp->GetBinLowEdge(binlow));
      
      TString dbline;
      dbfile << dbline.Format("bb.gem.m%d.threshold_clustersum = %12.5g", imod, thresh[0] ) << endl << endl;

      ///// Add this histogram to the pdf output
      TH1D *htemp = hADCtemp;
      htemp->SetName(Form("hclust_m%i",imod));
      htemp->SetTitle(Form("Module %i Cluster Sum ADC",imod));
      c->cd(3);
      DrawHist(htemp,thresh,1);
    }
    

    TCanvas *c2= new TCanvas("c2","",1600,1200);
    c2->Divide(2,3);

    TString hname;
    TH1D *htempU, *htempV;

    htempU = hmaxstrip_tU_mod->ProjectionY( "htempU", imod+1, imod+1 );
    htempV = hmaxstrip_tV_mod->ProjectionY( "htempV", imod+1, imod+1 );
    
    if( htempU->GetEntries() >= 300 ){
      
      FitGaus_FWHM( htempU, 0.3 );
      FitGaus_FWHM( htempV, 0.3 );
      
      TF1 *fitfuncU = (TF1*) (htempU->GetListOfFunctions()->FindObject("gaus"));
      TF1 *fitfuncV = (TF1*) (htempV->GetListOfFunctions()->FindObject("gaus"));
      
      maxstrip_t0[imod][0] = fitfuncU->GetParameter("Mean");
      maxstrip_t0[imod][1] = fitfuncV->GetParameter("Mean");
      maxstrip_tsigma[imod][0] = fitfuncU->GetParameter("Sigma");
      maxstrip_tsigma[imod][1] = fitfuncV->GetParameter("Sigma");

      TString dbline;
      dbfile_tcuts << dbline.Format("bb.gem.m%d.maxstrip_t0 = %10.4g %10.4g",imod,maxstrip_t0[imod][0], maxstrip_t0[imod][1] ) << endl;
      dbfile_tcuts << dbline.Format("bb.gem.m%d.maxstrip_tsigma = %10.4g %10.4g",imod,maxstrip_tsigma[imod][0], maxstrip_tsigma[imod][1] ) << endl;
      dbfile_tcuts << dbline.Format("bb.gem.m%d.maxstrip_tcut = %10.4g %10.4g", imod, nsigma, nsigma ) << endl << endl;
      
      ///// Now we create the histograms for the pdf results//////
      vector<double> resultsU, resultsV;
      resultsU.push_back(maxstrip_t0[imod][0]);
      resultsU.push_back(maxstrip_tsigma[imod][0]);
      resultsV.push_back(maxstrip_t0[imod][1]);
      resultsV.push_back(maxstrip_tsigma[imod][1]);

      TH1D *htU = htempU;
      htU->SetName(Form("htmaxstripU_m%i",imod));
      htU->SetTitle(Form("Module %i U Time;time (ns)",imod));
      c2->cd(1);
      DrawHist(htU,resultsU,2);
      fitfuncU->Draw("same");

      TH1D *htV = htempV;
      htV->SetName(Form("htmaxstripV_m%i",imod));
      htV->SetTitle(Form("Module %i V Time;time (ns)",imod));
      c2->cd(2);
      DrawHist(htV,resultsV,2);
      fitfuncV->Draw("same");
    }

    htempU = hmaxstrip_tU_deconv_mod->ProjectionY("htempU", imod+1,imod+1);
    htempV = hmaxstrip_tV_deconv_mod->ProjectionY("htempV", imod+1,imod+1);
    
    if( htempU->GetEntries() >= 300 ){
      
      FitGaus_FWHM( htempU, 0.3 );
      FitGaus_FWHM( htempV, 0.3 );
      
      TF1 *fitfuncU = (TF1*) (htempU->GetListOfFunctions()->FindObject("gaus"));
      TF1 *fitfuncV = (TF1*) (htempV->GetListOfFunctions()->FindObject("gaus"));
      
      maxstrip_t0_deconv[imod][0] = fitfuncU->GetParameter("Mean");
      maxstrip_t0_deconv[imod][1] = fitfuncV->GetParameter("Mean");
      maxstrip_tsigma_deconv[imod][0] = fitfuncU->GetParameter("Sigma");
      maxstrip_tsigma_deconv[imod][1] = fitfuncV->GetParameter("Sigma");

      TString dbline;
      dbfile_tcuts << dbline.Format("bb.gem.m%d.maxstrip_t0_deconv = %10.4g %10.4g",imod,maxstrip_t0_deconv[imod][0], maxstrip_t0_deconv[imod][1] ) << endl;
      dbfile_tcuts << dbline.Format("bb.gem.m%d.maxstrip_tsigma_deconv = %10.4g %10.4g",imod,maxstrip_tsigma_deconv[imod][0], maxstrip_tsigma_deconv[imod][1] ) << endl;
      dbfile_tcuts << dbline.Format("bb.gem.m%d.maxstrip_tcut_deconv = %10.4g %10.4g", imod, nsigma, nsigma ) << endl << endl;

      ///// Now we create the histograms for the pdf results//////
      vector<double> resultsU, resultsV;
      resultsU.push_back(maxstrip_t0_deconv[imod][0]);
      resultsU.push_back(maxstrip_tsigma_deconv[imod][0]);
      resultsV.push_back(maxstrip_t0_deconv[imod][1]);
      resultsV.push_back(maxstrip_tsigma_deconv[imod][1]);
      
      TH1D *htU = htempU;
      htU->SetName(Form("htstripdeconU_m%i",imod));
      htU->SetTitle(Form("Module %i U Deconvolution;deconv time (ns)",imod));
      c2->cd(3);
      DrawHist(htU,resultsU,2);
      fitfuncU->Draw("same");

      TH1D *htV = htempV;
      htV->SetName(Form("htstripdeconV_m%i",imod));
      htV->SetTitle(Form("Module %i V Deconvolution;deconv time (ns)",imod));
      c2->cd(4);
      DrawHist(htV,resultsV,2);
      fitfuncV->Draw("same");
    }
 
    htempU = hmaxstrip_tU_fit_mod->ProjectionY("htempU", imod+1,imod+1);
    htempV = hmaxstrip_tV_fit_mod->ProjectionY("htempV", imod+1,imod+1);
    
    if( htempU->GetEntries() >= 300 ){
      
      FitGaus_FWHM( htempU, 0.3 );
      FitGaus_FWHM( htempV, 0.3 );
      
      TF1 *fitfuncU = (TF1*) (htempU->GetListOfFunctions()->FindObject("gaus"));
      TF1 *fitfuncV = (TF1*) (htempV->GetListOfFunctions()->FindObject("gaus"));
      
      maxstrip_t0_fit[imod][0] = fitfuncU->GetParameter("Mean");
      maxstrip_t0_fit[imod][1] = fitfuncV->GetParameter("Mean");
      maxstrip_tsigma_fit[imod][0] = fitfuncU->GetParameter("Sigma");
      maxstrip_tsigma_fit[imod][1] = fitfuncV->GetParameter("Sigma");

      TString dbline;
      dbfile_tcuts << dbline.Format("bb.gem.m%d.maxstrip_t0_fit = %10.4g %10.4g",imod,maxstrip_t0_fit[imod][0], maxstrip_t0_fit[imod][1] ) << endl;
      dbfile_tcuts << dbline.Format("bb.gem.m%d.maxstrip_tsigma_fit = %10.4g %10.4g",imod,maxstrip_tsigma_fit[imod][0], maxstrip_tsigma_fit[imod][1] ) << endl;
      dbfile_tcuts << dbline.Format("bb.gem.m%d.maxstrip_tcut_fit = %10.4g %10.4g", imod, nsigma, nsigma ) << endl << endl;


      ///// Now we create the histograms for the pdf results//////
      vector<double> resultsU, resultsV;
      resultsU.push_back(maxstrip_t0_fit[imod][0]);
      resultsU.push_back(maxstrip_tsigma_fit[imod][0]);
      resultsV.push_back(maxstrip_t0_fit[imod][1]);
      resultsV.push_back(maxstrip_tsigma_fit[imod][1]);
      
      TH1D *htU = htempU;
      htU->SetName(Form("htstripfitU_m%i",imod));
      htU->SetTitle(Form("Module %i U Fit Time;deconv time (ns)",imod));
      c2->cd(5);
      DrawHist(htU,resultsU,2);
      fitfuncU->Draw("same");

      TH1D *htV = htempV;
      htV->SetName(Form("htstripfitV_m%i",imod));
      htV->SetTitle(Form("Module %i V Fit Time;deconv time (ns)",imod));
      c2->cd(6);
      DrawHist(htV,resultsV,2);
      fitfuncV->Draw("same");
    }

    TCanvas *c3= new TCanvas("c3","",1200,800);
    c3->Divide(2,2);


    TH1D *htemp;

    htemp = hdeltat_mod->ProjectionY("htemp", imod+1,imod+1 );

    if( htemp->GetEntries() >= 300 ){
      FitGaus_FWHM( htemp, 0.3 );

      TF1 *fitfunc = (TF1*) htemp->GetListOfFunctions()->FindObject("gaus");

      double mean = fitfunc->GetParameter("Mean");
      double sigma = fitfunc->GetParameter("Sigma");

      TString dbline;
      dbfile_tcuts << dbline.Format("bb.gem.m%d.deltat_sigma = %12.5g", imod, sigma ) << endl;
      dbfile_tcuts << dbline.Format("bb.gem.m%d.deltat_cut = %12.5g", imod, 10.0*sigma ) << endl;

      vector<double> results = {10.0*sigma};

      TH1D *hplot = htemp;
      hplot->SetName(Form("hdeltat_m%i",imod));
      hplot->SetTitle(Form("Module %i #Deltat; #Deltat (ns)",imod));
      c3->cd(1);
      DrawHist(hplot,results,2);
      fitfunc->Draw("same");
    }

    htemp = hdeltat_deconv_mod->ProjectionY("htemp", imod+1,imod+1 );

    if( htemp->GetEntries() >= 300 ){
      FitGaus_FWHM( htemp, 0.3 );

      TF1 *fitfunc = (TF1*) htemp->GetListOfFunctions()->FindObject("gaus");

      double mean = fitfunc->GetParameter("Mean");
      double sigma = fitfunc->GetParameter("Sigma");

      TString dbline;
      dbfile_tcuts << dbline.Format("bb.gem.m%d.deltat_sigma_deconv = %12.5g", imod, sigma ) << endl;
      dbfile_tcuts << dbline.Format("bb.gem.m%d.deltat_cut_deconv = %12.5g", imod, 10.0*sigma ) << endl;

      vector<double> results = {10.0*sigma};

      TH1D *hplot = htemp;
      hplot->SetName(Form("hdeltatdeconv_m%i",imod));
      hplot->SetTitle(Form("Module %i Deconv #Deltat;deconv #Deltat(ns)",imod));
      c3->cd(2);
      DrawHist(hplot,results,2);
      fitfunc->Draw("same");
    }

    htemp = hdeltat_fit_mod->ProjectionY("htemp", imod+1,imod+1 );

    if( htemp->GetEntries() >= 300 ){
      FitGaus_FWHM( htemp, 0.3 );

      TF1 *fitfunc = (TF1*) htemp->GetListOfFunctions()->FindObject("gaus");

      double mean = fitfunc->GetParameter("Mean");
      double sigma = fitfunc->GetParameter("Sigma");

      TString dbline;
      dbfile_tcuts << dbline.Format("bb.gem.m%d.deltat_sigma_fit = %12.5g", imod, sigma ) << endl;
      dbfile_tcuts << dbline.Format("bb.gem.m%d.deltat_cut_fit = %12.5g", imod, 10.0*sigma ) << endl << endl;

      vector<double> results = {10.0*sigma};

      TH1D *hplot = htemp;
      hplot->SetName(Form("hdeltatfit_m%i",imod));
      hplot->SetTitle(Form("Module %i Fit #Deltat;Fit #Deltat (ns)",imod));
      c3->cd(3);
      DrawHist(hplot,results,2);
      fitfunc->Draw("same");
    }


    TCanvas *c4= new TCanvas("c4","",1600,1200);
    c4->Divide(2,3);

    htempU = htU_mod->ProjectionY("htempU",imod+1,imod+1);
    htempV = htV_mod->ProjectionY("htempV",imod+1,imod+1);

    if( htempU->GetEntries() >= 300 ){
      FitGaus_FWHM( htempU, 0.3 );
      FitGaus_FWHM( htempV, 0.3 );
      TF1 *fitfuncU = (TF1*) (htempU->GetListOfFunctions()->FindObject("gaus"));
      TF1 *fitfuncV = (TF1*) (htempV->GetListOfFunctions()->FindObject("gaus"));

      tmean[imod][0] = fitfuncU->GetParameter("Mean");
      tmean[imod][1] = fitfuncV->GetParameter("Mean");
      tsigma[imod][0] = fitfuncU->GetParameter("Sigma");
      tsigma[imod][1] = fitfuncV->GetParameter("Sigma");

      TString dbline;

      dbfile_tcuts << dbline.Format("bb.gem.m%d.HitTimeMean = %12.5g %12.5g", imod, tmean[imod][0], tmean[imod][1] ) << endl;
      dbfile_tcuts << dbline.Format("bb.gem.m%d.HitTimeSigma = %12.5g %12.5g", imod, tsigma[imod][0], tsigma[imod][1] ) << endl;

      ///// Now we create the histograms for the pdf results//////
      vector<double> resultsU, resultsV;
      resultsU.push_back(tmean[imod][0]);
      resultsU.push_back(tsigma[imod][0]);
      resultsV.push_back(tmean[imod][1]);
      resultsV.push_back(tsigma[imod][1]);

      TH1D *htU = htempU;
      htU->SetName(Form("htclustU_m%i",imod));
      htU->SetTitle(Form("Module %i U Cluster Time;deconv time (ns)",imod));
      c4->cd(1);
      DrawHist(htU,resultsU,2);
      fitfuncU->Draw("same");

      TH1D *htV = htempV;
      htV->SetName(Form("htclustV_m%i",imod));
      htV->SetTitle(Form("Module %i V Cluster Time;deconv time (ns)",imod));
      c4->cd(2);
      DrawHist(htV,resultsV,2);
      fitfuncV->Draw("same");
      
    }

    htempU = htU_deconv_mod->ProjectionY("htempU",imod+1,imod+1);
    htempV = htV_deconv_mod->ProjectionY("htempV",imod+1,imod+1);

    if( htempU->GetEntries() >= 300 ){
      FitGaus_FWHM( htempU, 0.3 );
      FitGaus_FWHM( htempV, 0.3 );
      TF1 *fitfuncU = (TF1*) (htempU->GetListOfFunctions()->FindObject("gaus"));
      TF1 *fitfuncV = (TF1*) (htempV->GetListOfFunctions()->FindObject("gaus"));

      tmean[imod][2] = fitfuncU->GetParameter("Mean");
      tmean[imod][3] = fitfuncV->GetParameter("Mean");
      tsigma[imod][2] = fitfuncU->GetParameter("Sigma");
      tsigma[imod][3] = fitfuncV->GetParameter("Sigma");

      TString dbline;

      dbfile_tcuts << dbline.Format("bb.gem.m%d.HitTimeMeanDeconv = %12.5g %12.5g", imod, tmean[imod][2], tmean[imod][3] ) << endl;
      dbfile_tcuts << dbline.Format("bb.gem.m%d.HitTimeSigmaDeconv = %12.5g %12.5g", imod, tsigma[imod][2], tsigma[imod][3] ) << endl;

      ///// Now we create the histograms for the pdf results//////
      vector<double> resultsU, resultsV;
      resultsU.push_back(tmean[imod][0]);
      resultsU.push_back(tsigma[imod][0]);
      resultsV.push_back(tmean[imod][1]);
      resultsV.push_back(tsigma[imod][1]);

      TH1D *htU = htempU;
      htU->SetName(Form("htclustdeconU_m%i",imod));
      htU->SetTitle(Form("Module %i U Cluster Deconv Time;deconv time (ns)",imod));
      c4->cd(3);
      DrawHist(htU,resultsU,2);
      fitfuncU->Draw("same");

      TH1D *htV = htempV;
      htV->SetName(Form("htclustdeconV_m%i",imod));
      htV->SetTitle(Form("Module %i V Cluster Deconv Time;deconv time (ns)",imod));
      c4->cd(4);
      DrawHist(htV,resultsV,2);
      fitfuncV->Draw("same");
      
    }

    htempU = htU_fit_mod->ProjectionY("htempU",imod+1,imod+1);
    htempV = htV_fit_mod->ProjectionY("htempV",imod+1,imod+1);

    if( htempU->GetEntries() >= 300 ){
      FitGaus_FWHM( htempU, 0.3 );
      FitGaus_FWHM( htempV, 0.3 );
      TF1 *fitfuncU = (TF1*) (htempU->GetListOfFunctions()->FindObject("gaus"));
      TF1 *fitfuncV = (TF1*) (htempV->GetListOfFunctions()->FindObject("gaus"));

      tmean[imod][4] = fitfuncU->GetParameter("Mean");
      tmean[imod][5] = fitfuncV->GetParameter("Mean");
      tsigma[imod][4] = fitfuncU->GetParameter("Sigma");
      tsigma[imod][5] = fitfuncV->GetParameter("Sigma");

      TString dbline;

      dbfile_tcuts << dbline.Format("bb.gem.m%d.HitTimeMeanFit = %12.5g %12.5g", imod, tmean[imod][4], tmean[imod][5] ) << endl;
      dbfile_tcuts << dbline.Format("bb.gem.m%d.HitTimeSigmaFit = %12.5g %12.5g", imod, tsigma[imod][4], tsigma[imod][5] ) << endl << endl;

      ///// Now we create the histograms for the pdf results//////
      vector<double> resultsU, resultsV;
      resultsU.push_back(tmean[imod][0]);
      resultsU.push_back(tsigma[imod][0]);
      resultsV.push_back(tmean[imod][1]);
      resultsV.push_back(tsigma[imod][1]);

      TH1D *htU = htempU;
      htU->SetName(Form("htclustfitU_m%i",imod));
      htU->SetTitle(Form("Module %i U Cluster Fit Time;deconv time (ns)",imod));
      c4->cd(5);
      DrawHist(htU,resultsU,2);
      fitfuncU->Draw("same");

      TH1D *htV = htempV;
      htV->SetName(Form("htclustfitV_m%i",imod));
      htV->SetTitle(Form("Module %i V Cluster Fit Time;deconv time (ns)",imod));
      c4->cd(6);
      DrawHist(htV,resultsV,2);
      fitfuncV->Draw("same");
      
    }

    TCanvas *c5= new TCanvas("c5","",1200,800);
    c5->Divide(2,2);

    htemp = hADCasym_mod->ProjectionY("htemp",imod+1,imod+1);

    if( htemp->GetEntries() >= 300 ){
      FitGaus_FWHM( htemp, 0.4 );
      TF1 *fitfunc = (TF1*) (htemp->GetListOfFunctions()->FindObject("gaus"));
      double mean = fitfunc->GetParameter("Mean");
      double sigma = fitfunc->GetParameter("Sigma");

      TString dbline;
      dbfile_tcuts << dbline.Format("bb.gem.m%d.ADCasym_sigma = %12.5g", imod, sigma ) << endl;
      
      ///// Now we create the histograms for the pdf results//////
      vector<double> results = {mean,sigma};
      
      TH1D *hplot = htemp;
      hplot->SetName(Form("hADCasym_m%i",imod));
      hplot->SetTitle(Form("Module %i ADC Asymmetry;ADC Asym",imod));
      c5->cd(1);
      DrawHist(hplot,results,2);
      fitfunc->Draw("same");
      
    }

    htemp = hADCratio_mod->ProjectionY("htemp",imod+1,imod+1);

    if( htemp->GetEntries() >= 300 ){
      FitGaus_FWHM( htemp, 0.4 );

      TF1 *fitfunc = (TF1*) (htemp->GetListOfFunctions()->FindObject("gaus"));
      double mean = fitfunc->GetParameter("Mean");
      double sigma = fitfunc->GetParameter("Sigma");

      TString dbline;
      dbfile_tcuts << dbline.Format("bb.gem.m%d.ADCratio_sigma = %12.5g", imod, sigma ) << endl;
      
      ///// Now we create the histograms for the pdf results//////
      vector<double> results = {mean,sigma};
      
      TH1D *hplot = htemp;
      hplot->SetName(Form("hADCratio_m%i",imod));
      hplot->SetTitle(Form("Module %i ADC Ratio;Clust ADC V/U",imod));
      c5->cd(2);
      DrawHist(hplot,results,2);
      fitfunc->Draw("same");
    }

    htemp = hADCratio_deconv_mod->ProjectionY("htemp",imod+1,imod+1);

    if( htemp->GetEntries() >= 300 ){
      FitGaus_FWHM( htemp, 0.4 );

      TF1 *fitfunc = (TF1*) (htemp->GetListOfFunctions()->FindObject("gaus"));
      double mean = fitfunc->GetParameter("Mean");
      double sigma = fitfunc->GetParameter("Sigma");
      
      TString dbline;
      dbfile_tcuts << dbline.Format("bb.gem.m%d.ADCratio_sigma_deconv = %12.5g", imod, sigma ) << endl << endl;

      ///// Now we create the histograms for the pdf results//////
      vector<double> results = {mean,sigma};
      
      TH1D *hplot = htemp;
      hplot->SetName(Form("hADCratiodecon_m%i",imod));
      hplot->SetTitle(Form("Module %i Deconvolution ADC Ratio;Clust V ADC / Clust U ADC",imod));
      c5->cd(3);
      DrawHist(hplot,results,2);
      fitfunc->Draw("same");
      
    }

    htemp = hdt_corr_mod->ProjectionY("htemp",imod+1,imod+1);

    if( htemp->GetEntries() >= 300 ){
      FitGaus_FWHM( htemp, 0.3 );

      TF1 *fitfunc = (TF1*) (htemp->GetListOfFunctions()->FindObject("gaus"));
      double mean = fitfunc->GetParameter("Mean");
      double sigma = fitfunc->GetParameter("Sigma");

      TString dbline;
      dbfile_tcuts << dbline.Format("bb.gem.m%d.sigma_tcorr = %12.5g", imod, sigma ) << endl << endl;
      
      ///// Now we create the histograms for the pdf results//////
      vector<double> results = {mean,sigma};
      
      TH1D *hplot = htemp;
      hplot->SetName(Form("hADCratiodecon_m%i",imod));
      hplot->SetTitle(Form("Module %i Hit Corrected Time Difference;t_{corr}-t_{0}^{track} (ns);",imod));
      c5->cd(4);
      DrawHist(hplot,results,2);
      fitfunc->Draw("same");
    }
    
    ///// Now we will plot the ADC threshold histograms /////
    if(imod == nmodules - 1) c->Print(outfilepdf + ")");
    else c->Print(outfilepdf);

    if(imod == 0) {
      c2->Print(outfilepdf_time + "(");
      c3->Print(outfilepdf_time);
      c4->Print(outfilepdf_time);
      c5->Print(outfilepdf_time);
    }
    else if(imod == nmodules - 1){
      c2->Print(outfilepdf_time);
      c3->Print(outfilepdf_time);
      c4->Print(outfilepdf_time);
      c5->Print(outfilepdf_time + ")");
    }
    else{
      c2->Print(outfilepdf_time);
      c3->Print(outfilepdf_time);
      c4->Print(outfilepdf_time);
      c5->Print(outfilepdf_time);
    }
    /////////////////////////////////////////////////////////////
    
  }

  //Tracking constraints:
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
  
  //hdeltat_mod->Draw("colz");

  fout->Write();

  
}
