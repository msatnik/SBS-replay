#include "TH1D.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TCut.h"
#include "TLegend.h"

const double Mp = 0.938272;
const double mu_p = 2.793;
const double mu_n = -1.913;

TH1D *hdxMCp, *hdxMCn, *hdxdata, *hdxbg;

//Option 0: fit MC p + MC n + polynomial background of (up to) order 4:
double dxtotal_poly(double *x, double *par){
  double Np = par[0];
  double R = par[1];
  
  double dx = x[0];

  double bg = 0.0;
  for( int i=0; i<=4; i++ ){
    bg += par[i+2]*pow(dx,i);
  }

  double Nprotons = hdxMCp->Interpolate(dx);
  double Nneutrons = hdxMCn->Interpolate(dx);
  
  return bg + Np * (Nprotons + R*Nneutrons);
  
}

double dxtotal_gaus(double *x, double *par){
  double Np = par[0];
  double R = par[1];
  double dx = x[0];
  double bg = par[2]*exp(-0.5*pow((dx-par[3])/par[4],2));

  double Nprotons = hdxMCp->Interpolate(dx);
  double Nneutrons = hdxMCn->Interpolate(dx);
  return bg + Np * (Nprotons + R*Nneutrons);
}

//BG flag = 0: polynomial background fit of order (up to) 4
//BG flag = 1: Gaussian background fit
//BG flag = 2: background from data; must supply cut
//BG flag = 3: background from MC; must provide background MC
//for now we only implement options 0,1 (and perhaps 2):

//We can put in some crude GMn extraction too if we generate Q^2 and epsilon distributions; 

void fit_deltax_MCdata_comp(const char *fname_MCp, const char *fname_MCn, const char *fname_data, double xminfit=-3.0, double xmaxfit=2.0, int bgflag=0, int order_bg=4){

  TFile *fout = new TFile("temp.root","RECREATE");
  
  hdxMCp = new TH1D("hdxMCp","Monte Carlo d(e,e'p)n;#Deltax(m);Normalized events",250,xminfit,xmaxfit);
  hdxMCn = new TH1D("hdxMCn","Monte Carlo d(e,e'n)p;#Deltax(m);Normalized events",250,xminfit,xmaxfit);
  hdxdata = new TH1D("hdxdata","Data with cuts;#Deltax(m);Events",250,xminfit,xmaxfit);

  hdxbg = new TH1D("hdxbg","Data with (anti)cuts;#Deltax(m);Events",250,xminfit,xmaxfit);

  TFile *fMCp = new TFile(fname_MCp,"READ");
  TFile *fMCn = new TFile(fname_MCn,"READ");
  TFile *fdata = new TFile(fname_data,"READ");

  TCut trackquality = "abs(xfp-thfp*0.9)<0.3&&abs(vz)<0.065";
  TCut datacut = "W2<1.15&&W2>0&&abs(deltat_adc)<12.0&&EHCAL>0.01&&abs(deltay_4vect)<1.0";

  TCut MCcut = "W2<1.15&&W2>0&&EHCAL>0&&abs(deltay_4vect)<1.0";

  TCut fiducialcut = "abs(xHCAL_expect_4vect+0.75)<1.59&&abs(xHCAL_expect_4vect+0.75-protondeflection)<=1.59&&-0.45<=yHCAL_expect_4vect&&yHCAL_expect_4vect<=0.6";

  TTree *TMCp, *TMCn, *Tdata;

  fMCp->GetObject("Tout",TMCp);
  fMCn->GetObject("Tout",TMCn);
  fdata->GetObject("Tout",Tdata);

  fout->cd();
  
  TMCp->Project("hdxMCp","deltax_4vect",trackquality&&MCcut);
  TMCn->Project("hdxMCn","deltax_4vect",trackquality&&MCcut);
  Tdata->Project("hdxdata","deltax_4vect",trackquality&&datacut);

  hdxdata->Draw();

  hdxdata->SetLineWidth(3);

  hdxMCp->SetLineColor(2);
  hdxMCp->SetLineWidth(3);
  hdxMCp->SetFillColorAlpha(2,0.25);
  //hdxMCp->Draw("SAME");

  hdxMCn->SetLineColor(4);
  hdxMCn->SetLineWidth(3);
  hdxMCn->SetFillColorAlpha(4,0.25);
  //hdxMCn->Draw("SAME");

  TF1 *fitfunc_polybg = new TF1("fitfunc_polybg",dxtotal_poly, xminfit,xmaxfit,7);
  fitfunc_polybg->SetParName(0,"NormP");
  fitfunc_polybg->SetParName(1,"Rnp");
  fitfunc_polybg->SetParName(2,"bgp0");
  fitfunc_polybg->SetParName(3,"bgp1");
  fitfunc_polybg->SetParName(4,"bgp2");
  fitfunc_polybg->SetParName(5,"bgp3");
  fitfunc_polybg->SetParName(6,"bgp4");

  fitfunc_polybg->SetNpx(1000);

  fitfunc_polybg->SetParameter(1,1.0);
  fitfunc_polybg->SetParameter(0,hdxMCp->GetBinContent(hdxMCp->GetMaximumBin())/hdxdata->GetBinContent(hdxMCp->GetMaximumBin()));
  for( int ipar=2; ipar<7; ipar++ ){
    fitfunc_polybg->SetParameter(ipar,0.0);
    if( ipar-2 > order_bg ) fitfunc_polybg->FixParameter(ipar,0.0);
  }

  TF1 *fitfunc_gausbg = new TF1("fitfunc_gausbg",dxtotal_gaus,xminfit,xmaxfit,5);

  fitfunc_gausbg->SetParName(0,"NormP");
  fitfunc_gausbg->SetParName(1,"Rnp");
  fitfunc_gausbg->SetParName(2,"NormBG");
  fitfunc_gausbg->SetParName(3,"MeanBG");
  fitfunc_gausbg->SetParName(4,"SigmaBG");

  fitfunc_gausbg->SetParameter(1,1.0);
  fitfunc_gausbg->SetParameter(0,hdxMCp->GetBinContent(hdxMCp->GetMaximumBin())/hdxdata->GetBinContent(hdxMCp->GetMaximumBin()));
  fitfunc_gausbg->SetParameter(2,1.0);
  fitfunc_gausbg->SetParameter(3,hdxdata->GetMean());
  fitfunc_gausbg->SetParameter(4,1.5);
  fitfunc_gausbg->SetParLimits(3,xminfit,xmaxfit);
  fitfunc_gausbg->SetParLimits(4,0.1,xmaxfit-xminfit);

  double Rnp, dRnp;
  
  if( bgflag==0 ){
    hdxdata->Fit(fitfunc_polybg,"S0","",xminfit,xmaxfit);
    hdxMCp->Scale(fitfunc_polybg->GetParameter(0));
    hdxMCn->Scale(fitfunc_polybg->GetParameter(0)*fitfunc_polybg->GetParameter(1));
    hdxMCp->Draw("SAME");
    hdxMCn->Draw("SAME");

    TF1 *polbg = new TF1("polbg","[0]+[1]*x+[2]*pow(x,2)+[3]*pow(x,3)+[4]*pow(x,4)",xminfit,xmaxfit);

    for( int ipar=0; ipar<5; ipar++){
      polbg->SetParameter(ipar,fitfunc_polybg->GetParameter(ipar+2));
    }

    polbg->SetLineWidth(2);
    polbg->SetLineColor(kGreen+2);

    polbg->SetNpx(1000);
    polbg->Draw("CSAME");
    Rnp = fitfunc_polybg->GetParameter(1);
    dRnp = fitfunc_polybg->GetParError(1);
   
    
  } else if( bgflag == 1 ){
    hdxdata->Fit(fitfunc_gausbg,"S0","",xminfit,xmaxfit);
    hdxMCp->Scale(fitfunc_gausbg->GetParameter(0));
    hdxMCn->Scale(fitfunc_gausbg->GetParameter(0)*fitfunc_polybg->GetParameter(1));
    hdxMCp->Draw("SAME");
    hdxMCn->Draw("SAME");

    TF1 *GausBG = new TF1("GausBG","[0]*exp(-0.5*pow((x-[1])/[2],2))",xminfit,xmaxfit);

    for( int ipar=0; ipar<3; ipar++ ){
      GausBG->SetParameter(ipar,fitfunc_gausbg->GetParameter(ipar+2));
    }

    GausBG->SetLineWidth(2);
    GausBG->SetLineColor(kGreen+2);
    GausBG->SetNpx(1000);
    GausBG->Draw("CSAME");

    Rnp = fitfunc_gausbg->GetParameter(1);
    dRnp = fitfunc_gausbg->GetParError(1);
  }

  //Get Q^2 and epsilon from Monte Carlo distributions:
  double Q2avg_p, Q2avg_n, epsilon_p, epsilon_n;

  fout->cd();
  
  TH1D *hQ2p = new TH1D("hQ2p","Monte Carlo d(e,e'p)n;Q^{2} (GeV^{2});",2000,0.0,20.0);
  TH1D *hQ2n = new TH1D("hQ2n","Monte Carlo d(e,e'n)p;Q^{2} (GeV^{2})",2000,0.0,20.0);
  TH1D *hQ2data = new TH1D("hQ2data","Data with cuts;Q^{2} (GeV^{2});",2000,0.0,20.0);

  TH1D *heps_p = new TH1D("heps_p","Monte Carlo d(e,e'p)n;#varepsilon;",250,0.0,1.0);
  TH1D *heps_n = new TH1D("heps_n","Monte Carlo d(e,e'n)p;#varepsilon;",250,0.0,1.0);
  TH1D *heps_data = new TH1D("heps_data","Data with cuts;#varepsilon;",250,0.0,1.0);

  TMCp->Project("hQ2p","Q2",trackquality&&MCcut);
  TMCn->Project("hQ2n","Q2",trackquality&&MCcut);
  Tdata->Project("hQ2data","Q2",trackquality&&datacut);

  TMCp->Project("heps_p","epsilon",trackquality&&MCcut);
  TMCn->Project("heps_n","epsilon",trackquality&&MCcut);
  Tdata->Project("heps_data","epsilon",trackquality&&datacut);

  double Q2mean = 0.5*(hQ2p->GetMean()+hQ2n->GetMean());
  double epsmean = 0.5*(heps_p->GetMean()+heps_n->GetMean());

  cout << "(Q^2, epsilon)=(" << Q2mean << ", " << epsmean << "), R_{np} = " << Rnp << " +/- " << dRnp << endl;

  // To do a crude, simplistic GMn extraction, we need to calculate sigma_p and GEn using the same parametrization as in g4sbs:

  double tau = Q2mean/(4.0*pow(Mp,2));
  double GD = pow(1.0+Q2mean/0.71,-2);
  double GEn_g4sbs = (1.520*tau + 2.629*tau*tau + 3.055*tau*tau*tau)*GD/(1.0+5.222*tau+0.040*tau*tau+11.438*tau*tau*tau);
  double GMn_g4sbs = -1.913*(1.0+2.33*tau)/(1.0 + 14.72*tau + 24.20*tau*tau + 84.1*tau*tau*tau );

  double GEp_g4sbs = (1.0-0.24*tau)/(1.0 + 10.98*tau + 12.82*tau*tau + 21.97*tau*tau*tau );
  double GMp_g4sbs = 2.79*(1.0+0.12*tau)/(1.0 + 10.97*tau + 18.86*tau*tau + 6.55*tau*tau*tau );

  double Rnp_g4sbs = (epsmean * pow(GEn_g4sbs,2) + tau * pow(GMn_g4sbs,2))/(epsmean * pow(GEp_g4sbs,2)+tau * pow( GMp_g4sbs,2));

  Rnp *= Rnp_g4sbs;
  dRnp *= Rnp_g4sbs;

  double RS_GMp12 = 1.0 -0.46*tau + 0.12 * pow(tau,2);
  double GMp_GMp12 = mu_p * (1.0 + 0.072*tau)/(1.0 + 10.73*tau + 19.81*pow(tau,2) + 4.75 * pow(tau,3) );

  //sigmaR = eps*GE^2 + tau*GM^2 = tau*GM^2*(1.0 + epsilon/tau * pow(GE^2/GM^2));
  // RS = mu_p^2 GE^2 / GM^2
  // sigmaR = tau * GM2 * ( 1.0 + epsilon/mu_p^2*tau * RS );
  double sigmaR_p = tau*pow(GMp_GMp12,2)*(1.0 + epsmean/tau * RS_GMp12/pow(mu_p,2));

  double sigmaR_n = sigmaR_p * Rnp;
  double dsigmaR_n = sigmaR_p * dRnp;

  //sigmaR_n = epsilon * GEn^2 + tau*GMn^2
  double GMn = sqrt( (sigmaR_p * Rnp - epsmean * pow(GEn_g4sbs,2))/tau )/fabs(mu_n)/GD;
  double GMnplus = sqrt( (sigmaR_p * (Rnp + dRnp) - epsmean*pow(GEn_g4sbs,2))/tau )/fabs(mu_n)/GD;
  double GMnminus = sqrt( (sigmaR_p * (Rnp - dRnp) - epsmean*pow(GEn_g4sbs,2))/tau )/fabs(mu_n)/GD;

  double dGMn = 0.5*fabs( GMnplus-GMnminus );

  cout << "GMn/(mu_n GD) = " << GMn << " +/- " << dGMn << endl;
  

  //d (sqrt(x)) = 1/(2 sqrt(x))*dx = 1/(2y)*dx 
    
  // sigma_n = Rnp * sigma_p 
}
