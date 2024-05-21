#include "gmn_tree.C"
#include "TChain.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include <iostream>
#include <fstream>
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TCut.h"
#include "TTreeFormula.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TDecompSVD.h"
#include "TVector3.h"
#include "TLorentzVector.h"

const double Mp=0.938272;

void TOFcal_new(const char *inputfilename, const char *outputfilename="TOFcal_temp.root"){

  ifstream configfile(inputfilename);

  if( !configfile ) return;

  TFile *fout = new TFile(outputfilename,"RECREATE");
  
  TChain *C = new TChain("T");
  
  TString currentline;
  while( currentline.ReadLine(configfile) && !currentline.BeginsWith("endlist") ){
    if( !currentline.BeginsWith("#") ){
      C->Add(currentline.Data());
    }
  }

  TCut globalcut = "";
  
  while( currentline.ReadLine(configfile) && !currentline.BeginsWith("endcut") ){
    if( !currentline.BeginsWith("#") ){
      globalcut += currentline.Data();
    }
  }

  //Default reference PMTs for aligning all other PMTs:
  int hodorefID = 43; // reference BAR (not PMT); we'll set the LEFT PMT to the reference by convention:
  int hcalrefID = 199; //row 16 column 6

  double zhodo = 1.854454; //meters
  double Lbar_hodo = 0.6; //meters
  double etof0 = (1.96+3.0)/0.299792458; // "central" TOF value:

  double thetaHCAL = 29.4;
  double dHCAL = 11.0;
  double dSBS = 2.25;

  double W2min=0.6, W2max=1.2;
  
  double sbsmaxfield = 1.26;
  double sbsfield = 0.7;

  const double Dgap = 48.0*2.54/100.0;
  
  double Ebeam = 5.9826;
  
  while( currentline.ReadLine(configfile) && !currentline.BeginsWith("endconfig") ){
    if( !currentline.BeginsWith("#") ){
      TObjArray *currentline_tokens = currentline.Tokenize(" ");
      if( currentline_tokens->GetEntries() >= 2 ){
	TString skey = ( (TObjString*) (*currentline_tokens)[0] )->GetString();

	TString sval = ( (TObjString*) (*currentline_tokens)[1] )->GetString();
	if( skey.BeginsWith("hodorefID") ){
	  hodorefID = sval.Atoi();
	}

	if( skey.BeginsWith("hcalrefID") ){
	  hcalrefID = sval.Atoi();
	}

	if( skey.BeginsWith("zhodo") ){
	  zhodo = sval.Atof();
	}

	if( skey.BeginsWith("Lbar_hodo") ){
	  Lbar_hodo = sval.Atof();
	}

	if( skey.BeginsWith("etof0") ){
	  etof0 = sval.Atof();
	}
      }
    }
  }

  thetaHCAL *= TMath::Pi()/180.0;
  
  TVector3 zaxis_HCAL(-sin(thetaHCAL),0,cos(thetaHCAL));
  TVector3 xaxis_HCAL(0,-1,0);
  TVector3 yaxis_HCAL = (zaxis_HCAL.Cross(xaxis_HCAL)).Unit();

  TVector3 HCALorigin = dHCAL*zaxis_HCAL;
  
  //int nparams_hodo = 3*180; //offset, walk correction, propagation speed for left and right PMTs for 90 bars minus offset fixed at zero for reference PMT
  //actually, vscint should be determined per paddle;
  int nparams_hodo = 180+2*90; //We actually only really want one zero offset (paddle-specific) that applies to the mean time, so we need 90 zero offsets, 90 vscint values, and 180 walk correction slopes.
  int nparams_HCAL = 224*2; //offset and walk correction for 224 modules.

  //Ordering of params:
  // t0ileft, t0iright, walkileft, walkiright, 1/vileft, 1/viright
  
  TMatrixD Mhodo(nparams_hodo,nparams_hodo);
  TMatrixD Mhcal(nparams_HCAL,nparams_HCAL);
  TVectorD bhodo(nparams_hodo);
  TVectorD bhcal(nparams_HCAL);

  for( int ipar=0; ipar<nparams_hodo; ipar++ ){
    for( int jpar=0; jpar<nparams_hodo; jpar++ ){
      Mhodo(ipar,jpar) = 0.0;
    }
    bhodo(ipar) = 0.0;
  }

  for( int ipar=0; ipar<nparams_HCAL; ipar++ ){
    for( int jpar=0; jpar<nparams_HCAL; jpar++ ){
      Mhcal(ipar,jpar) = 0.0;
    }
    bhcal(ipar) = 0.0;
  }
  
  //HCAL and hodo TDCs both use BigBite trigger as reference time, so both have reference time subtracted:

  TTreeFormula *GlobalCut = new TTreeFormula("GlobalCut",globalcut,C);

  gmn_tree *T = new gmn_tree(C);

  int treenum=-1,oldtreenum=-1;

  long nevent=0;

  while( T->GetEntry(nevent) ){

    if( nevent % 1000 == 0 ){
      cout << "nevent = " << nevent << ", run number = " << T->g_runnum << endl;
    }
    
    treenum = C->GetTreeNumber();
    if( nevent == 0 || treenum != oldtreenum ){
      oldtreenum = treenum;
      GlobalCut->UpdateFormulaLeaves();
    }

    bool passed_global_cut = GlobalCut->EvalInstance(0) != 0;

    if( passed_global_cut ){ //do the things:

      //grab needed track parameters:
      double vz = T->bb_tr_vz[0];
      double pathl = T->bb_tr_pathl[0];
      double etof = pathl/0.299792458; //electron TOF from vertex to hodo.

      double tleft = T->bb_hodotdc_clus_tleft[0];
      double tright = T->bb_hodotdc_clus_tright[0];
      double totleft = T->bb_hodotdc_clus_totleft[0];
      double totright = T->bb_hodotdc_clus_totright[0];

      int ID = int(T->bb_hodotdc_clus_id[0]);

      double yhodo = T->bb_tr_y[0]+zhodo*T->bb_tr_ph[0];

      // how to define our chi2 statistic? We want to correct all hodo PMT times to tvertex = 0:
      // tPMT = etof + t0 - walk * TOT + d/vscint
      // chi2 = sum_{i=1}^N \sum_{j=1}^Nhit (tPMT-etof-t0+walk*TOT - d/vscint)^2
      // dchi2/dt0_k = 2*\sum_i,j (tPMT-etof-t0+walk*TOT-d/vscint)*-1
      // dchi2/dwalk_k = 2*sum_i,j (tPMT-etof-t0+walk*TOT-d/vscint)*TOT
      // dchi2/d(1/v)_k = 2*sum_i,j (tPMT-etof-t0+walk*TOT-d/vscint)*-d

      double dLEFT = std::min(Lbar_hodo,std::max(0.0,Lbar_hodo/2.0 - yhodo));
      double dRIGHT = std::min(Lbar_hodo,std::max(0.0,Lbar_hodo/2.0 + yhodo));

      // cout << "Best hodo hit (tL,tR, totL, totR, yhodo, dL, dR, pathl, etof)=("
      // 	   << tleft << ", " << tright << ", " << totleft << ", " << totright << ", " << yhodo << ", "
      // 	   << dLEFT << ", " << dRIGHT << ", " << pathl << ", " << etof << ")" << endl;

      // we want to minimize
      // chi2 = sum_{i=1}^Nevent (tleft_corr^2 + tright_corr^2)
      // dchi2/dt0_j = sum_i (2*tleft_corr_i*dtleft_corr_i/dt0_j + 2*tright_corr_i*dtright_corr_i/dt0_j) = 0
      //we COULD just use TMinuit but that's lazy.
      //let parameters 0-89 be the t0L, 90-179 be t0R, 180-269 be wL, 270-359 be wR, and 360-449 be 1/vscint:
      int ipar_t0 = ID;
      int ipar_vinv = ID+90;
      //int ipar_t0R = ID+90;
      int ipar_wL = ID+180;
      int ipar_wR = ID+270;
      

      

      // double coeff_t0L = 
      
      //sum_i (tL - etof - t0 + w*TOTL - dL/v) = 0  
      // tL-etof = t0L -w*TOTL + dL/v
      //there is a more elegant way to do this:
      //Left PMT:
      Mhodo(ipar_t0,ipar_t0) += 1.0;
      Mhodo(ipar_t0,ipar_wL) += -totleft;
      Mhodo(ipar_t0,ipar_vinv) += dLEFT;

      Mhodo(ipar_wL, ipar_t0) += -totleft;
      Mhodo(ipar_wL, ipar_wL) += pow(totleft,2);
      Mhodo(ipar_wL, ipar_vinv) += -totleft*dLEFT;

      Mhodo(ipar_vinv, ipar_t0) += dLEFT;
      Mhodo(ipar_vinv, ipar_wL) += -totleft*dLEFT;
      Mhodo(ipar_vinv, ipar_vinv) += pow(dLEFT,2);

      bhodo(ipar_t0) += (tleft - (etof-etof0))*1.0;
      bhodo(ipar_wL) += (tleft - (etof-etof0))*(-totleft);
      bhodo(ipar_vinv) += (tleft -(etof-etof0))*dLEFT;


      //Right PMT:
      Mhodo(ipar_t0,ipar_t0) += 1.0;
      Mhodo(ipar_t0,ipar_wR) += -totright;
      Mhodo(ipar_t0,ipar_vinv) += dRIGHT;

      Mhodo(ipar_wR, ipar_t0) += -totright;
      Mhodo(ipar_wR, ipar_wR) += pow(totright,2);
      Mhodo(ipar_wR, ipar_vinv) += -totright*dRIGHT;

      Mhodo(ipar_vinv, ipar_t0) += dRIGHT;
      Mhodo(ipar_vinv, ipar_wR) += -totright*dRIGHT;
      Mhodo(ipar_vinv, ipar_vinv) += pow(dRIGHT,2);

      bhodo(ipar_t0) += (tright - (etof-etof0))*1.0;
      bhodo(ipar_wR) += (tright - (etof-etof0))*(-totright);
      bhodo(ipar_vinv) += (tright -(etof-etof0))*dRIGHT;
     
      
    }
    nevent++;
  }

  TDecompSVD Ahodo(Mhodo);
  Ahodo.Solve(bhodo);

  bhodo.Print();

  double HODOt0[90], HODOwL[90], HODOwR[90], vscint[90];

  for( int i=0; i<90; i++ ){
    HODOt0[i] = bhodo[i];
    //HODOt0R[i] = bhodo[i+90];
    HODOwL[i] = bhodo[i+180];
    HODOwR[i] = bhodo[i+270];
    vscint[i] = 1.0/bhodo[i+90];
  }

  nevent = 0;

  treenum=-1,oldtreenum=-1;

  cout << "2nd loop, testing new hodoscope calibrations..." << endl;

  fout->cd();
  
  TH2D *htmean_vs_ID_old = new TH2D("htmean_vs_ID_old","OLD ;bar ID; bar t_{mean} (ns)", 90,-0.5,89.5, 150,-30,30);
  TH2D *htmean_vs_ID_new = new TH2D("htmean_vs_ID_new","NEW ;bar ID; bar t_{mean} (ns)", 90,-0.5,89.5, 150,-30,30);

  TH2D *htdiff_yhodo_old = new TH2D("htdiff_yhodo_old", "OLD; y_{hodo} (m); bar t_{L}-t_{R} (ns)", 150,-0.3,0.3,150,-15,15);
  TH2D *htdiff_yhodo_new = new TH2D("htdiff_yhodo_new", "NEW; y_{hodo} (m); bar t_{L}-t_{R} (ns)", 150,-0.3,0.3,150,-15,15);

  //Really want a large quantity of zero-field LH2 data for this purpose:
  TH2D *htdiffHCAL_vs_HCALID_new = new TH2D("htdiffHCAL_vs_HCALID_new", "NEW hodo, old HCAL; HCAL ID; t_{HCAL}-t_{HODO,corr} (ns)",224,-0.5,223.5,250,-50,50); 

  TH2D *hdxdy = new TH2D("hdxdy",";#Deltay (m); #Deltax (m)",250,-1.25,1.25,250,-1.25,1.25);
  
  while( T->GetEntry(nevent) ){

    if( nevent % 1000 == 0 ){
      cout << "nevent = " << nevent << ", run number = " << T->g_runnum << endl;
    }
    
    treenum = C->GetTreeNumber();
    if( nevent == 0 || treenum != oldtreenum ){
      oldtreenum = treenum;
      GlobalCut->UpdateFormulaLeaves();
    }

    bool passed_global_cut = GlobalCut->EvalInstance(0) != 0;

    if( passed_global_cut ){ //do the things:
      double vz = T->bb_tr_vz[0];
      double pathl = T->bb_tr_pathl[0];
      double etof = pathl/0.299792458;

      double tleft = T->bb_hodotdc_clus_tleft[0];
      double tright = T->bb_hodotdc_clus_tright[0];
      double totleft = T->bb_hodotdc_clus_totleft[0];
      double totright = T->bb_hodotdc_clus_totright[0];

      int ID = int(T->bb_hodotdc_clus_id[0]);

      double yhodo = T->bb_tr_y[0]+zhodo*T->bb_tr_ph[0];

      // how to define our chi2 statistic? We want to correct all hodo PMT times to tvertex = 0:
      // tPMT = etof + t0 - walk * TOT + d/vscint
      // chi2 = sum_{i=1}^N \sum_{j=1}^Nhit (tPMT-etof-t0+walk*TOT - d/vscint)^2
      // dchi2/dt0_k = 2*\sum_i,j (tPMT-etof-t0+walk*TOT-d/vscint)*-1
      // dchi2/dwalk_k = 2*sum_i,j (tPMT-etof-t0+walk*TOT-d/vscint)*TOT
      // dchi2/d(1/v)_k = 2*sum_i,j (tPMT-etof-t0+walk*TOT-d/vscint)*-d

      double dLEFT = std::min(Lbar_hodo,std::max(0.0,Lbar_hodo/2.0 - yhodo));
      double dRIGHT = std::min(Lbar_hodo,std::max(0.0,Lbar_hodo/2.0 + yhodo));

      double tleft_CORR = tleft - (etof-etof0) - HODOt0[ID] + HODOwL[ID]*totleft - dLEFT/vscint[ID];
      double tright_CORR = tright - (etof-etof0) - HODOt0[ID] + HODOwR[ID]*totright - dRIGHT/vscint[ID];

      double tmean_old = T->bb_hodotdc_clus_tmean[0];
      double tdiff_old = T->bb_hodotdc_clus_tdiff[0];

      //calculate a corrected time difference without the propagation correction:
      double tdiff_CORR = (tleft - HODOt0[ID]+ HODOwL[ID]*totleft)-
	(tright-HODOt0[ID]+HODOwR[ID]*totright);
      
      htmean_vs_ID_old->Fill( ID, tmean_old );
      htmean_vs_ID_new->Fill( ID, 0.5*(tleft_CORR+tright_CORR) );
      htdiff_yhodo_old->Fill( yhodo, tdiff_old );
      htdiff_yhodo_new->Fill( yhodo, tdiff_CORR );

      double tHCAL = T->sbs_hcal_tdctimeblk;
      double eblkHCAL = T->sbs_hcal_eblk;
      int idblkHCAL = int(T->sbs_hcal_idblk);
      double xHCAL = T->sbs_hcal_x;
      double yHCAL = T->sbs_hcal_y;
      
      //Eventually we'll need to develop parametrized TOF for protons/neutrons for HCAL for different field settings.

      TVector3 ep(T->bb_tr_px[0], T->bb_tr_py[0], T->bb_tr_pz[0]);
      TLorentzVector k(0,0,Ebeam,Ebeam);
      TLorentzVector P(0,0,0,Mp);

      TLorentzVector kprime(ep,ep.Mag());
      TLorentzVector q = k-kprime;
      TLorentzVector Pprime = P + q;

      double Q2 = q.M2();
      
      double etheta = ep.Theta();
      double ephi = ep.Phi();

      
      
      double Eprime_etheta = Ebeam/(1.0+Ebeam/Mp*(1.0-cos(etheta)));
      double Tp_etheta = Ebeam-Eprime_etheta;
      double pp_etheta = sqrt(pow(Tp_etheta,2)+2.0*Mp*Tp_etheta);
      double ptheta_etheta = acos( (Ebeam-Eprime_etheta*cos(etheta))/(pp_etheta));
      double pphi_ephi = ephi + TMath::Pi();

      
      double ptheta_4vect = q.Vect().Theta();
      double pphi_4vect = q.Vect().Phi();

      TVector3 pnhat_expect( sin(ptheta_etheta)*cos(pphi_ephi),
			     sin(ptheta_etheta)*sin(pphi_ephi),
			     cos(ptheta_etheta) );

      TVector3 vertex(0,0,vz);

      double sintersect = (HCALorigin-vertex).Dot(zaxis_HCAL)/(pnhat_expect.Dot(zaxis_HCAL));

      TVector3 HCAL_intersect = vertex + sintersect * pnhat_expect;

      double xHCAL_expect = (HCAL_intersect - HCALorigin).Dot( xaxis_HCAL );
      double yHCAL_expect = (HCAL_intersect - HCALorigin).Dot( yaxis_HCAL );

      double Lpath_HCAL = (HCAL_intersect - vertex).Mag();

      double beta_proton = pp_etheta/sqrt(pow(pp_etheta,2)+pow(Mp,2));
      
      double TOF_HCAL = Lpath_HCAL/(beta_proton*0.299792458);
      double HCALt0 = dHCAL/(beta_proton*0.299792458);

      double W2 = T->e_kine_W2;

      // thetabend = 0.3/p * BdL ;

      double BdL = Dgap * sbsmaxfield * sbsfield;
      double thetabend = 0.3/pp_etheta * BdL;
     
      double protondeflection = tan(thetabend)*(dHCAL-(dSBS+Dgap/2.0)); 
      
      //int idblkHCAL = T->sbs_hcal_idblk;

      double deltax = xHCAL-(xHCAL_expect-protondeflection);
      double deltay = yHCAL-(yHCAL_expect);

      hdxdy->Fill( deltay,deltax );
      
      //if( W2min < W2 && W2 < W2max && sqrt(pow(deltax,2)+pow(deltay,2))<=0.24 && eblkHCAL>0.02 ){
      if( eblkHCAL>0.02 && W2>W2min && W2<W2max && sqrt(pow(deltax,2)+pow(deltay,2))<=0.24 ){
	// cout << "(Eprime,etheta,ephi)=("
	//      << ep.Mag() << ", " << etheta*57.3 << ", "
	//      << ephi*57.3 << ")" << endl;
	// cout << "(ptheta_etheta,pp_etheta)=("
	//      << ptheta_etheta*57.3 << ", " << pp_etheta << ")" << endl;

	// cout << "(ephi,pphi)=(" << ephi*57.3 << ", " << pphi_ephi*57.3 << ")"
	//      << endl;
	
	htdiffHCAL_vs_HCALID_new->Fill( idblkHCAL, tHCAL-0.5*(tleft_CORR+tright_CORR)-(TOF_HCAL-HCALt0));
	//}
      }
    }
    nevent++;
  }

  fout->Write();
  
}
