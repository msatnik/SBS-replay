#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TObjArray.h"
#include "TCut.h"
#include "TObjString.h"
#include "gmn_elastic_tree.C"
#include "TTreeFormula.h"
#include <iostream>
#include <fstream>

#include "TStyle.h"

void HCAL_p_eff_uniformity(const char *setupfilename, const char *outfilename){
  TH1::SetDefaultSumw2();

  gStyle->SetMarkerStyle(20);
  gStyle->SetLineColor(1);
  gStyle->SetMarkerColor(1);
  gROOT->ForceStyle();
  
  ifstream setupfile(setupfilename);
  
  if( !setupfile ) return;

  TChain *C = new TChain("Tout");

  TString currentline;
  while(currentline.ReadLine(setupfile) && !currentline.BeginsWith("endlist")){
    if( !currentline.BeginsWith("#") ){
      C->Add(currentline.Data());
    }
  }

  // we need to define four sets of cuts here: 
  // one is a global/optics/track quality cut
  // two is a W cut
  // three is a "good hcal" cut (dx, dy, possibly also dt, EHCAL)
  // four is a set of "fiducial" cuts in x and y 
  // since cuts 2 and 4 are generally going to be rectangular in nature, we can specify lower and upper limits for these via command-line arguments:

  TCut globalcut = "";
  
  while(currentline.ReadLine(setupfile) && !currentline.BeginsWith("endglobalcut")){
    if( !currentline.BeginsWith("#")) globalcut += currentline.Data();
  }

  TCut hcalcut = "";

  while(currentline.ReadLine(setupfile) && !currentline.BeginsWith("endhcalcut")){
    if( !currentline.BeginsWith("#")) hcalcut += currentline.Data();
  }

  //provide default values for these that can be overridden by key-value pairs in the "config" section:
  //Note that the default "fiducial" limits correspond to the neutron case; i.e., no deflection:
  double W2min=0.6, W2max=1.0, fid_xmin=-2.5, fid_xmax=1.0, fid_ymin=-0.6, fid_ymax=0.6;

  int protondeflectionflag = 0; //1 = use protondeflection from tree, 2 = use mean from config file:

  double pdeflectmean=0.0;
  
  while( currentline.ReadLine(setupfile) && !currentline.BeginsWith("endconfig")){
    if( !currentline.BeginsWith("#") ){
      TObjArray *currentline_tokens = currentline.Tokenize(" ");

      int ntokens = currentline_tokens->GetEntries();

      if( ntokens >= 2 ){
	TString skey = ( (TObjString*) (*currentline_tokens)[0] )->GetString();
	TString sval = ( (TObjString*) (*currentline_tokens)[1] )->GetString();

	if( skey == "W2min" ){
	  W2min = sval.Atof();
	}

	if( skey == "W2max" ){
	  W2max = sval.Atof();
	}
	
	if( skey == "fid_xmin" ){
	  fid_xmin = sval.Atof();
	}

	if( skey == "fid_xmax" ){
	  fid_xmax = sval.Atof();
	}
	
	if( skey == "fid_ymin" ){
	  fid_ymin = sval.Atof();
	}

	if( skey == "fid_ymax" ){
	  fid_ymax = sval.Atof();
	}

	if( skey == "protondeflectionflag" ){
	  protondeflectionflag = sval.Atoi();
	}

	if( skey == "protondeflectionmean" ){
	  pdeflectmean = sval.Atof();
	}
	
	
      }

      currentline_tokens->Delete();
    }
  }

  TFile *fout = new TFile(outfilename,"RECREATE");

  //Histograms to fill:
  TH1D *hW2_globalcut = new TH1D("hW2_globalcut","Global cut and fiducial cut;W^{2} (GeV^{2};",300,-0.5,2.5);
  TH1D *hW2_hcalcut = new TH1D("hW2_hcalcut","Global+fiducial+HCAL;W^{2} (GeV^{2});",300,-0.5,2.5);
  TH1D *hW2_anticut = new TH1D("hW2_anticut","Global+fiducial+!HCAL;W^{2} (GeV^{2});", 300,-0.5,2.5);

  TH1D *hxHCAL_expect_all = new TH1D("hxHCAL_expect_all","Global cut; Expected x_{HCAL} (m);", 250,-3.25,1.75);
  TH1D *hxHCAL_expect_cut = new TH1D("hxHCAL_expect_cut","Global+HCAL cuts; Expected x_{HCAL} (m);", 250,-3.25,1.75);
  TH1D *hxHCAL_expect_acut = new TH1D("hxHCAL_expect_acut","Global + !HCAL cuts; Expected x_{HCA} (m);", 250,-3.25,1.75);
 
  TH1D *hyHCAL_expect_all = new TH1D("hyHCAL_expect_all","Global cut; Expected y_{HCAL} (m);", 250, -1.25,1.25);
  TH1D *hyHCAL_expect_cut = new TH1D("hyHCAL_expect_cut","Global+HCAL cuts; Expected y_{HCAL} (m);", 250,-1.25,1.25);
  TH1D *hyHCAL_expect_acut = new TH1D("hyHCAL_expect_acut","Global+HCAL cuts; Expected y_{HCAL} (m);", 250,-1.25,1.25);
  
  TH2D *hxyHCAL_expect_all = new TH2D("hxyHCAL_expect_all","Global cut; Expected y_{HCAL} (m); Expected x_{HCAL} (m)", 125,-1.25,1.25,125,-3.25,1.75);
  TH2D *hxyHCAL_expect_cut = new TH2D("hxyHCAL_expect_cut","Global+HCAL cuts; Expected y_{HCAL} (m); Expected x_{HCAL} (m)", 125,-1.25,1.25,125,-3.25,1.75);
  TH2D *hxyHCAL_expect_acut = new TH2D("hxyHCAL_expect_acut","Global+HCAL cuts; Expected y_{HCAL} (m); Expected x_{HCAL} (m)", 125,-1.25,1.25,125,-3.25,1.75);

  TH1D *hxHCAL_all = new TH1D("hxHCAL_all","Global cut; Observed x_{HCAL} (m);", 500,-3,2);
  TH1D *hyHCAL_all = new TH1D("hyHCAL_all","Global cut; Observed y_{HCAL} (m);", 250,-1.25,1.25);  
  TH2D *hxyHCAL_all = new TH2D("hxyHCAL_all", "Global cut; y_{HCAL} (m); x_{HCAL} (m)", 125,-1.25,1.25,250,-3,2);
 
  TH1D *hxHCAL_cut = new TH1D("hxHCAL_cut","Global+HCAL cut; Observed x_{HCAL} (m);", 500,-3,2);
  TH1D *hyHCAL_cut = new TH1D("hyHCAL_cut","Global+HCAL cut; Observed y_{HCAL} (m);", 250,-1.25,1.25);  
  TH2D *hxyHCAL_cut = new TH2D("hxyHCAL_cut", "Global+HCAL cut; y_{HCAL} (m); x_{HCAL} (m)", 125,-1.25,1.25,250,-3,2);

  TH1D *hxHCAL_acut = new TH1D("hxHCAL_acut","Global+!HCAL cut; Observed x_{HCAL} (m);", 500,-3,2);
  TH1D *hyHCAL_acut = new TH1D("hyHCAL_acut","Global+!HCAL cut; Observed y_{HCAL} (m);", 250,-1.25,1.25);  
  TH2D *hxyHCAL_acut = new TH2D("hxyHCAL_acut", "Global+!HCAL cut; y_{HCAL} (m); x_{HCAL} (m)", 125,-1.25,1.25,250,-3,2);

  TH1D *hrowHCAL_all = new TH1D("hrowHCAL_all", "Global cut; HCAL row; ", 24,-0.5,23.5);
  TH1D *hcolHCAL_all = new TH1D("hcolHCAL_all", "Global cut; HCAL col; ", 12,-0.5,11.5);
  TH2D *hrowcolHCAL_all = new TH2D("hrowcolHCAL_all", "Global cut; HCAL col; HCAL row", 12,-0.5,11.5,24,-0.5,23.5);

  TH1D *hrowHCAL_cut = new TH1D("hrowHCAL_cut", "Global+HCAL cut; HCAL row; ", 24,-0.5,23.5);
  TH1D *hcolHCAL_cut = new TH1D("hcolHCAL_cut", "Global+HCAL cut; HCAL col; ", 12,-0.5,11.5);
  TH2D *hrowcolHCAL_cut = new TH2D("hrowcolHCAL_cut", "Global+HCAL cut; HCAL col; HCAL row", 12,-0.5,11.5,24,-0.5,23.5);

  TH1D *hrowHCAL_acut = new TH1D("hrowHCAL_acut", "Global+!HCAL cut; HCAL row; ", 24,-0.5,23.5);
  TH1D *hcolHCAL_acut = new TH1D("hcolHCAL_acut", "Global+!HCAL cut; HCAL col; ", 12,-0.5,11.5);
  TH2D *hrowcolHCAL_acut = new TH2D("hrowcolHCAL_acut", "Global+!HCAL cut; HCAL col; HCAL row", 12,-0.5,11.5,24,-0.5,23.5);

  TH1D *hdxall = new TH1D("hdxall", "Global cut; #Deltax (m);", 250,-3,2);
  TH1D *hdxcut = new TH1D("hdxcut", "Global +HCAL cuts; #Deltax (m);", 250,-3,2);
  TH1D *hdxacut = new TH1D("hdxacut", "Global +!HCAL cuts; #Deltax (m);", 250,-3,2);

  TH1D *hdyall = new TH1D("hdyall", "Global cut; #Deltay (m);", 250,-1.25,1.25);
  TH1D *hdycut = new TH1D("hdycut", "Global +HCAL cuts; #Deltay (m);", 250,-1.25,1.25);
  TH1D *hdyacut = new TH1D("hdyacut", "Global +!HCAL cuts; #Deltay (m);", 250,-1.25,1.25);
  
  TH2D *hdxdyall = new TH2D("hdxdyall", "Global cut; #Delta y (m); #Delta x (m)", 250,-1.25,1.25,250,-3,2); 
  TH2D *hdxdycut = new TH2D("hdxdycut", "Global + HCAL cut; #Delta y (m); #Delta x (m)", 250,-1.25,1.25,250,-3,2); 
  TH2D *hdxdyacut = new TH2D("hdxdyacut", "Global + !HCAL cut; #Delta y (m); #Delta x (m)", 250,-1.25,1.25,250,-3,2); 

  TH1D *hdtall = new TH1D("hdtall", "Global cut; #Deltat_{ADC} (ns);", 250,-25,25);
  TH1D *hdtcut = new TH1D("hdtcut", "Global +HCAL cuts; #Deltat_{ADC} (ns);", 250,-25,25);
  TH1D *hdtacut = new TH1D("hdtacut", "Global +!HCAL cuts; #Deltat_{ADC} (ns);", 250,-25,25);

  TH1D *hEHCALall = new TH1D("hEHCALall", "Global cut; E_{HCAL} (GeV);", 500,0.0,2.5);
  TH1D *hEHCALcut = new TH1D("hEHCALcut", "Global +HCAL cuts; E_{HCAL} (GeV);", 500, 0.0, 2.5);
  TH1D *hEHCALacut = new TH1D("hEHCALacut", "Global +!HCAL cuts; E_{HCAL} (GeV);", 500, 0.0, 2.5);
  
  TH2D *hEHCAL_vs_idblk = new TH2D("hEHCAL_vs_idblk", "Global + HCAL cuts; HCAL block ID; HCAL energy", 288, -0.5,287.5, 500,0.0,2.5);

  gmn_elastic_tree *T = new gmn_elastic_tree(C);

  //If this doesn't work, we can use "C" instead of "T":
  TTreeFormula *GlobalCut = new TTreeFormula("GlobalCut",globalcut,C);
  TTreeFormula *HCalCut = new TTreeFormula("HCalCut",hcalcut,C);

  long nevent=0;
  long treenum=0, currenttreenum=0;

  
  while( C->GetEntry( nevent ) ){
    currenttreenum = C->GetTreeNumber();

    if( nevent > 0 && nevent % 10000 == 0 ) cout << "nevent = " << nevent << endl;
    
    if( currenttreenum != treenum || nevent == 0 ){

      treenum = currenttreenum;
      GlobalCut->UpdateFormulaLeaves();
      HCalCut->UpdateFormulaLeaves();


      cout << "switching to new file " << C->GetFile()->GetName() << endl;
    }

    bool passed_global_cut = GlobalCut->EvalInstance(0) != 0;
    bool passed_hcalcut = HCalCut->EvalInstance(0) != 0;

    double protondeflection = 0.0;

    if( protondeflectionflag == 1 ){
      protondeflection = T->protondeflection;
    } else if ( protondeflectionflag == 2 ){
      protondeflection = pdeflectmean;
    }
    
    bool fidcutx = T->xHCAL_expect - protondeflection >= fid_xmin && T->xHCAL_expect - protondeflection <= fid_xmax;
    bool fidcuty = T->yHCAL_expect >= fid_ymin && T->yHCAL_expect <= fid_ymax;

    if( passed_global_cut ){
  
      if( fidcutx && fidcuty ){
	hW2_globalcut->Fill( T->W2 );
	if( passed_hcalcut ){
	  hW2_hcalcut->Fill( T->W2 );
	} else {
	  hW2_anticut->Fill( T->W2 );
	}
      }

      if( T->W2 >= W2min && T->W2 <= W2max ){

	if( fidcuty ){	
	  hxHCAL_expect_all->Fill( T->xHCAL_expect - protondeflection );
	  hxHCAL_all->Fill( T->xHCAL );
	  hrowHCAL_all->Fill( T->rowblkHCAL[0] );
	  if( passed_hcalcut ){
	    hxHCAL_expect_cut->Fill( T->xHCAL_expect - protondeflection );
	    hxHCAL_cut->Fill( T->xHCAL );
	    hrowHCAL_cut->Fill( T->rowblkHCAL[0] );
	  } else {
	    hxHCAL_expect_acut->Fill( T->xHCAL_expect - protondeflection );
	    hxHCAL_acut->Fill( T->xHCAL );
	    hrowHCAL_acut->Fill( T->rowblkHCAL[0] );
	  }
	}
	
	if( fidcutx ){
	  hyHCAL_expect_all->Fill( T->yHCAL_expect );
	  hyHCAL_all->Fill( T->yHCAL );
	  hcolHCAL_all->Fill( T->colblkHCAL[0] );
	  if( passed_hcalcut ){
	    hyHCAL_expect_cut->Fill( T->yHCAL_expect );
	    hyHCAL_cut->Fill( T->yHCAL );
	    hcolHCAL_cut->Fill( T->colblkHCAL[0] );
	  } else {
	    hyHCAL_expect_acut->Fill( T->yHCAL_expect );
	    hyHCAL_acut->Fill( T->yHCAL );
	    hcolHCAL_acut->Fill( T->colblkHCAL[0] );
	  }
	}

	hxyHCAL_expect_all->Fill(T->yHCAL_expect, T->xHCAL_expect - protondeflection );
	if( passed_hcalcut ){
	  hxyHCAL_expect_cut->Fill( T->yHCAL_expect, T->xHCAL_expect - protondeflection );
	} else {
	  hxyHCAL_expect_acut->Fill( T->yHCAL_expect, T->xHCAL_expect - protondeflection );
	}
	//Fill xHCAL, yHCAL, and row and column histograms only for events passing fiducial cut:

	if( fidcutx && fidcuty ){
	  hdxall->Fill( T->deltax );
	  hdyall->Fill( T->deltay );
	  hdxdyall->Fill( T->deltay, T->deltax );

	  hxyHCAL_all->Fill( T->yHCAL, T->xHCAL );
	  hrowcolHCAL_all->Fill( T->colblkHCAL[0], T->rowblkHCAL[0] );

	  hdtall->Fill( T->deltat_adc );
	  hEHCALall->Fill( T->EHCAL );

	  if( passed_hcalcut ){
	  
	    hdtcut->Fill( T->deltat_adc );
	    hEHCALcut->Fill( T->EHCAL );

	    hdxcut->Fill( T->deltax );
	    hdycut->Fill( T->deltay );
	    hdxdycut->Fill( T->deltay, T->deltax );

	    hxyHCAL_cut->Fill( T->yHCAL, T->xHCAL );
	    hrowcolHCAL_cut->Fill( T->colblkHCAL[0], T->rowblkHCAL[0] );

	    hEHCAL_vs_idblk->Fill( T->idblkHCAL[0], T->EHCAL );

	  } else {
	  
	    hdtacut->Fill( T->deltat_adc );
	    hEHCALacut->Fill( T->EHCAL );

	    hdxacut->Fill( T->deltax );
	    hdyacut->Fill( T->deltay );
	    hdxdyacut->Fill( T->deltay, T->deltax );

	    hxyHCAL_acut->Fill( T->yHCAL, T->xHCAL );
	    hrowcolHCAL_acut->Fill( T->colblkHCAL[0], T->rowblkHCAL[0] );
	  }
	}
      }
    }

    
    nevent++;

  }  
  
  TH1D *heff_vs_xexpect = new TH1D( *hxHCAL_expect_cut );
  heff_vs_xexpect->SetName( "heff_vs_xexpect" );
  heff_vs_xexpect->Divide( hxHCAL_expect_cut, hxHCAL_expect_all );

  //Calculate correct error for efficiency determination:
  
  for( int i=1; i<=heff_vs_xexpect->GetNbinsX(); i++ ){
    double eff = heff_vs_xexpect->GetBinContent(i);
    double N = std::max(1.0,hxHCAL_expect_all->GetBinContent(i));
    heff_vs_xexpect->SetBinError(i,sqrt(eff*(1.0-eff)/N));
  }

  TH1D *heff_vs_yexpect = new TH1D( *hyHCAL_expect_cut );
  heff_vs_yexpect->SetName( "heff_vs_yexpect" );
  heff_vs_yexpect->Divide( hyHCAL_expect_cut, hyHCAL_expect_all );

  for( int i=1; i<=heff_vs_yexpect->GetNbinsX(); i++ ){
    double eff = heff_vs_yexpect->GetBinContent(i);
    //prevent divide-by-zero errors
    double N = std::max(1.0,hyHCAL_expect_all->GetBinContent(i));
    heff_vs_yexpect->SetBinError(i,sqrt(eff*(1.0-eff)/N));
  }
  
  TH2D *heff_vs_xyexpect = new TH2D( *hxyHCAL_expect_cut );
  heff_vs_xyexpect->SetName("heff_vs_xyexpect" );
  heff_vs_xyexpect->Divide( hxyHCAL_expect_cut, hxyHCAL_expect_all );

  for( int i=1; i<=heff_vs_xyexpect->GetNbinsX(); i++ ){
    for( int j=1; j<=heff_vs_xyexpect->GetNbinsY(); j++ ){

      int bin = heff_vs_xyexpect->GetBin(i,j);
      
      double eff = heff_vs_xyexpect->GetBinContent(bin);
      double N = std::max(1.0,hxyHCAL_expect_all->GetBinContent(bin));
      heff_vs_xyexpect->SetBinError(bin,sqrt(eff*(1.0-eff)/N));
    }
  }

  TH1D *heff_vs_x = new TH1D( *hxHCAL_cut );
  heff_vs_x->SetName("heff_vs_x");
  heff_vs_x->Divide( hxHCAL_cut, hxHCAL_all );

  TH1D *heff_vs_y = new TH1D( *hyHCAL_cut );
  heff_vs_y->SetName("heff_vs_y");
  heff_vs_y->Divide( hyHCAL_cut, hyHCAL_all );

  TH2D *heff_vs_xy = new TH2D( *hxyHCAL_cut );
  heff_vs_xy->SetName("heff_vs_xy" );
  heff_vs_xy->Divide( hxyHCAL_cut, hxyHCAL_all );
  
  TH1D *heff_vs_row = new TH1D( *hrowHCAL_cut );
  heff_vs_row->SetName("heff_vs_row");
  heff_vs_row->Divide( hrowHCAL_cut, hrowHCAL_all );

  TH1D *heff_vs_col = new TH1D( *hcolHCAL_cut );
  heff_vs_col->SetName("heff_vs_col");
  heff_vs_col->Divide( hcolHCAL_cut, hcolHCAL_all );
  
  TH2D *heff_vs_rowcol = new TH2D( *hrowcolHCAL_cut );
  heff_vs_rowcol->SetName( "heff_vs_rowcol" );
  heff_vs_rowcol->Divide( hrowcolHCAL_cut, hrowcolHCAL_all );


  TH1D *heff_vs_W2 = new TH1D( *hW2_hcalcut );
  heff_vs_W2->SetName("heff_vs_W2");
  heff_vs_W2->Divide( hW2_hcalcut, hW2_globalcut );

  for( int i=1; i<=heff_vs_W2->GetNbinsX(); i++ ){
    double eff = heff_vs_W2->GetBinContent(i);
    double N = std::max(1.0,hW2_globalcut->GetBinContent(i));
    heff_vs_W2->SetBinError(i,sqrt(eff*(1.0-eff)/N));
  }
  
  fout->Write();

  
}
