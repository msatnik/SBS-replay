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
#include "TStyle.h"
#include <iostream>
#include <fstream>

void HCAL_np_ratio_uniformity(const char *setupfilename, const char *outfilename){
  TH1::SetDefaultSumw2();

  gStyle->SetMarkerStyle(20);
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
  // three is a "good hcal" cut (dx, dy, possibly also dt, EHCAL) for proton and neutron:
  // four is a set of "fiducial" cuts in x and y 
  // since cuts 2 and 4 are generally going to be rectangular in nature, we can specify lower and upper limits for these via command-line arguments:

  TCut globalcut = "";
  
  while(currentline.ReadLine(setupfile) && !currentline.BeginsWith("endglobalcut")){
    if( !currentline.BeginsWith("#")) globalcut += currentline.Data();
  }

  TCut hcalcut_all = "";
  TCut hcalcut_proton = "";
  TCut hcalcut_neutron = "";

  //This cut is intended for HCAL energy threshold and (loose) coincidence time
  while(currentline.ReadLine(setupfile) && !currentline.BeginsWith("endhcalcut") ){
    if( !currentline.BeginsWith("#") ) hcalcut_all += currentline.Data();
  }
  
  while(currentline.ReadLine(setupfile) && !currentline.BeginsWith("endprotoncut")){
    if( !currentline.BeginsWith("#")) hcalcut_proton += currentline.Data();
  }

  while(currentline.ReadLine(setupfile) && !currentline.BeginsWith("endneutroncut")){
    if( !currentline.BeginsWith("#")) hcalcut_neutron += currentline.Data();
  }
  

  //provide default values for these that can be overridden by key-value pairs in the "config" section:
  //Note that the default "fiducial" limits correspond to the neutron case; i.e., no deflection:
  double W2min=0.65, W2max=0.95, fid_xmin=-2.5, fid_xmax=1.0, fid_ymin=-0.6, fid_ymax=0.6;
  double thetapq_cut_p=0.05, thetapq_cut_n=0.05;
  
  int protondeflectionflag = 0; //1 = use protondeflection from tree, 2 = use mean from config file:

  double pdeflectmean=0.0;

  double dxmin=-0.4,dxmax=0.4,dymin=-0.5,dymax=0.5;
  double dxshift_np = 0.0, dyshift_np = 0.0;
  
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

	if( skey == "thetapq_cut_p" ){
	  thetapq_cut_p = sval.Atof();
	}

	if( skey == "thetapq_cut_n" ){
	  thetapq_cut_n = sval.Atof();
	}

	if( skey == "dxmin" ){
	  dxmin = sval.Atof();
	}

	if( skey == "dxmax" ){
	  dxmax = sval.Atof();
	}

	if( skey == "dymin" ){
	  dymin = sval.Atof();
	}

	if( skey == "dymax" ){
	  dymax = sval.Atof();
	}

	if( skey == "dxshift_np" ){
	  dxshift_np = sval.Atof();
	}

	if( skey == "dyshift_np" ){
	  dyshift_np = sval.Atof();
	}
	
      }

      currentline_tokens->Delete();
    }
  }

  TFile *fout = new TFile(outfilename,"RECREATE");

  //Histograms to fill:
  TH1D *hW2_globalcut = new TH1D("hW2_globalcut","Global cut and fiducial cut;W^{2} (GeV^{2};",300,-0.5,2.5);
  TH1D *hW2_hcalcut_p = new TH1D("hW2_hcalcut_p","Global+fiducial+pcut;W^{2} (GeV^{2});",300,-0.5,2.5);
  TH1D *hW2_anticut_p = new TH1D("hW2_anticut_p","Global+fiducial+!pcut;W^{2} (GeV^{2});", 300,-0.5,2.5);
  TH1D *hW2_hcalcut_n = new TH1D("hW2_hcalcut_n","Global+fiducial+ncut;W^{2} (GeV^{2});",300,-0.5,2.5);
  TH1D *hW2_anticut_n = new TH1D("hW2_anticut_n","Global+fiducial+!ncut;W^{2} (GeV^{2});", 300,-0.5,2.5);

  TH1D *hxHCAL_expect_all_p = new TH1D("hxHCAL_expect_all_p","Global cut; Expected proton x_{HCAL} (m);", 250,-3.25,1.75);
  TH1D *hxHCAL_expect_cut_p = new TH1D("hxHCAL_expect_cut_p","Global+proton cuts; Expected x_{HCAL} (m);", 250,-3.25,1.75);
  TH1D *hxHCAL_expect_acut_p = new TH1D("hxHCAL_expect_acut_p","Global + !proton cuts; Expected x_{HCAL} (m);", 250,-3.25,1.75);

  TH1D *hxHCAL_nexpect_cut_p = new TH1D("hxHCAL_nexpect_cut_p", "Global + proton cuts; Expected neutron x_{HCAL} (m);", 250,-3.25,1.75);
  
  TH1D *hxHCAL_expect_all_n = new TH1D("hxHCAL_expect_all_n","Global cut; Expected neutron x_{HCAL} (m);", 250,-3.25,1.75);
  TH1D *hxHCAL_expect_cut_n = new TH1D("hxHCAL_expect_cut_n","Global+neutron cuts; Expected x_{HCAL} (m);", 250,-3.25,1.75);
  TH1D *hxHCAL_expect_acut_n = new TH1D("hxHCAL_expect_acut_n","Global + !neutron cuts; Expected x_{HCAL} (m);", 250,-3.25,1.75);
 
  TH1D *hyHCAL_expect_all = new TH1D("hyHCAL_expect_all","Global cut; Expected y_{HCAL} (m);", 250, -1.25,1.25);
  TH1D *hyHCAL_expect_cut_p = new TH1D("hyHCAL_expect_cut_p","Global+proton cuts; Expected y_{HCAL} (m);", 250,-1.25,1.25);
  TH1D *hyHCAL_expect_acut_p = new TH1D("hyHCAL_expect_acut_p","Global+!proton cuts; Expected y_{HCAL} (m);", 250,-1.25,1.25);

  TH1D *hyHCAL_expect_cut_n = new TH1D("hyHCAL_expect_cut_n","Global+neutron cuts; Expected y_{HCAL} (m);", 250,-1.25,1.25);
  TH1D *hyHCAL_expect_acut_n = new TH1D("hyHCAL_expect_acut_n","Global+!neutron cuts; Expected y_{HCAL} (m);", 250,-1.25,1.25);
  
  TH2D *hxyHCAL_expect_all_p = new TH2D("hxyHCAL_expect_all_p","Global cut; Expected y_{HCAL} (m); Expected proton x_{HCAL} (m)", 125,-1.25,1.25,125,-3.25,1.75);
  TH2D *hxyHCAL_expect_all_n = new TH2D("hxyHCAL_expect_all_n","Global cut; Expected y_{HCAL} (m); Expected neutron x_{HCAL} (m)", 125,-1.25,1.25,125,-3.25,1.75);
  TH2D *hxyHCAL_expect_cut_p = new TH2D("hxyHCAL_expect_cut_p","Global+proton cuts; Expected y_{HCAL} (m); Expected x_{HCAL} (m)", 125,-1.25,1.25,125,-3.25,1.75);
  TH2D *hxyHCAL_expect_acut_p = new TH2D("hxyHCAL_expect_acut_p","Global+!proton cuts; Expected y_{HCAL} (m); Expected x_{HCAL} (m)", 125,-1.25,1.25,125,-3.25,1.75);

  TH2D *hxyHCAL_nexpect_cut_p = new TH2D("hxyHCAL_nexpect_cut_p","Global + proton cuts; Expected y_{HCAL} (m); Expected x_{HCAL} (m)", 125,-1.25,1.25,125,-3.25,1.75);
  
  TH2D *hxyHCAL_expect_cut_n = new TH2D("hxyHCAL_expect_cut_n","Global+neutron cuts; Expected y_{HCAL} (m); Expected x_{HCAL} (m)", 125,-1.25,1.25,125,-3.25,1.75);
  TH2D *hxyHCAL_expect_acut_n = new TH2D("hxyHCAL_expect_acut_n","Global+!neutron cuts; Expected y_{HCAL} (m); Expected x_{HCAL} (m)", 125,-1.25,1.25,125,-3.25,1.75);

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
  TH1D *hdxcut_p = new TH1D("hdxcut_p", "Global +proton cuts; #Deltax (m);", 250,-3,2);
  TH1D *hdxacut_p = new TH1D("hdxacut_p", "Global +!proton cuts; #Deltax (m);", 250,-3,2);
  TH1D *hdxcut_n = new TH1D("hdxcut_n", "Global +neutron cuts; #Deltax (m);", 250,-3,2);
  TH1D *hdxacut_n = new TH1D("hdxacut_n", "Global +!neutron cuts; #Deltax (m);", 250,-3,2);

  TH1D *hdyall = new TH1D("hdyall", "Global cut; #Deltay (m);", 250,-1.25,1.25);
  TH1D *hdycut_p = new TH1D("hdycut_p", "Global +proton cuts; #Deltay (m);", 250,-1.25,1.25);
  TH1D *hdyacut_p = new TH1D("hdyacut_p", "Global +!proton cuts; #Deltay (m);", 250,-1.25,1.25);
  TH1D *hdycut_n = new TH1D("hdycut_n", "Global +neutron cuts; #Deltay (m);", 250,-1.25,1.25);
  TH1D *hdyacut_n = new TH1D("hdyacut_n", "Global +!neutron cuts; #Deltay (m);", 250,-1.25,1.25);
  
  TH2D *hdxdyall = new TH2D("hdxdyall", "Global cut; #Delta y (m); #Delta x (m)", 250,-1.25,1.25,250,-3,2); 
  TH2D *hdxdycut_p = new TH2D("hdxdycut_p", "Global + proton cut; #Delta y (m); #Delta x (m)", 250,-1.25,1.25,250,-3,2); 
  TH2D *hdxdyacut_p = new TH2D("hdxdyacut_p", "Global + !proton cut; #Delta y (m); #Delta x (m)", 250,-1.25,1.25,250,-3,2);
  TH2D *hdxdycut_n = new TH2D("hdxdycut_n", "Global + neutron cut; #Delta y (m); #Delta x (m)", 250,-1.25,1.25,250,-3,2); 
  TH2D *hdxdyacut_n = new TH2D("hdxdyacut_n", "Global + !neutron cut; #Delta y (m); #Delta x (m)", 250,-1.25,1.25,250,-3,2); 

  TH1D *hdtall = new TH1D("hdtall", "Global cut; #Deltat_{ADC} (ns);", 250,-25,25);
  TH1D *hdtcut_p = new TH1D("hdtcut_p", "Global +proton cuts; #Deltat_{ADC} (ns);", 250,-25,25);
  TH1D *hdtacut_p = new TH1D("hdtacut_p", "Global +!proton cuts; #Deltat_{ADC} (ns);", 250,-25,25);
  TH1D *hdtcut_n = new TH1D("hdtcut_n", "Global +neutron cuts; #Deltat_{ADC} (ns);", 250,-25,25);
  TH1D *hdtacut_n = new TH1D("hdtacut_n", "Global +!neutron cuts; #Deltat_{ADC} (ns);", 250,-25,25);

  TH1D *hEHCALall = new TH1D("hEHCALall", "Global cut; E_{HCAL} (GeV);", 500,0.0,2.5);
  TH1D *hEHCALcut_p = new TH1D("hEHCALcut_p", "Global +proton cuts; E_{HCAL} (GeV);", 500, 0.0, 2.5);
  TH1D *hEHCALacut_p = new TH1D("hEHCALacut_p", "Global +!proton cuts; E_{HCAL} (GeV);", 500, 0.0, 2.5);
  TH1D *hEHCALcut_n = new TH1D("hEHCALcut_n", "Global +neutron cuts; E_{HCAL} (GeV);", 500, 0.0, 2.5);
  TH1D *hEHCALacut_n = new TH1D("hEHCALacut_n", "Global +!neutron cuts; E_{HCAL} (GeV);", 500, 0.0, 2.5);
  
  TH2D *hEHCAL_vs_idblk = new TH2D("hEHCAL_vs_idblk", "Global + HCAL cuts; HCAL block ID; HCAL energy", 288, -0.5,287.5, 500,0.0,2.5);

  TH1D *hxnexpect_all_fid = new TH1D("hxnexpect_all_fid","Global+fiducial+W^{2} cuts; Expected x_{neutron} (m);",250,-3.25,1.75);
  TH1D *hxnexpect_cut_fid_p = new TH1D("hxnexpect_cut_fid_p","Global+fid+W^{2} + p cuts; Expected x_{neutron} (m);",250,-3.25,1.75);
  TH1D *hxnexpect_cut_fid_n = new TH1D("hxnexpect_cut_fid_n","Global+fid+W^{2} + n cuts; Expected x_{neutron} (m);",250,-3.25,1.75);
  
  TH1D *hynexpect_all_fid = new TH1D("hynexpect_all_fid","Global+fiducial+W^{2} cuts; Expected y_{neutron} (m);",250,-1.25,1.25);
  TH1D *hynexpect_cut_fid_p = new TH1D("hynexpect_cut_fid_p","Global+fid+W^{2} + p cuts; Expected y_{neutron} (m);",250,-1.25,1.25);
  TH1D *hynexpect_cut_fid_n = new TH1D("hynexpect_cut_fid_n","Global+fid+W^{2} + n cuts; Expected y_{neutron} (m);",250,-1.25,1.25);

  TH2D *hxynexpect_all_fid = new TH2D("hxynexpect_all_fid","Global+fid+W^{2} cuts; Expected y_{n} (m); Expected x_{n} (m)",125,-1.25,1.25,125,-3.25,1.75);
  TH2D *hxynexpect_cut_fid_p = new TH2D("hxynexpect_cut_fid_p","Global+fid+W^{2} + p cuts; Expected y_{n} (m); Expected x_{n} (m)",125,-1.25,1.25,125,-3.25,1.75);
  TH2D *hxynexpect_cut_fid_n = new TH2D("hxynexpect_cut_fid_n","Global+fid+W^{2} + n cuts; Expected y_{n} (m); Expected x_{n} (m)",125,-1.25,1.25,125,-3.25,1.75);

  TH2D *hW2_thetapq_all_n = new TH2D("hW2_thetapq_all_n","Global cuts; #theta_{pq}^{n} (rad);W^{2} (GeV^{2})",150,0,0.3,150,-0.5,2.5);

  TH2D *hW2_thetapq_all_p = new TH2D("hW2_thetapq_all_p","Global cuts; #theta_{pq}^{p} (rad);W^{2} (GeV^{2})",150,0,0.3,150,-0.5,2.5);

  TH2D *hW2_thetapq_fid_n = new TH2D("hW2_thetapq_fid_n","Global+fiducial cuts; #theta_{pq}^{n} (rad);W^{2} (GeV^{2})",150,0,0.3,150,-0.5,2.5);

  TH2D *hW2_thetapq_fid_p = new TH2D("hW2_thetapq_fid_p","Global+fiducial cuts; #theta_{pq}^{p} (rad);W^{2} (GeV^{2})",150,0,0.3,150,-0.5,2.5);

  // What other histograms do we want?
  // We want to compare thetapq for p and n (with fid. and W^2/etc cuts)
  // We also want to look at dx with all OTHER cuts (incl. dy) and dy with all OTHER cuts (including dx)
  // We also want to look at ratios for dx for neutron and dx+protondeflection for proton (with all other cuts)

  TH1D *hthetapq_n = new TH1D("hthetapq_n","Global + fid. + W^2 cuts; #theta_{pq} (n hypothesis)",150,0,0.3);
  TH1D *hthetapq_p = new TH1D("hthetapq_p","Global + fid. + W^2 cuts; #theta_{pq} (p hypothesis)",150,0,0.3);

  TH1D *hdxn_allothercuts = new TH1D("hdxn_allothercuts","Global + fid. + W^2 + #Deltay cuts; #Deltax (n hypothesis)",300,-3,3);
  TH1D *hdxp_allothercuts = new TH1D("hdxp_allothercuts","Global + fid. + W^2 + #Deltay cuts; #Deltax (p hypothesis)",300,-3,3);

  TH1D *hdyn_allothercuts = new TH1D("hdyn_allothercuts","Global + fid. + W^2 + #Deltax(neutron) cuts; #Deltay (m)",300,-1.5,1.5); 
  TH1D *hdyp_allothercuts = new TH1D("hdyp_allothercuts","Global + fid. + W^2 + #Deltax(proton) cuts; #Deltay (m)",300,-1.5,1.5); 
  
  gmn_elastic_tree *T = new gmn_elastic_tree(C);

  //If this doesn't work, we can use "C" instead of "T":
  TTreeFormula *GlobalCut = new TTreeFormula("GlobalCut",globalcut,C);
  TTreeFormula *HCalCut = new TTreeFormula("HCalCut", hcalcut_all,C);
  TTreeFormula *HCalCutP = new TTreeFormula("HCalCutP",hcalcut_proton,C);
  TTreeFormula *HCalCutN = new TTreeFormula("HCalCutN",hcalcut_neutron,C);
  
  long nevent=0;
  long treenum=0, currenttreenum=0;

  
  while( C->GetEntry( nevent ) ){
    currenttreenum = C->GetTreeNumber();

    if( nevent > 0 && nevent % 10000 == 0 ) cout << "nevent = " << nevent << endl;
    
    if( currenttreenum != treenum || nevent == 0 ){

      treenum = currenttreenum;
      GlobalCut->UpdateFormulaLeaves();
      HCalCut->UpdateFormulaLeaves();
      HCalCutP->UpdateFormulaLeaves();
      HCalCutN->UpdateFormulaLeaves();

      cout << "switching to new file " << C->GetFile()->GetName() << endl;
    }

    bool passed_global_cut = GlobalCut->EvalInstance(0) != 0;
    bool passed_hcalcut = HCalCut->EvalInstance(0) != 0;
    bool passed_hcalcut_p = HCalCutP->EvalInstance(0) != 0;
    bool passed_hcalcut_n = HCalCutN->EvalInstance(0) != 0;
    
    double protondeflection = 0.0;

    if( protondeflectionflag == 1 ){
      protondeflection = T->protondeflection;
    } else if ( protondeflectionflag == 2 ){
      protondeflection = pdeflectmean;
    }

    //Now here we must apply the fiducial cut for both proton and neutron simultaneously:
    
    bool fidcutx_p = T->xHCAL_expect_4vect - protondeflection >= fid_xmin && T->xHCAL_expect_4vect - protondeflection <= fid_xmax;
    bool fidcuty_p = T->yHCAL_expect_4vect >= fid_ymin && T->yHCAL_expect_4vect <= fid_ymax;

    bool fidcutx_n = T->xHCAL_expect_4vect >= fid_xmin && T->xHCAL_expect_4vect <= fid_xmax;
    bool fidcuty_n = fidcuty_p;

    bool fidcutx = fidcutx_p && fidcutx_n;
    bool fidcuty = fidcuty_p && fidcuty_n;

    
    
    if( passed_global_cut && passed_hcalcut ){

      hW2_thetapq_all_n->Fill( T->thetapq_n, T->W2 );
      hW2_thetapq_all_p->Fill( T->thetapq_p, T->W2 );
      
      if( fidcutx && fidcuty ){
	hW2_globalcut->Fill( T->W2 );
	if( passed_hcalcut_p ){
	  hW2_hcalcut_p->Fill( T->W2 );
	} else if( !passed_hcalcut_n ){
	  hW2_anticut_p->Fill( T->W2 );
	}

	hW2_thetapq_fid_n->Fill( T->thetapq_n, T->W2 );
	hW2_thetapq_fid_p->Fill( T->thetapq_p, T->W2 );
	
	if( passed_hcalcut_n ){
	  hW2_hcalcut_n->Fill( T->W2 );
	} else if( !passed_hcalcut_p ){
	  hW2_anticut_n->Fill( T->W2 );
	}

	if( T->W2 >= W2min && T->W2 <= W2max ){
	  hxnexpect_all_fid->Fill( T->xHCAL_expect_4vect );
	  hynexpect_all_fid->Fill( T->yHCAL_expect_4vect );
	  hxynexpect_all_fid->Fill( T->yHCAL_expect_4vect, T->xHCAL_expect_4vect );
	  if( passed_hcalcut_n ){
	    hxnexpect_cut_fid_n->Fill( T->xHCAL_expect_4vect );
	    hynexpect_cut_fid_n->Fill( T->yHCAL_expect_4vect );
	    hxynexpect_cut_fid_n->Fill( T->yHCAL_expect_4vect, T->xHCAL_expect_4vect );
	  }

	  if( passed_hcalcut_p ){
	    hxnexpect_cut_fid_p->Fill( T->xHCAL_expect_4vect );
	    hynexpect_cut_fid_p->Fill( T->yHCAL_expect_4vect );
	    hxynexpect_cut_fid_p->Fill( T->yHCAL_expect_4vect, T->xHCAL_expect_4vect );
	  }
	}
	
      }

      if( T->W2 >= W2min && T->W2 <= W2max ){

	if( fidcuty ){	
	  hxHCAL_expect_all_p->Fill( T->xHCAL_expect_4vect - protondeflection );
	  hxHCAL_expect_all_n->Fill( T->xHCAL_expect_4vect );
	  hxHCAL_all->Fill( T->xHCAL );
	  hrowHCAL_all->Fill( T->rowblkHCAL[0] );
	  if( passed_hcalcut_p ){
	    hxHCAL_expect_cut_p->Fill( T->xHCAL_expect_4vect - protondeflection );
	    hxHCAL_nexpect_cut_p->Fill( T->xHCAL_expect_4vect );
	    hxHCAL_cut->Fill( T->xHCAL );
	    hrowHCAL_cut->Fill( T->rowblkHCAL[0] );
	  } else if( !passed_hcalcut_n ){
	    hxHCAL_expect_acut_p->Fill( T->xHCAL_expect_4vect - protondeflection );
	    hxHCAL_acut->Fill( T->xHCAL );
	    hrowHCAL_acut->Fill( T->rowblkHCAL[0] );
	  }

	  if( passed_hcalcut_n ){
	    hxHCAL_expect_cut_n->Fill( T->xHCAL_expect_4vect );
	    hxHCAL_cut->Fill( T->xHCAL );
	    hrowHCAL_cut->Fill( T->rowblkHCAL[0] );
	  } else if( !passed_hcalcut_p ){
	    hxHCAL_expect_acut_n->Fill( T->xHCAL_expect_4vect );
	    hxHCAL_acut->Fill( T->xHCAL );
	    hrowHCAL_acut->Fill( T->rowblkHCAL[0] );
	  }
	}
	
	if( fidcutx ){
	  hyHCAL_expect_all->Fill( T->yHCAL_expect_4vect );
	  hyHCAL_all->Fill( T->yHCAL );
	  hcolHCAL_all->Fill( T->colblkHCAL[0] );
	  if( passed_hcalcut_p ){
	    hyHCAL_expect_cut_p->Fill( T->yHCAL_expect_4vect );
	    hyHCAL_cut->Fill( T->yHCAL );
	    hcolHCAL_cut->Fill( T->colblkHCAL[0] );
	  } else if( !passed_hcalcut_n ) {
	    hyHCAL_expect_acut_p->Fill( T->yHCAL_expect_4vect );
	    hyHCAL_acut->Fill( T->yHCAL );
	    hcolHCAL_acut->Fill( T->colblkHCAL[0] );
	  }

	  if( passed_hcalcut_n ){
	    hyHCAL_expect_cut_n->Fill( T->yHCAL_expect_4vect );
	    hyHCAL_cut->Fill( T->yHCAL );
	    hcolHCAL_cut->Fill( T->colblkHCAL[0] );
	  } else if( !passed_hcalcut_p ) {
	    hyHCAL_expect_acut_n->Fill( T->yHCAL_expect_4vect );
	    hyHCAL_acut->Fill( T->yHCAL );
	    hcolHCAL_acut->Fill( T->colblkHCAL[0] );
	  }
	}

	hxyHCAL_expect_all_p->Fill(T->yHCAL_expect_4vect, T->xHCAL_expect_4vect - protondeflection );
	hxyHCAL_expect_all_n->Fill(T->yHCAL_expect_4vect, T->xHCAL_expect_4vect );
	if( passed_hcalcut_p ){
	  hxyHCAL_expect_cut_p->Fill( T->yHCAL_expect_4vect, T->xHCAL_expect_4vect - protondeflection );
	  hxyHCAL_nexpect_cut_p->Fill( T->yHCAL_expect_4vect, T->xHCAL_expect_4vect );
	} else if( !passed_hcalcut_n ){
	  hxyHCAL_expect_acut_p->Fill( T->yHCAL_expect_4vect, T->xHCAL_expect_4vect - protondeflection );
	}

	if( passed_hcalcut_n ){
	  hxyHCAL_expect_cut_n->Fill( T->yHCAL_expect_4vect, T->xHCAL_expect_4vect );
	} else if( !passed_hcalcut_p ){
	  hxyHCAL_expect_acut_n->Fill( T->yHCAL_expect_4vect, T->xHCAL_expect_4vect );
	}
	
	//Fill xHCAL, yHCAL, and row and column histograms only for events passing fiducial cut:
	
	if( fidcutx && fidcuty ){

	  hthetapq_n->Fill( T->thetapq_n );
	  hthetapq_p->Fill( T->thetapq_p );

	  if( dymin <= T->deltay_4vect && T->deltay_4vect <= dymax ){
	    hdxn_allothercuts->Fill( T->deltax_4vect );
	    hdxp_allothercuts->Fill( T->deltax_4vect + protondeflection - dxshift_np );
	  }

	  if( dxmin <= T->deltax_4vect && T->deltax_4vect <= dxmax ){
	    hdyn_allothercuts->Fill( T->deltay_4vect );
	  }
	  if( dxmin <= T->deltax_4vect + protondeflection && T->deltax_4vect + protondeflection <= dxmax ){
	    hdyp_allothercuts->Fill( T->deltay_4vect );
	  }
	  
	  hdxall->Fill( T->deltax_4vect );
	  hdyall->Fill( T->deltay_4vect );
	  hdxdyall->Fill( T->deltay_4vect, T->deltax_4vect );

	  hxyHCAL_all->Fill( T->yHCAL, T->xHCAL );
	  hrowcolHCAL_all->Fill( T->colblkHCAL[0], T->rowblkHCAL[0] );

	  hdtall->Fill( T->deltat_adc );
	  hEHCALall->Fill( T->EHCAL );

	  if( passed_hcalcut_p ){
	  
	    hdtcut_p->Fill( T->deltat_adc );
	    hEHCALcut_p->Fill( T->EHCAL );

	    hdxcut_p->Fill( T->deltax_4vect );
	    hdycut_p->Fill( T->deltay_4vect );
	    hdxdycut_p->Fill( T->deltay_4vect, T->deltax_4vect );

	    hxyHCAL_cut->Fill( T->yHCAL, T->xHCAL );
	    hrowcolHCAL_cut->Fill( T->colblkHCAL[0], T->rowblkHCAL[0] );

	    hEHCAL_vs_idblk->Fill( T->idblkHCAL[0], T->EHCAL );

	  } else if( !passed_hcalcut_n ){
	  
	    hdtacut_p->Fill( T->deltat_adc );
	    hEHCALacut_p->Fill( T->EHCAL );

	    hdxacut_p->Fill( T->deltax_4vect );
	    hdyacut_p->Fill( T->deltay_4vect );
	    hdxdyacut_p->Fill( T->deltay_4vect, T->deltax_4vect );

	    hxyHCAL_acut->Fill( T->yHCAL, T->xHCAL );
	    hrowcolHCAL_acut->Fill( T->colblkHCAL[0], T->rowblkHCAL[0] );
	  }

	  if( passed_hcalcut_n ){
	    hdtcut_n->Fill( T->deltat_adc );
	    hEHCALcut_n->Fill( T->EHCAL );
	    hdxcut_n->Fill( T->deltax_4vect );
	    hdycut_n->Fill( T->deltay_4vect );
	    hdxdycut_n->Fill( T->deltay_4vect, T->deltax_4vect );

	    hxyHCAL_cut->Fill( T->yHCAL, T->xHCAL );
	    hrowcolHCAL_cut->Fill( T->colblkHCAL[0], T->rowblkHCAL[0] );
	    hEHCAL_vs_idblk->Fill( T->idblkHCAL[0], T->EHCAL );
	  } else if( !passed_hcalcut_p ){
	    hdtacut_n->Fill( T->deltat_adc );
	    hEHCALacut_n->Fill( T->EHCAL );
	    hdxacut_n->Fill( T->deltax_4vect );
	    hdyacut_n->Fill( T->deltay_4vect );
	    hdxdyacut_n->Fill( T->deltay_4vect, T->deltax_4vect );

	    hxyHCAL_acut->Fill( T->yHCAL, T->xHCAL );
	    hrowcolHCAL_acut->Fill( T->colblkHCAL[0], T->rowblkHCAL[0] );
	  }
	    
	}
      }
    }

    
    nevent++;

  }  

  TH1D *hnpratio_x = new TH1D( *hxnexpect_cut_fid_n );
  hnpratio_x->SetName( "hnpratio_x" );
  hnpratio_x->SetTitle("; x_{expect} (m); N/P ratio");
  hnpratio_x->Divide( hxnexpect_cut_fid_n, hxnexpect_cut_fid_p );

  TH1D *hnpratio_y = new TH1D( *hynexpect_cut_fid_n );
  hnpratio_y->SetName( "hnpratio_y" );
  hnpratio_y->SetTitle("; y_{expect} (m); N/P ratio");
  hnpratio_y->Divide( hynexpect_cut_fid_n, hynexpect_cut_fid_p );

  TH2D *hnpratio_xy = new TH2D( *hxynexpect_cut_fid_n );
  hnpratio_xy->SetName( "hnpratio_xy" );
  hnpratio_xy->SetTitle( "N/P ratio ; y_{expect} (m); x_{expect} (m)" );
  hnpratio_xy->Divide(hxynexpect_cut_fid_n, hxynexpect_cut_fid_p );

  TH1D *hnpratio_W2 = new TH1D( *hW2_hcalcut_n );
  hnpratio_W2->SetName("hnpratio_W2");
  hnpratio_W2->SetTitle("; W^{2} (GeV^{2});N/P ratio");
  hnpratio_W2->Divide( hW2_hcalcut_n, hW2_hcalcut_p );

  TH1D *hnpratio_EHCAL = new TH1D( *hEHCALcut_n );
  hnpratio_EHCAL->SetName("hnpratio_EHCAL");
  hnpratio_EHCAL->SetTitle("Not necessarily expected to be constant; E_{HCAL} (GeV);N/P ratio");
  hnpratio_EHCAL->Divide( hEHCALcut_n, hEHCALcut_p );
  
  TH1D *hnpratio_x_nofid = new TH1D( *hxHCAL_expect_cut_n );
  hnpratio_x_nofid->SetName("hnpratio_x_nofid" );
  hnpratio_x_nofid->SetTitle( "Global + y fiducial cuts; Expected x_{HCAL} (m); N/P ratio" );
  hnpratio_x_nofid->Divide( hxHCAL_expect_cut_n, hxHCAL_nexpect_cut_p );

  TH1D *hnpratio_y_nofid = new TH1D( *hyHCAL_expect_cut_n );
  hnpratio_y_nofid->SetName("hnpratio_y_nofid" );
  hnpratio_y_nofid->SetTitle( "Global + x fiducial cuts; Expected y_{HCAL} (m); N/P ratio" );
  hnpratio_y_nofid->Divide( hyHCAL_expect_cut_n, hyHCAL_expect_cut_p );

  TH2D *hnpratio_xy_nofid = new TH2D( *hxyHCAL_expect_cut_n );
  hnpratio_xy_nofid->SetName("hnpratio_xy_nofid" );
  hnpratio_xy_nofid->SetTitle( "N/P ratio; Expected y_{HCAL} (m); Expected x_{HCAL} (m)" );
  hnpratio_xy_nofid->Divide( hxyHCAL_expect_cut_n, hxyHCAL_nexpect_cut_p );

  TH1D *hnpratio_thpq = new TH1D(*hthetapq_n);
  hnpratio_thpq->SetName("hnpratio_thpq");
  hnpratio_thpq->SetTitle("Global+fid+W^2+HCAL E,t;#theta_{pq} (rad); N/P ratio");
  hnpratio_thpq->Divide( hthetapq_n, hthetapq_p );

  TH1D *hnpratio_dx = new TH1D(*hdxn_allothercuts);
  hnpratio_dx->SetName("hnpratio_dx");
  hnpratio_dx->SetTitle("all other cuts; #Deltax (m);N/P ratio");
  hnpratio_dx->Divide(hdxn_allothercuts, hdxp_allothercuts);

  TH1D *hnpratio_dy = new TH1D(*hdyn_allothercuts);
  hnpratio_dy->SetName("hnpratio_dy");
  hnpratio_dy->SetTitle("all other cuts; #Deltay (m);N/P ratio");
  hnpratio_dy->Divide(hdyn_allothercuts, hdyp_allothercuts);
  // TH1D *heff_vs_xexpect = new TH1D( *hxHCAL_expect_cut );
  // heff_vs_xexpect->SetName( "heff_vs_xexpect" );
  // heff_vs_xexpect->Divide( hxHCAL_expect_cut, hxHCAL_expect_all );

  // TH1D *heff_vs_yexpect = new TH1D( *hyHCAL_expect_cut );
  // heff_vs_yexpect->SetName( "heff_vs_yexpect" );
  // heff_vs_yexpect->Divide( hyHCAL_expect_cut, hyHCAL_expect_all );

  // TH2D *heff_vs_xyexpect = new TH2D( *hxyHCAL_expect_cut );
  // heff_vs_xyexpect->SetName("heff_vs_xyexpect" );
  // heff_vs_xyexpect->Divide( hxyHCAL_expect_cut, hxyHCAL_expect_all );


  // TH1D *heff_vs_x = new TH1D( *hxHCAL_cut );
  // heff_vs_x->SetName("heff_vs_x");
  // heff_vs_x->Divide( hxHCAL_cut, hxHCAL_all );

  // TH1D *heff_vs_y = new TH1D( *hyHCAL_cut );
  // heff_vs_y->SetName("heff_vs_y");
  // heff_vs_y->Divide( hyHCAL_cut, hyHCAL_all );

  // TH2D *heff_vs_xy = new TH2D( *hxyHCAL_cut );
  // heff_vs_xy->SetName("heff_vs_xy" );
  // heff_vs_xy->Divide( hxyHCAL_cut, hxyHCAL_all );
  
  // TH1D *heff_vs_row = new TH1D( *hrowHCAL_cut );
  // heff_vs_row->SetName("heff_vs_row");
  // heff_vs_row->Divide( hrowHCAL_cut, hrowHCAL_all );

  // TH1D *heff_vs_col = new TH1D( *hcolHCAL_cut );
  // heff_vs_col->SetName("heff_vs_col");
  // heff_vs_col->Divide( hcolHCAL_cut, hcolHCAL_all );
  
  // TH2D *heff_vs_rowcol = new TH2D( *hrowcolHCAL_cut );
  // heff_vs_rowcol->SetName( "heff_vs_rowcol" );
  // heff_vs_rowcol->Divide( hrowcolHCAL_cut, hrowcolHCAL_all );


  fout->Write();

  
}
