void coinBBDiff(){

  cout << "Processing coinBBDiff macro.." << endl;

  TStopwatch *st = new TStopwatch();
  st->Start(kTRUE);

  TPad *pads = (TPad*) gPad;
  pads->Divide(3,2);

  Int_t Ndata_bb_tdctrig_tdcelemID = 0;
  Double_t bb_sh_nclus = 0., sbs_hcal_nclus = 0.;
  Double_t bb_tdctrig_tdc[6] = {0.}, bb_tdctrig_tdcelemID[6] = {0.};
  Double_t hcal_reftdc[4];
  Double_t hcal_reftdc_elemID[4];

  Int_t Ndata_hcal_reftdc;
  
  TH1D *h1_BB = new TH1D("h1_BB","BB Trig time; ns", 1000, 0, 1000);
  TH1D *h1_HCAL = new TH1D("h1_HCAL","HCAL Trig time; ns", 1000, 0, 1000);
  TH1D *h1_calDiff = new TH1D("h1_calDiff","HCal/BBCal Difference; ns", 250, 400, 650);
  TH1D *h1_RF_time = new TH1D("h1_RF_time","Accelerator RF Time (BB); ns", 501, -0.5, 500.5);
  TH1D *h1_RF_time_hcal = new TH1D("h1_RF_time_hcal","Accelerator RF Time (HCAL); ns;",501,-0.5,500.5);
  TH2D *h2_RF_time_hcal_vs_BB = new TH2D("h2_RF_time_hcal_vs_BB","RF time HCAL vs BB; BB RF time (ns); HCAL RF time (ns)",501,-0.5,500.5,501,-0.5,500.5);

  double trigbits;
  
  // Declare branches
  TTree *T = (TTree*) gDirectory->Get("T");
  T->SetBranchStatus("*",0);
  T->SetBranchStatus("bb.tdctrig.tdc",1);
  T->SetBranchStatus("bb.sh.nclus",1);
  T->SetBranchStatus("sbs.hcal.nclus",1);
  T->SetBranchStatus("Ndata.bb.tdctrig.tdcelemID",1);
  T->SetBranchStatus("bb.tdctrig.tdcelemID",1);
  T->SetBranchStatus("sbs.hcal.Ref.tdc",1);
  T->SetBranchStatus("sbs.hcal.Ref.tdcelemID",1);
  T->SetBranchStatus("Ndata.sbs.hcal.Ref.tdcelemID",1);
  T->SetBranchStatus("g.trigbits",1);
  
  T->SetBranchAddress("bb.tdctrig.tdc", &bb_tdctrig_tdc);
  T->SetBranchAddress("bb.tdctrig.tdcelemID", &bb_tdctrig_tdcelemID);
  T->SetBranchAddress("Ndata.bb.tdctrig.tdcelemID", &Ndata_bb_tdctrig_tdcelemID);
  T->SetBranchAddress("bb.sh.nclus", &bb_sh_nclus);
  T->SetBranchAddress("sbs.hcal.nclus", &sbs_hcal_nclus);
  T->SetBranchAddress("sbs.hcal.Ref.tdc",hcal_reftdc);
  T->SetBranchAddress("sbs.hcal.Ref.tdcelemID",hcal_reftdc_elemID);
  T->SetBranchAddress("Ndata.sbs.hcal.Ref.tdcelemID",&Ndata_hcal_reftdc);
  T->SetBranchAddress("g.trigbits", &trigbits);
  
  // Loop through events
  Long64_t nevents = T->GetEntries();

  int nbadRF_BB = 0;
  int nbadRF_HCAL = 0;

  for(Long64_t nevent=0; nevent<nevents; nevent++){

    T->GetEntry(nevent);

    //if( bb_sh_nclus==0 || sbs_hcal_nclus==0 ) continue;

    Double_t bb_time=0.;
    Double_t coin_time=0., hcal_time=0, bbslow_time=0, rf_time=0;

    bool anyBBRF = false, anyHCALRF = false;
    
    for(Int_t ihit=0; ihit<Ndata_bb_tdctrig_tdcelemID; ihit++){
      if(bb_tdctrig_tdcelemID[ihit]==0) hcal_time=bb_tdctrig_tdc[ihit];
      if(bb_tdctrig_tdcelemID[ihit]==1) coin_time=bb_tdctrig_tdc[ihit];
      if(bb_tdctrig_tdcelemID[ihit]==2) bbslow_time=bb_tdctrig_tdc[ihit];
      if(bb_tdctrig_tdcelemID[ihit]==4) {
	rf_time=bb_tdctrig_tdc[ihit];
	anyBBRF = true;
      }
      if(bb_tdctrig_tdcelemID[ihit]==5) bb_time=bb_tdctrig_tdc[ihit];
    }

    Double_t hcal_rftime = 0.0, bbtime_hcal=0.0;
    for( int ihit=0; ihit<Ndata_hcal_reftdc; ihit++ ){
      if( hcal_reftdc_elemID[ihit]==2 ){ //BigBite trigger time:
	bbtime_hcal = hcal_reftdc[ihit];
      }
      if( hcal_reftdc_elemID[ihit]==3 ){
	hcal_rftime = hcal_reftdc[ihit];
	anyHCALRF = true;
      }
    }
    
    if( trigbits == 4 && !anyHCALRF ) nbadRF_HCAL++;
    if( trigbits == 4 && !anyBBRF ) nbadRF_BB++;
    
    //Take difference to eliminate 25ns jitter

    if( trigbits == 4 ){
      h1_BB->Fill(bbslow_time);
      h1_HCAL->Fill(hcal_time);
      h1_calDiff->Fill(hcal_time-bb_time);
      h1_RF_time->Fill(rf_time);
      
      h1_RF_time_hcal->Fill( hcal_rftime );
      
      h2_RF_time_hcal_vs_BB->Fill( rf_time, hcal_rftime );
    }
  }
  
  pads->cd(1);
  h1_BB->SetStats(0);
  h1_BB->Draw();

  pads->cd(2);
  h1_HCAL->SetStats(0);
  h1_HCAL->Draw();

  pads->cd(3);
  h1_calDiff->SetStats(0);
  h1_calDiff->Draw();

  TString stitle;
  
  pads->cd(4)->SetLogy();
  h1_RF_time->SetTitle(stitle.Format("RF time BB (N bad = %d); ns;events",nbadRF_BB));
  h1_RF_time->SetStats(0);
  h1_RF_time->Draw();

  pads->cd(5)->SetLogy();
  h1_RF_time_hcal->SetTitle(stitle.Format("RF time HCAL (N bad = %d);ns;events", nbadRF_HCAL));
  h1_RF_time_hcal->SetStats(0);
  h1_RF_time_hcal->Draw();

  pads->cd(6)->SetLogz();
  h2_RF_time_hcal_vs_BB->SetStats(0);
  h2_RF_time_hcal_vs_BB->Draw("colz");
  
  cout << "Processed macro with " << nevents << " entries." << endl;
  cout << "Note: a large spread indicates that HCal is setting the trigger time." << endl;

  st->Stop();
  cout << "CPU time = " << st->CpuTime() << " s " << " Real time = " << st->RealTime() << " s " << endl;
}
