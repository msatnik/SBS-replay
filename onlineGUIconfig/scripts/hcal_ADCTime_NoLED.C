void hcal_ADCTime_NoLED(){

  cout << "Processing hcal_ADCTime_NoLED.." << endl;

  TStopwatch *st = new TStopwatch();
  st->Start(kTRUE);

  TPad *pads = (TPad*) gPad;
  pads->Divide(1,2);

  Double_t sbs_hcal_clus_blk_e[290] = {0.}, sbs_hcal_clus_blk_atime[290] = {0.}, sbs_hcal_clus_blk_id[290]={0.}, sbs_hcal_nclus = 0., sbs_hcal_ledbit = -2.; 
  UInt_t fEvtHdr_fTrigBits = -1;

  TH1D *h1_hcal_adctime_noled = new TH1D("h1_hcal_adctime","HCal ADC Clus Time; ns",200,0.,200.);
  TH2D *h2_hcal_adctime_vs_elemID = new TH2D("h2_adctime_vs_elemID","HCAL ADC Clus Time vs ElemID; ElemID; adc time [ns]",289,0.,289.,200,0.,200.);

  // Declare branches
  TTree *T = (TTree*) gDirectory->Get("T");
  T->SetBranchStatus("*",0);
  T->SetBranchStatus("sbs.hcal.clus_blk.atime",1);
  T->SetBranchStatus("sbs.hcal.clus_blk.e",1);
  T->SetBranchStatus("sbs.hcal.clus_blk.id",1);
  T->SetBranchStatus("sbs.hcal.nclus",1);
  //T->SetBranchStatus("sbs.hcal.ledbit",1);
  T->SetBranchStatus("fEvtHdr.fTrigBits",1);

  T->SetBranchAddress("sbs.hcal.clus_blk.atime", &sbs_hcal_clus_blk_atime);
  T->SetBranchAddress("sbs.hcal.clus_blk.e", &sbs_hcal_clus_blk_e);
  T->SetBranchAddress("sbs.hcal.clus_blk.id", &sbs_hcal_clus_blk_id);
  T->SetBranchAddress("sbs.hcal.nclus", &sbs_hcal_nclus);
  //T->SetBranchAddress("sbs.hcal.ledbit", &sbs_hcal_ledbit);
  T->SetBranchAddress("fEvtHdr.fTrigBits", &fEvtHdr_fTrigBits);

  // Acquire the number of entries
  Long64_t nevents = T->GetEntries();
  
  for(Long64_t nevent=0; nevent<nevents; nevent++){

    T->GetEntry(nevent);
    
    //sbs_hcal_clus_blk_e[0]<0.35
    if( sbs_hcal_clus_blk_atime[0]<1 || fEvtHdr_fTrigBits==32 || sbs_hcal_nclus==0 ) continue;
    
    h1_hcal_adctime_noled->Fill(sbs_hcal_clus_blk_atime[0]);

    h2_hcal_adctime_vs_elemID->Fill(sbs_hcal_clus_blk_id[0],sbs_hcal_clus_blk_atime[0]);

  }
  
  pads->cd(1);
  h1_hcal_adctime_noled->SetMinimum(-0.1);
  h1_hcal_adctime_noled->SetStats(0);
  h1_hcal_adctime_noled->Draw();

  //TH2D *hATIME_HCAL_VS_ID; hATIME_HCAL_VS_ID = (TH2D*)gDirectory->Get("hATIME_HCAL_VS_ID");
  //hATIME_HCAL_VS_ID->Draw("colz");
  //hATIME_HCAL_VS_ID->SetStats(0);

  pads->cd(2);
  //gStyle->SetPalette(53);
  h2_hcal_adctime_vs_elemID->SetStats(0);
  h2_hcal_adctime_vs_elemID->Draw("colz");

  cout << "Processed macro with " << nevents << " entries." << endl;
  
  st->Stop();
  cout << "CPU time = " << st->CpuTime() << " s " << " Real time = " << st->RealTime() << " s " << endl;
}
