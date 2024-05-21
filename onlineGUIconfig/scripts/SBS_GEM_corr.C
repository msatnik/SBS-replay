void SBS_GEM_corr(){

  gStyle->SetPalette(1);


  TTree *t = (TTree*) gFile->Get("T");

  TH2F *hit2d = new TH2F("hh","HCAL, SBS GEM correlation;SBS GEM #theta (rad);HCal x (m)",100,-0.2,0.1,100,-3,2);

  
  t->Draw("sbs.hcal.x:sbs.tr.th[0]>>hh","","colz" );

  gStyle->SetOptStat(0);  
 
}
