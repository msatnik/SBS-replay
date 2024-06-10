void plotHCALnp( TH1D *heff, double xmin=-3, double xmax=1.1 ){
  //heff->GetXaxis()->SetRangeUser( xmin - 0.2*(xmax-xmin), 
  //				  xmax + 0.2*(xmax-xmin) );

  gStyle->SetOptFit(1111); 
  gStyle->SetOptStat(0);
  
  TFitResultPtr fitresult = heff->Fit( "pol0", "SQ", "", xmin, xmax );
  
  double eff = fitresult->Parameter(0);
  double efferr = (fitresult->GetErrors())[0];

  gPad->SetGrid();

  heff->SetLineColor(1);
  heff->SetMarkerColor(1);
  heff->SetMarkerStyle(20);
  heff->SetMarkerSize(0.75);

  heff->GetYaxis()->SetRangeUser(0,1.0);
  
  heff->Draw();

  TString newtitle; 

  newtitle.Form( "N/P ratio = %6.4g #pm %6.4g", eff, efferr );
  heff->SetTitle(newtitle.Data());

  heff->GetYaxis()->SetTitle( "N/P ratio" );


}
