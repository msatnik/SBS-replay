void overlay(TString hist1name, TString hist2name, TString legentry1="", TString legentry2="",TString samestr="",Bool_t golden=false){
  TH1F* H1;
  TH1F* H2;
  
  H1 = (TH1F*) gFile->Get(hist1name);
  H2 = (TH1F*) gFile->Get(hist2name);

  H1->SetStats(false);
  H2->SetStats(false);
  H1->SetTitle("");
  H2->SetTitle("");

  H1->SetMarkerStyle(20);
  H2->SetMarkerStyle(21);
  H1->SetMarkerColor(1);
  H1->SetLineColor(1);
  H2->SetMarkerColor(2);
  H2->SetLineColor(2);

  TString DrawOpt="HIST";
  
  if(golden){
    H1->SetLineColor(30);
    H1->SetFillColor(30);
    H1->SetFillStyle(3027);
    H2->SetLineColor(46);
    H2->SetFillColor(46);
    H2->SetFillStyle(3027);
  } else {
    H1->SetLineColor(1);
    H2->SetLineColor(2);
    
    H1->SetLineWidth(2);
    H2->SetLineWidth(2);
    H2->SetLineStyle(10);
    
    //H1->SetFillColorAlpha(4,0.2);
    //H2->SetFillColorAlpha(2,0.2);
    //H1->SetMarkerStyle(20);
    //H1->SetMarkerColor(4);
    //H1->SetMarkerSize(1.0);
    //H2->SetMarkerStyle(25);
    //H2->SetMarkerColor(2);
    //H2->SetMarkerSize(1.0);
  }

  if(golden){
    if(!samestr.Contains("sames")) samestr += "sames";
    DrawOpt += samestr;
    
    H1->Draw(DrawOpt);
    H2->Draw(DrawOpt);
  } else {
    if(!golden && !samestr.Contains("sames")){
      // pick which hist to draw first:
      if(H1->GetMaximum() > H2->GetMaximum()){
	H1->Draw(DrawOpt);
	if(!samestr.Contains("sames")) samestr += "sames";
	H2->Draw(DrawOpt+samestr);
      } else {
	H2->Draw(DrawOpt);
	if(!samestr.Contains("sames")) samestr += "sames";
	H1->Draw(DrawOpt+samestr);
      }
    } else {
      H1->Draw(DrawOpt);
      if(!samestr.Contains("sames")) samestr += "sames";
      H2->Draw(DrawOpt+samestr);
    }
  }

  if( legentry1 == "" ) legentry1 = H1->GetName();
  if( legentry2 == "" ) legentry2 = H2->GetName();
  
  //Now add a small legend in the upper right corner:
  TLegend *leg_thisplot = new TLegend( 0.6,0.55,0.9,0.88, "", "brNDC" );

  leg_thisplot->SetFillStyle(0);
  leg_thisplot->SetBorderSize(0);
  leg_thisplot->AddEntry(H1, legentry1, "lf" );
  leg_thisplot->AddEntry(H2, legentry2, "lf" );

  leg_thisplot->Draw();
  
}
