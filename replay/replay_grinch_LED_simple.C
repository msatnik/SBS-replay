
//#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <iostream>

#include "TSystem.h"
#include "THaGlobals.h"
#include "TString.h"
#include "TFile.h"
#include "TList.h"
#include "TObject.h"
#include "TClonesArray.h"

#include "THaEvData.h"
#include "THaEvent.h"
#include "THaRun.h"
#include "THaAnalyzer.h"
#include "THaVarList.h"
#include "THaInterface.h"
#include "THaGoldenTrack.h"
#include "THaPrimaryKine.h"
#include "THaDecData.h"

#include "SBSBigBite.h"
#include "SBSBBShower.h"
#include "SBSBBTotalShower.h"
#include "SBSGRINCH.h"
//#include "SBSCherenkovDetector.h"
#include "SBSEArm.h"
#include "SBSHCal.h"
#include "SBSTimingHodoscope.h"
#include "SBSGEMSpectrometerTracker.h"
#include "SBSGEMTrackerBase.h"
#include "SBSRasteredBeam.h"
#include "LHRSScalerEvtHandler.h"
#include "SBSScalerEvtHandler.h"
//#endif

void replay_grinch_LED_simple(UInt_t runnum=4744, Long_t nevents=50000, Long_t firstevent=0, const char *fname_prefix="e1209016", UInt_t firstsegment=0, UInt_t maxsegments=1, Int_t pedestalmode=0, Int_t cmplots=0)
{

  THaAnalyzer* analyzer = new THaAnalyzer;
  
  SBSBigBite* bigbite = new SBSBigBite("bb", "BigBite spectrometer" );
  SBSGRINCH *grinch_tdc = new SBSGRINCH("grinch_tdc","GRINCH TDC data");

  grinch_tdc->SetModeTDC(SBSModeTDC::kTDC);
  //grinch_tdc->SetModeTDC(SBSModeTDC::kCommonStartTDC);
  grinch_tdc->SetModeADC(SBSModeADC::kNone);
  grinch_tdc->SetStoreEmptyElements(kFALSE);
  grinch_tdc->SetStoreRawHits(kFALSE);// set to kTRUE to get the "hits" branch to look at all the hits in the multi-hit TDC
  grinch_tdc->SetDisableRefTDC(true);
  bigbite->AddDetector(grinch_tdc);
  gHaApps->Add(bigbite);
    

  // add decoder
  THaApparatus* decL = new THaDecData("DL","Misc. Decoder Data");
  gHaApps->Add( decL );
  
  gHaPhysics->Add( new THaGoldenTrack( "BB.gold", "BigBite golden track", "bb" ));
  gHaPhysics->Add( new THaPrimaryKine( "e.kine", "electron kinematics", "bb", 0.0, 0.938272 ));
  
  //THaInterface::SetDecoder( SBSSimDecoder::Class() );
  THaEvent* event = new THaEvent;

  TString prefix = gSystem->Getenv("DATA_DIR");
  
  bool segmentexists = true;
  int segment=firstsegment; 

  int lastsegment=firstsegment;
  
  TDatime now = TDatime();
  
  
  //EPAF: copied the following from replay_BBGEM.C, as this script seems to be thought to handle splits properly.
  int stream = 0;
 
  TClonesArray *filelist = new TClonesArray("THaRun",10);

  vector<TString> pathlist;
  pathlist.push_back( prefix );

  if( prefix != "/adaqeb1/data1" )
    pathlist.push_back( "/adaqeb1/data1" );

  if( prefix != "/adaqeb2/data1" )
    pathlist.push_back( "/adaqeb2/data1" );

  if( prefix != "/adaqeb3/data1" )
    pathlist.push_back( "/adaqeb3/data1" );

  if( prefix != "/adaq1/data1/sbs" )
    pathlist.push_back( "/adaq1/data1/sbs" );

  if( prefix != "/cache/mss/halla/sbs/raw" )
    pathlist.push_back( "/cache/mss/halla/sbs/raw" );

  if( prefix != "/cache/mss/halla/sbs/GEnII/raw" )
    pathlist.push_back( "/cache/mss/halla/sbs/GEnII/raw" );

  for( int i=0; i<pathlist.size(); i++ ){
    cout << "search paths = " << pathlist[i] << endl;
  }

  TDatime RunDate = TDatime(); 

  int max1 = maxsegments;

  int segcounter=0;
  
  
  //This loop adds all file segments found to the list of THaRuns to process:
  while( segcounter < max1 && segment - firstsegment < maxsegments ){
    
    TString codafilename;
    //codafilename.Form( "%s/bbgem_%d.evio.%d", prefix.Data(), runnum, segment );
    codafilename.Form("%s_%d.evio.%d.%d", fname_prefix, runnum, stream, segment );

    TString ftest(fname_prefix);
    if( ftest == "bbgem" || ftest == "e1209019_trigtest" ){
      codafilename.Form("%s_%d.evio.%d", fname_prefix, runnum, segment );
    }

    segmentexists = false;

    cout << "codafilename = " << codafilename << endl;

    for( int ipath=0; ipath<pathlist.size(); ipath++ ){
      TString searchname;
      searchname.Form( "%s/%s", pathlist[ipath].Data(), codafilename.Data() );
      
      if( !gSystem->AccessPathName( searchname.Data() ) ){
	segmentexists = true;
	break;
      }
    }
   
    if( segmentexists ){
      new( (*filelist)[segcounter] ) THaRun( pathlist, codafilename.Data(), "GMN run" );
      cout << "Added segment " << segment << ", CODA file name = " << codafilename << endl;

    } 
    if( segmentexists ){
      segcounter++;
      lastsegment = segment;
    }
    segment++;
  }

  cout << "n segments to analyze = " << segcounter << endl;

  prefix = gSystem->Getenv("OUT_DIR");

  TString outfilename, outfilebase;

  if( nevents > 0 ){ 

    outfilebase.Form( "%s_replayed_grinch_%d_stream%d_seg%d_%d_firstevent%d_nevent%d", fname_prefix, runnum,
		      stream, firstsegment, lastsegment, firstevent, nevents );
    outfilename.Form( "%s/%s.root", prefix.Data(), outfilebase.Data());
  } else {
    outfilebase.Form( "%s_fullreplay_grinch_%d_stream%d_seg%d_%d", fname_prefix, runnum,
		      stream, firstsegment, lastsegment );
    outfilename.Form( "%s/%s.root", prefix.Data(), outfilebase.Data());
  }
 

  analyzer->SetVerbosity(2);
  analyzer->SetMarkInterval(100);

  analyzer->EnableBenchmarks();
  
  // Define the analysis parameters
  analyzer->SetEvent( event );
  analyzer->SetOutFile( outfilename.Data() );
  // File to record cuts accounting information
  
  prefix = gSystem->Getenv("LOG_DIR");
  // analyzer->SetSummaryFile(Form("%s/replay_gmn_%d_stream%d_seg%d_%d.log", prefix.Data(), runnum, 
  // 				stream, firstsegment, lastsegment));
  analyzer->SetSummaryFile(Form("%s/%s.log", prefix.Data(), outfilebase.Data()));
  prefix = gSystem->Getenv("SBS_REPLAY");
  
  prefix += "/replay/";
  //cout << "Hello World, if you get to this point you are reading the right script!" << endl;
  int myrun = (int) runnum;	
  TString odef_filename;
  TString ftest(fname_prefix);

  odef_filename = "replay_grinch_LED.odef";


  odef_filename.Prepend(prefix);
  //odef_filename.Prepend("/work/halla/sbs/ewertz/SBS-replay/replay/");

  analyzer->SetOdefFile( odef_filename );
  
  //added cut list in order to have 
  TString cdef_filename = "replay_grinch_LED.cdef";
  
  cdef_filename.Prepend( prefix );
  
  analyzer->SetCutFile( cdef_filename );
  //analyzer->SetCompressionLevel(0); // turn off compression

  filelist->Compress();

  for( int iseg=0; iseg<filelist->GetEntries(); iseg++ ){
    THaRun *run = ( (THaRun*) (*filelist)[iseg] );
    if( nevents > 0 ) run->SetLastEvent(firstevent+nevents-1); //not sure if this will work as we want it to for multiple file segments chained together

    run->SetFirstEvent( firstevent );
    
    run->SetDataRequired(THaRunBase::kDate|THaRunBase::kRunNumber);

    run->Init();
    
    if( run->GetSegment() >= firstsegment && run->GetSegment() - firstsegment < maxsegments ){
      //std::cout 
      analyzer->Process(run);     // start the actual analysis
    }
  }
  
}

