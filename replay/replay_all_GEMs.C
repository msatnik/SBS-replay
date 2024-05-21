#include "TSystem.h"
#include "THaGlobals.h"
#include "TString.h"
#include "TFile.h"
#include "TList.h"
#include "TObject.h"
#include "TClonesArray.h"

#include "THaEvData.h"
#include "MultiFileRun.h"
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

void replay_all_GEMs(UInt_t runnum=10491, Long_t nevents=-1, Long_t firstevent=1, const char *fname_prefix="e1209016", UInt_t firstsegment=0, UInt_t maxsegments=1, Int_t maxstream=2, Int_t pedestalmode=0, Int_t cmplots=1, Int_t usesbsgems=1)
{

  //  gSystem->Load("libsbs.so");

  int stream = 0;

  SBSBigBite   *bb = new SBSBigBite("bb", "Generic apparatus");
  //SBSGEMStand *gems = new SBSGEMStand("gems", "Collection of GEMs in stand");
  SBSGEMSpectrometerTracker *bbgem = new SBSGEMSpectrometerTracker("gem", "BigBite Hall A GEM data");
  bool pm =  ( pedestalmode != 0 );
  //this will override the database setting:
  bbgem->SetPedestalMode( pm );
  bbgem->SetMakeCommonModePlots( cmplots );
  bb->AddDetector(bbgem);

  //Add trigger TDC info:

  SBSGenericDetector* tdctrig= new SBSGenericDetector("tdctrig","BigBite shower TDC trig");
  tdctrig->SetModeADC(SBSModeADC::kNone);
  tdctrig->SetModeTDC(SBSModeTDC::kTDC);
  tdctrig->SetStoreEmptyElements(kFALSE);
  bb->AddDetector( tdctrig );

  SBSBBTotalShower* ts= new SBSBBTotalShower("ts", "sh", "ps", "BigBite shower");
  ts->SetDataOutputLevel(0);
  bb->AddDetector( ts );
  ts->SetStoreEmptyElements(kFALSE);
  ts->GetShower()->SetStoreEmptyElements(kFALSE);
  ts->GetPreShower()->SetStoreEmptyElements(kFALSE);

  SBSGenericDetector* bbtrig= new SBSGenericDetector("bbtrig","BigBite shower ADC trig");
  bbtrig->SetModeADC(SBSModeADC::kADC);
  bbtrig->SetModeTDC(SBSModeTDC::kTDC);
  bbtrig->SetStoreEmptyElements(kFALSE);
  bb->AddDetector( bbtrig );

  //bool pm =  ( pedestalmode != 0 );
  //this will override the database setting:
  ( static_cast<SBSGEMTrackerBase *> (bbgem) )->SetPedestalMode( pedestalmode );
  

    SBSEArm *harm = new SBSEArm("sbs","Hadron Arm with HCal");
  SBSHCal* hcal =  new SBSHCal("hcal","HCAL");
  //hcal->SetStoreRawHits(kTRUE);
  hcal->SetStoreEmptyElements(kFALSE);
  harm->AddDetector(hcal);

  SBSGenericDetector* sbstrig= new SBSGenericDetector("trig","HCal trigs");
  sbstrig->SetModeADC(SBSModeADC::kWaveform);
  //sbstrig->SetStoreRawHits(kTRUE);
  sbstrig->SetStoreEmptyElements(kFALSE);
  harm->AddDetector( sbstrig );  

  SBSGEMSpectrometerTracker *sbsgem = new SBSGEMSpectrometerTracker("gem", "Super BigBite Hall A GEM data");

  sbsgem->SetPedestalMode( pm );
  sbsgem->SetMakeCommonModePlots( cmplots );  
  if (usesbsgems != 0 ) harm->AddDetector(sbsgem);

  gHaApps->Add(harm);
  
  // Set up the equipment to be analyzed.
  
  // add the two spectrometers with the "standard" configuration
  // (VDC planes, S1, and S2)
  // Collect information about a easily modified random set of channels
  // (see DB_DIR/*/db_D.dat)
  /*
    THaApparatus* DECDAT = new THaDecData("D","Misc. Decoder Data");
    gHaApps->Add( DECDAT );
  */
  

  // Set up the analyzer - we use the standard one,
  // but this could be an experiment-specific one as well.
  // The Analyzer controls the reading of the data, executes
  // tests/cuts, loops over Apparatus's and PhysicsModules,
  // and executes the output routines.
  THaAnalyzer* analyzer = new THaAnalyzer;
  
  gHaApps->Add(bb);

  // A simple event class to be output to the resulting tree.
  // Creating your own descendant of THaEvent is one way of
  // defining and controlling the output.
 TString prefix = gSystem->Getenv("DATA_DIR");

  vector<string> pathlist;
  if( !prefix.IsNull() )
    pathlist.push_back( prefix.Data() );

  if( prefix != "/adaqeb[1-3]/data1" )
    pathlist.push_back( "/adaqeb[1-3]/data1" );

  if( prefix != "/adaq1/data1/sbs" )
    pathlist.push_back( "/adaq1/data1/sbs" );

  // Do not use wildcard filenames with this directory
  if( prefix != "/cache/mss/halla/sbs/raw" )
    pathlist.push_back( "/cache/mss/halla/sbs/raw" );

  if( prefix != "/cache/mss/halla/sbs/raw" )
    pathlist.push_back( "/cache/mss/halla/sbs/GEnII/raw" );

  for( const auto& path: pathlist ) {
    cout << "search paths = " << path << endl;
  }

  
  TString ftest(fname_prefix);
  bool test_data = ( ftest == "bbgem" || ftest == "e1209019_trigtest" );

  // Build list of files we'd like to analyze. These are candidates. If some of them
  // do not exist, they will be quietly skipped by the MultiFileRun class. A list of
  // the actually found files will be printed during initialization.
  vector<string> filenames;
  Int_t actual_maxstream = test_data ? 0 : maxstream;
  for( Int_t istream = 0; istream <= actual_maxstream; ++istream ) {
    for( UInt_t iseg = 0; iseg < maxsegments; ++iseg ) {
      TString codafilename;
      if( test_data )
	codafilename.Form("%s_%u.evio.%u", fname_prefix, runnum, istream);
      else
	codafilename.Form("%s_%u.evio.%d.%u", fname_prefix, runnum, istream, iseg);
      cout << "codafilename = " << codafilename << endl;
      filenames.emplace_back(codafilename.Data());
    }
  }
  if( filenames.empty() ) {
    cout << "No files. Check your parameters." << endl;
    return;
  }

  auto* run = new Podd::MultiFileRun(pathlist, filenames);
  run->SetFirstSegment(firstsegment);  // Starting segment number
  run->SetMaxSegments(maxsegments);    // Number of segments
  run->SetMaxStreams(maxstream+1);  // Number of streams

  // Count physics events sequentially, i.e. always starting at 1, for each invocation of
  // analyzer->Process(). This count will be compared to the run's event range, so it
  // can be specified independently of absolute event numbers. This is relevant if
  // the analysis starts with a segment > 0. Actual (raw, CODA) physics event numbers
  // will then start with some largish number (the first event number in the continuation
  // segment), but the analyzer's count will still start at 1 for that analysis,
  // so you can give the run an event range of, say, 1-50000 to analyze 50k events,
  // regardless of what segments are being replayed. The event header written for each
  // entry in the output ROOT tree does contain the actual CODA event number, which is
  // what you usually want for analysis.

  analyzer->SetCountMode(THaAnalyzer::kCountPhysics);

  if( nevents > 0 ) run->SetLastEvent(firstevent+nevents-1);
  run->SetFirstEvent( firstevent );
  run->SetDataRequired(THaRunBase::kDate|THaRunBase::kRunNumber);
  Int_t st = run->Init();
  if( st != THaRunBase::READ_OK ) {
    cerr << "========= Error initializing run" << endl;
    return;
  }

  cout << "n segments to analyze = " << run->GetNFiles() << endl;

  prefix = gSystem->Getenv("OUT_DIR");

  TString outfilename;
  UInt_t lastsegment = run->GetLastSegment();
  if( nevents > 0 ){
    outfilename.Form( "%s/%s_replayed_%u_stream%d_%d_seg%u_%u_firstevent%ld_nevent%ld.root", prefix.Data(), fname_prefix, runnum,
		      0, maxstream, firstsegment, lastsegment, firstevent, nevents );
  } else {
    outfilename.Form( "%s/%s_fullreplay_%u_stream%d_%d_seg%u_%u.root", prefix.Data(), fname_prefix, runnum,
		      0, maxstream, firstsegment, lastsegment );
  }

  analyzer->SetVerbosity(2);
  analyzer->SetMarkInterval(100);

  analyzer->EnableBenchmarks();
  
  // Define the analysis parameters
    analyzer->SetOutFile( outfilename.Data() );
  // File to record cuts accounting information
  analyzer->SetSummaryFile("replay_BBGEM.log"); // optional

  prefix = gSystem->Getenv("SBS_REPLAY");
  prefix += "/replay/";

  TString odef_filename = "replay_all_GEMs.odef";
  
  odef_filename.Prepend( prefix );

  analyzer->SetOdefFile( odef_filename );
  
  //analyzer->SetCompressionLevel(0); // turn off compression
  
  st = analyzer->Init(run);
  if( st != 0 ) {
    cerr << "========= Error initializing analyzer" << endl;
    return;
  }
  analyzer->Process(run);     // start the actual analysis
}
