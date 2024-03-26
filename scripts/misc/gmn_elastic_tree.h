//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Feb  4 14:47:04 2024 by ROOT version 6.26/10
// from TTree Tout/Tree containing variables for momentum calibration
// found on file: /Users/puckett/SBSrootfiles3/pass2/SBS4/LH2/Elastic_SBS4_LH2_SBS0percent.root
//////////////////////////////////////////////////////////

#ifndef gmn_elastic_tree_h
#define gmn_elastic_tree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class gmn_elastic_tree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           runnum;
   Int_t           HCALcut;
   Int_t           BBcut;
   Double_t        Ebeam;
   Double_t        Q2;
   Double_t        etheta;
   Double_t        ephi;
   Double_t        ep_recon;
   Double_t        ep_elastic;
   Double_t        ep_incident;
   Double_t        thetabend;
   Double_t        dpel;
   Double_t        W2;
   Double_t        xfp;
   Double_t        yfp;
   Double_t        thfp;
   Double_t        phfp;
   Double_t        xfp0;
   Double_t        yfp0;
   Double_t        thfp0;
   Double_t        phfp0;
   Double_t        thtgt;
   Double_t        phtgt;
   Double_t        ytgt;
   Double_t        xtgt;
   Double_t        xfcp;
   Double_t        yfcp;
   Double_t        zfcp;
   Double_t        xbcp;
   Double_t        ybcp;
   Double_t        zbcp;
   Double_t        bbtrchi2;
   Double_t        bbtrchi2hits;
   Double_t        bbtrnhits;
   Double_t        bbtrngoodhits;
   Double_t        bbtrt0;
   Double_t        vx;
   Double_t        vy;
   Double_t        vz;
   Double_t        BBdist;
   Double_t        BBtheta;
   Double_t        HCALdist;
   Double_t        HCALtheta;
   Double_t        xHCAL;
   Double_t        yHCAL;
   Double_t        xHCAL_expect;
   Double_t        yHCAL_expect;
   Double_t        xHCAL_expect_4vect;
   Double_t        yHCAL_expect_4vect;
   Double_t        EHCAL;
   Double_t        THCAL_adc;
   Double_t        THCAL_tdc;
   Double_t        deltax;
   Double_t        deltay;
   Double_t        pp_expect;
   Double_t        ptheta_expect;
   Double_t        pphi_expect;
   Double_t        EPS;
   Double_t        ESH;
   Double_t        Etot;
   Double_t        xSH;
   Double_t        ySH;
   Double_t        xPS;
   Double_t        TSH;
   Double_t        TPS;
   Double_t        Q2_4vect;
   Double_t        W2_4vect;
   Double_t        epx;
   Double_t        epy;
   Double_t        epz;
   Double_t        qx;
   Double_t        qy;
   Double_t        qz;
   Double_t        qmag;
   Double_t        thetaq;
   Double_t        phiq;
   Double_t        thetarecon_n;
   Double_t        phirecon_n;
   Double_t        thetarecon_p;
   Double_t        phirecon_p;
   Double_t        thetapq_n;
   Double_t        thetapq_p;
   Double_t        deltat;
   Double_t        deltat_adc;
   Double_t        deltax_4vect;
   Double_t        deltay_4vect;
   Double_t        protondeflection;
   Int_t           ibest_HCAL;
   Int_t           nblkHCAL;
   Double_t        againblkHCAL[42];   //[nblkHCAL]
   Double_t        atimeblkHCAL[42];   //[nblkHCAL]
   Double_t        eblkHCAL[42];   //[nblkHCAL]
   Double_t        idblkHCAL[42];   //[nblkHCAL]
   Double_t        rowblkHCAL[42];   //[nblkHCAL]
   Double_t        colblkHCAL[42];   //[nblkHCAL]
   Double_t        xblkHCAL[42];   //[nblkHCAL]
   Double_t        yblkHCAL[42];   //[nblkHCAL]
   Int_t           grinch_clsize;
   Int_t           grinch_tridx;
   Double_t        grinch_tmean;
   Double_t        grinch_xmean;
   Double_t        grinch_ymean;
   Double_t        grinch_adc;

   // List of branches
   TBranch        *b_runnum;   //!
   TBranch        *b_HCALcut;   //!
   TBranch        *b_BBcut;   //!
   TBranch        *b_Ebeam;   //!
   TBranch        *b_Q2;   //!
   TBranch        *b_etheta;   //!
   TBranch        *b_ephi;   //!
   TBranch        *b_ep_recon;   //!
   TBranch        *b_ep_elastic;   //!
   TBranch        *b_ep_incident;   //!
   TBranch        *b_thetabend;   //!
   TBranch        *b_dpel;   //!
   TBranch        *b_W2;   //!
   TBranch        *b_xfp;   //!
   TBranch        *b_yfp;   //!
   TBranch        *b_thfp;   //!
   TBranch        *b_phfp;   //!
   TBranch        *b_xfp0;   //!
   TBranch        *b_yfp0;   //!
   TBranch        *b_thfp0;   //!
   TBranch        *b_phfp0;   //!
   TBranch        *b_thtgt;   //!
   TBranch        *b_phtgt;   //!
   TBranch        *b_ytgt;   //!
   TBranch        *b_xtgt;   //!
   TBranch        *b_xfcp;   //!
   TBranch        *b_yfcp;   //!
   TBranch        *b_zfcp;   //!
   TBranch        *b_xbcp;   //!
   TBranch        *b_ybcp;   //!
   TBranch        *b_zbcp;   //!
   TBranch        *b_bbtrchi2;   //!
   TBranch        *b_bbtrchi2hits;   //!
   TBranch        *b_bbtrnhits;   //!
   TBranch        *b_bbtrngoodhits;   //!
   TBranch        *b_bbtrt0;   //!
   TBranch        *b_vx;   //!
   TBranch        *b_vy;   //!
   TBranch        *b_vz;   //!
   TBranch        *b_BBdist;   //!
   TBranch        *b_BBtheta;   //!
   TBranch        *b_HCALdist;   //!
   TBranch        *b_HCALtheta;   //!
   TBranch        *b_xHCAL;   //!
   TBranch        *b_yHCAL;   //!
   TBranch        *b_xHCAL_expect;   //!
   TBranch        *b_yHCAL_expect;   //!
   TBranch        *b_xHCAL_expect_4vect;   //!
   TBranch        *b_yHCAL_expect_4vect;   //!
   TBranch        *b_EHCAL;   //!
   TBranch        *b_THCAL_adc;   //!
   TBranch        *b_THCAL_tdc;   //!
   TBranch        *b_deltax;   //!
   TBranch        *b_deltay;   //!
   TBranch        *b_pp_expect;   //!
   TBranch        *b_ptheta_expect;   //!
   TBranch        *b_pphi_expect;   //!
   TBranch        *b_EPS;   //!
   TBranch        *b_ESH;   //!
   TBranch        *b_Etot;   //!
   TBranch        *b_xSH;   //!
   TBranch        *b_ySH;   //!
   TBranch        *b_xPS;   //!
   TBranch        *b_TSH;   //!
   TBranch        *b_TPS;   //!
   TBranch        *b_Q2_4vect;   //!
   TBranch        *b_W2_4vect;   //!
   TBranch        *b_epx;   //!
   TBranch        *b_epy;   //!
   TBranch        *b_epz;   //!
   TBranch        *b_qx;   //!
   TBranch        *b_qy;   //!
   TBranch        *b_qz;   //!
   TBranch        *b_qmag;   //!
   TBranch        *b_thetaq;   //!
   TBranch        *b_phiq;   //!
   TBranch        *b_thetarecon_n;   //!
   TBranch        *b_phirecon_n;   //!
   TBranch        *b_thetarecon_p;   //!
   TBranch        *b_phirecon_p;   //!
   TBranch        *b_thetapq_n;   //!
   TBranch        *b_thetapp_n;   //!
   TBranch        *b_deltat;   //!
   TBranch        *b_deltat_ADC;   //!
   TBranch        *b_deltax_4vect;   //!
   TBranch        *b_deltay_4vect;   //!
   TBranch        *b_protondeflection;   //!
   TBranch        *b_ibest_HCAL;   //!
   TBranch        *b_nblkHCAL;   //!
   TBranch        *b_againblkHCAL;   //!
   TBranch        *b_atimeblkHCAL;   //!
   TBranch        *b_eblkHCAL;   //!
   TBranch        *b_idblkHCAL;   //!
   TBranch        *b_rowblkHCAL;   //!
   TBranch        *b_colblkHCAL;   //!
   TBranch        *b_xblkHCAL;   //!
   TBranch        *b_yblkHCAL;   //!
   TBranch        *b_grinch_clsize;   //!
   TBranch        *b_grinch_tridx;   //!
   TBranch        *b_grinch_tmean;   //!
   TBranch        *b_grinch_xmean;   //!
   TBranch        *b_grinch_ymean;   //!
   TBranch        *b_grinch_adc;   //!

   gmn_elastic_tree(TTree *tree=0);
   virtual ~gmn_elastic_tree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef gmn_elastic_tree_cxx
gmn_elastic_tree::gmn_elastic_tree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/Users/puckett/SBSrootfiles3/pass2/SBS4/LH2/Elastic_SBS4_LH2_SBS0percent.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/Users/puckett/SBSrootfiles3/pass2/SBS4/LH2/Elastic_SBS4_LH2_SBS0percent.root");
      }
      f->GetObject("Tout",tree);

   }
   Init(tree);
}

gmn_elastic_tree::~gmn_elastic_tree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t gmn_elastic_tree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t gmn_elastic_tree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void gmn_elastic_tree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("runnum", &runnum, &b_runnum);
   fChain->SetBranchAddress("HCALcut", &HCALcut, &b_HCALcut);
   fChain->SetBranchAddress("BBcut", &BBcut, &b_BBcut);
   fChain->SetBranchAddress("Ebeam", &Ebeam, &b_Ebeam);
   fChain->SetBranchAddress("Q2", &Q2, &b_Q2);
   fChain->SetBranchAddress("etheta", &etheta, &b_etheta);
   fChain->SetBranchAddress("ephi", &ephi, &b_ephi);
   fChain->SetBranchAddress("ep_recon", &ep_recon, &b_ep_recon);
   fChain->SetBranchAddress("ep_elastic", &ep_elastic, &b_ep_elastic);
   fChain->SetBranchAddress("ep_incident", &ep_incident, &b_ep_incident);
   fChain->SetBranchAddress("thetabend", &thetabend, &b_thetabend);
   fChain->SetBranchAddress("dpel", &dpel, &b_dpel);
   fChain->SetBranchAddress("W2", &W2, &b_W2);
   fChain->SetBranchAddress("xfp", &xfp, &b_xfp);
   fChain->SetBranchAddress("yfp", &yfp, &b_yfp);
   fChain->SetBranchAddress("thfp", &thfp, &b_thfp);
   fChain->SetBranchAddress("phfp", &phfp, &b_phfp);
   fChain->SetBranchAddress("xfp0", &xfp0, &b_xfp0);
   fChain->SetBranchAddress("yfp0", &yfp0, &b_yfp0);
   fChain->SetBranchAddress("thfp0", &thfp0, &b_thfp0);
   fChain->SetBranchAddress("phfp0", &phfp0, &b_phfp0);
   fChain->SetBranchAddress("thtgt", &thtgt, &b_thtgt);
   fChain->SetBranchAddress("phtgt", &phtgt, &b_phtgt);
   fChain->SetBranchAddress("ytgt", &ytgt, &b_ytgt);
   fChain->SetBranchAddress("xtgt", &xtgt, &b_xtgt);
   fChain->SetBranchAddress("xfcp", &xfcp, &b_xfcp);
   fChain->SetBranchAddress("yfcp", &yfcp, &b_yfcp);
   fChain->SetBranchAddress("zfcp", &zfcp, &b_zfcp);
   fChain->SetBranchAddress("xbcp", &xbcp, &b_xbcp);
   fChain->SetBranchAddress("ybcp", &ybcp, &b_ybcp);
   fChain->SetBranchAddress("zbcp", &zbcp, &b_zbcp);
   fChain->SetBranchAddress("bbtrchi2", &bbtrchi2, &b_bbtrchi2);
   fChain->SetBranchAddress("bbtrchi2hits", &bbtrchi2hits, &b_bbtrchi2hits);
   fChain->SetBranchAddress("bbtrnhits", &bbtrnhits, &b_bbtrnhits);
   fChain->SetBranchAddress("bbtrngoodhits", &bbtrngoodhits, &b_bbtrngoodhits);
   fChain->SetBranchAddress("bbtrt0", &bbtrt0, &b_bbtrt0);
   fChain->SetBranchAddress("vx", &vx, &b_vx);
   fChain->SetBranchAddress("vy", &vy, &b_vy);
   fChain->SetBranchAddress("vz", &vz, &b_vz);
   fChain->SetBranchAddress("BBdist", &BBdist, &b_BBdist);
   fChain->SetBranchAddress("BBtheta", &BBtheta, &b_BBtheta);
   fChain->SetBranchAddress("HCALdist", &HCALdist, &b_HCALdist);
   fChain->SetBranchAddress("HCALtheta", &HCALtheta, &b_HCALtheta);
   fChain->SetBranchAddress("xHCAL", &xHCAL, &b_xHCAL);
   fChain->SetBranchAddress("yHCAL", &yHCAL, &b_yHCAL);
   fChain->SetBranchAddress("xHCAL_expect", &xHCAL_expect, &b_xHCAL_expect);
   fChain->SetBranchAddress("yHCAL_expect", &yHCAL_expect, &b_yHCAL_expect);
   fChain->SetBranchAddress("xHCAL_expect_4vect", &xHCAL_expect_4vect, &b_xHCAL_expect_4vect);
   fChain->SetBranchAddress("yHCAL_expect_4vect", &yHCAL_expect_4vect, &b_yHCAL_expect_4vect);
   fChain->SetBranchAddress("EHCAL", &EHCAL, &b_EHCAL);
   fChain->SetBranchAddress("THCAL_adc", &THCAL_adc, &b_THCAL_adc);
   fChain->SetBranchAddress("THCAL_tdc", &THCAL_tdc, &b_THCAL_tdc);
   fChain->SetBranchAddress("deltax", &deltax, &b_deltax);
   fChain->SetBranchAddress("deltay", &deltay, &b_deltay);
   fChain->SetBranchAddress("pp_expect", &pp_expect, &b_pp_expect);
   fChain->SetBranchAddress("ptheta_expect", &ptheta_expect, &b_ptheta_expect);
   fChain->SetBranchAddress("pphi_expect", &pphi_expect, &b_pphi_expect);
   fChain->SetBranchAddress("EPS", &EPS, &b_EPS);
   fChain->SetBranchAddress("ESH", &ESH, &b_ESH);
   fChain->SetBranchAddress("Etot", &Etot, &b_Etot);
   fChain->SetBranchAddress("xSH", &xSH, &b_xSH);
   fChain->SetBranchAddress("ySH", &ySH, &b_ySH);
   fChain->SetBranchAddress("xPS", &xPS, &b_xPS);
   fChain->SetBranchAddress("TSH", &TSH, &b_TSH);
   fChain->SetBranchAddress("TPS", &TPS, &b_TPS);
   fChain->SetBranchAddress("Q2_4vect", &Q2_4vect, &b_Q2_4vect);
   fChain->SetBranchAddress("W2_4vect", &W2_4vect, &b_W2_4vect);
   fChain->SetBranchAddress("epx", &epx, &b_epx);
   fChain->SetBranchAddress("epy", &epy, &b_epy);
   fChain->SetBranchAddress("epz", &epz, &b_epz);
   fChain->SetBranchAddress("qx", &qx, &b_qx);
   fChain->SetBranchAddress("qy", &qy, &b_qy);
   fChain->SetBranchAddress("qz", &qz, &b_qz);
   fChain->SetBranchAddress("qmag", &qmag, &b_qmag);
   fChain->SetBranchAddress("thetaq", &thetaq, &b_thetaq);
   fChain->SetBranchAddress("phiq", &phiq, &b_phiq);
   fChain->SetBranchAddress("thetarecon_n", &thetarecon_n, &b_thetarecon_n);
   fChain->SetBranchAddress("phirecon_n", &phirecon_n, &b_phirecon_n);
   fChain->SetBranchAddress("thetarecon_p", &thetarecon_p, &b_thetarecon_p);
   fChain->SetBranchAddress("phirecon_p", &phirecon_p, &b_phirecon_p);
   fChain->SetBranchAddress("thetapq_n", &thetapq_n, &b_thetapq_n);
   fChain->SetBranchAddress("thetapq_p", &thetapq_p, &b_thetapp_n);
   fChain->SetBranchAddress("deltat", &deltat, &b_deltat);
   fChain->SetBranchAddress("deltat_adc", &deltat_adc, &b_deltat_ADC);
   fChain->SetBranchAddress("deltax_4vect", &deltax_4vect, &b_deltax_4vect);
   fChain->SetBranchAddress("deltay_4vect", &deltay_4vect, &b_deltay_4vect);
   fChain->SetBranchAddress("protondeflection", &protondeflection, &b_protondeflection);
   fChain->SetBranchAddress("ibest_HCAL", &ibest_HCAL, &b_ibest_HCAL);
   fChain->SetBranchAddress("nblkHCAL", &nblkHCAL, &b_nblkHCAL);
   fChain->SetBranchAddress("againblkHCAL", againblkHCAL, &b_againblkHCAL);
   fChain->SetBranchAddress("atimeblkHCAL", atimeblkHCAL, &b_atimeblkHCAL);
   fChain->SetBranchAddress("eblkHCAL", eblkHCAL, &b_eblkHCAL);
   fChain->SetBranchAddress("idblkHCAL", idblkHCAL, &b_idblkHCAL);
   fChain->SetBranchAddress("rowblkHCAL", rowblkHCAL, &b_rowblkHCAL);
   fChain->SetBranchAddress("colblkHCAL", colblkHCAL, &b_colblkHCAL);
   fChain->SetBranchAddress("xblkHCAL", xblkHCAL, &b_xblkHCAL);
   fChain->SetBranchAddress("yblkHCAL", yblkHCAL, &b_yblkHCAL);
   fChain->SetBranchAddress("grinch_clsize", &grinch_clsize, &b_grinch_clsize);
   fChain->SetBranchAddress("grinch_tridx", &grinch_tridx, &b_grinch_tridx);
   fChain->SetBranchAddress("grinch_tmean", &grinch_tmean, &b_grinch_tmean);
   fChain->SetBranchAddress("grinch_xmean", &grinch_xmean, &b_grinch_xmean);
   fChain->SetBranchAddress("grinch_ymean", &grinch_ymean, &b_grinch_ymean);
   fChain->SetBranchAddress("grinch_adc", &grinch_adc, &b_grinch_adc);
   Notify();
}

Bool_t gmn_elastic_tree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void gmn_elastic_tree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t gmn_elastic_tree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef gmn_elastic_tree_cxx
