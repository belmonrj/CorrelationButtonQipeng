#include "../NonFlowSubtractor.C"


void testSimple(TH1D*, TH1D*, TH1D*);

void run_v33_pAu()
{
  const int nmultbins = 4;
  const int nptbins = 6;

  TH1D* h_CNT_BBCN[nmultbins][nptbins];
  TH1D* h_CNT_BBCS[nmultbins][nptbins];
  TH1D* h_CNT_FVTN[nmultbins][nptbins];
  TH1D* h_CNT_FVTS[nmultbins][nptbins];

  TH1D* h_BBCN_BBCS[nmultbins];
  TH1D* h_BBCN_FVTN[nmultbins];
  TH1D* h_BBCN_FVTS[nmultbins];

  TH1D* h_FVTN_BBCS[nmultbins];
  TH1D* h_FVTN_FVTS[nmultbins];
  TH1D* h_FVTS_BBCS[nmultbins];

  TFile* fSeyoung = TFile::Open("seyoung_data.root");

  for ( int imult = 0; imult < nmultbins; ++imult )
    {
      // ---
      h_BBCN_BBCS[imult] = (TH1D*)fSeyoung->Get(Form("h_BBCN_BBCS_C%d",imult));
      h_BBCN_FVTN[imult] = (TH1D*)fSeyoung->Get(Form("h_BBCN_FVTN_C%d",imult));
      h_BBCN_FVTS[imult] = (TH1D*)fSeyoung->Get(Form("h_BBCN_FVTS_C%d",imult));
      if ( h_BBCN_BBCS[imult] == NULL ) cout << "uh oh BBCN_BBCS" << endl;
      if ( h_BBCN_FVTN[imult] == NULL ) cout << "uh oh BBCN_FVTN" << endl;
      if ( h_BBCN_FVTS[imult] == NULL ) cout << "uh oh BBCN_FVTS" << endl;
      // ---
      h_FVTN_BBCS[imult] = (TH1D*)fSeyoung->Get(Form("h_FVTN_BBCS_C%d",imult));
      h_FVTN_FVTS[imult] = (TH1D*)fSeyoung->Get(Form("h_FVTN_FVTS_C%d",imult));
      h_FVTS_BBCS[imult] = (TH1D*)fSeyoung->Get(Form("h_FVTS_BBCS_C%d",imult));
      if ( h_FVTN_BBCS[imult] == NULL ) cout << "uh oh FVTN_BBCS" << endl;
      if ( h_FVTN_FVTS[imult] == NULL ) cout << "uh oh FVTN_FVTS" << endl;
      if ( h_FVTS_BBCS[imult] == NULL ) cout << "uh oh FVTS_BBCS"  << endl;
      // ---
      for ( int ipt = 0; ipt < nptbins; ++ipt )
        {
          // ---
          h_CNT_BBCN[imult][ipt] = (TH1D*)fSeyoung->Get(Form("h_CNT_BBCN_C%d_pT%d",imult,ipt));
          if ( h_CNT_BBCN[imult][ipt] == NULL ) cout << "uh oh CNT_BBCN" << endl;
          // ---
          h_CNT_BBCS[imult][ipt] = (TH1D*)fSeyoung->Get(Form("h_CNT_BBCS_C%d_pT%d",imult,ipt));
          if ( h_CNT_BBCS[imult][ipt] == NULL ) cout << "uh oh CNT_BBCS" << endl;
          // ---
          h_CNT_FVTN[imult][ipt] = (TH1D*)fSeyoung->Get(Form("h_CNT_FVTN_C%d_pT%d",imult,ipt));
          if ( h_CNT_FVTN[imult][ipt] == NULL ) cout << "uh oh CNT_FVTN" << endl;
          // ---
          h_CNT_FVTS[imult][ipt] = (TH1D*)fSeyoung->Get(Form("h_CNT_FVTS_C%d_pT%d",imult,ipt));
          if ( h_CNT_FVTS[imult][ipt] == NULL ) cout << "uh oh CNT_FVTS" << endl;
        }
    }

  cout << "If you don't see an \"uh oh\" then all histograms have been pulled, now processing" << endl;

  cout << "BBCN_BBCS ";
  testSimple(h_BBCN_BBCS[3],h_BBCN_BBCS[0],h_BBCN_BBCS[2]);
  cout << "BBCN_FVTN ";
  testSimple(h_BBCN_FVTN[3],h_BBCN_FVTN[0],h_BBCN_FVTN[2]);
  cout << "BBCN_FVTS ";
  testSimple(h_BBCN_FVTS[3],h_BBCN_FVTS[0],h_BBCN_FVTS[2]);

  cout << "FVTN_BBCS ";
  testSimple(h_FVTN_BBCS[3],h_FVTN_BBCS[0],h_FVTN_BBCS[2]);
  cout << "FVTN_FVTS ";
  testSimple(h_FVTN_FVTS[3],h_FVTN_FVTS[0],h_FVTN_FVTS[2]);
  cout << "FVTS_BBCS ";
  testSimple(h_FVTS_BBCS[3],h_FVTS_BBCS[0],h_FVTS_BBCS[2]);

  cout << "Well that was fun, let's try the pT dependence now" << endl;

  for ( int ipt = 0; ipt < nptbins; ++ipt )
    {
      cout << "CNT_BBCN ";
      testSimple(h_CNT_BBCN[3][ipt],h_CNT_BBCN[0][ipt],h_CNT_BBCN[2][ipt]);
    }
  for ( int ipt = 0; ipt < nptbins; ++ipt )
    {
      cout << "CNT_BBCS ";
      testSimple(h_CNT_BBCS[3][ipt],h_CNT_BBCS[0][ipt],h_CNT_BBCS[2][ipt]);
    }
  for ( int ipt = 0; ipt < nptbins; ++ipt )
    {
      cout << "CNT_FVTN ";
      testSimple(h_CNT_FVTN[3][ipt],h_CNT_FVTN[0][ipt],h_CNT_FVTN[2][ipt]);
    }
  for ( int ipt = 0; ipt < nptbins; ++ipt )
    {
      cout << "CNT_FVTS ";
      testSimple(h_CNT_FVTS[3][ipt],h_CNT_FVTS[0][ipt],h_CNT_FVTS[2][ipt]);
    }

}

// --- nearly a direct copy from Qipeng
void testSimple(TH1D* h_correlation_LM, TH1D* h_correlation_HM, TH1D* h_correlation_LM2) {

    // TFile* fin = new TFile(Form("./correlationHistTest.root"),"READ");

    // TH1F* h_correlation_LM  = (TH1F*)fin->Get( Form("PTYRaw/dphi_Nch_ref") );
    // TH1F* h_correlation_HM  = (TH1F*)fin->Get( Form("PTYRaw/dphi_Nch5") );

    // // Correlation for the second lowest LM bin 
    // // will be used for correction
    // TH1F* h_correlation_LM2 = (TH1F*)fin->Get( Form("PTYRaw/dphi_Nch2") );

    NonFlowSubtractor subTool;
    NonFlowSubtractor subTool2;
    // change the fitter configure before call init() here
    //subTool.setAtlasFixedC3();
    //subTool.setAtlasFixedC4();
    //subTool.setDebug(); // Qipeng suggests to not use debug
    subTool.init();
    subTool2.setZYAM();
    subTool2.init();

    //--------------------------------------------------
    // ATLAS template fit
    //--------------------------------------------------
    //subResult theResult = subTool.templateFit(h_correlation_LM, h_correlation_HM, h_correlation_LM2); // this is for ATLAS improved method
    subResult theResult = subTool.templateFit(h_correlation_LM, h_correlation_HM); // no ZYAM
    subResult theResultR = subTool.referenceFit(h_correlation_LM, h_correlation_HM); // reference fit
    subResult theResult2 = subTool2.templateFit(h_correlation_LM, h_correlation_HM); // ATLAS with ZYAM
    subResult theResultRZ = subTool2.referenceFit(h_correlation_LM, h_correlation_HM); // reference fit with ZYAM (no change)
    //subResult theResult = subTool.templateHistFit(h_correlation_LM, h_correlation_HM);

    //--------------------------------------------------
    // Test for using symmetrized dphi correlation
    // OK from first look
    // need further confirmation
    //--------------------------------------------------
    //TH1* hsym_LM = Symmetrize(h_correlation_LM);
    //hsym_LM->SetName("hs_LM");
    //TH1* hsym_HM = Symmetrize(h_correlation_HM);
    //hsym_HM->SetName("hs_HM");
    //subResult theResult = subTool.templateFit(hsym_LM, hsym_HM);

    //--------------------------------------------------
    // Access the fitted results and plots
    //--------------------------------------------------
    //theResult.getCoeffRawValue(3)
    //cout << "v33 = " << theResult.getCoeffRawValue(3) << " +/- " << theResult.getCoeffRawError(3) << endl;
    cout << "v33 = " << theResult.getCoeffRawValue(3) << " +/- " << theResult.getCoeffRawError(3) << " and "; // raw
    cout << "v33subA = " << theResult.getCoeffSubValue(3) << " +/- " << theResult.getCoeffSubError(3) << " and "; // ATLAS no ZYAM sub
    cout << "v33subAZ = " << theResult2.getCoeffSubValue(3) << " +/- " << theResult2.getCoeffSubError(3) << " and "; // ATLAS with ZYAM sub
    cout << "v33subRZ = " << theResultRZ.getCoeffSubValue(3) << " +/- " << theResultRZ.getCoeffSubError(3) << " and "; // ATLAS with ZYAM sub
    cout << "v33subR = " << theResultR.getCoeffSubValue(3) << " +/- " << theResultR.getCoeffSubError(3) << endl; // Reference fitting method
    //cout << "Improved v33 = " << theResult.getCoeffSubImpValue(3) << " +/- " << theResult.getCoeffSubImpError(3) << endl; // 3 histos instead of 2, ATLAS improved method

    //TCanvas* c1 = new TCanvas("c1","scaling",50,50, 600,700);
    //TCanvas* c1 = new TCanvas("c1","scaling",50,50, 600,600);
    //subTool.plotAtlasHistHM(c1);
    //subTool.plotAtlasSubHM(c1);
    //c1->cd();
    // add addtional lengends

    //TCanvas* c2 = new TCanvas("c2","scaling",50,50, 600,600);
    //subTool.plotAtlasLM(c2);
    //c2->cd();
    // add addtional lengends

}
