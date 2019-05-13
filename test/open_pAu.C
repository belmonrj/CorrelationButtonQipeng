void open_pAu()
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

  cout << "If you don't see an \"uh oh\" then everything is copasetic" << endl;

}

