void gather_pAu_morebins()
{

  TFile* file_CNT_BBCN = TFile::Open("CorHisto_CNT_BBCN_morebins.root");
  TFile* file_CNT_BBCS = TFile::Open("CorHisto_CNT_BBCS_morebins.root");
  TFile* file_CNT_FVTN = TFile::Open("CorHisto_CNT_FVTN_morebins.root");
  TFile* file_CNT_FVTS = TFile::Open("CorHisto_CNT_FVTS_morebins.root");

  TFile* file_BBCN_BBCS = TFile::Open("CorHisto_BBCN_BBCS.root");
  TFile* file_BBCN_FVTN = TFile::Open("CorHisto_BBCN_FVTN.root");
  TFile* file_BBCN_FVTS = TFile::Open("CorHisto_BBCN_FVTS.root");

  TFile* file_FVTN_BBCS = TFile::Open("CorHisto_FVTN_BBCS.root");
  TFile* file_FVTN_FVTS = TFile::Open("CorHisto_FVTN_FVTS.root");
  TFile* file_FVTS_BBCS = TFile::Open("CorHisto_FVTS_BBCS.root");

  const int nmultbins = 4;
  const int nptbins = 14;

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

  TFile* fSeyoung = new TFile("seyoung_data_morebins.root","recreate");

  for ( int imult = 0; imult < nmultbins; ++imult )
    {
      // ---
      h_BBCN_BBCS[imult] = (TH1D*)file_BBCN_BBCS->Get(Form("SIG_BBCN_BBCS_C%d",imult));
      h_BBCN_FVTN[imult] = (TH1D*)file_BBCN_FVTN->Get(Form("SIG_COR_BBC_FVTN_C%d",imult));
      h_BBCN_FVTS[imult] = (TH1D*)file_BBCN_FVTS->Get(Form("SIG_COR_BBCN_FVTS_C%d",imult));
      if ( h_BBCN_BBCS[imult] == NULL ) cout << "uh oh BBCN_BBCS" << endl;
      if ( h_BBCN_FVTN[imult] == NULL ) cout << "uh oh BBCN_FVTN" << endl;
      if ( h_BBCN_FVTS[imult] == NULL ) cout << "uh oh BBCN_FVTS" << endl;
      // ---
      h_FVTN_BBCS[imult] = (TH1D*)file_FVTN_BBCS->Get(Form("SIG_COR_BBCS_FVTN_C%d",imult));
      h_FVTN_FVTS[imult] = (TH1D*)file_FVTN_FVTS->Get(Form("SIG_FVTN_FVTS_C%d",imult));
      h_FVTS_BBCS[imult] = (TH1D*)file_FVTS_BBCS->Get(Form("SIG_COR_BBC_FVTS_C%d",imult));
      if ( h_FVTN_BBCS[imult] == NULL ) cout << "uh oh FVTN_BBCS" << endl;
      if ( h_FVTN_FVTS[imult] == NULL ) cout << "uh oh FVTN_FVTS" << endl;
      if ( h_FVTS_BBCS[imult] == NULL ) cout << "uh oh FVTS_BBCS"  << endl;
      // ---
      h_BBCN_BBCS[imult]->SetName(Form("h_BBCN_BBCS_C%d",imult));
      h_BBCN_FVTN[imult]->SetName(Form("h_BBCN_FVTN_C%d",imult));
      h_BBCN_FVTS[imult]->SetName(Form("h_BBCN_FVTS_C%d",imult));
      h_BBCN_BBCS[imult]->Write();
      h_BBCN_FVTN[imult]->Write();
      h_BBCN_FVTS[imult]->Write();
      h_FVTN_BBCS[imult]->SetName(Form("h_FVTN_BBCS_C%d",imult));
      h_FVTN_FVTS[imult]->SetName(Form("h_FVTN_FVTS_C%d",imult));
      h_FVTS_BBCS[imult]->SetName(Form("h_FVTS_BBCS_C%d",imult));
      h_FVTN_BBCS[imult]->Write();
      h_FVTN_FVTS[imult]->Write();
      h_FVTS_BBCS[imult]->Write();
      // ---
      for ( int ipt = 0; ipt < nptbins; ++ipt )
        {
          // ---
          h_CNT_BBCN[imult][ipt] = (TH1D*)file_CNT_BBCN->Get(Form("SIG_COR_CNT_BBCN_C%d_pT%d",imult,ipt));
          if ( h_CNT_BBCN[imult][ipt] == NULL ) cout << "uh oh CNT_BBCN" << endl;
          // ---
          h_CNT_BBCS[imult][ipt] = (TH1D*)file_CNT_BBCS->Get(Form("SIG_COR_CNT_BBCS_C%d_pT%d",imult,ipt));
          if ( h_CNT_BBCS[imult][ipt] == NULL ) cout << "uh oh CNT_BBCS" << endl;
          // ---
          h_CNT_FVTN[imult][ipt] = (TH1D*)file_CNT_FVTN->Get(Form("SIG_COR_CNT_FVTN_C%d_pT%d",imult,ipt));
          if ( h_CNT_FVTN[imult][ipt] == NULL ) cout << "uh oh CNT_FVTN" << endl;
          // ---
          h_CNT_FVTS[imult][ipt] = (TH1D*)file_CNT_FVTS->Get(Form("SIG_COR_CNT_FVTS_C%d_pT%d",imult,ipt));
          if ( h_CNT_FVTS[imult][ipt] == NULL ) cout << "uh oh CNT_FVTS" << endl;
          // ---
          h_CNT_BBCN[imult][ipt]->SetName(Form("h_CNT_BBCN_C%d_pT%d",imult,ipt));
          h_CNT_BBCN[imult][ipt]->Write();
          h_CNT_BBCS[imult][ipt]->SetName(Form("h_CNT_BBCS_C%d_pT%d",imult,ipt));
          h_CNT_BBCS[imult][ipt]->Write();
          h_CNT_FVTN[imult][ipt]->SetName(Form("h_CNT_FVTN_C%d_pT%d",imult,ipt));
          h_CNT_FVTN[imult][ipt]->Write();
          h_CNT_FVTS[imult][ipt]->SetName(Form("h_CNT_FVTS_C%d_pT%d",imult,ipt));
          h_CNT_FVTS[imult][ipt]->Write();
        }
    }



}

