void plot_pAu()
{

  ifstream fin("full_data.txt");

  bool verbose = false;

  string system;
  string junk;
  double v22raw[30];
  double ev22raw[30];
  double v22sub[30];
  double ev22sub[30];
  for ( int i = 0; i < 30; ++i )
    {
      fin>>system>>junk>>junk>>v22raw[i]>>junk>>ev22raw[i]>>junk>>junk>>junk>>v22sub[i]>>junk>>ev22sub[i];
      if ( verbose ) cout << i << " " << system << " " << v22raw[i] << " " << ev22raw[i] << " " << v22sub[i] << " " << ev22sub[i] << endl;
    }

  const int index_BBCN_BBCS = 0;
  const int index_BBCN_FVTN = 1;
  const int index_BBCN_FVTS = 2;
  const int index_FVTN_BBCS = 3;
  const int index_FVTN_FVTS = 4;
  const int index_FVTS_BBCS = 5;
  const int nptbins = 6;
  const int index_CNT_BBCN = 6; // 6--11
  const int index_CNT_BBCS = 12; // 12--17
  const int index_CNT_FVTN = 18; // 18--23
  const int index_CNT_FVTS = 24; // 24--29

  double c2_CNT_FVTN_FVTS_raw[nptbins];
  double ec2_CNT_FVTN_FVTS_raw[nptbins];
  double c2_CNT_FVTN_FVTS_sub[nptbins];
  double ec2_CNT_FVTN_FVTS_sub[nptbins];
  double ptvalues[nptbins] = {0.5,1.0,1.5,2.0,2.5,3.0};
  for ( int i = 0; i < nptbins; ++i )
    {
      // ---
      double a = v22raw[index_CNT_FVTN+i];
      double b = v22raw[index_CNT_FVTS+i];
      double c = v22raw[index_FVTN_FVTS];
      double ea = v22raw[index_CNT_FVTN+i];
      double eb = v22raw[index_CNT_FVTS+i];
      double ec = v22raw[index_FVTN_FVTS];
      double c2 = a*b/c;
      double ec2 = c2*sqrt((ea/a)*(ea/a)+(eb/b)*(eb/b)+(ec/c)*(ec/c));
      c2_CNT_FVTN_FVTS_raw[i] = c2;
      ec2_CNT_FVTN_FVTS_raw[i] = ec2;
      // ---
      a = v22sub[index_CNT_FVTN+i];
      b = v22sub[index_CNT_FVTS+i];
      c = v22sub[index_FVTN_FVTS];
      ea = v22sub[index_CNT_FVTN+i];
      eb = v22sub[index_CNT_FVTS+i];
      ec = v22sub[index_FVTN_FVTS];
      c2 = a*b/c;
      ec2 = c2*sqrt((ea/a)*(ea/a)+(eb/b)*(eb/b)+(ec/c)*(ec/c));
      c2_CNT_FVTN_FVTS_sub[i] = c2;
      ec2_CNT_FVTN_FVTS_sub[i] = ec2;
    }

  TCanvas* c1 = new TCanvas("c1","");

  TGraphErrors* tge_c2_CNT_FVTN_FVTS_raw = new TGraphErrors(nptbins,ptvalues,c2_CNT_FVTN_FVTS_raw,0,ec2_CNT_FVTN_FVTS_raw);
  tge_c2_CNT_FVTN_FVTS_raw->SetMarkerColor(kBlack);
  tge_c2_CNT_FVTN_FVTS_raw->SetMarkerStyle(kOpenCircle);
  tge_c2_CNT_FVTN_FVTS_raw->Draw("ap");

  c1->Print("test_fig.png");

}
