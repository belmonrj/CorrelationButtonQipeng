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

void arguments(int,int,int, const char*);

void plot_pAu()
{
  arguments(index_CNT_FVTN,index_CNT_FVTS,index_FVTN_FVTS,"CNT_FVTN_FVTS");
  arguments(index_CNT_BBCS,index_CNT_FVTS,index_FVTS_BBCS,"CNT_BBCS_FVTS");
}

void arguments(int indexA, int indexB, int indexC, const char* name)
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

  double c2_raw[nptbins];
  double ec2_raw[nptbins];
  double c2_sub[nptbins];
  double ec2_sub[nptbins];
  double v2_raw[nptbins];
  double ev2_raw[nptbins];
  double v2_sub[nptbins];
  double ev2_sub[nptbins];
  double ptvalues[nptbins] = {0.35,0.7,1.25,1.75,2.35,3.5};
  for ( int i = 0; i < nptbins; ++i )
    {
      // ---
      double a = v22raw[indexA+i];
      double b = v22raw[indexB+i];
      double c = v22raw[indexC];
      double ea = ev22raw[indexA+i];
      double eb = ev22raw[indexB+i];
      double ec = ev22raw[indexC];
      double c2 = a*b/c;
      double ec2 = c2*sqrt((ea/a)*(ea/a)+(eb/b)*(eb/b)+(ec/c)*(ec/c));
      c2_raw[i] = c2;
      ec2_raw[i] = ec2;
      // ---
      a = v22sub[indexA+i];
      b = v22sub[indexB+i];
      c = v22sub[indexC];
      ea = ev22sub[indexA+i];
      eb = ev22sub[indexB+i];
      ec = ev22sub[indexC];
      c2 = a*b/c;
      ec2 = c2*sqrt((ea/a)*(ea/a)+(eb/b)*(eb/b)+(ec/c)*(ec/c));
      c2_sub[i] = c2;
      ec2_sub[i] = ec2;
      // ---
      if ( c2_raw[i] >= 0 ) v2_raw[i] = sqrt(c2_raw[i]);
      else v2_raw[i] = sqrt(-c2_raw[i]);
      ev2_raw[i] = ec2_raw[i]/(2*v2_raw[i]);
      // ---
      if ( c2_sub[i] >= 0 ) v2_sub[i] = sqrt(c2_sub[i]);
      else v2_sub[i] = sqrt(-c2_sub[i]);
      ev2_sub[i] = ec2_sub[i]/(2*v2_sub[i]);
    }

  TCanvas* c1 = new TCanvas("c1","");

  TGraphErrors* tge_c2_raw = new TGraphErrors(nptbins,ptvalues,c2_raw,0,ec2_raw);
  tge_c2_raw->SetMarkerColor(kBlack);
  tge_c2_raw->SetMarkerStyle(kFullCircle);
  tge_c2_raw->Draw("ap");
  if ( verbose ) c1->Print(Form("syh_c2_%s_raw.png",name));

  TGraphErrors* tge_c2_sub = new TGraphErrors(nptbins,ptvalues,c2_sub,0,ec2_sub);
  tge_c2_sub->SetMarkerColor(kBlack);
  tge_c2_sub->SetMarkerStyle(kOpenCircle);
  tge_c2_sub->Draw("ap");
  if ( verbose ) c1->Print(Form("syh_c2_%s_sub.png",name));

  TGraphErrors* tge_v2_raw = new TGraphErrors(nptbins,ptvalues,v2_raw,0,ev2_raw);
  tge_v2_raw->SetMarkerColor(kBlack);
  tge_v2_raw->SetMarkerStyle(kFullCircle);
  tge_v2_raw->Draw("ap");
  if ( verbose ) c1->Print(Form("syh_v2_%s_raw.png",name));

  TGraphErrors* tge_v2_sub = new TGraphErrors(nptbins,ptvalues,v2_sub,0,ev2_sub);
  tge_v2_sub->SetMarkerColor(kBlack);
  tge_v2_sub->SetMarkerStyle(kOpenCircle);
  tge_v2_sub->Draw("ap");
  if ( verbose ) c1->Print(Form("syh_v2_%s_sub.png",name));

  double xmin = 0.0;
  double xmax = 4.0;
  double ymin = -0.05;
  double ymax = 0.35;
  TH2D* hdummy = new TH2D("hdummy","",1,xmin,xmax,1,ymin,ymax);
  hdummy->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hdummy->GetYaxis()->SetTitle("v_{2}");
  hdummy->Draw();
  tge_v2_raw->Draw("p");
  tge_v2_sub->Draw("p");
  TLegend* leg = new TLegend(0.18,0.72,0.38,0.92);
  leg->SetHeader(name);
  leg->AddEntry(tge_v2_raw,"raw v_{2}","p");
  leg->AddEntry(tge_v2_sub,"sub v_{2} (ATLAS method)","p");
  leg->Draw();
  TLine* line = new TLine(xmin,0,xmax,0);
  line->SetLineWidth(2);
  line->SetLineStyle(2);
  line->Draw();
  c1->Print(Form("syh_v2_%s.png",name));

  xmin = 0.0;
  xmax = 4.0;
  ymin = -0.01;
  ymax = 0.1;
  delete hdummy;
  hdummy = new TH2D("hdummy","",1,xmin,xmax,1,ymin,ymax);
  hdummy->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hdummy->GetYaxis()->SetTitle("c_{2}");
  hdummy->Draw();
  tge_c2_raw->Draw("p");
  tge_c2_sub->Draw("p");
  delete leg;
  leg = new TLegend(0.18,0.72,0.38,0.92);
  leg->SetHeader(name);
  leg->AddEntry(tge_c2_raw,"raw c_{2}","p");
  leg->AddEntry(tge_c2_sub,"sub c_{2} (ATLAS method)","p");
  leg->Draw();
  line->Draw();
  c1->Print(Form("syh_c2_%s.png",name));

  delete c1;
  delete leg;
  delete hdummy;

}
