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

TFile* outfile;

void arguments(int,int,int, const char*);

void newplot_v33_pAu()
{

  outfile = new TFile("hdata_v33_pAu_pAulm.root","recreate");

  arguments(index_CNT_BBCN,index_CNT_BBCS,index_BBCN_BBCS,"CNT_BBCN_BBCS");
  arguments(index_CNT_FVTN,index_CNT_FVTS,index_FVTN_FVTS,"CNT_FVTN_FVTS");
  arguments(index_CNT_BBCS,index_CNT_FVTS,index_FVTS_BBCS,"CNT_BBCS_FVTS");
  arguments(index_CNT_BBCS,index_CNT_FVTN,index_FVTN_BBCS,"CNT_BBCS_FVTN");
  arguments(index_CNT_BBCN,index_CNT_FVTN,index_BBCN_FVTN,"CNT_BBCN_FVTN");
  arguments(index_CNT_BBCN,index_CNT_FVTS,index_BBCN_FVTS,"CNT_BBCN_FVTS");

  outfile->Write();
  outfile->Close();

}

void arguments(int indexA, int indexB, int indexC, const char* name)
{

  ifstream fin("data_v33_pAu.txt");

  bool verbose = false;

  string system;
  string junk;
  double v33raw[30];
  double ev33raw[30];
  double v33subA[30];
  double ev33subA[30];
  double v33subAZ[30];
  double ev33subAZ[30];
  double v33subRZ[30];
  double ev33subRZ[30];
  double v33subR[30];
  double ev33subR[30];
  for ( int i = 0; i < 30; ++i )
    {
      fin>>system>>junk>>junk>>v33raw[i]>>junk>>ev33raw[i]
         >>junk>>junk>>junk>>v33subA[i]>>junk>>ev33subA[i]
         >>junk>>junk>>junk>>v33subAZ[i]>>junk>>ev33subAZ[i]
         >>junk>>junk>>junk>>v33subRZ[i]>>junk>>ev33subRZ[i]
         >>junk>>junk>>junk>>v33subR[i]>>junk>>ev33subR[i];
      if ( verbose ) cout << i << " " << system << " "
                          << v33raw[i] << " " << ev33raw[i] << " "
                          << v33subA[i] << " " << ev33subA[i] << " "
                          << v33subAZ[i] << " " << ev33subAZ[i] << " "
                          << v33subRZ[i] << " " << ev33subRZ[i] << " "
                          << v33subR[i] << " " << ev33subR[i] << endl;
    }
  fin.close();

  double c3_raw[nptbins];
  double ec3_raw[nptbins];
  double v3_raw[nptbins];
  double ev3_raw[nptbins];
  double v3_subA[nptbins];
  double ev3_subA[nptbins];
  double c3_subA[nptbins];
  double ec3_subA[nptbins];
  double v3_subAZ[nptbins];
  double ev3_subAZ[nptbins];
  double c3_subAZ[nptbins];
  double ec3_subAZ[nptbins];
  double v3_subRZ[nptbins];
  double ev3_subRZ[nptbins];
  double c3_subRZ[nptbins];
  double ec3_subRZ[nptbins];
  double v3_subR[nptbins];
  double ev3_subR[nptbins];
  double c3_subR[nptbins];
  double ec3_subR[nptbins];
  double ptvalues[nptbins] = {0.35,0.7,1.2,1.7,2.35,3.5};
  double a,b,c,ea,eb,ec,c3,ec3;
  // ---
  char foutname[50];
  sprintf(foutname,"DataTextFiles/data_v3_pAu_%s.txt",name);
  cout << "foutname is " << foutname << endl;
  ofstream fout(foutname);
  // sprintf(foutname,"data_%s.root",name);
  // cout << "foutname is " << foutname << endl;
  // TFile* outfile = new TFile(foutname,"recreate");
  TH1D* h_v3_raw = new TH1D(Form("h_v3_raw_%s",name),"",nptbins,-0.5,nptbins-0.5);
  TH1D* h_v3_subA = new TH1D(Form("h_v3_subA_%s",name),"",nptbins,-0.5,nptbins-0.5);
  TH1D* h_v3_subAZ = new TH1D(Form("h_v3_subAZ_%s",name),"",nptbins,-0.5,nptbins-0.5);
  TH1D* h_v3_subRZ = new TH1D(Form("h_v3_subRZ_%s",name),"",nptbins,-0.5,nptbins-0.5);
  TH1D* h_v3_subR = new TH1D(Form("h_v3_subR_%s",name),"",nptbins,-0.5,nptbins-0.5);
  // ---
  for ( int i = 0; i < nptbins; ++i )
    {
      fout << ptvalues[i] << " ";
      // --- raw
      a = v33raw[indexA+i];
      b = v33raw[indexB+i];
      c = v33raw[indexC];
      ea = ev33raw[indexA+i];
      eb = ev33raw[indexB+i];
      ec = ev33raw[indexC];
      c3 = a*b/c;
      ec3 = c3*sqrt((ea/a)*(ea/a)+(eb/b)*(eb/b)+(ec/c)*(ec/c));
      c3_raw[i] = c3;
      ec3_raw[i] = ec3;
      // --- calculate the v3
      if ( c3_raw[i] >= 0 ) v3_raw[i] = sqrt(c3_raw[i]);
      else v3_raw[i] = -sqrt(-c3_raw[i]);
      ev3_raw[i] = ec3_raw[i]/(2*v3_raw[i]);
      if ( ev3_raw[i] < 0 ) ev3_raw[i] *= -1;
      if ( ev3_raw[i] > 1 ) ev3_raw[i] = 0;
      fout << v3_raw[i] << " " << ev3_raw[i] << " ";
      h_v3_raw->SetBinContent(i+1,v3_raw[i]);
      h_v3_raw->SetBinError(i+1,ev3_raw[i]);
      // --- subA
      a = v33subA[indexA+i];
      b = v33subA[indexB+i];
      c = v33subA[indexC];
      ea = ev33subA[indexA+i];
      eb = ev33subA[indexB+i];
      ec = ev33subA[indexC];
      c3 = a*b/c;
      ec3 = c3*sqrt((ea/a)*(ea/a)+(eb/b)*(eb/b)+(ec/c)*(ec/c));
      c3_subA[i] = c3;
      ec3_subA[i] = ec3;
      if ( c3_subA[i] >= 0 ) v3_subA[i] = sqrt(c3_subA[i]);
      else v3_subA[i] = -sqrt(-c3_subA[i]);
      ev3_subA[i] = ec3_subA[i]/(2*v3_subA[i]);
      if ( ev3_subA[i] < 0 ) ev3_subA[i] *= -1;
      if ( ev3_subA[i] > 1 ) ev3_subA[i] = 0;
      fout << v3_subA[i] << " " << ev3_subA[i] << " ";
      h_v3_subA->SetBinContent(i+1,v3_subA[i]);
      h_v3_subA->SetBinError(i+1,ev3_subA[i]);
      // --- subAZ
      a = v33subAZ[indexA+i];
      b = v33subAZ[indexB+i];
      c = v33subAZ[indexC];
      ea = ev33subAZ[indexA+i];
      eb = ev33subAZ[indexB+i];
      ec = ev33subAZ[indexC];
      c3 = a*b/c;
      ec3 = c3*sqrt((ea/a)*(ea/a)+(eb/b)*(eb/b)+(ec/c)*(ec/c));
      c3_subAZ[i] = c3;
      ec3_subAZ[i] = ec3;
      if ( c3_subAZ[i] >= 0 ) v3_subAZ[i] = sqrt(c3_subAZ[i]);
      else v3_subAZ[i] = -sqrt(-c3_subAZ[i]);
      ev3_subAZ[i] = ec3_subAZ[i]/(2*v3_subAZ[i]);
      if ( ev3_subAZ[i] < 0 ) ev3_subAZ[i] *= -1;
      if ( ev3_subAZ[i] > 1 ) ev3_subAZ[i] = 0;
      fout << v3_subAZ[i] << " " << ev3_subAZ[i] << " ";
      h_v3_subAZ->SetBinContent(i+1,v3_subAZ[i]);
      h_v3_subAZ->SetBinError(i+1,ev3_subAZ[i]);
      // --- subRZ
      a = v33subRZ[indexA+i];
      b = v33subRZ[indexB+i];
      c = v33subRZ[indexC];
      ea = ev33subRZ[indexA+i];
      eb = ev33subRZ[indexB+i];
      ec = ev33subRZ[indexC];
      c3 = a*b/c;
      ec3 = c3*sqrt((ea/a)*(ea/a)+(eb/b)*(eb/b)+(ec/c)*(ec/c));
      c3_subRZ[i] = c3;
      ec3_subRZ[i] = ec3;
      if ( c3_subRZ[i] >= 0 ) v3_subRZ[i] = sqrt(c3_subRZ[i]);
      else v3_subRZ[i] = sqrt(-c3_subRZ[i]);
      ev3_subRZ[i] = ec3_subRZ[i]/(2*v3_subRZ[i]);
      if ( ev3_subRZ[i] < 0 ) ev3_subRZ[i] *= -1;
      if ( ev3_subRZ[i] > 1 ) ev3_subRZ[i] = 0;
      fout << v3_subRZ[i] << " " << ev3_subRZ[i] << " ";
      h_v3_subRZ->SetBinContent(i+1,v3_subRZ[i]);
      h_v3_subRZ->SetBinError(i+1,ev3_subRZ[i]);
      // --- subR
      a = v33subR[indexA+i];
      b = v33subR[indexB+i];
      c = v33subR[indexC];
      ea = ev33subR[indexA+i];
      eb = ev33subR[indexB+i];
      ec = ev33subR[indexC];
      c3 = a*b/c;
      ec3 = c3*sqrt((ea/a)*(ea/a)+(eb/b)*(eb/b)+(ec/c)*(ec/c));
      c3_subR[i] = c3;
      ec3_subR[i] = ec3;
      if ( c3_subR[i] >= 0 ) v3_subR[i] = sqrt(c3_subR[i]);
      else v3_subR[i] = -sqrt(-c3_subR[i]);
      ev3_subR[i] = ec3_subR[i]/(2*v3_subR[i]);
      if ( ev3_subR[i] < 0 ) ev3_subR[i] *= -1;
      if ( ev3_subR[i] > 1 ) ev3_subR[i] = 0;
      fout << v3_subR[i] << " " << ev3_subR[i] << endl;
      h_v3_subR->SetBinContent(i+1,v3_subR[i]);
      h_v3_subR->SetBinError(i+1,ev3_subR[i]);
    }
  fout.close();
  outfile->cd();
  h_v3_raw->Write();
  h_v3_subA->Write();
  h_v3_subAZ->Write();
  h_v3_subRZ->Write();
  h_v3_subR->Write();
  delete  h_v3_raw;
  delete  h_v3_subA;
  delete  h_v3_subAZ;
  delete  h_v3_subRZ;
  delete  h_v3_subR;

  TCanvas* c1 = new TCanvas("c1","");

  TGraphErrors* tge_c3_raw = new TGraphErrors(nptbins,ptvalues,c3_raw,0,ec3_raw);
  tge_c3_raw->SetMarkerColor(kBlack);
  tge_c3_raw->SetMarkerStyle(kFullCircle);
  tge_c3_raw->Draw("ap");
  if ( verbose ) c1->Print(Form("PlotFigs/syh_pAu_c3_%s_raw.png",name));

  TGraphErrors* tge_v3_raw = new TGraphErrors(nptbins,ptvalues,v3_raw,0,ev3_raw);
  tge_v3_raw->SetMarkerColor(kBlack);
  tge_v3_raw->SetMarkerStyle(kFullCircle);
  tge_v3_raw->Draw("ap");
  if ( verbose ) c1->Print(Form("PlotFigs/syh_pAu_v3_%s_raw.png",name));

  TGraphErrors* tge_c3_subA = new TGraphErrors(nptbins,ptvalues,c3_subA,0,ec3_subA);
  tge_c3_subA->SetMarkerColor(kBlack);
  tge_c3_subA->SetMarkerStyle(kOpenCircle);
  tge_c3_subA->Draw("ap");
  if ( verbose ) c1->Print(Form("PlotFigs/syh_pAu_c3_%s_subA.png",name));

  TGraphErrors* tge_v3_subA = new TGraphErrors(nptbins,ptvalues,v3_subA,0,ev3_subA);
  tge_v3_subA->SetMarkerColor(kBlack);
  tge_v3_subA->SetMarkerStyle(kOpenCircle);
  tge_v3_subA->Draw("ap");
  if ( verbose ) c1->Print(Form("PlotFigs/syh_pAu_v3_%s_subA.png",name));

  TGraphErrors* tge_c3_subAZ = new TGraphErrors(nptbins,ptvalues,c3_subAZ,0,ec3_subAZ);
  tge_c3_subAZ->SetMarkerColor(kBlack);
  tge_c3_subAZ->SetMarkerStyle(kOpenSquare);
  tge_c3_subAZ->Draw("ap");
  if ( verbose ) c1->Print(Form("PlotFigs/syh_pAu_c3_%s_subAZ.png",name));

  TGraphErrors* tge_v3_subAZ = new TGraphErrors(nptbins,ptvalues,v3_subAZ,0,ev3_subAZ);
  tge_v3_subAZ->SetMarkerColor(kBlack);
  tge_v3_subAZ->SetMarkerStyle(kOpenSquare);
  tge_v3_subAZ->Draw("ap");
  if ( verbose ) c1->Print(Form("PlotFigs/syh_pAu_v3_%s_subAZ.png",name));

  TGraphErrors* tge_c3_subRZ = new TGraphErrors(nptbins,ptvalues,c3_subRZ,0,ec3_subRZ);
  tge_c3_subRZ->SetMarkerColor(kBlack);
  tge_c3_subRZ->SetMarkerStyle(kOpenCross);
  tge_c3_subRZ->Draw("ap");
  if ( verbose ) c1->Print(Form("PlotFigs/syh_pAu_c3_%s_subRZ.png",name));

  TGraphErrors* tge_v3_subRZ = new TGraphErrors(nptbins,ptvalues,v3_subRZ,0,ev3_subRZ);
  tge_v3_subRZ->SetMarkerColor(kBlack);
  tge_v3_subRZ->SetMarkerStyle(kOpenCross);
  tge_v3_subRZ->Draw("ap");
  if ( verbose ) c1->Print(Form("PlotFigs/syh_pAu_v3_%s_subRZ.png",name));

  TGraphErrors* tge_c3_subR = new TGraphErrors(nptbins,ptvalues,c3_subR,0,ec3_subR);
  tge_c3_subR->SetMarkerColor(kBlack);
  tge_c3_subR->SetMarkerStyle(kOpenDiamond);
  tge_c3_subR->SetMarkerSize(2.5);
  tge_c3_subR->Draw("ap");
  if ( verbose ) c1->Print(Form("PlotFigs/syh_pAu_c3_%s_subR.png",name));

  TGraphErrors* tge_v3_subR = new TGraphErrors(nptbins,ptvalues,v3_subR,0,ev3_subR);
  tge_v3_subR->SetMarkerColor(kBlack);
  tge_v3_subR->SetMarkerStyle(kOpenDiamond);
  tge_v3_subR->SetMarkerSize(2.5);
  tge_v3_subR->Draw("ap");
  if ( verbose ) c1->Print(Form("PlotFigs/syh_pAu_v3_%s_subR.png",name));

  double seyoungdata_CNT_BBCS_FVTS_raw[nptbins] = {0.01999,0.04210,0.07368,0.09789,0.12,0.14210};
  double seyoungdata_CNT_FVTN_FVTS_raw[nptbins] = {0.03368,0.06421,0.11157,0.15263,0.20210,0.28947};
  double seyoungdata_CNT_FVTN_FVTS_sub[nptbins] = {0.03052,0.05157,0.05684,0.06,0.01999,0};

  // double ppg191_pt[13] = {0.49594,0.69400,0.89383,1.09417,1.29436,1.49426,1.69450,1.89463,2.09479,2.29479,2.49508,2.69540,2.89522};
  // double ppg191_v3[13] = {0.029465,0.041451,0.054005,0.064461,0.076305,0.088801,0.098890,0.105504,0.114177,0.124270,0.130828,0.138556,0.137665};
  // double ppg191_su[13] = {0.0021042,0.0029602,0.0038567,0.0046034,0.0054492,0.0063416,0.0070621,0.0075344,0.0083902,0.0091319,0.0096138,0.0101817,0.0101163};
  // double ppg191_sd[13] = {0.0039422,0.0060338,0.0078612,0.0117858,0.0139513,0.0162361,0.0205304,0.0219035,0.0257567,0.0280336,0.0295129,0.0328588,0.0326475};

  // TGraphAsymmErrors* tgae_ppg191 = new TGraphAsymmErrors(13,ppg191_pt,ppg191_v3,0,0,ppg191_sd,ppg191_su);
  // tgae_ppg191->SetLineWidth(2);
  // tgae_ppg191->SetLineColor(kRed);
  // tgae_ppg191->SetFillColorAlpha(kRed,0.35);

  // --- PPG216 data table
  double p216_x[10];
  double p216_y[10];
  double p216_ey[10];
  double p216_esyl[10];
  double p216_esyh[10];
  // pT                  v3                     v3StatErr               v3PosSystErr                 v3NegSysErr
  p216_x[0] = 0.5;    p216_y[0] = 0.0033;    p216_ey[0] = 0.0012;    p216_esyh[0] = 0.0007;       p216_esyl[0] = 0.0006;
  p216_x[1] = 0.7;    p216_y[1] = 0.0067;    p216_ey[1] = 0.0013;    p216_esyh[1] = 0.0018;       p216_esyl[1] = 0.0012;
  p216_x[2] = 0.9;    p216_y[2] = 0.0081;    p216_ey[2] = 0.0017;    p216_esyh[2] = 0.0028;       p216_esyl[2] = 0.0015;
  p216_x[3] = 1.1;    p216_y[3] = 0.0098;    p216_ey[3] = 0.0022;    p216_esyh[3] = 0.0042;       p216_esyl[3] = 0.0018;
  p216_x[4] = 1.3;    p216_y[4] = 0.0130;    p216_ey[4] = 0.0027;    p216_esyh[4] = 0.0066;       p216_esyl[4] = 0.0023;
  p216_x[5] = 1.5;    p216_y[5] = 0.0158;    p216_ey[5] = 0.0034;    p216_esyh[5] = 0.0094;       p216_esyl[5] = 0.0028;
  p216_x[6] = 1.7;    p216_y[6] = 0.0102;    p216_ey[6] = 0.0043;    p216_esyh[6] = 0.0069;       p216_esyl[6] = 0.0018;
  p216_x[7] = 1.9;    p216_y[7] = 0.0133;    p216_ey[7] = 0.0053;    p216_esyh[7] = 0.0102;       p216_esyl[7] = 0.0024;
  p216_x[8] = 2.25;   p216_y[8] = 0.0193;    p216_ey[8] = 0.0046;    p216_esyh[8] = 0.0178;       p216_esyl[8] = 0.0035;
  p216_x[9] = 2.75;   p216_y[9] = 0.0211;    p216_ey[9] = 0.0094;    p216_esyh[9] = 0.0241;       p216_esyl[9] = 0.0038;
  TGraphAsymmErrors* tgae_p216_sys = new TGraphAsymmErrors(10,p216_x,p216_y,0,0,p216_esyl,p216_esyh);
  tgae_p216_sys->SetLineWidth(20);
  tgae_p216_sys->SetLineColor(kGray);
  tgae_p216_sys->SetLineStyle(1);
  tgae_p216_sys->SetFillColorAlpha(kBlack,0.35);
  //tgae_p216_sys->Draw("pz");
  TGraphErrors* tge_p216_stat = new TGraphErrors(10,p216_x,p216_y,0,p216_ey);
  tge_p216_stat->SetMarkerColor(kBlack);
  tge_p216_stat->SetMarkerStyle(kFullCircle);
  //tge_p216_stat->Draw("p");


  // TGraph* tg_sy_CNT_BBCS_FVTS_raw = new TGraph(nptbins,ptvalues,seyoungdata_CNT_BBCS_FVTS_raw);
  // tg_sy_CNT_BBCS_FVTS_raw->SetLineColor(kBlue);
  // TGraph* tg_sy_CNT_FVTN_FVTS_raw = new TGraph(nptbins,ptvalues,seyoungdata_CNT_FVTN_FVTS_raw);
  // tg_sy_CNT_FVTN_FVTS_raw->SetLineColor(kBlue);
  // //TGraph* tg_sy_CNT_FVTN_FVTS_sub = new TGraph(nptbins,ptvalues,seyoungdata_CNT_FVTN_FVTS_sub);
  // TGraph* tg_sy_CNT_FVTN_FVTS_sub = new TGraph(5,ptvalues,seyoungdata_CNT_FVTN_FVTS_sub);
  // tg_sy_CNT_FVTN_FVTS_sub->SetLineColor(kBlue);

  double xmin = 0.0;
  double xmax = 4.0;
  double ymin = -0.05;
  double ymax = 0.1;
  //arguments(index_CNT_BBCS,index_CNT_FVTS,index_FVTS_BBCS,"CNT_BBCS_FVTS");
  bool is_CNT_BBCS_FVTS = ( indexA == index_CNT_BBCS && indexB == index_CNT_FVTS && indexC == index_FVTS_BBCS );
  bool is_CNT_FVTN_FVTS = ( indexA == index_CNT_FVTN && indexB == index_CNT_FVTS && indexC == index_FVTN_FVTS );
  if ( is_CNT_BBCS_FVTS ) ymin = -0.01;
  if ( is_CNT_BBCS_FVTS ) ymax = 0.05;
  TH2D* hdummy = new TH2D("hdummy","",1,xmin,xmax,1,ymin,ymax);
  hdummy->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hdummy->GetYaxis()->SetTitle("v_{3}");
  hdummy->Draw();
  //tgae_ppg191->Draw("L3");
  tgae_p216_sys->Draw("pz");
  tge_p216_stat->Draw("p");
  tge_v3_raw->Draw("p");
  // if ( is_CNT_BBCS_FVTS ) tg_sy_CNT_BBCS_FVTS_raw->Draw("l");
  // if ( is_CNT_FVTN_FVTS ) tg_sy_CNT_FVTN_FVTS_raw->Draw("l");
  TLegend* leg = new TLegend(0.18,0.65,0.38,0.92);
  leg->SetHeader(name);
  // if ( is_CNT_BBCS_FVTS ) leg->AddEntry(tg_sy_CNT_BBCS_FVTS_raw,"raw v_{3} (Seyoung)","l");
  // if ( is_CNT_FVTN_FVTS ) leg->AddEntry(tg_sy_CNT_FVTN_FVTS_raw,"raw v_{3} (Seyoung)","l");
  leg->AddEntry(tge_v3_raw,"raw v_{3}","p");
  leg->Draw();
  TLine* line = new TLine(xmin,0,xmax,0);
  line->SetLineWidth(2);
  line->SetLineStyle(2);
  line->Draw();
  c1->Print(Form("PlotFigs/syh_pAu_v3_%s.png",name));
  tge_v3_subA->Draw("p");
  //  if ( is_CNT_FVTN_FVTS ) tg_sy_CNT_FVTN_FVTS_sub->Draw("l");
  leg->AddEntry(tge_v3_subA,"sub v_{3} (ATLAS)","p");
  c1->Print(Form("PlotFigs/syh_pAu_v3_%s_sub1.png",name));
  tge_v3_subAZ->Draw("p");
  leg->AddEntry(tge_v3_subAZ,"sub v_{3} (ATLAS, ZYAM)","p");
  c1->Print(Form("PlotFigs/syh_pAu_v3_%s_sub2.png",name));
  tge_v3_subR->Draw("p");
  leg->AddEntry(tge_v3_subR,"sub v_{3} (Reference)","p");
  c1->Print(Form("PlotFigs/syh_pAu_v3_%s_sub3.png",name));

  xmin = 0.0;
  xmax = 4.0;
  ymin = -0.01;
  ymax = 0.1;
  delete hdummy;
  hdummy = new TH2D("hdummy","",1,xmin,xmax,1,ymin,ymax);
  hdummy->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hdummy->GetYaxis()->SetTitle("c_{3}");
  hdummy->Draw();
  tge_c3_raw->Draw("p");
  tge_c3_subA->Draw("p");
  tge_c3_subAZ->Draw("p");
  tge_c3_subR->Draw("p");
  delete leg;
  leg = new TLegend(0.18,0.72,0.38,0.92);
  leg->SetHeader(name);
  leg->AddEntry(tge_c3_raw,"raw c_{3}","p");
  leg->AddEntry(tge_c3_subA,"sub c_{3} (ATLAS)","p");
  leg->AddEntry(tge_c3_subAZ,"sub c_{3} (ATLAS,ZYAM)","p");
  leg->AddEntry(tge_c3_subR,"sub c_{3} (Reference)","p");
  leg->Draw();
  line->Draw();
  c1->Print(Form("PlotFigs/syh_pAu_c3_%s.png",name));

  delete c1;
  delete leg;
  delete hdummy;

}

//subtracted cnt fvtn fvts
// 0.3530, 0.03052
// 0.7081, 0.05157
// 1.2102, 0.05684
// 1.7061, 0.06
// 2.3428, 0.01999

//unsubtracted cnt fvtn fvts
// 0.3530, 0.03368
// 0.7081, 0.06421
// 1.2040, 0.11157
// 1.7061, 0.15263
// 2.3428, 0.20210
// 3.5367, 0.28947

// unsubtracted cnt bbcs fvts
// 0.3591, 0.01999
// 0.7020, 0.04210
// 1.2040, 0.07368
// 1.7061, 0.09789
// 2.3428, 0.12
// 3.5367, 0.14210

// pT      v3       stat. unc. sys. up   sys. down
// 0.49594 0.029465 0.00044149 0.0021042 0.0039422
// 0.69400 0.041451 0.00054303 0.0029602 0.0060338
// 0.89383 0.054005 0.00069428 0.0038567 0.0078612
// 1.09417 0.064461 0.00088503 0.0046034 0.0117858
// 1.29436 0.076305 0.00111408 0.0054492 0.0139513
// 1.49426 0.088801 0.00139669 0.0063416 0.0162361
// 1.69450 0.098890 0.00174562 0.0070621 0.0205304
// 1.89463 0.105504 0.00217696 0.0075344 0.0219035
// 2.09479 0.114177 0.00268917 0.0083902 0.0257567
// 2.29479 0.124270 0.00332228 0.0091319 0.0280336
// 2.49508 0.130828 0.00408305 0.0096138 0.0295129
// 2.69540 0.138556 0.00496819 0.0101817 0.0328588
// 2.89522 0.137665 0.00601138 0.0101163 0.0326475
