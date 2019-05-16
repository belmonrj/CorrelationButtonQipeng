const int index_BBCN_BBCS = 0;
const int index_BBCN_FVTN = 1;
const int index_BBCN_FVTS = 2;
const int index_FVTN_BBCS = 3;
const int index_FVTN_FVTS = 4;
const int index_FVTS_BBCS = 5;
const int nptbins = 5;
const int index_CNT_BBCN = 6; // 6--10
const int index_CNT_BBCS = 11; // 11--15
const int index_CNT_FVTN = 16; // 16--20
const int index_CNT_FVTS = 21; // 21--25

TFile* outfile;

void arguments(int,int,int, const char*);

void newplot_pAu_pp()
{

  //outfile = new TFile("hdata_pAu_pplm.root","recreate");

  arguments(index_CNT_BBCN,index_CNT_BBCS,index_BBCN_BBCS,"CNT_BBCN_BBCS");
  arguments(index_CNT_FVTN,index_CNT_FVTS,index_FVTN_FVTS,"CNT_FVTN_FVTS");
  arguments(index_CNT_BBCS,index_CNT_FVTS,index_FVTS_BBCS,"CNT_BBCS_FVTS");
  arguments(index_CNT_BBCS,index_CNT_FVTN,index_FVTN_BBCS,"CNT_BBCS_FVTN");
  arguments(index_CNT_BBCN,index_CNT_FVTN,index_BBCN_FVTN,"CNT_BBCN_FVTN");
  arguments(index_CNT_BBCN,index_CNT_FVTS,index_BBCN_FVTS,"CNT_BBCN_FVTS");

  // outfile->Write();
  // outfile->Close();

}

void arguments(int indexA, int indexB, int indexC, const char* name)
{

  ifstream fin("new_data_pp.txt");

  bool verbose = false;

  string system;
  string junk;
  double v22raw[26];
  double ev22raw[26];
  double v22subA[26];
  double ev22subA[26];
  double v22subAZ[26];
  double ev22subAZ[26];
  double v22subRZ[26];
  double ev22subRZ[26];
  double v22subR[26];
  double ev22subR[26];
  for ( int i = 0; i < 26; ++i )
    {
      fin>>system>>junk>>junk>>v22raw[i]>>junk>>ev22raw[i]
         >>junk>>junk>>junk>>v22subA[i]>>junk>>ev22subA[i]
         >>junk>>junk>>junk>>v22subAZ[i]>>junk>>ev22subAZ[i]
         >>junk>>junk>>junk>>v22subRZ[i]>>junk>>ev22subRZ[i]
         >>junk>>junk>>junk>>v22subR[i]>>junk>>ev22subR[i];
      if ( verbose ) cout << i << " " << system << " "
                          << v22raw[i] << " " << ev22raw[i] << " "
                          << v22subA[i] << " " << ev22subA[i] << " "
                          << v22subAZ[i] << " " << ev22subAZ[i] << " "
                          << v22subRZ[i] << " " << ev22subRZ[i] << " "
                          << v22subR[i] << " " << ev22subR[i] << endl;
    }
  fin.close();

  double c2_raw[nptbins];
  double ec2_raw[nptbins];
  double v2_raw[nptbins];
  double ev2_raw[nptbins];
  double v2_subA[nptbins];
  double ev2_subA[nptbins];
  double c2_subA[nptbins];
  double ec2_subA[nptbins];
  double v2_subAZ[nptbins];
  double ev2_subAZ[nptbins];
  double c2_subAZ[nptbins];
  double ec2_subAZ[nptbins];
  double v2_subRZ[nptbins];
  double ev2_subRZ[nptbins];
  double c2_subRZ[nptbins];
  double ec2_subRZ[nptbins];
  double v2_subR[nptbins];
  double ev2_subR[nptbins];
  double c2_subR[nptbins];
  double ec2_subR[nptbins];
  double ptvalues[nptbins] = {0.35,0.7,1.2,1.7,2.35};
  double a,b,c,ea,eb,ec,c2,ec2;
  // ---
  char* foutname;
  sprintf(foutname,"data_ppref_%s.txt",name);
  cout << "foutname is " << foutname << endl;
  ofstream fout(foutname);
  // sprintf(foutname,"data_%s.root",name);
  // cout << "foutname is " << foutname << endl;
  //TFile* outfile = new TFile(foutname,"recreate");
  // TH1D* h_v2_raw = new TH1D(Form("h_v2_raw_%s",name),"",nptbins,-0.5,nptbins-0.5);
  // TH1D* h_v2_subA = new TH1D(Form("h_v2_subA_%s",name),"",nptbins,-0.5,nptbins-0.5);
  // TH1D* h_v2_subAZ = new TH1D(Form("h_v2_subAZ_%s",name),"",nptbins,-0.5,nptbins-0.5);
  // TH1D* h_v2_subRZ = new TH1D(Form("h_v2_subRZ_%s",name),"",nptbins,-0.5,nptbins-0.5);
  // TH1D* h_v2_subR = new TH1D(Form("h_v2_subR_%s",name),"",nptbins,-0.5,nptbins-0.5);
  // ---
  for ( int i = 0; i < nptbins; ++i )
    {
      fout << ptvalues[i] << " ";
      // --- raw
      a = v22raw[indexA+i];
      b = v22raw[indexB+i];
      c = v22raw[indexC];
      ea = ev22raw[indexA+i];
      eb = ev22raw[indexB+i];
      ec = ev22raw[indexC];
      c2 = a*b/c;
      ec2 = c2*sqrt((ea/a)*(ea/a)+(eb/b)*(eb/b)+(ec/c)*(ec/c));
      c2_raw[i] = c2;
      ec2_raw[i] = ec2;
      // --- calculate the v2
      if ( c2_raw[i] >= 0 ) v2_raw[i] = sqrt(c2_raw[i]);
      else v2_raw[i] = sqrt(-c2_raw[i]);
      ev2_raw[i] = ec2_raw[i]/(2*v2_raw[i]);
      if ( ev2_raw[i] < 0 ) ev2_raw[i] *= -1;
      if ( ev2_raw[i] > 1 ) ev2_raw[i] = 0;
      fout << v2_raw[i] << " " << ev2_raw[i] << " ";
      // h_v2_raw->SetBinContent(i+1,v2_raw[i]);
      // h_v2_raw->SetBinError(i+1,ev2_raw[i]);
      // --- subA
      a = v22subA[indexA+i];
      b = v22subA[indexB+i];
      c = v22subA[indexC];
      ea = ev22subA[indexA+i];
      eb = ev22subA[indexB+i];
      ec = ev22subA[indexC];
      c2 = a*b/c;
      ec2 = c2*sqrt((ea/a)*(ea/a)+(eb/b)*(eb/b)+(ec/c)*(ec/c));
      c2_subA[i] = c2;
      ec2_subA[i] = ec2;
      if ( c2_subA[i] >= 0 ) v2_subA[i] = sqrt(c2_subA[i]);
      else v2_subA[i] = sqrt(-c2_subA[i]);
      ev2_subA[i] = ec2_subA[i]/(2*v2_subA[i]);
      if ( ev2_subA[i] < 0 ) ev2_subA[i] *= -1;
      if ( ev2_subA[i] > 1 ) ev2_subA[i] = 0;
      fout << v2_subA[i] << " " << ev2_subA[i] << " ";
      // h_v2_subA->SetBinContent(i+1,v2_subA[i]);
      // h_v2_subA->SetBinError(i+1,ev2_subA[i]);
      // --- subAZ
      a = v22subAZ[indexA+i];
      b = v22subAZ[indexB+i];
      c = v22subAZ[indexC];
      ea = ev22subAZ[indexA+i];
      eb = ev22subAZ[indexB+i];
      ec = ev22subAZ[indexC];
      c2 = a*b/c;
      ec2 = c2*sqrt((ea/a)*(ea/a)+(eb/b)*(eb/b)+(ec/c)*(ec/c));
      c2_subAZ[i] = c2;
      ec2_subAZ[i] = ec2;
      if ( c2_subAZ[i] >= 0 ) v2_subAZ[i] = sqrt(c2_subAZ[i]);
      else v2_subAZ[i] = sqrt(-c2_subAZ[i]);
      ev2_subAZ[i] = ec2_subAZ[i]/(2*v2_subAZ[i]);
      if ( ev2_subAZ[i] < 0 ) ev2_subAZ[i] *= -1;
      if ( ev2_subAZ[i] > 1 ) ev2_subAZ[i] = 0;
      fout << v2_subAZ[i] << " " << ev2_subAZ[i] << " ";
      // h_v2_subAZ->SetBinContent(i+1,v2_subAZ[i]);
      // h_v2_subAZ->SetBinError(i+1,ev2_subAZ[i]);
      // --- subRZ
      a = v22subRZ[indexA+i];
      b = v22subRZ[indexB+i];
      c = v22subRZ[indexC];
      ea = ev22subRZ[indexA+i];
      eb = ev22subRZ[indexB+i];
      ec = ev22subRZ[indexC];
      c2 = a*b/c;
      ec2 = c2*sqrt((ea/a)*(ea/a)+(eb/b)*(eb/b)+(ec/c)*(ec/c));
      c2_subRZ[i] = c2;
      ec2_subRZ[i] = ec2;
      if ( c2_subRZ[i] >= 0 ) v2_subRZ[i] = sqrt(c2_subRZ[i]);
      else v2_subRZ[i] = sqrt(-c2_subRZ[i]);
      ev2_subRZ[i] = ec2_subRZ[i]/(2*v2_subRZ[i]);
      if ( ev2_subRZ[i] < 0 ) ev2_subRZ[i] *= -1;
      if ( ev2_subRZ[i] > 1 ) ev2_subRZ[i] = 0;
      fout << v2_subRZ[i] << " " << ev2_subRZ[i] << " ";
      // h_v2_subRZ->SetBinContent(i+1,v2_subRZ[i]);
      // h_v2_subRZ->SetBinError(i+1,ev2_subRZ[i]);
      // --- subR
      a = v22subR[indexA+i];
      b = v22subR[indexB+i];
      c = v22subR[indexC];
      ea = ev22subR[indexA+i];
      eb = ev22subR[indexB+i];
      ec = ev22subR[indexC];
      c2 = a*b/c;
      ec2 = c2*sqrt((ea/a)*(ea/a)+(eb/b)*(eb/b)+(ec/c)*(ec/c));
      c2_subR[i] = c2;
      ec2_subR[i] = ec2;
      if ( c2_subR[i] >= 0 ) v2_subR[i] = sqrt(c2_subR[i]);
      else v2_subR[i] = sqrt(-c2_subR[i]);
      ev2_subR[i] = ec2_subR[i]/(2*v2_subR[i]);
      if ( ev2_subR[i] < 0 ) ev2_subR[i] *= -1;
      if ( ev2_subR[i] > 1 ) ev2_subR[i] = 0;
      fout << v2_subR[i] << " " << ev2_subR[i] << endl;
      // h_v2_subR->SetBinContent(i+1,v2_subR[i]);
      // h_v2_subR->SetBinError(i+1,ev2_subR[i]);
    }
  fout.close();
  //outfile->cd();
  // h_v2_raw->Write();
  // h_v2_subA->Write();
  // h_v2_subAZ->Write();
  // h_v2_subRZ->Write();
  // h_v2_subR->Write();
  // delete  h_v2_raw;
  // delete  h_v2_subA;
  // delete  h_v2_subAZ;
  // delete  h_v2_subRZ;
  // delete  h_v2_subR;



  TCanvas* c1 = new TCanvas("c1","");

  TGraphErrors* tge_c2_raw = new TGraphErrors(nptbins,ptvalues,c2_raw,0,ec2_raw);
  tge_c2_raw->SetMarkerColor(kBlack);
  tge_c2_raw->SetMarkerStyle(kFullCircle);
  tge_c2_raw->Draw("ap");
  if ( verbose ) c1->Print(Form("PlotFigs/syh_ppref_c2_%s_raw.png",name));

  TGraphErrors* tge_v2_raw = new TGraphErrors(nptbins,ptvalues,v2_raw,0,ev2_raw);
  tge_v2_raw->SetMarkerColor(kBlack);
  tge_v2_raw->SetMarkerStyle(kFullCircle);
  tge_v2_raw->Draw("ap");
  if ( verbose ) c1->Print(Form("PlotFigs/syh_ppref_v2_%s_raw.png",name));

  TGraphErrors* tge_c2_subA = new TGraphErrors(nptbins,ptvalues,c2_subA,0,ec2_subA);
  tge_c2_subA->SetMarkerColor(kBlack);
  tge_c2_subA->SetMarkerStyle(kOpenCircle);
  tge_c2_subA->Draw("ap");
  if ( verbose ) c1->Print(Form("PlotFigs/syh_ppref_c2_%s_subA.png",name));

  TGraphErrors* tge_v2_subA = new TGraphErrors(nptbins,ptvalues,v2_subA,0,ev2_subA);
  tge_v2_subA->SetMarkerColor(kBlack);
  tge_v2_subA->SetMarkerStyle(kOpenCircle);
  tge_v2_subA->Draw("ap");
  if ( verbose ) c1->Print(Form("PlotFigs/syh_ppref_v2_%s_subA.png",name));

  TGraphErrors* tge_c2_subAZ = new TGraphErrors(nptbins,ptvalues,c2_subAZ,0,ec2_subAZ);
  tge_c2_subAZ->SetMarkerColor(kBlack);
  tge_c2_subAZ->SetMarkerStyle(kOpenSquare);
  tge_c2_subAZ->Draw("ap");
  if ( verbose ) c1->Print(Form("PlotFigs/syh_ppref_c2_%s_subAZ.png",name));

  TGraphErrors* tge_v2_subAZ = new TGraphErrors(nptbins,ptvalues,v2_subAZ,0,ev2_subAZ);
  tge_v2_subAZ->SetMarkerColor(kBlack);
  tge_v2_subAZ->SetMarkerStyle(kOpenSquare);
  tge_v2_subAZ->Draw("ap");
  if ( verbose ) c1->Print(Form("PlotFigs/syh_ppref_v2_%s_subAZ.png",name));

  TGraphErrors* tge_c2_subRZ = new TGraphErrors(nptbins,ptvalues,c2_subRZ,0,ec2_subRZ);
  tge_c2_subRZ->SetMarkerColor(kBlack);
  tge_c2_subRZ->SetMarkerStyle(kOpenCross);
  tge_c2_subRZ->Draw("ap");
  if ( verbose ) c1->Print(Form("PlotFigs/syh_ppref_c2_%s_subRZ.png",name));

  TGraphErrors* tge_v2_subRZ = new TGraphErrors(nptbins,ptvalues,v2_subRZ,0,ev2_subRZ);
  tge_v2_subRZ->SetMarkerColor(kBlack);
  tge_v2_subRZ->SetMarkerStyle(kOpenCross);
  tge_v2_subRZ->Draw("ap");
  if ( verbose ) c1->Print(Form("PlotFigs/syh_ppref_v2_%s_subRZ.png",name));

  TGraphErrors* tge_c2_subR = new TGraphErrors(nptbins,ptvalues,c2_subR,0,ec2_subR);
  tge_c2_subR->SetMarkerColor(kBlack);
  tge_c2_subR->SetMarkerStyle(kOpenDiamond);
  tge_c2_subR->SetMarkerSize(2.5);
  tge_c2_subR->Draw("ap");
  if ( verbose ) c1->Print(Form("PlotFigs/syh_ppref_c2_%s_subR.png",name));

  TGraphErrors* tge_v2_subR = new TGraphErrors(nptbins,ptvalues,v2_subR,0,ev2_subR);
  tge_v2_subR->SetMarkerColor(kBlack);
  tge_v2_subR->SetMarkerStyle(kOpenDiamond);
  tge_v2_subR->SetMarkerSize(2.5);
  tge_v2_subR->Draw("ap");
  if ( verbose ) c1->Print(Form("PlotFigs/syh_ppref_v2_%s_subR.png",name));

  double seyoungdata_CNT_BBCS_FVTS_raw[6] = {0.01999,0.04210,0.07368,0.09789,0.12,0.14210};
  double seyoungdata_CNT_FVTN_FVTS_raw[6] = {0.03368,0.06421,0.11157,0.15263,0.20210,0.28947};
  double seyoungdata_CNT_FVTN_FVTS_sub[6] = {0.03052,0.05157,0.05684,0.06,0.01999,0};

  double ppg191_pt[13] = {0.49594,0.69400,0.89383,1.09417,1.29436,1.49426,1.69450,1.89463,2.09479,2.29479,2.49508,2.69540,2.89522};
  double ppg191_v2[13] = {0.029465,0.041451,0.054005,0.064461,0.076305,0.088801,0.098890,0.105504,0.114177,0.124270,0.130828,0.138556,0.137665};
  double ppg191_su[13] = {0.0021042,0.0029602,0.0038567,0.0046034,0.0054492,0.0063416,0.0070621,0.0075344,0.0083902,0.0091319,0.0096138,0.0101817,0.0101163};
  double ppg191_sd[13] = {0.0039422,0.0060338,0.0078612,0.0117858,0.0139513,0.0162361,0.0205304,0.0219035,0.0257567,0.0280336,0.0295129,0.0328588,0.0326475};

  TGraphAsymmErrors* tgae_ppg191 = new TGraphAsymmErrors(13,ppg191_pt,ppg191_v2,0,0,ppg191_sd,ppg191_su);
  tgae_ppg191->SetLineColor(kRed);
  tgae_ppg191->SetFillColorAlpha(kRed,0.35);


  TGraph* tg_sy_CNT_BBCS_FVTS_raw = new TGraph(nptbins,ptvalues,seyoungdata_CNT_BBCS_FVTS_raw);
  tg_sy_CNT_BBCS_FVTS_raw->SetLineColor(kBlue);
  TGraph* tg_sy_CNT_FVTN_FVTS_raw = new TGraph(nptbins,ptvalues,seyoungdata_CNT_FVTN_FVTS_raw);
  tg_sy_CNT_FVTN_FVTS_raw->SetLineColor(kBlue);
  //TGraph* tg_sy_CNT_FVTN_FVTS_sub = new TGraph(nptbins,ptvalues,seyoungdata_CNT_FVTN_FVTS_sub);
  TGraph* tg_sy_CNT_FVTN_FVTS_sub = new TGraph(5,ptvalues,seyoungdata_CNT_FVTN_FVTS_sub);
  tg_sy_CNT_FVTN_FVTS_sub->SetLineColor(kBlue);

  double xmin = 0.0;
  double xmax = 4.0;
  double ymin = -0.05;
  double ymax = 0.35;
  //arguments(index_CNT_BBCS,index_CNT_FVTS,index_FVTS_BBCS,"CNT_BBCS_FVTS");
  bool is_CNT_BBCS_FVTS = ( indexA == index_CNT_BBCS && indexB == index_CNT_FVTS && indexC == index_FVTS_BBCS );
  bool is_CNT_FVTN_FVTS = ( indexA == index_CNT_FVTN && indexB == index_CNT_FVTS && indexC == index_FVTN_FVTS );
  if ( is_CNT_BBCS_FVTS ) ymin = -0.02;
  if ( is_CNT_BBCS_FVTS ) ymax = 0.2;
  TH2D* hdummy = new TH2D("hdummy","",1,xmin,xmax,1,ymin,ymax);
  hdummy->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hdummy->GetYaxis()->SetTitle("v_{2}");
  hdummy->Draw();
  tgae_ppg191->Draw("L3");
  tge_v2_raw->Draw("p");
  if ( is_CNT_BBCS_FVTS ) tg_sy_CNT_BBCS_FVTS_raw->Draw("l");
  if ( is_CNT_FVTN_FVTS ) tg_sy_CNT_FVTN_FVTS_raw->Draw("l");
  TLegend* leg = new TLegend(0.18,0.65,0.38,0.92);
  leg->SetHeader(name);
  if ( is_CNT_BBCS_FVTS ) leg->AddEntry(tg_sy_CNT_BBCS_FVTS_raw,"raw v_{2} (Seyoung)","l");
  if ( is_CNT_FVTN_FVTS ) leg->AddEntry(tg_sy_CNT_FVTN_FVTS_raw,"raw v_{2} (Seyoung)","l");
  leg->AddEntry(tge_v2_raw,"raw v_{2}","p");
  leg->Draw();
  TLine* line = new TLine(xmin,0,xmax,0);
  line->SetLineWidth(2);
  line->SetLineStyle(2);
  line->Draw();
  c1->Print(Form("PlotFigs/syh_ppref_v2_%s.png",name));
  tge_v2_subA->Draw("p");
  if ( is_CNT_FVTN_FVTS ) tg_sy_CNT_FVTN_FVTS_sub->Draw("l");
  leg->AddEntry(tge_v2_subA,"sub v_{2} (ATLAS)","p");
  c1->Print(Form("PlotFigs/syh_ppref_v2_%s_sub1.png",name));
  tge_v2_subAZ->Draw("p");
  leg->AddEntry(tge_v2_subAZ,"sub v_{2} (ATLAS, ZYAM)","p");
  c1->Print(Form("PlotFigs/syh_ppref_v2_%s_sub2.png",name));
  tge_v2_subR->Draw("p");
  leg->AddEntry(tge_v2_subR,"sub v_{2} (Reference)","p");
  c1->Print(Form("PlotFigs/syh_ppref_v2_%s_sub3.png",name));

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
  tge_c2_subA->Draw("p");
  tge_c2_subAZ->Draw("p");
  tge_c2_subR->Draw("p");
  delete leg;
  leg = new TLegend(0.18,0.72,0.38,0.92);
  leg->SetHeader(name);
  leg->AddEntry(tge_c2_raw,"raw c_{2}","p");
  leg->AddEntry(tge_c2_subA,"sub c_{2} (ATLAS)","p");
  leg->AddEntry(tge_c2_subAZ,"sub c_{2} (ATLAS,ZYAM)","p");
  leg->AddEntry(tge_c2_subR,"sub c_{2} (Reference)","p");
  leg->Draw();
  line->Draw();
  c1->Print(Form("PlotFigs/syh_ppref_c2_%s.png",name));

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

// pT      v2       stat. unc. sys. up   sys. down
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
