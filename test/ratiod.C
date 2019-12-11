void ratiod()
{

  const int ntypes = 6;
  const int nlowmult = 2;
  const int nptbins = 6;

  double ptvalues[nptbins] = {0.35,0.7,1.2,1.7,2.35,3.5};

  string handle[ntypes];
  handle[0] = "CNT_BBCN_BBCS";
  handle[1] = "CNT_BBCN_FVTN";
  handle[2] = "CNT_BBCN_FVTS";
  handle[3] = "CNT_BBCS_FVTN";
  handle[4] = "CNT_BBCS_FVTS";
  handle[5] = "CNT_FVTN_FVTS";

  //double pt[nlowmult][ntypes][nptbins];
  double v2_raw[nlowmult][ntypes][nptbins];
  double ev2_raw[nlowmult][ntypes][nptbins];
  double v2_subA[nlowmult][ntypes][nptbins];
  double ev2_subA[nlowmult][ntypes][nptbins];
  double v2_subAZ[nlowmult][ntypes][nptbins];
  double ev2_subAZ[nlowmult][ntypes][nptbins];
  double v2_subRZ[nlowmult][ntypes][nptbins];
  double ev2_subRZ[nlowmult][ntypes][nptbins];
  double v2_subR[nlowmult][ntypes][nptbins];
  double ev2_subR[nlowmult][ntypes][nptbins];

  ifstream fin;
  for ( int ilowmult = 0; ilowmult < nlowmult; ++ilowmult )
    {
      for ( int itype = 0; itype < ntypes; ++itype )
        {
          int ptstop = nptbins;
          char fin_name[50];
          sprintf ( fin_name, "DataTextFiles/data_pAu_%s.txt", handle[itype].c_str() );
          if ( ilowmult == 1 )
            {
              //ptstop = 5;
              sprintf ( fin_name, "DataTextFiles/data_ppref_%s.txt", handle[itype].c_str() );
            }
          cout << fin_name << endl;
          fin.open(fin_name);
          for ( int ipt = 0; ipt < ptstop; ++ipt )
            {
              //cout << ilowmult << " " << itype << " " << ipt << endl;
              double x;
              fin >> x
                  >> v2_raw[ilowmult][itype][ipt]   >> ev2_raw[ilowmult][itype][ipt]
                  >> v2_subA[ilowmult][itype][ipt]  >> ev2_subA[ilowmult][itype][ipt]
                  >> v2_subAZ[ilowmult][itype][ipt] >> ev2_subAZ[ilowmult][itype][ipt]
                  >> v2_subRZ[ilowmult][itype][ipt] >> ev2_subRZ[ilowmult][itype][ipt]
                  >> v2_subR[ilowmult][itype][ipt]  >> ev2_subR[ilowmult][itype][ipt];
              //cout << ptvalues[ipt] << " " << x << endl;
            } // loop over pt values
          fin.close();
        } // loop over types
    } // loop over low mult ref

  cout << "Got the data..." << endl;

  // ---

  TCanvas* c1 = new TCanvas("c1","",1000,600);
  //TH2D* hdummy = new TH2D("hdummy","",1,0.0,4.0,1,-10.0,10.0);
  TH2D* hdummy = new TH2D("hdummy","",1,0.0,4.0,1,-1.0,3.0);
  hdummy->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hdummy->GetYaxis()->SetTitle("Ratio v_{2} peripheral-pAu/pp");
  hdummy->Draw();
  TLine* line = new TLine(0.0,1.0,4.0,1.0);
  line->SetLineWidth(2);
  line->SetLineStyle(2);
  line->Draw();
  TLegend* leg = new TLegend(0.58,0.58,0.88,0.92);
  leg->SetFillStyle(0);

  int color[6] = {kBlack,kBlue,kRed,kGreen+2,kMagenta+2,kOrange};
  int marker[6] = {kFullCircle,kFullSquare,kFullTriangleUp,kFullTriangleDown,kFullDiamond,kFullCross};

  double lmratio_raw[ntypes][nptbins];
  double elmratio_raw[ntypes][nptbins];
  double lmratio_subR[ntypes][nptbins];
  double elmratio_subR[ntypes][nptbins];
  TGraphErrors* tge_lmratio_raw[ntypes];
  TGraphErrors* tge_lmratio_subR[ntypes];
  cout << "starting type loop" << endl;
  for ( int itype = 0; itype < ntypes; ++itype )
    {
      cout << "starting pt loop" << endl;
      for ( int ipt = 0; ipt < nptbins-1; ++ipt )
        {
          lmratio_raw[itype][ipt] = v2_raw[0][itype][ipt]/v2_raw[1][itype][ipt];
          lmratio_subR[itype][ipt] = v2_subR[0][itype][ipt]/v2_subR[1][itype][ipt];
          elmratio_raw[itype][ipt] = lmratio_raw[itype][ipt] *
            sqrt(
                 pow(ev2_raw[0][itype][ipt]/v2_raw[0][itype][ipt],2.0)
                 +
                 pow(ev2_raw[1][itype][ipt]/v2_raw[1][itype][ipt],2.0)
                 );
          elmratio_subR[itype][ipt] = lmratio_subR[itype][ipt] *
            sqrt(
                 pow(ev2_subR[0][itype][ipt]/v2_subR[0][itype][ipt],2.0)
                 +
                 pow(ev2_subR[1][itype][ipt]/v2_subR[1][itype][ipt],2.0)
                 );
          cout << "calculated ratio " << v2_raw[0][itype][ipt] << "/" << v2_raw[1][itype][ipt] << " = " << lmratio_raw[itype][ipt] << endl;
          cout << "calculated ratio " << v2_subR[0][itype][ipt] << "/" << v2_subR[1][itype][ipt] << " = " << lmratio_subR[itype][ipt] << endl;
        }
      tge_lmratio_raw[itype] = new TGraphErrors(nptbins-1,ptvalues,lmratio_raw[itype],0,elmratio_raw[itype]);
      tge_lmratio_raw[itype]->SetMarkerStyle(marker[itype]);
      tge_lmratio_raw[itype]->SetMarkerColor(color[itype]);
      //tge_lmratio_raw[itype]->Draw("p");
      tge_lmratio_subR[itype] = new TGraphErrors(nptbins-1,ptvalues,lmratio_subR[itype],0,elmratio_subR[itype]);
      tge_lmratio_subR[itype]->SetMarkerStyle(marker[itype]);
      tge_lmratio_subR[itype]->SetMarkerColor(color[itype]);
      tge_lmratio_subR[itype]->Draw("p");
      leg->AddEntry(tge_lmratio_subR[itype],handle[itype].c_str(),"p");
    }
  leg->Draw();
  c1->Print("PlotFigs/fig_pAu_lm_ratio.png");

  // --

  TCanvas* c2 = new TCanvas("c2","",1000,600);
  TH2D* hdummy1 = new TH2D("hdummy1","",1,0.0,4.0,1,-1.0,4.0);
  hdummy1->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hdummy1->GetYaxis()->SetTitle("Ratio v_{2} type/CNT_FVTN_FVTS");
  c1->cd();
  hdummy1->Draw();
  TH2D* hdummy2 = new TH2D("hdummy2","",1,0.0,4.0,1,-1.0,10.0);
  hdummy2->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hdummy2->GetYaxis()->SetTitle("Ratio v_{2} type/CNT_FVTN_FVTS");
  c2->cd();
  hdummy2->Draw();
  // TLine* line = new TLine(0.0,1.0,4.0,1.0);
  // line->SetLineWidth(2);
  // line->SetLineStyle(2);
  // line->Draw();
  delete leg;
  //leg = new TLegend(0.58,0.58,0.88,0.92);
  leg = new TLegend(0.18,0.58,0.48,0.92);

  // int color[6] = {kBlack,kBlue,kRed,kGreen+2,kMagenta+2,kOrange};
  // int marker[6] = {kFullCircle,kFullSquare,kFullTriangleUp,kFullTriangleDown,kFullDiamond,kFullCross};

  int type_anchor = 5; // can change as needed

  double typeratio_raw[ntypes][nptbins];
  double etyperatio_raw[ntypes][nptbins];
  double typeratio_subR[ntypes][nptbins];
  double etyperatio_subR[ntypes][nptbins];
  TGraphErrors* tge_typeratio_raw[ntypes];
  TGraphErrors* tge_typeratio_subR[ntypes];
  cout << "starting type loop" << endl;
  for ( int itype = 0; itype < ntypes; ++itype )
    {
      if ( itype == type_anchor ) continue;
      cout << "starting pt loop" << endl;
      for ( int ipt = 0; ipt < nptbins; ++ipt )
        {
          typeratio_raw[itype][ipt] = v2_raw[0][itype][ipt]/v2_raw[0][type_anchor][ipt];
          typeratio_subR[itype][ipt] = v2_subR[0][itype][ipt]/v2_subR[0][type_anchor][ipt];
          etyperatio_raw[itype][ipt] = 0; // do this later
          etyperatio_subR[itype][ipt] = 0; // do this later
          etyperatio_raw[itype][ipt] = typeratio_raw[itype][ipt] *
            sqrt(
                 pow(ev2_raw[0][itype][ipt]/v2_raw[0][type_anchor][ipt],2.0)
                 +
                 pow(ev2_raw[0][itype][ipt]/v2_raw[0][type_anchor][ipt],2.0)
                 );
          etyperatio_subR[itype][ipt] = typeratio_subR[itype][ipt] *
            sqrt(
                 pow(ev2_subR[0][itype][ipt]/v2_subR[0][type_anchor][ipt],2.0)
                 +
                 pow(ev2_subR[0][itype][ipt]/v2_subR[0][type_anchor][ipt],2.0)
                 );
          cout << "calculated ratio " << v2_raw[0][itype][ipt] << "/" << v2_raw[0][type_anchor][ipt] << " = " << typeratio_raw[itype][ipt] << endl;
          cout << "calculated ratio " << v2_subR[0][itype][ipt] << "/" << v2_subR[0][type_anchor][ipt] << " = " << typeratio_subR[itype][ipt] << endl;
        }
      c1->cd();
      tge_typeratio_raw[itype] = new TGraphErrors(nptbins,ptvalues,typeratio_raw[itype],0,etyperatio_raw[itype]);
      tge_typeratio_raw[itype]->SetMarkerStyle(marker[itype]);
      tge_typeratio_raw[itype]->SetMarkerColor(color[itype]);
      tge_typeratio_raw[itype]->Draw("p");
      c2->cd();
      tge_typeratio_subR[itype] = new TGraphErrors(nptbins,ptvalues,typeratio_subR[itype],0,etyperatio_subR[itype]);
      tge_typeratio_subR[itype]->SetMarkerStyle(marker[itype]);
      tge_typeratio_subR[itype]->SetMarkerColor(color[itype]);
      tge_typeratio_subR[itype]->Draw("p");
      leg->AddEntry(tge_typeratio_subR[itype],handle[itype].c_str(),"p");
    }
  // ---
  c1->cd();
  line->Draw();
  leg->Draw();
  c1->Print("PlotFigs/fig_pAu_raw_type_ratio.png");
  // ---
  c2->cd();
  line->Draw();
  leg->Draw();
  c2->Print("PlotFigs/fig_pAu_subR_type_ratio.png");

  TGraphErrors* tge_v2_raw[ntypes];
  TGraphErrors* tge_v2_subR[ntypes];
  cout << "starting type loop" << endl;
  delete hdummy1;
  hdummy1 = new TH2D("hdummy1","",1,0.0,4.0,1,-0.1,0.35);
  hdummy1->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hdummy1->GetYaxis()->SetTitle("v_{2}");
  c1->cd();
  hdummy1->Draw();
  delete hdummy2;
  hdummy2 = new TH2D("hdummy2","",1,0.0,4.0,1,-0.1,0.35);
  hdummy2->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hdummy2->GetYaxis()->SetTitle("v_{2}^{ref}");
  c2->cd();
  hdummy2->Draw();
  delete leg;
  //leg = new TLegend(0.58,0.58,0.88,0.92);
  leg = new TLegend(0.18,0.58,0.28,0.92);
  leg->SetFillStyle(0);
  delete line;
  line = new TLine(0.0,0.0,4.0,0.0);
  line->SetLineWidth(2);
  line->SetLineStyle(2);
  for ( int itype = 0; itype < ntypes; ++itype )
    {
      c1->cd();
      tge_v2_raw[itype] = new TGraphErrors(nptbins,ptvalues,v2_raw[0][itype],0,ev2_raw[0][itype]);
      tge_v2_raw[itype]->SetMarkerStyle(marker[itype]);
      tge_v2_raw[itype]->SetMarkerColor(color[itype]);
      tge_v2_raw[itype]->Draw("p");
      c2->cd();
      tge_v2_subR[itype] = new TGraphErrors(nptbins,ptvalues,v2_subR[0][itype],0,ev2_subR[0][itype]);
      tge_v2_subR[itype]->SetMarkerStyle(marker[itype]);
      tge_v2_subR[itype]->SetMarkerColor(color[itype]);
      tge_v2_subR[itype]->Draw("p");
      leg->AddEntry(tge_v2_subR[itype],handle[itype].c_str(),"p");
    }
  // ---
  c1->cd();
  line->Draw();
  leg->Draw();
  c1->Print("PlotFigs/fig_pAu_raw_v2_comparison.png");
  // ---
  c2->cd();
  line->Draw();
  leg->Draw();
  c2->Print("PlotFigs/fig_pAu_subR_v2_comparison.png");

}
