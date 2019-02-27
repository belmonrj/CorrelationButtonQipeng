#ifndef ANAHELPER_H
#define ANAHELPER_H

void plotText(Double_t x,Double_t y,Color_t color, const char *text, Double_t tsize=0.04);
void plotMarkerText(Double_t x,Double_t y,Int_t color,Int_t mstyle,const char *text,Float_t msize=2., Double_t tsize=0.04);
void plotMarkerLineText(     Double_t x, Double_t y,Float_t msize,Int_t mcolor,Int_t mstyle,Int_t lcolor,Int_t lstyle, const char *text, Double_t tsize=0.04, bool EX0 = false);

void PrintHistBinning(const TH1* _hist) {
    if (_hist) {
        int _Nbins_Nch = _hist->GetXaxis()->GetNbins();
        const double* _dNch_range = _hist->GetXaxis()->GetXbins()->GetArray();
        for (int _index = 0; _index < _Nbins_Nch + 1; _index++) {
            cout << _dNch_range[_index];
            if (_index < _Nbins_Nch) cout << ", ";
            else cout << endl;
        }
    }
}

TH1* Symmetrize(const TH1* h_dphi) {

    int nbins = h_dphi->GetNbinsX();
    int nbins_half = nbins/2.;
    
    TH1D* m_dphi = new TH1D("m_dphi", "Symmetrized correlation", nbins_half, 0,TMath::Pi());

    //cout << "total number of bins: " << nbins << endl;
    //cout << "half number of bins: " << nbins_half << endl;
    for (int i = 1; i < nbins_half+1; i++) {
        if (i < nbins_half/2+1) {
	    //cout << "combine " << nbins_half/2+i << " with " << nbins_half/2-i+1 << endl;
            m_dphi->SetBinContent(i,h_dphi->GetBinContent(nbins_half/2+i) + h_dphi->GetBinContent(nbins_half/2-i+1));
            float error = TMath::Sqrt( pow(h_dphi->GetBinError(nbins_half/2+i),2) + pow(h_dphi->GetBinError(nbins_half/2-i+1),2) );
            m_dphi->SetBinError(i,error);
        } else if (i > nbins_half/2) {
	    //cout << "combine " << nbins_half/2+i << " with " << nbins-(i-nbins_half/2-1) << endl;
            m_dphi->SetBinContent(i,h_dphi->GetBinContent(nbins_half/2+i) + h_dphi->GetBinContent(nbins-(i-nbins_half/2-1)));
            float error = TMath::Sqrt( pow(h_dphi->GetBinError(nbins_half/2+i),2) + pow(h_dphi->GetBinError(nbins-(i-nbins_half/2-1)),2) );
            m_dphi->SetBinError(i,error);
        }
    }
    m_dphi->Scale(0.5);
    return m_dphi;
}


void ZYAMF(TH1F* hist) {

    double min = hist->GetBinContent(hist->GetMinimumBin());
    for (int i=1; i<hist->GetNbinsX()+1; i++) {
        hist->SetBinContent(i, hist->GetBinContent(i) - min);
    }
}

void ZYAM(TH1* hist) {

    double min = hist->GetBinContent(hist->GetMinimumBin());
    for (int i=1; i<hist->GetNbinsX()+1; i++) {
        hist->SetBinContent(i, hist->GetBinContent(i) - min);
    }
}


void plotMarkerLineText(Double_t x, Double_t y,Float_t msize,Int_t mcolor,Int_t mstyle,Int_t lcolor,Int_t lstyle, const char *text, Double_t tsize, bool EX0)
{
  //tsize = 0.04;
  TMarker *marker = new TMarker(x-(0.55*tsize),y,8);
  marker->SetMarkerColor(mcolor);  marker->SetNDC();
  marker->SetMarkerStyle(mstyle);
  marker->SetMarkerSize(msize);

  TLatex l; l.SetTextAlign(12); l.SetTextSize(tsize);
  l.SetNDC();
  l.DrawLatex(x,y,text);

  Double_t y1=y-0.5*tsize;
  Double_t y2=y+0.5*tsize;
  Double_t x2=x-0.15*tsize;
  Double_t x1=x-0.95*tsize;

  TLine mline;
  mline.SetLineWidth(2);
  mline.SetLineColor(lcolor);
  mline.SetLineStyle(lstyle);
  Double_t y_new=(y1+y2)/2.;
  Double_t x_new=(x1+x2)/2.;

  double size_y = 0.3*tsize;
  double size_x = 0.28*tsize;
  //double size = msize;
  if (msize!=0) {
    if (!EX0) {
      mline.DrawLineNDC(x1,       y_new,x1+size_x,y_new);
      mline.DrawLineNDC(x2-size_x,y_new,x2,       y_new);
    }
    mline.DrawLineNDC(x_new,y2-size_y,x_new,y2);
    mline.DrawLineNDC(x_new,y1,x_new,y1+size_y);
  } else {
    mline.DrawLineNDC(x-0.95*tsize,y_new,x-0.15*tsize,y_new);
  }
  marker->Draw();
}


void plotMarkerText(Double_t x,Double_t y,Int_t color,Int_t mstyle, const char *text,Float_t msize, Double_t tsize)
{
  TMarker *marker = new TMarker(x-(0.55*tsize),y,8);
  marker->SetMarkerColor(color);  marker->SetNDC();
  marker->SetMarkerStyle(mstyle);
  marker->SetMarkerSize(msize);
  marker->Draw();

  TLatex l; l.SetTextAlign(12); l.SetTextSize(tsize);
  l.SetNDC();
  l.DrawLatex(x,y,text);

  Double_t y1=y-0.5*tsize;
  Double_t y2=y+0.5*tsize;
  Double_t x2=x-0.15*tsize;
  Double_t x1=x-0.95*tsize;

  TLine mline;
  mline.SetLineWidth(2);
  mline.SetLineColor(color);
  mline.SetLineStyle(1);
  Double_t y_new=(y1+y2)/2.;
  Double_t x_new=(x1+x2)/2.;

  double size_y = 0.3*tsize;
  double size_x = 0.28*tsize;
  //double size = msize;
  //mline.DrawLineNDC(x1,y_new,x1+size_x,y_new);
  //mline.DrawLineNDC(x2-size_x,y_new,x2,y_new);
  mline.DrawLineNDC(x_new,y2-size_y,x_new,y2);
  mline.DrawLineNDC(x_new,y1,x_new,y1+size_y);
}


void plotText(Double_t x,Double_t y,Color_t color, const char *text, Double_t tsize) {

  TLatex l; l.SetTextAlign(12); l.SetTextSize(tsize);
  l.SetNDC();
  l.SetTextColor(color);
  l.DrawLatex(x,y,text);
}
#endif
