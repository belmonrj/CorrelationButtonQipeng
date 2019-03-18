#ifndef CORRELATIONMAKER_H
#define CORRELATIONMAKER_H

#include "FlowAnaConfig.h"

class CorrelationMaker {
    
    const FlowAnaConfig* m_anaConfig;
    bool m_debug;
    int m_rebin;

    // index to show which trigger weight function to be applied;
    // 0 for no trigger weight
    // 1 for p+Pb HF analysis 
    // supports for other analysis to be added
    int m_weightCorrIndex; 

    float m_etaBinInterval;
    float m_etaBoundary;

    float m_etaBinIndexLow;
    float m_etaBinIndexHigh;
    float m_jetEtaIndexLow;
    float m_jetEtaIndexHigh;

    float m_exclude_low;
    float m_exclude_high;
    bool  m_setExclude;
    // in case one wants to check the 1D signal and 1D mix
    TH1* m_h1_sig;
    TH1* m_h1_mix;

  public:
    CorrelationMaker();
    virtual ~CorrelationMaker() {};
    void init();

    void setWeightIndex (int _n) { m_weightCorrIndex = _n; }
    int  getWeightIndex () { return m_weightCorrIndex; }
    void setRebin (int _n) { m_rebin = _n; }
    void setDebug (bool _debug = true) { m_debug = _debug; }
    void setConfig (const FlowAnaConfig* _config) { m_anaConfig = _config; }

    // interface for J/psi analysis to exclude signal region
    void setExclude(float _low, float _high) {m_exclude_low = _low; m_exclude_high = _high; m_setExclude = true;};

    const FlowAnaConfig* getConfig () const { return m_anaConfig; }

    void get1DSigHist() const {return m_h1_sig;}
    void get1DMixHist() const {return m_h1_mix;}

    void generateHist_bkgSub();

    float weight2016HFAna(int index_Nch);

    // main function for generating single 1D correlation histogram
    // method1 -- ATLAS mixed event normalization, projection first, then take ratio (default ATLAS method)
    // method2 -- ATLAS mixed event normalization, take 2D ratio, then 1D projection
    // method3 -- CMS mixed event normalization, projection first, then take ratio
    // method4 -- CMS mixed event normalization, take 2D ratio, then 1D projjection (default CMS method)
    // method5 -- no mixed event correction  (for systematic study of mixing)
    // by default doing analaysis in pt and Nch
    virtual TH1* MakeCorr(float _pt_low, float _pt_high, float _Nch_low, float _Nch_high, int method = 1);

    // for UPC analysis
    virtual TH1* MakeCorrUpc(float _pt_low, float _pt_high, float _Nch_low, float _Nch_high, int method = 1);

    // Additional dimensions
    // for D meson, invariant mass
    // for muon, momentum imbalance
    // for UPC, gap
    virtual TH1* MakeCorr(float _pt_low, float _pt_high, float _Nch_low, float _Nch_high, float _add_low, float _add_high, int method = 1);

    virtual TH2* Make2DCorrUpc(float _pt_low, float _pt_high, float _Nch_low, float _Nch_high, int method = 1);
    //virtual TH2* Make2DCorr(float _pt_low, float _pt_high, float _Nch_low, float _Nch_high);
    //virtual TH2* Make2DCorr(float _pt_low, float _pt_high, float _Nch_low, float _Nch_high, float _add_low, float _add_high);
    
    void Plot2DCorrelation(TH2* h1, TCanvas* c1);
};
#endif


CorrelationMaker::CorrelationMaker() {
    m_anaConfig = 0;
    m_debug = false;
    m_rebin = 1;
    m_weightCorrIndex = 0;

    m_etaBinInterval = 0.1;
    m_etaBoundary = 5.0;

    m_h1_sig = 0;
    m_h1_mix = 0;

    m_exclude_low = 0;
    m_exclude_high = 0;;
    m_setExclude = false;
}



void CorrelationMaker::init() {
    m_anaConfig->print();
    m_etaBoundary    = m_anaConfig->getCorrEtaBoundary();
    m_etaBinInterval = m_anaConfig->getCorrEtaInterval();

    // two index to exclude deta range due to upper boundary cut
    m_etaBinIndexLow  = (m_etaBoundary - m_anaConfig->getEtaRangeHigh() )/m_etaBinInterval + 1;
    m_etaBinIndexHigh = 2.0*m_etaBoundary/m_etaBinInterval - (m_etaBoundary - m_anaConfig->getEtaRangeHigh())/m_etaBinInterval;

    // two index to exclude deta range due to lower boundary cut
    m_jetEtaIndexLow = (m_anaConfig->getEtaRangeHigh() - m_anaConfig->getEtaRangeLow())/m_etaBinInterval + m_etaBinIndexLow;
    m_jetEtaIndexHigh = m_etaBinIndexHigh - (m_anaConfig->getEtaRangeHigh() - m_anaConfig->getEtaRangeLow())/m_etaBinInterval;
}



float CorrelationMaker::weight2016HFAna(int index_Nch) {
    // use 1./lumi for this analysis
    float _lumi = 1.;
    if (index_Nch<14) _lumi = 0.082;
    else if (index_Nch<20) _lumi = 0.687;
    else if (index_Nch<24) _lumi = 4.394;
    else _lumi = 55.163;
    return 1./_lumi;
}


// a little bit wired to have multiplicity cut as float instead of int
// however, just keep it this way for now
TH1* CorrelationMaker::MakeCorr(float _pt_low, float _pt_high, float _Nch_low, float _Nch_high, int method) {

    TH2* h1 = 0;
    TH2* h2 = 0;

    float nsig = 0;
    float nsigMix = 0;

    string path_sig = m_anaConfig->getCorrHistPathSame();
    string path_mix = m_anaConfig->getCorrHistPathMix();

    m_anaConfig->inputFile()->cd();
    TH2* h2_nsig    = (TH2*)gDirectory->Get(m_anaConfig->getTrigYieldHistName().c_str());
    TH2* h2_nsigMix = (TH2*)gDirectory->Get(m_anaConfig->getMixTrigYieldHistName().c_str());

    if (!h2_nsig) cout << "Cannot retrieve yield hitogram for PTY calculation" << endl;
    if (!h2_nsigMix) cout << "Cannot retrieve mix-event yield hitogram for PTY calculation" << endl;
    if (m_debug) {h2_nsig->Print(); h2_nsigMix->Print();}
    if (m_debug) cout << "Merging the following slices" << endl;

    int index_Nch_low = 0;
    int index_Nch_high = 0;
    for (int iNch=1; iNch<m_anaConfig->getInputMultiBinning()->GetXaxis()->GetNbins() + 1; iNch++){
        if (m_anaConfig->getInputMultiBinning()->GetXaxis()->GetBinLowEdge(iNch) == _Nch_low)  index_Nch_low  = iNch - 1; 
        if (m_anaConfig->getInputMultiBinning()->GetXaxis()->GetBinUpEdge(iNch)  == _Nch_high) index_Nch_high = iNch - 1; 
    }
    if (m_debug) cout << "index_Nch_low = " << index_Nch_low << ",\tindex_Nch_high = " << index_Nch_high << endl;


    int index_pt_low = 0;
    int index_pt_high = 0;
    for (int ipt=1; ipt<m_anaConfig->getInputPtBinning()->GetXaxis()->GetNbins() + 1; ipt++){
        if (m_anaConfig->getInputPtBinning()->GetXaxis()->GetBinLowEdge(ipt) == _pt_low)  index_pt_low  = ipt - 1; 
        if (m_anaConfig->getInputPtBinning()->GetXaxis()->GetBinUpEdge(ipt)  == _pt_high) index_pt_high = ipt - 1; 
    }

    for (int iNch=index_Nch_low; iNch<index_Nch_high+1; iNch++) {
        for (int ipt=index_pt_low; ipt<index_pt_high+1; ipt++) {

	        TH2* htemp_1 = (TH2*)gDirectory->Get(Form("%sNch%d_pt%d",path_sig.c_str(),iNch,ipt));
	        TH2* htemp_2 = (TH2*)gDirectory->Get(Form("%sNch%d_pt%d",path_mix.c_str(),iNch,ipt));

            float _weight = 1.;
            if (getWeightIndex() == 1) {
                _weight = weight2016HFAna(iNch);
            }
            // to be extend to support other analysis

            htemp_1->Scale(_weight);
            htemp_2->Scale(_weight);
            nsig += h2_nsig->GetBinContent(ipt+1, iNch+1) * _weight;
            nsigMix += h2_nsigMix->GetBinContent(ipt+1, iNch+1) * _weight;

	        if (iNch==index_Nch_low && ipt==index_pt_low) {
	            h1 = (TH2*) htemp_1->Clone(Form("h1_Nch%d_%d", index_Nch_low, index_Nch_high));
	            h2 = (TH2*) htemp_2->Clone(Form("h2_Nch%d_%d", index_Nch_low, index_Nch_high));
	        } else {
                h1->Add(htemp_1);
                h2->Add(htemp_2);
	        }
	        delete htemp_1;
	        delete htemp_2;
	        if (m_debug) cout << "pt index = " << ipt << ", Nch index = " << iNch << endl;

        }
    }
    if (m_debug) cout << "Merged same event entries: " << h1->GetEntries()/1e6 << "M" << endl;

    TH1* hphi;

    if (m_h1_sig) {delete m_h1_sig; m_h1_sig = 0;}
    if (m_h1_mix) {delete m_h1_mix; m_h1_mix = 0;}

    // nbins to project
    int nbins_low  = (m_jetEtaIndexLow-1) - m_etaBinIndexLow + 1;
    int nbins_high = m_etaBinIndexHigh - (m_jetEtaIndexHigh+1) + 1;
    if (nbins_low != nbins_high) cout << nbins_low << ",\t " << nbins_high << endl;

    if (method == 1) {
        // ATLAS style mixed event normalization
        // make 1D project of signal and mix first
        // then take ratio of two 1D hist

        // Make projections first
        m_h1_sig    = h1->ProjectionY(Form("h_sig_phi_Nch%d",  index_Nch_low), m_etaBinIndexLow,    m_jetEtaIndexLow-1, "e" );
        m_h1_sig->Add(h1->ProjectionY(Form("h_sig_phi2_Nch%d", index_Nch_low), m_jetEtaIndexHigh+1, m_etaBinIndexHigh,  "e" ) );
        m_h1_mix    = h2->ProjectionY(Form("h_mix_phi_Nch%d",  index_Nch_low), m_etaBinIndexLow,    m_jetEtaIndexLow-1, "e" );
        m_h1_mix->Add(h2->ProjectionY(Form("h_mix_phi2_Nch%d", index_Nch_low), m_jetEtaIndexHigh+1, m_etaBinIndexHigh,  "e" ) );

        float int_S = m_h1_sig->Integral();

        hphi = (TH1*)m_h1_sig->Clone(Form("_hphi"));
        hphi->Divide(m_h1_mix);

        float int_C = hphi->Integral();
        float K = int_S / int_C;
        hphi->Scale(K/nsig, "width"); //per trigger & per deltaPhi yield

    } else if (method == 2) {
        // ATLAS style mixed event normalization
        // take 2D ratio, then project to 1D
        
        // start with N_{pair}
        m_h1_sig    = h1->ProjectionY(Form("h_sig_phi_Nch%d",  index_Nch_low), m_etaBinIndexLow,    m_jetEtaIndexLow-1, "e" );
        m_h1_sig->Add(h1->ProjectionY(Form("h_sig_phi2_Nch%d", index_Nch_low), m_jetEtaIndexHigh+1, m_etaBinIndexHigh,  "e" ) );
        m_h1_mix    = h2->ProjectionY(Form("h_mix_phi_Nch%d",  index_Nch_low), m_etaBinIndexLow,    m_jetEtaIndexLow-1, "e" );
        m_h1_mix->Add(h2->ProjectionY(Form("h_mix_phi2_Nch%d", index_Nch_low), m_jetEtaIndexHigh+1, m_etaBinIndexHigh,  "e" ) );
        float int_S = m_h1_sig->Integral();

        h1->Scale(1./h1->Integral());
        h2->Scale(1./h2->Integral());
        TH2F* h_correlation = (TH2F*)h1->Clone();
        h_correlation->Divide(h2);

        hphi    = h_correlation->ProjectionY(Form("h_phi_Nch%d", index_Nch_low), m_etaBinIndexLow, m_jetEtaIndexLow-1, "e");
        hphi->Add(h_correlation->ProjectionY(Form("h_phi2_Nch%d",index_Nch_low), m_jetEtaIndexHigh+1, m_etaBinIndexHigh, "e"));

        float int_C = hphi->Integral();
        float K = int_S / int_C;
        hphi->Scale(K/nsig, "width"); //per trigger & per deltaPhi yield

    } else if (method == 3) { 
        // CMS style mixed event normalization
        // make 1D project of signal and mix first
        // then take ratio of two 1D hist
        h1->Scale(1./nsig,    "width");
        h2->Scale(1./nsigMix, "width");
        int bin_00 = h2->FindBin(0,0);
        float content_00 = h2->GetBinContent(bin_00);
        h2->Scale( 1./content_00 );

        // Make projections first
        m_h1_sig    = h1->ProjectionY(Form("h_sig_phi_Nch%d",  index_Nch_low), m_etaBinIndexLow,    m_jetEtaIndexLow-1, "e" );
        m_h1_sig->Add(h1->ProjectionY(Form("h_sig_phi2_Nch%d", index_Nch_low), m_jetEtaIndexHigh+1, m_etaBinIndexHigh,  "e" ) );
        m_h1_mix    = h2->ProjectionY(Form("h_mix_phi_Nch%d",  index_Nch_low), m_etaBinIndexLow,    m_jetEtaIndexLow-1, "e" );
        m_h1_mix->Add(h2->ProjectionY(Form("h_mix_phi2_Nch%d", index_Nch_low), m_jetEtaIndexHigh+1, m_etaBinIndexHigh,  "e" ) );

        m_h1_sig->Scale(1./(nbins_low + nbins_high));
        m_h1_mix->Scale(1./(nbins_low + nbins_high));

        hphi = (TH1*)m_h1_sig->Clone(Form("_hphi"));
        hphi->Divide(m_h1_mix);
    } else if (method == 4) {
        // CMS style mixed event normalization
        // take 2D ratio, then project to 1D
        h1->Scale(1./nsig,    "width");
        h2->Scale(1./nsigMix, "width");
        int bin_00 = h2 ->FindBin(0,0);
        float content_00 = h2->GetBinContent(bin_00);
        h2->Scale( 1./content_00 );
        TH2F* h_correlation = (TH2F*)h1->Clone();
        h_correlation->Divide(h2);

        hphi    = h_correlation->ProjectionY(Form("h_sig_phi_Nch%d", index_Nch_low), m_etaBinIndexLow, m_jetEtaIndexLow-1, "e");
        hphi->Add(h_correlation->ProjectionY(Form("h_sig_phi2_Nch%d",index_Nch_low), m_jetEtaIndexHigh+1, m_etaBinIndexHigh, "e"));
        hphi->Scale(1./(nbins_low + nbins_high)); 
        

    } else if (method == 5) {
        // no mixing applied
        // copied from method 1

        // Make projections first
        m_h1_sig    = h1->ProjectionY(Form("h_sig_phi_Nch%d",  index_Nch_low), m_etaBinIndexLow,    m_jetEtaIndexLow-1, "e" );
        m_h1_sig->Add(h1->ProjectionY(Form("h_sig_phi2_Nch%d", index_Nch_low), m_jetEtaIndexHigh+1, m_etaBinIndexHigh,  "e" ) );

        hphi = (TH1*)m_h1_sig->Clone(Form("_hphi"));
        hphi->Scale(1./nsig, "width"); //per trigger & per deltaPhi yield

    } else {
        cout << "method is not supported, please choose from 1~3" << endl;
    }

    if (m_debug) cout << "------------------------------------------" << endl;
    delete h1;
    delete h2;

    hphi->Rebin(m_rebin);
    hphi->Scale(1./m_rebin);
    hphi->SetName(Form("h_pty_dphi_gap%.1fto%.1f_Nch%.0fto%.0f_pt%.0fto%.0f",m_anaConfig->getEtaRangeLow(), m_anaConfig->getEtaRangeHigh(), _Nch_low,_Nch_high,_pt_low,_pt_high));

    return hphi;
}



TH1* CorrelationMaker::MakeCorrUpc(float _pt_low, float _pt_high, float _Nch_low, float _Nch_high, int method) {

    TH2D* h1 = 0;
    TH2D* h2 = 0;

    string path_sig = m_anaConfig->getCorrHistPathSame();
    string path_mix = m_anaConfig->getCorrHistPathMix();

    if (m_debug) cout << "Merging the following slices" << endl;
    int index_Nch_low = 0;
    int index_Nch_high = 0;
    for (int iNch=1; iNch<m_anaConfig->getInputMultiBinning()->GetXaxis()->GetNbins() + 1; iNch++){
        if (m_anaConfig->getInputMultiBinning()->GetXaxis()->GetBinLowEdge(iNch) == _Nch_low)  index_Nch_low  = iNch - 1; 
        if (m_anaConfig->getInputMultiBinning()->GetXaxis()->GetBinUpEdge(iNch)  == _Nch_high) index_Nch_high = iNch - 1; 
    }


    int index_pt_low = 0;
    int index_pt_high = 0;
    for (int ipt=1; ipt<m_anaConfig->getInputPtBinning()->GetXaxis()->GetNbins() + 1; ipt++){
        if (m_anaConfig->getInputPtBinning()->GetXaxis()->GetBinLowEdge(ipt) == _pt_low)  index_pt_low  = ipt - 1; 
        if (m_anaConfig->getInputPtBinning()->GetXaxis()->GetBinUpEdge(ipt)  == _pt_high) index_pt_high = ipt - 1; 
    }
    if (m_debug) {
        cout << "index_Nch_low = " << index_Nch_low << ",\tindex_Nch_high = " << index_Nch_high << endl;
        cout << "index_pt_low  = " << index_pt_low << ", \tindex_pt_high  = " << index_pt_high << endl;
    }

    m_anaConfig->inputFile()->cd();
    TH2F* h2_nsig = (TH2F*)gDirectory->Get(m_anaConfig->getTrigYieldHistName().c_str());
    if (!h2_nsig) cout << "Cannot retrieve yield hitogram for PTY calculation" << endl;
    if (m_debug) {
        cout << "pt binning of h2_nsig: " << endl;
        int ix = 1;
        for (; ix < h2_nsig->GetNbinsX()+1; ix++) {
            cout << h2_nsig->GetXaxis()->GetBinLowEdge(ix) << ", ";
        }
        cout << h2_nsig->GetXaxis()->GetBinUpEdge(ix-1) << endl;

        cout << "Multiplicity binning of h2_nsig: " << endl;
        ix = 1;
        for (; ix < h2_nsig->GetNbinsY()+1; ix++) {
            cout << h2_nsig->GetYaxis()->GetBinLowEdge(ix) << ", ";
        }
        cout << h2_nsig->GetYaxis()->GetBinUpEdge(ix-1) << endl;
    }

    m_anaConfig->inputFile()->cd();
    float nsig = 0;
    for (int iNch=index_Nch_low; iNch<index_Nch_high+1; iNch++) {
        for (int ipt=index_pt_low; ipt<index_pt_high+1; ipt++) {

            // example names
            //h_deta_dphi_same_Mq4_ptq2_Pq1;1	
            //h2_deta_dphi_mixed_Mq4_ptq2_Sq0_tq0_Pq1;1 
	        TH2D* htemp_1 = (TH2D*)gDirectory->Get(Form("%sMq%d_ptq%d_Pq1",path_sig.c_str(),iNch,ipt));
	        TH2D* htemp_2 = (TH2D*)gDirectory->Get(Form("%sMq%d_ptq%d_Sq0_tq0_Pq1",path_mix.c_str(),iNch,ipt));

            nsig += h2_nsig->GetBinContent(ipt+1, iNch+1);

	        if (iNch==index_Nch_low && ipt==index_pt_low) {
	            h1 = (TH2D*) htemp_1->Clone(Form("h1_Nch%d_%d_pt%d_pt%d", index_Nch_low, index_Nch_high, index_pt_low, index_pt_high));
	            h2 = (TH2D*) htemp_2->Clone(Form("h2_Nch%d_%d_pt%d_pt%d", index_Nch_low, index_Nch_high, index_pt_low, index_pt_high));
	        } else {
                h1->Add(htemp_1);
                h2->Add(htemp_2);
	        }

            //double _width1 = m_anaConfig->getCorrEtaInterval();
            //double _width2 = h1->GetXaxis()->GetBinWidth(1);
            //if (_width2 != _width1) {
            //    cout << "Input bin width = " << h1->GetXaxis()->GetBinWidth(1) << endl;
            //    cout << "tool configuration = " << m_anaConfig->getCorrEtaInterval() << endl;
            //}

	        delete htemp_1;
	        delete htemp_2;
            if (m_debug) {
                cout << "pt: " << h2_nsig->GetXaxis()->GetBinLowEdge(ipt+1) << " ~ " << h2_nsig->GetXaxis()->GetBinUpEdge(ipt+1) << endl;
                cout << "nch: " << h2_nsig->GetYaxis()->GetBinLowEdge(iNch+1) << " ~ " << h2_nsig->GetYaxis()->GetBinUpEdge(iNch+1) << endl;
	            cout << "pt index = " << ipt << ", Nch index = " << iNch << endl;
	            cout << "ntrig = " << h2_nsig->GetBinContent(ipt+1, iNch+1) << endl;
            }
        }
    }

    TH1* hphi;
    if (m_debug) cout << "Merged same event entries: " << h1->GetEntries()/1e6 << "M" << endl;
    if (m_h1_sig) {delete m_h1_sig; m_h1_sig = 0;}
    if (m_h1_mix) {delete m_h1_mix; m_h1_mix = 0;}

    // nbins to project
    int nbins_low  = (m_jetEtaIndexLow-1) - m_etaBinIndexLow + 1;
    int nbins_high = m_etaBinIndexHigh - (m_jetEtaIndexHigh+1) + 1;

    if (nbins_low != nbins_high) cout << nbins_low << ",\t " << nbins_high << endl;
    if (m_debug) {
        cout << endl;
        cout << "Debug the gap index determination: " << endl;
        cout << m_etaBinIndexLow << ",\t" << m_jetEtaIndexLow -1 
             << " <==> " << h1->GetXaxis()->GetBinLowEdge(m_etaBinIndexLow) << ",\t" << h1->GetXaxis()->GetBinUpEdge(m_jetEtaIndexLow -1) << endl;
        cout << m_jetEtaIndexHigh+1 << ",\t" << m_etaBinIndexHigh
             << " <==> " << h1->GetXaxis()->GetBinLowEdge(m_jetEtaIndexHigh+1) << ",\t" << h1->GetXaxis()->GetBinUpEdge(m_etaBinIndexHigh) << endl;
        cout << "\tMerged Nbins negative side =  " << nbins_low << ", positive side =  " << nbins_high << endl;
        cout << endl;
    }

    if (method == 1) {
        // ATLAS style mixed event normalization
        // make 1D project of signal and mix first
        // then take ratio of two 1D hist

        // Make projections first
        m_h1_sig    = h1->ProjectionY(Form("h_sig_phi_Nch%d",  index_Nch_low), m_etaBinIndexLow,    m_jetEtaIndexLow-1, "e" );
        m_h1_sig->Add(h1->ProjectionY(Form("h_sig_phi2_Nch%d", index_Nch_low), m_jetEtaIndexHigh+1, m_etaBinIndexHigh,  "e" ) );
        m_h1_mix    = h2->ProjectionY(Form("h_mix_phi_Nch%d",  index_Nch_low), m_etaBinIndexLow,    m_jetEtaIndexLow-1, "e" );
        m_h1_mix->Add(h2->ProjectionY(Form("h_mix_phi2_Nch%d", index_Nch_low), m_jetEtaIndexHigh+1, m_etaBinIndexHigh,  "e" ) );

        float int_S = m_h1_sig->Integral();

        hphi = (TH1*)m_h1_sig->Clone(Form("_hphi"));
        hphi->Divide(m_h1_mix);

        float int_C = hphi->Integral();
        float K = int_S / int_C;
        if (m_debug) {
            cout << "int_S = " << int_S << endl;
            cout << "int_C = " << int_C << endl;
        }
        hphi->Scale(K/nsig, "width"); //per trigger & per deltaPhi yield

    } else if (method == 2) {
        // ATLAS style mixed event normalization
        // take 2D ratio, then project to 1D
        
        // start with N_{pair}
        m_h1_sig    = h1->ProjectionY(Form("h_sig_phi_Nch%d",  index_Nch_low), m_etaBinIndexLow,    m_jetEtaIndexLow-1, "e" );
        m_h1_sig->Add(h1->ProjectionY(Form("h_sig_phi2_Nch%d", index_Nch_low), m_jetEtaIndexHigh+1, m_etaBinIndexHigh,  "e" ) );
        m_h1_mix    = h2->ProjectionY(Form("h_mix_phi_Nch%d",  index_Nch_low), m_etaBinIndexLow,    m_jetEtaIndexLow-1, "e" );
        m_h1_mix->Add(h2->ProjectionY(Form("h_mix_phi2_Nch%d", index_Nch_low), m_jetEtaIndexHigh+1, m_etaBinIndexHigh,  "e" ) );
        float int_S = m_h1_sig->Integral();

        h1->Scale(1./h1->Integral());
        h2->Scale(1./h2->Integral());
        TH2D* h_correlation = (TH2D*)h1->Clone();
        h_correlation->Divide(h2);

        hphi    = h_correlation->ProjectionY(Form("h_phi_Nch%d", index_Nch_low), m_etaBinIndexLow, m_jetEtaIndexLow-1, "e");
        hphi->Add(h_correlation->ProjectionY(Form("h_phi2_Nch%d",index_Nch_low), m_jetEtaIndexHigh+1, m_etaBinIndexHigh, "e"));

        float int_C = hphi->Integral();
        float K = int_S / int_C;
        hphi->Scale(K/nsig, "width"); //per trigger & per deltaPhi yield

    } else if (method == 3) { 
        // CMS style mixed event normalization
        // make 1D project of signal and mix first
        // then take ratio of two 1D hist
        h1->Scale(1./nsig, "width");
        h2->Scale(1., "width");
        int bin_00 = h2->FindBin(0,0);
        float content_00 = h2->GetBinContent(bin_00);
        h2->Scale( 1./content_00 );

        // Make projections first
        m_h1_sig    = h1->ProjectionY(Form("h_sig_phi_Nch%d",  index_Nch_low), m_etaBinIndexLow,    m_jetEtaIndexLow-1, "e" );
        m_h1_sig->Add(h1->ProjectionY(Form("h_sig_phi2_Nch%d", index_Nch_low), m_jetEtaIndexHigh+1, m_etaBinIndexHigh,  "e" ) );
        m_h1_mix    = h2->ProjectionY(Form("h_mix_phi_Nch%d",  index_Nch_low), m_etaBinIndexLow,    m_jetEtaIndexLow-1, "e" );
        m_h1_mix->Add(h2->ProjectionY(Form("h_mix_phi2_Nch%d", index_Nch_low), m_jetEtaIndexHigh+1, m_etaBinIndexHigh,  "e" ) );

        m_h1_sig->Scale(1./(nbins_low + nbins_high));
        m_h1_mix->Scale(1./(nbins_low + nbins_high));

        hphi = (TH1*)m_h1_sig->Clone(Form("_hphi"));
        hphi->Divide(m_h1_mix);
    } else if (method == 4) {
        // CMS style mixed event normalization
        // take 2D ratio, then project to 1D
        h1->Scale(1./nsig, "width");
        h2->Scale(1., "width");
        int bin_00 = h2 ->FindBin(0,0);
        float content_00 = h2->GetBinContent(bin_00);
        h2->Scale( 1./content_00 );
        TH2D* h_correlation = (TH2D*)h1->Clone();
        h_correlation->Divide(h2);

        hphi    = h_correlation->ProjectionY(Form("h_sig_phi_Nch%d", index_Nch_low), m_etaBinIndexLow, m_jetEtaIndexLow-1, "e");
        hphi->Add(h_correlation->ProjectionY(Form("h_sig_phi2_Nch%d",index_Nch_low), m_jetEtaIndexHigh+1, m_etaBinIndexHigh, "e"));
        hphi->Scale(1./(nbins_low + nbins_high)); 
        
    } else if (method == 5) {
        // no mixing applied
        // copied from method 1

        // Make projections first
        
        if (m_debug) { 
            cout << "Integral before projection: " << h1->Integral() << endl;
        }
        m_h1_sig    = h1->ProjectionY(Form("h_sig_phi_Nch%d",  index_Nch_low), m_etaBinIndexLow,    m_jetEtaIndexLow-1, "e" );
        m_h1_sig->Add(h1->ProjectionY(Form("h_sig_phi2_Nch%d", index_Nch_low), m_jetEtaIndexHigh+1, m_etaBinIndexHigh,  "e" ) );
        if (m_debug) { 
            cout << "Integral after projection: " << m_h1_sig->Integral() << endl;
        }

        hphi = (TH1*)m_h1_sig->Clone(Form("_hphi"));
        hphi->Scale(1./nsig, "width"); //per trigger & per deltaPhi yield

        if (m_debug) { 
            cout << "nsig = " << nsig << endl;
            cout << "per trigger yield = " << hphi->Integral("width") << endl;
        }

    } else {
        cout << "method is not supported, please choose from 1~3" << endl;
    }

    if (m_debug) cout << "------------------------------------------" << endl;
    delete h1;
    delete h2;

    hphi->Rebin(m_rebin);
    hphi->Scale(1./m_rebin);
    hphi->SetName(Form("h_pty_dphi_gap%.1fto%.1f_Nch%.0fto%.0f_pt%.0fto%.0f",m_anaConfig->getEtaRangeLow(), m_anaConfig->getEtaRangeHigh(), _Nch_low,_Nch_high,_pt_low,_pt_high));

    return hphi;
}



TH2* CorrelationMaker::Make2DCorrUpc(float _pt_low, float _pt_high, float _Nch_low, float _Nch_high, int method) {

    TH2* h1 = 0;
    TH2* h2 = 0;

    string path_sig = m_anaConfig->getCorrHistPathSame();
    string path_mix = m_anaConfig->getCorrHistPathMix();

    if (m_debug) cout << "Merging the following slices" << endl;
    int index_Nch_low = 0;
    int index_Nch_high = 0;
    for (int iNch=1; iNch<m_anaConfig->getInputMultiBinning()->GetXaxis()->GetNbins() + 1; iNch++){
        if (m_anaConfig->getInputMultiBinning()->GetXaxis()->GetBinLowEdge(iNch) == _Nch_low)  index_Nch_low  = iNch - 1; 
        if (m_anaConfig->getInputMultiBinning()->GetXaxis()->GetBinUpEdge(iNch)  == _Nch_high) index_Nch_high = iNch - 1; 
    }


    int index_pt_low = 0;
    int index_pt_high = 0;
    for (int ipt=1; ipt<m_anaConfig->getInputPtBinning()->GetXaxis()->GetNbins() + 1; ipt++){
        if (m_anaConfig->getInputPtBinning()->GetXaxis()->GetBinLowEdge(ipt) == _pt_low)  index_pt_low  = ipt - 1; 
        if (m_anaConfig->getInputPtBinning()->GetXaxis()->GetBinUpEdge(ipt)  == _pt_high) index_pt_high = ipt - 1; 
    }
    if (m_debug) {
        cout << "index_Nch_low = " << index_Nch_low << ",\tindex_Nch_high = " << index_Nch_high << endl;
        cout << "index_pt_low  = " << index_pt_low << ", \tindex_pt_high  = " << index_pt_high << endl;
    }

    m_anaConfig->inputFile()->cd();
    TH2F* h2_nsig = (TH2F*)gDirectory->Get(m_anaConfig->getTrigYieldHistName().c_str());
    if (!h2_nsig) cout << "Cannot retrieve yield hitogram for PTY calculation" << endl;

    float nsig = 0;
    for (int iNch=index_Nch_low; iNch<index_Nch_high+1; iNch++) {
        for (int ipt=index_pt_low; ipt<index_pt_high+1; ipt++) {

	        TH2* htemp_1 = (TH2*)gDirectory->Get(Form("%sMq%d_ptq%d_Pq1",path_sig.c_str(),iNch,ipt));
	        TH2* htemp_2 = (TH2*)gDirectory->Get(Form("%sMq%d_ptq%d_Sq0_tq0_Pq1",path_mix.c_str(),iNch,ipt));

            nsig += h2_nsig->GetBinContent(ipt+1, iNch+1);

	        if (iNch==index_Nch_low && ipt==index_pt_low) {
	            h1 = (TH2*) htemp_1->Clone(Form("h1_Nch%d_%d_pt%d_pt%d", index_Nch_low, index_Nch_high, index_pt_low, index_pt_high));
	            h2 = (TH2*) htemp_2->Clone(Form("h2_Nch%d_%d_pt%d_pt%d", index_Nch_low, index_Nch_high, index_pt_low, index_pt_high));
	        } else {
                h1->Add(htemp_1);
                h2->Add(htemp_2);
	        }
	        delete htemp_1;
	        delete htemp_2;
	        if (m_debug) cout << "pt index = " << ipt << ", Nch index = " << iNch << endl;
        }
    }

    TH2* h_correlation = (TH2*) h1->Clone(Form("h2_pty_dphi_gap%.1fto%.1f_Nch%.0fto%.0f_pt%.0fto%.0f",m_anaConfig->getEtaRangeLow(), m_anaConfig->getEtaRangeHigh(), _Nch_low,_Nch_high,_pt_low,_pt_high));

    if (method == 1) {
        // with mixing correction
        float int_S = h1->Integral();
        h1->Scale(1./h1->Integral());
        h2->Scale(1./h2->Integral());
        h_correlation->Divide(h2);

        float int_C = h_correlation->Integral();
        float K = int_S / int_C;
        if (m_debug) {
            cout << h_correlation->Integral() << endl;
            cout << "int_S = " << int_S << endl;
            cout << "int_C = " << int_C << endl;
            cout << "K = " << K << endl;
        }
        cout << h_correlation->Integral() << endl;
        h_correlation->Scale(K/nsig, "width"); //per trigger & per deltaPhi yield 
        cout << h_correlation->Integral() << endl;

    } else if (method == 2) {
        // no mixxing correction
        h_correlation->Scale(1./nsig, "width"); 
    } else {
        cout << "method is not supported, please choose from 1 or 2" << endl;
    }

    if (m_debug) cout << "------------------------------------------" << endl;

    delete h1;
    delete h2;

    return h_correlation;
}



// a little bit wired to have multiplicity cut as float instead of int
// however, just keep it this way for now
// used for 3D case
TH1* CorrelationMaker::MakeCorr(float _pt_low, float _pt_high, float _Nch_low, float _Nch_high, float _add_low, float _add_high, int method) {

    TH2* h1 = 0;
    TH2* h2 = 0;

    float nsig = 0;
    float nsigMix = 0;

    string path_sig = m_anaConfig->getCorrHistPathSame();
    string path_mix = m_anaConfig->getCorrHistPathMix();

    m_anaConfig->inputFile()->cd();
    TH3* h2_nsig    = (TH3*)gDirectory->Get(m_anaConfig->getTrigYieldHistName().c_str());
    TH3* h2_nsigMix = (TH3*)gDirectory->Get(m_anaConfig->getMixTrigYieldHistName().c_str());

    if (!h2_nsig) {
        cout << "Cannot retrieve yield hitogram for PTY calculation" << endl;
        cout << "object list in file:" << endl;
        m_anaConfig->inputFile()->ls();
        cout << "target histogram name:" << m_anaConfig->getTrigYieldHistName().c_str() << endl;
    }
    if (!h2_nsigMix) cout << "Cannot retrieve mix-event yield hitogram for PTY calculation" << endl;
    if (m_debug) {h2_nsig->Print(); h2_nsigMix->Print();}
    if (m_debug) cout << "Merging the following slices" << endl;

    int index_Nch_low = 0;
    int index_Nch_high = 0;
    for (int iNch=1; iNch<m_anaConfig->getInputMultiBinning()->GetXaxis()->GetNbins() + 1; iNch++){
        if (m_anaConfig->getInputMultiBinning()->GetXaxis()->GetBinLowEdge(iNch) == _Nch_low)  index_Nch_low  = iNch - 1; 
        if (m_anaConfig->getInputMultiBinning()->GetXaxis()->GetBinUpEdge(iNch)  == _Nch_high) index_Nch_high = iNch - 1; 
    }


    int index_pt_low = 0;
    int index_pt_high = 0;
    for (int ipt=1; ipt<m_anaConfig->getInputPtBinning()->GetXaxis()->GetNbins() + 1; ipt++){
        if (m_anaConfig->getInputPtBinning()->GetXaxis()->GetBinLowEdge(ipt) == _pt_low)  index_pt_low  = ipt - 1; 
        if (m_anaConfig->getInputPtBinning()->GetXaxis()->GetBinUpEdge(ipt)  == _pt_high) index_pt_high = ipt - 1; 
    }


    int index_add_low = 0;
    int index_add_high = 0;
    for (int iadd=1; iadd<m_anaConfig->getInputThirdBinning()->GetXaxis()->GetNbins() + 1; iadd++){
        if (m_anaConfig->getInputThirdBinning()->GetXaxis()->GetBinLowEdge(iadd) == _add_low)  index_add_low  = iadd - 1; 
        if (m_anaConfig->getInputThirdBinning()->GetXaxis()->GetBinUpEdge(iadd)  == _add_high) index_add_high = iadd - 1; 
    }

    int _exclude_index_add_low = 0;
    int _exclude_index_add_high = 0;
    if (m_setExclude) {
        for (int iadd=1; iadd<m_anaConfig->getInputThirdBinning()->GetXaxis()->GetNbins() + 1; iadd++){
            if (m_anaConfig->getInputThirdBinning()->GetXaxis()->GetBinLowEdge(iadd) == m_exclude_low)  _exclude_index_add_low  = iadd - 1; 
            if (m_anaConfig->getInputThirdBinning()->GetXaxis()->GetBinUpEdge(iadd)  == m_exclude_high) _exclude_index_add_high = iadd - 1; 
        }
    }

    if (m_debug) {
        cout << "index_Nch_low = " << index_Nch_low << ",\tindex_Nch_high = " << index_Nch_high << endl;
        cout << "index_pt_low = "  << index_pt_low  << ",\tindex_pt_high = "  << index_pt_high  << endl;
        cout << "index_3rd_low = " << index_add_low << ",\tindex_3rd_high = " << index_add_high << endl;
        cout << "exclude index = " << _exclude_index_add_low << ",\t" << _exclude_index_add_high << endl;
    }


    for (int iNch=index_Nch_low; iNch<index_Nch_high+1; iNch++) {
        for (int ipt=index_pt_low; ipt<index_pt_high+1; ipt++) {
            for (int iadd=index_add_low; iadd<index_add_high+1; iadd++) {
                if (m_setExclude && iadd >= _exclude_index_add_low && iadd <= _exclude_index_add_high) continue;

	            TH2* htemp_1 = (TH2*)gDirectory->Get(Form("%sNch%d_pt%d_%s%d",path_sig.c_str(),iNch,ipt, m_anaConfig->getThirdDimensionName().c_str(), iadd));
	            TH2* htemp_2 = (TH2*)gDirectory->Get(Form("%sNch%d_pt%d_%s%d",path_mix.c_str(),iNch,ipt, m_anaConfig->getThirdDimensionName().c_str(), iadd));

                float _weight = 1.;
                if (getWeightIndex() == 1) {
                    _weight = weight2016HFAna(iNch);
                }
                // to be extend to support other analysis

                htemp_1->Scale(_weight);
                htemp_2->Scale(_weight);
                nsig += h2_nsig->GetBinContent(ipt+1, iNch+1, iadd+1) * _weight;
                nsigMix += h2_nsigMix->GetBinContent(ipt+1, iNch+1, iadd+1) * _weight;

	            if (iNch==index_Nch_low && ipt==index_pt_low && iadd==index_add_low) {
	                h1 = (TH2*) htemp_1->Clone(Form("h1_Nch%d_%d", index_Nch_low, index_Nch_high));
	                h2 = (TH2*) htemp_2->Clone(Form("h2_Nch%d_%d", index_Nch_low, index_Nch_high));
	            } else {
                    h1->Add(htemp_1);
                    h2->Add(htemp_2);
	            }
	            delete htemp_1;
	            delete htemp_2;
	            if (m_debug) cout << "pt index = " << ipt 
                                  << ", Nch index = " << iNch 
                                  << ", " << m_anaConfig->getThirdDimensionName().c_str() << " index = " << iadd 
                                  << endl;
            }
        }
    }
    if (m_debug) cout << "Merged same event entries: " << h1->GetEntries()/1e6 << "M" << endl;

    TH1* hphi;

    if (m_h1_sig) {delete m_h1_sig; m_h1_sig = 0;}
    if (m_h1_mix) {delete m_h1_mix; m_h1_mix = 0;}

    // nbins to project
    int nbins_low  = (m_jetEtaIndexLow-1) - m_etaBinIndexLow + 1;
    int nbins_high = m_etaBinIndexHigh - (m_jetEtaIndexHigh+1) + 1;
    if (nbins_low != nbins_high) cout << nbins_low << ",\t " << nbins_high << endl;

    if (method == 1) {
        // ATLAS style mixed event normalization
        // make 1D project of signal and mix first
        // then take ratio of two 1D hist

        // Make projections first
        m_h1_sig    = h1->ProjectionY(Form("h_sig_phi_Nch%d",  index_Nch_low), m_etaBinIndexLow,    m_jetEtaIndexLow-1, "e" );
        m_h1_sig->Add(h1->ProjectionY(Form("h_sig_phi2_Nch%d", index_Nch_low), m_jetEtaIndexHigh+1, m_etaBinIndexHigh,  "e" ) );
        m_h1_mix    = h2->ProjectionY(Form("h_mix_phi_Nch%d",  index_Nch_low), m_etaBinIndexLow,    m_jetEtaIndexLow-1, "e" );
        m_h1_mix->Add(h2->ProjectionY(Form("h_mix_phi2_Nch%d", index_Nch_low), m_jetEtaIndexHigh+1, m_etaBinIndexHigh,  "e" ) );

        float int_S = m_h1_sig->Integral();

        hphi = (TH1*)m_h1_sig->Clone(Form("_hphi"));
        hphi->Divide(m_h1_mix);

        float int_C = hphi->Integral();
        float K = int_S / int_C;
        hphi->Scale(K/nsig, "width"); //per trigger & per deltaPhi yield

    } else if (method == 2) {
        // ATLAS style mixed event normalization
        // take 2D ratio, then project to 1D
        
        // start with N_{pair}
        m_h1_sig    = h1->ProjectionY(Form("h_sig_phi_Nch%d",  index_Nch_low), m_etaBinIndexLow,    m_jetEtaIndexLow-1, "e" );
        m_h1_sig->Add(h1->ProjectionY(Form("h_sig_phi2_Nch%d", index_Nch_low), m_jetEtaIndexHigh+1, m_etaBinIndexHigh,  "e" ) );
        m_h1_mix    = h2->ProjectionY(Form("h_mix_phi_Nch%d",  index_Nch_low), m_etaBinIndexLow,    m_jetEtaIndexLow-1, "e" );
        m_h1_mix->Add(h2->ProjectionY(Form("h_mix_phi2_Nch%d", index_Nch_low), m_jetEtaIndexHigh+1, m_etaBinIndexHigh,  "e" ) );
        float int_S = m_h1_sig->Integral();

        h1->Scale(1./h1->Integral());
        h2->Scale(1./h2->Integral());
        TH2F* h_correlation = (TH2F*)h1->Clone();
        h_correlation->Divide(h2);

        hphi    = h_correlation->ProjectionY(Form("h_phi_Nch%d", index_Nch_low), m_etaBinIndexLow, m_jetEtaIndexLow-1, "e");
        hphi->Add(h_correlation->ProjectionY(Form("h_phi2_Nch%d",index_Nch_low), m_jetEtaIndexHigh+1, m_etaBinIndexHigh, "e"));

        float int_C = hphi->Integral();
        float K = int_S / int_C;
        hphi->Scale(K/nsig, "width"); //per trigger & per deltaPhi yield

    } else if (method == 3) { 
        // CMS style mixed event normalization
        // make 1D project of signal and mix first
        // then take ratio of two 1D hist
        h1->Scale(1./nsig,    "width");
        h2->Scale(1./nsigMix, "width");
        int bin_00 = h2->FindBin(0,0);
        float content_00 = h2->GetBinContent(bin_00);
        h2->Scale( 1./content_00 );

        // Make projections first
        m_h1_sig    = h1->ProjectionY(Form("h_sig_phi_Nch%d",  index_Nch_low), m_etaBinIndexLow,    m_jetEtaIndexLow-1, "e" );
        m_h1_sig->Add(h1->ProjectionY(Form("h_sig_phi2_Nch%d", index_Nch_low), m_jetEtaIndexHigh+1, m_etaBinIndexHigh,  "e" ) );
        m_h1_mix    = h2->ProjectionY(Form("h_mix_phi_Nch%d",  index_Nch_low), m_etaBinIndexLow,    m_jetEtaIndexLow-1, "e" );
        m_h1_mix->Add(h2->ProjectionY(Form("h_mix_phi2_Nch%d", index_Nch_low), m_jetEtaIndexHigh+1, m_etaBinIndexHigh,  "e" ) );

        m_h1_sig->Scale(1./(nbins_low + nbins_high));
        m_h1_mix->Scale(1./(nbins_low + nbins_high));

        hphi = (TH1*)m_h1_sig->Clone(Form("_hphi"));
        hphi->Divide(m_h1_mix);
    } else if (method == 4) {
        // CMS style mixed event normalization
        // take 2D ratio, then project to 1D
        h1->Scale(1./nsig,    "width");
        h2->Scale(1./nsigMix, "width");
        int bin_00 = h2 ->FindBin(0,0);
        float content_00 = h2->GetBinContent(bin_00);
        h2->Scale( 1./content_00 );
        TH2F* h_correlation = (TH2F*)h1->Clone();
        h_correlation->Divide(h2);

        hphi    = h_correlation->ProjectionY(Form("h_sig_phi_Nch%d", index_Nch_low), m_etaBinIndexLow, m_jetEtaIndexLow-1, "e");
        hphi->Add(h_correlation->ProjectionY(Form("h_sig_phi2_Nch%d",index_Nch_low), m_jetEtaIndexHigh+1, m_etaBinIndexHigh, "e"));
        hphi->Scale(1./(nbins_low + nbins_high)); 
        

    } else if (method == 5) {
        // no mixing applied
        // copied from method 1

        // Make projections first
        m_h1_sig    = h1->ProjectionY(Form("h_sig_phi_Nch%d",  index_Nch_low), m_etaBinIndexLow,    m_jetEtaIndexLow-1, "e" );
        m_h1_sig->Add(h1->ProjectionY(Form("h_sig_phi2_Nch%d", index_Nch_low), m_jetEtaIndexHigh+1, m_etaBinIndexHigh,  "e" ) );

        hphi = (TH1*)m_h1_sig->Clone(Form("_hphi"));
        hphi->Scale(1./nsig, "width"); //per trigger & per deltaPhi yield

    } else {
        cout << "method is not supported, please choose from 1~3" << endl;
    }

    if (m_debug) cout << "------------------------------------------" << endl;
    delete h1;
    delete h2;

    hphi->Rebin(m_rebin);
    hphi->Scale(1./m_rebin);
    hphi->SetName(Form("h_pty_dphi_gap%.1fto%.1f_Nch%.0fto%.0f_pt%.0fto%.0f_%s%dto%d",
                         m_anaConfig->getEtaRangeLow(), m_anaConfig->getEtaRangeHigh(), 
                         _Nch_low,_Nch_high,
                         _pt_low,_pt_high,
                         m_anaConfig->getThirdDimensionName().c_str(),index_add_low,index_add_high
                       ));

    return hphi;
}


void CorrelationMaker::Plot2DCorrelation(TH2* h1, TCanvas* c1) {
    TPad* pad = (TPad*)c1->cd();
    h1->SetLineWidth(1);

    double min;
    double max;

    max =  h1->GetBinContent(h1->FindBin(0,3.14))*1.02;
    min = (h1->GetBinContent(h1->FindBin(3.5,0)) 
         + h1->GetBinContent(h1->FindBin(3.0,0))
         + h1->GetBinContent(h1->FindBin(-3.5,0)) 
         + h1->GetBinContent(h1->FindBin(-3.0,0))
           )/4.*.99;

    h1->GetXaxis()->SetTitle("#Delta#it{#eta}");
    h1->GetXaxis()->SetRangeUser(-4.5,4.5); // hard-coded for now
    h1->GetZaxis()->SetRangeUser(min,max);
    h1->GetZaxis()->SetNdivisions(505,kTRUE);
    h1->GetYaxis()->SetNdivisions(509,kTRUE);
    h1->GetXaxis()->SetNdivisions(505,kTRUE);
    h1->GetXaxis()->SetTitle("#Delta#it{#eta}");
    h1->GetYaxis()->SetTitle("#Delta#it{#phi}");
    h1->GetZaxis()->SetTitle("#it{Y}(#Delta#it{#phi}, #Delta#it{#eta})");
    h1->GetXaxis()->CenterTitle();
    h1->GetYaxis()->CenterTitle();
    h1->GetZaxis()->CenterTitle();
    h1->GetZaxis()->SetTitleOffset(1.20);
    h1->Draw("SURF1FB");
    ATLASLabel(0.02,0.94,"Internal");
    //myText(    0.02,0.89,1,"#it{p}+Pb #sqrt{#it{s}_{NN}} = 8.16 TeV");
    //myText(    0.02,0.83,1,"#it{h}-#it{h} Correlation");
    //myText(    0.02,0.12,1,"0.5 < #it{p}_{T}^{trig,assco} < 5 GeV");
    pad->SetTheta(60); // default is 30
    pad->SetPhi(45); // default is 30
    pad->Update();
}

/*
void CorrelationMaker::generateHist_multi() {

    fout->mkdir("PTYRaw");
    fout->mkdir("PTYRaw/SameEvtCorrelation");
    fout->mkdir("PTYRaw/MixedEvtCorrelation");

    const int _nBins =  m_anaConfig->getOutputMultiBinning()->GetXaxis()->GetNbins()
    const double* _multiOutput_range = m_anaConfig->getOutputMultiBinning()->GetXaxis()->GetXbins()->GetArray();

    for (int i=0; i<_nBins; i++) {

       	// index is the histogram index (0~Nbins-1) in the files
       	// the histogram index number is bin index - 1
       	float Nch_low  = _multiOutput_range[i];
       	float Nch_high = _multiOutput_range[i+1];

       	int index_Nch_low = 0;
       	int index_Nch_high = 0;
       	for (int iNch=1; iNch<hNch_binning->GetXaxis()->GetNbins() + 1; iNch++){
       	    if (hNch_binning->GetXaxis()->GetBinLowEdge(iNch) == Nch_low) index_Nch_low = iNch - 1; 
       	    if (hNch_binning->GetXaxis()->GetBinUpEdge(iNch) == Nch_high) index_Nch_high = iNch - 1; 
       	}


        {
            // targeting pt selection for trigger particles 
            // index is the histogram index (0~Nbins-1) in the files
            // the histogram index number is bin index - 1
            int index_pt_low = 0;
            int index_pt_high = 0;
            for (int ipt=1; ipt<hpt_binning->GetXaxis()->GetNbins() + 1; ipt++){
                if (hpt_binning->GetXaxis()->GetBinLowEdge(ipt) == ptCutLow_dNch) index_pt_low = ipt - 1; 
                if (hpt_binning->GetXaxis()->GetBinUpEdge(ipt)  == ptCutHigh_dNch) index_pt_high = ipt - 1; 
            }

            for (int method = 1; method < 3; method++) {
                for (int iwgt = 1; iwgt < 2; iwgt++) {
                    bool _isWeighted = true;
                    string path_wgt = "Wgt";
                    if (iwgt == 1) {
                        _isWeighted = false;
                        path_wgt = "Raw";
                    }
                    string path_method = "PTY";
                    if (method == 2) path_method = "Corr";

                    if (i==0) {
                        // peripheral reference
                        TH1* h_dphi_ref = MakeCorr(index_pt_low, index_pt_high, index_Nch_ref_low, index_Nch_ref_high, _isWeighted, method);
                        h_dphi_ref->SetNameTitle("dphi_Nch_ref","");
                        fout->cd(Form("%s%s",path_method.c_str(), path_wgt.c_str()));
                        h_dphi_ref->Write();
                    }

                    TH1* h_dphi = MakeCorr(index_pt_low, index_pt_high, index_Nch_low, index_Nch_high, _isWeighted, method);
                    h_dphi->SetNameTitle(Form("dphi_Nch%d",i),"");
                    fout->cd(Form("%s%s",path_method.c_str(), path_wgt.c_str()));
                    h_dphi->Write();

                    if (h1_sig) {
                        fout->cd(Form("%s%s/SameEvtCorrelation",path_method.c_str(), path_wgt.c_str()));
                        h1_sig->SetNameTitle(Form("dphi_sig_Nch%d",i),"");
                        h1_sig->Write();
                    }
                    if (h1_mix) {
                        fout->cd(Form("%s%s/MixedEvtCorrelation",path_method.c_str(), path_wgt.c_str()));
                        h1_mix->SetNameTitle(Form("dphi_mix_Nch%d",i),"");
                        h1_mix->Write();
                    }
                }
            }

        } // test histogram for pt slices
}

*/
