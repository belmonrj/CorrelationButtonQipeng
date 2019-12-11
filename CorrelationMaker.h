#ifndef CORRELATIONMAKER_H
#define CORRELATIONMAKER_H

#include "FlowAnaConfig.h"

class CorrelationMaker {
    
    const FlowAnaConfig* m_anaConfig;
    bool m_debug;
    int m_rebin;

    // Depth of the mixing pool
    // default is 5;
    int m_mixDepth;

    // index to show which trigger weight function to be applied;
    // 0 or > 4 for no trigger weight
    // 1 for p+Pb hh correlation analysis 
    // 2 for p+Pb Dh correlation analysis for L1_TE200
    // 3 for p+Pb high run 313295
    // 4 for p+Pb low run 313435
    // supports for other analysis to be added
    int m_weightCorrIndex; 

    // remove the double counted pair from the correlation function errors
    // false means no error correction
    // true means run error correction for correlation of asso_pt_low <= trig_pt <= asso_pt_high
    bool m_runStatCorr;

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

    void setStatCorrection (bool _flag) { m_runStatCorr = _flag; }
    bool getStatCorrection () { return m_runStatCorr; }

    void setRebin (int _n) { m_rebin = _n; }
    void setDebug (bool _debug = true) { m_debug = _debug; }
    void setConfig (const FlowAnaConfig* _config) { m_anaConfig = _config; }

    // interface for J/psi analysis to exclude signal region
    void setExclude(float _low, float _high) {m_exclude_low = _low; m_exclude_high = _high; m_setExclude = true;};

    const FlowAnaConfig* getConfig () const { return m_anaConfig; }

    TH1* get1DSigHist() const {return m_h1_sig;}
    TH1* get1DMixHist() const {return m_h1_mix;}

    void generateHist_bkgSub();

    float weight2016hhAna(int index_Nch);
    float weight2016DzeroAna(int index_Nch);
    float weight295(int index_Nch);
    float weight435(int index_Nch);
    float weight435AnaSel(int index_Nch);

    // main function for generating single 1D correlation histogram
    // method1 -- ATLAS mixed event normalization, projection first, then take ratio (default ATLAS method)
    // method2 -- ATLAS mixed event normalization, take 2D ratio, then 1D projection
    // method3 -- CMS mixed event normalization, projection first, then take ratio
    // method4 -- CMS mixed event normalization, take 2D ratio, then 1D projjection (default CMS method)
    // method5 -- no mixed event correction  (for systematic study of mixing)
    // by default doing analaysis in pt and Nch
    virtual TH1* MakeCorr(float _pt_low, float _pt_high, float _Nch_low, float _Nch_high, int method = 1);
    //virtual TH2* Make2DCorr(float _pt_low, float _pt_high, float _Nch_low, float _Nch_high);
    // Additional 3rd dimensions
    // for D meson, invariant mass
    // for muon, momentum imbalance
    virtual TH1* MakeCorr(float _pt_low, float _pt_high, float _Nch_low, float _Nch_high, float _add_low, float _add_high, int method = 1);

    // for UPC analysis
    virtual TH1*  MakeCorrUpc(float _pt_low, float _pt_high, float _Nch_low, float _Nch_high, int method = 1, int type = 1);
    virtual TH1*  MakeCorrUpc(float _pt_low, float _pt_high, float _Nch_low, float _Nch_high, float _add_low, float _add_high, int method = 1, int type = 1);
    virtual float GetPtMean  (float _pt_low, float _pt_high, float _Nch_low, float _Nch_high, int type = 1);
    virtual float GetMultMean(float _pt_low, float _pt_high, float _Nch_low, float _Nch_high, int type = 1);
    virtual float GetMultMean(float _pt_low, float _pt_high, float _Nch_low, float _Nch_high, float _add_low, float _add_high, int type = 1);

    virtual TH2* Make2DCorrUpc(float _pt_low, float _pt_high, float _Nch_low, float _Nch_high, int method = 1);
    virtual TH2* Make2DCorrUpc(float _pt_low, float _pt_high, float _Nch_low, float _Nch_high, float _add_low, float _add_high, int method = 1);
    //virtual TH2* Make2DCorr(float _pt_low, float _pt_high, float _Nch_low, float _Nch_high, float _add_low, float _add_high);
    
    void Plot2DCorrelation(TH2* h1, TCanvas* c1);
};
#endif


CorrelationMaker::CorrelationMaker() {
    m_anaConfig = 0;
    m_debug = false;
    m_rebin = 1;
    m_weightCorrIndex = 0;
    m_runStatCorr = true;

    m_etaBinInterval = 0.1;
    m_etaBoundary = 5.0;


    m_mixDepth = 5;

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



float CorrelationMaker::weight2016hhAna(int index_Nch) {
    // use 1./lumi for this analysis
    float _lumi = 1.;
    if (index_Nch<14) _lumi = 0.082;
    else if (index_Nch<20) _lumi = 0.687;
    else if (index_Nch<24) _lumi = 4.394;
    else _lumi = 55.163;
    return 1./_lumi;
}



float CorrelationMaker::weight2016DzeroAna(int index_Nch) {
    // use 1./lumi for this analysis
    float efficiency = 1.;
    if (index_Nch >= 28) efficiency = 0.85; // hard-coded hatch for D zero analysis systemaitcs for L1_TE200 ineffiency
    return 1./efficiency;
}


float CorrelationMaker::weight295(int index_Nch) {
    // luminosity for triggers used in run 313295 in p+Pb
    // use 1./lumi for this analysis
    float _lumi = 1.;
    if (index_Nch<14) _lumi = 0.004155;
    else if (index_Nch<20) _lumi = 0.060830;
    else if (index_Nch<24) _lumi = 0.376292;
    else _lumi = 2.719786;
    return 1./_lumi;
}


float CorrelationMaker::weight435(int index_Nch) {
    // luminosity for triggers used in run 313295 in p+Pb
    // use 1./lumi for this analysis
    float _lumi = 1.;
    if (index_Nch<20) _lumi = 0.016328;
    else _lumi = 0.262875;
    return 1./_lumi;
}



float CorrelationMaker::weight435AnaSel(int index_Nch) {
    // luminosity for triggers used in run 313295 in p+Pb
    // use 1./lumi for this analysis
    float _lumi = 1.;
    if (index_Nch<14) _lumi = 0.016328;
    else if (index_Nch<20) _lumi = 0.0020476;
    else if (index_Nch<24) _lumi = 0.262875;
    else _lumi = 0.262595;
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

    // Hard-coded accesse to raw 1D multiplicity distributions to support Erorr correction for hh correlation
    // BE CAREFULL HERE
    TH1F* h1_raw_Nch = 0;
    h1_raw_Nch = (TH1F*)gDirectory->Get("h_Nch_selected_raw");
    if (!h1_raw_Nch) cout << "Cannot retrieve 1D multiplicity hitogram for PTY calculation" << endl;

    if (!h2_nsig) cout << "Cannot retrieve yield hitogram for PTY calculation" << endl;
    if (!h2_nsigMix) cout << "Cannot retrieve mix-event yield hitogram for PTY calculation" << endl;
    if (m_debug) {h2_nsig->Print(); h2_nsigMix->Print();}
    if (m_debug) cout << "Merging the following slices" << endl;

    int index_Nch_low = 0;
    int index_Nch_high = 0;
    for (int iNch=1; iNch<m_anaConfig->getInputMultiBinning()->GetXaxis()->GetNbins() + 1; iNch++){
        if (fabs(m_anaConfig->getInputMultiBinning()->GetXaxis()->GetBinLowEdge(iNch) - _Nch_low)  < 0.001)  index_Nch_low  = iNch - 1; 
        if (fabs(m_anaConfig->getInputMultiBinning()->GetXaxis()->GetBinUpEdge(iNch)  - _Nch_high)  < 0.001) index_Nch_high = iNch - 1; 
    }
    if (m_debug) cout << "index_Nch_low = " << index_Nch_low << ",\tindex_Nch_high = " << index_Nch_high << endl;


    int index_pt_low = 0;
    int index_pt_high = 0;
    for (int ipt=1; ipt<m_anaConfig->getInputPtBinning()->GetXaxis()->GetNbins() + 1; ipt++){
        if ( fabs(m_anaConfig->getInputPtBinning()->GetXaxis()->GetBinLowEdge(ipt) - _pt_low)  < 0.001) index_pt_low  = ipt - 1; 
        if ( fabs(m_anaConfig->getInputPtBinning()->GetXaxis()->GetBinUpEdge(ipt)  - _pt_high) < 0.001) index_pt_high = ipt - 1; 
    }

    for (int iNch=index_Nch_low; iNch<index_Nch_high+1; iNch++) {
        for (int ipt=index_pt_low; ipt<index_pt_high+1; ipt++) {


	        TH2* htemp_1 = (TH2*)gDirectory->Get(Form("%sNch%d_pt%d",path_sig.c_str(),iNch,ipt));
	        TH2* htemp_2 = (TH2*)gDirectory->Get(Form("%sNch%d_pt%d",path_mix.c_str(),iNch,ipt));

            // Error correction
            // for double counting trig-asso pairs
            if (m_runStatCorr) {
                // only apply the correction in case the trigger particle collection is a subset of the associate particle collection
                if (  m_anaConfig->getInputPtBinning()->GetXaxis()->GetBinUpEdge(ipt)  <= m_anaConfig->getAssoPtHigh()
                   && m_anaConfig->getInputPtBinning()->GetXaxis()->GetBinLowEdge(ipt) >= m_anaConfig->getAssoPtLow() ) {

                    
                    // nEvent is calculated from raw 1D multiplicity distribution
                    // only consequence would be slightly underestimation <Ntrig> due to too large Nevt including events without particle passing trigger selection
                    // still better than no error correction
                    int _nEvt = 0;
                    if (h1_raw_Nch) {
                        int _Nch_bin_low  = h1_raw_Nch->FindBin(_Nch_low);
                        int _Nch_bin_high = h1_raw_Nch->FindBin(_Nch_high);
                        _nEvt = h1_raw_Nch->Integral(_Nch_bin_low, _Nch_bin_high);
                    }

                    // error correction for double counting based on effective sample size
                    float _Neff_trig = pow(h2_nsig->GetBinContent(ipt+1, iNch+1)/h2_nsig->GetBinError(ipt+1, iNch+1), 2)/_nEvt;
                    float _Neff_pair_same  = htemp_1->GetEffectiveEntries()/_nEvt;
                    float _Ncorr_pair_same = _Neff_pair_same - 0.5*_Neff_trig*(_Neff_trig-1);
                    float _errorScale_same = sqrt(_Neff_pair_same/_Ncorr_pair_same);

                    if (m_debug) {
                        cout << "Running the error correction: " << endl; 
                        cout << "-----------------------------------------" << endl; 
                        cout << "Raw effective size: " << _Neff_pair_same 
                             << ", trigger particle effective size: " << _Neff_trig
                             << ", Corrected effective size: " << _Ncorr_pair_same
                             << ", Error correction scale " << _errorScale_same << endl;
                    } 
                    // mixed event error corrections
                    // not sure how useful they are
                    float _errorScale_mix;
                    if (htemp_2) {
                        float _Neff_pair_mix  = htemp_1->GetEffectiveEntries()/_nEvt;
                        float _Ncorr_pair_mix = _Neff_pair_mix - 0.5*m_mixDepth*_Neff_trig*(_Neff_trig-1);
                        _errorScale_mix = sqrt(_Neff_pair_mix/_Ncorr_pair_mix);
                    }

                    // hard-coded scaling
                    //_errorScale_same = sqrt(2);
                    //_errorScale_mix  = sqrt(2);

                    for (int xbin=1; xbin<htemp_1->GetNbinsX()+1; xbin++) {
                        for (int ybin=1; ybin<htemp_1->GetNbinsY()+1; ybin++) {
                            htemp_1->SetBinError(xbin, ybin, htemp_1->GetBinError(xbin,ybin)*_errorScale_same);
                            if (htemp_2) htemp_2->SetBinError(xbin, ybin, htemp_2->GetBinError(xbin,ybin)*_errorScale_mix );
                        }
                    }
                }
            }

            float _weight = 1.;
            if (getWeightIndex() == 1) {
                _weight = weight2016hhAna(iNch);
            } else if (getWeightIndex() == 2) {
                _weight = weight2016DzeroAna(iNch);
            } else if (getWeightIndex() == 3) {
                _weight = weight295(iNch);
            } else if (getWeightIndex() == 4) {
                _weight = weight435(iNch);
            } else if (getWeightIndex() == 5) {
                _weight = weight435AnaSel(iNch);
            }
            // to be extend to support other analysis
            htemp_1->Scale(_weight);
            if (htemp_2) htemp_2->Scale(_weight);
            nsig    += h2_nsig   ->GetBinContent(ipt+1, iNch+1) * _weight;
            nsigMix += h2_nsigMix->GetBinContent(ipt+1, iNch+1) * _weight;

            // apply stat error correction here
            // getStatCorrection()

	        if (iNch==index_Nch_low && ipt==index_pt_low) {
	            h1 = (TH2*) htemp_1->Clone(Form("h1_Nch%d_%d", index_Nch_low, index_Nch_high));
	            if (htemp_2) h2 = (TH2*) htemp_2->Clone(Form("h2_Nch%d_%d", index_Nch_low, index_Nch_high));
	        } else {
                h1->Add(htemp_1);
                if (htemp_2) h2->Add(htemp_2);
	        }

	        if (m_debug) {
                cout << "pt index = " << ipt << ", Nch index = " << iNch << endl;
                cout << "\tintegrated evt countst, current hist = " << htemp_1->Integral() << endl;
                cout << "\tevt counts for test, current hist = " << htemp_1->GetBinContent(10,10) << endl;
                cout << "\tstat error for test, current hist = " << 100* htemp_1->GetBinError(10,10) / htemp_1->GetBinContent(10,10) << endl;
                cout << "-------------------" << endl;
                cout << "\tintegrated evt countst, merged hist = " << h1->Integral() << endl;
                cout << "\tevt counts for test, merged hist  = " << h1->GetBinContent(10,10) << endl;
                cout << "\tstat error for test, merged hist  = " << 100* h1->GetBinError(10,10) / h1->GetBinContent(10,10) << endl;
            }

	        //if (htemp_1) {delete htemp_1; htemp_1 =0;}
	        //if (htemp_2) {delete htemp_2; htemp_2 =0;}

        }
    }
    if (m_debug) {
        cout << "Merged same  event entries: " << h1->Integral()/1e6 << "M" << endl;
        cout << "Merged mixed event entries: " << h2->Integral()/1e6 << "M" << endl;
    }

    TH1* hphi;

    if (m_h1_sig) {delete m_h1_sig; m_h1_sig = 0;}
    if (m_h1_mix) {delete m_h1_mix; m_h1_mix = 0;}

    // nbins to project
    int nbins_low  = 0.00001 + (m_jetEtaIndexLow-1) - m_etaBinIndexLow + 1;
    int nbins_high = 0.00001 + m_etaBinIndexHigh - (m_jetEtaIndexHigh+1) + 1;
    if (nbins_low != nbins_high) cout << nbins_low << ",\t " << nbins_high << endl;

    if (method == 1) {
        // ATLAS style mixed event normalization
        // make 1D project of signal and mix first
        // then take ratio of two 1D hist

        // Make projections first
        m_h1_sig    = h1->ProjectionY(Form("h_sig_phi_Nch%d",  index_Nch_low), 0.00001 + m_etaBinIndexLow,    0.00001 + m_jetEtaIndexLow-1, "e" );
        m_h1_sig->Add(h1->ProjectionY(Form("h_sig_phi2_Nch%d", index_Nch_low), 0.00001 + m_jetEtaIndexHigh+1, 0.00001 + m_etaBinIndexHigh,  "e" ) );
        m_h1_mix    = h2->ProjectionY(Form("h_mix_phi_Nch%d",  index_Nch_low), 0.00001 + m_etaBinIndexLow,    0.00001 + m_jetEtaIndexLow-1, "e" );
        m_h1_mix->Add(h2->ProjectionY(Form("h_mix_phi2_Nch%d", index_Nch_low), 0.00001 + m_jetEtaIndexHigh+1, 0.00001 + m_etaBinIndexHigh,  "e" ) );

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
        m_h1_sig    = h1->ProjectionY(Form("h_sig_phi_Nch%d",  index_Nch_low), 0.00001 + m_etaBinIndexLow,    0.00001 + m_jetEtaIndexLow-1, "e" );
        m_h1_sig->Add(h1->ProjectionY(Form("h_sig_phi2_Nch%d", index_Nch_low), 0.00001 + m_jetEtaIndexHigh+1, 0.00001 + m_etaBinIndexHigh,  "e" ) );
        m_h1_mix    = h2->ProjectionY(Form("h_mix_phi_Nch%d",  index_Nch_low), 0.00001 + m_etaBinIndexLow,    0.00001 + m_jetEtaIndexLow-1, "e" );
        m_h1_mix->Add(h2->ProjectionY(Form("h_mix_phi2_Nch%d", index_Nch_low), 0.00001 + m_jetEtaIndexHigh+1, 0.00001 + m_etaBinIndexHigh,  "e" ) );
        float int_S = m_h1_sig->Integral();

        h1->Scale(1./h1->Integral());
        h2->Scale(1./h2->Integral());
        TH2F* h_correlation = (TH2F*)h1->Clone();
        h_correlation->Divide(h2);

        hphi    = h_correlation->ProjectionY(Form("h_phi_Nch%d", index_Nch_low), 0.00001 + m_etaBinIndexLow,    0.00001 + m_jetEtaIndexLow-1, "e");
        hphi->Add(h_correlation->ProjectionY(Form("h_phi2_Nch%d",index_Nch_low), 0.00001 + m_jetEtaIndexHigh+1, 0.00001 + m_etaBinIndexHigh, "e"));

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
        m_h1_sig    = h1->ProjectionY(Form("h_sig_phi_Nch%d",  index_Nch_low), 0.00001 + m_etaBinIndexLow,    0.00001 + m_jetEtaIndexLow-1, "e" );
        m_h1_sig->Add(h1->ProjectionY(Form("h_sig_phi2_Nch%d", index_Nch_low), 0.00001 + m_jetEtaIndexHigh+1, 0.00001 + m_etaBinIndexHigh,  "e" ) );
        m_h1_mix    = h2->ProjectionY(Form("h_mix_phi_Nch%d",  index_Nch_low), 0.00001 + m_etaBinIndexLow,    0.00001 + m_jetEtaIndexLow-1, "e" );
        m_h1_mix->Add(h2->ProjectionY(Form("h_mix_phi2_Nch%d", index_Nch_low), 0.00001 + m_jetEtaIndexHigh+1, 0.00001 + m_etaBinIndexHigh,  "e" ) );

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

        hphi    = h_correlation->ProjectionY(Form("h_sig_phi_Nch%d", index_Nch_low), 0.00001 + m_etaBinIndexLow,    0.00001 + m_jetEtaIndexLow-1, "e");
        hphi->Add(h_correlation->ProjectionY(Form("h_sig_phi2_Nch%d",index_Nch_low), 0.00001 + m_jetEtaIndexHigh+1, 0.00001 + m_etaBinIndexHigh, "e"));
        hphi->Scale(1./(nbins_low + nbins_high)); 
        

    } else if (method == 5) {
        // no mixing applied
        // copied from method 1

        // Make projections first
        m_h1_sig    = h1->ProjectionY(Form("h_sig_phi_Nch%d",  index_Nch_low), 0.00001 + m_etaBinIndexLow,    0.00001 + m_jetEtaIndexLow-1, "e" );
        m_h1_sig->Add(h1->ProjectionY(Form("h_sig_phi2_Nch%d", index_Nch_low), 0.00001 + m_jetEtaIndexHigh+1, 0.00001 + m_etaBinIndexHigh,  "e" ) );

        hphi = (TH1*)m_h1_sig->Clone(Form("_hphi"));
        hphi->Scale(1./nsig, "width"); //per trigger & per deltaPhi yield

    } else if (method == 6) {
        // 1D projection with out any normalization
        m_h1_sig    = h1->ProjectionY(Form("h_sig_phi_Nch%d",  index_Nch_low), 0.00001 + m_etaBinIndexLow,    0.00001 + m_jetEtaIndexLow-1, "e" );
        m_h1_sig->Add(h1->ProjectionY(Form("h_sig_phi2_Nch%d", index_Nch_low), 0.00001 + m_jetEtaIndexHigh+1, 0.00001 + m_etaBinIndexHigh,  "e" ) );
        hphi = (TH1*)m_h1_sig->Clone(Form("_hphi"));
        hphi->Scale(1., "width"); //per trigger & per deltaPhi yield

    } else {
        cout << "method is not supported, please choose from 1~3" << endl;
    }

    if (m_debug) cout << "------------------------------------------" << endl;
    delete h1;
    delete h2;

    hphi->Rebin(m_rebin);
    hphi->Scale(1./m_rebin);
    hphi->SetName(Form("h_pty_dphi_gap%.1fto%.1f_Nch%.0fto%.0f_pt%.0fto%.0f",m_anaConfig->getEtaRangeLow(), m_anaConfig->getEtaRangeHigh(), _Nch_low,_Nch_high,_pt_low,_pt_high));

    if (m_h1_sig) {
        m_h1_sig->Rebin(m_rebin);
        m_h1_sig->Scale(1./m_rebin);
        if (m_h1_mix) m_h1_sig->Scale(1./m_h1_mix->Integral());
        m_h1_sig->SetName(Form("h_sig_dphi_gap%.1fto%.1f_Nch%.0fto%.0f_pt%.0fto%.0f",m_anaConfig->getEtaRangeLow(), m_anaConfig->getEtaRangeHigh(), _Nch_low,_Nch_high,_pt_low,_pt_high));
    }
    
    if (m_h1_mix){
        m_h1_mix->Rebin(m_rebin);
        m_h1_mix->Scale(1./m_rebin);
        m_h1_mix->Scale(1./m_h1_mix->Integral());
        m_h1_mix->SetName(Form("h_mix_dphi_gap%.1fto%.1f_Nch%.0fto%.0f_pt%.0fto%.0f",m_anaConfig->getEtaRangeLow(), m_anaConfig->getEtaRangeHigh(), _Nch_low,_Nch_high,_pt_low,_pt_high));
    }

    return hphi;
}



float CorrelationMaker::GetPtMean(float _pt_low, float _pt_high, float _Nch_low, float _Nch_high, int type) {

    TH2D* h1 = 0;

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

    m_anaConfig->inputFile()->cd();

    for (int iNch=index_Nch_low; iNch<index_Nch_high+1; iNch++) {
        for (int ipt=index_pt_low; ipt<index_pt_high+1; ipt++) {

	        TH2F* htemp_1 = (TH2F*)gDirectory->Get(Form("h_pt_mult_trig_same_Mq%d_ptq%d_Sq0_tq0_Pq%d", iNch, ipt, type));

	        if (iNch==index_Nch_low && ipt==index_pt_low) {
	            h1 = (TH2D*) htemp_1->Clone(Form("h1_Nch%d_%d_pt%d_pt%d", index_Nch_low, index_Nch_high, index_pt_low, index_pt_high));
	        } else {
                h1->Add(htemp_1);
	        }
	        delete htemp_1;
        }
    }
    return h1->GetMean(1);
}



float CorrelationMaker::GetMultMean(float _pt_low, float _pt_high, float _Nch_low, float _Nch_high, int type) {

    TH2D* h1 = 0;

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

    m_anaConfig->inputFile()->cd();

    for (int iNch=index_Nch_low; iNch<index_Nch_high+1; iNch++) {
        for (int ipt=index_pt_low; ipt<index_pt_high+1; ipt++) {

	        TH2F* htemp_1 = (TH2F*)gDirectory->Get(Form("h_pt_mult_trig_same_Mq%d_ptq%d_Sq0_tq0_Pq%d", iNch, ipt, type));

	        if (iNch==index_Nch_low && ipt==index_pt_low) {
	            h1 = (TH2D*) htemp_1->Clone(Form("h1_Nch%d_%d_pt%d_pt%d", index_Nch_low, index_Nch_high, index_pt_low, index_pt_high));
	        } else {
                h1->Add(htemp_1);
	        }
	        delete htemp_1;
        }
    }
    return h1->GetMean(2);
}



float CorrelationMaker::GetMultMean(float _pt_low, float _pt_high, float _Nch_low, float _Nch_high, float _add_low, float _add_high, int type) {

    TH2D* h1 = 0;

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

    m_anaConfig->inputFile()->cd();

    for (int iNch=index_Nch_low; iNch<index_Nch_high+1; iNch++) {
        for (int ipt=index_pt_low; ipt<index_pt_high+1; ipt++) {
            for (int iadd=index_add_low; iadd<index_add_high+1; iadd++) {

	            TH2F* htemp_1 = (TH2F*)gDirectory->Get(Form("h_pt_mult_trig_same_Mq%d_ptq%d_Sq%d_tq0_Pq%d", iNch, ipt, iadd, type));

	            if (iNch==index_Nch_low && ipt==index_pt_low) {
	                h1 = (TH2D*) htemp_1->Clone(Form("h1_Nch%d_%d_pt%d_pt%d", index_Nch_low, index_Nch_high, index_pt_low, index_pt_high));
	            } else {
                    h1->Add(htemp_1);
	            }
	            delete htemp_1;
            }
        }
    }
    return h1->GetMean(2);
}




TH1* CorrelationMaker::MakeCorrUpc(float _pt_low, float _pt_high, float _Nch_low, float _Nch_high, int method, int type) {

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
    TH2F* h2_nevt = (TH2F*)gDirectory->Get("h_NeventTrig");

    if (!h2_nsig) cout << "Cannot retrieve yield hitogram for PTY calculation" << endl;
    if (!h2_nevt) cout << "Cannot retrieve Nevt hitogram for PTY calculation" << endl;
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
            //cout << iNch << "\t" << ipt << endl;

            // example histogram names from Blair's file
            //h_deta_dphi_same_Mq4_ptq2_Pq1	
            //h_deta_dphi_same_Mq1_ptq4_Pq0
            //h2_deta_dphi_mixed_Mq4_ptq2_Sq0_tq0_Pq1;1 
	        TH2D* htemp_1 = (TH2D*)gDirectory->Get(Form("%sMq%d_ptq%d_Pq%d",path_sig.c_str(),iNch, ipt, type));
	        TH2D* htemp_2 = (TH2D*)gDirectory->Get(Form("%sMq%d_ptq%d_Sq0_tq0_Pq%d",path_mix.c_str(),iNch, ipt, type));

            // event counts
            // example histogram names from Blair's file
            // h_pt_mult_trig_same_Mq8_ptq7_Sq0_tq0_Pq1
	        TH2F* htempt_Nevt = (TH2F*)gDirectory->Get(Form("h_pt_mult_trig_same_Mq%d_ptq%d_Sq0_tq0_Pq%d", iNch, ipt, type));

            // Error correction
            // for double counting trig-asso pairs
            if (m_runStatCorr) {
                // only apply error correction when trigger and associate particle selection overlaps
                if (  m_anaConfig->getInputPtBinning()->GetXaxis()->GetBinUpEdge(ipt)  <= m_anaConfig->getAssoPtHigh()
                   && m_anaConfig->getInputPtBinning()->GetXaxis()->GetBinLowEdge(ipt) >= m_anaConfig->getAssoPtLow() ) {

                    float _nEvt = h2_nevt->GetBinContent(ipt+1, iNch+1);

                    // error correction for double counting based on effective sample size
                    float _Neff_trig = pow(h2_nsig->GetBinContent(ipt+1, iNch+1)/h2_nsig->GetBinError(ipt+1, iNch+1), 2)/_nEvt;
                    float _Neff_pair_same  = htemp_1->GetEffectiveEntries()/_nEvt;
                    float _Ncorr_pair_same = _Neff_pair_same - 0.5*_Neff_trig*(_Neff_trig-1);
                    float _errorScale_same = sqrt(_Neff_pair_same/_Ncorr_pair_same);
                    if (m_debug) {
                        cout << "Running the error correction: " << endl; 
                        cout << "Raw effective size: " << _Neff_pair_same 
                             << ", trigger particle effective size: " << _Neff_trig
                             << ", Corrected effective size: " << _Ncorr_pair_same
                             << ", Error correction scale " << _errorScale_same << endl;
                    } 
                    // mixed event error corrections
                    // not sure how useful they are
                    float _Neff_pair_mix  = htemp_1->GetEffectiveEntries()/_nEvt;
                    float _Ncorr_pair_mix = _Neff_pair_mix - 0.5*m_mixDepth*_Neff_trig*(_Neff_trig-1);
                    float _errorScale_mix = sqrt(_Neff_pair_mix/_Ncorr_pair_mix);

                    for (int xbin=1; xbin<htemp_1->GetNbinsX()+1; xbin++) {
                        for (int ybin=1; ybin<htemp_1->GetNbinsY()+1; ybin++) {
                            htemp_1->SetBinError(xbin, ybin, htemp_1->GetBinError(xbin,ybin)*_errorScale_same);
                            htemp_2->SetBinError(xbin, ybin, htemp_2->GetBinError(xbin,ybin)*_errorScale_mix );
                        }
                    }
                }
            }

            nsig += h2_nsig->GetBinContent(ipt+1, iNch+1);
            //cout << Form("%sMq%d_ptq%d_Pq%d",path_sig.c_str(),iNch, ipt, type) << endl;
            //cout << "\tnsig = " <<  h2_nsig->GetBinContent(ipt+1, iNch+1) << endl;
            //htemp_1->Print();
            //cout << "-------------" << endl;

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
            delete htempt_Nevt;
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
    int nbins_low  = 0.00001 + (m_jetEtaIndexLow-1) - m_etaBinIndexLow + 1;
    int nbins_high = 0.00001 + m_etaBinIndexHigh - (m_jetEtaIndexHigh+1) + 1;

    if (nbins_low != nbins_high) cout << nbins_low << ",\t " << nbins_high << endl;
    if (m_debug) {
        cout << endl;
        cout << "Debug the gap index determination: " << endl;

        cout << m_etaBinIndexLow << ",\t" << m_jetEtaIndexLow -1 
             << " <==> " << h1->GetXaxis()->GetBinLowEdge(0.00001 + m_etaBinIndexLow) << ",\t" << h1->GetXaxis()->GetBinUpEdge(0.00001 + m_jetEtaIndexLow -1) << endl;
        cout << m_jetEtaIndexHigh+1 << ",\t" << m_etaBinIndexHigh
             << " <==> " << h1->GetXaxis()->GetBinLowEdge(0.00001 + m_jetEtaIndexHigh+1) << ",\t" << h1->GetXaxis()->GetBinUpEdge(0.00001 + m_etaBinIndexHigh) << endl;
        cout << "\tMerged Nbins negative side =  " << nbins_low << ", positive side =  " << nbins_high << endl;
        cout << endl;
    }

    if (method == 1) {
        // ATLAS style mixed event normalization
        // make 1D project of signal and mix first
        // then take ratio of two 1D hist

        // Make projections first
        // since we are converting float index to int index
        // +0.00001 to avoid fluctuations in C++
        m_h1_sig    = h1->ProjectionY(Form("h_sig_phi_Nch%d",  index_Nch_low), 0.00001 + m_etaBinIndexLow,    0.00001 + m_jetEtaIndexLow-1, "e" );
        m_h1_sig->Add(h1->ProjectionY(Form("h_sig_phi2_Nch%d", index_Nch_low), 0.00001 + m_jetEtaIndexHigh+1, 0.00001 + m_etaBinIndexHigh,  "e" ) );
        m_h1_mix    = h2->ProjectionY(Form("h_mix_phi_Nch%d",  index_Nch_low), 0.00001 + m_etaBinIndexLow,    0.00001 + m_jetEtaIndexLow-1, "e" );
        m_h1_mix->Add(h2->ProjectionY(Form("h_mix_phi2_Nch%d", index_Nch_low), 0.00001 + m_jetEtaIndexHigh+1, 0.00001 + m_etaBinIndexHigh,  "e" ) );

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
        m_h1_sig    = h1->ProjectionY(Form("h_sig_phi_Nch%d",  index_Nch_low), 0.00001 + m_etaBinIndexLow,    0.00001 + m_jetEtaIndexLow-1, "e" );
        m_h1_sig->Add(h1->ProjectionY(Form("h_sig_phi2_Nch%d", index_Nch_low), 0.00001 + m_jetEtaIndexHigh+1, 0.00001 + m_etaBinIndexHigh,  "e" ) );
        m_h1_mix    = h2->ProjectionY(Form("h_mix_phi_Nch%d",  index_Nch_low), 0.00001 + m_etaBinIndexLow,    0.00001 + m_jetEtaIndexLow-1, "e" );
        m_h1_mix->Add(h2->ProjectionY(Form("h_mix_phi2_Nch%d", index_Nch_low), 0.00001 + m_jetEtaIndexHigh+1, 0.00001 + m_etaBinIndexHigh,  "e" ) );
        float int_S = m_h1_sig->Integral();

        h1->Scale(1./h1->Integral());
        h2->Scale(1./h2->Integral());
        TH2D* h_correlation = (TH2D*)h1->Clone();
        h_correlation->Divide(h2);

        hphi    = h_correlation->ProjectionY(Form("h_phi_Nch%d", index_Nch_low), 0.00001 + m_etaBinIndexLow,    0.00001 + m_jetEtaIndexLow-1, "e");
        hphi->Add(h_correlation->ProjectionY(Form("h_phi2_Nch%d",index_Nch_low), 0.00001 + m_jetEtaIndexHigh+1, 0.00001 + m_etaBinIndexHigh, "e"));

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
        m_h1_sig    = h1->ProjectionY(Form("h_sig_phi_Nch%d",  index_Nch_low), 0.00001 + m_etaBinIndexLow,    0.00001 + m_jetEtaIndexLow-1, "e" );
        m_h1_sig->Add(h1->ProjectionY(Form("h_sig_phi2_Nch%d", index_Nch_low), 0.00001 + m_jetEtaIndexHigh+1, 0.00001 + m_etaBinIndexHigh,  "e" ) );
        m_h1_mix    = h2->ProjectionY(Form("h_mix_phi_Nch%d",  index_Nch_low), 0.00001 + m_etaBinIndexLow,    0.00001 + m_jetEtaIndexLow-1, "e" );
        m_h1_mix->Add(h2->ProjectionY(Form("h_mix_phi2_Nch%d", index_Nch_low), 0.00001 + m_jetEtaIndexHigh+1, 0.00001 + m_etaBinIndexHigh,  "e" ) );

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

        hphi    = h_correlation->ProjectionY(Form("h_sig_phi_Nch%d", index_Nch_low), 0.00001 + m_etaBinIndexLow,    0.00001 + m_jetEtaIndexLow-1, "e");
        hphi->Add(h_correlation->ProjectionY(Form("h_sig_phi2_Nch%d",index_Nch_low), 0.00001 + m_jetEtaIndexHigh+1, 0.00001 + m_etaBinIndexHigh, "e"));
        hphi->Scale(1./(nbins_low + nbins_high)); 
        
    } else if (method == 5) {
        // no mixing applied
        // copied from method 1

        // Make projections first
        
        if (m_debug) { 
            cout << "Integral before projection: " << h1->Integral() << endl;
        }
        m_h1_sig    = h1->ProjectionY(Form("h_sig_phi_Nch%d",  index_Nch_low), 0.00001 + m_etaBinIndexLow,    0.00001 + m_jetEtaIndexLow-1, "e" );
        m_h1_sig->Add(h1->ProjectionY(Form("h_sig_phi2_Nch%d", index_Nch_low), 0.00001 + m_jetEtaIndexHigh+1, 0.00001 + m_etaBinIndexHigh,  "e" ) );
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
    hphi->SetName(Form("h_pty_dphi_gap%.1fto%.1f_Nch%.0fto%.0f_pt%.1fto%.1f",m_anaConfig->getEtaRangeLow(), m_anaConfig->getEtaRangeHigh(), _Nch_low,_Nch_high,_pt_low,_pt_high));

    if (m_h1_sig) {
        m_h1_sig->Rebin(m_rebin);
        //m_h1_sig->Scale(1./m_rebin);
        m_h1_sig->Scale(1./m_h1_sig->Integral());
        m_h1_sig->SetName(Form("h_sig_dphi_gap%.1fto%.1f_Nch%.0fto%.0f_pt%.1fto%.1f",m_anaConfig->getEtaRangeLow(), m_anaConfig->getEtaRangeHigh(), _Nch_low,_Nch_high,_pt_low,_pt_high));
    }
    
    if (m_h1_mix){
        m_h1_mix->Rebin(m_rebin);
        //m_h1_mix->Scale(1./m_rebin);
        m_h1_mix->Scale(1./m_h1_mix->Integral());
        m_h1_mix->SetName(Form("h_mix_dphi_gap%.1fto%.1f_Nch%.0fto%.0f_pt%.1fto%.1f",m_anaConfig->getEtaRangeLow(), m_anaConfig->getEtaRangeHigh(), _Nch_low,_Nch_high,_pt_low,_pt_high));
    }

    return hphi;
}



TH1* CorrelationMaker::MakeCorrUpc(float _pt_low, float _pt_high, float _Nch_low, float _Nch_high, float _add_low, float _add_high, int method, int type) {

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


    int index_add_low = 0;
    int index_add_high = 0;
    for (int iadd=1; iadd<m_anaConfig->getInputThirdBinning()->GetXaxis()->GetNbins() + 1; iadd++){
        if (m_anaConfig->getInputThirdBinning()->GetXaxis()->GetBinLowEdge(iadd) == _add_low)  index_add_low  = iadd - 1; 
        if (m_anaConfig->getInputThirdBinning()->GetXaxis()->GetBinUpEdge(iadd)  == _add_high) index_add_high = iadd - 1; 
    }


    m_anaConfig->inputFile()->cd();
    TH3* h2_nsig = (TH3*)gDirectory->Get(m_anaConfig->getTrigYieldHistName().c_str());
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
            for (int iadd=index_add_low; iadd<index_add_high+1; iadd++) {

                // example names
                //h_deta_dphi_same_Mq8_ptq7_Sq15_Pq1
                //h2_deta_dphi_mixed_Mq8_ptq7_Sq15_tq0_Pq1
	            TH2D* htemp_1 = (TH2D*)gDirectory->Get(Form("%sMq%d_ptq%d_Sq%d_Pq%d", path_sig.c_str(), iNch, ipt, iadd, type));
	            TH2D* htemp_2 = (TH2D*)gDirectory->Get(Form("%sMq%d_ptq%d_Sq%d_tq0_Pq%d",path_mix.c_str(),iNch, ipt, iadd, type));

                if (m_debug) {
                    cout << "pt: " << h2_nsig->GetXaxis()->GetBinLowEdge(ipt+1) << " ~ " << h2_nsig->GetXaxis()->GetBinUpEdge(ipt+1) << endl;
                    cout << "nch: " << h2_nsig->GetYaxis()->GetBinLowEdge(iNch+1) << " ~ " << h2_nsig->GetYaxis()->GetBinUpEdge(iNch+1) << endl;
                    cout << "gap: " << h2_nsig->GetZaxis()->GetBinLowEdge(iadd+1) << " ~ " << h2_nsig->GetZaxis()->GetBinUpEdge(iadd+1) << endl;
	                cout << "pt index = " << ipt << ", Nch index = " << iNch << ", gap index = " << iadd << endl;
	                cout << "ntrig = " << h2_nsig->GetBinContent(ipt+1, iNch+1, iadd+1) << endl;
                    if (m_debug) cout << "same event entries:   " << htemp_1->GetEntries() << endl;
                    if (m_debug) cout << "same event integral:  " << htemp_1->Integral()   << endl;
                    if (m_debug) cout << "mixed event entries:  " << htemp_2->GetEntries() << endl;
                    if (m_debug) cout << "mixed event integral: " << htemp_2->Integral()   << endl;
                    cout << endl;
                }

                if (htemp_1->GetEntries() == 0) {
	                delete htemp_1;
	                delete htemp_2;
                    continue;
                }

                nsig += h2_nsig->GetBinContent(ipt+1, iNch+1, iadd+1);

	            if (iNch==index_Nch_low && ipt==index_pt_low && iadd==index_add_low) {
	                h1 = (TH2D*) htemp_1->Clone(Form("h1_Nch%d_%d_pt%d_pt%d", index_Nch_low, index_Nch_high, index_pt_low, index_pt_high));
	                h2 = (TH2D*) htemp_2->Clone(Form("h2_Nch%d_%d_pt%d_pt%d", index_Nch_low, index_Nch_high, index_pt_low, index_pt_high));
	            } else {
                    h1->Add(htemp_1);
                    h2->Add(htemp_2);
	            }

	            delete htemp_1;
	            delete htemp_2;
            }
        }
    }

    // Error correction
    // for double counting trig-asso pairs
    // here we simple assume the symmetric trig-asso selection will always be used
    // so the erorr correction factor is just sqrt(2) to histogram after merging
    // otherwise, the 2D-like error correction should be applied to each histogram before merging
    if (m_runStatCorr) {
        for (int xbin=1; xbin<h1->GetNbinsX()+1; xbin++) {
            for (int ybin=1; ybin<h1->GetNbinsY()+1; ybin++) {
                h1->SetBinError(xbin, ybin, h1->GetBinError(xbin,ybin)*sqrt(2));
                h2->SetBinError(xbin, ybin, h2->GetBinError(xbin,ybin)*sqrt(2));
            }
        }
    }

    TH1* hphi;
    if (m_debug) cout << "Merged same event entries: " << h1->GetEntries()/1e6 << "M" << endl;
    if (m_debug) cout << "Merged same event integral: " << h1->Integral()/1e6 << "M" << endl;
    if (m_debug) cout << "Merged mixed event entries: " << h2->GetEntries()/1e6 << "M" << endl;
    if (m_debug) cout << "Merged mixed event integral: " << h2->Integral()/1e6 << "M" << endl;

    if (m_h1_sig) {delete m_h1_sig; m_h1_sig = 0;}
    if (m_h1_mix) {delete m_h1_mix; m_h1_mix = 0;}

    // nbins to project
    int nbins_low  = 0.00001 + (m_jetEtaIndexLow-1) - m_etaBinIndexLow + 1;
    int nbins_high = 0.00001 + m_etaBinIndexHigh - (m_jetEtaIndexHigh+1) + 1;

    if (nbins_low != nbins_high) cout << nbins_low << ",\t " << nbins_high << endl;
    if (m_debug) {
        cout << endl;
        cout << "Debug the gap index determination: " << endl;

        cout << m_etaBinIndexLow << ",\t" << m_jetEtaIndexLow -1 
             << " <==> " << h1->GetXaxis()->GetBinLowEdge(0.00001 + m_etaBinIndexLow) << ",\t" << h1->GetXaxis()->GetBinUpEdge(0.00001 + m_jetEtaIndexLow -1) << endl;
        cout << m_jetEtaIndexHigh+1 << ",\t" << m_etaBinIndexHigh
             << " <==> " << h1->GetXaxis()->GetBinLowEdge(0.00001 + m_jetEtaIndexHigh+1) << ",\t" << h1->GetXaxis()->GetBinUpEdge(0.00001 + m_etaBinIndexHigh) << endl;
        cout << "\tMerged Nbins negative side =  " << nbins_low << ", positive side =  " << nbins_high << endl;
        cout << endl;
    }

    if (method == 1) {
        // ATLAS style mixed event normalization
        // make 1D project of signal and mix first
        // then take ratio of two 1D hist

        // Make projections first
        // since we are converting float index to int index
        // +0.00001 to avoid fluctuations in C++
        m_h1_sig    = h1->ProjectionY(Form("h_sig_phi_Nch%d",  index_Nch_low), 0.00001 + m_etaBinIndexLow,    0.00001 + m_jetEtaIndexLow-1, "e" );
        m_h1_sig->Add(h1->ProjectionY(Form("h_sig_phi2_Nch%d", index_Nch_low), 0.00001 + m_jetEtaIndexHigh+1, 0.00001 + m_etaBinIndexHigh,  "e" ) );
        m_h1_mix    = h2->ProjectionY(Form("h_mix_phi_Nch%d",  index_Nch_low), 0.00001 + m_etaBinIndexLow,    0.00001 + m_jetEtaIndexLow-1, "e" );
        m_h1_mix->Add(h2->ProjectionY(Form("h_mix_phi2_Nch%d", index_Nch_low), 0.00001 + m_jetEtaIndexHigh+1, 0.00001 + m_etaBinIndexHigh,  "e" ) );

        float int_S = m_h1_sig->Integral();

        hphi = (TH1*)m_h1_sig->Clone(Form("_hphi"));
        hphi->Divide(m_h1_mix);

        float int_C = hphi->Integral();
        float K = int_S / int_C;
        if (m_debug) {
            cout << "same integral = " << m_h1_sig->Integral() << endl;
            cout << "mix  integral = " << m_h1_mix->Integral() << endl;
            cout << "nsig = " << nsig << endl;
            cout << "int_S = " << int_S << endl;
            cout << "int_C = " << int_C << endl;
        }
        hphi->Scale(K/nsig, "width"); //per trigger & per deltaPhi yield

    } else if (method == 2) {
        // ATLAS style mixed event normalization
        // take 2D ratio, then project to 1D
        
        // start with N_{pair}
        m_h1_sig    = h1->ProjectionY(Form("h_sig_phi_Nch%d",  index_Nch_low), 0.00001 + m_etaBinIndexLow,    0.00001 + m_jetEtaIndexLow-1, "e" );
        m_h1_sig->Add(h1->ProjectionY(Form("h_sig_phi2_Nch%d", index_Nch_low), 0.00001 + m_jetEtaIndexHigh+1, 0.00001 + m_etaBinIndexHigh,  "e" ) );
        m_h1_mix    = h2->ProjectionY(Form("h_mix_phi_Nch%d",  index_Nch_low), 0.00001 + m_etaBinIndexLow,    0.00001 + m_jetEtaIndexLow-1, "e" );
        m_h1_mix->Add(h2->ProjectionY(Form("h_mix_phi2_Nch%d", index_Nch_low), 0.00001 + m_jetEtaIndexHigh+1, 0.00001 + m_etaBinIndexHigh,  "e" ) );
        float int_S = m_h1_sig->Integral();

        h1->Scale(1./h1->Integral());
        h2->Scale(1./h2->Integral());
        TH2D* h_correlation = (TH2D*)h1->Clone();
        h_correlation->Divide(h2);

        hphi    = h_correlation->ProjectionY(Form("h_phi_Nch%d", index_Nch_low), 0.00001 + m_etaBinIndexLow,    0.00001 + m_jetEtaIndexLow-1, "e");
        hphi->Add(h_correlation->ProjectionY(Form("h_phi2_Nch%d",index_Nch_low), 0.00001 + m_jetEtaIndexHigh+1, 0.00001 + m_etaBinIndexHigh, "e"));

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
        m_h1_sig    = h1->ProjectionY(Form("h_sig_phi_Nch%d",  index_Nch_low), 0.00001 + m_etaBinIndexLow,    0.00001 + m_jetEtaIndexLow-1, "e" );
        m_h1_sig->Add(h1->ProjectionY(Form("h_sig_phi2_Nch%d", index_Nch_low), 0.00001 + m_jetEtaIndexHigh+1, 0.00001 + m_etaBinIndexHigh,  "e" ) );
        m_h1_mix    = h2->ProjectionY(Form("h_mix_phi_Nch%d",  index_Nch_low), 0.00001 + m_etaBinIndexLow,    0.00001 + m_jetEtaIndexLow-1, "e" );
        m_h1_mix->Add(h2->ProjectionY(Form("h_mix_phi2_Nch%d", index_Nch_low), 0.00001 + m_jetEtaIndexHigh+1, 0.00001 + m_etaBinIndexHigh,  "e" ) );

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

        hphi    = h_correlation->ProjectionY(Form("h_sig_phi_Nch%d", index_Nch_low), 0.00001 + m_etaBinIndexLow,    0.00001 + m_jetEtaIndexLow-1, "e");
        hphi->Add(h_correlation->ProjectionY(Form("h_sig_phi2_Nch%d",index_Nch_low), 0.00001 + m_jetEtaIndexHigh+1, 0.00001 + m_etaBinIndexHigh, "e"));
        hphi->Scale(1./(nbins_low + nbins_high)); 
        
    } else if (method == 5) {
        // no mixing applied
        // copied from method 1

        // Make projections first
        
        if (m_debug) { 
            cout << "Integral before projection: " << h1->Integral() << endl;
        }
        m_h1_sig    = h1->ProjectionY(Form("h_sig_phi_Nch%d",  index_Nch_low), 0.00001 + m_etaBinIndexLow,    0.00001 + m_jetEtaIndexLow-1, "e" );
        m_h1_sig->Add(h1->ProjectionY(Form("h_sig_phi2_Nch%d", index_Nch_low), 0.00001 + m_jetEtaIndexHigh+1, 0.00001 + m_etaBinIndexHigh,  "e" ) );
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
    hphi->SetName(Form("h_pty_dphi_gap%.1fto%.1f_Nch%.0fto%.0f_pt%.1fto%.1f",m_anaConfig->getEtaRangeLow(), m_anaConfig->getEtaRangeHigh(), _Nch_low,_Nch_high,_pt_low,_pt_high));

    if (m_h1_sig) {
        m_h1_sig->Rebin(m_rebin);
        //m_h1_sig->Scale(1./m_rebin);
        m_h1_sig->Scale(1./m_h1_sig->Integral());
        m_h1_sig->SetName(Form("h_sig_dphi_gap%.1fto%.1f_Nch%.0fto%.0f_pt%.1fto%.1f",m_anaConfig->getEtaRangeLow(), m_anaConfig->getEtaRangeHigh(), _Nch_low,_Nch_high,_pt_low,_pt_high));
    }
    
    if (m_h1_mix){
        m_h1_mix->Rebin(m_rebin);
        //m_h1_mix->Scale(1./m_rebin);
        m_h1_mix->Scale(1./m_h1_mix->Integral());
        m_h1_mix->SetName(Form("h_mix_dphi_gap%.1fto%.1f_Nch%.0fto%.0f_pt%.1fto%.1f",m_anaConfig->getEtaRangeLow(), m_anaConfig->getEtaRangeHigh(), _Nch_low,_Nch_high,_pt_low,_pt_high));
    }

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

            //h2_deta_dphi_mixed_Mq7_ptq4_Sq0_tq0_Pq1
	        TH2* htemp_1 = (TH2*)gDirectory->Get(Form("%sMq%d_ptq%d_Pq1",path_sig.c_str(),iNch,ipt));
	        TH2* htemp_2 = (TH2*)gDirectory->Get(Form("%sMq%d_ptq%d_Sq0_tq0_Pq1",path_mix.c_str(),iNch,ipt));


            htemp_2->Scale(htemp_1->Integral()/htemp_2->Integral());

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
        h_correlation->Scale(K/nsig, "width"); //per trigger & per deltaPhi yield 

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



// Making 2D correction in pt, eta and gap
TH2* CorrelationMaker::Make2DCorrUpc(float _pt_low, float _pt_high, float _Nch_low, float _Nch_high, float _add_low, float _add_high, int method) {

    TH2* h1 = 0;
    TH2* h2 = 0;
    
    int type = 1; // UPC events

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

    int index_add_low = 0;
    int index_add_high = 0;
    for (int iadd=1; iadd<m_anaConfig->getInputThirdBinning()->GetXaxis()->GetNbins() + 1; iadd++){
        if (m_anaConfig->getInputThirdBinning()->GetXaxis()->GetBinLowEdge(iadd) == _add_low)  index_add_low  = iadd - 1;
        if (m_anaConfig->getInputThirdBinning()->GetXaxis()->GetBinUpEdge(iadd)  == _add_high) index_add_high = iadd - 1;
    }

    m_anaConfig->inputFile()->cd();
    TH2F* h2_nsig = (TH2F*)gDirectory->Get(m_anaConfig->getTrigYieldHistName().c_str());
    if (!h2_nsig) cout << "Cannot retrieve yield hitogram for PTY calculation" << endl;

    float nsig = 0;
    for (int iNch=index_Nch_low; iNch<index_Nch_high+1; iNch++) {
        for (int ipt=index_pt_low; ipt<index_pt_high+1; ipt++) {
            for (int iadd=index_add_low; iadd<index_add_high+1; iadd++) {

                // example names
                //h_deta_dphi_same_Mq8_ptq7_Sq15_Pq1
                //h2_deta_dphi_mixed_Mq8_ptq7_Sq15_tq0_Pq1
                TH2D* htemp_1 = (TH2D*)gDirectory->Get(Form("%sMq%d_ptq%d_Sq%d_Pq%d", path_sig.c_str(), iNch, ipt, iadd, type));
                TH2D* htemp_2 = (TH2D*)gDirectory->Get(Form("%sMq%d_ptq%d_Sq%d_tq0_Pq%d",path_mix.c_str(),iNch, ipt, iadd, type));

                if (htemp_1->GetEntries() == 0) {
	                delete htemp_1;
	                delete htemp_2;
                    continue;
                }

                htemp_2->Scale(htemp_1->Integral()/htemp_2->Integral());
                nsig += h2_nsig->GetBinContent(ipt+1, iNch+1, iadd+1);

                if (iNch==index_Nch_low && ipt==index_pt_low && iadd==index_add_low) {
                    h1 = (TH2D*) htemp_1->Clone(Form("h1_Nch%d_%d_pt%d_pt%d", index_Nch_low, index_Nch_high, index_pt_low, index_pt_high));
                    h2 = (TH2D*) htemp_2->Clone(Form("h2_Nch%d_%d_pt%d_pt%d", index_Nch_low, index_Nch_high, index_pt_low, index_pt_high));
                } else {
                    h1->Add(htemp_1);
                    h2->Add(htemp_2);
                }    

                delete htemp_1;
                delete htemp_2;
                if (m_debug) {
                    cout << "pt: " << h2_nsig->GetXaxis()->GetBinLowEdge(ipt+1) << " ~ " << h2_nsig->GetXaxis()->GetBinUpEdge(ipt+1) << endl;
                    cout << "nch: " << h2_nsig->GetYaxis()->GetBinLowEdge(iNch+1) << " ~ " << h2_nsig->GetYaxis()->GetBinUpEdge(iNch+1) << endl;
                    cout << "gap: " << h2_nsig->GetZaxis()->GetBinLowEdge(iadd+1) << " ~ " << h2_nsig->GetZaxis()->GetBinUpEdge(iadd+1) << endl;
                    cout << "pt index = " << ipt << ", Nch index = " << iNch << ", gap index = " << iadd << endl;
                    cout << "ntrig = " << h2_nsig->GetBinContent(ipt+1, iNch+1, iadd+1) << endl;
                    cout << endl;
                }    
            } 

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
        h_correlation->Scale(K/nsig, "width"); //per trigger & per deltaPhi yield 

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

                // error correction for double counting based on effective sample size
                float _Neff_trig = pow(h2_nsig->GetBinContent(ipt+1, iNch+1, iadd+1)/h2_nsig->GetBinError(ipt+1, iNch+1, iadd+1), 2); 
                float _Neff_pair_same  = htemp_1->GetEffectiveEntries();
                float _Ncorr_pair_same = _Neff_pair_same - 0.5*_Neff_trig*(_Neff_trig-1);
                float _errorScale_same = sqrt(_Neff_pair_same/_Ncorr_pair_same);
                // mixed event error correction
                // not sure how useful they are
                float _Neff_pair_mix  = htemp_1->GetEffectiveEntries();
                float _Ncorr_pair_mix = _Neff_pair_mix - 0.5*m_mixDepth*_Neff_trig*(_Neff_trig-1);
                float _errorScale_mix = sqrt(_Neff_pair_mix/_Ncorr_pair_mix);

                // Error correction
                // no need for D-h, psi-h correlation study
                // not implemented for here

                float _weight = 1.;
                if (getWeightIndex() == 1) {
                    _weight = weight2016hhAna(iNch);
                } else if (getWeightIndex() == 2) {
                    _weight = weight2016DzeroAna(iNch);
                } else if (getWeightIndex() == 3) {
                    _weight = weight295(iNch);
                } else if (getWeightIndex() == 4) {
                    _weight = weight435(iNch);
                } else if (getWeightIndex() == 5) {
                    _weight = weight435AnaSel(iNch);
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
    int nbins_low  = 0.00001 + (m_jetEtaIndexLow-1) - m_etaBinIndexLow + 1;
    int nbins_high = 0.00001 + m_etaBinIndexHigh - (m_jetEtaIndexHigh+1) + 1;
    if (nbins_low != nbins_high) cout << nbins_low << ",\t " << nbins_high << endl;

    if (method == 1) {
        // ATLAS style mixed event normalization
        // make 1D project of signal and mix first
        // then take ratio of two 1D hist

        // Make projections first
        m_h1_sig    = h1->ProjectionY(Form("h_sig_phi_Nch%d",  index_Nch_low), 0.00001 + m_etaBinIndexLow,    0.00001 + m_jetEtaIndexLow-1, "e" );
        m_h1_sig->Add(h1->ProjectionY(Form("h_sig_phi2_Nch%d", index_Nch_low), 0.00001 + m_jetEtaIndexHigh+1, 0.00001 + m_etaBinIndexHigh,  "e" ) );
        m_h1_mix    = h2->ProjectionY(Form("h_mix_phi_Nch%d",  index_Nch_low), 0.00001 + m_etaBinIndexLow,    0.00001 + m_jetEtaIndexLow-1, "e" );
        m_h1_mix->Add(h2->ProjectionY(Form("h_mix_phi2_Nch%d", index_Nch_low), 0.00001 + m_jetEtaIndexHigh+1, 0.00001 + m_etaBinIndexHigh,  "e" ) );

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
        m_h1_sig    = h1->ProjectionY(Form("h_sig_phi_Nch%d",  index_Nch_low), 0.00001 + m_etaBinIndexLow,    0.00001 + m_jetEtaIndexLow-1, "e" );
        m_h1_sig->Add(h1->ProjectionY(Form("h_sig_phi2_Nch%d", index_Nch_low), 0.00001 + m_jetEtaIndexHigh+1, 0.00001 + m_etaBinIndexHigh,  "e" ) );
        m_h1_mix    = h2->ProjectionY(Form("h_mix_phi_Nch%d",  index_Nch_low), 0.00001 + m_etaBinIndexLow,    0.00001 + m_jetEtaIndexLow-1, "e" );
        m_h1_mix->Add(h2->ProjectionY(Form("h_mix_phi2_Nch%d", index_Nch_low), 0.00001 + m_jetEtaIndexHigh+1, 0.00001 + m_etaBinIndexHigh,  "e" ) );
        float int_S = m_h1_sig->Integral();

        h1->Scale(1./h1->Integral());
        h2->Scale(1./h2->Integral());
        TH2F* h_correlation = (TH2F*)h1->Clone();
        h_correlation->Divide(h2);

        hphi    = h_correlation->ProjectionY(Form("h_phi_Nch%d", index_Nch_low), 0.00001 + m_etaBinIndexLow,    0.00001 + m_jetEtaIndexLow-1, "e");
        hphi->Add(h_correlation->ProjectionY(Form("h_phi2_Nch%d",index_Nch_low), 0.00001 + m_jetEtaIndexHigh+1, 0.00001 + m_etaBinIndexHigh, "e"));

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
        m_h1_sig    = h1->ProjectionY(Form("h_sig_phi_Nch%d",  index_Nch_low), 0.00001 + m_etaBinIndexLow,    0.00001 + m_jetEtaIndexLow-1, "e" );
        m_h1_sig->Add(h1->ProjectionY(Form("h_sig_phi2_Nch%d", index_Nch_low), 0.00001 + m_jetEtaIndexHigh+1, 0.00001 + m_etaBinIndexHigh,  "e" ) );
        m_h1_mix    = h2->ProjectionY(Form("h_mix_phi_Nch%d",  index_Nch_low), 0.00001 + m_etaBinIndexLow,    0.00001 + m_jetEtaIndexLow-1, "e" );
        m_h1_mix->Add(h2->ProjectionY(Form("h_mix_phi2_Nch%d", index_Nch_low), 0.00001 + m_jetEtaIndexHigh+1, 0.00001 + m_etaBinIndexHigh,  "e" ) );

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

        hphi    = h_correlation->ProjectionY(Form("h_sig_phi_Nch%d", index_Nch_low), 0.00001 + m_etaBinIndexLow,    0.00001 + m_jetEtaIndexLow-1, "e");
        hphi->Add(h_correlation->ProjectionY(Form("h_sig_phi2_Nch%d",index_Nch_low), 0.00001 + m_jetEtaIndexHigh+1, 0.00001 + m_etaBinIndexHigh, "e"));
        hphi->Scale(1./(nbins_low + nbins_high)); 
        

    } else if (method == 5) {
        // no mixing applied
        // copied from method 1

        // Make projections first
        m_h1_sig    = h1->ProjectionY(Form("h_sig_phi_Nch%d",  index_Nch_low), 0.00001 + m_etaBinIndexLow,    0.00001 + m_jetEtaIndexLow-1, "e" );
        m_h1_sig->Add(h1->ProjectionY(Form("h_sig_phi2_Nch%d", index_Nch_low), 0.00001 + m_jetEtaIndexHigh+1, 0.00001 + m_etaBinIndexHigh,  "e" ) );

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
    h1->GetZaxis()->SetTitle("#it{C}(#Delta#it{#phi}, #Delta#it{#eta})");
    h1->GetXaxis()->CenterTitle();
    h1->GetYaxis()->CenterTitle();
    h1->GetZaxis()->CenterTitle();
    h1->GetZaxis()->SetTitleOffset(1.60);
    h1->Draw("SURF1FB");
    //ATLASLabel(0.02,0.94,"Internal");
    //myText(    0.02,0.89,1,"#it{p}+Pb #sqrt{#it{s}_{NN}} = 8.16 TeV");
    //myText(    0.02,0.83,1,"#it{h}-#it{h} Correlation");
    //myText(    0.02,0.12,1,"0.5 < #it{p}_{T}^{trig,asso} < 5 GeV");
    pad->SetTheta(60); // default is 30
    pad->SetPhi(45); // default is 30
    pad->Update();
}
