#ifndef CORRELATIONMAKER_H
#define CORRELATIONMAKER_H

#include "FlowAnaConfig.h"

class CorrelationMaker {
    
    const FlowAnaConfig* m_anaConfig;
    bool m_debug;

    float m_nEtaBinInterval;
    float m_etaBoundary;
    float m_etaBinIndexLow;
    float m_etaBinIndexHigh;
    float m_jetEtaIndexLow;
    float m_jetEtaIndexHigh;

    // in case one wants to check the 1D signal and 1D mix
    TH1* m_h1_sig;
    TH1* m_h1_mix;

  public:
    CorrelationMaker();
    virtual ~CorrelationMaker() {};
    void init();

    void setDebug (bool _debug) { m_debug = _debug; }
    void setConfig (const FlowAnaConfig* _config) { m_anaConfig = _config; }
    const FlowAnaConfig* getConfig () const { return m_anaConfig; }

    void get1DSigHist() const {return m_h1_sig;}
    void get1DMixHist() const {return m_h1_mix;}

    void generateHist_multi();
    void generateHist_pt();
    // additional dimensions
    void generateHist_bkgSub();

    // main function for generating single 1D correlation histogram
    // method1 -- ATLAS mixed event normalization, projection first, then take ratio (default ATLAS method)
    // method2 -- ATLAS mixed event normalization, take 2D ratio, then 1D projection
    // method3 -- CMS mixed event normalization, projection first, then take ratio
    // method4 -- CMS mixed event normalization, take 2D ratio, then 1D projjection (default CMS method)
    // method5 -- no mixed event correction  (for systematic study of mixing)
    virtual TH1* MakeCorr(float _pt_low, float _pt_high, float _Nch_low, float _Nch_high, int method = 1);
};
#endif


CorrelationMaker::CorrelationMaker() {
    m_anaConfig = 0;
    m_debug = false;

    m_nEtaBinInterval = 0.1;
    m_etaBoundary = 5.0;

    m_h1_sig = 0;
    m_h1_mix = 0;
}



void CorrelationMaker::init() {
    m_anaConfig->print();

    // two index to exclude deta range due to upper boundary cut
    m_etaBinIndexLow  = (m_etaBoundary - m_anaConfig->getEtaRangeHigh() )/m_nEtaBinInterval + 1;
    m_etaBinIndexHigh = 2.0*m_etaBoundary/m_nEtaBinInterval - (m_etaBoundary - m_anaConfig->getEtaRangeHigh())/m_nEtaBinInterval;

    // two index to exclude deta range due to lower boundary cut
    m_jetEtaIndexLow = (m_anaConfig->getEtaRangeHigh() - m_anaConfig->getEtaRangeLow())/m_nEtaBinInterval + m_etaBinIndexLow;
    m_jetEtaIndexHigh = m_etaBinIndexHigh - (m_anaConfig->getEtaRangeHigh() - m_anaConfig->getEtaRangeLow())/m_nEtaBinInterval;

    //const TH1* _hist = MakeCorr( 1, 1, 1, 1);
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

            nsig += h2_nsig->GetBinContent(ipt+1, iNch+1);
            nsigMix += h2_nsigMix->GetBinContent(ipt+1, iNch+1);

	        TH2* htemp_1 = (TH2*)gDirectory->Get(Form("%sNch%d_pt%d",path_sig.c_str(),iNch,ipt));
	        TH2* htemp_2 = (TH2*)gDirectory->Get(Form("%sNch%d_pt%d",path_mix.c_str(),iNch,ipt));
            if (!htemp_1->GetDefaultSumw2()) {
                htemp_1->Sumw2(kTRUE);
                htemp_2->Sumw2(kTRUE);
            }

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

        hphi = (TH1*)m_h1_sig->Clone(Form("hphi_pt%d", index_pt_low));
        hphi->Divide(m_h1_mix);

        float int_C = hphi->Integral();
        float K = int_S / int_C;
        hphi->Scale(K/nsig, "width"); //per trigger & per deltaPhi yield

    } else if (method == 2) {
        // ATLAS style mixed event normalization
        // take 2D ratio, then project to 1D
        
        // start with d^{2}N_{pair}/(deta dphi)
        float int_S = h1->Integral();

        h1->Scale(1./h1->Integral());
        h2->Scale(1./h2->Integral());
        TH2F* h_correlation = (TH2F*)h1->Clone();
        h_correlation->Divide(h2);
        float int_C = h_correlation->Integral();
        float K = int_S / int_C;
        h_correlation->Scale(K/nsig, "width"); //per trigger & per deltaPhi yield

        hphi    = h_correlation->ProjectionY(Form("h_sig_phi_Nch%d", index_Nch_low), m_etaBinIndexLow, m_jetEtaIndexLow-1, "e");
        hphi->Add(h_correlation->ProjectionY(Form("h_sig_phi2_Nch%d",index_Nch_low), m_jetEtaIndexHigh+1, m_etaBinIndexHigh, "e"));
        hphi->Scale(1./(nbins_low + nbins_high));

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

        hphi = (TH1*)m_h1_sig->Clone(Form("hphi_pt%d", index_pt_low));
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
        // place holder for now

    } else {
        cout << "method is not supported, please choose from 1~3" << endl;
    }

    if (m_debug) cout << "------------------------------------------" << endl;
    delete h1;
    delete h2;

    hphi->SetName(Form("h_pty_dphi_Nch%.0fto%.0f_pt%.0fto%.0f",_Nch_low,_Nch_high,_pt_low,_pt_high));

    return hphi;
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