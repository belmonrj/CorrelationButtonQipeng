#ifndef NONFLOWSUBTRACTOR_H 
#define NONFLOWSUBTRACTOR_H

#include "FlowAnaConfig.h"
#include "SubtractResults.h"



// Create the GlobalCHi2 structure
struct GlobalChi2 {
    GlobalChi2(  ROOT::Math::IMultiGenFunction & f1,
                 ROOT::Math::IMultiGenFunction & f2,
                 int n) :
                 fChi2_1(&f1), fChi2_2(&f2), nhar(n) {}

    // parameter vector is first background (in common 1 and 2)
    // and then is signal (only in 2)
    double operator() (const double *par) const {
        double p1[nhar+1];
        double p2[2*nhar+3];

        for (int i = 0; i < nhar+1; ++i) p1[i] = par[i];
        for (int i = 0; i < 2*nhar+3; ++i) p2[i] = par[i];

        return (*fChi2_1)(p1) + (*fChi2_2)(p2);
    }
    // number of orders of harmonics
    // passed from the sutractor class
    const int nhar;
    const ROOT::Math::IMultiGenFunction * fChi2_1;
    const ROOT::Math::IMultiGenFunction * fChi2_2;
};



class NonFlowSubtractor {

    const FlowAnaConfig* m_anaConfig;
    bool m_debug;

    // control how many orders of harmonics should be used 
    // default value is 4
    // However, 5 is used in previous MinBias publication
    static int m_nhar;
  public:
    static int getNHar() {return m_nhar;}
    static void setNHar(int n) {m_nhar = n;}

  private:
    // total number of parameters of ATLAS fit for LM events
    // index 0 ~ m_nhar
    int m_nparLmAtlas; // should be m_nhar + 1
    std::vector<string> m_parName_LM;

    // total number of parameters of ATLAS fit for HM events
    // index 0 ~ m_nhar for LM parameters, must match with m_par_LM
    // m_nhar + 1 ~ 2*m_nhar + 3 for rest parameters
    int m_nparHmAtlas; // should be 2*m_nhar + 3
    std::vector<string> m_parName_HM;

    std::string _formula;
    double m_dphiRangeLow;
    double m_dphiRangeHigh;


  public:
    NonFlowSubtractor();
    virtual ~NonFlowSubtractor() {};

    void init();
    void setDebug (bool _debug = true) { m_debug = _debug; }

    // interface for using AnaConfig classess
    // may be useless for other users
    void setConfig (const FlowAnaConfig* _config) { m_anaConfig = _config; }
    const FlowAnaConfig* getConfig () const { return m_anaConfig; }

    //---------------------------------
    // Different subtraction method
    //---------------------------------
    // input should be 1D DeltaPhi per-trigger pair yield
    // One can use uncorrected templateFit to self-normalized correlation function C(DeltaPhi)
    // !!NEVER!! used improved methods for self-normalized correlation function C(DeltaPhi)
    
  public:
    //---------------------------------
    // Atlas template fit method, if hist_LM2 = 0, the 1 step correction (HION-2017-07) will be applied to cX_corr_value in subResult
    //---------------------------------
    // uncorrected version
    subResult templateFit(TH1* hist_LM, TH1* hist_HM);

    // corrected version with non-zeor hist_LM2
    subResult templateFit(TH1* hist_LM, TH1* hist_HM, TH1* hist_LM2);

    // corrected version with external cn_LM value and error
    subResult templateFit(TH1* hist_LM, TH1* hist_HM, float cn_LM, float cn_LM_error);

    //---------------------------------
    // peripheral subtraction used by CMS
    //---------------------------------
    // Corrected cX_corr_value are the subtract ones and cX_value are identical with Fourier fit
    // with option to run non-zero LM flow correction
    // orders of the histograms matter, be careful
    // default version used in CMS publications with low mulitplicty has zero flow
    subResult periphSub  (TH1* h_sr_lm, TH1* h_lr_lm, 
                          TH1* h_sr_hm, TH1* h_lr_hm);

    // orders of the histograms matter, be careful
    // extended version with non-zero low multiplicity flow
    subResult periphSub  (TH1* h_sr_lm,  TH1* h_lr_lm, 
                          TH1* h_sr_hm,  TH1* h_lr_hm, 
                          TH1* h_sr_lm2, TH1* h_lr_lm2);

    // corrected version with external cn_LM value and error
    subResult periphSub  (TH1* h_sr_lm,  TH1* h_lr_lm, 
                          TH1* h_sr_hm,  TH1* h_lr_hm, 
                          float cn_LM,  float cn_LM_error);


    //---------------------------------
    // Direct Fourier fit
    //---------------------------------
    // Extracted unsubtracted flow in p+Pb
    // Could be used to study the rho in Pythia if the second histogram is not zero
    subResult fourierFit (TH1* hist, TH1* hist_LM = 0);



    // Additional interface to configure the ATLAS Template fitting
    //---------------------------------
    // flag to control c1 term, set to true (fixed at zero) by default 
    // otherwise the fit procedure have problems with converging
    // detailed in p+Pb HF flow supporting note.
    void setAtlasFixedC1 (bool _isFixed = true) {m_fixC1 = _isFixed;}
    // Same for c3, false by default
    // not necessary actually
    void setAtlasFixedC3 (bool _isFixed = true) {m_fixC3 = _isFixed;}

    // plot utilities for default ATLAS template fitting (no improved correction applied yet)
    bool plotAtlasHM(TPad* thePad); // small figure
    void plotAtlasHMLabel(TPad* thePad); // small figure
    bool plotAtlasLM(TPad* thePad); // small figure
    bool plotAtlasHM(TCanvas* theCanvas); // big figure
    bool plotAtlasLM(TCanvas* theCanvas); // big figure


    // no support for ZYAM in Atlas method to void confusing for improved fit correction
    // becomes a support method for running CMS method using ATLAS PTY definition
    void setZYAM (bool _isApplied = true) {m_applyZYAM = _isApplied;}

  private:
    // private members for ATLAS tempalte fit
    static double f_periph(double *x, double *par);
    static double f_ridge (double *x, double *par);
    static double templ   (double *x, double *par);

    bool m_fixC1;
    bool m_fixC3;
    bool m_applyZYAM;

    TH1F* m_hist_LM;
    TH1F* m_hist_HM;
    TF1* f_LM;
    TF1* f_HM;

    TF1* f_show_periph;
    TF1* f_show_flow2; //c2
    TF1* f_show_flow3; //c3

};
#endif

int NonFlowSubtractor::m_nhar = 4;

NonFlowSubtractor :: NonFlowSubtractor() {
    m_debug = false;
    m_fixC1 = true;
    m_fixC3 = false;

    m_dphiRangeLow = -0.5*TMath::Pi();
    m_dphiRangeHigh = 1.5*TMath::Pi();
}

void NonFlowSubtractor :: init() {
    if (m_debug) cout << "Running Non-flow subtraction up to order " << getNHar() << endl;

    _formula = "[0] * ( 1 + 2*(";
    for (int ihar=0; ihar<getNHar(); ihar++){
        _formula += Form("[%d]*cos(%d*x)", ihar+1, ihar+1);
        if (ihar != getNHar()-1) _formula += (" + ");
    }
    _formula += (") )");

    if (m_debug) cout << "Fourier function: " << _formula.c_str() << endl << endl;

    // initial ATLAS method
    m_nparLmAtlas = getNHar() + 1;
    m_nparHmAtlas = 2*getNHar() + 3;

    m_parName_LM.clear();
    m_parName_HM.clear();
    for (int ipar=0; ipar < m_nparLmAtlas; ipar++) {
        if (ipar == 0) {
            m_parName_LM.push_back("G_LM");
            m_parName_HM.push_back("G_LM");
        } else {
            m_parName_LM.push_back(Form("a%d_LM",ipar));
            m_parName_HM.push_back(Form("a%d_LM",ipar));
        }
    }
    for (int ipar=1; ipar < m_nparLmAtlas; ipar++) {
        m_parName_HM.push_back(Form("c%d_HM",ipar));
    }

    m_parName_HM.push_back("F_temp");
    m_parName_HM.push_back("G_temp");

    if (m_debug) {
        cout << " ===================================================" << endl;
        cout << " ATLAS template fitting parameterization:" << endl;
        cout << "  number of LM event parameters:" << m_nparLmAtlas << endl;
        for (int ipar=0; ipar < m_parName_LM.size(); ipar++) {
            cout << "   parameter "  << ipar << " :" << m_parName_LM.at(ipar).c_str() << endl;       
        }
        cout << endl;
        cout << "  number of HM event parameters:" << m_nparHmAtlas << endl;
        for (int ipar=0; ipar < m_parName_HM.size(); ipar++) {
            cout << "   parameter "  << ipar << " :" << m_parName_HM.at(ipar).c_str() << endl;       
        }
        cout << endl;
    }


}


// G^{LM}*(1+2*_SIGMA(an*cos(n*dphi))), n=1~4
// G^{LM} + _SIGMA 2*G^{LM}*(a_n*cos(n*dphi))
// 5 parameters in total
double NonFlowSubtractor::f_periph(double *x, double *par) {
    double xx = x[0];
    double value = (par[0]);
    for (int ihar=0; ihar < getNHar(); ihar++) {
        value += 2.*(par[0])*(par[ihar+1]*TMath::Cos((ihar+1)*xx));
    }
    return value;
}


// 1+2*Sigma_n(c_n*Cos(n*dphi)), n=1~4
// 4 parameters in total
// c1, c2, c3, c4
double NonFlowSubtractor::f_ridge(double *x, double *par) {
    double xx = x[0];
    double value = 1;
    for (int ihar=0; ihar < getNHar(); ihar++) {
        value += 2.*(par[ihar]*TMath::Cos((ihar+1)*xx));
    }
    return value;
}


// F*f_periph + (G^{HM} - F*G^{LM})*f_ridge
// 11 parameters in total
double NonFlowSubtractor::templ(double *x, double *par) {
    double xx = x[0];
    double f = par[2*getNHar()+1]*f_periph(x, par) + par[2*getNHar()+2]*f_ridge(x, &par[getNHar()+1]);
    return f;
}
