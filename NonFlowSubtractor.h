#ifndef NONFLOWSUBTRACTOR_H 
#define NONFLOWSUBTRACTOR_H

#include "FlowAnaConfig.h"
#include "SubtractResults.h"



// Create the GlobalCHi2 structure
struct GlobalChi2 {
    GlobalChi2(  ROOT::Math::IMultiGenFunction & f1,
                 ROOT::Math::IMultiGenFunction & f2,
                 int n1, int n2) :
                 fChi2_1(&f1), fChi2_2(&f2), nhar_LM(n1), nhar_HM(n2) {}

    // parameter vector is first background (in common 1 and 2)
    // and then is signal (only in 2)
    double operator() (const double *par) const {
        double p1[nhar_LM + 1];
        double p2[nhar_LM + nhar_HM + 3];

        for (int i = 0; i < nhar_LM+1; ++i) p1[i] = par[i];
        for (int i = 0; i < nhar_LM + nhar_HM + 3; ++i) p2[i] = par[i];

        return (*fChi2_1)(p1) + (*fChi2_2)(p2);
    }
    // number of orders of harmonics
    // passed from the sutractor class
    const int nhar_LM;
    const int nhar_HM;
    const ROOT::Math::IMultiGenFunction * fChi2_1;
    const ROOT::Math::IMultiGenFunction * fChi2_2;
};



class NonFlowSubtractor {

    const FlowAnaConfig* m_anaConfig;
    bool m_debug;

    // control how many orders of harmonics should be used 
    // default value is 4
    // However, 5 is used in previous MinBias publication
    static int m_nhar_LM;
    static int m_nhar_HM;
    static TH1F* m_hist_LMtemp; // histogram tempalte for LM
    static TH1F* m_hist_HMtemp; // histogram tempalte for HM
    static TH1F* m_hist_LMtemp_bulk; // histogram tempalte for LM of bulk-bulk correlation

  public:
    // interface to control Nhar for HM in atlas fit
    // used for Fourier fit and CMS mehtod as well
    static int getNHar() {return m_nhar_HM;}
    static void setNHar(int n) {m_nhar_HM = n;}

    // interface to control LM fit number of harmonics 
    static int getNHarLM() {return m_nhar_LM;}
    static void setNHarLM(int n) {m_nhar_LM = n;}

    TH1F* getChi2Hist() {return h_chi2_c2;}

    // interface for getting and setting template fit parameters
    std::vector<float> getParVector() {return m_parVector;} 
    void setParVector(std::vector<float> _vec) { m_parVector = _vec;} 

    void setInteration(bool _flag) {m_doIteration = _flag;}

  private:
    // total number of parameters of ATLAS fit for LM events
    int m_nparLmAtlas; // should be m_nhar_LM + 1
    std::vector<string> m_parName_LM;

    // total number of parameters of ATLAS fit for HM events
    int m_nparHmAtlas; // should be m_nhar_LM + m_nhar_HM + 3
    std::vector<string> m_parName_HM;
    std::vector<float> m_parVector;

    bool m_doIteration;

    std::string _formula_LM;
    std::string _formula_LM_ZYAM;
    std::string _formula;
    double m_dphiRangeLow;
    double m_dphiRangeHigh;

    float _global_chi2;
    int   _global_ndof;

  public:
    NonFlowSubtractor();
    virtual ~NonFlowSubtractor() {};

    void init();
    void setDebug (bool _debug = true) { m_debug = _debug; }
    void setFixLM (bool _fix = true) { m_fixLM = _fix; }
    void setFitPrintLevel (int _level) { m_fitPrintLevel = _level; }

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
    subResult referenceFit(TH1* hist_LM, TH1* hist_HM); // duplication of template fit with modification to parameterization and names

    // corrected version with non-zeor hist_LM2
    subResult templateFit(TH1* hist_LM, TH1* hist_HM, TH1* hist_LM2);

    // corrected version with external cn_LM value and error
    subResult templateFit(TH1* hist_LM, TH1* hist_HM, float cn_LM, float cn_LM_error);

    // using hist_LM as template instead of fitting it
    subResult templateHistFit(TH1* hist_LM, TH1* hist_HM);
    subResult templateHistFit(TH1* hist_LM, TH1* hist_HM, TH1* hist_LM2);
    subResult templateHistFit(TH1* hist_LM, TH1* hist_HM, float cn_LM, float cn_LM_error);

    // old methods
    // might be needed for future development
    subResult templateHistFit2(TH1* hist_LM, TH1* hist_HM);
    subResult templateHistFit2(TH1* hist_LM, TH1* hist_LM_bulk, TH1* hist_HM);

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
    // supportin function for template fitting
    subResult fourierFitLM (TH1* hist, bool _zyam = false);


    TH1F* getPseudoHist_HM () {return m_hist_pseudo_HM;}
    TH1F* getPseudoHist_LM () {return m_hist_pseudo_LM;}
    void  setPseudoVaryScale (float _scale) { m_pseudoVaryScale = _scale;}
    float getPseudoVaryScale () {return m_pseudoVaryScale;}


    // Additional interface to configure the ATLAS Template fitting
    //---------------------------------
    // flag to control c1 term, set to true (fixed at zero) by default 
    // otherwise the fit procedure have problems with converging
    // detailed in p+Pb HF flow supporting note.
    void setAtlasFixedC1 (bool _isFixed = true) {m_fixC1 = _isFixed;}

    // Same for c3 and c4, false by default
    void setAtlasFixedC3 (bool _isFixed = true) {if (getNHar() > 2) {m_fixC3 = _isFixed;} else {m_fixC3 = true;}}
    void setAtlasFixedC4 (bool _isFixed = true) {if (getNHar() > 3) {m_fixC4 = _isFixed;} else {m_fixC4 = true;}}

    // plot utilities for default ATLAS template fitting (no improved correction applied yet)
    // major plotting style with HM template fit and ridge fit
    bool plotAtlasSubHM(TCanvas* theCanvas); // big figure
    bool plotAtlas3pHM(TCanvas* theCanvas); // big figure
    bool plotAtlasHistSubHM(TCanvas* theCanvas); // big figure

    bool plotAtlasHM(TPad* thePad); // small figure
    void plotAtlasHMLabel(TPad* thePad); // small figure
    bool plotAtlasLM(TPad* thePad); // small figure
    bool plotAtlasHM(TCanvas* theCanvas); // big figure
    bool plotAtlasLM(TCanvas* theCanvas); // big figure
    bool plotAtlasHMSimple(TCanvas* theCanvas); // big figure

    bool plotAtlasHistHM(TCanvas* theCanvas); // small figure
    TH1F* getPullHist() {return m_h_pull;} // get pull for Brain's check

    // by default no ZYAM is applied for any method
    // one could apply ZYAM to low multiplicity reference in the ATLAS method
    // one could apply ZYAM to LM and HM using CMS method, in case of using ATLAS PTY defintion (projection first, then take ratio)
    void setZYAM (bool _isApplied = true) {m_applyZYAM = _isApplied;}

  private:
    // private members for ATLAS tempalte fit
    static double f_periph(double *x, double *par);
    static double f_periph_zyam(double *x, double *par);
    static double f_ridge (double *x, double *par);
    static double templ   (double *x, double *par);
    static double templ_zyam   (double *x, double *par);

    // LM hist template version
    static double templ_hist   (double *x, double *par);
    static double templ_hist2  (double *x, double *par);
    static double f_periph_hist(double *x, double *par);
    static double f_periph_hist2(double *x, double *par);

    // function for Seyoung's reference fit
    // for reference, we can use f_periph for now
    //static double seyoung_periph (double *x, double *par);
    static double seyoung_templ  (double *x, double *par);

    // function for hist-based template fit
    static void hist_fcn(int &npar, double *gin, double &f, double *par, int flag);

    bool m_fixC1;
    bool m_fixC3;
    bool m_fixC4;
    bool m_applyZYAM;
    bool m_fixLM;
    bool m_plotBulkRef;

    // print level 
    // 0 = Quiet (Default)
    // 1 = just fit results (Default in this tool)
    // 2 = fit results with iterations
    // 3 = Verbose (additional details of iterations)
    int m_fitPrintLevel;

    // members for toy fit model check
    TH1F* m_hist_pseudo_HM;
    TH1F* m_hist_pseudo_LM;
    float m_pseudoVaryScale;

    float m_ridge_scaleFactor;

    TH1F* m_hist_LM;
    TH1F* m_hist_HM;
    TH1F* m_h_pull;
    TH1F* m_h_ridge;


    TF1* f_LM;
    TF1* f_HM;

    TF1* f_show_periph;
    TF1* f_show_flow2; //with pedstal
    TF1* f_show_flow3; //with pedstal
    TF1* f_show_ridge; //w/o pedstal
    TF1* f_show_ridge2; //w/o pedstal
    TF1* f_show_ridge3; //w/o pedstal

    TH1F* h_show_periph;
    TH1F* h_show_periph_bulk;
    TH1F* h_show_HM;
    TH1F* h_chi2_c2;

};
#endif



int NonFlowSubtractor::m_nhar_LM = 4;
int NonFlowSubtractor::m_nhar_HM = 4;
TH1F* NonFlowSubtractor::m_hist_LMtemp = 0;
TH1F* NonFlowSubtractor::m_hist_HMtemp = 0;
TH1F* NonFlowSubtractor::m_hist_LMtemp_bulk = 0;



NonFlowSubtractor :: NonFlowSubtractor() {
    m_debug = false;
    m_fixC1 = true;
    m_fixC3 = false;
    m_fixC4 = false;
    m_fixLM = false;
    m_h_pull = 0;
    m_h_ridge = 0;

    m_doIteration = true;

    m_plotBulkRef = false;
    m_ridge_scaleFactor = 1;

    m_pseudoVaryScale = 0.0;

    m_dphiRangeLow = -0.5*TMath::Pi();
    m_dphiRangeHigh = 1.5*TMath::Pi();

    m_applyZYAM = false;
    m_fitPrintLevel = 0;

    _global_chi2 = 0;;
    _global_ndof = 0;;

    m_parVector.clear();
}



void NonFlowSubtractor :: init() {
    if (m_debug) cout << "Running Non-flow subtraction up to order " << getNHar() << endl;

    if (getNHar() < 3) m_fixC3 = true;
    if (getNHar() < 4) m_fixC4 = true;

    _formula = "[0] * ( 1 + 2*(";
    for (int ihar=0; ihar<getNHar(); ihar++){
        _formula += Form("[%d]*cos(%d*x)", ihar+1, ihar+1);
        if (ihar != getNHar()-1) _formula += (" + ");
    }
    _formula += (") )");

    // for ATLAS method, in case nhar is different for LM and HM
    _formula_LM = "[0] * ( 1 + 2*(";
    _formula_LM_ZYAM = "[0] * ( 2*(";
    for (int ihar=0; ihar<getNHarLM(); ihar++){
        _formula_LM += Form("[%d]*cos(%d*x)", ihar+1, ihar+1);
        _formula_LM_ZYAM += Form("[%d]*cos(%d*x)", ihar+1, ihar+1);
        if (ihar != getNHarLM()-1) {
            _formula_LM += (" + ");
            _formula_LM_ZYAM += (" + ");
        }
    }
    _formula_LM += (") )");
    _formula_LM_ZYAM += (") )");


    if (m_debug) cout << "Fourier function: " << _formula.c_str() << endl << endl;
    if (m_debug) cout << "Fourier function (atlas LM): " << _formula_LM.c_str() << endl << endl;

    // initial ATLAS method
    m_nparLmAtlas = getNHarLM() + 1;
    m_nparHmAtlas = getNHarLM() + getNHar() + 3;

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
    for (int ipar=1; ipar < getNHar()+1; ipar++) {
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


// G^{LM}*(1+2*_SIGMA(an*cos(n*dphi)))
// G^{LM} + _SIGMA 2*G^{LM}*(a_n*cos(n*dphi))
double NonFlowSubtractor::f_periph(double *x, double *par) {
    double xx = x[0];
    double value = (par[0]);
    for (int ihar=0; ihar < getNHarLM(); ihar++) {
        value += 2.*(par[0])*(par[ihar+1]*TMath::Cos((ihar+1)*xx));
    }
    return value;
}



// G^{LM}*(2*_SIGMA(an*cos(n*dphi)))
// _SIGMA 2*G^{LM}*(a_n*cos(n*dphi))
double NonFlowSubtractor::f_periph_zyam(double *x, double *par) {
    double xx = x[0];
    double value = 0;
    for (int ihar=0; ihar < getNHarLM(); ihar++) {
        value += 2.*(par[0])*(par[ihar+1]*TMath::Cos((ihar+1)*xx));
    }
    return value;
}



// 1+2*Sigma_n(c_n*Cos(n*dphi))
double NonFlowSubtractor::f_ridge(double *x, double *par) {
    double xx = x[0];
    double value = 1;
    for (int ihar=0; ihar < getNHar(); ihar++) {
        value += 2.*(par[ihar]*TMath::Cos((ihar+1)*xx));
    }
    return value;
}



// F*f_periph + (G^{HM} - F*G^{LM})*f_ridge
double NonFlowSubtractor::templ(double *x, double *par) {
    double xx = x[0];
    double f = par[getNHar()+getNHarLM()+1]*f_periph(x, par) + par[getNHarLM()+getNHar()+2]*f_ridge(x, &par[getNHarLM()+1]);
    return f;
}



// F*f_periph + G^{HM}*f_ridge + (G^{HM} - F * G^{LM})
double NonFlowSubtractor::seyoung_templ(double *x, double *par) {
    double xx = x[0];
    double f = par[getNHarLM()+getNHar()+1]*f_periph(x, par)  //F
             + par[getNHarLM()+getNHar()+2]*f_ridge(x, &par[getNHarLM()+1]) // G{HM}
             //+ par[getNHarLM()+getNHar()+2] - par[getNHar()+getNHarLM()+1]*par[0];
             - par[getNHar()+getNHarLM()+1]*par[0];
    return f;
}



// F*f_periph_zyam + (G^{HM} - F*G^{LM})*f_ridge
double NonFlowSubtractor::templ_zyam(double *x, double *par) {
    double xx = x[0];
    double f = par[getNHarLM()+getNHar()+1]*f_periph_zyam(x, par)  //F
             + par[getNHarLM()+getNHar()+2]*f_ridge(x, &par[getNHarLM()+1]); // G{temp} = G{HM}
    return f;
}



double NonFlowSubtractor::templ_hist(double *x, double *par) {
    double xx = x[0];
    double _pty_LM = par[getNHar()] * m_hist_LMtemp->GetBinContent(m_hist_LMtemp->FindBin(xx));
    double f = _pty_LM + par[getNHar()+1]*f_ridge(x, &par[0]);
    return f;
}


// for plotting
double NonFlowSubtractor::f_periph_hist(double *x, double *par) {
    double xx = x[0];
    double f = par[0]*m_hist_LMtemp->GetBinContent(m_hist_LMtemp->FindBin(xx))+par[1] ;
    return f;
}


double NonFlowSubtractor::templ_hist2(double *x, double *par) {
    double xx = x[0];
    double _pty_LM = par[getNHar()] * ( par[getNHar()+2]     * m_hist_LMtemp     ->GetBinContent(m_hist_LMtemp->FindBin(xx)) 
                                      + (1-par[getNHar()+2]) * m_hist_LMtemp_bulk->GetBinContent(m_hist_LMtemp_bulk->FindBin(xx)) 
                                      );
    double f = _pty_LM + par[getNHar()+1]*f_ridge(x, &par[0]);
    return f;
}

// for plotting
double NonFlowSubtractor::f_periph_hist2(double *x, double *par) {
    double xx = x[0];
    double f = par[0] * ( par[2] * m_hist_LMtemp     ->GetBinContent(m_hist_LMtemp->FindBin(xx)) 
                    + (1-par[2]) * m_hist_LMtemp_bulk->GetBinContent(m_hist_LMtemp_bulk->FindBin(xx))
                       ) + par[1] ;
    return f;
}


void NonFlowSubtractor::hist_fcn(int &npar, double *gin, double &f, double *par, int flag) {
    double chisq = 0;
    double delta;
    for (int i=1; i<=m_hist_HMtemp->GetNbinsX(); i++) {
        Double_t xvalue = m_hist_HMtemp->GetBinCenter(i);
        delta  = m_hist_HMtemp->GetBinContent(i) - templ_hist(&xvalue,par);
        
        chisq += delta*delta /
                ( pow(m_hist_HMtemp->GetBinError(i), 2) + pow(par[getNHar()]*m_hist_LMtemp->GetBinError(i), 2) );
    }
    f=chisq;
}
