#ifndef FLOWANACONFIG_H
#define FLOWANACONFIG_H

#include "AnaHelper.h"

using namespace std;
using std::string;

class FlowAnaConfig {

    // reference bin Nch range 
    float m_multiRef_low;
    float m_multiRef_high;
    // second loweset LM Nch range for improved fit results
    float m_multiRef2_low;
    float m_multiRef2_high;
    
    float m_corrEtaBoundary;
    float m_corrEtaInterval;

    // gap cut for abs(DeltaEta)
    float m_etaRange_low;
    float m_etaRange_high;

    // information might be needed for ploting 
    float m_assoPt_low; // associate particle pT lower boundary
    float m_assoPt_high; // associate particle pT upper boundary

    string m_inputPath;
    string m_outputPath;
    string m_outputFigPath;

    string m_inputFileName;
    string m_outputFileName;
    TFile* m_inputFile;
    TFile* m_outputFile;

    // variable pointer to constant TH1F objects
    const TH1F* m_input_multiBinning;
    const TH1F* m_input_ptBinning;
    const TH1F* m_input_3rdBinning;

    string m_inputName_multiBinning;
    string m_inputName_ptBinning;
    // binning for dimension that used for background subtraction in hard-soft correlation
    // usually it's invariant mass, could be momemtum inbalance for mu-hadron correlation
    string m_inputName_3rdBinning; // name for input binning TH1F
    string m_Name_3rdBinning; // naame for dimension

    const TH1F* m_output_multiBinning;
    const TH1F* m_output_ptBinning;
    const TH1F* m_output_3rdBinning;

    string m_corrHistPathSame;
    string m_corrHistPathMix;
    string m_triggerYiledHistName;
    string m_mixTriggerYiledHistName;

  public:

    FlowAnaConfig();
    virtual ~FlowAnaConfig() {};
    void checkBinningConsistency() const;
    bool init();
    void print() const; 

    // inline member functions
    TFile* inputFile() const { return m_inputFile; }
    TFile* outputFile() const { return m_outputFile; }

    void  setCorrEtaBoundary(float _value) { m_corrEtaBoundary = _value; }
    void  setCorrEtaInterval(float _value) { m_corrEtaInterval = _value; }
    float getCorrEtaBoundary() const { return m_corrEtaBoundary; }
    float getCorrEtaInterval() const { return m_corrEtaInterval; }

    void setReference(float _low, float _high) { m_multiRef_low = _low; m_multiRef_high = _high; }
    void setReference2(float _low, float _high) { m_multiRef2_low = _low; m_multiRef2_high = _high; }
    void setEtaRange (float _low, float _high) { m_etaRange_low = _low; m_etaRange_high = _high; }
    void setAssoPt   (float _low, float _high) { m_assoPt_low   = _low; m_assoPt_high = _high; }
    float getReferenceLow()  const { return m_multiRef_low; }
    float getReferenceHigh() const { return m_multiRef_high; }
    float getReference2Low()  const { return m_multiRef2_low; }
    float getReference2High() const { return m_multiRef2_high; }
    float getEtaRangeLow()   const { return m_etaRange_low; }
    float getEtaRangeHigh()  const { return m_etaRange_high; }
    float getAssoPtLow()     const { return m_assoPt_low; }
    float getAssoPtHigh()    const { return m_assoPt_high; }

    void setInputPath     (string _path) { m_inputPath = _path; }
    void setOutputPath    (string _path) { m_outputPath = _path; }
    void setOutputFigPath (string _path) { m_outputFigPath = _path; }
    string getInputPath     () const { return m_inputPath; }
    string getOutputPath    () const { return m_outputPath; }
    string getOutputFigPath () const { return m_outputFigPath; }

    void setInputFileName  (string _name) { m_inputFileName = _name; }
    void setOutputFileName (string _name) { m_outputFileName = _name; }
    string getInputFileName  () const { return m_inputFileName; }
    string getOutputFileName () const { return m_outputFileName; }

    void setInputMultiBinningName  (string _name) { m_inputName_multiBinning = _name; }
    void setInputPtBinningName     (string _name) { m_inputName_ptBinning = _name; }
    void setInputThirdBinningName  (string _name) { m_inputName_3rdBinning = _name; }

    void setThirdDimensionName  (string _name) { m_Name_3rdBinning = _name; }
    string getThirdDimensionName  () const { return m_Name_3rdBinning; }

    void setOutputMultiBinning (const int _nbins, const float* _range ) { m_output_multiBinning = new TH1F("_output_multiBinning","",_nbins, _range); }
    void setOutputPtBinning    (const int _nbins, const float* _range ) { m_output_ptBinning    = new TH1F("_output_ptBinning","",_nbins, _range); }
    void setOutputThirdBinning    (const int _nbins, const float* _range ) { m_output_3rdBinning = new TH1F("_output_3rdBinning","",_nbins, _range); }


    const TH1F* getOutputMultiBinning() const {return m_output_multiBinning;}
    const TH1F* getOutputPtBinning() const {return m_output_ptBinning;}
    const TH1F* getOutputThirdBinning() const {return m_output_3rdBinning;}

    const TH1F* getInputMultiBinning() const {return m_input_multiBinning;}
    const TH1F* getInputPtBinning() const    {return m_input_ptBinning;}
    const TH1F* getInputThirdBinning() const {return m_input_3rdBinning;}

    void setCorrHistPathSame (string _path) { m_corrHistPathSame = _path; }
    void setCorrHistPathMix  (string _path) { m_corrHistPathMix = _path;  }
    string getCorrHistPathSame () const { return m_corrHistPathSame; }
    string getCorrHistPathMix  () const { return m_corrHistPathMix;  }

    void setTrigYieldHistName  (string _name) { m_triggerYiledHistName = _name;  }
    void setMixTrigYieldHistName  (string _name) { m_mixTriggerYiledHistName = _name;  }
    string getTrigYieldHistName  () const { return m_triggerYiledHistName;  }
    string getMixTrigYieldHistName  () const { return m_mixTriggerYiledHistName;  }
};


FlowAnaConfig::FlowAnaConfig() {

    // additional dimensions needs to be added
    // like invariant mass for hard-soft correlation
    m_input_multiBinning = 0;
    m_input_ptBinning = 0;
    m_input_3rdBinning = 0;

    m_output_multiBinning = 0;
    m_output_ptBinning = 0;
    m_output_3rdBinning = 0;

    // default value of paths and names
    // meaningless for different analyzer 
    m_inputPath = "../Selection/";
    m_outputPath = "../ButtonRootFiles/";
    m_outputFigPath = "../ButtonPlots/";
    m_inputName_multiBinning = "hNch_binning";
    m_inputName_ptBinning    = "hpt_binning";
    m_inputName_3rdBinning = "empty";

    m_corrHistPathSame = "Correlation_2D_raw/h2pc_sig_";
    m_corrHistPathMix  = "Correlation_2D_raw/h2pc_mix_";
    m_triggerYiledHistName = "h2_nsig_raw";
    m_mixTriggerYiledHistName = "h2_nsig_raw"; // to be changed for cms defintion

    // default ATLAS associated particle selection
    // it is not really controlled by analysis at this stage 
    // just for information
    m_assoPt_low = 0.5;
    m_assoPt_high = 5;

    m_etaRange_low = 2.0;
    m_etaRange_high = 5.0;

    m_corrEtaBoundary = 5.0; // only for hh, 4.5 for HF-h correlation
    m_corrEtaInterval = 0.1;

    m_multiRef_low = 10;
    m_multiRef_high = 20;
    m_multiRef2_low = 20;
    m_multiRef2_high = 30;
}


void FlowAnaConfig::checkBinningConsistency() const {
    if (m_input_multiBinning == 0 ) {
        cout << "cannot find input multiplicity binning histogram !" << endl;
    } else {
        if (m_output_multiBinning == 0) { 
            cout << "cannot find output multiplicity binning histogram !" << endl;
        } else {
            const double* _multiOutput_range = m_output_multiBinning->GetXaxis()->GetXbins()->GetArray();
            const double* _multiInput_range  =  m_input_multiBinning->GetXaxis()->GetXbins()->GetArray();
            for (int obin = 0; obin < m_output_multiBinning->GetXaxis()->GetNbins() + 1; obin++) {
                bool _hasInconsistency = true;
                for (int ibin = 0; ibin < m_input_multiBinning->GetXaxis()->GetNbins() + 1; ibin++) {
                    if (_multiOutput_range[obin] == _multiInput_range[ibin]) {
                        _hasInconsistency = false;
                        break;
                    }
                }
                if (_hasInconsistency) cout << "Find inconsistency between input and output multiplicity binning" << endl;
            }
        } 
    }


    if (m_input_ptBinning == 0 ) {
        cout << "cannot find input pt binning histogram !" << endl;
    } else {
        if (m_output_ptBinning == 0) { 
            cout << "cannot find output pt binning histogram !" << endl;
        } else {
            const double* _ptOutput_range = m_output_ptBinning->GetXaxis()->GetXbins()->GetArray();
            const double* _ptInput_range  =  m_input_ptBinning->GetXaxis()->GetXbins()->GetArray();
            for (int obin = 0; obin < m_output_ptBinning->GetXaxis()->GetNbins() + 1; obin++) {
                bool _hasInconsistency = true;
                for (int ibin = 0; ibin < m_input_ptBinning->GetXaxis()->GetNbins() + 1; ibin++) {
                    if ( fabs(_ptOutput_range[obin] - _ptInput_range[ibin]) < 0.01) {
                        _hasInconsistency = false;
                        break;
                    }
                }
                if (_hasInconsistency) {
                    cout << "Find inconsistency between input and output pt binning" << endl;
                }
            }
        } 
    }

    if (m_input_3rdBinning && m_output_3rdBinning) { 
        const double* _Output_range = m_output_3rdBinning->GetXaxis()->GetXbins()->GetArray();
        const double* _Input_range  =  m_input_3rdBinning->GetXaxis()->GetXbins()->GetArray();
        for (int obin = 0; obin < m_output_3rdBinning->GetXaxis()->GetNbins() + 1; obin++) {
            bool _hasInconsistency = true;
            for (int ibin = 0; ibin < m_input_3rdBinning->GetXaxis()->GetNbins() + 1; ibin++) {
                if (_Output_range[obin] == _Input_range[ibin]) {
                    _hasInconsistency = false;
                    break;
                }
            }
            if (_hasInconsistency) cout << "Find inconsistency between input and output 3rd binning" << endl;
        }
    } 


}



bool FlowAnaConfig::init() {
    bool _isAllRight = true;
    m_inputFile = new TFile(Form("%s/%s.root",m_inputPath.c_str(), m_inputFileName.c_str() ),"READ");
    if (m_inputFile->IsZombie()) {
        // Root will print out the error information
        return false;
    }
    m_input_multiBinning = (TH1F*) m_inputFile->Get(m_inputName_multiBinning.c_str());
    m_input_ptBinning = (TH1F*) m_inputFile->Get(m_inputName_ptBinning.c_str());
    if (m_inputName_3rdBinning.compare("empty") != 0) m_input_3rdBinning = (TH1F*) m_inputFile->Get(m_inputName_3rdBinning.c_str());

    m_outputFile = new TFile(Form("%s/%s.root",m_outputPath.c_str(), m_outputFileName.c_str() ),"RECREATE");
    if (m_outputFile->IsZombie()) {
        return false;
    }
    return _isAllRight;
}



// print information of the configuration
void FlowAnaConfig::print() const {
    cout << "----------------------------------------------------" << endl;
    if (m_inputFile) cout << "Open file: " << m_inputFile->GetName() << endl;
    else cout << "  Error:: Cannot read the input file !!" << endl;

    if (m_input_multiBinning) {
        cout << "----------------------------------------------------" << endl;
        cout << "Input multiplicity binning: " << endl;
        PrintHistBinning (m_input_multiBinning);
        if (m_output_multiBinning) {
            cout << "Output multiplicity binning: " << endl;
            PrintHistBinning (m_output_multiBinning);
        }
    }

    if (m_input_ptBinning) {
        cout << "----------------------------------------------------" << endl;
        cout << "Input pt binning: " << endl;
        PrintHistBinning (m_input_ptBinning);
        if (m_output_ptBinning) {
            cout << "Output pt binning: " << endl;
            PrintHistBinning (m_output_ptBinning);
        }
    }
    
    if (m_input_3rdBinning) {
        cout << "----------------------------------------------------" << endl;
        cout << "Input " << m_Name_3rdBinning.c_str() << " binning: " << endl;
        PrintHistBinning (m_input_3rdBinning);
        if (m_output_3rdBinning) {
            cout << "Output " << m_Name_3rdBinning.c_str() << "  binning: " << endl;
            PrintHistBinning (m_output_3rdBinning);
        }
    }

    checkBinningConsistency();

    cout << endl;
    cout << "----------------------------------------------------" << endl;
    cout << "Running the 2pc correlation analysis with the following setup: " << endl;

    cout << "Deta correlation boundary = " << getCorrEtaBoundary() << endl;
    cout << "Deta correlation Interval = " << getCorrEtaInterval() << endl;

    cout << "LM Reference binning: " << endl; 
    cout << "\tLM : " << getReferenceLow() << " <= Nch < " << getReferenceHigh() << endl;
    cout << "LM2 Reference binning: " << endl; 
    cout << "\tLM2: " << getReference2Low() << " <= Nch < " << getReference2High() << endl;
    cout << "dEta cut: " << endl;
    cout << "\t" << getEtaRangeLow() << " < |dEta| < " << getEtaRangeHigh() << endl;
    cout << "associate particle selection: " << endl;
    cout << "\t" << getAssoPtLow() << " < pt^{asso} < " << getAssoPtHigh() << " GeV" << endl;
    cout << "----------------------------------------------------" << endl;
    cout << endl;
}

#endif
