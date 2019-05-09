#include "NonFlowSubtractor.h"


subResult NonFlowSubtractor :: templateHistFit(TH1* hist_LM, TH1* hist_HM) {
    m_hist_LMtemp = (TH1F*)hist_LM->Clone("hist_LMtemp");

    f_HM = new TF1("f_HM", templ_hist, m_dphiRangeLow, m_dphiRangeHigh, getNHar()+2);

    for (int ipar=0; ipar<getNHar(); ipar++) {
        f_HM->SetParName(ipar,Form("c%d",ipar+1));
    }
    f_HM->SetParName(getNHar(),"F");
    f_HM->SetParName(getNHar()+1,"G");

    f_HM->FixParameter(0,0);

    hist_HM->Fit("f_HM", "0US");
    m_hist_HM = (TH1F*)hist_HM->Clone("__hist_HM");

    subResult theResult;
    vector<float> vec_value_sub;
    vector<float> vec_error_sub;

    for (int i=0; i<getNHar(); i++) {
        vec_value_sub.push_back(f_HM->GetParameter(i));
        vec_error_sub.push_back(f_HM->GetParError(i));

    }
    theResult.setCoeffSub(vec_value_sub, vec_error_sub);

    // construct stuffs for plotting
    f_show_periph = new TF1("f_show_periph", f_periph_hist, m_dphiRangeLow, m_dphiRangeHigh, 2);
    f_show_periph->SetParameter(0, f_HM->GetParameter(getNHar()));
    f_show_periph->SetParameter(1, f_HM->GetParameter(getNHar()+1));

    f_show_flow2 = new TF1("f_show_flow2","[0] + 2*[1]*[2]*cos(2*x)", m_dphiRangeLow, m_dphiRangeHigh);
    f_show_flow2->SetParameter(0, f_show_periph->Eval(0));
    f_show_flow2->SetParameter(1, f_HM->GetParameter(getNHar()+1) );
    f_show_flow2->SetParameter(2, f_HM->GetParameter(1)); // c2

    if (!m_fixC3) {  
        f_show_flow3 = new TF1("f_show_flow3","[0] + [1]*2*[2]*cos(3*x)", m_dphiRangeLow, m_dphiRangeHigh);
        f_show_flow3->SetParameter(0, f_show_periph->Eval(0));
        f_show_flow3->SetParameter(1, f_HM->GetParameter(getNHar()+1));
        f_show_flow3->SetParameter(2, f_HM->GetParameter(2)); // c3
    }

    h_show_periph = (TH1F*)hist_LM->Clone("h_show_periph");
    h_show_HM = (TH1F*)hist_LM->Clone("h_show_HM");
    h_show_periph->Reset();
    h_show_HM->Reset();

    for (int ibin = 1; ibin < h_show_periph->GetNbinsX()+1; ibin++) {
        float _xx = h_show_periph->GetXaxis()->GetBinCenter(ibin);
        h_show_periph->SetBinContent(ibin, f_show_periph->Eval(_xx));
        h_show_HM->SetBinContent(ibin, f_HM->Eval(_xx));
    }

    m_h_ridge = (TH1F*) m_hist_HM->Clone("h_ridge");
    m_h_ridge->Reset();
    for (int i=1; i<m_hist_HM->GetXaxis()->GetNbins()+1; i++){
        float _residual = m_hist_HM->GetBinContent(i) - h_show_periph->GetBinContent(i);
        _residual *= m_ridge_scaleFactor;
        float _residual_error = m_hist_HM->GetBinError(i)*m_ridge_scaleFactor;
        m_h_ridge->SetBinContent(i,_residual);
        m_h_ridge->SetBinError  (i,_residual_error);
    }

    f_show_ridge2 = new TF1("f_show_ridge2","[0]*2*[1]*cos(2*x)", m_dphiRangeLow, m_dphiRangeHigh);
    f_show_ridge2->SetParameter(0, f_HM->GetParameter(getNHar()+1)*m_ridge_scaleFactor);
    f_show_ridge2->SetParameter(1, f_HM->GetParameter(1)); // c2

    if (!m_fixC3) {  
        f_show_ridge3 = new TF1("f_show_ridge3","[0]*2*[1]*cos(3*x)", m_dphiRangeLow, m_dphiRangeHigh);
        f_show_ridge3->SetParameter(0, f_HM->GetParameter(getNHar()+1)*m_ridge_scaleFactor);
        f_show_ridge3->SetParameter(1, f_HM->GetParameter(2)); // c3
    }

    std::string _formula_ridge;
    _formula_ridge = "[0]*2*([1]*cos(2*x)";
    if (!m_fixC3) {
        _formula_ridge += (" + [2]*cos(3*x)");
    }
    if (!m_fixC4) {
        _formula_ridge += (" + [3]*cos(4*x)");
    }
    _formula_ridge += (")");

    f_show_ridge = new TF1("f_show_ridge",_formula_ridge.c_str(), m_dphiRangeLow, m_dphiRangeHigh);
    f_show_ridge->SetParameter(0, f_HM->GetParameter(getNHar()+1)*m_ridge_scaleFactor);
    f_show_ridge->SetParameter(1, f_HM->GetParameter(1)); // c2
    if (!m_fixC3) {
        f_show_ridge->SetParameter(2, f_HM->GetParameter(2)); // c3
    }
    if (!m_fixC4) {
        f_show_ridge->SetParameter(3, f_HM->GetParameter(3)); // c4
    }

    /*
    // calculate chi2 distribution vs. c2
    //
    TF1* f_chi2Copy_HM = (TF1*)f_HM->Clone("f_chi2Copy_HM");

    float _paraLow = 0;
    float _paraHigh = 0.01;
    int n_step = 20;
    float _step_size = (_paraHigh - _paraLow)/n_step;
    h_chi2_c2 = new TH1F("h_chi2_c2","",n_step,_paraLow,_paraHigh);

    for (int istep = 0; istep < n_step; istep++) {
        float _v22 = _paraLow + istep*_step_size;

        double* _paras = f_chi2Copy_HM->GetParameters();
        cout << f_chi2Copy_HM->GetParName(1) << ": " <<  f_chi2Copy_HM->GetParameter(1) << ",\t";
        //f_chi2Copy_HM->FixParameter(1, _v22);
        f_chi2Copy_HM->SetParameter(1, _v22);
        //f_chi2Copy_HM->Update();
        cout << f_chi2Copy_HM->GetParameter(1) << ",\t";

        double _x[1] = {2.};
        //cout << f_chi2Copy_HM->EvalPar(_x,_paras) << endl;
        cout << f_chi2Copy_HM->Eval(2.) << endl;

        float _chi2 = 0;
        for (int i=1; i<m_hist_HM->GetXaxis()->GetNbins()+1; i++){
            float _residual_HM = m_hist_HM->GetBinContent(i) - f_chi2Copy_HM->Eval(m_hist_HM->GetBinCenter(i));

            //cout << _residual_HM << ",\t " << _residual_LM << endl;
            //_chi2 += pow(_residual_HM/m_hist_HM->GetBinError(i),2) + pow(_residual_LM/m_hist_LM->GetBinError(i),2);
            _chi2 += pow(_residual_HM/m_hist_HM->GetBinError(i),2);
            //cout << _v22 << ",\t" << _residual_HM << endl;
        } 

        h_chi2_c2->SetBinContent(istep+1, _chi2);
        h_chi2_c2->SetBinError(istep+1, 0);
    }
    */

    return theResult;
}



subResult NonFlowSubtractor :: templateHistFit(TH1* hist_LM, TH1* hist_LM_bulk, TH1* hist_HM) {
    m_hist_LMtemp = (TH1F*)hist_LM->Clone("hist_LMtemp");
    m_hist_LMtemp_bulk = (TH1F*)hist_LM_bulk->Clone("hist_LMtemp_bulk");

    f_HM = new TF1("f_HM", templ_hist2, m_dphiRangeLow, m_dphiRangeHigh, getNHar()+3);

    for (int ipar=0; ipar<getNHar(); ipar++) {
        f_HM->SetParName(ipar,Form("c%d",ipar+1));
    }
    f_HM->SetParName(getNHar(),"F");
    f_HM->SetParName(getNHar()+1,"G");
    f_HM->SetParName(getNHar()+2,"p");

    f_HM->FixParameter(0,0);
    f_HM->SetParameter(getNHar()+2,0.15);
    //f_HM->SetParameter(getNHar()+3,0.75);

    hist_HM->Fit("f_HM", "0US");
    m_hist_HM = (TH1F*)hist_HM->Clone("__hist_HM");

    subResult theResult;
    vector<float> vec_value_sub;
    vector<float> vec_error_sub;

    for (int i=0; i<getNHar(); i++) {
        vec_value_sub.push_back(f_HM->GetParameter(i));
        vec_error_sub.push_back(f_HM->GetParError(i));

    }
    theResult.setCoeffSub(vec_value_sub, vec_error_sub);

    // construct stuffs for plotting
    f_show_periph = new TF1("f_show_periph", f_periph_hist2, m_dphiRangeLow, m_dphiRangeHigh, 3);
    f_show_periph->SetParameter(0, f_HM->GetParameter(getNHar()));
    f_show_periph->SetParameter(1, f_HM->GetParameter(getNHar()+1));
    f_show_periph->SetParameter(2, f_HM->GetParameter(getNHar()+2));
    //f_show_periph->SetParameter(3, f_HM->GetParameter(getNHar()+3));

    f_show_flow2 = new TF1("f_show_flow2","[0] + 2*[1]*[2]*cos(2*x)", m_dphiRangeLow, m_dphiRangeHigh);
    f_show_flow2->SetParameter(0, f_show_periph->Eval(0));
    f_show_flow2->SetParameter(1, f_HM->GetParameter(getNHar()+1) );
    f_show_flow2->SetParameter(2, f_HM->GetParameter(1)); // c2

    if (!m_fixC3) {  
        f_show_flow3 = new TF1("f_show_flow3","[0] + [1]*2*[2]*cos(3*x)", m_dphiRangeLow, m_dphiRangeHigh);
        f_show_flow3->SetParameter(0, f_show_periph->Eval(0));
        f_show_flow3->SetParameter(1, f_HM->GetParameter(getNHar()+1));
        f_show_flow3->SetParameter(2, f_HM->GetParameter(2)); // c3
    }

    h_show_periph      = (TH1F*)hist_LM     ->Clone("h_show_periph");
    h_show_periph_bulk = (TH1F*)hist_LM_bulk->Clone("h_show_periph_bulk");
    m_plotBulkRef = true;

    h_show_HM = (TH1F*)hist_LM->Clone("h_show_HM");
    h_show_HM->Reset();

    for (int ibin = 1; ibin < h_show_periph->GetNbinsX()+1; ibin++) {
        float _xx = h_show_periph->GetXaxis()->GetBinCenter(ibin);
        h_show_periph     ->SetBinContent(ibin, hist_LM->GetBinContent(ibin)      * f_HM->GetParameter(getNHar()+2)     * f_HM->GetParameter(getNHar()) + f_HM->GetParameter(getNHar()+1));
        h_show_periph_bulk->SetBinContent(ibin, hist_LM_bulk->GetBinContent(ibin) * (1-f_HM->GetParameter(getNHar()+2)) * f_HM->GetParameter(getNHar()) + f_HM->GetParameter(getNHar()+1));
        //h_show_periph     ->SetBinContent(ibin, hist_LM->GetBinContent(ibin)      * f_HM->GetParameter(getNHar()) + f_HM->GetParameter(getNHar()+1));
        //h_show_periph_bulk->SetBinContent(ibin, hist_LM_bulk->GetBinContent(ibin) * f_HM->GetParameter(getNHar()) + f_HM->GetParameter(getNHar()+1));
        //h_show_periph->SetBinContent(ibin, f_HM->Eval(_xx));
        h_show_HM->SetBinContent(ibin, f_HM->Eval(_xx));
    }

    m_h_ridge = (TH1F*) m_hist_HM->Clone("h_ridge");
    m_h_ridge->Reset();
    for (int i=1; i<m_hist_HM->GetXaxis()->GetNbins()+1; i++){
        float _xx = m_hist_HM->GetXaxis()->GetBinCenter(i);
        float _residual = m_hist_HM->GetBinContent(i) - f_show_periph->Eval(_xx);
        _residual *= m_ridge_scaleFactor;
        float _residual_error = m_hist_HM->GetBinError(i)*m_ridge_scaleFactor;
        m_h_ridge->SetBinContent(i,_residual);
        m_h_ridge->SetBinError  (i,_residual_error);
    }

    f_show_ridge2 = new TF1("f_show_ridge2","[0]*2*[1]*cos(2*x)", m_dphiRangeLow, m_dphiRangeHigh);
    f_show_ridge2->SetParameter(0, f_HM->GetParameter(getNHar()+1)*m_ridge_scaleFactor);
    f_show_ridge2->SetParameter(1, f_HM->GetParameter(1)); // c2

    if (!m_fixC3) {  
        f_show_ridge3 = new TF1("f_show_ridge3","[0]*2*[1]*cos(3*x)", m_dphiRangeLow, m_dphiRangeHigh);
        f_show_ridge3->SetParameter(0, f_HM->GetParameter(getNHar()+1)*m_ridge_scaleFactor);
        f_show_ridge3->SetParameter(1, f_HM->GetParameter(2)); // c3
    }

    std::string _formula_ridge;
    _formula_ridge = "[0]*2*([1]*cos(2*x)";
    if (!m_fixC3) {
        _formula_ridge += (" + [2]*cos(3*x)");
    }
    if (!m_fixC4) {
        _formula_ridge += (" + [3]*cos(4*x)");
    }
    _formula_ridge += (")");

    f_show_ridge = new TF1("f_show_ridge",_formula_ridge.c_str(), m_dphiRangeLow, m_dphiRangeHigh);
    f_show_ridge->SetParameter(0, f_HM->GetParameter(getNHar()+1)*m_ridge_scaleFactor);
    f_show_ridge->SetParameter(1, f_HM->GetParameter(1)); // c2
    if (!m_fixC3) {
        f_show_ridge->SetParameter(2, f_HM->GetParameter(2)); // c3
    }
    if (!m_fixC4) {
        f_show_ridge->SetParameter(3, f_HM->GetParameter(3)); // c4
    }

    return theResult;
}



subResult NonFlowSubtractor :: templateFit(TH1* hist_LM, TH1* hist_HM) {
    // procedure largely based on https://root.cern/doc/master/combinedFit_8C_source.html

    const int Npar(getNHarLM() + getNHar() + 3);
    const int Npar_LM(getNHarLM() + 1);

    const int _flowCoefIndex_begin = Npar_LM;
    const int _flowCoefIndex_end = Npar_LM + getNHar();

    subResult theResult;

    f_LM = new TF1("f_LM", f_periph, m_dphiRangeLow, m_dphiRangeHigh, Npar_LM, 1);
    f_LM->SetLineColor(kSpring-6);
    f_LM->SetLineStyle(3);

    f_HM = new TF1("f_HM", templ, m_dphiRangeLow, m_dphiRangeHigh, Npar, 1);
    f_HM->SetLineColor(2);

    if (m_debug) {
        cout << "-----" << endl;
        cout << "Npar = " << Npar << endl;
        f_HM->Print();
    }


    vector <double> parInitValues;
    parInitValues.clear();
    for (int i=0; i<Npar; i++) { 
        parInitValues.push_back(0);
    }

    //perform fourier fit first to get initial values
    subResult _init_LM = fourierFitLM(hist_LM);
    double _value_G_LM = _init_LM.getPedstalValue();
    double _error_G_LM = _init_LM.getPedstalError();
    parInitValues.at(0) = _init_LM.getPedstalValue();
    for (int i=1; i<getNHarLM()+1; i++) {
        parInitValues.at(i) = _init_LM.getCoeffRawValue(i);
    }

    // fit HM using LM paramerization
    // be carefull to the index
    subResult _init_HM = fourierFitLM(hist_HM);
    vector<float> vec_value_raw;
    vector<float> vec_error_raw;
    for (int i=getNHarLM()+1; i<getNHarLM()+getNHar()+1; i++) { 
        vec_value_raw.push_back(_init_HM.getCoeffRawValue(i-getNHarLM()));
        vec_error_raw.push_back(_init_HM.getCoeffRawError(i-getNHarLM()));
        parInitValues.at(i) = _init_HM.getCoeffRawValue(i-getNHarLM())/2.;
    }
    theResult.setCoeffRaw(vec_value_raw, vec_error_raw);
    double _value_G_HM = _init_HM.getPedstalValue();
    double _error_G_HM = _init_HM.getPedstalError();

    parInitValues.at(Npar-2) = 0.5;
    parInitValues.at(Npar-1) = 0.2;

    ROOT::Math::WrappedMultiTF1 wf_LM(*f_LM, 1);
    ROOT::Math::WrappedMultiTF1 wf_HM(*f_HM, 1);
    ROOT::Fit::DataOptions opt;
    ROOT::Fit::DataRange range;
    range.SetRange(m_dphiRangeLow,m_dphiRangeHigh);
    ROOT::Fit::BinData data_LM(opt,range);
    ROOT::Fit::BinData data_HM(opt,range);
    ROOT::Fit::FillData(data_LM, hist_LM);
    ROOT::Fit::FillData(data_HM, hist_HM);
    ROOT::Fit::Chi2Function chi2_LM(data_LM, wf_LM);
    ROOT::Fit::Chi2Function chi2_HM(data_HM, wf_HM);
    // combine the two pdf for LM and HM
    // common parameters are decleared
    GlobalChi2 globalChi2(chi2_LM, chi2_HM, getNHarLM(), getNHar());

    ROOT::Fit::Fitter fitter;
    fitter.Config().SetParamsSettings(Npar, &parInitValues.at(0));
    // print level 
    // 0 = Quiet (Default in ROOT)
    // 1 = just fit results (Default in this tool)
    // 2 = fit results with iterations
    // 3 = Verbose (additional details of iterations)
    fitter.Config().MinimizerOptions().SetPrintLevel(0);
    fitter.Config().SetMinimizer("Minuit2","Migrad");
    vector<unsigned int> minosErr_index;
    // only run minos errors for flow coefficents
    for (int i = _flowCoefIndex_begin; i < _flowCoefIndex_end; i++) {
        minosErr_index.push_back(i);
    }
    fitter.Config().SetMinosErrors( minosErr_index);

    if (m_fixC1) {
        fitter.Config().ParSettings(Npar_LM).SetValue(0);
        fitter.Config().ParSettings(Npar_LM).Fix();
    }
    if (m_fixC3 && getNHar() > 2) {
        fitter.Config().ParSettings(Npar_LM+2).SetValue(0);
        fitter.Config().ParSettings(Npar_LM+2).Fix();
    }
    if (m_fixC4 && getNHar() > 3) {
        fitter.Config().ParSettings(Npar_LM+3).SetValue(0);
        fitter.Config().ParSettings(Npar_LM+3).Fix();
    }

    // fix LM referece to initialize HM parameters
    for (int ipar = 1; ipar < getNHarLM()+1; ipar++) {
        fitter.Config().ParSettings(ipar).Fix();
    }

    for (int ipar = 0; ipar < Npar; ipar++) {
        fitter.Config().ParSettings(ipar).SetName(m_parName_HM.at(ipar).c_str());
        f_HM->SetParName(ipar, m_parName_HM.at(ipar).c_str());
    }
    
    if (m_debug) {
        cout << " ===================================================" << endl;
        cout << " Running ATLAS template fitting with the following initial values: " << endl;
        for (int ipar = 0; ipar < fitter.Config().NPar(); ipar++) {
            cout << "   parameter " << ipar << ":\t\tname = " << fitter.Config().ParSettings(ipar).Name().c_str() << ",\t\tvalue = " << fitter.Config().ParSettings(ipar).Value() << endl;
        }
        cout << endl;
    }


    //bool ROOT::Fit::Fitter::FitFCN (unsigned int npar, Function& fcn, const double* params = 0, unsigned int dataSize = 0, bool chi2fit = false)
    fitter.FitFCN(Npar, globalChi2, 0, data_LM.Size()+data_HM.Size(), true);
    
    if (!m_fixLM) {
        for (int ipar = 1; ipar < getNHarLM()+1; ipar++) {
            fitter.Config().ParSettings(ipar).Release();
        }
        // simultaneous fit with proper initial values
        fitter.FitFCN(Npar, globalChi2, 0, data_LM.Size()+data_HM.Size(), true);
    }

    ROOT::Fit::FitResult result = fitter.Result();
    if (m_debug) {
        cout << endl;
        cout << "****************************************" << endl;
        cout << "Fit results: " << endl;
        result.Print(std::cout);
        cout << endl;
        cout << "+===================================================" << endl;
        cout << "CovMatrix after fit: " << endl;
        result.PrintCovMatrix(std::cout);
        cout << endl;
    }

    double _value_F_temp = result.Value   (Npar-2);
    double _error_F_temp = result.ParError(Npar-2);

    double _value_rho_atlas = (_value_G_LM*_value_F_temp)/_value_G_HM;

    // Ignore the correlations between G_HM/G_LM (from Fourier fit) and F_temp
    // It's dominated by contribution from F_temp any way
    double _error_rho_atlas = _value_rho_atlas*sqrt(pow(_error_G_LM/_value_G_LM, 2) 
                                                  + pow(_error_G_HM/_value_G_HM, 2)
                                                  + pow(_error_F_temp/_value_F_temp,2));

    // fill the output result class
    theResult.setNHar(getNHar());
    theResult.setChi2(result.Chi2() / result.Ndf());
    theResult.setPedstalValue (_value_G_HM);
    theResult.setPedstalError (_error_G_HM);
    theResult.setRhoValue (_value_rho_atlas);
    theResult.setRhoError (_error_rho_atlas);

    vector<float> vec_value_sub;
    vector<float> vec_error_sub;
    vector<float> vec_correlation;
    vector<float> vec_minos_lower;
    vector<float> vec_minos_upper;
    for (int ihar=0; ihar<getNHar(); ihar++){
        vec_value_sub.push_back(result.Value   (Npar_LM+ihar) );
        vec_error_sub.push_back(result.ParError(Npar_LM+ihar) );
        // F-coefficient correlation is used for rho-coefficient correlation
        // used for improved fit error calculation
        vec_correlation.push_back(result.Correlation(Npar-2, Npar_LM+ihar)*(_value_G_LM/_value_G_HM) );

        vec_minos_lower.push_back(result.LowerError(Npar_LM+ihar));
        vec_minos_upper.push_back(result.UpperError(Npar_LM+ihar));
    }

    theResult.setCoeffSub(vec_value_sub, vec_error_sub);
    theResult.setRhoCorrelation(vec_correlation);
    theResult.setCoeffSubMinos(vec_minos_lower, vec_minos_upper);

    if (m_debug) {
        cout << endl;
        cout << endl;
        cout << "+===================================================" << endl;
        cout << " Check relative errors" << endl;
        cout << "   G_LM\tvalue = " << std::scientific << _value_G_LM 
             << "\terror = " << _error_G_LM 
             << "\trel_error = " << _error_G_LM/ _value_G_LM 
             << endl;

        cout << "   G_HM\tvalue = " << std::scientific << _value_G_HM 
             << "\terror = " << _error_G_HM 
             << "\trel_error = " << _error_G_HM/ _value_G_HM 
             << endl;

        cout << "   F_tm\tvalue = " << std::scientific << _value_F_temp 
             << "\terror = " << _error_F_temp 
             << "\trel_error = " << _error_F_temp/ _value_F_temp 
             << endl;

        cout << "   rho\tvalue = " << std::scientific << _value_rho_atlas 
             << "\terror = " << _error_rho_atlas 
             << "\trel_error = " << _error_rho_atlas/_value_rho_atlas 
             << endl << endl;

        cout << "+===================================================" << endl;
        cout << " Compare different error estimations:" << endl;
        // difference should indiate the non-linearity of the problem
        for (int i = _flowCoefIndex_begin; i < _flowCoefIndex_end; i++) {
            cout << "parameter: " << fitter.Config().ParSettings(i).Name() << endl;
            cout << "\tHESSE error = " << result.ParError(i) << endl;
            cout << "\tMinos error lower = " << result.LowerError(i) << ", upper = " << result.UpperError(i) << endl;
        }
    }

    // construct stuffs for plotting
    // -------------------------------------------------------------
    // Prepare histograms and function for plotting 
    // uncorrected cn parameters are used for visulization
    string _f_periph_forPlot = "[0] * ( 1 + 2*(";
    for (int ihar=0; ihar<getNHarLM(); ihar++){
        _f_periph_forPlot += Form("[%d]*cos(%d*x)", ihar+1, ihar+1);
        if (ihar != getNHarLM()-1) _f_periph_forPlot += (" + ");
    }
    _f_periph_forPlot += Form(") )*[%d] + [%d]", getNHarLM()+1, getNHarLM()+2);
    if (m_debug) cout << "Scaled peripheral function: " << _f_periph_forPlot.c_str() << endl;
    f_show_periph = new TF1("f_show_periph", _f_periph_forPlot.c_str(), m_dphiRangeLow, m_dphiRangeHigh);

    for (int ipar = 0; ipar < getNHarLM()+1; ipar++) {
        f_show_periph->SetParameter(ipar, result.Value(ipar));
    }

    f_show_periph->SetParameter(getNHarLM()+1, result.Value(Npar - 2));
    f_show_periph->SetParameter(getNHarLM()+2, result.Value(Npar - 1));

    f_show_flow2 = new TF1("f_show_flow2","[0] + [1]*(1+2*[2]*cos(2*x))", m_dphiRangeLow, m_dphiRangeHigh);
    f_show_flow2->SetParameter(0, f_LM->Eval(0)*result.Value(Npar - 2));
    f_show_flow2->SetParameter(1, result.Value(Npar - 1));
    f_show_flow2->SetParameter(2, result.Value(getNHarLM()+2)); // c2

    if (!m_fixC3) {  
        f_show_flow3 = new TF1("f_show_flow3","[0] + [1]*(1+2*[2]*cos(3*x) )", m_dphiRangeLow, m_dphiRangeHigh);
        f_show_flow3->SetParameter(0, f_LM->Eval(0)*result.Value(Npar - 2));
        f_show_flow3->SetParameter(1, result.Value(Npar - 1));
        f_show_flow3->SetParameter(2, result.Value(getNHarLM()+3)); // c3
    }

    m_hist_HM = (TH1F*)hist_HM->Clone("__hist_HM");
    m_hist_LM = (TH1F*)hist_LM->Clone("__hist_LM");

    m_h_ridge = (TH1F*) m_hist_HM->Clone("h_ridge");
    m_h_ridge->Reset();
    for (int i=1; i<m_hist_HM->GetXaxis()->GetNbins()+1; i++){
        float _residual = m_hist_HM->GetBinContent(i) - f_show_periph->Eval(m_hist_HM->GetBinCenter(i)) + result.Value(Npar - 1); // add back the pedestal
        _residual *= m_ridge_scaleFactor;
        float _residual_error = _residual * ( m_hist_HM->GetBinError(i) / m_hist_HM->GetBinContent(i) ) * m_ridge_scaleFactor;
        //float _residual_error = m_hist_HM->GetBinError(i)*m_ridge_scaleFactor;
        m_h_ridge->SetBinContent(i,_residual);
        m_h_ridge->SetBinError  (i,_residual_error);
    }

    f_show_ridge2 = new TF1("f_show_ridge2","[0]*(1+2*[1]*cos(2*x))", m_dphiRangeLow, m_dphiRangeHigh);
    f_show_ridge2->SetParameter(0, result.Value(Npar - 1)*m_ridge_scaleFactor);
    f_show_ridge2->SetParameter(1, result.Value(getNHarLM()+2)); // c2

    if (!m_fixC3) {  
        f_show_ridge3 = new TF1("f_show_ridge3","[0]*(1+2*[1]*cos(3*x))", m_dphiRangeLow, m_dphiRangeHigh);
        f_show_ridge3->SetParameter(0, result.Value(Npar - 1)*m_ridge_scaleFactor);
        f_show_ridge3->SetParameter(1, result.Value(getNHarLM()+3)); // c3
    }

    std::string _formula_ridge;
    _formula_ridge = "[0]*(1+ 2*[1]*cos(2*x)";
    if (!m_fixC3) {
        _formula_ridge += (" + 2*[2]*cos(3*x)");
    }
    if (!m_fixC4) {
        _formula_ridge += (" + 2*[3]*cos(4*x)");
    }
    _formula_ridge += (")");
    f_show_ridge = new TF1("f_show_ridge", _formula_ridge.c_str(), m_dphiRangeLow, m_dphiRangeHigh);
    f_show_ridge->SetParameter(0, result.Value(Npar - 1)*m_ridge_scaleFactor);
    f_show_ridge->SetParameter(1, result.Value(getNHarLM()+2)); // c3
    if (!m_fixC3) {
        f_show_ridge->SetParameter(2, result.Value(getNHarLM()+3)); // c3
    } else {
        f_show_ridge->SetParameter(2, 0); // c3
    }
    if (!m_fixC4) {
        f_show_ridge->SetParameter(3, result.Value(getNHarLM()+4)); // c3
    } else {
        f_show_ridge->SetParameter(3, 0); // c3
    }


    // for parameterization bias test
    // current only bias on v22 has been tested
    if (m_pesudoVaryScale != 0) {
        TF1* f_pesudo_HM = new TF1("f_pesudo_HM", templ, m_dphiRangeLow, m_dphiRangeHigh, Npar);

        for (int ipar = 0; ipar<f_HM->GetNpar(); ipar++) {
            float _value = f_HM->GetParameter(ipar);
            float _error = f_HM->GetParError(ipar);
            f_pesudo_HM->SetParameter(ipar, _value);
            f_pesudo_HM->SetParError (ipar, _error);

            if (ipar == Npar_LM+1) f_pesudo_HM->SetParameter(ipar, _value*(1.+m_pesudoVaryScale));
        }
        f_pesudo_HM->Update();
        m_hist_pesudo_HM = (TH1F*) hist_HM->Clone("m_hist_pesudo_HM");
        m_hist_pesudo_LM = (TH1F*) hist_LM->Clone("m_hist_pesudo_LM");

        for (int ibin=1; ibin < m_hist_pesudo_HM->GetNbinsX()+1; ibin++) {
            float _binCenter = m_hist_pesudo_HM->GetXaxis()->GetBinCenter(ibin);
            m_hist_pesudo_HM->SetBinContent(ibin, f_pesudo_HM->Eval(_binCenter));
        }
        delete f_pesudo_HM;
    }


    return theResult;
}



subResult NonFlowSubtractor::templateFit (TH1* hist_LM,  TH1* hist_HM, TH1* hist_LM2) {

    subResult result_hm  = templateFit(hist_LM, hist_HM);

    // for keeping the plot stuff to be the same as running unimproved fit
    // correction fit is done in different object 
    NonFlowSubtractor subLM2;
    subLM2.setNHar(NonFlowSubtractor::m_nhar_HM);
    subLM2.setNHarLM(NonFlowSubtractor::m_nhar_LM);
    subLM2.init();
    if (m_debug) {
        cout << endl;
        cout << " ===================================================" << endl;
        cout << " Running default template fitting to LM2 for obtaining correction" << endl;
    }
    subResult result_lm2 = subLM2.templateFit(hist_LM, hist_LM2);

    vector<float> vec_subImp_value;
    vector<float> vec_subImp_error;
    for (int ihar=0; ihar<getNHar(); ihar++){
        float _rho      = result_hm .getRhoValue();
        float _rho_error= result_hm .getRhoError();
        float _coeff_hm = result_hm .getCoeffSubValue(ihar+1);
        float _error_hm = result_hm .getCoeffSubError(ihar+1);
        float _coeff_lm = result_lm2.getCoeffSubValue(ihar+1);
        float _error_lm = result_lm2.getCoeffSubError(ihar+1);

        float _correlaiton = result_hm .getRhoCorrelation(ihar+1);

        //ignore error on rho
        float _coeff_corrected = _coeff_hm - _rho*(_coeff_hm - _coeff_lm);
        float _error = sqrt( pow((1-_rho)*_error_hm,2) 
                           + pow( _rho*_error_lm,2)
                           + pow( (_coeff_hm-_coeff_lm)*_rho_error,2)
                           - 2*(1-_rho)*(_coeff_hm-_coeff_lm)*(_rho_error*_error_hm*_correlaiton)
                           // ignore correlation between C_LM and rho/C_HM
                           ); // to be improved

        //special case for error estimation
        if (_coeff_hm == _coeff_lm) _error = _error_hm; // using HM as LM2, error should stay the same as before improving
        if (_coeff_hm == 0) _error = 0; // protection for n=1 when c1 == 0;

        vec_subImp_value.push_back(_coeff_corrected);
        vec_subImp_error.push_back(_error);
    }
    result_hm.setCoeffSubImp(vec_subImp_value, vec_subImp_error);

    if (m_debug) {
        cout << " ===================================================" << endl;
        cout << std::setprecision(5) << " Improved template fitting results: " << endl;
        cout << "  -unsubtracted:\t v22 value = " << result_hm.getV22RawValue() 
             << ",\terror = " << result_hm.getV22RawError() 
             << ",\trel_error = " << result_hm.getV22RawError()/result_hm.getV22RawValue() << endl;

        cout << "  -default method:\t v22 value = " << result_hm.getV22SubValue() 
             << ",\terror = " << result_hm.getV22SubError() 
             << ",\trel_error = " << result_hm.getV22SubError()/result_hm.getV22SubValue() << endl;

        cout << "  -improved method:\t v22 value = " << result_hm.getV22SubImpValue() 
             << ",\terror = " << result_hm.getV22SubImpError() 
             << ",\trel_error = " << result_hm.getV22SubImpError()/result_hm.getV22SubImpValue() << endl << endl;
    }

    return result_hm;
}



subResult NonFlowSubtractor::templateFit (TH1* hist_LM,  TH1* hist_HM, float cn_LM, float cn_LM_error) {

    subResult result_hm  = templateFit(hist_LM, hist_HM);
    if (m_debug) {
        cout << endl;
        cout << " ===================================================" << endl;
        cout << " Running default template fitting with external cn_LM input" << endl;
    }

    vector<float> vec_subImp_value;
    vector<float> vec_subImp_error;
    for (int ihar=0; ihar<getNHar(); ihar++){
        float _rho      = result_hm .getRhoValue();
        float _rho_error= result_hm .getRhoError();
        float _coeff_hm = result_hm .getCoeffSubValue(ihar+1);
        float _error_hm = result_hm .getCoeffSubError(ihar+1);
        float _coeff_lm = cn_LM;
        float _error_lm = cn_LM_error;

        float _correlaiton = result_hm .getRhoCorrelation(ihar+1);

        //ignore error on rho
        float _coeff_corrected = _coeff_hm - _rho*(_coeff_hm - _coeff_lm);
        float _error = sqrt( pow((1-_rho)*_error_hm,2) 
                           + pow( _rho*_error_lm,2)
                           + pow( (_coeff_hm-_coeff_lm)*_rho_error,2)
                           - 2*(1-_rho)*(_coeff_hm-_coeff_lm)*(_rho_error*_error_hm*_correlaiton)
                           // ignore correlation between C_LM and rho/C_HM
                           ); // to be improved

        //special case for error estimation
        if (_coeff_hm == _coeff_lm) _error = _error_hm; // using HM as LM2, error should stay the same as before improving
        if (_coeff_hm == 0) _error = 0; // protection for n=1 when c1 == 0;

        vec_subImp_value.push_back(_coeff_corrected);
        vec_subImp_error.push_back(_error);
    }
    result_hm.setCoeffSubImp(vec_subImp_value, vec_subImp_error);

    if (m_debug) {
        cout << " ===================================================" << endl;
        cout << std::setprecision(5) << " Improved template fitting results: " << endl;
        cout << "  -unsubtracted:\t v22 value = " << result_hm.getV22RawValue() 
             << ",\terror = " << result_hm.getV22RawError() 
             << ",\trel_error = " << result_hm.getV22RawError()/result_hm.getV22RawValue() << endl;

        cout << "  -default method:\t v22 value = " << result_hm.getV22SubValue() 
             << ",\terror = " << result_hm.getV22SubError() 
             << ",\trel_error = " << result_hm.getV22SubError()/result_hm.getV22SubValue() << endl;

        cout << "  -improved method:\t v22 value = " << result_hm.getV22SubImpValue() 
             << ",\terror = " << result_hm.getV22SubImpError() 
             << ",\trel_error = " << result_hm.getV22SubImpError()/result_hm.getV22SubImpValue() << endl << endl;
    }

    return result_hm;
}



subResult NonFlowSubtractor::fourierFit (TH1* hist, TH1* hist_LM) {

    TF1* f_fourier = new TF1("f_fourier", _formula_LM.c_str(), m_dphiRangeLow, m_dphiRangeHigh);

    hist->Fit("f_fourier", "0");

    vector<float> vec_value;
    vector<float> vec_error;
    for (int ihar=0; ihar<getNHar(); ihar++){
        vec_value.push_back(f_fourier->GetParameter(ihar+1));
        vec_error.push_back(f_fourier->GetParError (ihar+1));
    }

    subResult theResult;
    theResult.setNHar(getNHar());
    theResult.setCoeffRaw(vec_value, vec_error);
    theResult.setPedstalValue (f_fourier->GetParameter(0));
    theResult.setPedstalError (f_fourier->GetParError(0));

    // if hist_LM is not zero, rho for second order is calculated 
    // used for pythia study
    if (hist_LM) {
        TF1* f_fourier_LM = new TF1("f_fourier_LM", _formula.c_str(), m_dphiRangeLow, m_dphiRangeHigh);
        hist_LM->Fit("f_fourier_LM", "0");
        float _rho = f_fourier->GetParameter(2) / f_fourier_LM->GetParameter(2);
        // driving by LM
        float _rho_error = fabs(_rho)*fabs(f_fourier_LM->GetParError(2)/f_fourier_LM->GetParameter(2));
        theResult.setRhoValue(_rho);
        theResult.setRhoError(_rho_error);
        delete f_fourier_LM; f_fourier_LM = 0;
    }

    delete f_fourier; f_fourier = 0;

    return theResult;
}


subResult NonFlowSubtractor::fourierFitLM (TH1* hist) {

    TF1* f_fourier = new TF1("f_fourier", _formula_LM.c_str(), m_dphiRangeLow, m_dphiRangeHigh);

    // Sometimes there are problmes with the simple fourier fit due to not proper initialization
    f_fourier->SetParameter(0, hist->Integral());
    for (int ihar=1; ihar<getNHar(); ihar++){
        float _cosx = 0;
        float _count = 0;
        for (int ibin = 1; ibin < hist->GetNbinsX()+1; ibin++) {
            _cosx += TMath::Cos(ihar*hist->GetBinCenter(ibin)) * hist->GetBinContent(ibin);
            _count += hist->GetBinContent(ibin);
        }
        _cosx /= _count;
        f_fourier->SetParameter(ihar, _cosx);
    }
    
    string option = "0Q";
    if (m_debug) option = "0V";
    hist->Fit("f_fourier", Form("%s",option.c_str()) );

    vector<float> vec_value;
    vector<float> vec_error;
    for (int ihar=0; ihar<getNHarLM(); ihar++){
        vec_value.push_back(f_fourier->GetParameter(ihar+1));
        vec_error.push_back(f_fourier->GetParError (ihar+1));
    }

    subResult theResult;
    theResult.setNHar(getNHarLM());
    theResult.setCoeffRaw(vec_value, vec_error);
    theResult.setPedstalValue (f_fourier->GetParameter(0));
    theResult.setPedstalError (f_fourier->GetParError(0));

    delete f_fourier; f_fourier = 0;

    return theResult;
}



subResult NonFlowSubtractor::periphSub  ( TH1* h_sr_lm, TH1* h_lr_lm, 
                                          TH1* h_sr_hm, TH1* h_lr_hm) {
    // SR - LR for low multiplicity
    TH1F* h_sub_lm = (TH1F*)h_sr_lm->Clone("h_sub_lm");
    h_sub_lm->Add(h_lr_lm, -1);
    
    // SR - LR for high multiplicity
    TH1F* h_sub_hm = (TH1F*)h_sr_hm->Clone("h_sub_hm");
    h_sub_hm->Add(h_lr_hm, -1);

    // should be very careful here. 
    // if one follows the CMS paper to make PTY (2D ratio then projection to 1D dphi correlation), ZYAM is not necessary
    // on the other hand, if you follows the ATPAS paper, ZYAM must be applied
    if (m_applyZYAM) {
        ZYAM(h_sub_lm);
        ZYAM(h_sub_hm);
    }

    // hard coded here for now.
    // Could make them as data members
    int dphiBin_low  = h_sub_lm->FindBin(-1.2);
    int dphiBin_high = h_sub_lm->FindBin(1.2);
    
    double _factor_hm_error;
    float factor_hm = h_sub_hm->IntegralAndError(dphiBin_low, dphiBin_high, _factor_hm_error);
    double _factor_lm_error;
    float factor_lm = h_sub_lm->IntegralAndError(dphiBin_low, dphiBin_high, _factor_lm_error);
    
    // fit long-range correlation
    TF1* f_lm = new TF1("f_lm", _formula.c_str(), m_dphiRangeLow, m_dphiRangeHigh);
    TF1* f_hm = new TF1("f_hm", _formula.c_str(), m_dphiRangeLow, m_dphiRangeHigh);

    h_lr_lm->Fit("f_lm","0");
    h_lr_hm->Fit("f_hm","0");

    float G_lm_count = f_lm->GetParameter(0);
    float G_hm_count = f_hm->GetParameter(0);
    
    float rhoCMS = (G_lm_count * factor_hm) / (G_hm_count * factor_lm);
    // assuming error is driving by factor_hm and factor_lm
    float rhoCMS_error = rhoCMS * sqrt( pow(_factor_lm_error/factor_lm, 2) + pow(_factor_hm_error/factor_hm, 2));

    vector<float> vec_raw_value;
    vector<float> vec_raw_error;
    vector<float> vec_sub_value;
    vector<float> vec_sub_error;
    for (int ihar=0; ihar<getNHar(); ihar++){
        vec_raw_value.push_back(f_hm->GetParameter(ihar+1));
        vec_raw_error.push_back(f_hm->GetParError (ihar+1));

        float _coeff_sub = f_hm->GetParameter(ihar+1) - f_lm->GetParameter(ihar+1)*rhoCMS;

        vec_sub_value.push_back(_coeff_sub);
        // need to improve the error estimation
        vec_sub_error.push_back(f_hm->GetParError(ihar+1));
    }
    subResult theResult;
    theResult.setNHar(getNHar());
    theResult.setCoeffRaw(vec_raw_value, vec_raw_error);
    theResult.setCoeffSub(vec_sub_value, vec_sub_error);
    theResult.setPedstalValue (f_hm->GetParameter(0));
    theResult.setPedstalError (f_hm->GetParError(0));

    theResult.setRhoValue(rhoCMS);
    theResult.setRhoError(rhoCMS_error);

    delete f_hm; f_hm = 0;
    delete f_lm; f_lm = 0;
    return theResult;
}



subResult NonFlowSubtractor::periphSub  (TH1* h_sr_lm,  TH1* h_lr_lm, 
                                         TH1* h_sr_hm,  TH1* h_lr_hm, 
                                         TH1* h_sr_lm2, TH1* h_lr_lm2) {

    subResult result_hm  = periphSub(h_sr_lm, h_lr_lm, h_sr_hm,  h_lr_hm);

    NonFlowSubtractor subLM2;
    subLM2.setNHar(NonFlowSubtractor::m_nhar_HM);
    subLM2.init();
    subResult result_lm2 = subLM2.periphSub(h_sr_lm, h_lr_lm, h_sr_lm2, h_lr_lm2);

    // since rho is calculated from short range correlation, correlation between rho and flow coefficient is ignored (to be improved)
    vector<float> vec_subImp_value;
    vector<float> vec_subImp_error;
    for (int ihar=0; ihar<getNHar(); ihar++){
        float _rho      = result_hm .getRhoValue();
        float _rho_error= result_hm .getRhoError();
        float _coeff_hm = result_hm .getCoeffSubValue(ihar+1);
        float _error_hm = result_hm .getCoeffSubError(ihar+1);
        float _coeff_lm = result_lm2.getCoeffSubValue(ihar+1);
        float _error_lm = result_lm2.getCoeffSubError(ihar+1);

        float _coeff_corrected = _coeff_hm + _rho*(_coeff_lm);
        float _error = sqrt(pow(_error_hm, 2) + pow(_rho*_error_lm, 2) + pow(_coeff_lm*_rho_error, 2)); 
        vec_subImp_value.push_back(_coeff_corrected);
        vec_subImp_error.push_back(_error);
    }

    //update result_hm with imporoved fit results
    result_hm.setCoeffSubImp(vec_subImp_value, vec_subImp_error);

    if (m_debug) {
        cout << " ===================================================" << endl;
        cout << std::setprecision(5) << " Improved peripheral sbutraction results: " << endl;
        cout << "  -unsubtracted:\t v22 value = " << result_hm.getV22RawValue() 
             << ",\terror = " << result_hm.getV22RawError() 
             << ",\trel_error = " << result_hm.getV22RawError()/result_hm.getV22RawValue() << endl;

        cout << "  -default method:\t v22 value = " << result_hm.getV22SubValue() 
             << ",\terror = " << result_hm.getV22SubError() 
             << ",\trel_error = " << result_hm.getV22SubError()/result_hm.getV22SubValue() << endl;

        cout << "  -improved method:\t v22 value = " << result_hm.getV22SubImpValue() 
             << ",\terror = " << result_hm.getV22SubImpError() 
             << ",\trel_error = " << result_hm.getV22SubImpError()/result_hm.getV22SubImpValue() << endl << endl;
    }

    return result_hm;
}




subResult NonFlowSubtractor::periphSub  (TH1* h_sr_lm,  TH1* h_lr_lm, 
                                         TH1* h_sr_hm,  TH1* h_lr_hm, 
                                         float cn_LM,  float cn_LM_error) {

    subResult result_hm  = periphSub(h_sr_lm, h_lr_lm, h_sr_hm,  h_lr_hm);

    vector<float> vec_subImp_value;
    vector<float> vec_subImp_error;
    for (int ihar=0; ihar<getNHar(); ihar++){
        float _rho      = result_hm .getRhoValue();
        float _rho_error= result_hm .getRhoError();
        float _coeff_hm = result_hm .getCoeffSubValue(ihar+1);
        float _error_hm = result_hm .getCoeffSubError(ihar+1);
        float _coeff_lm = cn_LM;
        float _error_lm = cn_LM_error;

        float _coeff_corrected = _coeff_hm + _rho*(_coeff_lm);
        float _error = sqrt(pow(_error_hm, 2) + pow(_rho*_error_lm, 2) + pow(_coeff_lm*_rho_error, 2)); 
        vec_subImp_value.push_back(_coeff_corrected);
        vec_subImp_error.push_back(_error);
    }

    //update result_hm with imporoved fit results
    result_hm.setCoeffSubImp(vec_subImp_value, vec_subImp_error);

    if (m_debug) {
        cout << " ===================================================" << endl;
        cout << std::setprecision(5) << " Improved peripheral sbutraction results: " << endl;
        cout << "  -unsubtracted:\t v22 value = " << result_hm.getV22RawValue() 
             << ",\terror = " << result_hm.getV22RawError() 
             << ",\trel_error = " << result_hm.getV22RawError()/result_hm.getV22RawValue() << endl;

        cout << "  -default method:\t v22 value = " << result_hm.getV22SubValue() 
             << ",\terror = " << result_hm.getV22SubError() 
             << ",\trel_error = " << result_hm.getV22SubError()/result_hm.getV22SubValue() << endl;

        cout << "  -improved method:\t v22 value = " << result_hm.getV22SubImpValue() 
             << ",\terror = " << result_hm.getV22SubImpError() 
             << ",\trel_error = " << result_hm.getV22SubImpError()/result_hm.getV22SubImpValue() << endl << endl;
    }

    return result_hm;
}



// plot style for atlas template fitting using histogram as templates
bool NonFlowSubtractor :: plotAtlasHistHM (TCanvas* theCanvas) {
    if (!theCanvas) return false;
    theCanvas->cd();

    m_hist_HM->SetMarkerSize(0.8);
    m_hist_HM->SetXTitle("#Delta#it{#phi}");
    m_hist_HM->SetYTitle("#it{Y} (#Delta#it{#phi})");
    m_hist_HM->SetLineWidth(1);
    m_hist_HM->Draw("EX0SAME");
    h_show_HM->SetLineColor(2);
    h_show_HM->Draw("HISTSAME");
    f_show_flow2->SetLineColor(kBlue);
    f_show_flow2->SetLineStyle(2);
    f_show_flow2->SetLineWidth(1);
    f_show_flow2->Draw("same");

    h_show_periph->SetMarkerStyle(24);
    h_show_periph->SetMarkerSize(0.8);
    h_show_periph->Draw("psame");

    return true;
}




bool NonFlowSubtractor :: plotAtlasHM (TPad* thePad) {
    if (!thePad) return false;
    thePad->cd();

    m_hist_HM->SetMarkerSize(0.8);
    m_hist_HM->SetXTitle("#Delta#it{#phi}");
    m_hist_HM->SetYTitle("#it{Y} (#Delta#it{#phi})");
    m_hist_HM->SetLineWidth(1);
    f_HM->SetLineWidth(2);
    m_hist_HM->GetListOfFunctions()->Add(f_HM);
    m_hist_HM->Draw("EX0SAME");
    f_show_flow2->SetLineColor(kBlue);
    f_show_flow2->SetLineStyle(2);
    f_show_flow2->SetLineWidth(1);
    f_show_flow2->Draw("same");
    f_show_periph->SetLineWidth(2);
    f_show_periph->SetLineColor(kSpring+4);
    f_show_periph->SetLineStyle(2);
    f_show_periph->Draw("same");

    return true;
}


bool NonFlowSubtractor :: plotAtlasHistSubHM (TCanvas* theCanvas) {
    if (!theCanvas) return false;
    TH1F* h_pull = (TH1F*) m_hist_HM->Clone("h_Pull");
    theCanvas->cd();

    TPad *pad1 = new TPad("pad1","top pad",0,0.36,1,1);
    pad1->SetTopMargin(0.07);
    pad1->SetRightMargin(0.02);
    pad1->SetBottomMargin(0.);
    pad1->SetLeftMargin(0.13);
    theCanvas->cd();
    pad1->Draw();
    TPad *pad2 = new TPad("pad2","bottom pad",0,0,1,0.36);
    pad2->SetTopMargin(0);
    pad2->SetRightMargin(0.02);
    pad2->SetBottomMargin(0.2);
    pad2->SetLeftMargin(0.13);
    theCanvas->cd();
    pad2->Draw();

    pad1->cd();
    int MaxBin = m_hist_HM->GetMaximumBin();
    int MinBin = m_hist_HM->GetMinimumBin();
    float distance = m_hist_HM->GetBinContent(MaxBin) - m_hist_HM->GetBinContent(MinBin);
    if (m_plotBulkRef) {
        distance = m_hist_HM->GetBinContent(MaxBin) - h_show_periph->GetBinContent(h_show_periph->GetMinimumBin());
    }
    distance /= 4;
    double Max = m_hist_HM->GetBinContent(MaxBin) + distance;
    double Min = m_hist_HM->GetBinContent(MinBin) - distance/2.5;
    if (m_plotBulkRef) {
        Min = h_show_periph->GetBinContent(h_show_periph->GetMinimumBin()) - distance/2.5;
    }

    m_hist_HM->SetMarkerSize(1.3);
    m_hist_HM->SetYTitle("#it{Y}(#Delta#it{#phi})");
    m_hist_HM->GetYaxis()->SetRangeUser(Min, Max);
    m_hist_HM->GetYaxis()->SetNdivisions(508,kTRUE);
    m_hist_HM->GetXaxis()->SetNdivisions(509,kTRUE);
    m_hist_HM->GetYaxis()->SetLabelFont(43);
    m_hist_HM->GetYaxis()->SetLabelSize(25);
    m_hist_HM->GetYaxis()->SetTitleFont(43);
    m_hist_HM->GetYaxis()->SetTitleSize(27);
    m_hist_HM->GetYaxis()->SetTitleOffset(1.95);
    m_hist_HM->GetXaxis()->SetTickLength(0.067);
    m_hist_HM->Draw("EX0SAME");
    h_show_HM->SetLineColor(2);
    h_show_HM->SetLineWidth(3);
    h_show_HM->Draw("HISTSAME");
    f_show_flow2->SetLineColor(kBlue);
    f_show_flow2->SetLineStyle(2);
    f_show_flow2->SetLineWidth(3);
    f_show_flow2->Draw("same");
    if (!m_fixC3) {
        f_show_flow3->SetLineColor(kOrange+1);
        f_show_flow3->SetLineStyle(3);
        f_show_flow3->SetLineWidth(3);
        f_show_flow3->Draw("same");
    }
    h_show_periph->SetMarkerStyle(24);
    h_show_periph->SetMarkerSize(1.3);
    h_show_periph->Draw("psame");
    if (m_plotBulkRef) {
        h_show_periph_bulk->SetMarkerStyle(26);
        h_show_periph_bulk->SetMarkerColor(4);
        h_show_periph_bulk->SetMarkerSize(1.3);
        h_show_periph_bulk->Draw("psame");
    }

    plotMarkerLineText(0.22, 0.78, 1.0, 1, 20, 1, 1,"#bf{HM Data}", 0.05, true);
    plotMarkerLineText(0.22, 0.71, 0, 2, 1, 2, 1,"#bf{Fit}", 0.05);
    if (m_plotBulkRef) {
        plotMarkerLineText(0.22, 0.64, 1.0, 1, 24, 1, 1,"#bf{#it{G} + #it{Fp}#it{Y}^{LM}_{1}}", 0.05, true);
        plotMarkerLineText(0.22, 0.57, 1.0, 4, 26, 4, 1,"#bf{#it{G} + #it{F(1-p)}#it{Y}^{LM}_{2}}", 0.05, true);
    } else {
        plotMarkerLineText(0.22, 0.64, 1.0, 1, 24, 1, 1,"#bf{#it{G} + #it{F}#it{Y}^{LM}}", 0.05, true);
    }
    plotMarkerLineText(              0.42,0.78, 0, 4, 0, 4, 2,"#bf{#it{Y}_{2}^{ridge} + #it{F}#it{Y}^{LM}(0)}",0.05);
    if (!m_fixC3) plotMarkerLineText(0.42,0.71, 0, kOrange+1, 0, kOrange+1, 3,"#bf{#it{Y}_{3}^{ridge} + #it{F}#it{Y}^{LM}(0)}", 0.05);

    /*
    float _chi2 = 0;
    for (int i=1; i<h_pull->GetXaxis()->GetNbins()+1; i++){
        float _residual = m_hist_HM->GetBinContent(i) - f_HM->Eval(m_hist_HM->GetBinCenter(i));
        float _pull = _residual / m_hist_HM->GetBinError(i);
        h_pull->SetBinContent(i,_pull);
        h_pull->SetBinError(i,m_hist_HM->GetBinError(i)/m_hist_HM->GetBinError(i));
        _chi2 += pow(_residual/m_hist_HM->GetBinError(i),2);
    }
    _chi2 /= h_pull->GetXaxis()->GetNbins()-3;
    m_h_pull = (TH1F*) h_pull->Clone("m_h_Pull");
    */

    pad2->cd();
    int MaxBin_sub = m_h_ridge->GetMaximumBin();
    int MinBin_sub = m_h_ridge->GetMinimumBin();
    double Max_sub = TMath::Max( fabs(m_h_ridge->GetBinContent(MaxBin_sub)), fabs(m_h_ridge->GetBinContent(MinBin_sub)) )*2.;
    double Min_sub = -1*Max_sub;
    m_h_ridge->SetMarkerSize(1.3);
    //m_h_ridge->SetLineWidth(1);
    m_h_ridge->GetYaxis()->SetRangeUser(Min_sub, Max_sub);
    m_h_ridge->SetXTitle("#Delta#it{#phi}");
    if (m_ridge_scaleFactor != 1.0) {
        m_h_ridge->SetYTitle(Form("(#it{Y}(#Delta#it{#phi}) - #it{G} - #it{F}#it{Y}^{LM}(#Delta#it{#phi}))#times%.0f", m_ridge_scaleFactor));
    } else {
        m_h_ridge->SetYTitle(Form("#it{Y}(#Delta#it{#phi}) - #it{G} - #it{F}#it{Y}^{LM}(#Delta#it{#phi})"));
    }
    m_h_ridge->GetYaxis()->SetNdivisions(506,kTRUE);
    m_h_ridge->GetYaxis()->SetLabelFont(43);
    m_h_ridge->GetYaxis()->SetLabelSize(25);
    m_h_ridge->GetYaxis()->SetTitleFont(43);
    m_h_ridge->GetYaxis()->SetTitleSize(27);
    m_h_ridge->GetYaxis()->SetTitleOffset(1.95);
    m_h_ridge->GetYaxis()->CenterTitle(kTRUE);

    m_h_ridge->GetXaxis()->SetTitleFont(43);
    m_h_ridge->GetXaxis()->SetTitleSize(27);
    m_h_ridge->GetXaxis()->SetLabelFont(43);
    m_h_ridge->GetXaxis()->SetLabelSize(25);
    m_h_ridge->GetXaxis()->SetTitleOffset(1.9);
    m_h_ridge->GetXaxis()->SetTickLength(0.08);
    m_h_ridge->Draw("EX0SAME");

    f_show_ridge->SetLineStyle(1);
    f_show_ridge->SetLineWidth(3);
    f_show_ridge->SetLineColor(2);
    f_show_ridge->Draw("same");

    f_show_ridge2->SetLineColor(kBlue);
    f_show_ridge2->SetLineWidth(3);
    f_show_ridge2->SetLineStyle(2);
    f_show_ridge2->Draw("same");
    if (!m_fixC3) {
        f_show_ridge3->SetLineColor(kOrange+1);
        f_show_ridge3->SetLineStyle(3);
        f_show_ridge3->SetLineWidth(3);
        f_show_ridge3->Draw("same");
    }
    //plotMarkerLineText(0.50, 0.95, 1.0,1, 20, 1,1,"Data", 0.075, true);
    //plotMarkerLineText(0.60, 0.88, 0, 2, 1, 2, 1,"#it{Y}^{ridge}", 0.075);
    //plotMarkerLineText(              0.74,0.95, 0, 4, 0, 4, 2,"#it{Y}_{2}^{ridge}",0.075);
    //if (!m_fixC3) plotMarkerLineText(0.74,0.88, 0, kOrange+1, 0, kOrange+1, 3,"#it{Y}_{3}^{ridge}", 0.075);

    pad1->cd();

    return true;
}


void NonFlowSubtractor :: plotAtlasHMLabel (TPad* thePad) {
    thePad->cd();
    plotMarkerLineText(0.25, 0.88, 1.2, 1, 20, 1,1,"HM Data", 0.08, true);
    plotMarkerLineText(0.25, 0.78, 0,   2, 1, 2, 1,"Fit", 0.08);
    plotMarkerLineText(0.25, 0.68, 0, kSpring+4, 0, kSpring+4, 2,"#it{G} + #it{F}#it{Y}^{LM}", 0.08);
    plotMarkerLineText(0.25, 0.58, 0, 4, 0, 4, 2,"#it{Y}_{2}^{ridge} + #it{F}#it{Y}^{LM}(0)",0.08);
}




bool NonFlowSubtractor :: plotAtlasSubHM (TCanvas* theCanvas) {
    if (!theCanvas) return false;
    TH1F* h_pull = (TH1F*) m_hist_HM->Clone("h_Pull");
    theCanvas->cd();

    TPad *pad1 = new TPad("pad1","top pad",0,0.5,1,1);
    pad1->SetTopMargin(0.07);
    pad1->SetBottomMargin(0.);
    theCanvas->cd();
    pad1->Draw();
    TPad *pad2 = new TPad("pad2","bottom pad",0,0,1,0.5);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.2);
    theCanvas->cd();
    pad2->Draw();

    pad1->cd();
    int MaxBin = m_hist_HM->GetMaximumBin();
    int MinBin = m_hist_HM->GetMinimumBin();
    float distance = m_hist_HM->GetBinContent(MaxBin) - m_hist_HM->GetBinContent(MinBin);
    distance /= 2.0;
    double Max = m_hist_HM->GetBinContent(MaxBin) + distance;
    double Min = m_hist_HM->GetBinContent(MinBin) - distance/3.;

    m_hist_HM->SetMarkerSize(0.8);
    m_hist_HM->SetYTitle("#it{Y}(#Delta#it{#phi})");
    m_hist_HM->GetYaxis()->SetRangeUser(Min, Max);
    m_hist_HM->GetListOfFunctions()->Add(f_HM);
    m_hist_HM->GetYaxis()->SetNdivisions(508,kTRUE);
    m_hist_HM->GetXaxis()->SetNdivisions(509,kTRUE);
    m_hist_HM->GetYaxis()->SetLabelSize(0.06);
    m_hist_HM->GetYaxis()->SetTitleSize(0.06);
    m_hist_HM->GetYaxis()->SetTitleOffset(1.30);
    m_hist_HM->GetXaxis()->SetTickLength(0.067);
    m_hist_HM->Draw("EX0SAME");
    f_show_flow2->SetLineColor(kBlue);
    f_show_flow2->SetLineStyle(2);
    f_show_flow2->Draw("same");
    if (!m_fixC3) {
        f_show_flow3->SetLineColor(kOrange+1);
        f_show_flow3->SetLineStyle(3);
        f_show_flow3->Draw("same");
    }
    f_show_periph->SetLineColor(kSpring+4);
    f_show_periph->SetLineStyle(2);
    f_show_periph->Draw("same");

    plotMarkerLineText(0.55, 0.85, 1.0,1, 20, 1,1,"HM Data", 0.05, true);
    plotMarkerLineText(0.55, 0.78, 0, 2, 1, 2, 1,"Fit", 0.05);
    plotMarkerLineText(0.55, 0.71, 0, kSpring+4, 0, kSpring+4, 2,"#it{G} + #it{F}#it{Y}^{LM}", 0.05);
    plotMarkerLineText(              0.74,0.85, 0, 4, 0, 4, 2,"#it{Y}_{2}^{ridge} + #it{F}#it{Y}^{LM}(0)",0.05);
    if (!m_fixC3) plotMarkerLineText(0.74,0.78, 0, kOrange+1, 0, kOrange+1, 3,"#it{Y}_{3}^{ridge} + #it{F}#it{Y}^{LM}(0)", 0.05);

    float _chi2 = 0;
    for (int i=1; i<h_pull->GetXaxis()->GetNbins()+1; i++){
        float _residual = m_hist_HM->GetBinContent(i) - f_HM->Eval(m_hist_HM->GetBinCenter(i));
        float _pull = _residual / m_hist_HM->GetBinError(i);
        h_pull->SetBinContent(i,_pull);
        h_pull->SetBinError(i,m_hist_HM->GetBinError(i)/m_hist_HM->GetBinError(i));
        _chi2 += pow(_residual/m_hist_HM->GetBinError(i),2);
    }
    _chi2 /= h_pull->GetXaxis()->GetNbins()-3;
    m_h_pull = (TH1F*) h_pull->Clone("m_h_Pull");

    pad2->cd();
    int MaxBin_sub = m_h_ridge->GetMaximumBin();
    int MinBin_sub = m_h_ridge->GetMinimumBin();
    double Max_sub = m_h_ridge->GetBinContent(MaxBin_sub);
    double Min_sub = m_h_ridge->GetBinContent(MinBin_sub);

    double sub_distance = Max_sub - Min_sub;

    m_h_ridge->SetMarkerSize(0.8);
    //m_h_ridge->SetLineWidth(1);
    m_h_ridge->GetYaxis()->SetRangeUser(Min_sub - sub_distance, Max_sub + 2.*sub_distance);

    m_h_ridge->SetXTitle("#Delta#it{#phi}");
    if (m_ridge_scaleFactor != 1.0) {
        m_h_ridge->SetYTitle(Form("(#it{Y}(#Delta#it{#phi}) - #it{F}#it{Y}^{LM}(#Delta#it{#phi}))#times%.0f", m_ridge_scaleFactor));
    } else {
        m_h_ridge->SetYTitle(Form("#it{Y}(#Delta#it{#phi}) - #it{F}#it{Y}^{LM}(#Delta#it{#phi})"));
    }
    m_h_ridge->GetYaxis()->SetNdivisions(406,kTRUE);
    m_h_ridge->GetYaxis()->SetLabelSize(0.06);
    m_h_ridge->GetYaxis()->SetTitleSize(0.06);
    m_h_ridge->GetYaxis()->SetTitleOffset(1.3);
    m_h_ridge->GetYaxis()->CenterTitle(kTRUE);

    m_h_ridge->GetXaxis()->SetTitleSize(0.07);
    m_h_ridge->GetXaxis()->SetTitleOffset(1.0);
    m_h_ridge->GetXaxis()->SetLabelSize(0.07);
    m_h_ridge->GetXaxis()->SetTickLength(0.08);
    m_h_ridge->Draw("EX0SAME");

    f_show_ridge->SetLineStyle(1);
    f_show_ridge->SetLineWidth(2);
    f_show_ridge->SetLineColor(2);
    f_show_ridge->Draw("same");

    f_show_ridge2->SetLineColor(kBlue);
    f_show_ridge2->SetLineStyle(2);
    f_show_ridge2->Draw("same");
    if (!m_fixC3) {
        f_show_ridge3->SetLineColor(kOrange+1);
        f_show_ridge3->SetLineStyle(3);
        f_show_ridge3->Draw("same");
    }
    plotMarkerLineText(0.60, 0.85, 1.0,1, 20, 1,1,"Data", 0.05, true);
    plotMarkerLineText(0.60, 0.78, 0, 2, 1, 2, 1,"#it{Y}^{ridge}", 0.05);
    plotMarkerLineText(              0.77,0.85, 0, 4, 0, 4, 2,"#it{Y}_{2}^{ridge}",0.05);
    if (!m_fixC3) plotMarkerLineText(0.77,0.78, 0, kOrange+1, 0, kOrange+1, 3,"#it{Y}_{3}^{ridge}", 0.05);

    pad1->cd();

    return true;
}

bool NonFlowSubtractor :: plotAtlasHM (TCanvas* theCanvas) {
    if (!theCanvas) return false;
    TH1F* h_pull = (TH1F*) m_hist_HM->Clone("h_Pull");

    theCanvas->cd();

    TPad *pad1 = new TPad("pad1","top pad",0,0.3,1,1);
    pad1->SetTopMargin(0.07);
    pad1->SetBottomMargin(0.);
    theCanvas->cd();
    pad1->Draw();
    TPad *pad2 = new TPad("pad2","bottom pad",0,0,1,0.3);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.3);
    theCanvas->cd();
    pad2->Draw();

    pad1->cd();
    int MaxBin = m_hist_HM->GetMaximumBin();
    int MinBin = m_hist_HM->GetMinimumBin();
    float distance = m_hist_HM->GetBinContent(MaxBin) - m_hist_HM->GetBinContent(MinBin);
    distance /= 2.0;
    double Max = m_hist_HM->GetBinContent(MaxBin) + distance;
    double Min = m_hist_HM->GetBinContent(MinBin) - distance/3.;
    //double Min = m_hist_HM->GetBinContent(MinBin) - distance;
    m_hist_HM->SetYTitle("#it{Y}(#Delta#it{#phi})");
    m_hist_HM->GetYaxis()->SetRangeUser(Min, Max);
    m_hist_HM->GetListOfFunctions()->Add(f_HM);
    m_hist_HM->GetYaxis()->SetNdivisions(508,kTRUE);
    m_hist_HM->GetXaxis()->SetNdivisions(509,kTRUE);
    m_hist_HM->GetYaxis()->SetLabelSize(0.06);
    m_hist_HM->GetYaxis()->SetTitleSize(0.06);
    m_hist_HM->GetYaxis()->SetTitleOffset(1.00);
    m_hist_HM->GetXaxis()->SetTickLength(0.067);
    m_hist_HM->Draw("EX0SAME");
    f_show_flow2->SetLineColor(kBlue);
    f_show_flow2->SetLineStyle(2);
    f_show_flow2->Draw("same");
    if (!m_fixC3) {
        f_show_flow3->SetLineColor(kOrange+1);
        f_show_flow3->SetLineStyle(3);
        f_show_flow3->Draw("same");
    }
    f_show_periph->SetLineColor(kSpring+4);
    f_show_periph->SetLineStyle(2);
    f_show_periph->Draw("same");

    plotMarkerLineText(0.55, 0.85, 1.2,1, 20, 1,1,"HM Data", 0.05, true);
    plotMarkerLineText(0.55, 0.78, 0, 2, 1, 2, 1,"Fit", 0.05);
    plotMarkerLineText(0.55, 0.71, 0, kSpring+4, 0, kSpring+4, 2,"#it{G} + #it{F}#it{Y}^{LM}", 0.05);
    plotMarkerLineText(              0.74,0.85, 0, 4, 0, 4, 2,"#it{Y}_{2}^{ridge} + #it{F}#it{Y}^{LM}(0)",0.05);
    if (!m_fixC3) plotMarkerLineText(0.74,0.78, 0, kOrange+1, 0, kOrange+1, 3,"#it{Y}_{3}^{ridge} + #it{F}#it{Y}^{LM}(0)", 0.05);

    float _chi2 = 0;
    for (int i=1; i<h_pull->GetXaxis()->GetNbins()+1; i++){
        float _residual = m_hist_HM->GetBinContent(i) - f_HM->Eval(m_hist_HM->GetBinCenter(i));
        float _pull = _residual / m_hist_HM->GetBinError(i);
        h_pull->SetBinContent(i,_pull);
        h_pull->SetBinError(i,m_hist_HM->GetBinError(i)/m_hist_HM->GetBinError(i));
        _chi2 += pow(_residual/m_hist_HM->GetBinError(i),2);
    }
    _chi2 /= h_pull->GetXaxis()->GetNbins()-3;
    m_h_pull = (TH1F*) h_pull->Clone("m_h_Pull");

    pad2->cd();
    h_pull->SetMaximum(4.9);
    h_pull->SetMinimum(-4.9);
    //h_pull->SetMaximum(20.9);
    //h_pull->SetMinimum(-20.9);
    h_pull->GetYaxis()->SetNdivisions(406,kTRUE);
    h_pull->GetYaxis()->SetLabelSize(0.14);
    h_pull->GetYaxis()->SetTitleSize(0.14);
    h_pull->GetYaxis()->SetTitleOffset(0.40);
    h_pull->GetYaxis()->CenterTitle(kTRUE);
    h_pull->GetXaxis()->SetTitleSize(0.13);
    h_pull->GetXaxis()->SetTitleOffset(1.0);
    h_pull->GetXaxis()->SetLabelSize(0.13);
    h_pull->GetXaxis()->SetTickLength(0.10);
    h_pull->SetXTitle("#Delta#it{#phi}");
    h_pull->SetYTitle("Pull");
    h_pull->Draw("EX0");
    TGraph* line0 = new TGraph(2);
    line0->SetPoint(0,-5,0);
    line0->SetPoint(1,15,0);
    line0->SetLineStyle(1);
    line0->SetLineColor(2);
    line0->SetLineWidth(2);
    line0->Draw("SAME");
    TGraph* line_p2 = new TGraph(2);
    line_p2->SetPoint(0,-5, 2);
    line_p2->SetPoint(1,15,2);
    line_p2->SetLineStyle(2);
    line_p2->SetLineColor(2);
    line_p2->SetLineWidth(2);
    TGraph* line_m2 = new TGraph(2);
    line_m2->SetPoint(0,-5, -2);
    line_m2->SetPoint(1,15,-2);
    line_m2->SetLineStyle(2);
    line_m2->SetLineColor(2);
    line_m2->SetLineWidth(2);
    line_p2->Draw("SAME");
    line_m2->Draw("SAME");

    plotText( 0.65, 0.44, 1, Form("#it{#chi}^{2}/ndof = %.2f",_chi2), 0.12);

    pad1->cd();
    return true;
}



bool NonFlowSubtractor :: plotAtlasLM (TPad* thePad) {
    if (!thePad) return false;
    thePad->cd();

    int MaxBin_ref = m_hist_LM->GetMaximumBin();
    int MinBin_ref = m_hist_LM->GetMinimumBin();
    double Max_ref = m_hist_LM->GetBinContent(MaxBin_ref)*1.8;
    double Min_ref = m_hist_LM->GetBinContent(MinBin_ref)*0.2;
    m_hist_LM->SetMarkerSize(0.8);
    m_hist_LM->SetLineWidth(1);
    m_hist_LM->SetXTitle("#Delta#it{#phi}");
    m_hist_LM->SetYTitle("#it{Y}(#Delta#it{#phi})");
    m_hist_LM->GetYaxis()->SetRangeUser(Min_ref, Max_ref);
    m_hist_LM->Draw("EX0SAME");
    m_hist_LM->SetMarkerStyle(24);
    f_LM->SetLineWidth(2);
    m_hist_LM->GetListOfFunctions()->Add(f_LM);
    return true;
}



bool NonFlowSubtractor :: plotAtlasLM (TCanvas* theCanvas) {
    if (!theCanvas) return false;

    m_hist_LM->SetMarkerStyle(24);
    TH1F* ref_pull = (TH1F*) m_hist_LM->Clone("ref_pull");

    theCanvas->cd();

    TPad *pad1 = new TPad("pad1","top pad",0,0.3,1,1);
    pad1->SetTopMargin(0.07);
    pad1->SetBottomMargin(0.);
    theCanvas->cd();
    pad1->Draw();
    TPad *pad2 = new TPad("pad2","bottom pad",0,0,1,0.3);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.3);
    theCanvas->cd();
    pad2->Draw();

    pad1->cd();
    int MaxBin_ref = m_hist_LM->GetMaximumBin();
    int MinBin_ref = m_hist_LM->GetMinimumBin();
    float distance = m_hist_LM->GetBinContent(MaxBin_ref) - m_hist_LM->GetBinContent(MinBin_ref);
    distance /= 4.;
    double Max_ref = m_hist_LM->GetBinContent(MaxBin_ref) + distance;
    double Min_ref = m_hist_LM->GetBinContent(MinBin_ref) - distance;
    m_hist_LM->SetXTitle("#Delta#it{#phi}");
    m_hist_LM->SetYTitle("#it{Y}(#Delta#it{#phi})");
    m_hist_LM->GetYaxis()->SetRangeUser(Min_ref, Max_ref);
    m_hist_LM->GetListOfFunctions()->Add(f_LM);
    m_hist_LM->GetYaxis()->SetNdivisions(508,kTRUE);
    m_hist_LM->GetXaxis()->SetNdivisions(509,kTRUE);
    m_hist_LM->GetYaxis()->SetLabelSize(0.06);
    m_hist_LM->GetYaxis()->SetTitleSize(0.06);
    m_hist_LM->GetYaxis()->SetTitleOffset(1.00);
    m_hist_LM->GetXaxis()->SetTickLength(0.067);
    m_hist_LM->Draw("EX0SAME");
    plotMarkerLineText(0.25,0.42, 1.2,1, 24, 1,1,"LM Data", 0.05, true);
    plotMarkerLineText(0.25,0.36, 0, 2, 1, kSpring-6, 3,"LM Fourier Fit", 0.05);

    float _chi2 = 0;
    for (int i=1; i<ref_pull->GetXaxis()->GetNbins()+1; i++){
        float _residual = m_hist_LM->GetBinContent(i) - f_LM->Eval(m_hist_LM->GetBinCenter(i));
        float _pull = _residual / m_hist_LM->GetBinError(i);
        ref_pull->SetBinContent(i,_pull);
        ref_pull->SetBinError(i,m_hist_LM->GetBinError(i)/m_hist_LM->GetBinError(i));
        _chi2 += pow(_residual/m_hist_LM->GetBinError(i),2);
    }
    _chi2 /= ref_pull->GetXaxis()->GetNbins()-4;

    pad2->cd();
    ref_pull->SetMaximum(4.9);
    ref_pull->SetMinimum(-4.9);
    ref_pull->GetYaxis()->SetNdivisions(406,kTRUE);
    ref_pull->GetYaxis()->SetLabelSize(0.14);
    ref_pull->GetYaxis()->SetTitleSize(0.14);
    ref_pull->GetYaxis()->SetTitleOffset(0.40);
    ref_pull->GetYaxis()->CenterTitle(kTRUE);
    ref_pull->GetXaxis()->SetTitleSize(0.13);
    ref_pull->GetXaxis()->SetTitleOffset(1.0);
    ref_pull->GetXaxis()->SetLabelSize(0.13);
    ref_pull->GetXaxis()->SetTickLength(0.10);
    ref_pull->SetXTitle("#Delta#it{#phi}");
    ref_pull->SetYTitle("Pull");
    ref_pull->Draw("EX0");
    TGraph* line0 = new TGraph(2);
    line0->SetPoint(0,-5,0);
    line0->SetPoint(1,15,0);
    line0->SetLineStyle(1);
    line0->SetLineColor(2);
    line0->SetLineWidth(2);
    line0->Draw("SAME");
    TGraph* line_p2 = new TGraph(2);
    line_p2->SetPoint(0,-5, 2);
    line_p2->SetPoint(1,15,2);
    line_p2->SetLineStyle(2);
    line_p2->SetLineColor(2);
    line_p2->SetLineWidth(2);
    line_p2->Draw("SAME");
    TGraph* line_m2 = new TGraph(2);
    line_m2->SetPoint(0,-5, -2);
    line_m2->SetPoint(1,15,-2);
    line_m2->SetLineStyle(2);
    line_m2->SetLineColor(2);
    line_m2->SetLineWidth(2);
    line_m2->Draw("SAME");

    plotText( 0.65, 0.44, 1, Form("#it{#chi}^{2}/ndof = %.2f",_chi2), 0.12);
    return true;
}
