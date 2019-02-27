#include "NonFlowSubtractor.h"


subResult NonFlowSubtractor :: templateFit(TH1* hist_LM, TH1* hist_HM) {

    const int Npar(2*getNHar()+3);
    const int Npar_LM(getNHar()+1);

    subResult theResult;

    f_LM = new TF1("f_LM", f_periph, m_dphiRangeLow, m_dphiRangeHigh, Npar_LM);
    f_LM->SetLineColor(kSpring-6);
    f_LM->SetLineStyle(3);

    f_HM = new TF1("f_HM", templ, m_dphiRangeLow, m_dphiRangeHigh, Npar);
    f_HM->SetLineColor(2);


    vector <double> parInitValues;
    parInitValues.clear();
    for (int i=0; i<Npar; i++) { 
        parInitValues.push_back(0);
    }

    //perform fourier fit first to get initial values
    subResult _init_LM = fourierFit(hist_LM);
    double _value_G_LM = _init_LM.getPedstalValue();
    double _error_G_LM = _init_LM.getPedstalError();
    parInitValues.at(0) = _init_LM.getPedstalValue();
    for (int i=1; i<m_nparLmAtlas; i++) {
        parInitValues.at(i) = _init_LM.getCoeffRawValue(i);
    }

    subResult _init_HM = fourierFit(hist_HM);
    vector<float> vec_value_raw;
    vector<float> vec_error_raw;
    for (int i=Npar_LM; i<Npar_LM+getNHar(); i++) { 
        vec_value_raw.push_back(_init_HM.getCoeffRawValue(i-m_nparLmAtlas+1));
        vec_error_raw.push_back(_init_HM.getCoeffRawError(i-m_nparLmAtlas+1));
        parInitValues.at(i) = _init_HM.getCoeffRawValue(i-m_nparLmAtlas+1);
    }
    theResult.setCoeffRaw(vec_value_raw, vec_error_raw);
    double _value_G_HM = _init_HM.getPedstalValue();
    double _error_G_HM = _init_HM.getPedstalError();

    parInitValues.at(Npar-2) = 5.;
    parInitValues.at(Npar-1) = _init_HM.getPedstalValue();

    if (m_debug) {
        for (int i=0; i<Npar; i++) { 
            cout << i << ": " << parInitValues.at(i) << endl;
        }
    }

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
    GlobalChi2 globalChi2(chi2_LM, chi2_HM, getNHar());

    ROOT::Fit::Fitter fitter;
    fitter.Config().SetParamsSettings(Npar, &parInitValues.at(0));
    fitter.Config().MinimizerOptions().SetPrintLevel(0);
    fitter.Config().SetMinimizer("Minuit2","Migrad");

    if (m_fixC1) {
        fitter.Config().ParSettings(Npar_LM).SetValue(0);
        fitter.Config().ParSettings(Npar_LM).Fix();
    }

    for (int ipar = 0; ipar < Npar; ipar++) {
        fitter.Config().ParSettings(ipar).SetName(m_parName_HM.at(ipar).c_str());
    }
    
    if (m_debug) {
        cout << "Fit parameter configurations: " << endl;
        for (int ipar = 0; ipar < fitter.Config().NPar(); ipar++) {
            cout << "parameter " << ipar << ":\tname = " << fitter.Config().ParSettings(ipar).Name().c_str() << ",\tvalue = " << fitter.Config().ParSettings(ipar).Value() << endl;
        }
    }

    fitter.FitFCN(Npar, globalChi2, 0, data_LM.Size()+data_HM.Size(), true);

    ROOT::Fit::FitResult result = fitter.Result();
    result.Print(std::cout);
    //result.PrintCovMatrix(std::cout);

    vector<float> vec_value_sub;
    vector<float> vec_error_sub;
    for (int ihar=0; ihar<getNHar(); ihar++){
        vec_value_sub.push_back(result.Value   (Npar_LM+ihar) );
        vec_error_sub.push_back(result.ParError(Npar_LM+ihar) );
    }
    theResult.setCoeffSub(vec_value_sub, vec_error_sub);
    double _value_F_temp = result.Value   (Npar-2);
    double _error_F_temp = result.ParError(Npar-2);

    double _value_rho_atlas = (_value_G_LM*_value_F_temp)/_value_G_HM;
    // naive error estimation
    double _error_rho_atlas = _value_rho_atlas*sqrt(pow(_error_G_LM/_value_G_LM, 2) + pow(_error_G_HM/_value_G_HM, 2));

    theResult.setNHar(getNHar());
    theResult.setChi2(result.Chi2() / result.Ndf());
    theResult.setPedstalValue (_value_G_HM);
    theResult.setPedstalError (_error_G_HM);
    theResult.setRhoValue (_value_rho_atlas);
    theResult.setRhoError (_error_rho_atlas);

    if (m_debug) {
        cout << "_value_G_LM = " << _value_G_LM << endl;
        cout << "_value_G_HM = " << _value_G_HM << endl;
        cout << "_value_F_temp = " << _value_F_temp << endl;
        cout << "rho = " << _value_rho_atlas << endl;
    }

    // construct stuffs for plotting
    // -------------------------------------------------------------
    // Prepare histograms and function for plotting 
    // uncorrected cn parameters are used for visulization
    //if (m_applyZYAM) {
    //    string _f_periph_forPlot = "[0] * ( 1 + 2*(";
    //    for (int ihar=0; ihar<getNHar(); ihar++){
    //        _f_periph_forPlot += Form("[%d]*(cos(%d*x)-1)", ihar+1, ihar+1);
    //        if (ihar != getNHar()-1) _formula += (" + ");
    //    }
    //    _f_periph_forPlot += (") )*[%d] + [%d]", getNHar()+1, getNHar()+2);
    //    f_show_periph = new TF1("f_show_periph", _f_periph_forPlot.c_str(), m_dphiRangeLow, m_dphiRangeHigh);
    //} else {
    if (1) {
        string _f_periph_forPlot = "[0] * ( 1 + 2*(";
        for (int ihar=0; ihar<getNHar(); ihar++){
            _f_periph_forPlot += Form("[%d]*cos(%d*x)", ihar+1, ihar+1);
            if (ihar != getNHar()-1) _f_periph_forPlot += (" + ");
        }
        _f_periph_forPlot += Form(") )*[%d] + [%d]", getNHar()+1, getNHar()+2);
        if (m_debug) cout << "Scaled peripheral function: " << _f_periph_forPlot.c_str() << endl;
        f_show_periph = new TF1("f_show_periph", _f_periph_forPlot.c_str(), m_dphiRangeLow, m_dphiRangeHigh);
    }

    for (int ipar = 0; ipar < getNHar()+1; ipar++) {
        f_show_periph->SetParameter(ipar, result.Value(ipar));
    }

    f_show_periph->SetParameter(getNHar()+1, result.Value(Npar - 2));
    f_show_periph->SetParameter(getNHar()+2, result.Value(Npar - 1));

    f_show_flow2 = new TF1("f_show_flow2","[0] + [1]*(1+2*[2]*cos(2*x))", m_dphiRangeLow, m_dphiRangeHigh);
    f_show_flow2->SetParameter(0, f_LM->Eval(0)*result.Value(Npar - 2));
    f_show_flow2->SetParameter(1, result.Value(Npar - 1));
    f_show_flow2->SetParameter(2, result.Value(getNHar()+2)); // c2


    f_show_flow3 = new TF1("f_show_flow3","[0] + [1]*(1+2*[2]*cos(3*x) )", m_dphiRangeLow, m_dphiRangeHigh);
    f_show_flow3->SetParameter(0, f_LM->Eval(0)*result.Value(Npar - 2));
    f_show_flow3->SetParameter(1, result.Value(Npar - 1));
    f_show_flow3->SetParameter(2, result.Value(getNHar()+3)); // c3

    if (m_debug) {
        // double check the parameter assignment is corrected
        cout << "c2 = " << result.Value(getNHar()+2) << endl;
        cout << "c3 = " << result.Value(getNHar()+3) << endl;
        cout << "G_temp = " << result.Value(Npar - 1) << endl;
        cout << "F_temp = " << result.Value(Npar - 2) << endl;
    }

    m_hist_HM = (TH1F*)hist_HM->Clone("__hist_HM");
    m_hist_LM = (TH1F*)hist_LM->Clone("__hist_LM");

    //delete f_LM; f_LM = 0;
    //delete f_HM; f_HM = 0;
    return theResult;
}



subResult NonFlowSubtractor::templateFit (TH1* hist_LM,  TH1* hist_HM, TH1* hist_LM2) {

    subResult result_hm  = templateFit(hist_LM, hist_HM);

    // for keeping the plot stuff to be the same as running unimproved fit
    // correction fit is done in different object 
    NonFlowSubtractor subLM2;
    subLM2.setNHar(NonFlowSubtractor::m_nhar);
    subLM2.init();
    subResult result_lm2 = subLM2.templateFit(hist_LM, hist_LM2);

    vector<float> vec_subImp_value;
    vector<float> vec_subImp_error;
    for (int ihar=0; ihar<getNHar(); ihar++){
        float _rho      = result_hm .getRhoValue();
        float _coeff_hm = result_hm .getCoeffSubValue(ihar+1);
        float _error_hm = result_hm .getCoeffSubError(ihar+1);
        float _coeff_lm = result_lm2.getCoeffSubValue(ihar+1);
        float _error_lm = result_lm2.getCoeffSubError(ihar+1);

        //ignore error on rho
        float _coeff_corrected = _coeff_hm - _rho*(_coeff_hm - _coeff_lm);
        float _error = sqrt( pow((1-_rho)*_error_hm,2) + pow( _rho*_error_lm,2)); // to be improved
        if (_coeff_hm == 0) _error = 0; // protection for c1 = 0

        if (m_debug) cout << ihar << ",\t" << _coeff_hm << ",\t" << _coeff_lm << ",\t" << _coeff_lm << endl;

        vec_subImp_value.push_back(_coeff_corrected);
        vec_subImp_error.push_back(_error);
    }
    result_hm.setCoeffSubImp(vec_subImp_value, vec_subImp_error);

    return result_hm;
}



subResult NonFlowSubtractor::fourierFit (TH1* hist, TH1* hist_LM) {

    TF1* f_fourier = new TF1("f_fourier", _formula.c_str(), m_dphiRangeLow, m_dphiRangeHigh);
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



subResult NonFlowSubtractor::periphSub  ( TH1* h_sr_lm, TH1* h_lr_lm, 
                                          TH1* h_sr_hm, TH1* h_lr_hm) {
    // SR - LR for low multiplicity
    TH1F* h_sub_lm = (TH1F*)h_sr_lm->Clone("h_sub_lm");
    h_sub_lm->Add(h_lr_lm, -1);
    
    // SR - LR for high multiplicity
    TH1F* h_sub_hm = (TH1F*)h_sr_hm->Clone("h_sub_hm");
    h_sub_hm->Add(h_lr_hm, -1);

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
    subLM2.setNHar(NonFlowSubtractor::m_nhar);
    subLM2.init();
    subResult result_lm2 = subLM2.periphSub(h_sr_lm, h_lr_lm, h_sr_lm2, h_lr_lm2);

    vector<float> vec_subImp_value;
    vector<float> vec_subImp_error;
    for (int ihar=0; ihar<getNHar(); ihar++){
        float _rho      = result_hm .getRhoValue();
        float _coeff_hm = result_hm .getCoeffSubValue(ihar);
        float _error_hm = result_hm .getCoeffSubError(ihar);
        float _coeff_lm = result_lm2.getCoeffSubValue(ihar);
        float _error_lm = result_lm2.getCoeffSubError(ihar);

        float _coeff_corrected = _coeff_hm + _rho*(_coeff_lm);
        float _error = sqrt(pow(_error_hm,2) + pow(_rho*_error_lm,2)); // ignore error on rho, to be improved
        vec_subImp_value.push_back(_coeff_corrected);
        vec_subImp_error.push_back(_error);
    }

    //update result_hm with imporoved fit results
    result_hm.setCoeffSubImp(vec_subImp_value, vec_subImp_error);

    return result_hm;
}





bool NonFlowSubtractor :: plotAtlasHM (TPad* thePad) {
    if (!thePad) return false;
    thePad->cd();

    int MaxBin = m_hist_HM->GetMaximumBin();
    int MinBin = m_hist_HM->GetMinimumBin();
    double Max = 10;
    double Min = 0.;
    //m_hist_HM->GetYaxis()->SetRangeUser(Min, Max);
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
    if (!m_fixC3) {
        f_show_flow3->SetLineWidth(2);
        f_show_flow3->SetLineColor(kOrange+1);
        f_show_flow3->SetLineStyle(3);
        f_show_flow3->Draw("same");
    }
    f_show_periph->SetLineWidth(2);
    f_show_periph->SetLineColor(kSpring+4);
    f_show_periph->SetLineStyle(2);
    f_show_periph->Draw("same");
    plotMarkerText(    0.60,0.88, 1, 20, "Data", 0.7);
    plotMarkerLineText(0.60,0.83, 0, 2, 1, 2, 1,"Fit");
    plotMarkerLineText(0.60,0.78, 0, kSpring+4, 0, kSpring+4, 2,"#it{g} + #it{f}#it{Y}^{ periph}");
    plotMarkerLineText(              0.78,0.83, 0, 4, 0, 4, 2,"#it{c}_{2}");
    if (!m_fixC3) plotMarkerLineText(0.78,0.78, 0, kOrange+1, 0, kOrange+1, 3,"#it{c}_{3}");

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
    distance /= 4.;
    double Max = m_hist_HM->GetBinContent(MaxBin) + distance;
    double Min = m_hist_HM->GetBinContent(MinBin) - distance;
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

    plotMarkerLineText(0.25,0.84, 1.2,1, 20, 1,1,"HM Data", 0.05, true);
    plotMarkerLineText(0.25,0.77, 0, 2, 1, 2, 1,"Fit", 0.05);
    plotMarkerLineText(0.25,0.70, 0, kSpring+4, 0, kSpring+4, 2,"#it{G} + #it{F}#it{Y}^{LM}", 0.05);
    plotMarkerLineText(              0.25,0.63, 0, 4, 0, 4, 2,"#it{c}_{2}",0.05);
    if (!m_fixC3) plotMarkerLineText(0.25,0.56, 0, kOrange+1, 0, kOrange+1, 3,"#it{c}_{3}", 0.05);

    float _chi2 = 0;
    for (int i=1; i<h_pull->GetXaxis()->GetNbins()+1; i++){
        float _residual = m_hist_HM->GetBinContent(i) - f_HM->Eval(m_hist_HM->GetBinCenter(i));
        float _pull = _residual / m_hist_HM->GetBinError(i);
        h_pull->SetBinContent(i,_pull);
        h_pull->SetBinError(i,m_hist_HM->GetBinError(i)/m_hist_HM->GetBinError(i));
        _chi2 += pow(_residual/m_hist_HM->GetBinError(i),2);
    }
    _chi2 /= h_pull->GetXaxis()->GetNbins()-3;

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

    myText( 0.65, 0.44, 1, Form("#it{#chi}^{2}/ndof = %.2f",_chi2), 0.12);

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

    myText( 0.65, 0.44, 1, Form("#it{#chi}^{2}/ndof = %.2f",_chi2), 0.12);
    return true;
}
