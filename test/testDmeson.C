#include "../../CorrelationButton/CorrelationMaker.h"
#include "../../CorrelationButton/NonFlowSubtractor.C"

int testDmeson() {

    // example to show how to run the tool with the additional dimensions
    // copied from the D-h correlation analysis in p+Pb
    FlowAnaConfig* theConfig = new FlowAnaConfig();
    
    // configurate the analyss
    // input/output
    theConfig->setInputPath("../Selection/");
    theConfig->setInputFileName("correlation_oct");
    theConfig->setOutputPath("./RootFiles");
    theConfig->setOutputFileName("_result_nominal");
    theConfig->setOutputFigPath("./Plots");

    theConfig->setCorrHistPathSame("Correlation_2D_wgt/h2pc_sig_weighted_");
    theConfig->setCorrHistPathMix ("Correlation_2D_wgt/h2pc_mix_weighted_");
    theConfig->setTrigYieldHistName("h2_nsig_weighted");
    theConfig->setMixTrigYieldHistName("h2_nsig_weighted");
    
    // gap/referece/cut
    theConfig->setReference (0,  40); // LM reference Nch selection
    theConfig->setReference2(40, 80); // LM2 selection
    theConfig->setEtaRange (1.5, 4.5); // Gap

    theConfig->setCorrEtaBoundary(4.5);
    theConfig->setCorrEtaInterval(0.1);

    // addtional interface for configuring the whole analysis 
    // pt/Nch/... binning of outputs
    // inconsistency between input and ouput binning will be checked since input bins can only be merged to produce output
    const int Nbins_out_Nch = 3;
    float dNch_out_range[Nbins_out_Nch+1] = {40,180,240,300};
    theConfig->setOutputMultiBinning(Nbins_out_Nch, dNch_out_range);

    const int Nbins_out_pt = 3;
    float dpt_out_range[Nbins_out_pt+1] = {3, 4, 6, 10};
    theConfig->setOutputPtBinning(Nbins_out_pt, dpt_out_range);

    theConfig->setInputThirdBinningName("hmass_binning");
    theConfig->setThirdDimensionName("mass");
    const int Nbins_out_mass = 13;
    float dmass_out_range[Nbins_out_mass+1] = {1.75, 1.77, 1.799, 1.821, 1.843, 1.858, 1.872, 1.887, 1.909, 1.931, 1.96, 2.00, 2.05, 2.1};
    theConfig->setOutputThirdBinning(Nbins_out_mass, dmass_out_range);

    if (!theConfig->init()) {
        return 1;
    }

    CorrelationMaker* theMaker = new CorrelationMaker();
    theMaker->setConfig(theConfig);
    theMaker->setRebin(3); 
    theMaker->setDebug();
    theMaker->init();

    TCanvas* c0;

    // prepare histogram here
    // if the output binning is specified, loop over the intervals
    for (int index_Nch = 1; index_Nch < theConfig->getOutputMultiBinning()->GetNbinsX()+1; index_Nch++) {
        float _Nch_low  = theConfig->getOutputMultiBinning()->GetXaxis()->GetBinLowEdge(index_Nch);
        float _Nch_high = theConfig->getOutputMultiBinning()->GetXaxis()->GetBinUpEdge (index_Nch);

        for (int index_pt = 1; index_pt < theConfig->getOutputPtBinning()->GetNbinsX()+1; index_pt++) {
            float _pt_trig_low  = theConfig->getOutputPtBinning()->GetXaxis()->GetBinLowEdge(index_pt);
            float _pt_trig_high = theConfig->getOutputPtBinning()->GetXaxis()->GetBinUpEdge (index_pt);
            c0 = new TCanvas("c0","",1400,800);
            c0->Divide(5,3);
            int ipad = 1;

            NonFlowSubtractor subTool;
            subTool.init();
            if (theConfig->getOutputThirdBinning()) {
                TH1F* h_v22_raw = (TH1F*) theConfig->getOutputThirdBinning()->Clone(Form("h_v22_Nch%d_pt%d_raw_d%s", index_Nch-1, index_pt-1, theConfig->getThirdDimensionName().c_str()) );
                TH1F* h_v22     = (TH1F*) theConfig->getOutputThirdBinning()->Clone(Form("h_v22_Nch%d_pt%d_sub_d%s", index_Nch-1, index_pt-1, theConfig->getThirdDimensionName().c_str()) );
                TH1F* h_v22_cms = (TH1F*) theConfig->getOutputThirdBinning()->Clone(Form("h_v22_Nch%d_pt%d_cms_d%s", index_Nch-1, index_pt-1, theConfig->getThirdDimensionName().c_str()) );
                h_v22->Reset();
                h_v22_raw->Reset();
                h_v22_cms->Reset();

                for (int index = 1; index < theConfig->getOutputThirdBinning()->GetNbinsX()+1; index++) {
                    float _mass_low  = theConfig->getOutputThirdBinning()->GetXaxis()->GetBinLowEdge(index);
                    float _mass_high = theConfig->getOutputThirdBinning()->GetXaxis()->GetBinUpEdge (index);

                    TH1* hist_LM = theMaker->MakeCorr(_pt_trig_low, _pt_trig_high, theConfig->getReferenceLow(),theConfig->getReferenceHigh(), _mass_low, _mass_high);
                    TH1* hist_LM2= theMaker->MakeCorr(_pt_trig_low, _pt_trig_high, theConfig->getReference2Low(),theConfig->getReference2High(), _mass_low, _mass_high);

                    TH1* hist_HM  = theMaker->MakeCorr(_pt_trig_low, _pt_trig_high, _Nch_low, _Nch_high, _mass_low, _mass_high);

                    // apply subtraction here
                    subResult theResult = subTool.templateFit(hist_LM, hist_HM);
                    cout << "v22 = " << theResult.getV22SubValue() << " +/- " << theResult.getV22SubError() << endl;
                    cout << "Improved v22 = " << theResult.getV22SubImpValue() << " +/- " << theResult.getV22SubImpError() << endl;

                    h_v22_raw->SetBinContent(index, theResult.getV22RawValue());
                    h_v22_raw->SetBinError  (index, theResult.getV22RawError());
                    h_v22->SetBinContent(index, theResult.getV22SubValue());
                    h_v22->SetBinError  (index, theResult.getV22SubError());

                    TPad* thePad = (TPad*)c0->cd(ipad);
                    subTool.plotAtlasHM(thePad);
                    plotText(    0.2,0.88,1,Form("%.2f < #it{m} < %.2f GeV", _mass_low, _mass_high),0.06);
                    ipad++;
                    
                    // add cms-like results for comparison
                    NonFlowSubtractor subTool_cms;
                    subTool_cms.init();
                    subResult theResult_cms = subTool_cms.templateFit(hist_LM, hist_HM, 0., 0.);
                    h_v22_cms->SetBinContent(index, theResult_cms.getV22SubImpValue());
                    h_v22_cms->SetBinError  (index, theResult_cms.getV22SubImpError());
                }
                theConfig->outputFile()->cd();
                h_v22_raw->Write();
                h_v22->Write();
                h_v22_cms->Write();

            } // mass slices

            TPad* thePad = (TPad*)c0->cd(ipad);
            subTool.plotAtlasHMLabel(thePad);
            plotText( 0.2,0.38,1,"#it{D-h} Correlation", 0.08);
            ipad++;

            thePad = (TPad*)c0->cd(ipad);
            plotText( 0.2,0.88,1,"#font[72]{ATLAS} Internal",0.08);
            plotText( 0.2,0.78,1,Form("%.0f < #||{#Delta#it{#eta}} < %.0f", theConfig->getEtaRangeLow(), theConfig->getEtaRangeHigh()), 0.08);
            plotText( 0.2,0.68,1,Form("LM %.0f #leq #it{N}_{ch} < %.0f", 0., 40.), 0.08);
            plotText( 0.2,0.58,1,Form("HM %.0f #leq #it{N}_{ch} < %.0f", _Nch_low, _Nch_high), 0.08);
            plotText( 0.2,0.48,1,Form("%.1f #leq #it{p}_{T}^{trig} < %.1f GeV", _pt_trig_low, _pt_trig_high), 0.08);
            plotText( 0.2,0.36,1,Form("%.1f #leq #it{p}_{T}^{ass0} < %.1f GeV", 0.5, 5.0), 0.08);

            c0->SaveAs(Form("%s/Dh_nominalTempFit_Nch%d_pt%d.pdf",theConfig->getOutputFigPath().c_str(), index_Nch, index_pt));
            delete c0; c0=0;

        }
    }

    theConfig->outputFile()->cd();
    theConfig->getOutputMultiBinning()->Write();
    theConfig->getOutputPtBinning()->Write();
    theConfig->getOutputThirdBinning()->Write();

    return 0;
}
