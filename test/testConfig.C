#include "../CorrelationMaker.h"
#include "../NonFlowSubtractor.C"

int testConfig() {
    FlowAnaConfig* theConfig = new FlowAnaConfig();
    
    // configurate the analyss
    // input/output
    theConfig->setInputPath("./");
    theConfig->setInputFileName("condor_output");
    theConfig->setOutputPath("./");
    theConfig->setOutputFileName("_result_testConfig");
    theConfig->setOutputFigPath("./figures");
    
    // gap/referece/cut
    theConfig->setReference (10, 20); // LM reference Nch selection
    theConfig->setReference2(20, 30); // LM2 selection
    theConfig->setEtaRange (2.5, 5.0); // Gap

    // addtional interface for configuring the whole analysis 
    // pt/Nch/... binning of outputs
    // inconsistency between input and ouput binning will be checked since input bins can only be merged to produce output
    const int Nbins_out_Nch = 13;
    float dNch_out_range[Nbins_out_Nch+1] = {20,30,40,50,60,70,80,90,100,110,120,130,140,150};
    theConfig->setOutputMultiBinning(Nbins_out_Nch, dNch_out_range);

    const int Nbins_out_pt = 9;
    float dpt_out_range[Nbins_out_pt+1] = {0.5,1,2,3,4,5,6,7,8,10};
    theConfig->setOutputPtBinning(Nbins_out_pt, dpt_out_range);

    if (!theConfig->init()) {
        return 1;
    }

    CorrelationMaker* theMaker = new CorrelationMaker();
    theMaker->setConfig(theConfig);
    theMaker->init();
    //theMaker->setDebug(true);

    // prepare histogram here
    // if the output binning is specified, loop over the intervals
    if (theConfig->getOutputMultiBinning()) {
        float _pt_trig_low = 0.5;
        float _pt_trig_high = 5.;
        TH1F* h_v22_dNch = (TH1F*) theConfig->getOutputMultiBinning()->Clone("h_v22_dNch");
        h_v22_dNch->Reset();

        TH1* hist_LM  = theMaker->MakeCorr(_pt_trig_low, _pt_trig_high, theConfig->getReferenceLow(), theConfig->getReferenceHigh());
        TH1* hist_LM2 = theMaker->MakeCorr(_pt_trig_low, _pt_trig_high, theConfig->getReference2Low(),theConfig->getReference2High());

        for (int index = 1; index < theConfig->getOutputMultiBinning()->GetNbinsX()+1; index++) {
            float _Nch_low  = theConfig->getOutputMultiBinning()->GetXaxis()->GetBinLowEdge(index);
            float _Nch_high = theConfig->getOutputMultiBinning()->GetXaxis()->GetBinUpEdge (index);
            TH1* hist_HM  = theMaker->MakeCorr(_pt_trig_low, _pt_trig_high, _Nch_low, _Nch_high);

            // apply subtraction here
            NonFlowSubtractor subTool;
            subTool.init();
            subResult theResult = subTool.templateFit(hist_LM, hist_HM, hist_LM2);
            //cout << "v22 = " << theResult.getV22SubValue() << " +/- " << theResult.getV22SubError() << endl;
            //cout << "Improved v22 = " << theResult.getV22SubImpValue() << " +/- " << theResult.getV22SubImpError() << endl;
            h_v22_dNch->SetBinContent(index, theResult.getV22SubImpValue());
            h_v22_dNch->SetBinError  (index, theResult.getV22SubImpError());

            // make some plots here
            TCanvas* c1 = new TCanvas("c1","scaling",50,50, 600,600);
            subTool.plotAtlasHM(c1);
            c1->cd();
            // same as myText() in ATLAS plotstyle
            plotText( 0.5,0.88,1,"#font[72]{ATLAS} Internal");
            plotText( 0.5,0.82,1,Form("%.0f < #||{#Delta#it{#eta}} < %.0f", theConfig->getEtaRangeLow(), theConfig->getEtaRangeHigh()));
            plotText( 0.5,0.76,1,Form("%.0f #leq #it{N}_{ch} < %.0f", _Nch_low, _Nch_high));
            plotText( 0.2,0.20,1,"#it{h-h} Correlation");
            c1->SaveAs(Form("%s/testfig_dNch%d.pdf",theConfig->getOutputFigPath().c_str(), index));
            delete c1; c1 = 0;
            
            //TCanvas* c2 = new TCanvas("c2","scaling",50,50, 600,600);
            //subTool.plotAtlasLM(c2);
            //c2->cd();
            // add addtional lengends
        }
        theConfig->outputFile()->cd();
        h_v22_dNch->Write();

    }

    // test 2D correlation maker (pt_low, pt_high, Nch_low, Nch_high)
    

    return 0;
}
