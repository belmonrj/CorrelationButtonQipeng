#include "../NonFlowSubtractor.C"
void testJetSample() {

    //TFile* fin = new TFile(Form("./dphiProjections_t14_centrality.root"),"READ");
    TFile* fin = new TFile(Form("./histosForQipeng.root"),"READ");

    TH1F* h_correlation_LM  = (TH1F*)fin->Get( Form("h_mbPythia_peripheral") );
    TH1F* h_correlation_HM  = (TH1F*)fin->Get( Form("h_mbPythia_central") );
    //TH1F* h_correlation_LM  = (TH1F*)fin->Get( Form("h_mb_peripheral") );
    //TH1F* h_correlation_HM  = (TH1F*)fin->Get( Form("h_mb_central") );
    TH1F* h_correlation_LM_bulk  = (TH1F*)fin->Get( Form("h_mbPythia_rejectJet_peripheral") );

    TH1F* h_correlation_HM_bulk  = (TH1F*)fin->Get( Form("h_mbPythia_rejectJet_central") );

    NonFlowSubtractor subTool;
    subTool.setNHar(4);// fit with v2, v3 and v4;
    subTool.init();

    NonFlowSubtractor subTool2;
    subTool2.setNHar(4);// fit with v2, v3 and v4;
    subTool2.init();

    //--------------------------------------------------
    // ATLAS template fit
    //--------------------------------------------------
    subResult theResult = subTool.templateHistFit(h_correlation_LM, h_correlation_LM_bulk, h_correlation_HM);

    subResult theResult2 = subTool2.templateHistFit(h_correlation_LM, h_correlation_HM);
    theResult2 = subTool2.templateHistFit(h_correlation_LM, h_correlation_HM);
    cout << theResult2.getChi2() << endl;

    //--------------------------------------------------
    // Access the fitted results and plots
    //--------------------------------------------------
    //cout << "v22 = " << theResult.getV22SubValue() << " +/- " << theResult.getV22SubError() << endl;
    //cout << "Improved v22 = " << theResult.getV22SubImpValue() << " +/- " << theResult.getV22SubImpError() << endl;

    TCanvas* c1 = new TCanvas("c1","scaling",50,50, 600,700);
    subTool.plotAtlasHistSubHM(c1);
    //c1->cd();
    // add addtional lengends
    plotText( 0.2,0.89,1,"#font[72]{ATLAS} Internal");
    //plotText( 0.2,0.84,1,Form("%.1f < #||{#Delta#it{#eta}} < %.1f", theConfig->getEtaRangeLow(), theConfig->getEtaRangeHigh()));
    //plotText( 0.2,0.79,1,Form("HM %.0f #leq #it{N}_{ch} < %.0f", _Nch_low, _Nch_high));
    //plotText( 0.2,0.74,1,Form("LM %.0f #leq #it{N}_{ch} < %.0f", theConfig->getReferenceLow(), theConfig->getReferenceHigh()));
    plotText( 0.2,0.45,1,Form("#it{v}_{2,2}#times10^{3}= %.4f #pm %.4f", theResult.getV22SubValue()*1e3, theResult.getV22SubError()*1e3) );
    plotText( 0.2,0.40,1,Form("#it{v}_{3,3}#times10^{3}= %.4f #pm %.4f", theResult.getCoeffSubValue(3)*1e3, theResult.getCoeffSubError(3)*1e3) );
    //plotText( 0.2,0.15,1,Form("%.1f < #it{p}_{T}^{trig} < %.1f GeV", _pt_trig_low, _pt_trig_high));
    //c1->SaveAs("old.pdf");

    TCanvas* c2 = new TCanvas("c2","scaling",50,50, 600,700);
    subTool2.plotAtlasHistSubHM(c2);
    plotText( 0.2,0.89,1,"#font[72]{ATLAS} Internal");
    plotText( 0.2,0.45,1,Form("#it{v}_{2,2}#times10^{3}= %.4f #pm %.4f", theResult2.getV22SubValue()*1e3,    theResult2.getV22SubError()*1e3) );
    plotText( 0.2,0.40,1,Form("#it{v}_{3,3}#times10^{3}= %.4f #pm %.4f", theResult2.getCoeffSubValue(3)*1e3, theResult2.getCoeffSubError(3)*1e3) );
}
