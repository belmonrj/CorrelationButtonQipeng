#include "../NonFlowSubtractor.C"
void testErrors() {

    TFile* fin = new TFile(Form("./correlationHistTest.root"),"READ");

    TH1F* h_correlation_LM  = (TH1F*)fin->Get( Form("PTYRaw/dphi_Nch_ref") );
    TH1F* h_correlation_HM  = (TH1F*)fin->Get( Form("PTYRaw/dphi_Nch5") );

    NonFlowSubtractor subTool;
    subTool.setNHar(3);
    subTool.init();

    NonFlowSubtractor subTool_fixLM;
    subTool_fixLM.setNHar(3);
    subTool_fixLM.setFixLM();
    subTool_fixLM.init();

    //--------------------------------------------------
    // ATLAS template fit
    //--------------------------------------------------
    subResult theResult       = subTool.      templateFit(h_correlation_LM, h_correlation_HM);
    subResult theResult_fixLM = subTool_fixLM.templateFit(h_correlation_LM, h_correlation_HM);

    cout << "v22 = " << theResult.getV22SubValue() << " +/- " << theResult.getV22SubError() << endl;
    cout << "MINOS lower error = " << theResult.getCoeffSubLowerError(2) << ", upper error = " << theResult.getCoeffSubUpperError(2) << endl;
    cout << "Fit with fixed LM = " << theResult_fixLM.getV22SubValue() << " +/- " << theResult_fixLM.getV22SubError() << endl;
    if (theResult.getNHar() > 2) {
        cout << endl;
        cout << "v33 = " << theResult.getCoeffSubValue(3) << " +/- " << theResult.getCoeffSubError(3) << endl;
        cout << "MINOS lower error = " << theResult.getCoeffSubLowerError(3) << ", upper error = " << theResult.getCoeffSubUpperError(3) << endl;
        cout << "Fit with fixed LM = " << theResult_fixLM.getCoeffSubValue(3) << " +/- " << theResult_fixLM.getCoeffSubError(3) << endl;
    }

    TCanvas* c1 = new TCanvas("c1","scaling",50,50, 600,700);
    subTool.plotAtlasHM(c1);
    c1->cd();
}
