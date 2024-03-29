#include "../NonFlowSubtractor.C"
void testSimple() {

    TFile* fin = new TFile(Form("./correlationHistTest.root"),"READ");

    TH1F* h_correlation_LM  = (TH1F*)fin->Get( Form("PTYRaw/dphi_Nch_ref") );
    TH1F* h_correlation_HM  = (TH1F*)fin->Get( Form("PTYRaw/dphi_Nch6") );

    // Correlation for the second lowest LM bin 
    // will be used for correction
    TH1F* h_correlation_LM2 = (TH1F*)fin->Get( Form("PTYRaw/dphi_Nch2") );

    NonFlowSubtractor subTool;
    // change the fitter configure before call init() here
    //subTool.setAtlasFixedC3();
    //subTool.setAtlasFixedC4();
    subTool.setDebug();
    //subTool.setZYAM();
    subTool.init();

    //--------------------------------------------------
    // ATLAS template fit
    //--------------------------------------------------
    //subResult theResult = subTool.templateFit(h_correlation_LM, h_correlation_HM, h_correlation_LM2);
    subResult theResult = subTool.templateFit(h_correlation_LM, h_correlation_HM);
    //subResult theResult = subTool.referenceFit(h_correlation_LM, h_correlation_HM);
    //subResult theResult = subTool.templateHistFit(h_correlation_LM, h_correlation_HM, h_correlation_LM2);
    //subResult theResult = subTool.templateHistFit(h_correlation_LM, h_correlation_HM);

    //--------------------------------------------------
    // Test for using symmetrized dphi correlation
    // OK from first look
    // need further confirmation
    //--------------------------------------------------
    //TH1* hsym_LM = Symmetrize(h_correlation_LM);
    //hsym_LM->SetName("hs_LM");
    //TH1* hsym_HM = Symmetrize(h_correlation_HM);
    //hsym_HM->SetName("hs_HM");
    //subResult theResult = subTool.templateFit(hsym_LM, hsym_HM);

    /*
    */
    //--------------------------------------------------
    // Access the fitted results and plots
    //--------------------------------------------------
    cout << "v22 = " << theResult.getV22SubValue() << " +/- " << theResult.getV22SubError() << endl;
    cout << "Improved v22 = " << theResult.getV22SubImpValue() << " +/- " << theResult.getV22SubImpError() << endl;

    //TCanvas* c1 = new TCanvas("c1","scaling",50,50, 600,600);
    TCanvas* c1 = new TCanvas("c1","scaling",50,50, 600,800);
    subTool.plotAtlas3pHM(c1);
    //subTool.plotAtlasHistSubHM(c1);
    c1->cd();
    c1->SaveAs("test.pdf");
    // add addtional lengends
    
    TCanvas* c2 = new TCanvas("c2","scaling",50,50, 600,600);
    subTool.plotAtlasLM(c2);
    c2->cd();
    // add addtional lengends

}
