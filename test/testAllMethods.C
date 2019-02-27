#include "../NonFlowSubtractor.C"
void testAllMethods() {

    bool _runFourier = true;
    bool _runCms   = true;
    bool _runAtlas = true;

    TFile* fin = new TFile(Form("./correlationHistTest.root"),"READ");
    const int _nhar = 4;

    // retrieve the binning information from the input file
    TH1F* hNch_out_binning = (TH1F*) fin->Get("hNch_out_binning");
    TH1F* hpt_out_binning  = (TH1F*) fin->Get("hpt_out_binning");
    const int Nbins_out_Nch = hNch_out_binning->GetXaxis()->GetNbins();
    const int Nbins_out_pt  = hpt_out_binning->GetXaxis()->GetNbins();
    const double* dNch_out_range = hNch_out_binning->GetXaxis()->GetXbins()->GetArray();;
    const double* dpt_out_range  = hpt_out_binning->GetXaxis()->GetXbins()->GetArray();;

    TFile* fout = new TFile("_result_testAllMethods.root","RECREATE");

    TH1F* hNch_pedG = new TH1F(Form("hNch_pedG"),"", Nbins_out_Nch, dNch_out_range);
    TH1F* hNch_rho = new TH1F(Form("hNch_rho"),"", Nbins_out_Nch, dNch_out_range);
    TH1F* hNch_coef[_nhar];
    TH1F* hNch_coef_sub[_nhar];
    TH1F* hNch_coef_subImp[_nhar];
    for (int ihar=0; ihar<_nhar; ihar++){
        hNch_coef[ihar]     = new TH1F(Form("hNch_coef%d",ihar+1),"", Nbins_out_Nch, dNch_out_range);
        hNch_coef_sub[ihar] = new TH1F(Form("hNch_coef%d_sub",ihar+1),"", Nbins_out_Nch, dNch_out_range);
        hNch_coef_subImp[ihar] = new TH1F(Form("hNch_coef%d_subImp",ihar+1),"", Nbins_out_Nch, dNch_out_range);
    }

    // fourier method
    if (_runFourier) {
        TH1F* h_correlation_LM  = (TH1F*)fin->Get( Form("PTYRaw/dphi_Nch_ref") );
        for (int i = 2; i < Nbins_out_Nch; i++) {

            TH1F* h_correlation  = (TH1F*)fin->Get( Form("PTYRaw/dphi_Nch%d", i) );
            NonFlowSubtractor subF;
            subF.init();
            subF.setNHar(_nhar);
            subResult resultF = subF.fourierFit(h_correlation, h_correlation_LM);

            hNch_pedG->SetBinContent(i+1, resultF.getPedstalValue());
            hNch_pedG->SetBinError  (i+1, resultF.getPedstalError());
            hNch_rho->SetBinContent(i+1, resultF.getRhoValue());
            hNch_rho->SetBinError  (i+1, resultF.getRhoError());

            for (int ihar=0; ihar<_nhar; ihar++){
                hNch_coef[ihar]->SetBinContent(i+1, resultF.getCoeffRawValue(ihar+1));
                hNch_coef[ihar]->SetBinError  (i+1, resultF.getCoeffRawError(ihar+1));
            }
        }

        // save the output to a root file 
        fout->mkdir("fourier");
        fout->cd("fourier");
        hNch_pedG->Write();
        hNch_rho->Write();
        for (int ihar=0; ihar<_nhar; ihar++){
            hNch_coef[ihar]->Write();
        }
    }

    // CMS peripheral subtraction
    if (_runCms) {
        hNch_pedG->Reset();
        hNch_rho->Reset();
        for (int ihar=0; ihar<_nhar; ihar++){
            hNch_coef_sub[ihar]->Reset();
        }
        TFile* fin2 = new TFile(Form("./correlationHistTest_SR.root"),"READ");
        TH1F* h_correlation_LR_LM  = (TH1F*)fin->Get( Form("PTYRaw/dphi_Nch_ref") );
        TH1F* h_correlation_SR_LM  = (TH1F*)fin2->Get( Form("PTYRaw/dphi_Nch_ref") );
        h_correlation_LR_LM->Print();
        h_correlation_SR_LM->Print();
        for (int i = 2; i < Nbins_out_Nch; i++) {
            TH1F* h_correlation_LR  = (TH1F*)fin->Get( Form("PTYRaw/dphi_Nch%d", i) );
            TH1F* h_correlation_SR  = (TH1F*)fin2->Get( Form("PTYRaw/dphi_Nch%d", i) );

            NonFlowSubtractor subC;
            subC.setNHar(_nhar);
            subC.init();
            subResult resultC = subC.periphSub(h_correlation_SR_LM, h_correlation_LR_LM, h_correlation_SR, h_correlation_LR);

            hNch_pedG->SetBinContent(i+1, resultC.getPedstalValue());
            hNch_pedG->SetBinError  (i+1, resultC.getPedstalError());
            hNch_rho->SetBinContent(i+1, resultC.getRhoValue());
            hNch_rho->SetBinError  (i+1, resultC.getRhoError());

            for (int ihar=0; ihar<_nhar; ihar++){
                hNch_coef_sub[ihar]->SetBinContent(i+1, resultC.getCoeffSubValue(ihar+1));
                hNch_coef_sub[ihar]->SetBinError  (i+1, resultC.getCoeffSubError(ihar+1));
            }
        }

        // save the output to a root file 
        fout->mkdir("cms");
        fout->cd("cms");
        hNch_pedG->Write();
        hNch_rho->Write();
        for (int ihar=0; ihar<_nhar; ihar++){
            hNch_coef_sub[ihar]->Write();
        }
    }

    //TCanvas* c1 = new TCanvas("c1","",600,600);
    // ATLAS template fit
    if (_runAtlas) {
        hNch_pedG->Reset();
        hNch_rho->Reset();
        for (int ihar=0; ihar<_nhar; ihar++){
            hNch_coef_sub[ihar]->Reset();
        }

        TH1F* h_correlation_LM  = (TH1F*)fin->Get( Form("PTYRaw/dphi_Nch_ref") );
        TH1F* h_correlation_LM2  = (TH1F*)fin->Get( Form("PTYRaw/dphi_Nch2") );
        for (int i = 2; i < Nbins_out_Nch; i++) {
            TH1F* h_correlation_HM  = (TH1F*)fin->Get( Form("PTYRaw/dphi_Nch%d", i) );

            NonFlowSubtractor subA;
            subA.setNHar(_nhar);
            subA.init();
            subResult resultA = subA.templateFit(h_correlation_LM, h_correlation_HM, h_correlation_LM2);

            hNch_pedG->SetBinContent(i+1, resultA.getPedstalValue());
            hNch_pedG->SetBinError  (i+1, resultA.getPedstalError());
            hNch_rho->SetBinContent(i+1, resultA.getRhoValue());
            hNch_rho->SetBinError  (i+1, resultA.getRhoError());

            for (int ihar=0; ihar<_nhar; ihar++){
                hNch_coef_sub[ihar]->SetBinContent(i+1, resultA.getCoeffSubValue(ihar+1));
                hNch_coef_sub[ihar]->SetBinError  (i+1, resultA.getCoeffSubError(ihar+1));
                hNch_coef_subImp[ihar]->SetBinContent(i+1, resultA.getCoeffSubImpValue(ihar+1));
                hNch_coef_subImp[ihar]->SetBinError  (i+1, resultA.getCoeffSubImpError(ihar+1));
            }
            //subA.plotAtlasHM(c1);
        }

        // save the output to a root file 
        fout->mkdir("atlas");
        fout->cd("atlas");
        hNch_pedG->Write();
        hNch_rho->Write();
        for (int ihar=0; ihar<_nhar; ihar++){
            hNch_coef_sub[ihar]->Write();
            hNch_coef_subImp[ihar]->Write();
        }
    }
}
