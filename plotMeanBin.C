
// ~~ POSTPROCESSING ~~

// --------------------------------------------------------------------------------------
// Module for plotting Data in certain multiplicity bins
//
// --------------------------------------------------------------------------------------

// == Includes ==

#include "Plotting.h"
#include "postproc.h"

// == Namespace ==

using namespace std;

// -- ++ Helpful Functions ++ ----------------------------------------------------------

Double_t Gauss(Double_t *var, Double_t *par){

    Double_t x = var[0];
    return 1./(par[0]*TMath::Sqrt(TMath::TwoPi())) *
                                            exp(-0.5*TMath::Power((x-par[1])/par[0],2));

}

// --------------------------------------------------------------------------------------
// postprocessing
// --------------------------------------------------------------------------------------

void plotMeanBin() {

    // ------ Define Program Specific Parameters ----------------------------------------

    Double_t energy = 5.02;

    enum Plot {raw, fit, gauss};
    Plot plot = gauss;

    // --------- Read in all Files ------------------------------------------------------

    TFile* inData = new TFile("./Daten/ptSpectraUnfolded.root", "READ");

    TH2D* data = (TH2D*) inData->
            FindObjectAny(Form("multPtUnfolded_pp_%dTeV_eta08", (int)energy));

    if(!data) cout << "Error opening histogram data " << endl;

    // ------------- Get Correct Bins ----------------------------------------------------

    TH1D *data10 = data->ProjectionY("data10", 10, 10);
    TH1D *data25 = data->ProjectionY("data25", 25, 25);
    TH1D *data40 = data->ProjectionY("data40", 40, 40);

    for(int bin = 1; bin <= data10->GetXaxis()->GetNbins(); bin++){

            data10->SetBinContent(bin, data10->GetBinContent(bin) *
                                  data10->GetBinCenter(bin));

        }

    for(int bin = 1; bin <= data25->GetXaxis()->GetNbins(); bin++){

         data25->SetBinContent(bin, data25->GetBinContent(bin) *
                               data25->GetBinCenter(bin));

        }

    for(int bin = 1; bin <= data40->GetXaxis()->GetNbins(); bin++){

        data40->SetBinContent(bin, data40->GetBinContent(bin) *
                              data40->GetBinCenter(bin));

        }


    SetHistogramProperties(data10, "#it{p}_{T} (GeV/#it{c})",
                           "1/N_{ev} d#it{N}_{ch}/d#it{p}_{T}d#it{#eta}", 0);
    data10->GetYaxis()->SetRangeUser(0.0000001, 10);

    SetHistogramProperties(data25, "", "", 1);
    data25->GetYaxis()->SetRangeUser(0.0000001, 10);

    SetHistogramProperties(data40, "", "", 2);
    data40->GetYaxis()->SetRangeUser(0.0000001, 10);

    // ------------- Add Lines corresponding to mean values -----------------------------

    Double_t mean[3];
    Double_t error;

    GetMoment(data10, 1, mean[0], error, 0.15, 9.99, false);
    TLine* mean10 = new TLine(mean[0], 0.018, mean[0], 0.067);
    TString string10 = Form("#it{N}_{ch,ev} = 10, <#it{p}_{T}>_{10} = %.2lf GeV/#it{c}", mean[0]);

    GetMoment(data25, 1, mean[1], error, 0.15, 9.99, false);
    TLine* mean25 = new TLine(mean[1], 0.003, mean[1], 0.02);
    TString string25 = Form("#it{N}_{ch,ev} = 25, <#it{p}_{T}>_{25} = %.2lf GeV/#it{c}", mean[1]);

    GetMoment(data40, 1, mean[2], error, 0.15, 9.99, false);
    TLine* mean40 = new TLine(mean[2], 0.0005, mean[2], 0.0035);
    TString string40 = Form("#it{N}_{ch,ev} = 40, <#it{p}_{T}>_{40} = %.2lf GeV/#it{c}", mean[2]);

    mean10->SetLineWidth(3);
    mean10->SetLineColor(kRed+1);

    mean25->SetLineWidth(3);
    mean25->SetLineColor(kGreen-1);

    mean40->SetLineWidth(3);
    mean40->SetLineColor(kBlue+1);

    // ------------- Fit Exponential to Tails to evaluate slope -------------------------

    TF1 *slope = new TF1("slope", "[0]*x^(-[1])", 0, 10);
    slope->SetParameters(1, 4);

    data10->Fit("slope", "NM", "", 3, 10);
    Double_t slope10 = slope->GetParameter(1);
    Double_t amp10 = slope->GetParameter(0);

    data25->Fit("slope", "NM", "", 3, 10);
    Double_t slope25 = slope->GetParameter(1);
    Double_t amp25 = slope->GetParameter(0);

    data40->Fit("slope", "NM", "", 3, 10);
    Double_t slope40 = slope->GetParameter(1);
    Double_t amp40 = slope->GetParameter(0);


    TLine *fitRange = new TLine(3, 0.000000104, 10, 0.000000104);
    fitRange->SetLineWidth(5);
    fitRange->SetLineColor(kCyan+1);

    // ---------- Add Gaussian ----------------------------------------------------------

    Double_t variance[3];

    GetCentralMoment(data10, 2, variance[0], error, 0.15, 9.99, false);
    GetCentralMoment(data25, 2, variance[1], error, 0.15, 9.99, false);
    GetCentralMoment(data40, 2, variance[2], error, 0.15, 9.99, false);

    TF1* gaussian = new TF1("gaussian", &Gauss, 0, 10, 2);
    TF1 *gaussian10, *gaussian25, *gaussian40;

    gaussian->SetParameters(mean[0], variance[0]);


    //gaussian->SetLineWidth(3);
    gaussian->SetMarkerColor(kRed+1);
    gaussian->SetLineColor(kRed+1);
    gaussian->SetMarkerStyle(kDiamond);

    gaussian10 = (TF1*)gaussian->Clone("gaussian10");



    //gaussian->SetParameters(mean[1], variance[1]);
    //gaussian25 = (TF1*)gaussian->Clone("gaussian25");
    gaussian25 = 0;
    //gaussian25->SetLineWidth(3);
    //gaussian25->SetLineColor(kGreen+3);

    //gaussian->SetParameters(mean[2], variance[2]);
    //gaussian40 = (TF1*)gaussian->Clone("gaussian40");
    gaussian40 = 0;
    //gaussian40->SetLineWidth(3);
    //gaussian40->SetLineColor(kBlue+1);


    // ------------ Plotting ------------------------------------------------------------

        // == Legends ==

        TLegend *lInfo = new TLegend(0.125, 0.740165, 0.367574, 0.89387);

        lInfo->SetTextFont(42);
        lInfo->SetTextSize(0.03);
        lInfo->SetBorderSize(0);

        lInfo->AddEntry((TObject*)0x0,
               Form("pp collisions at #sqrt{s} = %.2lf TeV", energy), "");
        lInfo->AddEntry((TObject*)0x0, "charged particles", "");
        lInfo->AddEntry((TObject*)0x0, "-0.8 < #eta < 0.8", "");
        lInfo->AddEntry((TObject*)0x0, "0.15 (GeV/#it{c}) < #it{p}_{T} < 10 (GeV/#it{c})", "");


        TLegend *lBins = new TLegend(0.162129, 0.133577, 0.361386, 0.244282);

        lBins->SetTextFont(42);
        lBins->SetTextSize(0.03);
        lBins->SetBorderSize(0);

        lBins->AddEntry(data10, string10, "lp");
        lBins->AddEntry(data25, string25, "lp");
        lBins->AddEntry(data40, string40, "lp");

        TLegend *lFit = new TLegend(0.157178, 0.122598, 0.355198, 0.287283);

        lFit->SetTextFont(42);
        lFit->SetTextSize(0.03);
        lFit->SetBorderSize(0);

        lFit->AddEntry((TObject*)0x0, "Parametrisierung: a#upointx^{-b}", "");
        lFit->AddEntry(data10, Form("#it{N}_{ev,ch} = 10, a_{10} = %.3lf, b_{10} = %.2lf",
                                    amp10, slope10), "lp");
        lFit->AddEntry(data25, Form("#it{N}_{ev,ch} = 25, a_{25} = %.3lf, b_{25} = %.2lf",
                                    amp25, slope25), "lp");
        lFit->AddEntry(data40, Form("#it{N}_{ev,ch} = 40, a_{40} = %.3lf, b_{40} = %.2lf",
                                    amp40, slope40), "lp");
        lFit->AddEntry(fitRange, "Parametrisierungsbereich", "l");


        // == TObject Arrays ==

        TObjArray *hist = new TObjArray();
        TObjArray *hist1 = new TObjArray();
        TObjArray *hist2 = new TObjArray();
        TObjArray *hist3 = new TObjArray();

        TCanvas *output, *output1, *output2, *output3;

        if (plot == gauss){

            hist1->Add(data10);
            hist1->Add(mean10);

            hist1->Add(gaussian10);
            hist1->Add(lInfo);

            output1 = makeCanvas(hist1, 0, "NoTime LogX LogY", 0, 0);

            data10->GetXaxis()->SetTitleOffset(1.5);
            output1->SetBottomMargin(0.11);

            output1->Update();
            output1->SaveAs("Meeting/BinsGauss10.png");
            output1->SaveAs("Meeting/BinsGauss10.root");


            hist2->Add(data25);
            hist2->Add(mean25);

            hist2->Add(gaussian25);
            hist2->Add(lInfo);

            output2 = makeCanvas(hist2, 0, "NoTime LogX LogY", 0, 0);

            data25->GetXaxis()->SetTitleOffset(1.5);
            output2->SetBottomMargin(0.11);

            output2->Update();
            output2->SaveAs("Meeting/BinsGauss25.png");
            output2->SaveAs("Meeting/BinsGauss25.root");


            hist3->Add(data40);
            hist3->Add(mean40);

            hist3->Add(lInfo);
            hist3->Add(gaussian40);

            output3 = makeCanvas(hist3, 0, "NoTime LogX LogY", 0, 0);
            data40->GetXaxis()->SetTitleOffset(1.5);
            output3->SetBottomMargin(0.11);

            output3->Update();

            output3->SaveAs("Meeting/BinsGauss40.png");
            output3->SaveAs("Meeting/BinsGauss40.root");

        }

        else{

            hist->Add(data10);
            hist->Add(data25);
            hist->Add(data40);

            hist->Add(mean10);
            hist->Add(mean25);
            hist->Add(mean40);

            hist->Add(lInfo);

            switch(plot){

                case raw:   hist->Add(lBins); break;

                case fit:   hist->Add(lFit);
                            hist->Add(fitRange); break;

                case gauss: break;

            }


            output = makeCanvas(hist, 0, "NoTime LogX LogY", 0, 0);

            switch(plot){

                case raw:   output->SaveAs("Meeting/Bins.png");
                            output->SaveAs("Meeting/Bins.root");
                            break;

                case fit:   output->SaveAs("Meeting/BinsFit.png");
                            output->SaveAs("Meeting/BinsFit.root");
                            break;

                case gauss: break;

            }

        }



    // ----- Cleaning Up ----------------------------------------------------------------

    inData->Close();



}
