#include "function.C"
void drawGraph(){

    gStyle->SetOptTitle(0); 
    gStyle->SetTitleAlign(13);
    gStyle->SetTitleX(0.2);

    TH1D * h1Tau10Cut[10];
    TFile *file1 = new TFile("merge_output_t_hat_18on110_ePb.root","READ");
    h1Tau10Cut[0] = (TH1D*)file1->Get("h1_t_hatCut0Newpara");
    h1Tau10Cut[1] = (TH1D*)file1->Get("h1_t_hatCut_main");
    if(!file1) cout<<"There is no file1"<<endl;
    if(!h1Tau10Cut[1]) cout<<"There is no cut1"<<endl;   

    TFile *file2 = new TFile("merge_output_IncoherentVetoing_ePb.root","READ");
    h1Tau10Cut[2] = (TH1D*)file2->Get("h1_t_hat_Cut1");
    h1Tau10Cut[3] = (TH1D*)file2->Get("h1_t_hat_Cut2");
    h1Tau10Cut[4] = (TH1D*)file2->Get("h1_t_hat_Cut3");
    h1Tau10Cut[5] = (TH1D*)file2->Get("h1_t_hat_Cut4");
    h1Tau10Cut[6] = (TH1D*)file2->Get("h1_t_hat_Cut5");
    h1Tau10Cut[7] = (TH1D*)file2->Get("h1_t_hat_Cut6"); 
    //h1Tau10Cut[0]->Draw("l");
    //return;

    TH1D * h1Tau6Cut[10];
    TFile *file3 = new TFile("merge_output_t_hat_18on110_ePb_tau6.root","READ");
    h1Tau6Cut[0] = (TH1D*)file3->Get("h1_t_hatCut0Newpara");
    h1Tau6Cut[1] = (TH1D*)file3->Get("h1_t_hatCut_main");
    if(!file3) cout<<"There is no file3"<<endl;
    if(!h1Tau6Cut[1]) cout<<"There is no cut1"<<endl;   

    TFile *file4 = new TFile("merge_output_IncoherentVetoing_ePb_tau6.root","READ");
    h1Tau6Cut[2] = (TH1D*)file4->Get("h1_t_hat_Cut1");
    h1Tau6Cut[3] = (TH1D*)file4->Get("h1_t_hat_Cut2");
    h1Tau6Cut[4] = (TH1D*)file4->Get("h1_t_hat_Cut3");
    h1Tau6Cut[5] = (TH1D*)file4->Get("h1_t_hat_Cut4");
    h1Tau6Cut[6] = (TH1D*)file4->Get("h1_t_hat_Cut5");
    h1Tau6Cut[7] = (TH1D*)file4->Get("h1_t_hat_Cut6"); 
/*
    for(int i=0; i<8; i++){
        h1Tau6Cut[i]->Scale(h1Tau10Cut[0]->GetEntries() / h1Tau6Cut[0]->GetEntries(), "width");
    }
*/
    TH1D * h1Tau14Cut[10];
    TFile *file5 = new TFile("merge_output_t_hat_18on110_ePb_tau14.root","READ");
    h1Tau14Cut[0] = (TH1D*)file5->Get("h1_t_hatCut0Newpara");
    h1Tau14Cut[1] = (TH1D*)file5->Get("h1_t_hatCut_main");
    if(!file5) cout<<"There is no file5"<<endl;
    if(!h1Tau14Cut[1]) cout<<"There is no cut1"<<endl;   

    TFile *file6 = new TFile("merge_output_IncoherentVetoing_ePb_tau14.root","READ");
    h1Tau14Cut[2] = (TH1D*)file6->Get("h1_t_hat_Cut1");
    h1Tau14Cut[3] = (TH1D*)file6->Get("h1_t_hat_Cut2");
    h1Tau14Cut[4] = (TH1D*)file6->Get("h1_t_hat_Cut3");
    h1Tau14Cut[5] = (TH1D*)file6->Get("h1_t_hat_Cut4");
    h1Tau14Cut[6] = (TH1D*)file6->Get("h1_t_hat_Cut5");
    h1Tau14Cut[7] = (TH1D*)file6->Get("h1_t_hat_Cut6"); 
  
    double scale6 = h1Tau10Cut[0]->GetEntries() / h1Tau6Cut[0]->GetEntries();
    double scale14 = h1Tau10Cut[0]->GetEntries() / h1Tau14Cut[0]->GetEntries();
  
    for(int i=0; i<8; i++){
        h1Tau6Cut[i]->Scale(scale6,"width");
        h1Tau14Cut[i]->Scale(scale14, "width");
    }


   for(int i=0; i<8; i++){
        h1Tau10Cut[i]->Scale(1, "width");
        
   }

   //no Cut 
    const int nBins = h1Tau10Cut[0]->GetNbinsX();
    double x[8][nBins];
    double y[8][nBins];
    double exl[8][nBins];
    double exh[8][nBins];
    double eyl[8][nBins];
    double eyh[8][nBins];
    TGraphAsymmErrors *GraphCut[8];
    for(int j=0; j<8; j++){
        for(int i=0; i<nBins; i++){
            x[j][i] = h1Tau10Cut[j]->GetBinCenter(i+1);
            y[j][i] = h1Tau10Cut[j]->GetBinContent(i+1);
            exl[j][i] = h1Tau10Cut[j]->GetBinWidth(i+1)/2.;
            exh[j][i] = h1Tau10Cut[j]->GetBinWidth(i+1)/2.;
            eyl[j][i] = TMath::Abs(h1Tau6Cut[j]->GetBinContent(i+1) - h1Tau10Cut[j]->GetBinContent(i+1));
            eyh[j][i] = TMath::Abs(h1Tau14Cut[j]->GetBinContent(i+1) - h1Tau10Cut[j]->GetBinContent(i+1));
        }

        GraphCut[j] = new TGraphAsymmErrors(nBins, x[j], y[j], exl[j], exh[j], eyl[j], eyh[j]); 
    //    GraphCut[j] ->Print();
    }

    TCanvas *cEvent = new TCanvas("cEvent", "", 600, 800);
    cEvent->cd();

    //cEvent->SetLogy();
    TPad *c1_1 = new TPad("c1_1", "c1_1", 0., 0.2, 1., 1.);
    c1_1->SetTopMargin(0.1);
    c1_1->SetRightMargin(0.01);
    c1_1->SetLeftMargin(0.13);
    c1_1->SetFillColor(0);
    c1_1->SetBorderMode(0);
    c1_1->SetBorderSize(2);
    c1_1->SetLogy();
    c1_1->SetBottomMargin(0.0);
    c1_1->SetFrameBorderMode(0);
    c1_1->SetFrameBorderMode(0);
    c1_1->Draw();
    c1_1->cd();
    gPad->SetTicks();

     h1Tau10Cut[0]->SetStats(0);
        //h1Tau10Cut[0]->SetTitleOffset(0.1);
        h1Tau10Cut[0]->GetYaxis()->SetTitle("dN/d#it{t} (GeV^{-2})");
        h1Tau10Cut[0]->GetYaxis()->SetRangeUser(2e2, 2e7);
        h1Tau10Cut[0]->GetYaxis()->SetLabelOffset(0.002);
        h1Tau10Cut[0]->GetXaxis()->SetNdivisions(505);
        h1Tau10Cut[0]->GetXaxis()->SetTitle("|#it{t}| (GeV^{2})");
        h1Tau10Cut[0]->GetYaxis()->SetTitleSize(h1Tau10Cut[0]->GetYaxis()->GetTitleSize()*0.9);
        // h1Tau10Cut[0]->GetXaxis()->SetTitleSize(h1Tau10Cut[0]->GetXaxis()->GetTitleSize()*1.5);
        h1Tau10Cut[0]->GetYaxis()->SetTitleOffset(h1Tau10Cut[0]->GetYaxis()->GetTitleOffset()*1.1);
        h1Tau10Cut[0]->GetXaxis()->SetTitleOffset(h1Tau10Cut[0]->GetXaxis()->GetTitleOffset()*2.5);
        h1Tau10Cut[0]->Draw("chist");
        GraphCut[0]->SetFillStyle(1001);
        GraphCut[0]->Draw("e3");
        h1Tau10Cut[0]->Draw("same chist");


    for(Int_t i = 1; i <8 ; i++){
        GraphCut[i]->SetFillStyle(1001);
        h1Tau10Cut[i]->SetLineColor(i+1);
        h1Tau10Cut[i]->SetLineWidth(4);
        if(i!=1 && i!=3 && i!=5){

        GraphCut[i]->Draw("e3");
        h1Tau10Cut[i]->Draw("same chist");
        }
    }
        h1Tau10Cut[0]->SetTitleOffset(0.7,"x");
        h1Tau10Cut[1]->SetLineStyle(1);
        h1Tau10Cut[2]->SetLineStyle(2);
       // h1Tau10Cut[3]->SetLineStyle(4);
        h1Tau10Cut[4]->SetLineStyle(1);
        h1Tau10Cut[5]->SetLineStyle(6);
        h1Tau10Cut[6]->SetLineStyle(1);


        GraphCut[0]->SetFillColorAlpha(1, 0.5);
        GraphCut[0]->SetLineColor(1);
        h1Tau10Cut[1]->SetLineColor(46);
        GraphCut[1]->SetFillColorAlpha(46, 0.5);
        GraphCut[1]->SetLineColor(46);
        //h1Tau10Cut[1]->SetFillColorAlpha(46, 0.5);
        h1Tau10Cut[2]->SetLineColor(kAzure-1);
        GraphCut[2]->SetFillColorAlpha(kAzure-1, 0.5);
        GraphCut[2]->SetLineColor(kAzure-1);
        //h1Tau10Cut[2]->SetFillColorAlpha(kAzure-1, 0.5);
        h1Tau10Cut[3]->SetLineColor(kMagenta+2);       
        GraphCut[3]->SetFillColorAlpha(kMagenta+2, 0.5);
        GraphCut[3]->SetLineColor(kMagenta+2);
        //h1Tau10Cut[3]->SetFillColorAlpha(kMagenta+2, 0.5);
        h1Tau10Cut[4]->SetLineColor(kRed);
        GraphCut[4]->SetFillColorAlpha(kRed-2, 0.5);
        GraphCut[4]->SetLineColor(kRed);
        //h1Tau10Cut[4]->SetFillColorAlpha(kRed-2, 0.5);
        h1Tau10Cut[5]->SetLineColor(5+1);
        GraphCut[5]->SetFillColorAlpha(5+1, 0.5);
        GraphCut[5]->SetLineColor(5+1);
        //h1Tau10Cut[5]->SetFillColorAlpha(5+1, 0.5);
        h1Tau10Cut[6]->SetLineColor(38);
        GraphCut[6]->SetFillColorAlpha(38, 0.5);
        GraphCut[6]->SetLineColor(38);
        //h1Tau10Cut[6]->SetFillColorAlpha(38, 0.5);
        h1Tau10Cut[7]->SetLineColor(kGreen-2);
        GraphCut[7]->SetFillColorAlpha(kGreen-2, 0.5);
        GraphCut[7]->SetLineColor(kGreen-2);
        //h1Tau10Cut[7]->SetFillColorAlpha(kGreen-2, 0.5);

//        h1Tau10Cut[5]->SetMarkerColor(6);
        h1Tau10Cut[0]->SetLineWidth(4);
        h1Tau10Cut[2]->SetLineWidth(4);
        h1Tau10Cut[4]->SetLineWidth(4);
        h1Tau10Cut[5]->SetLineWidth(4);
        h1Tau10Cut[6]->SetLineWidth(4);
        h1Tau10Cut[7]->SetLineWidth(4);

       
        drawLine(0.015, 1205/0.0025, 0.019, 1205/0.0025,2, 1, 1); 
        drawLine(0.015, 954/0.0025, 0.019, 954/0.0025,2, 1, 1); 
        TArrow *ar1 = new TArrow(0.017, 975/0.0025, 0.017, 1150/0.0025,0.003,"<|>");
        ar1->SetLineWidth(2);
        ar1->SetFillColor(kRed);
        ar1->SetLineColor(kRed);
        ar1->Draw();
        drawLatex(0.011, 200000/*600/0.0025*/ , "#bf{1^{st} min.}", 42, 0.03, 1, 0); 

        drawLine(0.053, 90/0.0025, 0.057, 90/0.0025,2, 1, 1); 
        drawLine(0.053, 72/0.0025, 0.057, 72/0.0025,2, 1, 1); 
        TArrow *ar2 = new TArrow(0.055, 74/0.0025, 0.055, 87/0.0025,0.004,"<|>");
        ar2->SetLineWidth(2);
        ar2->SetFillColor(kRed);
        ar2->SetLineColor(kRed);
        ar2->Draw();
        drawLatex(0.049, 15500/*50/0.0025*/, "#bf{2^{nd} min.}", 42, 0.03, 1, 0); 

        drawLine(0.108, 20/0.0025, 0.112,20/0.0025,2, 1, 1); 
        drawLine(0.108, 13/0.0025, 0.112, 13/0.0025,2, 1, 1); 
        TArrow *ar3 = new TArrow(0.11, 13.5/0.0025, 0.11, 19.5/0.0025 ,0.005,"<|>");
        ar3->SetLineWidth(2);
        ar3->SetFillColor(kRed);
        ar3->SetLineColor(kRed);
        ar3->Draw();
        drawLatex(0.104, 3000/*9/0.0025*/, "#bf{3^{rd} min.}", 42, 0.03, 1, 0);  
    
    
        TLegend * leg = new TLegend(0.15,0.02,0.35,0.2,NULL,"brNDC");
    leg->AddEntry(GraphCut[0], "Total", "lf");
    leg->AddEntry(GraphCut[2], "Veto.2", "lf");
     leg->AddEntry(GraphCut[4], "Veto.4", "lf");
    
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.04);
    leg->Draw("same");

    TLegend * leg2 = new TLegend(0.35,0.05,0.55,0.17,NULL,"brNDC");
    //leg2->AddEntry(h1Tau10Cut[5], "veto.5", "l");
    leg2->AddEntry(GraphCut[6], "Veto.6", "lf");
    leg2->AddEntry(GraphCut[7], "Veto.7", "lf");
    
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->SetTextFont(42);
    leg2->SetTextSize(0.04);
    leg2->Draw("same");

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(28);
    latex.SetTextFont(43);
    latex.DrawLatex(0.6, 0.84, "e+Pb #rightarrow e'+J/#psi+X");

    TLatex latexb;
    latexb.SetNDC();
    latexb.SetTextSize(28);
    latexb.SetTextFont(43);
    latexb.DrawLatex(0.15, 0.92, "#bf{BeAGLE}");

    TLatex latexa;
    latexa.SetNDC();
    latexa.SetTextSize(28);
    latexa.SetTextFont(43);                               
    latexa.DrawLatex(0.7, 0.92, "18x110 GeV^{2}");

    cEvent->cd();

    TPad * c1_2 = new TPad("c1_2", "c1_2", 0., 0. , 1., 0.2);  
    c1_2->SetTopMargin(0.0);
    c1_2->SetRightMargin(0.01);
    c1_2->SetLeftMargin(0.13);
    c1_2->SetBottomMargin(0.4);
    c1_2->SetFillColor(0);
    c1_2->SetBorderMode(0);
    c1_2->SetBorderSize(2);
//    c1_2->SetGridx();
//    c1_2->SetGridy();
    c1_2->SetFrameBorderMode(0);
    c1_2->SetFrameBorderMode(0);
    c1_2->Draw();
    c1_2->cd();
    c1_2->SetLogy();
    
    TH1D *ratio[8]={} ;
    ratio[1] = (TH1D*)h1Tau10Cut[1]->Clone("ratio1");
    ratio[2] = (TH1D*)h1Tau10Cut[2]->Clone("ratio2");
    ratio[3] = (TH1D*)h1Tau10Cut[3]->Clone("ratio3");
    ratio[4] = (TH1D*)h1Tau10Cut[4]->Clone("ratio4");
    ratio[5] = (TH1D*)h1Tau10Cut[5]->Clone("ratio5");
    ratio[6] = (TH1D*)h1Tau10Cut[6]->Clone("ratio6");
    ratio[7] = (TH1D*)h1Tau10Cut[7]->Clone("ratio7");
        

        ratio[2]->GetYaxis()->SetTitle("ratio to total");
        ratio[2]->GetYaxis()->SetRangeUser(2e-3, 9);
        ratio[2]->GetYaxis()->SetNdivisions(2, 0, 4);
        ratio[2]->GetYaxis()->SetTitleSize(0.17);
        ratio[2]->GetYaxis()->SetTitleOffset(0.4);
        ratio[2]->GetYaxis()->SetLabelOffset(-0.002);
        ratio[2]->GetYaxis()->SetLabelSize(0.16);
        ratio[2]->GetXaxis()->SetTitle("|#it{t}| (GeV^{2})");
        ratio[2]->GetXaxis()->SetTitleSize(0.2);
        ratio[2]->GetXaxis()->SetLabelSize(0.18);
        //ratio[2]->GetXaxis()->SetLabelSize(ratio[1]->GetXaxis()->GetLabelSize()*1.2);
        ratio[2]->GetXaxis()->SetTitleOffset(0.85);
        ratio[2]->GetXaxis()->SetNdivisions(5,0,5);

        for(int i=1; i<8; i++){

        ratio[i]->Divide(h1Tau10Cut[0]);
            if(i!=1 && i!=3 && i!=5){
                ratio[i]->Draw("same chist");
            }
        }
        ratio[2]->SetStats(0);      
        ratio[2]->SetTitle(0);












}
