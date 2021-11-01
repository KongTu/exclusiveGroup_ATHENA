#include "RiceStyle.h"
using namespace std;

void plot_DVCS_ep(){

	TFile* file = new TFile("dvcs-ep_output.root");
	
	TH1D* h_t_MC = (TH1D*) file->Get("h_t_MC");
	TH1D* h_t_REC = (TH1D*) file->Get("h_t_REC");

	TH1D* h_Angle_gamma_MC = (TH1D*) file->Get("h_Angle_gamma_MC");

	TH1D* h_Eta_gamma_MC_match = (TH1D*) file->Get("h_Eta_gamma_MC_match");
	TH1D* h_Eta_gamma_MC = (TH1D*) file->Get("h_Eta_gamma_MC");

	TH1D* h_Eta_gamma_REC = (TH1D*) file->Get("h_Eta_gamma_REC");
	TH1D* h_Eta_gamma_REC_not_match = (TH1D*) file->Get("h_Eta_gamma_REC_not_match");

	TH1D* h_Eta_scatElec_MC_match = (TH1D*) file->Get("h_Eta_scatElec_MC_match");
	TH1D* h_Eta_scatElec_MC = (TH1D*) file->Get("h_Eta_scatElec_MC");

	//only plotting phi
	TCanvas* c11 = new TCanvas("c11","c11",1,1,600,600);
	gPad->SetLogy(1);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	gPad->SetRightMargin(0.01);
	TH1D* base11 = makeHist("base11", "", "|#it{t} | (GeV^{2})", "dN / d|#it{t} |  (GeV^{-2}) ", 100,0,2,kBlack);
	base11->GetYaxis()->SetRangeUser(4e-2, 1e5);
	base11->GetXaxis()->SetTitleColor(kBlack);
	// TGaxis::SetMaxDigits(3);
	fixedFontHist1D(base11,1.2,1.4);
	base11->GetYaxis()->SetTitleSize(base11->GetYaxis()->GetTitleSize()*1.3);
	base11->GetXaxis()->SetTitleSize(base11->GetXaxis()->GetTitleSize()*1.3);
	base11->GetYaxis()->SetLabelSize(base11->GetYaxis()->GetLabelSize()*1.5);
	base11->GetXaxis()->SetLabelSize(base11->GetXaxis()->GetLabelSize()*1.5);
	base11->GetXaxis()->SetNdivisions(4,4,0);
	base11->GetYaxis()->SetNdivisions(5,5,0);
	base11->Draw();

	h_t_MC->SetLineColor(kBlack);
	h_t_MC->Draw("hist same");

	h_t_REC->SetFillColorAlpha(kBlue,0.2);
    h_t_REC->SetFillStyle(1001);
	h_t_REC->SetMarkerStyle(24);
	h_t_REC->SetMarkerColor(kBlue);
	h_t_REC->Draw("PE e3 same");

	TLatex* r42 = new TLatex(0.18, 0.91, "ep 10x100 GeV^{2}");
	r42->SetNDC();
	r42->SetTextSize(22);
	r42->SetTextFont(43);
	r42->SetTextColor(kBlack);
	r42->Draw("same");

	TLatex* r43 = new TLatex(0.83,0.91, "ATHENA");
	r43->SetNDC();
	r43->SetTextSize(0.04);
	r43->Draw("same");

	TLatex* r44 = new TLatex(0.18, 0.84, "Q^{2} > 1 GeV^{2}");
	r44->SetNDC();
	r44->SetTextSize(20);
	r44->SetTextFont(43);
	r44->SetTextColor(kBlack);
	r44->Draw("same");

	TLatex* r44_2 = new TLatex(0.75, 0.84, "e+p #rightarrow e'+p'+#gamma" );
	r44_2->SetNDC();
	r44_2->SetTextSize(22);
	r44_2->SetTextFont(43);
	r44_2->SetTextColor(kBlack);
	r44_2->Draw("same");

	TLegend *w6 = new TLegend(0.23,0.23,0.5,0.33);
	w6->SetLineColor(kWhite);
	w6->SetFillColor(0);
	w6->SetTextSize(20);
	w6->SetTextFont(45);
	w6->AddEntry(h_t_MC, "EpIC truth ", "L");
	w6->AddEntry(h_t_REC, "ATHENA reco", "P");
	w6->Draw("same");


	//only plotting phi
	TCanvas* c22 = new TCanvas("c22","c22",1,1,600,600);
	gPad->SetLogy(1);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	gPad->SetRightMargin(0.01);
	TH1D* base22 = makeHist("base22", "", "3D photon angular resolution (rad)", "counts ", 100,0,0.18,kBlack);
	base22->GetYaxis()->SetRangeUser(1, 1e5);
	base22->GetXaxis()->SetTitleColor(kBlack);
	// TGaxis::SetMaxDigits(3);
	fixedFontHist1D(base22,1.2,1.4);
	base22->GetYaxis()->SetTitleSize(base22->GetYaxis()->GetTitleSize()*1.3);
	base22->GetXaxis()->SetTitleSize(base22->GetXaxis()->GetTitleSize()*1.3);
	base22->GetYaxis()->SetLabelSize(base22->GetYaxis()->GetLabelSize()*1.5);
	base22->GetXaxis()->SetLabelSize(base22->GetXaxis()->GetLabelSize()*1.5);
	base22->GetXaxis()->SetNdivisions(4,4,0);
	base22->GetYaxis()->SetNdivisions(5,5,0);
	base22->Draw();

	h_Angle_gamma_MC->Draw("PEsame");
	r42->Draw("same");
	r43->Draw("same");
	r44->Draw("same");
	r44_2->Draw("same");

	//only plotting phi
	TCanvas* c33 = new TCanvas("c33","c33",1,1,600,600);
	gPad->SetLogy(0);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	gPad->SetRightMargin(0.01);
	TH1D* base33 = makeHist("base33", "", "#eta (photon)", "Photon reco. ", 100,-4,4,kBlack);
	base33->GetYaxis()->SetRangeUser(0.001, 1.2);
	base33->GetXaxis()->SetTitleColor(kBlack);
	// TGaxis::SetMaxDigits(3);
	fixedFontHist1D(base33,1.2,1.4);
	base33->GetYaxis()->SetTitleSize(base33->GetYaxis()->GetTitleSize()*1.3);
	base33->GetXaxis()->SetTitleSize(base33->GetXaxis()->GetTitleSize()*1.3);
	base33->GetYaxis()->SetLabelSize(base33->GetYaxis()->GetLabelSize()*1.5);
	base33->GetXaxis()->SetLabelSize(base33->GetXaxis()->GetLabelSize()*1.5);
	base33->GetXaxis()->SetNdivisions(4,4,0);
	base33->GetYaxis()->SetNdivisions(5,5,0);
	base33->Draw();

	TGraphAsymmErrors * eff[5];
   	eff[0] = new TGraphAsymmErrors();
   	eff[0]->Divide(h_Eta_gamma_MC_match,h_Eta_gamma_MC);
   	eff[0]->SetMarkerStyle(20);
	eff[0]->Draw("PEsame");

	TH1D* ratio = (TH1D*) h_Eta_gamma_REC_not_match->Clone("ratio");
	for(int i=0;i<h_Eta_gamma_REC_not_match->GetNbinsX();i++){

		double value = h_Eta_gamma_REC_not_match->GetBinContent(i+1);
		double value2 = h_Eta_gamma_REC->GetBinContent(i+1);

		if(! ((value/value2) >0.0 && (value/value2) < 1.0) ){
			h_Eta_gamma_REC->SetBinContent(i+1,1e9);
		}
		value2 = h_Eta_gamma_REC->GetBinContent(i+1);
		ratio->SetBinContent(i+1, value/value2);
		ratio->SetBinError(i+1, h_Eta_gamma_REC_not_match->GetBinError(i+1)/value2);
 	}
	
   	ratio->SetMarkerStyle(24);
	ratio->Draw("PEsame");

	TLegend *w7 = new TLegend(0.43,0.23,0.6,0.33);
	w7->SetLineColor(kWhite);
	w7->SetFillColor(0);
	w7->SetTextSize(18);
	w7->SetTextFont(45);
	w7->AddEntry(eff[0], "Efficiency ", "P");
	w7->AddEntry(ratio, "Fake rate", "P");
	w7->Draw("same");


	TLine* l1 = new TLine(-4,1,4,1);
	l1->SetLineStyle(2);
	l1->SetLineWidth(2);
	l1->SetLineColor(kBlue);
	l1->Draw("same");

	r42->Draw("same");
	r43->Draw("same");
	r44->Draw("same");
	r44_2->Draw("same");

	TCanvas* c44 = new TCanvas("c44","c44",1,1,600,600);
	gPad->SetLogy(0);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	gPad->SetRightMargin(0.01);
	TH1D* base44 = makeHist("base44", "", "#eta (e')", "e' reco. ", 100,-4,4,kBlack);
	base44->GetYaxis()->SetRangeUser(0.001, 1.2);
	base44->GetXaxis()->SetTitleColor(kBlack);
	// TGaxis::SetMaxDigits(3);
	fixedFontHist1D(base44,1.2,1.4);
	base44->GetYaxis()->SetTitleSize(base44->GetYaxis()->GetTitleSize()*1.3);
	base44->GetXaxis()->SetTitleSize(base44->GetXaxis()->GetTitleSize()*1.3);
	base44->GetYaxis()->SetLabelSize(base44->GetYaxis()->GetLabelSize()*1.5);
	base44->GetXaxis()->SetLabelSize(base44->GetXaxis()->GetLabelSize()*1.5);
	base44->GetXaxis()->SetNdivisions(4,4,0);
	base44->GetYaxis()->SetNdivisions(5,5,0);
	base44->Draw();

	eff[1] = new TGraphAsymmErrors();
   	eff[1]->Divide(h_Eta_scatElec_MC_match,h_Eta_scatElec_MC);
   	eff[1]->SetMarkerStyle(20);
	eff[1]->Draw("PEsame");

	l1->Draw("same");
	r42->Draw("same");
	r43->Draw("same");
	r44->Draw("same");
	r44_2->Draw("same");

	c11->Print("figure_1.pdf");
	c22->Print("figure_2.pdf");
	c33->Print("figure_3.pdf");
	c44->Print("figure_4.pdf");


}