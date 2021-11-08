#include "RiceStyle.h"
using namespace std;

void plot_DVCS_ed(){

	TFile* file = new TFile("dvcs-d-neutron_output.root");
	TFile* file_base = new TFile("../dvcs-ep/dvcs-ep_output.root");
	
	TH1D* h_t_MC = (TH1D*) file->Get("h_t_MC");
	TH1D* h_t_REC = (TH1D*) file->Get("h_t_REC");
	TH1D* h_t_MC_ep = (TH1D*) file_base->Get("h_t_MC");
	TH1D* h_Q2_sim = (TH1D*) file->Get("h_Q2_sim");
	TH1D* h_Q2_sim_base = (TH1D*) file_base->Get("h_Q2_sim");

	int factor = h_Q2_sim_base->GetEntries()/h_Q2_sim->GetEntries();
	h_t_REC->Scale(factor);

	TH1D* h_Momentum_proton_MC = (TH1D*) file->Get("h_Momentum_proton_MC");
	TH1D* h_Momentum_proton_REC = (TH1D*) file->Get("h_Momentum_proton_REC");
	TH1D* h_Momentum_neutron_MC = (TH1D*) file->Get("h_Momentum_neutron_MC");
	TH1D* h_Momentum_neutron_REC = (TH1D*) file->Get("h_Momentum_neutron_REC");

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

	// for(int i=0;i<h_t_MC->GetNbinsX();i++){
	// 	double bincenter=h_t_MC->GetBinCenter(i+1);
	// 	if(bincenter>0.75){
	// 		h_t_MC->SetBinContent(i+1, 0.);
	// 	}
	// }
	h_t_MC_ep->SetLineColor(kBlack);
	h_t_MC_ep->Draw("hist same");

	h_t_REC->SetFillColorAlpha(kBlue,0.2);
    h_t_REC->SetFillStyle(1001);
	h_t_REC->SetMarkerStyle(24);
	h_t_REC->SetMarkerColor(kBlue);
	h_t_REC->Draw("PE e3 same");

	TLatex* r42 = new TLatex(0.18, 0.91, "ed 10x100 GeV^{2}");
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

	TLatex* r44_2 = new TLatex(0.7, 0.84, "e+d #rightarrow e'+n'+p'+#gamma" );
	r44_2->SetNDC();
	r44_2->SetTextSize(22);
	r44_2->SetTextFont(43);
	r44_2->SetTextColor(kBlack);
	r44_2->Draw("same");

	TLegend *w6 = new TLegend(0.63,0.6,0.8,0.73);
	w6->SetLineColor(kWhite);
	w6->SetFillColor(0);
	w6->SetTextSize(20);
	w6->SetTextFont(45);
	w6->AddEntry(h_t_MC_ep, "EpIC truth ", "L");
	w6->AddEntry(h_t_REC, "ATHENA reco. ", "P");
	w6->AddEntry(h_t_REC, "double tagging", "");
	w6->Draw("same");

	//only plotting phi
	TCanvas* c22 = new TCanvas("c22","c22",1,1,600,600);
	gPad->SetLogy(1);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	gPad->SetRightMargin(0.01);
	TH1D* base22 = makeHist("base22", "", "P_{miss} (GeV/c)", "Counts ", 100,0,0.56,kBlack);
	base22->GetYaxis()->SetRangeUser(4e-1, 1e3);
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

	h_Momentum_proton_MC->SetFillColor(kGray);
    h_Momentum_proton_MC->SetFillStyle(1001);
	h_Momentum_proton_MC->SetMarkerStyle(24);
	h_Momentum_proton_MC->SetMarkerColor(kBlack);
	h_Momentum_proton_MC->SetLineColor(kBlack);

	h_Momentum_proton_MC->Draw("hist F same");
	h_Momentum_proton_REC->Draw("PEsame");

	TLatex* r44_3 = new TLatex(0.7, 0.84, "Proton spectator" );
	r44_3->SetNDC();
	r44_3->SetTextSize(24);
	r44_3->SetTextFont(43);
	r44_3->SetTextColor(kBlack);
	r44_3->Draw("same");

	TLegend *w7 = new TLegend(0.73,0.64,0.95,0.73);
	w7->SetLineColor(kWhite);
	w7->SetFillColor(0);
	w7->SetTextSize(18);
	w7->SetTextFont(45);
	w7->AddEntry(h_Momentum_proton_MC, "EpIC truth ", "F");
	w7->AddEntry(h_Momentum_proton_REC, "ATHENA reco. ", "P");
	w7->Draw("same");
	r43->Draw("same");

	//only plotting phi
	TCanvas* c33 = new TCanvas("c33","c33",1,1,600,600);
	gPad->SetLogy(1);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	gPad->SetRightMargin(0.01);
	TH1D* base33 = makeHist("base33", "", "P_{miss} (GeV/c)", "Counts ", 100,0,1.,kBlack);
	base33->GetYaxis()->SetRangeUser(4e-1, 1e3);
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

	h_Momentum_neutron_MC->SetFillColor(kGray);
    h_Momentum_neutron_MC->SetFillStyle(1001);
	h_Momentum_neutron_MC->SetMarkerStyle(24);
	h_Momentum_neutron_MC->SetMarkerColor(kBlack);
	h_Momentum_neutron_MC->SetLineColor(kBlack);

	h_Momentum_neutron_MC->Draw("hist F same");
	h_Momentum_neutron_REC->Draw("PEsame");

	TLatex* r44_4 = new TLatex(0.73, 0.84, "neutron active" );
	r44_4->SetNDC();
	r44_4->SetTextSize(24);
	r44_4->SetTextFont(43);
	r44_4->SetTextColor(kBlack);
	r44_4->Draw("same");

	r43->Draw("same");

	TLegend *w8 = new TLegend(0.73,0.64,0.95,0.73);
	w8->SetLineColor(kWhite);
	w8->SetFillColor(0);
	w8->SetTextSize(18);
	w8->SetTextFont(45);
	w8->AddEntry(h_Momentum_neutron_MC, "EpIC truth ", "F");
	w8->AddEntry(h_Momentum_neutron_REC, "ATHENA reco. ", "P");
	w8->Draw("same");

	c11->Print("Figure_1.pdf");
	c22->Print("Figure_2.pdf");
	c33->Print("Figure_3.pdf");


}