#include "utility.h"
void plotdSigmadt(TString name="phi", bool veto_ = false, bool PHP_ = false){

	/* Beagle */

	if( (name=="rho" && PHP_)||
	 	(name=="phi" && PHP_)||
	 	(name=="jpsi" && PHP_)) {
	 	
	 	cout << "inconsistent settings between beagle and sartre! "<< endl;
	 	return;

	 }

	TString inputROOT="../rootfiles/beagle_allVMs_w_breakups.root";
	if(PHP_) inputROOT="../rootfiles/beagle_allVMs_w_breakups_PHP.root";
	if(veto_&&!PHP_) inputROOT="../rootfiles/beagle_allVMs_w_breakups_w_vetos.root";
	if(veto_&&PHP_) inputROOT="../rootfiles/beagle_allVMs_w_breakups_w_vetos_PHP.root";
	TFile* file_beagle = new TFile(inputROOT);
	TH1D* t_hat_all = (TH1D*) file_beagle->Get("h_trueT");
	TH1D* h_VM[2][3][5];
	TH1D* h_VM_daughter[2][3][5];

	for(int ibreak=0;ibreak<2;ibreak++){
		for(int ivm=0;ivm<3;ivm++){
			for(int ipro=0;ipro<5;ipro++){
				h_VM[ibreak][ivm][ipro] = (TH1D*) file_beagle->Get(Form("h_VM_%d_%d_%d",ibreak,ivm,ipro));
				h_VM_daughter[ibreak][ivm][ipro] = (TH1D*) file_beagle->Get(Form("h_VM_daughter_%d_%d_%d",ibreak,ivm,ipro));
			}
		}
	}

	/* Sartre */

	TFile* file_sartre = new TFile("../rootfiles/sartre_"+name+"_bnonsat.root");
	TH1D* h_coh_sartre = (TH1D*) file_sartre->Get("hist_t_coherent");
	TH1D* h_incoh_sartre = (TH1D*) file_sartre->Get("hist_t_incoherent");

	TCanvas* c1 = new TCanvas("c1","c1",1,1,600,600);
	gPad->SetLogy(1);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	TH1D* base1 = makeHist("base1", "", "|#it{t} | (GeV^{2})", "d#sigma/d|#it{t} | (nb/GeV^{2}) ", 100,0,0.18,kBlack);
	if(name=="rho"){base1->GetYaxis()->SetRangeUser(1e-1, 1e8);}
	if(name=="phi"){base1->GetYaxis()->SetRangeUser(1e-2, 1e7);}
	if(name=="jpsi"){base1->GetYaxis()->SetRangeUser(1e-3, 1e6);}
	if(name=="rho_photo"){base1->GetYaxis()->SetRangeUser(1e-1, 1e10);}
	if(name=="phi_photo"){base1->GetYaxis()->SetRangeUser(1e-1, 1e9);}
	if(name=="jpsi_photo"){base1->GetYaxis()->SetRangeUser(1e-3, 1e8);}

	base1->GetXaxis()->SetTitleColor(kBlack);
	fixedFontHist1D(base1,1.,1.1);
	base1->GetYaxis()->SetTitleSize(base1->GetYaxis()->GetTitleSize()*1.5);
	base1->GetXaxis()->SetTitleSize(base1->GetXaxis()->GetTitleSize()*1.5);
	base1->GetYaxis()->SetLabelSize(base1->GetYaxis()->GetLabelSize()*1.5);
	base1->GetXaxis()->SetLabelSize(base1->GetXaxis()->GetLabelSize()*1.5);
	base1->GetXaxis()->SetNdivisions(4,4,0);
	base1->GetYaxis()->SetNdivisions(5,5,0);
	base1->Draw();

	measureXsection(name, h_coh_sartre, 1);
	h_coh_sartre->SetLineColor(kBlue);
	h_coh_sartre->SetLineWidth(2);
	h_coh_sartre->Draw(" hist same");

	measureXsection(name, h_incoh_sartre, 1);
	h_incoh_sartre->SetLineColor(kRed);
	h_incoh_sartre->SetLineWidth(2);
	h_incoh_sartre->Draw(" hist same");

	TH1D* h_VM_background = (TH1D*) h_VM[0][vm_index][4]->Clone("h_VM_background");
	h_VM_background->Add(h_VM[1][vm_index][4],+1);
	//adding elastic and dissoc. together.
	measureXsection(name, h_VM_background, 0, t_hat_all->GetEntries(), PHP_);
	h_VM_background->Rebin(2);
	h_VM_background->Scale(1./2);
	h_VM_background->SetMarkerStyle(24);
	h_VM_background->Draw("P SAME");

	TLegend *w5 = new TLegend(0.53,0.74,0.73,0.87);
	w5->SetLineColor(kWhite);
	w5->SetFillColor(0);
	w5->SetTextSize(18);
	w5->SetTextFont(45);
	w5->AddEntry(h_coh_sartre, "Sartre "+legendName+" coherent  ", "PL");
	w5->AddEntry(h_incoh_sartre, "Sartre "+legendName+" incoherent  ", "PL");
	w5->AddEntry(h_VM_background, "BeAGLE "+legendName+" incoherent ", "P");
	w5->Draw("same");

	TLatex* r42 = new TLatex(0.18, 0.91, "eAu 18x110 GeV^{2}");
	r42->SetNDC();
	r42->SetTextSize(22);
	r42->SetTextFont(43);
	r42->SetTextColor(kBlack);
	r42->Draw("same");

	TLatex* r43 = new TLatex(0.6,0.91, "ATHENA #it{Internal}");
	r43->SetNDC();
	r43->SetTextSize(0.04);
	r43->Draw("same");

	TLatex* r44 = new TLatex(0.18, 0.84, "1<Q^{2}<20 GeV^{2}, |#eta_{dau}|<4.0");
	r44->SetNDC();
	r44->SetTextSize(20);
	r44->SetTextFont(43);
	r44->SetTextColor(kBlack);

	TLatex* r44_1 = new TLatex(0.18, 0.84, "Q^{2}<0.2 GeV^{2}, |#eta_{dau}|<4.0");
	r44_1->SetNDC();
	r44_1->SetTextSize(20);
	r44_1->SetTextFont(43);
	r44_1->SetTextColor(kBlack);
	if(PHP_) r44_1->Draw("same");
	else r44->Draw("same");

	if(veto_) c1->Print("../figures/dsigmadt/veto_dsigma_dt_"+name+".pdf");
	else c1->Print("../figures/dsigmadt/dsigma_dt_"+name+".pdf");

}