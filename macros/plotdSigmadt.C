#include "utility.h"
void plotdSigmadt(TString name="phi", int veto_ = 0, int PHP_ = 0, double minPt_=0.15){

	/* Beagle */

	if( ((name=="rho" && PHP_)||
	 	(name=="phi" && PHP_)||
	 	(name=="jpsi" && PHP_)) ||
	 	((name=="rho_photo" && !PHP_)||
	 	(name=="phi_photo" && !PHP_)||
	 	(name=="jpsi_photo" && !PHP_)) ) {
	 	
	 	cout << "inconsistent settings between beagle and sartre! "<< endl;
	 	return;

	 }

	TString inputROOT=Form("../rootfiles/beagle_output_PHP_%d_veto_%d_minPt_%.2f_smear_0.root",PHP_,veto_,minPt_);
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
	TH1D* h_t_reco[3][3][3][3];
	for(int ibreak=0;ibreak<3;ibreak++){
		for(int ivm=0;ivm<3;ivm++){
			for(int imethod=0;imethod<3;imethod++){
				for(int imass=0;imass<3;imass++){
					h_t_reco[ibreak][ivm][imethod][imass] = (TH1D*) file_beagle->Get(Form("h_t_reco_%d_%d_%d_%d",ibreak,ivm,imethod,imass));
				}
			}
		}
	}

	/* Sartre */

	TFile* file_sartre = new TFile(Form("../rootfiles/sartre_"+name+"_bnonsat_PID_0_minPt_%.2f_smear_0.root",minPt_));
	TH1D* h_coh_sartre = (TH1D*) file_sartre->Get("hist_t_coherent");
	TH1D* h_incoh_sartre = (TH1D*) file_sartre->Get("hist_t_incoherent");
	TH1D* hist_t_afterPhaseSpace_coherent_sartre = (TH1D*) file_sartre->Get("hist_t_afterPhaseSpace_coherent");
	TH1D* hist_t_afterPhaseSpace_incoherent_sartre = (TH1D*) file_sartre->Get("hist_t_afterPhaseSpace_incoherent");

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
	if(name=="phi_photo"){base1->GetYaxis()->SetRangeUser(1e-3, 1e9);}
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

	measureXsection(name, hist_t_afterPhaseSpace_coherent_sartre, 1);
	hist_t_afterPhaseSpace_coherent_sartre->SetLineColor(kBlue);
	hist_t_afterPhaseSpace_coherent_sartre->SetLineStyle(2);
	hist_t_afterPhaseSpace_coherent_sartre->SetLineWidth(2);
	hist_t_afterPhaseSpace_coherent_sartre->Draw(" hist same");

	measureXsection(name, h_incoh_sartre, 1);
	h_incoh_sartre->SetLineColor(kRed);
	h_incoh_sartre->SetLineWidth(2);
	h_incoh_sartre->Draw(" hist same");

	measureXsection(name, hist_t_afterPhaseSpace_incoherent_sartre, 1);
	hist_t_afterPhaseSpace_incoherent_sartre->SetLineColor(kRed);
	hist_t_afterPhaseSpace_incoherent_sartre->SetLineStyle(2);
	hist_t_afterPhaseSpace_incoherent_sartre->SetLineWidth(2);
	hist_t_afterPhaseSpace_incoherent_sartre->Draw(" hist same");

	TH1D* h_VM_background = (TH1D*) h_VM[0][vm_index][4]->Clone("h_VM_background");
	h_VM_background->Add(h_VM[1][vm_index][4],+1);
	//use Method A. only
	int method=0;
	if(PHP_) method=1;
	TH1D* h_VM_background_afterPhaseSpace = (TH1D*) h_t_reco[0][vm_index][method][0]->Clone("h_VM_background_afterPhaseSpace");
	h_VM_background_afterPhaseSpace->Add(h_t_reco[1][vm_index][method][0],+1);

	//adding elastic and dissoc. together.
	measureXsection(name, h_VM_background, 0, t_hat_all->GetEntries(), PHP_, 0);
	h_VM_background->Rebin(2);
	h_VM_background->Scale(1./2);
	h_VM_background->SetMarkerStyle(24);
	h_VM_background->Draw("P SAME");

	measureXsection(name, h_VM_background_afterPhaseSpace, 0, t_hat_all->GetEntries(), PHP_, 1);
	h_VM_background_afterPhaseSpace->SetMarkerStyle(25);
	h_VM_background_afterPhaseSpace->Draw("P SAME");

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

	TLatex* r44 = new TLatex(0.18, 0.84, "1<Q^{2}<20 GeV^{2}, |y_{VM}|<4.0");
	r44->SetNDC();
	r44->SetTextSize(20);
	r44->SetTextFont(43);
	r44->SetTextColor(kBlack);

	TLatex* r44_1 = new TLatex(0.18, 0.84, "Q^{2}<0.2 GeV^{2}, |y_{VM}|<4.0");
	r44_1->SetNDC();
	r44_1->SetTextSize(20);
	r44_1->SetTextFont(43);
	r44_1->SetTextColor(kBlack);
	if(PHP_) r44_1->Draw("same");
	else r44->Draw("same");

	TLegend *w6 = new TLegend(0.18,0.18,0.45,0.28);
	w6->SetLineColor(kWhite);
	w6->SetFillColor(0);
	w6->SetTextSize(15);
	w6->SetTextFont(45);
	w6->AddEntry(hist_t_afterPhaseSpace_coherent_sartre, "Sartre w. daug. cut "+legendName+" coherent  ", "PL");
	w6->AddEntry(hist_t_afterPhaseSpace_incoherent_sartre, "Sartre w. daug. cut "+legendName+" incoherent  ", "PL");
	w6->AddEntry(h_VM_background_afterPhaseSpace, "BeAGLE w. daug. cut "+legendName+" incoherent ", "P");
	w6->Draw("same");

	// cout << "beagle: " << h_VM_background->Integral(h_VM_background->FindBin(0.02),h_VM_background->FindBin(0.2),"width") << endl;
	// cout << "sartre: " << h_incoh_sartre->Integral(h_incoh_sartre->FindBin(0.02),h_incoh_sartre->FindBin(0.2),"width") << endl;

	// if(veto_) c1->Print("../figures/dsigmadt_2/veto_dsigma_dt_"+name+".pdf");
	// else c1->Print("../figures/dsigmadt_2/dsigma_dt_"+name+".pdf");

}