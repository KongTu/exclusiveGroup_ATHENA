#include "utility.h"
void plotVM(TString name="phi", bool veto_ = false, bool PHP_ = false){

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
	//VM histograms//
	/* first   index VM process, 91=0, 93=1*/
	/* second  index VM species, rho=0, phi=1, jpsi=2*/
	/* third   index VM property, pt=0, eta=1, phi=2, theta=3, reserved=4*/

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
	TH1D* h_VM_sartre[2][3];
	TH1D* h_VM_daughter_sartre[2][3];
	//VM histograms//
	/* first   index coh/incoh, coh=0, incoh=1*/
	/* second  index pt,eta,phi, pt=0, eta=1, phi=2*/
	for(int icoh=0;icoh<2;icoh++){
		for(int ipro=0;ipro<3;ipro++){
				h_VM_sartre[icoh][ipro] = (TH1D*) file_sartre->Get(Form("h_VM_%d_%d",icoh,ipro));
				h_VM_daughter_sartre[icoh][ipro] = (TH1D*) file_sartre->Get(Form("h_VM_daughter_%d_%d",icoh,ipro));
			
		}
	}

	TCanvas* c1 = new TCanvas("c1","c1",1,1,600,600);
	gPad->SetLogy(1);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	TGaxis::SetMaxDigits(2);
	TH1D* base1 = makeHist("base1", "", "#eta_{lab}", "d#sigma/d#eta (nb) ", 100,-8,8,kBlack);
	if(name=="rho"){base1->GetYaxis()->SetRangeUser(0.1, 9e5);vm_index=0;}
	if(name=="phi"){base1->GetYaxis()->SetRangeUser(0.1, 4e5);vm_index=1;}
	if(name=="jpsi"){base1->GetYaxis()->SetRangeUser(0.1, 2e3);vm_index=2;}//jpsi nodecay.
	if(name=="rho_photo"){base1->GetYaxis()->SetRangeUser(0.1, 9e7);vm_index=0;}
	if(name=="phi_photo"){base1->GetYaxis()->SetRangeUser(0.1, 4e7);vm_index=1;}
	if(name=="jpsi_photo"){base1->GetYaxis()->SetRangeUser(0.1, 2e5);vm_index=2;}//jpsi nodecay.

	base1->GetXaxis()->SetTitleColor(kBlack);
	fixedFontHist1D(base1,1.,1.1);
	base1->GetYaxis()->SetTitleSize(base1->GetYaxis()->GetTitleSize()*1.5);
	base1->GetXaxis()->SetTitleSize(base1->GetXaxis()->GetTitleSize()*1.5);
	base1->GetYaxis()->SetLabelSize(base1->GetYaxis()->GetLabelSize()*1.5);
	base1->GetXaxis()->SetLabelSize(base1->GetXaxis()->GetLabelSize()*1.5);
	base1->GetXaxis()->SetNdivisions(5,5,0);
	base1->GetYaxis()->SetNdivisions(3,2,0);
	base1->Draw();

	measureXsection(name,h_VM[0][vm_index][1],0,t_hat_all->GetEntries(),PHP_);
	h_VM[0][vm_index][1]->SetFillColorAlpha(kBlue,0.4);
    h_VM[0][vm_index][1]->SetFillStyle(1001);
	h_VM[0][vm_index][1]->SetMarkerStyle(24);
	h_VM[0][vm_index][1]->SetMarkerColor(kBlue-1);

	measureXsection(name,h_VM[1][vm_index][1],0,t_hat_all->GetEntries(),PHP_);
	h_VM[1][vm_index][1]->SetFillColorAlpha(kRed,0.4);
    h_VM[1][vm_index][1]->SetFillStyle(1001);
	h_VM[1][vm_index][1]->SetMarkerStyle(25);
	h_VM[1][vm_index][1]->SetMarkerColor(kRed);

	h_VM[0][vm_index][1]->Draw("P e3 same");
	h_VM[1][vm_index][1]->Draw("P e3 same");

	/* Sartre */
	measureXsection(name,h_VM_sartre[1][1],1);
	h_VM_sartre[1][1]->SetLineColor(kBlack);
	h_VM_sartre[1][1]->SetLineWidth(2);
	h_VM_sartre[1][1]->Draw("hist same");
	
	TLegend *w5 = new TLegend(0.45,0.74,0.6,0.87);
	w5->SetLineColor(kWhite);
	w5->SetFillColor(0);
	w5->SetTextSize(18);
	w5->SetTextFont(45);
	// w5->AddEntry(h_phi_coh_sartre, "Sartre "+legendName+" coherent  ", "PL");
	w5->AddEntry(h_VM_sartre[1][1], "Sartre "+legendName+" incoherent  ", "L");
	w5->AddEntry(h_VM[0][vm_index][1], "BeAGLE elas. "+legendName+" incoh. ", "PF");
	w5->AddEntry(h_VM[1][vm_index][1], "BeAGLE diss. "+legendName+" incoh. ", "PF");
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

	TLatex* r44 = new TLatex(0.18, 0.84, "1<Q^{2}<20 GeV^{2}");
	r44->SetNDC();
	r44->SetTextSize(20);
	r44->SetTextFont(43);
	r44->SetTextColor(kBlack);

	TLatex* r44_1 = new TLatex(0.18, 0.84, "Q^{2}<0.2 GeV^{2}");
	r44_1->SetNDC();
	r44_1->SetTextSize(20);
	r44_1->SetTextFont(43);
	r44_1->SetTextColor(kBlack);
	if(PHP_) r44_1->Draw("same");
	else r44->Draw("same");

	TCanvas* c2 = new TCanvas("c2","c2",1,1,600,600);
	gPad->SetLogy(1);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	TGaxis::SetMaxDigits(2);
	TH1D* base2 = makeHist("base2", "", "p_{T} (GeV)", "d#sigma/dp_{T} (nb/GeV) ", 100,0,5,kBlack);
	if(name=="rho"){base2->GetYaxis()->SetRangeUser(0.1, 2.5e6);vm_index=0;}
	if(name=="phi"){base2->GetYaxis()->SetRangeUser(0.1, 2e5);vm_index=1;}
	if(name=="jpsi"){base2->GetYaxis()->SetRangeUser(0.1, 1e4);vm_index=2;}//jpsi nodecay.
	if(name=="rho_photo"){base2->GetYaxis()->SetRangeUser(0.1, 2.5e8);vm_index=0;}
	if(name=="phi_photo"){base2->GetYaxis()->SetRangeUser(0.1, 2e7);vm_index=1;}
	if(name=="jpsi_photo"){base2->GetYaxis()->SetRangeUser(0.1, 1e6);vm_index=2;}//jpsi nodecay.
	base2->GetXaxis()->SetTitleColor(kBlack);
	fixedFontHist1D(base2,1.,1.1);
	base2->GetYaxis()->SetTitleSize(base2->GetYaxis()->GetTitleSize()*1.5);
	base2->GetXaxis()->SetTitleSize(base2->GetXaxis()->GetTitleSize()*1.5);
	base2->GetYaxis()->SetLabelSize(base2->GetYaxis()->GetLabelSize()*1.5);
	base2->GetXaxis()->SetLabelSize(base2->GetXaxis()->GetLabelSize()*1.5);
	base2->GetXaxis()->SetNdivisions(5,5,0);
	base2->GetYaxis()->SetNdivisions(3,2,0);
	base2->Draw();

	measureXsection(name,h_VM[0][vm_index][0],0,t_hat_all->GetEntries(),PHP_);
	h_VM[0][vm_index][0]->SetFillColorAlpha(kBlue,0.4);
    h_VM[0][vm_index][0]->SetFillStyle(1001);
	h_VM[0][vm_index][0]->SetMarkerStyle(24);
	h_VM[0][vm_index][0]->SetMarkerColor(kBlue-1);

	measureXsection(name,h_VM[1][vm_index][0],0,t_hat_all->GetEntries(),PHP_);
	h_VM[1][vm_index][0]->SetFillColorAlpha(kRed,0.4);
    h_VM[1][vm_index][0]->SetFillStyle(1001);
	h_VM[1][vm_index][0]->SetMarkerStyle(25);
	h_VM[1][vm_index][0]->SetMarkerColor(kRed);

	h_VM[0][vm_index][0]->Draw("P e3 same");
	h_VM[1][vm_index][0]->Draw("P e3 same");

	/* Sartre */
	measureXsection(name,h_VM_sartre[1][0],1);
	h_VM_sartre[1][0]->SetLineColor(kBlack);
	h_VM_sartre[1][0]->SetLineWidth(2);
	h_VM_sartre[1][0]->Draw("hist same");

	w5->Draw("same");
	r42->Draw("same");
	r43->Draw("same");
	if(PHP_) r44_1->Draw("same");
	else r44->Draw("same");

	TCanvas* c3 = new TCanvas("c3","c3",1,1,600,600);
	gPad->SetLogy(1);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	TGaxis::SetMaxDigits(2);
	TH1D* base3 = makeHist("base3", "", "p_{T} (GeV)", "d#sigma/dp_{T} (nb/GeV) ", 100,0,5,kBlack);
	if(name=="rho"){base3->GetYaxis()->SetRangeUser(0.1, 2.5e6);vm_index=0;}
	if(name=="phi"){base3->GetYaxis()->SetRangeUser(0.1, 2e5);vm_index=1;}
	if(name=="jpsi"){base3->GetYaxis()->SetRangeUser(0.1, 1e4);vm_index=2;}//jpsi nodecay.
	if(name=="rho_photo"){base3->GetYaxis()->SetRangeUser(0.1, 2.5e8);vm_index=0;}
	if(name=="phi_photo"){base3->GetYaxis()->SetRangeUser(0.1, 2e7);vm_index=1;}
	if(name=="jpsi_photo"){base3->GetYaxis()->SetRangeUser(0.1, 1e6);vm_index=2;}//jpsi nodecay.
	base3->GetXaxis()->SetTitleColor(kBlack);
	fixedFontHist1D(base3,1.,1.1);
	base3->GetYaxis()->SetTitleSize(base3->GetYaxis()->GetTitleSize()*1.5);
	base3->GetXaxis()->SetTitleSize(base3->GetXaxis()->GetTitleSize()*1.5);
	base3->GetYaxis()->SetLabelSize(base3->GetYaxis()->GetLabelSize()*1.5);
	base3->GetXaxis()->SetLabelSize(base3->GetXaxis()->GetLabelSize()*1.5);
	base3->GetXaxis()->SetNdivisions(5,5,0);
	base3->GetYaxis()->SetNdivisions(3,2,0);
	base3->Draw();

	/* Sartre */
	measureXsection(name,h_VM_sartre[0][0],1);
	h_VM_sartre[0][0]->SetLineColor(kBlack);
	h_VM_sartre[0][0]->SetLineWidth(2);
	h_VM_sartre[0][0]->Draw("hist same");

	r42->Draw("same");
	r43->Draw("same");
	if(PHP_) r44_1->Draw("same");
	else r44->Draw("same");

	TLegend *w6 = new TLegend(0.6,0.84,0.75,0.87);
	w6->SetLineColor(kWhite);
	w6->SetFillColor(0);
	w6->SetTextSize(18);
	w6->SetTextFont(45);
	// w6->AddEntry(h_phi_coh_sartre, "Sartre "+legendName+" coherent  ", "PL");
	w6->AddEntry(h_VM_sartre[0][0], "Sartre "+legendName+" coherent  ", "L");
	w6->Draw("same");

	TCanvas* c4 = new TCanvas("c4","c4",1,1,600,600);
	gPad->SetLogy(1);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	TGaxis::SetMaxDigits(2);
	TH1D* base4 = makeHist("base4", "", "#eta_{lab}", "d#sigma/d#eta (nb) ", 100,-8,8,kBlack);
	if(name=="rho"){base4->GetYaxis()->SetRangeUser(0.1, 9e5);vm_index=0;}
	if(name=="phi"){base4->GetYaxis()->SetRangeUser(0.1, 4e5);vm_index=1;}
	if(name=="jpsi"){base4->GetYaxis()->SetRangeUser(0.1, 2e3);vm_index=2;}//jpsi nodecay.
	if(name=="rho_photo"){base4->GetYaxis()->SetRangeUser(0.1, 9e7);vm_index=0;}
	if(name=="phi_photo"){base4->GetYaxis()->SetRangeUser(0.1, 4e7);vm_index=1;}
	if(name=="jpsi_photo"){base4->GetYaxis()->SetRangeUser(0.1, 2e5);vm_index=2;}//jpsi nodecay.
	
	base4->GetXaxis()->SetTitleColor(kBlack);
	fixedFontHist1D(base4,1.,1.1);
	base4->GetYaxis()->SetTitleSize(base4->GetYaxis()->GetTitleSize()*1.5);
	base4->GetXaxis()->SetTitleSize(base4->GetXaxis()->GetTitleSize()*1.5);
	base4->GetYaxis()->SetLabelSize(base4->GetYaxis()->GetLabelSize()*1.5);
	base4->GetXaxis()->SetLabelSize(base4->GetXaxis()->GetLabelSize()*1.5);
	base4->GetXaxis()->SetNdivisions(5,5,0);
	base4->GetYaxis()->SetNdivisions(3,2,0);
	base4->Draw();

	/* Sartre */
	measureXsection(name,h_VM_sartre[0][1],1);
	h_VM_sartre[0][1]->SetLineColor(kBlack);
	h_VM_sartre[0][1]->SetLineWidth(2);
	h_VM_sartre[0][1]->Draw("hist same");

	r42->Draw("same");
	r43->Draw("same");
	w6->Draw("same");
	if(PHP_) r44_1->Draw("same");
	else r44->Draw("same");

	if(veto_){
		c1->Print("../figures/kinematicsDistributions/veto_"+name+"_eta_incoh.pdf");
		c2->Print("../figures/kinematicsDistributions/veto_"+name+"_pt_incoh.pdf");
		c3->Print("../figures/kinematicsDistributions/veto_"+name+"_pt_coh.pdf");
		c4->Print("../figures/kinematicsDistributions/veto_"+name+"_eta_coh.pdf");
	}
	else{
		c1->Print("../figures/kinematicsDistributions/"+name+"_eta_incoh.pdf");
		c2->Print("../figures/kinematicsDistributions/"+name+"_pt_incoh.pdf");
		c3->Print("../figures/kinematicsDistributions/"+name+"_pt_coh.pdf");
		c4->Print("../figures/kinematicsDistributions/"+name+"_eta_coh.pdf");

	}
	

}