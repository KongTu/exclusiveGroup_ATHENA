#include "utility.h"
void plotVM_mass_and_t(TString name="phi", bool veto_=false, bool PHP_ = false, int method=0, int coh = 1){

	setVM(name);
	/* Beagle */
	if(name=="rho_photo"
	 	|| name=="phi_photo"
	 	  || name=="jpsi_photo") PHP_ = true;

	int sartre_vm_index=0;
	int beagle_vm_index=0;

	if(name=="phi"||name=="phi_photo") {
		beagle_vm_index=0;
		sartre_vm_index=1;
	}
	if(name=="jpsi"||name=="jpsi_photo") {
		beagle_vm_index=0;
		sartre_vm_index=2;
	}

	int processindex=1;//91,93

	TString inputROOT="../rootfiles/beagle_allVMs_w_breakups.root";
	if(PHP_) inputROOT="../rootfiles/beagle_allVMs_w_breakups_PHP.root";
	if(veto_&&!PHP_) inputROOT="../rootfiles/beagle_allVMs_w_breakups_w_vetos.root";
	if(veto_&&PHP_) inputROOT="../rootfiles/beagle_allVMs_w_breakups_w_vetos_PHP.root";
	TFile* file_beagle = new TFile(inputROOT);
	TH1D* t_hat_all = (TH1D*) file_beagle->Get("h_trueT");

	//beagle
	TH1D* h_t_reco[3][3][3][3];
	TH2D* h_VM_t_mass[3][3][3][3];
	for(int ibreak=0;ibreak<3;ibreak++){
		for(int ivm=0;ivm<3;ivm++){
			for(int imethod=0;imethod<3;imethod++){
				for(int imass=0;imass<3;imass++){
					h_t_reco[ibreak][ivm][imethod][imass] = (TH1D*) file_beagle->Get(Form("h_t_reco_%d_%d_%d_%d",ibreak,ivm,imethod,imass));
					h_VM_t_mass[ibreak][ivm][imethod][imass] = (TH2D*) file_beagle->Get(Form("h_VM_t_mass_%d_%d_%d_%d",ibreak,ivm,imethod,imass));
				}
			}
		}
	}
	TH1D* h_VM_mass[3][3][3];
	for(int ibreak=0;ibreak<3;ibreak++){
		for(int ivm=0;ivm<3;ivm++){
			for(int imass=0;imass<3;imass++){
				h_VM_mass[ibreak][ivm][imass] = (TH1D*) file_beagle->Get(Form("h_VM_mass_%d_%d_%d",ibreak,ivm,imass));
				h_VM_mass[ibreak][ivm][imass]->SetLineColor(kRed);
			}
		}
	}
	/* Sartre */
	TFile* file_sartre_all[3];
	if(PHP_){
		file_sartre_all[0] = new TFile("../rootfiles/sartre_rho_photo_bnonsat.root");
		file_sartre_all[1] = new TFile("../rootfiles/sartre_phi_photo_bnonsat.root");
		file_sartre_all[2] = new TFile("../rootfiles/sartre_jpsi_photo_bnonsat.root");
	}
	else{
		file_sartre_all[0] = new TFile("../rootfiles/sartre_rho_bnonsat.root");
		file_sartre_all[1] = new TFile("../rootfiles/sartre_phi_bnonsat.root");
		file_sartre_all[2] = new TFile("../rootfiles/sartre_jpsi_bnonsat.root");
	}
	//sartre wrong mass setting in different rootfiles unfortunately.
	TFile* file_sartre_all_wrongmass[3];
	if(PHP_){
		file_sartre_all_wrongmass[0] = new TFile("../rootfiles/sartre_rho_photo_bnonsat_extra.root");
		file_sartre_all_wrongmass[1] = new TFile("../rootfiles/sartre_phi_photo_bnonsat_extra.root");
		file_sartre_all_wrongmass[2] = new TFile("../rootfiles/sartre_jpsi_photo_bnonsat_extra.root");
	}
	else{
		file_sartre_all_wrongmass[0] = new TFile("../rootfiles/sartre_rho_bnonsat_extra.root");
		file_sartre_all_wrongmass[1] = new TFile("../rootfiles/sartre_phi_bnonsat_extra.root");
		file_sartre_all_wrongmass[2] = new TFile("../rootfiles/sartre_jpsi_bnonsat_extra.root");
	}
	

	TFile* file_sartre = new TFile("../rootfiles/sartre_"+name+"_bnonsat.root");
	//not use here
	TH1D* h_coh_sartre = (TH1D*) file_sartre->Get("hist_t_coherent");
	TH1D* h_incoh_sartre = (TH1D*) file_sartre->Get("hist_t_incoherent");
	TH1D* hist_t_afterPhaseSpace_coherent_sartre = (TH1D*) file_sartre->Get("hist_t_afterPhaseSpace_coherent");
	TH1D* hist_t_afterPhaseSpace_incoherent_sartre = (TH1D*) file_sartre->Get("hist_t_afterPhaseSpace_incoherent");
	//end not use here
	// sartre
    TH1D* h_t_reco_sartre[2][3][3];
    TH2D* h_VM_t_mass_sartre[2][3][3];
    TH2D* h_VM_t_mass_sartre_mixed[2][3][3][3];
    TH2D* h_VM_t_mass_sartre_wrong[2][3][3][3];
    for(int ibreak=0;ibreak<2;ibreak++){
        for(int imethod=0;imethod<3;imethod++){
        	for(int imass=0;imass<3;imass++){
				h_t_reco_sartre[ibreak][imethod][imass] = (TH1D*) file_sartre->Get(Form("h_t_reco_%d_%d_%d",ibreak,imethod,imass));
				h_VM_t_mass_sartre[ibreak][imethod][imass] = (TH2D*) file_sartre->Get(Form("h_VM_t_mass_%d_%d_%d",ibreak,imethod,imass));
        		for(int ifile=0;ifile<3;ifile++){
	    			h_VM_t_mass_sartre_mixed[ibreak][imethod][imass][ifile] = (TH2D*) file_sartre_all[ifile]->Get(Form("h_VM_t_mass_%d_%d_%d",ibreak,imethod,imass));
	    			h_VM_t_mass_sartre_wrong[ibreak][imethod][imass][ifile] = (TH2D*) file_sartre_all_wrongmass[ifile]->Get(Form("h_VM_t_mass_%d_%d_%d",ibreak,imethod,imass));
	    		}
        	}
        }
    }
    TH1D* h_VM_mass_sartre[2][3];
    TH1D* h_VM_mass_sartre_mixed[2][3][3];
    for(int ibreak=0;ibreak<2;ibreak++){
    	for(int imass=0;imass<3;imass++){
			h_VM_mass_sartre[ibreak][imass] = (TH1D*) file_sartre->Get(Form("h_VM_mass_%d_%d",ibreak,imass));
    		h_VM_mass_sartre[ibreak][imass]->SetLineColor(kBlue);
    		for(int ifile=0;ifile<3;ifile++){
    			h_VM_mass_sartre_mixed[ibreak][imass][ifile] = (TH1D*) file_sartre_all[ifile]->Get(Form("h_VM_mass_%d_%d",ibreak,imass));
	    		
	    		// if(PHP_) h_VM_mass_sartre_mixed[ibreak][imass][ifile]->Scale(sigma_sartre_photo[ifile]);
	    		// else h_VM_mass_sartre_mixed[ibreak][imass][ifile]->Scale(sigma_sartre_elect[ifile]);
	    		
	    		h_VM_mass_sartre_mixed[ibreak][imass][ifile]->SetLineColor(kBlue);
	    		h_VM_mass_sartre_mixed[ibreak][imass][ifile]->SetMarkerStyle(28);
    		
    		}
    	}
    }
    TCanvas* c1_1 = new TCanvas("c1_1","c1_1",1,1,600,600);
	gPad->SetLogy(1);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	TH1D* base1 = makeHist("base1", "", "Mass (GeV/c^{2})", "Normalized counts ", 100,0,3.5,kBlack);
	base1->GetYaxis()->SetRangeUser(1e-7, 3e0);
	base1->GetXaxis()->SetTitleColor(kBlack);
	TGaxis::SetMaxDigits(2);
	fixedFontHist1D(base1,1.,1.3);
	base1->GetYaxis()->SetTitleSize(base1->GetYaxis()->GetTitleSize()*1.5);
	base1->GetXaxis()->SetTitleSize(base1->GetXaxis()->GetTitleSize()*1.5);
	base1->GetYaxis()->SetLabelSize(base1->GetYaxis()->GetLabelSize()*1.5);
	base1->GetXaxis()->SetLabelSize(base1->GetXaxis()->GetLabelSize()*1.5);
	base1->GetXaxis()->SetNdivisions(4,4,0);
	base1->GetYaxis()->SetNdivisions(5,5,0);
	base1->Draw();

    h_VM_mass_sartre[coh][sartre_vm_index]->GetXaxis()->SetTitle("Mass (GeV/c^{2})");
    h_VM_mass_sartre[coh][sartre_vm_index]->GetYaxis()->SetTitle("Normalized Counts");
    h_VM_mass_sartre[coh][sartre_vm_index]->SetTitle("first t bin");

    h_VM_mass_sartre[coh][sartre_vm_index]->DrawNormalized("same");
	h_VM_mass[processindex][vm_index][beagle_vm_index]->DrawNormalized("same");

	TLegend *w5 = new TLegend(0.53,0.74,0.73,0.87);
	w5->SetLineColor(kWhite);
	w5->SetFillColor(0);
	w5->SetTextSize(18);
	w5->SetTextFont(45);
	w5->AddEntry(h_VM_mass_sartre[coh][sartre_vm_index], "Sartre "+legendName+" incoherent  ", "L");
	w5->AddEntry(h_VM_mass[processindex][vm_index][beagle_vm_index], "BeAGLE "+legendName+" incoherent  ", "L");
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
	c1_1->Print("../figures/mass/"+name+"_mass.pdf");

	if(name!="phi"&&name!="phi_photo") return;

	TCanvas* c2 = new TCanvas("c2","c2",1,1,600,600);
	gPad->SetLogy(1);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	TH1D* base2 = (TH1D*) base1->Clone("base2");
	base2->GetYaxis()->SetTitle("~ cross section");
	base2->GetYaxis()->SetRangeUser(1e0,1e7);
	base2->Draw();
	TString name_temp = "rho";
	if(PHP_) name_temp = "rho_photo";
	measureXsection(name_temp,h_VM_mass_sartre_mixed[coh][sartre_vm_index][0], 1);
	h_VM_mass_sartre_mixed[coh][sartre_vm_index][0]->Draw("hist same");
	name_temp = "phi";
	if(PHP_) name_temp = "phi_photo";
	measureXsection(name_temp,h_VM_mass_sartre_mixed[coh][sartre_vm_index][1], 1);
	h_VM_mass_sartre_mixed[coh][sartre_vm_index][1]->Draw("hist same");
	name_temp = "jpsi";
	if(PHP_) name_temp = "jpsi_photo";
	measureXsection(name_temp,h_VM_mass_sartre_mixed[coh][sartre_vm_index][2], 1);
	h_VM_mass_sartre_mixed[coh][sartre_vm_index][2]->Draw("hist same");
	setVM(name);

	if(PHP_) r44_1->Draw("same");
	else r44->Draw("same");
	TLegend *w4 = new TLegend(0.53,0.8,0.73,0.87);
	w4->SetLineColor(kWhite);
	w4->SetFillColor(0);
	w4->SetTextSize(18);
	w4->SetTextFont(45);
	w4->AddEntry(h_VM_mass_sartre_mixed[coh][sartre_vm_index][0], "Sartre "+legendName+" w. PID  ", "L");
	// w4->AddEntry(h_VM_mass[processindex][beagle_vm_index][1], "BeAGLE "+legendName+" w. PID  ", "L");
	w4->Draw("same");
	r42->Draw("same");
	r43->Draw("same");
	c2->Print("../figures/mass/sartre_mass_mixing_"+name+".pdf");

	TCanvas* c3 = new TCanvas("c3","c3",1,1,600,600);
	gPad->SetLogy(1);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	TH1D* base3 = (TH1D*) base1->Clone("base3");
	base3->GetYaxis()->SetRangeUser(1,1e6);
	base3->GetYaxis()->SetTitle("~ cross section");
	base3->Draw();
	beagle_vm_index=2; //all mixed.
	h_VM_mass[processindex][vm_index][1]->SetLineColor(kBlack);
	h_VM_mass[processindex][vm_index][1]->Draw("same");
	h_VM_mass[processindex][vm_index][2]->Draw("same");


	if(PHP_) r44_1->Draw("same");
	else r44->Draw("same");
	TLegend *w6 = new TLegend(0.53,0.74,0.73,0.87);
	w6->SetLineColor(kWhite);
	w6->SetFillColor(0);
	w6->SetTextSize(18);
	w6->SetTextFont(45);
	w6->AddEntry(h_VM_mass[processindex][vm_index][1], "BeAGLE "+legendName+" No PID  ", "L");
	w6->AddEntry(h_VM_mass[processindex][beagle_vm_index][1], "BeAGLE "+legendName+" w. PID  ", "L");
	w6->Draw("same");
	r42->Draw("same");
	r43->Draw("same");
	c3->Print("../figures/mass/beagle_mass_mixing_"+name+".pdf");

	
	//only plotting phi
	TCanvas* c11 = new TCanvas("c11","c11",1,1,600,600);
	gPad->SetLogy(1);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	TH1D* base11 = makeHist("base11", "", "|#it{t} | (GeV^{2})", "d#sigma/d|#it{t} | (nb/GeV^{2}) ", 100,0,0.18,kBlack);
	
	if(name=="phi"){base11->GetYaxis()->SetRangeUser(1e-2, 1e7);}
	if(name=="phi_photo"){base11->GetYaxis()->SetRangeUser(1e-1, 1e8);}

	base11->GetXaxis()->SetTitleColor(kBlack);
	TGaxis::SetMaxDigits(2);
	fixedFontHist1D(base11,1.,1.1);
	base11->GetYaxis()->SetTitleSize(base11->GetYaxis()->GetTitleSize()*1.5);
	base11->GetXaxis()->SetTitleSize(base11->GetXaxis()->GetTitleSize()*1.5);
	base11->GetYaxis()->SetLabelSize(base11->GetYaxis()->GetLabelSize()*1.5);
	base11->GetXaxis()->SetLabelSize(base11->GetXaxis()->GetLabelSize()*1.5);
	base11->GetXaxis()->SetNdivisions(4,4,0);
	base11->GetYaxis()->SetNdivisions(5,5,0);
	base11->Draw();

	//coherent
	TH1D* h_mass_for_binning = (TH1D*) h_VM_t_mass_sartre[0][method][1]->ProjectionX("h_mass_for_binning",1,1000);
	int bin_lower = h_mass_for_binning->FindBin(1.019-0.02);
	int bin_upper = h_mass_for_binning->FindBin(1.019+0.02);

	//true phi sample, with phi ->kk mass.
	TH1D* h_t_from_mass_coherent = (TH1D*) h_VM_t_mass_sartre[0][method][sartre_vm_index]->ProjectionY("h_t_from_mass_coherent",bin_lower,bin_upper);
	h_t_from_mass_coherent->SetLineColor(kBlue);
	measureXsection(name, h_t_from_mass_coherent, 1);
	h_t_from_mass_coherent->Draw("hist same");
	
	//use mixed
	TH1D* h_t_from_mass_coherent_mixed_rho = (TH1D*) h_VM_t_mass_sartre_mixed[0][method][sartre_vm_index][0]->ProjectionY("h_t_from_mass_coherent_mixed_rho",bin_lower,bin_upper);
	name_temp = "rho";
	if(PHP_) name_temp = "rho_photo";
	measureXsection(name_temp, h_t_from_mass_coherent_mixed_rho, 1);
	
	TH1D* h_t_from_mass_coherent_mixed_phi = (TH1D*) h_VM_t_mass_sartre_mixed[0][method][sartre_vm_index][1]->ProjectionY("h_t_from_mass_coherent_mixed_phi",bin_lower,bin_upper);
	name_temp = "phi";
	if(PHP_) name_temp = "phi_photo";
	measureXsection(name_temp, h_t_from_mass_coherent_mixed_phi, 1);

	h_t_from_mass_coherent_mixed_phi->SetLineColor(kBlack);
	h_t_from_mass_coherent_mixed_phi->SetMarkerStyle(25);
	h_t_from_mass_coherent_mixed_phi->Add(h_t_from_mass_coherent_mixed_rho,+1);
	h_t_from_mass_coherent_mixed_phi->Draw("P same");

	//using wrong.
	TH1D* h_t_from_mass_coherent_wrong_rho = (TH1D*) h_VM_t_mass_sartre_wrong[0][method][sartre_vm_index][0]->ProjectionY("h_t_from_mass_coherent_wrong_rho",bin_lower,bin_upper);
	name_temp = "rho";
	if(PHP_) name_temp = "rho_photo";
	measureXsection(name_temp, h_t_from_mass_coherent_wrong_rho, 1);

	TH1D* h_t_from_mass_coherent_wrong_phi = (TH1D*) h_VM_t_mass_sartre_wrong[0][method][sartre_vm_index][1]->ProjectionY("h_t_from_mass_coherent_wrong_phi",bin_lower,bin_upper);
	name_temp = "phi";
	if(PHP_) name_temp = "phi_photo";
	measureXsection(name_temp, h_t_from_mass_coherent_wrong_phi, 1);

	h_t_from_mass_coherent_wrong_phi->SetMarkerColor(kRed);
	h_t_from_mass_coherent_wrong_phi->SetLineColor(kRed);
	h_t_from_mass_coherent_wrong_phi->SetMarkerStyle(25);
	h_t_from_mass_coherent_wrong_phi->Add(h_t_from_mass_coherent_wrong_rho,+1);
	h_t_from_mass_coherent_wrong_phi->Draw("p same");

	setVM(name);
	
	//incoherent
	TH1D* h_t_from_truemass_incoherent_91 = (TH1D*) h_VM_t_mass[0][vm_index][method][0]->ProjectionY("h_t_from_truemass_incoherent_91",bin_lower,bin_upper);
	TH1D* h_t_from_truemass_incoherent_93 = (TH1D*) h_VM_t_mass[1][vm_index][method][0]->ProjectionY("h_t_from_truemass_incoherent_93",bin_lower,bin_upper);
	h_t_from_truemass_incoherent_91->Add(h_t_from_truemass_incoherent_93,+1);
	measureXsection(name, h_t_from_truemass_incoherent_91, 0, t_hat_all->GetEntries(), PHP_);
	h_t_from_truemass_incoherent_91->Rebin(2);
	h_t_from_truemass_incoherent_91->Scale(1./2);
	h_t_from_truemass_incoherent_91->SetMarkerStyle(24);
	h_t_from_truemass_incoherent_91->SetMarkerStyle(20);
	h_t_from_truemass_incoherent_91->SetLineStyle(2);
	h_t_from_truemass_incoherent_91->SetLineWidth(2);
	h_t_from_truemass_incoherent_91->SetLineColor(kRed);
	h_t_from_truemass_incoherent_91->Draw("hist same");

	//incoherent wrong mass
	TH1D* h_t_from_wrongmass_incoherent_91 = (TH1D*) h_VM_t_mass[0][vm_index][method][1]->ProjectionY("h_t_from_wrongmass_incoherent_91",bin_lower,bin_upper);
	TH1D* h_t_from_wrongmass_incoherent_93 = (TH1D*) h_VM_t_mass[1][vm_index][method][1]->ProjectionY("h_t_from_wrongmass_incoherent_93",bin_lower,bin_upper);
	h_t_from_wrongmass_incoherent_91->Add(h_t_from_wrongmass_incoherent_93,+1);
	measureXsection(name, h_t_from_wrongmass_incoherent_91, 0, t_hat_all->GetEntries(), PHP_);
	h_t_from_wrongmass_incoherent_91->Rebin(2);
	h_t_from_wrongmass_incoherent_91->Scale(1./2);
	h_t_from_wrongmass_incoherent_91->SetMarkerStyle(24);
	h_t_from_wrongmass_incoherent_91->SetMarkerColor(kRed);
	h_t_from_wrongmass_incoherent_91->Draw("P same");

	//incoherent pid mass
	TH1D* h_t_from_PIDmass_incoherent_91 = (TH1D*) h_VM_t_mass[0][vm_index][method][2]->ProjectionY("h_t_from_PIDmass_incoherent_91",bin_lower,bin_upper);
	TH1D* h_t_from_PIDmass_incoherent_93 = (TH1D*) h_VM_t_mass[1][vm_index][method][2]->ProjectionY("h_t_from_PIDmass_incoherent_93",bin_lower,bin_upper);
	h_t_from_PIDmass_incoherent_91->Add(h_t_from_PIDmass_incoherent_93,+1);
	measureXsection(name, h_t_from_PIDmass_incoherent_91, 0, t_hat_all->GetEntries(), PHP_);
	h_t_from_PIDmass_incoherent_91->Rebin(2);
	h_t_from_PIDmass_incoherent_91->Scale(1./2);
	h_t_from_PIDmass_incoherent_91->SetMarkerStyle(24);
	h_t_from_PIDmass_incoherent_91->SetMarkerColor(kBlue);
	h_t_from_PIDmass_incoherent_91->Draw("P same");

	TLegend *w7 = new TLegend(0.43,0.65,0.66,0.8);
	w7->SetLineColor(kWhite);
	w7->SetFillColor(0);
	w7->SetTextSize(18);
	w7->SetTextFont(45);
	w7->AddEntry(h_t_from_mass_coherent, "Sartre "+legendName+" coh. truth  ", "L");
	w7->AddEntry(h_t_from_mass_coherent_wrong_phi, "Sartre "+legendName+" coh. mixed NO PID  ", "PL");
	w7->AddEntry(h_t_from_mass_coherent_mixed_phi, "Sartre "+legendName+" coh. mixed but w. PID  ", "PL");
	w7->AddEntry(h_t_from_truemass_incoherent_91, "BeAGLE "+legendName+" incoh. truth  ", "L");
	w7->AddEntry(h_t_from_wrongmass_incoherent_91, "BeAGLE "+legendName+" incoh mixed NO PID  ", "PL");
	w7->AddEntry(h_t_from_PIDmass_incoherent_91, "BeAGLE "+legendName+" incoh mixed w. PID  ", "PL");
	w7->Draw("same");

	r42->Draw("same");
	r43->Draw("same");
	if(PHP_) r44_1->Draw("same");
	else r44->Draw("same");

	if(veto_) c11->Print("../figures/dsigmadt_mass/veto_dsigma_dt_"+name+".pdf");
	else c11->Print("../figures/dsigmadt_mass/dsigma_dt_"+name+".pdf");

	TCanvas* c22 = new TCanvas("c22","c22",1,1,600,600);
	gPad->SetLogy(0);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	TH1D* base22 = makeHist("base22", "", "|#it{t} | (GeV^{2})", "Ratio to truth ", 100,0,0.18,kBlack);
	if(name=="phi"){base22->GetYaxis()->SetRangeUser(0, 2);}
	if(name=="phi_photo"){base22->GetYaxis()->SetRangeUser(0, 2);}
	base22->GetXaxis()->SetTitleColor(kBlack);
	TGaxis::SetMaxDigits(2);
	fixedFontHist1D(base22,1.,1.1);
	base22->GetYaxis()->SetTitleSize(base22->GetYaxis()->GetTitleSize()*1.5);
	base22->GetXaxis()->SetTitleSize(base22->GetXaxis()->GetTitleSize()*1.5);
	base22->GetYaxis()->SetLabelSize(base22->GetYaxis()->GetLabelSize()*1.5);
	base22->GetXaxis()->SetLabelSize(base22->GetXaxis()->GetLabelSize()*1.5);
	base22->GetXaxis()->SetNdivisions(4,4,0);
	base22->GetYaxis()->SetNdivisions(5,5,0);
	base22->Draw();

	TH1D* ratio_1 = (TH1D*) h_t_from_mass_coherent_wrong_phi->Clone("ratio_1");
	TH1D* ratio_2 = (TH1D*) h_t_from_mass_coherent_mixed_phi->Clone("ratio_2");
	ratio_1->Divide(h_t_from_mass_coherent);
	ratio_2->Divide(h_t_from_mass_coherent);

	ratio_1->Draw("Psame");
	ratio_2->Draw("Psame");

	TLegend *w8 = new TLegend(0.43,0.32,0.66,0.4);
	w8->SetLineColor(kWhite);
	w8->SetFillColor(0);
	w8->SetTextSize(18);
	w8->SetTextFont(45);
	w8->AddEntry(ratio_1, "Sartre "+legendName+" coh. mixed NO PID  ", "PL");
	w8->AddEntry(ratio_2, "Sartre "+legendName+" coh. mixed w. PID  ", "PL");

	w8->Draw("same");
	if(veto_) c22->Print("../figures/dsigmadt_mass/ratio_veto_dsigma_dt_"+name+".pdf");
	else c22->Print("../figures/dsigmadt_mass/ratio_dsigma_dt_"+name+".pdf");


}