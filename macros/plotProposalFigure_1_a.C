#include "utility.h"
void plotProposalFigure_1_a(TString name="phi",int PHP_ = 0, int veto_ = 0, double minPt_=0.15, int method=2){

	TFile* beam_piple_file = new TFile("../analysis/include/beam_pipe.root");
	TH1D* h_beam_pipe = (TH1D*) beam_piple_file->Get("beam_pipe_factor");

	int sartre_vm_index=1;
	int beagle_vm_index=0;
	int processindex=1;//91,93
	double scale_factor = 0.822;
	if(name=="phi_photo") {scale_factor=0.335;PHP_=1;}

	TString inputROOT=Form("../rootfiles/beagle_output_PHP_%d_veto_%d_minPt_%.2f_forPID_new.root",PHP_,veto_,minPt_);
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
		file_sartre_all[0] = new TFile(Form("../rootfiles/sartre_rho_photo_bnonsat_PID_1_minPt_%.2f_forPID.root",minPt_));
		file_sartre_all[1] = new TFile(Form("../rootfiles/sartre_phi_photo_bnonsat_PID_1_minPt_%.2f_forPID.root",minPt_));
		file_sartre_all[2] = new TFile(Form("../rootfiles/sartre_jpsi_photo_bnonsat_PID_1_minPt_%.2f_forPID.root",minPt_));
	}
	else{
		file_sartre_all[0] = new TFile(Form("../rootfiles/sartre_rho_bnonsat_PID_1_minPt_%.2f_forPID.root",minPt_));//changed to _smear_0 to see full |Eta|<4.0 PID effect.
		file_sartre_all[1] = new TFile(Form("../rootfiles/sartre_phi_bnonsat_PID_1_minPt_%.2f_forPID.root",minPt_));
		file_sartre_all[2] = new TFile(Form("../rootfiles/sartre_jpsi_bnonsat_PID_1_minPt_%.2f_forPID.root",minPt_));
	}	
	//sartre wrong mass setting in different rootfiles unfortunately.
	TFile* file_sartre_all_wrongmass[3];
	if(PHP_){
		file_sartre_all_wrongmass[0] = new TFile(Form("../rootfiles/sartre_rho_photo_bnonsat_PID_0_minPt_%.2f_forPID.root",minPt_));
		file_sartre_all_wrongmass[1] = new TFile(Form("../rootfiles/sartre_phi_photo_bnonsat_PID_0_minPt_%.2f_forPID.root",minPt_));
		file_sartre_all_wrongmass[2] = new TFile(Form("../rootfiles/sartre_jpsi_photo_bnonsat_PID_0_minPt_%.2f_forPID.root",minPt_));
	}
	else{
		file_sartre_all_wrongmass[0] = new TFile(Form("../rootfiles/sartre_rho_bnonsat_PID_0_minPt_%.2f_forPID.root",minPt_));
		file_sartre_all_wrongmass[1] = new TFile(Form("../rootfiles/sartre_phi_bnonsat_PID_0_minPt_%.2f_forPID.root",minPt_));
		file_sartre_all_wrongmass[2] = new TFile(Form("../rootfiles/sartre_jpsi_bnonsat_PID_0_minPt_%.2f_forPID.root",minPt_));
	}

	TFile* file_sartre = new TFile(Form("../rootfiles/sartre_"+name+"_bnonsat_PID_1_minPt_%.2f_forPID.root",minPt_));
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
	    		
	    		h_VM_mass_sartre_mixed[ibreak][imass][ifile]->SetLineColor(kBlue);
	    		h_VM_mass_sartre_mixed[ibreak][imass][ifile]->SetMarkerStyle(28);
    		
    		}
    	}
    }
 
   //only plotting phi
	TCanvas* c11 = new TCanvas("c11","c11",1,1,600,600);
	gPad->SetLogy(1);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	TH1D* base11 = makeHist("base11", "", "|#it{t} | (GeV^{2})", "d#sigma/d|#it{t} | (nb/GeV^{2}) ", 100,0,0.18,kBlack);
	base11->GetYaxis()->SetRangeUser(1e-1, 1e8);
	if(name=="phi"){base11->GetYaxis()->SetRangeUser(1e-2, 1e7);}
	if(name=="phi_photo"){base11->GetYaxis()->SetRangeUser(1e-1, 1e9);}
	base11->GetXaxis()->SetTitleColor(kBlack);
	TGaxis::SetMaxDigits(3);
	fixedFontHist1D(base11,1.,1.1);
	base11->GetYaxis()->SetTitleSize(base11->GetYaxis()->GetTitleSize()*1.5);
	base11->GetXaxis()->SetTitleSize(base11->GetXaxis()->GetTitleSize()*1.5);
	base11->GetYaxis()->SetLabelSize(base11->GetYaxis()->GetLabelSize()*1.5);
	base11->GetXaxis()->SetLabelSize(base11->GetXaxis()->GetLabelSize()*1.5);
	base11->GetXaxis()->SetNdivisions(4,4,0);
	base11->GetYaxis()->SetNdivisions(5,5,0);
	base11->Draw();

	measureXsection(name, h_coh_sartre, 1);
	h_coh_sartre->SetLineColor(kBlue);
	h_coh_sartre->SetLineWidth(2);
	h_coh_sartre->Draw("hist same");

	measureXsection(name, hist_t_afterPhaseSpace_coherent_sartre, 1);
	hist_t_afterPhaseSpace_coherent_sartre->SetLineColor(kRed);
	// hist_t_afterPhaseSpace_coherent_sartre->Draw("hist same");

	//coherent
	TH1D* h_mass_for_binning = (TH1D*) h_VM_t_mass_sartre[0][method][1]->ProjectionX("h_mass_for_binning",1,1000);
	int bin_lower = h_mass_for_binning->FindBin(1.019-0.04);
	int bin_upper = h_mass_for_binning->FindBin(1.019+0.04);

	//use PID
	TH1D* h_t_from_mass_coherent_pid_rho = (TH1D*) h_VM_t_mass_sartre_mixed[0][method][sartre_vm_index][0]->ProjectionY("h_t_from_mass_coherent_pid_rho",bin_lower,bin_upper);
	TString name_temp = "rho";
	if(PHP_) name_temp = "rho_photo";
	measureXsection(name_temp, h_t_from_mass_coherent_pid_rho, 1);
	h_t_from_mass_coherent_pid_rho->SetLineWidth(2);
	h_t_from_mass_coherent_pid_rho->SetLineStyle(2);
	h_t_from_mass_coherent_pid_rho->SetLineColor(kBlack);
	// h_t_from_mass_coherent_pid_rho->Draw("hist same");

	// //using wrong.
	TH1D* h_t_from_mass_coherent_wrong_rho = (TH1D*) h_VM_t_mass_sartre_wrong[0][method][sartre_vm_index][0]->ProjectionY("h_t_from_mass_coherent_wrong_rho",bin_lower,bin_upper);
	name_temp = "rho";
	if(PHP_) name_temp = "rho_photo";
	measureXsection(name_temp, h_t_from_mass_coherent_wrong_rho, 1);
	h_t_from_mass_coherent_wrong_rho->SetLineWidth(2);
	h_t_from_mass_coherent_wrong_rho->SetLineStyle(3);
	h_t_from_mass_coherent_wrong_rho->SetLineColor(kRed);
	// h_t_from_mass_coherent_wrong_rho->Draw("hist same");

	/*
	Incoherent contributions, plot separately.
	*/

	//incoherent true mass - phi incoherent.
	TH1D* h_t_from_truemass_incoherent_91 = (TH1D*) h_VM_t_mass[0][1][method][0]->ProjectionY("h_t_from_truemass_incoherent_91",bin_lower,bin_upper);
	TH1D* h_t_from_truemass_incoherent_93 = (TH1D*) h_VM_t_mass[1][1][method][0]->ProjectionY("h_t_from_truemass_incoherent_93",bin_lower,bin_upper);
	h_t_from_truemass_incoherent_91->Add(h_t_from_truemass_incoherent_93,+1);
	measureXsection(name, h_t_from_truemass_incoherent_91, 0, t_hat_all->GetEntries(), PHP_);
	h_t_from_truemass_incoherent_91->SetLineStyle(2);
	h_t_from_truemass_incoherent_91->SetLineWidth(3);
	h_t_from_truemass_incoherent_91->SetLineColor(kBlack);
	
	h_t_from_truemass_incoherent_91->Scale(1./scale_factor);//number coming from integral ratio for t>0.0;
	//add beam pipe effect.
	for(int ibin=0;ibin<h_t_from_truemass_incoherent_91->GetNbinsX();ibin++){
		double bincenter = h_t_from_truemass_incoherent_91->GetBinCenter(ibin+1);
		double weight = h_beam_pipe->GetBinContent( h_beam_pipe->FindBin(bincenter) );
		if(veto_){
			h_t_from_truemass_incoherent_91->SetBinContent(ibin+1, h_t_from_truemass_incoherent_91->GetBinContent(ibin+1)*weight);
			h_t_from_truemass_incoherent_91->SetBinError(ibin+1, h_t_from_truemass_incoherent_91->GetBinError(ibin+1)*weight);
		}
	}
	TH1D* h_true = (TH1D*) h_t_from_truemass_incoherent_91->Clone("h_true");
	h_true->Rebin(2);
	h_true->Scale(0.5);
	h_true->Draw("hist same");

	//incoherent wrong mass
	TH1D* h_t_from_wrongmass_incoherent_91 = (TH1D*) h_VM_t_mass[0][1][method][1]->ProjectionY("h_t_from_wrongmass_incoherent_91",bin_lower,bin_upper);
	TH1D* h_t_from_wrongmass_incoherent_93 = (TH1D*) h_VM_t_mass[1][1][method][1]->ProjectionY("h_t_from_wrongmass_incoherent_93",bin_lower,bin_upper);
	h_t_from_wrongmass_incoherent_91->Add(h_t_from_wrongmass_incoherent_93,+1);
	measureXsection(name, h_t_from_wrongmass_incoherent_91, 0, t_hat_all->GetEntries(), PHP_);
	h_t_from_wrongmass_incoherent_91->Scale(1./scale_factor);//number coming from integral ratio for t>0.0;
	h_t_from_wrongmass_incoherent_91->SetMarkerStyle(26);
	h_t_from_wrongmass_incoherent_91->SetMarkerColor(kBlue);
	//add beam pipe effect.
	for(int ibin=0;ibin<h_t_from_wrongmass_incoherent_91->GetNbinsX();ibin++){
		double bincenter = h_t_from_wrongmass_incoherent_91->GetBinCenter(ibin+1);
		double weight = h_beam_pipe->GetBinContent( h_beam_pipe->FindBin(bincenter) );
		if(veto_){
			h_t_from_wrongmass_incoherent_91->SetBinContent(ibin+1, h_t_from_wrongmass_incoherent_91->GetBinContent(ibin+1)*weight);
			h_t_from_wrongmass_incoherent_91->SetBinError(ibin+1, h_t_from_wrongmass_incoherent_91->GetBinError(ibin+1)*weight);
		}
	}
	h_t_from_wrongmass_incoherent_91->SetMarkerColor(kRed);
	TH1D* h_wrong = histogramSubtraction(h_t_from_wrongmass_incoherent_91,h_t_from_truemass_incoherent_91);
	h_wrong->Add(h_t_from_mass_coherent_wrong_rho,+1);
	h_wrong->Rebin(2);
	h_wrong->Scale(0.5);
	h_wrong->Draw("PE same");

	//incoherent PID mass
	TH1D* h_t_from_PIDmass_incoherent_91 = (TH1D*) h_VM_t_mass[0][1][method][2]->ProjectionY("h_t_from_PIDmass_incoherent_91",bin_lower,bin_upper);
	TH1D* h_t_from_PIDmass_incoherent_93 = (TH1D*) h_VM_t_mass[1][1][method][2]->ProjectionY("h_t_from_PIDmass_incoherent_93",bin_lower,bin_upper);
	h_t_from_PIDmass_incoherent_91->Add(h_t_from_PIDmass_incoherent_93,+1);
	measureXsection(name, h_t_from_PIDmass_incoherent_91, 0, t_hat_all->GetEntries(), PHP_);
	h_t_from_PIDmass_incoherent_91->Scale(1./scale_factor);//number coming from integral ratio for t>0.0;
	h_t_from_PIDmass_incoherent_91->SetMarkerStyle(25);
	h_t_from_PIDmass_incoherent_91->SetMarkerColor(kBlue);
	//add beam pipe effect.
	for(int ibin=0;ibin<h_t_from_PIDmass_incoherent_91->GetNbinsX();ibin++){
		double bincenter = h_t_from_PIDmass_incoherent_91->GetBinCenter(ibin+1);
		double weight = h_beam_pipe->GetBinContent( h_beam_pipe->FindBin(bincenter) );
		if(veto_){
			h_t_from_PIDmass_incoherent_91->SetBinContent(ibin+1, h_t_from_PIDmass_incoherent_91->GetBinContent(ibin+1)*weight);
			h_t_from_PIDmass_incoherent_91->SetBinError(ibin+1, h_t_from_PIDmass_incoherent_91->GetBinError(ibin+1)*weight);
		}
		
	}
	h_t_from_PIDmass_incoherent_91->SetMarkerColor(kBlack);
	TH1D* h_pid = histogramSubtraction(h_t_from_PIDmass_incoherent_91, h_t_from_truemass_incoherent_91);
	h_pid->Add(h_t_from_mass_coherent_pid_rho,+1);
	h_pid->Rebin(2);
	h_pid->Scale(0.5);
	h_pid->Draw("PE same");

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

	TLatex* r44 = new TLatex(0.18, 0.84, "1 < Q^{2} < 20 GeV^{2}, |y_{VM}| < 4.0");
	r44->SetNDC();
	r44->SetTextSize(20);
	r44->SetTextFont(43);
	r44->SetTextColor(kBlack);

	TLatex* r44_1 = new TLatex(0.18, 0.84, "Q^{2} < 0.2 GeV^{2}, |y_{VM}| < 4.0");
	r44_1->SetNDC();
	r44_1->SetTextSize(20);
	r44_1->SetTextFont(43);
	r44_1->SetTextColor(kBlack);
	if(PHP_) r44_1->Draw("same");
	else r44->Draw("same");

	TLatex* r44_2 = new TLatex(0.18, 0.79, Form("|#eta_{daug.}| < 1.0, p_{T,daug.} > %.2f GeV/c",minPt_) );
	r44_2->SetNDC();
	r44_2->SetTextSize(20);
	r44_2->SetTextFont(43);
	r44_2->SetTextColor(kBlack);
	r44_2->Draw("same");
	
	TString method_name = "Method L";
	if(method==0){method_name = "Method E";}
	if(method==1){method_name = "Method A";}

	TLegend *w6 = new TLegend(0.5,0.58,0.72,0.73);
	w6->SetLineColor(kWhite);
	w6->SetFillColor(0);
	w6->SetTextSize(16);
	w6->SetTextFont(45);
	w6->AddEntry(h_coh_sartre, "Sartre coherent #phi truth ", "L");
	w6->AddEntry(h_true, "BeAGLE incoherent #phi truth ", "L");
	w6->AddEntry(h_wrong, "#rho^{0} contamination", "P");
	w6->AddEntry(h_pid, "#rho^{0} residue with PID |#eta|<1", "P");
	w6->Draw("same");

	TLegend *w7 = new TLegend(0.18,0.18,0.45,0.28);
	w7->SetLineColor(kWhite);
	w7->SetFillColor(0);
	w7->SetTextSize(15);
	w7->SetTextFont(45);
	w7->AddEntry(h_t_from_truemass_incoherent_91, "BeAGLE incoherent #phi truth ", "P");
	w7->AddEntry(h_wrong, "BeAGLE incoherent #rho^{0} contamination", "P");
	w7->AddEntry(h_pid, "BeAGLE incoherent #rho^{0} residue with PID", "P");


	c11->Print(Form("../figures/proposal_ideal/PID_proposal_figure_PHP_%d_veto_%d_minPt_%.2f.pdf",PHP_,veto_,minPt_));

}