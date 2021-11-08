#include "utility.h"
void plotProposalFigure_1_ratio_Bfield(TString name="phi", int smear_ = 0, int PHP_ = 0, int veto_ = 1, double minPt_=0.15, int method=2){

	TFile* beam_piple_file = new TFile("../analysis/include/beam_pipe.root");
	TH1D* h_beam_pipe = (TH1D*) beam_piple_file->Get("beam_pipe_factor");

	int sartre_vm_index=1;
	int beagle_vm_index=0;
	int processindex=1;//91,93
	double scale_factor = 0.822;
	if(name=="phi_photo") {scale_factor=0.335;PHP_=1;}

	TString inputROOT=Form("../rootfiles/beagle_output_PHP_%d_veto_%d_minPt_%.2f_smear_%d.root",PHP_,veto_,minPt_,smear_);
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
		file_sartre_all[0] = new TFile(Form("../rootfiles/sartre_rho_photo_bnonsat_PID_1_minPt_%.2f_smear_0.root",minPt_));
		file_sartre_all[1] = new TFile(Form("../rootfiles/sartre_phi_photo_bnonsat_PID_1_minPt_%.2f_smear_0.root",minPt_));
		file_sartre_all[2] = new TFile(Form("../rootfiles/sartre_jpsi_photo_bnonsat_PID_1_minPt_%.2f_smear_0.root",minPt_));
	}
	else{
		file_sartre_all[0] = new TFile(Form("../rootfiles/sartre_rho_bnonsat_PID_1_minPt_%.2f_smear_%d.root",minPt_,smear_));
		file_sartre_all[1] = new TFile(Form("../rootfiles/sartre_phi_bnonsat_PID_1_minPt_%.2f_smear_%d.root",minPt_,smear_));
		file_sartre_all[2] = new TFile(Form("../rootfiles/sartre_jpsi_bnonsat_PID_1_minPt_%.2f_smear_0.root",minPt_));
	}	

	TFile* file_sartre = new TFile(Form("../rootfiles/sartre_"+name+"_bnonsat_PID_1_minPt_%.2f_smear_0.root",minPt_));
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
    for(int ibreak=0;ibreak<2;ibreak++){
        for(int imethod=0;imethod<3;imethod++){
        	for(int imass=0;imass<3;imass++){
				h_t_reco_sartre[ibreak][imethod][imass] = (TH1D*) file_sartre->Get(Form("h_t_reco_%d_%d_%d",ibreak,imethod,imass));
				h_VM_t_mass_sartre[ibreak][imethod][imass] = (TH2D*) file_sartre->Get(Form("h_VM_t_mass_%d_%d_%d",ibreak,imethod,imass));
        		for(int ifile=0;ifile<3;ifile++){
	    			h_VM_t_mass_sartre_mixed[ibreak][imethod][imass][ifile] = (TH2D*) file_sartre_all[ifile]->Get(Form("h_VM_t_mass_%d_%d_%d",ibreak,imethod,imass));
	    		}
        	}
        }
    }
 

	TCanvas* c1 = new TCanvas("c1","c1",1,1,600,800);
	TPad* drawPad[2];
	drawPad[0] = new TPad("pad0","pad0",0.0,0.3,1.0,1.0);
	drawPad[1] = new TPad("pad1","pad1",0.0,0.0,1.0,0.29);
	
	drawPad[0]->SetLeftMargin(0.13);
	drawPad[0]->SetRightMargin(0.01);
	drawPad[0]->SetTopMargin(0.1);
	drawPad[0]->SetBottomMargin(0);
	drawPad[0]->Draw();

	drawPad[1]->SetLeftMargin(0.13);
	drawPad[1]->SetRightMargin(0.01);
	drawPad[1]->SetTopMargin(0.);
	drawPad[1]->SetBottomMargin(0.35);
	drawPad[1]->Draw();

	TPad* smallPad = new TPad("smallPad","smallPad",0.65,0.56,0.95,0.72);
	// smallPad->Draw();

	TH1D* base11 = makeHist("base11", "", "|#it{t} | (GeV^{2})", "d#sigma / d|#it{t} |  (nb/GeV^{2}) ", 100,0,0.102,kBlack);
	base11->GetYaxis()->SetRangeUser(4e-2, 1e8);
	if(name=="phi"){base11->GetYaxis()->SetRangeUser(2e-0, 5e6);}
	if(name=="phi_photo"){base11->GetYaxis()->SetRangeUser(4e-2, 1e9);}
	base11->GetXaxis()->SetTitleColor(kBlack);
	TGaxis::SetMaxDigits(3);
	fixedFontHist1D(base11,1.2,1.9);
	base11->GetYaxis()->SetTitleSize(base11->GetYaxis()->GetTitleSize()*1.3);
	base11->GetXaxis()->SetTitleSize(base11->GetXaxis()->GetTitleSize()*1.3);
	base11->GetYaxis()->SetLabelSize(base11->GetYaxis()->GetLabelSize()*1.5);
	base11->GetXaxis()->SetLabelSize(base11->GetXaxis()->GetLabelSize()*1.5);
	base11->GetXaxis()->SetNdivisions(5,5,0);
	base11->GetYaxis()->SetNdivisions(5,5,0);
	drawPad[0]->cd();
	gPad->SetLogy(1);
	base11->Draw();

	measureXsection(name, h_coh_sartre, 1);
	h_coh_sartre->SetLineColor(kBlack);
	h_coh_sartre->SetLineWidth(1);
	h_coh_sartre->Rebin(2);
	h_coh_sartre->Scale(0.5);
	h_coh_sartre->Draw("hist same");

	measureXsection(name, hist_t_afterPhaseSpace_coherent_sartre, 1);
	hist_t_afterPhaseSpace_coherent_sartre->SetLineColor(kBlue);
	// hist_t_afterPhaseSpace_coherent_sartre->Draw("hist same");

	//coherent
	TH1D* h_mass_for_binning = (TH1D*) h_VM_t_mass_sartre[0][method][1]->ProjectionX("h_mass_for_binning",1,1000);
	int bin_lower = h_mass_for_binning->FindBin(1.019-0.02);
	int bin_upper = h_mass_for_binning->FindBin(1.019+0.02);
	
	//use mixed
	//residue using kk mass on rho vm. 
	TH1D* h_t_from_mass_coherent_mixed_rho = (TH1D*) h_VM_t_mass_sartre_mixed[0][method][sartre_vm_index][0]->ProjectionY("h_t_from_mass_coherent_mixed_rho",bin_lower,bin_upper);
	TString name_temp = "rho";
	if(PHP_) name_temp = "rho_photo";
	measureXsection(name_temp, h_t_from_mass_coherent_mixed_rho, 1);
	// h_t_from_mass_coherent_mixed_rho->Draw("same");
	//true phi vm. 
	TH1D* h_t_from_mass_coherent_mixed_phi = (TH1D*) h_VM_t_mass_sartre_mixed[0][method][sartre_vm_index][1]->ProjectionY("h_t_from_mass_coherent_mixed_phi",bin_lower,bin_upper);
	name_temp = "phi";
	if(PHP_) name_temp = "phi_photo";
	measureXsection(name_temp, h_t_from_mass_coherent_mixed_phi, 1);
	h_t_from_mass_coherent_mixed_phi->SetLineColor(kBlack);
	h_t_from_mass_coherent_mixed_phi->SetMarkerStyle(25);
	h_t_from_mass_coherent_mixed_phi->Add(h_t_from_mass_coherent_mixed_rho,+1);//add both together
	// h_t_from_mass_coherent_mixed_phi->Draw();

	//incoherent pid mass
	TH1D* h_t_from_PIDmass_incoherent_91 = (TH1D*) h_VM_t_mass[0][vm_index][method][2]->ProjectionY("h_t_from_PIDmass_incoherent_91",bin_lower,bin_upper);
	TH1D* h_t_from_PIDmass_incoherent_93 = (TH1D*) h_VM_t_mass[1][vm_index][method][2]->ProjectionY("h_t_from_PIDmass_incoherent_93",bin_lower,bin_upper);
	h_t_from_PIDmass_incoherent_91->Add(h_t_from_PIDmass_incoherent_93,+1);
	measureXsection(name, h_t_from_PIDmass_incoherent_91, 0, t_hat_all->GetEntries(), PHP_);
	h_t_from_PIDmass_incoherent_91->SetMarkerStyle(24);
	h_t_from_PIDmass_incoherent_91->SetMarkerColor(kBlue);
	h_t_from_PIDmass_incoherent_91->Scale(1./scale_factor);//number coming from integral ratio for t>0.0;
	//add beam pipe effect.
	for(int ibin=0;ibin<h_t_from_PIDmass_incoherent_91->GetNbinsX();ibin++){
		double bincenter = h_t_from_PIDmass_incoherent_91->GetBinCenter(ibin+1);
		double weight = h_beam_pipe->GetBinContent( h_beam_pipe->FindBin(bincenter) );
		if(veto_ && smear_){
			h_t_from_PIDmass_incoherent_91->SetBinContent(ibin+1, h_t_from_PIDmass_incoherent_91->GetBinContent(ibin+1)*weight);
			h_t_from_PIDmass_incoherent_91->SetBinError(ibin+1, h_t_from_PIDmass_incoherent_91->GetBinError(ibin+1)*weight);
		}		
	}
	/*
	Start to play with subtraction. 
	*/
	TFile* file_temp = 0;
	if(smear_){
		file_temp = new TFile(Form("veto_extrapolation_PHP_%d_wBeamPipe.root",PHP_));
	}
	else{
		file_temp = new TFile(Form("veto_extrapolation_PHP_%d_wNoBeamPipe.root",PHP_));
	}
	
	TH1D* h_tagged = (TH1D*) file_temp->Get("ratio");
	TH1D* h_tagged_p10 = (TH1D*) h_tagged->Clone("ratio_p10");
	TH1D* h_tagged_m10 = (TH1D*) h_tagged->Clone("ratio_m10");
	
	h_tagged->Scale(h_t_from_PIDmass_incoherent_91->Integral()*1.0);
	h_tagged_p10->Scale(h_t_from_PIDmass_incoherent_91->Integral()*1.1);
	h_tagged_m10->Scale(h_t_from_PIDmass_incoherent_91->Integral()*0.9);

	TH1D* h_total = (TH1D*) h_t_from_mass_coherent_mixed_phi->Clone("h_total");
	h_total->Add(h_t_from_PIDmass_incoherent_91,+1);
	h_total->SetFillColorAlpha(kRed,0.4);
    h_total->SetFillStyle(1001);
	h_total->SetMarkerStyle(25);
	h_total->SetMarkerColor(kRed);
	// h_total->Draw("P E1 same");
	
	TH1D* h_total_subtract = (TH1D*) h_total->Clone("h_total_subtract");
	TH1D* h_total_subtract_p10 = (TH1D*) h_total->Clone("h_total_subtract_p10");
	TH1D* h_total_subtract_m10 = (TH1D*) h_total->Clone("h_total_subtract_m10");
	
	h_total_subtract->Add(h_tagged,-1);
	h_total_subtract_p10->Add(h_tagged_p10,-1);
	h_total_subtract_m10->Add(h_tagged_m10,-1);

	h_total_subtract->SetMarkerColor(kBlack);
	h_total_subtract->SetMarkerStyle(24);
	
	TGraphAsymmErrors* gr = new TGraphAsymmErrors();
	for(int ibin=0;ibin<h_total_subtract->GetNbinsX();ibin++){
		double value = h_total_subtract->GetBinContent(ibin+1);
		double center = h_total_subtract->GetBinCenter(ibin+1);
		double value_1 = h_total_subtract_p10->GetBinContent(ibin+1);
		double value_2 = h_total_subtract_m10->GetBinContent(ibin+1);
		gr->SetPoint(ibin, center, value);
		gr->SetPointEYhigh(ibin, value_2-value);
		gr->SetPointEYlow(ibin, value-value_1);
	}
	gr->SetFillColorAlpha(kBlue,0.4);
    gr->SetFillStyle(1001);
	gr->SetMarkerColor(kBlue);
	gr->SetMarkerStyle(24);
	gr->Draw("P C e3 same");
	h_total_subtract->Draw("P e1 same");

	TFile* input_file = new TFile("proposal_1.5T.root");
	TH1D* h_1p5T = (TH1D*) input_file->Get("h_1p5T");
	TGraphErrors* gr_1p5T = (TGraphErrors*) input_file->Get("gr_1p5T");
	gr_1p5T->SetFillColorAlpha(kGray+2,0.4);
    gr_1p5T->SetFillStyle(1001);
	gr_1p5T->SetMarkerColor(kGray+2);
	gr_1p5T->SetMarkerStyle(25);
	gr_1p5T->SetMarkerSize(1.);
	gr_1p5T->Draw("P e3 same");

	TFile* input_file2 = new TFile("proposal_0.75T.root");
	TH1D* h_0p75T = (TH1D*) input_file2->Get("h_0p75T");
	TGraphErrors* gr_0p75T = (TGraphErrors*) input_file2->Get("gr_0p75T");
	gr_0p75T->SetFillColorAlpha(kRed,0.4);
    gr_0p75T->SetFillStyle(1001);
	gr_0p75T->SetMarkerColor(kRed);
	gr_0p75T->SetMarkerStyle(28);
	gr_0p75T->SetMarkerSize(1.2);
	gr_0p75T->Draw("P e3 same");


	TLatex* r42 = new TLatex(0.18, 0.91, "eAu 18x110 GeV^{2}");
	r42->SetNDC();
	r42->SetTextSize(22);
	r42->SetTextFont(43);
	r42->SetTextColor(kBlack);
	r42->Draw("same");

	TLatex* r43 = new TLatex(0.7,0.91, "ATHENA #it{Internal}");
	r43->SetNDC();
	r43->SetTextSize(0.04);
	r43->Draw("same");

	TLatex* r44 = new TLatex(0.18, 0.84, "1 < Q^{2} < 20 GeV^{2}, |y_{#phi}| < 4.0");
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

	TLatex* r44_2 = new TLatex(0.18, 0.79, Form("|#eta_{kaon}| < 4.0, p_{T,kaon} > %.2f GeV/c",minPt_) );
	r44_2->SetNDC();
	r44_2->SetTextSize(20);
	r44_2->SetTextFont(43);
	r44_2->SetTextColor(kBlack);
	r44_2->Draw("same");
	
	TString method_name = "Method L";
	if(method==0){method_name = "Method E";}
	if(method==1){method_name = "Method A";}

	TLegend *w6 = new TLegend(0.53,0.52,0.9,0.67);
	w6->SetLineColor(kWhite);
	w6->SetFillColor(0);
	w6->SetTextSize(17);
	w6->SetTextFont(45);
	w6->AddEntry(h_coh_sartre, "Sartre coherent #phi truth ", "PL");
	w6->AddEntry(gr_0p75T, "0.75T: reco'. #phi "+method_name+"", "PE");
	w6->AddEntry(gr_1p5T, "1.5T: reco'. #phi "+method_name+"", "PE");
	w6->AddEntry(gr, "3.0T: reco'. #phi "+method_name+"" , "PE");
	w6->Draw("same");

	TBox* box1 = new TBox(0.088-0.038,2.E4,0.098-0.042,0.33E5);
    box1->SetFillColorAlpha(kGray+2,0.3);
    box1->SetFillStyle(1001);
	box1->SetLineWidth(0);
	box1->Draw("same");

	TBox* box2 = new TBox(0.088-0.038,0.96E4,0.098-0.042,1.8E4);
    box2->SetFillColorAlpha(kBlue,0.3);
    box2->SetFillStyle(1001);
	box2->SetLineWidth(0);
	box2->Draw("same");

	TBox* box3 = new TBox(0.088-0.038,3.5E4,0.098-0.042,6E4);
    box3->SetFillColorAlpha(kRed,0.3);
    box3->SetFillStyle(1001);
	box3->SetLineWidth(0);
	box3->Draw("same");

	TLatex* r44_3 = new TLatex(0.2, 0.22, "zoom in" );
	r44_3->SetNDC();
	r44_3->SetTextSize(14);
	r44_3->SetTextFont(43);
	r44_3->SetTextColor(kBlack);
	// r44_3->Draw("same");

	TH1D* base22 = makeHist("base22", "", "|#it{t} | (GeV^{2})", "reco' / truth ", 100,0,0.102,kBlack);
	base22->GetYaxis()->SetRangeUser(0.1, 80);
	base22->GetXaxis()->SetTitleColor(kBlack);
	TGaxis::SetMaxDigits(3);
	fixedFontHist1D(base22,4,1.9);
	base22->GetYaxis()->SetTitleSize(base22->GetYaxis()->GetTitleSize()*1.3);
	base22->GetXaxis()->SetTitleSize(base22->GetXaxis()->GetTitleSize()*1.3);
	base22->GetYaxis()->SetLabelSize(base22->GetYaxis()->GetLabelSize()*1.5);
	base22->GetXaxis()->SetLabelSize(base22->GetXaxis()->GetLabelSize()*1.5);
	base22->GetXaxis()->SetNdivisions(5,5,0);
	base22->GetYaxis()->SetNdivisions(5,5,0);

	drawPad[1]->cd();
	gPad->SetLogy(1);
	base22->Draw("");
	TH1D* h_ratio = (TH1D*) h_total_subtract->Clone("h_ratio");
	h_ratio->SetFillColorAlpha(kBlue,0.4);
    h_ratio->SetFillStyle(1001);
	h_ratio->SetMarkerStyle(25);
	h_ratio->SetMarkerColor(kBlue);
	h_ratio->SetLineColor(kBlue);
	h_ratio->SetMarkerStyle(24);
	h_ratio->Divide(h_coh_sartre);
	h_ratio->Draw("PE1 same");

	TH1D* h_ratio_1p5T = (TH1D*) h_1p5T->Clone("h_ratio");
	h_ratio_1p5T->SetFillColorAlpha(kGray+2,0.4);
    h_ratio_1p5T->SetFillStyle(1001);
	h_ratio_1p5T->SetMarkerColor(kGray+2);
	h_ratio_1p5T->SetLineColor(kGray+2);
	h_ratio_1p5T->SetMarkerStyle(25);
	h_ratio_1p5T->Divide(h_coh_sartre);
	h_ratio_1p5T->Draw("PE1 same");

	TH1D* h_ratio_0p75T = (TH1D*) h_0p75T->Clone("h_ratio");
	h_ratio_0p75T->SetFillColorAlpha(kRed,0.4);
    h_ratio_0p75T->SetFillStyle(1001);
	h_ratio_0p75T->SetMarkerColor(kRed);
	h_ratio_0p75T->SetLineColor(kRed);
	h_ratio_0p75T->SetMarkerStyle(28);
	h_ratio_0p75T->SetMarkerSize(1.2);
	h_ratio_0p75T->Divide(h_coh_sartre);
	h_ratio_0p75T->Draw("PE1 same");

	// TFile* output_file = new TFile("proposal_1.5T.root","RECREATE");
	// h_total_subtract->SetName("h_1p5T");
	// h_total_subtract->Write();
	// gr_copy->SetName("gr_1p5T");
	// gr_copy->Write();

	c1->Print(Form("../figures/proposal_ideal/proposal_figure_PHP_%d_veto_%d_minPt_%.2f_smear_%d_1.5T.pdf",PHP_,veto_,minPt_,smear_));

}