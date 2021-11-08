#include "utility.h"
void plotProposalFigure_1_b(TString name="phi",int PHP_ = 0, int veto_ = 1, double minPt_=0.15, int method=2){

	TFile* beam_piple_file = new TFile("../analysis/include/beam_pipe.root");
	TH1D* h_beam_pipe = (TH1D*) beam_piple_file->Get("beam_pipe_factor");

	int sartre_vm_index=1;
	int beagle_vm_index=0;
	int processindex=1;//91,93
	double scale_factor = 0.822;
	if(name=="phi_photo") {scale_factor=0.335;PHP_=1;}

	//spectial here:
	TString inputROOT_noVeto=Form("../rootfiles/beagle_output_PHP_%d_veto_0_minPt_%.2f_smear_0.root",PHP_,minPt_);
	TFile* file_beagle_noVeto = new TFile(inputROOT_noVeto);
	//end

	TString inputROOT=Form("../rootfiles/beagle_output_PHP_%d_veto_%d_minPt_%.2f_smear_0.root",PHP_,veto_,minPt_);
	TFile* file_beagle = new TFile(inputROOT);
	TH1D* t_hat_all = (TH1D*) file_beagle->Get("h_trueT");

	//beagle
	TH2D* h_VM_t_mass[3][3][3][3];
	TH2D* h_VM_t_mass_noVeto[3][3][3][3];
	for(int ibreak=0;ibreak<3;ibreak++){
		for(int ivm=0;ivm<3;ivm++){
			for(int imethod=0;imethod<3;imethod++){
				for(int imass=0;imass<3;imass++){
					h_VM_t_mass_noVeto[ibreak][ivm][imethod][imass] = (TH2D*) file_beagle_noVeto->Get(Form("h_VM_t_mass_%d_%d_%d_%d",ibreak,ivm,imethod,imass));
					h_VM_t_mass[ibreak][ivm][imethod][imass] = (TH2D*) file_beagle->Get(Form("h_VM_t_mass_%d_%d_%d_%d",ibreak,ivm,imethod,imass));
				}
			}
		}
	}

	/* Sartre */
	TFile* file_sartre = new TFile(Form("../rootfiles/sartre_"+name+"_bnonsat_PID_1_minPt_%.2f_smear_0.root",minPt_));
	//not use here
	TH1D* h_coh_sartre = (TH1D*) file_sartre->Get("hist_t_coherent");
	TH1D* h_incoh_sartre = (TH1D*) file_sartre->Get("hist_t_incoherent");
	TH1D* hist_t_afterPhaseSpace_coherent_sartre = (TH1D*) file_sartre->Get("hist_t_afterPhaseSpace_coherent");
	TH1D* hist_t_afterPhaseSpace_incoherent_sartre = (TH1D*) file_sartre->Get("hist_t_afterPhaseSpace_incoherent");
	//end not use here
	// sartre
    TH2D* h_VM_t_mass_sartre[2][3][3];
    for(int ibreak=0;ibreak<2;ibreak++){
        for(int imethod=0;imethod<3;imethod++){
        	for(int imass=0;imass<3;imass++){
				h_VM_t_mass_sartre[ibreak][imethod][imass] = (TH2D*) file_sartre->Get(Form("h_VM_t_mass_%d_%d_%d",ibreak,imethod,imass));
      
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
	base11->GetYaxis()->SetRangeUser(1e-1, 1e7);
	if(name=="phi"){base11->GetYaxis()->SetRangeUser(1e-1, 1e7);}
	if(name=="phi_photo"){base11->GetYaxis()->SetRangeUser(1e-1, 1e7);}
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
	int bin_lower = h_mass_for_binning->FindBin(1.019-0.02);
	int bin_upper = h_mass_for_binning->FindBin(1.019+0.02);

	/*
	Incoherent contributions, plot separately.
	*/

	//incoherent true mass - phi incoherent.
	TH1D* h_t_from_truemass_incoherent_91 = (TH1D*) h_VM_t_mass[0][1][method][0]->ProjectionY("h_t_from_truemass_incoherent_91",bin_lower,bin_upper);
	TH1D* h_t_from_truemass_incoherent_93 = (TH1D*) h_VM_t_mass[1][1][method][0]->ProjectionY("h_t_from_truemass_incoherent_93",bin_lower,bin_upper);
	h_t_from_truemass_incoherent_91->Add(h_t_from_truemass_incoherent_93,+1);
	measureXsection(name, h_t_from_truemass_incoherent_91, 0, t_hat_all->GetEntries(), PHP_);
	h_t_from_truemass_incoherent_91->SetMarkerStyle(24);
	h_t_from_truemass_incoherent_91->SetMarkerColor(kBlue);
	h_t_from_truemass_incoherent_91->Scale(1./scale_factor);//number coming from integral ratio for t>0.0;
	// for(int ibin=0;ibin<h_t_from_truemass_incoherent_91->GetNbinsX();ibin++){
	// 	double bincenter = h_t_from_truemass_incoherent_91->GetBinCenter(ibin+1);
	// 	double weight = h_beam_pipe->GetBinContent( h_beam_pipe->FindBin(bincenter) );
	// 	if(veto_){
	// 		h_t_from_truemass_incoherent_91->SetBinContent(ibin+1, h_t_from_truemass_incoherent_91->GetBinContent(ibin+1)*weight);
	// 		h_t_from_truemass_incoherent_91->SetBinError(ibin+1, h_t_from_truemass_incoherent_91->GetBinError(ibin+1)*weight);
	// 	}		
	// }
	h_t_from_truemass_incoherent_91->Draw("PE same");

	TH1D* h_noVeto_91 = (TH1D*) h_VM_t_mass_noVeto[0][1][method][0]->ProjectionY("h_noVeto_91",bin_lower,bin_upper);
	TH1D* h_noVeto_93 = (TH1D*) h_VM_t_mass_noVeto[1][1][method][0]->ProjectionY("h_noVeto_93",bin_lower,bin_upper);
	h_noVeto_91->Add(h_noVeto_93,+1);
	measureXsection(name, h_noVeto_91, 0, t_hat_all->GetEntries(), PHP_);
	h_noVeto_91->SetMarkerStyle(24);
	h_noVeto_91->SetMarkerColor(kBlue);
	// h_noVeto_91->Rebin(2);
	// h_noVeto_91->Scale(0.5);
	h_noVeto_91->Scale(1./scale_factor);//number coming from integral ratio for t>0.0;
	h_noVeto_91->SetLineWidth(2);
	h_noVeto_91->SetLineStyle(2);

	h_noVeto_91->Draw("hist same");
	TH1D* h_tagged = histogramSubtraction(h_noVeto_91,h_t_from_truemass_incoherent_91);
	h_tagged->SetMarkerColor(kRed);
	h_tagged->Draw("PE same");

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

	TLatex* r44_2 = new TLatex(0.18, 0.79, Form("|#eta_{daug.}| < 4.0, p_{T,daug.} > %.2f GeV/c",minPt_) );
	r44_2->SetNDC();
	r44_2->SetTextSize(20);
	r44_2->SetTextFont(43);
	r44_2->SetTextColor(kBlack);
	r44_2->Draw("same");
	
	TString method_name = "Method L";
	if(method==0){method_name = "Method E";}
	if(method==1){method_name = "Method A";}

	TLegend *w6 = new TLegend(0.38,0.58,0.75,0.73);
	w6->SetLineColor(kWhite);
	w6->SetFillColor(0);
	w6->SetTextSize(13);
	w6->SetTextFont(45);
	w6->AddEntry(h_coh_sartre, "Sartre coherent #phi truth ", "L");
	w6->AddEntry(h_noVeto_91, "BeAGLE incoherent #phi truth", "L");
	w6->AddEntry(h_tagged, "BeAGLE incoherent #phi tagged", "P");
	w6->AddEntry(h_t_from_truemass_incoherent_91, "BeAGLE incoherent #phi residue ", "P");
	w6->Draw("same");

	TCanvas*c2 = new TCanvas("c2","c2",1,1,600,600);
	gPad->SetLogy(1);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	TH1D* base2 = makeHist("base2", "", "|#it{t} | (GeV^{2})", "Self Normalized #it{t} ", 100,0,0.18,kBlack);
	base2->GetYaxis()->SetRangeUser(1e-3, 1e0);
	base2->GetXaxis()->SetTitleColor(kBlack);
	TGaxis::SetMaxDigits(3);
	fixedFontHist1D(base2,1.,1.2);
	base2->GetYaxis()->SetTitleSize(base2->GetYaxis()->GetTitleSize()*1.5);
	base2->GetXaxis()->SetTitleSize(base2->GetXaxis()->GetTitleSize()*1.5);
	base2->GetYaxis()->SetLabelSize(base2->GetYaxis()->GetLabelSize()*1.5);
	base2->GetXaxis()->SetLabelSize(base2->GetXaxis()->GetLabelSize()*1.5);
	base2->GetXaxis()->SetNdivisions(4,4,0);
	base2->GetYaxis()->SetNdivisions(5,5,0);
	base2->Draw();
	TH1D* ratio = (TH1D*) h_tagged->Clone("ratio");
	ratio->Scale( 1./ratio->Integral() );
	TH1D* ratio_d = (TH1D*) h_t_from_truemass_incoherent_91->Clone("ratio_d");
	ratio_d->Scale( 1./ratio_d->Integral() );
	// ratio->Divide(ratio_d);
	ratio_d->Draw("P same ");
	ratio->Draw("P same");

	TLegend *w7 = new TLegend(0.38,0.68,0.75,0.75);
	w7->SetLineColor(kWhite);
	w7->SetFillColor(0);
	w7->SetTextSize(15);
	w7->SetTextFont(45);
	w7->AddEntry(ratio, "BeAGLE incoherent #phi tagged  ", "P");
	w7->AddEntry(ratio_d, "BeAGLE incoherent #phi residue", "P");
	w7->Draw("same");

	TFile* output = new TFile(Form("veto_extrapolation_PHP_%d.root",PHP_),"RECREATE");
	ratio->Write();
	ratio_d->Write();

	// c11->Print(Form("../figures/proposal_ideal/PID_proposal_figure_PHP_%d_veto_%d_minPt_%.2f.pdf",PHP_,veto_,minPt_));

}