#include "utility.h"
void plotdSigmadt_RECO(TString name="phi", int method=0){

	/* Beagle */
	
	TFile* file_beagle = new TFile("../rootfiles/beagle_allVMs_w_breakups.root");
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
	//beagle
	TH1D* h_t_reco[3][3][3];
	TH2D* h_VM_t_mass[3][3][3];
	for(int ibreak=0;ibreak<3;ibreak++){
		for(int ivm=0;ivm<3;ivm++){
			for(int imethod=0;imethod<3;imethod++){
				h_t_reco[ibreak][ivm][imethod] = (TH1D*) file_beagle->Get(Form("h_t_reco_%d_%d_%d",ibreak,ivm,imethod));
				h_VM_t_mass[ibreak][ivm][imethod] = (TH2D*) file_beagle->Get(Form("h_VM_t_mass_%d_%d_%d",ibreak,ivm,imethod));
			}
		}
	}


	/* Sartre */

	TFile* file_sartre = new TFile("../rootfiles/sartre_"+name+"_bnonsat.root");
	TH1D* h_coh_sartre = (TH1D*) file_sartre->Get("hist_t_coherent");
	TH1D* h_incoh_sartre = (TH1D*) file_sartre->Get("hist_t_incoherent");
	// sartre
    TH1D* h_t_reco_sartre[2][3];
    TH2D* h_VM_t_mass_sartre[2][3];
    for(int ibreak=0;ibreak<2;ibreak++){
        for(int imethod=0;imethod<3;imethod++){
            h_t_reco_sartre[ibreak][imethod] = (TH1D*) file_sartre->Get(Form("h_t_reco_%d_%d",ibreak,imethod));
			h_VM_t_mass_sartre[ibreak][imethod] = (TH2D*) file_sartre->Get(Form("h_VM_t_mass_%d_%d",ibreak,imethod));
        }
    }

	TCanvas* c1 = new TCanvas("c1","c1",1,1,600,600);
	gPad->SetLogy(1);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	TH1D* base1 = makeHist("base1", "", "|#it{t} | (GeV^{2})", "d#sigma/d|#it{t} | (nb/GeV^{2}) ", 100,0,0.18,kBlack);
	if(name=="rho"){base1->GetYaxis()->SetRangeUser(1e-1, 1e8);}
	if(name=="phi"){base1->GetYaxis()->SetRangeUser(1e-2, 1e7);}
	if(name=="jpsi"){base1->GetYaxis()->SetRangeUser(1e-3, 1e6);}
	
	base1->GetXaxis()->SetTitleColor(kBlack);
	TGaxis::SetMaxDigits(2);
	fixedFontHist1D(base1,1.,1.1);
	base1->GetYaxis()->SetTitleSize(base1->GetYaxis()->GetTitleSize()*1.5);
	base1->GetXaxis()->SetTitleSize(base1->GetXaxis()->GetTitleSize()*1.5);
	base1->GetYaxis()->SetLabelSize(base1->GetYaxis()->GetLabelSize()*1.5);
	base1->GetXaxis()->SetLabelSize(base1->GetXaxis()->GetLabelSize()*1.5);
	base1->GetXaxis()->SetNdivisions(4,4,0);
	base1->GetYaxis()->SetNdivisions(5,5,0);
	base1->Draw();

	measureXsection(name, h_t_reco_sartre[0][method], 1);
	h_t_reco_sartre[0][method]->SetLineColor(kBlue);
	h_t_reco_sartre[0][method]->SetLineWidth(2);
	h_t_reco_sartre[0][method]->Draw(" hist same");

	measureXsection(name, h_t_reco_sartre[1][method], 1);
	h_t_reco_sartre[1][method]->SetLineColor(kRed);
	h_t_reco_sartre[1][method]->SetLineWidth(2);
	h_t_reco_sartre[1][method]->Draw(" hist same");

	TH1D* h_VM_background = (TH1D*) h_t_reco[0][vm_index][method]->Clone("h_VM_background");
	h_VM_background->Add(h_t_reco[1][vm_index][method],+1);
	//adding elastic and dissoc. together.
	measureXsection(name, h_VM_background, 0, t_hat_all->GetEntries(),0);
	h_VM_background->SetMarkerStyle(24);
	h_VM_background->Draw("P SAME");

	TLegend *w5 = new TLegend(0.53,0.74,0.73,0.87);
	w5->SetLineColor(kWhite);
	w5->SetFillColor(0);
	w5->SetTextSize(18);
	w5->SetTextFont(45);
	w5->AddEntry(h_t_reco_sartre[0][method], "Sartre "+legendName+" coherent  ", "PL");
	w5->AddEntry(h_t_reco_sartre[1][method], "Sartre "+legendName+" incoherent  ", "PL");
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
	r44->Draw("same");

	c1->Print("../figures/methods/dsigma_dt_"+name+"_Method"+std::to_string(method)+".pdf");

	TCanvas* c2 = new TCanvas("c2","c2",1,1,600,600);
	gPad->SetLogy(0);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	TH1D* base2 = (TH1D*) base1->Clone("base2");
	base2->GetYaxis()->SetRangeUser(1e-4,8e-2);
	base2->GetXaxis()->SetRangeUser(0,0.5);
	base1->GetXaxis()->SetTitleColor(kBlack);
	fixedFontHist1D(base2,1.,1.5);
	base2->GetYaxis()->SetTitleSize(base2->GetYaxis()->GetTitleSize()*1.5);
	base2->GetXaxis()->SetTitleSize(base2->GetXaxis()->GetTitleSize()*1.5);
	base2->GetYaxis()->SetLabelSize(base2->GetYaxis()->GetLabelSize()*1.5);
	base2->GetXaxis()->SetLabelSize(base2->GetXaxis()->GetLabelSize()*1.5);
	base2->GetXaxis()->SetNdivisions(4,4,0);
	base2->GetYaxis()->SetNdivisions(5,5,0);
	int processindex=0;
	if( processindex==0 ){
		base2->SetTitle("subprocess 91");
	}
	else{
		base2->SetTitle("subprocess 93");
	}
	base2->Draw();
	
	for(int imethod=0;imethod<3;imethod++){
		h_t_reco[processindex][vm_index][imethod]->SetMarkerStyle(24+imethod);
	}
	h_t_reco[processindex][vm_index][0]->SetMarkerColor(kBlack);
	h_t_reco[processindex][vm_index][1]->SetMarkerColor(kBlue);
	h_t_reco[processindex][vm_index][2]->SetMarkerColor(kRed);
	h_t_reco[processindex][vm_index][0]->DrawNormalized("PE same");
	h_t_reco[processindex][vm_index][1]->DrawNormalized("PE same");
	h_t_reco[processindex][vm_index][2]->DrawNormalized("PE same");

	TLegend *w6 = new TLegend(0.53,0.74,0.73,0.87);
	w6->SetLineColor(kWhite);
	w6->SetFillColor(0);
	w6->SetTextSize(18);
	w6->SetTextFont(45);
	w6->AddEntry(h_t_reco[processindex][vm_index][0], "BeAGLE "+legendName+" Method E  ", "PL");
	w6->AddEntry(h_t_reco[processindex][vm_index][1], "BeAGLE "+legendName+" Method A  ", "PL");
	w6->AddEntry(h_t_reco[processindex][vm_index][2], "BeAGLE "+legendName+" Method L ", "PL");
	w6->Draw("same");
	if(processindex==0){
		c2->Print("../figures/methods/beagle_"+name+"_process_91.pdf");
	}
	else{
		c2->Print("../figures/methods/beagle_"+name+"_process_93.pdf");
	}
	TCanvas* c3 = new TCanvas("c3","c3",1,1,600,600);
	gPad->SetLogy(1);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	TH1D* base3 = (TH1D*) base1->Clone("base3");
	base3->GetYaxis()->SetRangeUser(1e-7,4e-0);
	base3->GetXaxis()->SetRangeUser(0,0.5);
	base1->GetXaxis()->SetTitleColor(kBlack);
	fixedFontHist1D(base3,1.,1.5);
	base3->GetYaxis()->SetTitleSize(base3->GetYaxis()->GetTitleSize()*1.5);
	base3->GetXaxis()->SetTitleSize(base3->GetXaxis()->GetTitleSize()*1.5);
	base3->GetYaxis()->SetLabelSize(base3->GetYaxis()->GetLabelSize()*1.5);
	base3->GetXaxis()->SetLabelSize(base3->GetXaxis()->GetLabelSize()*1.5);
	base3->GetXaxis()->SetNdivisions(4,4,0);
	base3->GetYaxis()->SetNdivisions(5,5,0);
	int coh=1;
	if( coh==0 ){
		base3->SetTitle("coherent");
	}
	else{
		base3->SetTitle("incoherent");
	}
	base3->Draw();
	TH1D* hist_temp[3];
	for(int imethod=0;imethod<3;imethod++){
		hist_temp[imethod] = (TH1D*) h_t_reco_sartre[coh][imethod]->Clone(Form("hist_temp_%d",imethod));
		hist_temp[imethod]->SetLineStyle(1+imethod);
		hist_temp[imethod]->SetLineWidth(2);
	}
	hist_temp[0]->SetLineColor(kBlack);
	hist_temp[1]->SetLineColor(kBlue);
	hist_temp[2]->SetLineColor(kRed);
	hist_temp[0]->DrawNormalized("hist same");
	hist_temp[1]->DrawNormalized("hist same");
	hist_temp[2]->DrawNormalized("hist same");

	TLegend *w7 = new TLegend(0.53,0.74,0.73,0.87);
	w7->SetLineColor(kWhite);
	w7->SetFillColor(0);
	w7->SetTextSize(18);
	w7->SetTextFont(45);
	w7->AddEntry(hist_temp[0], "Sartre "+legendName+" Method E  ", "PL");
	w7->AddEntry(hist_temp[1], "Sartre "+legendName+" Method A  ", "PL");
	w7->AddEntry(hist_temp[2], "Sartre "+legendName+" Method L ", "PL");
	w7->Draw("same");
	if(coh==0){
		c3->Print("../figures/methods/sartre_"+name+"_coherent.pdf");
	}else{
		c3->Print("../figures/methods/sartre_"+name+"_incoherent.pdf");
	}
	


}