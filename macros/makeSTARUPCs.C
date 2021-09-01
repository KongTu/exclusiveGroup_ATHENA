#include "utility.h"

double photon_flux_from_average = 0.161227;// this is cut on W [8,47] GeV
// double photon_flux_from_average = 0.15124;//this is cut -1<y_J<4

void upcXsections(TH1D* hist){
	double binwidth=hist->GetBinWidth(1);
	double BR = 1.0;//ee & mumu channels.
	double beagle_lumi = 4e7/(1.971E3*197);//nanobarn
	double flux = photon_flux_from_average;
	double flux_UPCs = 11.78;
	double factor = (flux * binwidth * BR * beagle_lumi);
	hist->Scale(1./factor);
	hist->Scale(flux_UPCs); // star photon flux.
	hist->Scale(2.); //star rapidity coverage.
	hist->Scale(0.001);//to microbarn;
}

void makeSTARUPCs(TString name="jpsi"){

	TString inputROOT="../rootfiles/remakeNuclearBreakUps_PHP.root";
	TFile* file_beagle = new TFile(inputROOT);
	//histograms
	TH1D* h_VM_t[2][2][3];
	for(int iveto=0;iveto<2;iveto++){
		for(int iprocess=0;iprocess<2;iprocess++){
			for(int ivm=0;ivm<3;ivm++){
				h_VM_t[iveto][iprocess][ivm] = (TH1D*) file_beagle->Get(Form("h_VM_t_%d_%d_%d",iveto,iprocess,ivm));
			}
		}
	}
	TH1D* h_photon[3];
	TH1D* h_w[3];
	for(int ivm=0;ivm<3;ivm++){
		h_photon[ivm] = (TH1D*) file_beagle->Get(Form("h_photon_%d",ivm));
		h_w[ivm] = (TH1D*) file_beagle->Get(Form("h_w_%d",ivm));
	}


	int vm_ = 2;
	if(name=="phi") vm_ = 1; 
	if(name=="rho") vm_ = 0;

	TCanvas* c1 = new TCanvas("c1","c1",1,1,600,600);
	gPad->SetLogy(1);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	TH1D* base1 = makeHist("base1", "", "|#it{t} | (GeV^{2})", "d#sigma/d|#it{t} | (#mub/GeV^{2}) ", 100,0,0.5,kBlack);
	base1->GetYaxis()->SetRangeUser(1e-1, 1e5);
	base1->GetXaxis()->SetTitleColor(kBlack);
	fixedFontHist1D(base1,1.,1.1);
	base1->GetYaxis()->SetTitleSize(base1->GetYaxis()->GetTitleSize()*1.5);
	base1->GetXaxis()->SetTitleSize(base1->GetXaxis()->GetTitleSize()*1.5);
	base1->GetYaxis()->SetLabelSize(base1->GetYaxis()->GetLabelSize()*1.5);
	base1->GetXaxis()->SetLabelSize(base1->GetXaxis()->GetLabelSize()*1.5);
	base1->GetXaxis()->SetNdivisions(4,4,0);
	base1->GetYaxis()->SetNdivisions(5,5,0);
	base1->Draw();

	photon_flux_from_average = h_photon[vm_]->GetMean();

	h_VM_t[0][0][vm_]->Add(h_VM_t[0][1][vm_],+1);
	h_VM_t[1][0][vm_]->Add(h_VM_t[1][1][vm_],+1);
	TH1D* h_raw_all = (TH1D*) h_VM_t[0][0][vm_]->Clone("h_raw_all");
	TH1D* h_raw_neutron = (TH1D*) h_VM_t[1][0][vm_]->Clone("h_raw_neutron");
	
	upcXsections(h_VM_t[0][0][vm_]);
	upcXsections(h_VM_t[1][0][vm_]);

	h_VM_t[0][0][vm_]->Draw("hist same");
	h_VM_t[1][0][vm_]->SetLineColor(kRed);
	h_VM_t[1][0][vm_]->Draw("hist same");

	TLatex* r42 = new TLatex(0.18, 0.91, "eAu 18x110 GeV^{2}");
	r42->SetNDC();
	r42->SetTextSize(22);
	r42->SetTextFont(43);
	r42->SetTextColor(kBlack);
	// r42->Draw("same");

	TLatex* r43 = new TLatex(0.7,0.91, "BeAGLE");
	r43->SetNDC();
	r43->SetTextSize(0.04);
	r43->Draw("same");

	TLatex* r44_1 = new TLatex(0.18, 0.84, "Q^{2}<0.2 GeV^{2}, |y_{VM}|<1.0");
	r44_1->SetNDC();
	r44_1->SetTextSize(20);
	r44_1->SetTextFont(43);
	r44_1->SetTextColor(kBlack);
	r44_1->Draw("same");

	TLatex* r44_2 = new TLatex(0.18, 0.79, "Prediction for STAR AuAu J/#psi UPCs");
	r44_2->SetNDC();
	r44_2->SetTextSize(20);
	r44_2->SetTextFont(43);
	r44_2->SetTextColor(kBlack);
	r44_2->Draw("same");

	TLegend *w5 = new TLegend(0.2,0.65,0.53,0.74);
	w5->SetLineColor(kWhite);
	w5->SetFillColor(0);
	w5->SetTextSize(18);
	w5->SetTextFont(45);
	w5->AddEntry(h_VM_t[0][0][vm_], "Incoherent all ", "PL");
	w5->AddEntry(h_VM_t[1][0][vm_], "Incoherent ZDC neutron vetoed  ", "PL");
	w5->Draw("same");

	TFile * output = new TFile("../UPCs/forSTAR_AuAu200_Jpsi.root","RECREATE");
	h_VM_t[0][0][vm_]->SetName("h_correct_all");
	h_VM_t[0][0][vm_]->GetYaxis()->SetTitle("d#sigma/d|#it{t} | (#mub/GeV^{2})");
	h_VM_t[0][0][vm_]->GetXaxis()->SetTitle("|#it{t} | (GeV^{2})");
	h_VM_t[1][0][vm_]->SetName("h_correct_neutron");
	h_VM_t[1][0][vm_]->GetYaxis()->SetTitle("d#sigma/d|#it{t} | (#mub/GeV^{2})");
	h_VM_t[1][0][vm_]->GetXaxis()->SetTitle("|#it{t} | (GeV^{2})");
	h_VM_t[0][0][vm_]->Write();
	h_VM_t[1][0][vm_]->Write();
	h_raw_all->GetYaxis()->SetTitle("Raw counts");
	h_raw_neutron->GetYaxis()->SetTitle("Raw counts");
	h_raw_all->GetXaxis()->SetTitle("|#it{t} | (GeV^{2})");
	h_raw_neutron->GetXaxis()->SetTitle("|#it{t} | (GeV^{2})");
	h_raw_all->Write();
	h_raw_neutron->Write();
	c1->Print("../UPCs/forSTAR_AuAu200_Jpsi.pdf");

}