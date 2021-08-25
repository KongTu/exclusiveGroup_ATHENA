#include "utility.h"

double photon_flux_from_average = 0.1;
void upcXsections(TH1D* hist){
	double binwidth=hist->GetBinWidth(1);
	double BR = 1.0;//ee & mumu channels.
	double beagle_lumi = 2e7/(1.971E3*197);//nanobarn
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
	TH1D* base1 = makeHist("base1", "", "|#it{t} | (GeV^{2})", "d#sigma/d|#it{t} | (#mub/GeV^{2}) ", 100,0,0.18,kBlack);
	base1->GetYaxis()->SetRangeUser(30, 1e5);
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
	upcXsections(h_VM_t[0][0][vm_]);
	upcXsections(h_VM_t[1][0][vm_]);

	h_VM_t[0][0][vm_]->Draw("hist same");
	h_VM_t[1][0][vm_]->SetLineColor(kRed);
	h_VM_t[1][0][vm_]->Draw("hist same");



}