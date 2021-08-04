#include "RiceStyle.h"
using namespace std;
#define PI 3.1415926
#define MASS_PROTON   0.93827208816
#define MASS_NEUTRON  0.93956542052
#define MASS_NUCLEON  0.93891875 //.93891875
#define MASS_DEUTERON  1.8756129 //1.8756129 //1.8756134 (most precise)

void plotdSigmadt(TString name="phi"){

	/* Beagle */
	
	TFile* file_beagle = new TFile("../rootfiles/beagle_allVMs.root");
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

	double beagle_lumi = 59813./(34.4*118) + (100000-59813)/(34.4*79);//nanobarn
	double beagle_delta_t = h_VM[0][1][4]->GetBinWidth(1);
	double BR_phiTokk = 0.489;//branching ratio

	/* Sartre */

	TFile* file_sartre = new TFile("../rootfiles/sartre_"+name+"_bnonsat.root");
	TH1D* h_phi_coh_sartre = (TH1D*) file_sartre->Get("hist_t_coherent");
	TH1D* h_phi_incoh_sartre = (TH1D*) file_sartre->Get("hist_t_incoherent");

	double sigma=-1;
	int vm_index=-1;
	if(name=="rho"){
		sigma = 3.58E+4;
		vm_index=0;
	}
	else if(name=="phi"){
		sigma=4.72E+3;
		vm_index=1;
	}
	double sartre_lumi = 20000000./sigma;//nanbarn
	double sartre_delta_t = h_phi_coh_sartre->GetBinWidth(1); 

	TCanvas* c1 = new TCanvas("c1","c1",1,1,600,600);
	gPad->SetLogy(1);
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	TH1D* base1 = makeHist("base1", "", "|#it{t} | (GeV^{2})", "d#sigma/d|#it{t} | (nb/GeV^{2}) ", 100,0,0.18,kBlack);
	base1->GetYaxis()->SetRangeUser(1e-1, 1e8);
	base1->GetXaxis()->SetTitleColor(kBlack);
	fixedFontHist1D(base1,1.,1.1);
	base1->GetYaxis()->SetTitleSize(base1->GetYaxis()->GetTitleSize()*1.7);
	base1->GetXaxis()->SetTitleSize(base1->GetXaxis()->GetTitleSize()*1.7);
	base1->GetYaxis()->SetLabelSize(base1->GetYaxis()->GetLabelSize()*1.8);
	base1->GetXaxis()->SetLabelSize(base1->GetXaxis()->GetLabelSize()*1.7);
	base1->GetXaxis()->SetNdivisions(4,4,0);
	base1->GetYaxis()->SetNdivisions(5,5,0);
	base1->Draw();

	h_phi_coh_sartre->Scale(1./(sartre_lumi * sartre_delta_t * BR_phiTokk));
	h_phi_coh_sartre->SetLineColor(kBlue);
	h_phi_coh_sartre->Draw("HIST same");

	h_phi_incoh_sartre->Scale(1./(sartre_lumi * sartre_delta_t * BR_phiTokk));
	h_phi_incoh_sartre->SetLineColor(kRed);
	h_phi_incoh_sartre->Draw("HIST same");

	TH1D* h_VM_background = (TH1D*) h_VM[0][vm_index][4]->Clone("h_VM_background");
	h_VM_background->Add(h_VM[1][vm_index][4],+1);
	h_VM_background->Rebin(2);
	h_VM_background->Scale(1./2);
	h_VM_background->SetMarkerStyle(24);
	h_VM_background->Scale( 1./ (beagle_lumi * beagle_delta_t * BR_phiTokk) );
	h_VM_background->Draw("P SAME");

	TLegend *w5 = new TLegend(0.56,0.74,0.78,0.87);
	w5->SetLineColor(kWhite);
	w5->SetFillColor(0);
	w5->SetTextSize(19);
	w5->SetTextFont(45);
	w5->AddEntry(h_phi_coh_sartre, "Sartre #phi coherent  ", "PL");
	w5->AddEntry(h_phi_incoh_sartre, "Sartre #phi incoherent  ", "PL");
	w5->AddEntry(h_VM_background, "BeAGLE #phi incoherent ", "P");
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

	TLatex* r44 = new TLatex(0.18, 0.84, "1<Q^{2}<20 GeV^{2}, |#eta(K)|<4.0");
	r44->SetNDC();
	r44->SetTextSize(20);
	r44->SetTextFont(43);
	r44->SetTextColor(kBlack);
	r44->Draw("same");

}