#include "RiceStyle.h"
#include "TGaxis.h"
using namespace std;
#define PI 3.1415926
#define MASS_PROTON   0.93827208816
#define MASS_NEUTRON  0.93956542052
#define MASS_NUCLEON  0.93891875 //.93891875
#define MASS_DEUTERON  1.8756129 //1.8756129 //1.8756134 (most precise)

void plotVMdaughters(TString name="phi"){

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

	TCanvas* c1 = new TCanvas("c1","c1",1,1,600,600);
	gPad->SetLogy(0);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	TGaxis::SetMaxDigits(2);
	TH1D* base1 = makeHist("base1", "", "#eta_{lab}", "dN/d#eta ", 100,-8,8,kBlack);
	int vm_index=-1;
	if(name=="rho"){base1->GetYaxis()->SetRangeUser(0, 6e4);vm_index=0;}
	if(name=="phi"){base1->GetYaxis()->SetRangeUser(0, 6e3);vm_index=1;}
	if(name=="jpsi"){base1->GetYaxis()->SetRangeUser(0, 6e2);vm_index=2;}//jpsi nodecay.
	
	base1->GetXaxis()->SetTitleColor(kBlack);
	fixedFontHist1D(base1,1.,1.1);
	base1->GetYaxis()->SetTitleSize(base1->GetYaxis()->GetTitleSize()*1.5);
	base1->GetXaxis()->SetTitleSize(base1->GetXaxis()->GetTitleSize()*1.5);
	base1->GetYaxis()->SetLabelSize(base1->GetYaxis()->GetLabelSize()*1.5);
	base1->GetXaxis()->SetLabelSize(base1->GetXaxis()->GetLabelSize()*1.5);
	base1->GetXaxis()->SetNdivisions(5,5,0);
	base1->GetYaxis()->SetNdivisions(3,2,0);
	base1->Draw();

	h_VM_daughter[0][vm_index][1]->SetFillColorAlpha(kBlue,0.4);
    h_VM_daughter[0][vm_index][1]->SetFillStyle(1001);
	h_VM_daughter[0][vm_index][1]->SetMarkerStyle(24);
	h_VM_daughter[0][vm_index][1]->SetMarkerColor(kBlue-1);

	h_VM_daughter[1][vm_index][1]->SetFillColorAlpha(kRed,0.4);
    h_VM_daughter[1][vm_index][1]->SetFillStyle(1001);
	h_VM_daughter[1][vm_index][1]->SetMarkerStyle(25);
	h_VM_daughter[1][vm_index][1]->SetMarkerColor(kRed);

	h_VM_daughter[0][vm_index][1]->Draw("P e3 same");
	h_VM_daughter[1][vm_index][1]->Draw("P e3 same");

	TString legendName="";
	if(name=="rho"){
		vm_index=0;
		legendName="#rho^{0}";
	}
	else if(name=="phi"){
		vm_index=1;
		legendName="#phi";
	}
	else if(name=="jpsi"){
		vm_index=2;
		legendName="J/#psi";
	}
	TLegend *w5 = new TLegend(0.5,0.74,0.63,0.87);
	w5->SetLineColor(kWhite);
	w5->SetFillColor(0);
	w5->SetTextSize(18);
	w5->SetTextFont(45);
	// w5->AddEntry(h_phi_coh_sartre, "Sartre "+legendName+" coherent  ", "PL");
	// w5->AddEntry(h_phi_incoh_sartre, "Sartre "+legendName+" incoherent  ", "PL");
	w5->AddEntry(h_VM_daughter[0][vm_index][1], "BeAGLE elas. "+legendName+" daug. incoh. ", "P");
	w5->AddEntry(h_VM_daughter[1][vm_index][1], "BeAGLE diss. "+legendName+" daug. incoh. ", "P");
	w5->Draw("same");

	TLatex* r42 = new TLatex(0.18, 0.84, "eAu 18x110 GeV^{2}");
	r42->SetNDC();
	r42->SetTextSize(22);
	r42->SetTextFont(43);
	r42->SetTextColor(kBlack);
	r42->Draw("same");

	TLatex* r43 = new TLatex(0.6,0.91, "ATHENA #it{Internal}");
	r43->SetNDC();
	r43->SetTextSize(0.04);
	r43->Draw("same");


	TCanvas* c2 = new TCanvas("c2","c2",1,1,600,600);
	gPad->SetLogy(0);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	TGaxis::SetMaxDigits(2);
	TH1D* base2 = makeHist("base2", "", "p_{T} (GeV)", "dN/dp_{T} (GeV^{-1}) ", 100,0,3,kBlack);
	if(name=="rho"){base2->GetYaxis()->SetRangeUser(0, 8e4);vm_index=0;}
	if(name=="phi"){base2->GetYaxis()->SetRangeUser(0, 8e3);vm_index=1;}
	if(name=="jpsi"){base2->GetYaxis()->SetRangeUser(0, 6e2);vm_index=2;}//jpsi nodecay.
	base2->GetXaxis()->SetTitleColor(kBlack);
	fixedFontHist1D(base2,1.,1.1);
	base2->GetYaxis()->SetTitleSize(base2->GetYaxis()->GetTitleSize()*1.5);
	base2->GetXaxis()->SetTitleSize(base2->GetXaxis()->GetTitleSize()*1.5);
	base2->GetYaxis()->SetLabelSize(base2->GetYaxis()->GetLabelSize()*1.5);
	base2->GetXaxis()->SetLabelSize(base2->GetXaxis()->GetLabelSize()*1.5);
	base2->GetXaxis()->SetNdivisions(5,5,0);
	base2->GetYaxis()->SetNdivisions(3,2,0);
	base2->Draw();
	h_VM_daughter[0][vm_index][0]->SetFillColorAlpha(kBlue,0.4);
    h_VM_daughter[0][vm_index][0]->SetFillStyle(1001);
	h_VM_daughter[0][vm_index][0]->SetMarkerStyle(24);
	h_VM_daughter[0][vm_index][0]->SetMarkerColor(kBlue-1);

	h_VM_daughter[1][vm_index][0]->SetFillColorAlpha(kRed,0.4);
    h_VM_daughter[1][vm_index][0]->SetFillStyle(1001);
	h_VM_daughter[1][vm_index][0]->SetMarkerStyle(25);
	h_VM_daughter[1][vm_index][0]->SetMarkerColor(kRed);

	h_VM_daughter[0][vm_index][0]->Draw("P e3 same");
	h_VM_daughter[1][vm_index][0]->Draw("P e3 same");


	w5->Draw("same");
	r42->Draw("same");
	r43->Draw("same");

}