#include "RiceStyle.h"
using namespace std;
#define PI 3.1415926
#define MASS_PROTON   0.93827208816
#define MASS_NEUTRON  0.93956542052
#define MASS_NUCLEON  0.93891875 //.93891875
#define MASS_DEUTERON  1.8756129 //1.8756129 //1.8756134 (most precise)

void plotdSigmadt(){

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

	double beagle_lumi = (43242./(3.4E+1*82)) + (65815./(3.44E+1*(208-82)));//nanobarn
	double beagle_delta_t = h_VM[0][1][4]->GetBinWidth(1);
	double BR_phiTokk = 0.489;//branching ratio

	/* Sartre */

	TFile* file_sartre = new TFile("../rootfiles/sartre_phi.root");
	TH1D* h_phi_coh_sartre = (TH1D*) file_sartre->Get("hist_t_coherent");
	TH1D* h_phi_incoh_sartre = (TH1D*) file_sartre->Get("hist_t_incoherent");

	double sartre_lumi = (20000000./(4.72E+3));//nanbarn
	double sartre_delta_t = h_phi_coh_sartre->GetBinWidth(1); 

	h_phi_coh_sartre->Scale(1./(sartre_lumi * sartre_delta_t * BR_phiTokk));
	h_phi_coh_sartre->SetLineColor(kBlue);
	h_phi_coh_sartre->Draw("HIST ");

	h_phi_incoh_sartre->Scale(1./(sartre_lumi * sartre_delta_t * BR_phiTokk));
	h_phi_incoh_sartre->SetLineColor(kRed);
	h_phi_incoh_sartre->Draw("HIST same");

	TH1D* h_VM_background = (TH1D*) h_VM[0][1][4]->Clone("h_VM_background");
	h_VM_background->Add(h_VM[1][1][4],+1);
	h_VM_background->SetMarkerStyle(24);
	h_VM_background->Scale( 1./ (beagle_lumi * beagle_delta_t * BR_phiTokk) );
	h_VM_background->Draw("P SAME");

}