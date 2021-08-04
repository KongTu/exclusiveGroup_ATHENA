#include "RiceStyle.h"
using namespace std;
#define PI 3.1415926
#define MASS_PROTON   0.93827208816
#define MASS_NEUTRON  0.93956542052
#define MASS_NUCLEON  0.93891875 //.93891875
#define MASS_DEUTERON  1.8756129 //1.8756129 //1.8756134 (most precise)

void plotdSigmadt(){

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

	TFile* file_sartre = new TFile("../rootfiles/sartre_phi.root");
	TH1D* h_phi_coh_sartre = (TH1D*) file_sartre->Get("hist_t_coherent");
	TH1D* h_phi_incoh_sartre = (TH1D*) file_sartre->Get("hist_t_incoherent");



	// h_VM[0][1][4]->DrawNormalized();
	h_phi_incoh_sartre->DrawNormalized("");


}