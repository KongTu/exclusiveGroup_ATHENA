#include "utility.h"
void plotVM_mass_and_t(TString name="phi", int coh = 1, int method=0){

	setVM(name);
	TFile* file_beagle = new TFile("../rootfiles/beagle_allVMs_w_breakups.root");
	TH1D* t_hat_all = (TH1D*) file_beagle->Get("h_trueT");

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
    TH1D* h_t_reco_sartre[2][3];
    TH2D* h_VM_t_mass_sartre[2][3];
    for(int ibreak=0;ibreak<2;ibreak++){
        for(int imethod=0;imethod<3;imethod++){
            h_t_reco_sartre[ibreak][imethod] = (TH1D*) file_sartre->Get(Form("h_t_reco_%d_%d",ibreak,imethod));
			h_VM_t_mass_sartre[ibreak][imethod] = (TH2D*) file_sartre->Get(Form("h_VM_t_mass_%d_%d",ibreak,imethod));
        }
    }

    TCanvas* c1 = new TCanvas("c1","c1",1,1,1000,1000);
    c1->Divide(10,10);
    TH1D* h_sartre_mass[100];
    TH1D* h_beagle_mass[100];
    for(int it=0;it<100;it++){
    	h_sartre_mass[it] = (TH1D*) h_VM_t_mass_sartre[coh][method]->ProjectionX(Form("h_sartre_mass_%d",it),it+1,it+1);
    	h_sartre_mass[it]->SetLineColor(kRed);
    	h_beagle_mass[it] = (TH1D*) h_VM_t_mass[coh][vm_index][method]->ProjectionX(Form("h_beagle_mass_%d",it),it+1,it+1);
    	h_beagle_mass[it]->SetLineColor(kBlue);
    	c1->cd(it+1);
    	gPad->SetLogy(1);
    	h_sartre_mass[it]->DrawNormalized();
    	h_beagle_mass[it]->DrawNormalized("same");
    }
    TCanvas* c1_1 = new TCanvas("c1_1","c1_1",1,1,600,600);
    gPad->SetLogy(1);
    h_sartre_mass[0]->GetXaxis()->SetTitle("Mass (GeV/c^{2})");
    h_sartre_mass[0]->GetYaxis()->SetTitle("Normalized Counts");
    h_sartre_mass[0]->SetTitle("first t bin");

    h_sartre_mass[0]->DrawNormalized();
	h_beagle_mass[0]->DrawNormalized("same");




}