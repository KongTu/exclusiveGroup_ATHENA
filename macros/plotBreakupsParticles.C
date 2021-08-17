#include "utility.h"
void plotBreakupsParticles(TString name="phi", bool veto_ = false, bool PHP_ = false){

	setVM(name);

	/* Beagle */
	
	if( (name=="rho" && PHP_)||
	 	(name=="phi" && PHP_)||
	 	(name=="jpsi" && PHP_)) {
	 	
	 	cout << "inconsistent settings between beagle and sartre! "<< endl;
	 	return;

	 }

	TString inputROOT="../rootfiles/beagle_allVMs_w_breakups.root";
	if(PHP_) inputROOT="../rootfiles/beagle_allVMs_w_breakups_PHP.root";
	if(veto_&&!PHP_) inputROOT="../rootfiles/beagle_allVMs_w_breakups_w_vetos.root";
	if(veto_&&PHP_) inputROOT="../rootfiles/beagle_allVMs_w_breakups_w_vetos_PHP.root";
	TFile* file_beagle = new TFile(inputROOT);

	TH2D* h_part[2][3][8];
	// first index, process 91 or 93
	// second index, vm particles
	// thrid, different species, proton, neutron, gamma, pi, kaon, e, mu, A*
	for(int ibreak=0;ibreak<2;ibreak++){
		for(int ivm=0;ivm<3;ivm++){
			for(int ipid=0;ipid<8;ipid++){
				h_part[ibreak][ivm][ipid] = (TH2D*) file_beagle->Get(Form("h_part_%d_%d_%d",ibreak,ivm,ipid));
			}
		}
	}

	TCanvas* c1 = new TCanvas("c1","c1",1,1,600,600);
	gPad->SetLogz(1);
	gPad->SetLeftMargin(0.12);
	gPad->SetBottomMargin(0.12);
	gPad->SetTopMargin(0.12);
	gPad->SetRightMargin(0.12);
	TGaxis::SetMaxDigits(2);
	TH1D* base1 = makeHist("base1", "", "Log(p) (GeV)", "Log(#theta) (mrad) ", 100,-6,11,kBlack);
	base1->GetYaxis()->SetRangeUser(-7, 5);
	base1->GetXaxis()->SetTitleColor(kBlack);
	fixedFontHist1D(base1,1.,1.1);
	base1->GetYaxis()->SetTitleSize(base1->GetYaxis()->GetTitleSize()*1.3);
	base1->GetXaxis()->SetTitleSize(base1->GetXaxis()->GetTitleSize()*1.3);
	base1->GetYaxis()->SetLabelSize(base1->GetYaxis()->GetLabelSize()*1.5);
	base1->GetXaxis()->SetLabelSize(base1->GetXaxis()->GetLabelSize()*1.5);
	base1->GetXaxis()->SetNdivisions(5,5,0);
	base1->GetYaxis()->SetNdivisions(6,6,0);
	base1->GetZaxis()->SetNdivisions(2,2,0);
	
	TGaxis *newaxis2 = new TGaxis(-6,
	                            5,
	                            11,
	                            5,
	                            0.0024787622,
	                            59873.699,
	                            510,"G-");
	newaxis2->SetLabelOffset(0.02);
	newaxis2->SetLabelFont(42);
	newaxis2->SetLabelSize( 0.04 );
	newaxis2->SetTitle("p (GeV)");
	newaxis2->SetTitleOffset(1.3);
	newaxis2->CenterTitle();

	base1->Draw("SS");
	newaxis2->Draw("SS");
	gPad->Update();

	TGaxis *newaxis3 = new TGaxis(11,
	                            -7,
	                            11,
	                            5,
	                            0.00091188626,
	                            148.41266,
	                            510,"G+");
	newaxis3->SetLabelOffset(0.06);
	newaxis3->SetLabelFont(42);
	newaxis3->SetLabelSize( 0.04 );
	newaxis3->SetTitle("#theta (mrad)");
	newaxis3->SetTitleOffset(1.5);
	newaxis3->CenterTitle();
	newaxis3->Draw("SS");
	gPad->Update();

	for(int ipid=0;ipid<8;ipid++){
		h_part[0][vm_index][ipid]->Add(h_part[1][vm_index][ipid],+1);
		h_part[0][vm_index][ipid]->Draw("cont1 same");
		if(veto_)c1->Print(Form("../figures/breakup_particles/veto_"+name+"_breakup_stagged_%d.pdf",ipid) );
		else c1->Print(Form("../figures/breakup_particles/"+name+"_breakup_stagged_%d.pdf",ipid) );
		
	}

}