#include "../include/pleaseIncludeMe.h"
void testBreakUps(const TString filename="eA_TEST", const int nEvents = 40000, bool PHP_ = false, bool veto_ = true){

	TChain *tree = new TChain("EICTree");
	tree->Add( filename+".root" );
	
	EventBeagle* event(NULL);
	tree->SetBranchAddress("event", &event);

	TFile* output = 0;
	TString outputROOT="../../rootfiles/testBreakUps.root";
	if(PHP_) outputROOT="../../rootfiles/testBreakUps_PHP.root";
	if(veto_&&!PHP_) outputROOT="../../rootfiles/testBreakUps_w_vetos.root";
	if(veto_&&PHP_) outputROOT="../../rootfiles/testBreakUps_w_vetos_PHP.root";
	output = new TFile(outputROOT,"RECREATE");
	
	TH1D* hist_multiplicity = new TH1D("hist_multiplicity",";N",50,-0.5,49.5);
	TH1D* h_trueT = new TH1D("h_trueT",";-t (GeV^{2})", 100,0,0.2);
	TH1D* h_trueT_91 = new TH1D("h_trueT_91",";-t (GeV^{2})", 100,0,0.2);
	TH1D* h_trueT_91_after = new TH1D("h_trueT_91_after",";-t (GeV^{2})", 100,0,0.2);
	TH1D* h_trueT_93 = new TH1D("h_trueT_93",";-t (GeV^{2})", 100,0,0.2);
	TH1D* h_trueT_93_after = new TH1D("h_trueT_93_after",";-t (GeV^{2})", 100,0,0.2);
	//Save each step:
	TH1D* h_veto_step_91[6];
	TH1D* h_veto_step_93[6];
		for(int ihist=0;ihist<6;ihist++){
			h_veto_step_91[ihist] = new TH1D(Form("h_veto_step_91_%d",ihist),";-t (GeV^{2})", 100,0,0.2);
			h_veto_step_93[ihist] = new TH1D(Form("h_veto_step_93_%d",ihist),";-t (GeV^{2})", 100,0,0.2);
		}
	//
	for(int i(0); i < nEvents; ++i ) {
      
		// Read the next entry from the tree.
		tree->GetEntry(i);
		double eFrac = 100.*static_cast<double>(i)/nEvents;
        if (i%100000 == 0) cout << "Processing event " << i << " (" << eFrac << "%)" << endl;
		
		double pzlep = event->pzlep;
		double pztarg = event->pztarg;
		int Atarg = event->Atarg;
		double pztargA=pztarg*Atarg;
		int struck_nucleon = event->nucleon;
		double MASS_NUCLEON = MASS_PROTON;

		//event information:
		double trueQ2 = event->GetTrueQ2();
		double trueW2 = event->GetTrueW2();
		double trueX = event->GetTrueX();
		double trueY = event->GetTrueY();
		double trueNu = event->nu;
		double s_hat = event->GetHardS();
		double t_hat = event->t_hat;
		double u_hat = event->GetHardU();
		double photon_flux = event->GetPhotonFlux();
		int event_process = event->GetProcess();
		int nParticles = event->GetNTracks();
		double pxf = event->pxf;
		double pyf = event->pyf;
		double pzf = event->pzf;
		double pf3 = sqrt(pxf*pxf+pyf*pyf+pzf*pzf);

		double impact_parameter = event->b;
		double Tb = event->Thickness;
		double distance = event->d1st;
		int N_nevap = event->Nnevap;
		int N_pevap = event->Npevap;

		//event cuts
		int processindex=-1;
		if( event_process==91) processindex=0;
		else if( event_process==93) processindex=1;
		else processindex=2;

		if(PHP_){
			if( trueQ2 > 0.2 || trueQ2 < 0.05 ) continue;
		}else{
			if( trueQ2 < 1. || trueQ2 > 20. ) continue;
		}
		// if( trueY > 0.95 || trueY < 0.01 ) continue;
		if( trueW2<TMath::Power(1.95772,2)||trueW2>TMath::Power(88.9985,2)) continue;//to match Sartre
		//do analysis, or fill historgrams for event levels
		h_trueT->Fill(-t_hat);
		if(processindex==0) h_trueT_91->Fill(-t_hat);
		if(processindex==1) h_trueT_93->Fill(-t_hat);

		int multiplicity=0;
		for(int j(0); j < nParticles; ++j ) {
			const erhic::ParticleMC* particle = event->GetTrack(j);
			int status = particle->GetStatus();
			double pt = particle->GetPt();
			double eta = particle->GetEta();
			int charge = particle->eA->charge;
			if( status!= 1 ) continue;
			if(TMath::Abs(eta)<4.0 && pt>0.15 && charge!=0) multiplicity++;
		}
		hist_multiplicity->Fill(multiplicity);

		//veto by step
		for(int istep=0;istep<6;istep++){
			if( !veto_this_event(event, nParticles,istep)&&processindex==0 ) h_veto_step_91[istep]->Fill(-t_hat);
			if( !veto_this_event(event, nParticles,istep)&&processindex==1 ) h_veto_step_93[istep]->Fill(-t_hat);
		}
		//perform veto.
		if( veto_ ){
			if( veto_this_event(event, nParticles) ) continue;
		}
		if(processindex==0) h_trueT_91_after->Fill(-t_hat);
		if(processindex==1) h_trueT_93_after->Fill(-t_hat);
		

	}//end event loop

	output->Write();
	output->Close();


}
