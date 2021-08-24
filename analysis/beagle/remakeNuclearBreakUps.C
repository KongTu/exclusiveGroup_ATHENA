#include "../include/pleaseIncludeMe.h"
void remakeNuclearBreakUps(const TString filename="eA_TEST", const int nEvents = 40000, bool PHP_ = false){

	TChain *tree = new TChain("EICTree");
	tree->Add( filename+".root" );
	
	EventBeagle* event(NULL);
	tree->SetBranchAddress("event", &event);

	TFile* output = 0;
	TString outputROOT="../../rootfiles/remakeNuclearBreakUps.root";
	if(PHP_) outputROOT="../../rootfiles/remakeNuclearBreakUps_PHP.root";
	output = new TFile(outputROOT,"RECREATE");
	
	TH1D* h_VM_t[2][2][3];
	for(int iveto=0;iveto<2;iveto++){
		for(int iprocess=0;iprocess<2;iprocess++){
			for(int ivm=0;ivm<3;ivm++){
				h_VM_t[iveto][iprocess][ivm] = new TH1D(Form("h_VM_t_%d_%d_%d",iveto,iprocess,ivm),";-t (GeV^{2})",100,0.,0.2);
			}
		}
	}
	for(int i(0); i < nEvents; ++i ) {
      
		// Read the next entry from the tree.
		tree->GetEntry(i);
		double eFrac = 100.*static_cast<double>(i)/nEvents;
        if (i%100000 == 0) cout << "Processing event " << i << " (" << eFrac << "%)" << endl;
		
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

		//event cuts
		int processindex=-1;
		if( event_process==91) processindex=0;
		else if( event_process==93) processindex=1;
		else continue;

		if(PHP_){
			if( trueQ2 > 0.2 ) continue;
		}else{
			if( trueQ2 < 1. || trueQ2 > 20. ) continue;
		}
		if( trueW2<TMath::Power(1.95772,2)||trueW2>TMath::Power(88.9985,2)) continue;//to match Sartre
		/*some conditions & initialization*/
		int pdglist[]={113,333,443};
		int statuslist[]={2,2,2};
		int hasvm[3]={0,0,0};
		TLorentzVector vm_vect[3];
		for(int ivm=0;ivm<3;ivm++){vm_vect[ivm].SetPxPyPzE(0.,0.,0.,0.);}
		//begin particle loop;
		for(int j(0); j < nParticles; ++j ) {
			const erhic::ParticleMC* particle = event->GetTrack(j);
			int status = particle->GetStatus();
			int pdg = particle->GetPdgCode();		
			for(int ivm=0;ivm<3;ivm++){
				if(pdg!=pdglist[ivm]) continue;
				if(status!=statuslist[ivm]) continue;
				vm_vect[ivm]=particle->Get4Vector();
				if( fabs(vm_vect[ivm].Rapidity()) < 1.0 ) {
					hasvm[ivm]=1;//found vm.
				}
			}
		}
		//after particle loop;
		for(int ivm=0;ivm<3;ivm++){
			if(hasvm[ivm]){
				h_VM_t[0][processindex][ivm]->Fill( -t_hat );
				//perform veto.
				if( veto_this_event(event, nParticles) ) continue;
				h_VM_t[1][processindex][ivm]->Fill( -t_hat );
			}
		}

	}//end event loop

	output->Write();
	output->Close();


}
