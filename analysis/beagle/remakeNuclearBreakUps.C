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
				h_VM_t[iveto][iprocess][ivm] = new TH1D(Form("h_VM_t_%d_%d_%d",iveto,iprocess,ivm),";-t (GeV^{2})",100,0.,2);
			}
		}
	}
	TH1D* h_photon[3];
	TH1D* h_w[3];
	TH1D* h_vm_y[3];
	TH1D* h_vm_y_weight[3];
	TH2D* h_neutron[3];
	TH2D* h_w_y[3];
	for(int ivm=0;ivm<3;ivm++){
		h_photon[ivm] = new TH1D(Form("h_photon_%d",ivm),";flux",100,0,0.2);
		h_w[ivm] = new TH1D(Form("h_w_%d",ivm),";W (GeV)",100,1.5,90);
		h_vm_y[ivm] = new TH1D(Form("h_vm_y_%d",ivm),";rapidity",100,-3,5);
		h_vm_y_weight[ivm] = new TH1D(Form("h_vm_y_weight_%d",ivm),";rapidity",100,-3,5);
		h_neutron[ivm] = new TH2D(Form("h_neutron_%d",ivm),Form("h_neutron_%d",ivm),100,0,5,100,-PI,PI);
		h_w_y[ivm] = new TH2D(Form("h_w_y%d",ivm),Form("h_w_y%d",ivm),100,-3,5,100,1.5,90);
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
			if( trueQ2 > 1.0 ) continue;
		}else{
			if( trueQ2 < 1. || trueQ2 > 20. ) continue;
		}
		//this is the right way to cut on the phase space, NOT thur J/psi rapidity.
		if( trueW2<TMath::Power(8,2)||trueW2>TMath::Power(47,2)) continue;
		/*some conditions & initialization*/
		int pdglist[]={113,333,443};
		int statuslist[]={2,2,2};
		int hasvm[3]={0,0,0};
		TLorentzVector vm_vect[3];
		vector< TLorentzVector> neutron_vect_443;
		for(int ivm=0;ivm<3;ivm++){vm_vect[ivm].SetPxPyPzE(0.,0.,0.,0.);}
		//begin particle loop;
		for(int j(0); j < nParticles; ++j ) {
			const erhic::ParticleMC* particle = event->GetTrack(j);
			int status = particle->GetStatus();
			int pdg = particle->GetPdgCode();
			if(status==1&&pdg==2112){neutron_vect_443.push_back(particle->Get4Vector());}		
			for(int ivm=0;ivm<3;ivm++){
				if(pdg!=pdglist[ivm]) continue;
				if(status!=statuslist[ivm]) continue;
				vm_vect[ivm]=particle->Get4Vector();
				hasvm[ivm]=1;//found vm.

			}
		}
		//after particle loop;
		for(int ivm=0;ivm<3;ivm++){
			if(hasvm[ivm]){
				h_VM_t[0][processindex][ivm]->Fill( -t_hat );
				h_photon[ivm]->Fill(photon_flux);
				h_w[ivm]->Fill(sqrt(trueW2));
				h_vm_y[ivm]->Fill( vm_vect[ivm].Rapidity() );
				h_vm_y_weight[ivm]->Fill( vm_vect[ivm].Rapidity(), 1./photon_flux );
				h_w_y[ivm]->Fill(vm_vect[ivm].Rapidity(), sqrt(trueW2));
				for(unsigned in=0;in<neutron_vect_443.size();in++){
					h_neutron[ivm]->Fill(neutron_vect_443[in].Theta()*1000, neutron_vect_443[in].Phi());//mrad
				}
				//perform veto.
				if( veto_this_event(event, nParticles, 6) ) continue;
				h_VM_t[1][processindex][ivm]->Fill( -t_hat );
			}
		}

	}//end event loop

	output->Write();
	output->Close();


}
