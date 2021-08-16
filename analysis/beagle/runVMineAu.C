#include "../include/pleaseIncludeMe.h"
void runVMineAu(const TString filename="eA_TEST", const int nEvents = 40000, bool veto_ = true){

	TChain *tree = new TChain("EICTree");
	tree->Add( filename+".root" );
	
	EventBeagle* event(NULL);
	tree->SetBranchAddress("event", &event);

	TFile* output = 0;
	TString outputROOT="../../rootfiles/beagle_allVMs_w_breakups.root";
	if(veto_) outputROOT="../../rootfiles/beagle_allVMs_w_breakups_w_vetos.root";
	output = new TFile(outputROOT,"RECREATE");
	
	TH1D* h_trueT = new TH1D("h_trueT",";-t (GeV^{2})", 100,0,0.5);
	//VM histograms//
	/* first   index VM process, 91=0, 93=1, everything else=2*/
	/* second  index VM species, rho=0, phi=1, jpsi=2*/
	/* third   index VM property, pt=0, eta=1, phi=2, theta=3, reserved=4*/
	double bin_lower[]={0.,-8.,0.,0.,0.};
	double bin_upper[]={5.0,8.,6.5,4.,0.2};
	TH1D* h_VM[3][3][5];
	for(int ibreak=0;ibreak<3;ibreak++){
		for(int ivm=0;ivm<3;ivm++){
			for(int ipro=0;ipro<5;ipro++){
				h_VM[ibreak][ivm][ipro] = new TH1D(Form("h_VM_%d_%d_%d",ibreak,ivm,ipro),
					Form("h_VM_%d_%d_%d",ibreak,ivm,ipro),100,bin_lower[ipro],bin_upper[ipro] );
			}
		}
	}
	//END VM histograms//
	// .
	// .
	//VM daughter histograms//
	TH1D* h_VM_daughter[3][3][5];
	for(int ibreak=0;ibreak<3;ibreak++){
		for(int ivm=0;ivm<3;ivm++){
			for(int ipro=0;ipro<5;ipro++){
				h_VM_daughter[ibreak][ivm][ipro] = new TH1D(Form("h_VM_daughter_%d_%d_%d",ibreak,ivm,ipro),
					Form("h_VM_daughter_%d_%d_%d",ibreak,ivm,ipro),100,bin_lower[ipro],bin_upper[ipro] );
			}
		}
	}
	//END VM daughter histograms
	// .
	// .
	//Nuclear Breakup histograms
	// first index, process 91 or 93
	// second index, vm particles
	// thrid, different species, proton, neutron, gamma, pi, kaon, e, mu, A*
	TH2D* h_part[3][3][8];
	for(int ibreak=0;ibreak<3;ibreak++){
		for(int ivm=0;ivm<3;ivm++){
			for(int ipid=0;ipid<8;ipid++){
				h_part[ibreak][ivm][ipid] = new TH2D(Form("h_part_%d_%d_%d",ibreak,ivm,ipid),
					Form("h_part_%d_%d_%d",ibreak,ivm,ipid),100,TMath::Log(1e-3),TMath::Log(3e4),100,TMath::Log(1e-3),TMath::Log(50));
				
			}
		}
	}
	//END Nuclear Breakup histograms
	// .
	// .
	// t_reco histograms
	TH1D* h_t_reco[3][3][3];
	for(int ibreak=0;ibreak<3;ibreak++){
		for(int ivm=0;ivm<3;ivm++){
			for(int imethod=0;imethod<3;imethod++){
				h_t_reco[ibreak][ivm][imethod] = new TH1D(Form("h_t_reco_%d_%d_%d",ibreak,ivm,imethod),
					Form("h_t_reco_%d_%d_%d",ibreak,ivm,imethod),1000,0,2 );
			}
		}
	}
	// vm mass from daughters
	TH1D* h_VM_mass[3][3];
	for(int ibreak=0;ibreak<3;ibreak++){
		for(int ivm=0;ivm<3;ivm++){
			h_VM_mass[ibreak][ivm] = new TH1D(Form("h_VM_mass_%d_%d",ibreak,ivm),
				Form("h_VM_mass_%d_%d",ibreak,ivm),1000,0.,4);
		}
	}
	//nuclear remnant mass
	TH2D* h_Amass[3][3];
	for(int ibreak=0;ibreak<3;ibreak++){
		for(int ivm=0;ivm<3;ivm++){
			h_Amass[ibreak][ivm] = new TH2D(Form("h_Amass_%d_%d",ibreak,ivm),
				Form("h_Amass_%d_%d",ibreak,ivm), 100,0,0.2,100,-3,3);
		}
	}

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
		if( struck_nucleon==2112 ) MASS_NUCLEON = MASS_NEUTRON;

		TLorentzVector e_beam(0.,0.,pzlep,sqrt(pzlep*pzlep+MASS_ELECTRON*MASS_ELECTRON));
		TLorentzVector p_beam(0.,0.,pztarg,sqrt(pztarg*pztarg+MASS_NUCLEON*MASS_NUCLEON));
		TLorentzVector A_beam(0.,0.,pztargA,sqrt(pztargA*pztargA+MASS_AU197*MASS_AU197));
		TLorentzVector e_scattered(0.,0.,0.,0.);
		TLorentzVector vm_vect[3];
		for(int ivm=0;ivm<3;ivm++){
			vm_vect[ivm].SetPxPyPzE(0.,0.,0.,0.);
		}

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
		
		double impact_parameter = event->b;
		double Tb = event->Thickness;
		double distance = event->d1st;
		int N_nevap = event->Nnevap;
		int N_pevap = event->Npevap;

		//do analysis, or fill historgrams for event levels
		h_trueT->Fill(-t_hat);

		//event cuts
		int processindex=-1;
		if( event_process==91) processindex=0;
		else if( event_process==93) processindex=1;
		else processindex=2;
		if( trueQ2 < 1. || trueQ2 > 20. ) continue;
		// if( trueY > 0.95 || trueY < 0.01 ) continue;
		if( trueW2<TMath::Power(1.95772,2)||trueW2>TMath::Power(88.9985,2)) continue;//to match Sartre
		//perform veto.
		if( veto_ ){
			if( veto_this_event(event, nParticles) ) continue;
		}
		
		//particle loop
		//rho^0 = 113, decay->pipi
		//phi = 333, decay->kk
		//jpsi = 443, nodecay
		int pdglist[]={113,333,443};
		int statuslist[]={2,2,2};
		int acceptance[3]={1,1,1};
		int hasvm[3]={0,0,0};
		int pdgdecaylist[]={2212,2112,22,211,321,11,13,80000};
		TLorentzVector VM_particle_of_interest(0.,0.,0.,0.);
		for(int j(0); j < nParticles; ++j ) {

			const erhic::ParticleMC* particle = event->GetTrack(j);

			int pdg = particle->GetPdgCode();
			int status = particle->GetStatus();
			int index = particle->GetIndex();//index 1 and 2 are incoming particle electron and proton.
			double pt = particle->GetPt();
			double eta = particle->GetEta();
			double phi = particle->GetPhi();
			double rap = particle->GetRapidity();
			double mass = particle->GetM();
			double theta = particle->GetTheta(); 
			double mom = particle->GetP();
			int charge = particle->eA->charge;
			int NoBAM = particle->eA->NoBam;
			if( index==3 ) e_scattered = particle->Get4Vector();

			//do analysis track-by-track
			for(int ivm=0;ivm<3;ivm++){
				if(pdg!=pdglist[ivm]) continue;
				if(status!=statuslist[ivm]&&status!=statuslist[ivm]+1) continue;
				hasvm[ivm]=1;//found vm.
				vm_vect[ivm]=particle->Get4Vector();
				
				int daug1=particle->GetChild1Index()-1;//Beagle list index starts at 1.
				int daug2=particle->GetChildNIndex()-1;
				if(daug1==-1 || daug2==-1) continue;
			
				const erhic::ParticleMC* particle_daug1 = event->GetTrack(daug1);
				const erhic::ParticleMC* particle_daug2 = event->GetTrack(daug2);

				//stable particles for daughters.
				if(particle_daug1->GetStatus()!=1 ||
					particle_daug2->GetStatus()!=1 ) continue;

				//for Jpsi only no mumu pairs.
				if( TMath::Abs(particle_daug1->GetPdgCode()) == 13 ||
					TMath::Abs(particle_daug2->GetPdgCode()) == 13 ) {
					acceptance[ivm]=0;
					continue;
				}
				if(ivm==0){
					if(TMath::Abs(particle_daug1->GetPdgCode()) != 211 ||
					TMath::Abs(particle_daug2->GetPdgCode()) != 211 ){
					acceptance[ivm]=0;
					continue;
					}
				}
				if(ivm==1){//phi meson
					if(TMath::Abs(particle_daug1->GetPdgCode()) != 321 ||
					TMath::Abs(particle_daug2->GetPdgCode()) != 321 ){
					acceptance[ivm]=0;
					continue;
					}
				}

				h_VM[processindex][ivm][0]->Fill(pt);
				h_VM[processindex][ivm][1]->Fill(eta);
				h_VM[processindex][ivm][2]->Fill(phi);
				h_VM[processindex][ivm][3]->Fill(theta);

				if(TMath::Abs(particle_daug1->GetEta())>4.0||
					TMath::Abs(particle_daug2->GetEta())>4.0) acceptance[ivm]=0;

				h_VM_daughter[processindex][ivm][0]->Fill(particle_daug1->GetPt());
				h_VM_daughter[processindex][ivm][1]->Fill(particle_daug1->GetEta());
				h_VM_daughter[processindex][ivm][2]->Fill(particle_daug1->GetPhi());
				h_VM_daughter[processindex][ivm][3]->Fill(particle_daug1->GetTheta());

				h_VM_daughter[processindex][ivm][0]->Fill(particle_daug2->GetPt());
				h_VM_daughter[processindex][ivm][1]->Fill(particle_daug2->GetEta());
				h_VM_daughter[processindex][ivm][2]->Fill(particle_daug2->GetPhi());
				h_VM_daughter[processindex][ivm][3]->Fill(particle_daug2->GetTheta());

				VM_particle_of_interest = particle_daug1->Get4Vector() + particle_daug2->Get4Vector();
				h_VM_mass[processindex][ivm]->Fill( VM_particle_of_interest.M() );

			}

		} // end of particle loop
		
		//for each vm; do...
		for(int ivm=0;ivm<3;ivm++){
			if(acceptance[ivm]&&hasvm[ivm]) {
				h_VM[processindex][ivm][4]->Fill(-t_hat);
				for(int imethod=0;imethod<3;imethod++){
					double t_reco = giveMe_t(imethod,e_beam,e_scattered,A_beam,vm_vect[ivm]);
					h_t_reco[processindex][ivm][imethod]->Fill( t_reco );
					if(imethod==2)h_Amass[processindex][ivm]->Fill(t_reco,giveMe_Amass(e_beam,e_scattered,A_beam,vm_vect[ivm]));
				}
				//loop over particle again
				for(int j(0); j < nParticles; ++j ) {
					const erhic::ParticleMC* particle = event->GetTrack(j);
					int pdg = particle->GetPdgCode();
					int status = particle->GetStatus();
					int index = particle->GetIndex();
					double mom = particle->GetP();
					double theta = particle->GetTheta();theta=theta*1000.;
					//fill breakup particles
					for(int ipid=0;ipid<8;ipid++){
						if(ipid<7){
							if( TMath::Abs(pdg) == pdgdecaylist[ipid]
								&& status==1 ){
								h_part[processindex][ivm][ipid]->Fill(TMath::Log(mom), TMath::Log(theta) );
							}
						}
						else{
							if( TMath::Abs(pdg) > pdgdecaylist[ipid]
								&& status==1 ){
								h_part[processindex][ivm][ipid]->Fill(TMath::Log(mom), TMath::Log(theta) );
							}
						}
						
					}
					//end filling breakup particles, and it doesn't matter if they are VM decays.
				}

			}
		}
		//end each vm;

	}

	output->Write();
	output->Close();


}
