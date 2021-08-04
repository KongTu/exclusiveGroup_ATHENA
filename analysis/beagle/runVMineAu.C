#include <iostream>
#include <vector>
#include <sstream>
#include <string>

#include <eicsmear/erhic/EventBase.h>
#include <eicsmear/erhic/EventMC.h>
#include <eicsmear/erhic/EventPythia.h>
#include <eicsmear/erhic/Particle.h>
#include <eicsmear/erhic/ParticleMC.h>
#include <eicsmear/erhic/Pid.h>

#include "TString.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TAxis.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLorentzVector.h"
#include "TBranchElement.h"

#define PI            3.1415926

#define MASS_MUON     0.1056
#define MASS_ELECTRON 0.000511
#define MASS_JPSI 	  3.09688
#define MASS_PROTON   0.93827
#define MASS_NEUTRON  0.93957
#define MASS_DEUTERON 1.8756129
#define MASS_TRITON   2.7937167208086358
#define MASS_HE3      2.7937167208086358
#define MASS_ALPHA    3.7249556277448477
#define MASS_LI6      5.5874334416172715
#define MASS_C12      11.174866883234543
#define MASS_CA40     37.249556277448477
#define MASS_XE131    121.99229680864376
#define MASS_AU197    183.45406466643374
#define MASS_PB208    193.69769264273208

using namespace std;
using namespace erhic;


void runVMineAu(const TString filename="eA_TEST", const int nEvents = 40000){

	TChain *tree = new TChain("EICTree");
	tree->Add( filename+".root" );
	
	EventBeagle* event(NULL);
	tree->SetBranchAddress("event", &event);

	TFile* output = new TFile("../../rootfiles/beagle_phi.root","RECREATE");
	TH1D* h_trueT = new TH1D("h_trueT",";-t (GeV^{2})", 100,0,0.5);
	//VM histograms//
	/* first   index VM process, 91=0, 93=1, 91+93=2*/
	/* second  index VM species, rho=0, phi=1, jpsi=2*/
	/* third   index VM property, pt=0, eta=1, phi=2, theta=3, reserved=4*/
	double bin_lower[]={0.,-8.,0.,0.,0.};
	double bin_upper[]={5.0,8.,6.5,100.,1.2};
	TH1D* h_VM[2][3][5];
	for(int ibreak=0;ibreak<2;ibreak++){
		for(int ivm=0;ivm<3;ivm++){
			for(int ipro=0;ipro<5;ipro++){
				h_VM[ibreak][ivm][ipro] = new TH1D(Form("h_VM_%d_%d_%d",ibreak,ivm,ipro),
					Form("h_VM_%d_%d_%d",ibreak,ivm,ipro),100,bin_lower[ipro],bin_upper[ipro] );
			}
		}
	}
	//END VM histograms//

	//VM daughter histograms//
	TH1D* h_VM_daughter[2][3][5];
	for(int ibreak=0;ibreak<2;ibreak++){
		for(int ivm=0;ivm<3;ivm++){
			for(int ipro=0;ipro<5;ipro++){
				h_VM_daughter[ibreak][ivm][ipro] = new TH1D(Form("h_VM_daughter_%d_%d_%d",ibreak,ivm,ipro),
					Form("h_VM_daughter_%d_%d_%d",ibreak,ivm,ipro),100,bin_lower[ipro],bin_upper[ipro] );
			}
		}
	}
	//END VM daughter histograms
	for(int i(0); i < nEvents; ++i ) {
      
		// Read the next entry from the tree.
		tree->GetEntry(i);
		
		double pzlep = event->pzlep;
		double pztarg = event->pztarg;
		int struck_nucleon = event->nucleon;
		double MASS_NUCLEON = MASS_PROTON;
		if( struck_nucleon==2112 ) MASS_NUCLEON = MASS_NEUTRON;

		TLorentzVector e_beam(0.,0.,pzlep,sqrt(pzlep*pzlep+MASS_ELECTRON*MASS_ELECTRON));
		TLorentzVector p_beam(0.,0.,pztarg,sqrt(pztarg*pztarg+MASS_NUCLEON*MASS_NUCLEON));
		TLorentzVector e_scattered(0.,0.,0.,0.);

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

		//event cuts
		int processindex=-1;
		if( event_process==91) processindex=0;
		if( event_process==93) processindex=1;
		if( processindex<0 ) continue;
		if( trueQ2 < 1. || trueQ2 > 20. ) continue;
		if( trueY > 0.95 || trueY < 0.01 ) continue;

		//do analysis, or fill historgrams for event levels
		h_trueT->Fill(-t_hat);

		//particle loop
		//rho^0 = 113, decay->pipi
		//phi = 333, decay->kk
		//jpsi = 443, nodecay
		int pdglist[]={113,333,443};
		int statuslist[]={2,2,1};
		int multiplicity[3]={0,0,0};
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
			theta = theta*1000.0; //change to mrad;
			double mom = particle->GetP();
			int charge = particle->eA->charge;
			int NoBAM = particle->eA->NoBam;

			//do analysis track-by-track
			for(int ivm=0;ivm<3;ivm++){
				if(pdg!=pdglist[ivm]) continue;
				if(status!=statuslist[ivm]) continue;
				h_VM[processindex][ivm][0]->Fill(pt);
				h_VM[processindex][ivm][1]->Fill(rap);
				h_VM[processindex][ivm][2]->Fill(phi);
				h_VM[processindex][ivm][3]->Fill(theta);

				//rho and phi daughters:
				if(ivm<2){
					int daug1=particle->GetChild1Index();
					int daug2=particle->GetChildNIndex();
					if(daug1==0 || daug2==0) continue;
					cout << "daug1~"<<daug1<<endl;
					cout << "daug2~"<<daug2<<endl;

					const erhic::ParticleMC* particle_daug1 = event->GetTrack(daug1);
					const erhic::ParticleMC* particle_daug2 = event->GetTrack(daug2);

					h_VM_daughter[processindex][ivm][0]->Fill(particle_daug1->GetPt());
					h_VM_daughter[processindex][ivm][1]->Fill(particle_daug1->GetEta());
					h_VM_daughter[processindex][ivm][2]->Fill(particle_daug1->GetPhi());
					h_VM_daughter[processindex][ivm][3]->Fill(particle_daug1->GetTheta()*1000);

					h_VM_daughter[processindex][ivm][0]->Fill(particle_daug2->GetPt());
					h_VM_daughter[processindex][ivm][1]->Fill(particle_daug2->GetEta());
					h_VM_daughter[processindex][ivm][2]->Fill(particle_daug2->GetPhi());
					h_VM_daughter[processindex][ivm][3]->Fill(particle_daug2->GetTheta()*1000);

					cout << "pt daug1 " << particle_daug1->GetPt() << endl;
					cout << "pt daug2 " << particle_daug2->GetPt() << endl;

				}
			}

		} // end of particle loop
		for(int ivm=0;ivm<3;ivm++){h_VM[processindex][ivm][4]->Fill(-t_hat);}
		cout << "Event # " << i << endl;
	}

	output->Write();
	output->Close();


}
