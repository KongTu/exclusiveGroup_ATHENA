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


void runVMinePb(const TString filename="eA_TEST", const int nEvents = 40000){

	TChain *tree = new TChain("EICTree");
	tree->Add( filename+".root" );
	
	EventBeagle* event(NULL);
	tree->SetBranchAddress("event", &event);

	TFile* output = new TFile("output.root","RECREATE");
	TH1D* h_trueT = new TH1D("h_trueT",";-t (GeV^{2})", 100,0,1);
	TH1D* h_decomposition = new TH1D("h_decomposition",";",10,0,10);
	TH2D* h_thetaVsMom[3];
	TH2D* h_thetaVsMom_min[3];
	TH2D* h_thetaVsBam[3];
	TH2D* h_MomVsBam[3];
	for(int k=0;k<2;k++){
		h_thetaVsMom[k] = new TH2D(Form("h_thetaVsMom_%d",k),";p (GeV);#theta (mrad)",2500,0,250,10000,0,1000);
		h_thetaVsMom_min[k] = new TH2D(Form("h_thetaVsMom_min_%d",k),";p (GeV);#theta (mrad)",2500,0,250,10000,0,1000);
		h_thetaVsBam[k] = new TH2D(Form("h_thetaVsBam_%d",k),";NoBAM;#theta (mrad)",37,-1,36,10000,0,1000);
		h_MomVsBam[k] = new TH2D(Form("h_MomVsBam_%d",k),";NoBAM;p (GeV)",37,-1,36,2500,0,250);
	}
		h_thetaVsMom[2] = new TH2D(Form("h_thetaVsMom_%d",2),";p (GeV);#theta (mrad)",2500,0,20,10000,0,1000);
		h_thetaVsMom_min[2] = new TH2D(Form("h_thetaVsMom_min_%d",2),";p (GeV);#theta (mrad)",2500,0,20,10000,0,1000);
		h_thetaVsBam[2] = new TH2D(Form("h_thetaVsBam_%d",2),";NoBAM;#theta (mrad)",37,-1,36,10000,0,1000);
		h_MomVsBam[2] = new TH2D(Form("h_MomVsBam_%d",2),";NoBAM;p (GeV)",37,-1,36,2500,0,250);

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
		if( event_process != 91 && event_process != 93 ) continue;
		if( trueQ2 < 1. || trueQ2 > 100. ) continue;
		if( trueY > 0.95 || trueY < 0.01 ) continue;

		//do analysis, or fill historgrams for event levels

		//particle loop
		bool hasJpsi = false;
		bool hasNeutron = false;
		bool hasProton = false;
		bool hasPhoton = false;
		vector< double> angle_neutron, angle_proton, angle_photon;
		vector< double> momentum_neutron, momentum_proton, momentum_photon;
		vector< int> bam_neutron, bam_proton, bam_photon;
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
			int NoBAM = particle->eA->NoBam;

			//only stable particles or j/psi.
			if( status != 1 && pdg != 443 ) continue;

			if( pdg == 443 ){ //jpsi
				hasJpsi = true;
			}
			if( pdg == 2112 ){ // neutrons
				hasNeutron = true;
				angle_neutron.push_back( theta );
				momentum_neutron.push_back( mom );
				bam_neutron.push_back( NoBAM );
			}
			if( pdg == 2212 ){ // proton
				hasProton = true;
				angle_proton.push_back( theta );
				momentum_proton.push_back( mom );
				bam_proton.push_back( NoBAM );

			}
			if( pdg == 22 ){ // photon
				hasPhoton = true;
				angle_photon.push_back( theta );
				momentum_photon.push_back( mom );
				bam_photon.push_back( NoBAM );
			}
			//do analysis track-by-track

		} // end of particle loop


		if( hasJpsi ) {
			if( hasNeutron && !hasProton && !hasPhoton )h_decomposition->Fill(0);
			if( !hasNeutron && hasProton && !hasPhoton )h_decomposition->Fill(1);
			if( !hasNeutron && !hasProton && hasPhoton )h_decomposition->Fill(2);
			if( hasNeutron && hasProton && !hasPhoton )h_decomposition->Fill(3);
			if( hasNeutron && !hasProton && hasPhoton )h_decomposition->Fill(4);
			if( !hasNeutron && hasProton && hasPhoton )h_decomposition->Fill(5);
			if( hasNeutron && hasProton && hasPhoton )h_decomposition->Fill(6);

			h_trueT->Fill( -t_hat );
			int index_min=-1;
			double angle_min=999.;
			for(unsigned ipart=0;ipart<angle_neutron.size();ipart++){
				if( angle_neutron[ipart] < angle_min) {
					angle_min=angle_neutron[ipart];
					index_min=ipart;
				}
				h_thetaVsMom[0]->Fill(momentum_neutron[ipart],angle_neutron[ipart]);
				h_thetaVsBam[0]->Fill(bam_neutron[ipart], angle_neutron[ipart]);
				h_MomVsBam[0]->Fill(bam_neutron[ipart], momentum_neutron[ipart]);
			}
			if(index_min!=-1) h_thetaVsMom_min[0]->Fill( momentum_neutron[index_min],angle_neutron[index_min] );
			index_min = -1;
			angle_min=999.;
			for(unsigned ipart=0;ipart<angle_proton.size();ipart++){
				if( angle_proton[ipart] < angle_min) {
					angle_min=angle_proton[ipart];
					index_min=ipart;
				}
				h_thetaVsMom[1]->Fill(momentum_proton[ipart],angle_proton[ipart]);
				h_thetaVsBam[1]->Fill(bam_proton[ipart], angle_proton[ipart]);
				h_MomVsBam[1]->Fill(bam_proton[ipart], momentum_proton[ipart]);
			}
			if(index_min!=-1) h_thetaVsMom_min[1]->Fill( momentum_proton[index_min],angle_proton[index_min] );
			index_min = -1;
			angle_min=999.;
			for(unsigned ipart=0;ipart<angle_photon.size();ipart++){
				if( angle_photon[ipart] < angle_min) {
					angle_min=angle_photon[ipart];
					index_min=ipart;
				}
				h_thetaVsMom[2]->Fill(momentum_photon[ipart],angle_photon[ipart]);
				h_thetaVsBam[2]->Fill(bam_photon[ipart], angle_photon[ipart]);
				h_MomVsBam[2]->Fill(bam_photon[ipart], momentum_photon[ipart]);
			}
			if(index_min!=-1) h_thetaVsMom_min[2]->Fill( momentum_photon[index_min],angle_photon[index_min] );
		}
		//fill histograms
	}

	output->Write();
	output->Close();


}
