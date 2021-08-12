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

bool veto_this_event(EventBeagle* event, int nParticles){

	bool veto = false;
	int multiplicity=0;
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
		theta=theta*1000.; //mrad
		double mom = particle->GetP();
		int charge = particle->eA->charge;
		int NoBAM = particle->eA->NoBam;

		if( status!= 1 ) continue;
		if( charge==0 ){//neutral particles, pi0, photons, neutrons
			if(theta<4.5 && mom>0.05 ) veto=true;
			if(theta>5.5 && theta<20. && mom>0.05) veto=true;
		}
		else if( charge!=0 ){//charged particles
			double rigidity_ratio = (2.* mass)/(MASS_PROTON+MASS_NEUTRON) / TMath::Abs(charge);
			rigidity_ratio = rigidity_ratio / 2.5;//Au beam rigidity
			if( theta>0. && theta<5.0 
				&& rigidity_ratio>0.45 && rigidity_ratio<0.65 ) veto=true;
			if( theta>0. && theta<5.0 
				&& rigidity_ratio>0.7 && rigidity_ratio<0.95 ) veto=true;
			if( theta>5.5 && theta<20.0 ) veto=true;

			if(TMath::Abs(eta)<4.0) multiplicity++;
		}
		else{
			cout << "something's wrong about the charge" << endl;
		}
	}

	if(multiplicity>3) veto=true;

	return veto;
}




