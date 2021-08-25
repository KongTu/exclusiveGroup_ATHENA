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
#include "TRandom2.h"
#include "Math/Interpolator.h"
#include "Math/InterpolationTypes.h"

#define PI            3.1415926

#define MASS_PION     0.13957
#define MASS_KAON     0.493667
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
double daughtermasslist[]={MASS_PION,MASS_KAON,MASS_ELECTRON};

using namespace std;
using namespace erhic;
//from zhangbu
TFile* PIDinput = new TFile("../include/PIDchi2.root","READ");
TH2D* hist_pion = (TH2D*) PIDinput->Get("hist_pion");
TH2D* hist_kaon = (TH2D*) PIDinput->Get("hist_kaon");
bool veto_this_event(EventBeagle* event, int nParticles, int step_=-1){

	bool veto[] = {false,false,false,false,false,false,false};
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
		if( particle->GetParentIndex()==3 )continue;

		if( charge==0 ){//neutral particles, pi0, photons, neutrons
			if(theta<4.5 && pdg==2112 ) veto[6]=true;
			if(theta<4.5 && mom>0.05 ) veto[1]=true;
			if(theta>5.5 && theta<20. && mom>0.05) veto[2]=true;
		}
		else if( charge!=0 ){//charged particles
			double rigidity_ratio = (2.* mass)/(MASS_PROTON+MASS_NEUTRON) / TMath::Abs(charge);
			rigidity_ratio = rigidity_ratio / 2.5;//Au beam rigidity
			if( theta>0. && theta<5.0 
				&& rigidity_ratio>0.45 && rigidity_ratio<0.65 ) veto[3]=true;
			if( theta>0. && theta<5.0 
				&& rigidity_ratio>0.7 && rigidity_ratio<0.95 ) veto[4]=true;
			if( theta>5.5 && theta<20.0 ) veto[5]=true;

			if(TMath::Abs(eta)<4.0 && pt>0.15) multiplicity++;
		}
		else{
			cout << "something's wrong about the charge" << endl;
		}
	}

	if(multiplicity>4||multiplicity<2) veto[0]=true;//first veto

	bool anyVeto = (veto[0]||veto[1]||veto[2]||veto[3]||veto[4]||veto[5]);
	if( step_ == -1 ) return anyVeto;
	else return veto[step_];
}
//
double giveMe_t(int option, TLorentzVector e_beam, TLorentzVector e_scattered, TLorentzVector p_beam, TLorentzVector vm_vect){

	double method_E = (-vm_vect-e_scattered+e_beam).Mag2();
	
	double method_A = -99;
	TVector2 sum_pt(vm_vect.Px()+e_scattered.Px(), vm_vect.Py()+e_scattered.Py());
	method_A = sum_pt.Mod2();
	
	double method_L = -99.;
	TLorentzVector p_beam_scattered = p_beam-(vm_vect+e_scattered-e_beam);
	double p_Aplus = p_beam_scattered.E()+p_beam_scattered.Pz();
	double p_TAsquared = TMath::Power(p_beam_scattered.Pt(),2);
	double p_Aminus = (MASS_AU197*MASS_AU197 + p_TAsquared) / p_Aplus;
	TLorentzVector p_beam_scattered_corr; 
	p_beam_scattered_corr.SetPxPyPzE(p_beam_scattered.Px(),p_beam_scattered.Py(),(p_Aplus-p_Aminus)/2., (p_Aplus+p_Aminus)/2. );
	method_L = (p_beam_scattered_corr-p_beam).Mag2();

	if(option==0) return -method_E;
	else if(option==1) return method_A;
	else if(option==2) return -method_L;
	else return -99;
}
//
double giveMe_Amass(TLorentzVector e_beam, TLorentzVector e_scattered, TLorentzVector p_beam, TLorentzVector vm_vect){
	TLorentzVector p_beam_scattered = p_beam-(vm_vect+e_scattered-e_beam);
	return p_beam_scattered.M()/MASS_AU197;
}
//
double giveMe_PIDChi2(TLorentzVector v, TH2D* hist){

	double p = v.P();
	TH1D* h_total_projection = (TH1D*) hist->ProjectionX("h_total_projection",1,1e8);//total y bins
	int bin_of_interest = h_total_projection->FindBin(p);
	TH1D* h_normChi2_1D = (TH1D*) hist->ProjectionY("h_normChi2_1D",bin_of_interest,bin_of_interest);
	if(bin_of_interest<20) return -99;
		double PID = h_normChi2_1D->GetRandom();
	return PID;

}
//
void printSTABLE(EventBeagle* event, int nParticles){
	for(int j(0); j < nParticles; ++j ) {
		const erhic::ParticleMC* particle = event->GetTrack(j);
		int pdg = particle->GetPdgCode();
		int status = particle->GetStatus();
		int index = particle->GetIndex();//index 1 and 2 are incoming particle electron and proton.
		double pt = particle->GetPt();
		double eta = particle->GetEta();
		double phi = particle->GetPhi();
		if(status!=1) continue;
		cout << "PDG = " << pdg << " , pt = " << pt << " , eta = " << eta << " , phi = " << phi << endl;  
	}

}



