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
double pathLength(double pt, double p){
  double LGADTOF = 50.0; //cm
  double BField = 3.0; //Tesla 
  double sintheta=LGADTOF*0.003*BField/2.0/pt;
  if (sintheta>1.0) return 0.0;
  double arc = 2.0*p/0.003/BField*TMath::ASin(sintheta);
  return arc; 
}
double giveMe_PIDChi2(TLorentzVector v1, TLorentzVector v2, double mass){

	double tofRes = 25.0;//picoseconds 
    double startRes = 30.0;//scattered electron timing picoseconds 

    double tof1 = pathLength(v1.Pt(),v1.P())/v1.Beta()*1000.0/30.0;//picoseconds 
    double tof2 = pathLength(v2.Pt(),v2.P())/v2.Beta()*1000.0/30.0;//picoseconds 
	
	if(tof1==0.||tof2==0.) return -99.;

	if(fabs(mass-MASS_PION)<1E-2) mass=MASS_PION;
	if(fabs(mass-MASS_KAON)<1E-2) mass=MASS_KAON;
	
	TLorentzVector v_real1,v_real2;
	v_real1.SetVectM(v1.Vect(), mass);
	v_real2.SetVectM(v2.Vect(), mass);
	
	double tof3 = pathLength(v_real1.Pt(),v_real1.P())/v_real1.Beta()*1000.0/30.0;//picoseconds 
    double tof4 = pathLength(v_real2.Pt(),v_real2.P())/v_real2.Beta()*1000.0/30.0;//picoseconds
	
	double starttiming = gRandom->Gaus(0.0,startRes);
    double timesmear1 =  (gRandom->Gaus(0.0,tofRes)+starttiming);
    double timesmear2 =  (gRandom->Gaus(0.0,tofRes)+starttiming);
    tof3 += timesmear1;//picoseconds 
    tof4 += timesmear2;//picoseconds 
    // cout << "to1,tof2,tof3,tof4 = " << tof1 << " " << tof2 << " " << tof3 << " " << tof4 << endl;
    double chi2 = 1.0/pow(tofRes,2)*(pow(tof1-tof3,2)+pow(tof2-tof4,2)-1.0/pow(tofRes,2)*pow(tof1+tof2-tof3-tof4,2)/(2.0/pow(tofRes,2)+1.0/pow(startRes,2)));
	// cout << "chi2 ? " << chi2 << endl;
	return chi2;
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



