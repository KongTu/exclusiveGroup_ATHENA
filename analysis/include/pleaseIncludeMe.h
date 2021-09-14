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
#include "TGeoMatrix.h"
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

double minPt_=0.1;

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
			if(theta<2.5 && pdg==2112 ) veto[6]=true;
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
	if( e_scattered.Eta() < -4.0 ) method_E = -0.01;//photoproduction protection

	double method_A = -99;
	TVector2 sum_pt(vm_vect.Px()+e_scattered.Px(), vm_vect.Py()+e_scattered.Py());
	method_A = sum_pt.Mod2();
	if( e_scattered.Eta() < -4.0 ) method_A = vm_vect.Pt()*vm_vect.Pt();//photoproduction;
	
	double method_L = -99.;
	TLorentzVector p_beam_scattered = p_beam-(vm_vect+e_scattered-e_beam);
	double p_Aplus = p_beam_scattered.E()+p_beam_scattered.Pz();
	double p_TAsquared = TMath::Power(p_beam_scattered.Pt(),2);
	double p_Aminus = (MASS_AU197*MASS_AU197 + p_TAsquared) / p_Aplus;
	TLorentzVector p_beam_scattered_corr; 
	p_beam_scattered_corr.SetPxPyPzE(p_beam_scattered.Px(),p_beam_scattered.Py(),(p_Aplus-p_Aminus)/2., (p_Aplus+p_Aminus)/2. );
	method_L = (p_beam_scattered_corr-p_beam).Mag2();
	if( e_scattered.Eta() < -4.0 ) method_L = -0.05;//photoproduction protection

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
  double BField = 0.;//Tesla
  if(minPt_>0.3) BField=3.0;
  if(minPt_>=0.15 && minPt_<=0.3) BField=1.5;
  if(minPt_<0.15 && minPt_>=0.07) BField=0.5;
  double LGADTOF = 50.0; //cm
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

vector<TLorentzVector> letsMakeItReal(TLorentzVector e_beam, TLorentzVector e_scattered,
 TLorentzVector A_beam, TLorentzVector daug_1, TLorentzVector daug_2){

		/*Smearing includes:
		- crossing angle 25 mrad. [not yet]
		- angular divergence of the beams.
		- e, e', VM daughters (2-prone), Au.
		*/


		//1. angular divergence: CDR
		double momentum_resolution_e = 10.9E-4;// GeV/c
		double momentum_resolution_Au = 6.2E-4;// GeV/c
		double theta_resolution_e[2]={0.101E-3,0.037E-3};//x,y rad
		double theta_resolution_h[2]={0.218E-3,0.379E-3};//x,y rad, w. strong hadron cooling
	
		//e incoming beam
		double p = e_beam.Pz();
		double px = TMath::Sin( gRandom->Gaus(0.0,theta_resolution_e[0]) ) * p;
		double py = TMath::Sin( gRandom->Gaus(0.0,theta_resolution_e[1]) ) * p;
		double theta = TMath::ASin(sqrt(px*px+py*py)/p);
		double pz = p*TMath::Cos(theta);
		pz = (1.+gRandom->Gaus(0.,momentum_resolution_e))*pz;
		TLorentzVector e_beam_smear(px, py, pz, sqrt(px*px+py*py+pz*pz+MASS_ELECTRON*MASS_ELECTRON));
		e_beam = e_beam_smear;
		
		// // Boost the system around doesn't work?
		// TVector3 e_beam_boost = e_beam.BoostVector();
		// TVector3 e_beam_reverse_boost = e_beam_smear.BoostVector();
		// e_scattered.Boost(-e_beam_boost);
		// e_scattered.Boost(e_beam_reverse_boost);

		//Modify scattered e'
		double simComp[3]; 
		simComp[0] = e_scattered.Px(); 
		simComp[1] = e_scattered.Py(); 
		simComp[2] = e_scattered.Pz();		
	
		TGeoRotation *horizDiv = new TGeoRotation();
		TGeoRotation *vertDiv  = new TGeoRotation();
		
		double horizAngle = gRandom->Gaus(0.0, theta_resolution_e[0]);   
		double vertAngle  = gRandom->Gaus(0.0, theta_resolution_e[1]); 

		horizDiv->RotateY(horizAngle);
		vertDiv->RotateX(vertAngle);

		double angDivOutHoriz[3], angDivOutHorizAndVert[3];
		horizDiv->MasterToLocalVect(simComp, angDivOutHoriz);
		vertDiv->MasterToLocalVect(angDivOutHoriz, angDivOutHorizAndVert);

		TVector3 e_scattered_modified(angDivOutHorizAndVert);
    e_scattered.SetVectM(e_scattered_modified,MASS_ELECTRON);

    //A incoming beam.
		TVector3 A_beam_boost = A_beam.BoostVector();
		p = A_beam.Pz();
		px = TMath::Sin(gRandom->Gaus(0.0,theta_resolution_h[0])) * p;
		py = TMath::Sin(gRandom->Gaus(0.0,theta_resolution_h[1])) * p;
		theta = TMath::ASin(sqrt(px*px+py*py)/p);
		pz = p*TMath::Cos(theta);
		pz = (1.+gRandom->Gaus(0.,momentum_resolution_Au))*pz;
		TLorentzVector A_beam_smear(px, py, pz, sqrt(px*px+py*py+pz*pz+MASS_AU197*MASS_AU197));
		A_beam = A_beam_smear;

		//2. pt resolution. YR. page 351.
		// B=3T
		// double pt_resolution[]={0.0005,0.001,0.001};
		// double pt_resolution_constant[]={0.005,0.005,0.005};
		// B=1.5T
		// double pt_resolution[]={0.001,0.002,0.002};
		// double pt_resolution_constant[]={0.01,0.01,0.01};
		// B = 0.75T
		double pt_resolution[]={0.002,0.004,0.004};
		double pt_resolution_constant[]={0.02,0.02,0.02};

		double eta_bins[]={0.0,1.0,2.5,4.0};
		int pt_index_e = -1;
		int pt_index_daug_1 = -1;
		int pt_index_daug_2 = -1;
		for(int i=0;i<3;i++){
			if(fabs(e_scattered.Eta())>eta_bins[i] 
				&& fabs(e_scattered.Eta()) < eta_bins[i+1]) pt_index_e = i;

			if(fabs(daug_1.Eta())>eta_bins[i] 
				&& fabs(daug_1.Eta()) < eta_bins[i+1]) pt_index_daug_1 = i;

			if(fabs(daug_2.Eta())>eta_bins[i] 
				&& fabs(daug_2.Eta()) < eta_bins[i+1]) pt_index_daug_2 = i;
		}
		//e_scattered:
		double pt_e_scattered = -99.;
		if(pt_index_e>=0) {
			double width = sqrt(TMath::Power(e_scattered.Pt()*pt_resolution[pt_index_e],2)+TMath::Power(pt_resolution_constant[pt_index_e],2));
			double resolution = gRandom->Gaus(0.0,width);
			pt_e_scattered = (1.+resolution)*e_scattered.Pt();
			e_scattered.SetPtEtaPhiM(pt_e_scattered,e_scattered.Eta(),e_scattered.Phi(),e_scattered.M());
		}
		//daughter 1:
		double pt_daug_1 = -99.;
		if(pt_index_daug_1>=0) {
			double width = sqrt(TMath::Power(daug_1.Pt()*pt_resolution[pt_index_daug_1],2)+TMath::Power(pt_resolution_constant[pt_index_daug_1],2));
			double resolution = gRandom->Gaus(0.0,width);
			pt_daug_1 = (1.+resolution)*daug_1.Pt();
			daug_1.SetPtEtaPhiM(pt_daug_1,daug_1.Eta(),daug_1.Phi(),daug_1.M());
		}
		//daughter 2:
		double pt_daug_2 = -99.;
		if(pt_index_daug_2>=0) {
			double width = sqrt(TMath::Power(daug_2.Pt()*pt_resolution[pt_index_daug_2],2)+TMath::Power(pt_resolution_constant[pt_index_daug_2],2));
			double resolution = gRandom->Gaus(0.0,width);
			pt_daug_2 = (1.+resolution)*daug_2.Pt();
			daug_2.SetPtEtaPhiM(pt_daug_2,daug_2.Eta(),daug_2.Phi(),daug_2.M());
		}
		
		vector<TLorentzVector > update;
		update.push_back(e_beam);
		update.push_back(e_scattered);
		update.push_back(A_beam);
		update.push_back(daug_1);
		update.push_back(daug_2);

		return update;
}

