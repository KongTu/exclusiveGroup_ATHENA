// ==============================================================================
//  File: readerExample.C
//
//  Analyze ROOT tree produced by Sartre's standard
//  implementation 'sartreMain'.
//
//  Here:Template/Example only
//
//  Author: Thomas Ullrich
//  Last update: June 4, 2021
//
// ------------------------------------------------------------------------------
//
//
// ==============================================================================
#include <iostream>
#include <fstream>
#include <vector>
#include "TObject.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TTree.h"
#include "TMath.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TChain.h"
#include "TLorentzVector.h"
#include "Math/Interpolator.h"
#include "Math/InterpolationTypes.h"
#include "TRandom2.h"
#include "TText.h"
#include "TGraph.h"
//#define PR(x) cout << #x << " = " << (x) << endl;

using namespace std;

//
//  Event structure
//
struct rootSartreEvent { // Event structure as used in Sartre
    double t;
    double Q2;
    double x;
    double s;
    double y;
    double W;
    double xpom;
    int    iEvent;
    int    pol;      // 0=transverse or 1=longitudinal
    int    dmode;    // 0=coherent, 1=Incoherent
};

//
//  Helper functions and operators
//
ostream& operator<<(ostream& os, const TLorentzVector& v)
{
    os << v.Px() << '\t' << v.Py() << '\t'  << v.Pz() << '\t'  << v.E() << '\t';
    double m2 = v*v;
    if (m2 < 0)
        os << '(' << -sqrt(-m2) << ')';
    else
        os << '(' << sqrt(m2) << ')';

    return os;
}

double RhoSpectral(){//STAR PRC 2017 parameters
  double M0rho = 0.7765;
  double Mpi = 0.139;
  double W0rho = 0.156;
  double Arho = 1.538;
  double B2pi = -1.21;
  double M2pi = 0.278;
  double Mrho=0.278;
  double MaxSpectral = 1.1*(1.0/W0rho+pow(B2pi/Arho,2));
  Mrho = M0rho;
  double Mspectral =0.0;
  while (Mspectral<gRandom->Uniform(0.0,MaxSpectral)){
    Mrho = gRandom->Uniform(M2pi,1.5);
    double Wrho = W0rho*M0rho/Mrho*pow(((pow(Mrho,2)-4*pow(Mpi,2))/(pow(M0rho,2)-4*pow(Mpi,2))),1.5);
    double PoleFactor = pow(Mrho*M0rho*Wrho, 0.5)/(pow(Mrho*Mrho-M0rho*M0rho,2)+M0rho*M0rho*Wrho*Wrho);
    Mspectral = pow((PoleFactor*(Mrho*Mrho-M0rho*M0rho)+B2pi/Arho),2)+pow(PoleFactor*Mrho*Wrho,2);
  }
  return Mrho;
}


// boost a daughter vector according to the parent (lab) vector boost
void applyBoost( TLorentzVector &_parent_lv, TLorentzVector &_d_lv ){

    float betaX = _parent_lv.Px() / _parent_lv.E();
    float betaY = _parent_lv.Py() / _parent_lv.E();
    float betaZ = _parent_lv.Pz() / _parent_lv.E();

    _d_lv.Boost(betaX,betaY,betaZ);
}

/* Two Body Decay in the rest frame of the parent
 *
 * The decay is computed and then the daughters are boosted into the frame of the parent
 */
vector<TLorentzVector> twoBodyDecay( TLorentzVector _parent_lv, double M_d1, double M_d2 ){
    // cout<< "twoBodyDecay( lv, M_d1=" << M_d1 << ", M_d2=" << M_d2 << ")" << endl;

    double dm2 = pow(M_d1 - M_d2,2);
    double sm2 = pow(M_d1 + M_d2,2);

    double M_p = _parent_lv.M();
    double M_p2 = pow(_parent_lv.M(), 2);

    double E1 = (M_p2 + pow(M_d1,2) - pow(M_d2,2)) / ( 2 * M_p );
    double E2 = M_p - E1;

    // cout << "E1 = " << E1 << endl;
    // cout << "E2 = " << E2 << endl;

    double p1 = sqrt( pow(E1,2) - pow(M_d1,2) );
    double p2 = sqrt( pow(E2,2) - pow(M_d2,2) );

    

    double p = sqrt( (M_p2 - sm2) * (M_p2 - dm2) ) / ( 2. * _parent_lv.M() );

    // cout << "p1 = " << p1 << " and p2 = " << p2 << " vs. p = " << p << endl;


    // double theta = TMath::ACos( gRandom->Uniform(-1, 1) );
    double costheta = gRandom->Uniform(0,1.);

    double phi = gRandom->Uniform(0,TMath::Pi()*2);

    double pz = p1*costheta;
    double px = p1*sqrt(1.-costheta*costheta)*cos(phi);
    double py = p1*sqrt(1.-costheta*costheta)*sin(phi);

    TLorentzVector daughter1( px, py, pz, sqrt( p1*p1 + M_d1*M_d1 ));
    TLorentzVector daughter2( -px, -py, -pz, sqrt( p2*p2 + M_d2*M_d2 ) );

    TVector3 v = _parent_lv.Vect().Unit();
    daughter1.RotateUz( v );
    daughter2.RotateUz( v );

    applyBoost(_parent_lv,daughter1);
    applyBoost(_parent_lv,daughter2);

    vector<TLorentzVector> dlv;
    dlv.push_back( daughter1 );
    dlv.push_back( daughter2 );
    return dlv;
}

double pathLength(double pt, double p){
  double LGADTOF = 50.0; //cm
  double BField = 3.0; //Tesla 
  double sintheta=LGADTOF*0.003*BField/2.0/pt;
  if (sintheta>1.0) return 0.0;
  double arc = 2.0*p/0.003/BField*TMath::ASin(sintheta);
  //  cout<<pt<<"  "<<p<<"  "<<TMath::ASin(sintheta)<<"  pt arc sintheta "<<arc<<"  "<<sintheta<<endl;
  return arc; 
}
//===========================================================================
//  Main function
//===========================================================================

void readerRhoSpectral(double fractionOfEventsToRead = 1){
    //
    //  Setup filenames
    //
    string fnames[16];
    int nnames = 0;

    fnames[0] = "/star/u/ullrich/data/sartre/data/sartre_bnonsat_Au_rho_1.root";
    fnames[1] = "/star/u/ullrich/data/sartre/data/sartre_bnonsat_Au_rho_2.root";
    fnames[2] = "/star/u/ullrich/data/sartre/data/sartre_bnonsat_Au_rho_3.root";
    fnames[3] = "/star/u/ullrich/data/sartre/data/sartre_bnonsat_Au_rho_4.root";
    fnames[4] = "/star/u/ullrich/data/sartre/data/sartre_bnonsat_Au_rho_5.root";
    fnames[5] = "/star/u/ullrich/data/sartre/data/sartre_bnonsat_Au_rho_6.root";
    fnames[6] = "/star/u/ullrich/data/sartre/data/sartre_bnonsat_Au_rho_7.root";
    fnames[7] = "/star/u/ullrich/data/sartre/data/sartre_bnonsat_Au_rho_8.root";
    fnames[8] = "/star/u/ullrich/data/sartre/data/sartre_bnonsat_Au_rho_9.root";
    fnames[9] = "/star/u/ullrich/data/sartre/data/sartre_bnonsat_Au_rho_10.root";
    nnames = 10;
   
    //
    //   Histogram Booking (example)
    //
    TH1D *hist_t_coherent = new TH1D("hist_t_coherentRho2phi", "coherent", 72, 0, 0.18);
    TH1D *hist_t_incoherent = new TH1D("hist_t_incoherentRho2phi", "incoherent", 72, 0, 0.18);
    TH1D *hist_mass1_coherent = new TH1D("hist_mass_coherentRho", "Mass1Coherent", 1000, 0, 2.0);
    TH1D *hist_mass2_coherent = new TH1D("hist_mass_coherentRho2phi", "Mass2Coherent", 1000, 0.9, 2.9);
    TH1D *hist_tPID_coherent = new TH1D("hist_tPID_coherentRho2phi", "tPIDcoherent", 72, 0, 0.18);
    TH2D *hist_PIDChi2_coherent = new TH2D("hist_PIDChi2_coherentRho2phi", "PIDChi2coherent", 100,0.,3.0,500,0.0,100.0);
    TH2D *hist_PIDNormChi2_coherent = new TH2D("hist_PIDNormChi2_coherentRho2phi", "PIDNormChi2coherent", 100,0.,3.0,500,0.0,100.0);
//    TH1D *hist_tPID_coherent = new TH1D("hist_tPID_coherentRho2phi", "tPIDcoherent", 72, 0, 0.18);
//    TH1D *hist_tsmear2_coherent = new TH1D("hist_tsmear2_coherentRho2phi", "tsmear2coherent", 72, 0, 0.18);
//    TH1D *hist_tsmear3_coherent = new TH1D("hist_tsmear3_coherentRho2phi", "tsmear3coherent", 72, 0, 0.18);

    //
    //  Build chain
    //
    TFile *f = 0;
    TTree *tree = 0;
    TChain *chain = 0;

    cout << "Setup new chain" << endl;
    chain = new TChain("tree"); // name of the tree is the argument
    for (int iadd = 0; iadd < nnames; iadd++) {
         chain->Add(fnames[iadd].c_str());
         cout << "\tadded file '" << fnames[iadd].c_str() << "' to chain." << endl;
    }
    tree = chain;
        
    //
    //  Number of entries and events to process
    //
    double nEntries = tree->GetEntries();        
    cout << "Available Entries = " << nEntries << endl;
    nEntries *= fractionOfEventsToRead;
    cout << "Used Entries = " << nEntries << endl;
        
    //
    //  Set up tree
    //
    rootSartreEvent myEvent;
    TLorentzVector *eIn = 0;   // e beam
    TLorentzVector *pIn = 0;   // hadron beam
    TLorentzVector *vm = 0;    // vector meson
    TLorentzVector *vmd1 = 0;  // VM daughter 1
    TLorentzVector *vmd2 = 0;  // daughter 2
    TLorentzVector *eOut = 0;  // scattered electron
    TLorentzVector *pOut = 0;  // outgoing hadron
    TLorentzVector *gamma = 0; // virtual gamma
        
    tree->SetBranchAddress("event",&myEvent);
    tree->SetBranchAddress("eIn",&eIn);
    tree->SetBranchAddress("pIn",&pIn);
    tree->SetBranchAddress("vm",&vm);
    tree->SetBranchAddress("vmDaughter1",&vmd1);
    tree->SetBranchAddress("vmDaughter2",&vmd2);
    tree->SetBranchAddress("eOut",&eOut);
    tree->SetBranchAddress("pOut",&pOut);
    tree->SetBranchAddress("gamma",&gamma);
    
    tree->GetBranch("event");
    tree->GetBranch("eIn");
    tree->GetBranch("pIn");
    tree->GetBranch("vm");
    tree->GetBranch("vmDaughter1");
    tree->GetBranch("vmDaughter2");
    tree->GetBranch("eOut");
    tree->GetBranch("pOut");
    tree->GetBranch("gamma");
        
    double mkaon = 0.494;
    double mpion = 0.139;
    double mrho=0.0;
    double tofRes = 25.0;//picoseconds 
    double startRes = 30.0;//scattered electron timing picoseconds 
    double vmdpt = 0.0; 
    //
    //  Loop events
    //
    int oldinfo = gErrorIgnoreLevel;
    gErrorIgnoreLevel = kError;
    unsigned int acceptedEvents = 0;
    for(int ievent=0; ievent< nEntries; ievent++) {
      //    for(int ievent=0; ievent< 10000000; ievent++) {
        double eFrac = 100.*static_cast<double>(ievent)/nEntries;
        if (ievent%100000 == 0) cout << "Processing event " << ievent << " (" << eFrac << "%)" << endl;
        tree->GetEntry(ievent);
        TLorentzVector eInVec(*eIn);
        TLorentzVector pInVec(*pIn);
        TLorentzVector vmVec(*vm);
        TLorentzVector vmd1Vec(*vmd1);
        TLorentzVector vmd2Vec(*vmd2);
        TLorentzVector eOutVec(*eOut);
        TLorentzVector pOutVec(*pOut);
        TLorentzVector gammaVec(*gamma);
        
        //=================================================================
        //  ==> At this point all information of the tuple is available <==
        //=================================================================
        
        // Apply cuts (example) ...
        /*
        mrho= RhoSpectral();
        vmVec.SetVectM(vmVec.Vect(),mrho);
        vector<TLorentzVector> pions;
        pions = twoBodyDecay(vmVec,mpion,mpion);
        TLorentzVector vmd1Vec(pions[0]);
        TLorentzVector vmd2Vec(pions[1]);
        */
        mrho= vmVec.M();
        double tof1 = pathLength(vmd1Vec.Vect().Perp(),vmd1Vec.Vect().Mag())/vmd1Vec.Beta()*1000.0/30.0;//picoseconds 
        double tof2 = pathLength(vmd2Vec.Vect().Perp(),vmd2Vec.Vect().Mag())/vmd2Vec.Beta()*1000.0/30.0;//picoseconds 

        vmd1Vec.SetVectM(vmd1Vec.Vect(),mkaon);
        vmd2Vec.SetVectM(vmd2Vec.Vect(),mkaon);

        vmVec = vmd1Vec+vmd2Vec;
        bool accepted = true;
        //      hist_mass_coherent->Fill(vmVec.M());
        if (fabs(vmd1Vec.PseudoRapidity()) >1.0) accepted = false;
        if (fabs(vmd2Vec.PseudoRapidity()) > 1.0) accepted = false;
        hist_mass1_coherent->Fill(mrho);
        hist_mass2_coherent->Fill(vmVec.M());
        if(fabs(vmVec.M()-1.019)>0.020) accepted = false;
        if (!accepted) continue;
        acceptedEvents++;
        //      cout<<mrho<<"  "<<vmVec.M()<<"  "<<acceptedEvents<<endl;
        double sqrt2t = fabs(myEvent.t);
        //      mrho = pow(sqrt2t+mrho*mrho,2);
        //        double tsmear1 = fabs(gRandom->Gaus(sqrt2t,sqrt2t*sqrt2t/mrho));
        //      mrho = pow(sqrt2t*sqrt2t+0.776*0.776,2);
        //        double tsmear2 = fabs(gRandom->Gaus(sqrt2t,sqrt2t*sqrt2t/mrho));
        //              mrho = pow(sqrt2t*sqrt2t+1.019*1.019,2);
        //        double tsmear3 = fabs(gRandom->Gaus(sqrt2t,sqrt2t*sqrt2t/mrho));
        
        //
        // Histogramming (example) ...
        //

        double tof3 = pathLength(vmd1Vec.Vect().Perp(),vmd1Vec.Vect().Mag())/vmd1Vec.Beta()*1000.0/30.0;//picoseconds 
        double tof4 = pathLength(vmd2Vec.Vect().Perp(),vmd2Vec.Vect().Mag())/vmd2Vec.Beta()*1000.0/30.0;//picoseconds 
        double starttiming = gRandom->Gaus(0.0,startRes);
        double timesmear1 =  (gRandom->Gaus(0.0,tofRes)+starttiming);
        double timesmear2 =  (gRandom->Gaus(0.0,tofRes)+starttiming);
        tof1 += timesmear1;//picoseconds 
        tof2 += timesmear2;//picoseconds 

        double chi2 = 1.0/pow(tofRes,2)*(pow(tof1-tof3,2)+pow(tof2-tof4,2)-1.0/pow(tofRes,2)*pow(tof1+tof2-tof3-tof4,2)/(2.0/pow(tofRes,2)+1.0/pow(startRes,2)));
        double Normchi2 = 1.0/pow(tofRes,2)*(pow(timesmear1,2)+pow(timesmear2,2)-1.0/pow(tofRes,2)*pow(timesmear1+timesmear2,2)/(2.0/pow(tofRes,2)+1.0/pow(startRes,2)));
        if (vmd1Vec.Vect().Mag()>vmd2Vec.Vect().Mag())
          {
            vmdpt = vmd1Vec.Vect().Mag();
          }
        else vmdpt = vmd2Vec.Vect().Mag();

        if (myEvent.dmode < 0.5) { // coherent
            hist_t_coherent->Fill(fabs(myEvent.t), 1);
            hist_PIDChi2_coherent->Fill(vmdpt,chi2);
            hist_PIDNormChi2_coherent->Fill(vmdpt,Normchi2);
            if(chi2<4.6)
            hist_tPID_coherent->Fill(fabs(myEvent.t),1.);
            //      hist_tsmear1_coherent->Fill(tsmear1,1);
            //      hist_tsmear2_coherent->Fill(tsmear2,1);
            //      hist_tsmear2_coherent->Fill(tsmear2,1);
        }
        else {                    // incoherent
            hist_t_incoherent->Fill(fabs(myEvent.t), 1);
        }
        
        //
        //   Some print outs
        //
        /*
        cout << "\nEvent data:" << endl;
        PR(myEvent.t);
        PR(myEvent.Q2);
        PR(myEvent.x);
        PR(vmVec);
        PR(eOutVec);    */        
    } // end event loop
    
    
    int nAcceptedPercent = nEntries > 0 ? 100.*acceptedEvents/nEntries : 0;
    cout << "Accepted events: " << acceptedEvents << " (" << nEntries << ") or " << nAcceptedPercent << "%" << endl;
    gErrorIgnoreLevel = oldinfo;
    
    if (f) f->Close();
            
    //
    //  Save histos
    //
    TFile hfile("exclusiveRho2phiSpectralPIDProposal.root","RECREATE");
    //rho total cross section from log file is 35.8 ub while phi cross section is 4.72ub, a factor of 7.6 
    hist_t_coherent->Write();
    hist_t_incoherent->Write();
    hist_tPID_coherent->Write();
    hist_PIDChi2_coherent->Write();
    hist_PIDNormChi2_coherent->Write();
    hist_mass1_coherent->Write();
    hist_mass2_coherent->Write();
    //    hist_tsmear1_coherent->Write();
    //    hist_tsmear2_coherent->Write();
    //    hist_tsmear3_coherent->Write();
    hfile.Close();
    cout << "All histos stored in file 'exampleHistosRho2phiSpectral.root'." << endl;
}