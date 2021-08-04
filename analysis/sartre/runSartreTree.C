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
#define PR(x) cout << #x << " = " << (x) << endl;

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

//===========================================================================
//  Main function
//===========================================================================

void runSartreTree(double fractionOfEventsToRead = 1)
{
    //
    //  Setup filenames
    //
    string fnames[16];
    int nnames = 0;

    fnames[0] = "/gpfs02/eic/DATA/sartre/data/sartre_bnonsat_Au_jpsi_1.root";
    fnames[1] = "/gpfs02/eic/DATA/sartre/data/sartre_bnonsat_Au_jpsi_2.root";
    fnames[2] = "/gpfs02/eic/DATA/sartre/data/sartre_bnonsat_Au_jpsi_3.root";
    fnames[3] = "/gpfs02/eic/DATA/sartre/data/sartre_bnonsat_Au_jpsi_4.root";
    fnames[4] = "/gpfs02/eic/DATA/sartre/data/sartre_bnonsat_Au_jpsi_5.root";
    fnames[5] = "/gpfs02/eic/DATA/sartre/data/sartre_bnonsat_Au_jpsi_6.root";
    fnames[6] = "/gpfs02/eic/DATA/sartre/data/sartre_bnonsat_Au_jpsi_7.root";
    fnames[7] = "/gpfs02/eic/DATA/sartre/data/sartre_bnonsat_Au_jpsi_8.root";
    fnames[8] = "/gpfs02/eic/DATA/sartre/data/sartre_bnonsat_Au_jpsi_9.root";
    fnames[9] = "/gpfs02/eic/DATA/sartre/data/sartre_bnonsat_Au_jpsi_10.root";
    nnames = 10;
   
    //
    //   Histogram Booking (example)
    //
    TH1D *hist_t_coherent = new TH1D("hist_t_coherent", "coherent", 72, 0, 0.18);
    TH1D *hist_t_incoherent = new TH1D("hist_t_incoherent", "incoherent", 72, 0, 0.18);

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
        
    
    //
    //  Loop events
    //
    int oldinfo = gErrorIgnoreLevel;
    gErrorIgnoreLevel = kError;
    unsigned int acceptedEvents = 0;
    for(int ievent=0; ievent< nEntries; ievent++) {
        
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
    
        bool accepted = true;
        if (vmd1Vec.PseudoRapidity() < -4) accepted = false;
        if (vmd1Vec.PseudoRapidity() > 4) accepted = false;
        if (!accepted) continue;
	acceptedEvents++;
        
        //
        // Histogramming (example) ...
        //
        if (myEvent.dmode < 0.5) { // coherent
            hist_t_coherent->Fill(fabs(myEvent.t), 1);
        }
        else {                    // incoherent
            hist_t_incoherent->Fill(fabs(myEvent.t), 1);
        }
        
        //
        //   Some print outs
        //
        // cout << "\nEvent data:" << endl;
        // PR(myEvent.t);
        // PR(myEvent.Q2);
        // PR(myEvent.x);
        // PR(vmVec);
        // PR(eOutVec);
        
    } // end event loop
    
    
    int nAcceptedPercent = nEntries > 0 ? 100.*acceptedEvents/nEntries : 0;
    cout << "Accepted events: " << acceptedEvents << " (" << nEntries << ") or " << nAcceptedPercent << "%" << endl;
    gErrorIgnoreLevel = oldinfo;
    
    if (f) f->Close();
            
    //
    //  Save histos
    //
    TFile hfile("../../rootfiles/sartre_jpsi_bnonsat.root","RECREATE");
    hist_t_coherent->Write();
    hist_t_incoherent->Write();
    hfile.Close();
    cout << "All histos stored in file '../../rootfiles/sartre_jpsi_bnonsat.root'." << endl;
}
