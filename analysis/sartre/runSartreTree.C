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
//kong's include
#include "../include/pleaseIncludeMe.h"
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

void runSartreTree(double fractionOfEventsToRead = 1, TString vm_name="jpsi", int PID_=0, double setLowPt_=0.1, bool smear_=false)
{
    minPt_ = setLowPt_;
    TString name_PID=Form("%d",PID_);
    TString name_LowPt=Form("%.2f",setLowPt_);
    TString name_smear=Form("%d",smear_);

    //  Setup filenames
    //
    string fnames[16];
    int nnames = 0;

    fnames[0] = "/gpfs02/eic/DATA/sartre/data/sartre_bnonsat_Au_"+vm_name+"_1.root";
    fnames[1] = "/gpfs02/eic/DATA/sartre/data/sartre_bnonsat_Au_"+vm_name+"_2.root";
    fnames[2] = "/gpfs02/eic/DATA/sartre/data/sartre_bnonsat_Au_"+vm_name+"_3.root";
    fnames[3] = "/gpfs02/eic/DATA/sartre/data/sartre_bnonsat_Au_"+vm_name+"_4.root";
    fnames[4] = "/gpfs02/eic/DATA/sartre/data/sartre_bnonsat_Au_"+vm_name+"_5.root";
    fnames[5] = "/gpfs02/eic/DATA/sartre/data/sartre_bnonsat_Au_"+vm_name+"_6.root";
    fnames[6] = "/gpfs02/eic/DATA/sartre/data/sartre_bnonsat_Au_"+vm_name+"_7.root";
    fnames[7] = "/gpfs02/eic/DATA/sartre/data/sartre_bnonsat_Au_"+vm_name+"_8.root";
    fnames[8] = "/gpfs02/eic/DATA/sartre/data/sartre_bnonsat_Au_"+vm_name+"_9.root";
    fnames[9] = "/gpfs02/eic/DATA/sartre/data/sartre_bnonsat_Au_"+vm_name+"_10.root";
    nnames = 10;
    
    //output root files:
    TFile hfile("../../rootfiles/sartre_"+vm_name+"_bnonsat_PID_"+name_PID+"_minPt_"+name_LowPt+"_smear_"+name_smear+".root","RECREATE");
    //
    //   Histogram Booking (example)
    //
    TH1D *hist_t_coherent = new TH1D("hist_t_coherent", "coherent", 72, 0, 0.18);
    TH1D *hist_t_incoherent = new TH1D("hist_t_incoherent", "incoherent", 72, 0, 0.18);
    TH1D *hist_t_afterPhaseSpace_coherent = new TH1D("hist_t_afterPhaseSpace_coherent", "coherent", 72, 0, 0.18);
    TH1D *hist_t_afterPhaseSpace_incoherent = new TH1D("hist_t_afterPhaseSpace_incoherent", "incoherent", 72, 0, 0.18);

    TH1D *h_VM[2][3];
    TH1D *h_VM_daughter[2][3];
    double bin_lower[]={0.,-8.,0.,0.,0.};
    double bin_upper[]={5.0,8.,6.5,4.,0.2};
    for(int icoh=0;icoh<2;icoh++){
        for(int ipro=0;ipro<3;ipro++){
            h_VM[icoh][ipro]=new TH1D(Form("h_VM_%d_%d",icoh,ipro),
            Form("h_VM_%d_%d",icoh,ipro),100,bin_lower[ipro],bin_upper[ipro] );
            h_VM_daughter[icoh][ipro]=new TH1D(Form("h_VM_daughter_%d_%d",icoh,ipro),
            Form("h_VM_daughter_%d_%d",icoh,ipro),100,bin_lower[ipro],bin_upper[ipro] );
        }
    }
    // t_reco histograms
    TH1D* h_t_reco[2][3][3];
    TH2D* h_VM_t_mass[2][3][3];
    for(int ibreak=0;ibreak<2;ibreak++){
        for(int imethod=0;imethod<3;imethod++){
            for(int imass=0;imass<3;imass++){
                h_t_reco[ibreak][imethod][imass] = new TH1D(Form("h_t_reco_%d_%d_%d",ibreak,imethod,imass),
                    Form("h_t_reco_%d_%d_%d",ibreak,imethod,imass),1000,0,2 );
                h_VM_t_mass[ibreak][imethod][imass] = new TH2D(Form("h_VM_t_mass_%d_%d_%d",ibreak,imethod,imass),
                        Form("h_VM_t_mass_%d_%d_%d",ibreak,imethod,imass),1000,0.,4,100,0.,0.2);
            }
        }
    }
    // vm mass from daughters
    TH1D* h_VM_mass[2][3];
    //nuclear remnant mass
    TH2D* h_Amass[2][3];
    for(int ibreak=0;ibreak<2;ibreak++){
        for(int imass=0;imass<3;imass++){
            h_VM_mass[ibreak][imass] = new TH1D(Form("h_VM_mass_%d_%d",ibreak,imass),
                Form("h_VM_mass_%d_%d",ibreak,imass),1000,0.,4);
            h_Amass[ibreak][imass] = new TH2D(Form("h_Amass_%d_%d",ibreak,imass),
                Form("h_Amass_%d_%d",ibreak,imass), 100,0,0.2,100,-3,3);
        }
    }
    TH2D* h_PID=new TH2D("h_PID",";p;chi2",100,0,3,500,0,100);

    
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
        TLorentzVector aInVec(0.,0.,pInVec.Pz()*197, sqrt(pInVec.Pz()*197*pInVec.Pz()*197 + MASS_AU197*MASS_AU197) );
        
        
        //=================================================================
        //  ==> At this point all information of the tuple is available <==
        //=================================================================
       
         //Kong's edits to fill some histograms:
        int coh_index=-1;
        if(myEvent.dmode < 0.5) coh_index=0;
        else coh_index=1;
        //VM.
        h_VM[coh_index][0]->Fill(vmVec.Pt());
        h_VM[coh_index][1]->Fill(vmVec.Eta());
        h_VM[coh_index][2]->Fill(vmVec.Phi());
         //daug.1
        h_VM_daughter[coh_index][0]->Fill(vmd1Vec.Pt());
        h_VM_daughter[coh_index][1]->Fill(vmd1Vec.Eta());
        h_VM_daughter[coh_index][2]->Fill(vmd1Vec.Phi());
        //daug.2
        h_VM_daughter[coh_index][0]->Fill(vmd2Vec.Pt());
        h_VM_daughter[coh_index][1]->Fill(vmd2Vec.Eta());
        h_VM_daughter[coh_index][2]->Fill(vmd2Vec.Phi());
        

        bool accepted = true;
        if (TMath::Abs(vmVec.Rapidity())>4.) accepted = false;
        if (!accepted) continue;
        if (myEvent.dmode < 0.5) { // coherent
            hist_t_coherent->Fill(fabs(myEvent.t), 1);
        }
        else {                    // incoherent
            hist_t_incoherent->Fill(fabs(myEvent.t), 1);
        }
        if (TMath::Abs(vmd1Vec.PseudoRapidity()) > 1.) accepted = false;
        if (TMath::Abs(vmd2Vec.PseudoRapidity()) > 1.) accepted = false;
        if (vmd1Vec.Pt() < minPt_) accepted = false;
        if (vmd2Vec.Pt() < minPt_) accepted = false;
        if (!accepted) continue;
        if (myEvent.dmode < 0.5) { // coherent
            hist_t_afterPhaseSpace_coherent->Fill(fabs(myEvent.t), 1);
        }
        else {                    // incoherent
            hist_t_afterPhaseSpace_incoherent->Fill(fabs(myEvent.t), 1);
        }
        
        acceptedEvents++;
        //smearing
        vector< TLorentzVector> update;
        if(smear_){
            update = letsMakeItReal(eInVec,eOutVec,aInVec,vmd1Vec,vmd2Vec);
            eInVec=update[0];eOutVec=update[1];aInVec=update[2];vmd1Vec=update[3];vmd2Vec=update[4];  
        }
        //VM t
        for(int imass=0;imass<3;imass++){
            TLorentzVector vmd1Vec_new,vmd2Vec_new,vmVec_new;
            TVector3 temp_v1=vmd1Vec.Vect();
            TVector3 temp_v2=vmd2Vec.Vect();
            vmd1Vec_new.SetVectM(temp_v1,daughtermasslist[imass]);
            vmd2Vec_new.SetVectM(temp_v2,daughtermasslist[imass]);
            vmVec_new = vmd1Vec_new+vmd2Vec_new;
            
            if(PID_==1){
                double chi2=-99.;
                if(TMath::Abs(vmd1Vec_new.Eta())<1.0 
                            && TMath::Abs(vmd2Vec_new.Eta())<1.0
                                && imass==1){
                    if( vm_name=="rho"||vm_name=="rho_photo" ){
                        chi2 = giveMe_PIDChi2(vmd1Vec_new, vmd2Vec_new, MASS_PION);
                    }
                    if( vm_name=="phi"||vm_name=="phi_photo" ){
                        chi2 = giveMe_PIDChi2(vmd1Vec_new, vmd2Vec_new, MASS_KAON);
                    }
                    h_PID->Fill(vmd1Vec_new.P(), chi2);
                    if( chi2>4.6 ) continue;
                }
            }

            double mass = vmVec_new.M();
            for(int imethod=0;imethod<3;imethod++){
                double t_reco = giveMe_t(imethod,eInVec,eOutVec,aInVec,vmVec_new);
                h_t_reco[coh_index][imethod][imass]->Fill( t_reco );
                h_VM_t_mass[coh_index][imethod][imass]->Fill( mass, t_reco );
                if(imethod==2)h_Amass[coh_index][imass]->Fill(t_reco,giveMe_Amass(eInVec,eOutVec,aInVec,vmVec_new));
            }
            //invMASS
            h_VM_mass[coh_index][imass]->Fill( mass );
        }
       

        // Apply cuts (example) ..
    
        
        
        //
        // Histogramming (example) ...
        //
        
        
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
    
    hfile.Write();
    hfile.Close();
    cout << "All histos stored in file '../../rootfiles/sartre_"+vm_name+"_bnonsat_PID_"+name_PID+"_minPt_"+name_LowPt+"_smear_"+name_smear+".root'." << endl;
}
