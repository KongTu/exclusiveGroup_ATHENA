#include "HepMC3/GenEvent.h"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/Print.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/GenParticle.h"

#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TROOT.h>
#include <TFile.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <iostream>
#include "TLorentzVector.h"

#define MASS_AU197    183.45406466643374
#define MASS_PROTON   0.93827
using namespace HepMC3;

TLorentzVector getFourMomentum(std::shared_ptr<const HepMC3::GenParticle> p){
    return TLorentzVector(p->momentum().px(), p->momentum().py(), p->momentum().pz(), p->momentum().e());
}

//main
int main(int argc, char **argv) {

    //input file
    if(argc != 2){

        std::cout << __func__ << ": error: usage should be: " << argv[0] << " inputFile" << std::endl;
        exit(0);
    }

    //histograms

    TFile* output = new TFile("output.root","RECREATE");
    TH1D* h_phi_mass = new TH1D("h_phi_mass",";mass (GeV)",100,0.1,1.5);
    TH2D* hQ2vsX = new TH2D("hQ2vsX",";xbj;Q2 (GeV^{2})",1000,1e-11,1e-6,10000,1e-11,1e-7);
    TH1D* hQ2 = new TH1D("hQ2",";Q2 (GeV^{2})",100000,0,1);
    //open file
    ReaderAscii inputFile(argv[1]);
    
    //loop over event
    size_t iEvent = 0;

    while(! inputFile.failed()) {

        //event
        GenEvent evt(Units::GEV,Units::MM);

        //read event from input file
        inputFile.read_event(evt);

        //if reading failed - exit loop
        if(inputFile.failed() ) break;

        //beam particles:
        TLorentzVector eIn = getFourMomentum(evt.particles().at(0));
        //virtual photon
        TLorentzVector gammaStar = getFourMomentum(evt.particles().at(1)); 
         //out electron
        TLorentzVector eOut = getFourMomentum(evt.particles().at(2));
        
        TLorentzVector pIn = getFourMomentum(evt.particles().at(3));

        //in kaons
        TLorentzVector kaonPlus = getFourMomentum(evt.particles().at(4)); 
        TLorentzVector kaonMinus = getFourMomentum(evt.particles().at(5)); 

        //out Ion
        TLorentzVector AOut = getFourMomentum(evt.particles().at(6)); 

        TLorentzVector phi = kaonPlus+kaonMinus;
        h_phi_mass->Fill( phi.M() );
        // if(fabs(kaonMinus.Eta())>4.0 
        //     || fabs(kaonPlus.Eta())>4.0
        //       || kaonMinus.Pt() < 0.15 
        //         || kaonPlus.Pt() < 0.15 ) continue;

        // if( fabs(phi.Rapidity()) > 4.0 ) continue;

        TLorentzVector q=eIn-eOut;
        double Q2= -q.Mag2();
        double Q2_starlight = gammaStar.Mag2();
        cout << "Q2 kong " << Q2 <<  " with Q2 from starlight " << Q2_starlight << endl;
        hQ2->Fill(Q2_starlight);
        double xbj = Q2_starlight / (2*pIn.Dot(gammaStar));
        hQ2vsX->Fill(xbj,Q2_starlight);

        //id
        iEvent++;
    }
    
    //close file
    inputFile.close();

    //print 

    output->Write();
    output->Close();

    return 0;
}
