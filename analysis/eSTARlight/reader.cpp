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

#include "TFile.h"
#include "TLorentzVector.h"

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

        //virtual photon
        TLorentzVector gammaStar = getFourMomentum(evt.particles().at(0)); 

        //in kaons
        TLorentzVector kaonPlus = getFourMomentum(evt.particles().at(1)); 
        TLorentzVector kaonMinus = getFourMomentum(evt.particles().at(2)); 
        TLorentzVector phi = kaonPlus+kaonMinus;
        h_phi_mass->Fill( phi.M() );

        //out electron
        TLorentzVector eOut = getFourMomentum(evt.particles().at(3)); 

        //out proton
        TLorentzVector pOut = getFourMomentum(evt.particles().at(4)); 

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
