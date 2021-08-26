#include "HepMC3/GenEvent.h"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/Print.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/GenParticle.h"

#include <TH1.h>
#include <TROOT.h>
#include <TFile.h>
#include <TGraph.h>
#include <TCanvas.h>

#include "DVCSEvent.h"

#include <iostream>

using namespace HepMC3;

//main
int main(int argc, char **argv) {

    //input file
    if(argc != 2){

        std::cout << __func__ << ": error: usage should be: " << argv[0] << " inputFile" << std::endl;
        exit(0);
    }

    //histograms
    std::vector<TH1D*> h_Q2(10, nullptr);
    std::vector<TH1D*> h_t(10, nullptr);
    std::vector<TH1D*> h_xB(10, nullptr);
    std::vector<TH1D*> h_y(10, nullptr);
    std::vector<TH1D*> h_phi(10, nullptr);

    h_Q2[0] = new TH1D("h_Q2_00", "", 90, 1., 10.);
    h_t[0] = new TH1D("h_t_00", "", 50, 0.0, 1.0);
    h_xB[0] = new TH1D("h_xB_00", "", 50, 0.0, 1.0);
    h_y[0] = new TH1D("h_y_00", "", 50, 0.0, 1.0);
    h_phi[0] = new TH1D("h_phi_00", "", 100, 0.0, 6.2831);
    
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

        //DVCS event
        DVCSEvent dvcsEvent(evt);

        TLorentzVector eIn = getFourMomentum(evt.particles().at(0)); 
        cout << "eIn: " << eIn.Px() << " " << eIn.Py() << " " << eIn.Pz() << " " << eIn.E() << endl;

        //fill
        h_Q2[0]->Fill(dvcsEvent.getQ2());
        h_t[0]->Fill(-1 * dvcsEvent.getT());
        h_xB[0]->Fill(dvcsEvent.getXb());
        h_y[0]->Fill(dvcsEvent.getY());
        h_phi[0]->Fill(dvcsEvent.getPhi());

        //id
        iEvent++;
    }
    
    //run info
    std::shared_ptr<HepMC3::GenRunInfo> runInfo = inputFile.run_info();
    
    //run info
    std::shared_ptr<HepMC3::Attribute> att_csValue = runInfo->attributes().find("integrated_cross_section_value")->second;
    if(att_csValue){
    	double csValue = std::stod(att_csValue->unparsed_string());
    	std::cout << csValue << std::endl; 
    }else{
	std::cout << "HepMC3::Attribute pointer is NULL" << std::endl;
    }
    
    //run info
    std::shared_ptr<HepMC3::Attribute> att_eventsNumber = runInfo->attributes().find("generated_events_number")->second;
    if(att_eventsNumber){
        double eventsNumber = std::stod(att_eventsNumber->unparsed_string());
        std::cout << eventsNumber << std::endl;
    }else{
	std::cout << "HepMC3::Attribute pointer is NULL" << std::endl;
    }
    
    //close file
    inputFile.close();

    //print 
   
    TCanvas* can = new TCanvas();
    can->Divide(3, 2);

    can->cd(1);
    h_Q2[0]->Draw();
    can->cd(2);
    h_t[0]->Draw();
    can->cd(3);
    h_xB[0]->Draw();
    can->cd(4);
    h_y[0]->Draw();
    can->cd(5);
    h_phi[0]->Draw();
  
    can->Print("plots.pdf", "pdf");
    
    return 0;
}