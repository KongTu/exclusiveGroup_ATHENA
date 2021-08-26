#include "HepMC3/GenEvent.h"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/Print.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/GenParticle.h"

#include <TH1.h>
#include <TF1.h>
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

    TF1 *deutNk_beagle = new TF1("Deuteron n(k) in fm^{-1}",getdNdkDeut,0,10,0);
    TF1* cthetaFlat= new TF1("cthetaFlat","0.5",-1.,1.);
    TF1* phiFlat= new TF1("phiFlat","1",-PI,PI);


    TFile* output = new TFile("output.root","RECREATE");

    h_Q2[0] = new TH1D("h_Q2_00", "", 90, 1., 10.);
    h_t[0] = new TH1D("h_t_00", "", 50, 0.0, 1.0);
    h_xB[0] = new TH1D("h_xB_00", "", 50, 0.0, 1.0);
    h_y[0] = new TH1D("h_y_00", "", 50, 0.0, 1.0);
    h_phi[0] = new TH1D("h_phi_00", "", 100, 0.0, 6.2831);

    TH1D* h_alpha = new TH1D("h_alpha","",100,0,2);
    TH1D* h_p = new TH1D("h_p","",100,0,1);
    TH1D* h_pz = new TH1D("h_pz","",100,-1,1);

    
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
        //out electron
        TLorentzVector eOut = getFourMomentum(evt.particles().at(1)); 
        //virtual photon
        TLorentzVector gammaStar = getFourMomentum(evt.particles().at(2)); 
        //in proton
        TLorentzVector pIn = getFourMomentum(evt.particles().at(3)); 
        //out photon
        TLorentzVector gammaOut = getFourMomentum(evt.particles().at(4)); 
        //out proton
        TLorentzVector pOut = getFourMomentum(evt.particles().at(5)); 
   
        //deuteron light front wave fucntion:
        double k1 = deutNk_beagle->GetRandom()*0.197;
        double theta1=TMath::ACos(cthetaFlat->GetRandom());
        double phi1 = phiFlat->GetRandom();
        double kx1=k1*TMath::Sin(theta1)*TMath::Cos(phi1);
        double ky1=k1*TMath::Sin(theta1)*TMath::Sin(phi1);
        double kz1=k1*TMath::Cos(theta1);
        double Enn = sqrt(kx1*kx1+ky1*ky1+kz1*kz1+MASS_NUCLEON*MASS_NUCLEON);
        double alpha_SN = 1. - kz1 / Enn;
        double alpha_AN = 2 - alpha_SN;
        double ANMT2 = MASS_PROTON*MASS_PROTON+kx1*kx1+ky1*ky1;
        double pz1 = -(alpha_AN*MASS_DEUTERON)/4. + ANMT2/(alpha_AN*MASS_DEUTERON);
        double E1 = (alpha_AN*MASS_DEUTERON)/4. + ANMT2/(alpha_AN*MASS_DEUTERON);
        double SNMT2 = MASS_NEUTRON*MASS_NEUTRON + kx1*kx1 + ky1*ky1;
        double pz2 = -(alpha_SN*MASS_DEUTERON)/4. + SNMT2/(alpha_SN*MASS_DEUTERON);
        double E2 = (alpha_SN*MASS_DEUTERON)/4. + SNMT2/(alpha_SN*MASS_DEUTERON);

        h_alpha->Fill(alpha_SN);
        TLorentzVector pIn_d(kx1,ky1,pz1,E1);
        TLorentzVector nIn_d(-kx1,-ky1,pz2,E2);
        h_p->Fill(nIn_d.P());
        h_pz->Fill(nIn_d.Pz());

        //begin boost:
        cout << "starting to boost around. Step.1 boost to CM frame" << endl;
        TLorentzVector cm = eIn+pIn;
        TVector3 cm_boost = cm.BoostVector();
        gammaOut.Boost(-cm_boost);
        eOut.Boost(-cm_boost);
        pOut.Boost(-cm_boost);

        PRINT4VECTOR(cm,1);

        TVector3 restframe = pIn.BoostVector();
        cout << "new cm frame ~ " << endl;
        pIn_d.Boost(restframe);
        TLorentzVector cm_new = eIn+pIn_d;
        TVector3 cm_new_boost = cm_new.BoostVector();
        gammaOut.Boost(cm_new_boost);
        eOut.Boost(cm_new_boost);
        pOut.Boost(cm_new_boost);

        PRINT4VECTOR(cm_new,1);
        cout << "conservation ~ " << endl;

        TLorentzVector all = eIn+pIn_d-eOut-pOut-gammaOut;
        PRINT4VECTOR(all,1);
        cout << "all particles " << endl;
        PRINT4VECTOR(eIn,1);
        PRINT4VECTOR(eOut,1);
        PRINT4VECTOR(pIn_d,1);
        PRINT4VECTOR(gammaOut,1);
        PRINT4VECTOR(pOut,1);


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
    
    output->Write();
    output->Close();

    return 0;
}
