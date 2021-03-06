#include "HepMC3/GenEvent.h"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/WriterAscii.h"
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

#include "DVCSEvent.h"

#include <iostream>

using namespace HepMC3;

//main
int main(int argc, char **argv) {

    //input file
    if(argc != 3){

        std::cout << __func__ << ": error: usage should be: " << argv[0] << " inputFile" << argv[1] << " outputFile" << std::endl;
        exit(0);
    }

    //histograms

    TF1 *deutNk_beagle = new TF1("Deuteron n(k) in fm^{-1}",getdNdkDeut,0,10,0);
    TF1* cthetaFlat= new TF1("cthetaFlat","0.5",-1.,1.);
    TF1* phiFlat= new TF1("phiFlat","1",-PI,PI);

    //open file
    ReaderAscii inputFile(argv[1]);
    WriterAscii text_output(argv[2]);
    
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

        //*Kong starts here.

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
        //swap proton w. neutron
        double ANMT2 = MASS_NEUTRON*MASS_NEUTRON+kx1*kx1+ky1*ky1;
        double pz1 = -(alpha_AN*MASS_DEUTERON)/4. + ANMT2/(alpha_AN*MASS_DEUTERON);
        double E1 = (alpha_AN*MASS_DEUTERON)/4. + ANMT2/(alpha_AN*MASS_DEUTERON);
        double SNMT2 = MASS_PROTON*MASS_PROTON + kx1*kx1 + ky1*ky1;
        double pz2 = -(alpha_SN*MASS_DEUTERON)/4. + SNMT2/(alpha_SN*MASS_DEUTERON);
        double E2 = (alpha_SN*MASS_DEUTERON)/4. + SNMT2/(alpha_SN*MASS_DEUTERON);

        //swap proton w. neutron
        TLorentzVector nIn_d(kx1,ky1,pz1,E1);
        TLorentzVector pIn_d(-kx1,-ky1,pz2,E2);
 
        //Deuteron beam.
        TLorentzVector dIn(0.,0.,pIn.Pz()*2,sqrt(pIn.Pz()*2*pIn.Pz()*2+MASS_DEUTERON*MASS_DEUTERON));
        
        TVector3 d_rf = dIn.BoostVector();
        nIn_d.Boost(d_rf);

        TVector3 p_rf = pIn.BoostVector();
        pIn_d.Boost(p_rf);
        TVector3 p_rf_new = pIn_d.BoostVector();

        pOut.Boost(-p_rf);
        pOut.Boost(p_rf_new);
        gammaOut.Boost(-p_rf);
        gammaOut.Boost(p_rf_new);

        TLorentzVector nOut;
        nOut.SetPtEtaPhiM(pOut.Pt(),pOut.Eta(),pOut.Phi(),MASS_NEUTRON);
        //swap proton w. neutron
        TLorentzVector all = eIn+dIn-eOut-gammaOut-nOut-pIn_d;
        
        //correcting energy momentum. 
        gammaStar.Boost(-d_rf);
        pIn_d.Boost(-d_rf);
        gammaOut.Boost(-d_rf);
        nOut.Boost(-d_rf);

        double qzkz = gammaStar.Pz() - (pIn_d.Pz());//qz-kz
        double numn = gammaStar.E() - pIn_d.E();//sqrt( MASS_NEUTRON*MASS_NEUTRON + pxf*pxf+pyf*pyf+pzf*pzf )
        double jx = gammaOut.Px();
        double jy = gammaOut.Py();
        double jz = gammaOut.Pz();
        double px = nOut.Px();
        double py = nOut.Py();
        double pz = nOut.Pz();

        //ad hoc momentum conservation. move excess momentum energy to photon and struck nucleon
        jz = getCorrJz(qzkz,numn,jx,jy,px,py,MASS_NEUTRON);
        pz = getCorrPz(qzkz,numn,jx,jy,px,py,MASS_NEUTRON);

        gammaOut.SetPxPyPzE(jx,jy,jz,sqrt(jx*jx+jy*jy+jz*jz));
        nOut.SetPxPyPzE(px,py,pz,sqrt(px*px+py*py+pz*pz+MASS_NEUTRON*MASS_NEUTRON));

        gammaStar.Boost(d_rf);
        pIn_d.Boost(d_rf);
        gammaOut.Boost(d_rf);
        nOut.Boost(d_rf);

        all = eIn+dIn-eOut-gammaOut-nOut-pIn_d;

        //Hepmc3 output

        GenEvent evt_w(Units::GEV,Units::MM);
        evt_w.set_event_number(iEvent);
        std::shared_ptr<GenCrossSection> cross_section = std::make_shared<GenCrossSection>();
        evt_w.add_attribute("GenCrossSection",cross_section);
        cross_section->set_cross_section(1.0,0.0);//same as input

        //                                                               px      py        pz       e     pdgid status
        GenParticlePtr p1 = std::make_shared<GenParticle>( FourVector( eIn.Px(), eIn.Py(),  eIn.Pz(),  eIn.E() ),11,  4 );
        GenParticlePtr p2 = std::make_shared<GenParticle>( FourVector( eOut.Px(), eOut.Py(),  eOut.Pz(),  eOut.E()),11,  1 );
        GenParticlePtr p3 = std::make_shared<GenParticle>( FourVector( gammaStar.Px(), gammaStar.Py(),  gammaStar.Pz(),  gammaStar.E() ),22,  3 );
        GenParticlePtr p4 = std::make_shared<GenParticle>( FourVector( dIn.Px(), dIn.Py(),  dIn.Pz(),  dIn.E()), 1000010020,  4 );
        
        GenVertexPtr v1 = std::make_shared<GenVertex>();
        v1->add_particle_in (p1);
        v1->add_particle_out (p2);
        v1->add_particle_out(p3);
        evt_w.add_vertex(v1);
        
        GenVertexPtr v2 = std::make_shared<GenVertex>();
        v2->add_particle_in(p3);
        v2->add_particle_in(p4);
        evt_w.add_vertex(v2);

        GenParticlePtr p5 = std::make_shared<GenParticle>( FourVector( gammaOut.Px(), gammaOut.Py(),  gammaOut.Pz(),  gammaOut.E() ),22,  1 );
        GenParticlePtr p6 = std::make_shared<GenParticle>( FourVector( nOut.Px(), nOut.Py(),  nOut.Pz(),  nOut.E()), 2112,  1 );
        GenParticlePtr p7 = std::make_shared<GenParticle>( FourVector( pIn_d.Px(), pIn_d.Py(),  pIn_d.Pz(),  pIn_d.E()), 2212,  1 );

        v2->add_particle_out(p5);
        v2->add_particle_out(p6);
        v2->add_particle_out(p7);
        
        // Print::listing(evt_w);
        text_output.write_event(evt_w);

        //id
        iEvent++;
    }
    
    std::cout << "Total events = " << iEvent << " has been processed and written in hepmc format" << std::endl;
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
    text_output.close();


    return 0;
}
