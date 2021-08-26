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
    TH1D* h_pz_diff_b = new TH1D("h_pz_diff_b","",100,-10,10);
    TH1D* h_E_diff_b = new TH1D("h_E_diff_b","",100,-10,10);
    TH1D* h_pz_diff = new TH1D("h_pz_diff","",100,-10,10);
    TH1D* h_E_diff = new TH1D("h_E_diff","",100,-10,10);

    TH1D* h_new_t_A = new TH1D("h_new_t_A", "", 50, 0.0, 1.0);
    TH1D* h_new_t_D = new TH1D("h_new_t_D", "", 50, 0.0, 1.0);
    
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

        TLorentzVector dIn(0.,0.,200.,sqrt(200*200+MASS_DEUTERON*MASS_DEUTERON));
        
        TVector3 d_rf = dIn.BoostVector();
        nIn_d.Boost(d_rf);

        TVector3 p_rf = pIn.BoostVector();
        pIn_d.Boost(p_rf);
        TVector3 p_rf_new = pIn_d.BoostVector();

        pOut.Boost(-p_rf);
        pOut.Boost(p_rf_new);
        gammaOut.Boost(-p_rf);
        gammaOut.Boost(p_rf_new);

        TLorentzVector all = eIn+dIn-eOut-gammaOut-pOut-nIn_d;
        h_pz_diff_b->Fill(all.Pz());
        h_E_diff_b->Fill(all.E());
        
        //correcting
        gammaStar.Boost(-d_rf);
        nIn_d.Boost(-d_rf);
        gammaOut.Boost(-d_rf);
        pOut.Boost(-d_rf);

        double qzkz = gammaStar.Pz() - (nIn_d.Pz());//qz-kz
        double numn = gammaStar.E() - nIn_d.E();//sqrt( MASS_NEUTRON*MASS_NEUTRON + pxf*pxf+pyf*pyf+pzf*pzf )
        double jx = gammaOut.Px();
        double jy = gammaOut.Py();
        double jz = gammaOut.Pz();
        double px = pOut.Px();
        double py = pOut.Py();
        double pz = pOut.Pz();

        //ad hoc momentum conservation. move excess momentum energy to photon and struck nucleon
        jz = getCorrJz(qzkz,numn,jx,jy,px,py,MASS_PROTON);
        pz = getCorrPz(qzkz,numn,jx,jy,px,py,MASS_PROTON);

        gammaOut.SetPxPyPzE(jx,jy,jz,sqrt(jx*jx+jy*jy+jz*jz));
        pOut.SetPxPyPzE(px,py,pz,sqrt(px*px+py*py+pz*pz+MASS_PROTON*MASS_PROTON));

        gammaStar.Boost(d_rf);
        nIn_d.Boost(d_rf);
        gammaOut.Boost(d_rf);
        pOut.Boost(d_rf);

        all = eIn+dIn-eOut-gammaOut-pOut-nIn_d;

        h_pz_diff->Fill(all.Pz());
        h_E_diff->Fill(all.E());

        //Method A.
        TVector2 sum_pt(eOut.Px()+gammaOut.Px(), eOut.Py()+gammaOut.Py());
        h_new_t_A->Fill(sum_pt.Mod2());

        nIn_d.Boost(-d_rf);
        pOut.Boost(-d_rf);

        Double_t E_bInt = (alpha_AN*MASS_DEUTERON)/4. + (nIn_d.Px()*nIn_d.Px()+
            nIn_d.Py()*nIn_d.Py()+MASS_PROTON*MASS_PROTON)/(alpha_AN*MASS_DEUTERON);
        Double_t Pz_bInt = -(alpha_AN*MASS_DEUTERON)/4. + (nIn_d.Px()*nIn_d.Px()+
            nIn_d.Py()*nIn_d.Py()+MASS_PROTON*MASS_PROTON)/(alpha_AN*MASS_DEUTERON);
        //new 4 vector for struck nucleon before interaction;
        TLorentzVector n_primeprime;
        n_primeprime.SetPxPyPzE(-nIn_d.Px(),-nIn_d.Py(),
        Pz_bInt,E_bInt);
        double t_doubletagging = -1*(pOut - n_primeprime).Mag2();
        h_new_t_D->Fill(t_doubletagging);

        nIn_d.Boost(d_rf);
        pOut.Boost(d_rf);

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
