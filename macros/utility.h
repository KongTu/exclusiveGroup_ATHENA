#include "RiceStyle.h"
#include "TGaxis.h"
using namespace std;
#define PI 3.1415926
#define MASS_PROTON   0.93827208816
#define MASS_NEUTRON  0.93956542052
#define MASS_NUCLEON  0.93891875 //.93891875
#define MASS_DEUTERON  1.8756129 //1.8756129 //1.8756134 (most precise)

int vm_index=-1;
TString legendName="";
void setVM(TString name="rho"){
	if(name=="rho"||name=="rho_photo"){
		vm_index=0;
		legendName="#rho^{0}";
	}
	else if(name=="phi"||name=="phi_photo"){
		vm_index=1;
		legendName="#phi";
	}
	else if(name=="jpsi"||name=="jpsi_photo"){
		vm_index=2;
		legendName="J/#psi";
	}
}
double sigma_sartre_elect[]={3.58E+4,4.72E+3,199.};
double sigma_sartre_photo[]{3.86E+5,2.42E+5,458.};

void measureXsection(TString name="rho", TH1D* hist=0, int sample=0, int total_events=1e5, bool PHP_=false, bool cutOnDaug_=true){

	//branching ratios
	double BR_beagle_decay[] = {1.0,0.489,0.5};//branching ratio, depend on what we select in beagle
	if(!cutOnDaug_) {
		BR_beagle_decay[0]=1.0;BR_beagle_decay[1]={1.0};BR_beagle_decay[2]={1.0};
	}
	double BR_sartre_decay[] = {1.0,1.0,1.0};//branching ratio
	//beagle constants
	double beagle_lumi = 1e5/(34.4*197);//nanobarn
	if(PHP_) beagle_lumi = 1e5/(1.971E3*197);//nanobarn
	double lumi_factor = total_events/1e5;
	beagle_lumi=beagle_lumi*lumi_factor;

	double sigma=-1;
	if(name=="rho"){
		sigma = 3.58E+4;
		vm_index=0;
		legendName="#rho^{0}";
	}
	else if(name=="phi"){
		sigma=4.72E+3;
		vm_index=1;
		legendName="#phi";
	}
	else if(name=="jpsi"){
		sigma=199.;
		vm_index=2;
		legendName="J/#psi";
	}
	//photoproduction
	else if(name=="rho_photo"){
		sigma = 3.86E+5;
		vm_index=0;
		legendName="#rho^{0}";
	}
	else if(name=="phi_photo"){
		sigma=2.42E+5;
		vm_index=1;
		legendName="#phi";
	}
	else if(name=="jpsi_photo"){
		sigma=458.;
		vm_index=2;
		legendName="J/#psi";
	}
	//here we use 2E7 events
	double sartre_lumi = 20000000./sigma;//nanbarn

	double BR_to_apply = 1.;
	double lumi_to_apply=-99;
	if(sample==0){
		lumi_to_apply = beagle_lumi;
		BR_to_apply = BR_beagle_decay[vm_index];
	}
	else if(sample==1){
		lumi_to_apply = sartre_lumi;
		BR_to_apply = BR_sartre_decay[vm_index];
	}
	double binwidth = hist->GetBinWidth(1);

	hist->Scale(1./ (lumi_to_apply * BR_to_apply * binwidth) );

}