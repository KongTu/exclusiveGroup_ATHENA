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
void measureXsection(TString name="rho", TH1D* hist=0, int sample=0){

	//branching ratios
	double BR_decay[] = {1.0,0.489,1.0};//branching ratio
	//beagle constants
	double beagle_lumi = 1e5/(34.4*197);//nanobarn
	beagle_lumi=beagle_lumi*100.;

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
	//here we use 2E7 events
	double sartre_lumi = 20000000./sigma;//nanbarn

	double lumi_to_apply=-99;
	if(sample==0){
		lumi_to_apply = beagle_lumi;
	}
	else if(sample==1){
		lumi_to_apply = sartre_lumi;
	}
	double binwidth = hist->GetBinWidth(1);

	hist->Scale(1./ (lumi_to_apply * BR_decay[vm_index] * binwidth) );

}