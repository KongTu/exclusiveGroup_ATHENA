#include "TROOT.h"
void makeGausSmearTable(){

	double elec_reso=9.8E-4;
	double gold_reso=6.6E-4;
	
	for(int i=0;i<800;i++){
		//double value = gRandom->Gaus(0.0,elec_reso);
		double value2 = gRandom->Gaus(0.0,gold_reso);
		cout << 110.*(1.+value2) << endl;//"      " << 110.*(1.+value2)  << endl;
	}
	
	

}