/***********************************************************************
Toy MC spherical generator.

Description:
  
	Spherical decay of a single particle with 'tecm' center-of-mass mass into
	'Nop' particles of PDG codes in 'idOut' array. 
	4-momenta of final particles are in 'pf[0]',...,'pf[Nop-1]'.
	
	
Short tutorial:

TODO

Authors:
	
	Radoslaw Kycia

How to cite:
	The program is a simple derivative of [1] and [2]


	[1] R. A. Kycia, J. Turnau, J. J. Chwastowski, R. Staszewski, M. Trzebiński, 'The adaptive Monte Carlo toolbox for phase space integration and generation', Commun. Comput. Phys., 25 (2019), pp. 1547-1563; DOI: 10.4208/cicp.OA-2018-0028; https://arxiv.org/abs/1411.6035

	[2] R. A. Kycia, J. Chwastowski, R. Staszewski, J. Turnau, 'GenEx: A simple generator structure for exclusive processes in high energy collisions', Commun. Comput. Phys., 24 (2018), pp. 860-884, DOI: 10.4208/cicp.OA-2017-0105; https://arxiv.org/abs/1411.6035

License:

	This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>



***********************************************************************/

#include <iostream>
#include <assert.h>
#include <cmath>
#include <complex>
#include <string>
#include <fstream>
#include <vector>
#include <utility>
#include <queue> 

#include <TMath.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TRandom3.h>
#include <TFoam.h>
#include <TFoamIntegrand.h>
#include <TLorentzVector.h>
#include <TDatabasePDG.h>

#include "TGenDecay.h"

using namespace std;

TDatabasePDG * PDGDatabese = TDatabasePDG::Instance();


////////////////////////////////////////////////////////////////////////
//GLOBAL VARIABLES
////////////////////////////////////////////////////////////////////////





////////////////////////////////////////////////////////////////////////
//MAIN
////////////////////////////////////////////////////////////////////////

int main()
{
	
	//MC Generator
	TGenDecay generator;

	
	//center of mass energy
	const double mmu = 0.1057; // muon rest mass GeV
	const double tecm = mmu;   // set muon mass in as total energy in CM frame
	
	//CM 4-momentum - blob that initiates decay
	TLorentzVector pbCM;
	//set blob energy to tecm with zero momentum - rest frame of particle with mass tecm
	pbCM.SetPxPyPzE( 0.0, 0.0, 0.0, tecm );
	

	//outgoing particles
	const int Nop = 2;

	//PDGID outgoing paricles (masses are taken from PDG table)
	//Useful codes: pi0=111, pi+= 211, pi-= -211; full list: https://pdg.lbl.gov/2007/reviews/montecarlorpp.pdf
	int idOut[Nop] = {11,11};  //e, e


	//out particles:
	TLorentzVector pf[Nop];
	
	
	//masses of products - from PDGDatabase
	double mass[Nop];
	
	for(int i=0; i < Nop; i++)
	{
		mass[i] = PDGDatabese->GetParticle( idOut[i] )->Mass();
	}
 
	//set decay configuration
	generator.SetDecay(pbCM, Nop, mass);
	
	
	//phase space element for decay
	double wtdecay = 0.0;
	
	
	//PROLOG USER SECTION (definitions etc.):
	
	//Histograms
		TH1D * h_rapidity = new TH1D("h_rapidity", "h_rapidity; y; events",100,-2.0,2.0);
		//h_rapidity->SetBit(TH1::kCanRebin);
		h_rapidity->Sumw2();
	
		TH1D * h_CM = new TH1D("h_CM", "h_CM; CM[GeV]; events",200,0.0,1.0);
		//h_CM->SetBit(TH1::kCanRebin);
		h_CM->Sumw2();
	
		TH2D * h_etaphi = new TH2D("h_etaphi", "h_etaphi; y; #phi",100,-3.0,3.0,100,-M_PI,M_PI);
		//h_etaphi->SetBit(TH1::kCanRebin);
		h_etaphi->Sumw2();
	

		
		//open file for txt data pf[]
		ofstream myfile;
		myfile.open ("events.txt");
		myfile << "# pf0.px, pf0.py, pf0.pz, pf0.E,...., Weight" << endl;
		
	
	//END OF PROLOG USER SECTION
	
	//integrand
	double sumIntegrand = 0.0;
	double sumIntegrand2 = 0.0;
	
	long NevTot = 10e1; //Total number of events
	
	//GENRATION LOOP
	long   loop;
	for(loop=0; loop<NevTot; loop++)
	{
		
		//Generate event
		
		//make decay
			wtdecay =  generator.Generate();
		
		//get out particles
			for( int i = 0; i < Nop; i++ )
			{
				pf[i] = *(generator.GetDecay( i ));
			}
		
	
		
		//|Matrix Element|^2 - components
					
						 
			 //MATRIX ELEMENT HERE:
			 //PUT HERE SPECIFIC MATRIX ELEMENT. You can use: pf[0]....pf[Nop-1] - 4 momenta of outgoing particles
			 
			double integrand = 1.0;
			
			//update integrad
			sumIntegrand += wtdecay * integrand;
			//update integrand^2
			sumIntegrand2 += pow(wtdecay * integrand, 2);
    
			
		//filling histograms - CHANGE ACCORDINGLY
			h_rapidity->Fill( pf[0].Rapidity(), wtdecay);
			h_rapidity->Fill( pf[1].Rapidity(), wtdecay);
			h_rapidity->Fill( pf[2].Rapidity(), wtdecay);
			
			TLorentzVector pCM;
			
			for( int i = 0; i < Nop; i++)
			{
				pCM += pf[i];
			}
		
			h_CM->Fill( pCM.M(), wtdecay);
			
			for( int i = 0; i < Nop; i++)
			{
				h_etaphi->Fill( pf[i].Rapidity(), pf[i].Phi(), wtdecay);
			}
		
		
		// save data into file
		for( int i = 1; i < Nop; i++ )
		{
				myfile << pf[i].Px() << ", " << pf[i].Py() << ", " << pf[i].Pz() << ", " << pf[i].E() << ", ";
		}
		myfile << wtdecay << endl;
		
		

		
		if( loop % 10000 == 0)
			cout << "loop = " << loop << endl;
		
	}
	
	//END GENRATION LOOP


	cout << "====== Events generated, entering Finalize" << endl;
	
	double integral = sumIntegrand / (double) NevTot;
	double error = ( pow(sumIntegrand,2) - sumIntegrand2/ (double) NevTot) * sqrt(NevTot);

	cout << "Integral = " << integral << " +- " << error << endl;

	double R2 = M_PI * sqrt( ( pow(mmu,2) - pow( mass[0]+mass[1],2) )* (pow(mmu,2)-pow(mass[0]-mass[1],2)) )/ (2.0 * pow(mmu,2) * pow(2.0*M_PI,2));
	
	cout << "theoretical Integral = " << R2 << endl;
	
	cout << "====================================" << endl;
	
	//EPILOG USER SECTION (writing files etc.)
	
	//Saving histograms - root file
		TFile RootFile("histograms.root","RECREATE","histograms");
		RootFile.ls();
		h_rapidity->Write();
		h_CM->Write();
		h_etaphi->Write();
		RootFile.Close();
	
	
	//Saving histograms - pdf file
	
			
	//save histograms
		
		TCanvas* canv1 = new TCanvas("canv","plot");	
		canv1->cd(1);
		h_rapidity->Draw("h");
		canv1->Update();
		canv1->Print( "rapidity.eps", "" );
		delete canv1;
		
		
		TCanvas* canv2 = new TCanvas("canv","plot");	
		canv2->cd(1);
		h_CM->Draw("h");
		canv2->Update();
		canv2->Print( "CM.eps", "" );
		delete canv2;
		
		TCanvas* canv3 = new TCanvas("canv","plot");	
		canv3->cd(1);
		h_etaphi->Draw("h");
		canv3->Update();
		canv3->Print( "etaPhi.eps", "" );
		delete canv3;
	
	
	//save event tree
		string rootFilename("events.root");
		
		TFile file( rootFilename.c_str(), "Update" ); 
		file.SetCompressionLevel(1); 
		file.cd();
		file.Close();
	
	//closing txt file
		myfile.close();
	
	
	//END EPILOG SECTION
	
	
	//CLEANING
		delete h_rapidity;
		delete h_CM;
		delete h_etaphi;
	
	
	cout << "***** End of Spherical Generation Program  *****" << endl;


	return 0;
	
};
