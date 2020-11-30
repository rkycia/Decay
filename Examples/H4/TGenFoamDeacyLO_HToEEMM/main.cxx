/***********************************************************************
Toy MC spherical generator.

The code realizes decay H-> 2e 2mu of

S. Boselli, C. M. Carloni Calame, G. Montagna, O. Nicrosini, and F. Piccinini, Higgs boson
decay into four leptons at NLOPS electroweak accuracy, JHEP 06 (2015) 023, arXiv:1503.07394
[hep-ph]

A. Bredenstein, A. Denner, S. Dittmaier, and M. Weber, Precise predictions for the Higgs-boson
decay H → WW/ZZ → 4 leptons, Phys. Rev. D 74 (2006) 013004, arXiv:hep-ph/0604011.

Description:
  
	Spherical decay of a single particle with 'tecm' center-of-mass mass into
	'Nop' particles of PDG codes in 'idOut' array. 
	4-momenta of final particles are in 'pf[0]',...,'pf[Nop-1]'.
	
	
Short tutorial:

	1. Write your integrand expression in Generator::Integrand() method. Example is provided.
	2. Set tecm - mass in CM of decaying particle
	3. Set Nop - number of outhoing particles.
	4. Set IdOut - PDG identificators of outgoing particles, see: https://pdg.lbl.gov/2007/reviews/montecarlorpp.pdf
	5. Done. Type in in console: make run


Warning:
 	
	It makes spherical decay and is unsiutable for strong non-spherical cuts, e.g., 
	on p_{t} of some particles.

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
#include <TLorentzVector.h>
#include <TDatabasePDG.h>

#include "TGenFoamDecay.h"

using namespace std;

TDatabasePDG * PDGDatabese = TDatabasePDG::Instance();


////////////////////////////////////////////////////////////////////////
//GLOBAL VARIABLES
////////////////////////////////////////////////////////////////////////

class Generator: public TGenFoamDecay 
{
Double_t Integrand( int fNt, TLorentzVector * pf );	
	
};


Double_t Generator::Integrand( int fNt, TLorentzVector * pf )
{

	//Here put your expression to integrate over LIPS
	// example provided below
	
	//constants
	const double pi = M_PI;

    const double am_e = 0.000511;
    const double am_mu = 0.105658;

	const double am_H = 125.1;
    const double am_W = 80.379;
    const double am_Z = 91.1876;
    const double Gam_Z = 2.4952;

    const double GF = 1.166379*1.e-5;   //in GeV^{-2}
    const double alpha_G = sqrt(2.0)*GF*pow(am_W,2)/pi*(1.0-pow(am_W,2)/pow(am_Z,2));
    const double gemG2 = 4.0*pi*alpha_G;
    const double gemG = sqrt(gemG2);

	//coupling constants

    const double  aQl = -1.0;
    const double aI3W = -0.5;

    const double cW = am_W/am_Z;
    const double sW = sqrt(1.0-pow(cW,2));

    const double ag6 = -sW/cW * aQl;
    const double ag7 = ag6 + aI3W/(sW*cW);

    const double g_HZZ = am_Z/cW/sW;
      
	//invariant masses and scalar products
	double am12 = (pf[0] + pf[1]).M();
	double am34 = (pf[2] + pf[3]).M();
	
	double prod_p1p2 = pf[0]*pf[1];
	double prod_p1p3 = pf[0]*pf[2];
	double prod_p1p4 = pf[0]*pf[3];
	double prod_p2p3 = pf[1]*pf[2];
	double prod_p2p4 = pf[1]*pf[3];
	double prod_p3p4 = pf[2]*pf[3];
	
	//matrix element squared
	//H -> ZZ -> e+e-mu+mu- (LO)
	//one diagram
    
	double prop2_12 = 1.0/(pow((pow(am12,2)-pow(am_Z,2)),2)+pow((am_Z*Gam_Z),2));
	double prop2_34 = 1.0/(pow((pow(am34,2)-pow(am_Z,2)),2)+pow((am_Z*Gam_Z),2));
    
	double amp2 = 4.0*pow(ag6,2)*pow(ag7,2)*pow(am_e,2)*pow(am_mu,2) \
                  + ag6*ag7*(pow(ag6,2)+pow(ag7,2))*pow(am_mu,2)*prod_p1p2 \
                  + ag6*ag7*(pow(ag6,2)+pow(ag7,2))*pow(am_e,2)*prod_p3p4 \
                  + 2.0*pow(ag6,2)*pow(ag7,2)*prod_p1p4*prod_p2p3 \
                  + (pow(ag6,4)+pow(ag7,4))*prod_p1p3*prod_p2p4;
	
	double amat2 = pow((pow(gemG,3)*g_HZZ),2) * 16.0*amp2 * prop2_12 * prop2_34;
	
	double anorm = 1.0 / (2.0*am_H);
	
	amat2 = amat2*anorm;
	
	return amat2;
};	

////////////////////////////////////////////////////////////////////////
//MAIN
////////////////////////////////////////////////////////////////////////

int main()
{
	
	//MC Generator
	Generator generator;

	
	//center of mass energy
	const double am_H = 125.1;
	const double tecm = am_H; // mass in CM of decaying particle (in GeV)
	
	//CM 4-momentum - blob that initiates decay
	TLorentzVector pbCM;
	//set blob energy to tecm with zero momentum - rest frame of particle with mass tecm
	pbCM.SetPxPyPzE( 0.0, 0.0, 0.0, tecm );
	
	//outgoing particles
	const int Nop = 4;

	//PDGID outgoing paricles (masses are taken from PDG table)
	//Useful codes: pi0=111, pi+= 211, pi-= -211; full list: https://pdg.lbl.gov/2007/reviews/montecarlorpp.pdf
	int idOut[Nop] = {11,-11, 13,-13};  //e-, e+, mu-, mu+

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
		TH1D * h_Rapidity = new TH1D("h_Rapidity", "h_Rapidity; Rapidity; events",100,-2.0,2.0);
		//h_Rapidity->SetCanExtend(TH1::kXaxis); 
		h_Rapidity->Sumw2();
	
		TH1D * h_CM = new TH1D("h_CM", "h_CM; CM [GeV]; events",100,0.0,200.0);
		//h_CM->SetCanExtend(TH1::kXaxis); 
		h_CM->Sumw2();
	
		TH2D * h_RapidityPhi = new TH2D("h_RapidityPhi", "h_RapidityPhi; Rapidity; #phi",100,-3.0,3.0,100,-M_PI,M_PI);
		//h_RapidityPhi->SetCanExtend(TH2::kXaxis); 
		h_RapidityPhi->Sumw2();
        
		TH1D * h_M12 = new TH1D("h_M12", "h_M12; M_{12} [GeV]; events",100,0.0,100.0);
		//h_M12->SetCanExtend(TH1::kXaxis); 
		h_M12->Sumw2();
	
		TH1D * h_M34 = new TH1D("h_M34", "h_M34; M_{34} [GeV]; events",100,0.0,100.0);
		//h_M34->SetCanExtend(TH1::kXaxis); 
		h_M34->Sumw2();
        
        TH1D * h_M13 = new TH1D("h_M13", "h_M13; M_{13} [GeV]; events",100,0.0,100.0);
		//h_M13->SetCanExtend(TH1::kXaxis); 
		h_M13->Sumw2();
	
		TH1D * h_M24 = new TH1D("h_M24", "h_M24; M_{24} [GeV]; events",100,0.0,100.0);
		//h_M24->SetCanExtend(TH1::kXaxis); 
		h_M24->Sumw2();
        
		TH2D * h_M12M34 = new TH2D("h_M12M34", "h_M12M34; M_{12} [GeV]; M_{34} [GeV]",100,0.0,100,100,0.0,100.0);
        //h_M12M34->SetCanExtend(TH2::kXaxis); 
		h_M12M34->Sumw2();
		
		TH2D * h_M13M24 = new TH2D("h_M13M24", "h_M13M24; M_{13} [GeV]; M_{24} [GeV]",100,0.0,100,100,0.0,100.0);
        //h_M13M24->SetCanExtend(TH2::kXaxis); 
		h_M13M24->Sumw2();
	
		
		//open file for txt data pf[]
		ofstream myfile;
		myfile.open ("events.txt");
		myfile << "# pf0.px, pf0.py, pf0.pz, pf0.E, ..., Weight" << endl;
		
	
	//END OF PROLOG USER SECTION
	
	
	long NevTot = 1e6; // Total number of events to generate
	
	//GENRATION LOOP
	long   loop;
	for(loop=0; loop<NevTot; loop++)
	{
		
		//Generate event
		
		//make decay
			wtdecay =  generator.Generate();
			//cout << "wt " << wtdecay << endl;
		
		//get out particles
			for( int i = 0; i < Nop; i++ )
			{
				pf[i] = *(generator.GetDecay( i ));
			}
    
			
		//filling histograms - CHANGE ACCORDINGLY
			h_Rapidity->Fill( pf[0].Rapidity(), wtdecay);
			h_Rapidity->Fill( pf[1].Rapidity(), wtdecay);
			h_Rapidity->Fill( pf[2].Rapidity(), wtdecay);
			
			TLorentzVector pCM;
			
			for( int i = 0; i < Nop; i++)
			{
				pCM += pf[i];
			}
		
			h_CM->Fill( pCM.M(), wtdecay);
			
			for( int i = 0; i < Nop; i++)
			{
				h_RapidityPhi->Fill( pf[i].Rapidity(), pf[i].Phi(), wtdecay);
			}
			
            Double_t M12 = (pf[0] + pf[1]).M();
            Double_t M34 = (pf[2] + pf[3]).M();
            Double_t M13 = (pf[0] + pf[2]).M();
            Double_t M24 = (pf[1] + pf[3]).M();

			h_M12->Fill( M12, wtdecay);
			h_M34->Fill( M34, wtdecay);
			h_M13->Fill( M13, wtdecay);
			h_M24->Fill( M24, wtdecay);
			
            h_M12M34->Fill( M12, M34, wtdecay);
            h_M13M24->Fill( M13, M24, wtdecay);
            
		// save data into file
		for( int i = 0; i < Nop; i++ )
		{
				myfile << pf[i].Px() << ", " << pf[i].Py() << ", " << pf[i].Pz() << ", " << pf[i].E() << ", ";
		}
		myfile << wtdecay << endl;
        // myfile << wtdecay << setprecision(8) << endl;
		
		

		
		if( loop % 10000 == 0)
			cout << "loop = " << loop << endl;
	}
	
	//END GENRATION LOOP
	cout << "====== Events generated, entering Finalize" << endl;
	
	generator.Finalize();
	


	
	cout << "====================================" << endl;
	
	//EPILOG USER SECTION (writing files etc.)
	
	//Saving histograms - root file
		TFile RootFile("histograms.root","RECREATE","histograms");
		RootFile.ls();
		h_Rapidity->Write();
		h_CM->Write();
		h_RapidityPhi->Write();
		h_M12->Write();
		h_M34->Write();
		h_M13->Write();
		h_M24->Write();
		h_M12M34->Write();
		h_M13M24->Write();
		RootFile.Close();
	
	//Saving histograms - eps/pdf file
	
		TCanvas* canv1 = new TCanvas("canv","plot",500);	
		canv1->cd(1);
		h_Rapidity->Draw("h");
		canv1->Update();
		canv1->Print( "Rapidity.eps", "" );
		delete canv1;
		
		TCanvas* canv2 = new TCanvas("canv","plot",500);	
		canv2->cd(1);
		h_CM->Draw("h");
		canv2->Update();
		canv2->Print( "CM.eps", "" );
		delete canv2;
		
		TCanvas* canv3 = new TCanvas("canv","plot",500);	
		canv3->cd(1);
		h_RapidityPhi->Draw("Colz");
		canv3->Update();
		canv3->Print( "RapidityPhi.eps", "" );
		delete canv1;
	
        TCanvas* canv11 = new TCanvas("canv","plot",500);	
		canv11->cd(1);
		h_M13->Draw("h");
		canv11->Update();
		canv11->Print( "M13.eps", "" );
		delete canv2;
        
        TCanvas* canv12 = new TCanvas("canv","plot",500);	
		canv12->cd(1);
		h_M24->Draw("h");
		canv12->Update();
		canv12->Print( "M24.eps", "" );
		delete canv2;
        
        TCanvas* canv13 = new TCanvas("canv","plot",500);	
		canv13->cd(1);
		h_M12->Draw("h");
		canv13->Update();
		canv13->Print( "M12.eps", "" );
		delete canv2;
        
        TCanvas* canv14 = new TCanvas("canv","plot",500);	
		canv14->cd(1);
		h_M34->Draw("h");
		canv14->Update();
		canv14->Print( "M34.eps", "" );
		delete canv2;
        
		TCanvas* canv15 = new TCanvas("canv","plot",500,500);	
		canv15->cd(1);
		h_M12M34->Draw("Colz");
		canv15->Update();
		canv15->Print( "M12M34.eps", "" );
        
		TCanvas* canv16 = new TCanvas("canv","plot",500,500);	
		canv16->cd(1);
		h_M13M24->Draw("Colz");
		canv16->Update();
		canv16->Print( "M13M24.eps", "" );
	
	//closing txt file
	myfile.close();
	
	
	
	//END EPILOG SECTION
	
	
	//CLEANING
		delete h_Rapidity;
		delete h_CM;
		delete h_RapidityPhi;
		delete h_M12;
		delete h_M34;
		delete h_M13;
		delete h_M24;
		delete h_M12M34;
		delete h_M13M24;
	
	
	cout << "***** End of Spherical Generation Program  *****" << endl;


	return 0;
	
};
