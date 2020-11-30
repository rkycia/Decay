/***********************************************************************

 TGenFoamDecay class
 
 
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

#ifndef ROOT_TGenFoamDecay
#define ROOT_TGenFoamDecay

/*!
  @class TGenFoamDecay

  @brief  Provides adaptive MC simulator of decay and integrator over LIPS.
 
	This class is multipurpose MC generator for decay which also integrate functions over phase space of decay. It uses FOAM adaptive integrator.
	
	For default integrated function is 1. 
	
	The best way to use this class is to make your stub class that inherits after TGenFoamDecay and redefine Integrand() method filling it with a integrand that will be integrated over LIPS (Lorent Invariant Phase Space). The volume element of LIPS is :
	 
	\f$dLIPS = (2\pi)^4 \delta^{(4)}(P-\sum_{i=1}^{n}p_{i}) \prod_{i=1}^{n} \frac{d^{3}p_{i}}{(2\pi)^3 2E_{i}}\f$
	
	See: H.Pilkhun, The Interactions of Hadrons North-Holland 1967, M.D. Schwartz, Quantum Field Theory and the Standard Model, Cambridge UP 2013
	
	Generator implements interface TGenInterface and work with it is along this interface.
	
	The general way of working with the generator is as follows:
	1. Prepare 4-momenta of decaying particle P and mass[nt] array of final products.
	2. Initialize generator:  SetDecay(P, nt, mass);
	3. Generate decay: Generate(void);
	4. For each of particles enumerated by 0..nt-1 get their 4-momentum: pfi = GetDecay(i);
	5. Repeat 3 and 4 for another decay.
	
	If you prepared the integrand to be integrated over LIPS then you can print its value by calling Finalize() method or get the numerical value by GetIntegMC().
	
	@warning If the integrand has 'small' support in phase space and you get some nonsense result, then you probably want to increase nSampl and nCells. It will increase probability that adaptive MC will spot support of the function.
	
	Citation info and detailed description:
	
	TODO
	

*/

#include <TGenInterface.h>

#include <TLorentzVector.h>
#include <TFoam.h>
#include <TFoamIntegrand.h>
#include <TRandom3.h>

#include <stdlib.h>
#include <iostream>       // std::cin, std::cout
#include <queue>          // std::queue
#include <assert.h>

#include <TDecay.h>

using namespace std;


class TGenFoamDecay : public TFoamIntegrand, public TGenInterface  {
private:  
   Int_t        fNt;             // number of decay particles
   Double_t     fMass[18];       // masses of particles
   Double_t     fBeta[3];        // betas of decaying particle
   Double_t     fTeCmTm;         // total energy in the C.M. minus the total mass
   TLorentzVector  fDecPro[18];  //kinematics of the generated particles
   Int_t       seed; 			 //seed for pseudorandom generator

   TDecay _decay;				 //decay engine
   TFoam*  _foam;				 //adaptive integrator
   TRandom3 _pseRan;				 // pseudorandom number generator 
   
	//FOAM parameters
	Int_t  nCells;   // Number of Cells
	Int_t  nSampl;   // Number of MC events per cell in build-up
	Int_t  nBin;   // Number of bins in build-up
	Int_t  OptRej;   // Wted events for OptRej=0; wt=1 for OptRej=1 (default)
	Int_t  OptDrive;   // (D=2) Option, type of Drive =0,1,2 for TrueVol,Sigma,WtMax
	Int_t  EvPerBin;   // Maximum events (equiv.) per bin in buid-up
	Int_t  Chat;   // Chat level



	/// @returns weight of the process
	/// @param nDim  dimension of integration
	/// @param Xarg vector of probablilistic variables from [0;1] of dim nDim
	/// @warning it is required by Foam integrator
	Double_t Density(int nDim, Double_t *Xarg);
	


public:

   /// constructor
   TGenFoamDecay(): fNt(0), fMass(), fBeta(), fTeCmTm(0.), seed(4357), nCells(1000), nSampl(1000), nBin(8), OptRej(1), OptDrive(2), EvPerBin(25), Chat(1) {_foam = new TFoam("FoamX");};
   /// copy constructor
   TGenFoamDecay(const TGenFoamDecay &gen);
   /// desctructor
   virtual ~TGenFoamDecay() { delete _foam;}
   /// assignemt
   TGenFoamDecay& operator=(const TGenFoamDecay &gen);

   /// Sets up configuration of decay
   /// @param P 4-momentum of decaying particle (Momentum, Energy units are Gev/C, GeV)
   /// @param nt number of products of decay
   /// @param mass mass matrix of products of decay mass[nt]
   /// @returns kTRUE - the decay is permitted by kinematics; kFALSE - the decay is forbidden by kinematics
   /// @warning This should be first method to call since it sets up decay configuration.	
   /// @warning The method also initialize FOAM.
   Bool_t SetDecay(TLorentzVector &P, Int_t nt, const Double_t *mass);
   
   /// Generate a single decay
   Double_t Generate(void);
   
   /// Collect 4-vector of products of decay
   ///  @param n  number of final particle to get from range 1...nt
   /// @warning You should call Generate() in first place.
   TLorentzVector *GetDecay(Int_t n); 

   /// @returns 4-vector of n-th product of decay
   ///  @param n  number of final particle to get from range 1...nt
   Int_t    GetNt()      const { return fNt;};
   	
	/// @returns the function under integral over LIPS (Lorentz Invariant Phase Space)
	/// @param fNt  number of outgoing particles
	/// @param pf  array of TLorentzVectors of outgoing particles
	/// @warning It is set to 1.0. You should to redefine it (in your derived class, after inheritance) when you want to use full adaptation features.
	virtual Double_t Integrand( int fNt, TLorentzVector* pf );
	
	///sets seed for pseudorandom number generator
	virtual void setSeed( Int_t seed ) {seed = seed; _pseRan.SetSeed(seed);};
	
	///finalize Foam printing out MC integral and statistics
	virtual void Finalize( void );
	
	///@returns integral +- error of MC integration of function in Integrand() over LIPS
	virtual void GetIntegMC(Double_t & inetgral, Double_t & error);
	
	///Sets FOAM number of Cells
	/// @warning It should be done before Generate()
	void SetnCells(Int_t nc) {nCells = nc;};
	
	///Sets FOAM number of MC events per cell in build-up
	/// @warning It should be done before Generate()
	void SetnSampl(Int_t ns) {nSampl = ns;};
	
	///Sets FOAM number of bins in build-up
	/// @warning It should be done before Generate()
	void SetnBin(Int_t nb) {nBin = nb;};
	
	///Sets FOAM Weigh events for OptRej=0; wt=1 for OptRej=1 (default)
	/// @warning It determines if the events will be weighted of not (weight=1).
	/// @warning It should be done before Generate()
	void SetOptRej(Int_t OptR) {OptRej = OptR;};
	
	///Sets FOAM (D=2) option, type of Drive =0,1,2 for TrueVol,Sigma,WtMax
	/// @warning It should be done before Generate()
	void SetOptDrive(Int_t OptD) {OptDrive = OptD;};
	
	///Sets FOAM maximum events (equiv.) per bin in buid-up
	/// @warning It should be done before Generate()
	void SetEvPerBin(Int_t Ev) {EvPerBin = Ev;};
	
	///Sets FOAM chat level
	/// @warning It should be done before Generate()
	void SetChat(Int_t Ch) {Chat = Ch;};
	
};

#endif

