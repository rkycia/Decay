/***********************************************************************

 TGenDecay class


Brief: Adapter class for TDecay

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

#ifndef ROOT_TGenDecay
#define ROOT_TGenDecay

/*!
  @class TGenDecay

  @brief  Adapter class for TDecay. 
  
	The class can generate kinematics and the volume element of this event in LIPS (Lorentz Invariant Phase Space).
  
	The volume element of LIPS is :
	 
	\f$dLIPS = (2\pi)^4 \delta^{(4)}(P-\sum_{i=1}^{n}p_{i}) \prod_{i=1}^{n} \frac{d^{3}p_{i}}{(2\pi)^3 2E_{i}}\f$
	
	See: H.Pilkhun, The Interactions of Hadrons North-Holland 1967, M.D. Schwartz, Quantum Field Theory and the Standard Model, Cambridge UP 2013
	
	Generator implements interface TGenInterface and work with it is along this interface.
	
	The general way of working with the generator is as follows:
	1. Prepare 4-momenta of decaying particle P and mass[nt] array of final products.
	2. Initialize generator:  SetDecay(P, nt, mass);
	3. Generate decay: Generate(void);
	4. For each of particles enumerated by 0..nt-1 get their 4-momentum: pfi = GetDecay(i);
	5. Repeat 3 and 4 for another decay.

*/

#include <TGenInterface.h>
#include <TDecay.h>

#include <TLorentzVector.h>
#include <TRandom3.h>

#include <iostream>       // std::cin, std::cout
#include <queue>          // std::queue
#include <assert.h>


using namespace std;


class TGenDecay : public TObject, public TGenInterface {
private:  
   Int_t        fNt;             // number of decay particles
   Double_t     fMass[18];       // masses of particles
   Double_t     fBeta[3];        // betas of decaying particle
   Double_t     fTeCmTm;         // total energy in the C.M. minus the total mass
   TLorentzVector  fDecPro[18];  //kinematics of the generated particles 
   Int_t         seed;			 //seed for pseudorandom number generator

   TDecay _decay;				 //decay engine
   TRandom3 _pseRan;			 //pseudorandom numbers
   
    queue<double> rndQueue;		 //queue for random numbers

public:

   /// constructor
   TGenDecay(): fNt(0), fMass(), fBeta(), fTeCmTm(0.), seed(4357) {}
   /// copy constructor
   TGenDecay(const TGenDecay &gen);
   /// desctuctor
   virtual ~TGenDecay() {}
   /// assignment
   TGenDecay& operator=(const TGenDecay &gen);

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
   Int_t    GetNt()      const { return fNt;}
   
   ///sets seed for pseudorandom number generator
   void setSeed( Int_t seed ) {seed = seed; _pseRan.SetSeed(seed);};

};

#endif

