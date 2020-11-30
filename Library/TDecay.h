/***********************************************************************

TDeacy class

Brief: Original GENBOD algorithm with w few modifications. Adopted from TGenPhaseSpace ROOT package.
 


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

#ifndef ROOT_TDecay
#define ROOT_TDecay

/*!
  @class TDecay

  @brief  Original GENBOD algorithm with w few modifications.  Adopted from TGenPhaseSpace ROOT package.
 
  Tis is a helper class that can be base to construct MC generators that simulate decay.
  
  The volume element of LIPS is :
	 
	\f$(dLIPS = (2\pi)^4 \delta^{(4)}(P-\sum_{i=1}^{n}p_{i}) \prod_{i=1}^{n} \frac{d^{3}p_{i}}{(2\pi)^3 2E_{i}}\f$
  
  See: H.Pilkhun, The Interactions of Hadrons North-Holland 1967, M.D. Schwartz, Quantum Field Theory and the Standard Model, Cambridge UP 2013
  
  The class requires random numbers in Generate(). They are used to make an event.
  The scheme of use:
  1. Prepare 4-momenta of decaying particle P and mass[nt] array of final products.
  2. Initialize generator:  SetDecay(P, nt, mass);
  3. Generate decay: Generate(void);
  4. For each of particles enumerated by 0..nt-1 get their 4-momentum: pfi = GetDecay(i);
  5. Repeat 3 and 4 for another decay.
 
*/


#include "TLorentzVector.h"

#include <iostream>       // std::cin, std::cout
#include <queue>          // std::queue
#include <assert.h>

using namespace std;


class TDecay : public TObject {
private:  
   Int_t        fNt;             // number of decay particles
   Double_t     fMass[18];       // masses of particles
   Double_t     fBeta[3];        // betas of decaying particle
   Double_t     fTeCmTm;         // total energy in the C.M. minus the total mass
   TLorentzVector  fDecPro[18];  //kinematics of the generated particles 

   Double_t PDK(Double_t a, Double_t b, Double_t c);
   
   ///factorial function
   int factorial(int n);  

public:

   /// constructor
   TDecay(): fNt(0), fMass(), fBeta(), fTeCmTm(0.) {};
   
   /// copy constructor
   TDecay(const TDecay &gen);
   
   /// desctructor
   virtual ~TDecay() {};
   
   /// assignment
   TDecay& operator=(const TDecay &gen);

   /// Sets up configuration of decay
   /// @param P 4-momentum of decaying particle
   /// @param nt number of products of decay
   /// @param mass mass matrix of products of decay mass[nt]
   /// @returns kTRUE - the decay is permitted by kinematics; kFALSE - the decay is forbidden by kinematics
   /// @warning This should be first method to call since it sets up decay configuration.
   Bool_t SetDecay(TLorentzVector &P, Int_t nt, const Double_t *mass);
   
   /// Generate a single decay
   /// @param rnd a queue of 3*nt-4 random numbers from external source.
   Double_t        Generate(std::queue<double> &rnd);
   
   /// @returns 4-vector of n-th product of decay
   /// @param n  number of final particle to get from range 1...nt
   /// @warning You should call Generate() in first place.
   TLorentzVector *GetDecay(Int_t n); 

   /// @returns number of final particles
   Int_t    GetNt()      const { return fNt;}

};

#endif

