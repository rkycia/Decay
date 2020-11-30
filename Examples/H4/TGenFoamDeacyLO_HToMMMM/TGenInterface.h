/***********************************************************************

 TGenInterface interface
 
 Description: Minimal interface implemented by MC generators for decay
 
 Usage:
  
  The general concept of working with MC generators is as follows:
  1. Prepare 4-momenta of decaying particle P and mass[nt] array of final products.
  2. Initialize generator:  SetDecay(P, nt, mass);
  3. Generate decay: Generate(void);
  4. For each of particles enumerated by 0..nt-1 get their 4-momentum: pfi = GetDecay(i);
  5. Repeat 3 and 4 for another decay.
 
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

#ifndef ROOT_TGenInterface
#define ROOT_TGenInterface

/*!
  @class TGenInterface

  @brief  Minimal interface implemented by MC generators for decay
  
  The general concept of working with MC generators is as follows:
  1. Prepare 4-momenta of decaying particle P and mass[nt] array of final products.
  2. Initialize generator:  SetDecay(P, nt, mass);
  3. Generate decay: Generate(void);
  4. For each of particles enumerated by 0..nt-1 get their 4-momentum: pfi = GetDecay(i);
  5. Repeat 3 and 4 for another decay.
 
*/


#include <TLorentzVector.h>


using namespace std;


class TGenInterface  {
private:  

   Int_t  fNt;  // number of decay particles
   
public:

	///constructor
	TGenInterface(): fNt(0) {};
	
	///destructor
	virtual ~TGenInterface() {};

	/// Sets up configuration of decay
	/// @param P 4-momentum of decaying particle
	/// @param nt number of products of decay
	/// @param mass mass matrix of products of decay mass[nt]
	/// @warning This should be first method to call since it sets up decay configuration.
	virtual Bool_t SetDecay(TLorentzVector &P, Int_t nt, const Double_t *mass)=0;
	
	/// Generate a single decay
	virtual Double_t  Generate(void)=0;
	
	/// @returns 4-vector of n-th product of decay
	///  @param n  number of final particle to get from range 1...nt
	/// @warning You should call Generate() in first place.
	virtual TLorentzVector *GetDecay(Int_t n)=0; 

	/// @returns number of final particles
	virtual Int_t GetNt()const { return fNt;};
   
	
};

#endif

