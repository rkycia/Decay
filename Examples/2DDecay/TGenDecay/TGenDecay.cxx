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

#include "TGenDecay.h"


//const Int_t kMAXP = 18;


//__________________________________________________________________________________________________
TGenDecay::TGenDecay(const TGenDecay &gen) : TObject(gen)
{
   //copy constructor
   fNt      = gen.fNt;
   fTeCmTm  = gen.fTeCmTm;
   fBeta[0] = gen.fBeta[0];
   fBeta[1] = gen.fBeta[1];
   fBeta[2] = gen.fBeta[2];
   _decay   = gen._decay;
   _pseRan  = gen._pseRan;
   for (Int_t i=0;i<fNt;i++) 
   {
      fMass[i]   = gen.fMass[i];
      fDecPro[i] = gen.fDecPro[i];
   }
}

//__________________________________________________________________________________________________
TGenDecay& TGenDecay::operator=(const TGenDecay &gen)
{
   // Assignment operator
   TObject::operator=(gen);
   fNt      = gen.fNt;
   fTeCmTm  = gen.fTeCmTm;
   fBeta[0] = gen.fBeta[0];
   fBeta[1] = gen.fBeta[1];
   fBeta[2] = gen.fBeta[2];
   _decay   = gen._decay;
   _pseRan  = gen._pseRan;
   for (Int_t i=0;i<fNt;i++) 
   {
      fMass[i]   = gen.fMass[i];
      fDecPro[i] = gen.fDecPro[i];
   }
   return *this;
}
   
//__________________________________________________________________________________________________
Double_t TGenDecay::Generate(void)
{

	//clear queue - in case there was a previous run
	while (!rndQueue.empty())
	{
		rndQueue.pop();
	}

	//number of degrees of freedom in the decay
	int nDim = 3*fNt-4;
	
	//put rnd numbers into queue
	for( int i = 0; i < nDim; i++)
	{
		rndQueue.push( _pseRan.Rndm() );
	}  

	return _decay.Generate( rndQueue );

}

//__________________________________________________________________________________
TLorentzVector *TGenDecay::GetDecay(Int_t n) 
{ 
   //return Lorentz vector corresponding to decay n
   if (n>fNt) return 0;
   
   return _decay.GetDecay(n);
}

//_____________________________________________________________________________________
Bool_t TGenDecay::SetDecay(TLorentzVector &P, Int_t nt, const Double_t *mass) 
{

   kMAXP = nt;

   Int_t n;
   fNt = nt;
   if (fNt<2 || fNt>18) return kFALSE;  // no more then 18 particle

   //
   fTeCmTm = P.Mag();           // total energy in C.M. minus the sum of the masses
   for (n=0;n<fNt;n++) {
      fMass[n]  = mass[n];
      fTeCmTm  -= mass[n];
   }

   if (fTeCmTm<=0) return kFALSE;    // not enough energy for this decay
   
   _pseRan.SetSeed(seed);			 //set seed

   _decay.SetDecay(P, fNt, fMass);  //set decay to TDecay

   return kTRUE; 
}
