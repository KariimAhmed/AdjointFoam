/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019 PCOpt/NTUA
    Copyright (C) 2013-2019 FOSS GP
    Copyright (C) 2019-2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "objectiveThickness.H"
#include "createZeroField.H"
#include "IOmanip.H"
#include "addToRunTimeSelectionTable.H"
#include <iostream>
#include <algorithm>
#include <cmath>  
using namespace std;
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace objectives
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(objectiveThickness, 1);
addToRunTimeSelectionTable
(
    objectiveIncompressible,
    objectiveThickness,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

objectiveThickness::objectiveThickness
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& adjointSolverName,
    const word& primalSolverName
)
:
    objectiveIncompressible(mesh, dict, adjointSolverName, primalSolverName),
    initThick_(Zero),
    objectivePatches_
    (
        mesh_.boundaryMesh().patchSet
        (
            dict.get<wordRes>("patches")
        ).sortedToc()
    )
{
	int arrSize = sizeof(objectivePatches_);
 	int i ; 
	int m ;
	int n ;      
	double ymaxi[arrSize];   // get size of objective patches 

    // Read target thickness if present. Else use the current one as a target
     
	for (const label patchi : objectivePatches_)
        {

		const fvPatch& patch = mesh_.boundary()[patchi];
		const fvPatchVectorField& faceCentres = mesh_.C().boundaryField()[patchi]; // get faces of the objective patches
		int arrSize2 = sizeof(faceCentres[patchi]); 
		double y[arrSize2];
		
		for(m = 0;m < arrSize2; ++m) // scan all faces for each patch to find max y coordinate
		{
		const vector& c = faceCentres[m];// get coordinate for face centre
		double x = c[0];
		double y = abs(c[1]); // abs because on pressure side y is negative
		
		}
		int yarr=sizeof(y);
		 //double arr=y;		
		for(i = 1;i < yarr; ++i)
    		{
       
       		if(y[0] < y[i])
           	y[0] = y[i];
    		}
		  ymaxi[patchi]= y[0]; // max y for each patch so we have max y for suction patch and max y for pressure patch

        }
    
 for(n =0;n < arrSize; ++n)
    		{
       initThick_ += ymaxi[n]; // add max y coordinate for both patches so we get the maximum thickness
    		}	


    // Allocate boundary field pointers
    bdxdbDirectMultPtr_.reset(createZeroBoundaryPtr<vector>(mesh_));
    bdSdbMultPtr_.reset(createZeroBoundaryPtr<vector>(mesh_));
}
	int i =0; 
	int m =0;
	int n =0; 
        double x=0;
        double y=0;
	double c=0;
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


scalar objectiveThickness::J()
{
    J_ = Zero;
	int arrSize = sizeof(objectivePatches_);
 	int i ; 
	int m ;
	int n ;      
	double ymaxi[arrSize];   // get size of objective patches 


    for (const label patchi : objectivePatches_)
    {
        const fvPatch& patch = mesh_.boundary()[patchi];
	const fvPatchVectorField& faceCentres = mesh_.C().boundaryField()[patchi];
		int arrSize2 = sizeof(faceCentres[patchi]); 
		double y[arrSize2];
		
		for(m = 0;m < arrSize2; ++m) // scan all faces for each patch to find max y coordinate
		{
		const vector& c = faceCentres[m];// get coordinate for face centre
		double x = c[0];
		double y = abs(c[1]);
		
		}
		int yarr=sizeof(y);
		 //double arr=y;		
		for(i = 1;i < yarr; ++i)
    		{
       
       		if(y[0] < y[i])
           	y[0] = y[i];
    		}
		  ymaxi[patchi]= y[0]; // max y for each patch so we have max y for suction patch and max y for pressure patch		
    }
for(n =0;n < arrSize; ++n)
    		{
       J_ += ymaxi[n]; // add max y coordinate for both patches so we get the maximum thickness
    		}

    J_ -= initThick_;
    J_ /= initThick_;
    return J_;
}




//- dxdbDirectMultiplier is Term multiplying delta(x)/delta b at the boundary
         //- for objectives that directly depend on x
         //- Needed in both FI and SI computations
// but no need for it here as we deal with x and y only 

//void objectiveThickness::update_dxdbDirectMultiplier()
//{
   // const scalar oneThird(1.0/3.0);
   // for (const label patchi : objectivePatches_)
    //{
    //    const fvPatch& patch = mesh_.boundary()[patchi];
    //    tmp<vectorField> tnf = patch.nf();
    //    const vectorField& nf = tnf();
    //    bdxdbDirectMultPtr_()[patchi] = -oneThird*nf/initThick_; 
    //}
//}




bool objectiveThickness::write(const bool valid) const
{
    if (Pstream::master())
    {
        // File is opened only upon invocation of the write function
        // in order to avoid various instantiations of the same objective
        // opening the same file
        unsigned int width = IOstream::defaultPrecision() + 6;
        if (objFunctionFilePtr_.empty())
        {
            setObjectiveFilePtr();
            objFunctionFilePtr_() << setw(4)     << "#"               << " ";
            objFunctionFilePtr_() << setw(width) << "(Y - YInit)/YInit" << " ";
            objFunctionFilePtr_() << setw(width) << "YInit" << endl;
        }

        objFunctionFilePtr_() << setw(4)     << mesh_.time().value() << " ";
        objFunctionFilePtr_() << setw(width) << J_ << " ";
        objFunctionFilePtr_() << setw(width) << initThick_ << endl;
    }

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace objectives
} // End namespace Foam

// ************************************************************************* //
