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

#include "objectiveChord.H"
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

defineTypeNameAndDebug(objectiveChord, 1);
addToRunTimeSelectionTable
(
    objectiveIncompressible,
    objectiveChord,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

objectiveChord::objectiveChord
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
	int arrSize = sizeof(objectivePatches_)/8.;
 	int i =0; 
	int m ;
	int n ;      
//    		Info<< "MABITE   " << arrSize << endl;
//    		Info<< "MABITE   " << objectivePatches_ << endl;
	double ymaxi[arrSize];   // get size of objective patches 
	double ymini[arrSize];   // get size of objective patches 

    // Read target Chord if present. Else use the current one as a target
     
	for (const label patchi : objectivePatches_)
        {
		const fvPatch& patch = mesh_.boundary()[patchi];
		const fvPatchVectorField& faceCentres = mesh_.C().boundaryField()[patchi]; // get faces of the objective patches
		ymaxi[i] = gMax(faceCentres)[0];
		ymini[i] = gMin(faceCentres)[0];
	i+=1;
        }
double tempMax=0;
double tempMin=0;
 for(n =0;n < arrSize; ++n)
    		{
	if(ymaxi[n]>tempMax)
		tempMax=ymaxi[n];
	if(ymini[n]<tempMin)
		tempMin=ymini[n];
    		}	
       initThick_ = tempMax-tempMin; // add max y coordinate for both patches so we get the maximum Chord
    		Info<< "initThick   " << initThick_ << endl;
    // Allocate boundary field pointers
    bdxdbDirectMultPtr_.reset(createZeroBoundaryPtr<vector>(mesh_));
    bdSdbMultPtr_.reset(createZeroBoundaryPtr<vector>(mesh_));
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


scalar objectiveChord::J()
{
    J_ = Zero;
	int arrSize = sizeof(objectivePatches_)/8.;
 	int i = 0; 
	int m ;
	int n ;      
	double ymaxi[arrSize];   // get size of objective patches 
	double ymini[arrSize];   // get size of objective patches 

    // Read target Chord if present. Else use the current one as a target
     
	for (const label patchi : objectivePatches_)
        {
		const fvPatch& patch = mesh_.boundary()[patchi];
		const fvPatchVectorField& faceCentres = mesh_.C().boundaryField()[patchi]; // get faces of the objective patches
		ymaxi[i] = gMax(faceCentres)[0];
		ymini[i] = gMin(faceCentres)[0];
	i+=1;
        }
double tempMax=0;
double tempMin=0;
 for(n =0;n < arrSize; ++n)
    		{
	if(ymaxi[n]>tempMax)
		tempMax=ymaxi[n];
	if(ymini[n]<tempMin)
		tempMin=ymini[n];
    		}	
       J_ = tempMax-tempMin; // add max y coordinate for both patches so we get the maximum Chord
    J_ -= initThick_;
    J_ /= initThick_;
    		Info<< "J   " << J_  << endl;
    return J_;
}


void objectiveChord::update_dxdbDirectMultiplier()
{
    const scalar oneThird(1.0/3.0);
    for (const label patchi : objectivePatches_)
    {
        const fvPatch& patch = mesh_.boundary()[patchi];
        tmp<vectorField> tnf = patch.nf();
        const vectorField& nf = tnf();
        bdxdbDirectMultPtr_()[patchi] = -oneThird*nf/initThick_;
    }
}


void objectiveChord::update_dSdbMultiplier()
{
    const scalar oneThird(1.0/3.0);
    for (const label patchi : objectivePatches_)
    {
        const fvPatch& patch = mesh_.boundary()[patchi];
        bdSdbMultPtr_()[patchi] = -oneThird*patch.Cf()/initThick_;
    }
}


bool objectiveChord::write(const bool valid) const
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
