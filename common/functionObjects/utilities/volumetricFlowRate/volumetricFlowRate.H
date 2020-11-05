/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
     \\/     M anipulation  |
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

Class
    Foam::functionObjects::volumetricFlowRate

Description
    Evaluates and outputs volumetricFlowRate for rhombohedron structure.

    Example of function object specification:
    \verbatim
    volumetricFlowRate1
    {
        type            volumetricFlowRate;
        libs            ("twoPhaseFlowRhombohedronStructureLibutilityFunctionObjects.so");

        phi             \<phi\> / \<rhoPhi\>;
        markerFields    (alpha.\<phase1\> alpha.\<phase2\> alpha.\<phase3\>);
        patches         (\<inflowPatchesSet\>);
    }
    \endverbatim

    Output data is written to the dir postProcessing/volumetricFlowRate1/\<timeDir\>/.

See also
    Foam::functionObject
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::logFiles
    Foam::functionObjects::timeControl

SourceFiles
    volumetricFlowRate.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_volumetricFlowRate_H
#define functionObjects_volumetricFlowRate_H

#include "fvMeshFunctionObject.H"
#include "logFiles.H"
#include "volFields.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
            Class volumetricFlowRate Declaration
\*---------------------------------------------------------------------------*/

class volumetricFlowRate
:
    public fvMeshFunctionObject,
    public logFiles
{
protected:

    // Protected data
    // Read from dictionary

        //- Name of flux field (optional)
        word phiName_;

        //- Marker fields
        wordList phaseSet_;

        //- Patches for flow rate calculation
        labelHashSet patchSet_;


    // Protected Member Functions

        //- File header information
        virtual void writeFileHeader(const label i);


private:

        //- Disallow default bitwise copy construct
        volumetricFlowRate(const volumetricFlowRate&);

        //- Disallow default bitwise assignment
        void operator=(const volumetricFlowRate&);


public:

    //- Runtime type information
    TypeName("volumetricFlowRate");


    // Constructors

        //- Construct from Time and dictionary
        volumetricFlowRate
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~volumetricFlowRate();


    // Member Functions

        //- Read the controls
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Calculate and write the volumetricFlowRate
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
