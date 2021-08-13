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

\*---------------------------------------------------------------------------*/

#include "saturation.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(saturation, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        saturation,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::saturation::writeFileHeader(const label i)
{
//     label oneMinusSumSaturation(0);

    if (Pstream::master())
    {
        writeHeader(file(), "saturation");
        writeCommented(file(), "Time");

        forAll(phaseSet_, phasei)
        {
            const word& phaseName = phaseSet_[phasei];

            if (mesh_.foundObject<volScalarField>(word("alpha." + phaseName)))
                writeTabbed(file(), phaseName);
//             else
//                 oneMinusSumSaturation = phasei;
        }

//         const word& phaseName = phaseSet_[oneMinusSumSaturation];
// 
//         writeTabbed(file(), phaseName);

        file() << endl;
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::saturation::saturation
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    logFiles(obr_, name),
    phaseSet_()
{
    read(dict);
    resetName(typeName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::saturation::~saturation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::saturation::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    if (mesh_.foundObject<IOdictionary>(word("transportProperties")))
    {
        const dictionary&
            caseProperties = mesh_.lookupObject<IOdictionary>("transportProperties");

        caseProperties.lookup("phases") >> phaseSet_;
    }
    else if (mesh_.foundObject<IOdictionary>(word("thermophysicalProperties")))
    {
        const dictionary&
            caseProperties = mesh_.lookupObject<IOdictionary>("thermophysicalProperties");

        caseProperties.lookup("phases") >> phaseSet_;
    }
    else
        FatalErrorInFunction
            << "Can not find transportProperties/thermophysicalProperties file"
            << exit(FatalError);

    return true;
}


bool Foam::functionObjects::saturation::execute()
{
    return true;
}


bool Foam::functionObjects::saturation::write()
{
    logFiles::write();

    if (Pstream::master())
        writeTime(file());

    scalar
        saturationi(0);
//         sumSaturation(0);

    forAll(phaseSet_, phasei)
    {
        const word& phaseName = phaseSet_[phasei];

        if (mesh_.foundObject<volScalarField>(word("alpha." + phaseName)))
        {
            const volScalarField&
                alphai = mesh_.lookupObject<volScalarField>(word("alpha." + phaseName));

            const scalarField&
                V = mesh_.V();

            saturationi =
                    gSum(V * alphai) / gSum(V);

//             sumSaturation += saturationi;

            Info<< "saturation(" << phaseName << ") = "
                <<  saturationi
                << " []" << nl;
        }

        if (Pstream::master())
        {
            file()
                << tab << saturationi;
//                 << endl;
        }
    }

    if (Pstream::master())
    {
        file()
//             << tab << scalar(1.0) - sumSaturation
            << endl;
    }

    return true;
}


// ************************************************************************* //
