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

#include "volumetricFlowRate.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(volumetricFlowRate, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        volumetricFlowRate,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::volumetricFlowRate::writeFileHeader(const label i)
{
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

//     label oneMinusSumAlphaPhase(0);

    if (Pstream::master())
    {
        writeHeader(file(), "volumetricFlowRate");
        writeCommented(file(), "");

        forAllConstIter(labelHashSet, patchSet_, iter)
        {
            label patchi = iter.key();

            writeTabbed(file(), pbm[patchi].name());

            for(label phasei(0); phasei < phaseSet_.size() - 1; phasei++)
//             forAll(phaseSet_, phasei)
                writeTabbed(file(), "");
        }

        file() << endl;

        writeCommented(file(), "Time");

        forAllConstIter(labelHashSet, patchSet_, iter)
        {
            forAll(phaseSet_, phasei)
            {
                const word& phaseName = phaseSet_[phasei];

                if (mesh_.foundObject<volScalarField>(word("alpha." + phaseName)))
                    writeTabbed(file(), phaseName);
//                 else
//                     oneMinusSumAlphaPhase = phasei;
            }
        }

//         const word& phaseName = phaseSet_[oneMinusSumAlphaPhase];
// 
//         writeTabbed(file(), phaseName);

        file() << endl;
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::volumetricFlowRate::volumetricFlowRate
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    logFiles(obr_, name),
    phiName_(),
    phaseSet_(),
    patchSet_()
{
    read(dict);
    resetName(typeName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::volumetricFlowRate::~volumetricFlowRate()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::volumetricFlowRate::read(const dictionary& dict)
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

    phiName_ = dict.lookupOrDefault<word>("phi", "phi");

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    patchSet_ =
        mesh_.boundaryMesh().patchSet
        (
            wordReList(dict.lookupOrDefault("patches", wordReList()))
        );

    Info<< type() << " " << name() << ":" << nl;

    if (patchSet_.empty())
    {
        if (Pstream::master())
        {
            FatalErrorInFunction
                << "Unable to find patches"
                << exit(FatalError);
        }
    }

    Info<< nl;

    Info<< "    patches =";

    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        label patchi = iter.key();

        Info<< " " << pbm[patchi].name();
    }

    Info<< nl << nl;

    return true;
}


bool Foam::functionObjects::volumetricFlowRate::execute()
{
    return true;
}


bool Foam::functionObjects::volumetricFlowRate::write()
{
    logFiles::write();

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    const surfaceScalarField&
        phi = mesh_.lookupObject<surfaceScalarField>(phiName_);

    if (Pstream::master())
        writeTime(file());

    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        label patchi = iter.key();

        Info<< pbm[patchi].name() << ":" << nl;

        scalar
            alphaiVolumetricFlowRate(0);
//             sumVolumetricFlowRate(0);

        forAll(phaseSet_, phasei)
        {
            const word& phaseName = phaseSet_[phasei];

            if (mesh_.foundObject<volScalarField>(word("alpha." + phaseName)))
            {
                const volScalarField&
                    alphai = mesh_.lookupObject<volScalarField>(word("alpha." + phaseName));

                if (phi.dimensions() == dimMass/dimTime)
                {
                    const volScalarField&
                        rho = mesh_.lookupObject<volScalarField>("rho");

    //                 label
    //                     patchID = mesh_.boundaryMesh().findPatchID(pbm[patchi].name());

                    alphaiVolumetricFlowRate =
                        gSum(
                            phi.boundaryField()[patchi]
                          * alphai.boundaryField()[patchi]
                          / rho.boundaryField()[patchi]
                        );

//                     sumVolumetricFlowRate += alphaiVolumetricFlowRate;

                    Info<< "    volumetricFlowRate(" << phaseName << ") = "
                        << alphaiVolumetricFlowRate << " [m^3 / s]" << nl;
                }
                else if (phi.dimensions() == dimVolume/dimTime)
                {
                    alphaiVolumetricFlowRate =
                        gSum(
                            phi.boundaryField()[patchi]
                          * alphai.boundaryField()[patchi]
                        );

//                     sumVolumetricFlowRate += alphaiVolumetricFlowRate;

                    Info<< "    volumetricFlowRate(" << phaseName << ") = "
                        << alphaiVolumetricFlowRate << " [m^3 / s]" << nl;
                }
                else
                {
                    FatalErrorInFunction
                        << "Incompatible dimensions for phi: " << phi.dimensions() << nl
                        << "Dimensions should be " << dimMass/dimTime << " or "
                        << dimVolume/dimTime << exit(FatalError);
                }

                if (Pstream::master())
                {
                    file()
                        << tab << alphaiVolumetricFlowRate;
//                         << endl;
                }
            }
        }
    }

    if (Pstream::master())
    {
        file()
//             << tab << scalar(1.0) - sumMassFlowRate
            << endl;
    }

    return true;
}


// ************************************************************************* //
