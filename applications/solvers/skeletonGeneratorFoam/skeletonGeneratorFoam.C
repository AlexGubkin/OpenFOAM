/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

Application
    laplacianFoam

Description
    Solves a simple Laplace equation, e.g. for thermal diffusion in a solid.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvModels.H"
#include "fvConstraints.H"
// #include "simpleControl.H"
#include "pimpleControl.H"
#include "fvCellSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"

    #include "createTime.H"
    #include "createMesh.H"
    #include "createDyMControls.H"
//     #include "createControl.H"

//     simpleControl simple(mesh);

    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating S distribution\n" << endl;

    fvCellSet set(transportProperties, mesh);

    while (pimple.run(runTime))
    {
        const scalarField& V = mesh.V();

        skeletonMarker = (S < threshold) ? 1 : 0;

        const scalar m = gSum(V*skeletonMarker)/gSum(V);

        Info<< "m = " << m << nl
            << "|dm| = " << (Foam::mag(m0-m)).value() << nl
            << nl << endl;

        // Do any mesh changes
        mesh.update();

        runTime++;

        Info<< "Time = " << runTime.userTimeName() << nl << endl;

        fvModels.correct();

        while (pimple.correctNonOrthogonal())
        {
            fvScalarMatrix SEqn
            (
                fvm::ddt(S) - fvm::laplacian(DS, S)
             ==
                fvModels.source(S)
            );

            fvConstraints.constrain(SEqn);
            SEqn.solve();
            fvConstraints.constrain(S);
        }

        #include "write.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
