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

    while (pimple.run(runTime))
    {
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

        const scalarField& V = mesh.V();

        skeletonMarker = (S > threshold) ? 1 : 0;

        porosity =
                1.0 - gSum(V*skeletonMarker)/gSum(V);

        Info<< "porosity = " << porosity << nl
            << "|dm| = " << (Foam::mag(m-porosity)).value() << nl
            << nl << endl;

        #include "write.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;



//     Info<< "\nCalculating S distribution\n" << endl;
//
//     while
//     (
//         runTime.run()
// //      || ((Foam::mag(m - porosity)).value() > 1e-3)
//     )
//     {
//         runTime++;
//
//         Info<< "Time = " << runTime.timeName() << nl << endl;
//
//         // Do any mesh changes
//         mesh.update();
//
//         while (pimple.correctNonOrthogonal())
//         {
//             fvScalarMatrix SEqn
//             (
//                 fvm::ddt(S) - fvm::laplacian(DS, S)
//              ==
//                 fvModels.source(S)
//             );
//
//             fvModels.constrain(SEqn);
//             SEqn.solve();
//             fvModels.constrain(S);
//         }
//
//         #include "write.H"
//
//         const scalarField& V = mesh.V();
//
//         skeletonMarker = (S > threshold) ? 1 : 0;
//
//         porosity =
//                 1.0 - gSum(V*skeletonMarker)/gSum(V);
//
//         Info<< "porosity = " << porosity << nl
//             << "|dm| = " << (Foam::mag(m-porosity)).value() << nl
//             << nl << endl;
//
//         Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
//             << "  ClockTime = " << runTime.elapsedClockTime() << " s"
//             << nl << endl;
//     }
//
//     Info<< "End\n" << endl;
//
//     return 0;
}


// ************************************************************************* //
