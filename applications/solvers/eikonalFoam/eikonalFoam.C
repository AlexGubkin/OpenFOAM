/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    eikonalFoam

Description
    Eikonal equation solver calculates distance field d.

    Convergence is accelerated by first generating an approximate solution
    using one of the simpler methods, e.g. Poisson.

    References:
    \verbatim
        P.G. Tucker, C.L. Rumsey, R.E. Bartels, R.T. Biedron,
        "Transport equation based wall distance computations aimed at flows
        with time-dependent geometry",
        NASA/TM-2003-212680, December 2003.
    \endverbatim
    
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "wallPolyPatch.H"
#include "simpleControl.H"
#include "nonOrthogonalSolutionControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);
    nonOrthogonalSolutionControl eikonal(mesh, "eikonal");

    #include "createControls.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< nl << "Calculating distance field" << endl;

    // Since solver contains no time loop it would never execute
    // function objects so do it ourselves
    runTime.functionObjects().start();

    label dPsiRefCell = 0;
    scalar dPsiRefValue = 0.0;

    // Non-orthogonal velocity potential corrector loop
    while (simple.correctNonOrthogonal())
    {
        fvScalarMatrix dPsiEqn
        (
            fvm::laplacian(dPsi)
         ==
            dimensionedScalar(dimless, -1.0)
        );

//         dPsiEqn.relax();
        dPsiEqn.setReference(dPsiRefCell, dPsiRefValue);
        dPsiEqn.solve();

        if (simple.finalNonOrthogonalIter())
        {
            graddPsi = fvc::grad(dPsi);
//             volScalarField magGraddPsi(mag(graddPsi));

            d = sqrt(magSqr(graddPsi) + 2*dPsi) - mag(graddPsi);
        }
    }

    // Write dPsi and graddPsi
    dPsi.write();
    graddPsi.write();

    // Non-orthogonal eikonal corrector loop
    while (eikonal.correctNonOrthogonal())
    {
        gradd = fvc::grad(d);

        surfaceVectorField graddf(fvc::interpolate(gradd));

        surfaceScalarField dPhi("dPhi", graddf & mesh.Sf());

        fvScalarMatrix dEqn
        (
            fvm::div(dPhi, d)
          - fvm::Sp(fvc::div(dPhi), d)
          - epsilon*d*fvm::laplacian(d)
         ==
            dimensionedScalar(dimless, 1.0)
        );

        dEqn.relax();
        dEqn.solve();
    }

    // Write d and gradd
    d.write();
    gradd.write();

    // Non-orthogonal eikonal corrector loop
    while (eikonal.correctNonOrthogonal())
    {
        d2 = d2*pos(mag(gradd) - dimensionedScalar(dimless, 0.8));

        gradd2 = fvc::grad(d2);

        surfaceVectorField gradd2f(fvc::interpolate(gradd2));

        surfaceScalarField d2Phi("d2Phi", gradd2f & mesh.Sf());

        fvScalarMatrix d2Eqn
        (
            fvm::div(d2Phi, d2)
          - fvm::Sp(fvc::div(d2Phi), d2)
          - epsilon*d2*fvm::laplacian(d2)
         ==
            dimensionedScalar(dimless, 1.0)
        );

        d2Eqn.relax();
        d2Eqn.solve();
    }

//     d2 = d2*pos(mag(gradd) - dimensionedScalar(dimless, 0.8));
//     gradd2 = fvc::grad(d2);

    // Write d and gradd
    d2.write();
    gradd2.write();

//     int iter = 0;
//     scalar initialResidual = 0;
//     do
//     {
//         nd = fvc::grad(d);
//         nd /= (mag(nd) + small);
// 
//         surfaceVectorField nf(fvc::interpolate(nd));
//         nf /= (mag(nf) + small);
// 
//         surfaceScalarField dPhi("dPhi", nf & mesh.Sf());
// 
//         fvScalarMatrix dEqn
//         (
//             fvm::div(dPhi, d)
//           - fvm::Sp(fvc::div(dPhi), d)
//           - epsilon*d*fvm::laplacian(d)
//         ==
//             dimensionedScalar(dimless, 1.0)
//         );
// 
//         dEqn.relax();
//         initialResidual = dEqn.solve().initialResidual();
//  
//     } while (initialResidual > tolerance && ++iter < maxIter);
}


// ************************************************************************* //
