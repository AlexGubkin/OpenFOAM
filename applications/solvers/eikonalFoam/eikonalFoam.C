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
    potentialFoam

Description
    Potential flow solver which solves for the velocity potential, to
    calculate the flux-field, from which the velocity field is obtained by
    reconstructing the flux.

    This application is particularly useful to generate starting fields for
    Navier-Stokes codes.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
// #include "nonOrthogonalSolutionControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

//     nonOrthogonalSolutionControl eikonal(mesh, "eikonal");

    #include "createControls.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< nl << "Calculating distance field" << endl;

    // Since solver contains no time loop it would never execute
    // function objects so do it ourselves
    runTime.functionObjects().start();

    int iter = 0;
    scalar initialResidual = 0;

    do
    {
        nd = fvc::grad(d);
        nd /= (mag(nd) + small);

        surfaceVectorField nf(fvc::interpolate(nd));
        nf /= (mag(nf) + small);

        surfaceScalarField dPhi("dPhi", nf & mesh.Sf());

        fvScalarMatrix dEqn
        (
            fvm::div(dPhi, d)
          - fvm::Sp(fvc::div(dPhi), d)
          - epsilon*d*fvm::laplacian(d)
        ==
            dimensionedScalar(dimless, 1.0)
        );

        dEqn.relax();
        initialResidual = dEqn.solve().initialResidual();
 
    } while (initialResidual > tolerance && ++iter < maxIter);
}


// ************************************************************************* //
