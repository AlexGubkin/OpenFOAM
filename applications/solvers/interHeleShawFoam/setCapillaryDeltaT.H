/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

Global
    setDeltaT

Description
    Reset the timestep to maintain a constant maximum courant Number.
    Reduction of time-step is immediate, but increase is damped to avoid
    unstable oscillations.

\*---------------------------------------------------------------------------*/

if (adjustTimeStep && pimple.correct())
{
    const volScalarField nu(mixture.nu());
    const volScalarField sigma(sigmaPtr->sigma());
    const scalarField V(mesh.V());

    scalar VMin = gMin(V);
    scalar VAvg =gSum(V*V)/gSum(V);
    scalar rhoMax = gMax(rho);
    scalar rhoMin = gMin(rho);
    scalar nuMax = gMax(nu);
    scalar nuMin = gMin(nu);
    scalar rhoAvg = 0.5*(rhoMin + rhoMax);
    scalar nuAvg = 0.5*(nuMin + nuMax);
    scalar sigmaMin = gMin(sigma);
    scalar tauMu = nuAvg*rhoAvg*Foam::pow(VAvg, 1.0/3.0)/sigmaMin;
    scalar tauRho = Foam::sqrt(rhoAvg*VAvg/sigmaMin);
    scalar C1 = 0.01;
    scalar C2 = 10.0;
    scalar deltaTFact = 0.9;

    scalar deltaT1 =
        0.5
       *(
            C2*tauMu
          + Foam::sqrt
            (
                C2*tauMu*C2*tauMu
              + 4.0*C1*tauRho*tauRho
            )
        );

    scalar deltaT2 =
        max
        (
            C2*tauMu,
            10.0*C1*tauRho
        );

    Info<<
        "VMin = " << VMin << endl <<
        "VAvg = " << VAvg << endl <<
        "rhoMax = " << rhoMax << endl <<
        "rhoMin = " << rhoMin << endl <<
        "nuMax = " << nuMax << endl <<
        "nuMin = " << nuMin << endl <<
        "rhoAvg = " << rhoAvg << endl <<
        "nuAvg = " << nuAvg << endl <<
        "sigmaMin = " << sigmaMin << endl <<
        "tauMu = " << tauMu << endl <<
        "tauRho = " << tauRho << endl <<
        "deltaT1 = " << deltaT1 << endl <<
        "deltaT2 = " << deltaT2 << endl;
    
    runTime.setDeltaT
    (
        min
        (
            min
            (
                min
                (
                    runTime.deltaTValue(),
                    maxDeltaT
                ),
                deltaTFact*deltaT1
            ),
            deltaTFact*deltaT2
        )
    );

    Info<< "deltaT = " <<  runTime.deltaTValue() << endl;
}

// ************************************************************************* //
