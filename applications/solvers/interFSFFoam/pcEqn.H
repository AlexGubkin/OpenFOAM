label pcRefCell = 0;
scalar pcRefValue = 0.0;
const scalar cPhicFilterCoeff(alphaControls.lookupOrDefault<scalar>("cPhicFilterCoeff", 0.1));
const scalar cPhicThresholdCoeff(alphaControls.lookupOrDefault<scalar>("cPhicThresholdCoeff", 0.01));

phic.relax();

while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix pcEqn
    (
        fvm::laplacian(pc)
     == fvc::div(phic)
    );

    pcEqn.setReference(pcRefCell, pcRefValue);
    pcEqn.solve();

    if (pimple.finalNonOrthogonalIter())
    {
        gradPcf =
            fvc::interpolate(w1)*fvc::interpolate(fvc::grad(pc));

//         gradPcf =
//             mag(deltaf)/(mag(deltaf) + mixture.deltaN())*fvc::interpolate(fvc::grad(pc));

        phicFilt =
            phicFilt.oldTime()
          + cPhicFilterCoeff*((gradPcf - (gradPcf & nHatfv)*nHatfv) & mesh.Sf());

//         phicFilt =
//             cPhicFilterCoeff*((gradPcf - (gradPcf & nHatfv)*nHatfv) & mesh.Sf());

        phic -= phicFilt;

        const dimensionedScalar fcAvg
        (
            average(mag(phic/mesh.magSf()))
           /(average(pos(mag(phic))) + SMALL)
        );

        Info<< "fcAvg = " << fcAvg.value() << nl;

        phic =
            phic
          - max
            (
                min
                (
                    phic,
                    cPhicThresholdCoeff*fcAvg*mesh.magSf()
                ),
              - cPhicThresholdCoeff*fcAvg*mesh.magSf()
            );
    }
}
