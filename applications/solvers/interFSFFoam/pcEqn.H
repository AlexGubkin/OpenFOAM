fcf.storePrevIter();

const scalar cPhicFilterCoeff(alphaControls.lookupOrDefault<scalar>("cPhicFilterCoeff", 0.1));
const scalar cPhicThresholdCoeff(alphaControls.lookupOrDefault<scalar>("cPhicThresholdCoeff", 0.01));

label pcRefCell = 0;
scalar pcRefValue = 0.0;

// gradPcf =
//     fvc::interpolate(w1)*fvc::interpolate(fvc::grad(pc));

gradPcf =
    fvc::interpolate(fvc::grad(pc));

fcf =
    fvc::interpolate(sigmaPtr->sigma())
   *fvc::interpolate(w2*kappa)
   /fvc::interpolate(w2)*deltaf;

fc = fvc::reconstruct(fcf*mesh.magSf());

fcfFilt =
    mag(deltaf)/(mag(deltaf) + mixture.deltaN())
   *(
        fcf.oldTime()
      + cPhicFilterCoeff*((gradPcf - (gradPcf & nHatfv)*nHatfv) & mesh.Sf()/mesh.magSf())
    );

fcFilt = fvc::reconstruct(fcfFilt*mesh.magSf());

fcf -= fcfFilt;

fcf.relax();

while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix pcEqn
    (
        fvm::laplacian(pc)
     == fvc::div(fcf*mesh.magSf())
    );

    pcEqn.setReference(pcRefCell, pcRefValue);
    pcEqn.solve();

    if (pimple.finalNonOrthogonalIter())
    {
        const dimensionedScalar fcAvg
        (
            average(mag(fcf))
           /(average(pos(mag(fcf))) + SMALL)
        );

        Info<< "fcAvg = " << fcAvg.value() << nl;

//         phic = fcf*mesh.magSf() - pcEqn.flux();
        phic = (fcf - fvc::snGrad(pc))*mesh.magSf();
        phic -=
            max
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
