label pcRefCell = 0;
scalar pcRefValue = 0.0;
const scalar cCapillaryForceFluxCoeff(alphaControls.lookupOrDefault<scalar>("cCapillaryForceFluxCoeff", 0.01));

// pc.storeOldTime();
// pc.oldTime();

phic =
    fvc::interpolate(sigmaPtr->sigma())
   *fvc::interpolate(w2*kappa)
   /fvc::interpolate(w2)*deltaf*mesh.magSf();

while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix pcEqn
    (
        fvm::laplacian(pc)
        ==
        fvc::div(phic)
    );

    pcEqn.setReference(pcRefCell, pcRefValue);
    pcEqn.solve();

    if (pimple.finalNonOrthogonalIter())
    {
        gradPcf =
            deltaf/(mag(deltaf) + mixture.deltaN())*fvc::interpolate(fvc::grad(pc));

        phic -=
            (gradPcf - (gradPcf & nHatfv)*nHatfv) & mesh.Sf();

        const dimensionedScalar capillaryForceAvg
        (
            dimForce/dimVolume,
            gSum(phic)
           /gSum(static_cast<surfaceScalarField>(fvc::interpolate(alphaSh)*mesh.magSf()))
        );

        Info<< "capillaryForceAvg = " << capillaryForceAvg.value() << endl;

        gradSigmaf =
            fvc::interpolate(fvc::grad(sigmaPtr->sigma()));

        tangentGradSigmaf =
            mag(gradAlphaShf)*(gradSigmaf - (gradSigmaf & nHatfv)*nHatfv) & mesh.Sf()/mesh.magSf();

        surfaceTensionForce =
            (
                phic
              - max
                (
                    min
                    (
                        phic,
                        cCapillaryForceFluxCoeff*capillaryForceAvg*mesh.magSf()
                    ),
                  - cCapillaryForceFluxCoeff*capillaryForceAvg*mesh.magSf()
                )
            )/mesh.magSf()
          + tangentGradSigmaf;
    }
}
