label   pcRefCell = 0;
scalar pcRefValue = 0.0;
const scalar cCapillaryForceFluxCoeff(alphaControls.lookupOrDefault<scalar>("cCapillaryForceFluxCoeff", 0.01));

while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix pcEqn
    (
        fvm::laplacian(pc)
        ==
        fvc::div(capillaryForceFlux)
    );

    pcEqn.setReference(pcRefCell, pcRefValue);
    pcEqn.solve();
    
    gradPcf =
        deltaf/(mag(deltaf) + mixture.deltaN())*fvc::interpolate(fvc::grad(pc));

    capillaryForceFlux =
        fvc::interpolate(w2*sigma*kappa)/fvc::interpolate(w2)*deltaf*mesh.magSf();

      capillaryForceFlux -= (gradPcf - (gradPcf & nHatfv)*nHatfv) & mesh.Sf();
      
    if (pimple.finalNonOrthogonalIter())
    {
        const scalar capillaryForceAvg =
            gSum(capillaryForceFlux)
          / gSum(static_cast<surfaceScalarField>(fvc::interpolate(alphaSh)*mesh.magSf()));

        gradSigmaf =
            fvc::interpolate(fvc::grad(sigma));

        tangentGradSigmaf =
            mag(gradAlphaShf)*(gradSigmaf - (gradSigmaf & nHatfv)*nHatfv) & mesh.Sf()/mesh.magSf();

        surfaceTensionForce =
            (
                capillaryForceFlux
              - max
                (
                    min
                    (
                        capillaryForceFlux,
                        cCapillaryForceFluxCoeff*capillaryForceAvg*mesh.magSf()
                    ),
                    -cCapillaryForceFluxCoeff*capillaryForceAvg*mesh.magSf()
                )
            )/mesh.magSf()
          + tangentGradSigmaf;
    }
}
