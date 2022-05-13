label   pcRefCell = 0;
scalar pcRefValue = 0.0;

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
        deltaf*fvc::interpolate(w2*sigma*kappa)/fvc::interpolate(w2)*magSf
      - (gradPcf - (gradPcf & nHatfv)*nHatfv) & Sf

    if (pimple.finalNonOrthogonalIter())
    {
        gradSigmaf =
            fvc::interpolate(fvc::grad(sigma));

        tangentGradSigmaf =
            mag(gradAlphaShf)*(gradSigmaf - (gradSigmaf & nHatfv)*nHatfv) & Sf/magSf;

        surfaceTensionForce = capillaryForceFlux/magSf + tangentGradSigmaf
    }
}
