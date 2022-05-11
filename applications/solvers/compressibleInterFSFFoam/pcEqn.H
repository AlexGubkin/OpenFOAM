label   pcRefCell = 0;
scalar pcRefValue = 0.0;

while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix pcEqn
    (
        fvm::laplacian(pc)
        ==
        fvc::div(capillaryForceFluxFiltered)
    );

    pcEqn.setReference(pcRefCell, pcRefValue);
    pcEqn.solve();

    const surfaceVectorField gradPcf(fvc::interpolate(fvc::grad(pc)));

    capillaryForceFluxFiltered = capillaryForceFlux - cutFactorf*((gradPcf - (gradPcf & nHatfv)*nHatfv) & Sf);

//     if (pimple.finalNonOrthogonalIter())
//     {
//     }
}
