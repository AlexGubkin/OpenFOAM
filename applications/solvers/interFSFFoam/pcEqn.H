label pcRefCell = 0;
scalar pcRefValue = 0.0;

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
        phic -= pcEqn.flux();
    }
}
