const label nKappaSmoothers(alphaControls.lookupOrDefault<label>("nKappaSmoothers", 2));

// Smoother loop
{
    kappa0 = -fvc::div(nHatf);
    kappa = kappa0;

    for (int kappaSmoother=0; kappaSmoother<nKappaSmoothers; kappaSmoother++)
    {
        Info<< "Kappa smoother #" << kappaSmoother + 1 << nl;

        tmp<volScalarField> tKappaSm
        (
            2.0*w1*kappa0
          + (scalar(1.0) - 2.0*w1)*fvc::surfaceSum(fvc::interpolate(w2*kappa)*mesh.magSf())()
           /fvc::surfaceSum(fvc::interpolate(w2)*mesh.magSf())()
        );
        kappa = tKappaSm();
    }
}
