const label nKappaSmoothers(alphaControls.lookupOrDefault<label>("nKappaSmoothers", 2));

// Smoother loop
kappa0 = -fvc::div(nHatf);
kappa = kappa0;

for (int kappaSmoother=0; kappaSmoother<nKappaSmoothers; kappaSmoother++)
{
    Info<< "Kappa smoother #" << kappaSmoother + 1 << nl;

//     kappa =
//         2.0*w1*kappa0
//       + (1.0 - 2.0*w1)*fvc::average(fvc::interpolate(w2*kappa))
//        /fvc::average(fvc::interpolate(w2));

    tmp<volScalarField> tKappa
    (
        2.0*w1*kappa0
      + (1.0 - 2.0*w1)*fvc::average(fvc::interpolate(w2*kappa))
       /fvc::average(fvc::interpolate(w2))
    );
    kappa = tKappa();
    kappa.correctBoundaryConditions();
}
