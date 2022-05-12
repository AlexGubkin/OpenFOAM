const label nKSmoothers(alphaControls.lookupOrDefault<label>("nKSmoothers", 2));

// Smoother loop
{
    K0 = -fvc::div(nHatf_);
    K = K0;

    for (int KSmoother=0; KSmoother<nKSmoothers; KSmoother++)
    {
        Info<< "K smoother #" << KSmoother + 1 << nl;

        tmp<volScalarField> tKSmoothed
        (
            2*w1*K0
            +(scalar(1) - 2*w1)*fvc::surfaceSum(fvc::interpolate(w2*K)*magSf)()/fvc::surfaceSum(fvc::interpolate(w2)*magSf)()
        );
        K = tKSmoothed();
    }
}
