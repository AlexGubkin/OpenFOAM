const label nAlphaSmoothers(alphaControls.lookupOrDefault<label>("nAlphaSmoothers", 3));
const scalar cSmoothingCoeff(alphaControls.lookupOrDefault<scalar>("cSmoothingCoeff", 0.5));
const scalar cSharpeningCoeff(alphaControls.lookupOrDefault<scalar>("cSharpeningCoeff", 0.5));

// Simple smooth for CSF
// {
//     alphaSm = alpha1;
// 
//     for (int aSmoother=0; aSmoother<nAlphaSmoothers; aSmoother++)
//     {
//         Info<< "Alpha smoother #" << aSmoother + 1 << nl;
// 
//         tmp<volScalarField> talphaSmoothed
//         (
//             fvc::surfaceSum(fvc::interpolate(alphaSm)*magSf)()
//            /fvc::surfaceSum(magSf)
//         );
//         alphaSm = talphaSmoothed();
//     }
// }

// Smoother loop for SSF
{
    alphaSm = alpha1;

    for (int aSmoother=0; aSmoother<nAlphaSmoothers; aSmoother++)
    {
        Info<< "Alpha smoother #" << aSmoother + 1 << nl;

        tmp<volScalarField> talphaSmoothed
        (
            cSmoothingCoeff*fvc::surfaceSum(fvc::interpolate(alphaSm)*magSf)()/fvc::surfaceSum(magSf)
            +(scalar(1) - cSmoothingCoeff)*alphaSm
        );
        alphaSm = talphaSmoothed();
    }
}

// Sharpener for alpha
{
    alphaC = min(scalar(1), max(alpha1, scalar(0)));

    w1 = sqrt(alphaC*(scalar(1) - alphaC));
    w2 = sqrt(w1*w1 + 1e-6);

    alphaSh =
    (
        scalar(1)/(scalar(1) - cSharpeningCoeff)
        *(
           min
           (
                max(alpha1, 0.5*cSharpeningCoeff),
                scalar(1) - 0.5*cSharpeningCoeff
           )
           -0.5*cSharpeningCoeff
        )
    );
}

{
    // Cell gradient of alpha
    const volVectorField gradAlpha(fvc::grad(alphaSm));

    // Interpolated face-gradient of alpha
    surfaceVectorField gradAlphaf(fvc::interpolate(gradAlpha));

    // Face unit interface normal
    surfaceVectorField nHatfv(gradAlphaf/(mag(gradAlphaf) + deltaN_));

    correctContactAngle(nHatfv.boundaryFieldRef(), gradAlphaf.boundaryField());

    // Face unit interface normal flux
    nHatf = nHatfv & Sf;
}
