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
            +(scalar(1.0) - cSmoothingCoeff)*alphaSm
        );
        alphaSm = talphaSmoothed();
    }
}

// Sharpener for alpha
{
    alphaC = min(scalar(1.0), max(alpha1, scalar(0)));

    w1 = sqrt(alphaC*(scalar(1.0) - alphaC));
    w2 = sqrt(w1*w1 + SMALL);

    alphaSh =
    (
        scalar(1.0)/(scalar(1.0) - cSharpeningCoeff)
        *(
           min
           (
                max(alpha1, 0.5*cSharpeningCoeff),
                scalar(1.0) - 0.5*cSharpeningCoeff
           )
           -0.5*cSharpeningCoeff
        )
    );
}

{
    // Cell gradient of alpha
    gradAlphaSm = fvc::grad(alphaSm);

    // Interpolated face-gradient of alpha
    gradAlphaSmf = fvc::interpolate(gradAlphaSm);

    // Face unit interface normal
    nHatfv = gradAlphaSmf/(mag(gradAlphaSmf) + mixture.deltaN());

    // Face unit interface normal flux
    nHatf = nHatfv & Sf;
}

// Correction for the boundary condition on the unit normal nHat on
// walls to produce the correct contact angle.

// The dynamic contact angle is calculated from the component of the
// velocity on the direction of the interface, parallel to the wall.
{
    volScalarField::Boundary& a1bf = alpha1.boundaryFieldRef();
    volScalarField::Boundary& a2bf = alpha2.boundaryFieldRef();

    const fvBoundaryMesh& boundary = mesh.boundary();

    forAll(boundary, patchi)
    {
        if (isA<alphaContactAngleFvPatchScalarField>(a1bf[patchi]))
        {
            alphaContactAngleFvPatchScalarField& a1cap =
                refCast<alphaContactAngleFvPatchScalarField>
                (
                    a1bf[patchi]
                );

            fvsPatchVectorField& nHatp = nHatfv[patchi];
            const scalarField theta
            (
                degToRad(a1cap.theta(U.boundaryField()[patchi], nHatp))
            );

            const vectorField nf
            (
                boundary[patchi].nf()
            );

            // Reset nHatp to correspond to the contact angle

            const scalarField a12(nHatp & nf);
            const scalarField b1(cos(theta));

            scalarField b2(nHatp.size());
            forAll(b2, facei)
            {
                b2[facei] = cos(acos(a12[facei]) - theta[facei]);
            }

            const scalarField det(1.0 - a12*a12);

            scalarField a((b1 - a12*b2)/det);
            scalarField b((b2 - a12*b1)/det);

            nHatp = a*nf + b*nHatp;
            nHatp /= (mag(nHatp) + deltaN_.value());

            a1cap.gradient() = (nf & nHatp)*mag(gradAlphaSmf[patchi]);
            a1cap.evaluate();
            a2bf[patchi] = 1 - a1cap;
        }
    }
}
