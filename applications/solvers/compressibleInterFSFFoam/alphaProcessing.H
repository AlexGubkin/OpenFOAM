const label nAlphaSmoothers(alphaControls.lookupOrDefault<label>("nAlphaSmoothers", 3));
const scalar cSmCoeff(alphaControls.lookupOrDefault<scalar>("cSmCoeff", 0.5));
const scalar cShCoeff(alphaControls.lookupOrDefault<scalar>("cShCoeff", 0.98));

// Smoother loop for SSF
alphaSm = alpha1;

for (int aSmoother=0; aSmoother<nAlphaSmoothers; aSmoother++)
{
    Info<< "Alpha smoother #" << aSmoother + 1 << nl;

//     alphaSm =
//         cSmCoeff*fvc::average(fvc::interpolate(alphaSm))
//       + (1.0 - cSmCoeff)*alphaSm;

    tmp<volScalarField> talphaSm
    (
        cSmCoeff*fvc::average(fvc::interpolate(alphaSm))
      + (1.0 - cSmCoeff)*alphaSm
    );
    alphaSm = talphaSm();
    alphaSm.correctBoundaryConditions();
}

// Sharpener for alpha
alphaC = min(1.0, max(alpha1, 0.0));
alphaC.correctBoundaryConditions();

alphaSh =
    1.0/(1.0 - cShCoeff)
   *(
        min
        (
            max
            (
                alpha1,
                0.5*cShCoeff
            ),
            1.0 - 0.5*cShCoeff
        )
      - 0.5*cShCoeff
    );
alphaSh.correctBoundaryConditions();

w1 = sqrt(alphaC*(1.0 - alphaC));
w2 = sqrt(w1*w1 + SMALL);
w1.correctBoundaryConditions();
w2.correctBoundaryConditions();


// Interpolated face-gradient of alpha
gradAlphaSm = fvc::grad(alphaSm);
gradAlphaSmf = fvc::interpolate(gradAlphaSm);

deltaf = fvc::snGrad(alphaSh);

// Face unit interface normal
// nHatfv = gradAlphaSmf/(mag(gradAlphaSmf) + mixture.deltaN());
nHatv = gradAlphaSm/(mag(gradAlphaSm) + mixture.deltaN());
nHatfv = fvc::interpolate(nHatv);


// Correction for the boundary condition on the unit normal nHat on
// walls to produce the correct contact angle.

// The dynamic contact angle is calculated from the component of the
// velocity on the direction of the interface, parallel to the wall.
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

        fvsPatchVectorField& nHatp = nHatfv.boundaryFieldRef()[patchi];
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
        const scalarField b1(Foam::cos(theta));

        scalarField b2(nHatp.size());
        forAll(b2, facei)
        {
            b2[facei] = Foam::cos(Foam::acos(a12[facei]) - theta[facei]);
        }

        const scalarField det(1.0 - a12*a12);

        scalarField a((b1 - a12*b2)/det);
        scalarField b((b2 - a12*b1)/det);

        nHatp = a*nf + b*nHatp;
        nHatp /= (mag(nHatp) + mixture.deltaN().value());

        a1cap.gradient() = (nf & nHatp)*mag(gradAlphaSmf.boundaryField()[patchi]);
        a1cap.evaluate();
        a2bf[patchi] = 1 - a1cap;
    }
}

// Face unit interface normal flux
nHatf = nHatfv & mesh.Sf();
