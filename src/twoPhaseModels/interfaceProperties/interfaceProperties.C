/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "interfaceProperties.H"
#include "alphaContactAngleFvPatchScalarField.H"
#include "unitConversion.H"
#include "surfaceInterpolate.H"
// #include "fvCFD.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "fvcSurfaceIntegrate.H"
#include "fvcReconstruct.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Correction for the boundary condition on the unit normal nHat on
// walls to produce the correct contact angle.

// The dynamic contact angle is calculated from the component of the
// velocity on the direction of the interface, parallel to the wall.

void Foam::interfaceProperties::correctContactAngle
(
    surfaceVectorField::Boundary& nHatb,
    const surfaceVectorField::Boundary& gradAlphaf
)
{
    const fvMesh& mesh = alpha1_.mesh();
    volScalarField::Boundary& a1bf = alpha1_.boundaryFieldRef();
    volScalarField::Boundary& a2bf = alpha2_.boundaryFieldRef();

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

            fvsPatchVectorField& nHatp = nHatb[patchi];
            const scalarField theta
            (
                degToRad(a1cap.theta(U_.boundaryField()[patchi], nHatp))
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

            a1cap.gradient() = (nf & nHatp)*mag(gradAlphaf[patchi]);
            a1cap.evaluate();
            a2bf[patchi] = 1 - a1cap;
        }
    }
}

void Foam::interfaceProperties::alphaSmoother()
{
    const fvMesh& mesh = alpha1_.mesh();
    const surfaceScalarField& magSf = mesh.magSf();

    const dictionary& alphaControls = mesh.solverDict(alpha1_.name());
    const label nAlphaSmoothers(alphaControls.lookupOrDefault<label>("nAlphaSmoothers", 2));
    const scalar cSmoothingCoeff(alphaControls.lookupOrDefault<label>("cSmoothingCoeff", 0.5));

    // Simple smooth for CSF
//     alphaSmoothed_ = alpha1_;
// 
//     for (int aSmoother=0; aSmoother<nAlphaSmoothers; aSmoother++)
//     {
//         Info<< "Alpha smoother #" << aSmoother + 1 << nl;
// 
//         tmp<volScalarField> talphaSmoothed
//         (
//             fvc::surfaceSum(fvc::interpolate(alphaSmoothed_)*magSf)()
//            /fvc::surfaceSum(magSf)
//         );
//         alphaSmoothed_ = talphaSmoothed();
//     }

    // Smoother loop for SSF
    alphaSmoothed_ = alpha1_;

    for (int aSmoother=0; aSmoother<3; aSmoother++)
    {
        Info<< "Alpha smoother #" << aSmoother + 1 << nl;

        tmp<volScalarField> talphaSmoothed
        (
            cSmoothingCoeff*fvc::surfaceSum(fvc::interpolate(alphaSmoothed_)*magSf)()/fvc::surfaceSum(magSf)
            +(1.0 - cSmoothingCoeff)*alphaSmoothed_
        );
        alphaSmoothed_ = talphaSmoothed();
    }
}

void Foam::interfaceProperties::alphaSharpener()
{
    const fvMesh& mesh = alpha1_.mesh();
    const surfaceScalarField& magSf = mesh.magSf();

    const dictionary& alphaControls = mesh.solverDict(alpha1_.name());
    const scalar cSharpeningCoeff(alphaControls.lookupOrDefault<label>("cSharpeningCoeff", 0.5));

    alphaCut_ = min(scalar(1), max(alpha1_, scalar(0)));

    w1_ = sqrt(alphaCut_*(1.0 - alphaCut_));
    w2_ = sqrt(w1_*w1_ + 1e-6);

    alphaSharpened_ =
    (
        1.0/(1.0 - cSharpeningCoeff)
        *(
           min
           (
                max(alpha1_, 0.5*cSharpeningCoeff),
                1.0 - 0.5*cSharpeningCoeff
           )
           -0.5*cSharpeningCoeff
        )
    );
}

void Foam::interfaceProperties::calculateK()
{
    const fvMesh& mesh = alpha1_.mesh();
    const surfaceVectorField& Sf = mesh.Sf();
    const surfaceScalarField& magSf = mesh.magSf();

    const dictionary& alphaControls = mesh.solverDict(alpha1_.name());
    const label nKSmoothers(alphaControls.lookupOrDefault<label>("nKSmoothers", 2));

    alphaSmoother();

    alphaSharpener();

    // Cell gradient of alpha
    const volVectorField gradAlpha(fvc::grad(alphaSmoothed_, "nHat"));

    // Interpolated face-gradient of alpha
    surfaceVectorField gradAlphaf(fvc::interpolate(gradAlpha));

    // gradAlphaf -=
    //    (mesh.Sf()/mesh.magSf())
    //   *(fvc::snGrad(alpha1_) - (mesh.Sf() & gradAlphaf)/mesh.magSf());

    // Face unit interface normal
    surfaceVectorField nHatfv(gradAlphaf/(mag(gradAlphaf) + deltaN_));
    // surfaceVectorField nHatfv
    // (
    //     (gradAlphaf + deltaN_*vector(0, 0, 1)
    //    *sign(gradAlphaf.component(vector::Z)))/(mag(gradAlphaf) + deltaN_)
    // );
    correctContactAngle(nHatfv.boundaryFieldRef(), gradAlphaf.boundaryField());

    // Face unit interface normal flux
    nHatf_ = nHatfv & Sf;

    // Cell gradient of smoothed alpha
//     const volVectorField gradAlphaSmoothed(fvc::grad(alphaSmoothed_, "nHatSmoothed"));
// 
//     // Interpolated face-gradient of smoothed alpha
//     surfaceVectorField gradAlphaSmoothedf(fvc::interpolate(gradAlphaSmoothed));
// 
//     // Face unit interface normal calculated by smoothed alpha
//     surfaceVectorField nHatSmoothedfv(gradAlphaSmoothedf/(mag(gradAlphaSmoothedf) + deltaN_));
//     correctContactAngle(nHatSmoothedfv.boundaryFieldRef(), gradAlphaSmoothedf.boundaryField());

//     nHatSmoothedf_ = nHatSmoothedfv & Sf;

    // Simple expression for curvature
//     K_ = -fvc::div(nHatf_);
    // surface curvature for CSF model
//     K_ = -fvc::div(nHatSmoothedf_);
    // surface curvature for SSF model
//     const volScalarField w2Smoothed
//     (
//         fvc::surfaceSum(fvc::interpolate(w2_)*magSf)()/fvc::surfaceSum(magSf)
//     );

    // Smoother loop
//     K0_ = -fvc::div(nHatSmoothedf_);
    K0_ = -fvc::div(nHatf_);
    K_ = K0_;

    for (int KSmoother=0; KSmoother<nKSmoothers; KSmoother++)
    {
        Info<< "K smoother #" << KSmoother + 1 << nl;

        tmp<volScalarField> tKSmoothed
        (
            2.0*w1_*K0_
            +(1.0 - 2.0*w1_)*fvc::surfaceSum(fvc::interpolate(w2_*K_)*magSf)()/fvc::surfaceSum(fvc::interpolate(w2_)*magSf)()
        );
        K_ = tKSmoothed();
    }

//     Info<< "K smoother #" << 3 << nl;
// 
//     tmp<volScalarField> tKSmoothed
//     (
//         fvc::surfaceSum(fvc::interpolate(w2_*K_)*magSf)()
//        /fvc::surfaceSum(fvc::interpolate(w2_)*magSf)()
//     );
//     K_ = tKSmoothed();





//     const volScalarField KSmoothed1
//     (
//         2.0*w1_*K
//        +(1.0 - 2.0*w1_)*fvc::surfaceSum(fvc::interpolate(w2_*K)*magSf)()/fvc::surfaceSum(magSf)
//        /w2Smoothed
//     );
// 
//     const volScalarField KSmoothed2
//     (
//         2.0*w1_*K
//        +(1.0 - 2.0*w1_)*fvc::surfaceSum(fvc::interpolate(w2_*KSmoothed1*magSf)()/fvc::surfaceSum(magSf)
//        /w2Smoothed
//     );

//     K_ = fvc::surfaceSum(fvc::interpolate(w2_*KSmoothed2)*magSf)()/fvc::surfaceSum(fvc::interpolate(w2_)*magSf)();

    // Complex expression for curvature.
    // Correction is formally zero but numerically non-zero.
    /*
    volVectorField nHat(gradAlpha/(mag(gradAlpha) + deltaN_));
    forAll(nHat.boundaryField(), patchi)
    {
        nHat.boundaryField()[patchi] = nHatfv.boundaryField()[patchi];
    }

    K_ = -fvc::div(nHatf_) + (nHat & fvc::grad(nHatfv) & nHat);
    */
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceProperties::interfaceProperties
(
    volScalarField& alpha1,
    volScalarField& alpha2,
    const volVectorField& U,
    const IOdictionary& dict
)
:
    transportPropertiesDict_(dict),

    sigmaPtr_(surfaceTensionModel::New(dict, alpha1.mesh())),

    deltaN_
    (
        "deltaN",
        1e-8/pow(average(alpha1.mesh().V()), 1.0/3.0)
    ),

    alpha1_(alpha1),
    alpha2_(alpha2),
    U_(U),

    nHatf_
    (
        IOobject
        (
            "nHatf",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        dimensionedScalar(dimArea, 0)
    ),

    K0_
    (
        IOobject
        (
            "interfaceProperties:K0",
            alpha1_.time().timeName(),
            alpha1_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha1_.mesh(),
        dimensionedScalar(dimless/dimLength, 0)
    ),

    K_
    (
        IOobject
        (
            "interfaceProperties:K",
            alpha1_.time().timeName(),
            alpha1_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha1_.mesh(),
        dimensionedScalar(dimless/dimLength, 0)
    ),

    alphaSmoothed_
    (
        IOobject
        (
            "alphaSmoothed",
            alpha1_.time().timeName(),
            alpha1_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha1_.mesh(),
        dimensionedScalar(dimless, 0)
    ),

    alphaCut_
    (
        IOobject
        (
            "alphaCut",
            alpha1_.time().timeName(),
            alpha1_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha1_.mesh(),
        dimensionedScalar(dimless, 0)
    ),

    alphaSharpened_
    (
        IOobject
        (
            "alphaSharpened",
            alpha1_.time().timeName(),
            alpha1_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha1_.mesh(),
        dimensionedScalar(dimless, 0)
    ),

    w1_
    (
        IOobject
        (
            "w1",
            alpha1_.time().timeName(),
            alpha1_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha1_.mesh(),
        dimensionedScalar(dimless, 0)
    ),

    w2_
    (
        IOobject
        (
            "w2",
            alpha1_.time().timeName(),
            alpha1_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha1_.mesh(),
        dimensionedScalar(dimless, 0)
    ),

    p_c_
    (
        IOobject
        (
            "p_c",
            alpha1_.time().timeName(),
            alpha1_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha1_.mesh(),
        dimensionedScalar(dimPressure, 0),
        zeroGradientFvPatchScalarField::typeName
    )
{
    calculateK();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::interfaceProperties::sigmaK() const
{
    return sigmaPtr_->sigma()*K_;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::interfaceProperties::surfaceTensionForce() const
{
    /*Original source code*/
//     return fvc::interpolate(sigmaK())*fvc::snGrad(alpha1_);

    /*Filtred surface force (FSF) approach*/
    const fvMesh& mesh = alpha1_.mesh();
    const surfaceVectorField& Sf = mesh.Sf();
    const surfaceScalarField& magSf = mesh.magSf();

    // Cell gradient of alpha
    const volVectorField gradAlpha(fvc::grad(alpha1_, "nHat"));

    // Interpolated face-gradient of alpha
    const surfaceVectorField gradAlphaf(fvc::interpolate(gradAlpha));

    // Face unit interface normal
    const surfaceVectorField nHatfv(gradAlphaf/(mag(gradAlphaf) + deltaN_));

    // Cell gradient of smoothed alpha
//     const volVectorField gradAlphaSmoothed(fvc::grad(alphaSmoothed_, "nHatSmoothed"));
    const volVectorField gradAlphaSharpened(fvc::grad(alphaSharpened_));

    // Interpolated face-gradient of smoothed alpha
//     const surfaceVectorField gradAlphaSmoothedf(fvc::interpolate(gradAlphaSmoothed));

    // Face unit interface normal calculated by smoothed alpha
//     const surfaceVectorField nHatSmoothedfv(gradAlphaSmoothedf/(mag(gradAlphaSmoothedf) + deltaN_));

    // Sigma
    const volScalarField sigma(sigmaPtr_->sigma());

    // Cell gradient of sigma
    const volVectorField gradSigma(fvc::grad(sigma));

    // Interpolated face-gradient of sigma
    const surfaceVectorField gradSigmaf(fvc::interpolate(gradSigma));

    // Interpolated tangent face-gradient of sigma
    const surfaceScalarField tangentGradSigmaf((gradSigmaf - (gradSigmaf & nHatSmoothedfv)*nHatSmoothedfv) & Sf/magSf);

    label p_cRefCell = 0;
    scalar p_cRefValue = 0.0;

    while (pimple.correctNonOrthogonal())
    {
        const surfaceVectorField f_n
        (
            fvc::interpolate(sigma)*fvc::interpolate(w2_*K_)/fvc::interpolate(w2_)*fvc::snGrad(alphaSharpened_)
           -fvc::snGrad(p_c_) - (fvc::grad(p_c_) & nHatSmoothedfv)*nHatSmoothedfv
        );

        fvScalarMatrix p_cEqn
        (
            fvm::laplacian(p_c_)
         ==
            fvc::div(f_n)
        );

        p_cEqn.setReference(p_cRefCell, p_cRefValue);
        p_cEqn.solve();
    }

    return fvc::interpolate(sigmaK())*fvc::snGrad(alphaSharpened_) + tangentGradSigmaf*fvc::interpolate(mag(gradAlphaSharpened));
//     return fvc::interpolate(sigmaK())*fvc::snGrad(alpha1_) + tangentGradSigmaf*mag(gradAlphaf);
}


Foam::tmp<Foam::volScalarField>
Foam::interfaceProperties::nearInterface() const
{
    return pos0(alpha1_ - 0.01)*pos0(0.99 - alpha1_);
}


void Foam::interfaceProperties::correct()
{
    calculateK();
}


bool Foam::interfaceProperties::read()
{
    sigmaPtr_->readDict(transportPropertiesDict_);

    return true;
}


// ************************************************************************* //

    /*Sharp surface tension force (SSF) approach*/

//     const fvMesh& mesh = alpha1_.mesh();
//     const scalarField& V = mesh.V();
//     const vectorField& CC = mesh.C();
// 
//     const scalar epsilon(0.99);
//     const scalar alpha10(0.5);
// 
//     volScalarField delta
//     (
//         IOobject
//         (
//             "delta",
//             alpha1_.time().timeName(),
//             mesh
//         ),
//         mesh,
//         dimensionedScalar(dimless, 0)
//     );
// 
//     forAll(CC, CVCi)
//     {
//         delta[CVCi] =
//             (Foam::mag(alpha1_[CVCi] - alpha10) <= 0.5*epsilon)
//             ? 1.0/epsilon*(1.0 + Foam::cos(Foam::constant::mathematical::twoPi*(alpha1_[CVCi] - alpha10)/epsilon))
//             : 0;
//     }
// 
//     volScalarField V
//     (
//         IOobject
//         (
//             "V",
//             alpha1_.time().timeName(),
//             mesh
//         ),
//         mesh,
//         dimensionedScalar(dimVolume, SMALL)
//     );
// 
//     V.ref() = mesh.V();
// 
//     // Cell gradient of sigma
// //     const volVectorField gradSigma(fvc::grad(sigmaPtr_->sigma(), "sigma"));
//     const volVectorField gradSigmaByV(fvc::grad(sigmaPtr_->sigma())/V);
// 
//     // Interpolated face-gradient of sigma
//     surfaceVectorField gradSigmaByVf(fvc::interpolate(delta)*fvc::interpolate(gradSigmaByV));
// 
//     surfaceScalarField tangentForce((gradSigmaByVf - (gradSigmaByVf & nHatfv)*nHatfv) & mesh.Sf());
// 
//     return fvc::interpolate(sigmaK())*fvc::snGrad(alpha1_) + tangentForce;
