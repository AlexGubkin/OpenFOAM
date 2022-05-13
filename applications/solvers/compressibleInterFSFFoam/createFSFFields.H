#include "surfaceTensionModel.H"

Info<< "Creation field kappa0\n" << endl;
volScalarField kappa0
(
    IOobject
    (
        "interfaceProperties:kappa0",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimless/dimLength, 0)
)

Info<< "Creation field kappa\n" << endl;
volScalarField kappa
(
    IOobject
    (
        "interfaceProperties:kappa",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimless/dimLength, 0)
)

Info<< "Creation field alphaSm\n" << endl;
volScalarField alphaSm
(
    IOobject
    (
        "alphaSm",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimless, 0)
)

Info<< "Creation field alphaC\n" << endl;
volScalarField alphaC
(
    IOobject
    (
        "alphaC",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimless, 0)
)

Info<< "Creation field alphaSh\n" << endl;
volScalarField alphaSh
(
    IOobject
    (
        "alphaSh",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimless, 0)
)

Info<< "Creation field w1\n" << endl;
volScalarField w1
(
    IOobject
    (
        "w1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimless, 0)
)

Info<< "Creation field w2\n" << endl;
volScalarField w2
(
    IOobject
    (
        "w2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimless, 0)
)

volScalarField pc
(
    IOobject
    (
        "pc",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimPressure, 0),
    zeroGradientFvPatchField<scalar>::typeName
//     zeroGradientFvPatchScalarField::typeName
)

{
    const surfaceVectorField& Sf = mesh.Sf();
    const surfaceScalarField& magSf = mesh.magSf();

    const dictionary& alphaControls = mesh.solverDict(alpha1.name());

    // Cell gradient of alpha
    volScalarField gradAlphaSm("gradAlphaSm", fvc::grad(alphaSm));

    // Interpolated face-gradient of alpha
    surfaceVectorField gradAlphaSmf("gradAlphaSmf", fvc::interpolate(gradAlphaSm));

    // Face unit interface normal
    surfaceVectorField nHatfv("nHatfv", gradAlphaSmf/(mag(gradAlphaSmf) + mixture.deltaN()));

    // Face unit interface normal flux
    surfaceScalarField nHatf("nHatf", nHatfv & Sf);

    // Cell gradient of sharpened alpha
    const volVectorField gradAlphaSharpened(fvc::grad(alphaSharpened_));

    // Interpolated face-gradient of sharpened alpha
    const surfaceVectorField gradAlphaSharpenedf(fvc::interpolate(gradAlphaSharpened));

    const surfaceScalarField deltaf(fvc::snGrad(alphaSharpened_));

    const surfaceScalarField cutFactorf(deltaf/(mag(deltaf) + deltaN_));



    // Sigma
    const volScalarField sigma(sigmaPtr_->sigma());

    // Cell gradient of sigma
    const volVectorField gradSigma(fvc::grad(sigma));

    // Interpolated face-gradient of sigma
    const surfaceVectorField gradSigmaf(fvc::interpolate(gradSigma));

    // Interpolated tangent face-gradient of sigma
    const surfaceScalarField tangentGradSigmaf((gradSigmaf - (gradSigmaf & nHatfv)*nHatfv) & Sf/magSf);

    const surfaceScalarField capillaryForceFlux(fvc::interpolate(w2_*sigmaK())/fvc::interpolate(w2_)*deltaf*magSf);
    surfaceScalarField capillaryForceFluxFiltered(capillaryForceFlux);

    //- Surface tension
    autoPtr<surfaceTensionModel> sigmaPtr;
    // sigmaPtr->readDict(transportPropertiesDict);
    // Sigma
//     const volScalarField sigma(sigmaPtr->sigma());
}
