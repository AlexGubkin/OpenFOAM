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
);

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
);

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
    alpha1
);

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
    alpha1
);

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
    alpha1
);

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
    sqrt(alphaC*(scalar(1.0) - alphaC))
);

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
    sqrt(w1*w1 + SMALL)
);


// Face unit interface normal initialization
// Gradient of alpha
volVectorField gradAlphaSm
(
    "gradAlphaSm",
    fvc::grad(alphaSm)
);
// Interpolated gradient of alpha
surfaceVectorField gradAlphaSmf
(
    "gradAlphaSmf",
    fvc::interpolate(gradAlphaSm)
);

// Interpolated face-gradient of sharpened alpha
surfaceScalarField deltaf
(
    "deltaf",
    fvc::snGrad(alphaSh)
);

// Face unit interface normal
volVectorField nHatv
(
    "nHatv",
    gradAlphaSm/(mag(gradAlphaSm) + mixture.deltaN())
);

// Face unit interface normal
surfaceVectorField nHatfv
(
    "nHatfv",
    gradAlphaSmf/(mag(gradAlphaSmf) + mixture.deltaN())
);

// Face unit interface normal flux
surfaceScalarField nHatf
(
    "nHatf",
    nHatfv & mesh.Sf()
);


// Initialization fields for surface tension force calculation
//- Surface tension
const dictionary&
    phasePropertiesDict = mesh.lookupObject<IOdictionary>("phaseProperties");

autoPtr<surfaceTensionModel> sigmaPtr
(
    surfaceTensionModel::New(phasePropertiesDict, mesh)
);
sigmaPtr->readDict(phasePropertiesDict);


// Surface tension force flux
surfaceScalarField phic
(
    "phic",
    fvc::interpolate(sigmaPtr->sigma())
   *fvc::interpolate(w2*kappa)
   /fvc::interpolate(w2)*deltaf*mesh.magSf()
);
