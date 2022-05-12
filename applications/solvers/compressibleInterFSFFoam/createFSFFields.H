#include "surfaceTensionModel.H"
//- Surface tension
autoPtr<surfaceTensionModel> sigmaPtr;
// sigmaPtr->readDict(transportPropertiesDict);
// Sigma
const volScalarField sigma(sigmaPtr->sigma());

Info<< "Creation field K0\n" << endl;
volScalarField K0
(
    IOobject
    (
        "interfaceProperties:K0",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimless/dimLength, 0)
)

Info<< "Creation field K\n" << endl;
volScalarField K
(
    IOobject
    (
        "interfaceProperties:K",
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

const surfaceVectorField& Sf = mesh.Sf();
const surfaceScalarField& magSf = mesh.magSf();

const dictionary& alphaControls = mesh.solverDict(alpha1.name());
