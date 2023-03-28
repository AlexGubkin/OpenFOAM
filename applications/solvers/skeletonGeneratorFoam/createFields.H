Info<< "Reading field S\n" << endl;

volScalarField S
(
    IOobject
    (
        "S",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField skeletonMarker
(
    IOobject
    (
        "skeletonMarker",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(S.dimensions(), 0)
);

Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

volScalarField threshold
(
    IOobject
    (
        "threshold",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(transportProperties.lookup("threshold"))
);

Info<< "Reading diffusivity DS\n" << endl;

dimensionedScalar DS
(
    transportProperties.lookup("DS")
);

dimensionedScalar m0
(
    transportProperties.lookup("m0")
);

#include "createFvModels.H"
#include "createFvConstraints.H"
