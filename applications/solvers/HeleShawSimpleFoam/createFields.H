Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading Hele-Shaw model parameters\n" << endl;
IOdictionary HeleShawModel
(
    IOobject
    (
        "HeleShawModel",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

Info<< "Reading gap thickness b\n" << endl;
dimensionedScalar b(HeleShawModel.lookup("b"));

dimensionedScalar advectionFactor
(
    dimless,
    1.2
);

dimensionedScalar rb
(
    dimless/dimLength,
    1/b.value()
);

dimensionedScalar resistanceFactor
(
    dimless/(dimLength*dimLength),
    12.0*rb.value()*rb.value()
);

#include "createPhi.H"

pressureReference pressureReference(p, simple.dict());

mesh.schemes().setFluxRequired(p.name());

autoPtr<viscosityModel> viscosity(viscosityModel::New(mesh));

autoPtr<incompressible::momentumTransportModel> turbulence
(
    incompressible::momentumTransportModel::New(U, phi, viscosity)
);

#include "createMRF.H"
#include "createFvModels.H"
#include "createFvConstraints.H"
