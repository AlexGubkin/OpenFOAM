// Get patch ids of walls
const labelHashSet patchIDs(mesh.boundaryMesh().findPatchIDs<wallPolyPatch>());

wordList dTypes
(
    mesh.boundary().size(),
    zeroGradientFvPatchField<scalar>::typeName
);

wordList d2Types
(
    mesh.boundary().size(),
    zeroGradientFvPatchField<scalar>::typeName
);

forAllConstIter(labelHashSet, patchIDs, iter)
{
    dTypes[iter.key()] = fixedValueFvPatchField<scalar>::typeName;
}


Info<< "Creation field alpha\n" << endl;
volScalarField alpha
(
    IOobject
    (
        "alpha",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volVectorField gradAlpha
(
    IOobject
    (
        "gradAlpha",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::grad(alpha)
);

gradAlpha.write();


Info<< "Creation distance field d\n" << endl;
volScalarField d
(
    IOobject
    (
        "d",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar(dimLength, small),
    dTypes
);

d = dimensionedScalar(dimLength, 1.0);
forAll(d, celli)
{
    if (mag(gradAlpha[celli]) > small)
    {
        d[celli] = small;
    }
}

Info<< "Creation field gradd\n" << endl;
volVectorField gradd
(
    IOobject
    (
        "gradd",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector(dimless, Zero)
);

Info<< "Creation distance field d2\n" << endl;
volScalarField d2
(
    IOobject
    (
        "d2",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar(dimLength, small),
    d2Types
);

d2 = dimensionedScalar(dimLength, 1.0);

Info<< "Creation field gradd2\n" << endl;
volVectorField gradd2
(
    IOobject
    (
        "gradd2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector(dimless, Zero)
);

Info<< "Creation field dPsi\n" << endl;
volScalarField dPsi
(
    IOobject
    (
        "dPsi",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar(sqr(dimLength), 0),
    d.boundaryFieldRef().types()
);

Info<< "Creation field graddPsi\n" << endl;
volVectorField graddPsi
(
    IOobject
    (
        "graddPsi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector(dimLength, Zero)
);
