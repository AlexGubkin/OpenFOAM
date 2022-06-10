Info<< "Reading distance field d\n" << endl;
volScalarField d
(
    IOobject
    (
        "d",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

const polyBoundaryMesh& pbm = mesh.boundaryMesh();

// const labelHashSet patchIDs(mesh.boundary().getPatchIDs());
const labelHashSet patchIDs;

// Get patch ids of walls
// const labelHashSet wallPatchIDs(wallUtils.getPatchIDs<wallPolyPatch>());

wordList dTypes
(
    mesh.boundary().size(),
    zeroGradientFvPatchField<vector>::typeName
);

forAllConstIter(labelHashSet, patchIDs, iter)
{
    dTypes[iter.key()] = fixedValueFvPatchField<vector>::typeName;
}

Info<< "Creation field nd\n" << endl;
volVectorField nd
(
    IOobject
    (
        "nd",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector(dimless, Zero),
    dTypes
);

const fvPatchList& patches = mesh.boundary();
volVectorField::Boundary& ndbf = nd.boundaryFieldRef();

forAllConstIter(labelHashSet, patchIDs, iter)
{
    label patchi = iter.key();
    ndbf[patchi] == -patches[patchi].nf();
}
 

// U = Zero;

// surfaceScalarField phi
// (
//     IOobject
//     (
//         "phi",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::AUTO_WRITE
//     ),
//     fvc::flux(U)
// );
// 
// if (args.optionFound("initialiseUBCs"))
// {
//     U.correctBoundaryConditions();
//     phi = fvc::flux(U);
// }
// 
// 
// // Default name for the pressure field
// word pName("p");
// 
// // Update name of the pressure field from the command-line option
// args.optionReadIfPresent("pName", pName);
// 
// // Infer the pressure BCs from the velocity
// wordList pBCTypes
// (
//     U.boundaryField().size(),
//     fixedValueFvPatchScalarField::typeName
// );
// 
// forAll(U.boundaryField(), patchi)
// {
//     if (U.boundaryField()[patchi].fixesValue())
//     {
//         pBCTypes[patchi] = zeroGradientFvPatchScalarField::typeName;
//     }
// }
// 
// Info<< "Constructing pressure field " << pName << nl << endl;
// volScalarField p
// (
//     IOobject
//     (
//         pName,
//         runTime.timeName(),
//         mesh,
//         IOobject::READ_IF_PRESENT,
//         IOobject::NO_WRITE
//     ),
//     mesh,
//     dimensionedScalar(pName, sqr(dimVelocity), 0),
//     pBCTypes
// );
// 
// // Infer the velocity potential BCs from the pressure
// wordList PhiBCTypes
// (
//     p.boundaryField().size(),
//     zeroGradientFvPatchScalarField::typeName
// );
// 
// forAll(p.boundaryField(), patchi)
// {
//     if (p.boundaryField()[patchi].fixesValue())
//     {
//         PhiBCTypes[patchi] = fixedValueFvPatchScalarField::typeName;
//     }
// }
// 
// Info<< "Constructing velocity potential field Phi\n" << endl;
// volScalarField Phi
// (
//     IOobject
//     (
//         "Phi",
//         runTime.timeName(),
//         mesh,
//         IOobject::READ_IF_PRESENT,
//         IOobject::NO_WRITE
//     ),
//     mesh,
//     dimensionedScalar(dimLength*dimVelocity, 0),
//     PhiBCTypes
// );
// 
// label PhiRefCell = 0;
// scalar PhiRefValue = 0;
// setRefCell
// (
//     Phi,
//     potentialFlow.dict(),
//     PhiRefCell,
//     PhiRefValue
// );
// mesh.setFluxRequired(Phi.name());
// 
// #include "createMRF.H"
