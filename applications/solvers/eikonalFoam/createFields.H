// const polyBoundaryMesh& pbm = mesh.boundaryMesh();

// Get patch ids of walls
const labelHashSet patchIDs(mesh.boundaryMesh().findPatchIDs<wallPolyPatch>());

wordList dTypes
(
    mesh.boundary().size(),
    zeroGradientFvPatchField<scalar>::typeName
);

wordList ndTypes
(
    mesh.boundary().size(),
    zeroGradientFvPatchField<vector>::typeName
);

forAllConstIter(labelHashSet, patchIDs, iter)
{
    dTypes[iter.key()] = fixedValueFvPatchField<scalar>::typeName;
    ndTypes[iter.key()] = fixedValueFvPatchField<vector>::typeName;
}

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
    ndTypes
);

const fvPatchList& patches = mesh.boundary();
volVectorField::Boundary& ndbf = nd.boundaryFieldRef();

forAllConstIter(labelHashSet, patchIDs, iter)
{
    label patchi = iter.key();
    ndbf[patchi] == -patches[patchi].nf();
}




// // Infer the pressure BCs from the velocity
// wordList dTypes
// (
//     mesh.boundary().size(),
//     zeroGradientFvPatchScalarField::typeName
// );
// 
// forAll(U.boundaryField(), patchi)
// {
//     if (U.boundaryField()[patchi].fixesValue())
//     {
//         pBCTypes[patchi] = zeroGradientFvPatchScalarField::typeName;
//     }
// }
