    if
    (
        runTime.writeTime()
//      || ((Foam::mag(m - porosity)).value() < 1e-3)
    )
    {
        volVectorField gradS(fvc::grad(S));

        volScalarField gradSx
        (
            IOobject
            (
                "gradSx",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            gradS.component(vector::X)
        );

        volScalarField gradSy
        (
            IOobject
            (
                "gradSy",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            gradS.component(vector::Y)
        );

        volScalarField gradSz
        (
            IOobject
            (
                "gradSz",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            gradS.component(vector::Z)
        );


        runTime.write();
    }
