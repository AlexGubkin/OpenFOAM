const dictionary& eikonalDict
(
    mesh.solutionDict().subDict("eikonal")
);

// const int nNonOrthCorr
// (
//     eikonalDict.lookupOrDefault<int>("nNonOrthogonalCorrectors", 0)
// );

// Sub-dictionary of coefficients
const dictionary& coeffs
(
    eikonalDict.optionalSubDict("eikonalCoeffs")
);

// Diffusion coefficient multiplying d*laplacian(d)
scalar epsilon(coeffs.lookupOrDefault<scalar>("epsilon", 0.1));

// Convergence tolerance for the iterations of the advection-diffusion
// equation to correct the distance-to-patch and normal-to-patch fields
// scalar tolerance(coeffs.lookupOrDefault<scalar>("tolerance", 1e-3));

// Maximum number of iterations of the advection-diffusion equation
// to correct the distance-to-patch and normal-to-patch fields
// int maxIter(coeffs.lookupOrDefault<int>("maxIter", 10));
