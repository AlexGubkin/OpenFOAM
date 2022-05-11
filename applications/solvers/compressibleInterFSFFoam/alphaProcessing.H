const dictionary& alphaControls = mesh.solverDict(alpha1.name());

const label nAlphaSmoothers(alphaControls.lookupOrDefault<label>("nAlphaSmoothers", 3));
const scalar cSmoothingCoeff(alphaControls.lookupOrDefault<scalar>("cSmoothingCoeff", 0.5));
