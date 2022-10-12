phic =
    fvc::interpolate(sigmaPtr->sigma())
   *fvc::interpolate(w2*kappa)
   /fvc::interpolate(w2)*deltaf*mesh.magSf();

gradSigmaf =
    fvc::interpolate(fvc::grad(sigmaPtr->sigma()));

phict =
    deltaf*(gradSigmaf - (gradSigmaf & nHatfv)*nHatfv) & mesh.Sf();
