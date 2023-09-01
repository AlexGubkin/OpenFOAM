gradSigma = fvc::grad(sigmaPtr->sigma());

fct = mag(fvc::grad(alpha1))*(gradSigma - (gradSigma & nHatSmv)*nHatSmv);

phict = fvc::interpolate(fct) & mesh.Sf();
