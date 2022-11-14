gradSigma = fvc::grad(sigmaPtr->sigma());

fct = mag(fvc::grad(alpha1))*(gradSigma - (gradSigma & nHatv)*nHatv);

phict = fvc::interpolate(fct) & mesh.Sf();
