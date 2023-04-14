    MRF.correctBoundaryVelocity(U);

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
      + MRF.DDt(rho, U)
      + turbulence.divDevTau(rho, U)
     ==
       phaseChange.SU(rho, rhoPhi, U)
     + fvModels.source(rho, U)
    );

    UEqn.relax();

    fvConstraints.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                  + fvc::snGrad(pc)
                ) * mesh.magSf()
              + phic
            )
          + 0.2*fvc::div(rhoPhi, U)
          - 12.0*mixture.mu()*U
        );

        fvConstraints.constrain(U);
    }
