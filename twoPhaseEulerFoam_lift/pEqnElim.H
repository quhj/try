surfaceScalarField alpha1f("alpha1f", fvc::interpolate(alpha1));
surfaceScalarField alpha2f("alpha2f", scalar(1) - alpha1f);

rAU1 = 1.0/U1Eqn.A();
rAU2 = 1.0/U2Eqn.A();

surfaceScalarField alpharAU1f(fvc::interpolate(alpha1*rAU1));
surfaceScalarField alpharAU2f(fvc::interpolate(alpha2*rAU2));

// --- Pressure corrector loop
while (pimple.correct())
{
    volVectorField HbyA1
    (
        IOobject::groupName("HbyA", phase1.name()),
        U1
    );
    HbyA1 = rAU1*U1Eqn.H();

    volVectorField HbyA2
    (
        IOobject::groupName("HbyA", phase2.name()),
        U2
    );
    HbyA2 = rAU2*U2Eqn.H();

    // Face force fluxes
    tmp<surfaceScalarField> phiF1;
    tmp<surfaceScalarField> phiF2;

    // Turbulent diffusion, particle-pressure lift and wall-lubrication fluxes
    {
        volScalarField turbulentDiffusivity(fluid.turbulentDiffusivity());
        volVectorField liftForce(fluid.liftForce());
        volVectorField wallLubricationForce(fluid.wallLubricationForce());
        surfaceScalarField snGradAlpha1(fvc::snGrad(alpha1)*mesh.magSf());

        phiF1 =
        (
            fvc::interpolate
            (
                rAU1*(turbulentDiffusivity + phase1.turbulence().pPrime())
            )*snGradAlpha1

          + (
                fvc::interpolate(rAU1*(wallLubricationForce + liftForce))
              & mesh.Sf()
            )
        );

        phiF2 =
        (
          - fvc::interpolate
            (
                rAU2*(turbulentDiffusivity + phase2.turbulence().pPrime())
            )*snGradAlpha1

          - (
                fvc::interpolate(rAU2*(wallLubricationForce + liftForce))
              & mesh.Sf()
            )
        );
    }

    // Mean density for buoyancy force and p_rgh -> p
    volScalarField rho("rho", fluid.rho());

    // Add Buoyancy force
    {
        surfaceScalarField ghSnGradRho(ghf*fvc::snGrad(rho)*mesh.magSf());

        phiF1() +=
            alpharAU1f
           *(
                ghSnGradRho
              - alpha2f*fvc::interpolate(rho1 - rho2)*(g & mesh.Sf())
            );

        phiF2() +=
            alpharAU2f
           *(
                ghSnGradRho
              - alpha1f*fvc::interpolate(rho2 - rho1)*(g & mesh.Sf())
            );
    }

    // ddtPhiCorr filter -- only apply in pure(ish) phases
    surfaceScalarField alpha1fBar(fvc::interpolate(fvc::average(alpha1f)));
    surfaceScalarField phiCorrCoeff1(pos(alpha1fBar - 0.99));
    surfaceScalarField phiCorrCoeff2(pos(0.01 - alpha1fBar));

    // Set ddtPhiCorr to 0 on non-coupled boundaries
    forAll(mesh.boundary(), patchi)
    {
        if
        (
            !mesh.boundary()[patchi].coupled()
         || isA<cyclicAMIFvPatch>(mesh.boundary()[patchi])
        )
        {
            phiCorrCoeff1.boundaryField()[patchi] = 0;
            phiCorrCoeff2.boundaryField()[patchi] = 0;
        }
    }

    // Phase-1 predicted flux
    surfaceScalarField phiHbyA1
    (
        IOobject::groupName("phiHbyA", phase1.name()),
        (fvc::interpolate(HbyA1) & mesh.Sf())
      + phiCorrCoeff1*fvc::interpolate(alpha1.oldTime()*rho1.oldTime()*rAU1)
       *(
            mrfZones.absolute(phi1.oldTime())
          - (fvc::interpolate(U1.oldTime()) & mesh.Sf())
        )/runTime.deltaT()
      - phiF1()
    );

    // Phase-2 predicted flux
    surfaceScalarField phiHbyA2
    (
        IOobject::groupName("phiHbyA", phase2.name()),
        (fvc::interpolate(HbyA2) & mesh.Sf())
      + phiCorrCoeff2*fvc::interpolate(alpha2.oldTime()*rho2.oldTime()*rAU2)
       *(
            mrfZones.absolute(phi2.oldTime())
          - (fvc::interpolate(U2.oldTime()) & mesh.Sf())
        )/runTime.deltaT()
      - phiF2()
    );

    // Face-drag coefficients
    surfaceScalarField D1f(fvc::interpolate(rAU1*dragCoeff));
    surfaceScalarField D2f(fvc::interpolate(rAU2*dragCoeff));

    // Construct the mean predicted flux
    // including explicit drag contributions based on absolute fluxes
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        alpha1f*(phiHbyA1 + D1f*mrfZones.absolute(phi2))
      + alpha2f*(phiHbyA2 + D2f*mrfZones.absolute(phi1))
    );
    mrfZones.makeRelative(phiHbyA);

    // Construct pressure "diffusivity"
    surfaceScalarField rAUf
    (
        "rAUf",
        mag(alpha1f*alpharAU1f + alpha2f*alpharAU2f)
    );

    // Update the fixedFluxPressure BCs to ensure flux consistency
    setSnGrad<fixedFluxPressureFvPatchScalarField>
    (
        p_rgh.boundaryField(),
        (
            phiHbyA.boundaryField()
          - mrfZones.relative
            (
                alpha1f.boundaryField()
               *(mesh.Sf().boundaryField() & U1.boundaryField())
              + alpha2f.boundaryField()
               *(mesh.Sf().boundaryField() & U2.boundaryField())
            )
        )/(mesh.magSf().boundaryField()*rAUf.boundaryField())
    );

    tmp<fvScalarMatrix> pEqnComp1;
    tmp<fvScalarMatrix> pEqnComp2;

    // Construct the compressibility parts of the pressure equation
    if (pimple.transonic())
    {
        surfaceScalarField phid1
        (
            IOobject::groupName("phid", phase1.name()),
            fvc::interpolate(psi1)*phi1
        );
        surfaceScalarField phid2
        (
            IOobject::groupName("phid", phase2.name()),
            fvc::interpolate(psi2)*phi2
        );

        pEqnComp1 =
            (
                contErr1
              - fvc::Sp(fvc::ddt(alpha1) + fvc::div(alphaPhi1), rho1)
            )/rho1
          + (alpha1/rho1)*correction
            (
                psi1*fvm::ddt(p_rgh)
              + fvm::div(phid1, p_rgh) - fvm::Sp(fvc::div(phid1), p_rgh)
            );
        deleteDemandDrivenData(pEqnComp1().faceFluxCorrectionPtr());
        pEqnComp1().relax();

        pEqnComp2 =
            (
                contErr2
              - fvc::Sp(fvc::ddt(alpha2) + fvc::div(alphaPhi2), rho2)
            )/rho2
          + (alpha2/rho2)*correction
            (
                psi2*fvm::ddt(p_rgh)
              + fvm::div(phid2, p_rgh) - fvm::Sp(fvc::div(phid2), p_rgh)
            );
        deleteDemandDrivenData(pEqnComp2().faceFluxCorrectionPtr());
        pEqnComp2().relax();
    }
    else
    {
        pEqnComp1 =
            (
                contErr1
              - fvc::Sp(fvc::ddt(alpha1) + fvc::div(alphaPhi1), rho1)
            )/rho1
          + (alpha1*psi1/rho1)*correction(fvm::ddt(p_rgh));

        pEqnComp2 =
            (
                contErr2
              - fvc::Sp(fvc::ddt(alpha2) + fvc::div(alphaPhi2), rho2)
            )/rho2
          + (alpha2*psi2/rho2)*correction(fvm::ddt(p_rgh));
    }

    // Cache p prior to solve for density update
    volScalarField p_rgh_0(p_rgh);

    // Iterate over the pressure equation to correct for non-orthogonality
    while (pimple.correctNonOrthogonal())
    {
        // Construct the transport part of the pressure equation
        fvScalarMatrix pEqnIncomp
        (
            fvc::div(phiHbyA)
          - fvm::laplacian(rAUf, p_rgh)
        );

        solve
        (
            pEqnComp1() + pEqnComp2() + pEqnIncomp,
            mesh.solver(p_rgh.select(pimple.finalInnerIter()))
        );

        // Correct fluxes and velocities on last non-orthogonal iteration
        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA + pEqnIncomp.flux();

            surfaceScalarField mSfGradp("mSfGradp", pEqnIncomp.flux()/rAUf);

            // Partial-elimination phase-flux corrector
            {
                surfaceScalarField phi1s
                (
                    phiHbyA1 + alpharAU1f*mSfGradp
                );

                surfaceScalarField phi2s
                (
                    phiHbyA2 + alpharAU2f*mSfGradp
                );

                surfaceScalarField phir
                (
                    ((phi1s + D1f*phi2s) - (phi2s + D2f*phi1s))/(1 - D1f*D2f)
                );

                phi1.boundaryField() ==
                    mrfZones.relative
                    (
                        mesh.Sf().boundaryField() & U1.boundaryField()
                    );
                phi1 = phi + alpha2f*phir;

                phi2.boundaryField() ==
                    mrfZones.relative
                    (
                        mesh.Sf().boundaryField() & U2.boundaryField()
                    );
                phi2 = phi - alpha1f*phir;
            }

            // Compressibility correction for phase-fraction equations
            fluid.dgdt() =
            (
                alpha1*(pEqnComp2 & p_rgh)
              - alpha2*(pEqnComp1 & p_rgh)
            );

            // Optionally relax pressure for velocity correction
            p_rgh.relax();

            // Update the static pressure
            p = max(p_rgh + rho*gh, pMin);
            p_rgh = p - rho*gh;

            mSfGradp = pEqnIncomp.flux()/rAUf;

            // Partial-elimination phase-velocity corrector
            {
                volVectorField U1s
                (
                    HbyA1 + fvc::reconstruct(alpharAU1f*mSfGradp - phiF1)
                );

                volVectorField U2s
                (
                    HbyA2 + fvc::reconstruct(alpharAU2f*mSfGradp - phiF2)
                );

                volScalarField D1(rAU1*dragCoeff);
                volScalarField D2(rAU2*dragCoeff);

                U = alpha1*(U1s + D1*U2) + alpha2*(U2s + D2*U1);
                volVectorField Ur(((1 - D2)*U1s - (1 - D1)*U2s)/(1 - D1*D2));

                U1 = U + alpha2*Ur;
                U1.correctBoundaryConditions();
                fvOptions.correct(U1);

                U2 = U - alpha1*Ur;
                U2.correctBoundaryConditions();
                fvOptions.correct(U2);

                U = fluid.U();
            }
        }
    }

    // Update densities from change in p
    rho1 += psi1*(p_rgh - p_rgh_0);
    rho2 += psi2*(p_rgh - p_rgh_0);

    // Update the phase kinetic energies
    K1 = 0.5*magSqr(U1);
    K2 = 0.5*magSqr(U2);

    // Update the pressure time-derivative if required
    if (thermo1.dpdt() || thermo2.dpdt())
    {
        dpdt = fvc::ddt(p);
    }
}
