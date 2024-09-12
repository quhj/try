{
    volVectorField HbyA(constrainHbyA(rAUc*UcEqn.H(), Uc, pOverRho));

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (
           fvc::flux(HbyA)
         + alphacf*rAUcf*fvc::ddtCorr(Uc, phic)
        )
    );

    if (pOverRho.needReference())
    {
        adjustPhi(phiHbyA, Uc, pOverRho);
    }

    phiHbyA += phicForces;

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(pOverRho, Uc, phiHbyA, rAUcf);

    // Non-orthogonal pressure corrector loop
    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(alphacf*rAUcf, pOverRho)
         ==
            fvc::ddt(alphac) + fvc::div(alphacf*phiHbyA)
        );

        pEqn.setReference(pRefCell, pRefValue);

        pEqn.solve(mesh.solver(pOverRho.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phic = phiHbyA - pEqn.flux()/alphacf;

            pOverRho.relax();

            Uc = HbyA
              + rAUc*fvc::reconstruct((phicForces - pEqn.flux()/alphacf)/rAUcf);
            Uc.correctBoundaryConditions();
        }
    }
}


#include "continuityErrs.H"
