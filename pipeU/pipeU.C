#include "fvCFD.H"


int main(int argc, char *argv[])
{
    argList::noParallel();
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    volVectorField U
    (
        IOobject
        (
            "U.water",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

	forAll(U, cell)
	{
		const scalar& x = mesh.C()[cell].x();
		const scalar& y = mesh.C()[cell].y();

		U[cell].x() = 0.0;
		U[cell].y() = 0.7 - mag(x - 0.075);
	}

    U.correctBoundaryConditions();

	U.write();

	Info<< "end" << endl;

    return 0;
}


// ************************************************************************* //
