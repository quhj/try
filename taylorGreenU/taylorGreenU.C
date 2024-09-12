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
            "U",
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

		U[cell].x() = Foam::sin(x)*Foam::cos(y);
		U[cell].y() = -Foam::sin(y)*Foam::cos(x);
	}

    U.correctBoundaryConditions();

	U.write();

	Info<< "end" << endl;

    return 0;
}


// ************************************************************************* //
