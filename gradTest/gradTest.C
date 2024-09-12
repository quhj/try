#include "fvCFD.H"


int main(int argc, char *argv[])
{
    argList::noParallel();
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

	volScalarField p 
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
		mesh
    );

	forAll(p, cellI)
	{
		p[cellI] = cellI;
	}	

    volVectorField gradp 
    (
        IOobject
        (
            "gradp",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
		fvc::grad(p)
    );

	gradp.write();
	p.write();

	Info<< "end" << endl;

    return 0;
}


// ************************************************************************* //
