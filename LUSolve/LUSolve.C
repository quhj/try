//#include "scalar.H"
#include "scalarMatrices.H"

using namespace Foam;

int main(int argc, char *argv[])
{
	scalar* M_ = new scalar[5];

	M_[0] = (558292.0573648899);
	M_[1] = (2233.1682294595607);
	M_[2] = (9.133658058489596);
	M_[3] = (0.038197186342054885);
	M_[4] = (0.00016333577394690517);

	//M_[0] = 1.0;
	//M_[1] = 2.0; 
	//M_[2] = 3.0;
	//M_[3] = 4.0;
	//M_[4] = 5.0;

	scalar V(0.5);
	scalar U(0.0);

	SquareMatrix<scalar>  momentsMatrix(3); 

	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			momentsMatrix(i, j) = M_[i + j];
		}
	}

	scalarDiagonalMatrix rhsU(3, 0);
	scalarDiagonalMatrix rhsV(3, 0);

	for(int i = 0; i < 3; i++)
	{
		rhsU[i] = M_[i]*U;
		rhsV[i] = M_[i]*V;
	}

	Info<< "Creating SquareMatrix momentsMatrix(3)" << nl
		<< momentsMatrix << nl << nl;

	Info<< "Creating scalarDiagonalMatrix rhsU(3, 0)" << nl
		<< rhsU << nl << nl;

	Info<< "Creating scalarDiagonalMatrix rhsV(3, 0)" << nl
		<< rhsV << nl << nl;

	LUsolve(momentsMatrix, rhsV);

	Info<< "Results: " << rhsV << nl;

	delete [] M_;

	return 0;
}


// ************************************************************************* //
