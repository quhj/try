/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
    @randomNumber   An random number generator 

	http://www.cfd-china.com/topic/1393/%E7%94%9F%E6%88%90random%E9%9A%8F%E6%9C%BA%E6%95%B0-of%E4%B8%AD%E8%AF%A5%E6%80%8E%E6%A0%B7%E5%AE%9E%E7%8E%B0
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


#include "Random.H"


using namespace Foam;

int main(int argc, char *argv[])
{
	label seed(1);

	//- Constructor
	Random randomNumber(seed);

	for (int i = 0; i < 100; i++)
	{
		//scalar testNumber(randomNumber.scalar01());
		scalar testNumber(randomNumber.GaussNormal());
       
		Info<< "Random number = " << testNumber << endl;
	}

	return 0;
}

