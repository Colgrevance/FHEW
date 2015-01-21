#include <iostream>
#include <omp.h>

int main ()
{
	/*
	#pragma omp parallel sections
	{
		std::cout << "asdasdasdas\n";

		#pragma omp section
		{
			std::cout << "Yo mama dude\n";
			for (int i = 0; i<50; ++i)
			{
				std::cout << i;
			}
			std::cout << "\n";
		}
	}
	std::cout << "bababababababa\n\n";
	*/

	long long tuchis[20][20];
	long long blas = 1;
	long long teto = 0;
	
	int tid;

	/*
	//#pragma omp parallel for
	for (int i=0; i<15; ++i)
	{
		for (int j=0; j<15; ++j)
		{
			tuchis[i][j] = blas+teto;
			blas = blas * 2;
		}
		blas = 1;
		teto = teto + 5;
	}
	*/





	#pragma omp parallel sections
	{
		#pragma omp section
		{
			for (int i=0; i<8; ++i)
			{

				for (int j=0; j<8; ++j)
				{
					//std::cout <<  "(\033[1;31m" << 2*i << "," << 2*j << ")\033[0m";
					tid = omp_get_thread_num();
					tuchis[i][j] = tid;
					//std::cout << tid;
				}
				//std::cout << "\n";
			}
		}
		#pragma omp section
		{
			for (int i=0; i<8; ++i)
			{

				for (int j=0; j<8; ++j)
				{
					//std::cout <<  "(\033[1;31m" << 2*i + 1 << "," << 2*j << ")\033[0m";
					tid = omp_get_thread_num();
					tuchis[i+8][j] = tid;
					//std::cout << tid;
				}
				//std::cout << "\n";
			}
		}
		#pragma omp section
		{
			for (int i=0; i<8; ++i)
			{

				for (int j=0; j<8; ++j)
				{
					//std::cout <<  "(\033[1;31m" << 2*i << "," << 2*j +1 << ")\033[0m";
					tid = omp_get_thread_num();
					tuchis[i][j+8] = tid;
					//std::cout << tid;
				}
				//std::cout << "\n";
			}
		}
		#pragma omp section
		{
			for (int i=0; i<8; ++i)
			{

				for (int j=0; j<8; ++j)
				{
					//std::cout <<  "(\033[1;31m" << 2*i + 1 << "," << 2*j +1 << ")\033[0m";
					tid = omp_get_thread_num();
					tuchis[i+8][j+8] = tid;
					//std::cout << tid;
				}
				//std::cout << "\n";
			}
		}
	}		// end of pragma sections	

	
	//#pragma omp parallel for schedule(static) private(tid)
	for (int i=0; i<16; ++i)
	{

		for (int j=0; j<16; ++j)
		{
			//std::cout << tuchis[i][j];
			std::cout <<  "(\033[1;31m" << tuchis[i][j] << ")\033[0m";
			tid = omp_get_thread_num();
			std::cout << tid;
		}
		//std::cout << "         " << i << "\n";
		std::cout << "\n";
	}
	

}