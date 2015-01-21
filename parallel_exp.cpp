#include <iostream>					// To print things
#include <omp.h>					// For openmp tools
#include <string>					// Just to play with strings because why not...?
#include "time_profiler.cpp"		// To measure running times

int naive_exponentiation (int base, int exponent)
{
	int res = 1;
	for (int i = 0; i<exponent; ++i)
	{
		res = res * base;
	}
	return res;
}

std::string entry (int s1, int s2, int s3)		// Got lazy and didn't want to write this a hundred times...
{
	std::string ent, s1s, s2s, s3s;
	s1s = std::to_string(s1);
	s2s = std::to_string(s2);
	s3s = std::to_string(s3);

	ent = s1s + "," + s2s + "," + s3s;
	return ent;
}

void printMyTuchis (std::string tuchis[16][16])
{
	int tid;
	for (int i=0; i<16; ++i)
	{

		for (int j=0; j<16; ++j)
		{
			//std::cout << tuchis[i][j];
			std::cout <<  "\033[1;31m(" << tuchis[i][j] << ")\033[0m";
			tid = omp_get_thread_num();
			std::cout << tid;
		}
		//std::cout << "         " << i << "\n";
		std::cout << "\n";
	}
}


int main ()
{
	int repetitions = 200;

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


	std::string tuchis[16][16];
	long long blas = 1;
	long long teto = 0;
	
	int tid, resi, resj;						//tid is meant to be the id of the thread that is taking care of that entry, resi and resj are the results of naive exponentiations
	std::string entryij;

	long long t1, t2; 							// initial a final times


	long long rt[6];
	for (int k = 0; k<6; ++k)
	{
		rt[k] = 0;
	}

	for (int k = 0; k<repetitions; ++k){

	
	/*
		///////////////////// --- CONTROL --- \\\\\\\\\\\\\\\\\\\\\
		Here we are profiling measuring the running time of filing 
		a matrix using a single core calling a function that is 
		implemented on the same code and filling it up with smiles!
	*/


	t1 = profiler::rdtsc();
	for (int i=0; i<16; ++i)
	{
		for (int j=0; j<16; ++j)
		{
			tuchis[i][j] = " :) ";
		}
	}
	t2 = profiler::rdtsc();
	std::cout << "The time to fill up tuchis with smiles using a single core is: " << t2-t1 << "\n"; 
	rt[0] = rt[0] + t2 - t1;


	// Print it for fun
	//printMyTuchis(tuchis);

	t1 = profiler::rdtsc();
	for (int i=0; i<16; ++i)
	{
		for (int j=0; j<16; ++j)
		{
			tid = omp_get_thread_num();
			resi = naive_exponentiation(2, i);
			resj = naive_exponentiation(3, j);
			entryij = entry(tid, resi, resj);
			tuchis[i][j] = entryij;
		}
	}
	t2 = profiler::rdtsc();
	std::cout << "The time to fill up tuchis exponentiating things using a single core is: " << t2-t1 << "\n"; 
	rt[1] = rt[1] + t2 - t1;

	// Print it for fun
	//printMyTuchis(tuchis);



	/*
		///////////////// --- FILLING WITH PARALLEL-FOR --- \\\\\\\\\\\\\\\\\
		Just filling a matrix up with a constant using omp parallel for
	*/
	t1 = profiler::rdtsc();
	#pragma omp parallel for
	for (int i = 0; i< 16; ++i)
	{
		for (int j = 0; j<16; ++j)
		{
			tuchis[i][j] = " :) ";
		}
	}
	t2 = profiler::rdtsc();
	std::cout << "The running time of filling tuchis using omp parallel for is " << t2 - t1 << "\n";

	rt[2] = rt[2] + t2 - t1;

	// Print it for fun
	//printMyTuchis (tuchis);

	/*
	t1 = profiler::rdtsc();
	#pragma omp parallel for
	for (int i = 0; i< 16; ++i)
	{
		for (int j = 0; j<16; ++j)
		{
			tid = omp_get_thread_num();
			resi = naive_exponentiation(2, i);
			resj = naive_exponentiation(3, j);
			entryij = entry(tid, resi, resj);
			tuchis[i][j] = entryij;
		}
	}
	t2 = profiler::rdtsc();
	std::cout << "The running time of filling tuchis exponentiating things using omp parallel for is " << t2 - t1 << "\n";

	// Print it for fun
	printMyTuchis (tuchis);

	*/

	t1 = profiler::rdtsc();
	#pragma omp parallel for private(resi, resj, entryij)
	for (int i = 0; i< 16; ++i)
	{
		for (int j = 0; j<16; ++j)
		{
			tid = omp_get_thread_num();
			resi = 1;
			for (int l = 0; l<i; ++l)
			{
				resi = resi * 2;
			}
			resj = 1;
			for (int m = 0; m<j; ++m)
			{
				resj = resj * 3;
			}
			tuchis[i][j] = std::to_string(tid) + "," + std::to_string(resi) + "," + std::to_string(resj);
		}
	}
	t2 = profiler::rdtsc();
	std::cout << "The running time of filling tuchis exponentiating things using omp parallel for is " << t2 - t1 << "\n";

	rt[3] = rt[3] + t2 - t1;

	// Print it for fun
	//printMyTuchis (tuchis);




	/*
		///////////// --- FILLING WITH PARALLEL-SECTIONS --- \\\\\\\\\\\\\\\\
		Filling a matrix up with the same constant using parallel-sections
	*/


	t1 = profiler::rdtsc();
	#pragma omp parallel sections
	{
		#pragma omp section
		{
			for (int i=0; i<8; ++i)
			{

				for (int j=0; j<8; ++j)
				{
					//std::cout <<  "(\033[1;31m" << 2*i << "," << 2*j << ")\033[0m";
					//tid = omp_get_thread_num();
					tuchis[i][j] = " :) ";
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
					//tid = omp_get_thread_num();
					tuchis[i+8][j] = " :) ";
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
					// tid = omp_get_thread_num();
					tuchis[i][j+8] = " :) ";
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
					// tid = omp_get_thread_num();
					tuchis[i+8][j+8] = " :) ";
					//std::cout << tid;
				}
				//std::cout << "\n";
			}
		}
	}		// end of pragma sections	
	t2 = profiler::rdtsc();
	std::cout << "The time to fill up tuchis using parallel-sections is: " << t2-t1 << "\n"; 

	rt[4] = rt[4] + t2 - t1;


	//printMyTuchis(tuchis);

	/*
	t1 = profiler::rdtsc();
	#pragma omp parallel sections
	{
		#pragma omp section
		{
			for (int i=0; i<8; ++i)
			{

				for (int j=0; j<8; ++j)
				{
					tid = omp_get_thread_num();
					resi = naive_exponentiation(2, i);
					resj = naive_exponentiation(3, j);
					entryij = entry(tid, resi, resj);
					tuchis[i][j] = entryij;
				}
			}
		}
		#pragma omp section
		{
			for (int i=0; i<8; ++i)
			{

				for (int j=0; j<8; ++j)
				{
					tid = omp_get_thread_num();
					resi = naive_exponentiation(2, i + 8);
					resj = naive_exponentiation(3, j);
					entryij = entry(tid, resi, resj);
					tuchis[i][j] = entryij;
				}
			}
		}
		#pragma omp section
		{
			for (int i=0; i<8; ++i)
			{

				for (int j=0; j<8; ++j)
				{
					tid = omp_get_thread_num();
					resi = naive_exponentiation(2, i);
					resj = naive_exponentiation(3, j + 8);
					entryij = entry(tid, resi, resj);
					tuchis[i][j] = entryij;
				}
			}
		}
		#pragma omp section
		{
			for (int i=0; i<8; ++i)
			{
				for (int j=0; j<8; ++j)
				{
					tid = omp_get_thread_num();
					resi = naive_exponentiation(2, i + 8);
					resj = naive_exponentiation(3, j + 8);
					entryij = entry(tid, resi, resj);
					tuchis[i][j] = entryij;
				}
			}
		}
	}		// end of pragma sections	
	t2 = profiler::rdtsc();
	std::cout << "The time to fill up tuchis using parallel-sections is: " << t2-t1 << "\n";

	// print it for fun
	printMyTuchis(tuchis);
	*/




	t1 = profiler::rdtsc();
	#pragma omp parallel sections
	{
		#pragma omp section
		{
			for (int i=0; i<8; ++i)
			{
				for (int j=0; j<8; ++j)
				{
					tid = omp_get_thread_num();
					resi = 1;
					for (int l = 0; l<i; ++l)
					{
						resi = resi * 2;
					}
					resj = 1;
					for (int m = 0; m<j; ++m)
					{
						resj = resj * 3;
					}
					tuchis[i][j] = std::to_string(tid) + "," + std::to_string(resi) + "," + std::to_string(resj);
				}
			}
		}
		#pragma omp section
		{
			for (int i=8; i<16; ++i)
			{
				for (int j=0; j<8; ++j)
				{
					tid = omp_get_thread_num();
					resi = 1;
					for (int l = 8; l<i; ++l)
					{
						resi = resi * 2;
					}
					resj = 1;
					for (int m = 0; m<j; ++m)
					{
						resj = resj * 3;
					}
					tuchis[i][j] = std::to_string(tid) + "," + std::to_string(resi) + "," + std::to_string(resj);
				}
			}
		}
		#pragma omp section
		{
			for (int i=0; i<8; ++i)
			{
				for (int j=8; j<16; ++j)
				{
					tid = omp_get_thread_num();
					resi = 1;
					for (int l = 0; l<i; ++l)
					{
						resi = resi * 2;
					}
					resj = 1;
					for (int m = 8; m<j; ++m)
					{
						resj = resj * 3;
					}
					tuchis[i][j] = std::to_string(tid) + "," + std::to_string(resi) + "," + std::to_string(resj);
				}
			}
		}
		#pragma omp section
		{
			for (int i=8; i<16; ++i)
			{
				for (int j=8; j<16; ++j)
				{
					tid = omp_get_thread_num();
					resi = 1;
					for (int l = 8; l<i; ++l)
					{
						resi = resi * 2;
					}
					resj = 1;
					for (int m = 8; m<j; ++m)
					{
						resj = resj * 3;
					}
					tuchis[i][j] = std::to_string(tid) + "," + std::to_string(resi) + "," + std::to_string(resj);
				}
			}
		}
	}		// end of pragma sections	
	t2 = profiler::rdtsc();
	std::cout << "The time to fill up tuchis using parallel-sections is: " << t2-t1 << "\n";

	rt[5] = rt[5] + t2 - t1;

	// print it for fun
	//printMyTuchis(tuchis);


	/*
	//#pragma omp parallel for schedule(static) private(tid)
	for (int i=0; i<16; ++i)
	{

		for (int j=0; j<16; ++j)
		{
			//std::cout << tuchis[i][j];
			std::cout <<  "\033[1;31m(" << tuchis[i][j] << ")\033[0m";
			tid = omp_get_thread_num();
			std::cout << tid;
		}
		//std::cout << "         " << i << "\n";
		std::cout << "\n";
	}
	*/


	}							// End of repetitions for	


	for (int k = 0; k<6; ++k)
	{
		std::cout << "The average for the " << k << "th experiment is " << rt[k]/repetitions << "\n";
	}


}