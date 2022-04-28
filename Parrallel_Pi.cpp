/*/------------------- Lab 7-------------------\*\
 Programmers: Jarod Teague, Jake Shockey
 Date: 8/29/2018
 Summary:
	This program calcultaes pi using the Gregory
	Leibnitz series approximation for Pi. The proof
	Substituted for a For Loop within the program
	is as follows:
	
               n
     π  =  4 * Σ (-1)^k / 2k + 1
             k = 1
	
	
compiler command: g++ -fopenmp Parrallel_Pi.cpp Timer.cpp -o Pi.exe
\*\--------------------------------------------/*/

#include <omp.h>
#include <iomanip>
#include <iostream>
#include "Timer.h"// timer.h does not calculate time with an accuracy in milliseconds.

#pragma comment( lib, "winmm.lib")

using namespace std;

void calculate_pi(int thread_count, int total_terms, double& pi_agg, double& time_agg);

int main()
{
	int thread_count = 0;
	int total_terms = 0;
	double pi_agg = 0;
	double time_agg = 0;
	
	
	int run_cycles = 0;
	bool flag = true;
	
	while(flag)// all of this extra code in the front is here merely to allow for faster value testing, giving the user the option to specify the number of cores and level of accuracy
	{
		int a = 0,b = 0,c = 0;
		
		if(thread_count > 0)
		{
			a = run_cycles;
			b = thread_count;
			c = total_terms;
			
			run_cycles = 0;
			thread_count = 0;
			total_terms = 0;
		}
		
		
		while(run_cycles <= 0 || run_cycles > 5)
		{
			if(!flag)
			{
				cout << "ERROR: invalid number,";
				if(run_cycles > 0)
					cout <<  "calculations may not exceed 5" << endl; // just to prevent the user from accidentally entering a number, that will take too long to calculate.	
			}
			
			cout <<"(previous value: " << setprecision(0) << a << ")" << " specifiy the number of calculations: ";
			cin >> run_cycles;
			
			flag = false;
		}
		
		cout << endl;
		
		while(thread_count <= 0)
		{
			if(flag)
				cout << "ERROR: invalid number" << endl;
			
			cout <<"(previous value: " << b << ")" << " specifiy the number of threads to use: ";
			cin >> thread_count;
			
			flag = true;
		}
		
		cout << endl;
		
		while(total_terms <= 0)
		{
			if(!flag)
				cout << "ERROR: invalid number" << endl;
			
			cout <<"(previous value: " << c << ")" << " specifiy the number of terms for accuracy: ";
			cin >> total_terms;
			
			flag = false;
		}
		
		cout << endl;
		
		for(int i = 0; i < run_cycles; i++)
		{
			calculate_pi(thread_count, total_terms, pi_agg, time_agg);
			if(i > 1)
			{
				//pi_agg /= 2; // callculates average with every loop so it does not lose 
			}
		}
		
		pi_agg /= run_cycles;
		time_agg /= run_cycles;
		
		cout << run_cycles << " calculations completed." << endl;
		cout << "average PI: " << setprecision(12) << pi_agg << endl;
		cout << "average time elapsed: " << setprecision(12) << time_agg<< "s " << endl << endl;
		
		cout<< "continue?(1 = yes, 0 = no): ";
		pi_agg = 0;
		time_agg = 0;
		cin >> flag;
	}
	
	return 0;
}

void calculate_pi(int thread_count, int total_terms, double& pi_agg, double& time_agg) // this is where the actual calculation is taking place
{
	int my_rank = 0;
	double sum = 0.0;
	time_t myStart = getTime();
	
	 /***********************************************/
	/*************** Multi Threading ***************/
	#pragma omp parallel num_threads(thread_count)// beggining of the openmp statement, calls the specified number of threads
	{
		my_rank = omp_get_thread_num();// obtains the rank of each thread
		int num_parts = total_terms/thread_count;// divides the number of terms between the threads
		int my_first_i = num_parts * my_rank;//allows each thread to calculate their starting point
		int my_last_i;// allows each thread to calculate their ending point
		
		if (my_rank == thread_count - 1)// tells the first thread to start at the beggining
		{
			my_last_i = total_terms;
		} 
		
		else
		{ 
			my_last_i = my_first_i + num_parts;// all of the rest will use their augmented starts
		}

		double factor;// part of the leibnitz notation, causes the value to oscillate
		if (my_first_i % 2 == 0.0) 
		{
			factor = 1.0;
		}
		else 
		{
			factor = -1.0;
		}
		
		double my_sum = 0.0;// the summation itself as mentioned above
		for(int i = my_first_i ; i < my_last_i ; i++, factor = -factor)// the function provided was not the leibnitz notation for pi, then again it never said it was either...
		{
			my_sum += 4*((pow(factor,i))/(2.0 * i + 1.0));
		}
		
		#pragma omp critical// end of the open mp statement
		sum += my_sum;
	}
	
	time_t myEnd = getTime();
	
	 /*********************************************/
	/***************** Result ********************/
	double time = totalTime(myStart, myEnd);
	pi_agg += sum;
	time_agg += time;
}