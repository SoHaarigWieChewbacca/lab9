#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <algorithm>

using namespace std;


void writeToFile(const double* const u, const string s, const double dx, const double xmin, const int N);
void initialize(double* const u, const double dx, const double xmin, const int N);
void ftcs(double* const u0, double* const u1, const double dx, const double V, double dt, const int N);


int main(){

	const double tEnd = 5.0;
	const double V = 1.0; 				// Velocity

	const int N  = 256;
	const double xmin = -10;
	const double xmax =  10;

	const double dx = (xmax-xmin)/(N-1);
	double dt = dx/V;				// Number of steps per subinterval
	const int Na = 10;
	const int Nk = int(tEnd/Na/dt);

	double* u0 = new double[N];
	double* u1 = new double[N];
 

	stringstream strm;

	initialize(u0,dx, xmin,N);

	writeToFile(u0, "u_0", dx, xmin, N);

	for(int i=1; i<=Na; i++){
		for(int j=0; j<Nk; j++){

			ftcs(u0, u1, dx, V, dt, N);
			swap(u0, u1);
		}
		strm.str("");
		strm << "u_" << i;
		writeToFile(u0, strm.str(), dx, xmin, N);
	}


	delete[] u0;
	delete[] u1;
	return 0;
}



void initialize(double* const u, const double dx, const double xmin, const int N){

	for(int i=0; i<N; i++){
		double x = xmin + i*dx;

		if (fabs(x)<=1.0)
			u[i] = 1;
		else
			u[i] = 0;
	}

}


void writeToFile(const double* const u, const string s, const double dx, const double xmin, const int N){

	ofstream out(s.c_str());

	for(int i=0; i<N; i++){
		double x = xmin + i * dx;
		out << x << "\t" << u[i] << endl;
	}

	out.close();
}


void ftcs(double* const u0, double* const u1, const double dx, const double V, double dt, const int N){

	double const C = V*dt/dx;

	u1[0] = u0[0];

	for(int i=1; i<N; i++)
		u1[i] = - C/2.0 * (u0[i+1] - u0[i-1]) + u0[i];
}



