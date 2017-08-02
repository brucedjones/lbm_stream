// The MIT License (MIT)
// 
// Copyright (c) 2017 Dr Bruce Jones
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include <fstream>
#include <string>
#include <stdlib.h>
#include "time.h"

#include "lbm_stream.h"

// Comment these lines to toggle between 
//	+ Origin Streaming / Two Lattice
//  + Output Files / No Output Files
// X-velocity is output to a binary file in IJ order

// Lattice Configuration
#define DIM 2
#define Q 9

// Lattice Constants
const int Lx = 20000;
const int Ly = 1000;
const int maxT = 100;
const double tau = 1.0;
const double force[DIM] = {0.00001, 0.0};
const int output_frequency = 100;
const double omega[Q] = { 4.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0 };
const int e[DIM][Q] = { { 0, 1, 0, -1, 0, 1, -1, -1, 1 },
						{ 0, 0, 1, 0, -1, 1, 1, -1, -1 } };
const int opp[Q] = { 0, 3, 4, 1, 2, 7, 8, 5, 6 };

// Data Structures
double *f, *ux;

#ifdef TWO_LATTICE
	double *f_prev;
#endif

void writeOutput(int t){
	char fname [50];
	sprintf(fname,"output_%05d.lbm",t);
	std::ofstream myFile;
	myFile.open (fname, std::ios::out | std::ios::binary);
	myFile.write ((char*)ux, Lx*Ly*sizeof(double));
}

int main(int argc, char* argv[])
{
	printf("%d %d %d\n", Q, Lx, Ly);
	// Allocate Memory
	f = new double[Q*Lx*Ly];
	#ifdef TWO_LATTICE
		printf("worked\n");
		f_prev = new double[Q*Lx*Ly];
	#endif
	#ifdef OUTPUT
		ux = new double[Lx*Ly];
	#endif

	// initialise domain to rest
	int ix;
	for(int j = 0; j<Ly; j++) {
		for(int i = 0; i<Lx; i++) {
			for(int q = 0; q<Q; q++) {
				//ix = q+(i*Q)+(j*Q*Lx);
				ix = (i) + (j*Lx) + (q*Q*Lx);
				f[ix] = omega[q];
				#ifdef TWO_LATTICE
					f_prev[ix] = omega[q];
				#endif
			}
		}
	}

	// Start the timer
	clock_t start_time = clock();

	// Main Time Loop
	for(int t = 0; t<maxT;t++) {
		// Loop over Y
		#pragma omp parallel for
		for(int j = 0; j<Ly; j++) {
			// Loop over X
			for(int i = 0; i<Lx; i++) {

				double f_eq[Q], f_now[Q];

				#ifdef ORIGIN_STREAMING
				// Load PDF's
				int idx;
				for(int q = 0; q<Q; q++) {
					// Compute PDF origin ijk
					int i_o = (abs(i - (t*e[0][q])) % Lx);

					// Compute PDF origin ijk
					int j_o = (abs(j - (t*e[1][q])) % Ly);

					//idx = q + (i_o*Q) + (j_o*Q*Lx);
					idx = (i_o)+(j_o*Lx) + (q*Q*Lx);
					f_now[q] = f[idx];
				}
				#endif
				
				#ifdef TWO_LATTICE
					for(int q = 0; q<Q; q++) {
						//int idx = q+(i*Q)+(j*Q*Lx);
						int idx = (i)+(j*Lx) + (q*Q*Lx);
						f_now[q] = f_prev[idx];
 					}
				#endif

				// Compute Macroscopic values
				double rho = 0;
				double u[DIM] = {0.0,0.0};

				for(int q = 0; q<Q; q++) {
					rho += f_now[q];
					for(int d = 0; d<DIM; d++) {
						u[d] += e[d][q]*f_now[q];
					}
				}

				// Add force contribution to momentum and divide by density to get velocity
				for(int d = 0; d<DIM; d++) {
					u[d] = (u[d]+(force[d]/2.0))/rho;
				}

				#ifdef OUTPUT
					// Store Velocity
					if(t%output_frequency ==0) {
						int idu = i+(j*Lx);
						ux[idu] = u[0];
					}
				#endif
				
				// Check for node for boundary conditions
				if(j==0 || j == Ly-1)
				{
					// Bounceback Boundary
					for(int q=0;q<Q;q++)
					{	
						#ifdef ORIGIN_STREAMING
						// Compute PDF origin ijk
						int i_o = (abs(i - (t*e[0][q])) % Lx);

						// Compute PDF origin ijk
						int j_o = (abs(j - (t*e[1][q])) % Ly);

						//int idxo = opp[q] + (i_o*Q) + (j_o*Q*Lx);
						int idxo = (i_o)+(j_o*Lx) + (opp[q]*Q*Lx);
							f[idxo] = f_now[opp[q]];
						#endif

						#ifdef TWO_LATTICE
							int i_target = (i+e[0][q]);
							int j_target = (j+e[1][q]);
							if(i_target<0) i_target = Lx-1;
							if(i_target>(Lx-1)) i_target = 0;
							if(j_target<0) j_target = Ly-1;
							if(j_target>(Ly-1)) j_target = 0;
							//int idx = q+(i_target*Q)+(j_target*Q*Lx);
							int idx = (i_target)+(j_target*Lx) + (q*Q*Lx);
							f[idx] = f_now[opp[q]];
						#endif
					}
				} else // Do standard LBM with Guo Force
				{					
					double F_coeff[DIM], force_term[Q];

					// Compute first helper term for equilibrium distribution
					double u_square, eu;
					u_square = 0.0;
					eu = 0.0;
					for(int d = 0; d<DIM; d++)
					{
						u_square += (u[d]*u[d]);
					}
					u_square*=1.5;
					
					for(int q=0;q<Q;q++)
					{
						// Compute second helper term for equilibrium distribution
						eu = 0.0;
						for(int d = 0; d<DIM; d++)
						{
							eu += (e[d][q]*u[d]);
						}
						eu*=3.0;
					
						// Compute equilibrium distribution
						f_eq[q] = rho*omega[q]*(1.0+eu+(0.5*eu*eu)-u_square);
					
						// Compute force coefficient
						for(int d = 0; d<DIM; d++)
						{
							F_coeff[d] = omega[q]*(1-(1/(2*(tau))))*(((e[d][q]-u[d])*3.0)+(e[d][q]*9.0*((e[0][q]*u[0])+(e[1][q]*u[1]))));
						}
					
						// Compute force term
						force_term[q] = 0;
						for(int d = 0; d<DIM; d++)
						{
							force_term[q] += F_coeff[d]*force[d];
						}
					
						// Execute collision
						#ifdef ORIGIN_STREAMING
							// Compute PDF origin ijk
							int i_o = (abs(i - (t*e[0][q])) % Lx);

							// Compute PDF origin ijk
							int j_o = (abs(j - (t*e[1][q])) % Ly);

							//idx = q + (i_o*Q) + (j_o*Q*Lx);
							idx = (i_o)+(j_o*Lx) + (q*Q*Lx);
							f[idx] = f_now[q] - (1.0/(tau)) * (f_now[q]-f_eq[q]) + force_term[q];
						#endif

						#ifdef TWO_LATTICE
							int i_target = (i+e[0][q]);
							int j_target = (j+e[1][q]);
							if(i_target<0) i_target = Lx-1;
							if(i_target>(Lx-1)) i_target = 0;
							if(j_target<0) j_target = Ly-1;
							if(j_target>(Ly-1)) j_target = 0;
							//int idx = q+(i_target*Q)+(j_target*Q*Lx);
							int idx = (i_target)+(j_target*Lx) + (q*Q*Lx);
							f[idx] = f_now[q] - (1.0/(tau)) * (f_now[q]-f_eq[q]) + force_term[q];
						#endif
					}
				}
			}
		}

		#ifdef OUTPUT
			if(t%output_frequency ==0) {
				writeOutput(t);
			}
		#endif
	}

	// Stop the timer, claculate and print time elapsed
	clock_t end_time = clock();
    double total_time = (double)(end_time - start_time);
    printf("Total Time: %.15f\n",total_time/CLOCKS_PER_SEC);
	getchar();
	
	return 0;
}