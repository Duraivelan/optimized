#include <iostream>
#include <random>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <tuple>
#include <dirent.h>
#include <algorithm>
#include <functional>
#include <array>
# include "defs.h"
# include "rigid_force.h"
//#include</storage3/usr/people/duraivelan/Downloads/eigen-eigen-bdd17ee3b1b3/Eigen/Eigenvalues>
#include<Eigen/Eigenvalues>

using namespace Eigen;

using namespace std;

 struct RowSort {
        bool operator()(vector<int> a, vector<int>  b)
        {   
            return a[0] < b[0];
        }   
    } ;

// random numbers using random_device option with normal distribution
/*void createInitialPosition_N_particles(std::string fileName, int N) {
      std::random_device rd, rd1;
      std::mt19937 genx(rd()),geny(rd1());
      std::ofstream outFile(fileName);
      std::normal_distribution<> d(0,1);
      for(int i=0;i<N;i++) {
      outFile<<d(genx)<<'\t'<<d(geny)<<std::endl;
      }
      outFile.close();
}
*/

// random numbers using rand function
void createInitialPosition_N_particles(std::string fileName, int N, double Lx, double Ly, double Lz) {
	double x,y,z;
 	srand (time(NULL)); // initialize random seed
 	std::ofstream outFile(fileName);
 	for(int i=0;i<N;i++) {
 		x=((double) rand() / (RAND_MAX/Lx))-Lx/2.0;  // create particle position from -Lx/2 to Lx/2
		y=((double) rand() / (RAND_MAX/Ly))-Ly/2.0;
		z=((double) rand() / (RAND_MAX/Lz))-Lz/2.0;
		outFile<<x<<'\t'<<y<<'\t'<<z<<std::endl;
 	}
 	outFile.close();
}

std::vector<int> radialDistFunc(double XYZ[][3], double Lx,double Ly, double Lz, double dr, int N) {
    std::vector<int> rdf((int) floor(sqrt(pow(Lx/2,2)+pow(Ly/2,2)+pow(Lz/2,2)))/dr,0);
	double r;
	for(int j=0;j<N;j++) {
		r=sqrt(pow(XYZ[j][0],2)+pow(XYZ[j][1],2)+pow(XYZ[j][2],2));
	    rdf[(int) floor(r/dr)]+=1;                        // put each particle in a bin according to its position from origin ie. (0,0)
	}
	return rdf;
}

// forceUpdate fucntion included as force.h header file


// mobility calculation routine

extern "C" {
    // LU decomoposition of a general matrix
    void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);

    // generate inverse of a matrix given its LU decomposition
    void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);
}

void inverse(double* A, int N)
{
    int *IPIV = new int[N+1];
    int LWORK = N*N;
    double *WORK = new double[LWORK];
    int INFO;

    dgetrf_(&N,&N,A,&N,IPIV,&INFO);
    dgetri_(&N,A,&N,IPIV,WORK,&LWORK,&INFO);

    delete IPIV;
    delete WORK;
}

void mobility_calc(int NrParticles){
                 
vctr3D dR, dr2;
double R, r2;
vector<SubData>  bead(NrParticles);

// variables for mobility tensor calculation
double eta_0=eta;
vctr3D e_ab , e_ab_unit ;
double e_ab2, e_ab2_inv, temp, temp1, temp2, temp3, tau ;
//

std::string fileName="new_cluster.dat";

//read x,y positions from XY.dat
std::ifstream dataFile(fileName);
if(!dataFile.good()) {
	std::cerr<<"Given file is corrupt /n"<<std::endl;
}
else {
    std::string line;
    for (int i=0;i<NrParticles;i++) {
		std::getline(dataFile,line);
    	std::istringstream currentLine(line);    
        currentLine >> bead[i].pos.comp[0];
        currentLine >> bead[i].pos.comp[1];
        currentLine >> bead[i].pos.comp[2];
		currentLine >> bead[i].radius;
    }
}	
	
 	
      
	mtrx3D Mobility_Tnsr_tt;
	mtrx3D Mobility_Tnsr_tr;
	mtrx3D Mobility_Tnsr_rt;
	mtrx3D Mobility_Tnsr_rr;
	mtrx35D Mobility_Tnsr_td;
	mtrx35D Mobility_Tnsr_rd;
	mtrx53D Mobility_Tnsr_dt;
	mtrx53D Mobility_Tnsr_dr;
	mtrx55D Mobility_Tnsr_dd;

	mtrx3D Resistance_Tnsr_tt;
	mtrx3D Resistance_Tnsr_tr;
	mtrx3D Resistance_Tnsr_rt;
	mtrx3D Resistance_Tnsr_rr;	
	mtrx35D Resistance_Tnsr_rd;
	mtrx35D Resistance_Tnsr_td;
	mtrx53D Resistance_Tnsr_dt;
	mtrx53D Resistance_Tnsr_dr;
	mtrx55D Resistance_Tnsr_dd;
	
// Kronecker delta

	double kron_del[3][3] = {	
								{1.0,0.0,0.0},
								{0.0,1.0,0.0},
								{0.0,0.0,1.0}
							};
							

// Levi-Civita 

	double Levi_Civi[3][3][3] = {
							{{0.0,0.0,0.0},{0.0,0.0,1.0},{0.0,-1.0,0.0}},
							{{0.0,0.0,-1.0},{0.0,0.0,0.0},{1.0,0.0,0.0}},
							{{0.0,1.0,0.0},{-1.0,0.0,0.0},{0.0,0.0,0.0}}
							};
	
// three and four index mobility matrices	
	double g_ijk[3][3][3];
	double h_ijk[3][3][3];
	double m_ijkl[3][3][3][3];

// three and four index resistance matrices	
	double G_IJK[3][3][3];
	double H_IJK[3][3][3];
	double M_IJKL[3][3][3][3];
		
// the five base matrices for strain tensor // option 5 :  equation 419 wouter's tex version clusterdyn_110816_1556

	double e[5][3][3]= {
							{{1.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,-1.0}},
							{{0.0,1.0,0.0},{1.0,0.0,0.0},{0.0,0.0,0.0}},
							{{0.0,0.0,1.0},{0.0,0.0,0.0},{1.0,0.0,0.0}},
							{{0.0,0.0,0.0},{0.0,0.0,1.0},{0.0,1.0,0.0}},
							{{0.0,0.0,0.0},{0.0,1.0,0.0},{0.0,0.0,-1.0}}
						};


	double e_l[5][3][3]= {
							{{ 2.0/3.0,0.0,0.0},{0.0,-1.0/3.0,0.0},{0.0,0.0,-1.0/3.0}},
							{{0.0,0.5,0.0},{0.5,0.0,0.0},{0.0,0.0,0.0}},
							{{0.0,0.0,0.5},{0.0,0.0,0.0},{0.5,0.0,0.0}},
							{{0.0,0.0,0.0},{0.0,0.0,0.5},{0.0,0.5,0.0}},
							{{-1.0/3.0,0.0,0.0},{0.0, 2.0/3.0,0.0},{0.0,0.0,-1.0/3.0}}
						};
   
   double mu_11N[121*NrParticles*NrParticles] ;  		// grand mobility matrix
   double zeta_11N[121*NrParticles*NrParticles] ;  	// grand resistance matrix
   double rho_11N[121*NrParticles*NrParticles] ;  	// grand resistance matrix
   double xi_11x11[11*11] ; 							// generalized friction matrix
   
   for (int i=0; i<121; i++)
		{
			xi_11x11[i] = 0.0; 
		}        
		
 std::ofstream outFile1("data.dat");

// important all lengths have been normalized by particle radius as metioned in Page 46, Appendix A - Durlofsky, Louis, John F. Brady, and Georges Bossis. 
				// "Dynamic simulation of hydrodynamically interacting particles." Journal of fluid mechanics 180 (1987): 21-49.
				// for ease of programming. 

				double	a_norm = 1.0; // /(6.0*M_PI*eta_0*bead[0].radius);						// mobility matrix a non-dimensionalized by 6*pi*mu*r
				double	b_norm = 1.0; // /(6.0*M_PI*eta_0*bead[0].radius*bead[0].radius);						// mobility matrix a non-dimensionalized by 6*pi*mu*r2
				double	c_norm = 1.0; // /(6.0*M_PI*eta_0*bead[0].radius*bead[0].radius*bead[0].radius);						// mobility matrix a non-dimensionalized by 6*pi*mu*r3
				double	g_norm = 1.0; // /(6.0*M_PI*eta_0*bead[0].radius*bead[0].radius*bead[0].radius);						//		assuming correction factor of 6*pi*mu*r3	
				double	h_norm = 1.0; // /(6.0*M_PI*eta_0*bead[0].radius*bead[0].radius*bead[0].radius);						//		assuming correction factor of 6*pi*mu*r3				
				double	m_norm = 1.0; // /(6.0*M_PI*eta_0*bead[0].radius*bead[0].radius*bead[0].radius);						//		assuming correction factor of 6*pi*mu*r3	
		
for (int a=0; a<NrParticles; a++)
	{
		for (int b=a; b<NrParticles; b++)
			{
				e_ab=bead[a].pos-bead[b].pos;
				mtrx3D Pab(e_ab, e_ab) ;
				e_ab2=e_ab.norm2();
				vctr3D col1(0.0, -e_ab.comp[2], e_ab.comp[1]);
				vctr3D col2(e_ab.comp[2],0.0,-e_ab.comp[0]);
				vctr3D col3(-e_ab.comp[1],e_ab.comp[0],0.0);
				mtrx3D epsilon_e_ab(col1 , col2 , col3);
				e_ab2_inv=1.0/e_ab2;
				
				e_ab_unit = e_ab*sqrt(e_ab2_inv); 
				if(a==b) {	e_ab_unit = e_ab*0.0; }	
				tau = 1.0/(6.0*M_PI*eta_0*bead[a].radius);
			    double r 	= sqrt(e_ab2)/bead[a].radius;			// distance between particle vector 'r' magnitude |r| normalized by particle radius 'a' ;
			    double r_1 	= 1.0/(r);
			    double r_2 	= 1.0/(r*r);			    
			    double r_3 	= 1.0/(r*r*r);
			    double r_4 	= 1.0/(r*r*r*r);
			    double r_5 	= 1.0/(r*r*r*r*r); 

				// mobility scalar values - as defined in Page 46, Appendix A - Durlofsky, Louis, John F. Brady, and Georges Bossis. 
				// "Dynamic simulation of hydrodynamically interacting particles." Journal of fluid mechanics 180 (1987): 21-49.

				double x_a[2][2] = {{	1.0		,	3.0*r_1/2.0		-	1.0*r_3			},{	3.0*r_1/2.0		-	1.0*r_3			,	1.0		} }; 
				double y_a[2][2] = {{	1.0		,	(3.0*r_1/4.0)	+	(1.0*r_3/2.0)	},{	3.0*r_1/4.0		+	1.0*r_3/2.0		,	1.0		} }; 
				double y_b[2][2] = {{	0.0		,  -3.0*r_2/4.0							},{	3.0*r_2/4.0							,	0.0		} }; 
				double x_c[2][2] = {{	3.0/4.0	,  	3.0*r_3/4.0							},{	3.0*r_3/4.0							,  	3.0/4.0	} }; 
				double y_c[2][2] = {{	3.0/4.0	,  -3.0*r_3/8.0							},{-3.0*r_3/8.0							,  	3.0/4.0	} }; 

				double x_g[2][2] = {{	0.0		,	9.0*r_2/4.0		-	18.0*r_4/5.0	},{-9.0*r_2/4.0		+	18.0*r_4/5.0	,	0.0		} };
				double y_g[2][2] = {{	0.0		,	6.0*r_4/5.0							},{-6.0*r_4/5.0							,	0.0		} };
				double y_h[2][2] = {{	0.0		,  -9.0*r_3/8.0							},{-9.0*r_3/8.0							,	0.0		} };
				double x_m[2][2] = {{	9.0/10.0,  -9.0*r_3/2.0		+ 	54.0*r_5/5.0	},{-9.0*r_3/2.0		+ 	54.0*r_5/5.0	,	9.0/10.0} };
				double y_m[2][2] = {{	9.0/10.0,   9.0*r_3/4.0		- 	36.0*r_5/5.0	},{ 9.0*r_3/4.0		- 	36.0*r_5/5.0	,	9.0/10.0} };
				double z_m[2][2] = {{	9.0/10.0,  					 	 9.0*r_5/5.0	},{ 				 	 9.0*r_5/5.0	,	9.0/10.0} };					

				if(a==b) {

				Mobility_Tnsr_tr	= 		null33D ;
					
				for (int i=0; i<3; i++)
					{
					for (int j=0; j<3; j++)
						{
							double ep_ijk_e_k = 0.0;
							
						for (int k=0; k<3; k++)
							{	

								double ep_jkl_e_l = 0.0;
								double ep_ikl_e_l = 0.0;
								
								for (int l=0; l<3; l++)
									{
										ep_jkl_e_l	+=	Levi_Civi[j][k][l]*e_ab_unit.comp[l];
										ep_ikl_e_l	+=	Levi_Civi[i][k][l]*e_ab_unit.comp[l];
										
										m_ijkl[i][j][k][l]	=	m_norm*((3.0/2.0)*x_m[1][1]*(e_ab_unit.comp[i]*e_ab_unit.comp[j] 					-	(1.0/3.0)*kron_del[i][j])*(e_ab_unit.comp[k]*e_ab_unit.comp[l]	
																-(1.0/3.0)*kron_del[k][l])
																+(1.0/2.0)*y_m[1][1]*(e_ab_unit.comp[i]*kron_del[j][l]*e_ab_unit.comp[k]	+	e_ab_unit.comp[j]*kron_del[i][l]*e_ab_unit.comp[k]
																					+ e_ab_unit.comp[i]*kron_del[j][k]*e_ab_unit.comp[l]	+ 	e_ab_unit.comp[j]*kron_del[i][k]*e_ab_unit.comp[l]
																					- 4.0*e_ab_unit.comp[i]*e_ab_unit.comp[j]*e_ab_unit.comp[k]*e_ab_unit.comp[l]	)
																			
																+(1.0/2.0)*z_m[1][1]*(kron_del[i][k]*kron_del[j][l]		+ 	kron_del[j][k]*kron_del[i][l]	- 	kron_del[i][j]*kron_del[k][l] 
																+ e_ab_unit.comp[i]*e_ab_unit.comp[j]*kron_del[k][l]	+	kron_del[i][j]*e_ab_unit.comp[k]*e_ab_unit.comp[l]	
																+ e_ab_unit.comp[i]*e_ab_unit.comp[j]*e_ab_unit.comp[k]*e_ab_unit.comp[l]
																- e_ab_unit.comp[i]*kron_del[j][l]*e_ab_unit.comp[k]	- 	e_ab_unit.comp[j]*kron_del[i][l]*e_ab_unit.comp[k]
																- e_ab_unit.comp[i]*kron_del[j][k]*e_ab_unit.comp[l]	- 	e_ab_unit.comp[j]*kron_del[i][k]*e_ab_unit.comp[l]
																));
																																							
									}	// l
									
								ep_ijk_e_k					+=	Levi_Civi[i][j][k]*e_ab_unit.comp[k];
								
								g_ijk[i][j][k]				=	g_norm*(x_g[1][1]*(e_ab_unit.comp[i]*e_ab_unit.comp[j] 	-	(1.0/3.0)*kron_del[i][j])*e_ab_unit.comp[k]
																+ 		y_g[1][1]*(e_ab_unit.comp[i]*kron_del[j][k]		+ 	e_ab_unit.comp[j]*kron_del[i][k]	-	2.0*e_ab_unit.comp[i]*e_ab_unit.comp[j]*e_ab_unit.comp[k]	)	);
										
								h_ijk[i][j][k]				= 	h_norm*(y_h[1][1]*(e_ab_unit.comp[i]*ep_jkl_e_l			+	e_ab_unit.comp[j]*ep_ikl_e_l										)	);
							}	// k		

							Mobility_Tnsr_tt.comp[i][j]		=	a_norm*(x_a[1][1]*e_ab_unit.comp[i]*e_ab_unit.comp[j]	+ 	y_a[1][1]*(kron_del[i][j]	- e_ab_unit.comp[i]*e_ab_unit.comp[j]	)	);
							
							Mobility_Tnsr_rt.comp[i][j]		=	b_norm*(													y_b[1][1]*ep_ijk_e_k													);
						
							Mobility_Tnsr_rr.comp[i][j]		=	c_norm*(x_c[1][1]*e_ab_unit.comp[i]*e_ab_unit.comp[j]	+ 	y_c[1][1]*(kron_del[i][j]	- e_ab_unit.comp[i]*e_ab_unit.comp[j]	)	);

						//		cout << "m_ijkl"	<< endl;				
						//		cout << m_ijkl[0][0][0][0] << endl;				
						}	// j
					}	// i
		/*
					for (int i=0; i<3; i++)
					{
					for (int j=0; j<3; j++)
						{									
							cout << "matrix"	<< endl;				
							cout << Mobility_Tnsr_rt.comp[i][j] << endl;
						}
					}														
			
				Mobility_Tnsr_tt	=	 	Unit_diag * tau ;
											
				Mobility_Tnsr_rr	=		Unit_diag * temp2 ;
				Mobility_Tnsr_rt	= 		null33D ;
				Mobility_Tnsr_tr	= 		null33D ;
		*/

			 } else {
				 					
				for (int i=0; i<3; i++)
					{
					for (int j=0; j<3; j++)
						{
							double ep_ijk_e_k = 0.0;
							
						for (int k=0; k<3; k++)
							{	

								double ep_jkl_e_l = 0.0;
								double ep_ikl_e_l = 0.0;
								
								for (int l=0; l<3; l++)
									{
										ep_jkl_e_l	+=	Levi_Civi[j][k][l]*e_ab_unit.comp[l];
										ep_ikl_e_l	+=	Levi_Civi[i][k][l]*e_ab_unit.comp[l];
										
										m_ijkl[i][j][k][l]	=	 m_norm*((3.0/2.0)*x_m[1][0]*(e_ab_unit.comp[i]*e_ab_unit.comp[j] 			-	(1.0/3.0)*kron_del[i][j])*(e_ab_unit.comp[k]*e_ab_unit.comp[l]	
																-(1.0/3.0)*kron_del[k][l])
																+(1.0/2.0)*y_m[1][0]*(e_ab_unit.comp[i]*kron_del[j][l]*e_ab_unit.comp[k]	+	e_ab_unit.comp[j]*kron_del[i][l]*e_ab_unit.comp[k]
																					+ e_ab_unit.comp[i]*kron_del[j][k]*e_ab_unit.comp[l]	+ 	e_ab_unit.comp[j]*kron_del[i][k]*e_ab_unit.comp[l]
																					- 4.0*e_ab_unit.comp[i]*e_ab_unit.comp[j]*e_ab_unit.comp[k]*e_ab_unit.comp[l]	)
																			
																+(1.0/2.0)*z_m[1][0]*(kron_del[i][k]*kron_del[j][l]		+ 	kron_del[j][k]*kron_del[i][l]	- 	kron_del[i][j]*kron_del[k][l]
																+ e_ab_unit.comp[i]*e_ab_unit.comp[j]*kron_del[k][l]	+	kron_del[i][j]*e_ab_unit.comp[k]*e_ab_unit.comp[l]	
																+ e_ab_unit.comp[i]*e_ab_unit.comp[j]*e_ab_unit.comp[k]*e_ab_unit.comp[l]
																- e_ab_unit.comp[i]*kron_del[j][l]*e_ab_unit.comp[k]	- 	e_ab_unit.comp[j]*kron_del[i][l]*e_ab_unit.comp[k]
																- e_ab_unit.comp[i]*kron_del[j][k]*e_ab_unit.comp[l]	- 	e_ab_unit.comp[j]*kron_del[i][k]*e_ab_unit.comp[l]
																));
																																							
									}	// l
									
								ep_ijk_e_k					+=	Levi_Civi[i][j][k]*e_ab_unit.comp[k];
								
								g_ijk[i][j][k]				=	g_norm*(x_g[1][0]*(e_ab_unit.comp[i]*e_ab_unit.comp[j] 	-	(1.0/3.0)*kron_del[i][j])*e_ab_unit.comp[k]
																+ 		y_g[1][0]*(e_ab_unit.comp[i]*kron_del[j][k]		+ 	e_ab_unit.comp[j]*kron_del[i][k]	-	2.0*e_ab_unit.comp[i]*e_ab_unit.comp[j]*e_ab_unit.comp[k]	)	);
										
								h_ijk[i][j][k]				= 	h_norm*(y_h[1][0]*(e_ab_unit.comp[i]*ep_jkl_e_l			+	e_ab_unit.comp[j]*ep_ikl_e_l										)	);
							}	// k		

							Mobility_Tnsr_tt.comp[i][j]		=	a_norm*(x_a[1][0]*e_ab_unit.comp[i]*e_ab_unit.comp[j]	+ 	y_a[1][0]*(kron_del[i][j]	- e_ab_unit.comp[i]*e_ab_unit.comp[j]	)	);
							
							Mobility_Tnsr_rt.comp[i][j]		=	b_norm*(													y_b[1][0]*ep_ijk_e_k													);
						
							Mobility_Tnsr_rr.comp[i][j]		=	c_norm*(x_c[1][0]*e_ab_unit.comp[i]*e_ab_unit.comp[j]	+ 	y_c[1][0]*(kron_del[i][j]	- e_ab_unit.comp[i]*e_ab_unit.comp[j]	)	);
									
							Mobility_Tnsr_tr	= 	    Mobility_Tnsr_rt*(1.0);		
		//						cout << "expresssion"	<< endl;				
		//						cout << Mobility_Tnsr_rt.comp[i][j] << endl;		
						}	// j
					}	// i	
				}				

//	extract the reduced index mobility tesnors from g_ijk, h_ijk and m_ijkl

// based on the equation mu^dt_{\p\g} = (e^p)_{\a\b}*mu^dt_{\a\b\g}	, etc in wouter notes equation no : 422 , version : 110816_1556
	
	Mobility_Tnsr_dt	= 		null53D ;	
	Mobility_Tnsr_dr	= 		null53D ;	
	Mobility_Tnsr_dd	= 		null55D ;
		
	Resistance_Tnsr_dt	= 		null53D ;	
	Resistance_Tnsr_dr	= 		null53D ;	
	Resistance_Tnsr_dd	= 		null55D ;
	
	Mobility_Tnsr_td	= 		null35D ;	
	Mobility_Tnsr_rd	= 		null35D ;	
		
	Resistance_Tnsr_td	= 		null35D ;	
	Resistance_Tnsr_rd	= 		null35D ;	
	
	for (int p=0; p<5; p++)
		{
		for (int g=0; g<3; g++)
			{
			for (int a=0; a<3; a++)
				{
				for (int b=0; b<3; b++)
					{		
						Mobility_Tnsr_dt.comp[p][g]		+=		e[p][a][b]*g_ijk[a][b][g];	
						Mobility_Tnsr_dr.comp[p][g]		+=		e[p][a][b]*h_ijk[a][b][g];		
						Mobility_Tnsr_td.comp[g][p]		+=		e[p][a][b]*g_ijk[a][b][g];	// going from g_dt matrix to ~g_td matrix hence circulation of indices 
						Mobility_Tnsr_rd.comp[g][p]		+=		e[p][a][b]*h_ijk[a][b][g];
					}
				}				
			}
		for (int s=0; s<5; s++)
			{
			for (int a=0; a<3; a++)
				{
				for (int b=0; b<3; b++)
					{
					for (int g=0; g<3; g++)
						{
						for (int d=0; d<3; d++)
							{							
								Mobility_Tnsr_dd.comp[p][s]		+=		e[p][a][b]*m_ijkl[a][b][g][d]*e[s][g][d];			
							}
						}													
					}
				}
			}
		}		
				
			for (int l=0; l<3; l++)
				{
					for (int k=0; k<3; k++)
						{		
							// 11N column major format
							zeta_11N[k	+	11*NrParticles*l	+	3*a	+	33*NrParticles*b														] 	=	 Mobility_Tnsr_tt.comp[k][l];
							zeta_11N[k	+	11*NrParticles*l	+	3*a	+	33*NrParticles*b	+	33*NrParticles*NrParticles						] 	=	 Mobility_Tnsr_tr.comp[k][l];
							zeta_11N[k	+	11*NrParticles*l	+	3*a	+	33*NrParticles*b	+	3*NrParticles									] 	=	 Mobility_Tnsr_rt.comp[k][l];
							zeta_11N[k	+	11*NrParticles*l	+	3*a	+	33*NrParticles*b	+	33*NrParticles*NrParticles	+	3*NrParticles	] 	=	 Mobility_Tnsr_rr.comp[k][l];							
							zeta_11N[k	+	11*NrParticles*l	+	3*b	+	33*NrParticles*a														] 	=	 Mobility_Tnsr_tt.comp[k][l];
							zeta_11N[k	+	11*NrParticles*l	+	3*b	+	33*NrParticles*a	+	33*NrParticles*NrParticles						] 	=	-Mobility_Tnsr_rt.comp[k][l];
							zeta_11N[k	+	11*NrParticles*l	+	3*b	+	33*NrParticles*a	+	3*NrParticles									] 	=	-Mobility_Tnsr_rt.comp[k][l];
							zeta_11N[k	+	11*NrParticles*l	+	3*b	+	33*NrParticles*a	+	33*NrParticles*NrParticles	+	3*NrParticles	] 	=	 Mobility_Tnsr_rr.comp[k][l];
/*							// 6N column major format
							zeta_11N[k	+	6*NrParticles*l	+	3*a	+	18*NrParticles*b														] 	=	 Mobility_Tnsr_tt.comp[k][l];
							zeta_11N[k	+	6*NrParticles*l	+	3*a	+	18*NrParticles*b	+	18*NrParticles*NrParticles						] 	=	 Mobility_Tnsr_tr.comp[k][l];
							zeta_11N[k	+	6*NrParticles*l	+	3*a	+	18*NrParticles*b	+	3*NrParticles									] 	=	 Mobility_Tnsr_rt.comp[k][l];
							zeta_11N[k	+	6*NrParticles*l	+	3*a	+	18*NrParticles*b	+	18*NrParticles*NrParticles	+	3*NrParticles	] 	=	 Mobility_Tnsr_rr.comp[k][l];							
							zeta_11N[k	+	6*NrParticles*l	+	3*b	+	18*NrParticles*a														] 	=	 Mobility_Tnsr_tt.comp[k][l];
							zeta_11N[k	+	6*NrParticles*l	+	3*b	+	18*NrParticles*a	+	18*NrParticles*NrParticles						] 	=	-Mobility_Tnsr_rt.comp[k][l];
							zeta_11N[k	+	6*NrParticles*l	+	3*b	+	18*NrParticles*a	+	3*NrParticles									] 	=	-Mobility_Tnsr_rt.comp[k][l];
							zeta_11N[k	+	6*NrParticles*l	+	3*b	+	18*NrParticles*a	+	18*NrParticles*NrParticles	+	3*NrParticles	] 	=	 Mobility_Tnsr_rr.comp[k][l];
*/						}
				}
				
			for (int l=0; l<5; l++)
				{
				for (int k=0; k<3; k++)
					{				
						// column major format
						zeta_11N[k	+	11*NrParticles*l	+	3*a	+	55*NrParticles*b	+	66*NrParticles*NrParticles						] 	=	 -Mobility_Tnsr_td.comp[k][l];
						zeta_11N[k	+	11*NrParticles*l	+	3*a	+	55*NrParticles*b	+	66*NrParticles*NrParticles	+	3*NrParticles	] 	=	 Mobility_Tnsr_rd.comp[k][l];
						zeta_11N[k	+	11*NrParticles*l	+	3*b	+	55*NrParticles*a	+	66*NrParticles*NrParticles						] 	=	 Mobility_Tnsr_td.comp[k][l];
						zeta_11N[k	+	11*NrParticles*l	+	3*b	+	55*NrParticles*a	+	66*NrParticles*NrParticles	+	3*NrParticles	] 	=	 Mobility_Tnsr_rd.comp[k][l];
						
					}
				}					
			for (int l=0; l<3; l++)
				{
				for (int k=0; k<5; k++)
					{				
						// column major format
						zeta_11N[k	+	11*NrParticles*l	+	5*a	+	33*NrParticles*b	+	6*NrParticles									] 	=	 Mobility_Tnsr_td.comp[l][k];
						zeta_11N[k	+	11*NrParticles*l	+	5*a	+	33*NrParticles*b	+	33*NrParticles*NrParticles	+	6*NrParticles	] 	=	 Mobility_Tnsr_rd.comp[l][k];
						zeta_11N[k	+	11*NrParticles*l	+	5*b	+	33*NrParticles*a	+	6*NrParticles									] 	=	 -Mobility_Tnsr_td.comp[l][k];
						zeta_11N[k	+	11*NrParticles*l	+	5*b	+	33*NrParticles*a	+	33*NrParticles*NrParticles	+	6*NrParticles	] 	=	 Mobility_Tnsr_rd.comp[l][k];						
					}
				}
			for (int l=0; l<5; l++)
				{
				for (int k=0; k<5; k++)
					{				
						// column major format
						zeta_11N[k	+	11*NrParticles*l	+	5*a	+	55*NrParticles*b	+	66*NrParticles*NrParticles	+	6*NrParticles	] 	=	 Mobility_Tnsr_dd.comp[k][l];
						zeta_11N[k	+	11*NrParticles*l	+	5*b	+	55*NrParticles*a	+	66*NrParticles*NrParticles	+	6*NrParticles	] 	=	 Mobility_Tnsr_dd.comp[k][l];
					}
				}

				
			}	

		}

	mtrx3D Friction_Tnsr_tt(0.0,0.0,0.0);
	mtrx3D Friction_Tnsr_tr(0.0,0.0,0.0);
	mtrx3D Friction_Tnsr_rt(0.0,0.0,0.0);
	mtrx3D Friction_Tnsr_rr(0.0,0.0,0.0);
 	mtrx53D Friction_Tnsr_dt 	=	null53D;
	mtrx53D Friction_Tnsr_dr	= 	null53D;
	mtrx35D Friction_Tnsr_td	=	null35D;
	mtrx35D Friction_Tnsr_rd	=	null35D;
	mtrx55D Friction_Tnsr_dd	=	null55D;
	
 for (int i=0;i<121;i++)
       { 
               for(int j=0;j<121;j++)
                       {
               cout<< std::setprecision(5) <<zeta_11N[i+11*NrParticles*j]<<'\t' ;
                       }
               cout<<std::endl; 
       }	
               cout<<std::endl; 
               
	inverse ( zeta_11N ,11*NrParticles )	 ; 	
	
	 for (int i=0;i<121;i++)
       { 
               for(int j=0;j<121;j++)
                       {
               cout<< std::setprecision(5) <<zeta_11N[i+11*NrParticles*j]<<'\t' ;
                       }
               cout<<std::endl; 
       }	
               cout<<std::endl; 
	
	for (int i=0; i<NrParticles; i++)
		{
			vctr3D col1 (  0.0 						,	bead[i].pos.comp[2] , -bead[i].pos.comp[1]	);
			vctr3D col2 ( -bead[i].pos.comp[2] 	,   0.0						,  bead[i].pos.comp[0]	);
			vctr3D col3 (  bead[i].pos.comp[1] 	,  -bead[i].pos.comp[0] , 0.0						);
			mtrx3D	Ai(	col1 , col2 , col3 ) ;   
					
			for (int j=0; j<NrParticles; j++)
				{
					
					vctr3D col4 (  0.0 						,	bead[j].pos.comp[2] , -bead[j].pos.comp[1]	);
					vctr3D col5 ( -bead[j].pos.comp[2] 	,   0.0						,  bead[j].pos.comp[0]	);
					vctr3D col6 (  bead[j].pos.comp[1] 	,  -bead[j].pos.comp[0] , 0.0						);
					mtrx3D	Aj(	col4 , col5 , col6 ) ; 

// if	(del_j)_{\alpha,p}	*	E^{inf\tilde}_{p}	= 	E^inf_{\alpha\beta} *	r_\beta 

// then	(del_j)_{\alpha,p}	*	E^{inf\tilde}_{p}	= 	(e_p)_{\alpha\beta}	*	E^{inf\tilde}_{p}	*	r_\beta	// sicne E^{inf}_{\alpha\beta}	= (e_p)_{\alpha\beta}*E^{inf\tilde}_{\p}

//hence	(del_j)_{\alpha,p}							= 	(e_p)_{\alpha\beta}	*	r_\beta	

					mtrx35D	Delj;
					mtrx53D	Deli;
					mtrx53D	Unit_tnsr_Redc;
					
						for (int k=0; k<5; k++)
						{							
							for (int a=0; a<3; a++)
							{
								Delj.comp[a][k]	=	0.0	;
								Deli.comp[k][a]	=	0.0	;
								for (int b=0; b<3; b++)
								{
									Delj.comp[a][k]	+=	e_l[k][a][b]	*	bead[j].pos.comp[b]	;
									Deli.comp[k][a]	+=	e_l[k][a][b]	*	bead[i].pos.comp[b]	;
									for (int c=0; c<3; c++)
									{
										Unit_tnsr_Redc.comp[k][a]	+=	e_l[k][b][c]		* (b==c)	*	bead[i].pos.comp[a]	;	
									}
								}
							}		
						}
	
					
					for (int l=0; l<3; l++)
						{
							for (int k=0; k<3; k++)
								{
							

// 									11N format 
									Resistance_Tnsr_tt.comp[k][l] =	zeta_11N[k	+	11*NrParticles*l	+	3*i	+	33*NrParticles*j														];
									Resistance_Tnsr_tr.comp[k][l] =	zeta_11N[k	+	11*NrParticles*l	+	3*i	+	33*NrParticles*j	+	33*NrParticles*NrParticles						];
									Resistance_Tnsr_rt.comp[k][l] =	zeta_11N[k	+	11*NrParticles*l	+	3*i	+	33*NrParticles*j	+	3*NrParticles									];
									Resistance_Tnsr_rr.comp[k][l] = zeta_11N[k	+	11*NrParticles*l	+	3*i	+	33*NrParticles*j	+	33*NrParticles*NrParticles	+	3*NrParticles	];


/*
// 									6N format
									Resistance_Tnsr_tt.comp[k][l] = zeta_11N[k+6*NrParticles*l+3*i+18*NrParticles*j] ;
									Resistance_Tnsr_tr.comp[k][l] = zeta_11N[k+6*NrParticles*l+3*i+18*NrParticles*j+18*NrParticles*NrParticles]	;
									Resistance_Tnsr_rt.comp[k][l] = zeta_11N[k+6*NrParticles*l+3*i+18*NrParticles*j+3*NrParticles] ;
									Resistance_Tnsr_rr.comp[k][l] = zeta_11N[k+6*NrParticles*l+3*i+18*NrParticles*j+18*NrParticles*NrParticles+3*NrParticles] ;
*/								
							}
						}


					for (int l=0; l<5; l++)
						{
							for (int k=0; k<3; k++)
								{				
									// column major format
									Resistance_Tnsr_td.comp[k][l] =	zeta_11N[k	+	11*NrParticles*l	+	3*i	+	55*NrParticles*j	+	66*NrParticles*NrParticles						];
									Resistance_Tnsr_rd.comp[k][l] =	zeta_11N[k	+	11*NrParticles*l	+	3*i	+	55*NrParticles*j	+	66*NrParticles*NrParticles	+	3*NrParticles	];						
								}
						}					
					
					for (int l=0; l<3; l++)
						{
							for (int k=0; k<5; k++)
								{				
									// column major format
									Resistance_Tnsr_dt.comp[k][l] =	zeta_11N[k	+	11*NrParticles*l	+	5*i	+	33*NrParticles*j	+	6*NrParticles									];
									Resistance_Tnsr_dr.comp[k][l] =	zeta_11N[k	+	11*NrParticles*l	+	5*i	+	33*NrParticles*j	+	33*NrParticles*NrParticles	+	6*NrParticles	];				
								}
						}
					
					for (int l=0; l<5; l++)
						{
							for (int k=0; k<5; k++)
								{				
									// column major format
									Resistance_Tnsr_dd.comp[k][l] =	zeta_11N[k	+	11*NrParticles*l	+	5*i	+	55*NrParticles*j	+	66*NrParticles*NrParticles	+	6*NrParticles	];
								}
						}						
					
					Friction_Tnsr_tt	+=		Resistance_Tnsr_tt ;  
					Friction_Tnsr_tr 	+= 	( 	Resistance_Tnsr_tr		- 	(	Resistance_Tnsr_tt*Aj	)	)	;
				//	Friction_Tnsr_rt 	+= 	(	Ai*Resistance_Tnsr_tt	+		Resistance_Tnsr_rt    	)    	; 
					Friction_Tnsr_rr 	+= 	( 	Resistance_Tnsr_rr    	-	(	Resistance_Tnsr_rt*Aj 	) 	+ 			Ai*Resistance_Tnsr_tr	-	Ai*Resistance_Tnsr_tt*Aj	)	;
				
					for (int l=0; l<5; l++)
						{
						for (int k=0; k<3; k++)
							{
								for (int m=0; m<3; m++)
									{
										Friction_Tnsr_td.comp[k][l]	+=	Resistance_Tnsr_tt.comp[k][m]*Delj.comp[m][l];
										Friction_Tnsr_rd.comp[k][l]	+= 	Resistance_Tnsr_rt.comp[k][m]*Delj.comp[m][l];
										//		Friction_Tnsr_dr.comp[l][k]	-= 	Resistance_Tnsr_dt.comp[l][m]*Aj.comp[m][k];
										for (int n=0; n<3; n++)
											{
												Friction_Tnsr_rd.comp[k][l]	+=	Ai.comp[k][n]	*	( 	Resistance_Tnsr_tt.comp[n][m]	*	Delj.comp[m][l]	)	;
											}
										Friction_Tnsr_rd.comp[k][l]	+=	Ai.comp[k][m]	*	( Resistance_Tnsr_td.comp[m][l]	)	;									
									}
								Friction_Tnsr_td.comp[k][l]	+=	Resistance_Tnsr_td.comp[k][l]	;
								Friction_Tnsr_rd.comp[k][l]	+=	Resistance_Tnsr_rd.comp[k][l]	;								
								Friction_Tnsr_dt.comp[k][l]	+=	Resistance_Tnsr_dt.comp[l][k]	;								
								Friction_Tnsr_dr.comp[l][k]	+=	Resistance_Tnsr_dr.comp[l][k]	;								
							}
							
	// assuming Frc_td = Frc_dt;
	
						for (int k=0; k<3; k++)
							{
								Friction_Tnsr_dt.comp[l][k]	=	Friction_Tnsr_td.comp[k][l]	;								
								Friction_Tnsr_dr.comp[l][k]	=	Friction_Tnsr_rd.comp[k][l]	;	
							}
							

// droping the extra terms like transpose and identity terms
						for (int k=0; k<5; k++)
						{								
							for (int m=0; m<3; m++)
							{
								Friction_Tnsr_dd.comp[l][k]	+= 	Resistance_Tnsr_dt.comp[l][m]*Delj.comp[m][k];
								Friction_Tnsr_dd.comp[l][k]	+= 	Deli.comp[l][m]*Resistance_Tnsr_td.comp[m][k];
									for (int n=0; n<3; n++)
									{
										Friction_Tnsr_dd.comp[l][k]	+=	Deli.comp[l][m]*Resistance_Tnsr_tt.comp[m][n]*Delj.comp[n][k] ;		
									} 				
							}
						
							Friction_Tnsr_dd.comp[l][k]	+=	Resistance_Tnsr_dd.comp[l][k]	;		
						}						

						}
			 
				}
		}
					Friction_Tnsr_rt = ~Friction_Tnsr_tr;

/*
			for (int l=0; l<3; l++)
				{
					for (int k=0; k<3; k++)
						{				
							// column major format
							xi_11x11[k	+	11*l					] 	=	 Friction_Tnsr_tt.comp[k][l];
							xi_11x11[k	+	11*l	+	33			] 	=	 Friction_Tnsr_tr.comp[k][l];
							xi_11x11[k	+	11*l	+	3			] 	=	 Friction_Tnsr_rt.comp[k][l];
							xi_11x11[k	+	11*l	+	33	+	3	] 	=	 Friction_Tnsr_rr.comp[k][l];							
						}
				}
				
			for (int l=0; l<5; l++)
				{
				for (int k=0; k<3; k++)
					{				
						// column major format
						xi_11x11[k	+	11*l	+	66			] 	=	 Friction_Tnsr_td.comp[k][l];
						xi_11x11[k	+	11*l	+	66	+	3	] 	=	 Friction_Tnsr_rd.comp[k][l];						
					}
				}					
			for (int l=0; l<3; l++)
				{
				for (int k=0; k<5; k++)
					{				
						// column major format
						xi_11x11[k	+	11*l	+	6			] 	=	 Friction_Tnsr_td.comp[l][k];
						xi_11x11[k	+	11*l	+	33	+	6	] 	=	 Friction_Tnsr_rd.comp[l][k];					
					}
				}
			for (int l=0; l<5; l++)
				{
				for (int k=0; k<5; k++)
					{				
						// column major format
						xi_11x11[k	+	11*l	+	66	+	6	] 	=	 Friction_Tnsr_dd.comp[k][l];
					}
				}
*/
	
	 			// 6x6 format					
	 				// column major format

					xi_11x11[0] = Friction_Tnsr_tt.comp[0][0] ;  
					xi_11x11[1] = Friction_Tnsr_tt.comp[0][1] ;  
					xi_11x11[2] = Friction_Tnsr_tt.comp[0][2] ; 
					xi_11x11[6] = Friction_Tnsr_tt.comp[1][0] ; 
					xi_11x11[7] = Friction_Tnsr_tt.comp[1][1] ;  
					xi_11x11[8] = Friction_Tnsr_tt.comp[1][2] ;  
					xi_11x11[12] = Friction_Tnsr_tt.comp[2][0] ;   
					xi_11x11[13] = Friction_Tnsr_tt.comp[2][1] ; 
					xi_11x11[14] = Friction_Tnsr_tt.comp[2][2] ; 				

					xi_11x11[18] = Friction_Tnsr_rt.comp[0][0] ;  
					xi_11x11[19] = Friction_Tnsr_rt.comp[0][1] ;  
					xi_11x11[20] = Friction_Tnsr_rt.comp[0][2] ; 
					xi_11x11[24] = Friction_Tnsr_rt.comp[1][0] ; 
					xi_11x11[25] = Friction_Tnsr_rt.comp[1][1] ;  
					xi_11x11[26] = Friction_Tnsr_rt.comp[1][2] ;  
					xi_11x11[30] = Friction_Tnsr_rt.comp[2][0] ;   
					xi_11x11[31] = Friction_Tnsr_rt.comp[2][1] ; 
					xi_11x11[32] = Friction_Tnsr_rt.comp[2][2] ; 				
										
					xi_11x11[3] = Friction_Tnsr_tr.comp[0][0] ;  
					xi_11x11[4] = Friction_Tnsr_tr.comp[0][1] ;  
					xi_11x11[5] = Friction_Tnsr_tr.comp[0][2] ; 
					xi_11x11[9] = Friction_Tnsr_tr.comp[1][0] ; 
					xi_11x11[10] = Friction_Tnsr_tr.comp[1][1] ;  
					xi_11x11[11] = Friction_Tnsr_tr.comp[1][2] ;  
					xi_11x11[15] = Friction_Tnsr_tr.comp[2][0] ;   
					xi_11x11[16] = Friction_Tnsr_tr.comp[2][1] ; 
					xi_11x11[17] = Friction_Tnsr_tr.comp[2][2] ; 					
										
					xi_11x11[21] = Friction_Tnsr_rr.comp[0][0] ;  
					xi_11x11[22] = Friction_Tnsr_rr.comp[0][1] ;  
					xi_11x11[23] = Friction_Tnsr_rr.comp[0][2] ; 
					xi_11x11[27] = Friction_Tnsr_rr.comp[1][0] ; 
					xi_11x11[28] = Friction_Tnsr_rr.comp[1][1] ;  
					xi_11x11[29] = Friction_Tnsr_rr.comp[1][2] ;  
					xi_11x11[33] = Friction_Tnsr_rr.comp[2][0] ;   
					xi_11x11[34] = Friction_Tnsr_rr.comp[2][1] ; 
					xi_11x11[35] = Friction_Tnsr_rr.comp[2][2] ; 				
	
		double h_clst_ijk[3][3][3] = {{{0}}};

	for (int a=0; a<3; a++)
		{
		for (int b=0; b<3; b++)
			{
			for (int g=0; g<3; g++)
				{
					h_clst_ijk[g][a][b] = 0.0;
					
				for (int p=0; p<5; p++)
					{
							h_clst_ijk[g][a][b]	+=		e[p][a][b]*Friction_Tnsr_rd.comp[g][p];		

					}													
				}
			}
		}

		double g_clst_ijk[3][3][3] = {{{0}}};

	for (int a=0; a<3; a++)
		{
		for (int b=0; b<3; b++)
			{
			for (int g=0; g<3; g++)
				{
					g_clst_ijk[g][a][b] = 0.0;
					
				for (int p=0; p<5; p++)
					{
							g_clst_ijk[g][a][b]	+=		e[p][a][b]*Friction_Tnsr_td.comp[g][p];		

					}													
				}
			}
		}
cout<<"g_clst_ijk"<<endl;

for (int s=0; s<3; s++)
	{							
		cout << setw(10) << g_clst_ijk[s][0][0] << "  " << setw(10) << g_clst_ijk[s][0][1] << "  " << setw(10) << g_clst_ijk[s][0][2] << endl;
		cout << setw(10) << g_clst_ijk[s][1][0] << "  " << setw(10) << g_clst_ijk[s][1][1] << "  " << setw(10) << g_clst_ijk[s][1][2] << endl;
		cout << setw(10) << g_clst_ijk[s][2][0] << "  " << setw(10) << g_clst_ijk[s][2][1] << "  " << setw(10) << g_clst_ijk[s][2][2] << endl;
    }	
cout<<"h_clst_ijk"<<endl;

for (int s=0; s<3; s++)
	{							
		cout << setw(10) << h_clst_ijk[s][0][0] << "  " << setw(10) << h_clst_ijk[s][0][1] << "  " << setw(10) << h_clst_ijk[s][0][2] << endl;
		cout << setw(10) << h_clst_ijk[s][1][0] << "  " << setw(10) << h_clst_ijk[s][1][1] << "  " << setw(10) << h_clst_ijk[s][1][2] << endl;
		cout << setw(10) << h_clst_ijk[s][2][0] << "  " << setw(10) << h_clst_ijk[s][2][1] << "  " << setw(10) << h_clst_ijk[s][2][2] << endl;
    }
    	
		cout<<std::endl ;
		cout<<xi_11x11[0]<<'\t'<<xi_11x11[6]<<'\t'<<xi_11x11[12]<<'\t'<<xi_11x11[18]<<'\t'<<xi_11x11[24]<<'\t'<<xi_11x11[30]<<std::endl ;
		cout<<xi_11x11[1]<<'\t'<<xi_11x11[7]<<'\t'<<xi_11x11[13]<<'\t'<<xi_11x11[19]<<'\t'<<xi_11x11[25]<<'\t'<<xi_11x11[31]<<std::endl ;
		cout<<xi_11x11[2]<<'\t'<<xi_11x11[8]<<'\t'<<xi_11x11[14]<<'\t'<<xi_11x11[20]<<'\t'<<xi_11x11[26]<<'\t'<<xi_11x11[32]<<std::endl ;
		cout<<std::endl ;
		cout<<xi_11x11[3]<<'\t'<<xi_11x11[9]<<'\t'<<xi_11x11[15]<<'\t'<<xi_11x11[21]<<'\t'<<xi_11x11[27]<<'\t'<<xi_11x11[33]<<std::endl ;
		cout<<xi_11x11[4]<<'\t'<<xi_11x11[10]<<'\t'<<xi_11x11[16]<<'\t'<<xi_11x11[22]<<'\t'<<xi_11x11[28]<<'\t'<<xi_11x11[34]<<std::endl ;
		cout<<xi_11x11[5]<<'\t'<<xi_11x11[11]<<'\t'<<xi_11x11[17]<<'\t'<<xi_11x11[23]<<'\t'<<xi_11x11[29]<<'\t'<<xi_11x11[35]<<std::endl ;
		
	 for (int i=0;i<6;i++)
       { 
               for(int j=0;j<6;j++)
                       {
               cout<< std::setprecision(5) <<xi_11x11[i+6*j]<<'\t' ;
                       }
               cout<<std::endl; 
       }	

	 for (int i=0;i<5;i++)
       { 
               for(int j=0;j<5;j++)
                       {
               cout<< std::setprecision(5) <<Friction_Tnsr_dd.comp[i][j]<<'\t' ;
                       }
               cout<<std::endl; 
       }	
		
	inverse ( xi_11x11 , 6 )	 ; 			

	 for (int i=0;i<6;i++)
       { 
               for(int j=0;j<6;j++)
                       {
               cout<< std::setprecision(5) <<xi_11x11[i+6*j]<<'\t' ;
                       }
               cout<<std::endl; 
       }	
		
/*	for (int i=0; i<36; i++)
		{
			xi_11x11[i]*=4.1419e-14;	// multiply by kbT in erg K-1
		} 	

*/
		outFile1<<std::endl ;
		outFile1<<xi_11x11[0]<<'\t'<<xi_11x11[6]<<'\t'<<xi_11x11[12]<<'\t'<<xi_11x11[18]<<'\t'<<xi_11x11[24]<<'\t'<<xi_11x11[30]<<std::endl ;
		outFile1<<xi_11x11[1]<<'\t'<<xi_11x11[7]<<'\t'<<xi_11x11[13]<<'\t'<<xi_11x11[19]<<'\t'<<xi_11x11[25]<<'\t'<<xi_11x11[31]<<std::endl ;
		outFile1<<xi_11x11[2]<<'\t'<<xi_11x11[8]<<'\t'<<xi_11x11[14]<<'\t'<<xi_11x11[20]<<'\t'<<xi_11x11[26]<<'\t'<<xi_11x11[32]<<std::endl ;
		outFile1<<std::endl ;
		outFile1<<xi_11x11[3]<<'\t'<<xi_11x11[9]<<'\t'<<xi_11x11[15]<<'\t'<<xi_11x11[21]<<'\t'<<xi_11x11[27]<<'\t'<<xi_11x11[33]<<std::endl ;
		outFile1<<xi_11x11[4]<<'\t'<<xi_11x11[10]<<'\t'<<xi_11x11[16]<<'\t'<<xi_11x11[22]<<'\t'<<xi_11x11[28]<<'\t'<<xi_11x11[34]<<std::endl ;
		outFile1<<xi_11x11[5]<<'\t'<<xi_11x11[11]<<'\t'<<xi_11x11[17]<<'\t'<<xi_11x11[23]<<'\t'<<xi_11x11[29]<<'\t'<<xi_11x11[35]<<std::endl ;
/*
		outFile1<<std::endl ;
		outFile1<<xi_11x11[0]<<'\t'<<xi_11x11[11]<<'\t'<<xi_11x11[22]<<'\t'<<xi_11x11[33]<<'\t'<<xi_11x11[44]<<'\t'<<xi_11x11[55]<<std::endl ;
		outFile1<<xi_11x11[1]<<'\t'<<xi_11x11[12]<<'\t'<<xi_11x11[23]<<'\t'<<xi_11x11[34]<<'\t'<<xi_11x11[45]<<'\t'<<xi_11x11[56]<<std::endl ;
		outFile1<<xi_11x11[2]<<'\t'<<xi_11x11[13]<<'\t'<<xi_11x11[24]<<'\t'<<xi_11x11[35]<<'\t'<<xi_11x11[46]<<'\t'<<xi_11x11[57]<<std::endl ;
		outFile1<<std::endl ;
		outFile1<<xi_11x11[3]<<'\t'<<xi_11x11[14]<<'\t'<<xi_11x11[25]<<'\t'<<xi_11x11[36]<<'\t'<<xi_11x11[47]<<'\t'<<xi_11x11[58]<<std::endl ;
		outFile1<<xi_11x11[4]<<'\t'<<xi_11x11[15]<<'\t'<<xi_11x11[26]<<'\t'<<xi_11x11[37]<<'\t'<<xi_11x11[48]<<'\t'<<xi_11x11[59]<<std::endl ;
		outFile1<<xi_11x11[5]<<'\t'<<xi_11x11[16]<<'\t'<<xi_11x11[27]<<'\t'<<xi_11x11[38]<<'\t'<<xi_11x11[49]<<'\t'<<xi_11x11[60]<<std::endl ;
*/

// using the trick of matrix inversion by parts, since the Stresslet and flow-field switch going from FTS to FTE when doing dynamics of the aggregates
double mu_d[6][5];
double mu_dd[5][5];

			for (int l=0; l<6; l++)
				{
				for (int k=0; k<5; k++)
					{	
						mu_d[l][k] = 0.0;
				//		mu_d[l+3][k] = 0.0;
					for (int m=0; m<3; m++)
						{				
							// column major format
							mu_d[l][k]	-=	xi_11x11[l	+	6*m]*Friction_Tnsr_td.comp[m][k];
							mu_d[l][k]	-=	xi_11x11[l	+	6*(m+3)]*Friction_Tnsr_rd.comp[m][k];
						}
				//	mu_d[l][k] *= g_norm;
					}
				}
			for (int l=0; l<5; l++)
				{
				for (int k=0; k<5; k++)
					{	
						mu_dd[l][k] = Friction_Tnsr_dd.comp[l][k];
					for (int m=0; m<3; m++)
						{				
							// column major format
							mu_dd[l][k]	+=	Friction_Tnsr_dt.comp[l][m]*mu_d[m][k];
							mu_dd[l][k]	+=	Friction_Tnsr_dr.comp[l][m]*mu_d[m+3][k];
						}
					}
				}				
		cout << setw(10) << Friction_Tnsr_td.comp[0][0] << "  " << setw(10) << Friction_Tnsr_td.comp[0][1] << "  " << setw(10) << Friction_Tnsr_td.comp[0][2] << "  " << setw(10) << Friction_Tnsr_td.comp[0][3] << "  " << setw(10) << Friction_Tnsr_td.comp[0][4] << endl;
		cout << setw(10) << Friction_Tnsr_td.comp[1][0] << "  " << setw(10) << Friction_Tnsr_td.comp[1][1] << "  " << setw(10) << Friction_Tnsr_td.comp[1][2] << "  " << setw(10) << Friction_Tnsr_td.comp[1][3] << "  " << setw(10) << Friction_Tnsr_td.comp[1][4] << endl;
		cout << setw(10) << Friction_Tnsr_td.comp[2][0] << "  " << setw(10) << Friction_Tnsr_td.comp[2][1] << "  " << setw(10) << Friction_Tnsr_td.comp[2][2] << "  " << setw(10) << Friction_Tnsr_td.comp[2][3] << "  " << setw(10) << Friction_Tnsr_td.comp[2][4] << endl;
 
 		cout << setw(10) << Friction_Tnsr_rd.comp[0][0] << "  " << setw(10) << Friction_Tnsr_rd.comp[0][1] << "  " << setw(10) << Friction_Tnsr_rd.comp[0][2] << "  " << setw(10) << Friction_Tnsr_rd.comp[0][3] << "  " << setw(10) << Friction_Tnsr_rd.comp[0][4] << endl;
		cout << setw(10) << Friction_Tnsr_rd.comp[1][0] << "  " << setw(10) << Friction_Tnsr_rd.comp[1][1] << "  " << setw(10) << Friction_Tnsr_rd.comp[1][2] << "  " << setw(10) << Friction_Tnsr_rd.comp[1][3] << "  " << setw(10) << Friction_Tnsr_rd.comp[1][4] << endl;
		cout << setw(10) << Friction_Tnsr_rd.comp[2][0] << "  " << setw(10) << Friction_Tnsr_rd.comp[2][1] << "  " << setw(10) << Friction_Tnsr_rd.comp[2][2] << "  " << setw(10) << Friction_Tnsr_rd.comp[2][3] << "  " << setw(10) << Friction_Tnsr_rd.comp[2][4] << endl;
 

		cout<<std::endl ;
		cout<<mu_d[0][0]<<'\t'<<mu_d[0][1]<<'\t'<<mu_d[0][2]<<'\t'<<mu_d[0][3]<<'\t'<<mu_d[0][4]<<std::endl ;
		cout<<mu_d[1][0]<<'\t'<<mu_d[1][1]<<'\t'<<mu_d[1][2]<<'\t'<<mu_d[1][3]<<'\t'<<mu_d[1][4]<<std::endl ;
		cout<<mu_d[2][0]<<'\t'<<mu_d[2][1]<<'\t'<<mu_d[2][2]<<'\t'<<mu_d[2][3]<<'\t'<<mu_d[2][4]<<std::endl ;
		cout<<std::endl ;
		cout<<mu_d[3][0]<<'\t'<<mu_d[3][1]<<'\t'<<mu_d[3][2]<<'\t'<<mu_d[3][3]<<'\t'<<mu_d[3][4]<<std::endl ;
		cout<<mu_d[4][0]<<'\t'<<mu_d[4][1]<<'\t'<<mu_d[4][2]<<'\t'<<mu_d[4][3]<<'\t'<<mu_d[4][4]<<std::endl ;
		cout<<mu_d[5][0]<<'\t'<<mu_d[5][1]<<'\t'<<mu_d[5][2]<<'\t'<<mu_d[5][3]<<'\t'<<mu_d[5][4]<<std::endl ;			
		outFile1<<std::endl ;
		outFile1<<mu_d[0][0]<<'\t'<<mu_d[0][1]<<'\t'<<mu_d[0][2]<<'\t'<<mu_d[0][3]<<'\t'<<mu_d[0][4]<<std::endl ;
		outFile1<<mu_d[1][0]<<'\t'<<mu_d[1][1]<<'\t'<<mu_d[1][2]<<'\t'<<mu_d[1][3]<<'\t'<<mu_d[1][4]<<std::endl ;
		outFile1<<mu_d[2][0]<<'\t'<<mu_d[2][1]<<'\t'<<mu_d[2][2]<<'\t'<<mu_d[2][3]<<'\t'<<mu_d[2][4]<<std::endl ;
		outFile1<<std::endl ;
		outFile1<<mu_d[3][0]<<'\t'<<mu_d[3][1]<<'\t'<<mu_d[3][2]<<'\t'<<mu_d[3][3]<<'\t'<<mu_d[3][4]<<std::endl ;
		outFile1<<mu_d[4][0]<<'\t'<<mu_d[4][1]<<'\t'<<mu_d[4][2]<<'\t'<<mu_d[4][3]<<'\t'<<mu_d[4][4]<<std::endl ;
		outFile1<<mu_d[5][0]<<'\t'<<mu_d[5][1]<<'\t'<<mu_d[5][2]<<'\t'<<mu_d[5][3]<<'\t'<<mu_d[5][4]<<std::endl ;		
		cout<<std::endl ;
		cout<<mu_dd[0][0]<<'\t'<<mu_dd[0][1]<<'\t'<<mu_dd[0][2]<<'\t'<<mu_dd[0][3]<<'\t'<<mu_dd[0][4]<<std::endl ;
		cout<<mu_dd[1][0]<<'\t'<<mu_dd[1][1]<<'\t'<<mu_dd[1][2]<<'\t'<<mu_dd[1][3]<<'\t'<<mu_dd[1][4]<<std::endl ;
		cout<<mu_dd[2][0]<<'\t'<<mu_dd[2][1]<<'\t'<<mu_dd[2][2]<<'\t'<<mu_dd[2][3]<<'\t'<<mu_dd[2][4]<<std::endl ;
		cout<<mu_dd[3][0]<<'\t'<<mu_dd[3][1]<<'\t'<<mu_dd[3][2]<<'\t'<<mu_dd[3][3]<<'\t'<<mu_dd[3][4]<<std::endl ;
		cout<<mu_dd[4][0]<<'\t'<<mu_dd[4][1]<<'\t'<<mu_dd[4][2]<<'\t'<<mu_dd[4][3]<<'\t'<<mu_dd[4][4]<<std::endl ;
		outFile1<<std::endl ;
		outFile1<<mu_dd[0][0]<<'\t'<<mu_dd[0][1]<<'\t'<<mu_dd[0][2]<<'\t'<<mu_dd[0][3]<<'\t'<<mu_dd[0][4]<<std::endl ;
		outFile1<<mu_dd[1][0]<<'\t'<<mu_dd[1][1]<<'\t'<<mu_dd[1][2]<<'\t'<<mu_dd[1][3]<<'\t'<<mu_dd[1][4]<<std::endl ;
		outFile1<<mu_dd[2][0]<<'\t'<<mu_dd[2][1]<<'\t'<<mu_dd[2][2]<<'\t'<<mu_dd[2][3]<<'\t'<<mu_dd[2][4]<<std::endl ;
		outFile1<<mu_dd[3][0]<<'\t'<<mu_dd[3][1]<<'\t'<<mu_dd[3][2]<<'\t'<<mu_dd[3][3]<<'\t'<<mu_dd[3][4]<<std::endl ;
		outFile1<<mu_dd[4][0]<<'\t'<<mu_dd[4][1]<<'\t'<<mu_dd[4][2]<<'\t'<<mu_dd[4][3]<<'\t'<<mu_dd[4][4]<<std::endl ;
		for (int a=0; a<3; a++)
		{
		for (int b=0; b<3; b++)
			{
			for (int g=0; g<3; g++)
				{
					h_clst_ijk[g][a][b] = 0.0;
					
				for (int p=0; p<5; p++)
					{
							h_clst_ijk[g][a][b]	+=		e_l[p][a][b]*mu_d[g][p];		

					}													
				}
			}
		}


	for (int a=0; a<3; a++)
		{
		for (int b=0; b<3; b++)
			{
			for (int g=0; g<3; g++)
				{
					g_clst_ijk[g][a][b] = 0.0;
					
				for (int p=0; p<5; p++)
					{
							g_clst_ijk[g][a][b]	+=		e_l[p][a][b]*mu_d[g+3][p];		

					}													
				}
			}
		}
cout<<"g_clst_ijk"<<endl;

for (int s=0; s<3; s++)
	{							
		cout << setw(10) << g_clst_ijk[s][0][0] << "  " << setw(10) << g_clst_ijk[s][0][1] << "  " << setw(10) << g_clst_ijk[s][0][2] << endl;
		cout << setw(10) << g_clst_ijk[s][1][0] << "  " << setw(10) << g_clst_ijk[s][1][1] << "  " << setw(10) << g_clst_ijk[s][1][2] << endl;
		cout << setw(10) << g_clst_ijk[s][2][0] << "  " << setw(10) << g_clst_ijk[s][2][1] << "  " << setw(10) << g_clst_ijk[s][2][2] << endl;
    }	
cout<<"h_clst_ijk"<<endl;

for (int s=0; s<3; s++)
	{							
		cout << setw(10) << h_clst_ijk[s][0][0] << "  " << setw(10) << h_clst_ijk[s][0][1] << "  " << setw(10) << h_clst_ijk[s][0][2] << endl;
		cout << setw(10) << h_clst_ijk[s][1][0] << "  " << setw(10) << h_clst_ijk[s][1][1] << "  " << setw(10) << h_clst_ijk[s][1][2] << endl;
		cout << setw(10) << h_clst_ijk[s][2][0] << "  " << setw(10) << h_clst_ijk[s][2][1] << "  " << setw(10) << h_clst_ijk[s][2][2] << endl;
    }	

}
void Collision(vector<SubData>& particle, vector<ParticleData>& cluster, int i, int j,int *Max_Cluster_N, vctr3D box, vctr3D rbox ) {

	/*	int cluster_swap ;
		
		if (i > j) {
			cluster_swap = i;
			i = j;
			j = cluster_swap;
			}
	*/
		vctr3D L_after, L_before ;
		vctr3D old_pos= cluster[i].pos;
		vctr3D dr;
		double temp_r, r;
		
		L_before= ((cluster[i].pos^cluster[i].vel)*cluster[i].mass + cluster[i].Iner_tnsr*cluster[i].omega +
				   (cluster[j].pos.revPBC(old_pos,box,rbox)^cluster[j].vel)*cluster[j].mass + cluster[j].Iner_tnsr*cluster[j].omega );
	//	L_before[0] = ( mass[i]*(Cluster_XYZ[i][1]*Velocity[i][2]-Cluster_XYZ[i][2]*Velocity[i][2])+ I[i][0][0]*Ang_Velocity[i][0]+I[i][0][1]*Ang_Velocity[i][1]+I[i][0][2]*Ang_Velocity[i][2] +	mass[j]*(Cluster_XYZ[j][1]*Velocity[j][2]-Cluster_XYZ[j][2]*Velocity[j][1])+ I[j][0][0]*Ang_Velocity[j][0]+I[j][0][1]*Ang_Velocity[j][1]+I[j][0][2]*Ang_Velocity[j][2] );
	//	L_before[1] = ( mass[i]*(Cluster_XYZ[i][2]*Velocity[i][0]-Cluster_XYZ[i][0]*Velocity[i][2])+ I[i][1][0]*Ang_Velocity[i][0]+I[i][1][1]*Ang_Velocity[i][1]+I[i][1][2]*Ang_Velocity[i][2] +	mass[j]*(Cluster_XYZ[j][2]*Velocity[j][0]-Cluster_XYZ[j][0]*Velocity[j][2])+ I[j][1][0]*Ang_Velocity[j][0]+I[j][1][1]*Ang_Velocity[j][1]+I[j][1][2]*Ang_Velocity[j][2] );
	//	L_before[2] = ( mass[i]*(Cluster_XYZ[i][0]*Velocity[i][1]-Cluster_XYZ[i][1]*Velocity[i][0])+ I[i][2][0]*Ang_Velocity[i][0]+I[i][2][1]*Ang_Velocity[i][1]+I[i][2][2]*Ang_Velocity[i][2] +	mass[j]*(Cluster_XYZ[j][0]*Velocity[j][1]-Cluster_XYZ[j][1]*Velocity[j][0])+ I[j][2][0]*Ang_Velocity[j][0]+I[j][2][1]*Ang_Velocity[j][1]+I[j][2][2]*Ang_Velocity[j][2] );

/*	L[i][0] = (I[i][0][0]*Ang_Velocity[i][0]+I[i][0][1]*Ang_Velocity[i][1]+I[i][0][2]*Ang_Velocity[i][2] + I[j][0][0]*Ang_Velocity[j][0]+I[j][0][1]*Ang_Velocity[j][1]+I[j][0][2]*Ang_Velocity[j][2] );
		L[i][1] = (I[i][1][0]*Ang_Velocity[i][0]+I[i][1][1]*Ang_Velocity[i][1]+I[i][1][2]*Ang_Velocity[i][2] + I[j][1][0]*Ang_Velocity[j][0]+I[j][1][1]*Ang_Velocity[j][1]+I[j][1][2]*Ang_Velocity[j][2] );
		L[i][2] = (I[i][2][0]*Ang_Velocity[i][0]+I[i][2][1]*Ang_Velocity[i][1]+I[i][2][2]*Ang_Velocity[i][2] + I[j][2][0]*Ang_Velocity[j][0]+I[j][2][1]*Ang_Velocity[j][1]+I[j][2][2]*Ang_Velocity[j][2] );
*/			
		cluster[i].pos=(cluster[i].pos*cluster[i].mass + cluster[j].pos*cluster[j].mass ) * (1.0/(cluster[i].mass+cluster[j].mass));		

		cluster[i].mass=cluster[i].mass+cluster[j].mass;
		
		L_after = ((cluster[i].pos^cluster[i].vel)*cluster[i].mass) ;
		
		cluster[i].clicked = 1 ; 
		
		cluster[i].angmom=(L_before - L_after);

	//	std::cout<<i<<'\t'<<j<<'\t'<<cluster[i].mass<<std::endl;

 
	//	std::cout<<"inside_collision"<<'\t'<<i<<'\t'<<A[i][0][0]<<'\t'<<Ang_Velocity[i][0]<<'\t'<<Q[i][0]<<'\t'<<L[i][0]<<'\t'<<torque[i][0]<<std::endl;
		
		remove("cluster_ID.dat");

		std::ofstream outFile70("cluster_ID.dat");

		cluster[i].radii_gyr=0.0;		     	     	
		
		for (int  k=0; k<cluster[i].Sub_Length; k++) {

		particle[cluster[i].sub[k]].pos_bdyfxd	 =  particle[cluster[i].sub[k]].pos-cluster[i].pos;
		particle[cluster[i].sub[k]].pos_bdyfxd.PBC(box,rbox);
	    cluster[i].radii_gyr+=particle[cluster[i].sub[k]].pos_bdyfxd.norm2()/(cluster[i].Sub_Length+cluster[j].Sub_Length);				
		outFile70<<cluster[i].sub[k]<<'\t'<<i<<endl;
		
		} 
		
		for (int  k=cluster[i].Sub_Length; k<cluster[i].Sub_Length+cluster[j].Sub_Length; k++) {
		
	//	cluster[i].sub[k]	=	cluster[j].sub[k-cluster[i].Sub_Length];
			cluster[i].sub[k] = cluster[j].sub[k-cluster[i].Sub_Length];
			particle[cluster[i].sub[k]].pos_bdyfxd	 =  particle[cluster[i].sub[k]].pos-cluster[i].pos;
			particle[cluster[i].sub[k]].pos_bdyfxd.PBC(box,rbox);				
			cluster[i].radii_gyr+=particle[cluster[i].sub[k]].pos_bdyfxd.norm2()/(cluster[i].Sub_Length+cluster[j].Sub_Length);		
			outFile70<<cluster[i].sub[k]<<'\t'<<j<<endl;
			particle[cluster[i].sub[k]].cluster=i;
			
		} 
		outFile70.close();
		cluster[i].pos.PBC(box,rbox);	
		cluster[i].Sub_Length=cluster[i].Sub_Length+cluster[j].Sub_Length;

		cluster[i].radii=0;
		
		for ( int l = 0 ; l < cluster[i].Sub_Length ; l ++ )
			{
				for ( int k = l+1 ; k < cluster[i].Sub_Length ; k ++ )
					{
						dr=(particle[cluster[i].sub[l]].pos_bdyfxd -particle[cluster[i].sub[k]].pos_bdyfxd );
					//	dr.PBC(box,rbox);
						temp_r= dr.norm2();
						r	=	 0.56	+	sqrt(temp_r);
						if(r>cluster[i].radii)
							{
								cluster[i].radii	= r;
							} 
					}
			}
					
		if (cluster[i].radii>max_size)
			{
				cout <<"Radius"<<'\t'<<cluster[i].radii<< endl;
				cout <<"Cluster No:"<<'\t'<<i<<'\t'<<", Nr. of particles in the cluster"<<'\t'<<cluster[i].Sub_Length<< endl;
				cout << "*** cluster reached maximum allowed size " << endl;
				abort();
			}	
		
		cluster[i].radii_gyr=sqrt(cluster[i].radii_gyr + 0.15/cluster[i].Sub_Length);  	// volume correction term for single spheres from paper Improved Calculation of Rotational Diffusion and Intrinsic Viscosity of Bead Models for
																						// Macromolecules and Nanoparticles , J. Garcı´a de la TorreJ. Phys. Chem. B 2007, 111, 955-961 955


		// std::cout<<Length_cluster[i]<<'\t'<<Length_cluster[j]<<std::endl;
		for ( int k = j ; k < *Max_Cluster_N-1 ; k ++ )
			{
		//	 cluster[k]=cluster[k+1];
			cluster[k].clicked =cluster[k+1].clicked;
			cluster[k].pos=cluster[k+1].pos;
			cluster[k].frc=cluster[k+1].frc;
			cluster[k].omega=cluster[k+1].omega;
			cluster[k].rotmat=cluster[k+1].rotmat;
			cluster[k].Iner_tnsr=cluster[k+1].Iner_tnsr;
			cluster[k].mobility_tnsr=cluster[k+1].mobility_tnsr;
			cluster[k].mobility_tnsr_sqrt=cluster[k+1].mobility_tnsr_sqrt;
			cluster[k].rot_mobility_tnsr=cluster[k+1].rot_mobility_tnsr;
			cluster[k].rot_mobility_tnsr_sqrt=cluster[k+1].rot_mobility_tnsr_sqrt;
			cluster[k].trq=cluster[k+1].trq;
			cluster[k].angmom=cluster[k+1].angmom;
			cluster[k].quat=cluster[k+1].quat;
			cluster[k].radii=cluster[k+1].radii;
			cluster[k].mass=cluster[k+1].mass;
 			cluster[k].Sub_Length=cluster[k+1].Sub_Length; 
			cluster[k].radii_gyr=cluster[k+1].radii_gyr;
			for ( int l = 0 ; l < cluster[k].Sub_Length ; l ++ )
			{
			cluster[k].sub[l]			=	cluster[k+1].sub[l];
			particle[cluster[k].sub[l]].cluster=k;
			} 
		}
				
		// modulus of A is one; A inverse is jsut the cofactor of A 
	
		*Max_Cluster_N=*Max_Cluster_N-1;
		if (*Max_Cluster_N<=1) {
		*Max_Cluster_N=1;
		cout<<"all particles combined into Single Cluster "<<endl;
	//	abort();
		} 

	}
std::random_device seed;
std::mt19937 gen{seed()};
std::normal_distribution<> R1(0.0,1.0),R2(0.0,1.0),R3(0.0,1.0),R4(0.0,1.0),R5(0.0,1.0),R6(0.0,1.0);

void brownian( int step , vector<ParticleData>& cluster, vector<SubData>& particle, int *Max_Cluster_N , double *KE_rot, double vel_scale) {
double a, b , c, lambda;
vctr4D quat_old;
*KE_rot=0;
double e[5][3][3]= {
						{{1.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,-1.0}},
						{{0.0,1.0,0.0},{1.0,0.0,0.0},{0.0,0.0,0.0}},
						{{0.0,0.0,1.0},{0.0,0.0,0.0},{1.0,0.0,0.0}},
						{{0.0,0.0,0.0},{0.0,0.0,1.0},{0.0,1.0,0.0}},
						{{0.0,0.0,0.0},{0.0,1.0,0.0},{0.0,0.0,-1.0}}
					};
	
double e_l[5][3][3]= {
							{{ 2.0/3.0,0.0,0.0},{0.0,-1.0/3.0,0.0},{0.0,0.0,-1.0/3.0}},
							{{0.0,0.5,0.0},{0.5,0.0,0.0},{0.0,0.0,0.0}},
							{{0.0,0.0,0.5},{0.0,0.0,0.0},{0.5,0.0,0.0}},
							{{0.0,0.0,0.0},{0.0,0.0,0.5},{0.0,0.5,0.0}},
							{{-1.0/3.0,0.0,0.0},{0.0, 2.0/3.0,0.0},{0.0,0.0,-1.0/3.0}}
						};
						    						
for(int i=0;i<*Max_Cluster_N;i++) 
	{
		vctr3D rand(R1(gen), R2(gen), R3(gen));
		vctr3D rand1(R4(gen), R5(gen), R6(gen));
		vctr3D u_inf(0.0,shear_rate*cluster[i].pos.comp[0],0.0); 		// shear flow gradient in y-direction
		mtrx3D E_inf_b = (~cluster[i].rotmat)*E_inf*cluster[i].rotmat;
		vctr5D E_inf_bt;
		mtrx3D S_b ;
		mtrx3D S_s ;
		for(int j=0;j<5;j++) 
			{		
				E_inf_bt.comp[j] = 0.0;
				
				for(int k=0;k<3;k++) 
					{
						for(int l=0;l<3;l++) 
							{
								E_inf_bt.comp[j]	+=	 e[j][k][l]*E_inf_b.comp[k][l];
							}
					}
			}	

		 
		if (cluster[i].Sub_Length>0) 
			{
				cluster[i].pos+=cluster[i].rotmat*cluster[i].mobility_tnsr*(~cluster[i].rotmat)*(cluster[i].frc*dt) 
							//	+cluster[i].rotmat*cluster[i].mobility_tnsr_tr*(~cluster[i].rotmat)*(w_inf*10000.0*dt)
								/*+ cluster[i].rotmat*cluster[i].mobility_tnsr_sqrt*(rand*kbT_dt)*/
								+u_inf*dt-cluster[i].rotmat*(cluster[i].mobility_tnsr_td*E_inf_bt)*dt ;
				
				for(int m=0;m<5;m++) 
					{
						cluster[i].Stresslet.comp[m]=0.0;			
					for(int n=0;n<5;n++) 
						{	
						cluster[i].Stresslet.comp[m]-=	cluster[i].mobility_tnsr_dd.comp[m][n]*E_inf_bt.comp[n]	;
						}	
					}
					
				for(int k=0;k<3;k++) 
					{
						for(int l=0;l<3;l++) 
						{
								S_b.comp[k][l] = 0.0;
								
							for(int j=0;j<5;j++) 
							{						
								S_b.comp[k][l]	+=	 e[j][k][l]*cluster[i].Stresslet.comp[j];
							}
						}
					}
					
					S_s = (cluster[i].rotmat)*S_b*(~cluster[i].rotmat);		
					
					for(int m=0;m<5;m++) 
						{		
						cluster[i].Stresslet.comp[m]=0.0;			
							
							for(int k=0;k<3;k++) 
								{
									for(int l=0;l<3;l++) 
										{
											cluster[i].Stresslet.comp[m]	+=	 e_l[m][k][l]*S_s.comp[k][l];
										}
								}
						}										
	
				if(xx_rotation)	
				{
			// update Q
				quat_old=cluster[i].quat;
				// translate space-fixed w*dt (i.e. theta) (3 dimensions) into qdot (4 dimensions).
				// based on the Wotuer's paper on An elementary singularity-free Rotational Brownian Dynamics algorithm for anisotropic particles 
				// J. Chem. Phys. 142, 114103 (2015)
				
				cluster[i].theta   	= 	cluster[i].rot_mobility_tnsr*(~cluster[i].rotmat)*(cluster[i].trq*dt)
									//	+cluster[i].rot_mobility_tnsr_rt*(~cluster[i].rotmat)*(w_inf*1.0E4*dt)
									//	+  cluster[i].rot_mobility_tnsr_sqrt*(rand1*kbT_dt*4.25)
										-  (cluster[i].mobility_tnsr_rd*E_inf_bt)*dt; 	// body fixed omega
				cluster[i].omega	=	w_inf*dt;						// space-fixed omega
				cluster[i].quat		= cluster[i].theta2quat() + cluster[i].omega2qdot() ;
			//	cout<<cluster[i].theta.comp[0]<<cluster[i].theta.comp[1]<<cluster[i].theta.comp[2]<<endl;
			// lagragian normalization of quaternions; see your notes;
			// after quaternion update you get new quaternion (say ~q) which non-normalised, i.e. |~q|!=1; 
			// assuming qi(t+dt) = ~qi + lambda*qi(t);
			// hence 	|qi(t+dt)| = |~qi + lambda*qi(t)| =1;
				a=1.0;
				b=cluster[i].quat*quat_old*2.0;
				c=cluster[i].quat.norm2()-1.0;
				lambda = (-b+sqrt(b*b-4.0*a*c))/(2.0*a);
				cluster[i].quat=cluster[i].quat+quat_old*lambda;
		
			// update A matrix
				}
				cluster[i].quat2rotmat();
	
				for (int j=0; j<cluster[i].Sub_Length; j++) 
					{
					//	particle[cluster[i].sub[j]].pos = cluster[i].pos + particle[cluster[i].sub[j]].pos_bdyfxd;
						particle[cluster[i].sub[j]].pos = cluster[i].pos + cluster[i].rotmat*particle[cluster[i].sub[j]].pos_bdyfxd;
						particle[cluster[i].sub[j]].pos.PBC(box,rbox);
					}
				cluster[i].pos.PBC(box,rbox);
			} 
			else 
			{
				cluster[i].radii	=	0.56;//rmin*0.5 ;		// radii of single particle is sqrt(rmin_x^2+rmin_y^2+rmin_z^2)
				cluster[i].pos+=cluster[i].frc*mu*dt+rand*mu_sqrt*kbT_dt;
				cluster[i].pos.PBC(box,rbox);
				*KE_rot += 	(cluster[i].omega)*(cluster[i].angmom)*0.5;	
				for (int j=0; j< cluster[i].Sub_Length; j++) 
					{
						particle[cluster[i].sub[j]].pos=cluster[i].pos;
					}
			}
	}		
}

int main() {


if(xxcluster_restart) {
	std::ifstream fin("random_device_state_new.txt");
  	fin >> gen;
	}
// current date/time based on current system
   time_t now = time(0);
   struct tm *ltm = localtime(&now);
   cout << "start time"<< '\t'
   << (ltm->tm_year + 1900) << '-'
   << (ltm->tm_mon + 1) << '-'
   <<  ltm->tm_mday << "\t"
   << ltm->tm_hour << ":"
   << ltm->tm_min << ":"
   << ltm->tm_sec << endl;

   cout << "\t : current Git Hash - extended version";
   system(" git log --pretty=format:'%H' -n 1 ");

int if_create_particles = xxcreate, ifrestart=xxrestart;
double tauT=0.1;
int cluster_combine;
double Temp=T0;
// double shear_rate = 0.0; //shear rate
int ifshear = 0;// set equal to 1 for shear
std::string dataFileName="../xxx",dataFileName_new="../xxxnew" ;
int Max_Cluster_N=NrParticles;
double simu_time=dt;
int step=0, nSteps=10000, frame=100;
int restart_frame_offset=0;
double vel_scale;
int if_Periodic =1;
std::cout<<'\n'<<cellx<<'\t'<<celly<<'\t'<<cellz<<std::endl;
double  T_Energy, K_Energy, P_Energy, p_energy=0;
vctr3D dR, dr2 , dr_vec;
double R, r2;
double dr=0.05; // step size for RDF calculation
// std::vector<int> RDF((int)  floor(sqrt((Lx/2)*(Lx/2)+(Ly/2)*(Ly/2)+(Lz/2)*(Lz/2)))/dr,0), RDF1((int)  floor(sqrt(Lx*Lx+Ly*Ly))/dr,0);
double KE_rot=0;
int NrSubs=NrParticles;

vector<SubData>  particle(NrParticles);
vector<ParticleData>  cluster( NrParticles, ParticleData(NrSubs) );
int combine_now=0;
int combine[NrParticles][4];

if(ifrestart)	{
	if(!xxcluster_restart)	{
std::string fileName=dataFileName+"/End_positions.dat";

//read x,y positions from XY.dat
std::ifstream dataFile(fileName);
if(!dataFile.good()) {
	std::cerr<<"Given file is corrupt tyu/n"<<std::endl;
}
else {
	std::cout<<"Reading X, Y, Z co-ordinates"<<std::endl;
    std::string line;
    for (int i=0;i<NrParticles;i++) {
		std::getline(dataFile,line);
    	std::istringstream currentLine(line);    
        currentLine >> particle[i].pos.comp[0];
        currentLine >> particle[i].pos.comp[1];
        currentLine >> particle[i].pos.comp[2];
    }    
}	
}  else { 
std::string fileName=dataFileName+"/End_Position_Full.xyz";

//read x,y positions from XY.dat
std::ifstream dataFile(fileName);
if(!dataFile.good()) {
	std::cerr<<"Given file is corrupt /n"<<std::endl;
}
else {
	std::cout<<"Reading X, Y, Z co-ordinates"<<std::endl;
    std::string line0;
    std::string line;
    std::string line1;
    std::string line2;
    std::string line3;
    std::string line4;
    std::string line5;
    std::string line6;
    std::string line7;
    int n=0;
    std::getline(dataFile,line0);
    std::istringstream currentLine0(line0);    
    currentLine0 >> Max_Cluster_N;
    currentLine0 >> restart_frame_offset;
    for (int i=0;i<Max_Cluster_N;i++) 
		{
			std::getline(dataFile,line3);
			std::istringstream currentLine3(line3);    
			currentLine3 >> cluster[i].Sub_Length;
			cluster[i].mass=cluster[i].Sub_Length;
			currentLine3 >> cluster[i].radii_gyr;
			currentLine3 >> cluster[i].pos.comp[0];
			currentLine3 >> cluster[i].pos.comp[1];
			currentLine3 >> cluster[i].pos.comp[2];
			std::getline(dataFile,line7);
			std::istringstream currentLine7(line7);
			currentLine7 >> cluster[i].quat.comp[0];
			currentLine7 >> cluster[i].quat.comp[1];
			currentLine7 >> cluster[i].quat.comp[2];
			currentLine7 >> cluster[i].quat.comp[3];

			if (cluster[i].Sub_Length>1) 
				{
					std::getline(dataFile,line);
					std::istringstream currentLine(line); 
					currentLine >> cluster[i].mobility_tnsr.comp[0][0];
					currentLine >> cluster[i].mobility_tnsr.comp[0][1];
					currentLine >> cluster[i].mobility_tnsr.comp[0][2]; 
					currentLine >> cluster[i].mobility_tnsr.comp[1][0];
					currentLine >> cluster[i].mobility_tnsr.comp[1][1];
					currentLine >> cluster[i].mobility_tnsr.comp[1][2];
					currentLine >> cluster[i].mobility_tnsr.comp[2][0];
					currentLine >> cluster[i].mobility_tnsr.comp[2][1];
					currentLine >> cluster[i].mobility_tnsr.comp[2][2];					
					std::getline(dataFile,line4);
					std::istringstream currentLine4(line4); 
					currentLine4 >> cluster[i].mobility_tnsr_sqrt.comp[0][0];
					currentLine4 >> cluster[i].mobility_tnsr_sqrt.comp[0][1];
					currentLine4 >> cluster[i].mobility_tnsr_sqrt.comp[0][2]; 
					currentLine4 >> cluster[i].mobility_tnsr_sqrt.comp[1][0];
					currentLine4 >> cluster[i].mobility_tnsr_sqrt.comp[1][1];
					currentLine4 >> cluster[i].mobility_tnsr_sqrt.comp[1][2];
					currentLine4 >> cluster[i].mobility_tnsr_sqrt.comp[2][0];
					currentLine4 >> cluster[i].mobility_tnsr_sqrt.comp[2][1];
					currentLine4 >> cluster[i].mobility_tnsr_sqrt.comp[2][2];
					if(xx_rotation)	
					{
					std::getline(dataFile,line5);
					std::istringstream currentLine5(line5); 
					currentLine5 >> cluster[i].rot_mobility_tnsr.comp[0][0];
					currentLine5 >> cluster[i].rot_mobility_tnsr.comp[0][1];
					currentLine5 >> cluster[i].rot_mobility_tnsr.comp[0][2]; 
					currentLine5 >> cluster[i].rot_mobility_tnsr.comp[1][0];
					currentLine5 >> cluster[i].rot_mobility_tnsr.comp[1][1];
					currentLine5 >> cluster[i].rot_mobility_tnsr.comp[1][2];
					currentLine5 >> cluster[i].rot_mobility_tnsr.comp[2][0];
					currentLine5 >> cluster[i].rot_mobility_tnsr.comp[2][1];
					currentLine5 >> cluster[i].rot_mobility_tnsr.comp[2][2];					
					std::getline(dataFile,line6);
					std::istringstream currentLine6(line6); 
					currentLine6 >> cluster[i].rot_mobility_tnsr_sqrt.comp[0][0];
					currentLine6 >> cluster[i].rot_mobility_tnsr_sqrt.comp[0][1];
					currentLine6 >> cluster[i].rot_mobility_tnsr_sqrt.comp[0][2]; 
					currentLine6 >> cluster[i].rot_mobility_tnsr_sqrt.comp[1][0];
					currentLine6 >> cluster[i].rot_mobility_tnsr_sqrt.comp[1][1];
					currentLine6 >> cluster[i].rot_mobility_tnsr_sqrt.comp[1][2];
					currentLine6 >> cluster[i].rot_mobility_tnsr_sqrt.comp[2][0];
					currentLine6 >> cluster[i].rot_mobility_tnsr_sqrt.comp[2][1];
					currentLine6 >> cluster[i].rot_mobility_tnsr_sqrt.comp[2][2];	}
				}
		for (int  j = 0 ; j < cluster[i].Sub_Length ; j ++ )
				{
					std::getline(dataFile,line1);
					std::istringstream currentLine1(line1);    
					cluster[i].sub[j]=n;
					n+=1;
					currentLine1 >> particle[cluster[i].sub[j]].pos.comp[0];
					currentLine1 >> particle[cluster[i].sub[j]].pos.comp[1];
					currentLine1 >> particle[cluster[i].sub[j]].pos.comp[2];
					std::getline(dataFile,line2);
					std::istringstream currentLine2(line2);  
					currentLine2 >> particle[cluster[i].sub[j]].pos_bdyfxd.comp[0];
					currentLine2 >> particle[cluster[i].sub[j]].pos_bdyfxd.comp[1];
					currentLine2 >> particle[cluster[i].sub[j]].pos_bdyfxd.comp[2];
					particle[cluster[i].sub[j]].mass=1.0;
					particle[cluster[i].sub[j]].radius=0.5;					
					particle[cluster[i].sub[j]].cluster=i;					
				}
    }
}
}
} else {
	std::string fileName="XYZ.dat";
if (if_create_particles) {
createInitialPosition_N_particles(fileName,NrParticles,Lx,Ly,Lz);
}
//read x,y positions from XY.dat
std::ifstream dataFile(fileName);
if(!dataFile.good()) {
	std::cerr<<"Given file is corrupt56 /n"<<std::endl;
}
else {
	std::cout<<"Reading X, Y, Z co-ordinates"<<std::endl;
    std::string line;
   

    for (int i=0;i<NrParticles;i++) {
		std::getline(dataFile,line);
    	std::istringstream currentLine(line);
        currentLine >> particle[i].pos.comp[1];
        currentLine >> particle[i].pos.comp[0];
        currentLine >> particle[i].pos.comp[2];
       cout<< particle[i].pos.comp[0]<<endl;
	//	std::getline(dataFile,line);        
	//	particle[i].pos.comp[2]+=10.0;
	//	particle[i].pos.comp[2]=particle[i].pos.comp[2]*(-1.0);
	//	particle[i].pos.comp[1]=particle[i].pos.comp[1]*(-1.0);
    }
}
if (xx_diffusion) {			// calculate the diffusion tensor for the particles read-in ; if starting with paticles of varies shape then initializing diffusion tensors
	
	/* sort particles into cells */
for ( int i = 0 ; i < Max_Cluster_N; i ++ )
	{
		if (!xxcluster_restart) 
			{
				cluster[i].Sub_Length=1;		// initially each cluster has size one
				cluster[i].mass=1.0;
				cluster[i].vel={0.0,0.0,0.0};

				// intialize Q, A matrix

				cluster[i].quat={1.0,0.0,0.0,0.0};
				cluster[i].quat2rotmat();
		}
	for ( int j = 0 ; j < cluster[i].Sub_Length ; j ++ )
		{ 
		if (!xxcluster_restart) {
			cluster[i].sub[j]=i;
			particle[cluster[i].sub[j]].cluster=i;
			particle[cluster[i].sub[j]].mass=cluster[i].mass;
			particle[cluster[i].sub[j]].radius=0.5;
			cluster[i].radii=0.56;
			// particle[i].pos is the position of cluster, and particle[i].sub[i].pos is the spaced fixed position of particles in the cluster; initially all clusters have 1 paricle per cluster, and position of cluster is same as position of spaced fixed sub-particle 
			particle[cluster[i].sub[j]].vel=cluster[i].vel;
			particle[cluster[i].sub[j]].pos_bdyfxd={0.0,0.0,0.0};//cluster[i].sub[j].pos;
			cluster[i].pos=particle[cluster[i].sub[j]].pos;
			cluster[i].omega={0.0,0.0,0.0};//{((double) rand()/(RAND_MAX)-0.5),((double) rand()/(RAND_MAX)-0.5),((double) rand()/(RAND_MAX)-0.5)};
			cluster[i].angmom={((double) rand()/(RAND_MAX)-0.5),((double) rand()/(RAND_MAX)-0.5),((double) rand()/(RAND_MAX)-0.5)};
			cluster[i].pos=particle[cluster[i].sub[j]].pos;
		} 				
	}
}
/*
		cluster[0].quat={0.9659,  0.0, 0.0,  0.2588	};
 
		// update A matrix

        cluster[0].quat2rotmat();     
        
		for (int  k=0; k< Max_Cluster_N; k++) {
			
			particle[k].pos = cluster[0].rotmat*particle[k].pos;
			}
			*/
	forceUpdate( particle, &p_energy, &combine_now , combine, 0);
cout<< combine_now<<"combine123"<<endl;
	if (combine_now>0) 
		{	
		//	cout<<combine_now<<endl;
			vector<vector<int>> temp_combine(combine_now+1,vector<int> (4)) ;
			for (int pn = 1; pn<=combine_now ; pn++) 
				{ 		
					for (int j = 0; j< 4 ; j ++) 
						{
							temp_combine[pn][j]=combine[pn][j];
						}
					//	cout<<pn<<'\t'<<temp_combine[pn][0]<<'\t'<<temp_combine[pn][1]<<'\t'<<"insdide main beroe sort"<<endl;
				}			
		
	if(combine_now>1) {	sort (temp_combine.begin()+1,temp_combine.end(), RowSort()); }
	
			/*	for (int pn = 1; pn<=combine_now ; pn++) 
				{ 		
						cout<<temp_combine[pn][0]<<'\t'<<temp_combine[pn][1]<<"insdide main after sort"<<endl;
				}	*/

	if(combine_now>1) {	
	int count=1;
	do
		{
			int j=1;
		do
			{
			if ((temp_combine[count][0]==temp_combine[count+j][0]) && (temp_combine[count][1]==temp_combine[count+j][1]) )
				{
					temp_combine.erase( temp_combine.begin() + count+ j );
					combine_now-=1;
					j-=1;
				}
				j+=1;
			}	while (j<=(combine_now-count));
			count=count+1;			
		} while (count<combine_now);
	}	
			/*	for (int pn = 1; pn<=combine_now ; pn++) 
				{ 		
						cout<<temp_combine[pn][0]<<'\t'<<temp_combine[pn][1]<<'\t'<<"insdide mainf after unique"<<endl;
				} */
	// collision detection
	 	//		cout<<combine_now<<endl;

	for ( int pn = 1 ; pn <=combine_now; pn ++ )
		{
			if(particle[temp_combine[pn][2]].cluster!=particle[temp_combine[pn][3]].cluster) {
 			Collision(particle, cluster, temp_combine[pn][0], temp_combine[pn][1], &Max_Cluster_N, box, rbox );
				
				for ( int pp = pn+1 ; pp <=combine_now; pp ++ )
					{
					 
							if ( temp_combine[pp][0]==temp_combine[pn][1] ) {
								if (temp_combine[pp][1] > temp_combine[pn][0]) {
								
									temp_combine[pp][0]=temp_combine[pn][0] ;
									
								}
								else {
										temp_combine[pp][0]=temp_combine[pp][1] ;
										temp_combine[pp][1]=temp_combine[pn][0] ;
									} 
									
								}
									
							if ( temp_combine[pp][1]==temp_combine[pn][1] ) {
								if (temp_combine[pp][0] < temp_combine[pn][0]) {
								
									temp_combine[pp][1]=temp_combine[pn][0] ;
									
								}
								else {
										temp_combine[pp][1]=temp_combine[pp][0] ;
										temp_combine[pp][0]=temp_combine[pn][0] ;
									} 
									
								}
					
						if (temp_combine[pp][0]>=temp_combine[pn][1]) {temp_combine[pp][0]-=1; } 
									
						if (temp_combine[pp][1]>=temp_combine[pn][1]) {temp_combine[pp][1]-=1; } 
		
					}
					
				}	
					//	sort (temp_combine.begin()+pn+1,temp_combine.end(), RowSort());

				/*	for (int px = 1; px<=combine_now ; pn++) 
				{ 		
						cout<<temp_combine[px][0]<<'\t'<<temp_combine[px][1]<<'\t'<<" after combine"<<endl;
				} */
			//	cout << "Doner one combine " << '\t'<< pn <<  endl;
		}
	cout <<cluster[0].pos.comp[0] << '\t'<<cluster[0].pos.comp[1] << '\t'<<cluster[0].pos.comp[2] << '\t'<<endl;
// calculate new diffusion tensors	
	for ( int i = 0 ; i < Max_Cluster_N; i ++ )
		{
			if(cluster[i].clicked == 1 ) {
		remove("new_cluster.dat");
		remove("data.dat");

		std::ofstream outFile7("new_cluster.dat");
		
	//	 cluster[i].pos.comp[1]-=5.1861;
		/*	cluster[i].pos.comp[0]=0.0;
			cluster[i].pos.comp[1]=0.0;
			cluster[i].pos.comp[2]=0.0;	
		*/
		// cluster[i].pos={0.0,7.86223	,-3.76172e-06};	
		
	//	cluster[i].quat={0.9239, 0.0, 0.0, 0.3827	};
	//	cluster[i].quat={0.9659,  0.0, 0.0,  0.2588	};
 
		// update A matrix

        cluster[i].quat2rotmat();
        
       cluster[i].rotmat.echo();
        
        
		for (int  k=0; k<cluster[i].Sub_Length; k++) {
			
		//	particle[cluster[i].sub[k]].pos = cluster[i].rotmat*particle[cluster[i].sub[k]].pos;

		//	particle[cluster[i].sub[k]].pos_bdyfxd = particle[cluster[i].sub[k]].pos;
	//	particle[cluster[i].sub[k]].pos_bdyfxd.comp[1]+=5.1861;
	//	particle[cluster[i].sub[k]].pos_bdyfxd=particle[cluster[i].sub[k]].pos;

		outFile7<<particle[cluster[i].sub[k]].pos_bdyfxd.comp[0]<<'\t'<<particle[cluster[i].sub[k]].pos_bdyfxd.comp[1]<<'\t'<<particle[cluster[i].sub[k]].pos_bdyfxd.comp[2]<<'\t'<<"0.35"<<std::endl;
		
		} 
		/*
particle[cluster[0].sub[0]].pos_bdyfxd.comp[0]=0.0;
particle[cluster[0].sub[0]].pos_bdyfxd.comp[1]=-2.5147;
particle[cluster[0].sub[0]].pos_bdyfxd.comp[2]=0.0;
*/
// mobility calculation 


		mobility_calc(cluster[i].Sub_Length);

        vctr3D CoD;

/*	cluster[i].pos+=CoD;
    for (int  k=0; k<cluster[i].Sub_Length; k++) {
    particle[cluster[i].sub[k]].pos_bdyfxd-=CoD;
    }
*/    



std::ifstream dataFile("data.dat");
std::string tmp;
if(!dataFile.good()) {
	std::cerr<<"Given file is corrupt /n"<<std::endl;
}
else {
    std::string line;
		std::getline(dataFile,line);
    for (int n=0;n<3;n++) {
		std::getline(dataFile,line);
    	std::istringstream currentLine(line);    
        currentLine >> cluster[i].mobility_tnsr.comp[n][0];
        currentLine >> cluster[i].mobility_tnsr.comp[n][1];
        currentLine >> cluster[i].mobility_tnsr.comp[n][2];
        currentLine >> cluster[i].mobility_tnsr_tr.comp[n][0];
        currentLine >> cluster[i].mobility_tnsr_tr.comp[n][1];
        currentLine >> cluster[i].mobility_tnsr_tr.comp[n][2];
    }
		std::getline(dataFile,line);

    for (int n=0;n<3;n++) {
		std::getline(dataFile,line);
    	std::istringstream currentLine(line);    
        currentLine >> cluster[i].rot_mobility_tnsr_rt.comp[n][0];
        currentLine >> cluster[i].rot_mobility_tnsr_rt.comp[n][1];
        currentLine >> cluster[i].rot_mobility_tnsr_rt.comp[n][2];
        currentLine >> cluster[i].rot_mobility_tnsr.comp[n][0];
        currentLine >> cluster[i].rot_mobility_tnsr.comp[n][1];
        currentLine >> cluster[i].rot_mobility_tnsr.comp[n][2];
 
    }
		std::getline(dataFile,line);
    for (int n=0;n<3;n++) {
		std::getline(dataFile,line);
    	std::istringstream currentLine(line);    
        currentLine >> cluster[i].mobility_tnsr_td.comp[n][0];
        currentLine >> cluster[i].mobility_tnsr_td.comp[n][1];
        currentLine >> cluster[i].mobility_tnsr_td.comp[n][2];
        currentLine >> cluster[i].mobility_tnsr_td.comp[n][3];
        currentLine >> cluster[i].mobility_tnsr_td.comp[n][4];
    }
		std::getline(dataFile,line);

    for (int n=0;n<3;n++) {
		std::getline(dataFile,line);
    	std::istringstream currentLine(line);    
        currentLine >> cluster[i].mobility_tnsr_rd.comp[n][0];
        currentLine >> cluster[i].mobility_tnsr_rd.comp[n][1];
        currentLine >> cluster[i].mobility_tnsr_rd.comp[n][2];
        currentLine >> cluster[i].mobility_tnsr_rd.comp[n][3];
        currentLine >> cluster[i].mobility_tnsr_rd.comp[n][4];
 
    }
		std::getline(dataFile,line);

    for (int n=0;n<5;n++) {
		std::getline(dataFile,line);
    	std::istringstream currentLine(line);    
        currentLine >> cluster[i].mobility_tnsr_dd.comp[n][0];
        currentLine >> cluster[i].mobility_tnsr_dd.comp[n][1];
        currentLine >> cluster[i].mobility_tnsr_dd.comp[n][2];
        currentLine >> cluster[i].mobility_tnsr_dd.comp[n][3];
        currentLine >> cluster[i].mobility_tnsr_dd.comp[n][4];
 
    }    
}	 
     
      //  currentLine >> cluster[i].mobility_tnsr.comp[n][0];

      //  currentLine >> cluster[i].rot_mobility_tnsr.comp[n][0];


	//	cluster[i].mobility_tnsr=cluster[i].mobility_tnsr*(1*10*2414323832351.228);				// multiply by kBT (assuming kB in erg/K and T as 300 K ) correct for 1/kBT term included in the value 
	//	cluster[i].mobility_tnsr_tr=cluster[i].mobility_tnsr_tr*(1*10*2414323832351.228);				// multiply by kBT (assuming kB in erg/K and T as 300 K ) correct for 1/kBT term included in the value 
	//	cluster[i].rot_mobility_tnsr=cluster[i].rot_mobility_tnsr*(1*10*2414323832351.228);		// outputed by hydro++
	//	cluster[i].rot_mobility_tnsr_rt=cluster[i].rot_mobility_tnsr_rt*(1*10*2414323832351.228);		// outputed by hydro++
	//	cluster[i].mobility_tnsr_td=cluster[i].mobility_tnsr_td*(1*10*2414323832351.228);		// outputed by hydro++
	//	cluster[i].mobility_tnsr_rd=cluster[i].mobility_tnsr_rd*(1*10*2414323832351.228);		// outputed by hydro++
		cluster[i].mobility_tnsr_sqrt=null33D;
		MatrixXd temp(6,6), temp_sqrt(6,6);
		for (int k=0;k<3;k++) {
			for (int l=0;l<3;l++) {
				temp(k,l)=cluster[i].mobility_tnsr.comp[k][l];
				cout<<cluster[i].mobility_tnsr.comp[k][l]<<endl;

			}
		}
		for (int k=0;k<3;k++) {
			for (int l=3;l<6;l++) {
				temp(k,l)=cluster[i].mobility_tnsr_tr.comp[k][l-3];
				cout<<cluster[i].mobility_tnsr_tr.comp[k][l-3]<<endl;
			}
		}
		for (int k=3;k<6;k++) {
			for (int l=0;l<3;l++) {
				temp(k,l)=cluster[i].rot_mobility_tnsr_rt.comp[k-3][l];
				cout<<cluster[i].rot_mobility_tnsr_rt.comp[k-3][l]<<endl;
			}
		}
		for (int k=3;k<6;k++) {		
			for (int l=3;l<6;l++) {
				temp(k,l)=cluster[i].rot_mobility_tnsr.comp[k-3][l-3];
				cout<<cluster[i].rot_mobility_tnsr.comp[k-3][l-3]<<endl;
			}
		}
	Eigen::SelfAdjointEigenSolver<MatrixXd> TRANS_MOBL_MAT(temp);
	temp_sqrt = TRANS_MOBL_MAT.operatorSqrt();
				
		cout<<"mobility_tnsr_sqrt"<<endl;

		for (int k=0;k<3;k++) {
			for (int l=0;l<3;l++) {
				cluster[i].mobility_tnsr_sqrt.comp[k][l]=temp_sqrt(k,l);
				cout<<cluster[i].mobility_tnsr_sqrt.comp[k][l]<<endl;
			}
		}

		cout<<"mobility_tnsr_tr_sqrt"<<endl;
		
		for (int k=0;k<3;k++) {
			for (int l=3;l<6;l++) {
				cluster[i].mobility_tnsr_tr_sqrt.comp[k][l-3]=temp_sqrt(k,l);
				cout<<cluster[i].mobility_tnsr_tr_sqrt.comp[k][l-3]<<endl;
			}
		}

		cout<<"rot_mobility_tnsr_rt_sqrt"<<endl;

		for (int k=3;k<6;k++) {
			for (int l=0;l<3;l++) {
				cluster[i].rot_mobility_tnsr_rt_sqrt.comp[k-3][l]=temp_sqrt(k,l);
				cout<<cluster[i].rot_mobility_tnsr_rt_sqrt.comp[k-3][l]<<endl;
			}
		}

		cout<<"rot_mobility_tnsr_sqrt"<<endl;

		for (int k=3;k<6;k++) {
			for (int l=3;l<6;l++) {
				cluster[i].rot_mobility_tnsr_sqrt.comp[k-3][l-3]=temp_sqrt(k,l);
				cout<<cluster[i].rot_mobility_tnsr_sqrt.comp[k-3][l-3]<<endl;
			}
		}
				
		cluster[i].quat={1.0,0.0,0.0,0.0};
	//	cluster[i].quat={0.7133, -0.1100, -0.3141, -0.6168};

		// update A matrix

        cluster[i].quat2rotmat();
	}
		cluster[i].clicked = 0; 
	}
	
	}	
}
}
	for ( int i = 0 ; i < Max_Cluster_N; i ++ )
		{

std::ifstream dataFile("data.dat");
std::string tmp;
if(!dataFile.good()) {
	std::cerr<<"Given file is corrupt /n"<<std::endl;
}
else {

    std::string line;
		std::getline(dataFile,line);
    for (int n=0;n<3;n++) {
		std::getline(dataFile,line);
    	std::istringstream currentLine(line);    
        currentLine >> cluster[i].mobility_tnsr.comp[n][0];
        currentLine >> cluster[i].mobility_tnsr.comp[n][1];
        currentLine >> cluster[i].mobility_tnsr.comp[n][2];
        currentLine >> cluster[i].mobility_tnsr_tr.comp[n][0];
        currentLine >> cluster[i].mobility_tnsr_tr.comp[n][1];
        currentLine >> cluster[i].mobility_tnsr_tr.comp[n][2];
    }
		std::getline(dataFile,line);

    for (int n=0;n<3;n++) {
		std::getline(dataFile,line);
    	std::istringstream currentLine(line);    
        currentLine >> cluster[i].rot_mobility_tnsr_rt.comp[n][0];
        currentLine >> cluster[i].rot_mobility_tnsr_rt.comp[n][1];
        currentLine >> cluster[i].rot_mobility_tnsr_rt.comp[n][2];
        currentLine >> cluster[i].rot_mobility_tnsr.comp[n][0];
        currentLine >> cluster[i].rot_mobility_tnsr.comp[n][1];
        currentLine >> cluster[i].rot_mobility_tnsr.comp[n][2];
 
    }
		std::getline(dataFile,line);
    for (int n=0;n<3;n++) {
		std::getline(dataFile,line);
    	std::istringstream currentLine(line);    
        currentLine >> cluster[i].mobility_tnsr_td.comp[n][0];
        currentLine >> cluster[i].mobility_tnsr_td.comp[n][1];
        currentLine >> cluster[i].mobility_tnsr_td.comp[n][2];
        currentLine >> cluster[i].mobility_tnsr_td.comp[n][3];
        currentLine >> cluster[i].mobility_tnsr_td.comp[n][4];
    }
		std::getline(dataFile,line);

    for (int n=0;n<3;n++) {
		std::getline(dataFile,line);
    	std::istringstream currentLine(line);    
        currentLine >> cluster[i].mobility_tnsr_rd.comp[n][0];
        currentLine >> cluster[i].mobility_tnsr_rd.comp[n][1];
        currentLine >> cluster[i].mobility_tnsr_rd.comp[n][2];
        currentLine >> cluster[i].mobility_tnsr_rd.comp[n][3];
        currentLine >> cluster[i].mobility_tnsr_rd.comp[n][4];
 
    }
		std::getline(dataFile,line);

    for (int n=0;n<5;n++) {
		std::getline(dataFile,line);
    	std::istringstream currentLine(line);    
        currentLine >> cluster[i].mobility_tnsr_dd.comp[n][0];
        currentLine >> cluster[i].mobility_tnsr_dd.comp[n][1];
        currentLine >> cluster[i].mobility_tnsr_dd.comp[n][2];
        currentLine >> cluster[i].mobility_tnsr_dd.comp[n][3];
        currentLine >> cluster[i].mobility_tnsr_dd.comp[n][4];
 
    }    
}	 
     
      //  currentLine >> cluster[i].mobility_tnsr.comp[n][0];

      //  currentLine >> cluster[i].rot_mobility_tnsr.comp[n][0];


	//	cluster[i].mobility_tnsr=cluster[i].mobility_tnsr*(1*10*2414323832351.228);				// multiply by kBT (assuming kB in erg/K and T as 300 K ) correct for 1/kBT term included in the value 
	//	cluster[i].mobility_tnsr_tr=cluster[i].mobility_tnsr_tr*(1*10*2414323832351.228);				// multiply by kBT (assuming kB in erg/K and T as 300 K ) correct for 1/kBT term included in the value 
	//	cluster[i].rot_mobility_tnsr=cluster[i].rot_mobility_tnsr*(1*10*2414323832351.228);		// outputed by hydro++
	//	cluster[i].rot_mobility_tnsr_rt=cluster[i].rot_mobility_tnsr_rt*(1*10*2414323832351.228);		// outputed by hydro++
	//	cluster[i].mobility_tnsr_td=cluster[i].mobility_tnsr_td*(1*10*2414323832351.228);		// outputed by hydro++
	//	cluster[i].mobility_tnsr_rd=cluster[i].mobility_tnsr_rd*(1*10*2414323832351.228);		// outputed by hydro++
		cluster[i].mobility_tnsr_sqrt=null33D;
		MatrixXd temp(6,6), temp_sqrt(6,6);
		for (int k=0;k<3;k++) {
			for (int l=0;l<3;l++) {
				temp(k,l)=cluster[i].mobility_tnsr.comp[k][l];
				cout<<cluster[i].mobility_tnsr.comp[k][l]<<endl;

			}
		}
		for (int k=0;k<3;k++) {
			for (int l=3;l<6;l++) {
				temp(k,l)=cluster[i].mobility_tnsr_tr.comp[k][l-3];
				cout<<cluster[i].mobility_tnsr_tr.comp[k][l-3]<<endl;
			}
		}
		for (int k=3;k<6;k++) {
			for (int l=0;l<3;l++) {
				temp(k,l)=cluster[i].rot_mobility_tnsr_rt.comp[k-3][l];
				cout<<cluster[i].rot_mobility_tnsr_rt.comp[k-3][l]<<endl;
			}
		}
		for (int k=3;k<6;k++) {		
			for (int l=3;l<6;l++) {
				temp(k,l)=cluster[i].rot_mobility_tnsr.comp[k-3][l-3];
				cout<<cluster[i].rot_mobility_tnsr.comp[k-3][l-3]<<endl;
			}
		}
	Eigen::SelfAdjointEigenSolver<MatrixXd> TRANS_MOBL_MAT(temp);
	temp_sqrt = TRANS_MOBL_MAT.operatorSqrt();
				
		cout<<"mobility_tnsr_sqrt"<<endl;

		for (int k=0;k<3;k++) {
			for (int l=0;l<3;l++) {
				cluster[i].mobility_tnsr_sqrt.comp[k][l]=temp_sqrt(k,l);
				cout<<cluster[i].mobility_tnsr_sqrt.comp[k][l]<<endl;
			}
		}

		cout<<"mobility_tnsr_tr_sqrt"<<endl;
		
		for (int k=0;k<3;k++) {
			for (int l=3;l<6;l++) {
				cluster[i].mobility_tnsr_tr_sqrt.comp[k][l-3]=temp_sqrt(k,l);
				cout<<cluster[i].mobility_tnsr_tr_sqrt.comp[k][l-3]<<endl;
			}
		}

		cout<<"rot_mobility_tnsr_rt_sqrt"<<endl;

		for (int k=3;k<6;k++) {
			for (int l=0;l<3;l++) {
				cluster[i].rot_mobility_tnsr_rt_sqrt.comp[k-3][l]=temp_sqrt(k,l);
				cout<<cluster[i].rot_mobility_tnsr_rt_sqrt.comp[k-3][l]<<endl;
			}
		}

		cout<<"rot_mobility_tnsr_sqrt"<<endl;

		for (int k=3;k<6;k++) {
			for (int l=3;l<6;l++) {
				cluster[i].rot_mobility_tnsr_sqrt.comp[k-3][l-3]=temp_sqrt(k,l);
				cout<<cluster[i].rot_mobility_tnsr_sqrt.comp[k-3][l-3]<<endl;
			}
		}
}
//delete all files before writing data

// following snippet taken from stakcflow link  http://stackoverflow.com/questions/11007494/how-to-delete-all-files-in-a-folder-but-not-delete-the-folder-c-linux
if (xxcluster_restart) {
dataFileName=dataFileName_new;
}
const char *dataFileNamePointer = dataFileName.c_str();  // covnert the datafilename to a char pointer ans pass it to the snippet below which delete all files in that folder before running the simulation
if (!xxcluster_restart) {
struct dirent *next_file;
DIR *theFolder;
char filepath[256];
theFolder = opendir(dataFileNamePointer);
while (( next_file = readdir(theFolder)) )
	{
    // build the full path for each file in the folder
    sprintf(filepath, "%s/%s",dataFileNamePointer, next_file->d_name);
    if(strcmp(filepath,"../xxx/log")!=0)
		{
			remove(filepath);
		}
	}
//
}
if(!xx_diffusion) {
/* sort particles into cells */
for ( int i = 0 ; i < Max_Cluster_N; i ++ )
	{
		if (!xxcluster_restart) 
			{
				cluster[i].Sub_Length=1;		// initially each cluster has size one
				cluster[i].mass=1.0;
				cluster[i].vel={0.0,0.0,0.0};

				// intialize Q, A matrix

				cluster[i].quat={1.0,0.0,0.0,0.0};
				cluster[i].quat2rotmat();
		}
	for ( int j = 0 ; j < cluster[i].Sub_Length ; j ++ )
		{ 
		if (!xxcluster_restart) {
			cluster[i].sub[j]=i;
			particle[cluster[i].sub[j]].cluster=i;
			particle[cluster[i].sub[j]].mass=cluster[i].mass;
			particle[cluster[i].sub[j]].radius=0.5;
			cluster[i].radii=0.56;
			// particle[i].pos is the position of cluster, and particle[i].sub[i].pos is the spaced fixed position of particles in the cluster; initially all clusters have 1 paricle per cluster, and position of cluster is same as position of spaced fixed sub-particle 
			particle[cluster[i].sub[j]].vel=cluster[i].vel;
			particle[cluster[i].sub[j]].pos_bdyfxd={0.0,0.0,0.0};//cluster[i].sub[j].pos;
			cluster[i].pos=particle[cluster[i].sub[j]].pos;
			cluster[i].omega={0.0,0.0,0.0};//{((double) rand()/(RAND_MAX)-0.5),((double) rand()/(RAND_MAX)-0.5),((double) rand()/(RAND_MAX)-0.5)};
			cluster[i].angmom={((double) rand()/(RAND_MAX)-0.5),((double) rand()/(RAND_MAX)-0.5),((double) rand()/(RAND_MAX)-0.5)};
			cluster[i].pos=particle[cluster[i].sub[j]].pos;
		} 				
	}
}
}

std::ofstream outFile1(dataFileName+"/PE_energy.dat");
std::ofstream outFile10(dataFileName+"/End_positions.dat");
std::ofstream outFile11(dataFileName+"/no_of_clusters.dat");
std::ofstream outFile_com(dataFileName+"/com.dat");
std::ofstream outFile_orient(dataFileName+"/orient.dat");

// perfrom MD steps
/*	if (ifrestart) {
	simu_time =10.001;
	new_neighbor = TOPMAP(cellx,celly,cellz,if_Periodic,neighbor,box,(R_cut+R_shell),simu_time*shear_rate*Ly);
	std::tie(Pxx[0],Pyy[0],Pzz[0],Pxy[0],Pxz[0],Pyz[0])= forceUpdate(&p_energy ,XYZ, new_neighbor, cell, cellLength, box, (R_cut+R_shell), N, force, shear_rate, simu_time, ifshear , epsilon, sigma, rs);
}
*/

step = restart_frame_offset*frame+1;

forceUpdate( particle, &p_energy, &combine_now , combine, &step);

	// convert subforces into total generalized forces on particles 

  for ( int i = 0 ; i < Max_Cluster_N; i ++ )
  {
	cluster[i].frc=null3D;
	cluster[i].trq=null3D;
	cluster[i].Iner_tnsr=null33D;

    for (int  j = 0 ; j < cluster[i].Sub_Length ; j ++ )
    {
		dr_vec = particle[cluster[i].sub[j]].pos-cluster[i].pos;
		dr_vec.PBC(box,rbox);
		cluster[i].frc +=                                                  particle[cluster[i].sub[j]].frc;		
		cluster[i].trq +=                                               dr_vec^particle[cluster[i].sub[j]].frc;
		mtrx3D dr_mat(dr_vec*dr_vec.comp[0],dr_vec*dr_vec.comp[1],dr_vec*dr_vec.comp[2]);
		cluster[i].Iner_tnsr+=(I_sphere+Unit_diag*(dr_vec.norm2())-dr_mat)*particle[cluster[i].sub[j]].mass; 	//	refer following paper , page 3 equa. 3 for interia tensor formula
																												//	Modification of Numerical Model for Ellipsoidal Monomers by Erwin Gostomski
    }
  }



std::ofstream outFile8(dataFileName+"/logfile");

	outFile8 << "start time"<< '\t'
   	<< (ltm->tm_year + 1900) << '-'
   	<< (ltm->tm_mon + 1) << '-'
  	<<  ltm->tm_mday << "\t"
   	<< ltm->tm_hour << ":"
   	<< ltm->tm_min << ":"
   	<< ltm->tm_sec << endl;
	outFile8<<"Rotational Brownian On"<<'\t'<<xx_rotation<<std::endl;
	outFile8<<"NrParticles"<<'\t'<<NrParticles<<std::endl;
	outFile8<<"mass"<<'\t'<<m<<std::endl;
	outFile8<<"Volume fraction"<<'\t'<<vol_frac<<std::endl;
	outFile8<<"kb"<<'\t'<<kb<<std::endl;
	outFile8<<"Temperature (T0) ,"<<'\t'<<T0<<std::endl;
	outFile8<<"box size (abosute units)"<<'\t'<<box.comp[0]<<'\t'<<box.comp[1]<<'\t'<<box.comp[2]<<std::endl;
	outFile8<<"shear rate"<<'\t'<<shear_rate<<std::endl;
	outFile8<<"Cut-off for Interaction Potetntial , R_cut"<<'\t'<<r_cut<<std::endl;
	outFile8<<"Saturation point for Interaction Potetntial , rs"<<'\t'<<rs<<std::endl;
	outFile8<<"epsilon"<<'\t'<<epsilon<<std::endl;
	outFile8<<"sigma"<<'\t'<<sigma<<std::endl;
	outFile8<<"Timestep, dt"<<'\t'<<dt<<std::endl;
	outFile8<<"Viscosity, eta"<<'\t'<<eta<<std::endl;
	outFile8<<"Mobility , mu"<<'\t'<<mu<<std::endl;
	outFile8<<'\n'<<" Data Folder and Git Vesrion : "<<'\n';
	system(" echo >> logfile & git log --pretty=format:'%h' -n 1 >> logfile   & echo >> logfile  &  pwd >> logfile & ");
	outFile8.close();

	std::ofstream outFile_inter_cluster("inter_cluster_data.dat");
	std::ofstream Stresslet_data("Stresslet_data.dat");

simu_time =dt;
do {
	p_energy=0;	

	if(step == 60000){
		
		outFile_inter_cluster<<cluster[0].quat.comp[0]<<'\t'<<cluster[0].quat.comp[1]<<'\t'<<cluster[0].quat.comp[2]<<'\t'<<cluster[0].quat.comp[3]<<'\t'<<endl;
		outFile_inter_cluster<<cluster[0].pos.comp[0]<<'\t'<<cluster[0].pos.comp[1]<<'\t'<<cluster[0].pos.comp[2]<<'\t'<<endl;
		
		
	for (int  k=0; k<cluster[0].Sub_Length; k++) {

		outFile_inter_cluster<<particle[cluster[0].sub[k]].pos.comp[0]<<'\t'<<particle[cluster[0].sub[k]].pos.comp[1]<<'\t'<<particle[cluster[0].sub[k]].pos.comp[2]<<std::endl;
		
		} 
	}
	

	brownian(step, cluster, particle, &Max_Cluster_N , &KE_rot, vel_scale )	;
	combine_now=0;
 	forceUpdate( particle, &p_energy, &combine_now , combine, &step);
	if (xxclustering && combine_now>0) 
		{	
		//	cout<<combine_now<<endl;
			vector<vector<int>> temp_combine(combine_now+1,vector<int> (4)) ;
			for (int pn = 1; pn<=combine_now ; pn++) 
				{ 		
					for (int j = 0; j< 4 ; j ++) 
						{
							temp_combine[pn][j]=combine[pn][j];
						}
					//	cout<<pn<<'\t'<<temp_combine[pn][0]<<'\t'<<temp_combine[pn][1]<<'\t'<<"insdide main beroe sort"<<endl;
				}			
		
	if(combine_now>1) {	sort (temp_combine.begin()+1,temp_combine.end(), RowSort()); }
	
			/*	for (int pn = 1; pn<=combine_now ; pn++) 
				{ 		
						cout<<temp_combine[pn][0]<<'\t'<<temp_combine[pn][1]<<"insdide main after sort"<<endl;
				}	*/

	if(combine_now>1) {	
	int count=1;
	do
		{
			int j=1;
		do
			{
			if ((temp_combine[count][0]==temp_combine[count+j][0]) && (temp_combine[count][1]==temp_combine[count+j][1]) )
				{
					temp_combine.erase( temp_combine.begin() + count+ j );
					combine_now-=1;
					j-=1;
				}
				j+=1;
			}	while (j<=(combine_now-count));
			count=count+1;			
		} while (count<combine_now);
	}	
			/*	for (int pn = 1; pn<=combine_now ; pn++) 
				{ 		
						cout<<temp_combine[pn][0]<<'\t'<<temp_combine[pn][1]<<'\t'<<"insdide mainf after unique"<<endl;
				} */
	// collision detection
	 	//		cout<<combine_now<<endl;

	for ( int pn = 1 ; pn <=combine_now; pn ++ )
		{
			if(particle[temp_combine[pn][2]].cluster!=particle[temp_combine[pn][3]].cluster) {
 			Collision(particle, cluster, temp_combine[pn][0], temp_combine[pn][1], &Max_Cluster_N, box, rbox );
				
				for ( int pp = pn+1 ; pp <=combine_now; pp ++ )
					{
					 
							if ( temp_combine[pp][0]==temp_combine[pn][1] ) {
								if (temp_combine[pp][1] > temp_combine[pn][0]) {
								
									temp_combine[pp][0]=temp_combine[pn][0] ;
									
								}
								else {
										temp_combine[pp][0]=temp_combine[pp][1] ;
										temp_combine[pp][1]=temp_combine[pn][0] ;
									} 
									
								}
									
							if ( temp_combine[pp][1]==temp_combine[pn][1] ) {
								if (temp_combine[pp][0] < temp_combine[pn][0]) {
								
									temp_combine[pp][1]=temp_combine[pn][0] ;
									
								}
								else {
										temp_combine[pp][1]=temp_combine[pp][0] ;
										temp_combine[pp][0]=temp_combine[pn][0] ;
									} 
									
								}
					
						if (temp_combine[pp][0]>=temp_combine[pn][1]) {temp_combine[pp][0]-=1; } 
									
						if (temp_combine[pp][1]>=temp_combine[pn][1]) {temp_combine[pp][1]-=1; } 
		
					}
					
				}	
					//	sort (temp_combine.begin()+pn+1,temp_combine.end(), RowSort());

				/*	for (int px = 1; px<=combine_now ; pn++) 
				{ 		
						cout<<temp_combine[px][0]<<'\t'<<temp_combine[px][1]<<'\t'<<" after combine"<<endl;
				} */
			//	cout << "Doner one combine " << '\t'<< pn <<  endl;
		}
	
// calculate new diffusion tensors	
	for ( int i = 0 ; i < Max_Cluster_N; i ++ )
		{
			if(cluster[i].clicked == 1 ) {
		remove("new_cluster.dat");
		remove("data.dat");

		std::ofstream outFile7("new_cluster.dat");

		for (int  k=0; k<cluster[i].Sub_Length; k++) {

		outFile7<<particle[cluster[i].sub[k]].pos_bdyfxd.comp[0]<<'\t'<<particle[cluster[i].sub[k]].pos_bdyfxd.comp[1]<<'\t'<<particle[cluster[i].sub[k]].pos_bdyfxd.comp[2]<<'\t'<<"0.4"<<std::endl;
		
		} 

// mobility calculation 


		mobility_calc(cluster[i].Sub_Length);

        vctr3D CoD;

/*	cluster[i].pos+=CoD;
    for (int  k=0; k<cluster[i].Sub_Length; k++) {
    particle[cluster[i].sub[k]].pos_bdyfxd-=CoD;
    }
*/    



std::ifstream dataFile("data.dat");
std::string tmp;
if(!dataFile.good()) {
	std::cerr<<"Given file is corrupt /n"<<std::endl;
}
else {
    std::string line;
		std::getline(dataFile,line);
    for (int n=0;n<3;n++) {
		std::getline(dataFile,line);
    	std::istringstream currentLine(line);    
        currentLine >> cluster[i].mobility_tnsr.comp[n][0];
        currentLine >> cluster[i].mobility_tnsr.comp[n][1];
        currentLine >> cluster[i].mobility_tnsr.comp[n][2];
    }
		std::getline(dataFile,line);

    for (int n=0;n<3;n++) {
		std::getline(dataFile,line);
    	std::istringstream currentLine(line);    
        currentLine >> cluster[i].rot_mobility_tnsr.comp[n][0];
        currentLine >> cluster[i].rot_mobility_tnsr.comp[n][0];
        currentLine >> cluster[i].rot_mobility_tnsr.comp[n][0];
        currentLine >> cluster[i].rot_mobility_tnsr.comp[n][0];
        currentLine >> cluster[i].rot_mobility_tnsr.comp[n][1];
        currentLine >> cluster[i].rot_mobility_tnsr.comp[n][2];
 
    }
		std::getline(dataFile,line);
    for (int n=0;n<3;n++) {
		std::getline(dataFile,line);
    	std::istringstream currentLine(line);    
        currentLine >> cluster[i].mobility_tnsr_td.comp[n][0];
        currentLine >> cluster[i].mobility_tnsr_td.comp[n][1];
        currentLine >> cluster[i].mobility_tnsr_td.comp[n][2];
        currentLine >> cluster[i].mobility_tnsr_td.comp[n][3];
        currentLine >> cluster[i].mobility_tnsr_td.comp[n][4];
    }
		std::getline(dataFile,line);

    for (int n=0;n<3;n++) {
		std::getline(dataFile,line);
    	std::istringstream currentLine(line);    
        currentLine >> cluster[i].mobility_tnsr_rd.comp[n][0];
        currentLine >> cluster[i].mobility_tnsr_rd.comp[n][1];
        currentLine >> cluster[i].mobility_tnsr_rd.comp[n][2];
        currentLine >> cluster[i].mobility_tnsr_rd.comp[n][3];
        currentLine >> cluster[i].mobility_tnsr_rd.comp[n][4];
 
    }
}	 
     
      //  currentLine >> cluster[i].mobility_tnsr.comp[n][0];

      //  currentLine >> cluster[i].rot_mobility_tnsr.comp[n][0];


		cluster[i].mobility_tnsr=cluster[i].mobility_tnsr*(1*10*2414323832351.228);				// multiply by kBT (assuming kB in erg/K and T as 300 K ) correct for 1/kBT term included in the value 
		cluster[i].rot_mobility_tnsr=cluster[i].rot_mobility_tnsr*(1*10*2414323832351.228);		// outputed by hydro++
		cluster[i].mobility_tnsr_td=cluster[i].mobility_tnsr_td*(1*10*2414323832351.228);		// outputed by hydro++
		cluster[i].mobility_tnsr_rd=cluster[i].mobility_tnsr_rd*(1*10*2414323832351.228);		// outputed by hydro++
		cluster[i].mobility_tnsr_sqrt=null33D;
		MatrixXd temp(3,3), temp_sqrt(3,3);
		for (int k=0;k<3;k++) {
			for (int l=0;l<3;l++) {
				temp(k,l)=cluster[i].mobility_tnsr.comp[k][l];
				cout<<cluster[i].mobility_tnsr.comp[k][l]<<endl;

			}
		}
	Eigen::SelfAdjointEigenSolver<MatrixXd> TRANS_MOBL_MAT(temp);
	temp_sqrt = TRANS_MOBL_MAT.operatorSqrt();
	//	temp_sqrt=temp.sqrt();
		for (int k=0;k<3;k++) {
			for (int l=0;l<3;l++) {
				cluster[i].mobility_tnsr_sqrt.comp[k][l]=temp_sqrt(k,l);
		//		cluster[i].mobility_tnsr_sqrt.comp[k][k]=sqrt(cluster[i].mobility_tnsr.comp[k][k]);
				cout<<cluster[i].mobility_tnsr_sqrt.comp[k][k]<<endl;
			}
		}
		
		cluster[i].rot_mobility_tnsr_sqrt=null33D;
		for (int k=0;k<3;k++) {
			for (int l=0;l<3;l++) {
				temp(k,l)=cluster[i].rot_mobility_tnsr.comp[k][l];
				cout<<cluster[i].rot_mobility_tnsr.comp[k][l]<<endl;

			}
		}
		Eigen::SelfAdjointEigenSolver<MatrixXd> ROT_MOBL_MAT(temp);
		temp_sqrt = ROT_MOBL_MAT.operatorSqrt();
	//	temp_sqrt=temp.sqrt();
		for (int k=0;k<3;k++) {
			for (int l=0;l<3;l++) {
				cluster[i].rot_mobility_tnsr_sqrt.comp[k][l]=temp_sqrt(k,l);
			//	cluster[i].rot_mobility_tnsr_sqrt.comp[k][k]=sqrt(cluster[i].rot_mobility_tnsr.comp[k][k]);
				cout<<cluster[i].rot_mobility_tnsr_sqrt.comp[k][k]<<endl;
			}
		}
		
		cluster[i].quat={1.0,0.0,0.0,0.0};

		// update A matrix

        cluster[i].quat2rotmat();
	}
		cluster[i].clicked = 0; 
	}
	
	}

	// convert subforces into total generalized forces on particles 

  for ( int i = 0 ; i < Max_Cluster_N; i ++ )
  {
	cluster[i].frc=null3D;
	cluster[i].trq=null3D;
	cluster[i].Iner_tnsr=null33D;

    for (int  j = 0 ; j < cluster[i].Sub_Length ; j ++ )
    {
		dr_vec = particle[cluster[i].sub[j]].pos-cluster[i].pos;
		dr_vec.PBC(box,rbox);
		cluster[i].frc +=                                                  particle[cluster[i].sub[j]].frc;		
		cluster[i].trq +=                                               dr_vec^particle[cluster[i].sub[j]].frc;
		mtrx3D dr_mat(dr_vec*dr_vec.comp[0],dr_vec*dr_vec.comp[1],dr_vec*dr_vec.comp[2]);
		cluster[i].Iner_tnsr+=(I_sphere+Unit_diag*(dr_vec.norm2())-dr_mat)*particle[cluster[i].sub[j]].mass; 	//	refer following paper , page 3 equa. 3 for interia tensor formula
																												//	Modification of Numerical Model for Ellipsoidal Monomers by Erwin Gostomski
    } 
  } 
/*
if (step%(frame*10)==0)
        {
                cout<<step<<endl;
                std::ofstream outFile_inter_endfile(dataFileName+"/End_Position_Full_"+std::to_string(step/(frame*10))+".xyz");
                std::ofstream outFile_inter_rand_state(dataFileName+"/random_device_state_"+std::to_string(step/(frame*10))+".txt");

 outFile_inter_endfile<<'\t'<<Max_Cluster_N<<'\t'<<(int) (step/frame)<<endl;

for ( int i = 0 ; i < Max_Cluster_N; i ++ )
        {
                 outFile_inter_endfile<<cluster[i].Sub_Length<<'\t'<<cluster[i].radii_gyr<<'\t'<<cluster[i].pos.comp[0]<<'\t'<<cluster[i].pos.comp[1]<<'\t'<<cluster[i].pos.comp[2]<<std::endl;
                 outFile_inter_endfile<<cluster[i].quat.comp[0]<<'\t'<<cluster[i].quat.comp[1]<<'\t'<<cluster[i].quat.comp[2]<<'\t'<<cluster[i].quat.comp[3]<<std::endl;
                if (cluster[i].Sub_Length>1)
                        {
                                cluster[i].mobility_tnsr.writeToFile(outFile_inter_endfile);
                                cluster[i].mobility_tnsr_sqrt.writeToFile( outFile_inter_endfile);
                if(xx_rotation)
                        {
								cluster[i].rot_mobility_tnsr.writeToFile(outFile_inter_endfile);
								cluster[i].rot_mobility_tnsr_sqrt.writeToFile(outFile_inter_endfile);
                        }
                        }
            for (int  j = 0 ; j < cluster[i].Sub_Length ; j ++ )
                        {
                                 outFile_inter_endfile<<'\t'<<particle[cluster[i].sub[j]].pos.comp[0]<<'\t'<<particle[cluster[i].sub[j]].pos.comp[1]<<'\t'<<particle[cluster[i].sub[j]].pos.comp[2]<<std::endl;
                                 outFile_inter_endfile<<'\t'<<particle[cluster[i].sub[j]].pos_bdyfxd.comp[0]<<'\t'<<particle[cluster[i].sub[j]].pos_bdyfxd.comp[1]<<'\t'<<particle[cluster[i].sub[j]].pos_bdyfxd.comp[2]<<std::endl;
                        }
        }

        outFile_inter_rand_state << gen;
        outFile_inter_rand_state.close();
         outFile_inter_endfile.close();
        }
*/
if (step%frame==0) 
	{ 

      std::ofstream outFile5(dataFileName+"/XYZ"+ std::to_string(step/frame) +".xyz");   
   		outFile5<<NrParticles<<std::endl;
   		outFile5<<"X Y Z co-ordinates"<<std::endl;
		outFile11<<step<<'\t'<<Max_Cluster_N<<std::endl;
		// save position, Kinetic energy, Potential energy, Forces every 'frame' steps and also store radii of gyration info
		
	//	std::ofstream outFile9(dataFileName+"/Cluster_dist"+ std::to_string(step/frame) +".dat");

		K_Energy=0;

		for ( int i = 0 ; i < Max_Cluster_N; i ++ )
			{
				if(cluster[i].Sub_Length>0)
				{
				Stresslet_data<<cluster[i].Stresslet.comp[0]<<'\t'<<cluster[i].Stresslet.comp[1]<<'\t'<<cluster[i].Stresslet.comp[2]<<'\t'<<cluster[i].Stresslet.comp[3]<<'\t'<<cluster[i].Stresslet.comp[4]<<'\t'<<std::endl;	
	//			outFile9<<cluster[i].radii_gyr<<'\t'<<cluster[i].Sub_Length<<std::endl;
				outFile_com<<cluster[i].pos.comp[0]<<'\t'<<cluster[i].pos.comp[1]<<'\t'<<cluster[i].pos.comp[2]<<std::endl;
			/*	outFile_orient<<particle[cluster[i].sub[cluster[i].Sub_Length-1]].pos.comp[0]-particle[cluster[i].sub[0]].pos.comp[0]<<'\t'<<
				particle[cluster[i].sub[cluster[i].Sub_Length-1]].pos.comp[1]-particle[cluster[i].sub[0]].pos.comp[1]<<'\t'<<
				particle[cluster[i].sub[cluster[i].Sub_Length-1]].pos.comp[2]-particle[cluster[i].sub[0]].pos.comp[2]<<'\t'<<i<<std::endl; */
				outFile_orient<<cluster[i].rotmat.comp[0][1]<<'\t'<<cluster[i].rotmat.comp[1][1]<<'\t'<<cluster[i].rotmat.comp[2][1]<<'\t'<<std::endl;
				}
			    for (int  j = 0 ; j < cluster[i].Sub_Length ; j ++ )
					{
					
					outFile5<<'H'<<'\t'<<particle[cluster[i].sub[j]].pos.comp[0]<<'\t'<<particle[cluster[i].sub[j]].pos.comp[1]<<'\t'<<particle[cluster[i].sub[j]].pos.comp[2]<<'\t'<<i<<std::endl;
					}
			}


/*		for ( int i = 0 ; i < NrParticles; i ++ )
			{
						outFile5<<'H'<<'\t'<<particle[i].pos.comp[0]<<'\t'<<particle[i].pos.comp[1]<<'\t'<<particle[i].pos.comp[2]<<std::endl;
						K_Energy+=0.5*m*(particle[i].vel.comp[0]*particle[i].vel.comp[0]
									   + particle[i].vel.comp[1]*particle[i].vel.comp[1]
									   + particle[i].vel.comp[2]*particle[i].vel.comp[2]);
			}
 */
 
     	outFile5<<'\n'<<std::endl;
	//	outFile1<<p_energy<<std::endl;
		outFile5.close();
	//	outFile9.close();
/*		
		// store info to restart file End_Position_Full.xyz
		std::ofstream outFile7(dataFileName+"/End_Position_Full_new.xyz");
		std::ofstream outFile_rand_state(dataFileName+"/random_device_state_new.txt");

outFile7<<'\t'<<Max_Cluster_N<<'\t'<<(int) (step/frame)<<endl;

for ( int i = 0 ; i < Max_Cluster_N; i ++ )
	{
		outFile7<<cluster[i].Sub_Length<<'\t'<<cluster[i].radii_gyr<<'\t'<<cluster[i].pos.comp[0]<<'\t'<<cluster[i].pos.comp[1]<<'\t'<<cluster[i].pos.comp[2]<<std::endl;
		outFile7<<cluster[i].quat.comp[0]<<'\t'<<cluster[i].quat.comp[1]<<'\t'<<cluster[i].quat.comp[2]<<'\t'<<cluster[i].quat.comp[3]<<std::endl;
		if (cluster[i].Sub_Length>1) 
			{
				cluster[i].mobility_tnsr.writeToFile(outFile7);
				cluster[i].mobility_tnsr_sqrt.writeToFile(outFile7);
		if(xx_rotation)	
			{
				cluster[i].rot_mobility_tnsr.writeToFile(outFile7);
				cluster[i].rot_mobility_tnsr_sqrt.writeToFile(outFile7);
			}
			}
	    for (int  j = 0 ; j < cluster[i].Sub_Length ; j ++ )
			{
				outFile7<<'\t'<<particle[cluster[i].sub[j]].pos.comp[0]<<'\t'<<particle[cluster[i].sub[j]].pos.comp[1]<<'\t'<<particle[cluster[i].sub[j]].pos.comp[2]<<std::endl;
				outFile7<<'\t'<<particle[cluster[i].sub[j]].pos_bdyfxd.comp[0]<<'\t'<<particle[cluster[i].sub[j]].pos_bdyfxd.comp[1]<<'\t'<<particle[cluster[i].sub[j]].pos_bdyfxd.comp[2]<<std::endl;
			}
	}

	outFile_rand_state << gen;
	outFile_rand_state.close();

	outFile7.close();

 remove("random_device_state.txt");
  char oldname_rand_state[] ="random_device_state_new.txt";
  char newname_rand_state[] ="random_device_state.txt";
  rename( oldname_rand_state , newname_rand_state );

 remove("End_Position_Full.xyz");
  char oldname[] ="End_Position_Full_new.xyz";
  char newname[] ="End_Position_Full.xyz";
  rename( oldname , newname );
*/
	}
	
	simu_time+=dt;
	step+=1;

} while(xxnstep);

std::ofstream outFile7(dataFileName+"/End_Position_Full_new.xyz");
		std::ofstream outFile_rand_state(dataFileName+"/random_device_state_new.txt");
	outFile_rand_state << gen;
	outFile_rand_state.close();

outFile7<<'\t'<<Max_Cluster_N<<'\t'<<(int) (step/frame )<<endl;

for ( int i = 0 ; i < Max_Cluster_N; i ++ )
	{
		outFile7<<cluster[i].Sub_Length<<'\t'<<cluster[i].radii_gyr<<'\t'<<cluster[i].pos.comp[0]<<'\t'<<cluster[i].pos.comp[1]<<'\t'<<cluster[i].pos.comp[2]<<std::endl;
		outFile7<<cluster[i].quat.comp[0]<<'\t'<<cluster[i].quat.comp[1]<<'\t'<<cluster[i].quat.comp[2]<<'\t'<<cluster[i].quat.comp[3]<<std::endl;
		if (cluster[i].Sub_Length>1) 
			{
				cluster[i].mobility_tnsr.writeToFile(outFile7);
				cluster[i].mobility_tnsr_sqrt.writeToFile(outFile7);
		if(xx_rotation)	
			{
				cluster[i].rot_mobility_tnsr.writeToFile(outFile7);
				cluster[i].rot_mobility_tnsr_sqrt.writeToFile(outFile7);
			}
			}
	    for (int  j = 0 ; j < cluster[i].Sub_Length ; j ++ )
			{
				outFile7<<'\t'<<particle[cluster[i].sub[j]].pos.comp[0]<<'\t'<<particle[cluster[i].sub[j]].pos.comp[1]<<'\t'<<particle[cluster[i].sub[j]].pos.comp[2]<<std::endl;
				outFile7<<'\t'<<particle[cluster[i].sub[j]].pos_bdyfxd.comp[0]<<'\t'<<particle[cluster[i].sub[j]].pos_bdyfxd.comp[1]<<'\t'<<particle[cluster[i].sub[j]].pos_bdyfxd.comp[2]<<std::endl;
			}
	}
	for (int i=0;i<NrParticles;i++)
		{
			outFile10<<particle[i].pos.comp[0]<<'\t'<<particle[i].pos.comp[1]<<'\t'<<particle[i].pos.comp[2]<<std::endl;
		}
outFile1.close();
outFile7.close();
outFile10.close();
outFile11.close();
outFile_com.close();
outFile_orient.close();
  remove("End_Position_Full.xyz");
  char oldname[] ="End_Position_Full_new.xyz";
  char newname[] ="End_Position_Full.xyz";
  rename( oldname , newname );

     // get time now
	now = time(0);
	ltm = localtime(&now);
   cout << "end time"<< '\t'
   << (ltm->tm_year + 1900) << '-' 
   << (ltm->tm_mon + 1) << '-'
   <<  ltm->tm_mday << "\t"
   << ltm->tm_hour << ":"
   << ltm->tm_min << ":"
   << ltm->tm_sec << endl;

return 0;

}
