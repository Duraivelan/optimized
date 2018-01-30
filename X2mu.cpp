/*
 * Computes the 11x11 mobility matrix 
 * Input file : 
 * 				** init.dat **
 * 					viscoity 
 * 					radius
 * 					no. of beads 
 * 				** XYZ.dat **
 *					X,Y,Z particle positons
 *
 *  Output file : 
 * 				 ** data.dat **
 * 					contains the elements of the 11x11 mobility matrix in non-dimensional values; 
 * 					as seen in eqaution 326 of wouter notes; caution eta_s in his note equals 6*pi*eta_0 here  
 * 							
 * 			
 * Current version : assumes all particle radii to be same, also doesn't compute the lubrication forces
 * 
 * Based on the suggestions by following authors :
 * Excerpt from : 	Seto, Ryohei, et al. "Restructuring of colloidal aggregates in shear flow: Coupling interparticle contact models with Stokesian dynamics."
 * 					arXiv preprint arXiv:1204.5680 (2012). 		
 * 	" It must be noted that the lubrication correction of SD 
is not applied in this work. For suspensions where the interparticle
interaction is absent, the lubrication forces play
essential role for near contact particles [33, 34]. On the
other hand, for rigid clusters, i.e. if the relative velocities
between particles are zero due to strong cohesive forces,
the lubrication correction has no contribution. Thus, the
lubrication correction to the mobility matrix can safely be
omitted [35–37]. "	
 * 			
 *  */

# include "structure_definitions.h"
#include <sstream>
#include <math.h>       /* pow */

// llapack input variable structuring 

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
 
int main() {
	                 
vctr3D dR, dr2;
double R, r2;


// variables for mobility tensor calculation
double eta_0,eta_6pi;
vctr3D e_ab , e_ab_unit ;
double e_ab2, e_ab2_inv ;
double radius;
int NrParticles ; 
//

std::string fileName="init.dat";

//read viscoisty and x,y,z positions from new_cluster.dat
std::ifstream dataFile;
dataFile.open(fileName);
if(!dataFile.good()) {
	std::cerr<<"Given file is corrupt /n"<<std::endl;
}
else {
    std::string line0;
	std::getline(dataFile,line0);
   	std::istringstream currentLine0(line0);  
   	currentLine0 >> eta_0;
   	eta_6pi = eta_0*6.0*M_PI ; 
   	cout << eta_0 << endl;
	
	currentLine0.str("");
	currentLine0.clear(); // Clear state flags.
	std::getline(dataFile,line0);
	currentLine0.str(line0); 

   	currentLine0 >> radius;
   	cout << radius << endl;

	currentLine0.str("");
	currentLine0.clear(); // Clear state flags.
	std::getline(dataFile,line0);
	currentLine0.str(line0); 

	currentLine0 >> NrParticles;
   	cout << NrParticles << endl;

}

dataFile.close();  
dataFile.clear();

vector<SubData>  bead(NrParticles);

fileName="XYZ.dat";

//read viscoisty and x,y,z positions from new_cluster.dat

dataFile.open(fileName);

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
 //       bead[i].pos.comp[2] = bead[i].pos.comp[2]*(-1.0);
      //  bead[i].radius = radius ;
		bead[i].radius = 1.0 ; 	// bead radii hard-coded as 1.0 to avoid errors arising from non-dimensionalization procedure, but in principle it could be anything
	//	bead[i].pos = bead[i].pos - cntr ;
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
/*
	double e_g_S[5][3][3]= {
							{{1.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,-1.0}},
							{{0.0,1.0,0.0},{1.0,0.0,0.0},{0.0,0.0,0.0}},
							{{0.0,0.0,1.0},{0.0,0.0,0.0},{1.0,0.0,0.0}},
							{{0.0,0.0,0.0},{0.0,0.0,1.0},{0.0,1.0,0.0}},
							{{0.0,0.0,0.0},{0.0,1.0,0.0},{0.0,0.0,-1.0}}
						};


	double e_S_a[5][3][3]= {
							{{ 2.0/3.0,0.0,0.0},{0.0,-1.0/3.0,0.0},{0.0,0.0,-1.0/3.0}},
							{{0.0,0.5,0.0},{0.5,0.0,0.0},{0.0,0.0,0.0}},
							{{0.0,0.0,0.5},{0.0,0.0,0.0},{0.5,0.0,0.0}},
							{{0.0,0.0,0.0},{0.0,0.0,0.5},{0.0,0.5,0.0}},
							{{-1.0/3.0,0.0,0.0},{0.0, 2.0/3.0,0.0},{0.0,0.0,-1.0/3.0}}
						};
   

	double e_E_a[5][3][3]= {
							{{1.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,-1.0}},
							{{0.0,1.0,0.0},{1.0,0.0,0.0},{0.0,0.0,0.0}},
							{{0.0,0.0,1.0},{0.0,0.0,0.0},{1.0,0.0,0.0}},
							{{0.0,0.0,0.0},{0.0,0.0,1.0},{0.0,1.0,0.0}},
							{{0.0,0.0,0.0},{0.0,1.0,0.0},{0.0,0.0,-1.0}}
						};


	double e_g_E[5][3][3]= {
							{{ 2.0/3.0,0.0,0.0},{0.0,-1.0/3.0,0.0},{0.0,0.0,-1.0/3.0}},
							{{0.0,0.5,0.0},{0.5,0.0,0.0},{0.0,0.0,0.0}},
							{{0.0,0.0,0.5},{0.0,0.0,0.0},{0.5,0.0,0.0}},
							{{0.0,0.0,0.0},{0.0,0.0,0.5},{0.0,0.5,0.0}},
							{{-1.0/3.0,0.0,0.0},{0.0, 2.0/3.0,0.0},{0.0,0.0,-1.0/3.0}}
						};
   
*/
// basis set for stress normal differences
	double e_S_a[5][3][3]= {
							{{1.0,0.0,0.0},{0.0,-1.0,0.0},{0.0,0.0,0.0}},
							{{0.0,0.5,0.0},{0.5,0.0,0.0},{0.0,0.0,0.0}},
							{{0.0,0.0,0.5},{0.0,0.0,0.0},{0.5,0.0,0.0}},
							{{0.0,0.0,0.0},{0.0,0.0,0.5},{0.0,0.5,0.0}},
							{{0.0,0.0,0.0},{0.0,1.0,0.0},{0.0,0.0,-1.0}}
						};


	double e_g_S[5][3][3]= {
							{{ 2.0/3.0,0.0,0.0},{0.0,-1.0/3.0,0.0},{0.0,0.0,-1.0/3.0}},
							{{0.0,1.0,0.0},{1.0,0.0,0.0},{0.0,0.0,0.0}},
							{{0.0,0.0,1.0},{0.0,0.0,0.0},{1.0,0.0,0.0}},
							{{0.0,0.0,0.0},{0.0,0.0,1.0},{0.0,1.0,0.0}},
							{{1.0/3.0,0.0,0.0},{0.0,1.0/3.0,0.0},{0.0,0.0,-2.0/3.0}}
						};
   

	double e_E_a[5][3][3]= {
							{{ 2.0/3.0,0.0,0.0},{0.0,-1.0/3.0,0.0},{0.0,0.0,-1.0/3.0}},
							{{0.0,1.0,0.0},{1.0,0.0,0.0},{0.0,0.0,0.0}},
							{{0.0,0.0,1.0},{0.0,0.0,0.0},{1.0,0.0,0.0}},
							{{0.0,0.0,0.0},{0.0,0.0,1.0},{0.0,1.0,0.0}},
							{{1.0/3.0,0.0,0.0},{0.0,1.0/3.0,0.0},{0.0,0.0,-2.0/3.0}}
						};


	double e_g_E[5][3][3]= {
							{{1.0,0.0,0.0},{0.0,-1.0,0.0},{0.0,0.0,0.0}},
							{{0.0,0.5,0.0},{0.5,0.0,0.0},{0.0,0.0,0.0}},
							{{0.0,0.0,0.5},{0.0,0.0,0.0},{0.5,0.0,0.0}},
							{{0.0,0.0,0.0},{0.0,0.0,0.5},{0.0,0.5,0.0}},
							{{0.0,0.0,0.0},{0.0,1.0,0.0},{0.0,0.0,-1.0}}
						};
						
						   
   double mu_11N[121*NrParticles*NrParticles] ;  		// grand mobility matrix
   double zeta_11N[121*NrParticles*NrParticles] ;  	// grand resistance matrix
   double rho_11N[121*NrParticles*NrParticles] ;  	// grand resistance matrix
   double xi_11x11[11*11] ; 							// generalized friction matrix
   	
	mtrx3D Friction_Tnsr_tt(0.0,0.0,0.0);
	mtrx3D Friction_Tnsr_tr(0.0,0.0,0.0);
	mtrx3D Friction_Tnsr_rt(0.0,0.0,0.0);
	mtrx3D Friction_Tnsr_rr(0.0,0.0,0.0);
 	mtrx53D Friction_Tnsr_dt 	=	null53D;
	mtrx53D Friction_Tnsr_dr	= 	null53D;
	mtrx35D Friction_Tnsr_td	=	null35D;
	mtrx35D Friction_Tnsr_rd	=	null35D;
	mtrx55D Friction_Tnsr_dd	=	null55D;
   
   for (int i=0; i<121; i++)
		{
			xi_11x11[i] = 0.0; 
		}        

		
 std::ofstream outFile1("data.dat");
 outFile1.precision(17);
 
// important all lengths have been normalized by particle radius as metioned in Page 46, Appendix A - Durlofsky, Louis, John F. Brady, and Georges Bossis. 
				// "Dynamic simulation of hydrodynamically interacting particles." Journal of fluid mechanics 180 (1987): 21-49.
				// for ease of programming. 
		
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
				if(a==b) 
					{	
						e_ab_unit = e_ab*0.0; 
					}	
				else
					{
						if(	e_ab2 < (4.0*bead[a].radius*bead[a].radius)	)	
							{	
							
								cout << e_ab2 << " "<< a << " and "<<b<<" numbered particles are touching or overlapping " << endl;
								cout << bead[a].pos.comp[0] << "\t" << bead[a].pos.comp[1] << "\t" << bead[a].pos.comp[2] << "\t" <<  endl;
								cout << bead[b].pos.comp[0] << "\t" << bead[b].pos.comp[1] << "\t" << bead[b].pos.comp[2] << "\t" <<  endl;

								abort();
							}	
					}

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
										
										m_ijkl[i][j][k][l]	=	((3.0/2.0)*x_m[1][1]*(e_ab_unit.comp[i]*e_ab_unit.comp[j] 					-	(1.0/3.0)*kron_del[i][j])*(e_ab_unit.comp[k]*e_ab_unit.comp[l]	
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
								
								g_ijk[i][j][k]				=	(x_g[1][1]*(e_ab_unit.comp[i]*e_ab_unit.comp[j] 	-	(1.0/3.0)*kron_del[i][j])*e_ab_unit.comp[k]
																+ 		y_g[1][1]*(e_ab_unit.comp[i]*kron_del[j][k]		+ 	e_ab_unit.comp[j]*kron_del[i][k]	-	2.0*e_ab_unit.comp[i]*e_ab_unit.comp[j]*e_ab_unit.comp[k]	)	);
										
								h_ijk[i][j][k]				= 	(y_h[1][1]*(e_ab_unit.comp[i]*ep_jkl_e_l			+	e_ab_unit.comp[j]*ep_ikl_e_l										)	);
							}	// k		

							Mobility_Tnsr_tt.comp[i][j]		=	(x_a[1][1]*e_ab_unit.comp[i]*e_ab_unit.comp[j]	+ 	y_a[1][1]*(kron_del[i][j]	- e_ab_unit.comp[i]*e_ab_unit.comp[j]	)	);
							
							Mobility_Tnsr_rt.comp[i][j]		=	(													y_b[1][1]*ep_ijk_e_k													);
						
							Mobility_Tnsr_rr.comp[i][j]		=	(x_c[1][1]*e_ab_unit.comp[i]*e_ab_unit.comp[j]	+ 	y_c[1][1]*(kron_del[i][j]	- e_ab_unit.comp[i]*e_ab_unit.comp[j]	)	);
		
						}	// j
					}	// i

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
										
										m_ijkl[i][j][k][l]	=	 ((3.0/2.0)*x_m[1][0]*(e_ab_unit.comp[i]*e_ab_unit.comp[j] 			-	(1.0/3.0)*kron_del[i][j])*(e_ab_unit.comp[k]*e_ab_unit.comp[l]	
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
								
								g_ijk[i][j][k]				=	(x_g[1][0]*(e_ab_unit.comp[i]*e_ab_unit.comp[j] 	-	(1.0/3.0)*kron_del[i][j])*e_ab_unit.comp[k]
																+ 		y_g[1][0]*(e_ab_unit.comp[i]*kron_del[j][k]		+ 	e_ab_unit.comp[j]*kron_del[i][k]	-	2.0*e_ab_unit.comp[i]*e_ab_unit.comp[j]*e_ab_unit.comp[k]	)	);
										
								h_ijk[i][j][k]				= 	(y_h[1][0]*(e_ab_unit.comp[i]*ep_jkl_e_l			+	e_ab_unit.comp[j]*ep_ikl_e_l										)	);
							}	// k		

							Mobility_Tnsr_tt.comp[i][j]		=	(x_a[1][0]*e_ab_unit.comp[i]*e_ab_unit.comp[j]	+ 	y_a[1][0]*(kron_del[i][j]	- e_ab_unit.comp[i]*e_ab_unit.comp[j]	)	);
							
							Mobility_Tnsr_rt.comp[i][j]		=	(													y_b[1][0]*ep_ijk_e_k													);
						
							Mobility_Tnsr_rr.comp[i][j]		=	(x_c[1][0]*e_ab_unit.comp[i]*e_ab_unit.comp[j]	+ 	y_c[1][0]*(kron_del[i][j]	- e_ab_unit.comp[i]*e_ab_unit.comp[j]	)	);
									
							Mobility_Tnsr_tr	= 	    Mobility_Tnsr_rt*(1.0);		
	
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

						Mobility_Tnsr_dt.comp[p][g]		+=		e_E_a[p][a][b]*g_ijk[a][b][g];	
						Mobility_Tnsr_dr.comp[p][g]		+=		e_E_a[p][a][b]*h_ijk[a][b][g];		
						Mobility_Tnsr_td.comp[g][p]		+=		g_ijk[a][b][g]*e_g_S[p][a][b];	// going from g_dt matrix to ~g_td matrix hence circulation of indices 
						Mobility_Tnsr_rd.comp[g][p]		+=		h_ijk[a][b][g]*e_g_S[p][a][b];
								
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
								Mobility_Tnsr_dd.comp[p][s]		+=		e_E_a[p][a][b]*m_ijkl[a][b][g][d]*e_g_S[s][g][d];			
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
							zeta_11N[k	+	11*NrParticles*l	+	3*b	+	33*NrParticles*a	+	33*NrParticles*NrParticles						] 	=	-Mobility_Tnsr_rt.comp[k][l];	// because Mobility_Tnsr_tr(a,b) = - Mobility_Tnsr_tr(b,a);
							zeta_11N[k	+	11*NrParticles*l	+	3*b	+	33*NrParticles*a	+	3*NrParticles									] 	=	-Mobility_Tnsr_rt.comp[k][l];	// because Mobility_Tnsr_rt(a,b) = - Mobility_Tnsr_rt(b,a);
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
						zeta_11N[k	+	11*NrParticles*l	+	3*a	+	55*NrParticles*b	+	66*NrParticles*NrParticles						] 	=	 -Mobility_Tnsr_td.comp[k][l];		// because Mobility_Tnsr_td(a,b) = - Mobility_Tnsr_td(b,a);
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
               
	inverse ( zeta_11N ,11*NrParticles )	 ; 	
	
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
									Delj.comp[a][k]	+=	e_g_E[k][b][a]	*	bead[j].pos.comp[b]	;
									Deli.comp[k][a]	+=	e_g_E[k][a][b]	*	bead[i].pos.comp[b]	;
							/*		for (int c=0; c<3; c++)
									{
										Unit_tnsr_Redc.comp[k][a]	+=	e_l[k][b][c]		* (b==c)	*	bead[i].pos.comp[a]	;	
									}
							*/	}
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
/*
 // test the 3-indice form of the g,h (td, rd) mobility matrices

	
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
							h_clst_ijk[g][a][b]	+=		e_E_a[p][a][b]*Friction_Tnsr_rd.comp[g][p];		

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
							g_clst_ijk[g][a][b]	+=		e_E_a[p][a][b]*Friction_Tnsr_td.comp[g][p];		

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
    	
		cout<<NrParticles<<std::endl ;
		cout<<xi_11x11[0]<<'\t'<<xi_11x11[6]<<'\t'<<xi_11x11[12]<<'\t'<<xi_11x11[18]<<'\t'<<xi_11x11[24]<<'\t'<<xi_11x11[30]<<std::endl ;
		cout<<xi_11x11[1]<<'\t'<<xi_11x11[7]<<'\t'<<xi_11x11[13]<<'\t'<<xi_11x11[19]<<'\t'<<xi_11x11[25]<<'\t'<<xi_11x11[31]<<std::endl ;
		cout<<xi_11x11[2]<<'\t'<<xi_11x11[8]<<'\t'<<xi_11x11[14]<<'\t'<<xi_11x11[20]<<'\t'<<xi_11x11[26]<<'\t'<<xi_11x11[32]<<std::endl ;
		cout<<std::endl ;
		cout<<xi_11x11[3]<<'\t'<<xi_11x11[9]<<'\t'<<xi_11x11[15]<<'\t'<<xi_11x11[21]<<'\t'<<xi_11x11[27]<<'\t'<<xi_11x11[33]<<std::endl ;
		cout<<xi_11x11[4]<<'\t'<<xi_11x11[10]<<'\t'<<xi_11x11[16]<<'\t'<<xi_11x11[22]<<'\t'<<xi_11x11[28]<<'\t'<<xi_11x11[34]<<std::endl ;
		cout<<xi_11x11[5]<<'\t'<<xi_11x11[11]<<'\t'<<xi_11x11[17]<<'\t'<<xi_11x11[23]<<'\t'<<xi_11x11[29]<<'\t'<<xi_11x11[35]<<std::endl ;

		outFile1<<std::endl ;
		outFile1<<xi_11x11[0]<<'\t'<<xi_11x11[6]<<'\t'<<xi_11x11[12]<<'\t'<<xi_11x11[18]<<'\t'<<xi_11x11[24]<<'\t'<<xi_11x11[30]<<std::endl ;
		outFile1<<xi_11x11[1]<<'\t'<<xi_11x11[7]<<'\t'<<xi_11x11[13]<<'\t'<<xi_11x11[19]<<'\t'<<xi_11x11[25]<<'\t'<<xi_11x11[31]<<std::endl ;
		outFile1<<xi_11x11[2]<<'\t'<<xi_11x11[8]<<'\t'<<xi_11x11[14]<<'\t'<<xi_11x11[20]<<'\t'<<xi_11x11[26]<<'\t'<<xi_11x11[32]<<std::endl ;
		outFile1<<std::endl ;
		outFile1<<xi_11x11[3]<<'\t'<<xi_11x11[9]<<'\t'<<xi_11x11[15]<<'\t'<<xi_11x11[21]<<'\t'<<xi_11x11[27]<<'\t'<<xi_11x11[33]<<std::endl ;
		outFile1<<xi_11x11[4]<<'\t'<<xi_11x11[10]<<'\t'<<xi_11x11[16]<<'\t'<<xi_11x11[22]<<'\t'<<xi_11x11[28]<<'\t'<<xi_11x11[34]<<std::endl ;
		outFile1<<xi_11x11[5]<<'\t'<<xi_11x11[11]<<'\t'<<xi_11x11[17]<<'\t'<<xi_11x11[23]<<'\t'<<xi_11x11[29]<<'\t'<<xi_11x11[35]<<std::endl ;

*/

/*
// Ellipsoid mobilities from Kim and Karrila book ; Page 64

double c = 1.0;	// short axis

double a = 5.0*c;	// long axis

double a_bead = 5.0;	// long axis of bead ellipsoid

double a_bead3 = a_bead*a_bead*a_bead;	

double e = ( pow(( a*a - c*c ),0.5 )) / a ; 	// eccentricity

double e2 = e*e ;

double e3 = e*e*e ;

double e5 = e2*e3 ;

double L = log((1.0+e)/(1.0-e)); 

double XA = (8.0/3.0)*(e3)/(-2.0*e + ( 1.0+e2 )*L) ;

cout << "XA = " << XA << endl;
cout << "aspect ratio = " << a << endl;

double YA = (16.0/3.0)*(e3)/(2.0*e + ( 3.0*(e2) - 1.0 )*L) ;

double YB = 0.0 ; 

double XC = (4.0/3.0)*(e3)*(1.0-e2)/(2.0*e - ( 1.0-e2 )*L) ;

double YC = (4.0/3.0)*(e3)*(2.0-e2)/(-2.0*e + ( 1.0+e2 )*L) ;

double XG = 0.0 ;

double YG = 0.0 ;

double YH = (4.0/3.0)*(e5)/(-2.0*e + ( 1.0+e2 )*L) ;

double XM =  (8.0/15.0)*(e5)/(( 3.0-e2 )*L -6.0*e) ;

double YM_div1= 2.0*e*(2.0*(e2)-3.0) + 3.0*(1.0-e2)*L ; 

double YM_div2 = -2.0*e + (1.0+e2)*L ;

double YM = (4.0/5.0)*(e5)*(2.0*e*(1.0-2.0*(e2)) -(1.0-e2)*L)/(YM_div1*YM_div2) ;

double ZM = (16.0/5.0)*(e5)*(1.0-e2)/( 3.0*((1.0-e2)*(1.0-e2))*L - 2.0*e*(3.0-5.0*(e2)) ) ; 

mtrx3D Friction_Tnsr_tt_anl;
mtrx3D Friction_Tnsr_tr_anl;
mtrx3D Friction_Tnsr_rt_anl;
mtrx3D Friction_Tnsr_rr_anl;
double G_IJK_anl[3][3][3];
double H_IJK_anl[3][3][3];
double M_IJKL_anl[3][3][3][3];

mtrx53D Friction_Tnsr_dt_anl 	=	null53D;
mtrx53D Friction_Tnsr_dr_anl	= 	null53D;
mtrx35D Friction_Tnsr_td_anl	=	null35D;
mtrx35D Friction_Tnsr_rd_anl	=	null35D;
mtrx55D Friction_Tnsr_dd_anl	=	null55D;

 e_ab_unit = {0.0,0.0,1.0}; 		// symmetry axis of ellipsoid; here we take it to be along z-axis
 
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
										
										M_IJKL_anl[i][j][k][l]	=	 ((3.0/2.0)*XM*(e_ab_unit.comp[i]*e_ab_unit.comp[j] 			-	(1.0/3.0)*kron_del[i][j])*(e_ab_unit.comp[k]*e_ab_unit.comp[l]	
																-(1.0/3.0)*kron_del[k][l])
																+(1.0/2.0)*YM*(e_ab_unit.comp[i]*kron_del[j][l]*e_ab_unit.comp[k]	+	e_ab_unit.comp[j]*kron_del[i][l]*e_ab_unit.comp[k]
																					+ e_ab_unit.comp[i]*kron_del[j][k]*e_ab_unit.comp[l]	+ 	e_ab_unit.comp[j]*kron_del[i][k]*e_ab_unit.comp[l]
																					- 4.0*e_ab_unit.comp[i]*e_ab_unit.comp[j]*e_ab_unit.comp[k]*e_ab_unit.comp[l]	)
																			
																+(1.0/2.0)*ZM*(kron_del[i][k]*kron_del[j][l]		+ 	kron_del[j][k]*kron_del[i][l]	- 	kron_del[i][j]*kron_del[k][l]
																+ e_ab_unit.comp[i]*e_ab_unit.comp[j]*kron_del[k][l]	+	kron_del[i][j]*e_ab_unit.comp[k]*e_ab_unit.comp[l]	
																+ e_ab_unit.comp[i]*e_ab_unit.comp[j]*e_ab_unit.comp[k]*e_ab_unit.comp[l]
																- e_ab_unit.comp[i]*kron_del[j][l]*e_ab_unit.comp[k]	- 	e_ab_unit.comp[j]*kron_del[i][l]*e_ab_unit.comp[k]
																- e_ab_unit.comp[i]*kron_del[j][k]*e_ab_unit.comp[l]	- 	e_ab_unit.comp[j]*kron_del[i][k]*e_ab_unit.comp[l]
																));
																																							
									}	// l
									
								ep_ijk_e_k					+=	Levi_Civi[i][j][k]*e_ab_unit.comp[k];
								
								G_IJK_anl[i][j][k]	=	(XG*(e_ab_unit.comp[i]*e_ab_unit.comp[j] 	-	(1.0/3.0)*kron_del[i][j])*e_ab_unit.comp[k]
																+ 		YG*(e_ab_unit.comp[i]*kron_del[j][k]		+ 	e_ab_unit.comp[j]*kron_del[i][k]	-	2.0*e_ab_unit.comp[i]*e_ab_unit.comp[j]*e_ab_unit.comp[k]	)	);
										
								H_IJK_anl[i][j][k]	= 	(YH*(e_ab_unit.comp[i]*ep_jkl_e_l			+	e_ab_unit.comp[j]*ep_ikl_e_l										)	);
							}	// k		

							Friction_Tnsr_tt_anl.comp[i][j]		=	(XA*e_ab_unit.comp[i]*e_ab_unit.comp[j]	+ 	YA*(kron_del[i][j]	- e_ab_unit.comp[i]*e_ab_unit.comp[j]	)	);
							
							Friction_Tnsr_rt_anl.comp[i][j]		=	(													YB*ep_ijk_e_k													);
						
							Friction_Tnsr_rr_anl.comp[i][j]		=	(XC*e_ab_unit.comp[i]*e_ab_unit.comp[j]	+ 	YC*(kron_del[i][j]	- e_ab_unit.comp[i]*e_ab_unit.comp[j]	)	);
									
							Friction_Tnsr_tr_anl	= 	    Friction_Tnsr_rt_anl*(1.0);		
	
						}	// j
					}	// i	


	for (int p=0; p<5; p++)
		{
		for (int g=0; g<3; g++)
			{
			for (int a=0; a<3; a++)
				{
				for (int b=0; b<3; b++)
					{

						Friction_Tnsr_dt_anl.comp[p][g]		+=		e_S_a[p][a][b]*G_IJK_anl[a][b][g];	
						Friction_Tnsr_dr_anl.comp[p][g]		+=		e_S_a[p][a][b]*H_IJK_anl[a][b][g]*a_bead3*(4.0/6.0) ;	// corection for ellipsoid size	
						Friction_Tnsr_td_anl.comp[g][p]		+=		G_IJK_anl[a][b][g]*e_g_E[p][a][b];	// going from g_dt matrix to ~g_td matrix hence circulation of indices 
						Friction_Tnsr_rd_anl.comp[g][p]		+=		H_IJK_anl[a][b][g]*e_g_E[p][a][b]*a_bead3*(4.0/6.0) ;	// corection for ellipsoid size
								
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
								Friction_Tnsr_dd_anl.comp[p][s]		+=		e_S_a[p][a][b]*M_IJKL_anl[a][b][g][d]*e_g_E[s][g][d]*a_bead3*(20.0/18.0);	// corection for ellipsoid size		
							}
						}													
					}
				}
			}
		}			

	cout << "Friction matrix bead-calc and analytical echos" << endl;
	Friction_Tnsr_tt.echo();
	Friction_Tnsr_tt_anl.echo();
	Friction_Tnsr_tr.echo();
	Friction_Tnsr_tr_anl.echo();
	Friction_Tnsr_rt.echo();
	Friction_Tnsr_rt_anl.echo();
	Friction_Tnsr_rr.echo();
	Friction_Tnsr_rr_anl.echo();

	Friction_Tnsr_tt = Friction_Tnsr_tt_anl*a_bead ;		// corection for ellipsoid size
	Friction_Tnsr_rt = Friction_Tnsr_rt_anl ;
	Friction_Tnsr_tr = Friction_Tnsr_tr_anl ;
	Friction_Tnsr_rr = Friction_Tnsr_rr_anl*a_bead3*(8.0/6.0) ;	// corection for ellipsoid size
	Friction_Tnsr_dt = Friction_Tnsr_dt_anl ;
	Friction_Tnsr_dr = Friction_Tnsr_dr_anl ;
	Friction_Tnsr_td = Friction_Tnsr_td_anl ;
	Friction_Tnsr_rd = Friction_Tnsr_rd_anl ;
	Friction_Tnsr_dd = Friction_Tnsr_dd_anl ;

	

	for (int k=0; k<5; k++)
		{				
			// column major format
			cout << Friction_Tnsr_dt.comp[k][0]<< '\t'<< Friction_Tnsr_dt.comp[k][1]<< '\t'<< Friction_Tnsr_dt.comp[k][2]<< endl;
		}

	for (int k=0; k<5; k++)
		{				
			// column major format
				cout << Friction_Tnsr_dt_anl.comp[k][0]<< '\t'<< Friction_Tnsr_dt_anl.comp[k][1]<< '\t'<< Friction_Tnsr_dt_anl.comp[k][2]<< endl;
		}

	for (int k=0; k<5; k++)
		{				
			// column major format
			cout << Friction_Tnsr_dr.comp[k][0]<< '\t'<< Friction_Tnsr_dr.comp[k][1]<< '\t'<< Friction_Tnsr_dr.comp[k][2]<< endl;
		}

	for (int k=0; k<5; k++)
		{				
			// column major format
				cout << Friction_Tnsr_dr_anl.comp[k][0]<< '\t'<< Friction_Tnsr_dr_anl.comp[k][1]<< '\t'<< Friction_Tnsr_dr_anl.comp[k][2]<< endl;
		}		

	for (int k=0; k<5; k++)
		{				
			// column major format
			cout << Friction_Tnsr_rd.comp[0][k]<< '\t'<< Friction_Tnsr_rd.comp[1][k]<< '\t'<< Friction_Tnsr_rd.comp[2][k]<< endl;
		}

	for (int k=0; k<5; k++)
		{				
			// column major format
				cout << Friction_Tnsr_rd_anl.comp[0][k]<< '\t'<< Friction_Tnsr_rd_anl.comp[1][k]<< '\t'<< Friction_Tnsr_rd_anl.comp[2][k]<< endl;
		}		
	
	for (int k=0; k<5; k++)
		{				
		for (int k=0; k<5; k++)
			{			// column major format
			cout << Friction_Tnsr_dd.comp[k][0]<< '\t' ;
			}
			cout<< endl;
		}

	for (int k=0; k<5; k++)
		{				
		for (int k=0; k<5; k++)
			{			// column major format
			cout << Friction_Tnsr_dd_anl.comp[k][0]<< '\t' ; 
			}
			cout<< endl;
		}
		
		
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
*/	
			inverse ( xi_11x11 , 6 )	 ; 			

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
		
		outFile1<<std::endl ;
		outFile1<<mu_d[0][0]<<'\t'<<mu_d[0][1]<<'\t'<<mu_d[0][2]<<'\t'<<mu_d[0][3]<<'\t'<<mu_d[0][4]<<std::endl ;
		outFile1<<mu_d[1][0]<<'\t'<<mu_d[1][1]<<'\t'<<mu_d[1][2]<<'\t'<<mu_d[1][3]<<'\t'<<mu_d[1][4]<<std::endl ;
		outFile1<<mu_d[2][0]<<'\t'<<mu_d[2][1]<<'\t'<<mu_d[2][2]<<'\t'<<mu_d[2][3]<<'\t'<<mu_d[2][4]<<std::endl ;
		outFile1<<std::endl ;
		outFile1<<mu_d[3][0]<<'\t'<<mu_d[3][1]<<'\t'<<mu_d[3][2]<<'\t'<<mu_d[3][3]<<'\t'<<mu_d[3][4]<<std::endl ;
		outFile1<<mu_d[4][0]<<'\t'<<mu_d[4][1]<<'\t'<<mu_d[4][2]<<'\t'<<mu_d[4][3]<<'\t'<<mu_d[4][4]<<std::endl ;
		outFile1<<mu_d[5][0]<<'\t'<<mu_d[5][1]<<'\t'<<mu_d[5][2]<<'\t'<<mu_d[5][3]<<'\t'<<mu_d[5][4]<<std::endl ;		
		outFile1<<std::endl ;
		outFile1<<mu_dd[0][0]<<'\t'<<mu_dd[0][1]<<'\t'<<mu_dd[0][2]<<'\t'<<mu_dd[0][3]<<'\t'<<mu_dd[0][4]<<std::endl ;
		outFile1<<mu_dd[1][0]<<'\t'<<mu_dd[1][1]<<'\t'<<mu_dd[1][2]<<'\t'<<mu_dd[1][3]<<'\t'<<mu_dd[1][4]<<std::endl ;
		outFile1<<mu_dd[2][0]<<'\t'<<mu_dd[2][1]<<'\t'<<mu_dd[2][2]<<'\t'<<mu_dd[2][3]<<'\t'<<mu_dd[2][4]<<std::endl ;
		outFile1<<mu_dd[3][0]<<'\t'<<mu_dd[3][1]<<'\t'<<mu_dd[3][2]<<'\t'<<mu_dd[3][3]<<'\t'<<mu_dd[3][4]<<std::endl ;
		outFile1<<mu_dd[4][0]<<'\t'<<mu_dd[4][1]<<'\t'<<mu_dd[4][2]<<'\t'<<mu_dd[4][3]<<'\t'<<mu_dd[4][4]<<std::endl ;
		
// center of diffusion calculation based on "Hydrodynamic properties of rigid particles: comparison of different modeling and computational procedures." Biophysical journal 76.6 (1999): 3044-3057.			
// Page 3046 equation 13.
vctr3D ctr_diff ; 
double temp_mat[3*3];
temp_mat[0] =  xi_11x11[28] + xi_11x11[35]	;
temp_mat[1] = -xi_11x11[27] 				; 
temp_mat[2] = -xi_11x11[33]				; 
temp_mat[3] = -xi_11x11[27]				;
temp_mat[4] =  xi_11x11[21] + xi_11x11[35]	;
temp_mat[5] = -xi_11x11[34]				;
temp_mat[6] = -xi_11x11[33]				;
temp_mat[7] = -xi_11x11[34]				;
temp_mat[8] =  xi_11x11[28] + xi_11x11[21]	;

inverse ( temp_mat , 3 )	 ; 	

ctr_diff.comp[0] = temp_mat[0]*(xi_11x11[16] - xi_11x11[11]) + temp_mat[3]*(xi_11x11[5] - xi_11x11[15]) + temp_mat[6]*(xi_11x11[9] - xi_11x11[4]) ;
ctr_diff.comp[1] = temp_mat[1]*(xi_11x11[16] - xi_11x11[11]) + temp_mat[4]*(xi_11x11[5] - xi_11x11[15]) + temp_mat[7]*(xi_11x11[9] - xi_11x11[4]) ;
ctr_diff.comp[2] = temp_mat[2]*(xi_11x11[16] - xi_11x11[11]) + temp_mat[5]*(xi_11x11[5] - xi_11x11[15]) + temp_mat[8]*(xi_11x11[9] - xi_11x11[4]) ;

// ctr_diff.echo();
 cout.precision(17);
 cout<<"center of diffusion"<<'\t'<<ctr_diff.comp[0]<<'\t'<<ctr_diff.comp[1]<<'\t'<<ctr_diff.comp[2]<<std::endl ;

mtrx3D D_tt, D_tr, D_rt , D_rr, U_OD;
mtrx35D D_td;

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
				Delj.comp[a][k]	+=	e_g_E[k][b][a]	*	ctr_diff.comp[b]	;
				Deli.comp[k][a]	+=	e_g_E[k][a][b]	*	ctr_diff.comp[b]	;
		/*		for (int c=0; c<3; c++)
				{
					Unit_tnsr_Redc.comp[k][a]	+=	e_l[k][b][c]		* (b==c)	*	bead[i].pos.comp[a]	;	
				}
		*/	}
		}		
	}
 
D_tt.comp[0][0] = xi_11x11[0];
D_tt.comp[1][0] = xi_11x11[1];
D_tt.comp[2][0] = xi_11x11[2];
D_tt.comp[0][1] = xi_11x11[6];
D_tt.comp[1][1] = xi_11x11[7];
D_tt.comp[2][1] = xi_11x11[8];
D_tt.comp[0][2] = xi_11x11[12];
D_tt.comp[1][2] = xi_11x11[13];
D_tt.comp[2][2] = xi_11x11[14];

D_tr.comp[0][0] = xi_11x11[3];
D_tr.comp[1][0] = xi_11x11[4];
D_tr.comp[2][0] = xi_11x11[5];
D_tr.comp[0][1] = xi_11x11[9];
D_tr.comp[1][1] = xi_11x11[10];
D_tr.comp[2][1] = xi_11x11[11];
D_tr.comp[0][2] = xi_11x11[15];
D_tr.comp[1][2] = xi_11x11[16];
D_tr.comp[2][2] = xi_11x11[17];

D_rt.comp[0][0] = xi_11x11[18];
D_rt.comp[1][0] = xi_11x11[19];
D_rt.comp[2][0] = xi_11x11[20];
D_rt.comp[0][1] = xi_11x11[24];
D_rt.comp[1][1] = xi_11x11[25];
D_rt.comp[2][1] = xi_11x11[26];
D_rt.comp[0][2] = xi_11x11[30];
D_rt.comp[1][2] = xi_11x11[31];
D_rt.comp[2][2] = xi_11x11[32];

D_rr.comp[0][0] = xi_11x11[21];
D_rr.comp[1][0] = xi_11x11[22];
D_rr.comp[2][0] = xi_11x11[23];
D_rr.comp[0][1] = xi_11x11[27];
D_rr.comp[1][1] = xi_11x11[28];
D_rr.comp[2][1] = xi_11x11[29];
D_rr.comp[0][2] = xi_11x11[33];
D_rr.comp[1][2] = xi_11x11[34];
D_rr.comp[2][2] = xi_11x11[35];

U_OD.comp[0][0] =  0.0;
U_OD.comp[1][0] =  ctr_diff.comp[2];
U_OD.comp[2][0] = -ctr_diff.comp[1];
U_OD.comp[0][1] = -ctr_diff.comp[2];
U_OD.comp[1][1] =  0.0;
U_OD.comp[2][1] =  ctr_diff.comp[0];
U_OD.comp[0][2] =  ctr_diff.comp[1];
U_OD.comp[1][2] = -ctr_diff.comp[0];
U_OD.comp[2][2] =  0.0;

mtrx3D D_tt_CoD = D_tt -  U_OD*D_rr*U_OD + D_rt*U_OD - U_OD*D_tr ; 
mtrx3D D_tr_CoD = D_tr +  D_rr*U_OD ;  // based on equations 42 from Wouter's notes "clusterdyn"
mtrx3D D_rt_CoD = D_rt -  U_OD*D_rr ;  // based on equations 43 from Wouter's notes "clusterdyn"

xi_11x11[0]  = D_tt_CoD.comp[0][0];
xi_11x11[1]  = D_tt_CoD.comp[1][0];
xi_11x11[2]  = D_tt_CoD.comp[2][0];
xi_11x11[6]  = D_tt_CoD.comp[0][1];
xi_11x11[7]  = D_tt_CoD.comp[1][1];
xi_11x11[8]  = D_tt_CoD.comp[2][1];
xi_11x11[12] = D_tt_CoD.comp[0][2];
xi_11x11[13] = D_tt_CoD.comp[1][2];
xi_11x11[14] = D_tt_CoD.comp[2][2];

xi_11x11[3]  = D_tr_CoD.comp[0][0];
xi_11x11[4]  = D_tr_CoD.comp[1][0];
xi_11x11[5]  = D_tr_CoD.comp[2][0];
xi_11x11[9]  = D_tr_CoD.comp[0][1];
xi_11x11[10] = D_tr_CoD.comp[1][1];
xi_11x11[11] = D_tr_CoD.comp[2][1];
xi_11x11[15] = D_tr_CoD.comp[0][2];
xi_11x11[16] = D_tr_CoD.comp[1][2];
xi_11x11[17] = D_tr_CoD.comp[2][2];

xi_11x11[18]  = D_rt_CoD.comp[0][0];
xi_11x11[19]  = D_rt_CoD.comp[1][0];
xi_11x11[20]  = D_rt_CoD.comp[2][0];
xi_11x11[24]  = D_rt_CoD.comp[0][1];
xi_11x11[25]  = D_rt_CoD.comp[1][1];
xi_11x11[26]  = D_rt_CoD.comp[2][1];
xi_11x11[30]  = D_rt_CoD.comp[0][2];
xi_11x11[31]  = D_rt_CoD.comp[1][2];
xi_11x11[32]  = D_rt_CoD.comp[2][2];
	
		outFile1<<std::endl ;
		outFile1<<xi_11x11[0]<<'\t'<<xi_11x11[6]<<'\t'<<xi_11x11[12]<<'\t'<<xi_11x11[18]<<'\t'<<xi_11x11[24]<<'\t'<<xi_11x11[30]<<std::endl ;
		outFile1<<xi_11x11[1]<<'\t'<<xi_11x11[7]<<'\t'<<xi_11x11[13]<<'\t'<<xi_11x11[19]<<'\t'<<xi_11x11[25]<<'\t'<<xi_11x11[31]<<std::endl ;
		outFile1<<xi_11x11[2]<<'\t'<<xi_11x11[8]<<'\t'<<xi_11x11[14]<<'\t'<<xi_11x11[20]<<'\t'<<xi_11x11[26]<<'\t'<<xi_11x11[32]<<std::endl ;
		outFile1<<std::endl ;
		outFile1<<xi_11x11[3]<<'\t'<<xi_11x11[9]<<'\t'<<xi_11x11[15]<<'\t'<<xi_11x11[21]<<'\t'<<xi_11x11[27]<<'\t'<<xi_11x11[33]<<std::endl ;
		outFile1<<xi_11x11[4]<<'\t'<<xi_11x11[10]<<'\t'<<xi_11x11[16]<<'\t'<<xi_11x11[22]<<'\t'<<xi_11x11[28]<<'\t'<<xi_11x11[34]<<std::endl ;
		outFile1<<xi_11x11[5]<<'\t'<<xi_11x11[11]<<'\t'<<xi_11x11[17]<<'\t'<<xi_11x11[23]<<'\t'<<xi_11x11[29]<<'\t'<<xi_11x11[35]<<std::endl ;
		

mu_d[0][0] = mu_d[0][0] - ( - Delj.comp[0][0] +  U_OD.comp[0][0]*mu_d[3][0] + U_OD.comp[0][1]*mu_d[4][0] + U_OD.comp[0][2]*mu_d[5][0] ) ;
mu_d[0][1] = mu_d[0][1] - ( - Delj.comp[0][1] +  U_OD.comp[0][0]*mu_d[3][1] + U_OD.comp[0][1]*mu_d[4][1] + U_OD.comp[0][2]*mu_d[5][1] ) ;
mu_d[0][2] = mu_d[0][2] - ( - Delj.comp[0][2] +  U_OD.comp[0][0]*mu_d[3][2] + U_OD.comp[0][1]*mu_d[4][2] + U_OD.comp[0][2]*mu_d[5][2] ) ;
mu_d[0][3] = mu_d[0][3] - ( - Delj.comp[0][3] +  U_OD.comp[0][0]*mu_d[3][3] + U_OD.comp[0][1]*mu_d[4][3] + U_OD.comp[0][2]*mu_d[5][3] ) ;
mu_d[0][4] = mu_d[0][4] - ( - Delj.comp[0][4] +  U_OD.comp[0][0]*mu_d[3][4] + U_OD.comp[0][1]*mu_d[4][4] + U_OD.comp[0][2]*mu_d[5][4] ) ;
mu_d[1][0] = mu_d[1][0] - ( - Delj.comp[1][0] +  U_OD.comp[1][0]*mu_d[3][0] + U_OD.comp[1][1]*mu_d[4][0] + U_OD.comp[1][2]*mu_d[5][0] ) ;
mu_d[1][1] = mu_d[1][1] - ( - Delj.comp[1][1] +  U_OD.comp[1][0]*mu_d[3][1] + U_OD.comp[1][1]*mu_d[4][1] + U_OD.comp[1][2]*mu_d[5][1] ) ;
mu_d[1][2] = mu_d[1][2] - ( - Delj.comp[1][2] +  U_OD.comp[1][0]*mu_d[3][2] + U_OD.comp[1][1]*mu_d[4][2] + U_OD.comp[1][2]*mu_d[5][2] ) ;
mu_d[1][3] = mu_d[1][3] - ( - Delj.comp[1][3] +  U_OD.comp[1][0]*mu_d[3][3] + U_OD.comp[1][1]*mu_d[4][3] + U_OD.comp[1][2]*mu_d[5][3] ) ;
mu_d[1][4] = mu_d[1][4] - ( - Delj.comp[1][4] +  U_OD.comp[1][0]*mu_d[3][4] + U_OD.comp[1][1]*mu_d[4][4] + U_OD.comp[1][2]*mu_d[5][4] ) ;
mu_d[2][0] = mu_d[2][0] - ( - Delj.comp[2][0] +  U_OD.comp[2][0]*mu_d[3][0] + U_OD.comp[2][1]*mu_d[4][0] + U_OD.comp[2][2]*mu_d[5][0] ) ;
mu_d[2][1] = mu_d[2][1] - ( - Delj.comp[2][1] +  U_OD.comp[2][0]*mu_d[3][1] + U_OD.comp[2][1]*mu_d[4][1] + U_OD.comp[2][2]*mu_d[5][1] ) ;
mu_d[2][2] = mu_d[2][2] - ( - Delj.comp[2][2] +  U_OD.comp[2][0]*mu_d[3][2] + U_OD.comp[2][1]*mu_d[4][2] + U_OD.comp[2][2]*mu_d[5][2] ) ;
mu_d[2][3] = mu_d[2][3] - ( - Delj.comp[2][3] +  U_OD.comp[2][0]*mu_d[3][3] + U_OD.comp[2][1]*mu_d[4][3] + U_OD.comp[2][2]*mu_d[5][3] ) ;
mu_d[2][4] = mu_d[2][4] - ( - Delj.comp[2][4] +  U_OD.comp[2][0]*mu_d[3][4] + U_OD.comp[2][1]*mu_d[4][4] + U_OD.comp[2][2]*mu_d[5][4] ) ;    
	
		outFile1<<std::endl ;
		outFile1<<mu_d[0][0]<<'\t'<<mu_d[0][1]<<'\t'<<mu_d[0][2]<<'\t'<<mu_d[0][3]<<'\t'<<mu_d[0][4]<<std::endl ;
		outFile1<<mu_d[1][0]<<'\t'<<mu_d[1][1]<<'\t'<<mu_d[1][2]<<'\t'<<mu_d[1][3]<<'\t'<<mu_d[1][4]<<std::endl ;
		outFile1<<mu_d[2][0]<<'\t'<<mu_d[2][1]<<'\t'<<mu_d[2][2]<<'\t'<<mu_d[2][3]<<'\t'<<mu_d[2][4]<<std::endl ;
		outFile1<<std::endl ;
		outFile1<<mu_d[3][0]<<'\t'<<mu_d[3][1]<<'\t'<<mu_d[3][2]<<'\t'<<mu_d[3][3]<<'\t'<<mu_d[3][4]<<std::endl ;
		outFile1<<mu_d[4][0]<<'\t'<<mu_d[4][1]<<'\t'<<mu_d[4][2]<<'\t'<<mu_d[4][3]<<'\t'<<mu_d[4][4]<<std::endl ;
		outFile1<<mu_d[5][0]<<'\t'<<mu_d[5][1]<<'\t'<<mu_d[5][2]<<'\t'<<mu_d[5][3]<<'\t'<<mu_d[5][4]<<std::endl ;		
		outFile1<<std::endl ;
		outFile1<<mu_dd[0][0]<<'\t'<<mu_dd[0][1]<<'\t'<<mu_dd[0][2]<<'\t'<<mu_dd[0][3]<<'\t'<<mu_dd[0][4]<<std::endl ;
		outFile1<<mu_dd[1][0]<<'\t'<<mu_dd[1][1]<<'\t'<<mu_dd[1][2]<<'\t'<<mu_dd[1][3]<<'\t'<<mu_dd[1][4]<<std::endl ;
		outFile1<<mu_dd[2][0]<<'\t'<<mu_dd[2][1]<<'\t'<<mu_dd[2][2]<<'\t'<<mu_dd[2][3]<<'\t'<<mu_dd[2][4]<<std::endl ;
		outFile1<<mu_dd[3][0]<<'\t'<<mu_dd[3][1]<<'\t'<<mu_dd[3][2]<<'\t'<<mu_dd[3][3]<<'\t'<<mu_dd[3][4]<<std::endl ;
		outFile1<<mu_dd[4][0]<<'\t'<<mu_dd[4][1]<<'\t'<<mu_dd[4][2]<<'\t'<<mu_dd[4][3]<<'\t'<<mu_dd[4][4]<<std::endl ;
std::ofstream ofile("data_binary.bin", ios::out | ios::binary);

 // cout << "temp_mu" << '\t';

double temp_mu;

		for (int l=0; l<36; l++)
				{
					ofile.write((char*) &xi_11x11[l], sizeof(xi_11x11[l]));

				}
				
			for (int l=0; l<6; l++)
				{
				for (int k=0; k<5; k++)
					{
						temp_mu = mu_d[l][k] ; 
						
					//	cout << temp_mu << '\t';

						ofile.write((char*) &temp_mu, sizeof(temp_mu));
					
					}
				}

			for (int l=0; l<5; l++)
				{
				for (int k=0; k<5; k++)
					{
						
						ofile.write((char*) &mu_dd[l][k], sizeof(mu_dd[l][k]));
						// cout << temp_mu << '\t';

					}
				}		
					
			for (int l=0; l<3; l++)
				{						
						ofile.write((char*) &ctr_diff.comp[l], sizeof(ctr_diff.comp[l]));
						// cout << temp_mu << '\t';
				}	
								
			ofile.close();
/*			
		ifstream File ( "data_binary.bin" , ios::in | ios::binary );

		for (int l=0; l<36; l++)
				{
							File.read( (char*) &temp_mu     , sizeof(temp_mu     ) );
							cout << temp_mu << '\n';
				}

			for (int l=0; l<6; l++)
				{
				for (int k=0; k<5; k++)
					{
							File.read( (char*) &temp_mu     , sizeof(temp_mu     ) );
							cout << temp_mu << '\t';
					}
				}


 cout.precision(17);

			for (int l=0; l<5; l++)
				{
				for (int k=0; k<5; k++)
					{
						
							File.read( (char*) &temp_mu     , sizeof(temp_mu     ) );
							cout << temp_mu << '\n';

					}
				}
*/
				 	
/*
 * test the 3-indice form of the g,h matrices
		
		for (int a=0; a<3; a++)
		{
		for (int b=0; b<3; b++)
			{
			for (int g=0; g<3; g++)
				{
					g_clst_ijk[g][a][b] = 0.0;
					
				for (int p=0; p<5; p++)
					{
							g_clst_ijk[g][a][b]	+=		e_S_a[p][a][b]*mu_d[g][p];		

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
					h_clst_ijk[g][a][b] = 0.0;
					
				for (int p=0; p<5; p++)
					{
							h_clst_ijk[g][a][b]	+=		e_g_S[p][a][b]*mu_d[g+3][p];		

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
*/
	return 0 ;
}
