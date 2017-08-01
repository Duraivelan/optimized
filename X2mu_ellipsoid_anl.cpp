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
	
		
// the five base matrices for strain tensor // option 5 :  equation 419 wouter's tex version clusterdyn_110816_1556

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
  


	double	a_norm = 1.0; // /(6.0*M_PI*eta_0*bead[0].radius);						// mobility matrix a non-dimensionalized by 6*pi*mu*r
	double	b_norm = 1.0; // /(6.0*M_PI*eta_0*bead[0].radius*bead[0].radius);						// mobility matrix a non-dimensionalized by 6*pi*mu*r2
	double	c_norm = 1.0; // /(6.0*M_PI*eta_0*bead[0].radius*bead[0].radius*bead[0].radius);						// mobility matrix a non-dimensionalized by 6*pi*mu*r3
	double	g_norm = 1.0; // /(6.0*M_PI*eta_0*bead[0].radius*bead[0].radius*bead[0].radius);						//		assuming correction factor of 6*pi*mu*r3	
	double	h_norm = 1.0; // /(6.0*M_PI*eta_0*bead[0].radius*bead[0].radius*bead[0].radius);						//		assuming correction factor of 6*pi*mu*r3				
	double	m_norm = 1.0; // /(6.0*M_PI*eta_0*bead[0].radius*bead[0].radius*bead[0].radius);						//		assuming correction factor of 6*pi*mu*r3	
				
// Ellipsoid mobilities from Kim and Karrila book ; Page 64

double c = 1.0;	// short axis

double a = 6.0*c;	// long axis

double e = ( pow(( a*a - c*c ),0.5 )) / a ; 	// eccentricity

double e2 = e*e ;

double e3 = e*e*e ;

double e5 = e2*e3 ;

double L = log((1+e)/(1-e)); 

double XA = (8/3)*(e3)/(-2*e + ( 1+e2 )*L) ;

double YA = (16/3)*(e3)/(2*e + ( 3*(e2) - 1 )*L) ;

double YB = 0.0 ; 

double XC = (4/3)*(e3)*(1-e2)/(2*e - ( 1-e2 )*L) ;

double YC = (4/3)*(e3)*(2-e2)/(-2*e + ( 1+e2 )*L) ;

double XG = 0.0 ;

double YG = 0.0 ;

double YH = (4/3)*(e5)/(-2*e + ( 1+e2 )*L) ;

double XM =  (8/15)*(e5)/(( 3-e2 )*L -6*e) ;

double YM_div1= 2*e*(2*(e2)-3) + 3*(1-e2)*L ; 

double YM_div2 = -2*e + (1+e2)*L ;

double YM = (4/5)*(e5)*(2*e*(1-2*(e2)) -(1-e2)*L)/(YM_div1*YM_div2) ;

double ZM = (16/5)*(e5)*(1-e2)/( 3*((1-e2)*(1-e2))*L - 2*e*(3-5*(e2)) ) ; 

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

vctr3D e_ab_unit = null3D; 
	
				 					
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
										
										M_IJKL_anl[i][j][k][l]	=	 m_norm*((3.0/2.0)*XM*(e_ab_unit.comp[i]*e_ab_unit.comp[j] 			-	(1.0/3.0)*kron_del[i][j])*(e_ab_unit.comp[k]*e_ab_unit.comp[l]	
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
								
								G_IJK_anl[i][j][k]	=	g_norm*(XG*(e_ab_unit.comp[i]*e_ab_unit.comp[j] 	-	(1.0/3.0)*kron_del[i][j])*e_ab_unit.comp[k]
																+ 		YG*(e_ab_unit.comp[i]*kron_del[j][k]		+ 	e_ab_unit.comp[j]*kron_del[i][k]	-	2.0*e_ab_unit.comp[i]*e_ab_unit.comp[j]*e_ab_unit.comp[k]	)	);
										
								H_IJK_anl[i][j][k]	= 	h_norm*(YH*(e_ab_unit.comp[i]*ep_jkl_e_l			+	e_ab_unit.comp[j]*ep_ikl_e_l										)	);
							}	// k		

							Friction_Tnsr_tt_anl.comp[i][j]		=	a_norm*(XA*e_ab_unit.comp[i]*e_ab_unit.comp[j]	+ 	YA*(kron_del[i][j]	- e_ab_unit.comp[i]*e_ab_unit.comp[j]	)	);
							
							Friction_Tnsr_rt_anl.comp[i][j]		=	b_norm*(													YB*ep_ijk_e_k													);
						
							Friction_Tnsr_rr_anl.comp[i][j]		=	c_norm*(XC*e_ab_unit.comp[i]*e_ab_unit.comp[j]	+ 	YC*(kron_del[i][j]	- e_ab_unit.comp[i]*e_ab_unit.comp[j]	)	);
									
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
						Friction_Tnsr_dr_anl.comp[p][g]		+=		e_S_a[p][a][b]*H_IJK_anl[a][b][g];		
						Friction_Tnsr_td_anl.comp[g][p]		+=		G_IJK_anl[a][b][g]*e_g_E[p][a][b];	// going from g_dt matrix to ~g_td matrix hence circulation of indices 
						Friction_Tnsr_rd_anl.comp[g][p]		+=		H_IJK_anl[a][b][g]*e_g_E[p][a][b];
								
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
								Friction_Tnsr_dd_anl.comp[p][s]		+=		e_S_a[p][a][b]*M_IJKL_anl[a][b][g][d]*e_g_E[s][g][d];			
							}
						}													
					}
				}
			}
		}			   

return 0 ; 
}
