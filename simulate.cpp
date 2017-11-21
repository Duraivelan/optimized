/*
 * simuate particle dynamics using the 11x11 mobility matrix outputted by X2mu.cpp
 * Input file : 
 * 				** init.dat **
 * 					viscoity 
 * 					radius
 * 					no. of beads 
 * 					kb 
 * 					T
 * 					shear rate
 * 				** XYZ.dat **
 *					X,Y,Z particle positons
 *
 *  Output file :  							
 * 			
 * Current version : assumes all particle radii to be same, also doesn't compute the lubrication forces
 * 					 the equation motion follows wouter's eq. 327, where every term is non-dimentionalized accordingly	
 * 			
 *  */
 
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
//#include</storage3/usr/people/duraivelan/Downloads/eigen-eigen-bdd17ee3b1b3/Eigen/Eigenvalues>
#include<Eigen/Eigenvalues>

using namespace Eigen;

using namespace std;

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


std::random_device seed;
std::mt19937 gen{seed()};
std::normal_distribution<> R1(0.0,1.0),R2(0.0,1.0),R3(0.0,1.0),R4(0.0,1.0),R5(0.0,1.0),R6(0.0,1.0);
/*
double e_g_S[5][3][3]= {
						{{1.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,-1.0}},
						{{0.0,1.0,0.0},{1.0,0.0,0.0},{0.0,0.0,0.0}},
						{{0.0,0.0,1.0},{0.0,0.0,0.0},{1.0,0.0,0.0}},
						{{0.0,0.0,0.0},{0.0,0.0,1.0},{0.0,1.0,0.0}},
						{{0.0,0.0,0.0},{0.0,1.0,0.0},{0.0,0.0,-1.0}}
					};

double e_E_a[5][3][3]= {
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

					
void brownian( long long int step , vector<ParticleData>& cluster, vector<SubData>& particle, int *Max_Cluster_N , double vel_scale , 
	const double force_norm	,
	const double torque_norm,
	const double pos_norm	,
	const double vel_norm 	,
	const double stochas_norm,
	const double shear_rate,
	const double sqrt_2kbTdt,
	const vctr3D box,
	const vctr3D rbox) {
double a, b , c, lambda;
vctr4D quat_old;
						    						
for(int i=0;i<*Max_Cluster_N;i++) 
	{
		vctr3D rand(R1(gen), R2(gen), R3(gen));
		vctr3D rand1(R4(gen), R5(gen), R6(gen));
		
		// simple shear flow;  flow in x-direction, gradient in y-direction, vorticity in z-direction

		vctr3D u_inf(0.0,shear_rate*cluster[i].pos.comp[0],0.0); 		// shear flow gradient in y-direction

		const mtrx3D E_inf(	
							{0.0,shear_rate/2.0,0.0},
							{shear_rate/2.0,0.0,0.0},
							{0.0,0.0,0.0});

		const vctr3D w_inf(0.0,0.0,0.5*shear_rate);
		
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
								E_inf_bt.comp[j]	+=	 e_E_a[j][k][l]*E_inf_b.comp[k][l];
							//	cout.precision(17);
							//	cout<< E_inf_b.comp[k][l]<<'\t';
							}
					}
				//	cout.precision(17);
				//	cout<< E_inf_bt.comp[j]<<'\t';
			}	

		 
//		if (cluster[i].Sub_Length>0) 
//			{
				cluster[i].pos+=cluster[i].rotmat*cluster[i].mobility_tnsr*(~cluster[i].rotmat)*(cluster[i].frc*force_norm*dt) 
								+cluster[i].rotmat*cluster[i].mobility_tnsr_tr*(~cluster[i].rotmat)*(cluster[i].trq*torque_norm*dt)
								+ cluster[i].rotmat*cluster[i].mobility_tnsr_sqrt*(rand*stochas_norm)
								+u_inf*vel_norm*dt-cluster[i].rotmat*(cluster[i].mobility_tnsr_td*E_inf_bt)*dt ;
				cluster[i].pos = cluster[i].pos*pos_norm ; 
				for(int m=0;m<5;m++) 
					{
						cluster[i].Stresslet.comp[m]=0.0;			
					for(int n=0;n<5;n++) 
						{	
						cluster[i].Stresslet.comp[m]-=	cluster[i].mobility_tnsr_dd.comp[m][n]*E_inf_bt.comp[n]	;
						cout.precision(17);
					//	cout<< cluster[i].mobility_tnsr_dd.comp[m][n]<<'\t';
					//	cout<< cluster[i].Stresslet.comp[m]<<'\n';
						}	
					//	cout.precision(17);
					//	cout<< cluster[i].Stresslet.comp[m]<<'\t';
					}

				
				for(int k=0;k<3;k++) 
					{
						for(int l=0;l<3;l++) 
						{
								S_b.comp[k][l] = 0.0;
								
							for(int j=0;j<5;j++) 
							{						
								S_b.comp[k][l]	+=	 e_g_S[j][k][l]*cluster[i].Stresslet.comp[j];
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
											cluster[i].Stresslet.comp[m]	+=	 e_S_a[m][k][l]*S_s.comp[k][l];
										}
								}
						}										
	
	//			if(xx_rotation)	
	//			{
			// update Q
				quat_old=cluster[i].quat;
				// translate space-fixed w*dt (i.e. theta) (3 dimensions) into qdot (4 dimensions).
				// based on the Wotuer's paper on An elementary singularity-free Rotational Brownian Dynamics algorithm for anisotropic particles 
				// J. Chem. Phys. 142, 114103 (2015)
				
				cluster[i].theta   	= 	cluster[i].rot_mobility_tnsr*(~cluster[i].rotmat)*(cluster[i].trq*torque_norm*dt)
									//	+cluster[i].rot_mobility_tnsr_rt*(~cluster[i].rotmat)*(cluster[i].frc*force_norm*dt)
										+  cluster[i].rot_mobility_tnsr_sqrt*(rand1*stochas_norm)
										-  (cluster[i].mobility_tnsr_rd*E_inf_bt)*dt; 	// body fixed omega
				cluster[i].omega	=	w_inf*dt;						// space-fixed omega
				cluster[i].quat		= cluster[i].theta2quat() + cluster[i].omega2qdot() ;
			//	cout<<cluster[i].theta.comp[0]<<cluster[i].theta.comp[1]<<cluster[i].theta.comp[2]<<endl;
			// lagragian normalization of quaternions; see your notes;
			// after quaternion update you get new quaternion (say ~q) which non-normalised, i.e. |~q|!=1; 
			// assuming qi(t+dt) = ~qi + lambda*qi(t);
			// hence 	|qi(t+dt)| = |~qi + lambda*qi(t)| =1;
				a=quat_old.norm2();
				b=cluster[i].quat*quat_old*2.0;
				c=cluster[i].quat.norm2()-1.0;
				lambda = (-b+sqrt(b*b-4.0*a*c))/(2.0*a);
				cluster[i].quat=cluster[i].quat+quat_old*lambda;
		
			// update A matrix
	//			}
				cluster[i].quat2rotmat();
	
				for (int j=0; j<cluster[i].Sub_Length; j++) 
					{
						particle[cluster[i].sub[j]].pos = cluster[i].pos + cluster[i].rotmat*particle[cluster[i].sub[j]].pos_bdyfxd;
						particle[cluster[i].sub[j]].pos.PBC(box,rbox);
					}
							
		
			//	cluster[i].pos.PBC(box,rbox);
//			} 
/*			else 
			{
				cluster[i].radii	=	0.56;//rmin*0.5 ;		// radii of single particle is sqrt(rmin_x^2+rmin_y^2+rmin_z^2)
				cluster[i].pos+=cluster[i].frc*mu*dt+rand*mu_sqrt*sqrt_2kbTdt;
				cluster[i].pos.PBC(box,rbox);
				for (int j=0; j< cluster[i].Sub_Length; j++) 
					{
						particle[cluster[i].sub[j]].pos=cluster[i].pos;
					}
			}
*/		 
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


double eta_0_temp ; 
double eta_6pi_temp ; 
double bead_radii_temp ; 	// radius of bead used in the mobility calculation i.e. X2mu.cpp file; ideally kept = 1.0 to avoid non-dimensionaliztion errors
int NrParticles_temp;
double kb_temp;
double T0_temp;
double shear_rate_temp;

std::ifstream dataFile;
std::string fileName="init.dat";

//read viscoisty and x,y,z positions from new_cluster.dat
dataFile.open(fileName);
if(!dataFile.good()) {
	std::cerr<<"Given file is corrupt /n"<<std::endl;
}
else {	
   	std::string line0;
	std::getline(dataFile,line0);
   	std::istringstream currentLine0(line0);  
   	currentLine0 >> eta_0_temp;
   	eta_6pi_temp = eta_0_temp*6.0*M_PI ; 
   	cout << eta_0_temp << endl;
	
	currentLine0.str("");
	currentLine0.clear(); // Clear state flags.
	std::getline(dataFile,line0);
	currentLine0.str(line0); 

   	currentLine0 >> bead_radii_temp;
   	cout << bead_radii_temp << endl;

	currentLine0.str("");
	currentLine0.clear(); // Clear state flags.
	std::getline(dataFile,line0);
	currentLine0.str(line0); 

	currentLine0 >> NrParticles_temp;
   	cout << NrParticles_temp << endl;

	currentLine0.str("");
	currentLine0.clear(); // Clear state flags.
	std::getline(dataFile,line0);
	currentLine0.str(line0); 

	currentLine0 >> kb_temp;
   	cout << kb_temp << endl;
   	
	currentLine0.str("");
	currentLine0.clear(); // Clear state flags.
	std::getline(dataFile,line0);
	currentLine0.str(line0); 

	currentLine0 >> T0_temp;
   	cout << T0_temp << endl;
   	
	currentLine0.str("");
	currentLine0.clear(); // Clear state flags.
	std::getline(dataFile,line0);
	currentLine0.str(line0); 

	currentLine0 >> shear_rate_temp;
   	cout << shear_rate_temp << endl;
}

dataFile.close();  
dataFile.clear();

const double eta_0 = eta_0_temp; 
const double eta_6pi = eta_6pi_temp ;
 
const double bead_radii = bead_radii_temp ; 
const int NrParticles = NrParticles_temp;

const double kb = kb_temp;
const double T0 = T0_temp;
const double sqrt_2kbTdt= sqrt(2.0*kb*T0*dt) ;

const double shear_rate = shear_rate_temp ; 

const double mu = 1.0/(6.0*pi*eta_0); // mu - mobility, eta - viscosity, r-radius of particle suspensions
const double mu_sqrt=sqrt(mu);

const int  cubic = 1 ; 	// cubic box 
	
const double Lx = pow(NrParticles*(4.0/3.0)*M_PI*(bead_radii*bead_radii*bead_radii),1.0/3.0 ); 
const double Ly = Lx ; // assuming cubic 
const double Lz = Lx ; 
const double Volume =Lx*Ly*Lz;
const double Volume_inv = 1.0/ Volume;
const double Particle_radius = 0.5 ; // sigma/2.0;
const double Particle_vol = 4.0*pi*(Particle_radius*Particle_radius*Particle_radius)/3.0;
const double vol_frac = (double) NrParticles * Particle_vol * Volume_inv;
const int cellx=(int) ceil(Lx/r_cut);
const int celly=(int) ceil(Ly/r_cut);
const int cellz=(int) ceil(Lz/r_cut);
 // based on the equation motion following wouter's eq. 327
 
 const double force_norm 	= 1.0/(eta_6pi*bead_radii*bead_radii)				; // since eta_s in wouter's note is eta)*6*pi here 
 const double torque_norm 	= 1.0/(eta_6pi*bead_radii*bead_radii*bead_radii)	;
 const double pos_norm 		= 1.0/(bead_radii)		;
 const double vel_norm 		= 1.0/(bead_radii)		;
 const double stochas_norm	= sqrt(2.0*kb*T0*dt/(eta_6pi*bead_radii*bead_radii*bead_radii))		;

const vctr3D box(Lx, Ly, Lz);
const vctr3D rbox(1.0/Lx, 1.0/Ly, 1.0/Lz);

int if_create_particles = xxcreate, ifrestart=xxrestart;
double tauT=0.1;
double Temp=T0;
// double shear_rate = 0.0; //shear rate
int ifshear = 0;// set equal to 1 for shear
std::string dataFileName="../xxx",dataFileName_new="../xxxnew" ;
double simu_time=dt;
long long int step=0, nSteps=10000, frame=100;
double vel_scale;
int if_Periodic =1;
int Max_Cluster_N=NrParticles;
int NrSubs=NrParticles;
int restart_frame_offset=0;

vctr5D Stresslet_mean = null5D;

double max_cos=0.0,min_cos=0.0,min_tan=0.0,max_tan=0.0, cos_val=0.0,tan_val=0.0;


// variables for Electric field 
vctr3D dipole_b(0.0,0.0,1.0);
vctr3D dipole_s(0.0,0.0,1.0);
double cos_theta = dipole_s*Elec_fld;
int hist_x[50]={};
int hist_phi_x[50]={};
int hist_y[50]={};
int hist_phi_y[50]={};
int hist_z[50]={};
int hist_phi_z[50]={};
int hist_C[100]={}; // histogram of orbit constant of ellipsids

// variables for polar , azimuthal angle histogram
double nbins = 100.0; 
double bin_cos = M_PI/nbins;			
double bin_tan = 2.0*M_PI/nbins;

cout << bin_cos << '\t' << bin_tan << endl;

long long int hist_pol_azi[100][100]={};

vctr3D dr_vec;
vector<SubData>  particle(NrParticles);
vector<ParticleData>  cluster( 1, ParticleData(NrSubs) );;

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
        currentLine >> particle[i].pos.comp[0];
        currentLine >> particle[i].pos.comp[1];
        currentLine >> particle[i].pos.comp[2];
    //   cout<< particle[i].pos.comp[0]<<endl;
	//	std::getline(dataFile,line);        
	//	particle[i].pos.comp[2]+=10.0;
	//	particle[i].pos.comp[2]=particle[i].pos.comp[2]*(-1.0);
	//	particle[i].pos.comp[0]=particle[i].pos.comp[0]*(-1.0);
	//	particle[i].pos.comp[1]=particle[i].pos.comp[1]*(-1.0);
    }
}
if (xx_diffusion) {			// calculate the diffusion tensor for the particles read-in ; if starting with paticles of varies shape then initializing diffusion tensors
	
	/* sort particles into cells */
for ( int i = 0 ; i < 1; i ++ )
	{
		if (!xxcluster_restart) 
			{
				cluster[i].Sub_Length=Max_Cluster_N;		// initially each cluster has size one
				cluster[i].mass=Max_Cluster_N;
				cluster[i].vel={0.0,0.0,0.0};

				// intialize Q, A matrix

				cluster[i].quat={1.0,0.0,0.0,0.0};
				cluster[i].quat2rotmat();
		}
	for ( int j = 0 ; j < cluster[i].Sub_Length ; j ++ )
		{ 
		if (!xxcluster_restart) {
			cluster[i].sub[j]=j;
			particle[cluster[i].sub[j]].cluster=0;
			particle[cluster[i].sub[j]].mass=1.0;
			particle[cluster[i].sub[j]].radius=0.5;
			cluster[i].radii=0.56;
			// particle[i].pos is the position of cluster, and particle[i].sub[i].pos is the spaced fixed position of particles in the cluster; initially all clusters have 1 paricle per cluster, and position of cluster is same as position of spaced fixed sub-particle 
			particle[cluster[i].sub[j]].vel=cluster[i].vel;
			particle[cluster[i].sub[j]].pos_bdyfxd=particle[cluster[i].sub[j]].pos;//cluster[i].sub[j].pos;
			cluster[i].pos={0.0,0.0,0.0};
			cluster[i].omega={0.0,0.0,0.0};//{((double) rand()/(RAND_MAX)-0.5),((double) rand()/(RAND_MAX)-0.5),((double) rand()/(RAND_MAX)-0.5)};
			cluster[i].angmom={((double) rand()/(RAND_MAX)-0.5),((double) rand()/(RAND_MAX)-0.5),((double) rand()/(RAND_MAX)-0.5)};
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
			
//	 forceUpdate( particle, &p_energy, &combine_now , combine, 0, NrParticles, Lx, Ly, Lz );

// calculate new diffusion tensors	
	for ( int i = 0 ; i < 1; i ++ )
		{
        
        
ifstream File ( "data_binary.bin" , ios::in | ios::binary );
if(!File.good()) {
	std::cerr<<"Given file is corrupt /n"<<std::endl;
}
else {			
		double xi_11x11[36];
		double temp_mu;
		
		for (int l=0; l<36; l++)
				{
					File.read( (char*) &xi_11x11[l]     , sizeof(xi_11x11[l]     ) );
					cout << temp_mu << '\n';
				}
				

					cluster[i].mobility_tnsr.comp[0][0] = xi_11x11[0] ;  
					cluster[i].mobility_tnsr.comp[0][1] = xi_11x11[1]  ;  
					cluster[i].mobility_tnsr.comp[0][2] = xi_11x11[2]  ; 
					cluster[i].mobility_tnsr.comp[1][0] = xi_11x11[6]  ; 
					cluster[i].mobility_tnsr.comp[1][1] = xi_11x11[7]  ;  
					cluster[i].mobility_tnsr.comp[1][2] = xi_11x11[8]  ;  
					cluster[i].mobility_tnsr.comp[2][0] = xi_11x11[12]  ;   
					cluster[i].mobility_tnsr.comp[2][1] = xi_11x11[13]  ; 
					cluster[i].mobility_tnsr.comp[2][2] = xi_11x11[14]  ; 				

					cluster[i].rot_mobility_tnsr_rt.comp[0][0] = xi_11x11[18] ;  
					cluster[i].rot_mobility_tnsr_rt.comp[0][1] = xi_11x11[19]  ;  
					cluster[i].rot_mobility_tnsr_rt.comp[0][2] = xi_11x11[20]  ; 
					cluster[i].rot_mobility_tnsr_rt.comp[1][0] = xi_11x11[24]  ; 
					cluster[i].rot_mobility_tnsr_rt.comp[1][1] = xi_11x11[25]  ;  
					cluster[i].rot_mobility_tnsr_rt.comp[1][2] = xi_11x11[26]  ;  
					cluster[i].rot_mobility_tnsr_rt.comp[2][0] = xi_11x11[30]  ;   
					cluster[i].rot_mobility_tnsr_rt.comp[2][1] = xi_11x11[31]  ; 
					cluster[i].rot_mobility_tnsr_rt.comp[2][2] = xi_11x11[32]  ; 

					cluster[i].mobility_tnsr_tr.comp[0][0] = xi_11x11[3] ;  
					cluster[i].mobility_tnsr_tr.comp[0][1] = xi_11x11[4]  ;  
					cluster[i].mobility_tnsr_tr.comp[0][2] = xi_11x11[5]  ; 
					cluster[i].mobility_tnsr_tr.comp[1][0] = xi_11x11[9]  ; 
					cluster[i].mobility_tnsr_tr.comp[1][1] = xi_11x11[10]  ;  
					cluster[i].mobility_tnsr_tr.comp[1][2] = xi_11x11[11]  ;  
					cluster[i].mobility_tnsr_tr.comp[2][0] = xi_11x11[15]  ;   
					cluster[i].mobility_tnsr_tr.comp[2][1] = xi_11x11[16]  ; 
					cluster[i].mobility_tnsr_tr.comp[2][2] = xi_11x11[17]  ; 				
																					
					cluster[i].rot_mobility_tnsr.comp[0][0] = xi_11x11[21] ;  
					cluster[i].rot_mobility_tnsr.comp[0][1] = xi_11x11[22]  ;  
					cluster[i].rot_mobility_tnsr.comp[0][2] = xi_11x11[23]  ; 
					cluster[i].rot_mobility_tnsr.comp[1][0] = xi_11x11[27]  ; 
					cluster[i].rot_mobility_tnsr.comp[1][1] = xi_11x11[28]  ;  
					cluster[i].rot_mobility_tnsr.comp[1][2] = xi_11x11[29]  ;  
					cluster[i].rot_mobility_tnsr.comp[2][0] = xi_11x11[33]  ;   
					cluster[i].rot_mobility_tnsr.comp[2][1] = xi_11x11[34]  ; 
					cluster[i].rot_mobility_tnsr.comp[2][2] = xi_11x11[35]  ; 				
																			

		for (int l=0; l<3; l++)
			{
				for (int k=0; k<5; k++)
					{
							File.read( (char*) &temp_mu     , sizeof(temp_mu     ) );
							cluster[i].mobility_tnsr_td.comp[l][k] = temp_mu;
							cout << temp_mu << '\t';
					}
			}
			
		for (int l=3; l<6; l++)
			{
				for (int k=0; k<5; k++)
					{
							File.read( (char*) &temp_mu     , sizeof(temp_mu     ) );
							cluster[i].mobility_tnsr_rd.comp[l-3][k] = temp_mu;
							cout << temp_mu << '\t';
					}
			}
			
		for (int l=0; l<5; l++)
			{
				for (int k=0; k<5; k++)
					{
							File.read( (char*) &temp_mu     , sizeof(temp_mu     ) );
							cluster[i].mobility_tnsr_dd.comp[l][k] = temp_mu;
							cout << temp_mu << '\n';
					}
			}	
        
 }
/*        
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
*/
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
	//	cluster[i].quat={0.8467   , 0.5320    ,   0.0   ,      0.0 };
	//	cluster[i].quat={0.7071   , -0.7071     ,  0.0   ,      0.0 };
	//	cluster[i].quat={0.972369920397677,	0.233445363855905,	0.,	0.};
	//	cluster[i].quat={0.9239  ,  0.3827   ,      0.0     ,    0.0};
	//	cluster[i].quat={0.987688340595138,	0.156434465040231,	0.0,	0.0};	//   pi/10;
	//	cluster[i].quat={0.951056516295154,	0.309016994374947,	0.0,	0.0};	// 2*pi/10;
		cluster[i].quat={0.891006524188368,	0.453990499739547,	0.0,	0.0};	// 3*pi/10;
	//	cluster[i].quat={0.809016994374948,	0.587785252292473,	0.0,	0.0};	// 4*pi/10;
	//	cluster[i].quat={0.707106781186548,  0.707106781186548, 0.0,    0.0};	// 5*pi/10;

		// update A matrix

        cluster[i].quat2rotmat();
        cluster[i].rotmat.echo();
	}
	
}
}
Max_Cluster_N =1;

/*
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
*/

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

std::ofstream outFile10(dataFileName+"/End_positions.dat");
std::ofstream outFile_com(dataFileName+"/com.dat");
std::ofstream outFile_orient(dataFileName+"/orient.dat");

// forceUpdate( particle, &p_energy, &combine_now , combine, &step, NrParticles , Lx, Ly, Lz);


// convert subforces into total generalized forces on particles 

/*
// For electric field generated torque
  for ( int i = 0 ; i < 1; i ++ )
  {
	cluster[i].frc=null3D;
	cluster[i].trq=dipole_s^Elec_fld;
	cluster[i].Iner_tnsr=null33D;
  }
*/
// For electric field generated force with gives torque

particle[0].charge = 1.0;
particle[1].charge = -1.0;
  for ( int i = 0 ; i < Max_Cluster_N; i ++ )
   {
   for (int  j = 0 ; j < cluster[i].Sub_Length ; j ++ )
		{
               particle[cluster[i].sub[j]].frc = Elec_fld*particle[cluster[i].sub[j]].charge;                                                                                                                                                                                                                              //      Modification of Numerical Model for Ellipsoidal Monomers by Erwin Gostomski
		}
   }

/*
  for ( int i = 0 ; i < 1; i ++ )
  {
	cluster[i].frc=null3D;
	cluster[i].trq=null3D;
	cluster[i].Iner_tnsr=null33D;
  }
  */
  
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
               cluster[i].Iner_tnsr+=(I_sphere+Unit_diag*(dr_vec.norm2())-dr_mat)*particle[cluster[i].sub[j]].mass;    //      refer following paper , page 3 equa. 3 for interia tensor formula
                                                                                                                                                                                                                               //      Modification of Numerical Model for Ellipsoidal Monomers by Erwin Gostomski
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
	outFile8<<"Viscosity, eta"<<'\t'<<eta_0<<std::endl;
	outFile8<<"Mobility , mu"<<'\t'<<mu<<std::endl;
	outFile8<<'\n'<<" Data Folder and Git Vesrion : "<<'\n';
	system(" echo >> logfile & git log --pretty=format:'%h' -n 1 >> logfile   & echo >> logfile  &  pwd >> logfile & ");
	outFile8.close();

	std::ofstream outFile_inter_cluster("inter_cluster_data.dat");
	std::ofstream Stresslet_data("Stresslet_data.dat");
	cout<<step<<endl;

	
std::ofstream outFile12("vec1.dat");
std::ofstream outFile13("vec2.dat");
std::ofstream outFile14("vec3.dat");
vctr3D eig1(1.0 , 0.0 , 0.0 ); 
vctr3D eig2(0.0 , 1.0 , 0.0 );
vctr3D eig3(0.0 , 0.0 , 1.0 );
vctr3D vec1, vec2, vec3;


// create a grid and bin the points of geodesic dome

	

 const int dm27[26][3] = { 	{  1,  0,  0 },
							{ -1,  0,  0 },
							{  0,  1,  0 },
							{  0, -1,  0 },
							{  0,  0,  1 },
							{  0,  0, -1 },
							{  1,  1,  0 },
							{  1, -1,  0 },
							{ -1,  1,  0 },
							{ -1, -1,  0 },
							{  1, 0,  1 },
							{  1, 0, -1 },
							{ -1, 0,  1 },
							{ -1, 0, -1 },
							{  0, 1,  1 },
							{  0, 1, -1 },
							{  0,-1,  1 },
							{  0,-1, -1 },
							{  1,  1,  1 },
							{  1, -1,  1 },
							{ -1,  1,  1 },
							{ -1, -1,  1 },
							{  1,  1, -1 },
							{  1, -1, -1 },
							{ -1,  1, -1 },
							{ -1, -1, -1 },
 };

	double binsSize[3];
	double boxEdge = 1.0*2.0+2.0*0.3; 	// max cubic box size enveloping unit sphere (i.e. enveloping the geodesic) 
	double rboxEdge = 1.0/boxEdge ;
	vctr3D envbox	= 	{boxEdge,boxEdge,boxEdge};
	vctr3D envRbox  =	{rboxEdge,rboxEdge,rboxEdge};
	vctr3D havenvbox  =	{boxEdge/2.0,boxEdge/2.0,boxEdge/2.0};
	
	double max_sepr = 0.1345; 	// max distance between neighbhoring points on the geodesic
	int    Nrbins[3],MaxNrbins ; 
	double gridsScale[3];
	int NrPoints = 1002 ; 		// no. of points on geodesic
	int orientHist[NrPoints] = {};
	int gridUpdate ; 
	int MaxPrperCell = 10 ; 	// max no. of points per gridcell
	int    i,j;
	int    ii,jj;
	int    mi[3],m,mj[3];
	vctr3D dr,dR;
	double theta, maxtheta ; 
	
	// read in geodesic points 
	
	double GeodesicPt[NrPoints][3] = {} ; 
	
	fileName="GeodesicPt_1002.dat";
	dataFile;

	dataFile.open(fileName);

	if(!dataFile.good()) {
		std::cerr<<"Given file is corrupt : GeodesicPt_1002.dat  /n"<<std::endl;
	}
	else {
		std::string line;
		for (int i=0;i<NrPoints;i++) {
			std::getline(dataFile,line);
			std::istringstream currentLine(line);    
			currentLine >> GeodesicPt[i][0];
			currentLine >> GeodesicPt[i][1];
			currentLine >> GeodesicPt[i][2];
		}
	}	
	
	dataFile.close();  
	dataFile.clear();	
	
	for ( i = 0 ; i < 3 ; i++ )
	{
    Nrbins[i] = floor ( envbox.comp[i] / (max_sepr) ); // cellnr runs from 0 to NrCells-1
    gridsScale[i] = Nrbins[i] * envRbox.comp[i];
    if ( Nrbins[i] < 3 ) { cout << "*** Nrbins[" << i << "] = " << Nrbins[i] << endl ; abort(); }
  }

// periodic boundary conditions

  MaxNrbins = max( max( Nrbins[0], Nrbins[1] ), Nrbins[2]);
  
// generate grid list
 	int binGrid[Nrbins[0]][Nrbins[1]][Nrbins[2]][MaxPrperCell+1];

  for ( mi[0] = 0 ; mi[0] < Nrbins[0] ; mi[0]++ )
  {
    for ( mi[1] = 0 ; mi[1] < Nrbins[1] ; mi[1]++ )
    {
      for ( mi[2] = 0 ; mi[2] < Nrbins[2] ; mi[2]++ )
      {

        binGrid[mi[0]][mi[1]][mi[2]][0] = 0;
      } // miz
    } // miy
  } // mix

for ( int i = 0 ; i < NrPoints ; i ++ )
  {
	  
    mi[x] = int ( (GeodesicPt[i][0]+havenvbox.comp[0]) * gridsScale[0] );
    mi[y] = int ( (GeodesicPt[i][1]+havenvbox.comp[1]) * gridsScale[1] );
    mi[z] = int ( (GeodesicPt[i][2]+havenvbox.comp[2]) * gridsScale[2] );       

    
    if ( int (binGrid[mi[0]][mi[1]][mi[2]][0]) >= MaxPerCell-1 )
    {
      cout << "*** cell overfull" << endl;
      cout << mi[0] << "  " << mi[1] << "  " << mi[2] << endl;
      abort();
    }

    binGrid[mi[0]][mi[1]][mi[2]][0] ++ ;
    
    binGrid[mi[0]][mi[1]][mi[2]][ int (binGrid[mi[0]][mi[1]][mi[2]][0])] = i;

} // i
	//		 cout<<"exit grid "<<endl;
 
// end grid creation 
cout << "cehck_here" << endl;

simu_time =dt;
do {

	brownian(step, cluster, particle, &Max_Cluster_N , vel_scale, force_norm, torque_norm, pos_norm, vel_norm, stochas_norm , shear_rate , sqrt_2kbTdt, box, rbox)	;


  	dipole_s  = cluster[0].rotmat*dipole_b;			// rotate the body fixed dipole
	cos_theta = dipole_s*Elec_fld;					// calucate the angle between dipole and electric field, dot product gives the cosine of angle

	// hist[int (floor((cos_theta+5.0)/0.4))]+=1;

	double phi;
	
	cos_theta = cluster[0].rotmat.comp[0][2];
	phi = atan2(cluster[0].rotmat.comp[1][2],cluster[0].rotmat.comp[0][2]);
	
	hist_phi_z[int (floor((phi+M_PI)/0.1257))]+=1;
	hist_z[int (floor((cos_theta+1.0)/0.04))]+=1;

	cos_theta = cluster[0].rotmat.comp[0][1];
	phi = atan2(cluster[0].rotmat.comp[1][1],cluster[0].rotmat.comp[0][1]);
	
	hist_phi_y[int (floor((phi+M_PI)/0.1257))]+=1;
	hist_y[int (floor((cos_theta+1.0)/0.04))]+=1;
	
	cos_theta = cluster[0].rotmat.comp[0][0];
	phi = atan2(cluster[0].rotmat.comp[1][0],cluster[0].rotmat.comp[0][0]);
	
	hist_phi_x[int (floor((phi+M_PI)/0.1257))]+=1;
	hist_x[int (floor((cos_theta+1.0)/0.04))]+=1;


// 	forceUpdate( particle, &p_energy, &combine_now , combine, &step , NrParticles, Lx, Ly, Lz);


// convert subforces into total generalized forces on particles 

/*
//  For electric field generated torque
  for ( int i = 0 ; i < 1; i ++ )
  {
	cluster[i].frc=null3D;
	cluster[i].trq=dipole_s^Elec_fld;
	cluster[i].Iner_tnsr=null33D;
  }
  */
/*
  // For electric field generated force with gives torque

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
               cluster[i].Iner_tnsr+=(I_sphere+Unit_diag*(dr_vec.norm2())-dr_mat)*particle[cluster[i].sub[j]].mass;    //      refer following paper , page 3 equa. 3 for interia tensor formula
                                                                                                                                                                                                                               //      Modification of Numerical Model for Ellipsoidal Monomers by Erwin Gostomski
		}
   }
   
*/
/*
  for ( int i = 0 ; i < 1; i ++ )
  {
	cluster[i].frc=null3D;
	cluster[i].trq=null3D;
	cluster[i].Iner_tnsr=null33D;
  }


  // for rotational relaxation check
 
			vec1 = cluster[0].rotmat*eig1;
			vec2 = cluster[0].rotmat*eig2;
			vec3 = cluster[0].rotmat*eig3;
			
		outFile12<<vec1.comp[0] <<'\t'<< vec1.comp[1] << '\t'<< vec1.comp[2] <<  endl;      
		outFile13<<vec2.comp[0] <<'\t'<< vec2.comp[1] << '\t'<< vec2.comp[2] <<  endl;      
		outFile14<<vec3.comp[0] <<'\t'<< vec3.comp[1] << '\t'<< vec3.comp[2] <<  endl;      
						
		outFile_com<<cluster[0].pos.comp[0]<<'\t'<<cluster[0].pos.comp[1]<<'\t'<<cluster[0].pos.comp[2]<<'\t'<<std::endl;
*/

	Stresslet_mean += cluster[0].Stresslet;
	
	// binning of orbital constant of ellipsoid , C
			
	vctr3D director = {cluster[0].rotmat.comp[0][2],cluster[0].rotmat.comp[1][2],cluster[0].rotmat.comp[2][2]}; 
	
	double tan_phi =  (director.comp[1])/(director.comp[0]);
	double tan_phi2 = tan_phi*tan_phi;
	double cos_theta = director.comp[2] ; 
	double cos_theta2 = cos_theta*cos_theta;
	double cos_phi = sqrt(1.0/(tan_phi2+1.0));
	double tan_theta = sqrt(1.0/cos_theta2-1.0); 
	double ar = 5.0;
	double C = tan_theta*cos_phi*sqrt(ar*ar+tan_phi2);
	double max_C_lim = ceil(M_PI/2.0) ; 
	// outFile_com<< C<< '\t'<< atan(C) <<std::endl;

	int i = round(atan(C)/(max_C_lim/100.0)); // (10.0*M_PI/2.0)
	// i = min( max(i,0), max_C_lim ) ;
	hist_C[i]++; 
	
	// end of C binning 

/*
if (step%(frame)==0) 
	{ 


		// save position every 'frame' steps 
		
		for ( int i = 0 ; i < 1; i ++ )
			{
				if(cluster[i].Sub_Length>0)
				{
				Stresslet_data.precision(17);
				Stresslet_data<<cluster[i].Stresslet.comp[0]<<'\t'<<cluster[i].Stresslet.comp[1]<<'\t'<<cluster[i].Stresslet.comp[2]<<'\t'<<cluster[i].Stresslet.comp[3]<<'\t'<<cluster[i].Stresslet.comp[4]<<'\t'<<std::endl;	
	//			outFile_com<<cos_theta<<'\t'<<cos_theta<<'\t'<<cos_theta<<std::endl;
	//			vctr3D director = cluster[0].rotmat*particle[cluster[0].sub[0]].pos_bdyfxd ; 
				vctr3D director = {cluster[0].rotmat.comp[0][2],cluster[0].rotmat.comp[1][2],cluster[0].rotmat.comp[2][2]}; 
	//			outFile_orient<<cluster[i].rotmat.comp[0][1]<<'\t'<<cluster[i].rotmat.comp[1][1]<<'\t'<<cluster[i].rotmat.comp[2][1]<<'\t'<<std::endl;

					double t3 = 2.0 * (cluster[i].quat.comp[0] * cluster[i].quat.comp[3]+ cluster[i].quat.comp[1] * cluster[i].quat.comp[2]);
					double t4 = 1.0 - 2.0 * (cluster[i].quat.comp[2] * cluster[i].quat.comp[2] + cluster[i].quat.comp[3] * cluster[i].quat.comp[3]);  
					double yaw = std::atan2(t3, t4);	; // rotation about Z-axis


		// update bin 
		
		mi[0] = int ( (director.comp[0]+havenvbox.comp[0]) * gridsScale[0] );
		mi[1] = int ( (director.comp[1]+havenvbox.comp[1]) * gridsScale[1] );
		mi[2] = int ( (director.comp[2]+havenvbox.comp[2]) * gridsScale[2] ); 
          // particle j in same cell as i
          maxtheta = 0;

          for ( int jj =  1 ; jj <= binGrid[mi[0]][mi[1]][mi[2]][0] ; jj++ )
          {
			int j = binGrid[mi[0]][mi[1]][mi[2]][jj];
			theta =  director.comp[0]*GeodesicPt[j][0] + director.comp[1]*GeodesicPt[j][1] + director.comp[2]*GeodesicPt[j][2]; 
			if (theta > maxtheta)
			{
				gridUpdate = j ;
				maxtheta = theta ; 
			}
			
          } // jj

          // particle j in neighbour cell to i
          for ( m = 0 ; m < 26 ; m++ )
          {
            mj[0]      =  mi[0] + dm[m][0] ;
            mj[1]      =  mi[1] + dm[m][1] ;
            mj[2]      =  mi[2] + dm[m][2] ;
            if ( mj[0] >= Nrbins[0] ||  mj[0] < 0 || mj[1] >= Nrbins[1] ||  mj[1] < 0 || mj[2] >= Nrbins[2] ||  mj[2] < 0 ) continue ; 	// going out of envbox bounds
			
            for ( jj = 1 ; jj <= binGrid[mj[0]][mj[1]][mj[2]][0] ; jj++ )
            {
				j = binGrid[mj[0]][mj[1]][mj[2]][jj];

				theta =  director.comp[0]*GeodesicPt[j][0] + director.comp[1]*GeodesicPt[j][1] + director.comp[2]*GeodesicPt[j][2]; 

				if (theta > maxtheta)
				{
					gridUpdate = j ; 
					maxtheta = theta ; 
				}				
            } // jj
          } // m

		orientHist[gridUpdate]++ ;
		


			//	outFile_orient<<director.comp[0]<<'\t'<<director.comp[1]<<'\t'<<director.comp[2]<<'\t'<< yaw << std::endl;
				
			//	outFile_com<<cluster[0].pos.comp[0]<<'\t'<<cluster[0].pos.comp[1]<<'\t'<<cluster[0].pos.comp[2]<<'\t'<<std::endl;

				}
			}

	}
*/	
/*
if (step%(frame)==0) 
	{ 

		std::ofstream outFile5(dataFileName+"/XYZ"+ std::to_string(step/(frame)) +".xyz");   
   		outFile5<<NrParticles<<std::endl;
   		outFile5<<"X Y Z co-ordinates"<<std::endl;

		// save position every 'frame' steps 
		
		for ( int i = 0 ; i < 1; i ++ )
			{

			    for (int  j = 0 ; j < cluster[i].Sub_Length ; j ++ )
					{

						particle[cluster[i].sub[j]].pos = cluster[i].pos + cluster[i].rotmat*particle[cluster[i].sub[j]].pos_bdyfxd;
					//	particle[cluster[i].sub[j]].pos.PBC(box,rbox);
						
					outFile5<<'H'<<'\t'<<particle[cluster[i].sub[j]].pos.comp[0]<<'\t'<<particle[cluster[i].sub[j]].pos.comp[1]<<'\t'<<particle[cluster[i].sub[j]].pos.comp[2]<<'\t'<<i<<std::endl;
					}
			}
 
     	outFile5<<'\n'<<std::endl;
		outFile5.close();
	}
*/
	simu_time+=dt;
	step+=1;

} while(xxnstep);
cout << step << endl;
/*
// output the histogram of the cosine angles
for ( int i = 0 ; i < 50; i ++ )
	{
		cout << hist_x[i] << '\t' << hist_phi_x[i] << '\t' << hist_y[i] << '\t' << hist_phi_y[i] << '\t' << hist_z[i] << '\t' << hist_phi_z[i] << endl;
	}



 // output the histogram of the polar and azimuthal angles
for ( int i = 0 ; i < 100; i ++ )
	{
		for ( int j = 0 ; j < 100; j ++ )
		{
			cout<< std::setprecision(5) << hist_pol_azi[i][j] <<'\t';
		}
		cout << endl;
	}
*/

 // output the histogram of the geodesic
 			cout<< "histogram of the geodesic" << endl;

for ( int i = 0 ; i < 100; i ++ )
	{
			outFile_orient<< hist_C[i] << endl;
	}
	cout<< max_cos << '\t' << min_cos << '\t' << max_tan << '\t'  << min_tan << endl;

	cout<<Stresslet_mean.comp[0]<<'\t'<<Stresslet_mean.comp[1]<<'\t'<<Stresslet_mean.comp[2]<<'\t'<<Stresslet_mean.comp[3]<<'\t'<<Stresslet_mean.comp[4]<<'\t'<<std::endl;	
	
					
	std::ofstream outFile_rand_state(dataFileName+"/random_device_state_new.txt");
	outFile_rand_state << gen;
	outFile_rand_state.close();

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
