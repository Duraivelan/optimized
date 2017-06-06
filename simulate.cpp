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
		vctr3D u_inf(shear_rate*cluster[i].pos.comp[1],0.0,0.0); 		// shear flow gradient in y-direction
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
										+  cluster[i].rot_mobility_tnsr_sqrt*(rand1*kbT_dt)
										/*-  (cluster[i].mobility_tnsr_rd*E_inf_bt)*dt*/; 	// body fixed omega
				cluster[i].omega	=	w_inf*dt;						// space-fixed omega
				cluster[i].quat		= cluster[i].theta2quat() /*+ cluster[i].omega2qdot()*/ ;
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
int step=0, nSteps=10000, frame=100000;
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

vctr3D dipole_b(0.0,0.0,1.0);
vctr3D dipole_s(0.0,0.0,1.0);
double cos_theta = dipole_s*Elec_fld;
int hist[25]={};

vector<SubData>  particle(NrParticles);
vector<ParticleData>  cluster( 1, ParticleData(NrSubs) );
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
        currentLine >> particle[i].pos.comp[0];
        currentLine >> particle[i].pos.comp[1];
        currentLine >> particle[i].pos.comp[2];
       cout<< particle[i].pos.comp[0]<<endl;
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
//	 forceUpdate( particle, &p_energy, &combine_now , combine, 0);
// cout<< combine_now<<"combine123"<<endl;

// calculate new diffusion tensors	
	for ( int i = 0 ; i < 1; i ++ )
		{
		remove("new_cluster.dat");
	//	remove("data.dat");

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
        
      /*  cluster[i].rotmat = {
							{1.0,		0.0,	0.0},
							{0.0,	cos(0.33),	sin(0.33)},
							{0.0,	-sin(0.33),	cos(0.33)},
							};
        */
		cluster[i].rotmat.echo();
        
        
		for (int  k=0; k<cluster[i].Sub_Length; k++) {
			
	//		particle[cluster[i].sub[k]].pos = cluster[i].rotmat*particle[cluster[i].sub[k]].pos;

	//		particle[cluster[i].sub[k]].pos_bdyfxd = particle[cluster[i].sub[k]].pos;
	//	particle[cluster[i].sub[k]].pos_bdyfxd.comp[1]+=5.1861;
	//	particle[cluster[i].sub[k]].pos_bdyfxd=particle[cluster[i].sub[k]].pos;

		outFile7<<particle[cluster[i].sub[k]].pos_bdyfxd.comp[0]<<'\t'<<particle[cluster[i].sub[k]].pos_bdyfxd.comp[1]<<'\t'<<particle[cluster[i].sub[k]].pos_bdyfxd.comp[2]<<'\t'<<"0.4"<<std::endl;
		
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
	
}
}
Max_Cluster_N =1;
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

// forceUpdate( particle, &p_energy, &combine_now , combine, &step);

	// convert subforces into total generalized forces on particles 

  for ( int i = 0 ; i < 1; i ++ )
  {
	cluster[i].frc=null3D;
	cluster[i].trq=dipole_s^Elec_fld;
	cluster[i].Iner_tnsr=null33D;
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
	cout<<step<<endl;

simu_time =dt;
do {
	p_energy=0;	

	brownian(step, cluster, particle, &Max_Cluster_N , &KE_rot, vel_scale )	;
	dipole_s  = cluster[0].rotmat*dipole_b;
	cos_theta = dipole_s*Elec_fld;

	hist[int (floor((cos_theta+5.0)/0.4))]+=1;
//				outFile_com<<cos_theta<<std::endl;

// 	forceUpdate( particle, &p_energy, &combine_now , combine, &step);
//	cout<<step<<endl;
	// convert subforces into total generalized forces on particles 
  for ( int i = 0 ; i < 1; i ++ )
  {
	cluster[i].frc=null3D;
	cluster[i].trq=dipole_s^Elec_fld;
	cluster[i].Iner_tnsr=null33D;
  }

if (step%frame==0) 
	{ 

      std::ofstream outFile5(dataFileName+"/XYZ"+ std::to_string(step/frame) +".xyz");   
   		outFile5<<NrParticles<<std::endl;
   		outFile5<<"X Y Z co-ordinates"<<std::endl;
		outFile11<<step<<'\t'<<Max_Cluster_N<<std::endl;
		// save position, Kinetic energy, Potential energy, Forces every 'frame' steps and also store radii of gyration info
		
	//	std::ofstream outFile9(dataFileName+"/Cluster_dist"+ std::to_string(step/frame) +".dat");

		K_Energy=0;

		for ( int i = 0 ; i < 1; i ++ )
			{
				if(cluster[i].Sub_Length>0)
				{
				Stresslet_data<<cluster[i].Stresslet.comp[0]<<'\t'<<cluster[i].Stresslet.comp[1]<<'\t'<<cluster[i].Stresslet.comp[2]<<'\t'<<cluster[i].Stresslet.comp[3]<<'\t'<<cluster[i].Stresslet.comp[4]<<'\t'<<std::endl;	
	//			outFile9<<cluster[i].radii_gyr<<'\t'<<cluster[i].Sub_Length<<std::endl;
	//			outFile_com<<cos_theta<<'\t'<<cos_theta<<'\t'<<cos_theta<<std::endl;
			/*	outFile_orient<<particle[cluster[i].sub[cluster[i].Sub_Length-1]].pos.comp[0]-particle[cluster[i].sub[0]].pos.comp[0]<<'\t'<<
				particle[cluster[i].sub[cluster[i].Sub_Length-1]].pos.comp[1]-particle[cluster[i].sub[0]].pos.comp[1]<<'\t'<<
				particle[cluster[i].sub[cluster[i].Sub_Length-1]].pos.comp[2]-particle[cluster[i].sub[0]].pos.comp[2]<<'\t'<<i<<std::endl; */
				outFile_orient<<cluster[i].rotmat.comp[0][2]<<'\t'<<cluster[i].rotmat.comp[1][2]<<'\t'<<cluster[i].rotmat.comp[2][2]<<'\t'<<std::endl;
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

for ( int i = 0 ; i < 25; i ++ )
	{
		cout << hist[i] << endl;
	}
	

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
