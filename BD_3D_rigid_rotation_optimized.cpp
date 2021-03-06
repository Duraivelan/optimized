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
#include</home/duraivelan/Downloads/eigen-eigen-10219c95fe65/Eigen/Eigenvalues>
//#include<Eigen/Eigenvalues>

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


void Collision(vector<SubData>& particle, vector<ParticleData>& cluster, int i, int j,int *Max_Cluster_N, vctr3D box, vctr3D rbox ) {
	  
	    int temp_min, temp_max;
	    

		 while (cluster[i].linked != 0 ) 
		{ 	
			cout << i << '\t'<< "inside" <<endl;
			i = cluster[i].link_clstr_ptr ;
		}
		while (cluster[j].linked != 0 ) 
		{ 	
			cout << j << '\t'<< "inside" <<endl;
			j = cluster[j].link_clstr_ptr ;
		} 	

		temp_min = min(i,j);
		temp_max = max(i,j);
		i = temp_min;
		j = temp_max; 
		cout << i << '\t'<< j <<endl;
		vctr3D L_after, L_before ;
		vctr3D old_pos= cluster[i].pos;
		vctr3D dr;
		double temp_r, r;
	
		cout<<particle[cluster[i].sub[0]].pos.comp[0]<<'\t'<<particle[cluster[i].sub[0]].pos.comp[1]<<'\t'<<particle[cluster[i].sub[0]].pos.comp[2]<<'\t'<<std::endl;
		cout<<particle[cluster[j].sub[0]].pos.comp[0]<<'\t'<<particle[cluster[j].sub[0]].pos.comp[1]<<'\t'<<particle[cluster[j].sub[0]].pos.comp[2]<<'\t'<<std::endl;
			
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
		
	
		cluster[i].angmom=(L_before - L_after);

	//	std::cout<<i<<'\t'<<j<<'\t'<<cluster[i].mass<<std::endl;

 
	//	std::cout<<"inside_collision"<<'\t'<<i<<'\t'<<A[i][0][0]<<'\t'<<Ang_Velocity[i][0]<<'\t'<<Q[i][0]<<'\t'<<L[i][0]<<'\t'<<torque[i][0]<<std::endl;
		
		remove("new_cluster.dat");
		remove("hydro++input.txt");

		std::ofstream outFile7("new_cluster.dat");
		std::ofstream outFile4("hydro++input.txt");

		outFile7<<1<<",    !Unit of length for coordinates and radii, cm (10 A)"<<endl;
		outFile7<<cluster[i].Sub_Length+cluster[j].Sub_Length<<",        !Number of beads"<<endl;

		cout<<particle[cluster[i].sub[0]].pos.comp[0]<<'\t'<<particle[cluster[i].sub[0]].pos.comp[1]<<'\t'<<particle[cluster[i].sub[0]].pos.comp[2]<<'\t'<<std::endl;
		cout<<particle[cluster[j].sub[0]].pos.comp[0]<<'\t'<<particle[cluster[j].sub[0]].pos.comp[1]<<'\t'<<particle[cluster[j].sub[0]].pos.comp[2]<<'\t'<<std::endl;


		cluster[i].radii_gyr=0.0;		     	     	
		
		for (int  k=0; k<cluster[i].Sub_Length; k++) {

		particle[cluster[i].sub[k]].pos_bdyfxd	 =  particle[cluster[i].sub[k]].pos-cluster[i].pos;
		particle[cluster[i].sub[k]].pos_bdyfxd.PBC(box,rbox);
	    cluster[i].radii_gyr+=particle[cluster[i].sub[k]].pos_bdyfxd.norm2()/(cluster[i].Sub_Length+cluster[j].Sub_Length);				
		outFile7<<particle[cluster[i].sub[k]].pos_bdyfxd.comp[0]<<'\t'<<particle[cluster[i].sub[k]].pos_bdyfxd.comp[1]<<'\t'<<particle[cluster[i].sub[k]].pos_bdyfxd.comp[2]<<'\t'<<particle[cluster[i].sub[k]].radius<<std::endl;

		
		} 
		
		for (int  k=cluster[i].Sub_Length; k<cluster[i].Sub_Length+cluster[j].Sub_Length; k++) {
		
	//	cluster[i].sub[k]	=	cluster[j].sub[k-cluster[i].Sub_Length];
			cluster[i].sub[k] = cluster[j].sub[k-cluster[i].Sub_Length];
			particle[cluster[i].sub[k]].pos_bdyfxd	 =  particle[cluster[i].sub[k]].pos-cluster[i].pos;
			particle[cluster[i].sub[k]].pos_bdyfxd.PBC(box,rbox);				
			cluster[i].radii_gyr+=particle[cluster[i].sub[k]].pos_bdyfxd.norm2()/(cluster[i].Sub_Length+cluster[j].Sub_Length);		
			outFile7<<particle[cluster[i].sub[k]].pos_bdyfxd.comp[0]<<'\t'<<particle[cluster[i].sub[k]].pos_bdyfxd.comp[1]<<'\t'<<particle[cluster[i].sub[k]].pos_bdyfxd.comp[2]<<'\t'<<particle[cluster[i].sub[k]].radius<<std::endl;			
			particle[cluster[i].sub[k]].cluster= i ;
			
		} 
		outFile7.close();
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


		
		outFile4<<"Square Tetramer                 Title"<<endl;
		outFile4<<"12-cluster                  filename for output files"<<endl;
		outFile4<<"new_cluster.dat              Structural (bead coords) file"<<endl;
		outFile4<<"12                              ICASE"<<endl;
		outFile4<<26.8500<<"                             Temperature, centigrade"<<endl;
		outFile4<<eta<<"                           Solvent viscosity"<<endl;
		outFile4<<cluster[i].Sub_Length<<"                          Molecular weight"<<endl;
		outFile4<<0.01<<"                           Specific volume of macromolecule"<<endl;
		outFile4<<"1.0                             Solution density"<<endl;
		outFile4<<"1,                 Number of values of H"<<endl;
		outFile4<<"0,             HMAX"<<endl;
		outFile4<<"1,                 Number of intervals for the distance distribution"<<endl;
		outFile4<<"0              RMAX"<<endl;
		outFile4<<"0,             (ONLY IF ISCA IS NOT ZERO) NTRIALS"<<endl;
		outFile4<<"1                   IDIF=1 (yes) for full diffusion tensors"<<endl;
		outFile4<<"*           End of file"<<endl;
		outFile4.close();
		
		system("../diffusion_tensor/hydro++10-lnx.exe < ../diffusion_tensor/input.txt  > /dev/null ");

		// cout<<"Done hydro"<<endl;
		std::ifstream dataFile("12-cluster-res.txt");
if(!dataFile.good()) {
	std::cerr<<"Given file is corrupt /n"<<std::endl;
}
else {
    std::string line;
    for (int n=0;n<56;n++) {
		std::getline(dataFile,line);
	}
    for (int n=0;n<3;n++) {
		std::getline(dataFile,line);
    	std::istringstream currentLine(line);    
        currentLine >> cluster[i].mobility_tnsr.comp[n][0];
        currentLine >> cluster[i].mobility_tnsr.comp[n][1];
        currentLine >> cluster[i].mobility_tnsr.comp[n][2];
    }
   for (int n=0;n<2;n++) {
		std::getline(dataFile,line);
	}
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
}	 
		cluster[i].mobility_tnsr=cluster[i].mobility_tnsr*(1*10*2414323832351.228);				// multiply by kBT (assuming kB in erg/K and T as 300 K ) correct for 1/kBT term included in the value 
		cluster[i].rot_mobility_tnsr=cluster[i].rot_mobility_tnsr*(1*10*2414323832351.228);		// outputed by hydro++
		cluster[i].mobility_tnsr_sqrt=null33D;
		MatrixXd temp(3,3), temp_sqrt(3,3);
		for (int k=0;k<3;k++) {
			for (int l=0;l<3;l++) {
				temp(k,l)=cluster[i].mobility_tnsr.comp[k][l];
			//	cout<<cluster[i].mobility_tnsr.comp[k][l]<<endl;

			}
		}
	Eigen::SelfAdjointEigenSolver<MatrixXd> TRANS_MOBL_MAT(temp);
	temp_sqrt = TRANS_MOBL_MAT.operatorSqrt();
	//	temp_sqrt=temp.sqrt();
		for (int k=0;k<3;k++) {
			for (int l=0;l<3;l++) {
				cluster[i].mobility_tnsr_sqrt.comp[k][l]=temp_sqrt(k,l);
		//		cluster[i].mobility_tnsr_sqrt.comp[k][k]=sqrt(cluster[i].mobility_tnsr.comp[k][k]);
		//		cout<<cluster[i].mobility_tnsr_sqrt.comp[k][k]<<endl;
			}
		}
		
		cluster[i].rot_mobility_tnsr_sqrt=null33D;
		for (int k=0;k<3;k++) {
			for (int l=0;l<3;l++) {
				temp(k,l)=cluster[i].rot_mobility_tnsr.comp[k][l];
		//		cout<<cluster[i].rot_mobility_tnsr.comp[k][l]<<endl;

			}
		}
		Eigen::SelfAdjointEigenSolver<MatrixXd> ROT_MOBL_MAT(temp);
		temp_sqrt = ROT_MOBL_MAT.operatorSqrt();
	//	temp_sqrt=temp.sqrt();
		for (int k=0;k<3;k++) {
			for (int l=0;l<3;l++) {
				cluster[i].rot_mobility_tnsr_sqrt.comp[k][l]=temp_sqrt(k,l);
			//	cluster[i].rot_mobility_tnsr_sqrt.comp[k][k]=sqrt(cluster[i].rot_mobility_tnsr.comp[k][k]);
		//		cout<<cluster[i].rot_mobility_tnsr_sqrt.comp[k][k]<<endl;
			}
		}
		
		cluster[i].quat={1.0,0.0,0.0,0.0};

		// update A matrix

        cluster[i].quat2rotmat();

		// modulus of A is one; A inverse is jsut the cofactor of A 

		cluster[j].linked = 1; 
		cluster[j].link_clstr_ptr= i ; 
		
		}
		
std::random_device seed;
std::mt19937 gen{seed()};
std::normal_distribution<> R1(0.0,1.0),R2(0.0,1.0),R3(0.0,1.0),R4(0.0,1.0),R5(0.0,1.0),R6(0.0,1.0);

void brownian( int step , vector<ParticleData>& cluster, vector<SubData>& particle, int *Max_Cluster_N , double *KE_rot, double vel_scale) {
double a, b , c, lambda;
vctr4D quat_old;
*KE_rot=0;
for(int i=0;i<*Max_Cluster_N;i++) 
	{
		vctr3D rand(R1(gen), R2(gen), R3(gen));
		vctr3D rand1(R4(gen), R5(gen), R6(gen));
		if (cluster[i].Sub_Length>1) 
			{
				cluster[i].pos+=cluster[i].rotmat*cluster[i].mobility_tnsr*cluster[i].rotmat*(cluster[i].frc*dt) + cluster[i].rotmat*cluster[i].mobility_tnsr_sqrt*(rand*kbT_dt);
	
				if(xx_rotation)	
				{
			// update Q
				quat_old=cluster[i].quat;
				// translate space-fixed w*dt (i.e. theta) (3 dimensions) into qdot (4 dimensions).
				// based on the Wotuer's paper on An elementary singularity-free Rotational Brownian Dynamics algorithm for anisotropic particles 
				// J. Chem. Phys. 142, 114103 (2015)
				
				cluster[i].theta   	= 	cluster[i].rot_mobility_tnsr*cluster[i].rotmat*(cluster[i].trq*dt) + cluster[i].rot_mobility_tnsr_sqrt*(rand1*kbT_dt);
				cluster[i].quat		=	cluster[i].theta2quat();
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
	std::ifstream fin("random_device_state.txt");
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
double shear_rate = 0.0; //shear rate
int ifshear = 0;// set equal to 1 for shear
std::string dataFileName="../xxx",dataFileName_new="../xxxnew" ;
int Max_Cluster_N=NrParticles;
double simu_time=dt;
int step=0, nSteps=10000, frame=1000;
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
int combine[NrParticles][2];

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
					particle[cluster[i].sub[j]].radius=0.56;					
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
			particle[cluster[i].sub[j]].radius=0.56;
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

std::ofstream outFile1(dataFileName+"/PE_energy.dat");
std::ofstream outFile10(dataFileName+"/End_positions.dat");
std::ofstream outFile11(dataFileName+"/no_of_clusters.dat");

// perfrom MD steps
/*	if (ifrestart) {
	simu_time =10.001;
	new_neighbor = TOPMAP(cellx,celly,cellz,if_Periodic,neighbor,box,(R_cut+R_shell),simu_time*shear_rate*Ly);
	std::tie(Pxx[0],Pyy[0],Pzz[0],Pxy[0],Pxz[0],Pyz[0])= forceUpdate(&p_energy ,XYZ, new_neighbor, cell, cellLength, box, (R_cut+R_shell), N, force, shear_rate, simu_time, ifshear , epsilon, sigma, rs);
}
*/

step = restart_frame_offset*frame+1;
 	cout <<cluster[0].link_clstr_ptr<<endl;

 forceUpdate( particle, &p_energy, &combine_now , combine, &step);

 	cout <<cluster[0].link_clstr_ptr<<endl;

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


simu_time =dt;
do {
	p_energy=0;	

	brownian(step, cluster, particle, &Max_Cluster_N , &KE_rot, vel_scale )	;
	combine_now=0;
 	forceUpdate( particle, &p_energy, &combine_now , combine, &step);

	if (xxclustering && combine_now>1) 
		{	
		//	cout<<combine_now<<endl;
			vector<vector<int>> temp_combine(combine_now+1,vector<int> (2)) ;
			for (int pn = 1; pn<=combine_now ; pn++) 
				{ 		
					for (int j = 0; j< 2 ; j ++) 
						{
							temp_combine[pn][j]=combine[pn][j];
						}
				//		cout<<pn<<'\t'<<temp_combine[pn][0]<<'\t'<<temp_combine[pn][1]<<"insdide main beroe sort"<<endl;
				}			
		
	sort (temp_combine.begin()+1,temp_combine.end(), RowSort());
	
		/*		for (int pn = 1; pn<=combine_now ; pn++) 
				{ 		
						cout<<temp_combine[pn][0]<<'\t'<<temp_combine[pn][1]<<"insdide main after sort"<<endl;
				}	
*/
	if(combine_now>0) {	
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
		/*		for (int pn = 1; pn<=combine_now ; pn++) 
				{ 		
						cout<<temp_combine[pn][0]<<'\t'<<temp_combine[pn][1]<<"insdide mainf after unique"<<endl;
				}*/
	// collision detection
	 		//	cout<<combine_now<<endl;

	for ( int pn = 1 ; pn <=combine_now; pn ++ )
		{
			Collision(particle, cluster, temp_combine[pn][0], temp_combine[pn][1], &Max_Cluster_N, box, rbox );
		    
		}
			 	cout <<cluster[0].link_clstr_ptr<<endl;

		
			/*			for (int pn = 1; pn<=combine_now ; pn++) 
				{ 		
						cout<<temp_combine[pn][0]<<'\t'<<temp_combine[pn][1]<<" after combine"<<endl;
				}
*/
	
		
		// std::cout<<Length_cluster[i]<<'\t'<<Length_cluster[j]<<std::endl;
		int j =1;
		do
			{ 
				if (cluster[j].link_clstr_ptr != -1 ) 
				{	
                   for ( int k = j ; k < Max_Cluster_N-1 ; k ++ )
                       {
                       cluster[k].link_clstr_ptr=cluster[k+1].link_clstr_ptr;
                       cluster[k].linked=cluster[k+1].linked;
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
                       cluster[k].sub[l]                       =       cluster[k+1].sub[l];
                       particle[cluster[k].sub[l]].cluster = particle[cluster[k+1].sub[l]].cluster ;
                       } 
					}
					Max_Cluster_N=Max_Cluster_N-1;
					if (Max_Cluster_N<=1) 
					{
						Max_Cluster_N=1;
						cout<<"all particles combined into Single Cluster "<<endl;
						abort();
					} 
				 j-=1;
				}
				j+=1;
			} while (j < Max_Cluster_N)	;
			
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

if (step%frame==0) 
	{ 

        std::ofstream outFile5(dataFileName+"/XYZ"+ std::to_string(step/frame) +".xyz");   
		outFile5<<NrParticles<<std::endl;
		outFile5<<"X Y Z co-ordinates"<<std::endl;
		outFile11<<step<<'\t'<<Max_Cluster_N<<std::endl;
		// save position, Kinetic energy, Potential energy, Forces every 'frame' steps and also store radii of gyration info
		
		std::ofstream outFile9(dataFileName+"/Cluster_dist"+ std::to_string(step/frame) +".dat");

		K_Energy=0;

		for ( int i = 0 ; i < Max_Cluster_N; i ++ )
			{
				if(cluster[i].Sub_Length>1)
				{
				outFile9<<cluster[i].radii_gyr<<'\t'<<cluster[i].Sub_Length<<std::endl;
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
		outFile1<<p_energy<<std::endl;
		outFile5.close();
		outFile9.close();
		
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
				cluster[i].mobility_tnsr.writeToFile(outFile7);
				cluster[i].mobility_tnsr_sqrt.writeToFile(outFile7);
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
				cluster[i].mobility_tnsr.writeToFile(outFile7);
				cluster[i].mobility_tnsr_sqrt.writeToFile(outFile7);
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
