#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "vctr3D.h"
#include "vctr4D.h"

using namespace std;

// useful constants

int const x = 0;
int const y = 1;
int const z = 2;

int const ang2D = z;  // dimension used to store rotation angle in 2D simulation

double const pi     = 3.1415926535;
double const twopi  = 2.  * pi;
double const halfpi = 0.5 * pi;
double const recipi = 1. / pi;

vctr3D const null3D  ( 0., 0., 0. );        // why can I not use these in calculations?
mtrx3D const null33D ( 0., 0., 0. );
vctr4D const null4D  ( 0., 0., 0., 0. );

//std::string const in_bld      ("input.build");  // "proper" syntax
string const in_conf   = "input.Config";
string const in_ctrl   = "input.control";
string const in_fld    = "input.field";
string const out_dcd   = "output.dcd";
string const out_conf  = "output.Config";
string const in_bld    = "input.build";
string const out_bld   = "output.build";

//======================================================================

struct WallData
{
  bool   on;
  double mean;
  double ampl;
  double freq;
  double phase, dphse;
  double pos;
  int    dirctn; // = box.dimens - 1 ( if on )

  WallData () // default settings
  {
    on     = false;
    mean   = -10.;
    pos    = mean;
    ampl   =   0.;
    freq   =   0.;
    phase  =   0.;
    dphse  =   0.;
//  dirctn =   z;
  }
}; // WallData

//======================================================================

struct BoxData
{
  int       file;
  int       type;
  int       dimens;
  vctr3D    edge;
  vctr3D    halve;
  vctr3D    recip;
  WallData  wall;
  int       angle_start, angle_stop;
  int       omega_start, omega_stop;
  double    BckgrndFric;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void update() // update box-size related parameters
  {
    for ( int i = 0 ; i < 3 ; i ++ )
    {
      recip.comp[i] = 1.0 / edge.comp[i];
      halve.comp[i] = 0.5 * edge.comp[i];
    }
  } // update
  
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void SetAngRange() // range of angle-paramaters depending on box dimensionaliy
  {
    switch ( dimens )
    {
      case 2  : if ( ang2D > 2 ) { cout << "*** Fatal ang2D: " << ang2D << endl; abort(); }
                angle_start = ang2D;  // one angle
                angle_stop  = ang2D;
                omega_start = ang2D;  // one angular velocity
                omega_stop  = ang2D;
                break;
      case 3  : angle_start = 0;
                angle_stop  = 3;  // four quaternions
                omega_start = 0;
                omega_stop  = 2;  // three Cartesian (?) angular velocities
                break;
      default : cout << "*** incorrect number of dimensions" << endl; abort();
    }
  } // SetAngRange
  
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void PBC( vctr3D& inp )
  {
    inp.comp[0] -= round( inp.comp[0] * recip.comp[0] ) * edge.comp[0];
    inp.comp[1] -= round( inp.comp[1] * recip.comp[1] ) * edge.comp[1];
    inp.comp[2] -= round( inp.comp[2] * recip.comp[2] ) * edge.comp[2];
  }
}; // BoxData

//======================================================================

struct SubData
{
  vctr3D pos;
  vctr3D pos_bdyfxd;
  vctr3D cell;
  vctr3D vel;
  vctr3D frc;
  double mass;
  int PartNum;
  double radius;
  int cluster;

/*double sigma;
  double eps;
  double mass; */
}; // SubData

//======================================================================

struct ParticleData
{
  vctr3D pos;
  vctr3D vel;
  vctr3D frc;
  vector <int> sub;
  vctr3D omega;         // lab-fixed
  mtrx3D rotmat;        // apply for body-fixed-to-lab-fixed
  mtrx3D Iner_tnsr;
  mtrx3D mobility_tnsr;
  mtrx3D mobility_tnsr_sqrt; 
  mtrx3D rot_mobility_tnsr;
  mtrx3D rot_mobility_tnsr_sqrt;
  vctr3D trq;           // lab-fixed = cartesian
  vctr3D theta;
  vctr3D angmom;        // lab-fixed
  vctr4D quat;
  int Sub_Length;
  int clicked; 
  double radii;
  double radii_gyr;
  double mass;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


  ParticleData( const int NrSubs ) // constructor
  {
    sub.resize(NrSubs);
	mass;
    pos    = null3D;
    vel    = null3D;
    angmom = null3D;
    omega  = null3D;
    quat   = null4D;
    Iner_tnsr = null33D;
	Sub_Length=0;

  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void echo( const BoxData& box )
  {
    int i;
    cout << setiosflags(ios::fixed) << setprecision(4);
    for ( i =  0 ; i < box.dimens ; i ++ ) { cout << "  " << setw(10) << pos.comp[i]; }
    for ( i =  0 ; i < box.dimens ; i ++ ) { cout << "  " << setw(10) << vel.comp[i]; }
    for ( i = box.angle_start ; i <= box.angle_stop ; i ++ )
                                           { cout << "  " << setw(10) << quat.comp[i]; }
    for ( i = box.omega_start ; i <= box.omega_stop ; i ++ )
                                           { cout << "  " << setw(10) << omega.comp[i]; }
    cout << endl;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void phi2rotmat() // for 2D box: conversion from rotation angle to rotation matrix
  {
    double cosphi,sinphi;
    cosphi = cos(quat.comp[ang2D]);
    sinphi = sin(quat.comp[ang2D]);
    rotmat.comp[0][0] =   cosphi;
    rotmat.comp[0][1] = - sinphi;
    rotmat.comp[1][0] =   sinphi;
    rotmat.comp[1][1] =   cosphi;
    rotmat.comp[0][2] =       0.;
    rotmat.comp[1][2] =       0.;
    rotmat.comp[2][0] =       0.;
    rotmat.comp[2][1] =       0.;
    rotmat.comp[2][2] =       1.;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void quat2rotmat() // for 3D box: conversion from quaternions to rotation matrix
  {
    double q2[4],q3[4][4];
    q2[1]    = quat.comp[1] * quat.comp[1];
    q2[2]    = quat.comp[2] * quat.comp[2];
    q2[3]    = quat.comp[3] * quat.comp[3];
    q2[0]    = quat.comp[0] * quat.comp[0] - q2[1] - q2[2] - q2[3];
    q3[0][1] = quat.comp[0] * quat.comp[1] * 2.;
    q3[0][2] = quat.comp[0] * quat.comp[2] * 2.;
    q3[0][3] = quat.comp[0] * quat.comp[3] * 2.;
    q3[1][2] = quat.comp[1] * quat.comp[2] * 2.;
    q3[1][3] = quat.comp[1] * quat.comp[3] * 2.;
    q3[2][3] = quat.comp[2] * quat.comp[3] * 2.;
    rotmat.comp[0][0] = q2[0]    + q2[1]   * 2.;  // = q0^2 +   q1^2 -   q2^2 -   q3^2
                                                  // = 1             - 2 q2^2 - 2 q3^2
    rotmat.comp[1][1] = q2[0]    + q2[2]   * 2.;  // = q0^2 -   q1^2 +   q2^2 -   q3^2
                                                  // = 1    - 2 q1^2          - 2 q3^2
    rotmat.comp[2][2] = q2[0]    + q2[3]   * 2.;  // = qo^2 -   q1^2 -   q2^2 +   q3^2
                                                  // = 1    - 2 q1^2 - 2 q2^2
    rotmat.comp[0][1] = q3[1][2] - q3[0][3];      // = 2 ( q1 * q2 - q0 * q3 )
    rotmat.comp[0][2] = q3[1][3] + q3[0][2];      // = 2 ( q1 * q3 + q0 * q2 )
    rotmat.comp[1][0] = q3[1][2] + q3[0][3];      // note: inverse of A in Goldstein (4.67) / (4.115)
    rotmat.comp[1][2] = q3[2][3] - q3[0][1];
    rotmat.comp[2][0] = q3[1][3] - q3[0][2];
    rotmat.comp[2][1] = q3[2][3] + q3[0][1];

/*  for testing puroposes:
    double dot;
    int    i,j;
    bool   halt = false; 
    mtrx3D maal;
    dot = quat * quat;
    if ( abs( 1. - dot ) > 1.e-6 ) { cout << "*** quaternion wrong length *** " << dot << endl; halt = true; }
    for ( i = 0 ; i < 3 ; i++ )
    {
      for ( j = i; j < 3 ; j++ )
      {
        dot = rotmat.comp[i][0] * rotmat.comp[j][0] + rotmat.comp[i][1] * rotmat.comp[j][1] + rotmat.comp[i][2] * rotmat.comp[j][2];
        if ( i == j )
        {
          if ( abs( 1. - dot ) > 1.e-6 ) { cout << "*** rotation matrix not orthonormal *** " << dot << endl; halt = true; }
        }
        else // i != j
        {
          if ( abs( 0. - dot ) > 1.e-6 ) { cout << "*** rotation matrix not orthogonal  *** " << dot << endl; halt = true; }
        }
      } // j
    } // i
    if ( halt )
    {
      maal = ( ~ rotmat ) * rotmat - mtrx3D( 1., 1., 1.);
      quat.echo();   cout << endl;
      rotmat.echo(); cout << endl;
      maal.echo();   cout << endl; abort();
    }
*/
  } // quat2rotmat()

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  vctr4D omega2qdot() // for 3D box: translate space-fixed omega (3 dimensions) into qdot (4 dimensions).
  {
/*
       (q1)      ( q0 -q1 -q2 -q3 ) (   0  )
    d  (q2)    1 ( q1  q0  q3 -q2 ) ( ws_x )
    -- (q3)  = - ( q2 -q3  q0  q1 ) ( ws_y )  = C Ws
    dt (q4)    2 ( q3  q2 -q1  q0 ) ( ws_z )
*/
    vctr4D result;
    result.comp[0] = - quat.comp[1] * omega.comp[x] - quat.comp[2] * omega.comp[y] - quat.comp[3] * omega.comp[z]; // + quat.comp[0] * 0.
    result.comp[1] =   quat.comp[0] * omega.comp[x] + quat.comp[3] * omega.comp[y] - quat.comp[2] * omega.comp[z]; // + quat.comp[1] * 0.
    result.comp[2] = - quat.comp[3] * omega.comp[x] + quat.comp[0] * omega.comp[y] + quat.comp[1] * omega.comp[z]; // + quat.comp[2] * 0.
    result.comp[3] =   quat.comp[2] * omega.comp[x] - quat.comp[1] * omega.comp[y] + quat.comp[0] * omega.comp[z]; // + quat.comp[3] * 0.
    result *= 0.5;
    return result;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  vctr4D theta2quat() // for 3D box: translate body-fixed w*dt (i.e. theta) (3 dimensions) into qdot (4 dimensions).
					  // based on the Wotuer's paper on An elementary singularity-free Rotational Brownian Dynamics algorithm for anisotropic particles 
					  // J. Chem. Phys. 142, 114103 (2015)
  {
/*
			  ( q0 -q1 -q2 -q3 ) 
			1 ( q1  q0 -q3  q2 ) 
    B =     - ( q2  q3  q0 -q1 )   
            2 ( q3 -q2  q1  q0 ) 
*/
    vctr4D result;
    result.comp[0] = - quat.comp[1] * theta.comp[x] - quat.comp[2] * theta.comp[y] - quat.comp[3] * theta.comp[z]; // + quat.comp[0] * 0.
    result.comp[1] =   quat.comp[0] * theta.comp[x] - quat.comp[3] * theta.comp[y] + quat.comp[2] * theta.comp[z]; // + quat.comp[1] * 0.
    result.comp[2] =   quat.comp[3] * theta.comp[x] + quat.comp[0] * theta.comp[y] - quat.comp[1] * theta.comp[z]; // + quat.comp[2] * 0.
    result.comp[3] = - quat.comp[2] * theta.comp[x] + quat.comp[1] * theta.comp[y] + quat.comp[0] * theta.comp[z]; // + quat.comp[3] * 0.
    result *= 0.5;
    return result;
  }
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  vctr4D theta_space2quat(const vctr3D& inp) // for 3D box: translate space-fixed w*dt (i.e. theta) (3 dimensions) into qdot (4 dimensions).
					  // based on the Wotuer's paper on An elementary singularity-free Rotational Brownian Dynamics algorithm for anisotropic particles 
					  // J. Chem. Phys. 142, 114103 (2015)
  {
/*
			  ( q0 -q1 -q2 -q3 ) 
			1 ( q1  q0  q3 -q2 ) 
    C  =    - ( q2 -q3  q0  q1 )   
            2 ( q3  q2 -q1  q0 ) 
*/
    vctr4D result;
    result.comp[0] = - quat.comp[1] * inp.comp[x] - quat.comp[2] * inp.comp[y] - quat.comp[3] * inp.comp[z]; // + quat.comp[0] * 0.
    result.comp[1] =   quat.comp[0] * inp.comp[x] + quat.comp[3] * inp.comp[y] - quat.comp[2] * inp.comp[z]; // + quat.comp[1] * 0.
    result.comp[2] = - quat.comp[3] * inp.comp[x] + quat.comp[0] * inp.comp[y] + quat.comp[1] * inp.comp[z]; // + quat.comp[2] * 0.
    result.comp[3] =   quat.comp[2] * inp.comp[x] - quat.comp[1] * inp.comp[y] + quat.comp[0] * inp.comp[z]; // + quat.comp[3] * 0.
    result *= 0.5;
    return result;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void qdot2omega( const vctr4D& inp, double& chck ) // for 3D box: translate qdot (4 dimensions) into omega (3 dimensions).
  {
/*
    (   0  )     (  q0  q1  q2  q3 ) ( q0dot )
    ( wb_x )     ( -q1  q0  q3 -q2 ) ( q1dot )     -1      T
    ( wb_y ) = 2 ( -q2 -q3  q0  q1 ) ( q2dot )    B   = 4 B
    ( wb_z )     ( -q3  q2 -q1  q0 ) ( q3dot )
*/
    chck           = (   quat.comp[0] * inp.comp[0] + quat.comp[1] * inp.comp[1] + quat.comp[2] * inp.comp[2] + quat.comp[3] * inp.comp[3]) * 2.;
    omega.comp[x]  = ( - quat.comp[1] * inp.comp[0] + quat.comp[0] * inp.comp[1] + quat.comp[3] * inp.comp[2] - quat.comp[2] * inp.comp[3]) * 2.;
    omega.comp[y]  = ( - quat.comp[2] * inp.comp[0] - quat.comp[3] * inp.comp[1] + quat.comp[0] * inp.comp[2] + quat.comp[1] * inp.comp[3]) * 2.;
    omega.comp[z]  = ( - quat.comp[3] * inp.comp[0] + quat.comp[2] * inp.comp[1] - quat.comp[1] * inp.comp[2] + quat.comp[0] * inp.comp[3]) * 2.;
//  cout << "check " << check << endl;
  }

}; // ParticleData

/*
struct PairData
{
  double sigma2;
  double eps4;
};
*/

struct EneData
{
  double tot;
  double kin, trns, rot;
  double pot, fld, wall, grav;
  int    pairs;

  EneData ()
  {
    wall = 0.;
    grav = 0.;
  }
}; // EneData



struct BodyFixedData
{	vctr3D pos;
	vector <SubData> sub;


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  BodyFixedData( const int NrSubs ) // constructor
  {
    sub.resize(NrSubs);

    pos    = null3D;

    for ( int si = 0 ; si < NrSubs ; si ++ )
      { sub[si].pos = null3D;} // sub[si].vel = null3D; }
  }
};

/*
struct ShapeData
{
  double length;
  double h_length;
  double diamtr;
  double radius;
  double TotLngth;
};
*/

struct CtrlData
{
  int       NrSteps;
  double    dt, rdt;
//double    tmprtr;
  int       logintrvl, lognext;
  int       dcdintrvl, dcdnext;
  WallData  wall;
  double    grav;
  double    BckgrndFric;
};
