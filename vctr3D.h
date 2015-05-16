#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>

using namespace std;

struct mtrx3D;

struct vctr3D
{

// data

  double comp[3];

// constructors

  vctr3D() {} 

  vctr3D( const double inp0, const double inp1, const double inp2 )
  {
    comp[0] = inp0;
    comp[1] = inp1;
    comp[2] = inp2;
  }

// constants (?)
//vctr3D null() { vctr3D result ( 0., 0., 0. );  return (result); }
  void null() { comp[0] = 0. ; comp[1] = 0. ; comp[2] = 0. ; }
//void ex  () { comp[0] = 1. ; comp[1] = 0. ; comp[2] = 0. ; }
//void ey  () { comp[0] = 0. ; comp[1] = 1. ; comp[2] = 0. ; }
//void ez  () { comp[0] = 0. ; comp[1] = 0. ; comp[2] = 1. ; }

// member functions

  void echo() // echo(stream st) to make it possible to send output to st
  {
    cout << setiosflags(ios::fixed) << setprecision(4);
    cout << setw(10) << comp[0] << "  "
         << setw(10) << comp[1] << "  "
         << setw(10) << comp[2] << endl;
  }

  double norm2()
  {
    return comp[0] * comp[0] + comp[1] * comp[1] + comp[2] * comp[2];
  }

void PBC(const vctr3D& box, const vctr3D& rbox)
  {
    comp[0] -= round( comp[0] * rbox.comp[0] ) * box.comp[0];
    comp[1] -= round( comp[1] * rbox.comp[1] ) * box.comp[1];
    comp[2] -= round( comp[2] * rbox.comp[2] ) * box.comp[2];
  } 
vctr3D& revPBC(const vctr3D& inp,const vctr3D& box, const vctr3D& rbox)
  {
	  if(abs(comp[0]-inp.comp[0])>box.comp[0]*0.5) {
		comp[0] -= box.comp[0]*(comp[0]/abs(comp[0]));
		}
	  if(abs(comp[1]-inp.comp[1])>box.comp[1]*0.5) {		
    comp[1] -= box.comp[1]*(comp[1]/abs(comp[1]));
		}
	  if(abs(comp[2]-inp.comp[2])>box.comp[2]*0.5) {		
    comp[2] -= box.comp[2]*(comp[2]/abs(comp[2]));
		}
	return (*this);
  } 
// operators

//        copy
/*vctr3D& operator = ( const vctr3D& inp )
  {
    comp[0] = inp.comp[0];
    comp[1] = inp.comp[1];
    comp[2] = inp.comp[2];
    return (*this);
  }
*/

//        addition 1
  vctr3D& operator += ( const vctr3D& inp )
  {
    comp[0] += inp.comp[0];
    comp[1] += inp.comp[1];
    comp[2] += inp.comp[2];
    return (*this);
  }

//       addition 2 // trailing 'const': value of receiving object remains unchanged
  vctr3D operator + ( const vctr3D& inp ) const
  {
    vctr3D result = *this;
    result += inp;
    return result;
  }

//        subtraction 1
  vctr3D& operator -= ( const vctr3D& inp )
  {
    comp[0] -= inp.comp[0];
    comp[1] -= inp.comp[1];
    comp[2] -= inp.comp[2];
    return (*this);
  }

//       subtraction 2
  vctr3D operator - ( const vctr3D& inp ) const
  {
    vctr3D result = *this;
    result -= inp;
    return result;
  }

//        scalar multiplication 1
  vctr3D& operator *= ( const double inp )
  {
    comp[0] *= inp;
    comp[1] *= inp;
    comp[2] *= inp;
    return (*this);
  }

//       scalar multiplication 2
  vctr3D operator * ( const double inp ) const
  {
    vctr3D result = *this;
    result *= inp;
    return result;
  }

/*       vector multiplication component-wise
  vctr3D operator & ( const vctr3D& inp ) const
  {
    vctr3D result = *this;
    result.comp[0] *= inp.comp[0];
    result.comp[1] *= inp.comp[1];
    result.comp[2] *= inp.comp[2];
    return result;
  } */

//       dot product
  double operator * ( const vctr3D& inp ) const
  {
    return ( comp[0] * inp.comp[0]
           + comp[1] * inp.comp[1]
           + comp[2] * inp.comp[2] );
  }

//       cross product
  vctr3D operator ^ ( const vctr3D& inp ) const
  {
    vctr3D result( comp[1] * inp.comp[2] - comp[2] * inp.comp[1],
                   comp[2] * inp.comp[0] - comp[0] * inp.comp[2],
                   comp[0] * inp.comp[1] - comp[1] * inp.comp[0] );
    return result;
  }

}; // vctr3D


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


struct mtrx3D
{
  double comp[3][3];

// constructors

  mtrx3D() {};

//         diagonal matrix 1
  mtrx3D ( const vctr3D& inp )
  {
    comp[0][0] = inp.comp[0];
    comp[0][1] = 0.;
    comp[0][2] = 0.;
    comp[1][0] = 0.;
    comp[1][1] = inp.comp[1];
    comp[1][2] = 0.;
    comp[2][0] = 0.;
    comp[2][1] = 0.;
    comp[2][2] = inp.comp[2];
  }

//         diagonal matrix 2
  mtrx3D ( const double inp0, const double inp1, const double inp2 )
  {
    comp[0][0] = inp0;
    comp[0][1] = 0.;
    comp[0][2] = 0.;
    comp[1][0] = 0.;
    comp[1][1] = inp1;
    comp[1][2] = 0.;
    comp[2][0] = 0.;
    comp[2][1] = 0.;
    comp[2][2] = inp2;
  }


//         dyadic product 1
  mtrx3D ( const vctr3D& inp0, const vctr3D& inp1 )
  {
    comp[0][0] = inp0.comp[0] * inp1.comp[0];
    comp[0][1] = inp0.comp[0] * inp1.comp[1];
    comp[0][2] = inp0.comp[0] * inp1.comp[2];
    comp[1][0] = inp0.comp[1] * inp1.comp[0];
    comp[1][1] = inp0.comp[1] * inp1.comp[1];
    comp[1][2] = inp0.comp[1] * inp1.comp[2];
    comp[2][0] = inp0.comp[2] * inp1.comp[0];
    comp[2][1] = inp0.comp[2] * inp1.comp[1];
    comp[2][2] = inp0.comp[2] * inp1.comp[2];
  }

//         full definition, vectors = columns
  mtrx3D ( const vctr3D& inp0, const vctr3D& inp1, const vctr3D& inp2 )
  {
    comp[0][0] = inp0.comp[0];
    comp[0][1] = inp1.comp[0];
    comp[0][2] = inp2.comp[0];
    comp[1][0] = inp0.comp[1];
    comp[1][1] = inp1.comp[1];
    comp[1][2] = inp2.comp[1];
    comp[2][0] = inp0.comp[2];
    comp[2][1] = inp1.comp[2];
    comp[2][2] = inp2.comp[2];
  }

// constants (?)
  mtrx3D null() { vctr3D v ( 0., 0., 0. ); mtrx3D result(v); return(result); }
  mtrx3D one () { vctr3D v ( 1., 1., 1. ); mtrx3D result(v); return(result); }

// member functions

  void echo()
  {
    cout << setw(10) << comp[0][0] << "  " << setw(10) << comp[0][1] << "  " << setw(10) << comp[0][2] << endl;
    cout << setw(10) << comp[1][0] << "  " << setw(10) << comp[1][1] << "  " << setw(10) << comp[1][2] << endl;
    cout << setw(10) << comp[2][0] << "  " << setw(10) << comp[2][1] << "  " << setw(10) << comp[2][2] << endl;
  }

  double trace() const
  {
    return ( comp[0][0] + comp[1][1] + comp[2][2] );
  }

  double det()
  {
    return ( comp[0][0] * ( comp[1][1] * comp[2][2] - comp[1][2] * comp[2][1] )
           + comp[0][1] * ( comp[1][2] * comp[2][0] - comp[1][0] * comp[2][2] )
           + comp[0][2] * ( comp[1][0] * comp[1][1] - comp[2][1] * comp[2][0] )
           );
  }

// operators

//        matrix addition 1
  mtrx3D& operator += ( const mtrx3D& inp )
  {
    comp[0][0] += inp.comp[0][0];
    comp[0][1] += inp.comp[0][1];
    comp[0][2] += inp.comp[0][2];
    comp[1][0] += inp.comp[1][0];
    comp[1][1] += inp.comp[1][1];
    comp[1][2] += inp.comp[1][2];
    comp[2][0] += inp.comp[2][0];
    comp[2][1] += inp.comp[2][1];
    comp[2][2] += inp.comp[2][2];
    return (*this);
  }

//       matrix addition 2
  mtrx3D operator + ( const mtrx3D& inp ) const
  {
    mtrx3D result = *this;
    result += inp;
    return result;
  }

//        matrix subtraction 1
  mtrx3D& operator -= ( const mtrx3D& inp )
  {
    comp[0][0] -= inp.comp[0][0];
    comp[0][1] -= inp.comp[0][1];
    comp[0][2] -= inp.comp[0][2];
    comp[1][0] -= inp.comp[1][0];
    comp[1][1] -= inp.comp[1][1];
    comp[1][2] -= inp.comp[1][2];
    comp[2][0] -= inp.comp[2][0];
    comp[2][1] -= inp.comp[2][1];
    comp[2][2] -= inp.comp[2][2];
    return (*this);
  }

//       matrix subtraction 2
  mtrx3D operator - ( const mtrx3D& inp ) const
  {
    mtrx3D result = *this;
    result -= inp;
    return result;
  }

//       matrix - scalar multiplication 1
  mtrx3D& operator *= ( const double inp )
  {
    comp[0][0] *= inp;
    comp[0][1] *= inp;
    comp[0][2] *= inp;
    comp[1][0] *= inp;
    comp[1][1] *= inp;
    comp[1][2] *= inp;
    comp[2][0] *= inp;
    comp[2][1] *= inp;
    comp[2][2] *= inp;
    return (*this);
  }

//       matrix - scalar multiplication 2
  mtrx3D operator * ( const double inp ) const
  {
    mtrx3D result = *this;
    result *= inp;
    return result;
  }

//       matrix - vector multiplication
  vctr3D operator * ( const vctr3D& inp ) const
  {
    vctr3D result;
    result.comp[0] = comp[0][0] * inp.comp[0] + comp[0][1] * inp.comp[1] + comp[0][2] * inp.comp[2];
    result.comp[1] = comp[1][0] * inp.comp[0] + comp[1][1] * inp.comp[1] + comp[1][2] * inp.comp[2];
    result.comp[2] = comp[2][0] * inp.comp[0] + comp[2][1] * inp.comp[1] + comp[2][2] * inp.comp[2];
    return result;
  }

//        matrix - matrix multiplication 1
  mtrx3D& operator *= ( const mtrx3D& inp )
  {
    mtrx3D tmp = *this;
    comp[0][0] = tmp.comp[0][0] * inp.comp[0][0] + tmp.comp[0][1] * inp.comp[1][0] + tmp.comp[0][2] * inp.comp[2][0];
    comp[0][1] = tmp.comp[0][0] * inp.comp[0][1] + tmp.comp[0][1] * inp.comp[1][1] + tmp.comp[0][2] * inp.comp[2][1];
    comp[0][2] = tmp.comp[0][0] * inp.comp[0][2] + tmp.comp[0][1] * inp.comp[1][2] + tmp.comp[0][2] * inp.comp[2][2];
    comp[1][0] = tmp.comp[1][0] * inp.comp[0][0] + tmp.comp[1][1] * inp.comp[1][0] + tmp.comp[1][2] * inp.comp[2][0];
    comp[1][1] = tmp.comp[1][0] * inp.comp[0][1] + tmp.comp[1][1] * inp.comp[1][1] + tmp.comp[1][2] * inp.comp[2][1];
    comp[1][2] = tmp.comp[1][0] * inp.comp[0][2] + tmp.comp[1][1] * inp.comp[1][2] + tmp.comp[1][2] * inp.comp[2][2];
    comp[2][0] = tmp.comp[2][0] * inp.comp[0][0] + tmp.comp[2][1] * inp.comp[1][0] + tmp.comp[2][2] * inp.comp[2][0];
    comp[2][1] = tmp.comp[2][0] * inp.comp[0][1] + tmp.comp[2][1] * inp.comp[1][1] + tmp.comp[2][2] * inp.comp[2][1];
    comp[2][2] = tmp.comp[2][0] * inp.comp[0][2] + tmp.comp[2][1] * inp.comp[1][2] + tmp.comp[2][2] * inp.comp[2][2];
    return (*this);
  }

//       matrix - matrix multiplication 2
  mtrx3D operator * ( const mtrx3D& inp ) const
  {
    mtrx3D result = *this;
    result *= inp;
    return result;
  }

//       transpose
  mtrx3D operator ~ () const
  {
    mtrx3D tmp = *this;
    mtrx3D result;
    result.comp[0][0] = tmp.comp[0][0];
    result.comp[0][1] = tmp.comp[1][0];
    result.comp[0][2] = tmp.comp[2][0];
    result.comp[1][0] = tmp.comp[0][1];
    result.comp[1][1] = tmp.comp[1][1];
    result.comp[1][2] = tmp.comp[2][1];
    result.comp[2][0] = tmp.comp[0][2];
    result.comp[2][1] = tmp.comp[1][2];
    result.comp[2][2] = tmp.comp[2][2];
    return result;
  }

//       dyadic product 2
//friend mtrx3D operator & ( const vctr3D& inp0, const vctr3D& inp1 );

   mtrx3D dyadic ( const vctr3D& inp0, const vctr3D& inp1 )
   {
     mtrx3D result( inp0, inp1 );
     return result;
   }

   mtrx3D inverse ()
   {
    mtrx3D result;
    double rdet;
    rdet = 1. / ( comp[0][0] * ( comp[1][1] * comp[2][2] - comp[1][2] * comp[2][1] )
                + comp[0][1] * ( comp[1][2] * comp[2][0] - comp[1][0] * comp[2][2] )
                + comp[0][2] * ( comp[1][0] * comp[2][1] - comp[1][1] * comp[2][0] ) );
    result.comp[0][0] =   ( comp[1][1] * comp[2][2] - comp[1][2] * comp[2][1] ) * rdet;
    result.comp[0][1] = - ( comp[0][1] * comp[2][2] - comp[0][2] * comp[2][1] ) * rdet;
    result.comp[0][2] =   ( comp[0][1] * comp[1][2] - comp[0][2] * comp[1][1] ) * rdet;
    result.comp[1][0] = - ( comp[1][0] * comp[2][2] - comp[1][2] * comp[2][0] ) * rdet;
    result.comp[1][1] =   ( comp[0][0] * comp[2][2] - comp[0][2] * comp[2][0] ) * rdet;
    result.comp[1][2] = - ( comp[0][0] * comp[1][2] - comp[0][2] * comp[1][0] ) * rdet;
    result.comp[2][0] =   ( comp[1][0] * comp[2][1] - comp[1][1] * comp[2][0] ) * rdet;
    result.comp[2][1] = - ( comp[0][0] * comp[2][1] - comp[0][1] * comp[2][0] ) * rdet;
    result.comp[2][2] =   ( comp[0][0] * comp[1][1] - comp[0][1] * comp[1][0] ) * rdet;
    return result;
   }

}; // mtrx3D

/*       dyadic product 2
  mtrx3D operator & ( const vctr3D& inp0, const vctr3D& inp1 )
  {
    mtrx3D result( inp0, inp1 );
    return result;
  }
*/

