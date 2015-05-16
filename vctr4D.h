#include <cstdio>
#include <cstdlib>
#include <iostream>

using namespace std;

struct vctr4D
{

// data

  double comp[4];

// constructors

  vctr4D() {} 

  vctr4D( const double inp0, const double inp1, const double inp2, const double inp3 )
  {
    comp[0] = inp0;
    comp[1] = inp1;
    comp[2] = inp2;
    comp[3] = inp3;
  }

// constants (?)
  void null() { comp[0] = 0. ; comp[1] = 0. ; comp[2] = 0. ; comp[3] = 0. ; }

// member functions

  void echo()
  {
    cout << setiosflags(ios::fixed) << setprecision(4);
    cout << setw(10) << comp[0] << "  "
         << setw(10) << comp[1] << "  "
         << setw(10) << comp[2] << "  "
         << setw(10) << comp[3] << endl;
  }

  double norm2()
  {
    return comp[0] * comp[0] + comp[1] * comp[1] + comp[2] * comp[2] + comp[3] * comp[3];
  }

// operators

/*        copy
  vctr4D& operator = ( const vctr4D& inp )
  {
    comp[0] = inp.comp[0];
    comp[1] = inp.comp[1];
    comp[2] = inp.comp[2];
    comp[3] = inp.comp[2];
    return (*this);
  }
*/

//        addition 1
  vctr4D& operator += ( const vctr4D& inp )
  {
    comp[0] += inp.comp[0];
    comp[1] += inp.comp[1];
    comp[2] += inp.comp[2];
    comp[3] += inp.comp[3];
    return (*this);
  }

//       addition 2
  vctr4D operator + ( const vctr4D& inp ) const
  {
    vctr4D result = *this;
    result += inp;
    return result;
  }

//        subtraction 1
  vctr4D& operator -= ( const vctr4D& inp )
  {
    comp[0] -= inp.comp[0];
    comp[1] -= inp.comp[1];
    comp[2] -= inp.comp[2];
    comp[3] -= inp.comp[3];
    return (*this);
  }

//       subtraction 2
  vctr4D operator - ( const vctr4D& inp ) const
  {
    vctr4D result = *this;
    result -= inp;
    return result;
  }

//        scalar multiplication 1
  vctr4D& operator *= ( double inp )
  {
    comp[0] *= inp;
    comp[1] *= inp;
    comp[2] *= inp;
    comp[3] *= inp;
    return (*this);
  }

//       scalar multiplication 2
  vctr4D operator * ( const double inp ) const
  {
    vctr4D result = *this;
    result *= inp;
    return result;
  }

//        scalar division 1
  vctr4D& operator /= ( const double inp )
  {
    comp[0] /= inp;
    comp[1] /= inp;
    comp[2] /= inp;
    comp[3] /= inp;
    return (*this);
  }

//       scalar division 2
  vctr4D operator / ( const double inp ) const
  {
    vctr4D result = *this;
    result *= inp;
    return result;
  }

//       dot product
  double operator * ( vctr4D& inp ) const
  {
    return (  comp[0] * inp.comp[0]
            + comp[1] * inp.comp[1]
            + comp[2] * inp.comp[2]
            + comp[3] * inp.comp[3] );
  }

}; // vctr4D

struct mtrx4D
{
  double comp[4][4];

// constructors

// diagonal matrix 1
  mtrx4D ( const vctr4D& inp )
  {
    for ( int i = 0 ; i <= 3 ; i ++ )
    {
      for ( int j = 0 ; j <= 3 ; j ++ )
      {
        comp[i][j] = 0.;
      }
      comp[i][i] = inp.comp[i];
    }
  }

// diagonal matrix 2
  mtrx4D ( const double inp0, const double inp1, const double inp2, const double inp3 )
  {
    for ( int i = 0 ; i < 4 ; i ++ )
    {
      for ( int j = 0 ; j < 4 ; j ++ )
      {
        comp[i][j] = 0.;
      }
    }
    comp[0][0] = inp0;
    comp[1][1] = inp1;
    comp[2][2] = inp2;
    comp[3][3] = inp3;
  }

// full definition, columns = vectors
  mtrx4D ( const vctr4D& inp0, const vctr4D& inp1, const vctr4D& inp2, const vctr4D& inp3 )
  {
    for ( int i = 0 ; i < 4 ; i ++ )
    {
      comp[i][0] = inp0.comp[i];
      comp[i][1] = inp1.comp[i];
      comp[i][2] = inp2.comp[i];
      comp[i][3] = inp3.comp[i];
    }
  }

}; // mtrx4D
