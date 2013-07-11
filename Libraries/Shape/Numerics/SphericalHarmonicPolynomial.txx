/*=========================================================================

  Author: Christine Xu, Ipek Oguz, Martin Styner

=========================================================================*/
#ifndef _namicSphericalHarmonicPolynomial_txx
#define _namicSphericalHarmonicPolynomial_txx

#include "SphericalHarmonicPolynomial.h"
#include <ctime>
#include <cstdlib>

#include <math.h>
#include <iostream>
#include <fstream>

namespace neurolib
{
#ifndef M_PI
#define M_PI 3.1415926535897932
#endif
#ifndef M_PI_2
#define M_PI_2 1.5707963267948966
#endif

template <unsigned int TDimension>
SphericalHarmonicPolynomial<TDimension>::SphericalHarmonicPolynomial()
{
  // m_Dimension = 2;
  srand( time( 0 ) );
}

template <unsigned int TDimension>
SphericalHarmonicPolynomial<TDimension>::~SphericalHarmonicPolynomial()
{

}

template <unsigned int TDimension>
void SphericalHarmonicPolynomial<TDimension>::SetCoefs(CoefListType& coeflist)
{
  if( coeflist.empty() )
    {
    throw SphericalHarmonicPolynomialException(__FILE__, __LINE__, "The input coefficient list is empty!");
    }
  unsigned int d = coeflist[0].GetPointDimension();
  if( !(TDimension == d) )
    {
    throw SphericalHarmonicPolynomialException(
            __FILE__, __LINE__, "The dimension of the input coefficient list does not match the specified dimension!");
    }
  m_Coefs = coeflist;
}

template <unsigned int TDimension>
void SphericalHarmonicPolynomial<TDimension>::Evaluate(unsigned int from_l,
                                                       unsigned int to_l, double phi, double theta,
                                                       double* sum)
{
  // if(m_Dimension != 2)
  //  throw SphericalHarmonicPolynomialException(__FILE__, __LINE__, "Only Spherical Harmonics of dimension two is
  // currently implemented.");
  // if(size > TDimension)
  //   throw SphericalHarmonicPolynomialException(__FILE__, __LINE__, "Dimension of the output sum exceeds the dimension
  // of the coefficients.");
  if( m_Coefs.empty() )
    {
    throw SphericalHarmonicPolynomialException(__FILE__, __LINE__, "Coefficients must be specified.");
    }
  if( from_l > to_l )
    {
    throw SphericalHarmonicPolynomialException(__FILE__, __LINE__,
                                               "The starting degree should be smaller or equal to the ending degree.");
    }
  if( to_l > m_Degree )
    {
    throw SphericalHarmonicPolynomialException(__FILE__, __LINE__,
                                               "The evalueated degree mustn't exceed the size of the coefficients.");
    }
  if( m_Coefs.size() < (to_l + 1) * (to_l + 1) )
    {
    throw SphericalHarmonicPolynomialException(__FILE__, __LINE__, "Coefficients size mismatch.");
    }

  unsigned int l, m, coord;
  for( coord = 0; coord < TDimension; coord++ )
    {
    sum[coord] = 0.0;
    for( l = from_l; l <= to_l; l++ ) // for the from_l degree to to_l degree
      {
      double* plm = new double[l + 1];
      for( m = 0; m <= l; m++ )
        {
        double coef = sqrt( (2 * l + 1.0) / 4.0 / M_PI / fac_quot(l + m, l - m) );
        double row = plgndr_row(l, m, cos(theta) );
        plm[m] = coef * row;
        }

      sum[coord] += m_Coefs[l * l][coord] * plm[0];
      for( m = 1; m <= l; m++ )
        {
        double sin_t = sin(m * phi);
        double cos_t = cos(m * phi);
        sum[coord] += ( cos_t * m_Coefs[l * l + 2 * m - 1][coord] // real part
                        + sin_t * m_Coefs[l * l + 2 * m][coord] ) // imaginary part
          * plm[m];
        }
      delete [] plm;
      }
    }
}

/*
output:
p_ptr:  the value of specified list of Associated Legendre polynomials (vertically, columnly)

input:
m:      the order of the polynomials
l:      the maximal degree of the polynomials
x:      the input parameter of the Associated Legendre polynomials

Legendre Polynomials come from the Sturm-Liouville Boundary Value Problem:
(1-x^2)y'' - 2xy' + n(n+1)y = 0

The Legendre Polynomials have the recurrence relation of:
(n+1)P_n+1(x) = (2n+1)xP_n(x) - n P_n-1(x)

The associate Legendre Polynomials are defined using Legendre Polynomials,
and they also obey the following recurrence relations
(l-m)P^m_l(x) = x(2l-1)P^m_l-1(x) - (l+m-1)P^m_l-2(x)

More about Legendre Polynomials:
http://mathworld.wolfram.com/LegendrePolynomial.html
*/
template <unsigned int TDimension>
double SphericalHarmonicPolynomial<TDimension>::plgndr_row(int l, int m, double x)
{
  double fact, pll, pmm, pmmp1, somx2;
  int    i, ll;

  pll = 0;
  if( m < 0 || m > l || fabs(x) > 1.0 )
    {
    std::cout << "Bad arguments in routine PLGNDR" << std::endl;
    exit(1);
    }
  pmm = 1.0;
  if( m > 0 )
    {
    somx2 = sqrt( (1.0 - x) * (1.0 + x) );
    fact = 1.0;
    for( i = 1; i <= m; i++ )
      {
      pmm *= -fact * somx2;
      fact += 2.0;
      }
    }
  // *p_ptr++= pmm;
  if( l == m )
    {
    return pmm;
    }
  else
    {
    pmmp1 = x * (2 * m + 1) * pmm;
    // *p_ptr++= pmmp1;
    if( l == (m + 1) )
      {
      return pmmp1;
      }
    else
      {
      for( ll = (m + 2); ll <= l; ll++ )
        {
        pll = (x * (2 * ll - 1) * pmmp1 - (ll + m - 1) * pmm) / (ll - m);
        // *p_ptr++= pll;
        pmm = pmmp1;
        pmmp1 = pll;
        }
      return pll;
      }
    }
}

/* factorial quotient a!/b!*/
template <unsigned int TDimension>
double SphericalHarmonicPolynomial<TDimension>::fac_quot(int a, int b)
{
  double res = 1;

  while( a > b )
    {
    res *= a--;
    }

  return res;
}

template <unsigned int TDimension>
void SphericalHarmonicPolynomial<TDimension>::EvaluateAt(unsigned int l, unsigned int m, double phi, double theta,
                                                         double* result)
{
  if( m_Coefs.empty() )
    {
    throw SphericalHarmonicPolynomialException(__FILE__, __LINE__, "Coefficients must be specified.");
    }
  if( l > m_Degree )
    {
    throw SphericalHarmonicPolynomialException(__FILE__, __LINE__,
                                               "The evalueated degree mustn't exceed the size of the coefficients.");
    }

  if( l < m )
    {
    result[0] = result[1] = 0;
    return;
    }

  double coef = sqrt( (2 * l + 1.0) / 4.0 / M_PI / fac_quot(l + m, l - m) );
  double plm = plgndr_row(l, m, cos(theta) );
  double factor = coef * plm;

  double sin_t = sin(m * phi);
  double cos_t = cos(m * phi);

  result[0] = cos_t * factor;
  result[1] = sin_t * factor;
}

template <unsigned int TDimension>
void SphericalHarmonicPolynomial<TDimension>::EvaluateDTheta(unsigned int from_l,
                                                             unsigned int to_l, double phi, double theta,
                                                             double* sum)
{
  if( m_Coefs.empty() )
    {
    throw SphericalHarmonicPolynomialException(__FILE__, __LINE__, "Coefficients must be specified.");
    }
  if( from_l > to_l )
    {
    throw SphericalHarmonicPolynomialException(__FILE__, __LINE__,
                                               "The starting degree should be smaller or equal to the ending degree.");
    }
  if( to_l > m_Degree )
    {
    throw SphericalHarmonicPolynomialException(__FILE__, __LINE__,
                                               "The evaluated degree mustn't exceed the size of the coefficients.");
    }
  if( m_Coefs.size() < (to_l + 1) * (to_l + 1) )
    {
    throw SphericalHarmonicPolynomialException(__FILE__, __LINE__, "Coefficients size mismatch.");
    }

  int l, m, coord;

  double cotTheta = cos(theta) / sin(theta);
  double cosPhi = cos(phi);
  double sinPhi = sin(phi);

  double realPart, imPart, temp;
  double eval0[2], eval1[2], coef[2];
  for( coord = 0; coord < (int) TDimension; coord++ )
    {
    sum[coord] = 0.;
    }
  for( l = from_l; l <= (int) to_l; l++ ) // for the from_l degree to to_l degree
    {
    for( m = 0; m <= l; m++ )
      {
      this->EvaluateAt(l, m, phi, theta, eval0);
      this->EvaluateAt(l, m + 1, phi, theta, eval1);

      temp = sqrt(l + l * l - m - m * m);
      for( coord = 0; coord < (int) TDimension; coord++ )
        {
        this->GetCoefAt(l, m, coord, coef);
        realPart = m * cotTheta * eval0[0] + temp * (cosPhi * eval1[0] + sinPhi * eval1[1]);
        imPart =   m * cotTheta * eval0[1] + temp * (cosPhi * eval1[1] - sinPhi * eval1[0]);

        sum[coord] += realPart * coef[0] + imPart * coef[1];
        }
      }
    }
}

template <unsigned int TDimension>
void SphericalHarmonicPolynomial<TDimension>::EvaluateDPhi(unsigned int from_l,
                                                           unsigned int to_l, double phi, double theta,
                                                           double* sum)
{
  if( m_Coefs.empty() )
    {
    throw SphericalHarmonicPolynomialException(__FILE__, __LINE__, "Coefficients must be specified.");
    }
  if( from_l > to_l )
    {
    throw SphericalHarmonicPolynomialException(__FILE__, __LINE__,
                                               "The starting degree should be smaller or equal to the ending degree.");
    }
  if( to_l > m_Degree )
    {
    throw SphericalHarmonicPolynomialException(__FILE__, __LINE__,
                                               "The evaluated degree mustn't exceed the size of the coefficients.");
    }
  if( m_Coefs.size() < (to_l + 1) * (to_l + 1) )
    {
    throw SphericalHarmonicPolynomialException(__FILE__, __LINE__, "Coefficients size mismatch.");
    }

  int l, m, coord;

  double eval[2];
  double coefs[2];
  double realPart;
  double imPart;
  for( coord = 0; coord < (int) TDimension; coord++ )
    {
    sum[coord] = 0.;
    }
  for( l = from_l; l <= (int) to_l; l++ ) // for the from_l degree to to_l degree
    {
    for( m = 0; m <= l; m++ )
      {
      this->EvaluateAt(l, m, phi, theta, eval);
      for( coord = 0; coord < (int) TDimension; coord++ )
        {
        this->GetCoefAt(l, m, coord, coefs);
        realPart = -eval[1];
        imPart = eval[0];

        sum[coord] += ( realPart * coefs[0] + imPart * coefs[1] ) * m;
        }
      }
    }
}

template <unsigned int TDimension>
void SphericalHarmonicPolynomial<TDimension>::EvaluateDDThetaTheta(unsigned int from_l,
                                                                   unsigned int to_l, double phi, double theta,
                                                                   double* sum)
{
  if( m_Coefs.empty() )
    {
    throw SphericalHarmonicPolynomialException(__FILE__, __LINE__, "Coefficients must be specified.");
    }
  if( from_l > to_l )
    {
    throw SphericalHarmonicPolynomialException(__FILE__, __LINE__,
                                               "The starting degree should be smaller or equal to the ending degree.");
    }
  if( to_l > m_Degree )
    {
    throw SphericalHarmonicPolynomialException(__FILE__, __LINE__,
                                               "The evaluated degree mustn't exceed the size of the coefficients.");
    }
  if( m_Coefs.size() < (to_l + 1) * (to_l + 1) )
    {
    throw SphericalHarmonicPolynomialException(__FILE__, __LINE__, "Coefficients size mismatch.");
    }

  double cotTheta = 1. / tan(theta);
  double cotTheta2 = cotTheta * cotTheta;
  double cscTheta = 1. / sin(theta);
  double cscTheta2 = cscTheta * cscTheta;
  double cosPhi = cos(phi);
  double sinPhi = sin(phi);
  double cos2Phi = cos(2 * phi);
  double sin2Phi = sin(2 * phi);

  int l, m, coord;
  for( coord = 0; coord < (int) TDimension; coord++ )
    {
    sum[coord] = 0.;
    }

  double eval0[2], eval1[2], eval2[2];
  double coefs[2];
  double realPart, imPart;
  double factor0, factor1, factor2;
  double temp1, temp2;
  double realLastTerm, imLastTerm, swap;
  for( l = (int) from_l; l <= (int) to_l; l++ ) // for the from_l degree to to_l degree
    {
    for( m = 0; m <= l; m++ )
      {
      this->EvaluateAt(l, m,   phi, theta, eval0);
      this->EvaluateAt(l, m + 1, phi, theta, eval1);
      this->EvaluateAt(l, m + 2, phi, theta, eval2);

      realLastTerm = eval2[0] * cos2Phi + eval2[1] * sin2Phi;
      imLastTerm = eval2[1] * cos2Phi - eval2[0] * sin2Phi;

      factor0 = sqrt(l + l * l - m * (m + 1) );
      factor1 = l + l * l - (m + 1) * (m + 2);
      if( factor1 > 0 )
        {
        factor2 = factor0 * sqrt(factor1);
        }
      else
        {
        factor2 = factor0 * sqrt(-factor1);
        swap = realLastTerm;
        realLastTerm = -imLastTerm;
        imLastTerm = swap;
        }
      temp1 = m * (-cscTheta2 + m * cotTheta2);
      temp2 = cotTheta * factor0 * (2 * m + 1);
      for( coord = 0; coord < (int) TDimension; coord++ )
        {
        this->GetCoefAt(l, m, coord, coefs);

        realPart = eval0[0] * temp1
          + (eval1[0] * cosPhi + eval1[1] * sinPhi) * temp2
          + realLastTerm * factor2;

        imPart = eval0[1] * temp1
          + (eval1[1] * cosPhi - eval1[0] * sinPhi) * temp2
          + imLastTerm * factor2;

        sum[coord] += ( realPart * coefs[0] + imPart * coefs[1] );
        }

      }
    }

}

template <unsigned int TDimension>
void SphericalHarmonicPolynomial<TDimension>::EvaluateDDThetaPhi(unsigned int from_l,
                                                                 unsigned int to_l, double phi, double theta,
                                                                 double* sum)
{
  if( m_Coefs.empty() )
    {
    throw SphericalHarmonicPolynomialException(__FILE__, __LINE__, "Coefficients must be specified.");
    }
  if( from_l > to_l )
    {
    throw SphericalHarmonicPolynomialException(__FILE__, __LINE__,
                                               "The starting degree should be smaller or equal to the ending degree.");
    }
  if( to_l > m_Degree )
    {
    throw SphericalHarmonicPolynomialException(__FILE__, __LINE__,
                                               "The evaluated degree mustn't exceed the size of the coefficients.");
    }
  if( m_Coefs.size() < (to_l + 1) * (to_l + 1) )
    {
    throw SphericalHarmonicPolynomialException(__FILE__, __LINE__, "Coefficients size mismatch.");
    }

  int l, m, coord;

  double eval0[2], eval1[2];
  double coefs[2];
  double realPart;
  double imPart;
  double factor;

  double cotTheta = 1. / tan(theta);
  double cosPhi = cos(phi);
  double sinPhi = sin(phi);
  for( coord = 0; coord < (int) TDimension; coord++ )
    {
    sum[coord] = 0.;
    }
  for( l = from_l; l <= (int) to_l; l++ ) // for the from_l degree to to_l degree
    {
    for( m = 0; m <= l; m++ )
      {
      this->EvaluateAt(l, m, phi, theta, eval0);
      this->EvaluateAt(l, m + 1, phi, theta, eval1);
      factor = sqrt(l + l * l - m - m * m);
      for( coord = 0; coord < (int) TDimension; coord++ )
        {
        this->GetCoefAt(l, m, coord, coefs);
        realPart = m * cotTheta * eval0[1] + factor * (cosPhi * eval1[1] - sinPhi * eval1[0]);
        imPart = m * cotTheta * eval0[0] + factor * (cosPhi * eval1[0] + sinPhi * eval1[1]);

        sum[coord] += ( -realPart * coefs[0] + imPart * coefs[1] ) * m;
        }
      }
    }

}

template <unsigned int TDimension>
void SphericalHarmonicPolynomial<TDimension>::EvaluateDDPhiPhi(unsigned int from_l,
                                                               unsigned int to_l, double phi, double theta,
                                                               double* sum)
{
  if( m_Coefs.empty() )
    {
    throw SphericalHarmonicPolynomialException(__FILE__, __LINE__, "Coefficients must be specified.");
    }
  if( from_l > to_l )
    {
    throw SphericalHarmonicPolynomialException(__FILE__, __LINE__,
                                               "The starting degree should be smaller or equal to the ending degree.");
    }
  if( to_l > m_Degree )
    {
    throw SphericalHarmonicPolynomialException(__FILE__, __LINE__,
                                               "The evaluated degree mustn't exceed the size of the coefficients.");
    }
  if( m_Coefs.size() < (to_l + 1) * (to_l + 1) )
    {
    throw SphericalHarmonicPolynomialException(__FILE__, __LINE__, "Coefficients size mismatch.");
    }

  int l, m, coord;

  double eval[2];
  double coefs[2];
  double realPart;
  double imPart;
  for( coord = 0; coord < (int) TDimension; coord++ )
    {
    sum[coord] = 0.;
    }
  for( l = (int) from_l; l <= (int) to_l; l++ ) // for the from_l degree to to_l degree
    {
    for( m = 0; m <= l; m++ )
      {
      this->EvaluateAt(l, m, phi, theta, eval);
      for( coord = 0; coord < (int) TDimension; coord++ )
        {
        this->GetCoefAt(l, m, coord, coefs);
        realPart = eval[0];
        imPart = eval[1];

        sum[coord] -= ( realPart * coefs[0] + imPart * coefs[1] ) * m * m;
        }
      }
    }
}

template <unsigned int TDimension>
void SphericalHarmonicPolynomial<TDimension>::EvaluateDDDThetaThetaTheta
  (unsigned int from_l, unsigned int to_l, double phi, double theta, double* sum)
{
  if( m_Coefs.empty() )
    {
    throw SphericalHarmonicPolynomialException(__FILE__, __LINE__, "Coefficients must be specified.");
    }
  if( from_l > to_l )
    {
    throw SphericalHarmonicPolynomialException(__FILE__, __LINE__,
                                               "The starting degree should be smaller or equal to the ending degree.");
    }
  if( to_l > m_Degree )
    {
    throw SphericalHarmonicPolynomialException(__FILE__, __LINE__,
                                               "The evaluated degree mustn't exceed the size of the coefficients.");
    }
  if( m_Coefs.size() < (to_l + 1) * (to_l + 1) )
    {
    throw SphericalHarmonicPolynomialException(__FILE__, __LINE__, "Coefficients size mismatch.");
    }

  double cotTheta = 1. / tan(theta);
  double cotTheta2 = cotTheta * cotTheta;
  double cscTheta = 1. / sin(theta);
  double cscTheta2 = cscTheta * cscTheta;
  double cosPhi = cos(phi);
  double sinPhi = sin(phi);
  double cos2Phi = cos(2 * phi);
  double sin2Phi = sin(2 * phi);
  double cos3Phi = cos(3 * phi);
  double sin3Phi = sin(3 * phi);

  int l, m, coord;
  for( coord = 0; coord < (int) TDimension; coord++ )
    {
    sum[coord] = 0.;
    }

  double eval0[2], eval1[2], eval2[2], eval3[2];
  double coefs[2];
  double realPart, imPart;
  double factor0, factor1, factor2;
  double temp;
  double term1[2], term2[2], term3[2], term4[2];
  for( l = (int) from_l; l <= (int) to_l; l++ ) // for the from_l degree to to_l degree
    {
    for( m = 0; m <= l; m++ )
      {
      this->EvaluateAt(l, m,   phi, theta, eval0);
      this->EvaluateAt(l, m + 1, phi, theta, eval1);
      this->EvaluateAt(l, m + 2, phi, theta, eval2);
      this->EvaluateAt(l, m + 3, phi, theta, eval3);

      factor0 = sqrt(l + l * l - m * (m + 1) );
      factor1 = l + l * l - (m + 1) * (m + 2);
      factor2 = l + l * l - (m + 2) * (m + 3);
      if( factor2 > 0 )
        {
        factor2 = sqrt(factor2);
        term4[0] = eval3[0] * factor2;
        term4[1] = eval3[1] * factor2;
        }
      else
        {
        factor2 = sqrt(-factor2);
        term4[0] = -eval3[1] * factor2;
        term4[1] = eval3[0] * factor2;
        }

      temp = 3 * cotTheta * (m + 1);
      term3[0] = temp * ( eval2[0] * cosPhi - eval2[1] * sinPhi ) + term4[0];
      term3[1] = temp * ( eval2[1] * cosPhi + eval2[0] * sinPhi ) + term4[1];

      if( factor1 > 0 )
        {
        factor1 = sqrt(factor1);
        term3[0] = term3[0] * factor1;
        term3[1] = term3[1] * factor1;
        }
      else
        {
        factor1 = sqrt(-factor1);
        term3[0] = -term3[1] * factor1;
        term3[1] = term3[0] * factor1;
        }

      temp = cotTheta2 * ( 1 + 3 * m + 3 * m * m ) - cscTheta2 * ( 1 + 3 * m );
      term2[0] = factor0 * temp * ( cos2Phi * eval1[0] - sin2Phi * eval1[1] ) + term3[0];
      term2[1] = factor0 * temp * ( cos2Phi * eval1[1] + sin2Phi * eval1[0] ) + term3[1];

      term1[0] = cos3Phi * term2[0] + sin3Phi * term2[1];
      term1[1] = cos3Phi * term2[1] - sin3Phi * term2[0];

      temp = m * cotTheta * ( m * m * cotTheta2 + (2 - 3 * m) * cscTheta2 );
      for( coord = 0; coord < (int) TDimension; coord++ )
        {
        this->GetCoefAt(l, m, coord, coefs);

        realPart = temp * eval0[0] + term1[0];

        imPart = temp * eval0[1] + term1[1];

        sum[coord] += ( realPart * coefs[0] + imPart * coefs[1] );
        }
      }
    }
}

template <unsigned int TDimension>
void SphericalHarmonicPolynomial<TDimension>::EvaluateDDDThetaThetaPhi
  (unsigned int from_l, unsigned int to_l, double phi, double theta, double* sum)
{
  if( m_Coefs.empty() )
    {
    throw SphericalHarmonicPolynomialException(__FILE__, __LINE__, "Coefficients must be specified.");
    }
  if( from_l > to_l )
    {
    throw SphericalHarmonicPolynomialException(__FILE__, __LINE__,
                                               "The starting degree should be smaller or equal to the ending degree.");
    }
  if( to_l > m_Degree )
    {
    throw SphericalHarmonicPolynomialException(__FILE__, __LINE__,
                                               "The evaluated degree mustn't exceed the size of the coefficients.");
    }
  if( m_Coefs.size() < (to_l + 1) * (to_l + 1) )
    {
    throw SphericalHarmonicPolynomialException(__FILE__, __LINE__, "Coefficients size mismatch.");
    }

  double cotTheta = 1. / tan(theta);
  double cotTheta2 = cotTheta * cotTheta;
  double cscTheta = 1. / sin(theta);
  double cscTheta2 = cscTheta * cscTheta;
  double cosPhi = cos(phi);
  double sinPhi = sin(phi);
  double cos2Phi = cos(2 * phi);
  double sin2Phi = sin(2 * phi);

  int l, m, coord;
  for( coord = 0; coord < (int) TDimension; coord++ )
    {
    sum[coord] = 0.;
    }

  double eval0[2], eval1[2], eval2[2];
  double coefs[2];
  double realPart, imPart;
  double factor0, factor1, factor2;
  double temp;
  double term1[2], term2[2], term3[2];
  for( l = (int) from_l; l <= (int) to_l; l++ ) // for the from_l degree to to_l degree
    {
    for( m = 0; m <= l; m++ )
      {
      this->EvaluateAt(l, m,   phi, theta, eval0);
      this->EvaluateAt(l, m + 1, phi, theta, eval1);
      this->EvaluateAt(l, m + 2, phi, theta, eval2);

      factor0 = sqrt(l + l * l - m * (m + 1) );
      factor1 = l + l * l - (m + 1) * (m + 2);
      if( factor1 > 0 )
        {
        factor2 = sqrt(factor1);
        term3[0] = eval2[0] * factor2;
        term3[1] = eval2[1] * factor2;
        }
      else
        {
        factor2 = sqrt(-factor1);
        term3[0] = -eval2[1] * factor2;
        term3[1] = eval2[0] * factor2;
        }

      temp = cotTheta * (2 * m + 1);
      term2[0] = temp * ( eval1[0] * cosPhi - eval1[1] * sinPhi ) + term3[0];
      term2[1] = temp * ( eval1[1] * cosPhi + eval1[0] * sinPhi ) + term3[1];

      temp = m * (-cscTheta2 + m * cotTheta2);
      term1[0] = temp * ( eval0[0] * cos2Phi - eval0[1] * sin2Phi ) + factor0 * term2[0];
      term1[1] = temp * ( eval0[1] * cos2Phi + eval0[0] * sin2Phi ) + factor0 * term2[1];
      for( coord = 0; coord < (int) TDimension; coord++ )
        {
        this->GetCoefAt(l, m, coord, coefs);

        realPart = term1[1] * cos2Phi - term1[0] * sin2Phi;

        imPart = term1[0] * cos2Phi + term1[1] * sin2Phi;

        sum[coord] += ( -realPart * coefs[0] + imPart * coefs[1] ) * m;
        }
      }
    }
}

template <unsigned int TDimension>
void SphericalHarmonicPolynomial<TDimension>::EvaluateDDDThetaPhiPhi
  (unsigned int from_l, unsigned int to_l, double phi, double theta, double* sum)
{
  if( m_Coefs.empty() )
    {
    throw SphericalHarmonicPolynomialException(__FILE__, __LINE__, "Coefficients must be specified.");
    }
  if( from_l > to_l )
    {
    throw SphericalHarmonicPolynomialException(__FILE__, __LINE__,
                                               "The starting degree should be smaller or equal to the ending degree.");
    }
  if( to_l > m_Degree )
    {
    throw SphericalHarmonicPolynomialException(__FILE__, __LINE__,
                                               "The evaluated degree mustn't exceed the size of the coefficients.");
    }
  if( m_Coefs.size() < (to_l + 1) * (to_l + 1) )
    {
    throw SphericalHarmonicPolynomialException(__FILE__, __LINE__, "Coefficients size mismatch.");
    }

  int l, m, coord;

  double eval0[2], eval1[2];
  double coefs[2];
  double realPart;
  double imPart;
  double factor;

  double cotTheta = 1. / tan(theta);
  double cosPhi = cos(phi);
  double sinPhi = sin(phi);
  for( coord = 0; coord < (int) TDimension; coord++ )
    {
    sum[coord] = 0.;
    }
  for( l = from_l; l <= (int) to_l; l++ ) // for the from_l degree to to_l degree
    {
    for( m = 0; m <= l; m++ )
      {
      this->EvaluateAt(l, m, phi, theta, eval0);
      this->EvaluateAt(l, m + 1, phi, theta, eval1);
      factor = sqrt(l + l * l - m - m * m);
      for( coord = 0; coord < (int) TDimension; coord++ )
        {
        this->GetCoefAt(l, m, coord, coefs);
        realPart = m * cotTheta * eval0[0] + factor * (cosPhi * eval1[0] + sinPhi * eval1[1]);
        imPart = m * cotTheta * eval0[1] + factor * (cosPhi * eval1[1] - sinPhi * eval1[0]);

        sum[coord] += ( realPart * coefs[0] + imPart * coefs[1] ) * m * m * (-1);
        }
      }
    }
}

template <unsigned int TDimension>
void SphericalHarmonicPolynomial<TDimension>::EvaluateDDDPhiPhiPhi
  (unsigned int from_l, unsigned int to_l, double phi, double theta, double* sum)
{
  if( m_Coefs.empty() )
    {
    throw SphericalHarmonicPolynomialException(__FILE__, __LINE__, "Coefficients must be specified.");
    }
  if( from_l > to_l )
    {
    throw SphericalHarmonicPolynomialException(__FILE__, __LINE__,
                                               "The starting degree should be smaller or equal to the ending degree.");
    }
  if( to_l > m_Degree )
    {
    throw SphericalHarmonicPolynomialException(__FILE__, __LINE__,
                                               "The evaluated degree mustn't exceed the size of the coefficients.");
    }
  if( m_Coefs.size() < (to_l + 1) * (to_l + 1) )
    {
    throw SphericalHarmonicPolynomialException(__FILE__, __LINE__, "Coefficients size mismatch.");
    }

  int l, m, coord;

  double eval[2];
  double coefs[2];
  double realPart;
  double imPart;
  for( coord = 0; coord < (int) TDimension; coord++ )
    {
    sum[coord] = 0.;
    }
  for( l = (int) from_l; l <= (int) to_l; l++ ) // for the from_l degree to to_l degree
    {
    for( m = 0; m <= l; m++ )
      {
      this->EvaluateAt(l, m, phi, theta, eval);
      for( coord = 0; coord < (int) TDimension; coord++ )
        {
        this->GetCoefAt(l, m, coord, coefs);
        realPart = eval[1];
        imPart = -eval[0];

        sum[coord] -= ( realPart * coefs[0] + imPart * coefs[1] ) * m * m * m;
        }
      }
    }

}

template <unsigned int TDimension>
void SphericalHarmonicPolynomial<TDimension>::GetCoefAt(unsigned int l, unsigned int m, unsigned int coord,
                                                        double *coef)
{
  if( m == 0 )
    {
    coef[0] = m_Coefs[l * l][coord];
    coef[1] = 0;
    }
  else
    {
    coef[0] = m_Coefs[l * l + 2 * m - 1][coord];
    coef[1] = m_Coefs[l * l + 2 * m][coord];
    }
}

template <unsigned int TDimension>
void SphericalHarmonicPolynomial<TDimension>::GetUnitNormal(unsigned int from_l, unsigned int to_l, double phi,
                                                            double theta,
                                                            double *n)
{
  double dtheta[3], dphi[3];
  double length;

  this->EvaluateDTheta(from_l, to_l, phi, theta, dtheta);
  this->EvaluateDPhi(from_l, to_l, phi, theta, dphi);

  cross(dtheta, dphi, n);

  length = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);

  if( length != 0 )
    {
    length = 1. / length;
    n[0] *= length;
    n[1] *= length;
    n[2] *= length;
    }
}

template <unsigned int TDimension>
void SphericalHarmonicPolynomial<TDimension>::EvaluateFirstFundamentalForm(unsigned int from_l, unsigned int to_l,
                                                                           double phi, double theta,
                                                                           double *result)
{
  double dtheta[3], dphi[3];

  this->EvaluateDTheta(from_l, to_l, phi, theta, dtheta);
  this->EvaluateDPhi(from_l, to_l, phi, theta, dphi);

  result[0] = dot(dtheta, dtheta);
  result[1] = dot(dtheta, dphi);
  result[2] = dot(dphi, dphi);
}

template <unsigned int TDimension>
void SphericalHarmonicPolynomial<TDimension>::EvaluateSecondFundamentalForm(unsigned int from_l, unsigned int to_l,
                                                                            double phi, double theta,
                                                                            double *result)
{
  double dThetaTheta[3], dThetaPhi[3], dPhiPhi[3], normal[3];

  this->EvaluateDDThetaTheta(from_l, to_l, phi, theta, dThetaTheta);
  this->EvaluateDDThetaPhi(from_l, to_l, phi, theta, dThetaPhi);
  this->EvaluateDDPhiPhi(from_l, to_l, phi, theta, dPhiPhi);
  this->GetUnitNormal(from_l, to_l, phi, theta, normal);

  result[0] = dot(dThetaTheta, normal);
  result[1] = dot(dThetaPhi, normal);
  result[2] = dot(dPhiPhi, normal);
}

template <unsigned int TDimension>
void SphericalHarmonicPolynomial<TDimension>::GetPrincipalCurvatures(unsigned int from_l, unsigned int to_l, double phi,
                                                                     double theta,
                                                                     double *kappa)
{
  double first[3], second[3];
  double a, b, c, delta, sqrtDelta;
  double lambda1, lambda2;

  this->EvaluateFirstFundamentalForm(from_l, to_l, phi, theta, first);
  this->EvaluateSecondFundamentalForm(from_l, to_l, phi, theta, second);

  // solve for lambda1, lambda2
  a = first[1] * second[2] - first[2] * second[1];
  b = first[0] * second[2] - first[2] * second[0];
  c = first[0] * second[1] - first[1] * second[0];

  if( a == 0 )
    {
    if( b != 0 )
      {
      lambda1 = lambda2 = -c / b;
      }
    else
      {
      lambda1 = lambda2 = 0;
      }
    }
  else
    {
    delta = b * b - 4 * a * c;
    if( delta < 0 )
      {
      throw SphericalHarmonicPolynomialException(__FILE__, __LINE__, "Complex roots to lambda.");
      }

    sqrtDelta = sqrt(delta);
    lambda1 = (-b - sqrtDelta) / (2 * a);
    lambda2 = (-b + sqrtDelta) / (2 * a);
    }

  // compute the principal curvatures
  double quotient;

  quotient = first[0] + 2 * first[1] * lambda1 + first[2] * lambda1 * lambda1;
  if( quotient != 0 )
    {
    kappa[0] = (second[0] + 2 * second[1] * lambda1 + second[2] * lambda1 * lambda1) / quotient;
    }
  else
    {
    kappa[0] = 0;
    }

  quotient = first[0] + 2 * first[1] * lambda2 + first[2] * lambda2 * lambda2;
  if( quotient != 0 )
    {
    kappa[1] = (second[0] + 2 * second[1] * lambda2 + second[2] * lambda2 * lambda2) / quotient;
    }
  else
    {
    kappa[1] = 0;
    }

  // make sure kappa0 is larger than kappa1
  if( kappa[0] < kappa[1] )
    {
    double temp = kappa[0];
    kappa[0] = kappa[1];
    kappa[1] = temp;
    std::cout << "swap" << std::endl;
    }
}

template <unsigned int TDimension>
void SphericalHarmonicPolynomial<TDimension>::GetPrincipalDirectionsUV(unsigned int from_l, unsigned int to_l,
                                                                       double phi, double theta,
                                                                       double *uv)
{
  double first[3], second[3];

  this->EvaluateFirstFundamentalForm(from_l, to_l, phi, theta, first);
  this->EvaluateSecondFundamentalForm(from_l, to_l, phi, theta, second);

  // compute the matrix A = M_I^-1 * M_II = ( (p q) (r s) )
  double p, q, r, s;
  p = first[2] * second[0] - first[1] * second[1];
  q = first[2] * second[1] - first[1] * second[2];
  r = first[0] * second[1] - first[1] * second[0];
  s = first[0] * second[2] - first[1] * second[1];

  // the characteristic equation of the matrix A is lambda2 + b lambda + c = 0 (a=1)
  double b, c, sqrtdelta, lambda[2];
  b = -p - s;
  c = p * s - q * r;
  sqrtdelta = sqrt( b * b - 4 * c );
  lambda[0] = ( -b + sqrtdelta ) / 2;
  lambda[1] = ( -b - sqrtdelta ) / 2;

  // then, the eigenvectors of A, which are the principal directions in the parameter space, are:
  int i;
  for( i = 0; i < 2; i++ )
    {
    uv[2 * i] = 1;
    if( q != 0 )
      {
      uv[2 * i + 1] = ( -p + lambda[i] ) / q;
      }
    else if( lambda[i] - s != 0 )
      {
      uv[2 * i + 1] = r / (lambda[i] - s);
      }
    else
      {
      uv[2 * i] = uv[2 * i + 1] = 0;
      }
    double length = sqrt( uv[2 * i + 1] * uv[2 * i + 1] + 1 );

    uv[2 * i] /= length;
    uv[2 * i + 1] /= length;
    }

  /*double kappas[2] ;
  this->GetPrincipalCurvatures (from_l, to_l, phi, theta, kappas) ;

  for ( int i = 0 ; i < 2 ; i++ )
  {
  double kappa = kappas[i] ;
  double a, b ;
  double temp, usquare, u, v ;
  a = second[0] - second[2] - kappa * first[0] + kappa * first[2] ;
  b = second[1] - kappa * first[1] ;

  temp = 4 * b * b + a * a ;
  if ( temp != 0 )
    usquare = ( 1 + a / sqrt ( temp ) ) / 2.0 ;
  else
    usquare = 0 ;

  u = sqrt ( usquare ) ;
  v = sqrt ( 1 - usquare ) ;
  uv[2*i] = u ;
  uv[2*i+1] = v ;
  }
  */
}

template <unsigned int TDimension>
void SphericalHarmonicPolynomial<TDimension>::GetPrincipalDirections(unsigned int from_l, unsigned int to_l, double phi,
                                                                     double theta,
                                                                     double *dir)
{
  double uv[4];

  this->GetPrincipalDirectionsUV(from_l, to_l, phi, theta, uv);

  // pull back to the xyz space
  double dtheta[3], dphi[3];
  this->EvaluateDTheta(from_l, to_l, phi, theta, dtheta);
  this->EvaluateDPhi(from_l, to_l, phi, theta, dphi);

  dir[0] = uv[0] * dtheta[0] + uv[1] * dphi[0];
  dir[1] = uv[0] * dtheta[1] + uv[1] * dphi[1];
  dir[2] = uv[0] * dtheta[2] + uv[1] * dphi[2];

  dir[3] = uv[2] * dtheta[0] + uv[3] * dphi[0];
  dir[4] = uv[2] * dtheta[1] + uv[3] * dphi[1];
  dir[5] = uv[2] * dtheta[2] + uv[3] * dphi[2];

  // finally, make them unit vectors
  double l = sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);
  if( l != 0 )
    {
    l = 1. / l;
    dir[0] *= l;
    dir[1] *= l;
    dir[2] *= l;
    }

  l = sqrt(dir[3] * dir[3] + dir[4] * dir[4] + dir[5] * dir[5]);
  if( l != 0 )
    {
    l = 1. / l;
    dir[3] *= l;
    dir[4] *= l;
    dir[5] *= l;
    }
}

template <unsigned int TDimension>
void SphericalHarmonicPolynomial<TDimension>::ComputePrincipalCurve(unsigned int from_l, unsigned int to_l,
                                                                    double startPhi, double startTheta, double *curve,
                                                                    double stepSize,
                                                                    int curveLength)
{
  int    i, curveIndex;
  double x[3];
  double phi, theta;
  double dir[4], prevDir[2];

  // double orientation ;

  phi = startPhi;
  theta = startTheta;

  prevDir[0] = prevDir[1] = 1;
  for( i = curveIndex = 0; i < curveLength; i++, curveIndex += 3 )
    {
    this->Evaluate(from_l, to_l, phi, theta, x);
    this->GetPrincipalDirectionsUV(from_l, to_l, phi, theta, dir);
    if( dir[0] * prevDir[0] + dir[1] * prevDir[1] < 0 )
      {
      dir[0] *= -1;
      dir[1] *= -1;
      }
    phi += dir[1] * stepSize;
    theta += dir[0] * stepSize;
    curve[curveIndex] = x[0];
    curve[curveIndex + 1] = x[1];
    curve[curveIndex + 2] = x[2];
    prevDir[0] = dir[0];
    prevDir[1] = dir[1];
    }
}

template <unsigned int TDimension>
double SphericalHarmonicPolynomial<TDimension>::ComputeC(unsigned int from_l, unsigned int to_l, double phi,
                                                         double theta)
{
  double kappa[2];
  double r;
  double c;

  this->GetPrincipalCurvatures(from_l, to_l, phi, theta, kappa);
  r = sqrt( (kappa[0] * kappa[0] + kappa[1] * kappa[1]) / 2 );
  if( r > 0 )
    {
    c = 2 * log(r) / M_PI;
    }
  else // -"infinity"
    {
    c = -20;
    }

  return c;
}

template <unsigned int TDimension>
double SphericalHarmonicPolynomial<TDimension>::ComputeS(unsigned int from_l, unsigned int to_l, double phi,
                                                         double theta)
{
  double kappa[2];
  double s;

  this->GetPrincipalCurvatures(from_l, to_l, phi, theta, kappa);
  if( kappa[0] != kappa[1] )
    {
    s = 2 * atan( ( kappa[0] + kappa[1] ) / ( kappa[1] - kappa[0] ) ) / M_PI;
    }
  else
    {
    if( kappa[0] > 0 )
      {
      s = 1;
      }
    else
      {
      s = -1;
      }
    }

  return s;
}

template <unsigned int TDimension>
double SphericalHarmonicPolynomial<TDimension>::ComputeGaussianCurvature(unsigned int from_l, unsigned int to_l,
                                                                         double phi,
                                                                         double theta)
{
  double kappa[2];

  this->GetPrincipalCurvatures(from_l, to_l, phi, theta, kappa);

  return kappa[0] * kappa[1];
}

template <unsigned int TDimension>
double SphericalHarmonicPolynomial<TDimension>::ComputeMeanCurvature(unsigned int from_l, unsigned int to_l, double phi,
                                                                     double theta)
{
  double kappa[2];

  this->GetPrincipalCurvatures(from_l, to_l, phi, theta, kappa);
  return ( kappa[0] + kappa[1] ) / 2.0;
}

template <unsigned int TDimension>
void SphericalHarmonicPolynomial<TDimension>::cross( double *u, double *v, double *result )
{
  result[0] = u[1] * v[2] - u[2] * v[1];
  result[1] = u[2] * v[0] - u[0] * v[2];
  result[2] = u[0] * v[1] - u[1] * v[0];
}

template <unsigned int TDimension>
double SphericalHarmonicPolynomial<TDimension>::dot( double *u, double *v )
{
  double result;

  result = u[0] * v[0] + u[1] * v[1] + u[2] * v[2];

  return result;
}

template <unsigned int TDimension>
void SphericalHarmonicPolynomial<TDimension>::ComputePrincipalCurveUV(unsigned int from_l, unsigned int to_l,
                                                                      double startPhi, double startTheta, double *curve,
                                                                      double stepSize,
                                                                      int curveLength)
{
  int    i, curveIndex;
  double phi, theta;
  double dir[4], prevDir[2];

  // double orientation ;

  phi = startPhi;
  theta = startTheta;

  prevDir[0] = prevDir[1] = 1;
  for( i = curveIndex = 0; i < curveLength; i++, curveIndex += 2 )
    {
    this->GetPrincipalDirectionsUV(from_l, to_l, phi, theta, dir);
    if( dir[0] * prevDir[0] + dir[1] * prevDir[1] < 0 )
      {
      dir[0] *= -1;
      dir[1] *= -1;
      }
    curve[curveIndex] = phi;
    curve[curveIndex + 1] = theta;

    phi += dir[1] * stepSize;
    theta += dir[0] * stepSize;

    prevDir[0] = dir[0];
    prevDir[1] = dir[1];
    }
}

template <unsigned int TDimension>
void SphericalHarmonicPolynomial<TDimension>::GetOneRidgePoint(unsigned int from_l, unsigned int to_l,
                                                               unsigned int family, double *uv,
                                                               double *point)
{
  // find a random starting point on the sphere
  double u, v;

  u = rand() / (double) RAND_MAX * M_PI;
  v = rand() / (double) RAND_MAX * M_PI * 2;

  // start a search for a ridge point on the principal curve
  while( true )
    {
    const int LENGTH = 10000;
    // compute a segment of the principal curve starting at this point
    double curve[LENGTH * 2];
    this->ComputePrincipalCurveUV(from_l, to_l, u, v, curve, 0.0005, LENGTH);
    double values[LENGTH];
    double kappa[2];
    for( int i = 0; i < LENGTH; i++ )
      {
      this->GetPrincipalCurvatures(from_l, to_l, curve[i * 2], curve[i * 2 + 1], kappa);
      values[i] = kappa[family];
      }

    // search for a local extremum of curvature on the principal curve
    double prev, current, next, currentMax;
    current = values[0];
    next = values[1];
    currentMax = current;
    bool          max;
    bool          flag = false;
    double        x[3];
    std::ofstream out;
    out.open("test.txt" );
    out << "NUMBER_OF_POINTS=9999" << std::endl << "DIMENSION=3" << std::endl << "TYPE=Curve" << std::endl;
    for( int i = 1; i < LENGTH; i++ )
      {
      prev = current;
      current = next;
      next = values[i + 1];

      this->Evaluate( from_l, to_l, curve[i * 2], curve[i * 2 + 1], x );
      out << x[0] << " " << x[1] << " " << x[2] << std::endl;
      // min = ( current < next ) && ( current < prev ) ;
      max = ( current > next ) && ( current > prev );
      if( max && current > currentMax )
        {
        uv[0] = curve[i * 2];
        uv[1] = curve[i * 2 + 1];
        this->Evaluate( from_l, to_l, uv[0], uv[1], point );
        currentMax = current;
        flag = true;
        // return ;
        }
      }
    if( flag )
      {
      return;
      }
    // this segment does not contain a local extremum - continue with the next segment
    u = curve[LENGTH * 2 - 2];
    v = curve[LENGTH * 2 - 1];
    }
}

template <unsigned int TDimension>
void SphericalHarmonicPolynomial<TDimension>::GradKappa(unsigned int from_l, unsigned int to_l, unsigned int family,
                                                        double *point,
                                                        double *gradKappa)
{
  double first[3], second[3], dTheta[3], dPhi[3], dThetaTheta[3], dThetaPhi[3], dPhiPhi[3], n[3];
  double dThetaThetaTheta[3], dThetaThetaPhi[3], dThetaPhiPhi[3], dPhiPhiPhi[3];
  double a, b, c, delta, sqrtDelta;
  double lambda;

  this->EvaluateFirstFundamentalForm(from_l, to_l, point[0], point[1], first);
  this->EvaluateSecondFundamentalForm(from_l, to_l, point[0], point[1], second);
  this->EvaluateDTheta(from_l, to_l, point[0], point[1], dTheta);
  this->EvaluateDPhi(from_l, to_l, point[0], point[1], dPhi);
  this->EvaluateDDThetaTheta(from_l, to_l, point[0], point[1], dThetaTheta);
  this->EvaluateDDThetaPhi(from_l, to_l, point[0], point[1], dThetaPhi);
  this->EvaluateDDPhiPhi(from_l, to_l, point[0], point[1], dPhiPhi);
  this->EvaluateDDDThetaThetaTheta(from_l, to_l, point[0], point[1], dThetaThetaTheta);
  this->EvaluateDDDThetaThetaPhi(from_l, to_l, point[0], point[1], dThetaThetaPhi  );
  this->EvaluateDDDThetaPhiPhi(from_l, to_l, point[0], point[1], dThetaPhiPhi    );
  this->EvaluateDDDPhiPhiPhi(from_l, to_l, point[0], point[1], dPhiPhiPhi      );
  this->GetUnitNormal(from_l, to_l, point[0], point[1], n);

  // solve for lambda1, lambda2
  a = first[1] * second[2] - first[2] * second[1];
  b = first[0] * second[2] - first[2] * second[0];
  c = first[0] * second[1] - first[1] * second[0];

  if( a == 0 )
    {
    if( b != 0 )
      {
      lambda = -c / b;
      }
    else
      {
      lambda = 0;
      }
    }
  else
    {
    delta = b * b - 4 * a * c;
    if( delta < 0 )
      {
      throw SphericalHarmonicPolynomialException(__FILE__, __LINE__, "Complex roots to lambda.");
      }
    if( delta < 0.0001 )
      {
      std::cout << "Umbilic!" << std::endl;
      gradKappa[0] = this->lastRidgeDir[0];
      gradKappa[1] = this->lastRidgeDir[1];

      return;
      }
    sqrtDelta = sqrt(delta);
    if( family == 0 )
      {
      lambda = (-b - sqrtDelta) / (2 * a);
      }
    else
      {
      lambda = (-b + sqrtDelta) / (2 * a);
      }
    }

  if( b == lambda * a * 2 )
    {
    gradKappa[0] =  this->lastRidgeDir[0];
    gradKappa[1] =  this->lastRidgeDir[1];
    return;
    }

  double EDTheta, FDTheta, GDTheta;
  double EDPhi, FDPhi, GDPhi;
  double eDTheta, fDTheta, gDTheta;
  double eDPhi, fDPhi, gDPhi;
  double lambdaDTheta, lambdaDPhi;
  double nDTheta[3], nDPhi[3];
  double aDTheta, bDTheta, cDTheta;
  double aDPhi, bDPhi, cDPhi;

  EDTheta = 2 * dot(dThetaTheta, dTheta);
  EDPhi = 2 * dot(dThetaPhi, dTheta);
  FDTheta = dot(dTheta, dThetaPhi) + dot(dPhi, dThetaTheta);
  FDPhi = dot(dTheta, dPhiPhi) + dot(dThetaPhi, dPhi);
  GDTheta = 2 * dot(dPhi, dThetaPhi);
  GDPhi = 2 * dot(dPhi, dPhiPhi);

  double temp1[3], temp2[3];
  cross( dTheta, dThetaPhi, temp1 );
  cross( dThetaTheta, dPhi, temp2 );
  nDTheta[0] = temp1[0] + temp2[0];
  nDTheta[1] = temp1[1] + temp2[1];
  nDTheta[2] = temp1[2] + temp2[2];

  cross( dTheta, dPhiPhi, temp1 );
  cross( dThetaPhi, dPhi, temp2 );
  nDPhi[0] = temp1[0] + temp2[0];
  nDPhi[1] = temp1[1] + temp2[1];
  nDPhi[2] = temp1[2] + temp2[2];

  eDTheta = dot( dThetaTheta, nDTheta ) + dot( dThetaThetaTheta, n );
  eDPhi = dot( dThetaTheta, nDPhi ) + dot( dThetaThetaPhi, n );
  fDTheta = dot( dThetaPhi, nDTheta ) + dot( dThetaThetaPhi, n );
  fDPhi = dot( dThetaPhi, nDPhi ) + dot( dThetaPhiPhi, n );
  gDTheta = dot( dPhiPhi, nDTheta ) + dot( dThetaPhiPhi, n );
  gDPhi = dot( dPhiPhi, nDPhi ) + dot( dPhiPhiPhi, n );

  aDTheta = second[2] * FDTheta + first[1] * gDTheta - second[1] * GDTheta - first[2] * fDTheta;
  aDPhi   = second[2] * FDPhi   + first[1] * gDPhi   - second[1] * GDPhi   - first[2] * fDPhi;
  bDTheta = second[2] * EDTheta + first[0] * gDTheta - second[0] * GDTheta - first[2] * eDTheta;
  bDPhi   = second[2] * EDPhi   + first[0] * gDPhi   - second[0] * GDPhi   - first[2] * eDPhi;
  cDTheta = second[1] * EDTheta + first[0] * fDTheta - second[0] * FDTheta - first[1] * eDTheta;
  cDPhi   = second[1] * EDPhi   + first[0] * fDPhi   - second[0] * FDPhi   - first[1] * eDPhi;

  lambdaDTheta = ( cDTheta + lambda * bDTheta + lambda * lambda * aDTheta ) / ( -2 * a * lambda + b );
  lambdaDPhi = ( cDPhi + lambda * bDPhi + lambda * lambda * aDPhi ) / ( -2 * a * lambda + b );

  double top = second[0] + 2 * second[1] * lambda + second[2] * lambda * lambda;
  double bottom = first[0] + 2 * first[1] * lambda + first[2] * lambda * lambda;

  double TopDTheta, BottomDTheta;
  double TopDPhi,   BottomDPhi;
  TopDTheta = eDTheta + 2 * lambda * fDTheta + 2 * second[1] * lambdaDTheta + lambda * lambda * gDTheta + 2
    * second[2] * lambda * lambdaDTheta;
  TopDPhi = eDPhi + 2 * lambda * fDPhi + 2 * second[1] * lambdaDPhi + lambda * lambda * gDPhi + 2 * second[2]
    * lambda * lambdaDPhi;
  BottomDTheta = EDTheta + 2 * lambda * FDTheta + 2 * first[1] * lambdaDTheta + lambda * lambda * GDTheta + 2
    * first[2] * lambda * lambdaDTheta;
  BottomDPhi = EDPhi + 2 * lambda * FDPhi + 2 * first[1] * lambdaDPhi + lambda * lambda * GDPhi + 2 * first[2]
    * lambda * lambdaDPhi;

  double KappaDTheta, KappaDPhi;
  KappaDTheta = ( TopDTheta * bottom - top * BottomDTheta );
  KappaDPhi = ( TopDPhi * bottom - top * BottomDPhi );

  double temp = sqrt( KappaDTheta * KappaDTheta + KappaDPhi * KappaDPhi );
  if( temp == 0 )
    {
    temp = 1;
    }
  this->lastRidgeDir[0] = gradKappa[0] = KappaDTheta / temp;
  this->lastRidgeDir[1] = gradKappa[1] = KappaDPhi / temp;
}

template <unsigned int TDimension>
void SphericalHarmonicPolynomial<TDimension>::GetGradKappa(unsigned int from_l, unsigned int to_l, unsigned int family,
                                                           double phi, double theta,
                                                           double *dir)
{
  double uv[2], point[2];

  point[0] = phi;
  point[1] = theta;
  this->GradKappa(from_l, to_l, family, point, uv);

  // pull back to the xyz space
  double dtheta[3], dphi[3];
  this->EvaluateDTheta(from_l, to_l, phi, theta, dtheta);
  this->EvaluateDPhi(from_l, to_l, phi, theta, dphi);

  dir[0] = uv[1] * dtheta[0] + uv[0] * dphi[0];
  dir[1] = uv[1] * dtheta[1] + uv[0] * dphi[1];
  dir[2] = uv[1] * dtheta[2] + uv[0] * dphi[2];

  // finally, make them unit vectors
  double l = sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);
  if( l != 0 )
    {
    l = 1. / l;
    dir[0] *= l;
    dir[1] *= l;
    dir[2] *= l;
    }

}

template <unsigned int TDimension>
void SphericalHarmonicPolynomial<TDimension>::FollowRidge(unsigned int from_l, unsigned int to_l, unsigned int family,
                                                          double *startPoint, double *uv, double *point,
                                                          double stepSize)
{
  double _GradKappa[2];

  this->GradKappa(from_l, to_l, 1 - family, startPoint, _GradKappa);
  uv[0] = startPoint[0] + stepSize * _GradKappa[0];
  uv[1] = startPoint[1] + stepSize * _GradKappa[1];

  this->Evaluate( from_l, to_l, uv[0], uv[1], point );
}

template <unsigned int TDimension>
void SphericalHarmonicPolynomial<TDimension>::GetCaustic(unsigned int from_l, unsigned int to_l, int family, double phi,
                                                         double theta,
                                                         double *point)
{
  double kappa[2], normal[3], x[3];

  this->GetPrincipalCurvatures(from_l, to_l, phi, theta, kappa);
  this->GetUnitNormal(from_l, to_l, phi, theta, normal);
  this->Evaluate( from_l, to_l, phi, theta, x );

  double dist;

  if( kappa[1 - family] != 0 )
    {
    dist = 1.0 / kappa[1 - family];
    }
  else
    {
    dist = 0;
    }
  for( int i = 0; i < 3; i++ )
    {
    point[i] = x[i] + normal[i] * dist;
    }

}

} // end namespace neurolib

#endif
