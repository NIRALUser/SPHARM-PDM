#ifndef _namicSphericalHarmonicPolynomial_h
#define _namicSphericalHarmonicPolynomial_h

/** \class SphericalHarmonicPolynomial
 *
 *  \brief This class provides 2D(phi, theta) Spherical Harmonics basis functions.
 *
 *  \author Christine Xu
 */

#include <vector>

#include <itkMacro.h>
#include <itkPoint.h>

namespace neurolib
{

class SphericalHarmonicPolynomialException : public itk::ExceptionObject
{
public:

  /** Run-time information. */
  itkTypeMacro( SphericalHarmonicPolynomialException, ExceptionObject );

  /** Constructor. */
  SphericalHarmonicPolynomialException(const char *file, unsigned int line,
                                       const char* message = "Error in evaluate Spherical Harmonic polynomials") :
    ExceptionObject(file, line)
  {
    SetDescription(message);
  }

  /** Constructor. */
  SphericalHarmonicPolynomialException(const std::string & file, unsigned int line,
                                       const char* message = "Error in evaluate Spherical Harmonic polynomials") :
    ExceptionObject(file, line)
  {
    SetDescription(message);
  }

};

template <unsigned int TDimension = 3>
class SphericalHarmonicPolynomial
{

public:
  typedef SphericalHarmonicPolynomial Self;

  typedef double                             ScalarType;
  typedef itk::Point<ScalarType, TDimension> CoefType;
  typedef std::vector<CoefType>              CoefListType;

  /*typedef SphericalHarmonicSpatialObject::ScalarType ScalarType;
  typedef SphericalHarmonicSpatialObject::CoefType CoefType;
  typedef SphericalHarmonicSpatialObject::CoefListType CoefListType;*/

  /*unsigned int GetDimension(void) const
  {return m_Dimension;}
  void SetDimension(unsigned int d)
  {m_Dimension = d;}*/

  unsigned int GetDegree(void) const
  {
    return m_Degree;
  }

  void SetDegree(unsigned int d)
  {
    m_Degree = d;
  }

  void SetCoefs(CoefListType& coeflist);

  void GetCoefAt(unsigned int l, unsigned int m, unsigned int coord, double *coef);

  void Evaluate(unsigned int from_l, unsigned int to_l, double phi, double theta, double* sum);

  void EvaluateAt(unsigned int l, unsigned int m, double phi, double theta, double* result);

  // First order derivatives
  void EvaluateDTheta(unsigned int from_l, unsigned int to_l, double phi, double theta, double* sum);

  void EvaluateDPhi(unsigned int from_l, unsigned int to_l, double phi, double theta, double* sum);

  // Second order derivatives
  void EvaluateDDThetaTheta(unsigned int from_l, unsigned int to_l, double phi, double theta, double* sum);

  void EvaluateDDThetaPhi(unsigned int from_l, unsigned int to_l, double phi, double theta, double* sum);

  void EvaluateDDPhiPhi(unsigned int from_l, unsigned int to_l, double phi, double theta, double* sum);

  // Third order derivatives
  void EvaluateDDDThetaThetaTheta(unsigned int from_l, unsigned int to_l, double phi, double theta, double* sum);

  void EvaluateDDDThetaThetaPhi(unsigned int from_l, unsigned int to_l, double phi, double theta, double* sum);

  void EvaluateDDDThetaPhiPhi(unsigned int from_l, unsigned int to_l, double phi, double theta, double* sum);

  void EvaluateDDDPhiPhiPhi(unsigned int from_l, unsigned int to_l, double phi, double theta, double* sum);

  // Unit normal vector to the surface
  void GetUnitNormal(unsigned int from_l, unsigned int to_l, double phi, double theta, double *n);

  // First and second fundamental forms
  // CHECK
  void EvaluateFirstFundamentalForm(unsigned int from_l, unsigned int to_l, double phi, double theta, double *result);

  void EvaluateSecondFundamentalForm(unsigned int from_l, unsigned int to_l, double phi, double theta, double *result);

  // CHECK
  void GetPrincipalCurvatures(unsigned int from_l, unsigned int to_l, double phi, double theta, double *kappa);

  void GetPrincipalDirections(unsigned int from_l, unsigned int to_l, double phi, double theta, double *dir);

  void GetPrincipalDirectionsUV(unsigned int from_l, unsigned int to_l, double phi, double theta, double *uv);

  void ComputePrincipalCurve(unsigned int from_l, unsigned int to_l, double startPhi, double startTheta, double *curve,
                             double stepSize,
                             int curveLength);

  // CHECK
  double ComputeC(unsigned int from_l, unsigned int to_l, double phi, double theta);

  double ComputeS(unsigned int from_l, unsigned int to_l, double phi, double theta);

  double ComputeGaussianCurvature(unsigned int from_l, unsigned int to_l, double phi, double theta);

  double ComputeMeanCurvature(unsigned int from_l, unsigned int to_l, double phi, double theta);

  void GetGradKappa(unsigned int from_l, unsigned int to_l, unsigned int family, double phi, double theta,
                    double *gradKappa);

  void GetOneRidgePoint(unsigned int from_l, unsigned int to_l, unsigned int family, double *uv, double *point);

  void FollowRidge(unsigned int from_l, unsigned int to_l, unsigned int family, double *startPoint, double *uv,
                   double *point,
                   double stepSize);

  void GetCaustic(unsigned int from_l, unsigned int to_l, int family, double phi, double theta, double *point);

  SphericalHarmonicPolynomial();

  ~SphericalHarmonicPolynomial();
protected:

  double plgndr_row(int l, int m, double x);

  double fac_quot(int a, int b);

  void ComputePrincipalCurveUV(unsigned int from_l, unsigned int to_l, double startPhi, double startTheta,
                               double *curve, double stepSize,
                               int curveLength);

  void GradKappa(unsigned int from_l, unsigned int to_l, unsigned int family, double *point, double *gradKappa);

  double lastRidgeDir[2];
private:

  // unsigned int m_Dimension;

  unsigned int m_Degree; // m_Degree >= 0

  CoefListType m_Coefs;

  void cross( double *u, double *v, double *result );

  double dot( double *u, double *v );

};

} // end namespace neurolib

#ifndef ITK_MANUAL_INSTANTIATION
#include "SphericalHarmonicPolynomial.txx"
#endif

#endif // _namicSphericalHarmonicPolynomial_h
