//# Copyright (C) 1996-2010
//# Associated Universities, Inc. Washington DC, USA.
//#
//# This library is free software; you can redistribute it and/or modify it
//# under the terms of the GNU Library General Public License as published by
//# the Free Software Foundation; either version 2 of the License, or (at your
//# option) any later version.
//#
//# This library is distributed in the hope that it will be useful, but WITHOUT
//# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
//# License for more details.
//#
//# You should have received a copy of the GNU Library General Public License
//# along with this library; if not, write to the Free Software Foundation,
//# Inc., 675 Massachusetts Ave, Cambridge, MA 02139, USA.
//#
//# Correspondence concerning AIPS++ should be addressed as follows:
//#        Internet email: casa-feedback@nrao.edu.
//#        Postal address: AIPS++ Project Office
//#                        National Radio Astronomy Observatory
//#                        520 Edgemont Road
//#                        Charlottesville, VA 22903-2475 USA
//#
//# $Id:  $AspClean2026.cc
// Same include list as in MatrixCleaner.cc
#include <casacore/casa/Arrays/Matrix.h>
#include <casacore/casa/Arrays/Cube.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/MatrixMath.h>
#include <casacore/casa/BasicMath/Math.h>
#include <casacore/casa/BasicSL/Complex.h>
#include <casacore/casa/Logging/LogIO.h>
#include <casacore/casa/OS/File.h>
#include <casacore/casa/Containers/Record.h>

#include <casacore/lattices/LRegions/LCBox.h>
#include <casacore/casa/Arrays/Slicer.h>
#include <casacore/scimath/Mathematics/FFTServer.h>
#include <casacore/casa/OS/HostInfo.h>
#include <casacore/casa/Arrays/ArrayError.h>
#include <casacore/casa/Arrays/ArrayIter.h>
#include <casacore/casa/Arrays/VectorIter.h>

#include <casacore/casa/Utilities/GenSort.h>
#include <casacore/casa/BasicSL/String.h>
#include <casacore/casa/Utilities/Assert.h>
#include <casacore/casa/Utilities/Fallible.h>

#include <casacore/casa/BasicSL/Constants.h>
#include <casacore/casa/Logging/LogSink.h>
#include <casacore/casa/Logging/LogMessage.h>

#include <synthesis/MeasurementEquations/MatrixCleaner.h>
#include <synthesis/TransformMachines/StokesImageUtil.h>
#include <synthesis/TransformMachines2/Utils.h>
#include <casacore/coordinates/Coordinates/TabularCoordinate.h>

#ifdef _OPENMP
#include <omp.h>
#endif

// Additional include files
#include <algorithm>
#include <cmath>
#include <limits>
#include <functional>
#include <vector>
#include <iostream>
#include <iomanip>

#include <synthesis/MeasurementEquations/AspClean2026.h>

// for alglib
#include <synthesis/MeasurementEquations/objfunc_alglib.h>

using namespace alglib;

using namespace casacore;
using namespace std;

namespace casa {

namespace {

// ALGLIB VJ interface: residual vector and analytic Jacobian.
// The lambdas are constructed inside optimizeComponent(), so no header changes
// or additional externally visible class members are required.
struct Asp2026LMContext {
  std::function<void(const real_1d_array&, real_1d_array&,
                     real_2d_array*)> evaluate;
};

void asp2026Residuals(const real_1d_array& x,
                      real_1d_array& fi, void* ptr)
{
  static_cast<Asp2026LMContext*>(ptr)->evaluate(x, fi, nullptr);
}

void asp2026ResidualsJacobian(const real_1d_array& x,
                              real_1d_array& fi,
                              real_2d_array& jac, void* ptr)
{
  static_cast<Asp2026LMContext*>(ptr)->evaluate(x, fi, &jac);
}

struct AspClean2026LBFGSContext { AspClean2026* cleaner; };

constexpr Float fracMin = 0.65f;
constexpr Float fracMax = 0.99f;
constexpr Float minEpsFrac = 1.0e-6f;
constexpr double maxTransform = 20.0;

} // namespace

AspClean2026::AspClean2026()
    : //itsResidualComponentSupport(2, 1),
      itsFrac(0.65f),
      itsEpsFrac(1.0f),
      itsPeakValue(0.0f),
      itsPeakSign(1.0f)
{}

AspClean2026::~AspClean2026() = default;

void AspClean2026::setInitScaleXfrs(const Float /*width*/)
{
  // Keep one initial scale so the AspMatrixCleaner base code
  // does not break. AspClean2026 does not use this scale to define its model.
  destroyAspScales();
  itsInitScaleSizes.resize(1, false);
  itsInitScaleSizes[0] = 0.0f;
  itsNInitScales = 1;
  itsNscales = 1;

  // do not initialize the following for asp2026 b/c
  // it causes seg fault as well as increasing memory usage
  /*itsInitScales.resize(1, false);
  itsInitScaleXfrs.resize(1, false);
  itsInitScales[0] = Matrix<Float>(psfShape_p);
  makeInitScaleImage(itsInitScales[0], 0.0f);*/
  fft = FFTServer<Float, Complex>(psfShape_p);
  //fft.fft0(itsInitScaleXfrs[0], itsInitScales[0]);
}

Bool AspClean2026::setInitScaleMasks(
    const Matrix<Float>& mask,
    const Float& maskThreshold)
{
  destroyMasks();
  destroyInitMasks();

  itsMask = new Matrix<Float>(mask.shape());
  itsMask->assign(mask);
  itsMaskThreshold = maskThreshold;
  noClean_p = max(*itsMask) < itsMaskThreshold;

  if (noClean_p)
    return false;

  // For AspClean2026 this is the CLEAN window Omega, not a scale mask.
  itsInitScaleMasks.resize(1);
  itsInitScaleMasks[0] = mask;

  blcDirty = IPosition(mask.shape().nelements(), 0);
  trcDirty = mask.shape() - 1;

  return true;
}

std::vector<Float> AspClean2026::getActiveSetAspen(const float /*peakres*/)
{
  // NO-OP
  // The AspClean2026 component is selected and optimized in optimizeComponent().
  itsGoodAspActiveSet.clear();
  itsGoodAspAmplitude.clear();
  itsGoodAspCenter.clear();
  return {};
}

Bool AspClean2026::insideCleanWindow(Int i, Int j) const
{
    if (itsMask.null())
        return true;

    if (itsMaskThreshold < 0.0f)
        return (*itsMask)(i, j) > 0.0f;

    return (*itsMask)(i, j) >= itsMaskThreshold;
}


void AspClean2026::extractResidualPeak()
{
  AlwaysAssert(!itsDirty.null(), AipsError);
  AlwaysAssert(itsInitScaleMasks.nelements() == 1, AipsError);

  const IPosition shape = itsDirty->shape();
  Float peak = 0.0f;
  IPosition peakPosition(shape.nelements(), 0);

  for (Int j = blcDirty(1); j <= trcDirty(1); ++j) 
  {
    for (Int i = blcDirty(0); i <= trcDirty(0); ++i) 
    {
      if (!insideCleanWindow(i, j))
        continue;

      const Float value = (*itsDirty)(i, j);
      //if (std::abs(value) > std::abs(peak)) // mode=absolute
      if (value > peak) // mode=positive
      {
        peak = value;
        peakPosition = IPosition(2, i, j);
      }
    }
  }

  itsPositionOptimum = peakPosition;
  itsPeakValue = peak;
  itsPeakResidual = std::abs(peak);
  itsPeakSign = (peak >= 0.0f) ? 1.0f : -1.0f;

  // These values keep the legacy stopping/bookkeeping code well-defined,
  // but they are not used to scale the 2026 component update.
  itsOptimumScale = 0;
  itsOptimumScaleSize = -1.0f;
}

Matrix<Float> AspClean2026::makeTanhComponent(
    const Float frac,
    const Float epsFrac) const
{
  const IPosition shape = itsDirty->shape();
  Matrix<Float> component(shape);
  component = 0.0f;

  if (itsPeakValue == 0.0f)
    return component;

  const Float peakAmplitude = std::abs(itsPeakValue);
  const Float threshold = frac * peakAmplitude;
  const Float epsilon =
      epsFrac * peakAmplitude + 1.0e-12f;

  // Eq. (20) and (27): the component is constructed from residual pixels
  // inside the clean window using the signed peak and the adaptive threshold.  
  for (Int j = blcDirty(1); j <= trcDirty(1); ++j) 
  {
    for (Int i = blcDirty(0); i <= trcDirty(0); ++i) 
    {
      if (!insideCleanWindow(i, j))
        continue;

      const Float residual = (*itsDirty)(i, j);
      const Float signedResidual = itsPeakSign * residual;

      if (signedResidual >= threshold)
      {
        component(i, j) = residual *
            std::tanh((signedResidual - threshold) / epsilon);
      }
    }
  }

  return component;
}

/*Matrix<Float> AspClean2026::convolveComponent(
    const Matrix<Float>& component)
{
  Matrix<Complex> componentXfr;
  fft.fft0(componentXfr, component);

  Matrix<Complex> work((*itsXfr) * componentXfr);
  Matrix<Float> convolved(component.shape());

  fft.fft0(convolved, work, false);
  fft.flip(convolved, false, false);

  // genie do I need shift correction like the Asp base class?
  /*const Int nx = convolved.shape()(0);
  const Int ny = convolved.shape()(1);
  IPosition nullnull(2, 0);
  IPosition sup(2, nx, ny);
  Matrix<Float> shifted(convolved.shape());

  if (itsdimensionsareeven) {
      Matrix<Float> sub = convolved(nullnull + 1, sup - 1);
      sub.assign_conforming(shifted(nullnull, sup - 2));
  } else {
      Matrix<Float> sub = convolved(nullnull + 2, sup - 1);
      sub.assign_conforming(shifted(nullnull, sup - 3));
  }* /

  return convolved;
}*/


Matrix<Float> AspClean2026::convolveComponent(
    const Matrix<Float>& component)
{
  // Forward FFT of the component.
  Matrix<Complex> componentXfr;
  fft.fft0(componentXfr, component);

  // Multiply by the PSF Fourier transform.
  Matrix<Complex> work((*itsXfr) * componentXfr);

  // Inverse FFT.
  Matrix<Float> convolved(component.shape());
  fft.fft0(convolved, work, false);
  fft.flip(convolved, false, false);

  // ---------------------------------------------------------
  // Pixel-origin correction used by AspMatrixCleaner.cc
  // (Hendrik's pixel-shifting fix).
  //
  // Even-sized image: shift by 1 pixel.
  // Odd-sized image:  shift by 2 pixels.
  //
  // Use the actual matrix dimensions to determine parity.
  // This avoids relying on itsdimensionsareeven being set
  // in the legacy AspMatrixCleaner::optimizeComponent().
  // ---------------------------------------------------------

  const IPosition support = convolved.shape();
  const IPosition nullnull(2, 0);

  const Int nx = support[0];
  const Int ny = support[1];

  const Bool evenDimensions =
      (nx % 2 == 0) && (ny % 2 == 0);

  const Int offset = evenDimensions ? 1 : 2;

  // Keep a copy of the original convolution.
  Matrix<Float> shifted(support);
  shifted.assign_conforming(convolved);

  // Initialize the output so that uncovered boundary
  // pixels do not retain the unshifted convolution.
  convolved = 0.0f;

  if (nx > offset && ny > offset)
  {
    const IPosition shiftOffset(2, offset, offset);

    Matrix<Float> destination =
        convolved(nullnull + shiftOffset, support - 1);

    Matrix<Float> source =
        shifted(nullnull, support - shiftOffset - 1);

    destination.assign_conforming(source);
  }

  return convolved;
}



double AspClean2026::objectiveAndGradient(
    const Float frac,
    const Float epsFrac,
    double& fracGradient,
    double& epsFracGradient) 
{
  const Matrix<Float> localComponent =
      makeTanhComponent(frac, epsFrac);
  const Matrix<Float> localConvolved =
      convolveComponent(localComponent);

  /*if (component != nullptr)
    *component = localComponent;
  if (convolved != nullptr)
    *convolved = localConvolved;

  Float chi2 = 0.0f;
  const IPosition shape = itsDirty->shape();*/
  double chi2 = 0.0;
  fracGradient = 0.0;
  epsFracGradient = 0.0;

  // Eq. 28: evaluate the objective only over the clean window.
  for (Int j = blcDirty(1); j <= trcDirty(1); ++j) 
  {
    for (Int i = blcDirty(0); i <= trcDirty(0); ++i) 
    {
      if (!insideCleanWindow(i, j))
        continue;

      //const Float difference =
      //    (*itsDirty)(i, j) - localConvolved(i, j);
      const double difference =
          double((*itsDirty)(i, j)) - double(localConvolved(i, j));

      chi2 += difference * difference;
    }
  }

  // Differentiate the smooth part of the thresholded component. Pixels
  // exactly on the moving threshold are non-differentiable and have measure
  // zero, so use the derivative of the active branch.
  Matrix<Float> componentDerivative(itsDirty->shape(), 0.0f);
  const Float peakAmplitude = std::abs(itsPeakValue);
  const Float threshold = frac * peakAmplitude;
  const Float epsilon = epsFrac * peakAmplitude + 1.0e-12f;

  for (Int j = blcDirty(1); j <= trcDirty(1); ++j)
  {
    for (Int i = blcDirty(0); i <= trcDirty(0); ++i)
    {
      if (!insideCleanWindow(i, j))
        continue;

      const Float residual = (*itsDirty)(i, j);
      const Float signedResidual = itsPeakSign * residual;
      if (signedResidual >= threshold)
      {
        const Float z = (signedResidual - threshold) / epsilon;
        const Float tanhZ = std::tanh(z);
        const Float sechSquared = 1.0f - tanhZ * tanhZ;
        componentDerivative(i, j) =
            -residual * sechSquared / epsFrac;
      }
    }
  }

  Matrix<Float> convolvedDerivative =
      convolveComponent(componentDerivative);
  for (Int j = blcDirty(1); j <= trcDirty(1); ++j)
    for (Int i = blcDirty(0); i <= trcDirty(0); ++i)
      if (insideCleanWindow(i, j))
        fracGradient += -2.0 *
            (double((*itsDirty)(i, j)) - double(localConvolved(i, j))) *
            double(convolvedDerivative(i, j));

  componentDerivative = 0.0f;
  for (Int j = blcDirty(1); j <= trcDirty(1); ++j)
  {
    for (Int i = blcDirty(0); i <= trcDirty(0); ++i)
    {
      if (!insideCleanWindow(i, j))
        continue;

      const Float residual = (*itsDirty)(i, j);
      const Float signedResidual = itsPeakSign * residual;
      if (signedResidual >= threshold)
      {
        const Float z = (signedResidual - threshold) / epsilon;
        const Float tanhZ = std::tanh(z);
        const Float sechSquared = 1.0f - tanhZ * tanhZ;
        componentDerivative(i, j) =
            -residual * sechSquared * z / epsFrac;
      }
    }
  }

  convolvedDerivative = convolveComponent(componentDerivative);
  for (Int j = blcDirty(1); j <= trcDirty(1); ++j)
    for (Int i = blcDirty(0); i <= trcDirty(0); ++i)
      if (insideCleanWindow(i, j))
        epsFracGradient += -2.0 *
            (double((*itsDirty)(i, j)) - double(localConvolved(i, j))) *
            double(convolvedDerivative(i, j));


  return chi2;
}

Float AspClean2026::transformedFrac(const double x)
{
  const double y = std::clamp(x, -maxTransform, maxTransform);
  const double logistic = 1.0 / (1.0 + std::exp(-y));
  return Float(fracMin + (fracMax - fracMin) * logistic);
}

Float AspClean2026::transformedEpsFrac(const double x)
{
  const double y = std::clamp(x, -maxTransform, maxTransform);
  return Float(std::exp(y));
}

double AspClean2026::inverseFrac(const Float frac)
{
  const double normalized = std::clamp(
      (double(frac) - fracMin) / (fracMax - fracMin),
      1.0e-8,
      1.0 - 1.0e-8);
  return std::log(normalized / (1.0 - normalized));
}

double AspClean2026::inverseEpsFrac(const Float epsFrac)
{
  return std::log(std::max(double(epsFrac), double(minEpsFrac)));
}

void AspClean2026::lbfgsObjective(
    const alglib::real_1d_array& x,
    double& func,
    alglib::real_1d_array& gradient,
    void* ptr)
{
  auto* context = static_cast<AspClean2026LBFGSContext*>(ptr);
  AspClean2026* cleaner = context->cleaner;

  const Float frac = transformedFrac(x[0]);
  const Float epsFrac = transformedEpsFrac(x[1]);

  //func = double(cleaner->objective(frac, epsFrac));
  double fracGradient = 0.0;
  double epsFracGradient = 0.0;
  func = cleaner->objectiveAndGradient(
      frac, epsFrac, fracGradient, epsFracGradient);

  const double fracRange = double(fracMax - fracMin);
  gradient[0] = std::abs(x[0]) >= maxTransform ? 0.0 : fracGradient *
      (double(frac - fracMin) * double(fracMax - frac) / fracRange);
  gradient[1] = std::abs(x[1]) >= maxTransform
      ? 0.0 : epsFracGradient * double(epsFrac);
}

/*void AspClean2026::optimizeComponent()
{
  extractResidualPeak();

  itsResidualComponent.resize(itsDirty->shape());
  itsResidualComponent = 0.0f;
  itsConvolvedComponent.resize(itsDirty->shape());
  itsConvolvedComponent = 0.0f;

  if (itsPeakValue == 0.0f)
  {
    itsStrengthOptimum = 0.0f;
    return;
  }

  // Author: lmfit(method="least_squares") -> scipy.optimize.least_squares,
  // using bounded TRF. ALGLIB minlm VJ is a bounded least-squares
  // alternative, not an identical implementation of SciPy TRF.
  // Supply analytic residual Jacobians to avoid finite-difference FFTs.
  constexpr int parameterCount = 2;
  real_1d_array x, lowerBounds, upperBounds, scale;
  x.setlength(parameterCount);
  lowerBounds.setlength(parameterCount);
  upperBounds.setlength(parameterCount);
  scale.setlength(parameterCount);

  x[0] = std::clamp(itsFrac, fracMin, fracMax);
  x[1] = std::max(itsEpsFrac, minEpsFrac);
  lowerBounds[0] = fracMin;
  upperBounds[0] = fracMax;
  lowerBounds[1] = minEpsFrac;
  upperBounds[1] = std::numeric_limits<double>::infinity();
  scale[0] = 1.0;
  scale[1] = 1.0;

  // The Python residual vector consists of pixels where window != 0.
  // Store their coordinates once per minor-cycle iteration.
  std::vector<IPosition> pixels;
  for (Int j = blcDirty(1); j <= trcDirty(1); ++j)
    for (Int i = blcDirty(0); i <= trcDirty(0); ++i)
      if (insideCleanWindow(i, j))
        pixels.emplace_back(2, i, j);

  if (pixels.empty()) return;
  const ae_int_t m = static_cast<ae_int_t>(pixels.size());

  std::size_t diagnosticEval = 0;

  Asp2026LMContext context;
  context.evaluate = [this, &pixels, &diagnosticEval](
                         const real_1d_array& trial,
                         real_1d_array& fi,
                         real_2d_array* jac) {
    const Float f = static_cast<Float>(trial[0]);
    const Float ef = static_cast<Float>(trial[1]);
    const Matrix<Float> trialComponent = makeTanhComponent(f, ef);
    const Matrix<Float> trialConvolved = convolveComponent(trialComponent);

    for (ae_int_t k = 0; k < static_cast<ae_int_t>(pixels.size()); ++k) {
      const Int i = pixels[k][0], j = pixels[k][1];
      fi[k] = double((*itsDirty)(i,j)) - double(trialConvolved(i,j));
    }
    if (!jac) return;

    const Float peak = std::abs(itsPeakValue);
    const Float threshold = f * peak;
    // Python: eps_frac * abs(max_res) + 1e-12, not max(...,1e-12).
    const Float epsilon = ef * peak + 1.0e-12f;

    // Analytic derivative of the component wrt frac.
    Matrix<Float> componentDerivative(itsDirty->shape(), 0.0f);
    for (Int j = blcDirty(1); j <= trcDirty(1); ++j)
      for (Int i = blcDirty(0); i <= trcDirty(0); ++i) {
        if (!insideCleanWindow(i,j)) continue;
        const Float r = (*itsDirty)(i,j);
        const Float znum = itsPeakSign*r - threshold;
        if (znum < 0.0f) continue;
        const Float z = znum/epsilon;
        const Float t = std::tanh(z);
        const Float sech2 = 1.0f-t*t;
        componentDerivative(i,j) = -r*sech2*peak/epsilon;
      }
    Matrix<Float> convolvedDerivative = convolveComponent(componentDerivative);
    for (ae_int_t k = 0; k < static_cast<ae_int_t>(pixels.size()); ++k)
      (*jac)[k][0] = -double(convolvedDerivative(pixels[k][0],pixels[k][1]));

    // Analytic derivative wrt eps_frac.
    componentDerivative = 0.0f;
    for (Int j = blcDirty(1); j <= trcDirty(1); ++j)
      for (Int i = blcDirty(0); i <= trcDirty(0); ++i) {
        if (!insideCleanWindow(i,j)) continue;
        const Float r = (*itsDirty)(i,j);
        const Float znum = itsPeakSign*r - threshold;
        if (znum < 0.0f) continue;
        const Float z = znum/epsilon;
        const Float t = std::tanh(z);
        componentDerivative(i,j) = -r*(1.0f-t*t)*z*peak/epsilon;
      }
    convolvedDerivative = convolveComponent(componentDerivative);
    for (ae_int_t k = 0; k < static_cast<ae_int_t>(pixels.size()); ++k)
      (*jac)[k][1] = -double(convolvedDerivative(pixels[k][0],pixels[k][1]));
  
    // debug info: objective and gradient of the scalar least-squares
    // objective F = sum(fi[k]^2).
    //
    // ALGLIB receives the residual vector fi and its Jacobian J.
    // Therefore, grad(F) = 2 * J^T * fi.

    /*double objective = 0.0;
    double gradFrac = 0.0;
    double gradEpsFrac = 0.0;

    for (ae_int_t k = 0;
         k < static_cast<ae_int_t>(pixels.size()); ++k)
    {
      const double residual = fi[k];

      objective += residual * residual;
      gradFrac += 2.0 * residual * (*jac)[k][0];
      gradEpsFrac += 2.0 * residual * (*jac)[k][1];
    }

    ++diagnosticEval;

    std::cerr << std::scientific << std::setprecision(8)
              << "[AspClean2026 optimizer]"
              << " eval=" << diagnosticEval
              << " frac=" << f
              << " epsFrac=" << ef
              << " objective=" << objective
              << " gradFrac=" << gradFrac
              << " gradEpsFrac=" << gradEpsFrac
              << " absGradFrac=" << std::abs(gradFrac)
              << " absGradEpsFrac=" << std::abs(gradEpsFrac)
              << '\n';* /
  };

  minlmstate state;
  minlmcreatevj(m, x, state);
  minlmsetbc(state, lowerBounds, upperBounds);
  minlmsetscale(state, scale);
  // This ALGLIB version uses (state, epsx, maxits), not the older
  // (state, epsg, epsf, epsx, maxits) signature. Not SciPy TRF tolerances.
  minlmsetcond(state, 1.0e-3, 10);
  
  minlmreport report;

  minlmoptimize(state, asp2026Residuals, asp2026ResidualsJacobian,
                nullptr, &context);
  minlmresults(state, x, report);

  /*std::cerr << std::scientific << std::setprecision(8)
            << "[AspClean2026 optimizer final]"
            << " terminationtype=" << report.terminationtype
            << " frac=" << x[0]
            << " epsFrac=" << x[1]
            << '\n';* /

  itsFrac = static_cast<Float>(x[0]);
  itsEpsFrac = static_cast<Float>(x[1]);

  itsResidualComponent = makeTanhComponent(itsFrac, itsEpsFrac);
  itsConvolvedComponent = convolveComponent(itsResidualComponent);
  
  // -------- BEGIN DIAGNOSTICS --------

  double dotRC = 0.0;
  double normC2 = 0.0;
  double normR2 = 0.0;

  const auto shape = itsDirty->shape();
  const int nx = shape[0];
  const int ny = shape[1];

  for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {

          if (!insideCleanWindow(i, j)) continue;

          const double r = (*itsDirty)(i, j);
          const double c = itsConvolvedComponent(i, j);

          dotRC += r * c;
          normC2 += c * c;
          normR2 += r * r;
      }
  }

  const double deltaF = -2.0 * dotRC + normC2;

  const double deltaF_gain =
      -2.0 * itsGain * dotRC
      + itsGain * itsGain * normC2;

  std::cerr << "[AspClean2026 component diagnostics]"
            << " frac=" << itsFrac
            << " epsFrac=" << itsEpsFrac
            << " dotRC=" << dotRC
            << " normC2=" << normC2
            << " F0=" << normR2
            << " deltaF=" << deltaF
            << " deltaF_gain=" << deltaF_gain
            << std::endl;

  // -------- END DIAGNOSTICS --------

  // The component is already amplitude-valued. Do not multiply it by
  // itsFrac again during the model/residual update.
  itsStrengthOptimum = max(abs(itsResidualComponent));

}*/


void AspClean2026::optimizeComponent()
{
  extractResidualPeak();
  std::cerr << "[AspClean2026 peak]"
          << " peakValue=" << itsPeakValue
          << " peakSign=" << itsPeakSign
          << " peakPosition=" << itsPositionOptimum
          << std::endl;
  

  itsResidualComponent.resize(itsDirty->shape());
  itsResidualComponent = 0.0f;
  itsConvolvedComponent.resize(itsDirty->shape());
  itsConvolvedComponent = 0.0f;

  if (itsPeakValue == 0.0f)
  {
    itsStrengthOptimum = 0.0f;
    return;
  }

  constexpr int parameterCount = 2;
  real_1d_array x, lowerBounds, upperBounds, scale;
  x.setlength(parameterCount);
  lowerBounds.setlength(parameterCount);
  upperBounds.setlength(parameterCount);
  scale.setlength(parameterCount);

  x[0] = std::clamp(itsFrac, fracMin, fracMax);
  x[1] = std::max(itsEpsFrac, minEpsFrac);

  lowerBounds[0] = fracMin;
  upperBounds[0] = fracMax;
  lowerBounds[1] = minEpsFrac;
  upperBounds[1] = std::numeric_limits<double>::infinity();

  scale[0] = 1.0;
  scale[1] = 1.0;

  // The Python residual vector consists of pixels where window != 0.
  // Store their coordinates once per minor-cycle iteration.
  std::vector<IPosition> pixels;

  for (Int j = blcDirty(1); j <= trcDirty(1); ++j)
    for (Int i = blcDirty(0); i <= trcDirty(0); ++i)
      if (insideCleanWindow(i, j))
        pixels.emplace_back(2, i, j);

  if (pixels.empty())
  {
    itsStrengthOptimum = 0.0f;
    return;
  }

  const ae_int_t m = static_cast<ae_int_t>(pixels.size());

  // ------------------------------------------------------------
  // Component diagnostics.
  //
  // All diagnostics use the same residual image and pixel list
  // as the optimizer.
  //
  // F0          = sum(R^2)
  // dotRC       = sum(R * convolvedComponent)
  // normC2      = sum(convolvedComponent^2)
  //
  // deltaF      = -2*dotRC + normC2
  // deltaF_gain = -2*gain*dotRC + gain^2*normC2
  //
  // Negative deltaF means the squared residual is reduced.
  // ------------------------------------------------------------

  auto printComponentDiagnostics =
      [this, &pixels](const char* label,
                      double frac,
                      double epsFrac,
                      const Matrix<Float>& convolvedComponent)
  {
    double dotRC = 0.0;
    double normC2 = 0.0;
    double normR2 = 0.0;

    for (const auto& pixel : pixels)
    {
      const Int i = pixel[0];
      const Int j = pixel[1];

      const double r = static_cast<double>((*itsDirty)(i, j));
      const double c =
          static_cast<double>(convolvedComponent(i, j));

      dotRC += r * c;
      normC2 += c * c;
      normR2 += r * r;
    }

    const double gain = static_cast<double>(itsGain);

    const double deltaF =
        -2.0 * dotRC + normC2;

    const double deltaF_gain =
        -2.0 * gain * dotRC
        + gain * gain * normC2;

    const double F_full = normR2 + deltaF;
    const double F_gain = normR2 + deltaF_gain;

    const double corr =
        (normR2 > 0.0 && normC2 > 0.0)
            ? dotRC / std::sqrt(normR2 * normC2)
            : 0.0;

    const double optimalGain =
        (normC2 > 0.0)
            ? dotRC / normC2
            : 0.0;

    std::cerr << std::scientific
              << std::setprecision(8)
              << "[AspClean2026 " << label << "]"
              << " frac=" << frac
              << " epsFrac=" << epsFrac
              << " dotRC=" << dotRC
              << " normC2=" << normC2
              << " corr=" << corr
              << " optimalGain=" << optimalGain
              << " F0=" << normR2
              << " deltaF=" << deltaF
              << " deltaF_gain=" << deltaF_gain
              << " F_full=" << F_full
              << " F_gain=" << F_gain
              << std::endl;
  };

  // ------------------------------------------------------------
  // INITIAL DIAGNOSTICS
  //
  // Evaluate the component produced by the optimizer's starting
  // parameters, before ALGLIB changes those parameters.
  // ------------------------------------------------------------

  const double initialFrac = x[0];
  const double initialEpsFrac = x[1];

  {
    const Matrix<Float> initialComponent =
        makeTanhComponent(
            static_cast<Float>(initialFrac),
            static_cast<Float>(initialEpsFrac));

    const Matrix<Float> initialConvolved =
        convolveComponent(initialComponent);

    const Int pi = itsPositionOptimum[0];
    const Int pj = itsPositionOptimum[1];

    std::cerr << "[AspClean2026 INITIAL peak response]"
              << " peakSign=" << itsPeakSign
              << " residual=" << (*itsDirty)(pi, pj)
              << " component=" << initialComponent(pi, pj)
              << " convolved=" << initialConvolved(pi, pj)
              << std::endl;

    printComponentDiagnostics(
        "INITIAL",
        initialFrac,
        initialEpsFrac,
        initialConvolved);
  }

  // ------------------------------------------------------------
  // Original optimizer context and analytic Jacobian.
  // ------------------------------------------------------------

  std::size_t diagnosticEval = 0;

  Asp2026LMContext context;

  context.evaluate =
      [this, &pixels, &diagnosticEval](
          const real_1d_array& trial,
          real_1d_array& fi,
          real_2d_array* jac)
  {
    const Float f = static_cast<Float>(trial[0]);
    const Float ef = static_cast<Float>(trial[1]);

    const Matrix<Float> trialComponent =
        makeTanhComponent(f, ef);

    const Matrix<Float> trialConvolved =
        convolveComponent(trialComponent);

    for (ae_int_t k = 0;
         k < static_cast<ae_int_t>(pixels.size()); ++k)
    {
      const Int i = pixels[k][0];
      const Int j = pixels[k][1];

      fi[k] =
          double((*itsDirty)(i, j))
          - double(trialConvolved(i, j));
    }

    if (!jac)
      return;

    const Float peak = std::abs(itsPeakValue);
    const Float threshold = f * peak;

    // Python:
    // eps_frac * abs(max_res) + 1e-12
    const Float epsilon = ef * peak + 1.0e-12f;

    // ----------------------------------------------------------
    // Analytic derivative wrt frac.
    // ----------------------------------------------------------

    Matrix<Float> componentDerivative(
        itsDirty->shape(), 0.0f);

    for (Int j = blcDirty(1); j <= trcDirty(1); ++j)
      for (Int i = blcDirty(0); i <= trcDirty(0); ++i)
      {
        if (!insideCleanWindow(i, j))
          continue;

        const Float r = (*itsDirty)(i, j);
        const Float znum =
            itsPeakSign * r - threshold;

        if (znum < 0.0f)
          continue;

        const Float z = znum / epsilon;
        const Float t = std::tanh(z);
        const Float sech2 = 1.0f - t * t;

        componentDerivative(i, j) =
            -r * sech2 * peak / epsilon;
      }

    Matrix<Float> convolvedDerivative =
        convolveComponent(componentDerivative);

    for (ae_int_t k = 0;
         k < static_cast<ae_int_t>(pixels.size()); ++k)
    {
      (*jac)[k][0] =
          -double(convolvedDerivative(
              pixels[k][0], pixels[k][1]));
    }

    // ----------------------------------------------------------
    // Analytic derivative wrt eps_frac.
    // ----------------------------------------------------------

    componentDerivative = 0.0f;

    for (Int j = blcDirty(1); j <= trcDirty(1); ++j)
      for (Int i = blcDirty(0); i <= trcDirty(0); ++i)
      {
        if (!insideCleanWindow(i, j))
          continue;

        const Float r = (*itsDirty)(i, j);
        const Float znum =
            itsPeakSign * r - threshold;

        if (znum < 0.0f)
          continue;

        const Float z = znum / epsilon;
        const Float t = std::tanh(z);

        componentDerivative(i, j) =
            -r * (1.0f - t * t)
            * z * peak / epsilon;
      }

    convolvedDerivative =
        convolveComponent(componentDerivative);

    for (ae_int_t k = 0;
         k < static_cast<ae_int_t>(pixels.size()); ++k)
    {
      (*jac)[k][1] =
          -double(convolvedDerivative(
              pixels[k][0], pixels[k][1]));
    }
  };

  // ------------------------------------------------------------
  // Original ALGLIB optimization.
  // ------------------------------------------------------------

  minlmstate state;
  minlmcreatevj(m, x, state);
  minlmsetbc(state, lowerBounds, upperBounds);
  minlmsetscale(state, scale);
  minlmsetcond(state, 1.0e-3, 10);

  minlmreport report;

  minlmoptimize(
      state,
      asp2026Residuals,
      asp2026ResidualsJacobian,
      nullptr,
      &context);

  minlmresults(state, x, report);

  itsFrac = static_cast<Float>(x[0]);
  itsEpsFrac = static_cast<Float>(x[1]);

  // ------------------------------------------------------------
  // FINAL COMPONENT
  // ------------------------------------------------------------

  itsResidualComponent =
      makeTanhComponent(itsFrac, itsEpsFrac);

  itsConvolvedComponent =
      convolveComponent(itsResidualComponent);

  // debug output
  const Int pi = itsPositionOptimum[0];
  const Int pj = itsPositionOptimum[1];

  std::cerr << "[AspClean2026 FINAL peak response]"
            << " residual=" << (*itsDirty)(pi,pj)
            << " component=" << itsResidualComponent(pi,pj)
            << " convolved=" << itsConvolvedComponent(pi,pj)
            << std::endl;

  // ------------------------------------------------------------
  // FINAL DIAGNOSTICS
  //
  // Evaluate the component that will actually be used by
  // updateModelAndResidual().
  // ------------------------------------------------------------

  printComponentDiagnostics(
      "FINAL",
      itsFrac,
      itsEpsFrac,
      itsConvolvedComponent);

  std::cerr << "[AspClean2026 optimizer final]"
            << " terminationtype=" << report.terminationtype
            << " frac=" << itsFrac
            << " epsFrac=" << itsEpsFrac
            << std::endl;

  itsStrengthOptimum =
      max(abs(itsResidualComponent));
}


void AspClean2026::updateModelAndResidual(
    Matrix<Float>& model,
    const IPosition& /*support*/,
    const Matrix<Float>& /*scale*/,
    const Matrix<Complex>& /*scaleXfr*/)
{
  if (itsPeakValue == 0.0f)
    return;

  // Eq. 30 and 31: the optimized component already contains 
  // amplitude, so only the gain is applied here.
  model += itsGain * itsResidualComponent;
  (*itsDirty) -= itsGain * itsConvolvedComponent;
}

void AspClean2026::makeInitScaleImage(
    Matrix<Float>& iscale,
    const Float& /*scaleSize*/)
{
  iscale = 0.0f;
}

void AspClean2026::makeScaleImage(
    Matrix<Float>& iscale,
    const Float& /*scaleSize*/,
    const Float& /*amp*/,
    const IPosition& /*center*/)
{
  // for AspClean2026 this is simply the optimized residual component.
  AlwaysAssert(itsResidualComponent.shape() == iscale.shape(), AipsError);
  iscale = itsResidualComponent;
}

IPosition AspClean2026::componentSupport(
    const IPosition& imageShape) const
{
  // The support is the adaptive residual component itself, not a Gaussian
  // support derived from itsPsfWidth or itsOptimumScaleSize.
  return imageShape;
}

Bool AspClean2026::useLegacyStrengthLogic() const
{ return false; }

Bool AspClean2026::useHogbomFallback() const
{ return false; }

Bool AspClean2026::useLegacyScaleImagePath() const
{ return false; }

} // namespace casa
