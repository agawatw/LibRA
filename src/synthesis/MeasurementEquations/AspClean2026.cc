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

#include <synthesis/MeasurementEquations/AspClean2026.h>

// for alglib
#include <synthesis/MeasurementEquations/objfunc_alglib.h>

using namespace alglib;

using namespace casacore;
using namespace std;

namespace casa {

namespace {

template <typename Objective>
void runAspClean2026LBFGS(
    minlbfgsstate& state,
    real_1d_array& x,
    minlbfgsreport& rep,
    Objective objective,
    void* ptrParam)
{
  alglib::minlbfgsoptimize(state, objective, NULL, ptrParam);
  minlbfgsresults(state, x, rep);
}

struct AspClean2026LBFGSContext
{
  AspClean2026* cleaner;
};

constexpr Float fracMin = 0.60f;
constexpr Float fracMax = 0.99f;
constexpr Float minEpsFrac = 1.0e-6f;
constexpr double maxTransform = 20.0;

} // namespace

AspClean2026::AspClean2026()
    : itsFrac(0.8f),
      itsEpsFrac(0.1f),
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

  // The 2026 path constructs its component directly from the residual and
  // never consumes the legacy initial-scale image or its transform.  Do not
  // allocate and FFT an image-sized zero matrix here: on large Hummbee images
  // those unused allocations can exhaust memory during initialization.
  fft = FFTServer<Float, Complex>(psfShape_p);
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
      if (std::abs(value) > std::abs(peak))
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
      std::max(epsFrac * peakAmplitude, 1.0e-12f);

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

Matrix<Float> AspClean2026::convolveComponent(
    const Matrix<Float>& component)
{
  Matrix<Complex> componentXfr;
  fft.fft0(componentXfr, component);

  Matrix<Complex> work((*itsXfr) * componentXfr);
  Matrix<Float> convolved(component.shape());

  fft.fft0(convolved, work, false);
  fft.flip(convolved, false, false);

  // Match the pixel-origin correction used by the legacy Asp update. Without
  // this shift the minor cycle subtracts a PSF response displaced from the
  // component written into the model. A following major cycle then computes
  // the true response from that model and restores the apparent peak.
  IPosition nullnull(2, 0);
  const IPosition support(convolved.shape());
  Matrix<Float> unshifted(convolved.shape());
  unshifted.assign_conforming(convolved);
  if (itsdimensionsareeven)
  {
    Matrix<Float> destination = convolved(nullnull + 1, support - 1);
    destination.assign_conforming(unshifted(nullnull, support - 2));
  }
  else
  {
    Matrix<Float> destination = convolved(nullnull + 2, support - 1);
    destination.assign_conforming(unshifted(nullnull, support - 3));
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
  const Float epsilon = std::max(epsFrac * peakAmplitude, 1.0e-12f);

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

void AspClean2026::optimizeComponent()
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

  // Reuse the existing lbfgs in the Asp base class. The Asp base
  // optimizes the Gaussian amplitude and scale. The AspClean2026 uses
  // the same optimizer through a new callback.
  constexpr int parameterCount = 2;
  constexpr int historySize = 2;

  real_1d_array x;
  x.setlength(parameterCount);
  x[0] = inverseFrac(std::clamp(itsFrac, fracMin + 1.0e-5f,
                                fracMax - 1.0e-5f));
  x[1] = inverseEpsFrac(std::max(itsEpsFrac, minEpsFrac));

  real_1d_array scale;
  scale.setlength(parameterCount);
  scale[0] = 1.0;
  scale[1] = 1.0;

  minlbfgsstate state;
  // The callback supplies the objective and its analytic gradient. ALGLIB
  // requires the history size M not to exceed the problem dimension N.
  minlbfgscreate(parameterCount, historySize, x, state);
  minlbfgssetcond(state, 1.0e-5, 1.0e-6, 1.0e-5, 25);
  minlbfgssetscale(state, scale);

  minlbfgsreport report;
  AspClean2026LBFGSContext context{this};

  runAspClean2026LBFGS(state, x, report, &AspClean2026::lbfgsObjective, &context);

  itsFrac = transformedFrac(x[0]);
  itsEpsFrac = transformedEpsFrac(x[1]);

  itsResidualComponent =
      makeTanhComponent(itsFrac, itsEpsFrac);
  itsConvolvedComponent =
      convolveComponent(itsResidualComponent);

  // The component is already amplitude-valued. Do not multiply it by
  // itsFrac again during the model/residual update.
  itsStrengthOptimum = max(abs(itsResidualComponent));
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
{
  return false;
}

Bool AspClean2026::useHogbomFallback() const
{
  return false;
}

Bool AspClean2026::useLegacyScaleImagePath() const
{
  return false;
}

} // namespace casa
