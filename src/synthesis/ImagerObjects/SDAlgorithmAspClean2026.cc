//# SDAlgorithmAspClean2026.cc: Implementation of SDAlgorithmAspClean2026

#include <synthesis/ImagerObjects/SDAlgorithmAspClean2026.h>

#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/Matrix.h>
#include <casacore/casa/Exceptions/Error.h>
#include <casacore/casa/Logging/LogIO.h>
#include <casacore/casa/Quanta/Quantum.h>
#include <casacore/casa/Utilities/Assert.h>

#include <synthesis/ImagerObjects/SIImageStore.h>

#include <vector>

using namespace casacore;

namespace casa {

SDAlgorithmAspClean2026::SDAlgorithmAspClean2026(bool isSingle)
    : SDAlgorithmBase(),
      itsMatPsf(),
      itsMatResidual(),
      itsMatModel(),
      itsMatMask(),
      itsCleaner(),
      itsMCsetup(true),
      itsIsSingle(isSingle)
{
  itsAlgorithmName = String("asp2026");
}

SDAlgorithmAspClean2026::~SDAlgorithmAspClean2026() = default;

void SDAlgorithmAspClean2026::initializeDeconvolver()
{
  LogIO os(LogOrigin("SDAlgorithmAspClean2026", "initializeDeconvolver", WHERE));
  AlwaysAssert(static_cast<bool>(itsImages), AipsError);

  itsImages->residual()->get(itsMatResidual, true);
  itsImages->model()->get(itsMatModel, true);
  itsImages->psf()->get(itsMatPsf, true);
  itsImages->mask()->get(itsMatMask, true);

  // A cube can present a different PSF for each channel, whereas an MFS
  // minor cycle can reuse this setup.
  if (itsMCsetup)
  {
    Matrix<Float> psf(itsMatPsf);
    itsCleaner.setPsf(psf);
    itsCleaner.setInitScaleXfrs(0.0f);
    itsCleaner.ignoreCenterBox(true);

    if (itsIsSingle)
      itsMCsetup = false;
  }

  Matrix<Float> mask(itsMatMask);
  itsCleaner.setInitScaleMasks(mask, 0.99f);
  itsCleaner.setaspcontrol(0, 0, 0, Quantity(0.0, "%"));

  Matrix<Float> residual;
  residual.reference(itsMatResidual);
  itsCleaner.setDirty(residual);

  // AspClean2026 builds its component directly from the residual, but the
  // shared Asp minor-cycle driver still requires a valid scale list.
  std::vector<Float> scales{0.0f};
  itsCleaner.defineAspScales(scales);
}

void SDAlgorithmAspClean2026::takeOneStep(
    Float loopgain,
    Int cycleNiter,
    Float cycleThreshold,
    Float& peakresidual,
    Float& modelflux,
    Int& iterdone)
{
  LogIO os(LogOrigin("SDAlgorithmAspClean2026", "takeOneStep", WHERE));

  itsCleaner.setaspcontrol(
      cycleNiter, loopgain, Quantity(cycleThreshold, "Jy"), Quantity(0.0, "%"));

  Matrix<Float> model;
  model.reference(itsMatModel);

  itsCleaner.startingIteration(0);
  const Int result = itsCleaner.aspclean(model);
  iterdone = itsCleaner.numberIterations();

  if (result == -1)
    os << LogIO::WARN << "AspClean2026 minor cycle stopped in point mode"
       << LogIO::POST;
  else if (result == -2)
    os << LogIO::WARN << "AspClean2026 minor cycle stopped after a large-scale divergence"
       << LogIO::POST;
  else if (result == -3)
    os << LogIO::WARN << "AspClean2026 minor cycle stopped because it is diverging"
       << LogIO::POST;

  itsMatResidual = itsCleaner.getterResidual();
  peakresidual = itsCleaner.getterPeakResidual();
  modelflux = sum(itsMatModel);
}

void SDAlgorithmAspClean2026::finalizeDeconvolver()
{
  itsImages->residual()->put(itsMatResidual);
  itsImages->model()->put(itsMatModel);
}

} // namespace casa
