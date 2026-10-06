//# SDAlgorithmAspClean2026.h: Definition for SDAlgorithmAspClean2026

#ifndef SYNTHESIS_SDALGORITHMASPCLEAN2026_H
#define SYNTHESIS_SDALGORITHMASPCLEAN2026_H

#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/BasicSL/String.h>

#include <synthesis/ImagerObjects/SDAlgorithmBase.h>
#include <synthesis/MeasurementEquations/AspClean2026.h>

namespace casa {

class SDAlgorithmAspClean2026 : public SDAlgorithmBase
{
public:
  explicit SDAlgorithmAspClean2026(bool isSingle = true);
  ~SDAlgorithmAspClean2026() override;

protected:
  void takeOneStep(
      casacore::Float loopgain,
      casacore::Int cycleNiter,
      casacore::Float cycleThreshold,
      casacore::Float& peakresidual,
      casacore::Float& modelflux,
      casacore::Int& iterdone) override;

  void initializeDeconvolver() override;
  void finalizeDeconvolver() override;

private:
  casacore::Array<casacore::Float> itsMatPsf;
  casacore::Array<casacore::Float> itsMatResidual;
  casacore::Array<casacore::Float> itsMatModel;
  casacore::Array<casacore::Float> itsMatMask;

  AspClean2026 itsCleaner;
  casacore::Bool itsMCsetup;
  bool itsIsSingle;
};

} // namespace casa

#endif
