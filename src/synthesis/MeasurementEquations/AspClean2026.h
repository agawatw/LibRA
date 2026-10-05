//# AspClean2026.h: Residual-patch AspClean2026 minor cycle
//# Copyright (C) 1996,1997,1998,1999,2000,2001,2002,2003
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
//#        Internet email: aips2-request@nrao.edu.
//#        Postal address: AIPS++ Project Office
//#                        National Radio Astronomy Observatory
//#                        520 Edgemont Road
//#                        Charlottesville, VA 22903-2475 USA
//#
//#
//# $Id: AspClean2026.h Genie H. 2020-04-06 <mhsieh@nrao.edu $

#ifndef SYNTHESIS_ASPCLEAN2026_H
#define SYNTHESIS_ASPCLEAN2026_H

#include <synthesis/MeasurementEquations/AspMatrixCleaner.h>

namespace casa {

// AspClean2026 only has one component from the residual around its peak. 
// It doesn't require initial scales. The tanh-mapped residual patch, after 
// optimization, is the component used for both the model and residual updates.

class AspClean2026 : public AspMatrixCleaner
{
public:
  AspClean2026();
  ~AspClean2026() override;

protected:
  void setInitScaleXfrs(const casacore::Float width) override;

  casacore::Bool setInitScaleMasks(
      const casacore::Matrix<casacore::Float>& mask,
      const casacore::Float& maskThreshold) override;

  std::vector<casacore::Float> getActiveSetAspen(
      const float peakres) override;

  void makeInitScaleImage(
      casacore::Matrix<casacore::Float>& iscale,
      const casacore::Float& scaleSize) override;

  void makeScaleImage(
      casacore::Matrix<casacore::Float>& iscale,
      const casacore::Float& scaleSize,
      const casacore::Float& amp,
      const casacore::IPosition& center) override;

  void optimizeComponent() override;

  void updateModelAndResidual(
    casacore::Matrix<casacore::Float>& model,
    const casacore::IPosition& support,
    const casacore::Matrix<casacore::Float>& scale,
    const casacore::Matrix<casacore::Complex>& scaleXfr) override;

  casacore::IPosition componentSupport(
      const casacore::IPosition& imageShape) const override;

  casacore::Bool useLegacyStrengthLogic() const override;
  casacore::Bool useHogbomFallback() const override;
  casacore::Bool useLegacyScaleImagePath() const override;

private:
  using AspMatrixCleaner::itsSwitchedToHogbom;
  
  casacore::Bool insideCleanWindow(casacore::Int i, casacore::Int j) const;

  void extractResidualPeak();

  casacore::Matrix<casacore::Float> makeTanhComponent(
      casacore::Float frac,
      casacore::Float epsFrac) const;

  casacore::Matrix<casacore::Float> convolveComponent(
      const casacore::Matrix<casacore::Float>& component) const;

  casacore::Float objective(
      casacore::Float frac,
      casacore::Float epsFrac,
      casacore::Matrix<casacore::Float>* component = nullptr,
      casacore::Matrix<casacore::Float>* convolved = nullptr) const;

  static void lbfgsObjective(
      const alglib::real_1d_array& x,
      double& func,
      void* ptr);

  static casacore::Float transformedFrac(double x);
  static casacore::Float transformedEpsFrac(double x);
  static double inverseFrac(casacore::Float frac);
  static double inverseEpsFrac(casacore::Float epsFrac);

  //casacore::IPosition itsResidualComponentSupport;
  casacore::Matrix<casacore::Float> itsResidualComponent;
  casacore::Matrix<casacore::Float> itsConvolvedComponent;

  casacore::Float itsFrac;
  casacore::Float itsEpsFrac;
  casacore::Float itsPeakValue;
  casacore::Float itsPeakSign;
};

} // namespace casa

#endif

