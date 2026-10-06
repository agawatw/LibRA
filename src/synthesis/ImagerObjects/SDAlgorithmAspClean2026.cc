//# SDAlgorithmAspClean2026.cc: Implementation of SDAlgorithmAspClean2026
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
//# $Id$


#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/OS/HostInfo.h>

#include <components/ComponentModels/SkyComponent.h>
#include <components/ComponentModels/ComponentList.h>
#include <casacore/images/Images/TempImage.h>
#include <casacore/images/Images/SubImage.h>
#include <casacore/images/Regions/ImageRegion.h>
#include <casacore/casa/OS/File.h>
#include <casacore/lattices/LEL/LatticeExpr.h>
#include <casacore/lattices/Lattices/TiledLineStepper.h>
#include <casacore/lattices/Lattices/LatticeStepper.h>
#include <casacore/lattices/Lattices/LatticeIterator.h>
#include <synthesis/TransformMachines/StokesImageUtil.h>
#include <casacore/coordinates/Coordinates/StokesCoordinate.h>
#include <casacore/casa/Exceptions/Error.h>
#include <casacore/casa/BasicSL/String.h>
#include <casacore/casa/Utilities/Assert.h>
#include <casacore/casa/OS/Directory.h>
#include <casacore/tables/Tables/TableLock.h>

#include<synthesis/ImagerObjects/SIMinorCycleController.h>
#include <sstream>

#include <casacore/casa/Logging/LogMessage.h>
#include <casacore/casa/Logging/LogIO.h>
#include <casacore/casa/Logging/LogSink.h>

#include <casacore/casa/System/Choice.h>
#include <msvis/MSVis/StokesVector.h>
#include <synthesis/ImagerObjects/SDAlgorithmAspClean2026.h>
#include <casacore/casa/Arrays/Matrix.h>
#include <casacore/casa/Quanta/Quantum.h>
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
  itsCleaner.setInitScaleMasks(mask);
  itsCleaner.setaspcontrol(0, 0, 0, Quantity(0.0, "%"));

  Matrix<Float> residual;
  residual.reference(itsMatResidual);
  itsCleaner.setDirty(residual);

  // AspClean2026 builds its component directly from the residual, but the
  // shared Asp minor-cycle code still requires a valid scale list.
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