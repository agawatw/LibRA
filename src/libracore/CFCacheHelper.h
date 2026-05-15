// -*- C++ -*-
//# CFCacheHelper.h: Definition for CFCacheHelper functions
//# Copyright (C) 2021
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
//# Correspondence concerning this should be addressed as follows:
//#        Postal address: National Radio Astronomy Observatory
//#                        1003 Lopezville Road,
//#                        Socorro, NM - 87801, USA
//#
//# $Id$
#ifndef LIBRACORE_CFCACHEHELPER_H
#define LIBRACORE_CFCACHEHELPER_H
#include <casacore/measures/Measures/Stokes.h>
#include <casacore/coordinates/Coordinates/CoordinateSystem.h>
#include <casacore/images/Images/ImageInterface.h>
#include <casacore/images/Images/PagedImage.h>
#include <casacore/images/Images/TempImage.h>
#include <casacore/tables/TaQL/ExprNode.h>

#include <synthesis/TransformMachines2/FTMachine.h>
#include <synthesis/TransformMachines2/AWVisResampler.h>
#include <synthesis/TransformMachines2/CFStore2.h>
#include <synthesis/TransformMachines2/AWVisResamplerHPG.h>
#include <synthesis/TransformMachines2/AWProjectWBFTHPG.h>
#include <synthesis/TransformMachines2/MakeCFArray.h>
#include <synthesis/TransformMachines2/ThreadCoordinator.h>

void makeCFS_inmemory(CountedPtr<casa::refim::CFStore2> cfs2_l,
		      CountedPtr<casa::refim::CFStore2> cfswt2_l,
		      const bool& psTerm,
		      const bool& aTerm,
		      const bool& conjBeams);

void fillCFC_inmemory(DataBase& db,
		      //CountedPtr<refim::ConvolutionFunction>& awcf_l
		      refim::ConvolutionFunction>& awcf_l,
		      const TempImage<Complex>& cgrid,
		      const Vector<double>& uvScale,
		      const Vector<double>& uvOffset);
#endif 
