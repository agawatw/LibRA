// -*- C++ -*-
//# CFCacheHeler.cc: Helper functions for making and managing CFCache
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

#include <CFCacheHelper.h>
//
//--------------------------------------------------------------------------
//
CountedPtr<refim::PolOuterProduct> setPOP(vi::VisBuffer2 &vb2,
					  Vector<casacore::Stokes::StokesTypes> visPolMap,
					  Vector<int> polMap,
					  std::string &stokes, std::string &mType)
{
  CountedPtr<refim::PolOuterProduct> pop_l = new PolOuterProduct;
  
  //------------------------a mess----------------------------------------------------
  Vector<Int> intpolmap(visPolMap.nelements());
  for (uInt kk=0; kk < intpolmap.nelements(); ++kk){
    intpolmap[kk]=Int(visPolMap[kk]);
  }
  pop_l->initCFMaps(intpolmap, polMap);
  
  PolMapType polMat, polIndexMat, conjPolMat, conjPolIndexMat;
  Vector<Int> visPol(vb2.correlationTypes());
  polMat = pop_l->makePolMat(visPol,polMap);
  polIndexMat = pop_l->makePol2CFMat(visPol,polMap);
  
  conjPolMat = pop_l->makeConjPolMat(visPol,polMap);
  conjPolIndexMat = pop_l->makeConjPol2CFMat(visPol,polMap);
  
  return pop_l;
}

void makeCFS_inmemory(CountedPtr<casa::refim::CFStore2> cfs2_l,
		      CountedPtr<casa::refim::CFStore2> cfswt2_l,
		      const bool& psTerm,
		      const bool& aTerm,
		      const bool& conjBeams)
{
  //
  // mode="fillcf" case.  The list of CFs in the CFC are
  // expected to be "blank CFs" with all the necessary meta
  // information to fill them.  Nothing else other than the
  // CFC and the list of CFs is necessary.  An already filled
  // CF (the IsFilled=1 entry in CFS*/miscInfo.rec) will be
  // left untouched.
  //
  Vector<double> dummyUVScale,uvOffset;
  Matrix<double> dummyvbFreqSel;
  // !!!!!!!!!!!!!!!! What should the cfCacheName be?!!!!!!!!!!!!!
  string cfCacheName; 
  AWConvFunc::makeConvFunction2(cfCacheName,
				dummyUVScale, uvOffset,	dummyvbFreqSel,
				*cfs2_l,*cfswt2_l,
				psTerm,	aTerm, conjBeams);


  // Report some stats.
  Double memUsed=cfs2_l->memUsage();
  String unit(" KB");
  memUsed = (Int)(memUsed/1024.0+0.5);
  if (memUsed > 1024) {memUsed /=1024; unit=" MB";}
}

void fillCFC_inmemory(DataBase& db,
		      //CountedPtr<refim::ConvolutionFunction>& awcf_l
		      CountedPtr<casa::refim::CFStore2> cfs2_l,
		      CountedPtr<casa::refim::CFStore2> cfswt2_l,
		      refim::ConvolutionFunction& awcf_l,
		      const TempImage<Complex>& cgrid,
		      int& nW, float& pa, float& dpa,
		      const Vector<double>& uvScale,
		      const Vector<double>& uvOffset,
		      std::string mType,
		      std::string stokes)
{
  //-------------------------------------------------------------------------------------------------
  // Instantiate the PolOuterProduce object which encapsulates the
  // polarization maps from VB to Muller (which is CFs) to the
  // image-plane polarizations requested.  Set the necessary maps in
  // the ConvolutionFunction object (awcf_l variable).
  //
  //Replace this with the actual polarization parser (this commented
  //existed earlier and I (SB) don't know what it may mean).
  Vector<int> polMap;
  Vector<casacore::Stokes::StokesTypes> visPolMap;//{0,1,2,3};

  refim::SynthesisUtils::matchPol(*(db.vb_l),cgrid.coordinates(),polMap,visPolMap);
  // Initialize pop to have the right values

  CountedPtr<refim::PolOuterProduct> pop_p = setPOP(*(db.vb_l), visPolMap, polMap, stokes, mType);

  // Vector<int> spwidList, fieldidList;
  // Vector<double> spwRefFreqList;
  // spwidList      = db.spwidList;
  // fieldidList    = db.fieldidList;
  // spwRefFreqList = db.spwRefFreqList;
  Matrix<double> mssFreqSel;
  mssFreqSel  = db.msSelection.getChanFreqList(NULL,true);
  awcf_l.setPolMap(polMap);
  //  awcf_l->setSpwSelection(spwChanSelFlag_p);
  // Replace mssFreqSel with a filtered version of it, filtered
  // for the SPW IDs in the stl::vector db.spwidList.tovector().
  mssFreqSel.assign(filterByFirstColumn(mssFreqSel,db.spwidList));
  awcf_l.setSpwFreqSelection(mssFreqSel);

  // Get the PA from the MS/VB if UI setting is outside the valid
  // range for PA [-180, +180].
  if (abs(pa) > 180.0) pa=getPA(*(db.vb_l));
  constexpr double D2R=M_PI/180.0;
  pa *= D2R; dpa *= D2R;
  //
  // This currently makes both CF and WTCF.  It will also fill the
  // CFs (serially) for fillCF=true.  For fillCF=false (default
  // value here), it will produce an "empty CFC" which has
  // single-pixel CFs with the necessary information in the
  // header/MiscInfo to fill them later. This allows
  // parallelization of the compute intensive step (filling the
  // CFs) which is highly parallelizable and scales well.
  //
  Matrix<Double> vbFreqSelection ;
  bool fillCF=true;
  awcf_l.makeConvFunction(cgrid , *(db.vb_l), nW,
			  pop_p, pa, dpa, uvScale, uvOffset,
			  vbFreqSelection,
			  *cfs2_l, *cfswt2_l,
			  fillCF);
  //
  // AWConvFunc::makeConvFunction() does not make the memory
  // model (CFStore) persistent.  So save the contents of the
  // CFStore on the disk.
  //
  // [07Jan2024] In the dryrun mode, only the meta info is
  // written as casacore::Records converted to
  // casacore::Tables.  Writing these with multi-threadings
  // seems to be work.  The bool parameter is therefore set to
  // true (it is false in the default interface).
  // cfs2_l->makePersistent(cfCacheObj_l->getCacheDir().c_str(),"","", Quantity(pa,"rad"),Quantity(dpa,"rad"),0,0,true);
  // cfswt2_l->makePersistent(cfCacheObj_l->getCacheDir().c_str(),"","WT",Quantity(pa,"rad"),Quantity(dpa,"rad"),0,0,true);

  // cfs2_l->makePersistent(cfCacheObj_l->getCacheDir().c_str(),"","", true);
  // cfswt2_l->makePersistent(cfCacheObj_l->getCacheDir().c_str(),"","WT",true);
}
