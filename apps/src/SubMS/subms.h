// # Copyright (C) 2021
// # Associated Universities, Inc. Washington DC, USA.
// #
// # This library is free software; you can redistribute it and/or modify it
// # under the terms of the GNU Library General Public License as published by
// # the Free Software Foundation; either version 2 of the License, or (at your
// # option) any later version.
// #
// # This library is distributed in the hope that it will be useful, but WITHOUT
// # ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// # FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
// # License for more details.is
// #
// # You should have received a copy of the GNU Library General Public License
// # along with this library; if not, write to the Free Software Foundation,
// # Inc., 675 Massachusetts Ave, Cambridge, MA 02139, USA.
// #
// # Correspondence concerning this should be addressed as follows:
// #        Postal address: National Radio Astronomy Observatory
// #                        1003 Lopezville Road,
// #                        Socorro, NM - 87801, USA
// #
// # $Id$

#ifndef SUBMS_H
#define SUBMS_H

#include <casacore/casa/aips.h>
#include <casacore/ms/MSSel/MSSelection.h>
#include <casacore/ms/MSSel/MSSelectionError.h>
#include <casacore/ms/MSSel/MSSelectionTools.h>
#include <msvis/MSVis/SubMS.h>

using namespace std;
using namespace casa;
using namespace casacore;

/**
 * @brief Creates a subset of a MeasurementSet using the CASA SubMS engine.
 *
 * Selection is driven by the SPW string syntax (e.g. "0:5~58" for SPW 0
 * channels 5-58), so no separate nchan/start vectors are needed.
 * Channel averaging is controlled by chanStep: step=1 is no averaging,
 * step=N averages every N channels across all selected SPWs.
 * Time averaging is controlled by integ (seconds); pass -1 to disable.
 * The combineStr parameter is forwarded to makeSubMS and controls whether
 * scans or SPWs are combined in the output (e.g. "scan,spw").
 * If overwrite=true and OutMSBuf already exists on disk it will be removed
 * before writing. If overwrite=false and OutMSBuf exists an error is thrown.
 * If verbose=true a summary of all selection parameters is printed before
 * the split runs.
 *
 * @param MSNBuf       Input Measurement Set name (required).
 * @param OutMSBuf     Output Measurement Set name (required).
 * @param WhichColStr  Data column to write: "data", "model", "corrected", "all".
 * @param fieldStr     Field selection (e.g. "0", "NGC*", "*").
 * @param timeStr      Time range selection (e.g. "2010/01/01/00:00~01:00").
 * @param spwStr       SPW selection with optional channel range (e.g. "0:5~58").
 * @param baselineStr  Baseline/antenna selection (e.g. "VA01&VA02").
 * @param scanStr      Scan selection (e.g. "1~5").
 * @param arrayStr     Sub-array selection.
 * @param uvdistStr    UV distance selection (e.g. ">1klambda").
 * @param taqlStr      Arbitrary TaQL expression for additional row selection.
 * @param corrStr      Correlation/polarisation selection (e.g. "RR,LL").
 * @param intentStr    Scan intent selection (e.g. "CALIBRATE_*").
 * @param obsStr       Observation ID selection.
 * @param combineStr   Combine scans/spws in output (e.g. "scan,spw", or "").
 * @param integ        Time bin size in seconds for averaging (-1 = no averaging).
 * @param chanStep     Channel averaging step applied to all selected SPWs (default 1 = no averaging).
 * @param overwrite    Remove OutMSBuf if it already exists (default false).
 * @param verbose      Print selection summary before running (default false).
 */
void SubMS_func(const string& MSNBuf="",      const string& OutMSBuf="",
                const string& WhichColStr="data",
                const string& fieldStr="*",   const string& timeStr="",
                const string& spwStr="*",     const string& baselineStr="",
                const string& scanStr="",     const string& arrayStr="",
                const string& uvdistStr="",   const string& taqlStr="",
                const string& corrStr="",     const string& intentStr="",
                const string& obsStr="",      const string& combineStr="",
                const double  integ=-1.0,
                const int     chanStep=1,
                const bool    overwrite=false,
                const bool    verbose=false);

/**
 * @brief Parafeed CLI user interface for subms.
 *
 * @param restart      Whether to re-enter the UI loop.
 * @param argc         Argument count.
 * @param argv         Argument vector.
 * @param interactive  Whether to run interactively.
 * @param MSNBuf       Input MS name.
 * @param OutMSBuf     Output MS name.
 * @param WhichColStr  Data column name.
 * @param fieldStr     Field selection string.
 * @param timeStr      Time range string.
 * @param spwStr       SPW selection string.
 * @param baselineStr  Baseline selection string.
 * @param scanStr      Scan selection string.
 * @param arrayStr     Sub-array selection string.
 * @param uvdistStr    UV distance selection string.
 * @param taqlStr      TaQL expression string.
 * @param corrStr      Correlation selection string.
 * @param intentStr    Scan intent selection string.
 * @param obsStr       Observation ID selection string.
 * @param combineStr   Combine string forwarded to makeSubMS.
 * @param integ        Time bin size in seconds (-1 = no averaging).
 * @param chanStep     Channel averaging step (1 = no averaging).
 * @param overwrite    Overwrite existing output MS.
 * @param verbose      Print selection summary before running.
 */
void UI(bool restart, int argc, char **argv, bool interactive,
        string& MSNBuf,       string& OutMSBuf,
        string& WhichColStr,
        string& fieldStr,     string& timeStr,
        string& spwStr,       string& baselineStr,
        string& scanStr,      string& arrayStr,
        string& uvdistStr,    string& taqlStr,
        string& corrStr,      string& intentStr,
        string& obsStr,       string& combineStr,
        double& integ,
        int&    chanStep,
        bool&   overwrite,
        bool&   verbose);

#endif
