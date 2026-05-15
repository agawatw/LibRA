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

#include <subms.h>
#include <filesystem>

void SubMS_func(const string& MSNBuf,      const string& OutMSBuf,
                const string& WhichColStr,
                const string& fieldStr,    const string& timeStr,
                const string& spwStr,      const string& baselineStr,
                const string& scanStr,     const string& arrayStr,
                const string& uvdistStr,   const string& taqlStr,
                const string& corrStr,     const string& intentStr,
                const string& obsStr,      const string& combineStr,
                const double  integ,
                const int     chanStep,
                const bool    overwrite,
                const bool    verbose)
{
    if (MSNBuf.empty())
        throw AipsError("SubMS: input MS name not set");
    if (OutMSBuf.empty())
        throw AipsError("SubMS: output MS name not set");
    if (chanStep < 1)
        throw AipsError("SubMS: chanstep must be >= 1");

    // Overwrite guard
    if (std::filesystem::exists(OutMSBuf)) {
        if (!overwrite)
            throw AipsError("SubMS: output MS already exists: " + OutMSBuf
                            + " -- set overwrite=true to replace it");
        std::filesystem::remove_all(OutMSBuf);
    }

    if (verbose) {
        LogIO os(LogOrigin("SubMS", "SubMS_func"));
        os << LogIO::NORMAL
           << "SubMS selection summary:" << endl
           << "  vis        = " << MSNBuf       << endl
           << "  outms      = " << OutMSBuf     << endl
           << "  datacolumn = " << WhichColStr  << endl
           << "  field      = " << fieldStr     << endl
           << "  spw        = " << spwStr       << endl
           << "  time       = " << timeStr      << endl
           << "  baseline   = " << baselineStr  << endl
           << "  scan       = " << scanStr      << endl
           << "  array      = " << arrayStr     << endl
           << "  uvdist     = " << uvdistStr    << endl
           << "  taql       = " << taqlStr      << endl
           << "  correlation= " << corrStr      << endl
           << "  intent     = " << intentStr    << endl
           << "  observation= " << obsStr       << endl
           << "  combine    = " << combineStr   << endl
           << "  integ      = " << integ        << endl
           << "  chanstep   = " << chanStep     << endl
           << LogIO::POST;
    }

    MeasurementSet ms(MSNBuf, TableLock(TableLock::AutoNoReadLocking));

    SubMS splitter(ms);

    String CspwStr(spwStr),       CfieldStr(fieldStr),
           CbaselineStr(baselineStr), CscanStr(scanStr),
           CuvdistStr(uvdistStr), CtaqlStr(taqlStr),
           CarrayStr(arrayStr),   CcorrStr(corrStr),
           CintentStr(intentStr), CobsStr(obsStr);

    Bool ok = splitter.setmsselect(CspwStr, CfieldStr, CbaselineStr, CscanStr,
                                   CuvdistStr, CtaqlStr,
                                   Vector<Int>(1, chanStep),
                                   CarrayStr, CcorrStr, CintentStr, CobsStr);

    if (!ok)
        throw AipsError("SubMS: setmsselect failed -- check selection strings");

    String CtimeStr(timeStr);
    splitter.selectTime(static_cast<Double>(integ), CtimeStr);

    String OutMSName(OutMSBuf), WhichCol(WhichColStr);
    String CcombineStr(combineStr);

    if (!splitter.makeSubMS(OutMSName, WhichCol, Vector<Int>(1, 0), CcombineStr))
        throw AipsError("SubMS: makeSubMS failed writing " + OutMSBuf);
}
