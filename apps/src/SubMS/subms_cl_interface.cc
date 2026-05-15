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

#include <cl.h>
#include <clinteract.h>
#include <subms.h>

#define RestartUI(Label)  {if(clIsInteractive()) {goto Label;}}

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
        bool&   verbose)
{
    clSetPrompt(interactive);

    if (!restart) {
        BeginCL(argc, argv);
        clInteractive(0);
    }

    try {
        int i;
        i=1; clgetSValp("ms",           MSNBuf,      i);
        i=1; clgetSValp("outms",        OutMSBuf,    i);
             clgetFullValp("datacolumn", WhichColStr);
             clSetOptions("datacolumn",  {"data", "model", "corrected", "all"});
             clgetFullValp("field",       fieldStr);
             clgetFullValp("time",        timeStr);
        i=1; clgetValp("integ",          integ,       i);
             clgetFullValp("spw",         spwStr);
        i=1; clgetIValp("chanstep",      chanStep,    i);
             clgetFullValp("baseline",    baselineStr);
             clgetFullValp("scan",        scanStr);
             clgetFullValp("array",       arrayStr);
             clgetFullValp("uvdist",      uvdistStr);
             dbgclgetFullValp("taql",     taqlStr);
             clgetFullValp("correlation", corrStr);
             clgetFullValp("intent",      intentStr);
             clgetFullValp("observation", obsStr);
             clgetFullValp("combine",     combineStr);
        i=1; clgetBValp("overwrite",     overwrite,   i);
        i=1; clgetBValp("verbose",       verbose,     i);
        EndCL();

        string mesgs;
        if (MSNBuf.empty())   mesgs += "Input MS name not set.\n";
        if (OutMSBuf.empty()) mesgs += "Output MS name not set.\n";
        if (chanStep < 1)     mesgs += "chanstep must be >= 1.\n";
        if (!mesgs.empty())   clThrowUp(mesgs, "###Fatal", CL_FATAL);
    }
    catch (clError& x) {
        x << x << endl;
        if (x.Severity() == CL_FATAL) { throw; }
    }
}

#ifndef SUBMS_LIBRARY_BUILD
int main(int argc, char **argv)
{
    string MSNBuf, OutMSBuf, WhichColStr="data",
           fieldStr="*", timeStr, spwStr="*",
           baselineStr, uvdistStr, taqlStr,
           scanStr, arrayStr, corrStr,
           intentStr, obsStr, combineStr;
    double integ     = -1.0;
    int    chanStep  = 1;
    bool   overwrite = false;
    bool   verbose   = false;
    bool   restartUI = false;

 RENTER:
    try {
        MSNBuf = OutMSBuf = timeStr = baselineStr = uvdistStr =
            taqlStr = scanStr = arrayStr = corrStr =
            intentStr = obsStr = combineStr = "";
        WhichColStr = "data";
        fieldStr    = "*";
        spwStr      = "*";
        integ       = -1.0;
        chanStep    = 1;
        overwrite   = false;
        verbose     = false;

        bool interactive = true;

        UI(restartUI, argc, argv, interactive,
           MSNBuf, OutMSBuf, WhichColStr,
           fieldStr, timeStr, spwStr, baselineStr,
           scanStr, arrayStr, uvdistStr, taqlStr,
           corrStr, intentStr, obsStr, combineStr,
           integ, chanStep, overwrite, verbose);
        restartUI = false;

        SubMS_func(MSNBuf, OutMSBuf, WhichColStr,
                   fieldStr, timeStr, spwStr, baselineStr,
                   scanStr, arrayStr, uvdistStr, taqlStr,
                   corrStr, intentStr, obsStr, combineStr,
                   integ, chanStep, overwrite, verbose);
    }
    catch (clError& x) {
        x << x.what() << endl;
        restartUI = true;
    }
    catch (MSSelectionError& x) {
        cerr << "###MSSelectionError: " << x.getMesg() << endl;
        restartUI = true;
    }
    catch (AipsError& x) {
        cerr << "###AipsError: " << x.getMesg() << endl;
        restartUI = true;
    }
    if (restartUI) RestartUI(RENTER);
}
#endif // SUBMS_LIBRARY_BUILD
