# Copyright (C) 2021, 2026
# Associated Universities, Inc. Washington DC, USA.
#
# This library is free software; you can redistribute it and/or modify it
# under the terms of the GNU Library General Public License as published by
# the Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This library is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
# License for more details.is
#
# You should have received a copy of the GNU Library General Public License
# along with this library; if not, write to the Free Software Foundation,
# Inc., 675 Massachusetts Ave, Cambridge, MA 02139, USA.
#
# Correspondence concerning this should be addressed as follows:
#        Postal address: National Radio Astronomy Observatory
#                        1003 Lopezville Road,
#                        Socorro, NM - 87801, USA
#
# $Id$


#!/bin/bash

input_file=$1
EXT_PREFIX="wis"

OMP_NUM_THREADS=1

# Default parameters
OSGjob=false
imagename="wistest0"
ncycle=4
LIBRAHOME=/home/dhruva/disk1/sanjay/Packages/LibRA/libra
logdir=${PWD}/wisimaging0.log
RUNAPP=../scripts/runapp.sh
libraBIN=${LIBRAHOME}/install_test/bin

MODE_PSF="wispsf"
MODE_WEIGHT="wisweight"
MODE_RESIDUAL="wisresidual"
MODE_PREDICT="predict"
mkdir -p ${logdir}

function runapp_l() {
    # The "!" at the end of the input .def file does an exlusive load
    # of parameter from the .def file that aren't set (in an
    # interactive session of via the commandline
    $@ help=def,${input_file}"!" 
    }


griddingAPP=${libraBIN}/roadrunner
deconvolutionAPP=${libraBIN}/hummbee
normalizationAPP=${libraBIN}/dale


# restart=1 will pick up files from an existing htclean run and continue CLEANing
# set start_index to the current number of imaging cycles + 1 (next imaging cycle after restart)
# set GPUENGINE=0 to run casatools based gridding (same as serialized htclean)
restart=0 
start_index=1

# Begin execution block

echo "runapp.sh" {$RUNAPP}
echo "Using gridding application:      "${griddingAPP}
echo "Using deconvolution application: "${deconvolutionAPP}
echo "Using normalization application: "${normalizationAPP}
echo ""

if [ "$restart" -eq "0" ]
then
    # Initialize the MODEL_DATA Column
    # roadrunner mode=predict datacolumn=model modelimage=""
    echo "Initializing the MODEL_DATA ..."
    runapp_l ${griddingAPP} mode=${MODE_PREDICT} datacolumn=model imagename= modeimagename= \
             2>| ${logdir}/predict0.log

    # generate weight
    echo "Making the initial ${MODE_WEIGHT} ..."
    runapp_l ${griddingAPP} mode=${MODE_WEIGHT} imagename=${imagename}.weight complexgrid=${imagename}_gridcf.vis0 \
             2>| ${logdir}/${MODE_WEIGHT}0.log

    # # generare residual
    echo "Making the initial ${MODE_RESIDUAL} ..."
    runapp_l ${griddingAPP} mode=${MODE_RESIDUAL} imagename=${imagename}.residual complexgrid=${imagename}_gridres.vis0 \
             2>| ${logdir}/${MODE_RESIDUAL}0.log

    # # generate psf
    echo "Making the initial ${MODE_PSF} ..."
    runapp_l ${griddingAPP} mode=${MODE_PSF} imagename=${imagename}.psf complexgrid=${imagename}_gridpsf.vis0 \
             2>| ${logdir}/${MODE_PSF}0.log
    
    # # normalize the residual
    echo "Normalizing the initial ${MODE_RESIDUAL} ..."
    runapp_l ${normalizationAPP} imtype=residual imagename=${imagename} \
             2>| ${logdir}/norm_res0.log
    cp -r ${imagename}.residual ${imagename}.dirty.cycle0

    # echo ${RUNAPP} ${normalizationAPP} normalize ${imagename} ${input_file} ${logdir} -t residual -c 0
    # ${RUNAPP} ${normalizationAPP} normalize ${imagename} ${input_file} ${logdir} -t residual -c 0
    # cp -r ${imagename}.residual ${imagename}.dirty.cycle0

    # # normalize the PSF and make primary beam
    echo "Normalizing the initial ${MODE_PSF} and PB ..."
    runapp_l ${normalizationAPP} imtype=psf imagename=${imagename} computepb=1\
             2>| ${logdir}/norm_psf0.log
    cp -r ${imagename}.psf ${imagename}.psf.cycle0
    cp -r ${imagename}.weight ${imagename}.weight.cycle0

    # ${RUNAPP} ${normalizationAPP} normalize ${imagename} ${input_file} ${logdir} -t psf
    # cp -r ${imagename}.psf ${imagename}.psf.cycle0
    # cp -r ${imagename}.weight ${imagename}.weight.cycle0
exit 0;

else
    echo "Doing only the update step..."
fi

i=$start_index

while [ ! -f stopIMCycles ] && [ "${i}" -lt "${ncycle}" ]
do
    echo "Cycle: "${i}

    # run hummbee for updateModel deconvolution iterations
    #------------------------------------------------------------------------------------
    echo "Making model update ${i} ..."
    runapp_l ${deconvolutionAPP} imagename=${imagename} mode=deconvolve\
             2>| ${logdir}/model_update${i}.log
    # ${RUNAPP} ${deconvolutionAPP} deconvolve ${imagename} ${input_file} ${logdir} -c ${i}
    #------------------------------------------------------------------------------------

    # work around to fix the NOOP in dale
    #------------------------------------------------------------------------------------
    sed -i 's/SubType =  normalized/SubType =/' ${imagename}.model/table.info

    # run dale to divide model by weights
    runapp_l ${normalizationAPP}  imtype=model imagename=${imagename} \
             2>| ${logdir}/model_normalize${i}.log
    #------------------------------------------------------------------------------------
             
    # ${RUNAPP} ${normalizationAPP} normalize ${imagename} ${input_file} ${logdir} -t model -c ${i}

    # save model at each cycle for debugging and remove unwanted images
    cp -r ${imagename}.divmodel ${imagename}.model.cycle${i}
    rm -rf ${imagename}.psf
    rm -rf ${imagename}.residual
    rm -rf ${imagename}.weight



    #------------------------------------------------------------------------------------
    # Fill the MODEL_DATA Column
    # roadrunner mode=predict datacolumn=model modelimage=<LATEST MODEL FROM HUMMBEE>
    echo "Making ${i} MODEL_DATA ..."
    runapp_l ${griddingAPP} mode=${MODE_PREDICT} datacolumn=model imagename=${imagename}.divmodel \
             2>| ${logdir}/predict${i}.log
    # echo -e "load psf_predict.def \ngo" | ${griddingAPP} help=dbg 2>| model_predict.log
    #eval "${griddingApp} help=def,psf_predict.def" 2>| model_predict.log
    #------------------------------------------------------------------------------------

    #------------------------------------------------------------------------------------
    # Make weight and save gridded vis
    echo "Making the ${i} ${MODE_WEIGHT} ..."
    runapp_l ${griddingAPP} mode=${MODE_WEIGHT} imagename=${imagename}.weight \
             complexgrid=${imagename}_gridcf.vis${i} \
             2>| ${logdir}/${MODE_WEIGHT}${i}.log

    # ${RUNAPP} ${griddingAPP} ${MODE_WEIGHT} ${imagename} ${input_file} ${logdir} -c ${i}
    # mv ${imagename}.${MODE_WEIGHT} ${imagename}.weight
    # cp -r ${imagename}_gridv.vis ${imagename}_gridcf.vis${i}
    #------------------------------------------------------------------------------------


    #------------------------------------------------------------------------------------
    # Make Redisual and save gridded vis
    echo "Making the ${i} ${MODE_RESIDUAL} ..."
    runapp_l ${griddingAPP} mode=${MODE_RESIDUAL} imagename=${imagename}.residual modelimagename=${imagename}.divmodel \
             complexgrid=${imagename}_gridres.vis${i} \
             2>| ${logdir}/${MODE_RESIDUAL}${i}.log

#     ${RUNAPP} ${griddingAPP} ${MODE_RESIDUAL} ${imagename} ${input_file} ${logdir} -m ${imagename}.divmodel -c ${i}
# #    ${RUNAPP} ${griddingAPP} ${MODE_RESIDUAL} ${imagename} ${input_file} ${logdir} -c ${i}
#     mv ${imagename}.${MODE_RESIDUAL} ${imagename}.residual
#     cp -r ${imagename}_gridv.vis ${imagename}_gridres.vis${i}

    # normalize the new residual 
    echo "Normalizing the ${i} ${MODE_RESIDUAL} ..."
    runapp_l ${normalizationAPP} imtype=residual imagename=${imagename} \
             2>| ${logdir}/norm_res${i}.log
    # ${RUNAPP} ${normalizationAPP} normalize ${imagename} ${input_file} ${logdir} -t residual -c ${i}
    #------------------------------------------------------------------------------------


    #------------------------------------------------------------------------------------
    # Make PSF and save gridded vis
    echo "Making the ${i} ${MODE_PSF} ..."
    runapp_l ${griddingAPP} mode=${MODE_PSF} imagename=${imagename}.psf complexgrid=${imagename}_gridpsf.vis${i} \
             2>| ${logdir}/${MODE_PSF}${i}.log
    
    # ${RUNAPP} ${griddingAPP} ${MODE_PSF} ${imagename} ${input_file} ${logdir} -c ${i}
    # mv ${imagename}.${MODE_PSF} ${imagename}.psf
    # cp -r ${imagename}_gridv.vis ${imagename}_gridpsf.vis{i}

    # normalize the new PSF
    echo "Normalizing the ${i} ${MODE_PSF} and PB ..."
    runapp_l ${normalizationAPP} imtype=psf imagename=${imagename} computepb=1\
             2>| ${logdir}/norm_psf0.log
    # ${RUNAPP} ${normalizationAPP} normalize ${imagename} ${input_file} ${logdir} -t psf -c ${i}
    #------------------------------------------------------------------------------------

    # save residual at each cycle for debugging
    cp -r ${imagename}.residual ${imagename}.dirty.cycle${i}
    cp -r ${imagename}.psf ${imagename}.psf.cycle${i}
    cp -r ${imagename}.weight ${imagename}.weight.cycle${i}

    i=$((i+1))
done

# run hummbee for restore
#${RUNAPP} ${deconvolutionAPP} restore ${imagename} ${input_file} ${logdir}

if [ "${OSGjob}" = "True" ]
then
    createTar weight sumwt psf mask model residual image image.pbcor
    lsdir
fi
