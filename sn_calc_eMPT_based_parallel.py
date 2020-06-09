#Written by MVM on 28 August 2018

#Want a function that will return the emission line S/N
#as a function of slit position, size, PA, all that jazz

#Edited by Emma Curtis Lake October 2018
#To work with command line arguments, with multiple processes and produce outputs for multiple lines


import numpy as np
#import matplotlib.pyplot as plt
from pandeia.engine.calc_utils import build_default_calc
from pandeia.engine.perform_calculation import perform_calculation
import pandeia.engine
from pathos.multiprocessing import ProcessingPool
import itertools
from astropy.io import fits
import time
from collections import OrderedDict, defaultdict
from astropy.table import Table
import sys
import argparse
import pickle
import os

#We are only using g140m/f070lp but pandeia only provides simulated spectra for this setup up to 1.2 um - perhaps because zeroth order contaminates spectra beyond this?
#I don't want a gap in the S/N coverage so I'll also include g140m/g100lp spectra
nConfig = 4
grating = ['g140m','g140m','g235m','g395m']
filter = ['f070lp','f100lp','f170lp','f290lp']

filterWlDict = {'clear':{'low':0.7,'high':5.1},\
                'f070lp':{'low':0.7,'high':1.2},\
                'f100lp':{'low':1.1,'high':1.8},\
                'f170lp':{'low':1.7,'high':3.1},\
                'f290lp':{'low':3.0,'high':5.1}}


#for the Sersic profile, since the ETC prefers not to use the half-light radius
#def bn(n):
#    return 2.*n-1./3.+4./(405.*n)+46./(25515.*n**2.)+131./(1148175.*n**3.)-2194697./(30690717750.*n**4.)
#def majorminor(n,re,ellip):
#    #ellipticity is defined as (major-minor)/major
#    scale_length=re/(bn(n)**n)
#    #scale length is the circularized radius, i.e. r_scale = sqrt(a*b)
#    major_axis=scale_length/np.sqrt(1.-ellip)
#    minor_axis=scale_length*np.sqrt(1.-ellip)
#    return (major_axis,minor_axis)





def sn_array(inputs,disperser='prism',filt='clear', ngroup=19, nint=2, nexp=36):
    #E - This function doesn't take account of variablility accross the fov of the detector/msa, and doesn't account for the placement of the detector gap within the spectrum.
    #E - This function also assumes that the line flux follows the surface brightness profile of the continuum, and that both have simple Sersic profiles
    xoff=inputs['xoff']
    yoff=inputs['yoff']
    configuration=build_default_calc('jwst','nirspec','msa')
    configuration['configuration']['instrument']['disperser']=disperser
    configuration['configuration']['instrument']['filter']=filt
    #E - I *think* shutter location is just so that the detector gap is placed in the correct place
    configuration['configuration']['instrument']['shutter_location']='q3_345_20'#'q4_345_20'
    #exposure specifications from DEEP APT file
    configuration['configuration']['detector']['ngroup']=ngroup
    configuration['configuration']['detector']['nint']=nint
    #PRISM
    configuration['configuration']['detector']['nexp']=nexp
    configuration['configuration']['detector']['readmode']='nrsirs2'
    configuration['strategy']['dithers'][0]['on_source'] = inputs['onSource']
    configuration['configuration']['instrument']['slitlet_shape'] = inputs['slitletShape']

    #default configuration['configuration'] has 1x3 shutter config
    scene = {}
    if (inputs['r_circ'] > 0):
        size=inputs['r_circ']
        sersic=inputs['sersic_n']
        ellip=1.-inputs['axis_ratio']
        pa=inputs['position_angle']
        #pandeia now takes half-light radii as input
        #major_axis,minor_axis=majorminor(sersic,size,ellip)
        major_axis = inputs['re_maj']
        minor_axis = inputs['re_maj']*inputs['axis_ratio']
        scene['position'] = {'x_offset':xoff, 'y_offset': yoff, 'orientation': pa, 'position_parameters':['x_offset','y_offset','orientation']}
        scene['shape']={'geometry':'sersic','major': major_axis,'minor':minor_axis,'sersic_index':sersic}
    else:
        pa=inputs['position_angle']
        #this is the dummy trigger to go for a point source
        scene['position'] = {'x_offset':xoff, 'y_offset': yoff, 'orientation': pa, 'position_parameters':['x_offset','y_offset','orientation']}
        scene['shape']={'geometry':'point'}



    scene['spectrum']={'name':'Flat Source','spectrum_parameters':['sed','normalization']}
    scene['spectrum']['sed'] = {'sed_type':'flat','unit':'fnu'}
    #normalized to have NO CONTINUUM
    scene['spectrum']['normalization'] = {'type': 'at_lambda','norm_wave': 1.4, 'norm_waveunit':'um', 'norm_flux': 0, 'norm_fluxunit':'mjy'}
    
    emission_line_array = []
    for i,f in enumerate(inputs['flux']):
        flux=f
        wave=inputs['wave'][i]
        if wave > filterWlDict[filt]['low'] and wave < filterWlDict[filt]['high'] and flux > 0:
            emission_line={}
            emission_line['emission_or_absorption']='emission'
            emission_line['center']=wave
            print 'wl: ', wave, flux
            #TODO: check units...
            emission_line['strength']=flux
            emission_line['profile']='gaussian'
            #assume the line is basically unresolved, i.e. 50 km/s
            emission_line['width']=50.#50.
            emission_line_array.append(emission_line)
    
#        print emission_line_array
#        sys.exit()
    if len(emission_line_array) > 0:
        print len(emission_line_array)
        scene['spectrum']['lines']=emission_line_array
        print emission_line_array
        print '** here'
        configuration['scene'][0]=scene
        report=perform_calculation(configuration)
        print 'success'
#            plt.figure()
        #plt.plot(report['1d']['extracted_flux'][0],report['1d']['extracted_flux'][1],color='blue')
#        plt.plot(report['1d']['sn'][0]/(1+inputs['redshift']),report['1d']['sn'][1],color='red')
#        plt.plot(report['1d']['total_flux'][0],report['1d']['total_flux'][1],'--',color='black')
#        plt.plot(report['1d']['extracted_flux_plus_bg'][0],report['1d']['extracted_flux_plus_bg'][1],'--',color='blue')
#        plt.plot(report['1d']['extracted_noise'][0],report['1d']['extracted_noise'][1],'--',color='orange')
#        plt.show()
        #integrate over the whole line
        #(well, just the part that has some significant S/N, i.e. > 0.1)
#        print '***', np.max(report['1d']['sn'][1])
#        print 'extracted_flux: ', report['1d']['extracted_flux']
#        print report['1d'].keys()
#        print 'wave: ', wave, flux
#        print report.keys()
#        sys.exit()

        #If just one position is provided for all exposures, the S/N can be output, but if you want to combine different exposures with different source offsets
        #If a sequence of dithers is provided as input, the integrated flux and noise estimates are needed to combine them at the end.
    integrated_sn=[]
    integrated_flux = []
    integrated_noise = []
#        plt.figure()
#        plt.plot(report['1d']['extracted_flux'][0]/(inputs['redshift']+1.),report['1d']['sn'][1],color='red')
#        plt.scatter(report['1d']['extracted_flux'][0]/(inputs['redshift']+1.),report['1d']['sn'][1],color='red')
    for i,w in enumerate(inputs['wave']):
        if w > filterWlDict[filt]['low'] and w < filterWlDict[filt]['high'] and inputs['flux'][i] > 0.: #we don't need to add in lines outside wavelength range of filter
            doublet = inputs['lineDoublet'][i]
            if doublet > 0:
                doubletIdx = np.where(inputs['lineDoublet'] == doublet)[0]
            else:
                doubletIdx = [i]
            if doubletIdx[0] == 0:
                line_lower_w = 0.7
            else:
                line_lower_w = np.sum(inputs['wave'][doubletIdx[0]-1:doubletIdx[0]+1])/2.
            if np.max(doubletIdx) == len(inputs['wave'])-1:
                line_upper_w = 7
            else:
                line_upper_w = np.sum(inputs['wave'][np.max(doubletIdx):np.max(doubletIdx)+2])/2.
            #need to isolate pixels contributing to the line
            #print w, line_lower_w, line_upper_w, inputs['lineDoublet']
            tempIdx = np.where((report['1d']['sn'][0] >= line_lower_w) & (report['1d']['sn'][0] <= line_upper_w) & \
                               (report['1d']['sn'][1] > 5e-3))[0]
            tempIdxHigh = np.where(report['1d']['sn'][0] > w)[0]
            idxCtr = tempIdxHigh[0]-1
                
            if len(tempIdx) == 0:
                integrated_sn.append(0.)
                integrated_flux.append(0.)
                integrated_noise.append(-99)
            else:
                #make sure not losing flux when lines start to become blended, always want a minimum of 3 pixels either side of the pixel containing the centre of the line
                if np.min(tempIdx) > idxCtr-3:
                    newIdx = list(np.fromiter((x for x in range(idxCtr-3,np.min(tempIdx))),np.int))+list(tempIdx)
                    tempIdx = newIdx
                if np.max(tempIdx) < idxCtr+3:
                    newIdx = list(tempIdx)+list(np.fromiter((x for x in range(np.max(tempIdx)+1,idxCtr+4)),np.int))
                    tempIdx = newIdx
#                    print 'plotting', w, report['1d']['extracted_flux'][0][tempIdx], report['1d']['extracted_flux'][0][tempIdx]
#                    plt.plot(report['1d']['sn'][0][tempIdx]/(1+inputs['redshift']),report['1d']['sn'][1][tempIdx],color='black')
                temp_sn = np.sum(report['1d']['sn'][1][tempIdx]/np.sqrt(len(tempIdx)))
                temp_flux = np.sum(report['1d']['extracted_flux'][1][tempIdx]) #this is not currently optimal, as you are more likely to fit a Gaussian, but still
                temp_noise2 = np.sum(np.power(report['1d']['extracted_noise'][1][tempIdx],2))
                if not np.isfinite(temp_sn):
                    integrated_sn.append(0.)
                    integrated_flux.append(0.)
                    integrated_noise.append(-99)
                else:
                    integrated_sn.append(temp_sn)
                    integrated_flux.append(temp_flux)
                    integrated_noise.append(np.sqrt(temp_noise2))
        else:
            #wavelength of line outside of wavelength range of filter
            integrated_sn.append(0.)
            integrated_flux.append(0.)
            integrated_noise.append(-99)
    #np.sum(report['1d']['sn'][1])/np.sqrt(np.sum(report['1d']['sn'][1] > 1e-1))
    
        #plt.show()
#        sys.exit()

    return integrated_sn,integrated_flux,integrated_noise

def add_to_outputDict(sn,n,f,outputDict, linesInOrder, lineDoublet, writeOutput, colStr, aa, temp_flux = None, temp_sn = None, temp_noise = None):
    for i,line in enumerate(linesInOrder):
        if writeOutput[i] == 1:
            if lineDoublet[i] > 0:
                tempIdx = np.where(lineDoublet == lineDoublet[i])[0]
                if len(tempIdx) > 2:
                    print "error, only accounting for doublets, shouldn't have any more than two lines blended, ", len(tempIdx)
                    sys.exit()
                outputLable = linesInOrder[tempIdx[0]]+'_'+linesInOrder[tempIdx[1]]
                if temp_sn is not None and temp_flux is not None and temp_noise is not None:
                    temp_sn[outputLable].append(sn[i])
                    temp_flux[outputLable].append(f[i])
                    temp_noise[outputLable].append(n[i])
                outputDict[outputLable+'_sn'+colStr].append(sn[i])
                outputDict[outputLable+'_f'+colStr].append(f[i])
                outputDict[outputLable+'_n'+colStr].append(n[i])
            else:
                if temp_sn is not None and temp_flux is not None and temp_noise is not None:
                    temp_sn[line].append(sn[i])
                    temp_flux[line].append(f[i])
                    temp_noise[line].append(n[i])
                outputDict[line+'_sn'+colStr].append(sn[i])
                outputDict[line+'_f'+colStr].append(f[i])
                outputDict[line+'_n'+colStr].append(n[i])


def sn_user_spec(inputs, disperser = 'prism', filt = 'clear', ngroup = 19, nint = 2, nexp = 36):
    wl = inputs['wl'] #in microns
    spec = inputs['spec'] #in mJy
    xoff = inputs['xoff']
    yoff = inputs['yoff']
    
    configuration=build_default_calc('jwst','nirspec','msa')
    configuration['configuration']['instrument']['disperser']=disperser
    configuration['configuration']['instrument']['filter']=filt
    #E - I *think* shutter location is just so that the detector gap is placed in the correct place
    configuration['configuration']['instrument']['shutter_location']='q3_345_20'#'q4_345_20'
    #exposure specifications from DEEP APT file
    configuration['configuration']['detector']['ngroup']=ngroup
    configuration['configuration']['detector']['nint']=nint
    #PRISM
    configuration['configuration']['detector']['nexp']=nexp
    configuration['configuration']['detector']['readmode']='nrsirs2'
    configuration['strategy']['dithers'][0]['on_source'] = inputs['onSource']
    configuration['configuration']['instrument']['slitlet_shape'] = inputs['slitletShape']
    
    #default configuration['configuration'] has 1x3 shutter config
    scene = {}
    if (inputs['r_circ'] > 0):
        size=inputs['r_circ']
        sersic=inputs['sersic_n']
        ellip=1.-inputs['axis_ratio']
        pa=inputs['position_angle']
        major_axis = inputs['re_maj']
        minor_axis = inputs['re_maj']*inputs['axis_ratio']
        #major_axis,minor_axis=majorminor(sersic,size,ellip)
        scene['position'] = {'x_offset':xoff, 'y_offset': yoff, 'orientation': pa, 'position_parameters':['x_offset','y_offset','orientation']}
        scene['shape']={'geometry':'sersic','major': major_axis,'minor':minor_axis,'sersic_index':sersic}
    else:
        pa=inputs['position_angle']
        #this is the dummy trigger to go for a point source
        scene['position'] = {'x_offset':xoff, 'y_offset': yoff, 'orientation': pa, 'position_parameters':['x_offset','y_offset','orientation']}
        scene['shape']={'geometry':'point'}
    
    scene['spectrum'] = {}
    scene['spectrum']['name'] = "continuum_spectrum"
    scene['spectrum']['redshift'] = 0 #because otherwise it shifts the wavelength array...
    #it doesn't seem to do anything with the normalization of the
    #source spectrum, however?!
    tempIdx = np.where((wl >= filterWlDict[filt]['low']) & (wl <= filterWlDict[filt]['high']))[0]
    scene['spectrum']['sed'] = {'sed_type': 'input', 'spectrum': [wl[tempIdx], spec[tempIdx]], 'unit':'mJy'}
    scene['spectrum']['normalization'] = {}
    scene['spectrum']['normalization'] = {'type': 'none'}
    
    configuration['scene'][0]=scene
    report=perform_calculation(configuration)
    return report

def add_cont_sn_to_outputDict(report,outputDict,keyOutStr, regions, redshift):
    for key in regions.keys():
        #calculate average pixel S/N in regions defined in the input dictionary
        region_wl_low = regions[key]['low']*(1+redshift)
        region_wl_high = regions[key]['high']*(1+redshift)
        tempIdx = np.where((report['1d']['sn'][0] >= region_wl_low) & (report['1d']['sn'][0] <= region_wl_high))[0]
        if len(tempIdx) > 0:
            temp_sn = np.mean(report['1d']['sn'][1][tempIdx])
            if np.isfinite(temp_sn):
                outputDict[key+'_avg_sn_pp'+keyOutStr] = temp_sn
            else:
                outputDict[key+'_avg_sn_pp'+keyOutStr] = 0.
        else:
            outputDict[key+'_avg_sn_pp'+keyOutStr] = 0.

def obj_sn(inputDict):
    print 'here I am'
    outputDict = defaultdict(list)
    idx = inputDict['idx']
    id_cat = inputDict['id_cat']
    cat = inputDict['cat']
    ditherDict = inputDict['ditherDict']
    linesInOrder = inputDict['linesInOrder']
    writeOutput_prism = inputDict['writeOutput_prism']
    writeOutput_R1000 = inputDict['writeOutput_R1000']
    lineDoublet_prism = inputDict['lineDoublet_prism']
    lineDoublet_R1000 = inputDict['lineDoublet_R1000']
    args = inputDict['args']
    print args.R100lines, args.R1000lines

    calcArr = [args.ps, args.psOff, args.es, args.esOff]
    outputStr = ['_ps', '_ps_off', '_es', '_es_off']
    combineArr = [False, False, False, False]
    if not args.singleDither:
        combineArr = [False, True, False, True]
    #does not produce desired effect, not sure how to set this correctly...
    #onSourceArr = [[False,True,False],[False,False,True],[True,False,False]]
    onSourceArr = [[False,True,False],[False,False,True],[True,False,False]]
    slitletShapeArr = [[[0,-1],[0,0],[0,1]],[[0,-2],[0,-1],[0,0]],[[0,0],[0,1],[0,2]]]
    outputFile = args.resultsDir+str(id_cat[idx])+".p"
    fileTest = os.path.isfile(outputFile)
    if fileTest == False:
        time_start = time.time()
        #random x and y offset in arcsec
        aa=np.where(cat['ID'] == id_cat[idx])[0]
        outputDict['ID'].append(id_cat[idx])
        outputDict['Re_circ'].append(cat['Re_circ'][aa[0]])
        outputDict['Re_maj'].append(cat['Re_maj'][aa[0]])
        outputDict['n'].append(cat['sersic_n'][aa[0]])
        outputDict['pa'].append(cat['position_angle'][aa[0]])
        #The script can work with dithers, but for that the centered estimates won't change
        outputDict['offx_0'].append(ditherDict['rx0'][idx])
        outputDict['offy_0'].append(ditherDict['ry0'][idx])
        outputDict['offx_1'].append(ditherDict['rx1'][idx])
        outputDict['offy_1'].append(ditherDict['ry1'][idx])
        outputDict['offx_2'].append(ditherDict['rx2'][idx])
        outputDict['offy_2'].append(ditherDict['ry2'][idx])
        if len(aa) == 1:
            lineFluxArr = []
            lineWaveArr = []
            for line in linesInOrder:
                #for line in ['HBaA_6563']:
                print line
                if args.BEAGLEcatFile:
                    lineFluxArr.append(cat[line+'_flux'][aa][0])
                else:
                    lineFluxArr.append(10**cat[line+'_flux'][aa][0])
                lineWaveArr.append((cat['redshift'][aa][0]+1.)*lineDict[line])
       
            lineDoublet = lineDoublet_prism #make sure something is assigned to lineDoublet for allocating in inputs, but change it based on the desired task in the following lines 
            for i,line in enumerate(linesInOrder):
                writeOutputTest = 0
                if args.R100lines:
                    writeOutputTest = writeOutput_prism[i]
                    lineDoublet = lineDoublet_prism
                if args.R1000lines:
                    writeOutputTest = writeOutput_R1000[i]
                    lineDoublet = lineDoublet_R1000
                if writeOutputTest == 1:
                    if lineDoublet[i] > 0:
                        tempIdx = np.where(lineDoublet == lineDoublet[i])[0]
                        if len(tempIdx) > 2:
                            print "error, only accounting for doublets, shouldn't have any more than two lines blended, ", len(tempIdx)
                            sys.exit()
                        outputLable = linesInOrder[tempIdx[0]]+'_'+linesInOrder[tempIdx[1]]
                        if args.BEAGLEcatFile:
                            outputFlux = cat[linesInOrder[tempIdx[0]]+'_flux'][aa[0]]+cat[linesInOrder[tempIdx[1]]+'_flux'][aa[0]]
                        else:
                            outputFlux = 10**cat[linesInOrder[tempIdx[0]]+'_flux'][aa[0]]+10**cat[linesInOrder[tempIdx[1]]+'_flux'][aa[0]]
                        outputDict[outputLable+'_flux'].append(outputFlux)
                    else:
                        outputLable = line
                        if args.BEAGLEcatFile:
                            outputDict[line+'_flux'].append(cat[line+'_flux'][aa[0]])
                        else:
                            outputDict[line+'_flux'].append(10.**cat[line+'_flux'][aa[0]])
                    print i, line, outputLable, outputDict[outputLable+'_flux'], aa
        

            rCirc = cat['Re_circ'][aa[0]]
            reMaj = cat['Re_maj'][aa[0]]
            xOff = ditherDict['rx0'][idx]
            yOff = ditherDict['ry0'][idx]
            sersic = cat['sersic_n'][aa[0]]
            axisRatio = cat['axis_ratio'][aa[0]]
            pa = cat['position_angle'][aa[0]]
            redshift = cat['redshift'][aa][0]
            
            inputs=[{'xoff':0,'yoff':0,'r_circ':-1,'re_maj':-1,'sersic_n':sersic,'axis_ratio':axisRatio,'position_angle':pa,'flux':lineFluxArr,'wave':lineWaveArr, 'redshift':redshift,     'lineDoublet':lineDoublet},\
                    {'xoff':xOff,'yoff':yOff,'r_circ':-1,'re_maj':-1,'sersic_n':sersic,'axis_ratio':axisRatio,'position_angle':pa,'flux':lineFluxArr,'wave':lineWaveArr, 'redshift':redshift, 'lineDoublet':lineDoublet},\
                    {'xoff':0,'yoff':0,'r_circ':rCirc,'re_maj':reMaj,'sersic_n':sersic,'axis_ratio':axisRatio,'position_angle':pa,'flux':lineFluxArr,'wave':lineWaveArr, 'redshift':redshift, 'lineDoublet':lineDoublet},\
                    {'xoff':xOff,'yoff':yOff,'r_circ':rCirc,'re_maj':reMaj,'sersic_n':sersic,'axis_ratio':axisRatio,'position_angle':pa,'flux':lineFluxArr,'wave':lineWaveArr, 'redshift':redshift, 'lineDoublet':lineDoublet}]#[ps,psOff,es,esOff]
            if args.R100lines:
                for i in range(len(calcArr)):
                    if calcArr[i] == True:
                        temp_sn = defaultdict(list)
                        temp_flux = defaultdict(list)
                        temp_noise = defaultdict(list)
                        if combineArr[i] == False:
                            nExp = args.nExp
                            ditherArr = [[inputs[i]['xoff'],inputs[i]['yoff']]]
                            onSource = [[False,True,False]]
                            slitletShape = [[[0,-1],[0,0],[0,1]]]
                            tempOutputStr = [outputStr[i]]
                        else:
                            nExp = args.nExp/3
                            ditherArr = []
                            tempOutputStr = []
                            onSource = onSourceArr
                            slitletShape = slitletShapeArr
                            for dither in range(3):
                                xKey = 'rx'+str(dither)
                                yKey = 'ry'+str(dither)
                                ditherArr.append([ditherDict[xKey][idx],ditherDict[yKey][idx]])
                                tempOutputStr.append(outputStr[i]+'_'+str(dither))
                        for j,dither in enumerate(ditherArr):
                            tempInputs = inputs[i]
                            tempInputs['xoff'] = dither[0]
                            tempInputs['yoff'] = dither[1]
                            tempInputs['onSource'] = onSource[j]
                            tempInputs['slitletShape'] = slitletShape[j]
                            sn,f,n = sn_array(tempInputs, ngroup=args.nGroup, nint=args.nInt, nexp=nExp)
                            add_to_outputDict(sn,n,f,outputDict, linesInOrder, lineDoublet_prism, writeOutput_prism, '_R100'+tempOutputStr[j], aa, temp_sn = temp_sn, temp_noise = temp_noise, temp_flux = temp_flux)
                        if combineArr[i] == True:
                            #combine info for dithers and save in outputDict
                            for key in temp_sn.keys():
                                tempIdx = np.where(np.array(temp_flux[key]) > 0)[0]
                                if len(tempIdx > 0):
                                    outputDict[key+'_f_R100'+outputStr[i]+'_comb'].append(np.mean(np.array(temp_flux[key])[tempIdx]))
                                    outputDict[key+'_n_R100'+outputStr[i]+'_comb'].append(np.sqrt(np.sum(np.power(np.array(temp_noise[key])[tempIdx],2)))/len(tempIdx))
                                    outputDict[key+'_sn_R100'+outputStr[i]+'_comb'].append(np.mean(np.array(temp_flux[key])[tempIdx])/(np.sqrt(np.sum(np.power(np.array(temp_noise[key])[tempIdx],2)))/len(tempIdx)))
                                else:
                                    outputDict[key+'_f_R100'+outputStr[i]+'_comb'].append(0.)
                                    outputDict[key+'_sn_R100'+outputStr[i]+'_comb'].append(0.)
                                    outputDict[key+'_n_R100'+outputStr[i]+'_comb'].append(-99)

            if args.R1000lines:
                for i in range(len(calcArr)):
                    if calcArr[i] == True:
                        temp_sn = defaultdict(list)
                        temp_flux = defaultdict(list)
                        temp_noise = defaultdict(list)
                        if combineArr[i] == False:
                            nExp = args.nExp
                            ditherArr = [[inputs[i]['xoff'],inputs[i]['yoff']]]
                            tempOutputStr = [outputStr[i]]
                            onSource = [[False,True,False]]
                            slitletShape = [[[0,-1],[0,0],[0,1]]]
          
                        else:
                            nExp = args.nExp/3
                            ditherArr = []
                            tempOutputStr = []
                            onSource = onSourceArr
                            slitletShape = slitletShapeArr
                            for dither in range(3):
                                xKey = 'rx'+str(dither)
                                yKey = 'ry'+str(dither)
                                ditherArr.append([ditherDict[xKey][idx],ditherDict[yKey][idx]])
                                tempOutputStr.append(outputStr[i]+'_'+str(dither))
                        for configIdx in range(nConfig):
                            for j,dither in enumerate(ditherArr):
                                tempInputs = inputs[i]
                                tempInputs['xoff'] = dither[0]
                                tempInputs['yoff'] = dither[1]
                                tempInputs['lineDoublet'] = lineDoublet_R1000
                                tempInputs['onSource'] = onSource[j]
                                tempInputs['slitletShape'] = slitletShape[j]
                                sn,f,n = sn_array(tempInputs, disperser=grating[configIdx], filt=filter[configIdx], ngroup = args.nGroup, nint = args.nInt, nexp = nExp)
                                add_to_outputDict(sn,n,f,outputDict, linesInOrder, lineDoublet_R1000, writeOutput_R1000, '_R1000_'+grating[configIdx]+'_'+filter[configIdx]+tempOutputStr[j], aa, temp_sn = temp_sn, temp_noise = temp_noise, temp_flux = temp_flux)
                            if combineArr[i] == True:
                                #combine info for dithers and save in outputDict
                                for key in temp_sn.keys():
                                    tempIdx = np.where(np.array(temp_flux[key]) > 0)[0]
                                    if len(tempIdx > 0):
                                        outputDict[key+'_f_R1000_'+grating[configIdx]+'_'+filter[configIdx]+outputStr[i]+'_comb'].append(np.mean(np.array(temp_flux[key])[tempIdx]))
                                        outputDict[key+'_n_R1000_'+grating[configIdx]+'_'+filter[configIdx]+outputStr[i]+'_comb'].append(np.sqrt(np.sum(np.power(np.array(temp_noise[key])[tempIdx],2)))/len(tempIdx))
                                        outputDict[key+'_sn_R1000_'+grating[configIdx]+'_'+filter[configIdx]+outputStr[i]+'_comb'].append(np.mean(np.array(temp_flux[key])[tempIdx])/(np.sqrt(np.sum(np.power(np.array(temp_noise[key])[tempIdx],2)))/len(tempIdx)))
                                    else:
                                        outputDict[key+'_f_R1000_'+grating[configIdx]+'_'+filter[configIdx]+outputStr[i]+'_comb'].append(0.)
                                        outputDict[key+'_sn_R1000_'+grating[configIdx]+'_'+filter[configIdx]+outputStr[i]+'_comb'].append(0.)
                                        outputDict[key+'_n_R1000_'+grating[configIdx]+'_'+filter[configIdx]+outputStr[i]+'_comb'].append(-99)

            if args.R100continuum:
                for i in range(len(calcArr)):
                    reportArr = []
                    if calcArr[i] == True:
                        if combineArr[i] == False:
                            nExp = args.nExp
                            ditherArr = [[inputs[i]['xoff'],inputs[i]['yoff']]]
                            tempOutputStr = [outputStr[i]]
                            onSource = [[False,True,False]]
                            slitletShape = [[[0,-1],[0,0],[0,1]]]
                        else:
                            nExp = args.nExp/3
                            ditherArr = []
                            tempOutputStr = []
                            onSource = onSourceArr
                            slitletShape = slitletShapeArr
                            for dither in range(3):
                                xKey = 'rx'+str(dither)
                                yKey = 'ry'+str(dither)
                                ditherArr.append([ditherDict[xKey][idx],ditherDict[yKey][idx]])
                                tempOutputStr.append(outputStr[i]+'_'+str(dither))
                        tempInputs = inputs[i]
                        tempInputs['wl'] = inputDict['wl']
                        tempInputs['spec'] = inputDict['spec']
                        for j,dither in enumerate(ditherArr):
                            tempInputs['xoff'] = dither[0]
                            tempInputs['yoff'] = dither[1]
                            tempInputs['onSource'] = onSource[j]
                            tempInputs['slitletShape'] = slitletShape[j]
                            report = sn_user_spec(tempInputs, ngroup=args.nGroup, nint=args.nInt, nexp=nExp)
                            tempOutputFile = args.resultsDir+str(id_cat[idx])+'_R100.p'
                            pickle.dump(report,open(tempOutputFile,'w'))
                            add_cont_sn_to_outputDict(report,outputDict, '_R100'+tempOutputStr[j], inputDict['regions'], tempInputs['redshift'])
                            reportArr.append(report)
                        if combineArr[i] == True:
                            tempReport = reportArr[0]
                            tempReport['1d']['extracted_noise'][1] = np.power(tempReport['1d']['extracted_noise'][1],2)
                            for j in range(1,3):
                                tempReport['1d']['extracted_flux'][1] = tempReport['1d']['extracted_flux'][1]+reportArr[j]['1d']['extracted_flux'][1]
                                tempReport['1d']['extracted_noise'][1] = tempReport['1d']['extracted_noise'][1] + np.power(reportArr[j]['1d']['extracted_noise'][1],2)
                            tempReport['1d']['extracted_noise'][1] = np.sqrt(tempReport['1d']['extracted_noise'][1])
                            add_cont_sn_to_outputDict(tempReport,outputDict, '_R100'+outputStr[i]+'_comb', inputDict['regions'], tempInputs['redshift'])


            if args.R1000continuum:
                for i in range(len(calcArr)):
                    if calcArr[i] == True:
                        if combineArr[i] == False:
                            nExp = args.nExp
                            ditherArr = [[inputs[i]['xoff'],inputs[i]['yoff']]]
                            tempOutputStr = [outputStr[i]]
                            onSource = [[False,True,False]]
                            slitletShape = [[[0,-1],[0,0],[0,1]]]
                        else:
                            nExp = args.nExp/3
                            ditherArr = []
                            tempOutputStr = []
                            onSource = onSourceArr
                            slitletShape = slitletShapeArr
                            for dither in range(3):
                                xKey = 'rx'+str(dither)
                                yKey = 'ry'+str(dither)
                                ditherArr.append([ditherDict[xKey][idx],ditherDict[yKey][idx]])
                                tempOutputStr.append(outputStr[i]+'_'+str(dither))
                        tempInputs = inputs[i]
                        tempInputs['wl'] = inputDict['wl']
                        tempInputs['spec'] = inputDict['spec']
                        for configIdx in range(nConfig):
                            reportArr = []
                            for j,dither in enumerate(ditherArr):
                                tempInputs['xoff'] = dither[0]
                                tempInputs['yoff'] = dither[1]
                                tempInputs['onSource'] = onSource[j]
                                tempInputs['slitletShape'] = slitletShape[j]
                                report = sn_user_spec(tempInputs, disperser=grating[configIdx], filt=filter[configIdx], ngroup=args.nGroup, nint=args.nInt, nexp=nExp)
                                tempOutputFile = args.resultsDir+str(id_cat[idx])+'_R1000_'+grating[configIdx]+'_'+filter[configIdx]+".p"
                                pickle.dump(report,open(tempOutputFile,'w'))
                                add_cont_sn_to_outputDict(report,outputDict, '_R1000_'+grating[configIdx]+'_'+filter[configIdx]+tempOutputStr[j], inputDict['regions'], tempInputs['redshift'])
                                reportArr.append(report)
                            if combineArr[i] == True:
                                tempReport = reportArr[0]
                                tempReport['1d']['extracted_noise'][1] = np.power(tempReport['1d']['extracted_noise'][1],2)
                                for j in range(1,3):
                                    tempReport['1d']['extracted_flux'][1] = tempReport['1d']['extracted_flux'][1]+reportArr[j]['1d']['extracted_flux'][1]
                                    tempReport['1d']['extracted_noise'][1] = tempReport['1d']['extracted_noise'][1] + np.power(reportArr[j]['1d']['extracted_noise'][1],2)
                                tempReport['1d']['extracted_noise'][1] = np.sqrt(tempReport['1d']['extracted_noise'][1])
                                add_cont_sn_to_outputDict(tempReport,outputDict, '_R1000_'+grating[configIdx]+'_'+filter[configIdx]+outputStr[i]+'_comb', inputDict['regions'], tempInputs['redshift'])


        else:
            for i,line in enumerate(linesInOrder):
                if args.R100:
                    if writeOutput_prism[i] == 1:
                        outputLable = line
                    if lineDoublet_prism[i] > 0:
                        tempIdx = np.where(lineDoublet_prism == lineDoublet_prism[i])[0]
                        if len(tempIdx) > 2:
                            print "error, only accounting for doublets, shouldn't have any more than two lines blended, ", len(tempIdx)
                            sys.exit()
                        outputLable = linesInOrder[tempIdx[0]]+'_'+linesInOrder[tempIdx[1]]
                if args.R1000:
                    if writeOutput_R1000[i] == 1:
                        outputLable = line
                    if lineDoublet_R1000[1] > 0:
                        tempIdx = np.where(lineDoublet_R1000 == lineDoublet_R1000[i])[0]
                        if len(tempIdx) > 2:
                            print "error, only accounting for doublets, shouldn't have any more than two lines blended, ", len(tempIdx)
                            sys.exit()
                        outputLable = linesInOrder[tempIdx[0]]+'_'+linesInOrder[tempIdx[1]]

                outputConfigurations = []
                if args.R100lines:
                    outputConfigurations.append(['R100'])
                if args.R1000lines:
                    for configIdx in range(nConfig ):
                        outputConfigurations.append(['R1000_'+grating[configIdx]+'_'+filter[configIdx]])
                outputDict[outputLable+'_flux'].append(0.)
                for config in outputConfigurations:
                    if args.ps:
                        outputDict[outputLable+'_sn_'+config+'_ps'].append(0.)
                        outputDict[outputLable+'_f_'+config+'_ps'].append(0.)
                        outputDict[outputLable+'_n_'+config+'_ps'].append(-99)
                    if args.es:
                        outputDict[outputLable+'_sn_'+config+'_es'].append(0.)
                        outputDict[outputLable+'_f_'+config+'_es'].append(0.)
                        outputDict[outputLable+'_n_'+config+'_es'].append(0.)
                    if args.singleDither:
                        if args.psOff:
                            outputDict[outputLable+'_sn_'+config+'_ps_off'].append(0.)
                            outputDict[outputLable+'_f_'+config+'_ps_off'].append(0.)
                            outputDict[outputLable+'_n_'+config+'_ps_off'].append(-99)
                        if args.esOff:
                            outputDict[outputLable+'_sn_'+config+'_es_off'].append(0.)
                            outputDict[outputLable+'_f_'+config+'_es_off'].append(0.)
                            outputDict[outputLable+'_n_'+config+'_es_off'].append(-99)
                    else:
                        for i in range(3):
                            if args.psOff:
                                outputDict[outputLable+'_sn_'+config+'_ps_off_'+str(i)].append(0.)
                                outputDict[outputLable+'_f_'+config+'_ps_off_'+str(i)].append(0.)
                                outputDict[outputLable+'_n_'+config+'_ps_off_'+str(i)].append(-99)
                            if args.esOff:
                                outputDict[outputLable+'_sn_'+config+'_es_off_'+str(i)].append(0.)
                                outputDict[outputLable+'_f_'+config+'_es_off_'+str(i)].append(0.)
                                outputDict[outputLable+'_n_'+config+'_es_off_'+str(i)].append(-99)
                        if args.psOff:
                            outputDict[outputLable+'_sn_'+config+'_ps_off_comb'].append(0.)
                            outputDict[outputLable+'_f_'+config+'_ps_off_comb'].append(0.)
                            outputDict[outputLable+'_n_'+config+'_ps_off_comb'].append(-99)
                        if args.esOff:
                            outputDict[outputLable+'_sn_'+config+'_es_off_comb'].append(0.)
                            outputDict[outputLable+'_f_'+config+'_es_off_comb'].append(0.)
                            outputDict[outputLable+'_n_'+config+'_es_off_comb'].append(-99)
            for region in inputDict['regions'].keys():
                outputConfigurations = []
                if args.R100continuum:
                    outputConfigurations.append(['R100'])
                if args.R1000continuum:
                    outputConfigurations.append(['R1000_'+grating[configIdx]+'_'+filter[configIdx]])
                for config in outputConfigurations:
                    if args.ps:
                        outputDict[region+'_avg_sn_pp'+config+'_ps']
                    if args.es:
                        outputDict[region+'_avg_sn_pp'+config+'_es']
                    if args.singleDither:
                        if args.psOff:
                            outputDict[region+'_avg_sn_pp'+config+'_ps_off']
                        if args.esOff:
                            outputDict[region+'_avg_sn_pp'+config+'_es_off']
                    else:
                        if args.psOff:
                            for i in range(3):
                                outputDict[region+'_avg_sn_pp'+config+'_ps_off_'+str(i)]
                            outputDict[region+'_avg_sn_pp'+config+'_ps_off_comb']
                        if args.esOff:
                            for i in range(3):
                                outputDict[region+'_avg_sn_pp'+config+'_es_off'+str(i)]
                            outputDict[region+'_avg_sn_pp'+config+'_es_off_comb']
        time_end = time.time()
        print 'time elapsed (min): ', (time_end-time_start)/60.

        print outputDict.keys()
        print outputDict
        for key in outputDict.keys():
            if 'flux' in key:
                print key, outputDict[key]
        pickle.dump(outputDict,open(outputFile,'w'))




parser = argparse.ArgumentParser()
parser.add_argument(
                    '--eMPT-out',
                    help="eMPT tile output file containing offset information",
                    action="store",
                    type=str,
                    dest="eMPTout",
                    required=True
                    )

#parser.add_argument(
#                    '--eMPT-out-list',
#                    help="if a list of eMPT tile output files are required to be read in, then --eMPT-out contains the filename of this list, there will be one output file per list and total S/N can be determine after from extracted flux and extracted noise columns",
#                    action="store_true",
#                    default=false,
#                    dest="eMPToutList"
#                    )

parser.add_argument(
                    '--single-dither',
                    help="if there is only a single dither included in the eMPT output file",
                    action="store_true",
                    default=False,
                    dest="singleDither"
                    )

parser.add_argument(
                    '-c','--cat-file',
                    help="JAGUAR catalogue file containing input parameters",
                    action="store",
                    type=str,
                    dest="catFile",
                    required=True
                    )

parser.add_argument(
                    '-B','--BEAGLE-cat-file',
                    help="The catalogue given is in the format of a BEAGLE output file",
                    action="store_true",
                    default=False,
                    dest="BEAGLEcatFile"
                    )

parser.add_argument(
                    '-s','--size-file',
                    help="If catalogue given is in the format of a BEAGLE output file, there is a separate file containing the size info",
                    action="store",
                    type=str,
                    dest="sizeFile"
                    )

parser.add_argument(
                    '-o','--output-file',
                    help="output fits file name",
                    action="store",
                    type=str,
                    dest="outputFile",
                    required=True
                    )

parser.add_argument(
                    '-r','--results-dir',
                    help="results for individual objects will be saved as pickle files in this directory",
                    action="store",
                    type=str,
                    dest="resultsDir",
                    required=True
                    )

parser.add_argument(
                    '--n-group',
                    help="number of groups",
                    action="store",
                    type=int,
                    default=19,
                    required=False,
                    dest="nGroup"
                    )

parser.add_argument(
                    '--n-int',
                    help="number of integrations, or resets within a single exposure",
                    action="store",
                    type=int,
                    default=2,
                    required=False,
                    dest="nInt"
                    )

parser.add_argument(
                    '--n-exp',
                    help="number of exposures",
                    action="store",
                    type=int,
                    default=36,
                    required=False,
                    dest="nExp"
                    )

parser.add_argument(
                    '--ps',
                    help="produce S/N for a centred point source?",
                    action="store_true",
                    default=False,
                    required=False,
                    dest="ps"
                    )

parser.add_argument(
                    '--ps-off',
                    help="produce S/N for point source at actual position in slit?",
                    action="store_true",
                    default=False,
                    required=False,
                    dest="psOff"
                    )

parser.add_argument(
                    '--es',
                    help="produce S/N for a centred extended source?",
                    action="store_true",
                    default=False,
                    required=False,
                    dest="es"
                    )

parser.add_argument(
                    '--es-off',
                    help="produce S/N for an extended source at actual position in slit?",
                    action="store_true",
                    default=False,
                    required=False,
                    dest="esOff"
                    )

parser.add_argument(
                    '--R100-lines',
                    help="produce S/N for emission lines in R100 spectra?",
                    action="store_true",
                    default=False,
                    required=False,
                    dest="R100lines"
                    )

parser.add_argument(
                    '--R1000-lines',
                    help="produce S/N for emission lines in R100 spectra?",
                    action="store_true",
                    default=False,
                    required=False,
                    dest="R1000lines"
                    )

parser.add_argument(
                    '--R100-continuum',
                    help="produce S/N for continuum regions in R100 spectra?",
                    action="store_true",
                    default=False,
                    required=False,
                    dest="R100continuum"
                    )

parser.add_argument(
                    '--R1000-continuum',
                    help="produce S/N for continuum regions in R1000 spectra?",
                    action="store_true",
                    default=False,
                    required=False,
                    dest="R1000continuum"
                    )


parser.add_argument(
                    '--make-output-fits',
                    help="do you just want to make the output fits file?  It will be made automatically at the end of the calculation, but you can also make it separatelly with this flag.  Make sure the es/ps etc. flags correspond to the data stored in the pickle files.",
                    action="store_true",
                    default=False,
                    required=False,
                    dest="makeOutputFits"
                    )

parser.add_argument(
                    '--n-proc',
                    help="number of processes to set running in the pool",
                    action="store",
                    type=int,
                    default=10,
                    required=False,
                    dest="nProc"
                    )

args = parser.parse_args()

#read in the offset information
a=np.genfromtxt(args.eMPTout, dtype=None, names=True)
id_cat=a['ID_cat']
clas=a['cla']
ditherDict={'rx0':a['rx0']*267.9e-3,\
            'ry0':a['ry0']*528.9e-3,\
            'rx1':a['rx1']*267.9e-3,\
            'ry1':a['ry1']*528.9e-3,\
            'rx2':a['rx2']*267.9e-3,\
            'ry2':a['ry2']*528.9e-3}


lineDict = {}
lineDict['C4_1548']=0.1548
lineDict['C4_1551']=0.1551
lineDict['O3_1661']=0.1661
lineDict['O3_1666']=0.1666
lineDict['C3_1907']=0.1907
lineDict['O2_3726']=0.3726
lineDict['O2_3729']=0.3729
lineDict['Ne3_3869']=0.3869
lineDict['O3_4363']=0.4363
lineDict['HBaB_4861']=0.4861
lineDict['O3_4959']=0.4959
lineDict['O3_5007']=0.5007
lineDict['HBaA_6563']=0.6563
lineDict['N2_6584']=0.6584
lineDict['S2_6716']=0.6716
lineDict['S2_6731']=0.6731

continuumDict = {}
#start off with some regions around lines - chunk lines together if their doublets
continuumDict['1550'] = {'low':0.1540, 'high':0.1560}
continuumDict['1660'] = {'low':0.1650, 'high':0.1670}
continuumDict['1910'] = {'low':0.1900, 'high':0.1920}
continuumDict['3730'] = {'low':0.3720, 'high':0.3740}
continuumDict['3870'] = {'low':0.3860, 'high':0.3880}
continuumDict['4360'] = {'low':0.4350, 'high':0.4370}
continuumDict['4860'] = {'low':0.4850, 'high':0.4870}
continuumDict['4960'] = {'low':0.4950, 'high':0.4970}
continuumDict['5010'] = {'low':0.5000, 'high':0.5020}
continuumDict['6560'] = {'low':0.6550, 'high':0.6570}
continuumDict['6720'] = {'low':0.6710, 'high':0.6730}
continuumDict['D4000_1'] = {'low':0.3750, 'high':0.3950} #following Poggianti 1997 defn of D4000
continuumDict['D4000_2'] = {'low':0.4050, 'high':0.4250}
continuumDict['Calzetti_1'] = {'low':0.1268, 'high':0.1284}
continuumDict['Calzetti_2'] = {'low':0.1309, 'high':0.1316}
continuumDict['Calzetti_3'] = {'low':0.1342, 'high':0.1371}
continuumDict['Calzetti_4'] = {'low':0.1407, 'high':0.1515}
continuumDict['Calzetti_5'] = {'low':0.1562, 'high':0.1583}
continuumDict['Calzetti_6'] = {'low':0.1677, 'high':0.1740}
continuumDict['Calzetti_7'] = {'low':0.1760, 'high':0.1833}
continuumDict['Calzetti_8'] = {'low':0.1866, 'high':0.1890}
continuumDict['Calzetti_9'] = {'low':0.1930, 'high':0.1950}
continuumDict['Calzetti_10'] = {'low':0.2400, 'high':0.2580}

#I'm not adding the following for now - I think it would be overkill
#index_CN1 4081 4285
#index_Ca4227 4212 4252
#index_G4300 4267 4336
#index_Fe4383 4360 4456
#index_Ca4455 4447 4493
#index_Fe4531 4505 4580
#index_Fe4668 4612 4757
#index_Fe5015 4946 5065
#index_Mg1 4897 5366
#index_Mg2 4897 5366
#index_Fe5270 5235 5319
#index_Fe5335 5307 5364
#index_Fe5406 5376 5425
#index_Fe5709 5674 5738
#index_Fe5782 5767 5813
#index_NaD 5863 5949
#index_TiO1 5819 6104
#index_TiO2 6069 6416


#load in the JAGUAR catalog
a=fits.open(args.catFile)
if 'C3_1909_flux' not in a['HII EMISSION'].data.dtype.names:
    lineDict['C3_1910']=0.1909
else:
    lineDict['C3_1909']=0.1909

if args.BEAGLEcatFile:
    s = fits.open(args.sizeFile)
    cat = {}
    cat['ID'] = a['IDs'].data['ID_cat']
    cat['Re_circ'] = s[1].data['Re_circ_arcsec']
    cat['Re_maj'] = s[1].data['Re_maj']
    if 'sersic_n' in s[1].data.dtype.names:
        cat['sersic_n'] = s[1].data['sersic_n']
        cat['position_angle'] = s[1].data['position_angle']
    else:
        if args.es or args.esOff:
            print s[1].data.dtype.names
            print 'ERROR: you sizes file does not contain sersic or position angle information!'
            sys.exit()
        cat['sersic_n'] = np.zeros_like(s[1].data['Re_maj'])+1.
        cat['position_angle'] = np.zeros_like(s[1].data['Re_maj'])
    cat['axis_ratio'] = s[1].data['axis_ratio']
    cat['redshift'] = a['GALAXY PROPERTIES'].data['redshift']
    for line in lineDict:
        cat[line+'_flux'] = a['HII EMISSION'].data[line+'_flux']
        print line, cat[line+'_flux'][:2]
    if args.R100continuum or args.R1000continuum:
        wl = a['CONTINUUM SED WL'].data[0][0]
        specArr = a['CONTINUUM SED'].data
    s.close()
else:
    cat={}
    for col in a[1].data.dtype.names():
        cat[col]=a[1].data[col]

if args.R1000continuum or args.R100continuum:
    if not args.BEAGLEcatFile:
        print 'error, need a BEAGLE input file containing continuum spectra to be able to calculate continuum S/N values'
        sys.exit()

a.close()
#pick some random ones
snarr_ps=[]
snarr_ps_off=[]
snarr_es=[]
snarr_es_off=[]
fluxarr=[]
sizearr=[]



lineDoublet_prism = np.array([1,1,2,2,3,3,4,4,0,0,0,5,5,0,0,6,6])
#write doublets as a single entry in output file - the next array says when to write the output so as not to duplicate information
writeOutput_prism = [1,0,1,0,1,0,1,0,1,1,1,1,0,1,1,1,0]

lineDoublet_R1000 = np.array([0,0,0,0,1,1,2,2,0,0,0,0,0,0,0,0,0])
writeOutput_R1000 = [1,1,1,1,1,0,1,0,1,1,1,1,1,1,1,1,1]



linesInOrder = ['C4_1548','C4_1551','O3_1661','O3_1666','C3_1907','C3_1909','O2_3726','O2_3729',\
                    'Ne3_3869','O3_4363','HBaB_4861','O3_4959','O3_5007','HBaA_6563','N2_6584',\
                    'S2_6716','S2_6731']


#Check through catalogue - sometimes C3_1909 is labelled C3_1910
if 'C3_1909_flux' not in cat.keys():
    linesInOrder[5] = 'C3_1910'

c_light = 2.99792e+18 # Ang/s
if not args.makeOutputFits:
    print 'here'
    inputArray = []
    for i in range(0,len(id_cat)):
        inputs = {'idx':i, 'id_cat':id_cat, 'cat':cat, 'ditherDict':ditherDict, 'linesInOrder':linesInOrder, 'writeOutput_prism':writeOutput_prism, 'lineDoublet_prism':lineDoublet_prism, 'args':args, \
            'lineDoublet_R1000':lineDoublet_R1000, 'writeOutput_R1000':writeOutput_R1000}
        if args.R100continuum or args.R1000continuum:
            aa=np.where(cat['ID'] == id_cat[1])[0]
            tempWl = wl*(1+cat['redshift'][aa[0]])
            tempSpec = specArr[aa,:]/(1+cat['redshift'][aa[0]])
            tempSpec = tempSpec[0]
            tempSpec = tempWl**2/c_light*tempSpec * 1E23 * 1.E3 #in mJy
            tempWl = tempWl/1.E4 #in microns
            deduped = {w: s for w, s in reversed(zip(list(tempWl), list(tempSpec)))}
            tempWl = np.asarray(sorted(deduped.keys()))
            tempSpec = np.asarray([deduped[k] for k in tempWl])
            inputs['wl'] = tempWl
            inputs['spec'] = tempSpec 
            inputs['regions'] = continuumDict
        inputArray.append(inputs)
    

    #pool = ProcessingPool(24)
    pool = ProcessingPool(nodes=args.nProc)
    pool.map(obj_sn, inputArray)
    #for i in range(len(inputArray)):
    #  obj_sn(inputArray[i])

outputColumnOrder=['ID','Re_circ','Re_maj','n','pa','offx_0','offy_0','offx_1','offy_1','offx_2','offy_2']

outputLinesInOrder_R100 = ['C4_1548_C4_1551','O3_1661_O3_1666','C3_1907_C3_1909','O2_3726_O2_3729',\
                           'Ne3_3869','O3_4363','HBaB_4861','O3_4959_O3_5007','HBaA_6563','N2_6584',\
                           'S2_6716_S2_6731']

outputLinesInOrder_R1000 = ['C4_1548','C4_1551','O3_1661','O3_1666','C3_1907_C3_1909','O2_3726_O2_3729',\
                      'Ne3_3869','O3_4363','HBaB_4861','O3_4959','O3_5007','HBaA_6563','N2_6584',\
                      'S2_6716','S2_6731']

outputContinuumInOrder= ['1550','1660','1910','3730','3870','4360','4860','4960','5010','6560','6720','D4000_1','D4000_2',\
                         'Calzetti_1','Calzetti_2','Calzetti_3','Calzetti_4','Calzetti_5','Calzetti_6','Calzetti_7',\
                         'Calzetti_8','Calzetti_9','Calzetti_10']

if 'C3_1909_flux' not in cat.keys():
    outputLinesInOrder_R100[2] = 'C3_1907_C3_1910'
    outputLinesInOrder_R1000[4] = 'C3_1907_C3_1910'

#this part needs to be made more flexible to be able to handle different resolution modes!
if args.R100lines:
    for line in outputLinesInOrder_R100:
        outputColumnOrder.append(line+'_flux')
        if args.ps:
            outputColumnOrder.append(line+'_sn_R100_ps')
            outputColumnOrder.append(line+'_f_R100_ps')
            outputColumnOrder.append(line+'_n_R100_ps')
        if args.psOff:
            if args.singleDither:
                outputColumnOrder.append(line+'_sn_R100_ps_off')
                outputColumnOrder.append(line+'_f_R100_ps_off')
                outputColumnOrder.append(line+'_n_R100_ps_off')
            else:
                for i in range(3):
                    outputColumnOrder.append(line+'_sn_R100_ps_off_'+str(i))
                    outputColumnOrder.append(line+'_f_R100_ps_off_'+str(i))
                    outputColumnOrder.append(line+'_n_R100_ps_off_'+str(i))
                outputColumnOrder.append(line+'_sn_R100_ps_off_comb')
                outputColumnOrder.append(line+'_f_R100_ps_off_comb')
                outputColumnOrder.append(line+'_n_R100_ps_off_comb')
        if args.es:
            outputColumnOrder.append(line+'_sn_R100_es')
            outputColumnOrder.append(line+'_f_R100_es')
            outputColumnOrder.append(line+'_n_R100_es')
        if args.esOff:
            if args.singleDither:
                outputColumnOrder.append(line+'_sn_R100_es_off')
                outputColumnOrder.append(line+'_f_R100_es_off')
                outputColumnOrder.append(line+'_n_R100_es_off')
            else:
                for i in range(3):
                    outputColumnOrder.append(line+'_sn_R100_es_off_'+str(i))
                    outputColumnOrder.append(line+'_f_R100_es_off_'+str(i))
                    outputColumnOrder.append(line+'_n_R100_es_off_'+str(i))
                outputColumnOrder.append(line+'_sn_R100_es_off_comb')
                outputColumnOrder.append(line+'_f_R100_es_off_comb')
                outputColumnOrder.append(line+'_n_R100_es_off_comb')
if args.R1000lines:
    for line in outputLinesInOrder_R1000:
        outputColumnOrder.append(line+'_flux')
        for configIdx in range(nConfig):
            if args.ps:
                outputColumnOrder.append(line+'_sn_R1000_'+grating[configIdx]+'_'+filter[configIdx]+'_ps')
                outputColumnOrder.append(line+'_f_R1000_'+grating[configIdx]+'_'+filter[configIdx]+'_ps')
                outputColumnOrder.append(line+'_n_R1000_'+grating[configIdx]+'_'+filter[configIdx]+'_ps')
            if args.psOff:
                if args.singleDither:
                    outputColumnOrder.append(line+'_sn_R1000_'+grating[configIdx]+'_'+filter[configIdx]+'_ps_off')
                    outputColumnOrder.append(line+'_f_R1000_'+grating[configIdx]+'_'+filter[configIdx]+'_ps_off')
                    outputColumnOrder.append(line+'_n_R1000_'+grating[configIdx]+'_'+filter[configIdx]+'_ps_off')
                else:
                    for i in range(3):
                        outputColumnOrder.append(line+'_sn_R1000_'+grating[configIdx]+'_'+filter[configIdx]+'_ps_off_'+str(i))
                        outputColumnOrder.append(line+'_f_R1000_'+grating[configIdx]+'_'+filter[configIdx]+'_ps_off_'+str(i))
                        outputColumnOrder.append(line+'_n_R1000_'+grating[configIdx]+'_'+filter[configIdx]+'_ps_off_'+str(i))
                    outputColumnOrder.append(line+'_sn_R1000_'+grating[configIdx]+'_'+filter[configIdx]+'_ps_off_comb')
                    outputColumnOrder.append(line+'_f_R1000_'+grating[configIdx]+'_'+filter[configIdx]+'_ps_off_comb')
                    outputColumnOrder.append(line+'_n_R1000_'+grating[configIdx]+'_'+filter[configIdx]+'_ps_off_comb')
            if args.es:
                outputColumnOrder.append(line+'_sn_R1000_'+grating[configIdx]+'_'+filter[configIdx]+'_es')
                outputColumnOrder.append(line+'_f_R1000_'+grating[configIdx]+'_'+filter[configIdx]+'_es')
                outputColumnOrder.append(line+'_n_R1000_'+grating[configIdx]+'_'+filter[configIdx]+'_es')
            if args.esOff:
                if args.singleDither:
                    outputColumnOrder.append(line+'_sn_R1000_'+grating[configIdx]+'_'+filter[configIdx]+'_es_off')
                    outputColumnOrder.append(line+'_f_R1000_'+grating[configIdx]+'_'+filter[configIdx]+'_es_off')
                    outputColumnOrder.append(line+'_n_R1000_'+grating[configIdx]+'_'+filter[configIdx]+'_es_off')
                else:
                    for i in range(3):
                        outputColumnOrder.append(line+'_sn_R1000_'+grating[configIdx]+'_'+filter[configIdx]+'_es_off_'+str(i))
                        outputColumnOrder.append(line+'_f_R1000_'+grating[configIdx]+'_'+filter[configIdx]+'_es_off_'+str(i))
                        outputColumnOrder.append(line+'_n_R1000_'+grating[configIdx]+'_'+filter[configIdx]+'_es_off_'+str(i))
                    outputColumnOrder.append(line+'_sn_R1000_'+grating[configIdx]+'_'+filter[configIdx]+'_es_off_comb')
                    outputColumnOrder.append(line+'_f_R1000_'+grating[configIdx]+'_'+filter[configIdx]+'_es_off_comb')
                    outputColumnOrder.append(line+'_n_R1000_'+grating[configIdx]+'_'+filter[configIdx]+'_es_off_comb')
if args.R100continuum:
    for region in outputContinuumInOrder:
        if args.ps:
            outputColumnOrder.append(region+'_avg_sn_pp_R100_ps')
        if args.psOff:
            if args.singleDither:
                outputColumnOrder.append(region+'_avg_sn_pp_R100_ps_off')
            else:
                for i in range(3):
                    outputColumnOrder.append(region+'_avg_sn_pp_R100_ps_off_'+str(i))
                outputColumnOrder.append(region+'_avg_sn_pp_R100_ps_off_comb')
        if args.es:
            outputColumnOrder.append(region+'_avg_sn_pp_R100_es')
        if args.esOff:
            if args.singleDither:
                outputColumnOrder.append(region+'_avg_sn_pp_R100_es_off')
            else:
                for i in range(3):
                    outputColumnOrder.append(region+'_avg_sn_pp_R100_es_off_'+str(i))
                outputColumnOrder.append(region+'_avg_sn_pp_R100_es_off_comb')
if args.R1000continuum:
    for region in outputContinuumInOrder:
        for configIdx in range(nConfig):
            if args.ps:
                outputColumnOrder.append(region+'_avg_sn_pp_R1000_'+grating[configIdx]+'_'+filter[configIdx]+'_ps')
            if args.psOff:
                if args.singleDither:
                    outputColumnOrder.append(region+'_avg_sn_pp_R1000_'+grating[configIdx]+'_'+filter[configIdx]+'_ps_off')
                else:
                    for i in range(3):
                        outputColumnOrder.append(region+'_avg_sn_pp_R1000_'+grating[configIdx]+'_'+filter[configIdx]+'_ps_off_'+str(i))
                    outputColumnOrder.append(region+'_avg_sn_pp_R1000_'+grating[configIdx]+'_'+filter[configIdx]+'_ps_off_comb')
            if args.es:
                outputColumnOrder.append(region+'_avg_sn_pp_R1000_'+grating[configIdx]+'_'+filter[configIdx]+'_es')
            if args.esOff:
                if args.singleDither:
                    outputColumnOrder.append(region+'_avg_sn_pp_R1000_'+grating[configIdx]+'_'+filter[configIdx]+'_es_off')
                else:
                    for i in range(3):
                        outputColumnOrder.append(region+'_avg_sn_pp_R1000_'+grating[configIdx]+'_'+filter[configIdx]+'_es_off_'+str(i))
                    outputColumnOrder.append(region+'_avg_sn_pp_R1000_'+grating[configIdx]+'_'+filter[configIdx]+'_es_off_comb')


outputDictOrdered = OrderedDict()
outputDict = defaultdict(list)

for id in id_cat:
    fileTest = os.path.isfile(args.resultsDir+str(id)+'.p')
    if fileTest is True:
        print args.resultsDir+str(id)+'.p'
        data = pickle.load(open(args.resultsDir+str(id)+'.p','r'))
#        print id, data['ID'], data['O3_4363_flux']
        tempIdx = np.where(data['ID'] == id)[0]
#        print tempIdx
        for col in outputColumnOrder:
            col2 = col
#            if 'comb' not in col:
#                if 'es' in col:
#                    col2 = col.replace('_es_off','__es_off')
#                if 'ps' in col:
#                    col2 = col.replace('_ps_off','__ps_off')
#            print col, col2, data[col2]
#            print isinstance(data[col2][0],list), isinstance(data[col2][0],np.ndarray)
#            print col2, data[col2]
            print col, data[col], hasattr(data[col], '__len__')
            #print data.keys()
            if hasattr(data[col],'__len__'):
                outputDict[col].append(data[col][0])
            else:
                outputDict[col].append(data[col])
#            print outputDict[col]
#        sys.exit()
#        anything = raw_input("prompt")


print '***'
for col in outputColumnOrder:
    outputDictOrdered[col]=np.array(outputDict[col])
    print col, outputDictOrdered[col].shape
tmp = args.outputFile
outputFile = tmp.replace(".fits",".p")
pickle.dump(outputDictOrdered,open(outputFile,"w"))
print outputFile
sys.exit()

#print outputDictOrdered
outputTable = Table(outputDictOrdered)
outputTable.write(args.outputFile, overwrite=True)

