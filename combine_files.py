import numpy as np
from astropy.io import fits
import argparse
from astropy.table import Table
from collections import OrderedDict, defaultdict
import pickle

def pickFilter_R1000(wl):
    filterWlDict = {'f070lp':{'low':0.7,'high':1.1},\
                    'f100lp':{'low':1.1,'high':1.75},\
                    'f170lp':{'low':1.75,'high':3.0},\
                    'f290lp':{'low':3.0,'high':5.1}}
    gratingDict = {'f070lp':'g140m',\
                    'f100lp':'g140m',\
                    'f170lp':'g235m',\
                    'f290lp':'g395m'}
    filter = None
    grating = None
    for filt in filterWlDict.keys():
        if wl >= filterWlDict[filt]['low'] and wl < filterWlDict[filt]['high']:
            filter = filt
            grating = gratingDict[filt]

    return filter,grating

parser = argparse.ArgumentParser()
parser.add_argument(
                    '-f1, --file1',
                    help="first SN file",
                    action="store",
                    type=str,
                    dest="file1",
                    required=True
                    )

parser.add_argument(
                    '-f2, --file2',
                    help="second SN file",
                    action="store",
                    type=str,
                    dest="file2",
                    default = None,
                    required=False
                    )

parser.add_argument(
                    '-f3, --file3',
                    help="third SN file",
                    action="store",
                    type=str,
                    dest="file3",
                    default = None,
                    required=False
                    )

parser.add_argument(
                    '--folder1',
                    help="folder containing output spectra for first dither",
                    action="store",
                    type=str,
                    dest="folder1",
                    required=False,
                    default=None
                    )

parser.add_argument(
                    '--folder2',
                    help="folder containing output spectra for second dither",
                    action="store",
                    type=str,
                    dest="folder2",
                    required=False,
                    default=None
                    )

parser.add_argument(
                    '--folder3',
                    help="folder containing output spectra for third dither",
                    action="store",
                    type=str,
                    dest="folder3",
                    required=False,
                    default=None 
                    )

parser.add_argument(
                    '-i, --inputFile',
                    help="inputCatFile",
                    action="store",
                    type=str,
                    dest="inputCatFile",
                    required=True
                    )

parser.add_argument(
                    '-o, --outputFile',
                    help="outputCatFile",
                    action="store",
                    type=str,
                    dest="outputFile",
                    required=True
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
                    '--1910',
                    help="1910 rather than 1909 label",
                    action="store_true",
                    default=False,
                    required=False,
                    dest="l1910"
                    )

parser.add_argument(
                    '--es',
                    help="if extended source",
                    action="store_true",
                    default=False,
                    required=False,
                    dest="es"
                    )

parser.add_argument(
                    '--fits',
                    help="set this flag if reading in fits files rather than pickle files",
                    action="store_true",
                    default=False,
                    required=False,
                    dest="fits"
                    )

args = parser.parse_args()

nConfig = 4
grating = ['g140m','g140m','g235m','g395m']
filter = ['f070lp','f100lp','f170lp','f290lp']







lineDict = {}
lineDict['C4_1548']=0.1548
lineDict['C4_1551']=0.1551
lineDict['O3_1661']=0.1661
lineDict['O3_1666']=0.1666
lineDict['C3_1907']=0.1907
if args.l1910:
    lineDict['C3_1910'] = 0.1909
else:
    lineDict['C3_1909'] = 0.1909
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

lineDoublet_prism = np.array([1,1,2,2,3,3,4,4,0,0,0,0,0,0,0,5,5])
#write doublets as a single entry in output file - the next array says when to write the output so as not to duplicate information
writeOutput_prism = [1,0,1,0,1,0,1,0,1,1,1,1,1,1,1,1,0]

lineDoublet_R1000 = np.array([0,0,0,0,1,1,2,2,0,0,0,0,0,0,0,0,0])
writeOutput_R1000 = [1,1,1,1,1,0,1,0,1,1,1,1,1,1,1,1,1]



linesInOrder = ['C4_1548','C4_1551','O3_1661','O3_1666','C3_1907','C3_1909','O2_3726','O2_3729',\
                'Ne3_3869','O3_4363','HBaB_4861','O3_4959','O3_5007','HBaA_6563','N2_6584',\
                'S2_6716','S2_6731']

if args.l1910:
    linesInOrder[5] = 'C3_1910'

outputColumnOrder=['ID','Re_circ','n','pa']

outputLinesInOrder_R100 = ['C4_1548_C4_1551','O3_1661_O3_1666','C3_1907_C3_1909','O2_3726_O2_3729',\
                           'Ne3_3869','O3_4363','HBaB_4861','O3_4959_O3_5007','HBaA_6563','N2_6584',\
                           'S2_6716_S2_6731']
if args.l1910:
    outputLinesInOrder_R100[2] = 'C3_1907_C3_1910'

outputLinesInOrder_R1000 = ['C4_1548','C4_1551','O3_1661','O3_1666','C3_1907_C3_1909','O2_3726_O2_3729',\
                            'Ne3_3869','O3_4363','HBaB_4861','O3_4959','O3_5007','HBaA_6563','N2_6584',\
                            'S2_6716','S2_6731']
if args.l1910:
    outputLinesInOrder_R1000[4] = 'C3_1907_C3_1910'

outputLineLambdaInOrder_R1000 = [0.1548,0.1551,0.1661,0.1666,0.1907,0.3726,0.3869,\
                                 0.4363,0.4861,0.4959,0.5007,0.6563,0.6584,0.6716,0.6731]

outputContinuumInOrder= ['1550','1660','1910','3730','3870','4360','4860','4960','5010','6560','6720','D4000_1','D4000_2',\
                         'Calzetti_1','Calzetti_2','Calzetti_3','Calzetti_4','Calzetti_5','Calzetti_6','Calzetti_7',\
                         'Calzetti_8','Calzetti_9','Calzetti_10']

input = fits.open(args.inputCatFile)

if args.fits:
    dither1 = fits.open(args.file1)
    if args.file2 is not None and args.file3 is not None:
        dither2 = fits.open(args.file2)
        dither3 = fits.open(args.file3)
        dithers = [dither1[1].data,dither2[1].data,dither3[1].data]
    else:
        dithers = [dither1[1].data]
else:
    dither1 = pickle.load(open(args.file1,"r"))
    if args.file2 is not None and args.file3 is not None:
        dither2 = pickle.load(open(args.file2,"r"))
        dither3 = pickle.load(open(args.file3,"r"))
        dithers = [dither1,dither2,dither3]
    else:
        dithers = [dither1]


#first go through dithers and figure out which objects in which dither
ditherDict = defaultdict(list)
nObj = len(input['IDs'].data['ID_cat'])
for i in range(nObj):
    ditherDict['ID'].append(input['IDs'].data['ID_cat'][i])
    ditherDict['redshift'].append(input['GALAXY PROPERTIES'].data['redshift'][i])
    ditherArr = []
    idxArr = []
    for d in range(len(dithers)):
        tempIdx = np.where(dithers[d]['ID'] == input['IDs'].data['ID_cat'][i])[0]
        if len(tempIdx) > 0:
            ditherArr.append(d)
            idxArr.append(tempIdx)
    ditherDict['dithers'].append(ditherArr)
    ditherDict['idxArr'].append(idxArr)



outputDict = OrderedDict()
outputDict['ID'] = ditherDict['ID']
outputDict['redshift'] = ditherDict['redshift']
#for col in outputColumnOrder:
#    outputDict[col] = dither1[col]

if args.es:
    psStr = 'es'
else:
    psStr = 'ps'
if args.R100lines:
    for col in outputLinesInOrder_R100:
        fluxArr = []
        noiseArr = []
        snArr = []
        nDitherArr = []
        inputFluxArr = []
        for i in range(nObj):
            tempFlux = 0.
            tempNoise = 0.
            inputFluxArr.append(dithers[ditherDict['dithers'][i][0]][col+'_flux'][ditherDict['idxArr'][i][0][0]])
            for j,d in enumerate(ditherDict['dithers'][i]):
                tempFlux = tempFlux + dithers[d][col+'_f_R100_'+psStr+'_off_comb'][ditherDict['idxArr'][i][j]][0]
                if dithers[d][col+'_n_R100_'+psStr+'_off_comb'][ditherDict['idxArr'][i][j]] > -90:
                    tempNoise = tempNoise + dithers[d][col+'_n_R100_'+psStr+'_off_comb'][ditherDict['idxArr'][i][j]][0]**2
            if tempNoise > 0:
                tempSn = tempFlux/np.sqrt(tempNoise)
            else:
                tempSn = 0.
            fluxArr.append(tempFlux)
            noiseArr.append(np.sqrt(tempNoise))
            snArr.append(tempSn)
            nDitherArr.append(len(ditherDict['dithers'][i]))
                                
        outputDict[col+'_flux'] = inputFluxArr
        outputDict[col+'_R100_sn'] = snArr
        outputDict[col+'_R100_f'] = fluxArr
        outputDict[col+'_R100_n'] = noiseArr
        outputDict[col+'_n_dithers'] = nDitherArr
                                

if args.R1000lines:
    for wlIdx,col in enumerate(outputLinesInOrder_R1000):
        fluxArr = []
        noiseArr = []
        snArr = []
        nDitherArr = []
        inputFluxArr = []
        for i in range(nObj):
            redshift = ditherDict['redshift'][i]
            tempFlux = 0.
            tempNoise = 0.
            inputFluxArr.append(dithers[ditherDict['dithers'][i][0]][col+'_flux'][ditherDict['idxArr'][i][0]][0])
            #pick grating/filter configuration
            print outputLineLambdaInOrder_R1000[wlIdx]*(1+redshift)
            filter,grating = pickFilter_R1000(outputLineLambdaInOrder_R1000[wlIdx]*(1+redshift))
            if filter is not None:
                for j,d in enumerate(ditherDict['dithers'][i]):
                    print col+'_f_R1000_'+grating+'_'+filter+'_'+psStr+'_off_comb'
                    tempFlux = tempFlux + dithers[d][col+'_f_R1000_'+grating+'_'+filter+'_'+psStr+'_off_comb'][ditherDict['idxArr'][i][j]]
                    tempNoise = tempNoise + dithers[d][col+'_n_R1000_'+grating+'_'+filter+'_'+psStr+'_off_comb'][ditherDict['idxArr'][i][j]]**2
                tempSn = tempFlux/np.sqrt(tempNoise)
                fluxArr.append(tempFlux[0])
                noiseArr.append(np.sqrt(tempNoise[0]))
                snArr.append(tempSn[0])
            else:
                fluxArr.append(-99)
                noiseArr.append(-99)
                snArr.append(-99)
            nDitherArr.append(len(ditherDict['dithers'][i]))
        outputDict[col+'_flux'] = inputFluxArr
        outputDict[col+'_R1000_sn'] = snArr
        outputDict[col+'_R1000_f'] = fluxArr
        outputDict[col+'_R1000_n'] = noiseArr
        outputDict[col+'_n_dithers'] = nDitherArr


if args.R100continuum:
    outputDict['redshift'] = ditherDict['redshift']
    continuumOutputDict = defaultdict(list)
    if args.folder1 is None:
        print 'must define folders containing the saved spectra ', args.folder1, args.folder2, args.folder3
        sys.exit()
    else:
        folderArr = [args.folder1,args.folder2,args.folder3]
    for i in range(nObj):
        wlArr = None
        continuumOutputDict['n_dithers'].append(len(ditherDict['dithers'][i]))
        for d in ditherDict['dithers'][i]:
            pickleFile = folderArr[d]+str(ditherDict['ID'][i])+'_R100.p'
            data = pickle.load(open(pickleFile,'r'))
            if wlArr is None:
                wlArr = data['1d']['extracted_flux'][0]
                fluxArr = data['1d']['extracted_flux'][1]
                noiseArr = np.power(data['1d']['extracted_noise'][1],2)
            else:
                fluxArr = fluxArr + data['1d']['extracted_flux'][1]
                noiseArr = noiseArr + np.power(data['1d']['extracted_noise'][1],2)
        for cont in outputContinuumInOrder:
            inputWl = data['input']['scene'][0]['spectrum']['sed']['spectrum'][0]
            inputSpec = data['input']['scene'][0]['spectrum']['sed']['spectrum'][1]
            tempIdx = np.where((inputWl >= continuumDict[cont]['low']*(1+ditherDict['redshift'][i])) & (inputWl <= continuumDict[cont]['high']*(1+ditherDict['redshift'][i])))[0]
            if len(tempIdx) > 1:
                inputFlux = np.mean(inputSpec[tempIdx])
            else:
                inputFlux = 0.
            continuumOutputDict[cont+'_av_fnu_mjy'].append(inputFlux)

            tempIdx = np.where((wlArr >= continuumDict[cont]['low']*(1+ditherDict['redshift'][i])) & (wlArr <= continuumDict[cont]['high']*(1+ditherDict['redshift'][i])))[0]
            if len(tempIdx) > 1:
                tempsn = np.mean(fluxArr[tempIdx]/np.sqrt(noiseArr[tempIdx]))
            else:
                tempsn = 0.
            continuumOutputDict[cont+'_av_sn_pp'].append(tempsn)

    for cont in outputContinuumInOrder:
        outputDict[cont+'_av_fnu_mjy'] = continuumOutputDict[cont+'_av_fnu_mjy']
        outputDict[cont+'_av_sn_pp'] = continuumOutputDict[cont+'_av_sn_pp']
            
                
if args.R1000continuum:
    outputDict['redshift'] = ditherDict['redshift']
    continuumOutputDict = defaultdict(list)
    if args.folder1 is None:
        print 'must define folders containing the saved spectra ', args.folder1, args.folder2, args.folder3
        sys.exit()
    else:
        folderArr = [args.folder1,args.folder2,args.folder3]
    for cont in outputContinuumInOrder:
        for i in range(nObj):
            wlArr = None
            continuumOutputDict['n_dithers'].append(len(ditherDict['dithers'][i]))
            filter,grating = pickFilter_R1000(continuumDict[cont]['low']*(1+ditherDict['redshift'][i]))
            if filter is not None and grating is not None:
                for d in ditherDict['dithers'][i]:
                    pickleFile = folderArr[d]+str(ditherDict['ID'][i])+'_R1000_'+grating+'_'+filter+'.p'
                    data = pickle.load(open(pickleFile,'r'))
                    if wlArr is None:
                        wlArr = data['1d']['extracted_flux'][0]
                        fluxArr = data['1d']['extracted_flux'][1]
                        noiseArr = np.power(data['1d']['extracted_noise'][1],2)
                        inputWl = data['input']['scene'][0]['spectrum']['sed']['spectrum'][0]
                        inputSpec = data['input']['scene'][0]['spectrum']['sed']['spectrum'][1]
                        tempIdx = np.where((inputWl >= continuumDict[cont]['low']*(1+ditherDict['redshift'][i])) & (inputWl <= continuumDict[cont]['high']*(1+ditherDict['redshift'][i])))[0]
                        if len(tempIdx) > 1:
                            inputFlux = np.mean(inputSpec[tempIdx])
                        else:
                            inputFlux = 0.
                        continuumOutputDict[cont+'_av_fnu_mjy'].append(inputFlux)
                    else:
                        fluxArr = fluxArr + data['1d']['extracted_flux'][1]
                        noiseArr = noiseArr + np.power(data['1d']['extracted_noise'][1],2)

                tempIdx = np.where((wlArr >= continuumDict[cont]['low']*(1+ditherDict['redshift'][i])) & (wlArr <= continuumDict[cont]['high']*(1+ditherDict['redshift'][i])))[0]
                if len(tempIdx) > 1:
                    tempsn = np.mean(fluxArr[tempIdx]/np.sqrt(noiseArr[tempIdx]))
                else:
                    tempsn = 0.
            else:
                continuumOutputDict[cont+'_av_fnu_mjy'].append(0.)
                tempsn = 0.
            continuumOutputDict[cont+'_av_sn_pp'].append(tempsn)

    for cont in outputContinuumInOrder:
        outputDict[cont+'_av_fnu_mjy'] = continuumOutputDict[cont+'_av_fnu_mjy']
        outputDict[cont+'_av_sn_pp'] = continuumOutputDict[cont+'_av_sn_pp']
#        fluxArr = []
#        noiseArr = []
#        snArr = []
#        nDitherArr = []
#        inputFluxArr = []
#        for i in range(nObj):
#            tempFlux = 0.
#            tempNoise = 0.
#            inputFluxArr.append(dithers[ditherDict['dithers'][i][0]][col+'_flux'][ditherDict]['idxArr'][i][0])
#            for j,d in enumerate(ditherDict['dithers']):
#                tempFlux = tempFlux + dithers[d][col+'_f_R100_'+psStr+'_off_comb'][ditherDict]['idxArr'][j]
#                tempNoise = tempNoise + dithers[d][col+'_n_R100_'+psStr+'_off_comb'][ditherDict]['idxArr'][j]**2
#            tempSn = tempFlux/np.sqrt(tempNoise)
#            fluxArr.append(tempFlux)
#            noiseArr.append(np.sqrt(tempNoise))
#            snArr.append(tempSn)
#            nDitherArr.append(len(ditherDict['dithers'][i]))
#        
#        outputDict[col+'_flux'] = inputFluxArr
#        outputDict[col+'_R100_sn'] = snArr
#        outputDict[col+'_R100_f'] = fluxArr
#        outputDict[col+'_R100_n'] = noiseArr
#        outputDict[col_'_n_dithers'] = nDitherArr


outputTable = Table(outputDict)
outputTable.write(args.outputFile, overwrite=True)
