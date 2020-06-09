from astropy.io import fits
import numpy as np
from pandeia.engine.calc_utils import build_default_calc
from pandeia.engine.perform_calculation import perform_calculation
import pandeia.engine
from astropy.table import Table
from collections import OrderedDict
import copy
import pylab
import argparse
from scipy.interpolate import interp1d
import os

parser = argparse.ArgumentParser()

parser.add_argument(
                    '-b', '--beagle-output-file',
                    help="the beagle output file containing simulated spectra",
                    action="store",
                    type=str,
                    dest="beagleOutputFile",
                    required=True
                    )

parser.add_argument(
                    '-s', '--sizes-file',
                    help="file containing size info.  Must have same number of rows as BEAGLE output file.  If not supplied, objects will be simulated as point source.",
                    action="store",
                    type=str,
                    dest="sizesFile",
                    required=False,
                    default=None
                    )

parser.add_argument(
                    '--sn-off',
                    help="just output all spectra",
                    action="store_true",
                    dest="snOff",
                    required=False,
                    default=False
                    )

parser.add_argument(
                    '--interpolate-on',
                    help="DANGER - if BEAGLE resolution too low, will interpolate report target output - this is not correct and should not be done for simulating spectra in general",
                    action="store_true",
                    dest="interpolate",
                    required=False,
                    default=False
                    )

parser.add_argument(
                    '-o', '--output-folder',
                    help="output folder where fits files containing noisy spectra are written",
                    action="store",
                    type=str,
                    dest="outputFolder",
                    required=False,
                    default=None
                    )



parser.add_argument(
                    '-n',
                    help="number of objects in file to make spectra of",
                    action="store",
                    type=str,
                    dest="n",
                    required=False,
                    default=None
                    )


parser.add_argument(
                    '--id-file',
                    help=".txt file containing list of IDs to create output spectra for",
                    action="store",
                    type=str,
                    dest="idFile",
                    required=False,
                    default=None
                    )




parser.add_argument(
                    '--add-lines-separately',
                    help="read in the continuum spectrum and add the lines after",
                    action="store_true",
                    dest="addLines",
                    required=False,
                    default=False
                    )

args = parser.parse_args()

#This script was adapted from a short script written by Michael Maseda demonstrating how to
#set up emission line S/N calculations in pandeia.

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
    if (inputs['re_circ'] > 0):
        sersic=inputs['sersic_n']
        re = inputs['re_maj']
        ellip=1.-inputs['axis_ratio']
        pa=inputs['position_angle']
        major_axis = re
        minor_axis = re*inputs['axis_ratio']
        #pandeia used to use scale lengths, but now is fine with half-light radii
        scene['shape']={'geometry':'sersic','major': major_axis,'minor':minor_axis,'sersic_index':sersic}
#        #pandiea wants scale lengths not half-light radii
#        major_axis,minor_axis=majorminor(sersic,rc,ellip)
#        scene['position'] = {'x_offset':xoff, 'y_offset': yoff, 'orientation': pa, 'position_parameters':['x_offset','y_offset','orientation']}
#        scene['shape']={'geometry':'sersic','major': major_axis,'minor':minor_axis,'sersic_index':sersic}
#        print scene['shape']
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
    
    if args.addLines:
      emission_line_array = []
      for i,f in enumerate(inputs['flux']):
        flux=f
        wave=inputs['wave'][i]
        if wave > filterWlDict[filt]['low'] and wave < filterWlDict[filt]['high']:
          emission_line={}
          emission_line['emission_or_absorption']='emission'
          emission_line['center']=wave
          #TODO: check units...
          emission_line['strength']=flux
          emission_line['profile']='gaussian'
          #assume the line is basically unresolved, i.e. 50 km/s (FWHM)
          emission_line['width']=50.#50.
          emission_line_array.append(emission_line)
      scene['spectrum']['lines']=emission_line_array

    configuration['scene'][0]=scene
    report=perform_calculation(configuration)
    return report


#
exposureDict = {'DEEP':{'clear':{'ngroup':19,'nint':2,'nexp':36},'f070lp':{'ngroup':19,'nint':2,'nexp':9},'f100lp':{'ngroup':19,'nint':2,'nexp':9},'f170lp':{'ngroup':19,'nint':2,'nexp':9},'f290lp':{'ngroup':19,'nint':2,'nexp':9}},\
                'MEDIUM':{'clear':{'ngroup':13,'nint':1,'nexp':9},'f070lp':{'ngroup':13,'nint':1,'nexp':9},'f100lp':{'ngroup':13,'nint':1,'nexp':9},'f170lp':{'ngroup':13,'nint':1,'nexp':9},'f290lp':{'ngroup':13,'nint':1,'nexp':9}},\
                'MEDIUM_HST':{'clear':{'ngroup':16,'nint':1,'nexp':6},'f070lp':{'ngroup':13,'nint':1,'nexp':6},'f100lp':{'ngroup':13,'nint':1,'nexp':6},'f170lp':{'ngroup':13,'nint':1,'nexp':6},'f290lp':{'ngroup':16,'nint':1,'nexp':6}},\
                'DEEP_WORST_CASE':{'clear':{'ngroup':19,'nint':2,'nexp':12}},\
                'MEDIUM_WORST_CASE':{'clear':{'ngroup':13,'nint':1,'nexp':3}},\
                'MEDIUM_HST_WORST_CASE':{'clear':{'ngroup':16,'nint':1,'nexp':3}}}
#exposureDict = {'DEEP':{'clear':{'ngroup':19,'nint':2,'nexp':36}},\
#                'MEDIUM':{'clear':{'ngroup':13,'nint':1,'nexp':9}},\
#                'MEDIUM_HST':{'clear':{'ngroup':16,'nint':1,'nexp':6}},\
#                'DEEP_WORST_CASE':{'clear':{'ngroup':19,'nint':2,'nexp':12}},\
#                'MEDIUM_WORST_CASE':{'clear':{'ngroup':13,'nint':1,'nexp':3}},\
#                'MEDIUM_HST_WORST_CASE':{'clear':{'ngroup':16,'nint':1,'nexp':2}}}
exposureDict = {'DEEP':{'clear':{'ngroup':19,'nint':2,'nexp':36}},\
                'MEDIUM':{'clear':{'ngroup':13,'nint':1,'nexp':9}}}
#exposureDict = {'MEDIUM_HST':{'clear':{'ngroup':16,'nint':1,'nexp':6}},\
#                'DEEP_WORST_CASE':{'clear':{'ngroup':19,'nint':2,'nexp':12}},\
#                'MEDIUM_WORST_CASE':{'clear':{'ngroup':13,'nint':1,'nexp':3}},\
#                'MEDIUM_HST_WORST_CASE':{'clear':{'ngroup':16,'nint':1,'nexp':3}}}
#exposureDict = {'DEEP':{'clear':{'ngroup':19,'nint':2,'nexp':36},'f070lp':{'ngroup':19,'nint':2,'nexp':9},'f100lp':{'ngroup':19,'nint':2,'nexp':9},'f170lp':{'ngroup':19,'nint':2,'nexp':9},'f290lp':{'ngroup':19,'nint':2,'nexp':9}}}
#exposureDict = {'DEEP':{'f100lp':{'ngroup':19,'nint':2,'nexp':9},'f170lp':{'ngroup':19,'nint':2,'nexp':9},'f290lp':{'ngroup':19,'nint':2,'nexp':9}}}
exposureDict = {'DEEP':{'clear':{'ngroup':19,'nint':2,'nexp':36}}}
exposureDict = {'DEEP':{'f070lp':{'ngroup':19,'nint':2,'nexp':9},\
                'f100lp':{'ngroup':19,'nint':2,'nexp':9},\
                'f170lp':{'ngroup':19,'nint':2,'nexp':9},\
                'f290lp':{'ngroup':19,'nint':2,'nexp':9}}}

cat = fits.open(args.beagleOutputFile)
sizesSupplied = False
if args.sizesFile is not None:
  input = fits.open(args.sizesFile)
  sizesSupplied = True

grating = ['g140m','g140m','g235m','g395m']
filter = ['f070lp','f100lp','f170lp','f290lp']
filterDict = {'clear':'prism',\
                'f070lp':'g140m',\
                'f100lp':'g140m',\
                'f170lp':'g235m',\
                'f290lp':'g395m'}

filterWlDict = {'clear':{'low':0.7,'high':5.1},\
                'f070lp':{'low':0.7,'high':1.2},\
                'f100lp':{'low':1.1,'high':1.8},\
                'f170lp':{'low':1.7,'high':3.1},\
                'f290lp':{'low':3.0,'high':5.1}}

if args.addLines:
  wl = cat['CONTINUUM SED WL'].data[0][0]
  specArr = cat['CONTINUUM SED'].data
  lineWlList = []
  lineLabelList = []
  for col in cat['HII EMISSION'].data.dtype.names:
    if 'flux' in col:
      print(col)
      tmp = col.split("_")
      lineWlList.append(np.float(tmp[1].split("_")[0]))
      lineLabelList.append(col)
      print lineWlList[-1], lineLabelList[-1]
else:
  wl = cat['FULL SED WL'].data[0][0]
  specArr = cat['FULL SED'].data
c_light = 2.99792e+18 # Ang/s
nObj = len(cat['GALAXY PROPERTIES'].data['redshift'])

idArr = np.fromiter((x for x in range(nObj)),np.int)
if args.idFile is not None:
    temp = np.genfromtxt(args.idFile,dtype=None, names=True)
    idArr = temp['ID']
    print idArr
elif args.n is not None:
    if args.n < nObj:
        idArr = fromiter((x for x in np.range(args.n)),np.int)

print 'idArr: ', idArr
for i in idArr:
    z = cat['GALAXY PROPERTIES'].data['redshift'][i]
    tempWl = wl*(1+z)
    tempSpec = specArr[i,:]/(1+z)
    tempSpec = tempWl**2/c_light*tempSpec*1E23*1E3 #in mJy
    tempWl = tempWl/1E4 #in microns
    deduped = {w: s for w, s in reversed(zip(list(tempWl), list(tempSpec)))}
    tempWl = np.asarray(sorted(deduped.keys()))
    tempSpec = np.asarray([deduped[k] for k in tempWl])
    inputs = {}
    inputs['wl'] = tempWl
    inputs['spec'] = tempSpec
    inputs['xoff'] = 0.
    inputs['yoff'] = 0.
#    pylab.plot(tempWl, tempSpec)
#    pylab.xlim(0.7,5)
#    pylab.xlabel("wl/micron")
#    pylab.ylabel("fnu/Jy")
#    pylab.tight_layout()
#    pylab.show()
    if sizesSupplied:
      inputs['axis_ratio'] = input[1].data['axis_ratio'][i]
      inputs['sersic_n'] = input[1].data['sersic_n'][i]
      inputs['position_angle'] = input[1].data['position_angle'][i]
      inputs['re_circ'] = input[1].data['Re_maj'][i]*np.sqrt(input[1].data['axis_ratio'][i])
      inputs['re_maj'] = input[1].data['Re_maj'][i]
    else:
      inputs['axis_ratio'] = 1.
      inputs['sersic_n'] = -99
      inputs['position_angle'] = 0.
      inputs['re_circ'] = -99
    inputs['onSource'] = [False,True,False]
    inputs['slitletShape'] = [[0,-1],[0,0],[0,1]]
    if args.addLines:
      lineFluxArr = []
      for line in lineLabelList:
        lineFluxArr.append(cat['HII EMISSION'].data[line][i])
      inputs['wave'] = np.array(lineWlList)*(1+z)/1.E4
      inputs['flux'] = lineFluxArr

    #produce a mock spectrum for extended and point source, for each of the filter/grating configurations
    for key in exposureDict.keys():
        for filt in exposureDict[key].keys():

            if filt == 'clear':
                folder = args.outputFolder+key+"_R100/"
            else:
                folder = args.outputFolder+key+"_R1000/"
            idStr = str(i)
            if sizesSupplied:
              if 'ID' in input[1].data.dtype.names:
                idStr = str(input[1].data['ID'][i])
            if sizesSupplied:
                outputFile = folder+'/'+idStr+'_'+filt+'_'+filterDict[filt]+'_extended.fits'
            else:
                outputFile = folder+'/'+idStr+'_'+filt+'_'+filterDict[filt]+'.fits'
            print filt
            
            #test = os.path.exists(outputFile)
            test = False
            if not test:
              print inputs
              report = sn_user_spec(inputs, disperser = filterDict[filt], filt = filt, ngroup = exposureDict[key][filt]['ngroup'], nint = exposureDict[key][filt]['nint'], nexp = exposureDict[key][filt]['nexp'])
              outputDict = OrderedDict()
              outputDict['wl'] = report['1d']['extracted_flux'][0]
  #            pylab.figure()
  #            pylab.plot(report['1d']['target'][0],report['1d']['sn'][0])
  #            pylab.show()
              print report['1d']['target'][0][:10]
              print report['1d']['sn'][0][:10]
#              pylab.plot(report['1d']['sn'][0], report['1d']['sn'][1])
#              pylab.xlabel("wl/micron")
#              pylab.ylabel("pandeia S/N")
#              pylab.tight_layout()
#              pylab.show()
#              pylab.plot(report['1d']['extracted_contamination'][0], report['1d']['extracted_contamination'][1], label="extracted contamination")
#              pylab.plot(report['1d']['extracted_bg_only'][0], report['1d']['extracted_bg_only'][1], c='r', label="extracted bg only")
#              pylab.plot(report['1d']['extracted_bg_total'][0], report['1d']['extracted_bg_total'][1], c='g', label="extracted bg total")
#              pylab.plot(report['1d']['extracted_flux_plus_bg'][0], report['1d']['extracted_flux_plus_bg'][1], c='orange', label="extracted flux plus bg")
#              pylab.xlabel("wl/micron")
#              pylab.ylabel("extracted flux, pandeia units")
#              pylab.legend()
#              pylab.tight_layout()
#              pylab.show()
              if len(report['1d']['target'][0]) != len(report['1d']['sn'][0]):
                  if args.interpolate:
                      f = interp1d(report['1d']['target'][0],report['1d']['target'][1])
  #                    pylab.figure()
  #                    pylab.plot(report['1d']['target'][0],report['1d']['target'][1])
                      targetNew = f(report['1d']['sn'][0])
                      report['1d']['target'][0] = report['1d']['sn'][0]
                      report['1d']['target'][1] = targetNew
  #                    pylab.plot(report['1d']['target'][0],report['1d']['target'][1], c='r')
  #                    pylab.show()
                  else:
                      print 'error: length of target and sn report arrays differ, check input resolution!'
                      sys.exit()
          
              print len(report['1d']['target'][1]), len(report['1d']['sn'][1])
  #            sys.exit()
              noise = report['1d']['target'][1]/report['1d']['sn'][1]
              noisy_spectrum = report['1d']['target'][1]+np.random.normal(size=len(report['1d']['target'][1]))*noise
              outputDict['fnu'] = noisy_spectrum
              outputDict['fnu_err'] = noise
              outputDict['sn'] = report['1d']['sn'][1]
              outputDict['pand_extracted_flux'] = report['1d']['extracted_flux'][1]
              outputDict['pand_extracted_noise'] = report['1d']['extracted_noise'][1]
              outputDict['fnu_noiseless'] = report['1d']['target'][1]
              tempIdx = np.where((report['1d']['sn'][0] >= 0.6450*(1+z)) & (report['1d']['sn'][0] <= 0.6650*(1+z)))[0]
              if args.snOff:
                writeOutput = True
              else:
                writeOutput = False
                print tempIdx
                if len(tempIdx) == 0: #use region around OII 3727
                  print 0.6450*(1+z), 0.3727*(1+z)
                  tempIdx = np.where((report['1d']['sn'][0] >= 0.3600*(1+z)) & (report['1d']['sn'][0] <= 0.3800*(1+z)))[0]
                  print tempIdx
                
                sn = np.max(report['1d']['sn'][1][tempIdx])
                if sn > 3:
                  writOutput=True
          
              print '*** ', filt, writeOutput
              if writeOutput:
                outputTable = Table(outputDict)
                print outputFile
                outputTable.write(outputFile,overwrite=True)
#            inputs2 = copy.deepcopy(inputs)
#            inputs2['re_maj'],inputs2['re_min'] = -1,-1
#            report = sn_user_spec(inputs2, disperser = filterDict[filt], filt = filt, ngroup = exposureDict[key][filt]['ngroup'], nint = exposureDict[key][filt]['nint'], nexp = exposureDict[key][filt]['nexp'])
#            outputDict = OrderedDict()
#            outputDict['wl'] = report['1d']['extracted_flux'][0]
#            outputDict['extracted_flux'] = report['1d']['extracted_flux'][1]
#            outputDict['extracted_noise'] = report['1d']['extracted_noise'][1]
#            outputDict['sn'] = report['1d']['sn'][1]
#            if filt == 'clear':
#                folder = "/Users/user/Documents/JWST/data-to-science/"+key+"_R100/"
#            else:
#                folder = "/Users/user/Documents/JWST/data-to-science/"+key+"_R1000/"
#            outputFile = folder+str(input[1].data['ID'][i])+'_'+filt+'_'+filterDict[filt]+'_ps.fits'
#            outputTable = Table(outputDict)
#            outputTable.write(outputFile,overwrite=True)
