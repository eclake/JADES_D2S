from astropy.io import fits
import numpy as np
from astropy.table import Table
from collections import OrderedDict

data = fits.open("catalogue.fits")
snData = fits.open("../line_fluxes/sn_R1000_lines_combined.fits")

lineArr = ["C4_1548","C4_1551","O3_1661","O3_1666","C3_1907_C3_1910","O2_3726_O2_3729","Ne3_3869","O3_4363","HBaB_4861","O3_4959","O3_5007","HBaA_6563","N2_6584","S2_6716","S2_6731"]

eMPTout1 = np.genfromtxt("output_dither_MDC_00.txt", dtype=None, names=True)
eMPTout2 = np.genfromtxt("output_dither_MDC_01.txt", dtype=None, names=True)
eMPTout3 = np.genfromtxt("output_dither_MDC_02.txt", dtype=None, names=True)
outputDict = OrderedDict()
outputDict["ID"] = snData[1].data["ID"]
id_for_peter = np.zeros_like(snData[1].data["ID"])-1
for i in range(len(id_for_peter)):
  tempIdx = np.where(eMPTout1["ID_cat"] == snData[1].data["ID"][i])[0]
  if len(tempIdx) > 0:
    id_for_peter[i] = eMPTout1["ID_emma"][tempIdx]-1
  else:
    tempIdx = np.where(eMPTout2["ID_cat"] == snData[1].data["ID"][i])[0]
    if len(tempIdx) > 0:
      id_for_peter[i] = eMPTout2["ID_emma"][tempIdx]-1
    else:
      tempIdx = np.where(eMPTout3["ID_cat"] == snData[1].data["ID"][i])[0]
      if len(tempIdx) > 0:
        id_for_peter[i] = eMPTout3["ID_emma"][tempIdx]-1
outputDict["ID_for_peter"] = id_for_peter
outputDict["redshift"] = data['GALAXY PROPERTIES'].data['redshift']

for line in lineArr:
  flux = data['SPECTRAL INDICES'].data[line]
  flux[tempIdx] = -99
  sn = snData[1].data[line+"_R1000_sn"]
  outputDict[line+"_noiseless_flux"] = flux
  print flux.shape, sn.shape
  tempIdx = np.where(sn == 0)[0]
  err = np.abs(flux/sn)
  err[tempIdx] = -99
  outputDict[line+"_flux"] = flux + np.random.normal(size=len(err))*err
  outputDict[line+"_flux"][tempIdx] = -99
  outputDict[line+"_err"] = err
  


outputTable = Table(outputDict)
outputTable.write("../line_fluxes/perturbedFluxes.fits",overwrite=True)
