import numpy as np
import aolUtil

offline_source = 'exp=amoi0314:run=15'

cbSourceString = "DetInfo(AmoETOF.0:Acqiris.0)"
tMin_us =  1.48
tMax_us =  1.6
baselineSubtraction = 'early'
baselineEnd_us = 1.5

# Define the basic configuration
basic_tof_config = {
    "acqCh": 10, 
    "baselineSubtraction":"early", # 'early', 'roi', 'none' 
    "baselineEnd_us": 1.5, 
    "calibFile": ("tof_calib_default.json"), 
    "filterMethod": "none", # wavelet, average, wienerDeconv 
    "filterWaveletLevels": 10, 
    "filterWinerSNR": 1,
    "filterWinerResponse":1,
    "filterAverageNumPoints":4,
    "detectorSource": "DetInfo(AmoETOF.0:Acqiris.0)", 
    "tMax_us": 1.6, 
    "tMin_us": 1.48, 
    "tSlice": True 
}

lclsConfig = {
	'lcls_photonEnergyA':1,
	'lcls_photonEnergyB':1}
lclsConfig = aolUtil.struct(lclsConfig)

retardationPV = 'AMO:R14:IOC:10:VHS0:CH0:VoltageMeasure'

# Acqiris channel asignment
acqCh = {
        0:0,
        1:1,
        2:2,
        3:3,
        4:4,
        5:5,
        6:6,
        7:7,
        8:8,
        9:9,
        10:10,
        11:11,
        12:12,
        13:13,
        14:14,
        15:15
        }

timeRoi0_us_common = [1.522, 1.538]	#red
timeRoi0_us = [timeRoi0_us_common]*16

timeRoi0Bg_us_common = [1.5, 1.51]
timeRoi0Bg_us = [timeRoi0Bg_us_common]*16

timeRoi1_us_common = [1.515, 1.522]	#green
timeRoi1_us = [timeRoi1_us_common]*16

energyRoi0_eV_common = [40, 60]
energyRoi0_eV = [energyRoi0_eV_common]*16



# Make copies of the basic configuration for each of the detectors
tofConfigList = [None] * 16

def makeTofConfigList(online=True):
    global tofConfigList
    for i in range(16):
        tofConfigList[i] = basic_tof_config.copy()
	tofConfigList[i]['calibFile'] = ('/reg/neh/operator/amoopr/'
		    + 'amoi0114/psana/tofCalibs/tof{}Calib.json'.format(i+1)) 
	if online:
            tofConfigList[i]['acqCh'] = acqCh[i]

makeTofConfigList(online=True)

minE_eV = 50
maxE_eV = 1000
nEnergyBins = 256

energyScaleBinLimits = np.linspace(minE_eV, maxE_eV, nEnergyBins + 1)


fitMask = np.array([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15])
#fitMask = np.array([0,1,2,3,4,5,6, 8,9,10,11, 13,14,15])
#fitMask = np.array([0,1,2,3,5,6,7,8,9,10,11,13,14,15])

boolFitMask = np.array([i in fitMask for i in range(16)])
nanFitMask = np.array([1 if b else np.nan for b in boolFitMask])


# For CookieBox class debugging
domainToDisplay = 'Time'


# Stuff below are used in the debugging of the tofData class

dataSource = 'exp=amoc8114:run=31'
nEvents = 10
useZmq = False
