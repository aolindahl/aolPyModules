import numpy as np


# Define the basic configuration
basicTofConfig = {
    "acqCh": 1, 
    "baselineSubtraction":'roi',
    "baselineEnd_us": 1.5, 
    "calibFile": "/reg/neh/home/alindahl/pythonModules/tofCalibDefault.json", 
    "filterMethod": "none", # wavelet, average, wienerDeconv 
    "filterWaveletLevels": 100, 
    "filterWinerSNR": 1,
    "filterWinerResponse":1,
    "filterAverageNumPoints":100,
    "detectorSource": "DetInfo(AmoETOF.0:Acqiris.0)", 
    #"tMax_us": 2.5, 
    #"tMin_us": 1, 
    "tMax_us": 1.57, 
    "tMin_us": 1.48, 
    "tSlice": True 
}


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

timeRoi0_us_common = [1.533, 1.539]

if timeRoi0_us_common is not None:
        timeRoi0_us = {}
        for i in range(16):
                timeRoi0_us[i] = timeRoi0_us_common


timeRoi1_us_common = [1.515, 1.521]

if timeRoi1_us_common is not None:
	timeRoi1_us = {}
	for i in range(16):
		timeRoi1_us[i] = timeRoi1_us_common




# Make copies of the basic configuration for each of the detectors
tofConfigList = [None] * 16

def makeTofConfigList(online=True):
    global tofConfigList
    for i in range(16):
        tofConfigList[i] = basicTofConfig.copy()
        tofConfigList[i]['timeRoi0_us'] = timeRoi0_us[i]
        tofConfigList[i]['timeRoi1_us'] = timeRoi1_us[i]
        if online:
            tofConfigList[i]['acqCh'] = acqCh[i]

makeTofConfigList(online=True)

minE_eV = 50
maxE_eV = 1000
nEnergyBins = 1024

energyScaleBinLimits = np.linspace(minE_eV, maxE_eV, nEnergyBins+1)

fitMask = np.array([0,1, 2, 3, 4,5,6,7,8,9,10,11,12,13,14,15])
#fitMask = np.array([0,1, 2, 3, 4,5,6,7,9,10,11,12,13,14,15])


# For CookieBox class debugging
domainToDisplay = 'Time'


# Stuff below are used in the debugging of the tofData class

dataSource = 'exp=amoc8114:run=31'
nEvents = 10
useZmq = False
