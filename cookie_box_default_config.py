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
    "calib_file": ("tof_calib_default.json"), 
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

time_roi_0_us_common = [1.522, 1.538]	#red
time_roi_0_us = [time_roi_0_us_common]*16

time_roi_0Bg_us_common = [1.5, 1.51]
time_roi_0Bg_us = [time_roi_0Bg_us_common]*16

time_roi_1_us_common = [1.515, 1.522]	#green
time_roi_1_us = [time_roi_1_us_common]*16

energy_roi_0_eV_common = [40, 60]
energy_roi_0_eV = [energy_roi_0_eV_common]*16



# Make copies of the basic configuration for each of the detectors
tof_config_list = [None] * 16

def makeTofConfigList():
    global tof_config_list
    for i in range(16):
        tof_config_list[i] = basic_tof_config.copy()
	    # tof_config_list[i]['calib_file'] = ('/reg/neh/operator/amoopr/'
	    #	    + 'amoi0114/psana/tofCalibs/tof{}Calib.json'.format(i+1)) 
        tof_config_list[i]['calib_filte'] = 'tof_calib_default.json'
        tof_config_list[i]['acqCh'] = acqCh[i]

makeTofConfigList()

minE_eV = 50
maxE_eV = 1000
n_energy_bins= 256

energy_scale_eV = np.linspace(minE_eV, maxE_eV, 2*n_energy_bins + 1)[1::2]


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
