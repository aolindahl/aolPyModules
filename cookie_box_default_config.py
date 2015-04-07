from numpy import array, linspace
from os.path import join, dirname
import aolUtil

offline_source = 'exp=amoi0314:run=15'

#bSourceString = 'DetInfo(AmoETOF.0:Acqiris.0)'
#tMin_us =  0.48
#tMax_us =  0.6
#baselineSubtraction = 'early'
#baselineEnd_us = 0.5

# Define the basic configuration
basic_tof_config = {
    "acqCh": 10, 
    "baselineSubtraction":"early", # 'early', 'roi', 'none' 
    "baselineEnd_us": 0.5, 
    "calib_file": join(dirname(__file__), "tof_calib_default.json"), 
    "filterMethod": "none", # wavelet, average, wienerDeconv 
    "filterWaveletLevels": 10, 
    "filterWinerSNR": 1,
    "filterWinerResponse":1,
    "filterAverageNumPoints":4,
    "detectorSource": 'DetInfo(AmoETOF.0:Acqiris.0)', 
    "tMax_us": 0.6, 
    "tMin_us": 0.48, 
    "tSlice": True 
}

lclsConfig = {
	'lcls_photonEnergyA':1,
	'lcls_photonEnergyB':1}
lclsConfig = aolUtil.struct(lclsConfig)

retardationPV = 'AMO:R14:IOC:10:VHS0:CH0:VoltageMeasure'

# Acqiris channel asignment
acqiris_setup = {
        0:['DetInfo(AmoETOF.0:Acqiris.0)', 0],
        1:['DetInfo(AmoETOF.0:Acqiris.0)', 1],
        2:['DetInfo(AmoETOF.0:Acqiris.0)', 2],
        3:['DetInfo(AmoETOF.0:Acqiris.0)', 3],
        4:['DetInfo(AmoETOF.0:Acqiris.0)', 4],
        5:['DetInfo(AmoETOF.0:Acqiris.0)', 5],
        6:['DetInfo(AmoETOF.0:Acqiris.0)', 6],
        7:['DetInfo(AmoETOF.0:Acqiris.0)', 7],
        8:['DetInfo(AmoETOF.0:Acqiris.0)', 8],
        9:['DetInfo(AmoETOF.0:Acqiris.0)', 9],
        10:['DetInfo(AmoETOF.0:Acqiris.0)', 10],
        11:['DetInfo(AmoETOF.0:Acqiris.0)', 11],
        12:['DetInfo(AmoETOF.0:Acqiris.0)', 12],
        13:['DetInfo(AmoETOF.0:Acqiris.0)', 13],
        14:['DetInfo(AmoETOF.0:Acqiris.0)', 14],
        15:['DetInfo(AmoETOF.0:Acqiris.0)', 15],
        }

time_roi_0_us_common = [0.522, 0.538]	#red
time_roi_0_us = [time_roi_0_us_common]*16

time_roi_0_bg_us_common = [0.5, 0.51]
time_roi_0_bg_us = [time_roi_0_bg_us_common]*16

time_roi_1_us_common = [0.515, 0.522]	#green
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
        tof_config_list[i]['detectorSource']
        tof_config_list[i]['acqCh'] = acqiris_setup[i][1]

        tof_config_list[i]['time_roi_0_us'] = time_roi_0_us[i]
        tof_config_list[i]['time_roi_0_bg_us'] = time_roi_0_bg_us[i]
        tof_config_list[i]['time_roi_1_us'] = time_roi_1_us[i]

makeTofConfigList()

minE_eV = 50
maxE_eV = 1000
n_energy_bins= 256

energy_scale_eV = linspace(minE_eV, maxE_eV, 2*n_energy_bins + 1)[1::2]


fitMask = array([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15])
#fitMask = array([0,1,2,3,4,5,6, 8,9,10,11, 13,14,15])
#fitMask = array([0,1,2,3,5,6,7,8,9,10,11,13,14,15])

boolFitMask = array([i in fitMask for i in range(16)])
nanFitMask = array([1 if b else np.nan for b in boolFitMask])


# For CookieBox class debugging
domainToDisplay = 'Time'


# Stuff below are used in the debugging of the tofData class

dataSource = 'exp=amoc8114:run=31'
nEvents = 10
useZmq = False
