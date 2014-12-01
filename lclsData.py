# -*- coding: utf-8 -*-
"""
Created on Tue May  6 17:56:30 2014

@author: alindahl
"""

#from setupEnvironment import *
import psana
import numpy as np

EBeamTypeList = (psana.Bld.BldDataEBeamV0,
                 psana.Bld.BldDataEBeamV1,
                 psana.Bld.BldDataEBeamV2, 
                 psana.Bld.BldDataEBeamV3,
                 psana.Bld.BldDataEBeamV4,
                 psana.Bld.BldDataEBeamV5,
                 psana.Bld.BldDataEBeamV6)
                 


class LCLSdata:
    def __init__(self, config, quiet=True):
        self._evt = None

        # EBeam
        self._EBeamSource = psana.Source('BldInfo(EBeam)')
        self._EBeamType = None
        self._phaseCavitySource = psana.Source('BldInfo(PhaseCavity)')
        
        # Evr
        self._EvrSource = psana.Source('DetInfo(NoDetector.0:Evr.0)')
        self._EvrType = psana.EvrData.DataV3

        # photon beam
        #photonCalib = loadConfig(calibFile, quiet=quiet)
        #self._photonA = photonCalib.A
        #self._photonB = photonCalib.B
        self._photonA = config.lcls_photonEnergyA
        self._photonB = config.lcls_photonEnergyB

        self._feeSource = psana.Source('BldInfo(FEEGasDetEnergy)')
        
        self._idData = None


    def setEvent(self, evt):
        self._evt = evt
        self._eBeamData = None
        self._evrData = None
        self._idData = None


    ###########################################
    # EBeam stuff
    
    def getEBeamEnergyL3_MeV(self):
        EBeamObject = self.getEBeamObject()
        if EBeamObject == None:
            return np.nan
        return EBeamObject.ebeamL3Energy()
        
    def getEBeamEnergyBC2_MeV(self):
        EBeamObject = self.getEBeamObject()
        if EBeamObject == None:
            return np.nan
        # The BC2 energy is calculated using the dispersion. The vispersion
        # value for BC2 is -3.647 mm at the nominal beam energy of 5 GeV [Email
        # from Timothy Maxwell to Anton Lindahl on June 2 2014]
        return ( EBeamObject.ebeamEnergyBC2() / -364.7 + 1 ) * 5e3
        
    def getEBeamCharge_nC(self):
        EBeamObject = self.getEBeamObject()
        if EBeamObject == None:
            return np.nan
        return EBeamObject.ebeamCharge()
        
    def getEBeamPkCurrentBC2_A(self):
        EBeamObject = self.getEBeamObject()
        if EBeamObject == None:
            return np.nan
        return EBeamObject.ebeamPkCurrBC2()
        
        
    def _determineEBeamType(self):
        if self._EBeamType is None:
            for type in EBeamTypeList:
                data = self._evt.get(type, self._EBeamSource)
                if data is not None:
                    self._EBeamType = type
                    break

                
    def getEBeamObject(self):
        # Initialize the EBeam type
        if self._EBeamType is None:
            self._determineEBeamType()
        if self._eBeamData is None and self._EBeamType != None:
            self._eBeamData = self._evt.get(self._EBeamType, self._EBeamSource)
        return self._eBeamData
        
    ###########################################
    # Phase cavity

    def getPhaseCavityTimes_ps(self):
        pc = self._evt.get(psana.Bld.BldDataPhaseCavity,
                self._phaseCavitySource)
        if pc is None:
            return [np.nan] * 2
        return [pc.fitTime1(), pc.fitTime2()]

        
    ###########################################
    # Evr stuff
        
    def getEvrObject(self):
        if self._evrData is None:
            self._evrData = self._evt.get(self._EvrType, self._EvrSource)
        return self._evrData
    
    def getEvrCodes(self):
        if self.getEvrObject() == None:
            return None
        return [fifo.eventCode() for fifo in self._evrData.fifoEvents()]
        
    def evrCodeInEvent(self, evr):
        return evr in self.getEvrFifoEventCodes(evt)
        
        
    
    ###########################################
    # x-ray stuff

    def getPhotonEnergy(self):
        E = self.getEBeamEnergyL3_MeV()
        I = self.getEBeamCharge_nC()
        return self._photonA * (E  + self._photonB * I)**2

    def getPulseEnergy_mJ(self):
        fee = self._evt.get(psana.Bld.BldDataFEEGasDetEnergyV1, self._feeSource)
        if fee is None:
            return [np.nan for i in range(4)]
        return [ fee.f_11_ENRC(), fee.f_12_ENRC(), fee.f_21_ENRC(),
                fee.f_22_ENRC()]

    
    
    ###########################################
    # info stuff
    
    def getEventTime(self):
        if self._idData == None:
            self._idData = self._evt.get(psana.EventId)
        if self._idData == None:
            return np.nan
        timeStamp = self._idData.time()
        return timeStamp[0] + timeStamp[1] * 1e-9

    def getEventFiducial(self):
        if self._idData == None:
            self._idData = self._evt.get(psana.EventId)
        if self._idData == None:
            return 0
        return self._idData.fiducials()
    

if __name__ == '__main__':
    from configuration import loadConfiguration

    print 'Testing the LCLSdata class.'

    config =\
            loadConfiguration('/reg/neh/home/alindahl/amoc8114/configFiles/config5-27.json')
    lcls = LCLSdata(config)
    
    ds = psana.DataSource("/reg/neh/home/alindahl/wolfiExp/xtc/e24-r0284-s01-c00.xtc")
    for num, evt in enumerate(ds.events()):
        lcls.setEvent(evt)
        if num >= 10:
            EBeam = lcls.getEBeamObject()
            break
        print 'event {}, {} MeV, {} pC, {} A,\nE_photon {}'.format(
            num,
            lcls.getEBeamEnergyL3_MeV(),
            lcls.getEBeamCharge_nC()*1e3,
            lcls.getEBeamPkCurrentBC2_A(),
            lcls.getPhotonEnergy())
        print 'Fifo event codes: {}'.format(lcls.getEvrCodes())
        print 'FEE:', lcls.getPulseEnergy_mJ() 
    
    
