try:
    import psana
except:
    print 'Using the module "lcls" without psana capabilities.'
    psana = None
import numpy as np

_EBeamType = None
_feeType = None
_evt = None
_currentFiducial = None

if psana is not None:
    EBeamTypeList = (
            psana.Bld.BldDataEBeamV0,
            psana.Bld.BldDataEBeamV1,
            psana.Bld.BldDataEBeamV2,
            psana.Bld.BldDataEBeamV3,
            psana.Bld.BldDataEBeamV4,
            psana.Bld.BldDataEBeamV5,
            psana.Bld.BldDataEBeamV6,
            psana.Bld.BldDataEBeamV7)
            
    _EBeamSource = psana.Source('BldInfo(EBeam)')

    feeTypeList = (psana.Bld.BldDataFEEGasDetEnergy,
                   psana.Bld.BldDataFEEGasDetEnergyV1)

    _feeSource = psana.Source('BldInfo(FEEGasDetEnergy)')

    _feeFuncs = ['f_11_ENRC', 'f_12_ENRC',
                 'f_21_ENRC', 'f_22_ENRC',
                 'f_63_ENRC', 'f_63_ENRC']

    _evrSource = psana.Source('DetInfo(NoDetector.0:Evr.0)')
    _evrType = psana.EvrData.DataV3
    evrTypeList = (psana.EvrData.DataV3)
else:
    EBeamTypeList = None
    _EBeamSource = None
    feeTypeList = None
    _feeSoutce = None
    _evrSource = None
    _evrType = None
    evrTypeList = None



def setEvent(evt, verbose=False):
    global _evt
    global _currentFiducial

    if verbose:
        print 'Trying to set event.'

    # If no event is given
    if evt is None:
        # And ther is no event set already
        if _evt is None:
            if verbose:
                print 'No event set.'
            return False
        #Use old event if there is one
        else:
            if verbose:
                print 'No new event given, using old event.'
            return True

    # If a new event is given check if it is the same as the old one using the
    # fiducial
    newFid = evt.get(psana.EventId).fiducials()
    # If it is not the same as the old one
    if newFid != _currentFiducial:
        if verbose:
            print 'Updating the internal event.'
        # Update the infomation
        _evt = evt
        _currentFiducial = newFid
        return True
    elif verbose:
        print 'Event already up to date.'

    return True

# e-beam data
def getEBeamEnergyL3_MeV(evt=None, verbose=False):
    EBeamObject = getEBeamObject(evt=evt, verbose=verbose)
    if EBeamObject == None:
        return np.nan
    return EBeamObject.ebeamL3Energy()

def getEBeamEnergyBC2_MeV(evt=None, verbose=False):
    if not setEvent(evt):
        return np.nan
    EBeamObject = getEBeamObject(verbose=verbose)
    if EBeamObject == None:
        return np.nan
    # The BC2 energy is calculated using the dispersion. The dispersion
    # value for BC2 is -0.3647 m at the nominal beam energy of 5 GeV [Email
    # from Timothy Maxwell to Anton Lindahl on June 2 2014]
    return (getEBeamPosOffsetBC2_mm / -364.7 + 1) * 5e3

def getEBeamPosOffsetBC2_mm(evt=None, verbose=False):
    EBeamObject = getEBeamObject(verbose=verbose)
    if EBeamObject == None:
        return np.nan
    return EBeamObject.ebeamEnergyBC2()


def getEBeamCharge_nC(evt=None, verbose=False):
    if not setEvent(evt):
        return np.nan
    EBeamObject = getEBeamObject(verbose=verbose)
    if EBeamObject == None:
        return np.nan
    return EBeamObject.ebeamCharge()

def getEBeamPkCurrentBC2_A(evt=None, verbose=False):
    if not setEvent(evt):
        return np.nan
    EBeamObject = getEBeamObject(verbose=verbose)
    if EBeamObject == None:
        return np.nan
    return EBeamObject.ebeamPkCurrBC2()

# e-beam setup
def getEBeamObject(evt=None, verbose=False):
    if not setEvent(evt, verbose):
        return None
    # Initialize the EBeam type
    if _EBeamType is None:
        if verbose:
            print 'Initializing the EBeam type.'
        _determineEBeamType(verbose=verbose)
    try:
        return _evt.get(_EBeamType, _EBeamSource)
    except:
        return None

def _determineEBeamType(evt=None, verbose=False):
    global _EBeamType
    if not setEvent(evt, verbose):
        return None
    if verbose:
        print 'Find the correct EBeam type.'
    for type in EBeamTypeList:
        if verbose:
            print 'Trying {};'.format(type),
        data = _evt.get(type, _EBeamSource)
        if data is not None:
            _EBeamType = type
            if verbose:
                print ' correct.'
            break
        elif verbose:
            print ' wrong one.'

    return _EBeamType

# fee
def getPulseEnergy_mJ(evt=None, nValues=4, verbose=False):
    if not setEvent(evt):
        return np.nan
    fee = getFeeObject(_evt, verbose)
    if fee is None:
        return np.array([np.nan for i in range(nValues)])
    return np.array([getattr(fee, func)() if hasattr(fee, func) else 0.0
                     for func in _feeFuncs[:nValues]])

def getFeeObject(evt, verbose=False):
    if _feeType is None:
        _determineFeeType(evt, verbose)
    try:
        return evt.get(_feeType, _feeSource)
    except:
        return None

def _determineFeeType(evt, verbose=False):
    global _feeType
    if verbose:
        print 'Find the correct fee type.'
    for type in feeTypeList:
        if verbose:
            print 'Trying {};'.format(type),
        data = evt.get(type, _feeSource)
        if data is not None:
            _feeType = type
            if verbose:
                print ' correct.'
            break
        elif verbose:
            print ' wrong one.'


# Event id
def getIdObject(evt, verbose=False):
    return evt.get(psana.EventId)

def getEventFiducial(evt=None, verbose=False):
    if not setEvent(evt, verbose):
        return 0
    id = getIdObject(_evt, verbose)
    if id == None:
        return 0
    return id.fiducials()

def getEventTime(evt=None, offset=0, verbose=False):
    if not setEvent(evt, verbose):
        return np.nan
    id = getIdObject(_evt, verbose)
    if id == None:
        return np.nan
    time =  id.time()
    return (time[0] - offset) + time[1] * 1e-9


##############################
# EVR functionality

def getEvrObject(evt=None, evrSource=_evrSource, verbose=False):
    if verbose:
        print 'Ger the evr object.'
    if not setEvent(evt, verbose):
        if verbose:
            print 'No velid event.'
        return None
    if _evrType is None:
        _determineEvrTypa(_evt, verbose)
    return _evt.get(_evrType, evrSource)

def _determineEvrType(evt, verbose=False):
    global _evrType
    if verbose:
        print 'Finding the correct EVR type.'
    for type in evrTypeList:
        if verbose:
            print 'Trying {};'.format(type),
        data = evt.get(type, _evrSource)
        if data is not None:
            _evrType = type
            if verbose:
                print ' correct.'
                break
        elif verbose:
            print ' wrong one.'


def getEvrCodes(evt=None, verbose=False):
    evrData = getEvrObject(evt, verbose=verbose)
    if evrData is None:
        if verbose:
            print 'No evr object returned.'
        return []
    return [fifo.eventCode() for fifo in evrData.fifoEvents()]

def evrCodeInEvent(evr, evt=None, verbose=False):
    if verbse:
        print 'Checking for EVR codes.'
    return evr in getEvrCodes(evt, verbose=verbose)


if __name__ == '__main__':
    print 'Connecting to data source.'
    ds = psana.DataSource('exp=amoi0314:run=14')
    for i in range(10):
        print 'Geting an event.'
        evt = ds.events().next()
        print 'Set the event'
        setEvent(evt, verbose=True)
        print 'E at L3 is {} MeV'.format(getEBeamEnergyL3_MeV(verbose=True))
        print 'E at BC2 is {} MeV'.format(getEBeamEnergyBC2_MeV(verbose=True))
        print 'Q is {} nC'.format(getEBeamCharge_nC(verbose=True))
        print 'I is {} A'.format(getEBeamPkCurrentBC2_A(verbose=True))
        print 'fee is {} mJ'.format(getPulseEnergy_mJ(verbose=True))
        print 'Evr codes in event: {}'.format(getEvrCodes(evt, verbose=True))
