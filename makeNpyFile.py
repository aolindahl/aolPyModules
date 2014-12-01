import numpy as np


def dataName(txtName, ext=False):
    if txtName[-4:] != '.txt':
        return None
    return txtName[:-4] + '_data.npy'

def indexName(txtName, ext=False):
    if txtName[-4:] != '.txt':
        return None
    return txtName[:-4] + '_colnamesDict.npy'
    
def npyFromTxt(fileName):
    try:
        with open(fileName) as file:
            cols = file.readline().split()
            data = np.loadtxt(file)
    except:
        print 'Failure reading the file "{}".'.format(fileName)
        raise
        sys.exit()
    
    ind = {}
    for i, k in enumerate(cols):
        ind[k] = i
    
    np.save(indexName(fileName), ind)
    np.save(dataName(fileName), data)


def loadFromNpy(fileName):
    try:
        return [np.load(indexName(fileName)).tolist(), np.load(dataName(fileName))]
    except:
        print 'Could not lode npy files. Trying to make them from "{}".'.format(
                fileName)
        npyFromTxt(fileName)
        return [np.load(indexName(fileName)).tolist(),
                np.load(dataName(fileName))]

if __name__ == '__main__':
    import sys

    args = sys.argv

    if len(args) < 2:
        print 'Needs a file.'
        sys.exit()

    loadFromNpy(args[1])
