from configuration import loadConfiguration


def prepConfig(fileName):
    return loadConfiguration(fileName)





if __name__ == '__main__':
    config = prepConfig('cookieBoxDefaultConfig.json')
