"""
YAML Parser. Takes in .yml config file and creates a config object used
in other programs. 

(C) 2020 Sean Cummins, Drogheda, Ireland
Realeased under the GNU Public License (GPL)
email seancummins16@gmail.com
"""

import yaml

class ConfigurationError(Exception):
    pass

class YAML_Configurator():

    def __init__(self, filename):

        self.filename = filename
        self.Detectors = {}

    def readfile(self):

        # load config file from Yaml file
        with open(self.filename) as file_:
            self.configDict = yaml.safe_load(file_)

    def loadParams(self):

        self.readfile()

        for DetectorName in self.configDict["Detectors"]:

            Detector = DetectorConfig()
            Detector.loadParams(self.configDict["Detectors"][DetectorName])

            self.Detectors[DetectorName] = Detector


# Parameters that can be had by any configuration object
class ConfigObj(object):

    def warnAndExit(self, param):
        message = "{0} not set!".format(param)
        raise ConfigurationError(message)

    def warnAndDefault(self, param, newValue):
        message = "{0} not set, default to {1}".format(param, newValue)
        self.__setattr__(param, newValue)
 
    def initParams(self):
        for param in self.requiredParams:
            self.__setattr__(param, None)

    def loadParams(self, configDict):

        for param in self.requiredParams:
            try:
                self.__setattr__(param, configDict[param])
            except KeyError:
                self.warnAndExit(param)
            except IndexError:
                raise ConfigurationError("Not enough values for {0}".format(param))
            except:
                raise ConfigurationError("Failed while loading {0}. Check config file.".format(param))

        for param in self.optionalParams:
            try:
                self.__setattr__(param[0], configDict[param[0]])
            except KeyError:
                self.warnAndDefault(param[0], param[1])
            except IndexError:
                raise ConfigurationError("Not enough values for {0}".format(param))
            except:
                raise ConfigurationError("Failed while loading {0}. Check config file.".format(param))


        self.calcParams()

    def calcParams(self):
        """
        Dummy method to be overidden if required
        """
        pass

    def __iter__(self):
        for param in self.requiredParams:
            yield param, self.__dict__[param]
        for param in self.optionalParams:
            yield param[0], self.__dict__[param[0]]

    def __len__(self):
        return len(self.requiredParams)+len(self.optionalParams)

    def __setattr__(self, name, value):
        if name in self.allowedAttrs:
            self.__dict__[name] = value
        else:
            raise ConfigurationError("'{}' Attribute not a configuration parameter".format(name))

    def __repr__(self):
        return str(dict(self))


class DetectorConfig(ConfigObj):
    """
    Configuration parameters relavent for the Detectors.

    Required:
        =============                  ===================
        **Parameter**                  **Description**
        -------------                  -------------------
        ``Sources``
        ``Geometry`
        =============                  ===================

    Optional:
        =============                  ===================
        **Parameter**                  **Description**
        -------------                  -------------------
        ``Resolution``
        ``Efficiency``
        ``Channels``
        ``FWHM``
        ``Area``
        =============                  ===================

    """
    requiredParams = ['Sources', 'Geometry']
       
    optionalParams = [('Resolution', 1),
                      ('Efficiency', 1),
                      ('Channels', 1),
                      ('FWHM', 1),
                      ('Area',1)]

    calculatedParams = []

    allowedAttrs = requiredParams + calculatedParams 

    for p in optionalParams:
        allowedAttrs.append(p[0])