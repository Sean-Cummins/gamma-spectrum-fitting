"""
Main program used for processing and fitting spectral data.

(C) 2020 Sean Cummins, Drogheda, Ireland
Realeased under the GNU Public License (GPL)
email seancummins16@gmail.com
"""

#!/usr/bin/env python

import itertools
import re
import os
import logging
import logging.config
import csv
#import pprint

import numpy
import pandas
import lmfit
import matplotlib.pyplot as plt
import tabulate
import yaml
import asdf

import confParse

#Dynamic nested dictionaries
class NestDict(dict):
    def __getitem__(self,item):
        try:
            return dict.__getitem__(self,item)
        except:
            value = self[item] = type(self)()
            return value

"""
Main Class used in spectrum fitting. Loads, Fits and Saves Spectrum Data

Features:
-Initialised with a YAML config file and YAML log file
-Takes in .SPE or .MCA file. Extracts and stores spectrum, channel and other data.
-Builds fitting model using config file and spectrum 
-Fits spectrum data with model
-Saves spectrum fit parameters to CSV file.
-Saves plot of spectrum data with best fit to /plots

Other:
-Pretty printing of spectrum fit parameters
"""
class SpectrumFitter():

    def __init__(self):

        self.Spectra = NestDict()
        self.Fits = NestDict()

    def LoadConfigFile(self, configfile):

        self.config = confParse.YAML_Configurator(configfile)
        self.config.readfile()
        self.config.loadParams()

    def LoadLogConfigFile(self, configfile):

        with open(configfile, 'rt') as f:
            config = yaml.safe_load(f.read())
            logging.config.dictConfig(config)

        self.logger = logging.getLogger('dev')
        self.logger.setLevel(logging.INFO)

    def LoadSpectrum(self, filepath):

        File = os.path.basename(filepath)
        Check = File.split('_')

        if Check and len(Check) is 3:

            detector, source, angle_file = Check
            angle, extension = angle_file.split('.')

            if extension == "SPE" or extension == "Spe" or extension == "MCA" or extension == "mca":

                self.logger.info("File Extension: {extension} - Valid!".format(extension=extension))

                if detector in self.config.Detectors:

                    self.logger.info("Dectector Input: {detector} - Valid!".format(detector=detector))

                    if source in self.config.Detectors[detector].Sources:

                        self.logger.info("Source Input: {source} - Valid!".format(source=source))

                        if float(angle) in range(-180, 180):
                            self.logger.info("Angle Input: {angle} degs - Valid".format(angle = angle))
                        else:
                            self.logger.info("Angle Input: {angle} degs - Invalid!".format(angle = angle))   
                                      
                    else: 

                        self.logger.info("Source Input: {source} - Invalid!".format(source=source))
                        return
                else: 

                    self.logger.info("Dectector Input: {detector} - Invalid!".format(detector=detector))
                    return

            else:

                self.logger.info("File Extension: {extension} - Invalid!".format(extension=extension))
                return

        if extension == "SPE" or extension == "Spe":
        
            Channels, Spectrum, ChannelNumber, t_meas, t_elapsed = self.ReadSPEFile(filepath)
        
            self.Spectra[detector][source][angle]['Channels'] = Channels
            self.Spectra[detector][source][angle]['Spectrum'] = Spectrum
            self.Spectra[detector][source][angle]['LiveTime'] = t_meas
            self.Spectra[detector][source][angle]['Time'] = t_elapsed
            self.Spectra[detector][source][angle]['DeadTime'] = t_elapsed - t_meas

            self.config.Detectors[detector].Channels = ChannelNumber
        
        if extension == "MCA" or extension == "mca":

            Channels, Spectrum, ChannelNumber, t_meas, t_elapsed  = self.ReadMCAFile(filepath)
        
            self.Spectra[detector][source][angle]['Channels'] = Channels
            self.Spectra[detector][source][angle]['Spectrum'] = Spectrum
            self.Spectra[detector][source][angle]['LiveTime'] = t_meas
            self.Spectra[detector][source][angle]['Time'] = t_elapsed
            self.Spectra[detector][source][angle]['DeadTime'] = t_elapsed - t_meas

            self.config.Detectors[detector].Channels = ChannelNumber      

        return  

    def ReadSPEFile(self, filepath):

        Data = pandas.read_csv(filepath) # Read .SPE file into pandas dataframe                                                      
        Data = numpy.asarray(Data.iloc[0:-1]) # Lines as elements in numpy array
        Data = Data.flatten() # Collapse unnecessary dimension
        
        Spectrum = []             
                           
        flag = bool(False) # Set flag for extracting counts data only
    
        for i, counts_in_channel in enumerate(Data):   

            if (counts_in_channel=='$MEAS_TIM:'):

                t_meas, t_elapsed = numpy.asarray(Data[i+1].split(), dtype = float)

            if (counts_in_channel == '$ROI:'): # After this - no more channels/counts to collect so lock

                flag = bool(False)   

            elif (counts_in_channel == '$DATA:'): # After this - Channel Range then channels/counts we want to collect so unlock
                
                flag = bool(True)   
            
            elif (flag == bool(True)): # While unlocked collect channels/counts
                
                Spectrum.append(counts_in_channel)
    
        ChannelRange = Spectrum[0].split(' ')
        ChannelNumber = ChannelRange[1]

        Spectrum.pop(0)

        Spectrum = numpy.asarray(Spectrum)
        Spectrum = Spectrum.astype(float)#numbers still strings so cast to float    
        
        Channels = numpy.arange(0, Spectrum.shape[0], 1)

        return Channels, Spectrum, ChannelNumber, t_meas, t_elapsed 

    def ReadMCAFile(self, filepath):

        Data = pandas.read_csv(filepath, encoding= 'unicode_escape') # Read .MCA file into pandas dataframe, encoding deals with <<>>                                                     
        Data = numpy.asarray(Data.iloc[0:-1]) # Lines as elements in numpy array
        Data = Data.flatten() # Collapse unnecessary dimension
        
        Spectrum = []             
                           
        flag = bool(False) # Set flag for extracting counts data only
    
        for i, counts_in_channel in enumerate(Data):   

            if ('LIVE_TIME' in counts_in_channel): # Filter recorded live 

                time, value = Data[i].split(' - ')
                t_meas = float(value)

            if ('REAL_TIME' in counts_in_channel): # Filter recorded live and dead time

                time, value = Data[i].split(' - ')
                t_elapsed = float(value)

            if ('MCAC' in counts_in_channel): # #Channels MCAC=8192;    MCA/MCS Channels

                value, comment = Data[i].split(';')
                name, value = value.split('=')

                ChannelNumber = float(value)

            if (counts_in_channel == '<<END>>'): # After this - no more channels/counts to collect so lock

                flag = bool(False)   

            elif (counts_in_channel == '<<DATA>>'): # After this - Channel Range then channels/counts we want to collect so unlock
                
                flag = bool(True)   
            
            elif (flag == bool(True)): # While unlocked collect channels/counts
                
                Spectrum.append(counts_in_channel)

        Spectrum = numpy.asarray(Spectrum)
        Spectrum = Spectrum.astype(float)#numbers still strings so cast to float    
        
        Channels = numpy.arange(0, Spectrum.shape[0], 1)

        return Channels, Spectrum, ChannelNumber, t_meas, t_elapsed     

    def BuildSpectrumModel(self, channels, counts, **kwargs):
        
        if (kwargs['Model'] == 'Gaussian'):

            xmin, xmax = kwargs['ChannelRange']

            filteredchannels = channels[xmin:xmax]
            filteredcounts = counts[xmin:xmax] 

            Gauss1 = lmfit.models.GaussianModel(prefix='p1_') #p=Peak
            Line = lmfit.models.LinearModel(prefix='bg_') #bg=Background
            Model = Gauss1 + Line

            Pars = Line.make_params(c = numpy.median(filteredcounts))
            Pars += Gauss1.guess(filteredcounts, x = filteredchannels)

            Pars['p1_amplitude'].set(min=0)

        if (kwargs['Model'] == 'DoubleGaussian'):

            xmin1, xmax1 = kwargs['ChannelRange'][0]
            xmin2, xmax2 = kwargs['ChannelRange'][1]

            filteredchannels = channels[xmin1:xmax2]
            filteredcounts = counts[xmin1:xmax2]

            Gauss1 = lmfit.models.GaussianModel(prefix='p1_') 
            Gauss2 = lmfit.models.GaussianModel(prefix='p2_') 
            Line = lmfit.models.LinearModel(prefix='bg_')
            Model = Gauss1 + Gauss2 + Line

            Pars = Line.make_params(c = numpy.median(filteredcounts))
            Pars += Gauss1.guess(counts[xmin1:xmax1], x = channels[xmin1:xmax1])
            Pars += Gauss2.guess(counts[xmin2:xmax2], x = channels[xmin2:xmax2])
        
            Pars['p1_amplitude'].set(min=0)
            Pars['p2_amplitude'].set(min=0)    

        return filteredchannels, filteredcounts, Model, Pars


    def BuildBackgroundModel(self, counts, **kwargs):

            BackgroundModel = lmfit.models.LinearModel(prefix='bg_') #bg=Background
            BackgroundPars = BackgroundModel.make_params(c = numpy.median(counts))

            return BackgroundModel, BackgroundPars

    def BuildPeakModel(self, channels, counts, peakname, **kwargs):

        if (kwargs['Model'] == 'Gaussian'):
            PeakModel = lmfit.models.GaussianModel(prefix='p{peak}_'.format(peak = peakname)) #p=Peak
            PeakPars = PeakModel.guess(counts, x = channels)
            PeakPars['p{peak}_amplitude'.format(peak=peakname)].set(min=0)            

        return  PeakModel, PeakPars

    def FitSpectrum(self, filepath, savefig=True):

        File = os.path.basename(filepath)
        CheckFileName = File.split('_')
        
        if CheckFileName and len(CheckFileName) is 3:

            Detector, Source, AngleExtension = CheckFileName
            Angle, Extension = AngleExtension.split('.')

            try: 
                
                CheckLoaded = self.Spectra[Detector][Source][Angle] # is spectrum loaded?

                self.logger.info("Spectrum: {spectrum} Loaded".format(spectrum=File))

                Channels = self.Spectra[Detector][Source][Angle]['Channels']
                Spectrum = self.Spectra[Detector][Source][Angle]['Spectrum']

                Peaks = self.config.Detectors[Detector].Sources[Source]['Peaks']

                for Peak in Peaks:

                    self.logger.info("Fitting Peak {peak}".format(peak=Peak))
                    
                    FittingSpecs = Peaks[Peak]['Fitting']

                    xmin, xmax = FittingSpecs['ChannelRange']
                    FilteredChannels = Channels[xmin:xmax]
                    FilteredCounts = Spectrum[xmin:xmax] 

                    BackgroundModel, BackgroundPars = self.BuildBackgroundModel(FilteredCounts)
                    PeakModel, PeakPars = self.BuildPeakModel(FilteredChannels, FilteredCounts, Peak, **FittingSpecs)

                    Model = BackgroundModel + PeakModel
                    Pars = BackgroundPars + PeakPars

                    if 'FitWith' in FittingSpecs.keys(): # if FitWith parameter enabled

                        self.logger.info("Parameter: FitWith ENABLED")

                        if isinstance(FittingSpecs['FitWith'], list): # if user has written parameter in yaml file correctly

                            self.logger.info("Parameter: FitWith TYPE VALID")
                            
                            check = FittingSpecs['FitWith'][0]

                            # check if peak has already been fitted by another in the composite (no need to repeat)
                            if check not in self.Fits[Detector][Source][Angle].keys():

                                self.logger.info("Parameter: FitWith FITTING COMPOSITE MODEL")
                                                                
                                for FitWithPeak in FittingSpecs['FitWith']:

                                    FitWithPeak_FittingSpecs = Peaks[FitWithPeak]['Fitting']
                                    
                                    xmin_new, xmax_new = Peaks[FitWithPeak]['Fitting']['ChannelRange']
                                    
                                    FitWithPeak_FilteredChannels = Channels[xmin_new:xmax_new]
                                    FitWithPeak_FilteredCounts = Spectrum[xmin_new:xmax_new] 
                                    
                                    FitWith_PeakModel, FitWith_PeakPars = self.BuildPeakModel(FitWithPeak_FilteredChannels,\
                                                                                            FitWithPeak_FilteredCounts,\
                                                                                            FitWithPeak, **FitWithPeak_FittingSpecs)

                                    Model += FitWith_PeakModel
                                    Pars += FitWith_PeakPars
                       
                                    #Get xmin and xmax - WARNING: Only works if in order                                                     
                                    FilteredChannels = Channels[xmin:xmax_new]
                                    FilteredCounts = Spectrum[xmin:xmax_new]

                            else:
                                self.logger.info("Peak: {check} already fitted. Skipping over".format(check=check))
                                #we've already fitted this peak when we were fitting it in composite with another
                                return
                        else: 
                            self.logger.info("FittingSpec DataType is: {DataType}, Should be: list".format(Datatype = type(FittingSpecs['FitWith'])))
                            #user has written in yaml file correctly
                            return
                    
                    Result = Model.fit(FilteredCounts, Pars, x = FilteredChannels, method = FittingSpecs['Method'])

                    if savefig:
                    
                        if not os.path.exists('plots'):
                            os.mkdir('plots')
                        if not os.path.exists('plots/{detector}'.format(detector=Detector)):
                            os.mkdir('plots/{detector}'.format(detector=Detector))

                        plt.figure()
                        plt.plot(FilteredChannels, FilteredCounts, 'b+' )

                        Components = Result.eval_components()

                        for Prefix, Component in Components.items():
                            plt.plot(FilteredChannels, Component, '--', label = Prefix)

                        plt.plot(FilteredChannels, Result.best_fit,'r-')
                        plt.xlabel('Channels')
                        plt.ylabel('Counts')
                        plt.legend()
                        plt.savefig('plots/{detector}/{detector}_{source}_{angle}_{peak}.png'.format(detector = Detector,\
                                    source = Source, angle = Angle, peak = Peak))
                        plt.close()

                    #How to store fitted parameters - Need better way
                    # For now analysis is easier if time is stored in fits
                    for Name, Param in Result.params.items():
                        prefix, parameter = Name.split('_')
                        try:# peak pars
                            prefix_search = re.search(r"\d", prefix)
                            peak = prefix_search.group()
                            self.Fits[Detector][Source][Angle][peak][parameter] = Param.value 
                            self.Fits[Detector][Source][Angle][peak]['Time'] = self.Spectra[Detector][Source][Angle]['LiveTime']
                        except:#background pars
                            self.Fits[Detector][Source][Angle][Peak][Name] = Param.value
                                                        

            except KeyError:
                self.logger.info("Spectrum: {spectrum} Not Loaded".format(spectrum=File))

    def SaveFitsCSV(self, detector, path=None):

        filename = "{detector}_fits.csv".format(detector = detector)
        # User specified path
        if path is not None:
            filepath = os.path.join(path,filename)
        else:
            filepath = filename
        #Delete file if already there - stops appending
        try:
            os.remove(filepath)
        except OSError:
            pass
        #https://stackoverflow.com/questions/13575090/construct-pandas-dataframe-from-items-in-nested-dictionary
        df = pandas.DataFrame.from_dict({(i,j,k): self.Fits[detector][i][j][k] \
                                            for i in self.Fits[detector].keys() \
                                            for j in self.Fits[detector][i].keys() \
                                            for k in self.Fits[detector][i][j].keys()}, orient='index')
        df.to_csv(filepath, mode='a', index=True, header=True)

    #WARNING - jsonschema version incompatibility - Fix needed
    def SaveASDF(self, Data, filename):
        af = asdf.AsdfFile(Data)
        af.write_to('{filename}.asdf'.format(filename=filename))
        self.logger.info("Data written to {filename}".format(filename=filename))

    def PrintTable(self, detector, source):

        rowheader = self.Fits[detector][source]['1'].keys()
        ordered_rowheader = sorted(rowheader, key=lambda x: float(x))

        colheader = self.Fits[detector][source]['1']['0'].keys()

        for Peak in self.Fits[detector][source]:

            rows = [self.Fits[detector][source][Peak][angle].values() for angle in ordered_rowheader]   

            print("DETECTOR: {detector} | SOURCE: {source} | PEAK: {peak}".format(detector=detector, source=source, peak=Peak))
            print(tabulate.tabulate(rows, colheader, tablefmt='presto', showindex=ordered_rowheader))

    def PlotSpectrum(self, detector, source, angle):
        
        plt.figure()
        plt.plot(self.Spectra[detector][source][angle]['Channels'], self.Spectra[detector][source][angle]['Spectrum'])
        plt.show()


#Get file paths
def SearchPath(root):
    kw_files = []
    kws = ['NaI','CdTe', 'BGO', 'HPGe']
    #kws = ['HPGe_Co']
    print('Searching in Dataset for files with keywords :', kws)
    for path, subdirs, files in os.walk(root):
        for name in files:
            if any(kw in name for kw in kws):
                print('Found :', name)
                kw_files.append(os.path.join(path, name)) 
            else:
                print('Found : None')

    return kw_files

if __name__ == "__main__":

    Fitter = SpectrumFitter()
    Fitter.LoadLogConfigFile("logconf.yml")
    Fitter.LoadConfigFile("conf.yml")
    
    filepaths = SearchPath('Data')

    for filepath in filepaths:
    
        Fitter.LoadSpectrum(filepath)
        Fitter.FitSpectrum(filepath)

    for detector, sources in Fitter.config.Detectors.items():
    #    for source in sources.Sources.keys():
    #        Fitter.PrintTable(detector, source)

        Fitter.SaveFitsCSV(detector)









