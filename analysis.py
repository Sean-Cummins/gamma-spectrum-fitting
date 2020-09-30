"""
Program used for analysing fitted spectral data.

(C) 2020 Sean Cummins, Drogheda, Ireland
Realeased under the GNU Public License (GPL)
email seancummins16@gmail.com
"""

import os
import logging
import logging.config
import datetime

import numpy
import pandas
import yaml
import matplotlib.pyplot as plt
import scipy
import tabulate

from main import NestDict, SpectrumFitter
import confParse

class DetectorAnalyser():

    def __init__(self):

        self.Calibration = NestDict()
        self.Resolution = NestDict()
        self.PeakEfficiency = NestDict()

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

    def LoadFitCSV(self, filepath):

        self.Fits = pandas.read_csv(filepath, index_col=[0,1,2]) 

    def CalibrateDetector(self, detector, savefig = True):
        
        Indices = list(self.Fits.index.values) 
        OnAxisIndices = [(Source, Angle, Peak) for (Source, Angle, Peak) in Indices if Angle is 0]

        PeakEnergies = numpy.zeros(len(OnAxisIndices))
        PeakChannels = numpy.zeros(len(OnAxisIndices))

        for i, index in enumerate(OnAxisIndices):

            Source, Angle, Peak = index
            PeakEnergy = self.config.Detectors[detector].Sources[Source]['Peaks'][str(Peak)]['Energy']
            PeakChannel = self.Fits.loc[index, 'center']
            
            PeakEnergies[i] = PeakEnergy
            PeakChannels[i] = PeakChannel

        #Now do regression
        slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(PeakChannels, PeakEnergies)

        if savefig:

            if not os.path.exists('analysis'):
                os.mkdir('analysis')
            if not os.path.exists('analysis/plots'):
                os.mkdir('analysis/plots')

            plt.figure()
            plt.plot(PeakChannels, PeakEnergies, '+')
            plt.plot(PeakChannels, intercept + slope * PeakChannels, '.')
            plt.xlabel('Channel #')
            plt.ylabel('Energy [keV]')
            plt.savefig('analysis/plots/{detector}_calibration.png'.format(detector = detector))
            plt.close()

    def OnAxisResolution(self, detector, savefig = True):

        Indices = list(self.Fits.index.values) # get hashable tuple list for dataframe
        OnAxisIndices = [(Source, Angle, Peak) for (Source, Angle, Peak) in Indices if Angle is 0]

        PeakEnergies = numpy.zeros(len(OnAxisIndices))
        Resolution = numpy.zeros(len(OnAxisIndices))

        for i, index in enumerate(OnAxisIndices):

            Source, Angle, Peak = index
            
            PeakEnergy = self.config.Detectors[detector].Sources[Source]['Peaks'][str(Peak)]['Energy']
            PeakFWHM = self.Fits.loc[index, 'fwhm'] 

            PeakEnergies[i] = PeakEnergy
            Resolution[i] = PeakFWHM / PeakEnergy

        if savefig:

            if not os.path.exists('analysis'):
                os.mkdir('analysis')
            if not os.path.exists('analysis/plots'):
                os.mkdir('analysis/plots')

            plt.figure()
            plt.plot(PeakEnergies, Resolution, '.') 
            plt.xlabel('Energy (keV)')
            plt.ylabel('Resolution')
            plt.savefig('analysis/plots/{detector}_on_axis_resolution.png'.format(detector = detector))
            plt.close()

    def DetectorEfficiency(self, Detector, savefig = True, printtable = True):
        
        Indices = list(self.Fits.index.values) # get hashable list of tuples for dataframe

        for index in Indices:

            Source, Angle, Peak = index

            Geometry = self.config.Detectors[Detector].Geometry
            CalibrationActivity = self.config.Detectors[Detector].Sources[Source]['CalibrationActivity']
            CalibrationTime = self.config.Detectors[Detector].Sources[Source]['CalibrationTime']
            HalfLife = self.config.Detectors[Detector].Sources[Source]['HalfLife']
            SourceDetectorSeparation = self.config.Detectors[Detector].Sources[Source]['SourceDetectorSeparation']
            DecayFraction = self.config.Detectors[Detector].Sources[Source]['Peaks'][str(Peak)]['DecayFraction']
  
            Counts = self.Fits.loc[index, 'height']
            Time = self.Fits.loc[index, 'Time']

            PeakEnergy = self.config.Detectors[Detector].Sources[Source]['Peaks'][str(Peak)]['Energy']

            Activity = self.ActivityNow(CalibrationActivity, CalibrationTime, HalfLife)
            Flux, Events = self.CalculateEvents(Activity, Time, DecayFraction)
            GeometricFactor = self.CalculateGeometricFactor(SourceDetectorSeparation, Angle, Geometry)

            Efficiency = self.CalculatePeakEfficiency(Counts, Time, Activity, DecayFraction, GeometricFactor)

            self.PeakEfficiency[Detector][PeakEnergy][Angle]['Counts'] = Counts
            self.PeakEfficiency[Detector][PeakEnergy][Angle]['Time'] = Time
            self.PeakEfficiency[Detector][PeakEnergy][Angle]['CountRate'] = Counts / Time
            self.PeakEfficiency[Detector][PeakEnergy][Angle]['Activity'] = Activity
            self.PeakEfficiency[Detector][PeakEnergy][Angle]['DecayFraction'] = DecayFraction
            self.PeakEfficiency[Detector][PeakEnergy][Angle]['GeometricFactor'] = GeometricFactor  
            self.PeakEfficiency[Detector][PeakEnergy][Angle]['Efficiency'] = Efficiency                                                                   

        if savefig:

            #plt.figure()
            
            for PeakEnergy, Value in self.PeakEfficiency[Detector].items():

                Angle = []
                Efficiencies = []

                for key, value in Value.items():

                    Angle.append(key)
                    Efficiencies.append(value['Efficiency'])

                    #print(key, value)
        
                Angle, Efficiencies = (list(t) for t in zip(*sorted(zip(Angle, Efficiencies)))) # sort before line plotting
                
                plt.figure()
                plt.plot(Angle, Efficiencies, '.-', label = (str(PeakEnergy) + str(' keV')))
                #plt.polar(numpy.deg2rad(Angle), Efficiencies, '.-', label = (str(PeakEnergy) + str(' keV')))

                plt.xlabel('Angle (degs)')
                plt.ylabel('Intrinsic Efficiency')
                #plt.yscale('log')
                plt.grid()
                plt.legend()
            plt.show()

        if printtable:
            for peakenergy, efficiencyvalues in self.PeakEfficiency[Detector].items():
                self.PrintTable(Detector, peakenergy)


    def ActivityNow(self, CalibrationActivity, CalibrationTime, HalfLife, Now = '2018-11-1-12'):

        HalfLife = HalfLife * 365.25 * 24 * 60 * 60 #seconds

        t1 = datetime.datetime.strptime(CalibrationTime, '%Y-%m-%d-%H') 
        t2 = datetime.datetime.strptime(Now, '%Y-%m-%d-%H')
        dt = (t2-t1).total_seconds()
    
        DecayConstant = numpy.log(2) / HalfLife

        Activity = CalibrationActivity * numpy.exp(-DecayConstant * dt)

        return Activity

    #-----------------------------------------------------------------------------#
    """Function to calculate events and fluences from known live times, activites and
    decay fractions"""
    def CalculateEvents(self, Activity, Time, DecayFraction):
        
        Flux = Activity * 1e-6 * 3.7e10 * DecayFraction
        Events = Flux * Time

        return Flux, Events

    #-----------------------------------------------------------------------------#
    """Function to calculate Geometric Factor from cartesian co-ordinates and
    known detector geometries"""
    def CalculateGeometricFactor(self, SourceDetectorSeparation, Angle, Geometry):

        Angle = numpy.deg2rad(Angle)
        
        if (Geometry['Shape'] == 'Cylinder'):

            Diameter = Geometry['Dimensions']['Diameter'] * 0.01 # cm to m
            Height = Geometry['Dimensions']['Height'] * 0.01
            Radius = Diameter / 2
            
            Face = numpy.pi * Radius * Radius
            Side = Diameter * Height

        if (Geometry['Shape'] == 'Cuboid'):

            Length = Geometry['Dimensions']['Length'] * 0.01
            Width = Geometry['Dimensions']['Width'] * 0.01
            Depth = Geometry['Dimensions']['Depth'] * 0.01
            
            Face = Length * Width
            Side = Length * Depth      
        
        Area = Face * pow(numpy.cos(Angle), 2) + Side * pow(numpy.sin(Angle), 2)

        GeometricFactor = Area / (4 * numpy.pi * pow(SourceDetectorSeparation, 2))

        return GeometricFactor
           
    #-----------------------------------------------------------------------------#
    """Function to calculate peak effiency"""
    def CalculatePeakEfficiency(self, Counts, Time, Activity, DecayFraction, GeometricFactor):

        Activity = Activity * 1e-6 * 3.7e10

        CountRate = Counts/Time

        PeakEfficiency = CountRate * (1 / (Activity * DecayFraction * GeometricFactor))

        return PeakEfficiency

    def unnest(self, d, keys=[]):
        result = []
        for k, v in d.items():
            if isinstance(v, dict):
                result.extend(self.unnest(v, keys + [k]))
            else:
                result.append(tuple(keys + [k, v]))
        return result

    def PrintTable(self, detector, peakenergy):

        rowheader = self.PeakEfficiency[detector][peakenergy].keys()
        ordered_rowheader = sorted(rowheader, key=lambda x: float(x))

        colheader = self.PeakEfficiency[detector][peakenergy][0].keys()

        rows = [self.PeakEfficiency[detector][peakenergy][angle].values() for angle in ordered_rowheader]   

        print("DETECTOR: {detector} | PEAK: {peak}".format(detector=detector, peak=peakenergy))
        print(tabulate.tabulate(rows, colheader, tablefmt='presto', showindex=ordered_rowheader))

if __name__ == "__main__":

    Analyser = DetectorAnalyser()

    Analyser.LoadConfigFile("conf.yml")
    Analyser.LoadLogConfigFile("logconf.yml")

    Detector = 'HPGe'
    
    Analyser.LoadFitCSV('{detector}_fits.csv'.format(detector=Detector))
    Analyser.CalibrateDetector(Detector)
    Analyser.OnAxisResolution(Detector)

    PeakEfficiency = Analyser.DetectorEfficiency(Detector)
