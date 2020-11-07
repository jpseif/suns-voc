import numpy as np
import matplotlib.pyplot as plt
import os
from semiconductor.recombination.intrinsic import Radiative
from semiconductor.material.bandgap_narrowing import BandGapNarrowing as BGN
from semiconductor.material.intrinsic_carrier_density import IntrinsicCarrierDensity as ICD
from semiconductor.electrical.mobility import Mobility
import scipy.constants as C
from scipy.optimize import minimize_scalar
from scipy.stats import linregress
import re
import platform
import pandas as pd

# Last review: 18.03.2019 by J.Seif
# a = [0.633841, 0.632358, 0.630932, 0.629482, 0.628043]
# suns = [1.062788474, 1.008779061, 0.953958622, 0.905329773, 0.858237162]
# power = np.array([0.010401910, 0.010407574, 0.010406887, 0.010408312, 0.010407651, 0.010405968, 0.010402238])
class Suns_Voc_Measurement:
    """
    This class contains a collection of functions to calculate different cell
    parameters from Suns-Voc measurements. Among them are the following:
    Voc, pFF, J01 and J02 [not implemented yet], as well as the
    local ideality factor [m(V)].

    The analysis is based on the code of Ron Sinton (Suns-Voc spreadsheet v4.4)
    and the Semiconductor package of Mattias Juhl.

    Inputs:
        1. Sample file (*.smpl) containing the following information:
            - Wafer thickness [cm]
            - Base resistivity [Ohm.cm]
            - Base type [n-type or p-type]
            - Calibration constant [V/suns]
            - Short circuit current density (Jsc) [A/cm2]
        2. Suns-Voc(T) or Suns-PL(T) data for the custom software, developed by
            Eric Okawa from BTi (eric.okawa@btimaging.com).

    Output:
        1. TXT file with all the calculated values, similar to the Sinton spreadsheet.
    """

    # def __init__(self, f, Na=1e16, Nd=1, W=0.02, binn=1, bkgnd_corr=0.95, bkgnd_loc='end', T=300, crop_start=0, crop_end=1, sinton_consts=(1,2,3)):
    def __init__(
        self,
        sp,
        f,
        fn,
        findex=0,
        jsc=3.7E-2,
        wThick=0.0180,
        rBase=2.8,
        idealityPts=4,
        BGNon=1,
        TdepON=1,
        Na=0,
        Nd=0,
        T=25.0,
        VocLO=0,
        VocOneSunLO=0,
        SunsLO=0,
        SunsDataLO=0,
        illumSetpoints=np.array([]),
        dType='n-type',
        refCal=0.115273,
        aMode='QSS'
        ):
        # file-related
        self.sp = sp # path to all the measurement files
        self.f = f # file path
        self.fn = fn # file name
        self.sn = re.findall('\[([^]]+)\]', self.fn)[1] # sample name
        self.findex = findex # this is the index of the file (the 3rd file that is imported), used for jsc
        # measurement-related
        self.refCal = refCal # reference calibration constant [V/suns]
        self.VocLO = VocLO
        self.VocOneSunLO = VocOneSunLO
        self.SunsLO = SunsLO
        self.indClosestVoc = 4
        self.HIShifted = 0
        self.names = ['t', 'ref_volt', 'suns', 'volt'] # for the columns
        self.mpp = 0
        self.idealityPts = idealityPts
        # exterior parameters
        self.T = T # sample temperature [°C]
        if TdepON == 1:
            self.TK = self.T + 273.15 # sample temperature [K]
        else:
            self.TK = 300
        # constants
        self.const_q = 1.602176487E-19 # elementary charge [C]
        self.const_k = 1.38064504E-23 # boltzmann constant [J/K]
        self.const_pi = 3.1415926535898 # Pi
        self.const_h = 6.62606896E-34 # planck constant [Js]
        self.const_e_Si = 11.7 # dielectric constant of Si (ref. Kittel)
        # analysis related
        self.BGNon = BGNon
        # cell parameters
        self.jsc = jsc # Jsc from I-V [A/cm2], e.g. 3.7E-2
        self.wThick = wThick # wafer thickness [cm], e.g. 0.0180
        self.rBase = rBase # base resistivity [Ohm.cm]
        self.dType = dType # doping type
        if self.rBase == 0:
            self.Na = Na
            self.Nd = Nd
            self.doping = np.abs(self.Na - self.Nd)
        else:
            self.Na = 0
            self.Nd = 0
            self.doping = self.calc_doping() # doping
        # self.tempCo = tempCo # Fs start value, necessary to calculate Vcorrected
        # might be necessary to implement again to correct for the Temperature
        # difference between the temperature at which the cell was measured and
        # 25°C (as Sinton does it): Vcorr = V - (self.T - 25)*tempCo
        # analysis-related
        self.aMode = aMode # analysis mode, e.g. 'QSS' or 'GEN'
        # data-related
        self.load_data()
        # calculated Values
        self.ni_initial = self.calc_ni_initial()
        self.dVdt = self.calc_dVdt()
        self.Dn = self.calc_Dn(ni_init = 1)

        if self.BGNon == 1:
            self.ni_eff = self.calc_ni_eff() # formerly calc_niBGN
            self.Dn = self.calc_Dn(ni_init = 0) # first iteration after ni_eff
        else:
            self.ni_eff = self.ni_initial

        self.dndt = self.calc_dndt()
        self.netSuns = self.calc_netSuns()
        self.effSuns = self.calc_effSuns()

        # Match LO and HI measurements and recalclate the suns values
        # Needs to come before all the calculation of the other parameters
        # to assure that the shifted suns are taken into account.
        if not VocLO == 0:
            self.ShiftSuns(VocLO, SunsLO, SunsDataLO)

        self.tauEff = self.calc_tauEff()
        self.jEq = self.calc_jEq()
        self.pDensity = self.calc_pDensity()
        self.Voc = self.calc_Voc()
        self.Vmpp = self.calc_Vmpp()
        self.Jmpp = self.calc_Jmpp()
        self.pFF = self.calc_pFF()
        self.pEff = self.Vmpp * self.Jmpp / 0.1
        self.localIdeality = self.calc_localIdeality(idealityPts)
        self.dsuns = 0.1
        self.mppSuns = self.effSuns[self.mpp] # to get the suns value at mpp
        self.mppIdeality = self.localIdeality[self.mpp] # to get the local ideality factor at mpp

        if not np.size(illumSetpoints) == 0:
            self.getVocVSIllum(illumSetpoints)

    # Last review: 11.10.2019 by J.Seif / WORKING
    def load_data(self):
        """
        Function to load the data from the text file specified by the file path.
        Note that to access the different values you use self.raw['volt'] or any of the other column names.
        """
        self.raw = np.genfromtxt(self.f, delimiter='\t', skip_header=11, names=self.names, usecols=(0,1,2,3))
        
        if self.raw['volt'][0]<0:
            print("The cell voltage data seems to need inverting.")
            self.raw['volt'] = -self.raw['volt']
        if self.raw['ref_volt'][-1]<0:
            print("The reference voltage and suns values rest at a value smaller than zero, this has been corrected.")
          
            # correct by the last millisecond. Use median to avoid outlier effects
            ref_volt_offset = np.median(self.raw['ref_volt'][-10::])
            self.raw['ref_volt'] += abs(ref_volt_offset)
            
            suns_offset = np.median(self.raw['suns'][-10::])
            self.raw['suns'] += abs(suns_offset)
            
    def DistAverageSuns(self, avgLO, avgHI, calConst):

        dsuns = avgLO - calConst * avgHI

    def ShiftSuns(self, VocLO, SunsLO, SunsDataLO):
        """
        This function is used to shift the HI Suns-Voc measurements to the LO
        measurements.

        This is done by finding the measured Voc value in the HI dataset that is closest to Xth value
        in the LO dataset and then minimizing the difference in suns by changing the refCal value
        iteratively until the difference is smaller than a given threshold.

        Inputs:
            - VocLO: is the voltage value of the LO data at a specific point, e.g. 30th point from the start.
                float
            - SunsLO: is the suns value of VocLO
                float
            - SunsDataLO: is the entire SunsData of a specific measurement
                object
        """

        """CHANGE VocLO to input parameter, i.e. point."""

        """Minimize dsuns"""
        # get the index in the HI dataset of the value closest to the Xth Voc in the LO dataset
        # self.VocLO: is the Voc value to which we want to match.
        dVoc = [abs(Val - self.VocLO) for Val in self.raw['volt']]
        indClosestVoc = dVoc.index(min(dVoc))
        self.indClosestVoc = indClosestVoc

        # Linear regression of LO around the top part of the Suns-Voc curve
        VocLOx = SunsDataLO.raw['volt'][0:10]
        SunsLOy = SunsDataLO.netSuns[0:10]
        # LogSunsLOy = [np.log(i) for i in SunsLOy]
        slope, intercept, r_value, p_value, std_err = linregress(VocLOx,SunsLOy)
        # suns value at the highest LO data point
        highSunsLO = intercept + slope *  SunsDataLO.raw['volt'][0]

        # Linear regression of HI around the indClosestVoc
        VocHIx = self.raw['volt'][self.indClosestVoc - 4:self.indClosestVoc + 4]
        SunsHIy = self.netSuns[self.indClosestVoc - 4:self.indClosestVoc + 4]
        slope, intercept, r_value, p_value, std_err = linregress(VocHIx,SunsHIy)
        # suns value at the same voltage as the highest LO data point
        highSunsHI = intercept + slope *  SunsDataLO.raw['volt'][0]

        # difference of the suns values
        dsuns = (highSunsLO - highSunsHI) / highSunsLO
        # set the increment to change the refCal
        # increment  < 0 => increases the suns value
        # increment  > 0 => decreases the suns value
        increment = -0.0001

        # If dsuns is smaller than 0 from the start, change the initial direction
        # of the incrementation.
        if dsuns < 0:
            increment = increment * (-1)

        dsuns_sign = np.sign(dsuns)

        count = 0
        while np.abs(dsuns) > 0.00001:

            # Change the refCal value by a certain increment until dVoc < 0.
            # This means that the shift was too far and the sign of the increment
            # has to be changed.
            count += 1
            # decrease the refCal value
            self.refCal = self.refCal + increment
            # recalculate the netSuns and effSuns
            self.netSuns = self.calc_netSuns()
            # calculate the next dsuns value
            # to get the the [Voc, suns] tuple that falls on the HI data just below the [VocLO, SunsLO]
            # SunsAverage = np.average(self.effSuns[indClosestVoc-3:indClosestVoc+3])
            SunsHIy = self.netSuns[self.indClosestVoc - 4:self.indClosestVoc + 4]
            # LogSunsHIy = [np.log(i) for i in SunsHIy]
            slope, intercept, r_value, p_value, std_err = linregress(VocHIx,SunsHIy)
            highSunsHI = intercept + slope *  SunsDataLO.raw['volt'][0]

            # recalclate the dsuns value
            dsuns = (highSunsLO - highSunsHI) / highSunsLO

            if dsuns > 100:

                increment = increment * (-1)
                # raise Exception("There has been an error when shifting the HI to the LO data.")

            # Change the sign of the increment and the step width.
            if np.sign(dsuns) * dsuns_sign == -1:

                increment = (-1) * increment / 20
                dsuns_sign = np.sign(dsuns)

            # Stop if the counts go beyond X
            if count > 10000:
                break

        # Recalculate the effective suns values
        self.highSunsHI = highSunsHI
        self.effSuns = self.calc_effSuns()
        self.dsuns = dsuns
        self.HIShifted = 1

    def MergeData(self, SunsDataListLO, SunsDataListHI, index):

        """
        This function merges two data sets, LO and HI, after the HI data has been shifted.

        Inputs:
            - LO data set (including all the calculated data), all the data
            - HI data set (including all the calculated data), continues on above LO

        Output:
            - numpy array containing the merged data
        """

        time = np.concatenate((SunsDataListHI[index].raw['t'][0:SunsDataListHI[index].indClosestVoc],
            SunsDataListLO[index].raw['t']),
            axis=0)
        refVolt = np.concatenate((SunsDataListHI[index].raw['ref_volt'][0:SunsDataListHI[index].indClosestVoc],
            SunsDataListLO[index].raw['ref_volt']),
            axis=0)
        volt = np.concatenate((SunsDataListHI[index].raw['volt'][0:SunsDataListHI[index].indClosestVoc],
            SunsDataListLO[index].raw['volt'][:]),
            axis=0)
        dVdt = np.concatenate((SunsDataListHI[index].dVdt[0:SunsDataListHI[index].indClosestVoc],
            SunsDataListLO[index].dVdt[:]),
            axis=0)
        Dn = np.concatenate((SunsDataListHI[index].Dn[0:SunsDataListHI[index].indClosestVoc],
            SunsDataListLO[index].Dn[:]),
            axis=0)
        netSuns = np.concatenate((SunsDataListHI[index].netSuns[0:SunsDataListHI[index].indClosestVoc],
            SunsDataListLO[index].netSuns[:]),
            axis=0)
        effSuns = np.concatenate((SunsDataListHI[index].effSuns[0:SunsDataListHI[index].indClosestVoc],
            SunsDataListLO[index].effSuns[:]),
            axis=0)
        tauEff = np.concatenate((SunsDataListHI[index].tauEff[0:SunsDataListHI[index].indClosestVoc],
            SunsDataListLO[index].tauEff[:]),
            axis=0)
        jEq = np.concatenate((SunsDataListHI[index].jEq[0:SunsDataListHI[index].indClosestVoc],
            SunsDataListLO[index].jEq[:]),
            axis=0)
        pDensity = np.concatenate((SunsDataListHI[index].pDensity[0:SunsDataListHI[index].indClosestVoc],
            SunsDataListLO[index].pDensity[:]),
            axis=0)
        normVoc = volt / SunsDataListLO[index].Voc

        localIdeality = np.concatenate((SunsDataListHI[index].localIdeality[0:SunsDataListHI[index].indClosestVoc],
            SunsDataListLO[index].localIdeality[:]),
            axis=0)

        # RECALCULATING m(V) generates an artefact. Not sure why, but probably due to the stitching of HI and LO.
        # xData = np.log(effSuns)
        # yData = volt
        #
        # localIdeality = self.calc_slope(xData, yData, numpts=self.idealityPts, calcmVoc=0) / (self.const_k * self.TK / self.const_q)

        MergedData = np.array([time,
            refVolt,
            volt,
            dVdt,
            Dn,
            netSuns,
            effSuns,
            tauEff,
            jEq,
            pDensity,
            localIdeality,
            normVoc]
        )

        return MergedData

    def calc_slope(self, xData, yData, numpts=2, calcmVoc=0):
        """
        Function to calculate the slope of two data sets for the n-th and the (n+1)-st data point.

        Inputs:
            - x data, i.e. xData=[any array]
                array
            - y data, i.e. yData=[any array]
                array
            - number of points, i.e. numpts=2 (std.), defines the number of points
                that should be included in the calculation of the slope.

        Output:
            - linear regression between a given number of points (numpts)

        Option:
            - calcmVoc == 0: will output the m value at a given point, using the
                linear fit of a given number of points (numpts)
            - calcmVoc == 1: will calculate the m value at the Voc
        """
        # numData = xData.shape[0] - (numpts - 1) # get the size of the dataset, #rows-1
        # slope = np.zeros(xData.shape[0] - (numpts - 1))
        slope_array = np.array([])

        if calcmVoc == 0:

            for n in range(xData.shape[0]):

                if n + numpts > xData.shape[0]:

                    slope_array = np.append(slope_array, 0)

                else:

                    xn = xData[n:n + numpts]
                    yn = yData[n:n + numpts]

                    slope, intercept, r_value, p_value, std_err = linregress(xn,yn)
                    slope_array = np.append(slope_array, slope)

            return slope_array

        elif calcmVoc == 1:

            # mask the data to get only the points around Voc or 1 sun
            mask = (self.effSuns > 0.8) & (self.effSuns < 1.2)
            xn = np.log(self.effSuns[mask])
            yn = self.raw['volt'][mask]

            slope, intercept, r_value, p_value, std_err = linregress(xn,yn)

            self.mVoc = slope / (self.const_k * self.TK / self.const_q)

    def calc_dVdt(self, **kwargs):
        """
        Function to calculate the variation of the raw voltage with respect to time.

        Inputs:
            - time data
                array
            - cell voltage data
                array
        Output:
            - derivative of the cell voltage with time
                array
        """
        dVdt = self.calc_slope(xData=self.raw['t'], yData=self.raw['volt'], numpts=2, calcmVoc=0)

        return dVdt

    # Last review: 11.10.2019 by J.Seif / WORKING
    def calc_doping(self, **kwargs):
        """
        Function to calculate the doping. Based on the Sinton Instruments Suns-Voc spreadsheet v4.4.

        Inputs:
            - base resistivity (self.rBase) [Ohm.cm]
                float
            - doping type (self.dType) ['p-type' or 'n-type']
                string

        Output:
            - doping density [cm-3]
                float
        """
        if self.dType == 'p-type':
            doping = 10**(-0.0006543 * np.log10(self.rBase)**6
            + 0.000754055 * np.log10(self.rBase)**5 + 0.0093332 * np.log10(self.rBase)**4
            - 0.03469 * np.log10(self.rBase)**3 + 0.06473 * np.log10(self.rBase)**2
            - 1.08286 * np.log10(self.rBase) + 16.17944)
            self.Na = doping

        if self.dType == 'n-type':
            doping = 10**(-0.000634661 * np.log10(self.rBase)**6
            + 0.000820326 * np.log10(self.rBase)**5 + 0.01243 * np.log10(self.rBase)**4
            - 0.04571 * np.log10(self.rBase)**3 + 0.07246 * np.log10(self.rBase)**2
            - 1.07969 * np.log10(self.rBase) + 15.69691)
            self.Nd = doping

        return doping

    # Last review: 11.10.2019 by J.Seif / WORKING
    def calc_ni_initial(self,**kwargs):
        """
        Function to calculate an initial value for n_i for a given temperature.
        This calculation is based on the Semiconductor package of Mattias Juhl.

        Inputs:
            - material, e.g. 'Si'
            - temperature [K], e.g. 300

        Output:
            - intrinsic carrier density n_i [cm-3], number
                Note this is not the effective intrinsic carrier concentration
        """
        ni = ICD(material='Si',temp=self.TK).update()[0]

        return ni

    # Last review: 21.10.2019 by J.Seif
    def calc_ni_eff(self):
        """
        Function to calculate the effective intrinsic carrier density, i.e.
        including the band gap narrowing (BGN).

        Inputs:
            - doping densities donors and acceptors (Nd and Na) [cm-3]
                Note - needs to be the same dimension as the excess carrier
                density.
            - temperature [K]
                float
            - excess carrier density [cm-3]
        Output:
            - n_i,eff data
        """
        nxc = self.Dn
        Na_array = np.array([self.Na] * nxc.shape[0])
        Nd_array = np.array([self.Nd] * nxc.shape[0])
        if self.BGNon == 1:
            ni_eff = BGN(
                material='Si',
                author='Schenk_1988fer',
                temp=self.TK,
                nxc=nxc,
                Na=Na_array,
                Nd=Nd_array
                ).ni_eff(self.ni_initial)
        else:
            ni_eff = self.ni_initial

        return ni_eff

    # Last review: 21.10.2019 by J.Seif / WORKING
    def calc_Dn(self, ni_init):
        """
        Function to calculate the excess carrier density Delta n.

        Inputs:
            - raw voltage measured on the cell [V] taken from the data file
            - temperature [K]
            - initial n_i value [cm-3]
            - doping density [cm-3]

        Output:
            - excess carrier density Dn [cm-3], array
        """
        if ni_init == 1:

            Dn = np.array([])

            # This calculates the Dn values using an initial ni value
            Dn = ((self.calc_doping() ** 2 + 4 * self.calc_ni_initial() ** 2
                * np.exp((self.const_q / (self.const_k * self.TK))
                * self.raw['volt'])) ** 0.5 - self.calc_doping()) / 2

        elif ni_init == 0:

            # This calculates the Dn values using ni_eff
            i = -1
            Dn = np.array([])

            for value in self.ni_eff:

                i += 1
                Dn_calc = ((self.doping ** 2 + 4 * self.ni_eff[i] ** 2
                    * np.exp((self.const_q / (self.const_k * self.TK))
                    * self.raw['volt'][i])) ** 0.5 - self.doping) / 2
                Dn = np.append(Dn, Dn_calc)

        return Dn

    # Last review: 21.10.2019 by J.Seif / WORKING
    def calc_dndt(self):
        """
        Function to calculated dn/dt. The first fromulat is with BGN, the second
        one is without. BGN is the default.

        Inputs:
            - effective ni (including BGN)
                array
            - raw voltage measured on the cell [V]
                array
            - dVdt (calculated before)
                array

        Output:
            - derivative of the excess carrier density with time (dndt)
        """
        dndt = (
            self.ni_eff ** 2 * self.const_q / (self.const_k * self.TK)) * \
            (np.exp(self.raw['volt'] * self.const_q / (self.const_k * self.TK)) / \
            (self.doping ** 2 + 4 * self.ni_eff ** 2 * np.exp(self.raw['volt'] * \
            self.const_q / (self.const_k * self.TK))) ** 0.5) * self.dVdt

        return dndt

    # Last review: 21.10.2019 by J.Seif / WORKING
    def calc_netSuns(self):
        """
        Function to calculate the net suns values according to Sinton.

        Inputs:
            - raw suns, i.e. the suns directly converted from the reference
                photovoltage using the calibration constant [V/suns] supplied
                by Sinton.
            - wafer thickness [cm]
            - dndt (calculated before)
            - short circuit current density jsc [A/cm2]
            - analysis mode ['GEN' or 'QSS']

        Output:
            - net suns
        """
        if self.aMode == "GEN":

            netSuns = self.raw['ref_volt'] / self.refCal - self.const_q * self.wThick * self.dndt / self.jsc

        elif self.aMode == "QSS":

            netSuns = self.raw['ref_volt'] / self.refCal

        return netSuns

    # Last review: 21.10.2019 by J.Seif / WORKING
    def calc_effSuns(self):
        """
        Function to calculate the effective suns values according to Sinton.
        The values are the same as for the net suns, however everything below
        0.001 suns is discaded and set to 0.001 suns.
        """
        effSuns = np.array([])
        for i in np.linspace(0, self.netSuns.shape[0] - 1, self.netSuns.shape[0]):

            if self.netSuns[int(i)] < 0.001:

                effSuns = np.append(effSuns, 0.001)

            else:

                effSuns = np.append(effSuns, self.netSuns[int(i)])

        return effSuns

    # Last review: 21.10.2019 by J.Seif / WORKING
    def calc_tauEff(self):
        """
        Function to calculate the effective lifetime according to Sinton.

        Inputs:
            - excess carrier density (Dn) [cm-3]
            - short circuit current density (jsc) [A/cm2]
            - wafer thickness [cm]

        Output:
            - effective minority carrier density
        """
        tauEff = self.Dn / ((self.jsc * self.effSuns) / (self.const_q * self.wThick))

        return tauEff

    def calc_jEq(self):
        """
        Function to calculate the equivalent current density according to Sinton.

        Inputs:
            - effective suns values (calculated before)
            - short circuit current density (jsc) [A/cm2]

        Output:
            - equivalent current (for pseudo I-V curve) [A/cm2]
        """
        jEq = self.jsc * (1 - self.effSuns)

        return jEq

    def calc_pDensity(self):
        """
        Function to calculate the power density according to Sinton.

        Inputs:
            - equivalent current (for pseudo I-V curve) [A/cm2]
            - raw voltage measured on cell [V]

        Output:
            - power density [W/cm2]
        """
        pDensity = self.jEq * self.raw['volt']

        return pDensity

    def calc_Voc(self):
        """
        Function to calculate the Voc (voltage at 1 sun illumination).

        Inputs:
            - raw volate measured on cell [V]
            - effective suns (calculated before)

        Output:
            - Voc at 1 sun [V]
        """
        fVocSuns = np.array([])
        fVoc = np.array([])
        for i in np.linspace(0, self.raw['volt'].shape[0] - 1, self.raw['volt'].shape[0]):

            if self.effSuns[int(i)] <= 1.2 and self.effSuns[int(i)] >= 0.8:

                fVocSuns = np.append(fVocSuns, self.effSuns[int(i)])
                fVoc = np.append(fVoc, self.raw['volt'][int(i)])

                if self.effSuns[int(i)] < 0.8:

                    break

        # Check if the range around 1 sun is available or not.
        if fVocSuns.size == 0 or fVoc.size == 0:

            # if no suns data is available between 1.2 and 0.8 suns, use the first 20 points
            fVocSuns = self.effSuns[0:25]
            fVoc = self.raw['volt'][0:25]

        slope, intercept, r_value, p_value, std_err = linregress(np.log10(fVocSuns), fVoc)

        Voc = intercept

        return Voc

    def calc_Vmpp(self):
        """
        Function to calculate the voltage at maximum power.

        Inputs:
            - power density [W/cm2]
            - raw voltage measured on cell [V]
        Output:
            - voltage at maxium power point (Vmpp) [V]
        """
        ind_Vmpp = np.argmax(self.pDensity)
        Vmpp = self.raw['volt'][ind_Vmpp]
        self.mpp = ind_Vmpp

        return Vmpp

    def calc_Jmpp(self):
        """
        Function to calculate the current at maximum power.

        Inputs:
            - power density [W/cm2]
            - equivalent current [A/cm2]
        Output:
            - current at maxium power point (Jmpp) [A/cm2]
        """
        ind_Jmpp = np.argmax(self.pDensity)
        Jmpp = self.jEq[ind_Jmpp]

        return Jmpp

    def calc_pFF(self):
        """
        Function to calculate the pseudo FF.

        Inputs:
            - voltage and current at maximum power point
        Output:
            - pseudo fill factor (pFF)
        """
        pFF = self.Vmpp * self.Jmpp / (self.Voc * self.jsc)

        return pFF

    def getVocVSIllum(self, setpoints):
        """
        Gets the Voc values for user-defined illumination levels.

        Inputs:
            - setpoints: numpy array containing the setpoint values for the illumination
            - illum: numpy array containing illumination data
            - voc: numpy array containing Voc data
        Outputs:
            - numpy array containing the Voc values for each setpoint illum.
        """

        setpoints = np.asarray(setpoints)
        setpointslog = np.log(setpoints)
        vocIllum = np.array([])
        mask = np.array([])

        for value in setpoints:
            # print("Sun setpoint: " + str(value))

            idx1 = (np.abs(self.netSuns - value)).argmin()
            # print("idx1: " + str(idx1))
            # print("self.effSuns[idx1]: " + str(self.netSuns[idx1]))
            if self.netSuns[idx1] >= value:
                idx2 = idx1 + 1
                mask = (self.netSuns <= self.netSuns[idx1]) & (self.netSuns >= self.netSuns[idx2])
                # print("idx1: " + str(idx1))
                # print("idx2: " + str(idx2))
            elif self.netSuns[idx1] <= value:
                idx2 = idx1 - 1
                mask = (self.netSuns >= self.netSuns[idx1]) & (self.netSuns <= self.netSuns[idx2])
                # print("idx1: " + str(idx1))
                # print("idx2: " + str(idx2))

            # print(mask)
            xn = np.log(self.netSuns[mask])
            # print(xn)
            yn = self.raw['volt'][mask]

            if np.size(xn) != 0:
                slope, intercept, r_value, p_value, std_err = linregress(xn,yn)
                vocIllum = np.append(vocIllum, np.log(value)*slope + intercept)
            else:
                vocIllum = np.append(vocIllum, 0)

        self.vocIllum = vocIllum

    def calc_localIdeality(self, numpts):
        """
        Function to calculate the local ideality factor.

        Inputs:
            - effective suns
            - raw voltage measured on cell [V]
        Output:
            - local ideality factor m(V), slope of the Voc(log10(suns))
                curve divided by kT/q
        """
        xData = np.log(self.effSuns)
        yData = self.raw['volt']

        m = self.calc_slope(xData, yData, numpts=numpts, calcmVoc=0) / (self.const_k * self.TK / self.const_q)
        self.calc_slope(xData, yData, numpts=numpts, calcmVoc=1)

        # return 1 / ((self.const_k * self.TK / self.const_q) * np.gradient(xData, yData, edge_order=2))

        return m

    def Summarize_Voc(self, Merge):
        """
        This function summarizes all the data (raw and calculated) obtained from the
        Suns-Voc measurements.
        """

        """To determine which operating system is running and which path-separator to use."""
        if platform.system() == 'Windows':

            slash = '\\'

        elif platform.system() == 'Darwin':

            slash = '/'

        check = 0
        if not os.path.isdir(self.sp + slash + "_calcData"):
            os.mkdir(self.sp + slash + "_calcData")

        f = open(os.path.join(self.sp + slash + "_calcData", self.sn + "_T-" + str(np.round(self.T)) + "_SV.txt"), "w")

        for i in range(len(self.raw['t'])):

            if check == 0:
                f.write("{0}\t {1}\t {2}\t {3}\t {4}\t {5}\t {6}\t {7}\t {8}\t {9}\t {10}"
                    .format(
                        "time",
                        "ref V (raw)",
                        "cell V (raw)",
                        "dVdt",
                        "Dn",
                        "net suns",
                        "eff. suns",
                        "tau_eff",
                        "J_eq",
                        "p density",
                        "ideality"
                    )
                )

                f.write("\n")

                f.write("{0}\t {1}\t {2}\t {3}\t {4}\t {5}\t {6}\t {7}\t {8}\t {9}\t {10}"
                    .format(
                        "[s]",
                        "[V]",
                        "[V]",
                        "[V/s]",
                        "[cm-3]",
                        "[suns]",
                        "[suns]",
                        "[s]",
                        "[A/cm2]",
                        "[W/cm2]",
                        "[]"
                    )
                )

                f.write("\n")

                f.write("{0}\t {1}\t {2}\t {3}\t {4}\t {5}\t {6}\t {7}\t {8}\t {9}\t {10}"
                    .format(
                        str(np.round(self.T)) + "°C",
                        str(np.round(self.T)) + "°C",
                        str(np.round(self.T)) + "°C",
                        str(np.round(self.T)) + "°C",
                        str(np.round(self.T)) + "°C",
                        str(np.round(self.T)) + "°C",
                        str(np.round(self.T)) + "°C",
                        str(np.round(self.T)) + "°C",
                        str(np.round(self.T)) + "°C",
                        str(np.round(self.T)) + "°C",
                        str(np.round(self.T)) + "°C"
                    )
                )

                f.write("\n")

                f.write("\n")

                check += 1

            else:

                f.write("{0}\t {1}\t {2}\t {3}\t {4}\t {5}\t {6}\t {7}\t {8}\t {9}\t {10}"
                    .format(
                        self.raw['t'][i],
                        self.raw['ref_volt'][i],
                        self.raw['volt'][i],
                        self.dVdt[i],
                        self.Dn[i],
                        self.netSuns[i],
                        self.effSuns[i],
                        self.tauEff[i],
                        self.jEq[i],
                        self.pDensity[i],
                        self.localIdeality[i]
                    )
                )

                f.write("\n")

        f.close()

    check = 0
