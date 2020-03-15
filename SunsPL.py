import numpy as np
import matplotlib.pyplot as plt
from semiconductor.recombination.intrinsic import Radiative
from semiconductor.material.bandgap_narrowing import BandGapNarrowing
from semiconductor.electrical.mobility import Mobility
import scipy.constants as C
from scipy.optimize import minimize_scalar
from scipy.stats import linregress
import os
import re

def implied_J(Jsc, suns):
    return Jsc * (1 - suns)

class Voltage_Measurement:

    def __init__(
        self,
        sp,
        f,
        fn,
        W=0.02,
        rBase=2.8,
        dType='n-type',
        jsc=0.037,
        refCal=0.115273,
        VocSunsLo=0,
        T=25,
        aMode='QSS',
        binn=1,
        bkgnd_corr=0.90,
        bkgnd_loc='end',
        crop_start=0,
        crop_end=1,
        Ai_overwrite=0,
        highestVocPLLO=0
        ):
        """
        This class is used to analyze PL data taken with an InGaAs detector
        at the Suns-PL setup.
        """
        # data-related
        self.sp = sp # path to all the measurement files
        self.f = f # file path
        self.fn = fn # file name
        self.sn = re.findall('\[([^]]+)\]', self.fn)[1] # sample name
        self.names = ['t', 'ref', 'suns', 'PL'] # names for the columns in the data file
        self.load_data()
        self.highestVocPLLO = highestVocPLLO
        # sample-related
        self.Na = 0
        self.Nd = 0
        self.rBase = rBase
        self.dType = dType
        self.doping = self.calc_doping()
        self.W = W # width of the wafers
        self.jsc = jsc # cell current
        # measurement-related
        self.T = T
        self.TK = self.T + 273.15 # sample temperature [K]
        # for the LO measurements
        self.VocSunsLo = VocSunsLo
        self.Ai = self.CalcStartAi("POLY3") # Ai start value based on model LIN or POLY3
        self.Ai_increment = self.Ai / 2 # start increment
        self.Fs = 1e20 # Fs start value
        self.R = 0.
        self.refCal = refCal
        # analysis-related
        self.binn = binn # Binning value
        self.bkgnd_corr = bkgnd_corr # background correction
        self.bkgnd_loc = bkgnd_loc # location of background correctione
        self.crop_start = crop_start # where to start cropping
        self.crop_end = crop_end # where to end cropping the data
        self.aMode = aMode
        self.data()
        # calculated values
        self.nxc_PL = self.nxc_from_PL()
        self.iVoc = self.iVoc_from_nxc(self.nxc_PL)
        # self.tauEff = self.Tau_eff(self.nxc_PL, method='transient')
        self.ideality = self.calc_local_ideality_factor(self.iVoc, self.raw['suns'])
        self.dndt = self.dndt(self.nxc_PL, self.raw['t'])
        self.netSuns = self.calc_netSuns()
        self.effSuns = self.calc_effSuns()
        self.iVoc_1sun = self.calc_iVoc_1sun_PL(corr=0) # iVoc at 1 sun
        self.iVoc_1sun_corr = self.iVoc_1sun

        # Correct the iVoc data from the PL using the Voc at 1 sun from the LO
        # measurement (VocSunsLo).
        dVoc_start = (self.VocSunsLo - self.iVoc_1sun) / self.VocSunsLo
        self.iVocCorr(self.Ai_increment, dVoc_start)

    def load_data(self):
        """
        This function loads all the data from a specified txt file.

        Inputs:
            - data file path
        Output:
            - data in instance of this class
        """
        self.raw = np.genfromtxt(self.f, delimiter='\t', skip_header=11, names=self.names, usecols = (0,1,2,3))

    def data(self, **kwargs):
        """
        Function that conditions the raw data and returns the processed data.
        Processing order: background correction -> binning -> cropping
        """

        # do bagkground correction for reference and PL signal
        corr_index = int(self.raw.shape[0]*self.bkgnd_corr)

        if self.bkgnd_loc == 'start':
            self.raw['ref'] -= np.mean(self.raw['ref'][:corr_index])
            self.raw['PL'] -= np.mean(self.raw['PL'][:corr_index])
        elif self.bkgnd_loc == 'end':
            self.raw['ref'] -= np.mean(self.raw['ref'][corr_index:])
            self.raw['PL'] -= np.mean(self.raw['PL'][corr_index:])

        data = self.crop_data(self.binn_data(self.raw, self.binn), self.crop_start, self.crop_end)
        return data

    def CalcStartAi(self, model):
        """
        Function to determine the start Ai value for the shift of the PL data.

        Inputs:
            - temperature of the measurement

        Outputs:
            - Ai value
        """

        if model == "LIN":

            """Model based on a linear fit of Log10Ai(T)."""
            # start value (LOG10(Ai)) as function of temperature / determined empirically
            slope = -0.042891707
            intercept = 18.42823024
            Log10Ai = self.T * slope + intercept

            Ai = 10 ** Log10Ai

        elif model == "POLY3":

            """Model based on a 3rd order polynomial and Log10Ai(Voc of Suns-Voc)."""
            # Parameters for dVoc = 0.001
            A = 62.414
            B = -66.792
            C = 36.506
            D = 3.6494
            # Parameters for dVoc = 0.05
            # A = 59.907
            # B = -64.289
            # C = 34.974
            # D = 3.8034

            # LOG10Ai when dVoc == 0.05
            x = self.VocSunsLo
            Log10Ai = A * x**3 + B * x**2 + C * x**1 + D

            Ai = 10 ** Log10Ai

        return Ai

    def iVocCorr(self, increment, dVoc_start):
        """
        Function to shift the iVoc vector to the LO Suns-Voc data using the Voc at 1 sun from
        the LO Suns-Voc data and changing the Ai value to shift the iVoc data of the PL.

        Inputs:
            - Ai value (that's the parameter that changes)

        OutputS:
            - iVoc vector
        """
        # If dsuns is smaller than 0 from the start, change the initial direction
        # of the incrementation.
        dVoc = dVoc_start

        if dVoc < 0:
            increment = increment * (-1)

        dVoc_sign = np.sign(dVoc)

        # print("Ai\tdVoc")
        # print("start iVoc PL: " + str(self.iVoc_1sun))
        # print("Voc SV: " + str(self.VocSunsLo))
        count = 0
        # to check the slope of the dVoc(Ai) curve
        # AiList = np.array([])
        # AiList = np.append(AiList, self.Ai)
        # dVocList = np.array([])
        # dVocList = np.append(dVocList, dVoc)

        while np.abs(dVoc) > 0.0001:

            count += 1

            self.Ai = self.Ai + increment
            # AiList = np.append(AiList, self.Ai)

            iVoc_corr = self.iVoc_from_nxc(self.nxc_from_PL())
            self.iVoc_corr = iVoc_corr
            self.iVoc_1sun_corr = self.calc_iVoc_1sun_PL(corr=1)
            
            dVoc = (self.VocSunsLo - self.iVoc_1sun_corr) / self.VocSunsLo
            # dVocList = np.append(dVocList, dVoc)

            # s, i, r, p, std_err = linregress(AiList[len(AiList)-2:len(AiList)], dVocList[len(dVocList)-2:len(dVocList)])

            if dVoc > 100:

                increment = increment * (-1)

            # DEBUG:
            # print(str(self.Ai) + "\t" + str(dVoc))
            # Change the sign of the increment and the step width.
            if np.sign(dVoc) * dVoc_sign == -1:
                increment = (-1) * increment / 10
                dVoc_sign = np.sign(dVoc)

            # Stop if the counts go beyond X
            if count > 10000:
                raise ConvergenceError

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

            netSuns = self.raw['ref'] / self.refCal - C.e * self.W * self.dndt / self.jsc

        elif self.aMode == "QSS":

            netSuns = self.raw['ref'] / self.refCal

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

    def nxc_from_PL(self, filt=slice(None)):
        """
        This function calculates the excess carrier density (nxc).

        Inputs:
            - voltage signal from the PL detector [V]
            - Ai calibration constant, needs to be adapted to fit the Suns-Voc measurement
            - doping [cm-3]
        Output:
            - nxc
        """
        PL = self.data()['PL'][filt]

        def nxc_PL(x):

            return (np.sqrt(self.doping**2 + 4 * self.Ai * PL / self.B_rad(nxc=x)) - self.doping) / 2

        guess_nxc = (np.sqrt(self.doping**2 + 4 * self.Ai * PL / self.B_rad(nxc=np.ones_like(PL), constant=True)) - self.doping) / 2

        nxc = self.find_iteratively(guess_nxc, nxc_PL, verbose=False)

        return nxc

    def iVoc_from_nxc(self, nxc):
        """
        This function calculates the iVoc data from the excess carrier density (nxc).

        Inputs:
            - excess carrier density (nxc) [cm-3]
            - doping [cm-3]
            - ni_eff
        Output:
            - iVoc
        """

        iVoc = self.V_T() * np.log(nxc * (self.doping + nxc) / self.ni_eff(nxc=nxc)**2)
        return iVoc

    def calc_iVoc_1sun_PL(self, corr):
        """
        Function to calculate the Voc (voltage at 1 sun illumination).

        Inputs:
            - raw voltage measured on cell [V]
            - effective suns (calculated before)
            - corr just indicates whether the non-corrected or corrected iVoc data is to be used

        Output:
            - Voc at 1 sun [V]
        """
        fVocSuns = np.array([])
        fVoc = np.array([])

        if corr == 0:

            iVoc = self.iVoc
            range = [int(i) for i in np.linspace(0, self.iVoc.shape[0] - 1, self.iVoc.shape[0])]

        elif corr == 1:

            iVoc = self.iVoc_corr
            range = [int(i) for i in np.linspace(0, self.iVoc_corr.shape[0] - 1, self.iVoc_corr.shape[0])]

        for i in range:

            if self.effSuns[i] <= 1.2 and self.effSuns[i] >= 0.8:

                fVocSuns = np.append(fVocSuns, self.effSuns[i])
                fVoc = np.append(fVoc, iVoc[i])

                if self.effSuns[i] < 0.8:

                    break

        # Check if the range around 1 sun is available or not.
        if fVocSuns.size == 0 or fVoc.size == 0:

            # if no suns data is available between 1.2 and 0.8 suns, use the first 20 points
            fVocSuns = self.effSuns[0:25]

            if corr == 0:

                fVoc = self.iVoc[0:25]

            elif corr == 1:

                fVoc = self.iVoc_corr[0:25]

        slope, intercept, r_value, p_value, std_err = linregress(np.log10(fVocSuns), fVoc)
        iVoc = intercept

        return iVoc

    def Tau_eff(self, nxc, method='general'):
        """
        This function calculates the effective minority carrier lifetime.

        Inputs:
            - time [s]
            - excess carrier density [cm-3]
        output:
            - tau_eff [s]
        """

        if method == 'general':
            tau = nxc / self.gen_net(self.dndt(nxc, self.data()['t']))
        elif method == 'steady-state':
            tau = nxc / self.gen_av()
        elif method == 'transient':
            tau = -1 * nxc / self.dndt(nxc, self.data()['t'])
        else:
            raise ValueError('tau_eff analysis method not recognized.')

        return tau

    def calc_local_ideality_factor(self, Voc, G):
        """
        This function calculates the local ideality factor.

        Inputs:
            - Voc [V]
            - generation [suns]
            - thermal voltage
        output:
            - ideality factor m
        """
        return 1 / (self.V_T()*np.gradient(np.log(G), Voc, edge_order=2))

    def calc_slope(self, xData, yData, numpts=2):
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
        """

        # numData = xData.shape[0] - (numpts - 1) # get the size of the dataset, #rows-1
        # slope = np.zeros(xData.shape[0] - (numpts - 1))
        slope_array = np.array([])

        for n in range(xData.shape[0]):

            if n + numpts > xData.shape[0]:

                slope_array = np.append(slope_array, 0)

            else:

                xn = xData[n:n + numpts]
                yn = yData[n:n + numpts]

                slope, intercept, r_value, p_value, std_err = linregress(xn,yn)
                slope_array = np.append(slope_array, slope)

        return slope_array

    def find_iteratively(self, guess, fn, e=0.00001, verbose=False):
        """
        This function finds values iteratively based on an initial guess.
        """
        diff = 1
        count = 0
        x_i = guess
        while np.mean(diff) > e:
            x_j = fn(x_i)
            diff = np.abs((x_j - x_i) / x_i)
            x_i = x_j
            if verbose:
                print(count, np.mean(diff) )
            count += 1
            if count > 100:
                raise ConvergenceError
        return x_i

    def dndt(self, n, t):
        """
        This function calculates dndt

        Inputs:
            - nxc (calculated)
            - time [s]
        Output:
            - dndt [cm-3/s]
        """
        return np.gradient(n, t, edge_order=2)

    def gen_av(self):
        ref = self.data()['ref']
        return self.Fs * ref * (1 - self.R) / self.W

    def gen_net(self, dndt, suns=True):

        return self.gen_av() - dndt

    def V_T(self):
        """
        Calculates the thermal voltage.

        Inputs:
            - temperature [K]
            - constants k and e
        Output:
            - V_T
        """
        return C.k * self.TK / C.e

    def ni_eff(self, nxc, constant=False):

        if constant == False:
            author = 'Schenk_1988fer'
            ni = 9.65e9
            ni_eff = ni*BandGapNarrowing(nxc=nxc,Nd=self.Nd,Na=self.Na,temp=self.TK,author=author).ni_multiplier()
        else:
            ni_eff = np.ones_like(nxc) * constant

        return ni_eff

    def B_rad(self, nxc, constant=False):

        if constant == False:
            rad = Radiative(nxc=nxc, temp=self.TK, Na=self.Na, Nd=self.Nd, author='Altermatt_2005')
            B = rad.get_B(nxc=nxc)
        elif constant == True:
            B = np.ones_like(nxc) * 4.73e-15
        return B

    def _SSdiffs_Ai(self, Ai, filt):
        """
        Returns the sum of squared differences between the Voc and iVoc curves
        within a range given by the filter. Ai is taken as an input paramter.
        Strictly for Ai optimization routines. Watch out for accidentally changing
        Ai when you don't mean to.
        """
        self.Ai = Ai
        return np.sum(np.power((self.data()['Voc'][filt] - self.iVoc_from_nxc(nxc=self.nxc_from_PL(filt=filt)))/self.data()['Voc'][filt], 2))

    def filter_data(self, data, field, data_range):
        """
        Returns a boolean array where values of the data within the data_range
        are True and values outside the data_range are False.

        data = a numpy array
        field = A string corresponding to a field of a structured array to do the
                filtering on. Can be an empty slice if data is not a structured
                array.
        data_range = a tuple of form (min_val, max_val)
        """

        filt = (data[field] >= data_range[0]) * (data[field] <= data_range[1])
        return filt

    def MergeData(self, SunsPLDataListLO, SunsPLDataListHI, index):

        """
        This function merges two data sets, LO and HI, after the HI data has been shifted.

        Inputs:
            - LO data set (including all the calculated data), all the data
            - HI data set (including all the calculated data), continues on above LO

        Output:
            - numpy array containing the merged data
        """

        highestSunsLO = SunsPLDataListLO[0].effSuns[4]
        dSuns = []
        for i in range(len(SunsPLDataListHI[0].effSuns)):
            dSuns.append(abs(SunsPLDataListHI[0].effSuns[i] - highestSunsLO))
        indSuns = dSuns.index(min(dSuns))

        time = np.concatenate((SunsPLDataListHI[index].raw['t'][0:indSuns],
            SunsPLDataListLO[index].raw['t']),
            axis=0)
        refVolt = np.concatenate((SunsPLDataListHI[index].raw['ref'][0:indSuns],
            SunsPLDataListLO[index].raw['ref']),
            axis=0)
        effSuns = np.concatenate((SunsPLDataListHI[index].effSuns[0:indSuns],
            SunsPLDataListLO[index].effSuns[:]),
            axis=0)
        volt = np.concatenate((SunsPLDataListHI[index].raw['PL'][0:indSuns],
            SunsPLDataListLO[index].raw['PL'][:]),
            axis=0)
        nxc_PL = np.concatenate((SunsPLDataListHI[index].nxc_PL[0:indSuns],
            SunsPLDataListLO[index].nxc_PL[:]),
            axis=0)
        iVoc = np.concatenate((SunsPLDataListHI[index].iVoc[0:indSuns],
            SunsPLDataListLO[index].iVoc[:]),
            axis=0)
        iVoc_corr = np.concatenate((SunsPLDataListHI[index].iVoc_corr[0:indSuns],
            SunsPLDataListLO[index].iVoc_corr[:]),
            axis=0)
        dndt = np.concatenate((SunsPLDataListHI[index].dndt[0:indSuns],
            SunsPLDataListLO[index].dndt[:]),
            axis=0)

        xData = np.log(effSuns)
        yData = iVoc_corr

        localIdeality = self.calc_slope(xData, yData, 13) / (C.k * self.TK / C.e)

        MergedData = np.array([time,
            refVolt,
            effSuns,
            volt,
            nxc_PL,
            iVoc,
            iVoc_corr,
            dndt,
            localIdeality]
        )

        return MergedData

    def Summarize_PL(self):
        """
        This function summarizes all the data (raw and calculated) obtained from the
        Suns-PL measurements.
        """
        if not os.path.isdir(self.sp + "\_calcData"):
            os.mkdir(self.sp + "\_calcData")

        f = open(os.path.join(self.sp + "\_calcData", self.sn + "_T-" + str(np.round(self.T)) + "_PL.txt"), "w")
        check = 0
        for i in range(len(self.raw['t'])):

            if check == 0:
                # f.write("%s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t"
                f.write("{0}\t {1}\t {2}\t {3}\t {4}\t {5}\t {6}\t {7}\t {8}\t {9}"
                        .format(
                            "time",
                            "ref V (raw)",
                            "eff. suns",
                            "PL V (raw)",
                            "Dn",
                            "iVoc",
                            "iVoc_corr",
                            "ideality factor",
                            "dndt",
                            "iVoc@1sun"
                        )
                       )
                f.write("\n")
                f.write("{0}\t {1}\t {2}\t {3}\t {4}\t {5}\t {6}\t {7}\t {8}\t {9}"
                        .format(
                            "[s]",
                            "[V]",
                            "[suns]",
                            "[V]",
                            "[cm-3]",
                            "[V]",
                            "[V]",
                            "[]",
                            "[cm-3/s]",
                            "[V]"
                        )
                       )
                f.write("\n")
                f.write("{0}\t {0}\t {0}\t {0}\t {0}\t {0}\t {0}\t {0}\t {0}\t {1}"
                        .format(
                            str(round(self.T,2)) + "°C",
                            self.iVoc_1sun_corr
                        )
                       )
                f.write("\n")
                check += 1

            f.write("{0}\t {1}\t {2}\t {3}\t {4}\t {5}\t {6}\t {7}\t {8}\t"
                    .format(
                        self.raw['t'][i],
                        self.raw['ref'][i],
                        self.effSuns[i],
                        self.raw['PL'][i],
                        self.nxc_PL[i],
                        self.iVoc[i],
                        self.iVoc_corr[i],
                        self.ideality[i],
                        self.dndt[i]
                    )
                   )
            f.write("\n")

        f.close()

##### LEGACY

    def binn_data(self, data, binn):
        """
        This function bins the data to reduce the noise in the measurement.
        X data points (bin parameter) are averaged over and this average will
        be the new point.

        Inputs:
            - binning parameter
            - data
        Output:
            - binned data
        """
        if binn == 1:
            return data

        binned = np.zeros(data.shape[0] // binn, dtype=data.dtype)

        for name in data.dtype.names:
            for i in range(data.shape[0] // binn):
                binned[name][i] = np.mean(
                    data[name][i * binn:(i + 1) * binn], axis=0)

        return binned

    def crop_data(self, data, start, end):
        """
        This function crops the data.

        Inputs:
            - cropping parameters
            - data
        output:
            - cropped data
        """
        start_index = round(data.shape[0] * start)
        end_index = round(data.shape[0] * end)
        return data[start_index:end_index]
