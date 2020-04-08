# %time for i in range(10): print(i)
# %whos
# ctl+shift+D
# ctrl+D
# ctrl+/
# julia = [i**2 if i%2==0 else i for i in range(10)]
# julia2 = [i**2+j if i%2==0 else j for i,j in zip([0,1,2],[3,4,5])]

#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#---    Initialization
#///////////////////////////////////////////
# %%--  Imports
import numpy as np
from scipy import stats
from SunsVoc import Suns_Voc_Measurement as sv
import SunsPL as pl
from glob import glob
from littleBigHelpers import walklevel, getTemp, printProgressBar
from scipy.stats import linregress

import os
import sys
import platform # to determine the platforme Windows or MacOS
import time
import configparser # to read the *.smpl file
import shutil

import warnings
warnings.filterwarnings("ignore")
# %%-

#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#---    Setup
#///////////////////////////////////////////

"""
This code was developed to load and analyse data from the Suns-Voc(T) and Suns-PL(T) setup.
The way it works is the following:

    1)  The user defines a path to the folder where all the data for one sample
        is stored including a UNIQUE *.smpl file that contains some parameters
        about the sample, e.g. wafer thickness, doping type, resistivity etc.
        Furthermore, the user defines the number of points idealPts that should be
        averaged over to calculate the local ideality factor.

    2)  By executing the code below, the program will automatically load all the
        data files and do the calculations.

        ***Note*** There is a filename CONVENTION
        The standard way to measure should be without the white box, i.e. low
        illumination intensity. This is indicated by '\_lo' in the sample name.
        The high illumination intensity measurements, with the white box and the
        flash lowered to the box edge is to be indicated with a '\_hi' in the
        sample name. For the PL data files we don't make a distinction between
        'lo' and 'hi'. The 'lo' Suns-Voc measurements are used to get a value
        for the Voc@1sun which is used to shift the 'hi' and PL measurements to
        the right position. If no 'lo' neither 'hi' files are found in the sample
        folder the code will load any other Suns-Voc file.

    3)  The program will output all the calculated data as *.txt files in a
        folder '_calcData' and '_calcData-merged'. These files may then be
        imported in any other software (e.g. Origin, Excel, etc.)
"""

# %%--  USER INPUTS
# Number of points that are used to calculate the local ideality factor, i.e.
# DataFolder = the folder in which all the data is stored of one sample is.
# smplFilePath = path to the *.smpl file that contains some information about the sample and measurement.
DataFolder = 'C:\\Users\\z3525973\\Desktop\\Data\\200306_Al-BSF\\'
smplFilePath = ''
for file in os.listdir(DataFolder):
    if file.endswith(".smpl"):
        smplFilePath = os.path.join(DataFolder, file)
# the local slope of the Suns-Voc curve.
idealPts = 8
# If set to True all the filenames of the imported files will be printed.
# This can be used to debug the code and see which data works and which doesn't.
PrintFileNames = False
# %%-

# %%--  HOuSEKEEPING
# Exit if no *.smpl file was found.
if smplFilePath == '':
    print('')
    print('ERROR: No *.smpl file was found. Please create one before running the code.')
    print('')
    exit()

# If a _calcData or _calcData-merged folder exist, delete them
# _calcData = folder that contains the individual files with calculated data
dirpath = os.path.join(DataFolder, "_calcData")
if os.path.exists(dirpath):
    shutil.rmtree(dirpath)
# _calcData-merged = folder that contains the merged LO and HI data
dirpath = os.path.join(DataFolder, "_calcData-merged")
if os.path.exists(dirpath):
    shutil.rmtree(dirpath)
del(dirpath)

# determine the OS (Windows or Mac) and define the which slash to use when saving the data
if platform.system() == 'Windows':
    slash = '\\'
# elif platform.system() == 'Darwin':
#    slash = '/'
# %%-

# %%-- IMPORT *.smpl FILE
with open(smplFilePath) as fp:

    config = configparser.ConfigParser(allow_no_value=True)
    config.read_file(fp)
    waferThickness = float(config['Sample']['waferThickness'])
    baseResistivity = float(config['Sample']['baseResistivity'])
    jsc = config['Sample']['jsc'].split(',') # takes single values of jsc or a list separated by ','
    jsc = [float(i) for i in jsc] # jsc list with float values
    TCjsc = float(config['Sample']['TCjsc']) # temperature coefficient Jsc
    JscT = float(config['Sample']['JscT']) # measurement temperature at which the only Jsc value was measured
    JscTS = config['Sample']['JscTS'].split(',') # takes a list of T values separated by ','
    JscTS = [float(i) for i in JscTS] # T list with float values
    TS = config['Sample']['TS'].split(',') # takes single values of TS or a list separated by ','
    TS = [float(i) for i in TS] # TS list with float values
    baseType = config['Sample']['baseType']
    calConstantLO = float(config['Measurement']['calConstantLO'])

######################################################################################################
# SUNS-VOC DATA ANALYSIS #############################################################################
######################################################################################################

# SunFiles is the array that contains all the names of the Suns-Voc files in the selected folder.
# These files will be sorted by low (lo) and high (hi) intensity measurements.
sunFilesLO = []
sunFilesHI = []
sunFilesSTD = []
NoFilesFound_lo = 0
NoFilesFound_hi = 0
NoFilesFound_std = 0
SunsDataListLO = []
SunsDataListHI = []
SunsDataListSTD = []

# Find all the lo, hi or standard Suns-Voc files
numFileList = [0, 0, 0]
for root, dirs, files in walklevel(DataFolder, level=0):
    for file in files:

        if file.find('_Voc') > -1 and file.find('.sun') > -1:
            if file.find('_lo') > -1:
                sunFilesLO.append(file)
                numFileList[0] += 1
            elif file.find('_hi') > -1:
                sunFilesHI.append(file)
                numFileList[1] += 1
            elif file.find('_lo') == -1 and file.find('_hi') == -1:
                sunFilesSTD.append(file)
                numFileList[2] += 1

if len(sunFilesLO) == 0:
    NoFilesFound_lo = 1

if len(sunFilesHI) == 0:
    NoFilesFound_hi = 1

if len(sunFilesSTD) == 0:
    NoFilesFound_std = 1

# sunPaths is the array that contains all the paths to the Suns-Voc files.
sunPathsLO = []
sunPathsHI = []
sunPathsSTD = []

if NoFilesFound_lo == 0:

    for i in sunFilesLO:
        sunPathsLO.append(os.path.join(DataFolder, i))

if NoFilesFound_hi == 0:

    for i in sunFilesHI:
        sunPathsHI.append(os.path.join(DataFolder, i))

if NoFilesFound_lo == 1 and NoFilesFound_hi == 1:

    for i in sunFilesSTD:
        sunPathsSTD.append(os.path.join(DataFolder, i))

# check if the number jsc values given in the smpl file corresponds to the number of data files
# if there's only one value extrapolate the rest of the values

check = 0
checkAgain = 0
# get the number of files
numFileLists = 3 - NoFilesFound_lo - NoFilesFound_hi - NoFilesFound_std
numFiles = int((len(sunFilesLO) + len(sunFilesHI) + len(sunFilesSTD)) / numFileLists)

if NoFilesFound_lo == 0 or NoFilesFound_hi == 0 or NoFilesFound_std == 0:
    if not (len(sunFilesLO) == len(jsc)) or not (len(sunFilesHI) == len(jsc)) or not (len(sunFilesSTD) == len(jsc)):

        if len(jsc) == 1:
            if TCjsc == 0:
                checkAgain = 1
            else:
                check = 1
                njsc = []
                # calculate the Jsc values for all the temperatures at which Suns-Voc was measured
                for i in range(numFiles):
                    njsc.append(jsc[0] + (TS[i] - JscT) * TCjsc)

                # replace the old jsc by the new one that was calculated
                jsc = njsc

        elif len(jsc) > 1:
            check = 1
            slope, intercept, r_value, p_value, std_err = linregress(JscTS, jsc)
            TCjsc = slope

            # get the start value
            njsc = []
            jsc_start = jsc[0]
            # calculate the Jsc values for all the temperatures at which Suns-Voc was measured
            for i in range(numFiles):
                njsc.append(jsc_start + (TS[i] - JscTS[0]) * slope)

            # replace the old jsc by the new one that was calculated
            jsc = njsc

if check == 1:
    print("Please NOTE that the Jsc values were extrapolated using the following values.\n")
    print("\tJsc [A/cm2]: %s" % str('{:.3e}'.format(jsc[0])))
    print("\tJsc measured at T [°C]: %s" % str(JscT))
    print("\tTCjsc [A/cm2/°C]: %s\n" % str('{:.3e}'.format(TCjsc)))

print("The following Suns-Voc files have been found. LO: %s, HI: %s and STD: %s\n"
      % (str(numFileList[0]), str(numFileList[1]), str(numFileList[2])))

if checkAgain == 1:
    print("No TCjsc value was defined in the *.smpl file. Please correct that and start over.")

else:
    # Create the directory
    if not os.path.isdir(DataFolder + slash + "_calcData"):
        os.mkdir(DataFolder + slash + "_calcData")

    a = smplFilePath.split(slash)[len(smplFilePath.split(slash))-1]
    b = a[:a.find(".smpl")]
    fcellParameters = open(os.path.join(DataFolder + slash + "_calcData", b + "Cell_Parameters.txt"), "w")
    fcellParameters.write("{0}\t {1}\t {2}\t {3}\t {4}\t {5}\t {6}\t {7}\t {8}\t {9}\t {10}\t {11}\n"
        .format(
            "sample",
            "T",
            "Voc",
            "Jsc",
            "Vmpp",
            "Jmpp",
            "pFF",
            "pEff",
            "suns@mpp",
            "m(Vmpp)",
            "m(Voc)",
            "delta Suns(LO to HI)",
            "cal const"
        )
    )
    fcellParameters.write("{0}\t {1}\t {2}\t {3}\t {4}\t {5}\t {6}\t {7}\t {8}\t {9}\t {10}\t {11}\n"
        .format(
            " ",
            "[°C]",
            "[V]",
            "[A/cm2]",
            "[V]",
            "[A/cm2]",
            "[%]",
            "[%]",
            "[suns]",
            " ",
            " ",
            "[V]",
            " "
        )
    )
    fcellParameters.write("\n")

    # If LO files have been found
    if NoFilesFound_lo == 0:

        # Calculate the Suns-Voc LO data
        print("Calculating Suns-Voc LO data... please wait.")

        # Collect Voc values. Used to shift the PL data.
        cVoc = np.array([])
        cVocOneSun = np.array([])
        cSuns = np.array([])

        j = -1
        for pathDataLO in sunPathsLO:

            j += 1

            # Extract the parameters from the data file.
            mTemp = getTemp(pathDataLO)

            # Create new Suns_Voc_Measurement objects, i.e. import the data, do the calculations and
            # export the calcuated values.

            SunsData = sv(
                sp = DataFolder, # path to all the measurement files
                f = pathDataLO, # file path
                fn = sunFilesLO[j], # file name
                jsc = jsc[j], # short-circuit current density [A/cm2]
                wThick = waferThickness, # wafer thickness [cm]
                idealityPts=idealPts, # number of points to calculate the local ideality factor
                rBase = baseResistivity, # base resistivity [Ohm.cm]
                T = mTemp, # sample temperature [°C]
                VocLO = 0, # Voc values from the lo measurments
                VocOneSunLO = 0,
                SunsLO = 0,
                SunsDataLO = 0,
                dType = baseType, # n-type or p-type
                refCal = calConstantLO, # calibration constant for the reference photodiode [V/suns]
                aMode = 'GEN'
            )

            fcellParameters.write("{0}\t {1}\t {2}\t {3}\t {4}\t {5}\t {6}\t {7}\t {8}\t {9}\t {10}\n"
                .format(
                    SunsData.sn,
                    np.round(SunsData.T),
                    np.round(SunsData.Voc,4),
                    np.round(SunsData.jsc,5),
                    np.round(SunsData.Vmpp,4),
                    np.round(SunsData.Jmpp,4),
                    np.round(SunsData.pFF / 0.01,1),
                    np.round(SunsData.pEff / 0.01,1),
                    np.round(SunsData.mppSuns,5),
                    np.round(SunsData.mppIdeality,3),
                    np.round(SunsData.mVoc,3)
                )
            )

            # cVoc is used to shift the HI data to the LO data. This is done by finding the
            # measured Voc value in the HI dataset that is closest to SunsData.raw['volt'][9]
            # value in the LO dataset and then minimizing the difference in suns by changing
            # the refCal value iteratively until the difference is smaller than a given threshold.
            cVoc = np.append(cVoc, SunsData.raw['volt'][0])
            cVocOneSun = np.append(cVocOneSun, SunsData.Voc)
            cSuns = np.append(cSuns, SunsData.raw['suns'][0])
            SunsDataListLO.append(SunsData)
            SunsData.Summarize_Voc(0)

            if PrintFileNames:

                print("\tData imported. " + SunsData.fn)

            else:

                printProgressBar(j+1, len(sunPathsLO))

            # get the doping for the PL measurements
            Na = SunsData.Na
            Nd = SunsData.Nd

        print('\n')

    if NoFilesFound_hi == 0:

        # Calculate the Suns-Voc HI data
        print("Calculating Suns-Voc HI data... please wait.")

        crefCal = np.array([])
        j = -1
        for pathDataHI in sunPathsHI:

            j += 1

            # Extract the parameters from the data file.
            mTemp = getTemp(pathDataHI)


            # Create new Suns_Voc_Measurement objects, i.e. import the data, do the calculations and
            # export the calcuated values.
            SunsData = sv(
                sp = DataFolder, # path to all the measurement files
                f = pathDataHI, # file path
                fn = sunFilesHI[j], # file name
                jsc = jsc[j], # short-circuit current density [A/cm2]
                wThick = waferThickness, # wafer thickness [cm]
                rBase = baseResistivity, # base resistivity [Ohm.cm]
                idealityPts=idealPts, # number of points to calculate the local ideality factor
                T = mTemp, # sample temperature [°C]
                VocLO = cVoc[j], # Voc values of the LO measurments to which the closest value in the HI data will be found
                VocOneSunLO = cVocOneSun[j],
                SunsLO = cSuns[j], # Suns value to which the HI measurement will be matched
                SunsDataLO = SunsDataListLO[j], # the entire raw and calculated data of a LO measurement at the same temperature
                dType = baseType, # n-type or p-type
                refCal = calConstantLO, # calibration constant for the reference photodiode [V/suns]
                aMode = 'GEN'
            )

            fcellParameters.write("{0}\t {1}\t {2}\t {3}\t {4}\t {5}\t {6}\t {7}\t {8}\t {9}\t {10}\t {11}\t {12}\n"
                .format(
                    SunsData.sn,
                    np.round(SunsData.T),
                    np.round(SunsData.Voc,4),
                    np.round(SunsData.jsc,5),
                    np.round(SunsData.Vmpp,4),
                    np.round(SunsData.Jmpp,4),
                    np.round(SunsData.pFF / 0.01,1),
                    np.round(SunsData.pEff / 0.01,1),
                    np.round(SunsData.mppSuns,5),
                    np.round(SunsData.mppIdeality,3),
                    np.round(SunsData.mVoc,3),
                    np.round(SunsData.dsuns, 4),
                    np.round(SunsData.refCal, 4)
                )
            )

            crefCal = np.append(crefCal, SunsData.refCal)
            SunsDataListHI.append(SunsData)
            SunsData.Summarize_Voc(0)

            if PrintFileNames:

                print("\tData imported. " + SunsData.fn)

            else:

                printProgressBar(j+1, len(sunPathsHI))

            # get the doping for the PL measurements
            Na = SunsData.Na
            Nd = SunsData.Nd

        print('\n')

    if NoFilesFound_std == 0:

        # Calculate the Suns-Voc STD data.
        print("Calculating Suns-Voc STD data... please wait.")
        # Collect Voc values. Used to shift the PL data.
        j = -1
        for pathDataSTD in sunPathsSTD:

            j += 1

            # Extract the parameters from the data file.
            mTemp = getTemp(pathDataSTD)

            # Create new Suns_Voc_Measurement objects, i.e. import the data, do the calculations and
            # export the calcuated values.

            SunsData = sv(
                sp = DataFolder, # path to all the measurement files
                f = pathDataSTD, # file path
                fn = sunFilesSTD[j], # file name
                jsc = jsc[j], # short-circuit current density [A/cm2]
                wThick = waferThickness, # wafer thickness [cm]
                rBase = baseResistivity, # base resistivity [Ohm.cm]
                idealityPts=idealPts, # number of points to calculate the local ideality factor
                T = mTemp, # sample temperature [°C]
                VocLO = 0, # Voc values from the lo measurments
                VocOneSunLO = 0,
                SunsLO = 0,
                SunsDataLO = 0,
                dType = baseType, # n-type or p-type
                refCal = calConstantLO, # calibration constant for the reference photodiode [V/suns]
                aMode = 'GEN'
            )

            fcellParameters.write("{0}\t {1}\t {2}\t {3}\t {4}\t {5}\t {6}\t {7}\t {8}\t {9}\t {10}\n"
                .format(
                    SunsData.sn,
                    np.round(SunsData.T),
                    np.round(SunsData.Voc,4),
                    np.round(SunsData.jsc,5),
                    np.round(SunsData.Vmpp,4),
                    np.round(SunsData.Jmpp,4),
                    np.round(SunsData.pFF / 0.01,1),
                    np.round(SunsData.pEff / 0.01,1),
                    np.round(SunsData.mppSuns,5),
                    np.round(SunsData.mppIdeality,3),
                    np.round(SunsData.mVoc,3)
                )
            )

            SunsDataListSTD.append(SunsData)
            SunsData.Summarize_Voc(0)

            # get the doping for the PL measurements
            Na = SunsData.Na
            Nd = SunsData.Nd

            if PrintFileNames:

                print("\tData imported. " + SunsData.fn)

            else:

                printProgressBar(j+1, len(sunPathsSTD))

        print('\n')

    fcellParameters.close()

    ######################################################################################################
    # MERGE LO AND HI DATA  ##############################################################################
    ######################################################################################################
    if not len(SunsDataListLO) == 0 and not len(SunsDataListHI) == 0:

        if not len(SunsDataListLO) == len(SunsDataListHI):
            print("The number of LO and HI is not the same. Check your data.")

        else:

            print("Merging the Suns-Voc LO and HI data... please wait.")

            for j in range(len(SunsDataListLO)):

                MergedData = SunsDataListLO[j].MergeData(SunsDataListLO, SunsDataListHI, j)

                if not os.path.isdir(SunsDataListLO[j].sp + slash + "_calcData-merged"):
                    os.mkdir(SunsDataListLO[j].sp + slash + "_calcData-merged")

                if not os.path.isdir(SunsDataListLO[j].sp + slash + "_calcData-merged" + slash + "SV"):
                    os.mkdir(SunsDataListLO[j].sp + slash + "_calcData-merged" + slash + "SV")

                f2 = open(os.path.join(SunsDataListLO[j].sp +
                                       slash +
                                       "_calcData-merged" + slash + "SV",
                                       SunsDataListLO[j].sn[:SunsDataListLO[j].sn.find("_lo")] +
                                       "_T-" +
                                       str(np.round(SunsDataListLO[j].T)) +
                                       "_SV.txt"), "w")

                check = 0

                for i in range(len(MergedData[0])):

                    if check == 0:
                        f2.write("{0}\t {1}\t {2}\t {3}\t {4}\t {5}\t {6}\t {7}\t {8}\t {9}\t {10}\t {11}"
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
                                "ideality",
                                "cell V/Voc(1 sun)"
                            )
                        )

                        f2.write("\n")

                        f2.write("{0}\t {1}\t {2}\t {3}\t {4}\t {5}\t {6}\t {7}\t {8}\t {9}\t {10}\t {11}"
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
                                "[]",
                                "[]"
                            )
                        )

                        f2.write("\n")

                        f2.write("{0}\t {1}\t {2}\t {3}\t {4}\t {5}\t {6}\t {7}\t {8}\t {9}\t {10}\t {11}"
                            .format(
                                str(np.round(SunsDataListLO[j].T)) + "°C",
                                str(np.round(SunsDataListLO[j].T)) + "°C",
                                str(np.round(SunsDataListLO[j].T)) + "°C",
                                str(np.round(SunsDataListLO[j].T)) + "°C",
                                str(np.round(SunsDataListLO[j].T)) + "°C",
                                str(np.round(SunsDataListLO[j].T)) + "°C",
                                str(np.round(SunsDataListLO[j].T)) + "°C",
                                str(np.round(SunsDataListLO[j].T)) + "°C",
                                str(np.round(SunsDataListLO[j].T)) + "°C",
                                str(np.round(SunsDataListLO[j].T)) + "°C",
                                str(np.round(SunsDataListLO[j].T)) + "°C",
                                str(np.round(SunsDataListLO[j].T)) + "°C"
                            )
                        )

                        f2.write("\n")

                        check += 1

                    else:

                        f2.write("{0}\t {1}\t {2}\t {3}\t {4}\t {5}\t {6}\t {7}\t {8}\t {9}\t {10}\t {11}"
                            .format(
                                MergedData[0][i],
                                MergedData[1][i],
                                MergedData[2][i],
                                MergedData[3][i],
                                MergedData[4][i],
                                MergedData[5][i],
                                MergedData[6][i],
                                MergedData[7][i],
                                MergedData[8][i],
                                MergedData[9][i],
                                MergedData[10][i],
                                MergedData[11][i]
                            )
                        )

                        f2.write("\n")

                f2.close()

                printProgressBar(j+1, len(SunsDataListLO))

            print('\n')

######################################################################################################
# PL DATA ANALYSIS ###################################################################################
######################################################################################################

# plFiles is the array that contains all the names of the PL data files in the selected folder.
plFilesLO = []
plFilesHI = []
plFilesSTD = []
NoFilesFound_lo_PL = 0
NoFilesFound_hi_PL = 0
NoFilesFound_std_PL = 0
numFileList = [0, 0, 0]
PLDataListLO = []
PLDataListHI = []
PLDataListSTD = []
AiValsLO = []
highestVocPLLO = []

for root, dirs, files in walklevel(DataFolder, level=0):
    for file in files:
        if file.find('_PL') > -1:
            if file.find('_lo') > -1:
                plFilesLO.append(file)
                numFileList[0] += 1
            elif file.find('_hi') > -1:
                plFilesHI.append(file)
                numFileList[1] += 1
            elif file.find('_lo') == -1 and file.find('_hi') == -1:
                plFilesSTD.append(file)
                numFileList[2] += 1

if len(plFilesLO) == 0:
    NoFilesFound_lo_PL = 1

if len(plFilesHI) == 0:
    NoFilesFound_hi_PL = 1

if len(plFilesSTD) == 0:
    NoFilesFound_std_PL = 1

print("The following Suns-PL files have been found. LO: %s, HI: %s and STD: %s\n"
      % (str(int(numFileList[0])), str(int(numFileList[1])), str(int(numFileList[2]))))

if NoFilesFound_lo_PL == 0:

    # plPaths_LO is the array that contains all the paths to the PL data files.
    plPaths_LO = []
    for i in plFilesLO:
        plPaths_LO.append(os.path.join(DataFolder, i))

    # Calculate the Suns-PL data
    print("Calculating Suns-PL LO data... please wait.")

    # Collect Voc values of the PL LO data. Used to shift the PL HI data.
    cVocPL = np.array([])

    j = -1
    for pathData in plPaths_LO:

        j += 1
        # Extract the parameters from the data file.
        mTemp = getTemp(pathData)

        # Create new PL_Measurement objects, i.e. import the data, do the calculations and
        # export the calcuated values.

        PLData = pl.Voltage_Measurement(
            sp = DataFolder, # path to all the measurement files
            f = pathData, # file path
            fn = plFilesLO[j], # file name
            W = waferThickness,
            rBase = baseResistivity,
            dType = baseType,
            jsc = jsc[j],
            refCal = calConstantLO,
            bkgnd_corr=0.9,
            VocSunsLo = cVocOneSun[j], # this is used to shift the PL data to the Suns-Voc data
            T = mTemp, # temperature in °C
            aMode = 'GEN',
            Ai_overwrite=0,
            highestVocPLLO=0
        )

        cVocPL = np.append(cVocPL, PLData.iVoc_1sun_corr)
        AiValsLO.append(PLData.Ai)
        highestVocPLLO.append(PLData.iVoc_corr[0])
        PLDataListLO.append(PLData)
        PLData.Summarize_PL()

        if PrintFileNames:

            print("\tData imported. " + PLData.fn)

        else:

            printProgressBar(j+1, len(plPaths_LO))

    print('\n')

if NoFilesFound_hi_PL == 0:

    # plPaths_HI is the array that contains all the paths to the PL data files.
    plPaths_HI = []
    for i in plFilesHI:
        plPaths_HI.append(os.path.join(DataFolder, i))

    # Calculate the Suns-PL data
    print("Calculating Suns-PL HI data... please wait.")

    j = -1
    for pathData in plPaths_HI:

        j += 1
        # Extract the parameters from the data file.
        mTemp = getTemp(pathData)

        # Create new PL_Measurement objects, i.e. import the data, do the calculations and
        # export the calcuated values.

        PLData = pl.Voltage_Measurement(
            sp = DataFolder, # path to all the measurement files
            f = pathData, # file path
            fn = plFilesHI[j], # file name
            W = waferThickness,
            rBase = baseResistivity,
            dType = baseType,
            jsc = jsc[j],
            refCal = crefCal[j], # this will ensure that the same refCal is used fro the HI Voc and PL data.
            bkgnd_corr=0.9,
            VocSunsLo = cVocPL[j], # this is used to shift the HI PL data to the LO PL data
            T = mTemp, # temperature in °C
            aMode = 'GEN',
            Ai_overwrite=AiValsLO[j],
            highestVocPLLO=highestVocPLLO[j]
        )

        PLDataListHI.append(PLData)
        PLData.Summarize_PL()

        if PrintFileNames:

            print("\tData imported. " + PLData.fn)

        else:

            printProgressBar(j+1, len(plPaths_HI))

    print('\n')

if NoFilesFound_std_PL == 0:

    # plPaths_STD is the array that contains all the paths to the PL data files.
    plPaths_STD = []
    for i in plFilesSTD:
        plPaths_STD.append(os.path.join(DataFolder, i))

    # Calculate the Suns-PL data
    print("Calculating Suns-PL HI data... please wait.")

    j = -1
    for pathData in plPaths_STD:

        j += 1
        # Extract the parameters from the data file.
        mTemp = getTemp(pathData)

        # Create new PL_Measurement objects, i.e. import the data, do the calculations and
        # export the calcuated values.

        PLData = pl.Voltage_Measurement(
            sp = DataFolder, # path to all the measurement files
            f = pathData, # file path
            fn = plFilesSTD[j], # file name
            W = waferThickness,
            rBase = baseResistivity,
            dType = baseType,
            jsc = jsc[j],
            refCal = calConstantLO,
            bkgnd_corr=0.9,
            VocSunsLo = cVoc[j], # this is used to shift the STD PL data to the STD Voc data
            T = mTemp, # temperature in °C
            aMode = 'GEN',
            Ai_overwrite=0,
            highestVocPLLO=0
        )

        PLDataListSTD.append(PLData)
        PLData.Summarize_PL()

        if PrintFileNames:

            print("\tData imported. " + PLData.fn)

        else:

            printProgressBar(j+1, len(plPaths_STD))

    print('\n')

else:

    print("No PL files have been found.")

######################################################################################################
# MERGE LO AND HI DATA  ##############################################################################
######################################################################################################
if not len(PLDataListLO) == 0 and not len(PLDataListHI) == 0:

    if not len(PLDataListLO) == len(PLDataListHI):
        print("The number of LO and HI is not the same. Check your data.")

    else:

        print("Merging the Suns-PL LO and HI data... please wait.")

        for j in range(len(PLDataListLO)):

            MergedData = PLDataListLO[j].MergeData(PLDataListLO, PLDataListHI, j)

            if not os.path.isdir(PLDataListLO[j].sp + slash + "_calcData-merged"):
                os.mkdir(PLDataListLO[j].sp + slash + "_calcData-merged")

            if not os.path.isdir(PLDataListLO[j].sp + slash + "_calcData-merged" + slash + "PL"):
                os.mkdir(PLDataListLO[j].sp + slash + "_calcData-merged" + slash + "PL")

            f2 = open(os.path.join(PLDataListLO[j].sp +
                                   slash +
                                   "_calcData-merged" + slash + "PL",
                                   PLDataListLO[j].sn[:PLDataListLO[j].sn.find("_lo")] +
                                   "_T-" +
                                   str(np.round(PLDataListLO[j].T)) +
                                   "_PL.txt"), "w")

            check = 0

            for i in range(len(MergedData[0])):

                if check == 0:
                    f2.write("{0}\t {1}\t {2}\t {3}\t {4}\t {5}\t {6}\t {7}\t {8}"
                        .format(
                            "time",
                            "ref V (raw)",
                            "eff. suns",
                            "PL V (raw)",
                            "Dn",
                            "iVoc",
                            "iVoc_corr",
                            "dndt",
                            "ideality"
                        )
                    )

                    f2.write("\n")

                    f2.write("{0}\t {1}\t {2}\t {3}\t {4}\t {5}\t {6}\t {7}\t {8}"
                        .format(
                            "[s]",
                            "[V]",
                            "[suns]",
                            "[V]",
                            "[cm-3]",
                            "[V]",
                            "[V]",
                            "[cm-3/s]",
                            "[]"
                        )
                    )

                    f2.write("\n")

                    f2.write("{0}\t {1}\t {2}\t {3}\t {4}\t {5}\t {6}\t {7}\t {8}"
                        .format(
                            str(np.round(SunsDataListLO[j].T)) + "°C",
                            str(np.round(SunsDataListLO[j].T)) + "°C",
                            str(np.round(SunsDataListLO[j].T)) + "°C",
                            str(np.round(SunsDataListLO[j].T)) + "°C",
                            str(np.round(SunsDataListLO[j].T)) + "°C",
                            str(np.round(SunsDataListLO[j].T)) + "°C",
                            str(np.round(SunsDataListLO[j].T)) + "°C",
                            str(np.round(SunsDataListLO[j].T)) + "°C",
                            str(np.round(SunsDataListLO[j].T)) + "°C"
                        )
                    )

                    f2.write("\n")

                    check += 1

                else:

                    f2.write("{0}\t {1}\t {2}\t {3}\t {4}\t {5}\t {6}\t {7}\t {8}"
                        .format(
                            MergedData[0][i],
                            MergedData[1][i],
                            MergedData[2][i],
                            MergedData[3][i],
                            MergedData[4][i],
                            MergedData[5][i],
                            MergedData[6][i],
                            MergedData[7][i],
                            MergedData[8][i]
                        )
                    )

                    f2.write("\n")

            f2.close()

            printProgressBar(j+1, len(PLDataListLO))

print("\nAll calculations DONE and data EXPORTED!")
