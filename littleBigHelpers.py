# some smaller functions
import os
import codecs # to ignore non-utf-8 characters in the data files
# # for calculations
# import numpy as np
# from scipy import stats
# from SunsVoc_v1_6 import Suns_Voc_Measurement as sv
# import SunsPL_v1_3 as pl
# from glob import glob
#
# # for jupyter
# from PyQt5.QtWidgets import QFileDialog, QApplication
# import ipywidgets as widgets
# from ipywidgets import Layout, Button, Box, VBox, HBox, IntProgress
# from IPython.display import clear_output, display
# import warnings
# warnings.filterwarnings('ignore')
# import matplotlib.pyplot as plt
#
# #other
# import matplotlib.pyplot as plt
# from scipy.stats import linregress
# import littleBigHelpers as lBH
# # import littleBigHelpers import walklevel, getTemp, mainRoutine
#
# # for file handling
# from ipyfilechooser import FileChooser # to choose the *.smpl file
# import os
# import sys
# import platform # to determine the platforme Windows or MacOS
# import time
# import configparser # to read the *.smpl file
# import shutil

def walklevel(some_dir, level=1):
    """
    This function allows to search for files down to a certain
    level of sub-folders.
    """
    some_dir = some_dir.rstrip(os.path.sep)
    assert os.path.isdir(some_dir)
    num_sep = some_dir.count(os.path.sep)
    for root, dirs, files in os.walk(some_dir):
        yield root, dirs, files
        num_sep_this = root.count(os.path.sep)
        if num_sep + level <= num_sep_this:
            del dirs[:]

def getTemp(path):
    """
    The purpose of this function is two-fold:

        1) get the measurement temperature from the data file
        2) remove all non-utf-8 characters from the data file
    """
    global mTemp

    mTemp = 0

    with codecs.open(path,'r',encoding='utf-8',errors='ignore') as fp:
        i = 0
        while i < 5:
            i += 1
            fp.readline()
            if i == 5:
                mTemp = float(fp.readline().split('\t')[1])
    fp.close()

    with codecs.open(path,'r',encoding='utf-8',errors='ignore') as fp:
        entireFile = fp.read()

    fp.close()

    with codecs.open(path, 'w', encoding='utf8') as f:
        f.write(entireFile)

    f.close()

    return mTemp

# def mainRoutine(smpl):
#
#     """Get the sample folder path, i.e. the folder in which all the data is stored."""
#     """Get the path to the *.smpl file that contains some information about the sample and measurement."""
#     sampleFolder = smpl.selected_path
#     smplPath = smpl.selected
#
#     # remove _calcData
#     dirpath = os.path.join(sampleFolder, "_calcData")
#     if os.path.exists(dirpath):
#         shutil.rmtree(dirpath)
#
#     # remove _calcData-merged
#     dirpath = os.path.join(sampleFolder, "_calcData-merged")
#     if os.path.exists(dirpath):
#         shutil.rmtree(dirpath)
#
#     """To determine which operating system is running and which path-separator to use."""
#     if platform.system() == 'Windows':
#
#         slash = '\\'
#
#     elif platform.system() == 'Darwin':
#
#         slash = '/'
#
#     """Extract the parameters from the selected *.smpl file."""
#     with open(smplPath) as fp:
#
#         config = configparser.ConfigParser(allow_no_value=True)
#         config.readfp(fp)
#         waferThickness = float(config['Sample']['waferThickness'])
#         baseResistivity = float(config['Sample']['baseResistivity'])
#         jsc = config['Sample']['jsc'].split(',') # takes single values of jsc or a list separated by ','
#         jsc = [float(i) for i in jsc] # jsc list with float values
#         TCjsc = float(config['Sample']['TCjsc']) # temperature coefficient Jsc
#         JscT = float(config['Sample']['JscT']) # measurement temperature at which the only Jsc value was measured
#         JscTS = config['Sample']['JscTS'].split(',') # takes a list of T values separated by ','
#         JscTS = [float(i) for i in JscTS] # T list with float values
#         TS = config['Sample']['TS'].split(',') # takes single values of TS or a list separated by ','
#         TS = [float(i) for i in TS] # TS list with float values
#         baseType = config['Sample']['baseType']
#         calConstantLO = float(config['Measurement']['calConstantLO'])
#
#     ######################################################################################################
#     # SUNS-VOC DATA ANALYSIS #############################################################################
#     ######################################################################################################
#
#     """sunFiles is the array that contains all the names of the Suns-Voc files in the selected folder."""
#     """These files will be sorted by low (lo) and high (hi) intensity measurements."""
#     sunFilesLO = []
#     sunFilesHI = []
#     sunFilesSTD = []
#     NoFilesFound_lo = 0
#     NoFilesFound_hi = 0
#     NoFilesFound_std = 0
#     SunsDataListLO = []
#     SunsDataListHI = []
#     SunsDataListSTD = []
#
#     numFileList = [0, 0, 0]
#     for root, dirs, files in walklevel(sampleFolder, level=0):
#         for file in files:
#
#             if file.find('_Voc') > -1 and file.find('.sun') > -1:
#                 if file.find('_lo') > -1:
#                     sunFilesLO.append(file)
#                     numFileList[0] += 1
#                 elif file.find('_hi') > -1:
#                     sunFilesHI.append(file)
#                     numFileList[1] += 1
#                 elif file.find('_lo') == -1 and file.find('_hi') == -1:
#                     sunFilesSTD.append(file)
#                     numFileList[2] += 1
#
#     if len(sunFilesLO) == 0:
#         NoFilesFound_lo = 1
#
#     if len(sunFilesHI) == 0:
#         NoFilesFound_hi = 1
#
#     if len(sunFilesSTD) == 0:
#         NoFilesFound_std = 1
#
#     """sunPaths is the array that contains all the paths to the Suns-Voc files."""
#     sunPathsLO = []
#     sunPathsHI = []
#     sunPathsSTD = []
#
#     if NoFilesFound_lo == 0:
#
#         for i in sunFilesLO:
#             sunPathsLO.append(os.path.join(sampleFolder, i))
#
#     if NoFilesFound_hi == 0:
#
#         for i in sunFilesHI:
#             sunPathsHI.append(os.path.join(sampleFolder, i))
#
#     if NoFilesFound_lo == 1 and NoFilesFound_hi == 1:
#
#         for i in sunFilesSTD:
#             sunPathsSTD.append(os.path.join(sampleFolder, i))
#
#     # check if the number jsc values given in the smpl file corresponds to the number of data files
#     # if there's only one value extrapolate the rest of the values
#
#     check = 0
#     checkAgain = 0
#     # get the number of files
#     numFileLists = 3 - NoFilesFound_lo - NoFilesFound_hi - NoFilesFound_std
#     numFiles = int((len(sunFilesLO) + len(sunFilesHI) + len(sunFilesSTD)) / numFileLists)
#
#     if NoFilesFound_lo == 0 or NoFilesFound_hi == 0 or NoFilesFound_std == 0:
#         if not (len(sunFilesLO) == len(jsc)) or not (len(sunFilesHI) == len(jsc)) or not (len(sunFilesSTD) == len(jsc)):
#
#             if len(jsc) == 1:
#                 if TCjsc == 0:
#                     checkAgain = 1
#                 else:
#                     check = 1
#                     njsc = []
#                     # calculate the Jsc values for all the temperatures at which Suns-Voc was measured
#                     for i in range(numFiles):
#                         njsc.append(jsc[0] + (TS[i] - JscT) * TCjsc)
#
#                     # replace the old jsc by the new one that was calculated
#                     jsc = njsc
#
#             elif len(jsc) > 1:
#                 check = 1
#                 slope, intercept, r_value, p_value, std_err = linregress(JscTS, jsc)
#
#                 # get the start value
#                 njsc = []
#                 jsc_start = jsc[0]
#                 # calculate the Jsc values for all the temperatures at which Suns-Voc was measured
#                 for i in range(numFiles):
#                     njsc.append(jsc_start + (TS[i] - JscTS[0]) * slope)
#
#                 # replace the old jsc by the new one that was calculated
#                 jsc = njsc
#
#     if check == 1:
#         print("Please NOTE that the Jsc values were extrapolated using the following values.\n")
#         print("\tJsc [A/cm2]: %s" % str(jsc[0]))
#         print("\tJsc measured at T [°C]: %s" % str(JscT))
#         print("\tTCjsc [A/cm2/°C]: %s\n" % str(TCjsc))
#
#     print("The following Suns-Voc files have been found. LO: %s, HI: %s and STD: %s\n"
#           % (str(numFileList[0]), str(numFileList[1]), str(numFileList[2])))
#
#     if checkAgain == 1:
#         print("No TCjsc value was defined in the *.smpl file. Please correct that and start over.")
#
#     else:
#         # Create the directory
#         if not os.path.isdir(sampleFolder + slash + "_calcData"):
#             os.mkdir(sampleFolder + slash + "_calcData")
#
#         fcellParameters = open(os.path.join(sampleFolder + slash + "_calcData", "Cell_Parameters.txt"), "w")
#         fcellParameters.write("{0}\t {1}\t {2}\t {3}\t {4}\t {5}\t {6}\t {7}\t {8}\t {9}\t {10}\t {11}\n"
#             .format(
#                 "sample",
#                 "T",
#                 "Voc",
#                 "Jsc",
#                 "Vmpp",
#                 "Jmpp",
#                 "pFF",
#                 "pEff",
#                 "suns@mpp",
#                 "m(Vmpp)",
#                 "delta Suns(LO to HI)",
#                 "cal const"
#             )
#         )
#         fcellParameters.write("{0}\t {1}\t {2}\t {3}\t {4}\t {5}\t {6}\t {7}\t {8}\t {9}\t {10}\t {11}\n"
#             .format(
#                 " ",
#                 "[°C]",
#                 "[V]",
#                 "[A/cm2]",
#                 "[V]",
#                 "[A/cm2]",
#                 "[%]",
#                 "[%]",
#                 "[suns]",
#                 " ",
#                 "[V]",
#                 " "
#             )
#         )
#         fcellParameters.write("\n")
#
#         # If LO files have been found
#         if NoFilesFound_lo == 0:
#
#             """Calculate the Suns-Voc LO data"""
#             progSuns = IntProgress(min=0, max=len(sunFilesLO))
#             print("Calculating Suns-Voc LO data... please wait.")
#             display(progSuns)
#
#             # Collect Voc values. Used to shift the PL data.
#             cVoc = np.array([])
#             cVocOneSun = np.array([])
#             cSuns = np.array([])
#
#             j = -1
#             for pathDataLO in sunPathsLO:
#
#                 j += 1
#
#                 """Extract the parameters from the data file."""
#                 mTemp = getTemp(pathDataLO)
#
#                 """
#                 Create new Suns_Voc_Measurement objects, i.e. import the data, do the calculations and
#                 export the calcuated values.
#                 """
#                 SunsData = sv(
#                     sp = sampleFolder, # path to all the measurement files
#                     f = pathDataLO, # file path
#                     fn = sunFilesLO[j], # file name
#                     jsc = jsc[j], # short-circuit current density [A/cm2]
#                     wThick = waferThickness, # wafer thickness [cm]
#                     rBase = baseResistivity, # base resistivity [Ohm.cm]
#                     T = mTemp, # sample temperature [°C]
#                     VocLO = 0, # Voc values from the lo measurments
#                     VocOneSunLO = 0,
#                     SunsLO = 0,
#                     SunsDataLO = 0,
#                     dType = baseType, # n-type or p-type
#                     refCal = calConstantLO, # calibration constant for the reference photodiode [V/suns]
#                     aMode = 'GEN'
#                 )
#
#                 fcellParameters.write("{0}\t {1}\t {2}\t {3}\t {4}\t {5}\t {6}\t {7}\t {8}\t {9}\n"
#                     .format(
#                         SunsData.sn,
#                         np.round(SunsData.T),
#                         np.round(SunsData.Voc,4),
#                         np.round(SunsData.jsc,5),
#                         np.round(SunsData.Vmpp,4),
#                         np.round(SunsData.Jmpp,4),
#                         np.round(SunsData.pFF / 0.01,1),
#                         np.round(SunsData.pEff / 0.01,1),
#                         np.round(SunsData.mppSuns,5),
#                         np.round(SunsData.mppIdeality,3)
#                     )
#                 )
#
#                 # cVoc is used to shift the HI data to the LO data. This is done by finding the
#                 # measured Voc value in the HI dataset that is closest to SunsData.raw['volt'][9]
#                 # value in the LO dataset and then minimizing the difference in suns by changing
#                 # the refCal value iteratively until the difference is smaller than a given threshold.
#                 cVoc = np.append(cVoc, SunsData.raw['volt'][0])
#                 cVocOneSun = np.append(cVocOneSun, SunsData.Voc)
#                 cSuns = np.append(cSuns, SunsData.raw['suns'][0])
#                 SunsDataListLO.append(SunsData)
#                 SunsData.Summarize_Voc(0)
#                 progSuns.value += 1
#
#                 # get the doping for the PL measurements
#                 Na = SunsData.Na
#                 Nd = SunsData.Nd
#
#         if NoFilesFound_hi == 0:
#
#             """Calculate the Suns-Voc HI data"""
#             progSunsHI = IntProgress(min=0, max=len(sunFilesHI))
#             print("Calculating Suns-Voc HI data... please wait.")
#             display(progSunsHI)
#
#             crefCal = np.array([])
#             j = -1
#             for pathDataHI in sunPathsHI:
#
#                 j += 1
#
#                 """Extract the parameters from the data file."""
#                 mTemp = getTemp(pathDataHI)
#
#                 """
#                 Create new Suns_Voc_Measurement objects, i.e. import the data, do the calculations and
#                 export the calcuated values.
#                 """
#                 SunsData = sv(
#                     sp = sampleFolder, # path to all the measurement files
#                     f = pathDataHI, # file path
#                     fn = sunFilesHI[j], # file name
#                     jsc = jsc[j], # short-circuit current density [A/cm2]
#                     wThick = waferThickness, # wafer thickness [cm]
#                     rBase = baseResistivity, # base resistivity [Ohm.cm]
#                     T = mTemp, # sample temperature [°C]
#                     VocLO = cVoc[j], # Voc values of the LO measurments to which the closest value in the HI data will be found
#                     VocOneSunLO = cVocOneSun[j],
#                     SunsLO = cSuns[j], # Suns value to which the HI measurement will be matched
#                     SunsDataLO = SunsDataListLO[j], # the entire raw and calculated data of a LO measurement at the same temperature
#                     dType = baseType, # n-type or p-type
#                     refCal = calConstantLO, # calibration constant for the reference photodiode [V/suns]
#                     aMode = 'GEN'
#                 )
#
#                 fcellParameters.write("{0}\t {1}\t {2}\t {3}\t {4}\t {5}\t {6}\t {7}\t {8}\t {9}\t {10}\t {11}\n"
#                     .format(
#                         SunsData.sn,
#                         np.round(SunsData.T),
#                         np.round(SunsData.Voc,4),
#                         np.round(SunsData.jsc,5),
#                         np.round(SunsData.Vmpp,4),
#                         np.round(SunsData.Jmpp,4),
#                         np.round(SunsData.pFF / 0.01,1),
#                         np.round(SunsData.pEff / 0.01,1),
#                         np.round(SunsData.mppSuns,5),
#                         np.round(SunsData.mppIdeality,3),
#                         np.round(SunsData.dsuns, 4),
#                         np.round(SunsData.refCal, 4)
#                     )
#                 )
#
#                 crefCal = np.append(crefCal, SunsData.refCal)
#                 SunsDataListHI.append(SunsData)
#                 SunsData.Summarize_Voc(0)
#                 progSunsHI.value += 1
#
#                 # get the doping for the PL measurements
#                 Na = SunsData.Na
#                 Nd = SunsData.Nd
#
#         if NoFilesFound_std == 0:
#
#             """Calculate the Suns-Voc STD data"""
#             progSunsSTD = IntProgress(min=0, max=len(sunFilesSTD))
#             print("Calculating Suns-Voc STD data... please wait.")
#             display(progSunsSTD)
#
#             # Collect Voc values. Used to shift the PL data.
#             cVoc = np.array([])
#
#             j = -1
#             for pathDataSTD in sunPathsSTD:
#
#                 j += 1
#
#                 """Extract the parameters from the data file."""
#                 mTemp = getTemp(pathDataSTD)
#
#                 """
#                 Create new Suns_Voc_Measurement objects, i.e. import the data, do the calculations and
#                 export the calcuated values.
#                 """
#                 SunsData = sv(
#                     sp = sampleFolder, # path to all the measurement files
#                     f = pathDataSTD, # file path
#                     fn = sunFilesSTD[j], # file name
#                     jsc = jsc[j], # short-circuit current density [A/cm2]
#                     wThick = waferThickness, # wafer thickness [cm]
#                     rBase = baseResistivity, # base resistivity [Ohm.cm]
#                     T = mTemp, # sample temperature [°C]
#                     VocLO = 0, # Voc values from the lo measurments
#                     VocOneSunLO = 0,
#                     SunsLO = 0,
#                     SunsDataLO = 0,
#                     dType = baseType, # n-type or p-type
#                     refCal = calConstantLO, # calibration constant for the reference photodiode [V/suns]
#                     aMode = 'GEN'
#                 )
#
#                 fcellParameters.write("{0}\t {1}\t {2}\t {3}\t {4}\t {5}\t {6}\t {7}\t {8}\t {9}\n"
#                     .format(
#                         SunsData.sn,
#                         np.round(SunsData.T),
#                         np.round(SunsData.Voc,4),
#                         np.round(SunsData.jsc,5),
#                         np.round(SunsData.Vmpp,4),
#                         np.round(SunsData.Jmpp,4),
#                         np.round(SunsData.pFF / 0.01,1),
#                         np.round(SunsData.pEff / 0.01,1),
#                         np.round(SunsData.mppSuns,5),
#                         np.round(SunsData.mppIdeality,3)
#                     )
#                 )
#
#                 cVoc = np.append(cVoc, SunsData.Voc)
#                 SunsDataListSTD.append(SunsData)
#                 SunsData.Summarize_Voc(0)
#                 progSunsSTD.value += 1
#
#                 # get the doping for the PL measurements
#                 Na = SunsData.Na
#                 Nd = SunsData.Nd
#
#         fcellParameters.close()
#
#         ######################################################################################################
#         # MERGE LO AND HI DATA  ##############################################################################
#         ######################################################################################################
#         if not len(SunsDataListLO) == 0 and not len(SunsDataListHI) == 0:
#
#             if not len(SunsDataListLO) == len(SunsDataListHI):
#                 print("The number of LO and HI is not the same. Check your data.")
#
#             else:
#
#                 progSunsMerge = IntProgress(min=0, max=len(sunFilesLO))
#                 print("Merging the Suns-Voc LO and HI data... please wait.")
#                 display(progSunsMerge)
#
#                 for j in range(len(SunsDataListLO)):
#
#                     MergedData = SunsDataListLO[j].MergeData(SunsDataListLO, SunsDataListHI, j)
#
#                     if not os.path.isdir(SunsDataListLO[j].sp + slash + "_calcData-merged"):
#                         os.mkdir(SunsDataListLO[j].sp + slash + "_calcData-merged")
#
#                     f2 = open(os.path.join(SunsDataListLO[j].sp +
#                                            slash +
#                                            "_calcData-merged",
#                                            SunsDataListLO[j].sn +
#                                            "_LO-T-" +
#                                            str(np.round(SunsDataListLO[j].T)) +
#                                            "_HI-T-" +
#                                            str(np.round(SunsDataListHI[j].T)) +
#                                            "_SV.txt"), "w")
#
#                     check = 0
#
#                     for i in range(len(MergedData[0])):
#
#                         if check == 0:
#                             f2.write("{0}\t {1}\t {2}\t {3}\t {4}\t {5}\t {6}\t {7}\t {8}\t {9}\t {10}"
#                                 .format(
#                                     "time",
#                                     "ref V (raw)",
#                                     "cell V (raw)",
#                                     "dVdt",
#                                     "Dn",
#                                     "net suns",
#                                     "eff. suns",
#                                     "tau_eff",
#                                     "J_eq",
#                                     "p density",
#                                     "ideality"
#                                 )
#                             )
#
#                             f2.write("\n")
#
#                             f2.write("{0}\t {1}\t {2}\t {3}\t {4}\t {5}\t {6}\t {7}\t {8}\t {9}\t {10}"
#                                 .format(
#                                     "[s]",
#                                     "[V]",
#                                     "[V]",
#                                     "[V/s]",
#                                     "[cm-3]",
#                                     "[suns]",
#                                     "[suns]",
#                                     "[s]",
#                                     "[A/cm2]",
#                                     "[W/cm2]",
#                                     "[]"
#                                 )
#                             )
#
#                             f2.write("\n")
#
#                             f2.write("{0}\t {1}\t {2}\t {3}\t {4}\t {5}\t {6}\t {7}\t {8}\t {9}\t {10}"
#                                 .format(
#                                     str(np.round(SunsDataListLO[j].T)) + "°C",
#                                     str(np.round(SunsDataListLO[j].T)) + "°C",
#                                     str(np.round(SunsDataListLO[j].T)) + "°C",
#                                     str(np.round(SunsDataListLO[j].T)) + "°C",
#                                     str(np.round(SunsDataListLO[j].T)) + "°C",
#                                     str(np.round(SunsDataListLO[j].T)) + "°C",
#                                     str(np.round(SunsDataListLO[j].T)) + "°C",
#                                     str(np.round(SunsDataListLO[j].T)) + "°C",
#                                     str(np.round(SunsDataListLO[j].T)) + "°C",
#                                     str(np.round(SunsDataListLO[j].T)) + "°C",
#                                     str(np.round(SunsDataListLO[j].T)) + "°C"
#                                 )
#                             )
#
#                             f2.write("\n")
#
#                             check += 1
#
#                         else:
#
#                             f2.write("{0}\t {1}\t {2}\t {3}\t {4}\t {5}\t {6}\t {7}\t {8}\t {9}\t {10}"
#                                 .format(
#                                     MergedData[0][i],
#                                     MergedData[1][i],
#                                     MergedData[2][i],
#                                     MergedData[3][i],
#                                     MergedData[4][i],
#                                     MergedData[5][i],
#                                     MergedData[6][i],
#                                     MergedData[7][i],
#                                     MergedData[8][i],
#                                     MergedData[9][i],
#                                     MergedData[10][i]
#                                 )
#                             )
#
#                             f2.write("\n")
#
#                     progSunsMerge.value += 1
#
#                     f2.close()
#
#         ######################################################################################################
#         # PL DATA ANALYSIS ###################################################################################
#         ######################################################################################################
#
#     #     for root, dirs, files in walklevel(sampleFolder, level=0):
#     #         for file in files:
#     #             if file.find('_Voc') > -1:
#     #                 if file.find('_lo') > -1:
#     #                     sunFilesLO.append(file)
#     #                     numFileList[0] += 1
#     #                 elif file.find('_hi') > -1:
#     #                     sunFilesHI.append(file)
#     #                     numFileList[1] += 1
#     #                 elif file.find('_lo') == -1 and file.find('_hi') == -1:
#     #                     sunFilesSTD.append(file)
#     #                     numFileList[2] += 1
#
#         """plFiles is the array that contains all the names of the PL data files in the selected folder."""
#         plFilesLO = []
#         plFilesHI = []
#         plFilesSTD = []
#         NoFilesFound_lo_PL = 0
#         NoFilesFound_hi_PL = 0
#         NoFilesFound_std_PL = 0
#         numFileList = [0, 0, 0]
#         PLDataListLO = []
#         PLDataListHI = []
#         PLDataListSTD = []
#
#         for root, dirs, files in walklevel(sampleFolder, level=0):
#             for file in files:
#                 if file.find('_PL') > -1:
#                     if file.find('_lo') > -1:
#                         plFilesLO.append(file)
#                         numFileList[0] += 1
#                     elif file.find('_hi') > -1:
#                         plFilesHI.append(file)
#                         numFileList[1] += 1
#                     elif file.find('_lo') == -1 and file.find('_hi') == -1:
#                         plFilesSTD.append(file)
#                         numFileList[2] += 1
#
#         if len(plFilesLO) == 0:
#             NoFilesFound_lo_PL = 1
#
#         if len(plFilesHI) == 0:
#             NoFilesFound_hi_PL = 1
#
#         if len(plFilesSTD) == 0:
#             NoFilesFound_std_PL = 1
#
#         print("The following Suns-PL files have been found. LO: %s, HI: %s and STD: %s\n"
#               % (str(int(numFileList[0])), str(int(numFileList[1])), str(int(numFileList[2]))))
#
#         if NoFilesFound_lo_PL == 0:
#
#             """plPaths_LO is the array that contains all the paths to the PL data files."""
#             plPaths_LO = []
#             for i in plFilesLO:
#                 plPaths_LO.append(os.path.join(sampleFolder, i))
#
#             """Calculate the Suns-PL data"""
#             progPL = IntProgress(min=0, max=len(plFilesLO))
#             print("Calculating Suns-PL LO data... please wait.")
#             display(progPL)
#
#             # Collect Voc values of the PL LO data. Used to shift the PL HI data.
#             cVocPL = np.array([])
#
#             j = -1
#             for pathData in plPaths_LO:
#
#                 j += 1
#                 """Extract the parameters from the data file."""
#                 mTemp = getTemp(pathData)
#
#                 """
#                 Create new PL_Measurement objects, i.e. import the data, do the calculations and
#                 export the calcuated values.
#                 """
#                 PLData = pl.Voltage_Measurement(
#                     sp = sampleFolder, # path to all the measurement files
#                     f = pathData, # file path
#                     fn = plFilesLO[j], # file name
#                     W = waferThickness,
#                     rBase = baseResistivity,
#                     dType = baseType,
#                     jsc = jsc[j],
#                     refCal = calConstantLO,
#                     bkgnd_corr=0.9,
#                     VocSunsLo = cVoc[j], # this is used to shift the PL data to the Suns-Voc data
#                     T = mTemp, # temperature in °C
#                     aMode = 'GEN'
#                 )
#
#                 cVocPL = np.append(cVocPL, PLData.iVoc_1sun)
#                 PLDataListLO.append(PLData)
#                 PLData.Summarize_PL()
#                 progPL.value += 1
#
#         if NoFilesFound_hi_PL == 0:
#
#             """plPaths_HI is the array that contains all the paths to the PL data files."""
#             plPaths_HI = []
#             for i in plFilesHI:
#                 plPaths_HI.append(os.path.join(sampleFolder, i))
#
#             """Calculate the Suns-PL data"""
#             progPL = IntProgress(min=0, max=len(plFilesHI))
#             print("Calculating Suns-PL HI data... please wait.")
#             display(progPL)
#
#             j = -1
#             for pathData in plPaths_HI:
#
#                 j += 1
#                 """Extract the parameters from the data file."""
#                 mTemp = getTemp(pathData)
#
#                 """
#                 Create new PL_Measurement objects, i.e. import the data, do the calculations and
#                 export the calcuated values.
#                 """
#                 PLData = pl.Voltage_Measurement(
#                     sp = sampleFolder, # path to all the measurement files
#                     f = pathData, # file path
#                     fn = plFilesHI[j], # file name
#                     W = waferThickness,
#                     rBase = baseResistivity,
#                     dType = baseType,
#                     jsc = jsc[j],
#                     refCal = crefCal[j], # this will ensure that the same refCal is used fro the HI Voc and PL data.
#                     bkgnd_corr=0.9,
#                     VocSunsLo = cVocPL[j], # this is used to shift the HI PL data to the LO PL data
#                     T = mTemp, # temperature in °C
#                     aMode = 'GEN'
#                 )
#
#                 PLDataListHI.append(PLData)
#                 PLData.Summarize_PL()
#                 progPL.value += 1
#
#         if NoFilesFound_std_PL == 0:
#
#             """plPaths_STD is the array that contains all the paths to the PL data files."""
#             plPaths_STD = []
#             for i in plFilesSTD:
#                 plPaths_STD.append(os.path.join(sampleFolder, i))
#
#             """Calculate the Suns-PL data"""
#             progPL = IntProgress(min=0, max=len(plFilesSTD))
#             print("Calculating Suns-PL HI data... please wait.")
#             display(progPL)
#
#             j = -1
#             for pathData in plPaths_STD:
#
#                 j += 1
#                 """Extract the parameters from the data file."""
#                 mTemp = getTemp(pathData)
#
#                 """
#                 Create new PL_Measurement objects, i.e. import the data, do the calculations and
#                 export the calcuated values.
#                 """
#                 PLData = pl.Voltage_Measurement(
#                     sp = sampleFolder, # path to all the measurement files
#                     f = pathData, # file path
#                     fn = plFilesSTD[j], # file name
#                     W = waferThickness,
#                     rBase = baseResistivity,
#                     dType = baseType,
#                     jsc = jsc[j],
#                     refCal = calConstantLO,
#                     bkgnd_corr=0.9,
#                     VocSunsLo = cVoc[j], # this is used to shift the STD PL data to the STD Voc data
#                     T = mTemp, # temperature in °C
#                     aMode = 'GEN'
#                 )
#
#                 PLDataListSTD.append(PLData)
#                 PLData.Summarize_PL()
#                 progPL.value += 1
#
#         else:
#
#             print("No PL files have been found.")
#
#         print("All calculations DONE! Check the _calcData folder.")
