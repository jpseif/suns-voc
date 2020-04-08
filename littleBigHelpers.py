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

# Print iterations progress
def printProgressBar(iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = printEnd)
    # Print New Line on Complete
    if iteration == total:
        print()
