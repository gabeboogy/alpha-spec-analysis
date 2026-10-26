from scipy.optimize import curve_fit
from scipy.stats import linregress
import scipy
#%matplotlib ipympl
import matplotlib.pyplot as plt
import pandas as pd
import csv
import numpy as np
import lcdatabase

file_input = "Berkley Target 5 Point 3 without Bismuth.spe"
major_peaks_input = 3
bins = 4096
