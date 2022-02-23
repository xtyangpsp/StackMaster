import scipy
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import hilbert
from scipy.fftpack import fft,ifft,next_fast_len
from stockwell import st
from tslearn.utils import to_time_series_dataset
from tslearn.clustering import TimeSeriesKMeans

###################
#####Utility functions#####
################
def rms(d):
    return np.sqrt(np.mean(d**2))
