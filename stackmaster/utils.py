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
def power2pad(data):
	"""Zero pad data such that its length is a power of 2"""
	N=int(2**np.ceil(np.log2(len(data))))
	pad_end=np.zeros(int(N-len(data)))

	return np.concatenate((data,pad_end))
