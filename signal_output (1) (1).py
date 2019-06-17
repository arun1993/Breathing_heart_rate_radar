#!/usr/bin/env python
# coding: utf-8

# In[11]:


import csv
import matlab.engine
import math
import numpy as np
from scipy.fftpack import fft
with open('Desktop/rawdata.csv', 'r') as f:
    reader = csv.reader(f)
    signal_list = list(reader)
#print(signal_list)


# In[1]:


from pylive import live_plotter
import numpy as np

size = 100
x_vec = np.linspace(0,1,size+1)[0:-1]
y_vec = np.random.randn(len(x_vec))
line1 = []
while True:
    rand_val = np.random.randn(1)
    y_vec[-1] = rand_val
    line1 = live_plotter(x_vec,y_vec,line1)
    y_vec = np.append(y_vec[1:],0.0)


# In[2]:


time_sec = float(signal_list[2][1]) - float(signal_list[1][1])
freq = round(float(signal_list[-1][1])/time_sec)
actual_freq = round(freq/2)
print(actual_freq)


# In[3]:


signal_list = signal_list[1:]
#print(signal_list)
signal = [float(x[0]) for x in signal_list]
signal_np = np.asarray(signal)


# In[4]:


from scipy.signal import butter, lfilter, filtfilt
def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a
def butter_bandpass_filter(data, lowcut, highcut, fs, order=1):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = filtfilt(b, a, data)
    return y


# In[5]:


FFT_length = math.pow(2,16)
fft_signal = abs(fft(signal_np,int(FFT_length)))
fft_signal = np.round(fft_signal/max(fft_signal),4)
print(np.around(fft_signal,4))


# In[6]:


y1 = butter_bandpass_filter(fft_signal,1,100,freq)
f_y1 = abs(fft(y1,int(FFT_length)))
f_y1 = f_y1/max(f_y1)
print(y1)
print(f_y1)


# In[7]:


def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

def butter_lowpass_filter(data, cutoff, fs, order=5):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = lfilter(b, a, data)
    return y


# In[8]:


#passband = [1, 100]
#b, a = butter(5, passband, 'bandpass')
#y2= filtfilt(b, a, y1)
y2= butter_lowpass_filter(y1, 3, freq);
f_y2 = abs(fft(y2, FFT_length));
f_y2 = f_y2/max(f_y2);
print(y2)
print(f_y2)


# In[9]:


xdft = f_y2[:actual_freq]
print((xdft))
I = max(abs(xdft))
print(I)
print(actual_freq)
print((freq/FFT_length))
final_freq = np.arange(0,actual_freq,float(freq/FFT_length))
print(60*final_freq[19])

