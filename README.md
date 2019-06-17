# Breathing_heart_rate_radar
Breathing, Heart rate detection using Walabot Radar

The signal reflection from the human body is used to measure the heart rate. The phase of this reflected signal is proportional to the distance traveled by the signal. The distance measurement experiment was done in the last report. Also from the below figure we can see variations due to inhaling, exhaling and heart-beats. A person’s heart beats cause miniature movements on human body which is captured in those small variations seen in the graph. These small variations are periodic in nature due to breathing and heart rate being periodic activities. 
By performing an FFT on the distance variations, breathing and heart rate can be detected. Once FFT is performed on the phase of the signal over a period of time, the peak of the signal output gives an estimate of the person’s breathing rate. 

To obtain a precise measurement, FFT output is filtered to have only peak and its two adjacent bins. This filtering approach eliminates the noise by other human movements. An inverse FFT is performed to convert the signal back to time-domain. By measuring the slope of this time-domain signal, breathing rate can be accurately measured. 
Heartbeat are minute variations and its magnitude is very less compared to the breathing signal. This can’t be detected using a FFT as the higher amplitude signal masks this minute variations. After FFT is performed over a shorter time on the phase of the reflected signal, a band pass filter is used to get the beats in human heartbeat range. Breathing rate is less than 10 and noise is more than 300. By having the band pass filter, we can measure the heart rate. 

To summarize this section:
#Removal of Background noise:To remove background noise, samples are collected in the environment over a period of time and average is found Navg. When samples are collected for identifying heart rate, the noise value is subtracted from that.  S = S - Navg. This step can be avoided if calibration is done.
#Periodicity detection: To detect the periodicity, FFT operation is performed to convert the signals from time domain to  frequency domain. This will make the non-periodic signal to have less frequency. 
#Filtering collected samples: Using a butterworth bandpass filter, filtering is done to remove the samples in unwanted frequency range. 
#Heart Rate detection: By finding peaks, we can detect the heart rate of the person. 
