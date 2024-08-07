# -*- coding: utf-8 -*-
"""
 Copyright (c) 2021 CSAN_LiU

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <https://www.gnu.org/licenses/>.
 """


"""
Created on Fri Dec 18 11:34:19 2020

@author: ilosz01

Many of the analysis ideas are inspired on code provided by the following sources:
Tucker-Davis Technologies ( https://www.tdt.com/support/python-sdk/offline-analysis-examples/fiber-photometry-epoch-averaging-example/)
Dr. David Barker (pMAT; Bruno C.A. et al. 2021, Pharm BioChem Behav 201, https://doi.org/10.1016/j.pbb.2020.173093)
Dr. Patrik Mulholland (Braunsheidel K.M. et al. 2019, J Neurosci 39(46), https://doi.org/10.1523/JNEUROSCI.1674-19.2019)

Many analysis ideas come from TDT
https://www.tdt.com/support/python-sdk/offline-analysis-examples/fiber-photometry-epoch-averaging-example/

"""
# https://github.com/LABSN/tdtpy
# https://www.tdt.com/support/python-sdk/offline-analysis-examples/introduction-to-python-tdt-package/
import tdt
from tdt.TDTfilter import combine_time
from tdt.TDTfilter import get_valid_ind
# https://github.com/LABSN/tdtpy
# https://www.tdt.com/support/python-sdk/offline-analysis-examples/introduction-to-python-tdt-package/
import numpy as np
from numpy.polynomial import Polynomial
from numpy.polynomial.polynomial import polyval
import math
# https://pypi.org/project/lowess/
#import lowess
import pandas as pd
# from scipy.ndimage import gaussian_filter
from scipy import stats
from sklearn.metrics import auc
#from loess.loess_1d import loess_1d
#import statsmodels.api as sm
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.signal import find_peaks, filtfilt
from scipy.interpolate import interp1d
import os
import copy
import warnings


''' get all data from the recording (reads from multiple files)
    read_block returns a structured object. It is a Python dictionary 
    
    epocs	[struct]
    snips	[struct]
    streams	[struct]
    scalars	[struct]
    info	[struct]
    time_ranges:	array([[ 0.],
       [inf]])

    Streams is the recording (names: x465A and  x405A)
    epocs[PrtA] has the codes for events (i.e tone)
    
    examples of how to access data:
        print(data.streams.Wav1.fs) # dot syntax
        print(data['streams']['Wav1']['fs'])
        print(data['streams'].Wav1['fs']) # mix of dot syntax and dict keys
'''

#################################
# PLOT GENERAL SETTINGS
##############################
# location of all plots legend
MY_LEGEND_LOC='upper left'
MY_LEGENG_POS = (0.96, 1.1)
MY_LEFT = 0.06
MY_RIGHT = 0.9
MY_TOP = 0.9
MY_BOTTOM = 0.1
MY_WSPACE = 0.2
MY_HSPACE = 0.5

AXIS_LABEL_FONTSIZE = 14
TITLE_FONT_SIZE = 15
FIGURE_TITLE_FONT_SIZE = 16

SIGNAL_COLOR_RGB = '#1A36F2'
CONTROL_COLOR_RGB = '#FC2908'

DPI4SVG = 1200
################################
# return raw data structure
def get_raw_data(path):
    '''Returns a raw data extracted by tdt from all recording files
    '''
    try:
        raw_data = tdt.read_block(path) 
    except:
        raw_data = None
    return raw_data

def get_channel_names(raw_data):
    '''
    Returns a list of unique channel names in streams from file
    i.e. _405A, _465A
    '''
    all_channel_names = [key for key in raw_data.streams.keys()]
    return all_channel_names

def get_events(raw_data):
    '''
    Returns a list of unique events from file that
    don't start with cam or tick
    '''
    # tdt StructTypes are dictionaries
    all_epocs_list = [key for key in raw_data.epocs.keys()]
    my_event_list = []
    # iterate over all events that start with Ptr
    for evt in all_epocs_list:
        # exclude cam or tickfrom events
        if (evt.lower().startswith("cam") == False and evt.lower().startswith("tick") == False) or evt.lower().startswith("cam1") == True:
            if evt.lower().startswith("cam1") == False:
                # create set of unique events
                unique_events = set(raw_data.epocs[evt].data)
                for el in unique_events:
                    evt_name = evt + " " + str(int(el))
                    my_event_list.append(evt_name)
            else:
                # check if include special case for cam1 (only if it has notes)
                data = raw_data.epocs[evt]
                match = False
                for k in data.keys():
                    if k == 'notes':
                        match = True
                if match == True:
                    my_event_list.append(evt)
    return my_event_list

# returns true if there are any events in the data
def check_events(path):
    events_present = False
    data = get_raw_data(path)
    if len(get_events(data)) > 0:
        events_present = True
    return events_present

def get_event_on_off(raw_data, event):
    ''' Returns onset and offset times of the tones
        First element of this array are onset times, 
        Second el are offset times
    '''
    event_name_split = event.split(" ")
    on_off = [[],[]]
    if event.lower().startswith("cam1") == False:
        # find tone onsets and offsets
        try:    # some data was not readable by tdt.epoc_filter
            evt_on_off= tdt.epoc_filter(raw_data, event_name_split[0], values=[int(event_name_split[1])])
            on_off = evt_on_off.time_ranges
        except: 
            try:
                print("Problem getting on off event data by tdt.epoc_filter")
                print(event_name_split[0]) # some names end with _ and that gets replaced in data tank
                if event_name_split[0][-1] == "_":
                    adjusted_name = event_name_split[0][:-1]+"/"
                    print(adjusted_name)
                    evt_on_off= tdt.epoc_filter(raw_data, adjusted_name, values=[int(event_name_split[1])])
                    on_off = evt_on_off.time_ranges
            except:
                print("After readjusting event name still problem getting on off event data by tdt.epoc_filter")
    else: # special case for Cam1
        data = raw_data.epocs[event]
        onsets = []
        try:
            onsets = data.notes.ts
        except:
            onsets = []
        on_off[0] = onsets
    if len(on_off[0]) == 0:
        print("Could not get onsets any of known ways")
    return on_off


def get_single_channel(raw_data,chnl_name):
#    print('channel 1:', len(raw_data.streams[chnl_name].data))
    chnl_data = raw_data.streams[chnl_name].data
    # get that in seconds 
    times = create_timestamps(raw_data,chnl_name,chnl_data)
    # round to whole seconds
    full_total_sec = int(times[-1])
    # get only recording up to that time value
    full_times = [el for el in times if el < full_total_sec]
    chnl_data = chnl_data[:len(full_times)]
    return chnl_data

def create_timestamps(raw_data,chnl_name,channel_data):
    num_samples = len(channel_data)
    times = np.linspace(1, num_samples, num_samples) / raw_data.streams[chnl_name].fs
    return times

def get_frequency(raw_data,chnl_name):
    return raw_data.streams[chnl_name].fs

def trim_raw_data(raw_data,signal_name,control_name,beginning_sec,ending_sec):
    # return numpy array for all GCaMP channel raw data
    GCaMP_data = get_single_channel(raw_data,signal_name)
    # return numpy arrayfor all control raw data
    control_data = get_single_channel(raw_data,control_name)
    ts_raw = create_timestamps(raw_data,signal_name,GCaMP_data)
    #########
    # check
    # make sure both arrays are the same size
    len_all_GCaMP = GCaMP_data.size
    len_all_control = control_data.size
    # if GCaMP recording is longer, trim it to the control size
    if len_all_GCaMP > len_all_control:
        GCaMP_data = GCaMP_data[:len_all_control]
    # if control data is longer, trim it to GCaMP size
    elif len_all_GCaMP < len_all_control:
        control_data = control_data[:len_all_GCaMP]
        
#    print(len(GCaMP_data),len(control_data))
    # calculate the first beginning_sec seconds from a signal channel
    t0 = int(beginning_sec * raw_data.streams[signal_name].fs) # int rounds it to the nearest integer
    t1 = int(ending_sec * raw_data.streams[signal_name].fs)
#    print(t0,t1)
    if t1 == 0:
        ts = ts_raw[t0:]
        signal_trimmed = GCaMP_data[t0:]
        control_trimmed = control_data[t0:]
    else:
        ts = ts_raw[t0:-t1]
        signal_trimmed = GCaMP_data[t0:-t1]
        control_trimmed = control_data[t0:-t1]
#    print(len(signal_trimmed),len(control_trimmed))
    # create a dictionary
    trimmed_dict = {"ts":ts,"signal":signal_trimmed,"control":control_trimmed}
    return trimmed_dict


# create path to single data tank (data folder+subject name+experiment name)
# return a dictionary subject name: path
def create_list_of_paths(main_path,subject_list,experiment_name):
    paths_dict = {}
    for el in subject_list:
        to_subject = os.path.join(main_path,el)
        to_data_tank = os.path.join(to_subject,experiment_name)
        paths_dict[el] = to_data_tank
    return paths_dict

# create path to single data tank (data folder+experiment name+subject name)
# return a dictionary subject name: path
def create_list_of_paths_experiment_subjects(main_path,subject_list,experiment_name):
    paths_dict = {}
    # keep only subjects in the experiment folder
    valid_subjects = []
    to_experiment = os.path.join(main_path,experiment_name)
    for el in subject_list:
        to_data_tank = os.path.join(to_experiment,el)
        # check if path exists, because subjects come from all experiments
        if os.path.exists(to_data_tank):
            paths_dict[el] = to_data_tank
            valid_subjects.append(el) 
    return valid_subjects, paths_dict

# downsample by using median instead of mean
def downsample_tdt(signal_dict,downsample_n):
    GCaMP_data = signal_dict["signal"]
    control_data = signal_dict["control"]
    ts = signal_dict["ts"]
    GCaMP_data_means = []
    control_data_means = []
    for i in range(0, len(GCaMP_data), downsample_n):
        GCaMP_data_means.append(np.mean(GCaMP_data[i:i+downsample_n-1])) # This is the moving window median
        control_data_means.append(np.mean(control_data[i:i+downsample_n-1]))
    # adjust timestamps every nth
    ts_adjusted = ts[:-1:downsample_n]
    
    # check for equal lengths
    if len(GCaMP_data_means) < len(ts_adjusted):
        ts_adjusted = ts_adjusted[:len(GCaMP_data_means)]
    elif len(GCaMP_data_means) > len(ts_adjusted):
        GCaMP_data_means = GCaMP_data_means[:len(ts_adjusted)]
        control_data_means = control_data_means[:len(ts_adjusted)]
    
    return {"ts":ts_adjusted,"signal":GCaMP_data_means,"control":control_data_means}

# use 1D linear interpolation to downsample precisely in Hz
def downsample(signal_dict,target_Hz):
    GCaMP_data = signal_dict["signal"]
    control_data = signal_dict["control"]
    ts_original = signal_dict["ts"]
    ts = signal_dict["ts"]
    if ts[0] > 1: # if data was trimmed go back from zero
        ts = [el-ts[0] for el in ts]
    # round last time to full integer
    last = math.floor(ts[-1])
    ts_adjusted = np.linspace(1/target_Hz,last,last*target_Hz)
    GCaMP_interpolated = interp1d(ts, GCaMP_data)
    GCaMP_resampled_data = GCaMP_interpolated(ts_adjusted)
    control_interpolated = interp1d(ts, control_data)
    control_resampled_data = control_interpolated(ts_adjusted)
    # reconstruct time
    if ts_original[0] > 1: # if data was trimmed go back to start like the original
        ts_adjusted = [el+ts_original[0] for el in ts_adjusted]
    # print(ts_original[0],ts_adjusted[0])
    return {"ts":ts_adjusted,"signal":GCaMP_resampled_data,"control":control_resampled_data}

# modified polynomial
def normalize_dff(signal_dict,show_as,smooth,smooth_window):
    #####################
    # change lists to numpy array for calculations
    # create new timestamps for the data (from zero)
    ts_arr = np.asarray(signal_dict["ts"])
    signal_arr =np.asarray(signal_dict["signal"])
    control_arr = np.asarray(signal_dict["control"])
    if smooth == True:
        print("Start smoothing",smooth_window)
        a = 1
        b = np.divide(np.ones((smooth_window,)), smooth_window)
        control_arr = filtfilt(b, a, control_arr)
        signal_arr = filtfilt(b, a, signal_arr)
        print("Done smoothing")
        
    # Before filtering out +/- 2xstdev
#     # fit time axis to the 465nm stream  
# #    bls_Ca = np.polyfit(ts_arr,signal_arr,1) # deprecieted
#     bls_Ca = np.polynomial.polynomial.Polynomial.fit(ts_arr,signal_arr,1)
# #    print("bls_Ca",bls_Ca.convert().coef[::-1])
# #    F0Ca = np.polyval(bls_Ca,ts_arr) # deprecieted
#     F0Ca = polyval(ts_arr,bls_Ca.convert().coef)
#     # dF/F for the 465 channel
#     dFFCa = (signal_arr - F0Ca)/F0Ca *100
#     # fit time axis the 405nm stream
# #    bls_ref = np.polyfit(ts_arr,control_arr,1) # depreciated
#     bls_ref = np.polynomial.polynomial.Polynomial.fit(ts_arr,control_arr,1)
# #    F0Ref = np.polyval(bls_ref,ts_arr) # deprecieted
#     F0Ref = polyval(ts_arr,bls_ref.convert().coef)
#     # dF/F for the 405 channel
#     dFFRef = (control_arr - F0Ref)/F0Ref *100
# #    print(dFFRef)
#     dFFnorm = dFFCa - dFFRef

    ############################################################################################
    # 23/01 filter out signal values that are below or above 2 standard deviations from the signal mean 
    mean_signal = np.mean(signal_arr)
    stdev_signal = np.std(signal_arr)
    mean_control = np.mean(control_arr)
    stdev_control = np.std(control_arr)
    indexes_signal = np.where((signal_arr<mean_signal+2*stdev_signal) & (signal_arr>mean_signal-2*stdev_signal))
    selected_signal = signal_arr[indexes_signal]
    selected_ts_signal_arr = ts_arr[indexes_signal]
    indexes_control = np.where((control_arr<mean_control+2*stdev_control) & (control_arr>mean_control-2*stdev_control))
    selected_control = control_arr[indexes_control]
    selected_ts_control_arr = ts_arr[indexes_control]
    # fit time axis to the 465nm stream  
    bls_Ca = np.polynomial.polynomial.Polynomial.fit(selected_ts_signal_arr,selected_signal,1)
    bls_ref = np.polynomial.polynomial.Polynomial.fit(selected_ts_control_arr,selected_control,1)
    F0Ca = polyval(ts_arr,bls_Ca.convert().coef)
    # dF/F for the 465 channel
    dFFCa = (signal_arr - F0Ca)/F0Ca *100
    F0Ref = polyval(ts_arr,bls_ref.convert().coef)
    # dF/F for the 405 channel
    dFFRef = (control_arr - F0Ref)/F0Ref *100
    #    print(dFFRef)
    dFFnorm = dFFCa - dFFRef
    ######################################################################################################

    # find all values of the normalized DF/F that are negative so you can next shift up the curve 
    # to make 0 the mean value for DF/F
    negative = dFFnorm[dFFnorm<0]
    dFF=dFFnorm-np.mean(negative)

    if show_as == "Z-Score":
        median_all = np.median(dFF)
        mad = stats.median_abs_deviation(dFF)
        dFF = (dFF - median_all)/mad

    return {"ts":ts_arr,"normalized_signal":dFF}

# modified polynomial that uses custom baseline
def normalize_dff_baseline(signal_dict,baseline_dict,show_as,smooth,smooth_window):
    #####################
    # change lists to numpy array for calculations
    # create new timestamps for the data (from zero)
    ts_arr = np.asarray(signal_dict["ts"])
    signal_arr =np.asarray(signal_dict["signal"])
    control_arr = np.asarray(signal_dict["control"])
    print(f"All trace: first: {ts_arr[0]}, last: {ts_arr[-1]}")
    # change baseline lists to numpy array for calculations
    ts_baseline_arr = np.asarray(baseline_dict["ts"])
    signal_baseline_arr =np.asarray(baseline_dict["signal"])
    control_baseline_arr = np.asarray(baseline_dict["control"])
    print(f"Baseline: first: {ts_baseline_arr[0]}, last: {ts_baseline_arr[-1]}")
    if smooth == True:
        print("Start smoothing",smooth_window)
        a = 1
        b = np.divide(np.ones((smooth_window,)), smooth_window)
        control_arr = filtfilt(b, a, control_arr)
        signal_arr = filtfilt(b, a, signal_arr)
        control_baseline_arr = filtfilt(b, a, control_baseline_arr)
        signal_baseline_arr = filtfilt(b, a, signal_baseline_arr)
        print("Done smoothing")
    ############################################################################################
    # 23/01 filter out signal values that are below or above 2 standard deviations from the signal mean 
    mean_baseline_signal = np.mean(signal_baseline_arr)
    stdev_baseline_signal = np.std(signal_baseline_arr)
    mean_baseline_control = np.mean(control_baseline_arr)
    stdev_baseline_control = np.std(control_baseline_arr)
    indexes_signal = np.where((signal_baseline_arr<mean_baseline_signal+2*stdev_baseline_signal) & (signal_baseline_arr>mean_baseline_signal-2*stdev_baseline_signal))
    selected_signal = signal_baseline_arr[indexes_signal]
    selected_ts_baseline_signal_arr = ts_baseline_arr[indexes_signal]
    indexes_control = np.where((control_baseline_arr<mean_baseline_control+2*stdev_baseline_control) & (control_baseline_arr>mean_baseline_control-2*stdev_baseline_control))
    selected_control = control_baseline_arr[indexes_control]
    selected_ts_baseline_control_arr = ts_baseline_arr[indexes_control]
    # fit time axis to the 465nm stream  
    bls_Ca = np.polynomial.polynomial.Polynomial.fit(selected_ts_baseline_signal_arr,selected_signal,1)
    bls_ref = np.polynomial.polynomial.Polynomial.fit(selected_ts_baseline_control_arr,selected_control,1)
    F0Ca = polyval(ts_arr,bls_Ca.convert().coef)
    # dF/F for the 465 channel
    dFFCa = (signal_arr - F0Ca)/F0Ca *100
    F0Ref = polyval(ts_arr,bls_ref.convert().coef)
    # dF/F for the 405 channel
    dFFRef = (control_arr - F0Ref)/F0Ref *100
    #    print(dFFRef)
    dFFnorm = dFFCa - dFFRef
    ######################################################################################################
#     # fit time axis to the 465nm stream  
#     bls_Ca = np.polynomial.polynomial.Polynomial.fit(ts_baseline_arr,signal_baseline_arr,1)
# #    print("bls_Ca",bls_Ca.convert().coef[::-1])
#     F0Ca = polyval(ts_arr,bls_Ca.convert().coef)
#     # dF/F for the 465 channel
#     dFFCa = (signal_arr - F0Ca)/F0Ca *100
#     # fit time axis the 405nm stream
#     bls_ref = np.polynomial.polynomial.Polynomial.fit(ts_baseline_arr,control_baseline_arr,1)
#     F0Ref = polyval(ts_arr,bls_ref.convert().coef)
#     # dF/F for the 405 channel
#     dFFRef = (control_arr - F0Ref)/F0Ref *100
# #    print(dFFRef)
#     dFFnorm = dFFCa - dFFRef
#     # find all values of the normalized DF/F that are negative so you can next shift up the curve 
#     # to make 0 the mean value for DF/F
#     # negative = dFFnorm[dFFnorm<0]
#     # dFF=dFFnorm-np.mean(negative)
    ###########################################
    # fix the shift towards zero
    # get fragment of length baseline
    # print(f"Baseline signal mean: {mean_baseline_signal}, Baseline control mean: {mean_baseline_control}")
    # print(f"Median before shift: {np.median(dFFnorm)}, MAD before shit: {stats.median_abs_deviation(dFFnorm)}")
    baseline_normalized = dFFnorm[:len(ts_baseline_arr)]
    negative_baseline = baseline_normalized[baseline_normalized<0]
    dFF = dFFnorm - np.mean(negative_baseline)
    ###########################################

    if show_as == "Z-Score":
        median_all = np.median(dFF)
        mad = stats.median_abs_deviation(dFF)
        dFF = (dFF - median_all)/mad
        # print(f"Median: {median_all}, MAD: {mad}")

    return {"ts":ts_arr,"normalized_signal":dFF}


# https://github.com/djamesbarker/pMAT
def normalize_pMat(signal_dict,show_as,smooth,smooth_window):
    # change lists to numpy array for calculations
    ts_arr = np.asarray(signal_dict["ts"])
    signal_arr =np.asarray(signal_dict["signal"])
    control_arr = np.asarray(signal_dict["control"])
    if smooth == True:
        print("Start smoothing",smooth_window)
        a = 1
        b = np.divide(np.ones((smooth_window,)), smooth_window)
        control_arr = filtfilt(b, a, control_arr)
        signal_arr = filtfilt(b, a, signal_arr)
        print("Done smoothing")
        
        
    # Before
#    bls = np.polyfit(control_arr, signal_arr, 1)
#    fit_line = np.multiply(bls[0], control_arr) + bls[1]
#    dff = (signal_arr - fit_line)/fit_line * 100
        
    # Before filtering out +/- 2xstdev
#     # https://stackoverflow.com/questions/45338872/matlab-polyval-function-with-three-outputs-equivalent-in-python-numpy
#     mu = np.mean(control_arr)
#     std = np.std(control_arr, ddof=0)
#     # Call np.polyfit(), using the shifted and scaled version of control_arr
# #    cscaled = np.polyfit((control_arr - mu)/std, signal_arr, 1) # depreciated
#     cscaled = np.polynomial.polynomial.Polynomial.fit((control_arr - mu)/std, signal_arr, 1)
#     # Create a poly1d object that can be called
# #    pscaled = np.poly1d(cscaled) # old obsolete function
#     # https://numpy.org/doc/stable/reference/routines.polynomials.html
#     pscaled = Polynomial(cscaled.convert().coef)
#     # Inputs to pscaled must be shifted and scaled using mu and std
#     F0 = pscaled((control_arr - mu)/std)
# #    print("F0?",F0[:20])
#     dffnorm = (signal_arr - F0)/F0 * 100

    #############################################################################################
    # filter out signal values that are below or above 2 standard deviations from the signal mean 
    mean_signal = np.mean(signal_arr)
    stdev_signal = np.std(signal_arr)
    indexes = np.where((signal_arr<mean_signal+2*stdev_signal) & (signal_arr>mean_signal-2*stdev_signal))
    selected_signal = signal_arr[indexes]
    selected_control = control_arr[indexes]
    mu = np.mean(selected_control)
    std = np.std(selected_control)
    cscaled = np.polynomial.polynomial.Polynomial.fit((selected_control - mu)/std, selected_signal, 1)
    pscaled = Polynomial(cscaled.convert().coef)
    F0 = pscaled((control_arr - mu)/std)
    dffnorm = (signal_arr - F0)/F0 * 100
    #################################################################################################
    # find all values of the normalized DF/F that are negative so you can next shift up the curve 
    # to make 0 the mean value for DF/F
    negative = dffnorm[dffnorm<0]
    dff=dffnorm-np.mean(negative)

    if show_as == "Z-Score":
        median_all = np.median(dff)
        mad = stats.median_abs_deviation(dff)
        dff = (dff - median_all)/mad

    return {"ts":ts_arr,"normalized_signal":dff}

# function uses a custom baseline range to fit all of the data later
def normalize_pMat_custom_baseline(signal_dict,baseline_dict,show_as,smooth,smooth_window):
    # change lists to numpy array for calculations
    ts_arr = np.asarray(signal_dict["ts"])
    signal_arr =np.asarray(signal_dict["signal"])
    control_arr = np.asarray(signal_dict["control"])
    print(f"All trace: first: {ts_arr[0]}, last: {ts_arr[-1]}")
    # change baseline lists to numpy array for calculations
    ts_baseline_arr = np.asarray(baseline_dict["ts"])
    signal_baseline_arr =np.asarray(baseline_dict["signal"])
    control_baseline_arr = np.asarray(baseline_dict["control"])
    print(f"Baseline: first: {ts_baseline_arr[0]}, last: {ts_baseline_arr[-1]}")
    if smooth == True:
        print("Start smoothing",smooth_window)
        a = 1
        b = np.divide(np.ones((smooth_window,)), smooth_window)
        control_arr = filtfilt(b, a, control_arr)
        signal_arr = filtfilt(b, a, signal_arr)
        control_baseline_arr = filtfilt(b, a, control_baseline_arr)
        signal_baseline_arr = filtfilt(b, a, signal_baseline_arr)
        print("Done smoothing")
        
  
    #############################################################################################
    # filter out signal values that are below or above 2 standard deviations from the signal mean  
    mean_baseline_signal = np.mean(signal_baseline_arr)
    stdev_baseline_signal = np.std(signal_baseline_arr)
    indexes = np.where((signal_baseline_arr<mean_baseline_signal+2*stdev_baseline_signal) & (signal_baseline_arr>mean_baseline_signal-2*stdev_baseline_signal))
    selected_signal = signal_baseline_arr[indexes]
    selected_control = control_baseline_arr[indexes]
    mu = np.mean(selected_control)
    std = np.std(selected_control)
    cscaled = np.polynomial.polynomial.Polynomial.fit((selected_control - mu)/std, selected_signal, 1)
    pscaled = Polynomial(cscaled.convert().coef)
    F0 = pscaled((control_arr - mu)/std)
    dffnorm = (signal_arr - F0)/F0 * 100
    ######################################################################################################
    # # https://stackoverflow.com/questions/45338872/matlab-polyval-function-with-three-outputs-equivalent-in-python-numpy
    # mu = np.mean(control_baseline_arr)
    # std = np.std(control_baseline_arr, ddof=0)
    # # Call np.polyfit(), using the shifted and scaled version of control_arr
    # cscaled = np.polynomial.polynomial.Polynomial.fit((control_baseline_arr - mu)/std, signal_baseline_arr, 1)
    # # Create a poly1d object that can be called
    # # https://numpy.org/doc/stable/reference/routines.polynomials.html
    # pscaled = Polynomial(cscaled.convert().coef)
    # # Inputs to pscaled must be shifted and scaled using mu and std
    # F0 = pscaled((control_arr - mu)/std)
    # dffnorm = (signal_arr - F0)/F0 * 100
    # # find all values of the normalized DF/F that are negative so you can next shift up the curve 
    # # to make 0 the mean value for DF/F
    # # negative = dffnorm[dffnorm<0]
    # # dff=dffnorm-np.mean(negative)

    ###########################################
    # fix the shift towards zero
    # get fragment of length baseline
    baseline_normalized = dffnorm[:len(ts_baseline_arr)]
    negative_baseline = baseline_normalized[baseline_normalized<0]
    dff = dffnorm - np.mean(negative_baseline)
    ###########################################

    if show_as == "Z-Score":
        median_all = np.median(dff)
        mad = stats.median_abs_deviation(dff)
        dff = (dff - median_all)/mad

    return {"ts":ts_arr,"normalized_signal":dff}

# modified tdt function to handle special Cam1 case with notes
def my_epoc_filter(data, epoc, *, values=None, t=None, tref=False, keepdata=True):
    """TDT tank data filter. Extract data around epoc events.
    data = epoc_filter(data, epoc) where data is the output of read_block,
    epoc is the name of the epoc to filter on, and parameter value pairs
    define the filtering conditions.
    
    If no parameters are specified, then the time range of the epoc event
    is used as a time filter.
    
    Also creates data.filter, a string that describes the filter applied.
    Optional keyword arguments:
        values      specify array of allowed values
                      ex: tempdata = epoc_filter(data, 'Freq', values=[9000, 10000])
                        > retrieves data when Freq = 9000 or Freq = 10000
        t           specify onset/offset pairs relative to epoc onsets. If the
                      offset is not provided, the epoc offset is used.
                      ex: tempdata = epoc_filter(data, 'Freq', t=[-0.1, 0.5])
                        > retrieves data from 0.1 seconds before Freq onset to 0.4
                          seconds after Freq onset. Negative time ranges are discarded.
        tref        use the epoc event onset as a time reference. All timestamps for
                      epoc, snippet, and scalar events are then relative to epoc onsets.
                      ex: tempdata = epoc_filter(data, 'Freq', tref=True)
                        > sets snippet timestamps relative to Freq onset
        keepdata    keep the original stream data array and add a field called
                      'filtered' that holds the data from each valid time range. 
                      Defaults to True.
    
    IMPORTANT! Use a time filter (t argument) only after all value filters have been set.
    """
    
    data = copy.deepcopy(data)
    
    filter_string = ''

    if not hasattr(data, 'epocs'):
        raise Exception('no epocs found')
    elif len(data.epocs.keys()) == 0:
        raise Exception('no epocs found')

    fff = data.epocs.keys()
    match = ''
    all_names = []
    for k in data.epocs.keys():
        all_names.append(data.epocs[k].name)
        if data.epocs[k].name == epoc:
            match = k

    if not hasattr(data.epocs, match):
        raise Exception('epoc {0} is not a valid epoc event, valid events are: {1}'.format(epoc, all_names))

    if t:
        try:
            if len(t) > 2:
                raise Exception('{0} t vector must have 1 or 2 elements only'.format(repr(t)))
        except:
            raise Exception('{0} t vector must have 1 or 2 elements only'.format(repr(t)))
    
    ddd = data.epocs[match]

    time_ranges = None
    
    # VALUE FILTER, only use time ranges where epoc value is in filter array
    valid = [False for el in ddd.data]
    if values:
        # find valid time ranges
        # use this for numbers
        for i in range(len(values)):
            for j in range(len(ddd.data)):
                if ddd.data[j] == values[i]:
                    valid[j] = True
                    break
        
        time_ranges = np.vstack((ddd.onset[valid], ddd.offset[valid]))
        # print(ddd.onset[valid],ddd.offset[valid])
        # print(f"time ranges: {time_ranges}")
        # create filter string
        filter_string = '{0}:VALUE in [{1}];'.format(epoc, ','.join([str(v) for v in values]))

        # AND time_range with existing time ranges
        if hasattr(data, 'time_ranges'):
            time_ranges = combine_time(time_ranges, data.time_ranges, 'AND')
            data.time_ranges = time_ranges
    
    #TIME FILTER
    t1 = 0
    t2 = np.nan
    if t:
        t1 = t[0]
        if len(t) == 2:
            t2 = t[1]
        
        if not isinstance(time_ranges, np.ndarray):
            # preallocate
            time_ranges = np.zeros((2, len(ddd.onset)))
            for j in range(len(ddd.onset)):
                time_ranges[:, j] = [ddd.onset[j], ddd.offset[j]]
        else:
            time_ranges = data.time_ranges
        
        # find valid time ranges
        for j in range(time_ranges.shape[1]): #= 1:size(time_ranges,2)
            if np.isnan(t2):
                time_ranges[:, j] = [time_ranges[0,j]+t1, time_ranges[1,j]]
            else:
                time_ranges[:, j] = [time_ranges[0,j]+t1, time_ranges[0,j]+t1+t2]
        
        # throw away negative time ranges
        if np.all(~np.isnan(time_ranges)):
            time_ranges = time_ranges[:,time_ranges[0,:]>0]
        
        # create filter string
        if np.isnan(t2):
            filter_string = 'TIME:{0} [{1}];'.format(epoc, np.round(t1,2))
        else:
            filter_string = 'TIME:{0} [{1}:{2}];'.format(epoc, np.round(t1,2), np.round(t2,2))
        data.time_ranges = time_ranges
        data.time_ref = [t1, t2]
    
    # TREF FILTER
    if tref:
        filter_string += 'REF:{0}'.format(epoc)
        if hasattr(tref, "__len__"):
            if len(tref) > 1:
                t1 = tref[0]
    
    if values or t or tref:
        pass
    else:
        # no filter specified, use EPOC time ranges
        time_ranges = np.vstack((ddd.onset, ddd.offset))
        # AND time_range with existing time ranges
        if hasattr(data, 'time_ranges'):
            time_ranges = combine_time(time_ranges, data.time_ranges, 'AND')
            data.time_ranges = time_ranges
        filter_string = 'EPOC:{0};'.format(epoc)
        
    # set filter string
    if hasattr(data, 'filter'):
        data.filter += filter_string
    else:
        data.filter = filter_string
               
    time_ranges = data.time_ranges
    # print(f"data.time_ranges: {time_ranges}")
    # FILTER ALL EXISTING DATA ON THESE TIME RANGES
    # filter streams
    if data.streams:
        for k in data.streams.keys():
            fs = data.streams[k].fs
            sf = 1/(2.56e-6*fs)
            td_sample = np.uint64(data.streams[k].start_time/2.56e-6)
            filtered = []
            max_ind = max(data.streams[k].data.shape)
            for j in range(time_ranges.shape[1]):
                tlo_sample = np.uint64(time_ranges[0,j]/2.56e-6)
                thi_sample = np.uint64(time_ranges[1,j]/2.56e-6)
                
                onset = np.uint64(np.maximum(np.round((tlo_sample - td_sample)/sf),0))
                # throw it away if onset or offset extends beyond recording window
                if np.isinf(time_ranges[1,j]):
                    if onset <= max_ind and onset > -1:
                        if data.streams[k].data.ndim == 1:
                            filtered.append(data.streams[k].data[onset:])
                        else:
                            filtered.append(data.streams[k].data[:,onset:])
                        break
                else:
                    offset = np.uint64(np.maximum(np.round((thi_sample-td_sample)/sf),0))
                    if offset <= max_ind and offset > -1 and onset <= max_ind and onset > -1:
                        if data.streams[k].data.ndim == 1:
                            filtered.append(data.streams[k].data[onset:offset])
                        else:
                            filtered.append(data.streams[k].data[:,onset:offset])
            if keepdata:
                data.streams[k].filtered = filtered
            else:
                data.streams[k].data = filtered
    
    # filter snips
    if data.snips:
        warning_value = -1
        for k in data.snips.keys():
            ts = data.snips[k].ts
            
            # preallocate
            keep = np.zeros(len(ts), dtype=np.uint64)
            diffs = np.zeros(len(ts)) # for relative timestamps
            keep_ind = 0
            
            for j in range(len(ts)):
                ind1 = ts[j] > time_ranges[0,:]
                ind2 = ts[j] < time_ranges[1,:]
                ts_ind = np.where(ind1 & ind2)[0]
                if len(ts_ind) > 0:
                    if len(ts_ind) > 1:
                        min_diff = np.min(np.abs(time_ranges[1, ts_ind[0]]-time_ranges[1, ts_ind[1]]), np.abs(time_ranges[1, ts_ind[0]]-time_ranges[1, ts_ind[1]]))
                        warning_value = min_diff
                        continue
                    keep[keep_ind] = np.uint64(j)
                    diffs[keep_ind] = ts[j] - time_ranges[0, ts_ind] + t1; # relative ts
                    keep_ind += 1
            
            if warning_value > 0:
                warnings.warn('time range overlap, consider a maximum time range of %.2fs'.format(warning_value), Warning)
            
            # truncate
            keep = keep[:keep_ind]
            diffs = diffs[:keep_ind]
            if hasattr(data.snips[k], 'data'):
                if len(data.snips[k].data) > 0:
                    data.snips[k].data = data.snips[k].data[keep]
            
            if tref:
                data.snips[k].ts = diffs
            else:
                data.snips[k].ts = data.snips[k].ts[keep]
            
            # if there are any extra fields, keep those
            for kk in data.snips[k].keys():
                if kk in ['ts', 'name', 'data', 'sortname', 'fs', 'code', 'size', 'type', 'dform']:
                    continue
                if len(data.snips[k][kk]) >= np.max(keep):
                    data.snips[k][kk] = data.snips[k][kk][keep]
    
    # filter scalars, include if timestamp falls in valid time range
    if data.scalars:
        for k in data.scalars.keys():
            ts = data.scalars[k].ts
            keep = get_valid_ind(ts, time_ranges)
            if len(keep) > 0:
                # scalars can have multiple rows
                if data.scalars[k].data.ndim > 1:
                    data.scalars[k].data = data.scalars[k].data[:,keep]
                else:
                    data.scalars[k].data = data.scalars[k].data[keep]
                data.scalars[k].ts = data.scalars[k].ts[keep]
            else:
                data.scalars[k].data = []
                data.scalars[k].ts = []

    # filter epocs, include if onset falls in valid time range
    if data.epocs:
        for k in data.epocs.keys():
            ts = data.epocs[k].onset
            keep = get_valid_ind(ts, time_ranges)
            if len(keep) > 0:
                data.epocs[k].data = data.epocs[k].data[keep]
                data.epocs[k].onset = data.epocs[k].onset[keep]
                if hasattr(data.epocs[k], 'notes'):
                    if hasattr(data.epocs[k].notes, 'ts'):
                        keep2 = get_valid_ind(data.epocs[k].notes.ts, time_ranges)
                        data.epocs[k].notes.ts = data.epocs[k].notes.ts[keep2]
                        data.epocs[k].notes.index = data.epocs[k].notes.index[keep2]
                        data.epocs[k].notes.notes = data.epocs[k].notes.notes[keep2]
                    else:
                        data.epocs[k].notes = data.epocs[k].notes[keep]
                if hasattr(data.epocs[k], 'offset'):
                    data.epocs[k].offset = data.epocs[k].offset[keep]
            else:
                data.epocs[k].data = []
                data.epocs[k].onset = []
                if hasattr(data.epocs[k], 'offset'):
                    data.epocs[k].offset = []
                if hasattr(data.epocs[k], 'notes'):
                    data.epocs[k].notes = []
    
    return data

#https://www.tdt.com/support/python-sdk/offline-analysis-examples/fiber-photometry-epoch-averaging-example/
def filter_data_around_event(raw_data,perievent_options_dict,settings_dict,signal_name,control_name):
    event_name_split = perievent_options_dict["event"].split(" ")
    before = -perievent_options_dict["sec_before"]
    till = perievent_options_dict["sec_before"]+perievent_options_dict["sec_after"]+0.1
    trange = [before,till]
    if event_name_split[0] == 'Cam1': # special case
        # find values = indexes of timestamps
        my_values = []
        data = raw_data.epocs[event_name_split[0]]
        my_onsets = []
        try:
            my_onsets = data.notes.ts
        except:
            my_onsets = []
        if len(my_onsets) > 0:
            all_values = data.data
            all_ts = data.onset
            indexes = []
            for el in my_onsets:
                try:
                    idx = np.where(all_ts == el)
                    indexes.append(int(idx[0][0]))
                except:
                    print("Onset not found")
            if len(indexes) > 0:
                for el in indexes:
                    my_values.append(all_values[el])
                print(f"my_values {my_values}")
                modified_data= my_epoc_filter(raw_data, event_name_split[0], values = my_values, t=trange, tref=True)
    try:    # some data was not readable by tdt.epoc_filter
        modified_data= tdt.epoc_filter(raw_data, event_name_split[0], t=trange, values=[int(event_name_split[1])], tref=True)
    except: 
        print("Problem getting on off event data by tdt.epoc_filter")
        print(event_name_split[0]) # some names end with _ and that gets replaced in data tank
        if event_name_split[0][-1] == "_":
            adjusted_name = event_name_split[0][:-1]+"/"
            print(adjusted_name)
            modified_data= tdt.epoc_filter(raw_data, adjusted_name, t=trange, values=[int(event_name_split[1])], tref=True)
    all_control_filtered = modified_data.streams[control_name].filtered
    all_signal_filtered = modified_data.streams[signal_name].filtered
    modified_data.streams[signal_name].filtered = all_signal_filtered
    modified_data.streams[control_name].filtered = all_control_filtered
    '''Applying a time filter to a uniformly sampled signal means that the length 
        of each segment could vary by one sample. Let's find the minimum length 
        so we can trim the excess off before calculating the median.
        '''
    if len(all_signal_filtered) > 0: # sometimes the first and only event gave empty list after filtering
        min1 = np.min([np.size(x) for x in modified_data.streams[signal_name].filtered])
        min2 = np.min([np.size(x) for x in modified_data.streams[control_name].filtered])
        modified_data.streams[signal_name].filtered = [x[1:min1] for x in modified_data.streams[signal_name].filtered]
        modified_data.streams[control_name].filtered = [x[1:min2] for x in modified_data.streams[control_name].filtered]
    else:
        print("Not enough data that satisfy your request. Try changing event or times around event.")
        return modified_data
    # downsample data as well
    N = settings_dict[0]["downsample"] 
    
    control = []
    signal = []
    fs = get_frequency(raw_data,signal_name)
    try:
        for lst in modified_data.streams[control_name].filtered:
            ts = np.linspace(1, len(lst), len(lst))/fs
            # round last time to full integer
            last = math.floor(ts[-1])
            # make new times till total expected
            ts_adjusted = np.linspace(1/N,last,last*N)
            control_interpolated = interp1d(ts, lst)
            resampled_data = control_interpolated(ts_adjusted)
            control.append(resampled_data)
        for lst in modified_data.streams[signal_name].filtered:
            ts = np.linspace(1, len(lst), len(lst))/fs
            # round last time to full integer
            last = math.floor(ts[-1])
            ts_adjusted = np.linspace(1/N,last,last*N)
            signal_interpolated = interp1d(ts, lst)
            resampled_data = signal_interpolated(ts_adjusted)
            signal.append(resampled_data)
        modified_data.streams[signal_name].filtered_downsampled = signal
        modified_data.streams[control_name].filtered_downsampled = control   
    except:
        print("Not enough data that satisfy your request. Try changing event or times around event.")
    
    return modified_data

def analyze_perievent_data(data,current_trials,perievent_options_dict,settings_dict,signal_name,control_name):
    # create a dictionary with analysed data for plotting
    analyzed_perievent_dict = {}
    # get downsampled data from around event 
    if len(signal_name) > 0 and len(control_name) > 0:
        GCaMP_perievent_data = data.streams[signal_name].filtered_downsampled
        control_perievent_data = data.streams[control_name].filtered_downsampled
    else:
        GCaMP_perievent_data = data["signal"]
        control_perievent_data = data["control"]

    # only selected trials
    if len(current_trials) > 0:
        GCaMP_perievent_data = [GCaMP_perievent_data[trial-1] for trial in current_trials]
        control_perievent_data = [control_perievent_data[trial-1] for trial in current_trials]

    
    # Create a mean signal, standard error(median absolute deviation in pMat) of signal, and DC offset
    mean_signal = np.mean(GCaMP_perievent_data, axis=0)
    std_signal = np.std(GCaMP_perievent_data, axis=0) / np.sqrt(len(GCaMP_perievent_data))
    dc_signal = np.mean(mean_signal)
    mean_control = np.mean(control_perievent_data, axis=0)
    std_control = np.std(control_perievent_data, axis=0) / np.sqrt(len(control_perievent_data))
    dc_control = np.mean(mean_control)
    
    # Create the time vector for each stream store
    N = settings_dict[0]["downsample"]
    tFrom0 = np.linspace(1,len(mean_signal),len(mean_signal))/N
    ts_signal4average = -perievent_options_dict["sec_before"] + tFrom0
    ts_control4average = ts_signal4average
    
    # Subtract DC offset to get signals on top of one another
    mean_signal = mean_signal - dc_signal
    mean_control = mean_control - dc_control
    
    # create dictionary with data for average plot
    avg_data2plot_dict = {}
    avg_data2plot_dict["ts_signal"] = ts_signal4average
    avg_data2plot_dict["ts_control"] = ts_control4average
    avg_data2plot_dict["signal"] = mean_signal
    avg_data2plot_dict["control"] = mean_control
    avg_data2plot_dict["std_signal"] = std_signal
    avg_data2plot_dict["std_control"] = std_control
    # add that dict to all analysed data dict
    analyzed_perievent_dict["average"] = avg_data2plot_dict
    
    if settings_dict[0]["filter"] == True:
        print("Start smoothing",settings_dict[0]["filter_window"])
        # smooth data using lowess filter
        temp_signal_smoothed = []
        temp_control_smoothed = []
        a = 1
        b = np.divide(np.ones((settings_dict[0]["filter_window"],)), settings_dict[0]["filter_window"])
        for i in range(len(GCaMP_perievent_data)):
            single_evt = filtfilt(b, a, GCaMP_perievent_data[i])
            temp_signal_smoothed.append(single_evt)
        for i in range(len(control_perievent_data)):
            single_evt = filtfilt(b, a, control_perievent_data[i])
            temp_control_smoothed.append(single_evt)
        if len(temp_signal_smoothed) > 0 and len(temp_control_smoothed) > 0:
            GCaMP_perievent_data = temp_signal_smoothed
            control_perievent_data = temp_control_smoothed
        print("Done smoothing")
        
    # normalize
    y_dff_all = []
    # find out how to normalize
    if settings_dict[0]["normalization"] == 'Standard Polynomial Fitting':
        # https://github.com/djamesbarker/pMAT
        for x, y in zip(control_perievent_data, GCaMP_perievent_data):
            x = np.array(x)
            y = np.array(y)           
            ########################################################################################################################
            # Before removing +/-2xstdev
        #     # https://stackoverflow.com/questions/45338872/matlab-polyval-function-with-three-outputs-equivalent-in-python-numpy
        #     mu = np.mean(x)
        #     std = np.std(x, ddof=0)
        #     # Call np.polyfit(), using the shifted and scaled version of control_arr
        #     cscaled = np.polynomial.polynomial.Polynomial.fit((x - mu)/std, y, 1)
        #     # Create a poly1d object that can be called
        #     # https://numpy.org/doc/stable/reference/routines.polynomials.html
        #     pscaled = Polynomial(cscaled.convert().coef)
        #     # Inputs to pscaled must be shifted and scaled using mu and std
        #     F0 = pscaled((x - mu)/std)
        # #    print("F0?",F0[:20])
        #     dffnorm = (y - F0)/F0 * 100
            #############################################################################################
            # filter out signal values that are below or above 2 standard deviations from the signal mean 
            mean_signal = np.mean(y)
            stdev_signal = np.std(y)
            indexes = np.where((y<mean_signal+2*stdev_signal) & (y>mean_signal-2*stdev_signal))
            selected_signal = y[indexes]
            selected_control = x[indexes]
            mu = np.mean(selected_control)
            std = np.std(selected_control)
            cscaled = np.polynomial.polynomial.Polynomial.fit((selected_control - mu)/std, selected_signal, 1)
            pscaled = Polynomial(cscaled.convert().coef)
            F0 = pscaled((x - mu)/std)
            dffnorm = (y - F0)/F0 * 100
            #############################################################################
            # find all values of the normalized DF/F that are negative so you can next shift up the curve 
            # to make 0 the mean value for DF/F
            negative = dffnorm[dffnorm<0]
            dff=dffnorm-np.mean(negative)
            y_dff_all.append(dff)
            
    elif settings_dict[0]["normalization"] == 'Modified Polynomial Fitting':
        for x, y in zip(control_perievent_data, GCaMP_perievent_data):
            x = np.array(x)
            y = np.array(y)
            ##################################################################################
            # Before removing +/-2xstdev
            # bls_signal = np.polynomial.polynomial.Polynomial.fit(ts_signal4average, y, 1)
            # F0_signal = polyval(ts_signal4average,bls_signal.convert().coef)
            # dFF_signal = (y - F0_signal)/F0_signal *100
            
            # bls_control = np.polynomial.polynomial.Polynomial.fit(ts_control4average,x,1)
            # F0_control = polyval(ts_control4average,bls_control.convert().coef)
            # dFF_control = (x - F0_control)/F0_control *100
            # dFFnorm = dFF_signal - dFF_control
            ############################################################################################
            # 23/01 filter out signal values that are below or above 2 standard deviations from the signal mean 
            mean_signal = np.mean(y)
            stdev_signal = np.std(y)
            mean_control = np.mean(x)
            stdev_control = np.std(x)
            indexes_signal = np.where((y<mean_signal+2*stdev_signal) & (y>mean_signal-2*stdev_signal))
            selected_signal = y[indexes_signal]
            selected_ts_signal_arr = ts_signal4average[indexes_signal]
            indexes_control = np.where((x<mean_control+2*stdev_control) & (x>mean_control-2*stdev_control))
            selected_control = x[indexes_control]
            selected_ts_control_arr = ts_control4average[indexes_control]
            # fit time axis to the 465nm stream  
            bls_signal = np.polynomial.polynomial.Polynomial.fit(selected_ts_signal_arr,selected_signal,1)
            bls_control = np.polynomial.polynomial.Polynomial.fit(selected_ts_control_arr,selected_control,1)
            F0_signal = polyval(ts_signal4average,bls_signal.convert().coef)
            F0_control = polyval(ts_signal4average,bls_control.convert().coef)
            # dF/F for the 465 channel
            dFF_signal = (y - F0_signal)/F0_signal *100
            # dF/F for the 405 channel
            dFF_control = (x - F0_control)/F0_control *100
            dFFnorm = dFF_signal - dFF_control
            #############################################################################################
            # find all values of the normalized DF/F that are negative so you can next shift up the curve 
            # to make 0 the mean value for DF/F
            negative = dFFnorm[dFFnorm<0]
            dFF = dFFnorm-np.mean(negative)
            y_dff_all.append(dFF)
    
                   
    # get the z-score and standard error(median absolute deviation in pMat)
    zscore_all = []
    for dF in y_dff_all:
       ind = np.where((np.array(ts_control4average)<perievent_options_dict["baseline_to"]) & (np.array(ts_control4average)>perievent_options_dict["baseline_from"]))
       zb = np.median(dF[ind])
       mad = stats.median_abs_deviation(dF[ind])
       zscore_all.append((dF - zb)/mad)
    zerror = np.std(zscore_all, axis=0)/np.sqrt(np.size(zscore_all, axis=0))
    
#    for dF in y_dff_all:
#        ind = np.where((np.array(ts_control4average)<perievent_options_dict["baseline_to"]) & (np.array(ts_control4average)>perievent_options_dict["baseline_from"]))
#        zb = np.mean(dF[ind])
#        zsd = np.std(dF[ind])
#        zscore_all.append((dF - zb)/zsd)
#
#    zerror = np.std(zscore_all, axis=0)/np.sqrt(np.size(zscore_all, axis=0))
    
    # create dictionary with data for z-score plot
    zscore_data2plot_dict = {}
    zscore_data2plot_dict["ts"] = ts_control4average
    zscore_data2plot_dict["zscored"] = zscore_all
    zscore_data2plot_dict["zerror"] = zerror
    # add that dict to all analysed data dict
    analyzed_perievent_dict["zscore"] = zscore_data2plot_dict
    
    # Quantify changes as an area under the curve
    pre_ind = np.where((np.array(ts_control4average)<perievent_options_dict["auc_pre_to"]) & (np.array(ts_control4average)>perievent_options_dict["auc_pre_from"]))
    AUC_pre = auc(ts_control4average[pre_ind], np.mean(zscore_all, axis=0)[pre_ind])
    post_ind = np.where((np.array(ts_control4average)>perievent_options_dict["auc_post_from"]) & (np.array(ts_control4average)<perievent_options_dict["auc_post_to"]))
    AUC_post= auc(ts_control4average[post_ind], np.mean(zscore_all, axis=0)[post_ind])
    AUC = [AUC_pre, AUC_post]  
    # run a two-sample T-test
    t_stat,p_val = stats.ttest_ind(np.mean(zscore_all, axis=0)[pre_ind],
                               np.mean(zscore_all, axis=0)[post_ind], equal_var=False)

    aucs_pre_by_trial =[]
    aucs_post_by_trial = []
    # get single trial aucs
    for i in range(len(zscore_all)):
        pre = auc(ts_control4average[pre_ind], zscore_all[i][pre_ind])
        aucs_pre_by_trial.append(pre)
        post = auc(ts_control4average[post_ind], zscore_all[i][post_ind])
        aucs_post_by_trial.append(post)

    # print("pre mean",np.mean(np.asarray(aucs_pre_by_trial)),AUC_pre)
    # print("post mean",np.mean(np.asarray(aucs_post_by_trial)),AUC_post)
    AUC_pre_err = np.std(np.asarray(aucs_pre_by_trial)/np.sqrt(len(aucs_pre_by_trial)))
    AUC_post_err = np.std(np.asarray(aucs_post_by_trial)/np.sqrt(len(aucs_post_by_trial)))
    # print("errors",AUC_pre_err,AUC_post_err)

    
    # bin auc data in one second bins
    # create indexes
    aucs = [[] for i in range(len(zscore_all))]
    start = -perievent_options_dict["sec_before"]
    end = start + 1
    while end < perievent_options_dict["sec_after"]:
        ind = np.where((np.array(ts_control4average)<end) & (np.array(ts_control4average)>start))
        for i in range(len(zscore_all)):
            my_auc = auc(ts_control4average[ind], zscore_all[i][ind])
            aucs[i].append(my_auc)
        start = end
        end = start+1
        
    # create dictionary with data for AUC plot
    auc_data2plot_dict = {}
    auc_data2plot_dict['auc_data'] = AUC
    auc_data2plot_dict['p'] = p_val
    auc_data2plot_dict['auc_err_pre'] = AUC_pre_err
    auc_data2plot_dict['auc_err_post'] = AUC_post_err
    # add that dict to all analysed data dict
    analyzed_perievent_dict["auc"] = auc_data2plot_dict
    analyzed_perievent_dict["auc_by_sec"] = aucs
    analyzed_perievent_dict["aucs_pre_by_trial"] = aucs_pre_by_trial
    analyzed_perievent_dict["aucs_post_by_trial"] = aucs_post_by_trial
        
    return analyzed_perievent_dict
    
# create settings dataframe for csv
#  headers   
def get_settings_df(settings_dict):
    downsample_no = settings_dict[0]["downsample"]
    downsample_rate = settings_dict[0]["entered_downsample"]
    normalization = settings_dict[0]["normalization"]
    filter_on = settings_dict[0]["filter"]
    filter_window = settings_dict[0]["filter_window"]
    normalization_as = settings_dict[0]["show_norm_as"]
    subjects = settings_dict[0]["subject"]
    group = settings_dict[0]["subject_group_name"]
    my_df = pd.DataFrame({"rate after downsampling (Hz)":[downsample_rate],
                          "number of averaged samples for downsampling":[downsample_no],
                          "normalization":[normalization],
                          "normalization as": normalization_as,
                          "filter":[filter_on],
                          "filter window":[filter_window],
                          "subjects":subjects,
                          "subject group name":group})
    if "baseline_from_sec" in settings_dict[0]: # if that is in settings dictionary, all perievent settings are there
        my_df["baseline_from_sec"] = [settings_dict[0]["baseline_from_sec"]]
        my_df["baseline_to_sec"] = [settings_dict[0]["baseline_to_sec"]]
        my_df["auc_pre_from_sec"] = [settings_dict[0]["auc_pre_from"]]
        my_df["auc_pre_to_sec"] = [settings_dict[0]["auc_pre_to"]]
        my_df["auc_post_from_sec"] = [settings_dict[0]["auc_post_from"]]
        my_df["auc_post_to_sec"] = [settings_dict[0]["auc_post_to"]]
    my_df_transposed = my_df.transpose()
    return my_df_transposed

def get_last_timestamp(raw_data,signal_name):
    # return numpy array for all GCaMP channel raw data
    GCaMP_data = get_single_channel(raw_data,signal_name)
    ###############################
    # create new timestamps for the data (from zero)
    ts = create_timestamps(raw_data,signal_name,GCaMP_data)

    return ts[-1]


#################################################################################################
# PLOTS

def plot_raw(canvas,subject,raw_data,signal_name,control_name,export,save_plots,export_loc_data):
    # a tuple path+file beginning
    dump_path,file_beginning = export_loc_data
    # return numpy array for all GCaMP channel raw data
    GCaMP_data = get_single_channel(raw_data,signal_name)
    # return numpy arrayfor all control raw data
    control_data = get_single_channel(raw_data,control_name)
    #########
    # check
    # make sure both arrays are the same size
    len_all_GCaMP = GCaMP_data.size
    len_all_control = control_data.size
    # if GCaMP recording is longer, trim it to the control size
    if len_all_GCaMP > len_all_control:
        GCaMP_data = GCaMP_data[:len_all_control]
    # if control data is longer, trim it to GCaMP size
    elif len_all_GCaMP < len_all_control:
        control_data = control_data[:len_all_GCaMP]
    ###############################
    # create new timestamps for the data (from zero)
    ts = create_timestamps(raw_data,signal_name,GCaMP_data)
    
    # plot
    canvas.fig.clf()
    ax = canvas.fig.add_subplot(111)
    ax.cla()  # Clear the canvas
    # plot the lines
    ax.plot(ts, GCaMP_data,
             color=SIGNAL_COLOR_RGB,
             linewidth=1,
             label = "signal")
    ax.plot(ts,control_data,
             color=CONTROL_COLOR_RGB,
             linewidth=1,
             label = "control")

    # create title, axis labels, and legend
    my_title = 'Raw data: ' + subject
    canvas.fig.suptitle(my_title, fontsize=FIGURE_TITLE_FONT_SIZE)
    ax.set_xlabel('Time (sec)', fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_ylabel('mV', fontsize=AXIS_LABEL_FONTSIZE)
    # hide top and right border
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend(loc=MY_LEGEND_LOC, bbox_to_anchor=MY_LEGENG_POS)
    
#    ax.set_xlim([min(ts),max(ts)])
#    # Get current tick locations and append last to this array
#    x_ticks = np.append(ax.get_xticks(), max(ts))
#    # Set xtick locations to the values of the array `x_ticks`
#    ax.set_xticks(x_ticks)
    # set the spacing between subplots 
    canvas.fig.subplots_adjust(left=MY_LEFT, 
                    bottom=MY_BOTTOM,  
                    right=MY_RIGHT,  
                    top=MY_TOP,  
                    wspace=MY_WSPACE,  
                    hspace=MY_HSPACE) 
    
    # create subfolder in current subject folder
    if export == True:
        file_name = file_beginning+"_raw.csv"
        dump_file_path = os.path.join(dump_path,file_name)
        sig = "Signal_"+signal_name
        cont = "Control_"+control_name
        raw_df= pd.DataFrame({"Time (sec)":ts,
                                 sig +' (mV)' : GCaMP_data, 
                                 cont +' (mV)': control_data})
        raw_df.to_csv(dump_file_path, index=False)  
            
    if export == True and save_plots == True:
        plot_file_name = file_beginning+"_raw.png"
        dump_plot_file_path = os.path.join(dump_path,plot_file_name)
        canvas.fig.savefig(dump_plot_file_path)
        # also as svg
        plot_file_name = file_beginning+"_raw.svg"
        dump_plot_file_path = os.path.join(dump_path,plot_file_name)
        canvas.fig.savefig(dump_plot_file_path, format='svg', dpi=DPI4SVG)
    else:
        canvas.draw() 
    
def plot_trimmed(canvas,subject,signal_dict,export,save_plots,export_loc_data,signal_name,control_name):
    # a tuple path+file beginning
    dump_path,file_beginning = export_loc_data
    GCaMP_data = signal_dict["signal"]
    control_data = signal_dict["control"]
    ts = signal_dict["ts"]

    total_seconds = ts[-1]-ts[0]
    # start time from zero
    ts_reset = [i*total_seconds/len(ts) for i in range(len(ts))]
    
    
    # create subfolder in current subject folder
    if export == True:
        file_name = file_beginning + "_raw_trimmed.csv"
        dump_file_path = os.path.join(dump_path,file_name)
        sig = "Signal_"+signal_name
        cont = "Control_"+control_name
        raw_df= pd.DataFrame({"Time (sec)":ts_reset,
                                 sig +' (mV)': GCaMP_data, 
                                 cont +' (mV)': control_data})
        raw_df.to_csv(dump_file_path, index=False)  
    
    
#    print(len(GCaMP_data),len(control_data),len(ts))
    # plot
    # clear previous figure
    canvas.fig.clf()
    ax = canvas.fig.add_subplot(111)
   
    # plot the lines
    ax.plot(ts_reset, GCaMP_data,
             color=SIGNAL_COLOR_RGB,
             linewidth=1,
             label = "signal")
    ax.plot(ts_reset,control_data,
             color=CONTROL_COLOR_RGB,
             linewidth=1,
             label = "control")
    
    # create title, axis labels, and legend
    my_title = 'Raw data: ' + subject
    canvas.fig.suptitle(my_title, fontsize=FIGURE_TITLE_FONT_SIZE)
    ax.set_xlabel('Time (sec)', fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_ylabel('mV', fontsize=AXIS_LABEL_FONTSIZE)
    # hide top and right border
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend(loc=MY_LEGEND_LOC, bbox_to_anchor=MY_LEGENG_POS)
    
#    ax.set_xlim([min(ts),max(ts)])
#    # Get current tick locations and append last to this array
#    x_ticks = np.append(ax.get_xticks(), max(ts))
#    # Set xtick locations to the values of the array `x_ticks`
#    ax.set_xticks(x_ticks)
#    
    # set the spacing between subplots 
    canvas.fig.subplots_adjust(left=MY_LEFT, 
                    bottom=MY_BOTTOM,  
                    right=MY_RIGHT,  
                    top=MY_TOP,  
                    wspace=MY_WSPACE,  
                    hspace=MY_HSPACE) 
    
    if export == True and save_plots == True:
        plot_file_name = file_beginning + "_raw_trimmed.png"
        dump_plot_file_path = os.path.join(dump_path,plot_file_name)
        canvas.fig.savefig(dump_plot_file_path)
        # also as svg
        plot_file_name = file_beginning+"_raw_trimmed.svg"
        dump_plot_file_path = os.path.join(dump_path,plot_file_name)
        canvas.fig.savefig(dump_plot_file_path, format='svg', dpi=DPI4SVG)
    else:
        canvas.draw()
    

def plot_with_event(canvas,subject,signal_dict,event_name,event2_name,event_data,export,save_plots,export_loc_data,signal_name,control_name):
    # a tuple path+file beginning
    dump_path,file_beginning = export_loc_data
    try:
        GCaMP_data = signal_dict["signal"]
        control_data = signal_dict["control"]
        ts = signal_dict["ts"]
    except:
        # return numpy array for all GCaMP channel raw data
        GCaMP_data = get_single_channel(signal_dict,signal_name)
        # return numpy arrayfor all control raw data
        control_data = get_single_channel(signal_dict,control_name)
        #########
        # check
        # make sure both arrays are the same size
        len_all_GCaMP = GCaMP_data.size
        len_all_control = control_data.size
        # if GCaMP recording is longer, trim it to the control size
        if len_all_GCaMP > len_all_control:
            GCaMP_data = GCaMP_data[:len_all_control]
        # if control data is longer, trim it to GCaMP size
        elif len_all_GCaMP < len_all_control:
            control_data = control_data[:len_all_GCaMP]
        ###############################
        # create new timestamps for the data (from zero)
        ts = create_timestamps(signal_dict,signal_name,GCaMP_data)
        
    total_seconds = ts[-1]-ts[0]
    # start time from zero
    ts_reset = [i*total_seconds/len(ts) for i in range(len(ts))]
    
    # plot
    # clear previous figure
    canvas.fig.clf()
    ax = canvas.fig.add_subplot(111)
    
    # plot the lines
    ax.plot(ts_reset, GCaMP_data,
             color=SIGNAL_COLOR_RGB,
             linewidth=1,
             label = "signal")
    ax.plot(ts_reset,control_data,
             color=CONTROL_COLOR_RGB,
             linewidth=1,
             label = "control")
    
    # get accurate event data timing
    first_evt_onset = event_data[0][0]
    first_evt_onset_translated = [el-ts[0] for el in first_evt_onset]
    first_evt_offset = event_data[0][1]
    first_evt_offset_translated = [el-ts[0] for el in first_evt_offset]
    if len(event_data) > 1 and len(event_data[1][0]) > 0:
        second_evt_onset = event_data[1][0]
        second_evt_onset_translated = [el-ts[0] for el in second_evt_onset]
        second_evt_offset = event_data[1][1]
        second_evt_offset_translated = [el-ts[0] for el in second_evt_offset]
        
    if len(first_evt_onset) > 0:
        # first event
        # plot tone onsets as vertical lines
        # first element in tone_on_off is a list with onset times
        ctr = 0
        for i in range(len(first_evt_onset_translated)):
            if first_evt_onset_translated[i] >= ts_reset[0] and first_evt_onset_translated[i] <= ts_reset[-1]:
                if ctr == 0:  # label only first one for the legend
                    ax.axvline(x=first_evt_onset_translated[i],color='k',alpha=0.5,linewidth=1,label=event_name)
                else:
                    ax.axvline(x=first_evt_onset_translated[i],color='k',alpha=0.5,linewidth=1)
                ctr+=1
        # if there is more than one event, plot the second one as well
        if len(event_data) > 1 and len(event_data[1][0]) > 0:
            ctr = 0
            for i in range(len(second_evt_onset_translated)):
                if second_evt_onset_translated[i] >= ts_reset[0] and second_evt_onset_translated[i] <= ts_reset[-1]:
                    if ctr == 0:  # label only first one for the legend
                        ax.axvline(x=second_evt_onset_translated[i],linestyle='dashed',color='k',alpha=0.5,linewidth=1,label=event2_name)
                    else:
                        ax.axvline(x=second_evt_onset_translated[i],linestyle='dashed',color='k',alpha=0.5,linewidth=1)
                    ctr+=1
    
    # create title, axis labels, and legend
    my_title = 'Raw data: ' + subject
    canvas.fig.suptitle(my_title, fontsize=FIGURE_TITLE_FONT_SIZE)
    ax.set_xlabel('Time (sec)', fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_ylabel('mV', fontsize=AXIS_LABEL_FONTSIZE)
    # hide top and right border
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend(loc=MY_LEGEND_LOC, bbox_to_anchor=MY_LEGENG_POS)
    
#    ax.set_xlim([min(ts),max(ts)])
#    # Get current tick locations and append last to this array
#    x_ticks = np.append(ax.get_xticks(), max(ts))
#    # Set xtick locations to the values of the array `x_ticks`
#    ax.set_xticks(x_ticks)
    # set the spacing between subplots 
    canvas.fig.subplots_adjust(left=MY_LEFT, 
                    bottom=MY_BOTTOM,  
                    right=MY_RIGHT,  
                    top=MY_TOP,  
                    wspace=MY_WSPACE,  
                    hspace=MY_HSPACE) 
    
    # create subfolder in current subject folder
    if export == True:
        file_name = file_beginning+"_raw.csv"
        dump_file_path = os.path.join(dump_path,file_name)
        sig = "Signal_"+signal_name
        cont = "Control_"+control_name
        raw_df= pd.DataFrame({"Time (sec)":ts_reset,
                                 sig +' (mV)': GCaMP_data, 
                                 cont +' (mV)': control_data})
        raw_df.to_csv(dump_file_path, index=False)  
        
        event1_file = file_beginning+"_"+event_name+".csv"
        dump_file_path = os.path.join(dump_path,event1_file)
        evt1_df = pd.DataFrame({"onset (s)":first_evt_onset_translated,
                                    "offset (s)" : first_evt_offset_translated})
        evt1_df.to_csv(dump_file_path, index=False)
        if len(event_data) > 1 and len(event_data[1][0]) > 0:
            event2_file = file_beginning+"_"+event2_name+".csv"
            dump_file_path = os.path.join(dump_path,event2_file)
            evt2_df = pd.DataFrame({"onset (s)":second_evt_onset_translated,
                                       "offset (s)" : second_evt_offset_translated})
            evt2_df.to_csv(dump_file_path, index=False)
            
    if export == True and save_plots == True:
        if event2_name == "---":
            plot_file_name = file_beginning+"_raw_"+event_name+".png"
            plot_file_name2 = file_beginning+"_raw_"+event_name+".svg"
        else:
            plot_file_name = file_beginning+"_raw_"+event_name+"_"+event2_name+".png"
            plot_file_name2 = file_beginning+"_raw_"+event_name+".svg"
        dump_plot_file_path = os.path.join(dump_path,plot_file_name)
        canvas.fig.savefig(dump_plot_file_path)
        # also as svg
        dump_plot_file_path = os.path.join(dump_path,plot_file_name2)
        canvas.fig.savefig(dump_plot_file_path, format='svg', dpi=DPI4SVG)
        
    else:
        canvas.draw() 
    

def plot_downsampled_alone(canvas,options_dict,downsampled_signal_dict,export,save_plots,export_loc_data,settings_dict,signal_name,control_name):
    settings_dict[0]["subject"] = options_dict["subject"]
    settings_dict[0]["subject_group_name"] = options_dict["subject_group_name"]
    # a tuple path+file beginning
    dump_path,file_beginning = export_loc_data
    # downsampled data
    GCaMP_downsampled_data = downsampled_signal_dict["signal"]
    control_downsampled_data = downsampled_signal_dict["control"]
    ts_downsampled = downsampled_signal_dict["ts"]
    total_seconds = ts_downsampled[-1]-ts_downsampled[0]
    # start time from zero
    ts_reset = [i*total_seconds/len(ts_downsampled) for i in range(len(ts_downsampled))]
    # plot
    # clear previous figure
    canvas.fig.clf()
    ax = canvas.fig.add_subplot(111)

    # plot downsampled
    ax.plot(ts_reset,GCaMP_downsampled_data,
             color=SIGNAL_COLOR_RGB,
             linewidth=1,
             label = "signal")
    ax.plot(ts_reset,control_downsampled_data,
             color=CONTROL_COLOR_RGB,
             linewidth=1,
             label = "control")
    # create title, axis labels, and legend
    my_title = 'Subject: ' + options_dict["subject"]
    canvas.fig.suptitle(my_title, fontsize=FIGURE_TITLE_FONT_SIZE)
    ax.set_title('Downsampled',fontdict={'fontsize': TITLE_FONT_SIZE, 'fontweight': 'medium'})
    ax.set_xlabel('Time (sec)', fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_ylabel('mV', fontsize=AXIS_LABEL_FONTSIZE)
    ax.legend(loc=MY_LEGEND_LOC, bbox_to_anchor=MY_LEGENG_POS)
    # hide top and right border
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # set the spacing between subplots 
    canvas.fig.subplots_adjust(left=MY_LEFT, 
                    bottom=MY_BOTTOM,  
                    right=MY_RIGHT,  
                    top=MY_TOP,  
                    wspace=MY_WSPACE,  
                    hspace=MY_HSPACE) 
    
    # create subfolder in current subject folder
    if export == True:
        file_name = file_beginning+"_downsampled_trimmed.csv"
        dump_file_path = os.path.join(dump_path,file_name)
        settings_df = get_settings_df(settings_dict)
        sig = "Signal_"+signal_name
        cont = "Control_"+control_name
        raw_df= pd.DataFrame({"Time (seconds)":ts_reset,
                                 sig +' (mV)': GCaMP_downsampled_data, 
                                 cont +' (mV)': control_downsampled_data})
        settings_df.to_csv(dump_file_path,header=False)            
        with open(dump_file_path,'a') as f:
            f.write("\n")
        raw_df.to_csv(dump_file_path, index=False,mode='a')
              
                 
    if export == True and save_plots == True:
        plot_file_name = file_beginning+"_downsampled_trimmed.png"
        dump_plot_file_path = os.path.join(dump_path,plot_file_name)
        canvas.fig.savefig(dump_plot_file_path)
        # also as svg
        plot_file_name = file_beginning+"_downsampled_trimmed.svg"
        dump_plot_file_path = os.path.join(dump_path,plot_file_name)
        canvas.fig.savefig(dump_plot_file_path, format='svg', dpi=DPI4SVG)
    else:
        canvas.draw()

    
def plot_downsampled_alone_with_event(canvas,options_dict,downsampled_signal_dict,event_name,event2_name,event_data,export,save_plots,export_loc_data,settings_dict,signal_name,control_name):
    settings_dict[0]["subject"] = options_dict["subject"]
    settings_dict[0]["subject_group_name"] = options_dict["subject_group_name"]
    # a tuple path+file beginning
    dump_path,file_beginning = export_loc_data
    # downsampled data
    GCaMP_downsampled_data = downsampled_signal_dict["signal"]
    control_downsampled_data = downsampled_signal_dict["control"]
    ts_downsampled = downsampled_signal_dict["ts"]
    total_seconds = ts_downsampled[-1]-ts_downsampled[0]
    # start time from zero
    ts_reset = [i*total_seconds/len(ts_downsampled) for i in range(len(ts_downsampled))]
    # plot
    # clear previous figure
    canvas.fig.clf()
    ax = canvas.fig.add_subplot(111)

    # plot downsampled
    ax.plot(ts_reset,GCaMP_downsampled_data,
             color=SIGNAL_COLOR_RGB,
             linewidth=1,
             label = "signal")
    ax.plot(ts_reset,control_downsampled_data,
             color=CONTROL_COLOR_RGB,
             linewidth=1,
             label = "control")
    
    # get accurate event data timing
    first_evt_onset = event_data[0][0]
    first_evt_onset_translated = [el-ts_downsampled[0] for el in first_evt_onset]
    first_evt_offset = event_data[0][1]
    first_evt_offset_translated = [el-ts_downsampled[0] for el in first_evt_offset]
    if len(event_data) > 1 and len(event_data[1][0]) > 0:
        second_evt_onset = event_data[1][0]
        second_evt_onset_translated = [el-ts_downsampled[0] for el in second_evt_onset]
        second_evt_offset = event_data[1][1]
        second_evt_offset_translated = [el-ts_downsampled[0] for el in second_evt_offset]
    
    if len(first_evt_onset) > 0:
        # first event
        # plot tone onsets as vertical lines
        # first element in tone_on_off is a list with onset times
        ctr = 0
        for i in range(len(first_evt_onset_translated)):
            if first_evt_onset_translated[i] >= ts_reset[0] and first_evt_onset_translated[i] <= ts_reset[-1]:
                if ctr == 0:  # label only first one for the legend
                    ax.axvline(x=first_evt_onset_translated[i],color='k',alpha=0.5,linewidth=1,label=event_name)
                else:
                    ax.axvline(x=first_evt_onset_translated[i],color='k',alpha=0.5,linewidth=1)
                ctr+=1
        # if there is more than one event, plot the second one as well
        if len(event_data) > 1 and len(event_data[1][0]) > 0:
            ctr = 0
            for i in range(len(second_evt_onset_translated)):
                if second_evt_onset_translated[i] >= ts_reset[0] and second_evt_onset_translated[i] <= ts_reset[-1]:
                    if ctr == 0:  # label only first one for the legend
                        ax.axvline(x=second_evt_onset_translated[i],linestyle='dashed',color='k',alpha=0.5,linewidth=1,label=event2_name)
                    else:
                        ax.axvline(x=second_evt_onset_translated[i],linestyle='dashed',color='k',alpha=0.5,linewidth=1)
                    ctr+=1
    
    # create title, axis labels, and legend
    my_title = 'Subject: ' + options_dict["subject"]
    canvas.fig.suptitle(my_title, fontsize=FIGURE_TITLE_FONT_SIZE)
    ax.set_title('Downsampled',fontdict={'fontsize': TITLE_FONT_SIZE, 'fontweight': 'medium'})
    ax.set_xlabel('Time (sec)', fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_ylabel('mV', fontsize=AXIS_LABEL_FONTSIZE)
    ax.legend(loc=MY_LEGEND_LOC, bbox_to_anchor=MY_LEGENG_POS)
    # hide top and right border
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # set the spacing between subplots 
    canvas.fig.subplots_adjust(left=MY_LEFT, 
                    bottom=MY_BOTTOM,  
                    right=MY_RIGHT,  
                    top=MY_TOP,  
                    wspace=MY_WSPACE,  
                    hspace=MY_HSPACE) 
    
    # create subfolder in current subject folder
    if export == True:
        file_name = file_beginning+"_downsampled_trimmed.csv"
        dump_file_path = os.path.join(dump_path,file_name)
        settings_df = get_settings_df(settings_dict)
        sig = "Signal_"+signal_name
        cont = "Control_"+control_name
        raw_df= pd.DataFrame({"Time (sec)":ts_reset,
                                 sig +' (mV)': GCaMP_downsampled_data, 
                                 cont +' (mV)': control_downsampled_data})
        settings_df.to_csv(dump_file_path,header=False)            
        with open(dump_file_path,'a') as f:
            f.write("\n")
        raw_df.to_csv(dump_file_path, index=False,mode='a')
            
        event1_file = file_beginning+"_"+event_name+".csv"
        dump_file_path = os.path.join(dump_path,event1_file)
        evt1_df = pd.DataFrame({"onset (s)":first_evt_onset_translated,
                                    "offset (s)" : first_evt_offset_translated})
        evt1_df.to_csv(dump_file_path, index=False)
        if len(event_data) > 1 and len(event_data[1][0]) > 0:
            event2_file = file_beginning+"_"+event2_name+".csv"
            dump_file_path = os.path.join(dump_path,event2_file)
            evt2_df = pd.DataFrame({"onset (s)":second_evt_onset_translated,
                                       "offset (s)" : second_evt_offset_translated})
            evt2_df.to_csv(dump_file_path, index=False)
            
    if export == True and save_plots == True:
        if event2_name == "---":
            plot_file_name = file_beginning+"_downsampled_trimmed_"+event_name+".png"
            plot_file_name2 = file_beginning+"_downsampled_trimmed_"+event_name+".svg"
        else:
            plot_file_name = file_beginning+"_downsampled_trimmed_"+event_name+"_"+event2_name+".png"
            plot_file_name2 = file_beginning+"_downsampled_trimmed_"+event_name+"_"+event2_name+".svg"
        dump_plot_file_path = os.path.join(dump_path,plot_file_name)
        canvas.fig.savefig(dump_plot_file_path)
        # also as svg
        dump_plot_file_path = os.path.join(dump_path,plot_file_name2)
        canvas.fig.savefig(dump_plot_file_path, format='svg', dpi=DPI4SVG)
    else:
        canvas.draw() 

    
def plot_normalized_alone(canvas,options_dict,normalized_dict,export,save_plots,export_loc_data,settings_dict):
    settings_dict[0]["subject"] = options_dict["subject"]
    settings_dict[0]["subject_group_name"] = options_dict["subject_group_name"]
    show_norm_as = settings_dict[0]["show_norm_as"]
    # a tuple path+file beginning
    dump_path,file_beginning = export_loc_data
    # normalized data
    normalized_data = normalized_dict["normalized_signal"]
    ts_normalized = normalized_dict["ts"]
    # plot
    # clear previous figure
    canvas.fig.clf()
    ax = canvas.fig.add_subplot(111)
    
    total_seconds = ts_normalized[-1]-ts_normalized[0]
    # start time from zero
    ts_reset = [i*total_seconds/len(ts_normalized) for i in range(len(ts_normalized))]
    
    # plot downsampled that starts from zero
    ax.plot(ts_reset,normalized_data,
             color=SIGNAL_COLOR_RGB,
             linewidth=1,
             label = "normalized")
    # create title, axis labels, and legend
    my_title = 'Subject: ' + options_dict["subject"]
    canvas.fig.suptitle(my_title, fontsize=FIGURE_TITLE_FONT_SIZE)
    ax.set_title('Normalized',fontdict={'fontsize': TITLE_FONT_SIZE, 'fontweight': 'medium'})
    ax.set_xlabel('Time (sec)', fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_ylabel('dF/F (%)', fontsize=AXIS_LABEL_FONTSIZE)
    if show_norm_as == "Z-Score":
        ax.set_ylabel("Z-Score",fontsize=AXIS_LABEL_FONTSIZE)
    ax.legend(loc=MY_LEGEND_LOC, bbox_to_anchor=MY_LEGENG_POS)
    # hide top and right border
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # set the spacing between subplots 
    canvas.fig.subplots_adjust(left=MY_LEFT, 
                    bottom=MY_BOTTOM,  
                    right=MY_RIGHT,  
                    top=MY_TOP,  
                    wspace=MY_WSPACE,  
                    hspace=MY_HSPACE) 

    
    # create subfolder in current subject folder
    if export == True:
        file_name = file_beginning+"_normalized.csv"
        if file_beginning != options_dict["subject"]:
            file_name = file_beginning+"_"+options_dict["subject"]+"_normalized.csv"
        dump_file_path = os.path.join(dump_path,file_name)
        settings_df = get_settings_df(settings_dict)
        if show_norm_as == 'Z-Score':
            normalized_header = "signal (z-score)"
        else:

            normalized_header = "signal (%df/F)"
        raw_df= pd.DataFrame({"Time (sec)":ts_reset,
                                 normalized_header: normalized_data})
        settings_df.to_csv(dump_file_path,header=False)            
        with open(dump_file_path,'a') as f:
            f.write("\n")
        raw_df.to_csv(dump_file_path, index=False,mode='a') 
    if export == True and save_plots == True:
        plot_file_name = file_beginning+"_normalized.png"
        plot_file_name2 = file_beginning+"_normalized.svg"
        if file_beginning != options_dict["subject"]:
            plot_file_name = file_beginning+"_"+options_dict["subject"]+"_normalized.png"
            plot_file_name2 = file_beginning+"_"+options_dict["subject"]+"_normalized.svg"
        dump_plot_file_path = os.path.join(dump_path,plot_file_name)
        canvas.fig.savefig(dump_plot_file_path)
        # also as svg
        dump_plot_file_path = os.path.join(dump_path,plot_file_name2)
        canvas.fig.savefig(dump_plot_file_path, format='svg', dpi=DPI4SVG)
    else:
        canvas.draw()
    
def plot_normalized_alone_with_event(canvas,options_dict,normalized_dict,event_name,event2_name,event_data,export,save_plots,export_loc_data,settings_dict):
    settings_dict[0]["subject"] = options_dict["subject"]
    settings_dict[0]["subject_group_name"] = options_dict["subject_group_name"]
    show_norm_as = settings_dict[0]["show_norm_as"]
    # a tuple path+file beginning
    dump_path,file_beginning = export_loc_data
    # normalized data
    normalized_data = normalized_dict["normalized_signal"]
    ts_normalized = normalized_dict["ts"]
    # plot
    # clear previous figure
    canvas.fig.clf()
    ax = canvas.fig.add_subplot(111)
    
    total_seconds = ts_normalized[-1]-ts_normalized[0]
    # start time from zero
    ts_reset = [i*total_seconds/len(ts_normalized) for i in range(len(ts_normalized))]
    
    # plot downsampled
    ax.plot(ts_reset,normalized_data,
             color=SIGNAL_COLOR_RGB,
             linewidth=1,
             label = "normalized")

#    # plot downsampled
#    ax.plot(ts_normalized,normalized_data,
#             color='green',
#             linewidth=1,
#             label = "normalized")
    
    # get accurate event data timing
    first_evt_onset = event_data[0][0]
    first_evt_onset_translated = [el-ts_normalized[0] for el in first_evt_onset]
    first_evt_offset = event_data[0][1]
    first_evt_offset_translated = [el-ts_normalized[0] for el in first_evt_offset]
    if len(event_data) > 1 and len(event_data[1][0]) > 0:
        second_evt_onset = event_data[1][0]
        second_evt_onset_translated = [el-ts_normalized[0] for el in second_evt_onset]
        second_evt_offset = event_data[1][1]
        second_evt_offset_translated = [el-ts_normalized[0] for el in second_evt_offset]
    
    if len(first_evt_onset) > 0:
        # first event
        # plot tone onsets as vertical lines
        # first element in tone_on_off is a list with onset times
        ctr = 0
        for i in range(len(first_evt_onset_translated)):
            if first_evt_onset_translated[i] >= ts_reset[0] and first_evt_onset_translated[i] <= ts_reset[-1]:
                if ctr == 0:  # label only first one for the legend
                    ax.axvline(x=first_evt_onset_translated[i],color='k',alpha=0.5,linewidth=1,label=event_name)
                else:
                    ax.axvline(x=first_evt_onset_translated[i],color='k',alpha=0.5,linewidth=1)
                ctr+=1
        # if there is more than one event, plot the second one as well
        if len(event_data) > 1 and len(event_data[1][0]) > 0:
            ctr = 0
            for i in range(len(second_evt_onset_translated)):
                if second_evt_onset_translated[i] >= ts_reset[0] and second_evt_onset_translated[i] <= ts_reset[-1]:
                    if ctr == 0:  # label only first one for the legend
                        ax.axvline(x=second_evt_onset_translated[i],linestyle='dashed',color='k',alpha=0.5,linewidth=1,label=event2_name)
                    else:
                        ax.axvline(x=second_evt_onset_translated[i],linestyle='dashed',color='k',alpha=0.5,linewidth=1)
                    ctr+=1
    
    # create title, axis labels, and legend
    my_title = 'Subject: ' + options_dict["subject"]
    canvas.fig.suptitle(my_title, fontsize=FIGURE_TITLE_FONT_SIZE)
    ax.set_title('Normalized',fontdict={'fontsize': TITLE_FONT_SIZE, 'fontweight': 'medium'})
    ax.set_xlabel('Time (sec)', fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_ylabel('dF/F (%)', fontsize=AXIS_LABEL_FONTSIZE)
    if show_norm_as == "Z-Score":
        ax.set_ylabel("Z-Score",fontsize=AXIS_LABEL_FONTSIZE)
    ax.legend(loc=MY_LEGEND_LOC, bbox_to_anchor=MY_LEGENG_POS)
    # hide top and right border
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # set the spacing between subplots 
    canvas.fig.subplots_adjust(left=MY_LEFT, 
                    bottom=MY_BOTTOM,  
                    right=MY_RIGHT,  
                    top=MY_TOP,  
                    wspace=MY_WSPACE,  
                    hspace=MY_HSPACE) 
    
    # create subfolder in current subject folder
    if export == True:
        file_name = file_beginning+"_normalized.csv"
        if file_beginning != options_dict["subject"]:
            file_name = file_beginning+"_"+options_dict["subject"]+"_normalized.csv"
        dump_file_path = os.path.join(dump_path,file_name)
        settings_df = get_settings_df(settings_dict)
        if show_norm_as == 'Z-Score':
            normalized_header = "signal (z-score)"
        else:

            normalized_header = "signal (%df/F)"
        raw_df= pd.DataFrame({"Time (sec)":ts_reset,
                                 normalized_header : normalized_data})
        settings_df.to_csv(dump_file_path,header=False)            
        with open(dump_file_path,'a') as f:
                f.write("\n")
        raw_df.to_csv(dump_file_path, index=False,mode='a')
            
        event1_file = file_beginning+"_"+event_name+".csv"
        dump_file_path = os.path.join(dump_path,event1_file)
        evt1_df = pd.DataFrame({"onset (s)":first_evt_onset_translated,
                                    "offset (s)" : first_evt_offset_translated})
        evt1_df.to_csv(dump_file_path, index=False)
        if len(event_data) > 1 and len(event_data[1][0]) > 0:
            event2_file = file_beginning+"_"+event2_name+".csv"
            dump_file_path = os.path.join(dump_path,event2_file)
            evt2_df = pd.DataFrame({"onset (s)":second_evt_onset_translated,
                                       "offset (s)" : second_evt_offset_translated})
            evt2_df.to_csv(dump_file_path, index=False)
            
    if export == True and save_plots == True:
        if event2_name == "---":
            plot_file_name = file_beginning+"_normalized_"+event_name+".png"
            plot_file_name2 = file_beginning+"_normalized_"+event_name+".svg"
            if file_beginning != options_dict["subject"]:
                plot_file_name = file_beginning+"_"+options_dict["subject"]+"_normalized_"+event_name+".png"
                plot_file_name2 = file_beginning+"_"+options_dict["subject"]+"_normalized_"+event_name+".svg"
        else:
            plot_file_name = file_beginning+"_normalized_"+event_name+"_"+event2_name+".png"
            plot_file_name2 = file_beginning+"_normalized_"+event_name+"_"+event2_name+".svg"
            if file_beginning != options_dict["subject"]:
                plot_file_name = file_beginning+"_"+options_dict["subject"]+"_normalized_"+event_name+"_"+event2_name+".png"
                plot_file_name2 = file_beginning+"_"+options_dict["subject"]+"_normalized_"+event_name+"_"+event2_name+".svg"
        dump_plot_file_path = os.path.join(dump_path,plot_file_name)
        canvas.fig.savefig(dump_plot_file_path)
        # also as svg
        dump_plot_file_path = os.path.join(dump_path,plot_file_name2)
        canvas.fig.savefig(dump_plot_file_path, format='svg', dpi=DPI4SVG)
    else:
        canvas.draw() 
    
    
def plot_downsampled_and_normalized_alone(canvas,subject,downsampled_signal_dict,normalized_dict,show_norm_as):
    # downsampled data
    GCaMP_downsampled_data = downsampled_signal_dict["signal"]
    control_downsampled_data = downsampled_signal_dict["control"]
    ts_downsampled = downsampled_signal_dict["ts"]
    # normalized data
    normalized_data = normalized_dict["normalized_signal"]
    ts_normalized = normalized_dict["ts"]
    
    total_seconds_downsampled = ts_downsampled[-1]-ts_downsampled[0]
    # start time from zero
    ts_reset_downsampled = [i*total_seconds_downsampled/len(ts_downsampled) for i in range(len(ts_downsampled))]
    total_seconds_normalized = ts_normalized[-1]-ts_normalized[0]
    # start time from zero
    ts_reset_normalized = [i*total_seconds_normalized/len(ts_normalized) for i in range(len(ts_normalized))]
    
    # plot
    # clear previous figure
    canvas.fig.clf()
    ax = canvas.fig.add_subplot(211)
    ax2 = canvas.fig.add_subplot(212,sharex=ax)
    
    # plot downsampled
    ax.plot(ts_reset_downsampled,GCaMP_downsampled_data,
             color=SIGNAL_COLOR_RGB,
             linewidth=1,
             label = "signal")
    ax.plot(ts_reset_downsampled,control_downsampled_data,
             color=CONTROL_COLOR_RGB,
             linewidth=1,
             label = "control")
    ax2.plot(ts_reset_normalized,normalized_data,
             color=SIGNAL_COLOR_RGB,
             linewidth=1,
             label = "normalized")
    # create title, axis labels, and legend
    my_title = 'Subject: ' + subject
    canvas.fig.suptitle(my_title, fontsize=FIGURE_TITLE_FONT_SIZE)
    ax.set_title('Downsampled',fontdict={'fontsize': TITLE_FONT_SIZE, 'fontweight': 'medium'})
    ax.set_xlabel('Time (sec)', fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_ylabel('mV', fontsize=AXIS_LABEL_FONTSIZE)
    ax.legend(loc=MY_LEGEND_LOC, bbox_to_anchor=MY_LEGENG_POS)
    # hide top and right border
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax2.set_title('Normalized',fontdict={'fontsize': TITLE_FONT_SIZE, 'fontweight': 'medium'})
    ax2.set_xlabel('Time (sec)', fontsize=AXIS_LABEL_FONTSIZE)
    ax2.set_ylabel('dF/F (%)', fontsize=AXIS_LABEL_FONTSIZE)
    if show_norm_as == "Z-Score":
        ax2.set_ylabel("Z-Score",fontsize=AXIS_LABEL_FONTSIZE)
    ax2.legend(loc=MY_LEGEND_LOC, bbox_to_anchor=MY_LEGENG_POS)
    # hide top and right border
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    # set the spacing between subplots 
    canvas.fig.subplots_adjust(left=MY_LEFT, 
                    bottom=MY_BOTTOM,  
                    right=MY_RIGHT,  
                    top=MY_TOP,  
                    wspace=MY_WSPACE,  
                    hspace=MY_HSPACE) 
    canvas.draw()
    
def plot_downsampled_and_normalized_with_event(canvas,subject,downsampled_signal_dict,normalized_dict,show_norm_as,event_name,event2_name,event_data):
    # downsampled data
    GCaMP_downsampled_data = downsampled_signal_dict["signal"]
    control_downsampled_data = downsampled_signal_dict["control"]
    ts_downsampled = downsampled_signal_dict["ts"]
    # normalized data
    normalized_data = normalized_dict["normalized_signal"]
    ts_normalized = normalized_dict["ts"]
    
    total_seconds_downsampled = ts_downsampled[-1]-ts_downsampled[0]
    # start time from zero
    ts_reset_downsampled = [i*total_seconds_downsampled/len(ts_downsampled) for i in range(len(ts_downsampled))]
    total_seconds_normalized = ts_normalized[-1]-ts_normalized[0]
    # start time from zero
    ts_reset_normalized = [i*total_seconds_normalized/len(ts_normalized) for i in range(len(ts_normalized))]
    
    # plot
    # clear previous figure
    canvas.fig.clf()
    ax = canvas.fig.add_subplot(211)
    ax2 = canvas.fig.add_subplot(212,sharex=ax)
    
    # plot downsampled
    ax.plot(ts_reset_downsampled,GCaMP_downsampled_data,
             color=SIGNAL_COLOR_RGB,
             linewidth=1,
             label = "signal")
    ax.plot(ts_reset_downsampled,control_downsampled_data,
             color=CONTROL_COLOR_RGB,
             linewidth=1,
             label = "control")
    ax2.plot(ts_reset_normalized,normalized_data,
             color=SIGNAL_COLOR_RGB,
             linewidth=1,
             label = "normalized")
    
    # get accurate event data timing
    first_evt_onset = event_data[0][0]
    first_evt_onset_translated = [el-ts_downsampled[0] for el in first_evt_onset]
    first_evt_offset = event_data[0][1]
    first_evt_offset_translated = [el-ts_downsampled[0] for el in first_evt_offset]
    if len(event_data) > 1 and len(event_data[1][0]) > 0:
        second_evt_onset = event_data[1][0]
        second_evt_onset_translated = [el-ts_downsampled[0] for el in second_evt_onset]
        second_evt_offset = event_data[1][1]
        second_evt_offset_translated = [el-ts_downsampled[0] for el in second_evt_offset]
        
    if len(first_evt_onset) > 0:
        # first event
        # plot tone onsets as vertical lines
        # first element in tone_on_off is a list with onset times
        ctr = 0
        ctr1 = 0
        for i in range(len(first_evt_onset_translated)):
            if first_evt_onset_translated[i] >= ts_reset_downsampled[0] and first_evt_onset_translated[i] <= ts_reset_downsampled[-1]:
                if ctr == 0:  # label only first one for the legend
                    ax.axvline(x=first_evt_onset_translated[i],color='k',alpha=0.5,linewidth=1,label=event_name)
                else:
                    ax.axvline(x=first_evt_onset_translated[i],color='k',alpha=0.5,linewidth=1)
                ctr+=1
            if first_evt_onset_translated[i] >= ts_reset_normalized[0] and first_evt_onset_translated[i] <= ts_reset_normalized[-1]:
                if ctr1 == 0:  # label only first one for the legend
                    ax2.axvline(x=first_evt_onset_translated[i],color='k',alpha=0.5,linewidth=1,label=event_name)
                else:
                    ax2.axvline(x=first_evt_onset_translated[i],color='k',alpha=0.5,linewidth=1)
                ctr1+=1
        # if there is more than one event, plot the second one as well
        if len(event_data) > 1 and len(event_data[1][0]) > 0:
            ctr = 0
            ctr1 = 0
            for i in range(len(second_evt_onset_translated)):
                if second_evt_onset_translated[i] >= ts_reset_downsampled[0] and second_evt_onset_translated[i] <= ts_reset_downsampled[-1]:
                    if ctr == 0:  # label only first one for the legend
                        ax.axvline(x=second_evt_onset_translated[i],linestyle='dashed',color='k',alpha=0.5,linewidth=1,label=event2_name)
                    else:
                        ax.axvline(x=second_evt_onset_translated[i],linestyle='dashed',color='k',alpha=0.5,linewidth=1)
                    ctr+=1
                if second_evt_onset_translated[i] >= ts_reset_normalized[0] and second_evt_onset_translated[i] <= ts_reset_normalized[-1]:
                    if ctr1 == 0:  # label only first one for the legend
                        ax2.axvline(x=second_evt_onset_translated[i],linestyle='dashed',color='k',alpha=0.5,linewidth=1,label=event2_name)
                    else:
                        ax2.axvline(x=second_evt_onset_translated[i],linestyle='dashed',color='k',alpha=0.5,linewidth=1)
                    ctr1+=1
    
    # create title, axis labels, and legend
    my_title = 'Subject: ' + subject
    canvas.fig.suptitle(my_title, fontsize=FIGURE_TITLE_FONT_SIZE)
    ax.set_title('Downsampled',fontdict={'fontsize': TITLE_FONT_SIZE, 'fontweight': 'medium'})
    ax.set_xlabel('Time (sec)', fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_ylabel('mV', fontsize=AXIS_LABEL_FONTSIZE)
    ax.legend(loc=MY_LEGEND_LOC, bbox_to_anchor=MY_LEGENG_POS)
    # hide top and right border
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax2.set_title('Normalized',fontdict={'fontsize': TITLE_FONT_SIZE, 'fontweight': 'medium'})
    ax2.set_xlabel('Time (sec)', fontsize=AXIS_LABEL_FONTSIZE)
    ax2.set_ylabel('dF/F (%)', fontsize=AXIS_LABEL_FONTSIZE)
    if show_norm_as == "Z-Score":
        ax2.set_ylabel("Z-Score",fontsize=AXIS_LABEL_FONTSIZE)
    ax2.legend(loc=MY_LEGEND_LOC, bbox_to_anchor=MY_LEGENG_POS)
    # hide top and right border
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    # set the spacing between subplots 
    canvas.fig.subplots_adjust(left=MY_LEFT, 
                    bottom=MY_BOTTOM,  
                    right=MY_RIGHT,  
                    top=MY_TOP,  
                    wspace=MY_WSPACE,  
                    hspace=MY_HSPACE) 
    canvas.draw()
    
def plot_separate_only(canvas,subject,trimmed_signal_dict,downsampled_signal_dict):
    # if there is downsampled data plot that
    if len(downsampled_signal_dict) > 0:
        GCaMP_data = downsampled_signal_dict["signal"]
        control_data = downsampled_signal_dict["control"]
        ts = downsampled_signal_dict["ts"]
    else:
        # trimmed raw data
        GCaMP_data = trimmed_signal_dict["signal"]
        control_data = trimmed_signal_dict["control"]
        ts = trimmed_signal_dict["ts"]
        
    total_seconds = ts[-1]-ts[0]
    # start time from zero
    ts_reset = [i*total_seconds/len(ts) for i in range(len(ts))]
        
    # plot
    # clear previous figure
    canvas.fig.clf()
    ax = canvas.fig.add_subplot(211)
    ax2 = canvas.fig.add_subplot(212,sharex=ax)

    # plot the trimmed 
    ax.plot(ts_reset, GCaMP_data,
             color=SIGNAL_COLOR_RGB,
             linewidth=1,
             label = "signal")
    ax2.plot(ts_reset,control_data,
             color=CONTROL_COLOR_RGB,
             linewidth=1,
             label = "control")
    
    # create title, axis labels, and legend
    my_title = 'Subject: ' + subject
    canvas.fig.suptitle(my_title, fontsize=FIGURE_TITLE_FONT_SIZE)
    ax.set_xlabel('Time (sec)', fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_ylabel('mV', fontsize=AXIS_LABEL_FONTSIZE)
    ax.legend(loc=MY_LEGEND_LOC, bbox_to_anchor=MY_LEGENG_POS)
    # hide top and right border
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax2.set_xlabel('Time (sec)', fontsize=AXIS_LABEL_FONTSIZE)
    ax2.set_ylabel('mV', fontsize=AXIS_LABEL_FONTSIZE)
    ax2.legend(loc=MY_LEGEND_LOC, bbox_to_anchor=MY_LEGENG_POS)
    # hide top and right border
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    # show only integer ticks
    ax2.yaxis.set_major_locator(MaxNLocator(integer=True))
    # set the spacing between subplots 
    canvas.fig.subplots_adjust(left=MY_LEFT, 
                    bottom=MY_BOTTOM,  
                    right=MY_RIGHT,  
                    top=MY_TOP,  
                    wspace=MY_WSPACE,  
                    hspace=MY_HSPACE) 
    canvas.draw()
    
def plot_separate_with_event(canvas,subject,trimmed_signal_dict,downsampled_signal_dict,event_name,event2_name,event_data):
    # if there is downsampled data plot that
    if len(downsampled_signal_dict) > 0:
        GCaMP_data = downsampled_signal_dict["signal"]
        control_data = downsampled_signal_dict["control"]
        ts = downsampled_signal_dict["ts"]
    else:
        # trimmed raw data
        GCaMP_data = trimmed_signal_dict["signal"]
        control_data = trimmed_signal_dict["control"]
        ts = trimmed_signal_dict["ts"]
    total_seconds = ts[-1]-ts[0]
    # start time from zero
    ts_reset = [i*total_seconds/len(ts) for i in range(len(ts))]
    
    # plot
    # clear previous figure
    canvas.fig.clf()
    ax = canvas.fig.add_subplot(211)
    ax2 = canvas.fig.add_subplot(212,sharex=ax)

    # plot the trimmed 
    ax.plot(ts_reset, GCaMP_data,
             color=SIGNAL_COLOR_RGB,
             linewidth=1,
             label = "signal")
    ax2.plot(ts_reset,control_data,
             color=CONTROL_COLOR_RGB,
             linewidth=1,
             label = "control")
    
    # get accurate event data timing
    first_evt_onset = event_data[0][0]
    first_evt_onset_translated = [el-ts[0] for el in first_evt_onset]
    first_evt_offset = event_data[0][1]
    first_evt_offset_translated = [el-ts[0] for el in first_evt_offset]
    if len(event_data) > 1 and len(event_data[1][0]) > 0:
        second_evt_onset = event_data[1][0]
        second_evt_onset_translated = [el-ts[0] for el in second_evt_onset]
        second_evt_offset = event_data[1][1]
        second_evt_offset_translated = [el-ts[0] for el in second_evt_offset]
        
    if len(first_evt_onset) > 0:
        # first event
        # plot tone onsets as vertical lines
        # first element in tone_on_off is a list with onset times
        ctr = 0
        for i in range(len(first_evt_onset_translated)):
            if first_evt_onset_translated[i] >= ts_reset[0] and first_evt_onset_translated[i] <= ts_reset[-1]:
                if ctr == 0:  # label only first one for the legend
                    ax.axvline(x=first_evt_onset_translated[i],color='k',alpha=0.5,linewidth=1,label=event_name)
                    ax2.axvline(x=first_evt_onset_translated[i],color='k',alpha=0.5,linewidth=1,label=event_name)
                else:
                    ax.axvline(x=first_evt_onset_translated[i],color='k',alpha=0.5,linewidth=1)
                    ax2.axvline(x=first_evt_onset_translated[i],color='k',alpha=0.5,linewidth=1)
                ctr+=1
        # if there is more than one event, plot the second one as well
        if len(event_data) > 1 and len(event_data[1][0]) > 0:
            ctr = 0
            for i in range(len(second_evt_onset_translated)):
                if second_evt_onset_translated[i] >= ts_reset[0] and second_evt_onset_translated[i] <= ts_reset[-1]:
                    if ctr == 0:  # label only first one for the legend
                        ax.axvline(x=second_evt_onset_translated[i],linestyle='dashed',color='k',alpha=0.5,linewidth=1,label=event2_name)
                        ax2.axvline(x=second_evt_onset_translated[i],linestyle='dashed',color='k',alpha=0.5,linewidth=1,label=event2_name)
                    else:
                        ax.axvline(x=second_evt_onset_translated[i],linestyle='dashed',color='k',alpha=0.5,linewidth=1)
                        ax2.axvline(x=second_evt_onset_translated[i],linestyle='dashed',color='k',alpha=0.5,linewidth=1)
                    ctr+=1
                           
    # create title, axis labels, and legend
    my_title = 'Subject: ' + subject
    canvas.fig.suptitle(my_title, fontsize=FIGURE_TITLE_FONT_SIZE)
    ax.set_xlabel('Time (sec)', fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_ylabel('mV', fontsize=AXIS_LABEL_FONTSIZE)
    ax.legend(loc=MY_LEGEND_LOC, bbox_to_anchor=MY_LEGENG_POS)
    # hide top and right border
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax2.set_xlabel('Time (sec)', fontsize=AXIS_LABEL_FONTSIZE)
    ax2.set_ylabel('mV', fontsize=AXIS_LABEL_FONTSIZE)
    ax2.legend(loc=MY_LEGEND_LOC, bbox_to_anchor=MY_LEGENG_POS)
    # hide top and right border
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    # show only integer ticks
    ax2.yaxis.set_major_locator(MaxNLocator(integer=True))
    # set the spacing between subplots 
    canvas.fig.subplots_adjust(left=MY_LEFT, 
                    bottom=MY_BOTTOM,  
                    right=MY_RIGHT,  
                    top=MY_TOP,  
                    wspace=MY_WSPACE,  
                    hspace=MY_HSPACE) 
    canvas.draw()
    
def plot_separate_with_normalized(canvas,subject,trimmed_signal_dict,downsampled_signal_dict,normalized_dict,show_norm_as):
    # if there is downsampled data plot that
    if len(downsampled_signal_dict) > 0:
        GCaMP_data = downsampled_signal_dict["signal"]
        control_data = downsampled_signal_dict["control"]
        ts = downsampled_signal_dict["ts"]
    else:
        # trimmed raw data
        GCaMP_data = trimmed_signal_dict["signal"]
        control_data = trimmed_signal_dict["control"]
        ts = trimmed_signal_dict["ts"]
    # get normalized data to plot
    normalized_data = normalized_dict["normalized_signal"]
    ts_normalized = normalized_dict["ts"]
    
    total_seconds = ts[-1]-ts[0]
    # start time from zero
    ts_reset = [i*total_seconds/len(ts) for i in range(len(ts))]
    total_seconds_normalized = ts_normalized[-1]-ts_normalized[0]
    # start time from zero
    ts_reset_normalized = [i*total_seconds_normalized/len(ts_normalized) for i in range(len(ts_normalized))]
        
    # plot
    # clear previous figure
    canvas.fig.clf()
    ax = canvas.fig.add_subplot(311)
    ax2 = canvas.fig.add_subplot(312,sharex=ax)
    ax3 = canvas.fig.add_subplot(313,sharex=ax)

    ax.plot(ts_reset, GCaMP_data,
             color=SIGNAL_COLOR_RGB,
             linewidth=1,
             label = "signal")
    ax2.plot(ts_reset,control_data,
             color=CONTROL_COLOR_RGB,
             linewidth=1,
             label = "control")
    ax3.plot(ts_reset_normalized,normalized_data,
             color=SIGNAL_COLOR_RGB,
             linewidth=1,
             label = "normalized")
    # create title, axis labels, and legend
    my_title = 'Subject: ' + subject
    canvas.fig.suptitle(my_title, fontsize=FIGURE_TITLE_FONT_SIZE)
    ax.set_xlabel('Time (sec)', fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_ylabel('mV', fontsize=AXIS_LABEL_FONTSIZE)
    ax.legend(loc=MY_LEGEND_LOC, bbox_to_anchor=MY_LEGENG_POS)
    # hide top and right border
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax2.set_xlabel('Time (sec)', fontsize=AXIS_LABEL_FONTSIZE)
    ax2.set_ylabel('mV', fontsize=AXIS_LABEL_FONTSIZE)
    ax2.legend(loc=MY_LEGEND_LOC, bbox_to_anchor=MY_LEGENG_POS)
    # hide top and right border
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax3.set_ylabel('dF/F (%)', fontsize=AXIS_LABEL_FONTSIZE)
    if show_norm_as == "Z-Score":
        ax3.set_ylabel("Z-Score",fontsize=AXIS_LABEL_FONTSIZE)
    ax3.set_xlabel('Time (sec)', fontsize=AXIS_LABEL_FONTSIZE)
    ax3.legend(loc=MY_LEGEND_LOC, bbox_to_anchor=MY_LEGENG_POS)
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    # set the spacing between subplots 
    canvas.fig.subplots_adjust(left=MY_LEFT, 
                    bottom=MY_BOTTOM,  
                    right=MY_RIGHT,  
                    top=MY_TOP,  
                    wspace=MY_WSPACE,  
                    hspace=MY_HSPACE) 
    canvas.draw()
    
def plot_separate_with_normalized_with_event(canvas,subject,trimmed_signal_dict,downsampled_signal_dict,normalized_dict,show_norm_as,event_name,event2_name,event_data):
    # if there is downsampled data plot that
    if len(downsampled_signal_dict) > 0:
        GCaMP_data = downsampled_signal_dict["signal"]
        control_data = downsampled_signal_dict["control"]
        ts = downsampled_signal_dict["ts"]
    else:
        # trimmed raw data
        GCaMP_data = trimmed_signal_dict["signal"]
        control_data = trimmed_signal_dict["control"]
        ts = trimmed_signal_dict["ts"]
    # get normalized data to plot
    normalized_data = normalized_dict["normalized_signal"]
    ts_normalized = normalized_dict["ts"]
    
    total_seconds = ts[-1]-ts[0]
    # start time from zero
    ts_reset = [i*total_seconds/len(ts) for i in range(len(ts))]
    total_seconds_normalized = ts_normalized[-1]-ts_normalized[0]
    # start time from zero
    ts_reset_normalized = [i*total_seconds_normalized/len(ts_normalized) for i in range(len(ts_normalized))]
        
    # plot
    # clear previous figure
    canvas.fig.clf()
    ax = canvas.fig.add_subplot(311)
    ax2 = canvas.fig.add_subplot(312,sharex=ax)
    ax3 = canvas.fig.add_subplot(313,sharex=ax)

    ax.plot(ts_reset, GCaMP_data,
             color=SIGNAL_COLOR_RGB,
             linewidth=1,
             label = "signal")
    ax2.plot(ts_reset,control_data,
             color=CONTROL_COLOR_RGB,
             linewidth=1,
             label = "control")
    ax3.plot(ts_reset_normalized,normalized_data,
             color=SIGNAL_COLOR_RGB,
             linewidth=1,
             label = "normalized")
    
    # get accurate event data timing
    first_evt_onset = event_data[0][0]
    first_evt_onset_translated = [el-ts[0] for el in first_evt_onset]
    first_evt_offset = event_data[0][1]
    first_evt_offset_translated = [el-ts[0] for el in first_evt_offset]
    if len(event_data) > 1 and len(event_data[1][0]) > 0:
        second_evt_onset = event_data[1][0]
        second_evt_onset_translated = [el-ts[0] for el in second_evt_onset]
        second_evt_offset = event_data[1][1]
        second_evt_offset_translated = [el-ts[0] for el in second_evt_offset]
        
    if len(first_evt_onset) > 0:
        # first event
        # plot tone onsets as vertical lines
        # first element in tone_on_off is a list with onset times
        ctr = 0
        ctr1 = 0
        for i in range(len(first_evt_onset_translated)):
            if first_evt_onset_translated[i] >= ts_reset[0] and first_evt_onset_translated[i] <= ts_reset[-1]:
                if ctr == 0:  # label only first one for the legend
                    ax.axvline(x=first_evt_onset_translated[i],color='k',alpha=0.5,linewidth=1,label=event_name)
                    ax2.axvline(x=first_evt_onset_translated[i],color='k',alpha=0.5,linewidth=1,label=event_name)
                else:
                    ax.axvline(x=first_evt_onset_translated[i],color='k',alpha=0.5,linewidth=1)
                    ax2.axvline(x=first_evt_onset_translated[i],color='k',alpha=0.5,linewidth=1)
                ctr+=1
            if first_evt_onset_translated[i] >= ts_reset_normalized[0] and first_evt_onset_translated[i] <= ts_reset_normalized[-1]:
                if ctr1 == 0:  # label only first one for the legend
                    ax3.axvline(x=first_evt_onset_translated[i],color='k',alpha=0.5,linewidth=1,label=event_name)
                else:
                    ax3.axvline(x=first_evt_onset_translated[i],color='k',alpha=0.5,linewidth=1)
                ctr1+=1
        # if there is more than one event, plot the second one as well
        if len(event_data) > 1 and len(event_data[1][0]) > 0:
            ctr = 0
            ctr1 = 0
            for i in range(len(second_evt_onset_translated)):
                if second_evt_onset_translated[i] >= ts_reset[0] and second_evt_onset_translated[i] <= ts_reset[-1]:
                    if ctr == 0:  # label only first one for the legend
                        ax.axvline(x=second_evt_onset_translated[i],linestyle='dashed',color='k',alpha=0.5,linewidth=1,label=event2_name)
                        ax2.axvline(x=second_evt_onset_translated[i],linestyle='dashed',color='k',alpha=0.5,linewidth=1,label=event2_name)
                    else:
                        ax.axvline(x=second_evt_onset_translated[i],linestyle='dashed',color='k',alpha=0.5,linewidth=1)
                        ax2.axvline(x=second_evt_onset_translated[i],linestyle='dashed',color='k',alpha=0.5,linewidth=1)
                    ctr+=1
                if second_evt_onset_translated[i] >= ts_reset_normalized[0] and second_evt_onset_translated[i] <= ts_reset_normalized[-1]:
                    if ctr1 == 0:  # label only first one for the legend
                        ax3.axvline(x=second_evt_onset_translated[i],linestyle='dashed',color='k',alpha=0.5,linewidth=1,label=event2_name)
                    else:
                        ax3.axvline(x=second_evt_onset_translated[i],linestyle='dashed',color='k',alpha=0.5,linewidth=1)
                    ctr1+=1
    
    # create title, axis labels, and legend
    my_title = 'Subject: ' + subject
    canvas.fig.suptitle(my_title, fontsize=FIGURE_TITLE_FONT_SIZE)
    ax.set_xlabel('Time (sec)', fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_ylabel('mV', fontsize=AXIS_LABEL_FONTSIZE)
    ax.legend(loc=MY_LEGEND_LOC, bbox_to_anchor=MY_LEGENG_POS)
    # hide top and right border
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax2.set_xlabel('Time (sec)', fontsize=AXIS_LABEL_FONTSIZE)
    ax2.set_ylabel('mV', fontsize=AXIS_LABEL_FONTSIZE)
    ax2.legend(loc=MY_LEGEND_LOC, bbox_to_anchor=MY_LEGENG_POS)
    # hide top and right border
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax3.set_ylabel('dF/F (%)', fontsize=AXIS_LABEL_FONTSIZE)
    if show_norm_as == "Z-Score":
        ax3.set_ylabel("Z-Score",fontsize=AXIS_LABEL_FONTSIZE)
    ax3.set_xlabel('Time (sec)', fontsize=AXIS_LABEL_FONTSIZE)
    ax3.legend(loc=MY_LEGEND_LOC, bbox_to_anchor=MY_LEGENG_POS)
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    # set the spacing between subplots 
    canvas.fig.subplots_adjust(left=MY_LEFT, 
                    bottom=MY_BOTTOM,  
                    right=MY_RIGHT,  
                    top=MY_TOP,  
                    wspace=MY_WSPACE,  
                    hspace=MY_HSPACE) 
    canvas.draw()
    
# plot raw but downsampled perievents
def plot_raw_perievents(canvas,subject,modified_data,current_trials,perievent_options_dict,settings_dict,signal_name,control_name,export,save_plots,group_name,export_loc_data):
    settings_dict[0]["subject"] = subject
    settings_dict[0]["subject_group_name"] = group_name
    settings_dict[0]["baseline_from_sec"] =perievent_options_dict["baseline_from"]
    settings_dict[0]["baseline_to_sec"] =perievent_options_dict["baseline_to"]
    settings_dict[0]["auc_pre_from"] =perievent_options_dict["auc_pre_from"]
    settings_dict[0]["auc_pre_to"] =perievent_options_dict["auc_pre_to"]
    settings_dict[0]["auc_post_from"] =perievent_options_dict["auc_post_from"]
    settings_dict[0]["auc_post_to"] =perievent_options_dict["auc_post_to"]
    show_norm_as = settings_dict[0]["show_norm_as"]
    dump_path,file_beginning = export_loc_data
    raw_df = None
    event_name = perievent_options_dict["event_name"] if len(perievent_options_dict["event_name"]) > 0 else perievent_options_dict["event"]
    GCaMP_perievent_data = modified_data.streams[signal_name].filtered_downsampled
    control_perievent_data = modified_data.streams[control_name].filtered_downsampled
    if len(current_trials) > 0:
        # only those trials that are selected
        # print(current_trials)
        GCaMP_perievent_data = [GCaMP_perievent_data[trial-1] for trial in current_trials]
        control_perievent_data = [control_perievent_data[trial-1] for trial in current_trials]
    
    N = settings_dict[0]["downsample"]
    # !!! Create times in sec first to get last time value not total samples!!!
    t1From0 = np.linspace(1/N,(perievent_options_dict["sec_before"]+perievent_options_dict["sec_after"]),(perievent_options_dict["sec_before"]+perievent_options_dict["sec_after"])*N)
    # print(t1From0)
    ts1 = -perievent_options_dict["sec_before"] + t1From0
    # print(ts1)
    ts2 = ts1

    # ts1 = -perievent_options_dict["sec_before"] + np.linspace(1, len(GCaMP_perievent_data[0]), len(GCaMP_perievent_data[0]))/modified_data.streams[signal_name].fs*N
    # ts2 = -perievent_options_dict["sec_before"] + np.linspace(1, len(control_perievent_data[0]), len(control_perievent_data[0]))/modified_data.streams[control_name].fs*N
    total_plots = len(GCaMP_perievent_data) # total number of events
    # allow to show only up till 9 plots
    if total_plots > 9:
        total_plots = 9
    
    if settings_dict[0]["filter"] == True:
        print("Start smoothing",settings_dict[0]["filter_window"])
        # smooth data using filter
        temp_signal_smoothed = []
        temp_control_smoothed = []
        a = 1
        b = np.divide(np.ones((settings_dict[0]["filter_window"],)), settings_dict[0]["filter_window"])
        for i in range(len(GCaMP_perievent_data)):
            single_evt = filtfilt(b, a, GCaMP_perievent_data[i])
            # single_evt = gaussian_filter(GCaMP_perievent_data[i], sigma=settings_dict[0]["filter_fraction"])
#            single_evt = lowess.lowess(pd.Series(ts_signal4average), pd.Series(GCaMP_perievent_data[i]), 
#                                       bandwidth=settings_dict[0]["filter_fraction"], polynomialDegree=1)
            temp_signal_smoothed.append(single_evt)
        for i in range(len(control_perievent_data)):
            single_evt = filtfilt(b, a, control_perievent_data[i])
            # single_evt = gaussian_filter(control_perievent_data[i], sigma=settings_dict[0]["filter_fraction"])
#            single_evt = lowess.lowess(pd.Series(ts_control4average), pd.Series(control_perievent_data[i]), 
#                                       bandwidth=settings_dict[0]["filter_fraction"], polynomialDegree=1)
            temp_control_smoothed.append(single_evt)
        if len(temp_signal_smoothed) > 0 and len(temp_control_smoothed) > 0:
            GCaMP_perievent_data = temp_signal_smoothed
            control_perievent_data = temp_control_smoothed
        print("Done smoothing")
        
    # normalize
    y_dff_all = []
    # find out how to normalize
    if settings_dict[0]["normalization"] == 'Standard Polynomial Fitting':
        # https://github.com/djamesbarker/pMAT
        for x, y in zip(control_perievent_data, GCaMP_perievent_data):
            x = np.array(x)
            y = np.array(y)
            # Before removing +/- 2xstdev
        #     # https://stackoverflow.com/questions/45338872/matlab-polyval-function-with-three-outputs-equivalent-in-python-numpy
        #     mu = np.mean(x)
        #     std = np.std(x, ddof=0)
        #     # Call np.polyfit(), using the shifted and scaled version of control_arr
        #     cscaled = np.polynomial.polynomial.Polynomial.fit((x - mu)/std, y, 1)
        #     # Create a poly1d object that can be called
        #     # https://numpy.org/doc/stable/reference/routines.polynomials.html
        #     pscaled = Polynomial(cscaled.convert().coef)
        #     # Inputs to pscaled must be shifted and scaled using mu and std
        #     F0 = pscaled((x - mu)/std)
        # #    print("F0?",F0[:20])
        #     dffnorm = (y - F0)/F0 * 100
        #############################################################################################
            # filter out signal values that are below or above 2 standard deviations from the signal mean 
            mean_signal = np.mean(y)
            stdev_signal = np.std(y)
            indexes = np.where((y<mean_signal+2*stdev_signal) & (y>mean_signal-2*stdev_signal))
            selected_signal = y[indexes]
            selected_control = x[indexes]
            mu = np.mean(selected_control)
            std = np.std(selected_control)
            cscaled = np.polynomial.polynomial.Polynomial.fit((selected_control - mu)/std, selected_signal, 1)
            pscaled = Polynomial(cscaled.convert().coef)
            F0 = pscaled((x - mu)/std)
            dffnorm = (y - F0)/F0 * 100
        #############################################################################################################################
            # find all values of the normalized DF/F that are negative so you can next shift up the curve 
            # to make 0 the mean value for DF/F
            negative = dffnorm[dffnorm<0]
            dff=dffnorm-np.mean(negative)

            if show_norm_as == "Z-Score":
                median_all = np.median(dff)
                mad = stats.median_abs_deviation(dff)
                dff = (dff - median_all)/mad

            y_dff_all.append(dff)
            
    elif settings_dict[0]["normalization"] == 'Modified Polynomial Fitting':
        for x, y in zip(control_perievent_data, GCaMP_perievent_data):
            x = np.array(x)
            y = np.array(y)
            # Before removing +/- 2xstdev
            # bls_signal = np.polynomial.polynomial.Polynomial.fit(ts1, y, 1)
            # F0_signal = polyval(ts1,bls_signal.convert().coef)
            # dFF_signal = (y - F0_signal)/F0_signal *100
            
            # bls_control = np.polynomial.polynomial.Polynomial.fit(ts2,x,1)
            # F0_control = polyval(ts2,bls_control.convert().coef)
            # dFF_control = (x - F0_control)/F0_control *100
            # dFFnorm = dFF_signal - dFF_control
            ############################################################################################
            # 23/01 filter out signal values that are below or above 2 standard deviations from the signal mean 
            mean_signal = np.mean(y)
            stdev_signal = np.std(y)
            mean_control = np.mean(x)
            stdev_control = np.std(x)
            indexes_signal = np.where((y<mean_signal+2*stdev_signal) & (y>mean_signal-2*stdev_signal))
            selected_signal = y[indexes_signal]
            selected_ts_signal_arr = ts1[indexes_signal]
            indexes_control = np.where((x<mean_control+2*stdev_control) & (x>mean_control-2*stdev_control))
            selected_control = x[indexes_control]
            selected_ts_control_arr = ts2[indexes_control]
            # fit time axis to the 465nm stream  
            bls_signal = np.polynomial.polynomial.Polynomial.fit(selected_ts_signal_arr,selected_signal,1)
            bls_control = np.polynomial.polynomial.Polynomial.fit(selected_ts_control_arr,selected_control,1)
            F0_signal = polyval(ts1,bls_signal.convert().coef)
            F0_control = polyval(ts2,bls_control.convert().coef)
            # dF/F for the 465 channel
            dFF_signal = (y - F0_signal)/F0_signal *100
            # dF/F for the 405 channel
            dFF_control = (x - F0_control)/F0_control *100
            dFFnorm = dFF_signal - dFF_control
            ######################################################################################
            # find all values of the normalized DF/F that are negative so you can next shift up the curve 
            # to make 0 the mean value for DF/F
            negative = dFFnorm[dFFnorm<0]
            dff = dFFnorm-np.mean(negative)

            if show_norm_as == "Z-Score":
                median_all = np.median(dff)
                mad = stats.median_abs_deviation(dff)
                dff = (dff - median_all)/mad

            y_dff_all.append(dff)
    # print(f"How many trials in preview plot perievent: {len(y_dff_all)}")
    # plot
    # clear previous figure
    canvas.fig.clf()
    
    my_title = 'Subject: ' + subject
    canvas.fig.suptitle(my_title, fontsize=FIGURE_TITLE_FONT_SIZE)
    
    # plot 3 in the same row. Idea for layout from:
    # https://towardsdatascience.com/customizing-multiple-subplots-in-matplotlib-a3e1c2e099bc
    coord = []# create coord array from 231, 232, 233, ..., 236 
    column = 3
    if total_plots == 1: # if there is only single event, plot in one column
        column = 1
    if total_plots == 2: # if there are just two, plot in two columns
        column = 2
    for i in range(1, total_plots+1): 
        row = math.ceil(total_plots/3)
        coord.append(str(row)+str(column)+str(i))
    first_ax = None
    # create subplots
    for i in range(len(coord)):
        if i == 0: # remember for sharing axis
            ax = canvas.fig.add_subplot(int(coord[i]))
            first_ax = ax
        else:
            ax = canvas.fig.add_subplot(int(coord[i]),sharex=first_ax, sharey=first_ax)

        ax.plot(ts1,y_dff_all[i],color=SIGNAL_COLOR_RGB,linewidth=1,label = "normalized\nsignal")
        # plot event as vertical line
        p3=ax.axvline(x=0, linewidth=2, color='k', alpha=0.5,label=event_name)
        # hide top and right border
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        if len(current_trials) > 0:
            ax.set_title('Trial: '+str(current_trials[i]),fontdict={'fontsize': TITLE_FONT_SIZE, 'fontweight': 'medium'})
        else:
            ax.set_title('Trial: '+str(i+1),fontdict={'fontsize': TITLE_FONT_SIZE, 'fontweight': 'medium'})
        ax.set_xlabel('Time (sec)', fontsize=AXIS_LABEL_FONTSIZE)
        ax.set_ylabel('dF/F (%)', fontsize=AXIS_LABEL_FONTSIZE)
        if show_norm_as == "Z-Score":
            ax.set_ylabel('Z-Score', fontsize=AXIS_LABEL_FONTSIZE)
    ax.legend(loc=MY_LEGEND_LOC, bbox_to_anchor=MY_LEGENG_POS)
        
    # set the spacing between subplots 
    canvas.fig.subplots_adjust(left=MY_LEFT, 
                    bottom=MY_BOTTOM,  
                    right=MY_RIGHT,  
                    top=MY_TOP,  
                    wspace=MY_WSPACE,  
                    hspace=MY_HSPACE)   
    # create a list of dataframes to join later
    dfs = []
    time_df = pd.DataFrame({"Time (sec)":ts1})
    dfs.append(time_df)
    if show_norm_as == 'Z-Score':
        norm_type = " (z-score)"
    else:
        norm_type = " (%df/F)"
    # for i in range(total_plots):
    for i in range(len(y_dff_all)):
        # header_signal = "Trial"+str(i+1)+norm_type
        header_signal = "Trial"+str(current_trials[i])+norm_type
        signal_data = y_dff_all[i]
        df = pd.DataFrame({header_signal:signal_data})
        dfs.append(df)
    raw_df= pd.concat(dfs,axis=1)
    # create subfolder in current subject folder
    if export == True:
        file_name = file_beginning+"_"+event_name+"_all_perievent_normalized.csv"
        if file_beginning != subject:
            file_name = file_beginning+"_"+subject+ "_"+event_name+"_all_perievent_normalized.csv"
        dump_file_path = os.path.join(dump_path,file_name)
        settings_df = get_settings_df(settings_dict)
        settings_df.to_csv(dump_file_path,header=False)            
        with open(dump_file_path,'a') as f:
            f.write("\n")
        raw_df.to_csv(dump_file_path, index=False,mode='a')
  
    if export == True and save_plots == True:
        plot_file_name = file_beginning+"_"+ event_name + "_all_perievent_normalized.png"
        plot_file_name2 = file_beginning+"_"+ event_name + "_all_perievent_normalized.svg"
        if file_beginning != subject:
           plot_file_name = file_beginning+"_"+ subject+"_"+event_name + "_all_perievent_normalized.png" 
           plot_file_name2 = file_beginning+"_"+ subject+"_"+event_name + "_all_perievent_normalized.svg" 
        dump_plot_file_path = os.path.join(dump_path,plot_file_name)
        canvas.fig.savefig(dump_plot_file_path)
        # also as svg
        dump_plot_file_path = os.path.join(dump_path,plot_file_name2)
        canvas.fig.savefig(dump_plot_file_path, format='svg', dpi=DPI4SVG)
    else:
        canvas.draw()
    return raw_df
        
    

def plot_perievent_average_alone(canvas,subject,current_trials,perievent_options_dict,analyzed_perievent_dict,export,save_plots,group_name,settings_dict,signal_name,control_name,export_loc_data):
    settings_dict[0]["subject"] = subject
    settings_dict[0]["subject_group_name"] = group_name
    settings_dict[0]["baseline_from_sec"] =perievent_options_dict["baseline_from"]
    settings_dict[0]["baseline_to_sec"] =perievent_options_dict["baseline_to"]
    settings_dict[0]["auc_pre_from"] =perievent_options_dict["auc_pre_from"]
    settings_dict[0]["auc_pre_to"] =perievent_options_dict["auc_pre_to"]
    settings_dict[0]["auc_post_from"] =perievent_options_dict["auc_post_from"]
    settings_dict[0]["auc_post_to"] =perievent_options_dict["auc_post_to"]
    dump_path,file_beginning = export_loc_data
    event_name = perievent_options_dict["event_name"] if len(perievent_options_dict["event_name"]) > 0 else perievent_options_dict["event"]
    # get data to plot
    ts_signal = analyzed_perievent_dict["average"]["ts_signal"]
    ts_control = analyzed_perievent_dict["average"]["ts_control"]
    signal = analyzed_perievent_dict["average"]["signal"]
    control = analyzed_perievent_dict["average"]["control"]
    std_signal = analyzed_perievent_dict["average"]["std_signal"]
    std_control = analyzed_perievent_dict["average"]["std_control"]
    
    # plot
    # clear previous figure
    canvas.fig.clf()
    ax = canvas.fig.add_subplot(111)
    ax.plot(ts_signal, signal,
             color=SIGNAL_COLOR_RGB,
             linewidth=1,
             label = "signal")
    ax.plot(ts_control, control,
             color=CONTROL_COLOR_RGB,
             linewidth=1,
             label = "control")
    # plot a vertical line at t=0
    ax.axvline(x=0, linewidth=2, color='k', alpha=0.3,label=event_name)
    # plot standard error bands
    ax.fill_between(ts_signal, signal+std_signal, signal-std_signal,
                      facecolor=SIGNAL_COLOR_RGB, alpha=0.1)
    ax.fill_between(ts_control, control+std_control, control-std_control,
                      facecolor=CONTROL_COLOR_RGB, alpha=0.1)
    
    my_title = 'Average From '+ str(len(current_trials))+' Trials. Subject: ' + subject
    canvas.fig.suptitle(my_title, fontsize=FIGURE_TITLE_FONT_SIZE)
    # hide top and right border
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel('Time (sec)', fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_ylabel('mV', fontsize=AXIS_LABEL_FONTSIZE)
    ax.legend(loc=MY_LEGEND_LOC, bbox_to_anchor=MY_LEGENG_POS)
    
    # create subfolder in current subject folder
    if export == True:
        file_name = file_beginning+"_"+ event_name+"_average_perievent.csv"
        if file_beginning != subject:
            file_name = file_beginning+"_"+ subject+"_"+event_name+"_average_perievent.csv"
        dump_file_path = os.path.join(dump_path,file_name)
        sig = "Signal_"+signal_name
        cont = "Control_"+control_name
        raw_df= pd.DataFrame({"Time (sec)":ts_signal,
                                 sig +' (mV)': signal, 
                                 cont +' (mV)': control,
                                 "Error_signal": std_signal,
                                 "Error_control": std_control})
        settings_df = get_settings_df(settings_dict)
        settings_df.to_csv(dump_file_path,header=False)            
        with open(dump_file_path,'a') as f:
            f.write("\n")
        raw_df.to_csv(dump_file_path, index=False,mode='a')  
    if export == True and save_plots == True:
        plot_file_name = file_beginning+"_"+ event_name+"_average_perievent.png"
        plot_file_name2 = file_beginning+"_"+ event_name+"_average_perievent.svg"
        if file_beginning != subject:
            plot_file_name = file_beginning+"_"+ subject+"_"+event_name+"_average_perievent.png"
            plot_file_name2 = file_beginning+"_"+ subject+"_"+event_name+"_average_perievent.svg"
        dump_plot_file_path = os.path.join(dump_path,plot_file_name)
        canvas.fig.savefig(dump_plot_file_path)
        # also as svg
        dump_plot_file_path = os.path.join(dump_path,plot_file_name2)
        canvas.fig.savefig(dump_plot_file_path, format='svg', dpi=DPI4SVG)
    else:
        canvas.draw()
    
    
def plot_perievent_zscore_alone(canvas,subject, current_trials, perievent_options_dict,analyzed_perievent_dict,export,save_plots,group_name,settings_dict,export_loc_data):
    # print(f"Trials for auc from z-score: {current_trials}")
    settings_dict[0]["subject"] = subject
    settings_dict[0]["subject_group_name"] = group_name
    settings_dict[0]["baseline_from_sec"] =perievent_options_dict["baseline_from"]
    settings_dict[0]["baseline_to_sec"] =perievent_options_dict["baseline_to"]
    settings_dict[0]["auc_pre_from"] =perievent_options_dict["auc_pre_from"]
    settings_dict[0]["auc_pre_to"] =perievent_options_dict["auc_pre_to"]
    settings_dict[0]["auc_post_from"] =perievent_options_dict["auc_post_from"]
    settings_dict[0]["auc_post_to"] =perievent_options_dict["auc_post_to"]
    dump_path,file_beginning = export_loc_data    
    event_name = perievent_options_dict["event_name"] if len(perievent_options_dict["event_name"]) > 0 else perievent_options_dict["event"]
    raw_df = None
    # get data to plot
    ts = analyzed_perievent_dict["zscore"]["ts"]
    zscore_all = analyzed_perievent_dict["zscore"]["zscored"]
    zerror = analyzed_perievent_dict["zscore"]["zerror"]
    
    # plot
    # clear previous figure
    canvas.fig.clf()
    ax = canvas.fig.add_subplot(211)
    ax2 = canvas.fig.add_subplot(212,sharex=ax)

    # Heat Map based on z score of control fit subtracted signal
    cs = ax.imshow(zscore_all, cmap=plt.cm.bwr, interpolation='none', aspect="auto",
                    extent=[-perievent_options_dict['sec_before'], perievent_options_dict['sec_after'], 
                            len(zscore_all),0])
    cb = canvas.fig.colorbar(cs, ax=ax,pad=0.01, fraction=0.02)
    # plot the z-score trace for the signal with std error bands
    mean_zscore = np.mean(zscore_all, axis=0)
    ax2.plot(ts, mean_zscore, linewidth=2, color=SIGNAL_COLOR_RGB, label='Mean Z-Score')
    ax2.fill_between(ts, mean_zscore+zerror,
                          mean_zscore-zerror, facecolor=SIGNAL_COLOR_RGB, alpha=0.2, label="Standard\nerror")
    ax2.axvline(x=0, linewidth=2, color='k', alpha=0.3, label=event_name)
  
    my_title = 'Subject: ' + subject
    canvas.fig.suptitle(my_title, fontsize=FIGURE_TITLE_FONT_SIZE)
    ax.set_title('Individual z-Score Traces',fontsize=TITLE_FONT_SIZE)
    ax.set_ylabel('Trials',fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_xlabel('Seconds from Event Onset',fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_yticks(np.arange(0.5,len(zscore_all), 1))
    ax.set_yticklabels(np.arange(1, len(zscore_all)+1, 1))
    # hide top and right border
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.set_ylabel('z-Score',fontsize=AXIS_LABEL_FONTSIZE)
    ax2.set_xlabel('Seconds',fontsize=AXIS_LABEL_FONTSIZE)
    ax2.set_xlim([-perievent_options_dict['sec_before'], perievent_options_dict['sec_after']])
    ax2.legend(loc=MY_LEGEND_LOC, bbox_to_anchor=MY_LEGENG_POS)

    # bounds returns (x0,y0,width,height)
    # print("heatmap:")
    x_heatmap,y_heatmap,width_heatmap,height_heatmap = ax.get_position().bounds
    # print(ax.get_position().bounds)
    # print("zscore:")
    x_zscore,y_zscore,width_zscore,height_zscore = ax2.get_position().bounds
    # print(ax2.get_position().bounds)
    # readjust the position of zscore to match the sice of heatmap
    ax2.set_position([x_zscore,y_zscore,width_heatmap,height_zscore])
    
    # set the spacing between subplots 
    # canvas.fig.subplots_adjust(left=MY_LEFT, 
    #                 bottom=MY_BOTTOM,  
    #                 right=MY_RIGHT,  
    #                 top=MY_TOP,  
    #                 wspace=MY_WSPACE,  
    #                 hspace=MY_HSPACE) 
    
    # create a list of dataframes to join later
    dfs = []
    # add time
    time_df = pd.DataFrame({"Time (sec)":ts})
    dfs.append(time_df)
    for i in range(len(current_trials)):
        header = "Trial"+str(current_trials[i])+"_zscore"
        data = zscore_all[i]
        df = pd.DataFrame({header:data})
        dfs.append(df)
    # add mean
    mean_df = pd.DataFrame({"Mean_zscore":mean_zscore,
                                "Error":zerror})
    dfs.append(mean_df)
    raw_df= pd.concat(dfs,axis=1)
#    print("In zscore")
#    print(raw_df)
    # create subfolder in current subject folder
    if export == True:
        file_name = file_beginning+"_"+ event_name + "_perievent_zscore.csv"
        if file_beginning != subject:
            file_name = file_beginning+"_"+ subject+"_"+event_name + "_perievent_zscore.csv"
        dump_file_path = os.path.join(dump_path,file_name)
        settings_df = get_settings_df(settings_dict)
        settings_df.to_csv(dump_file_path,header=False)            
        with open(dump_file_path,'a') as f:
            f.write("\n")
        raw_df.to_csv(dump_file_path, index=False,mode='a') 
            
    if export == True and save_plots == True:
        plot_file_name = file_beginning+"_"+ event_name + "_perievent_zscore.png"
        plot_file_name2 = file_beginning+"_"+ event_name + "_perievent_zscore.svg"
        if file_beginning != subject:
            plot_file_name = file_beginning+"_"+ subject+"_"+event_name + "_perievent_zscore.png"
            plot_file_name2 = file_beginning+"_"+ subject+"_"+event_name + "_perievent_zscore.svg"
        dump_plot_file_path = os.path.join(dump_path,plot_file_name)
        canvas.fig.savefig(dump_plot_file_path)
        # also as svg
        dump_plot_file_path = os.path.join(dump_path,plot_file_name2)
        canvas.fig.savefig(dump_plot_file_path, format='svg', dpi=DPI4SVG)
    else:
        canvas.draw()
    return raw_df

def plot_perievent_zscore_with_trials_alone(canvas,subject, current_trials, perievent_options_dict,analyzed_perievent_dict,export,save_plots,group_name,settings_dict,export_loc_data):
    # print(f"Trials for auc from z-score: {current_trials}")
    settings_dict[0]["subject"] = subject
    settings_dict[0]["subject_group_name"] = group_name
    settings_dict[0]["baseline_from_sec"] =perievent_options_dict["baseline_from"]
    settings_dict[0]["baseline_to_sec"] =perievent_options_dict["baseline_to"]
    settings_dict[0]["auc_pre_from"] =perievent_options_dict["auc_pre_from"]
    settings_dict[0]["auc_pre_to"] =perievent_options_dict["auc_pre_to"]
    settings_dict[0]["auc_post_from"] =perievent_options_dict["auc_post_from"]
    settings_dict[0]["auc_post_to"] =perievent_options_dict["auc_post_to"]
    dump_path,file_beginning = export_loc_data    
    event_name = perievent_options_dict["event_name"] if len(perievent_options_dict["event_name"]) > 0 else perievent_options_dict["event"]
    raw_df = None
    # get data to plot
    ts = analyzed_perievent_dict["zscore"]["ts"]
    zscore_all = analyzed_perievent_dict["zscore"]["zscored"]
    zerror = analyzed_perievent_dict["zscore"]["zerror"]
    
    # plot
    # clear previous figure
    canvas.fig.clf()
    ax = canvas.fig.add_subplot(211)
    ax2 = canvas.fig.add_subplot(212,sharex=ax)
   
    # Heat Map based on z score of control fit subtracted signal
    cs = ax.imshow(zscore_all, cmap=plt.cm.bwr, interpolation='none', aspect="auto",
                    extent=[-perievent_options_dict['sec_before'], perievent_options_dict['sec_after'], 
                            len(zscore_all),0])
    canvas.fig.colorbar(cs, ax=ax,pad=0.01, fraction=0.02)
    # plot the z-score trace for the signal with std error bands
    mean_zscore = np.mean(zscore_all, axis=0)
    ax2.plot(ts, mean_zscore, linewidth=2, color=SIGNAL_COLOR_RGB, label='Mean Z-Score')
    # plot all trials
    for i in range(len(zscore_all)):
        if i == 0:
            ax2.plot(ts, zscore_all[i], linewidth=0.5, color=SIGNAL_COLOR_RGB, alpha=0.2, label='Trials')
        else:
            ax2.plot(ts, zscore_all[i], linewidth=0.5, color=SIGNAL_COLOR_RGB, alpha=0.2)

    ax2.axvline(x=0, linewidth=2, color='k', alpha=0.3, label=event_name)
  
    my_title = 'Subject: ' + subject
    canvas.fig.suptitle(my_title, fontsize=FIGURE_TITLE_FONT_SIZE)
    ax.set_title('Individual z-Score Traces',fontsize=TITLE_FONT_SIZE)
    ax.set_ylabel('Trials',fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_xlabel('Seconds from Event Onset',fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_yticks(np.arange(0.5,len(zscore_all), 1))
    ax.set_yticklabels(np.arange(1, len(zscore_all)+1, 1))
    # hide top and right border
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.set_ylabel('z-Score',fontsize=AXIS_LABEL_FONTSIZE)
    ax2.set_xlabel('Seconds',fontsize=AXIS_LABEL_FONTSIZE)
    ax2.legend(loc=MY_LEGEND_LOC, bbox_to_anchor=MY_LEGENG_POS)
    ax2.set_xlim([-perievent_options_dict['sec_before'], perievent_options_dict['sec_after']])

    # bounds returns (x0,y0,width,height)
    # print("heatmap:")
    x_heatmap,y_heatmap,width_heatmap,height_heatmap = ax.get_position().bounds
    # print(ax.get_position().bounds)
    # print("zscore:")
    x_zscore,y_zscore,width_zscore,height_zscore = ax2.get_position().bounds
    # print(ax2.get_position().bounds)
    # readjust the position of zscore to match the sice of heatmap
    ax2.set_position([x_zscore,y_zscore,width_heatmap,height_zscore])
    
    # # set the spacing between subplots 
    # canvas.fig.subplots_adjust(left=MY_LEFT, 
    #                 bottom=MY_BOTTOM,  
    #                 right=MY_RIGHT,  
    #                 top=MY_TOP,  
    #                 wspace=MY_WSPACE,  
    #                 hspace=MY_HSPACE) 
    # create a list of dataframes to join later
    dfs = []
    # add time
    time_df = pd.DataFrame({"Time (sec)":ts})
    dfs.append(time_df)
    for i in range(len(current_trials)):
        header = "Trial"+str(current_trials[i])+"_zscore"
        data = zscore_all[i]
        df = pd.DataFrame({header:data})
        dfs.append(df)
    # add mean
    mean_df = pd.DataFrame({"Mean_zscore":mean_zscore,
                                "Error":zerror})
    dfs.append(mean_df)
    raw_df= pd.concat(dfs,axis=1)
#    print("In zscore")
#    print(raw_df)
    # create subfolder in current subject folder
    if export == True:
        file_name = file_beginning+"_"+ event_name + "_perievent_zscore.csv"
        if file_beginning != subject:
            file_name = file_beginning+"_"+ subject+"_"+event_name + "_perievent_zscore.csv"
        dump_file_path = os.path.join(dump_path,file_name)
        settings_df = get_settings_df(settings_dict)
        settings_df.to_csv(dump_file_path,header=False)            
        with open(dump_file_path,'a') as f:
            f.write("\n")
        raw_df.to_csv(dump_file_path, index=False,mode='a') 
            
    if export == True and save_plots == True:
        plot_file_name = file_beginning+"_"+ event_name + "_perievent_zscore.png"
        plot_file_name2 = file_beginning+"_"+ event_name + "_perievent_zscore.svg"
        if file_beginning != subject:
            plot_file_name = file_beginning+"_"+ subject+"_"+event_name + "_perievent_zscore.png"
            plot_file_name2 = file_beginning+"_"+ subject+"_"+event_name + "_perievent_zscore.svg"
        dump_plot_file_path = os.path.join(dump_path,plot_file_name)
        canvas.fig.savefig(dump_plot_file_path)
        # also as svg
        dump_plot_file_path = os.path.join(dump_path,plot_file_name2)
        canvas.fig.savefig(dump_plot_file_path, format='svg', dpi=DPI4SVG)
    else:
        canvas.draw()
    return raw_df    
    
def plot_perievent_avg_zscore(canvas,subject,current_trials, perievent_options_dict,analyzed_perievent_dict):
    event_name = perievent_options_dict["event_name"] if len(perievent_options_dict["event_name"]) > 0 else perievent_options_dict["event"]
    
    # get data to plot avg
    ts_signal = analyzed_perievent_dict["average"]["ts_signal"]
    ts_control = analyzed_perievent_dict["average"]["ts_control"]
    signal = analyzed_perievent_dict["average"]["signal"]
    control = analyzed_perievent_dict["average"]["control"]
    std_signal = analyzed_perievent_dict["average"]["std_signal"]
    std_control = analyzed_perievent_dict["average"]["std_control"]
    # zscore
    ts = analyzed_perievent_dict["zscore"]["ts"]
    zscore_all = analyzed_perievent_dict["zscore"]["zscored"]
    zerror = analyzed_perievent_dict["zscore"]["zerror"]
    
    # plot
    # clear previous figure
    canvas.fig.clf()
    ax = canvas.fig.add_subplot(311)
    ax2 = canvas.fig.add_subplot(312,sharex=ax)
    ax3 = canvas.fig.add_subplot(313,sharex=ax)
    
    # average
    ax.plot(ts_signal, signal,
             color=SIGNAL_COLOR_RGB,
             linewidth=1,
             label = "signal")
    ax.plot(ts_control, control,
             color=CONTROL_COLOR_RGB,
             linewidth=1,
             label = "control")
    # plot a vertical line at t=0
    ax.axvline(x=0, linewidth=2, color='k', alpha=0.3,label=event_name)
    # plot standard error bands
    ax.fill_between(ts_signal, signal+std_signal, signal-std_signal,
                      facecolor=SIGNAL_COLOR_RGB, alpha=0.1)
    ax.fill_between(ts_control, control+std_control, control-std_control,
                      facecolor=CONTROL_COLOR_RGB, alpha=0.1)
    
    # zscore
    # Heat Map based on z score of control fit subtracted signal
    cs = ax2.imshow(zscore_all, cmap=plt.cm.bwr, interpolation='none', aspect="auto",
                    extent=[-perievent_options_dict['sec_before'], perievent_options_dict['sec_after'], 
                            len(zscore_all),0])
    canvas.fig.colorbar(cs, ax=ax2,pad=0.01, fraction=0.02)
    # plot the z-score trace for the signal with std error bands
    ax3.plot(ts, np.mean(zscore_all, axis=0), linewidth=2, color=SIGNAL_COLOR_RGB, label='Mean Z-Score')
    ax3.fill_between(ts, np.mean(zscore_all, axis=0)+zerror
                          ,np.mean(zscore_all, axis=0)-zerror, facecolor=SIGNAL_COLOR_RGB, alpha=0.2, label="Standard\nerror")
    ax3.axvline(x=0, linewidth=2, color='k', alpha=0.3,label=event_name)
    
    my_title = 'Subject: ' + subject
    canvas.fig.suptitle(my_title, fontsize=FIGURE_TITLE_FONT_SIZE)
    ax.set_title('Average From '+str(len(current_trials))+' Trials',fontsize=TITLE_FONT_SIZE)
    # hide top and right border
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel('Time (sec)', fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_ylabel('mV', fontsize=AXIS_LABEL_FONTSIZE)
    ax.legend(loc=MY_LEGEND_LOC, bbox_to_anchor=MY_LEGENG_POS)
    ax.set_xlim([-perievent_options_dict['sec_before'], perievent_options_dict['sec_after']])
    ax2.set_title('Individual z-Score Traces',fontsize=TITLE_FONT_SIZE)
    ax2.set_ylabel('Trials',fontsize=AXIS_LABEL_FONTSIZE)
    ax2.set_xlabel('Seconds from Event Onset',fontsize=AXIS_LABEL_FONTSIZE)
    ax2.set_yticks(np.arange(0.5,len(zscore_all), 1))
    ax2.set_yticklabels(np.arange(1, len(zscore_all)+1, 1))
    # hide top and right border
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    ax3.set_ylabel('z-Score',fontsize=AXIS_LABEL_FONTSIZE)
    ax3.set_xlabel('Seconds',fontsize=AXIS_LABEL_FONTSIZE)
    ax3.legend(loc=MY_LEGEND_LOC, bbox_to_anchor=MY_LEGENG_POS)
    ax3.set_xlim([-perievent_options_dict['sec_before'], perievent_options_dict['sec_after']])

    # bounds returns (x0,y0,width,height)
    # print("heatmap:")
    x_heatmap,y_heatmap,width_heatmap,height_heatmap = ax2.get_position().bounds
    # print(ax.get_position().bounds)
    # Average
    x_average,y_average,width_average,height_average = ax.get_position().bounds
    # print(ax.get_position().bounds)
    # readjust the position to match the sice of heatmap
    ax.set_position([x_average,y_average,width_heatmap,height_average])
    x_zscore,y_zscore,width_zscore,height_zscore = ax3.get_position().bounds
    # readjust the position of zscore to match the sice of heatmap
    ax3.set_position([x_zscore,y_zscore,width_heatmap,height_zscore])
    
    # # set the spacing between subplots 
    # canvas.fig.subplots_adjust(left=MY_LEFT, 
    #                 bottom=MY_BOTTOM,  
    #                 right=MY_RIGHT,  
    #                 top=MY_TOP,  
    #                 wspace=MY_WSPACE,  
    #                 hspace=MY_HSPACE) 
    canvas.draw()
    
def plot_perievent_avg_zscore_trials(canvas,subject,current_trials, perievent_options_dict,analyzed_perievent_dict):
    event_name = perievent_options_dict["event_name"] if len(perievent_options_dict["event_name"]) > 0 else perievent_options_dict["event"]
    
    # get data to plot avg
    ts_signal = analyzed_perievent_dict["average"]["ts_signal"]
    ts_control = analyzed_perievent_dict["average"]["ts_control"]
    signal = analyzed_perievent_dict["average"]["signal"]
    control = analyzed_perievent_dict["average"]["control"]
    std_signal = analyzed_perievent_dict["average"]["std_signal"]
    std_control = analyzed_perievent_dict["average"]["std_control"]
    # zscore
    ts = analyzed_perievent_dict["zscore"]["ts"]
    zscore_all = analyzed_perievent_dict["zscore"]["zscored"]
    zerror = analyzed_perievent_dict["zscore"]["zerror"]
    
    # plot
    # clear previous figure
    canvas.fig.clf()
    ax = canvas.fig.add_subplot(311)
    ax2 = canvas.fig.add_subplot(312,sharex=ax)
    ax3 = canvas.fig.add_subplot(313,sharex=ax)
    
    # average
    ax.plot(ts_signal, signal,
             color=SIGNAL_COLOR_RGB,
             linewidth=1,
             label = "signal")
    ax.plot(ts_control, control,
             color=CONTROL_COLOR_RGB,
             linewidth=1,
             label = "control")
    # plot a vertical line at t=0
    ax.axvline(x=0, linewidth=2, color='k', alpha=0.3,label=event_name)
    # plot standard error bands
    ax.fill_between(ts_signal, signal+std_signal, signal-std_signal,
                      facecolor=SIGNAL_COLOR_RGB, alpha=0.1)
    ax.fill_between(ts_control, control+std_control, control-std_control,
                      facecolor=CONTROL_COLOR_RGB, alpha=0.1)
    
    # zscore
    # Heat Map based on z score of control fit subtracted signal
    cs = ax2.imshow(zscore_all, cmap=plt.cm.bwr, interpolation='none', aspect="auto",
                    extent=[-perievent_options_dict['sec_before'], perievent_options_dict['sec_after'], 
                            len(zscore_all),0])
    canvas.fig.colorbar(cs, ax=ax2,pad=0.01, fraction=0.02)
    # plot the z-score trace for the signal with std error bands
    ax3.plot(ts, np.mean(zscore_all, axis=0), linewidth=2, color=SIGNAL_COLOR_RGB, label='Mean Z-Score')
    # plot all trials
    for i in range(len(zscore_all)):
        if i == 0:
            ax3.plot(ts, zscore_all[i], linewidth=0.5, color=SIGNAL_COLOR_RGB, alpha=0.2, label='Trials')
        else:
            ax3.plot(ts, zscore_all[i], linewidth=0.5, color=SIGNAL_COLOR_RGB, alpha=0.2)
    ax3.axvline(x=0, linewidth=2, color='k', alpha=0.3,label=event_name)
    
    my_title = 'Subject: ' + subject
    canvas.fig.suptitle(my_title, fontsize=FIGURE_TITLE_FONT_SIZE)
    ax.set_title('Average From '+str(len(current_trials))+' Trials',fontsize=TITLE_FONT_SIZE)
    # hide top and right border
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel('Time (sec)', fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_ylabel('mV', fontsize=AXIS_LABEL_FONTSIZE)
    ax.legend(loc=MY_LEGEND_LOC, bbox_to_anchor=MY_LEGENG_POS)
    ax.set_xlim([-perievent_options_dict['sec_before'], perievent_options_dict['sec_after']])
    ax2.set_title('Individual z-Score Traces',fontsize=TITLE_FONT_SIZE)
    ax2.set_ylabel('Trials',fontsize=AXIS_LABEL_FONTSIZE)
    ax2.set_xlabel('Seconds from Event Onset',fontsize=AXIS_LABEL_FONTSIZE)
    ax2.set_yticks(np.arange(0.5,len(zscore_all), 1))
    ax2.set_yticklabels(np.arange(1, len(zscore_all)+1, 1))
    # hide top and right border
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    ax3.set_ylabel('z-Score',fontsize=AXIS_LABEL_FONTSIZE)
    ax3.set_xlabel('Seconds',fontsize=AXIS_LABEL_FONTSIZE)
    ax3.legend(loc=MY_LEGEND_LOC, bbox_to_anchor=MY_LEGENG_POS)
    ax3.set_xlim([-perievent_options_dict['sec_before'], perievent_options_dict['sec_after']])

    # bounds returns (x0,y0,width,height)
    # print("heatmap:")
    x_heatmap,y_heatmap,width_heatmap,height_heatmap = ax2.get_position().bounds
    # print(ax.get_position().bounds)
    # Average
    x_average,y_average,width_average,height_average = ax.get_position().bounds
    # print(ax.get_position().bounds)
    # readjust the position to match the sice of heatmap
    ax.set_position([x_average,y_average,width_heatmap,height_average])
    x_zscore,y_zscore,width_zscore,height_zscore = ax3.get_position().bounds
    # readjust the position of zscore to match the sice of heatmap
    ax3.set_position([x_zscore,y_zscore,width_heatmap,height_zscore])
    
    # # set the spacing between subplots 
    # canvas.fig.subplots_adjust(left=MY_LEFT, 
    #                 bottom=MY_BOTTOM,  
    #                 right=MY_RIGHT,  
    #                 top=MY_TOP,  
    #                 wspace=MY_WSPACE,  
    #                 hspace=MY_HSPACE) 
    canvas.draw()
    
def plot_perievent_auc_alone(canvas,subject,perievent_options_dict,analyzed_perievent_dict,export,save_plots,group_name,settings_dict,export_loc_data):   
    settings_dict[0]["subject"] = subject
    settings_dict[0]["subject_group_name"] = group_name
    settings_dict[0]["baseline_from_sec"] =perievent_options_dict["baseline_from"]
    settings_dict[0]["baseline_to_sec"] =perievent_options_dict["baseline_to"]
    settings_dict[0]["auc_pre_from"] =perievent_options_dict["auc_pre_from"]
    settings_dict[0]["auc_pre_to"] =perievent_options_dict["auc_pre_to"]
    settings_dict[0]["auc_post_from"] =perievent_options_dict["auc_post_from"]
    settings_dict[0]["auc_post_to"] =perievent_options_dict["auc_post_to"]
    dump_path,file_beginning = export_loc_data 
    event_name = perievent_options_dict["event_name"] if len(perievent_options_dict["event_name"]) > 0 else perievent_options_dict["event"]
    # get data to plot
    AUC = analyzed_perievent_dict["auc"]["auc_data"]
    p = analyzed_perievent_dict["auc"]["p"]
    AUC_err_pre = analyzed_perievent_dict["auc"]["auc_err_pre"]
    AUC_err_post = analyzed_perievent_dict["auc"]["auc_err_post"]
    
    # plot
    # clear previous figure
    canvas.fig.clf()
    ax = canvas.fig.add_subplot(111)
    # polt bars with error bars
    ax.bar(np.arange(len(AUC)), AUC,
           color=[.8, .8, .8], align='center', alpha=0.5)
    ax.axhline(y=0, linewidth=0.5, color='k')
    ax.errorbar(np.arange(len(AUC)), AUC, yerr=[AUC_err_pre,AUC_err_post], fmt='none',color='k',capsize=4)
    x1, x2 = 0, 1 # columns indices for labels
    # y, h, col = max(AUC) + 2, 2, 'k'
    y, h, col = max(AUC) + max([AUC_err_pre,AUC_err_post]), 2, 'k'
    # # plot significance bar if p < .05
    # if p < 0.05:
    #     ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    #     ax.text((x1+x2)*.5, y+h, "*", ha='center', va='bottom', color=col)
        
    my_title = 'Subject: ' + subject
    canvas.fig.suptitle(my_title, fontsize=14)
    # hide top and right border
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    if min(AUC) < 0:
        ax.set_ylim(min(AUC)-2*h, y+2*h)
    else:
        ax.set_ylim(0-2, y+2*h)
    ax.set_ylabel('AUC',fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_title('Pre vs Post',fontsize=TITLE_FONT_SIZE)
    ax.set_xticks(np.arange(-1, len(AUC)+1))
    ax.set_xticklabels(['', 'PRE', 'POST', ''])
    ax.set_yticks([round(min(AUC), 2),round(max(AUC), 2)])
    
    # create subfolder in current subject folder
    if export == True:
        file_name = file_beginning+"_"+ event_name+"_area_under_curve_perievent.csv"
        if file_beginning != subject:
            file_name = file_beginning+"_"+ subject+"_"+event_name+"_area_under_curve_perievent.csv"
        dump_file_path = os.path.join(dump_path,file_name)
        # raw_df= pd.DataFrame({"pre_AUC":[AUC[0]],
        #                          "post_AUC" : [AUC[1]]})
         # save each trial's auc instead
        raw_df= pd.DataFrame({"Trial":[i+1 for i in range(len(analyzed_perievent_dict["aucs_pre_by_trial"]))],
                                    "Pre_AUC" : analyzed_perievent_dict["aucs_pre_by_trial"],
                                    "Post_AUC":analyzed_perievent_dict["aucs_post_by_trial"]})
        settings_df = get_settings_df(settings_dict)
        settings_df.to_csv(dump_file_path,header=False)            
        with open(dump_file_path,'a') as f:
            f.write("\n")
        raw_df.to_csv(dump_file_path, index=False,mode='a') 
        # export also auc binned in one second bins
        sec_auc_dfs = []
        for i in range(len(analyzed_perievent_dict["auc_by_sec"])):
            d = pd.DataFrame({"Trial_"+str(i+1):analyzed_perievent_dict["auc_by_sec"][i]})
            sec_auc_dfs.append(d)
        sec_auc_df = pd.concat(sec_auc_dfs,axis = 1)
      
        auc_file_name = file_beginning+"_"+ subject+"_"+event_name+"_auc_one_second_bins.csv"
        if file_beginning == subject:
            auc_file_name = file_beginning+"_"+event_name+"_auc_one_second_bins.csv"
        dump_file_path = os.path.join(dump_path,auc_file_name)
        settings_df.to_csv(dump_file_path,header=False)            
        with open(dump_file_path,'a') as f:
            f.write("\n")
        sec_auc_df.to_csv(dump_file_path, index=False,mode='a') 
             
    if export == True and save_plots == True:
        plot_file_name = file_beginning+"_"+ event_name+"_area_under_curve_perievent.png"
        plot_file_name2 = file_beginning+"_"+ event_name+"_area_under_curve_perievent.svg"
        if file_beginning != subject:
            plot_file_name = file_beginning+"_"+ subject+"_"+ event_name+"_area_under_curve_perievent.png"
            plot_file_name2 = file_beginning+"_"+ subject+"_"+ event_name+"_area_under_curve_perievent.svg"
        dump_plot_file_path = os.path.join(dump_path,plot_file_name)
        canvas.fig.savefig(dump_plot_file_path)
        # also as svg
        dump_plot_file_path = os.path.join(dump_path,plot_file_name2)
        canvas.fig.savefig(dump_plot_file_path, format='svg', dpi=DPI4SVG)
    else:
        canvas.draw()
    
    
def plot_perievent_avg_auc(canvas,subject,current_trials,perievent_options_dict,analyzed_perievent_dict):       
    event_name = perievent_options_dict["event_name"] if len(perievent_options_dict["event_name"]) > 0 else perievent_options_dict["event"]
    # get avg data to plot
    ts_signal = analyzed_perievent_dict["average"]["ts_signal"]
    ts_control = analyzed_perievent_dict["average"]["ts_control"]
    signal = analyzed_perievent_dict["average"]["signal"]
    control = analyzed_perievent_dict["average"]["control"]
    std_signal = analyzed_perievent_dict["average"]["std_signal"]
    std_control = analyzed_perievent_dict["average"]["std_control"]
    # AUC data
    AUC = analyzed_perievent_dict["auc"]["auc_data"]
    p = analyzed_perievent_dict["auc"]["p"]
    AUC_err_pre = analyzed_perievent_dict["auc"]["auc_err_pre"]
    AUC_err_post = analyzed_perievent_dict["auc"]["auc_err_post"]
    
    # plot
    # clear previous figure
    canvas.fig.clf()
    ax = canvas.fig.add_subplot(211)
    ax2 = canvas.fig.add_subplot(212)
    
    ax.plot(ts_signal, signal,
             color=SIGNAL_COLOR_RGB,
             linewidth=1,
             label = "signal")
    ax.plot(ts_control, control,
             color=CONTROL_COLOR_RGB,
             linewidth=1,
             label = "control")
    # plot a vertical line at t=0
    ax.axvline(x=0, linewidth=2, color='k', alpha=0.3,label=event_name)
    # plot standard error bands
    ax.fill_between(ts_signal, signal+std_signal, signal-std_signal,
                      facecolor=SIGNAL_COLOR_RGB, alpha=0.1)
    ax.fill_between(ts_control, control+std_control, control-std_control,
                      facecolor=CONTROL_COLOR_RGB, alpha=0.1)
    # polt bars with error bars
    ax2.bar(np.arange(len(AUC)), AUC, color=[.8, .8, .8], align='center', alpha=0.5)
    ax2.axhline(y=0, linewidth=0.5, color='k')
    ax2.errorbar(np.arange(len(AUC)), AUC, yerr=[AUC_err_pre,AUC_err_post], fmt='none',color='k',capsize=4)
    x1, x2 = 0, 1 # columns indices for labels
    # y, h, col = max(AUC) + 2, 2, 'k'
    y, h, col = max(AUC) + max([AUC_err_pre,AUC_err_post]), 2, 'k'
    # # plot significance bar if p < .05
    # if p < 0.05:
    #     ax2.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    #     ax2.text((x1+x2)*.5, y+h, "*", ha='center', va='bottom', color=col)
        
    my_title = 'Subject: ' + subject
    canvas.fig.suptitle(my_title, fontsize=FIGURE_TITLE_FONT_SIZE)
    ax.set_title('Average From '+str(len(current_trials))+' Trials',fontsize=TITLE_FONT_SIZE)
    # hide top and right border
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel('Time (sec)', fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_ylabel('mV', fontsize=AXIS_LABEL_FONTSIZE)
    ax.legend(loc=MY_LEGEND_LOC, bbox_to_anchor=MY_LEGENG_POS)
    
    # hide top and right border
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
    if min(AUC) < 0:
        ax2.set_ylim(min(AUC)-2*h, y+2*h)
    else:
        ax2.set_ylim(0-2, y+2*h)
    ax2.set_ylabel('AUC',fontsize=AXIS_LABEL_FONTSIZE)
    ax2.set_title('Pre vs Post',fontsize=TITLE_FONT_SIZE)
    ax2.set_xticks(np.arange(-1, len(AUC)+1))
    ax2.set_xticklabels(['', 'PRE', 'POST', ''])
    ax2.set_yticks([round(min(AUC), 2),round(max(AUC), 2)])
    
    canvas.draw()
    
def plot_perievent_zscore_auc(canvas,subject, perievent_options_dict,analyzed_perievent_dict):
    event_name = perievent_options_dict["event_name"] if len(perievent_options_dict["event_name"]) > 0 else perievent_options_dict["event"]
    # get data to plot
    ts = analyzed_perievent_dict["zscore"]["ts"]
    zscore_all = analyzed_perievent_dict["zscore"]["zscored"]
    zerror = analyzed_perievent_dict["zscore"]["zerror"]
    # AUC data
    AUC = analyzed_perievent_dict["auc"]["auc_data"]
    p = analyzed_perievent_dict["auc"]["p"]
    AUC_err_pre = analyzed_perievent_dict["auc"]["auc_err_pre"]
    AUC_err_post = analyzed_perievent_dict["auc"]["auc_err_post"]
    
    # plot
    # clear previous figure
    canvas.fig.clf()
    ax = canvas.fig.add_subplot(311)
    ax2 = canvas.fig.add_subplot(312,sharex=ax)
    ax3 = canvas.fig.add_subplot(313)
    
    # Heat Map based on z score of control fit subtracted signal
    cs = ax.imshow(zscore_all, cmap=plt.cm.bwr, interpolation='none', aspect="auto",
                    extent=[-perievent_options_dict['sec_before'], perievent_options_dict['sec_after'], 
                            len(zscore_all),0])
    canvas.fig.colorbar(cs, ax=ax,pad=0.01, fraction=0.02)
    # plot the z-score trace for the signal with std error bands
    ax2.plot(ts, np.mean(zscore_all, axis=0), linewidth=2, color=SIGNAL_COLOR_RGB,label='Mean Z-Score')
    ax2.fill_between(ts, np.mean(zscore_all, axis=0)+zerror
                          ,np.mean(zscore_all, axis=0)-zerror, facecolor=SIGNAL_COLOR_RGB, alpha=0.2,label="Standard\nerror")
    ax2.axvline(x=0, linewidth=2, color='k', alpha=0.3, label=event_name)
    
    # AUC
    # polt bars with error bars
    ax3.bar(np.arange(len(AUC)), AUC, color=[.8, .8, .8], align='center', alpha=0.5)
    ax3.axhline(y=0,linewidth=0.5, color='k')
    ax3.errorbar(np.arange(len(AUC)), AUC, yerr=[AUC_err_pre,AUC_err_post], fmt='none',color='k',capsize=4)
    x1, x2 = 0, 1 # columns indices for labels
    # y, h, col = max(AUC) + 2, 2, 'k'
    y, h, col = max(AUC) + max([AUC_err_pre,AUC_err_post]), 2, 'k'
    # # plot significance bar if p < .05
    # if p < 0.05:
    #     ax3.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    #     ax3.text((x1+x2)*.5, y+h, "*", ha='center', va='bottom', color=col)
  
    my_title = 'Subject: ' + subject
    canvas.fig.suptitle(my_title, fontsize=FIGURE_TITLE_FONT_SIZE)
    # hide top and right border
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax.set_title('Individual z-Score Traces',fontsize=TITLE_FONT_SIZE)
    ax.set_ylabel('Trials',fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_xlabel('Seconds from Event Onset',fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_yticks(np.arange(0.5,len(zscore_all), 1))
    ax.set_yticklabels(np.arange(1, len(zscore_all)+1, 1))
    ax2.set_ylabel('z-Score',fontsize=AXIS_LABEL_FONTSIZE)
    ax2.set_xlabel('Seconds',fontsize=AXIS_LABEL_FONTSIZE)
    ax2.legend(loc=MY_LEGEND_LOC, bbox_to_anchor=MY_LEGENG_POS)
    ax2.set_xlim([-perievent_options_dict['sec_before'], perievent_options_dict['sec_after']])
    # hide top and right border
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    ax3.spines['bottom'].set_visible(False)
    if min(AUC) < 0:
        ax3.set_ylim(min(AUC)-2*h, y+2*h)
    else:
        ax3.set_ylim(0-2, y+2*h)
    ax3.set_ylabel('AUC',fontsize=AXIS_LABEL_FONTSIZE)
    ax3.set_title('Pre vs Post',fontsize=TITLE_FONT_SIZE)
    ax3.set_xticks(np.arange(-1, len(AUC)+1))
    ax3.set_xticklabels(['', 'PRE', 'POST', ''])
    ax3.set_yticks([round(min(AUC), 2),round(max(AUC), 2)])

    # bounds returns (x0,y0,width,height)
    # print("heatmap:")
    x_heatmap,y_heatmap,width_heatmap,height_heatmap = ax.get_position().bounds
    x_zscore,y_zscore,width_zscore,height_zscore = ax2.get_position().bounds
    # readjust the position of zscore to match the sice of heatmap
    ax2.set_position([x_zscore,y_zscore,width_heatmap,height_zscore])
    x_auc,y_auc,width_auc,height_auc = ax3.get_position().bounds
    # readjust the position of zscore to match the sice of heatmap
    ax3.set_position([x_auc,y_auc,width_heatmap,height_auc])
    
    # # set the spacing between subplots 
    # canvas.fig.subplots_adjust(left=MY_LEFT, 
    #                 bottom=MY_BOTTOM,  
    #                 right=MY_RIGHT,  
    #                 top=MY_TOP,  
    #                 wspace=MY_WSPACE,  
    #                 hspace=MY_HSPACE) 
    canvas.draw()
    
def plot_perievent_zscore_trials_auc(canvas,subject, perievent_options_dict,analyzed_perievent_dict):
    event_name = perievent_options_dict["event_name"] if len(perievent_options_dict["event_name"]) > 0 else perievent_options_dict["event"]
    # get data to plot
    ts = analyzed_perievent_dict["zscore"]["ts"]
    zscore_all = analyzed_perievent_dict["zscore"]["zscored"]
    zerror = analyzed_perievent_dict["zscore"]["zerror"]
    # AUC data
    AUC = analyzed_perievent_dict["auc"]["auc_data"]
    p = analyzed_perievent_dict["auc"]["p"]
    AUC_err_pre = analyzed_perievent_dict["auc"]["auc_err_pre"]
    AUC_err_post = analyzed_perievent_dict["auc"]["auc_err_post"]
    
    # plot
    # clear previous figure
    canvas.fig.clf()
    ax = canvas.fig.add_subplot(311)
    ax2 = canvas.fig.add_subplot(312,sharex=ax)
    ax3 = canvas.fig.add_subplot(313)
    
    # Heat Map based on z score of control fit subtracted signal
    cs = ax.imshow(zscore_all, cmap=plt.cm.bwr, interpolation='none', aspect="auto",
                    extent=[-perievent_options_dict['sec_before'], perievent_options_dict['sec_after'], 
                            len(zscore_all),0])
    canvas.fig.colorbar(cs, ax=ax,pad=0.01, fraction=0.02)
    # plot the z-score trace for the signal with std error bands
    ax2.plot(ts, np.mean(zscore_all, axis=0), linewidth=2, color=SIGNAL_COLOR_RGB,label='Mean Z-Score')
    # plot all trials
    for i in range(len(zscore_all)):
        if i == 0:
            ax2.plot(ts, zscore_all[i], linewidth=0.5, color=SIGNAL_COLOR_RGB, alpha=0.2, label='Trials')
        else:
            ax2.plot(ts, zscore_all[i], linewidth=0.5, color=SIGNAL_COLOR_RGB, alpha=0.2)
    ax2.axvline(x=0, linewidth=2, color='k', alpha=0.3, label=event_name)
    
    # AUC
    ax3.bar(np.arange(len(AUC)), AUC, color=[.8, .8, .8], align='center', alpha=0.5)
    ax3.axhline(y=0,linewidth=0.5,color='k')
    ax3.errorbar(np.arange(len(AUC)), AUC, yerr=[AUC_err_pre,AUC_err_post], fmt='none',color='k',capsize=4)
    x1, x2 = 0, 1 # columns indices for labels
    # y, h, col = max(AUC) + 2, 2, 'k'
    y, h, col = max(AUC) + max([AUC_err_pre,AUC_err_post]), 2, 'k'
    # # plot significance bar if p < .05
    # if p < 0.05:
    #     ax3.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    #     ax3.text((x1+x2)*.5, y+h, "*", ha='center', va='bottom', color=col)
  
    my_title = 'Subject: ' + subject
    canvas.fig.suptitle(my_title, fontsize=FIGURE_TITLE_FONT_SIZE)
    ax.set_title('Individual z-Score Traces',fontsize=TITLE_FONT_SIZE)
    ax.set_ylabel('Trials',fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_xlabel('Seconds from Event Onset',fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_yticks(np.arange(0.5,len(zscore_all), 1))
    ax.set_yticklabels(np.arange(1, len(zscore_all)+1, 1))
    # hide top and right border
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.set_ylabel('z-Score',fontsize=AXIS_LABEL_FONTSIZE)
    ax2.set_xlabel('Seconds',fontsize=AXIS_LABEL_FONTSIZE)
    ax2.legend(loc=MY_LEGEND_LOC, bbox_to_anchor=MY_LEGENG_POS)
    ax2.set_xlim([-perievent_options_dict['sec_before'], perievent_options_dict['sec_after']])
    # hide top and right border
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    ax3.spines['bottom'].set_visible(False)
    if min(AUC) < 0:
        ax3.set_ylim(min(AUC)-2*h, y+2*h)
    else:
        ax3.set_ylim(0-2, y+2*h)
    ax3.set_ylabel('AUC',fontsize=AXIS_LABEL_FONTSIZE)
    ax3.set_title('Pre vs Post',fontsize=TITLE_FONT_SIZE)
    ax3.set_xticks(np.arange(-1, len(AUC)+1))
    ax3.set_xticklabels(['', 'PRE', 'POST', ''])
    ax3.set_yticks([round(min(AUC), 2),round(max(AUC), 2)])

    # bounds returns (x0,y0,width,height)
    # print("heatmap:")
    x_heatmap,y_heatmap,width_heatmap,height_heatmap = ax.get_position().bounds
    x_zscore,y_zscore,width_zscore,height_zscore = ax2.get_position().bounds
    # readjust the position of zscore to match the sice of heatmap
    ax2.set_position([x_zscore,y_zscore,width_heatmap,height_zscore])
    x_auc,y_auc,width_auc,height_auc = ax3.get_position().bounds
    # readjust the position of zscore to match the sice of heatmap
    ax3.set_position([x_auc,y_auc,width_heatmap,height_auc])
    
    # # set the spacing between subplots 
    # canvas.fig.subplots_adjust(left=MY_LEFT, 
    #                 bottom=MY_BOTTOM,  
    #                 right=MY_RIGHT,  
    #                 top=MY_TOP,  
    #                 wspace=MY_WSPACE,  
    #                 hspace=MY_HSPACE) 
    canvas.draw()
    
def plot_all_perievent(canvas,subject,current_trials, perievent_options_dict,analyzed_perievent_dict):
    event_name = perievent_options_dict["event_name"] if len(perievent_options_dict["event_name"]) > 0 else perievent_options_dict["event"]
    # get avg data to plot
    ts_signal = analyzed_perievent_dict["average"]["ts_signal"]
    ts_control = analyzed_perievent_dict["average"]["ts_control"]
    signal = analyzed_perievent_dict["average"]["signal"]
    control = analyzed_perievent_dict["average"]["control"]
    std_signal = analyzed_perievent_dict["average"]["std_signal"]
    std_control = analyzed_perievent_dict["average"]["std_control"]
    # get zscore data to plot
    ts = analyzed_perievent_dict["zscore"]["ts"]
    zscore_all = analyzed_perievent_dict["zscore"]["zscored"]
    zerror = analyzed_perievent_dict["zscore"]["zerror"]
    # AUC data
    AUC = analyzed_perievent_dict["auc"]["auc_data"]
    p = analyzed_perievent_dict["auc"]["p"]
    AUC_err_pre = analyzed_perievent_dict["auc"]["auc_err_pre"]
    AUC_err_post = analyzed_perievent_dict["auc"]["auc_err_post"]
    
    # plot
    # clear previous figure
    canvas.fig.clf()
    ax = canvas.fig.add_subplot(411)
    ax2 = canvas.fig.add_subplot(412,sharex=ax)
    ax3 = canvas.fig.add_subplot(413,sharex=ax)
    ax4 = canvas.fig.add_subplot(414)
    
    ax.plot(ts_signal, signal,
             color=SIGNAL_COLOR_RGB,
             linewidth=1,
             label = "signal")
    ax.plot(ts_control, control,
             color=CONTROL_COLOR_RGB,
             linewidth=1,
             label = "control")
    # plot a vertical line at t=0
    ax.axvline(x=0, linewidth=2, color='k', alpha=0.3,label=event_name)
    # plot standard error bands
    ax.fill_between(ts_signal, signal+std_signal, signal-std_signal,
                      facecolor=SIGNAL_COLOR_RGB, alpha=0.1)
    ax.fill_between(ts_control, control+std_control, control-std_control,
                      facecolor=CONTROL_COLOR_RGB, alpha=0.1)
    
    # Heat Map based on z score of control fit subtracted signal
    cs = ax2.imshow(zscore_all, cmap=plt.cm.bwr, interpolation='none', aspect="auto",
                    extent=[-perievent_options_dict['sec_before'], perievent_options_dict['sec_after'], 
                            len(zscore_all),0])
    canvas.fig.colorbar(cs, ax=ax2,pad=0.01, fraction=0.02)
    # plot the z-score trace for the signal with std error bands
    ax3.plot(ts, np.mean(zscore_all, axis=0), linewidth=2, color=SIGNAL_COLOR_RGB, label='Mean Z-Score')
    ax3.fill_between(ts, np.mean(zscore_all, axis=0)+zerror
                          ,np.mean(zscore_all, axis=0)-zerror, facecolor=SIGNAL_COLOR_RGB, alpha=0.2,label="Standard\nerror")
    ax3.axvline(x=0, linewidth=2, color='k', alpha=0.3, label=event_name)
    
    # AUC
    ax4.bar(np.arange(len(AUC)), AUC, color=[.8, .8, .8], align='center', alpha=0.5)
    ax4.axhline(y=0,linewidth=0.5,color='k')
    ax4.errorbar(np.arange(len(AUC)), AUC, yerr=[AUC_err_pre,AUC_err_post], fmt='none',color='k',capsize=4)
    x1, x2 = 0, 1 # columns indices for labels
    # y, h, col = max(AUC) + 2, 2, 'k'
    y, h, col = max(AUC) + max([AUC_err_pre,AUC_err_post]), 2, 'k'
    # # plot significance bar if p < .05
    # if p < 0.05:
    #     ax4.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    #     ax4.text((x1+x2)*.5, y+h, "*", ha='center', va='bottom', color=col)
        
    
    my_title = 'Subject: ' + subject
    canvas.fig.suptitle(my_title, fontsize=FIGURE_TITLE_FONT_SIZE)
    ax.set_title('Average From '+str(len(current_trials))+' Trials',fontsize=TITLE_FONT_SIZE)
    # hide top and right border
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
#    ax.set_xlabel('Time (sec)', fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_ylabel('mV', fontsize=AXIS_LABEL_FONTSIZE)
    ax.legend(loc=MY_LEGEND_LOC, bbox_to_anchor=MY_LEGENG_POS)
    ax.set_xlim([-perievent_options_dict['sec_before'], perievent_options_dict['sec_after']])
    ax2.set_title('Individual z-Score Traces',fontsize=TITLE_FONT_SIZE)
    ax2.set_ylabel('Trials',fontsize=AXIS_LABEL_FONTSIZE)
    ax2.set_yticks(np.arange(0.5,len(zscore_all), 1))
    ax2.set_yticklabels(np.arange(1, len(zscore_all)+1, 1))
#    ax2.set_xlabel('Seconds from Event Onset')
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    ax3.set_ylabel('z-Score',fontsize=AXIS_LABEL_FONTSIZE)
    ax3.legend(loc=MY_LEGEND_LOC, bbox_to_anchor=MY_LEGENG_POS)
    ax3.set_xlim([-perievent_options_dict['sec_before'], perievent_options_dict['sec_after']])
#    ax3.set_xlabel('Seconds')
    # hide top and right border
    ax4.spines['top'].set_visible(False)
    ax4.spines['right'].set_visible(False)
    ax4.spines['bottom'].set_visible(False)
    if min(AUC) < 0:
        ax4.set_ylim(min(AUC)-2*h, y+2*h)
    else:
        ax4.set_ylim(0-2, y+2*h)
    ax4.set_ylabel('AUC',fontsize=AXIS_LABEL_FONTSIZE)
    ax4.set_title('Pre vs Post',fontsize=TITLE_FONT_SIZE)
    ax4.set_xticks(np.arange(-1, len(AUC)+1))
    ax4.set_xticklabels(['', 'PRE', 'POST', ''])
    ax4.set_yticks([round(min(AUC), 2),round(max(AUC), 2)])

    # bounds returns (x0,y0,width,height)
    # print("heatmap:")
    x_heatmap,y_heatmap,width_heatmap,height_heatmap = ax2.get_position().bounds
    # Average
    x_average,y_average,width_average,height_average = ax.get_position().bounds
    # readjust the position to match the sice of heatmap
    ax.set_position([x_average,y_average,width_heatmap,height_average])
    x_zscore,y_zscore,width_zscore,height_zscore = ax3.get_position().bounds
    # readjust the position of zscore to match the sice of heatmap
    ax3.set_position([x_zscore,y_zscore,width_heatmap,height_zscore])
    x_auc,y_auc,width_auc,height_auc = ax4.get_position().bounds
    # readjust the position of zscore to match the sice of heatmap
    ax4.set_position([x_auc,y_auc,width_heatmap,height_auc])

    # # set the spacing between subplots 
    # canvas.fig.subplots_adjust(left=MY_LEFT, 
    #                 bottom=MY_BOTTOM,  
    #                 right=MY_RIGHT,  
    #                 top=MY_TOP,  
    #                 wspace=MY_WSPACE,  
    #                 hspace=MY_HSPACE) 
    canvas.draw()

def plot_all_perievent_zscore_trials(canvas,subject,current_trials, perievent_options_dict,analyzed_perievent_dict):
    event_name = perievent_options_dict["event_name"] if len(perievent_options_dict["event_name"]) > 0 else perievent_options_dict["event"]
    # get avg data to plot
    ts_signal = analyzed_perievent_dict["average"]["ts_signal"]
    ts_control = analyzed_perievent_dict["average"]["ts_control"]
    signal = analyzed_perievent_dict["average"]["signal"]
    control = analyzed_perievent_dict["average"]["control"]
    std_signal = analyzed_perievent_dict["average"]["std_signal"]
    std_control = analyzed_perievent_dict["average"]["std_control"]
    # get zscore data to plot
    ts = analyzed_perievent_dict["zscore"]["ts"]
    zscore_all = analyzed_perievent_dict["zscore"]["zscored"]
    zerror = analyzed_perievent_dict["zscore"]["zerror"]
    # AUC data
    AUC = analyzed_perievent_dict["auc"]["auc_data"]
    p = analyzed_perievent_dict["auc"]["p"]
    AUC_err_pre = analyzed_perievent_dict["auc"]["auc_err_pre"]
    AUC_err_post = analyzed_perievent_dict["auc"]["auc_err_post"]
    # plot
    # clear previous figure
    canvas.fig.clf()
    ax = canvas.fig.add_subplot(411)
    ax2 = canvas.fig.add_subplot(412,sharex=ax)
    ax3 = canvas.fig.add_subplot(413,sharex=ax)
    ax4 = canvas.fig.add_subplot(414)
    
    ax.plot(ts_signal, signal,
             color=SIGNAL_COLOR_RGB,
             linewidth=1,
             label = "signal")
    ax.plot(ts_control, control,
             color=CONTROL_COLOR_RGB,
             linewidth=1,
             label = "control")
    # plot a vertical line at t=0
    ax.axvline(x=0, linewidth=2, color='k', alpha=0.3,label=event_name)
    # plot standard error bands
    ax.fill_between(ts_signal, signal+std_signal, signal-std_signal,
                      facecolor=SIGNAL_COLOR_RGB, alpha=0.1)
    ax.fill_between(ts_control, control+std_control, control-std_control,
                      facecolor=CONTROL_COLOR_RGB, alpha=0.1)
    
    # Heat Map based on z score of control fit subtracted signal
    cs = ax2.imshow(zscore_all, cmap=plt.cm.bwr, interpolation='none', aspect="auto",
                    extent=[-perievent_options_dict['sec_before'], perievent_options_dict['sec_after'], 
                            len(zscore_all),0])
    canvas.fig.colorbar(cs, ax=ax2,pad=0.01, fraction=0.02)
    # plot the z-score trace for the signal with std error bands
    ax3.plot(ts, np.mean(zscore_all, axis=0), linewidth=2, color=SIGNAL_COLOR_RGB, label='Mean Z-Score')
    # plot all trials
    for i in range(len(zscore_all)):
        if i == 0:
            ax3.plot(ts, zscore_all[i], linewidth=0.5, color=SIGNAL_COLOR_RGB, alpha=0.2, label='Trials')
        else:
            ax3.plot(ts, zscore_all[i], linewidth=0.5, color=SIGNAL_COLOR_RGB, alpha=0.2)
    ax3.axvline(x=0, linewidth=2, color='k', alpha=0.3, label=event_name)
    
    # AUC
    ax4.bar(np.arange(len(AUC)), AUC, color=[.8, .8, .8], align='center', alpha=0.5)
    ax4.axhline(y=0,linewidth=0.5,color='k')
    ax4.errorbar(np.arange(len(AUC)), AUC, yerr=[AUC_err_pre,AUC_err_post], fmt='none',color='k',capsize=4)
    x1, x2 = 0, 1 # columns indices for labels
    # y, h, col = max(AUC) + 2, 2, 'k'
    y, h, col = max(AUC) + max([AUC_err_pre,AUC_err_post]), 2, 'k'
    # # plot significance bar if p < .05
    # if p < 0.05:
    #     ax4.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    #     ax4.text((x1+x2)*.5, y+h, "*", ha='center', va='bottom', color=col)
        
    
    my_title = 'Subject: ' + subject
    canvas.fig.suptitle(my_title, fontsize=FIGURE_TITLE_FONT_SIZE)
    ax.set_title('Average From '+str(len(current_trials))+' Trials',fontsize=TITLE_FONT_SIZE)
    # hide top and right border
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
#    ax.set_xlabel('Time (sec)', fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_ylabel('mV', fontsize=AXIS_LABEL_FONTSIZE)
    ax.legend(loc=MY_LEGEND_LOC, bbox_to_anchor=MY_LEGENG_POS)
    ax.set_xlim([-perievent_options_dict['sec_before'], perievent_options_dict['sec_after']])
    ax2.set_title('Individual z-Score Traces',fontsize=TITLE_FONT_SIZE)
    ax2.set_ylabel('Trials',fontsize=AXIS_LABEL_FONTSIZE)
    ax2.set_yticks(np.arange(0.5,len(zscore_all), 1))
    ax2.set_yticklabels(np.arange(1, len(zscore_all)+1, 1))
#    ax2.set_xlabel('Seconds from Event Onset')
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    ax3.set_ylabel('z-Score',fontsize=AXIS_LABEL_FONTSIZE)
    ax3.legend(loc=MY_LEGEND_LOC, bbox_to_anchor=MY_LEGENG_POS)
    ax3.set_xlim([-perievent_options_dict['sec_before'], perievent_options_dict['sec_after']])
#    ax3.set_xlabel('Seconds')
    # hide top and right border
    ax4.spines['top'].set_visible(False)
    ax4.spines['right'].set_visible(False)
    ax4.spines['bottom'].set_visible(False)
    if min(AUC) < 0:
        ax4.set_ylim(min(AUC)-2*h, y+2*h)
    else:
        ax4.set_ylim(0-2, y+2*h)
    ax4.set_ylabel('AUC',fontsize=AXIS_LABEL_FONTSIZE)
    ax4.set_title('Pre vs Post',fontsize=TITLE_FONT_SIZE)
    ax4.set_xticks(np.arange(-1, len(AUC)+1))
    ax4.set_xticklabels(['', 'PRE', 'POST', ''])
    ax4.set_yticks([round(min(AUC), 2),round(max(AUC), 2)])

    # bounds returns (x0,y0,width,height)
    # print("heatmap:")
    x_heatmap,y_heatmap,width_heatmap,height_heatmap = ax2.get_position().bounds
    # Average
    x_average,y_average,width_average,height_average = ax.get_position().bounds
    # readjust the position to match the sice of heatmap
    ax.set_position([x_average,y_average,width_heatmap,height_average])
    x_zscore,y_zscore,width_zscore,height_zscore = ax3.get_position().bounds
    # readjust the position of zscore to match the sice of heatmap
    ax3.set_position([x_zscore,y_zscore,width_heatmap,height_zscore])
    x_auc,y_auc,width_auc,height_auc = ax4.get_position().bounds
    # readjust the position of zscore to match the sice of heatmap
    ax4.set_position([x_auc,y_auc,width_heatmap,height_auc])

    # # set the spacing between subplots 
    # canvas.fig.subplots_adjust(left=MY_LEFT, 
    #                 bottom=MY_BOTTOM,  
    #                 right=MY_RIGHT,  
    #                 top=MY_TOP,  
    #                 wspace=MY_WSPACE,  
    #                 hspace=MY_HSPACE) 
    canvas.draw()
    
def plot_peaks(canvas,subject,data,options,blocks_windows_values,export,save_plots,group_name,settings_dict,export_loc_data):
    settings_dict[0]["subject"] = subject
    settings_dict[0]["subject_group_name"] = group_name
    show_norm_as = settings_dict[0]["show_norm_as"]
    dump_path,file_beginning = export_loc_data 
    # a list of three tuples. Each tuple is seconds from, seconds till
    time_windows_values = blocks_windows_values
    block1_from = time_windows_values[0][0]
    block1_till = time_windows_values[0][1]
    block2_from = time_windows_values[1][0]
    block2_till = time_windows_values[1][1]
    block3_from = time_windows_values[2][0]
    block3_till = time_windows_values[2][1]
    # check if there were any blocks
    blocks = False
    block1_spike_ts = []
    block1_spike_val = []
    block2_spike_ts = []
    block2_spike_val = []
    block3_spike_ts = []
    block3_spike_val = []
    # get peak params from gui
     # get peak params from gui
    height = None
    threshold = None
    distance = None
    prominence = None
    width = None
    wlen = None
    relHeight = 0.5
    plateau = None
    # if export was false, check also last parameter from options
    if export == False:
        export = options[-1]
        if export == True: # assume that user wants to plot as well
            save_plots = True
    
    if len(options[0]) > 0 and options[0] !="None":
        try:
            height = float(options[0]) 
        except:
            print("height problem")
    if len(options[1]) > 0 and options[1] !="None":
        try:
            threshold = float(options[1]) 
        except:
            print("threshold problem")
    if len(options[2]) > 0 and options[2] !="None":
        try:
            distance = float(options[2]) 
        except:
            print("distance problem")
    if len(options[3]) > 0 and options[3] !="None":
        try:
            prominence = float(options[3]) 
        except:
            print("prominence problem")
    if len(options[4]) > 0 and options[4] !="None":
        try:
            width = float(options[4]) 
            print("width",width)
        except:
            print("width problem")
    if len(options[5]) > 0 and options[5] !="None":
        try:
            wlen = float(options[5]) 
        except:
            print("wlen problem")
    if len(options[6]) > 0 and options[6] !="None":
        try:
            relHeight = float(options[6]) 
        except:
            print("relHeight problem")
    if len(options[7]) > 0 and options[7] !="None":
        try:
            plateau = float(options[7]) 
        except:
            print("plateau problem")
#    print(height,threshold,distance,prominence,width,wlen,relHeight,plateau)
    
    ts = data["ts"]
    total_seconds = ts[-1]-ts[0]
    normalized_signal = data["normalized_signal"]
#    print("total seconds",total_seconds)
#    print("total samples",len(normalized_signal))
#    print("1 sample",total_seconds/len(normalized_signal))
    # start time from zero
    ts_reset = [i*total_seconds/len(ts) for i in range(len(ts))]
    # distance from seconds to samples
    if distance != None and distance != "None":
        # find how many samples are per second now
        # since time is in second find when one starts
        samples_in_second = 0
        for i in range(len(ts_reset)):
            if ts_reset[i] > 1:
                samples_in_second = len(ts_reset[:i])
                break
        distance = distance*samples_in_second
    #safety check if user entered 0 they probably mean no distance
    if distance == 0:
        distance = None
#    print(ts_reset[-1])
    peaks, _ = find_peaks(normalized_signal, height= height, threshold = threshold, distance = distance,
                          prominence = prominence, width = width, wlen = wlen, rel_height = relHeight, plateau_size = plateau)
#    print("peaks\n",peaks)
    # adjust peaks results to the existing timestamps
    translated_peaks = [el*total_seconds/len(normalized_signal) for el in peaks]
    total_peaks = len(translated_peaks)
#    print(translated_peaks)
    
    # plot
    # clear previous figure
    canvas.fig.clf()
    ax = canvas.fig.add_subplot(3,1,(1, 2))
    ax2 = canvas.fig.add_subplot(3,1,3,sharex=ax)
    
    ax.plot(ts_reset,normalized_signal,
             color=SIGNAL_COLOR_RGB,
             linewidth=1,
             label = "signal")

#    ax.plot(peaks, normalized_signal[peaks], "xr")
    ax.plot(translated_peaks, normalized_signal[peaks], "xr")
    
    total_label = str(total_peaks)+ " Total Spikes"
    ax2.eventplot(translated_peaks, orientation='horizontal', linelengths=0.5, color='red',label = total_label)
    
    # get y limits of the top plot to place the block text under it
    bottom, top = ax.get_ylim()
    block_text_y_pos = bottom-0.6
    # bbox_props = dict(boxstyle='round', facecolor='k', alpha=0.1, hatch='X')
    bbox_props = dict(boxstyle='round', facecolor='k', alpha=0.1)
    # Plot windows if there are any
    if block1_from != "None" and block1_till != "None":
        print(f"From: {block1_from}, till: {block1_till}")
        blocks = True
        ax.axvline(x=block1_from,linestyle='dashed',color='k',alpha=0.5,linewidth=1)
        ax.axvline(x=block1_till,linestyle='dashed',color='k',alpha=0.5,linewidth=1)
        # ax.axvspan(block1_from,block1_till, alpha=0.1,color='k',hatch='X')
        ax.axvspan(block1_from,block1_till, alpha=0.1,color='k')
        # count spikes for that block
        for i in range(len(translated_peaks)):
            if translated_peaks[i] >= block1_from and translated_peaks[i] <= block1_till:
                block1_spike_ts.append(translated_peaks[i])
                block1_spike_val.append(normalized_signal[peaks][i])
        block1_text = str(len(block1_spike_ts))+" Spikes"
        canvas.fig.text(block1_from+1, block_text_y_pos, block1_text, ha='left', va='top', fontsize=AXIS_LABEL_FONTSIZE-2, bbox=bbox_props, transform=ax.transData)
    if block2_from != "None" and block2_till != "None":
        print(f"From: {block2_from}, till: {block2_till}")
        blocks = True
        ax.axvline(x=block2_from,linestyle='dashed',color='k',alpha=0.5,linewidth=1)
        ax.axvline(x=block2_till,linestyle='dashed',color='k',alpha=0.5,linewidth=1)
        # ax.axvspan(block2_from,block2_till, alpha=0.1,color='k',hatch='X')
        ax.axvspan(block2_from,block2_till, alpha=0.1,color='k')
        # count spikes for that block
        for i in range(len(translated_peaks)):
            if translated_peaks[i] >= block2_from and translated_peaks[i] <= block2_till:
                block2_spike_ts.append(translated_peaks[i])
                block2_spike_val.append(normalized_signal[peaks][i])
        block2_text = str(len(block2_spike_ts))+" Spikes"
        canvas.fig.text(block2_from+1, block_text_y_pos, block2_text, ha='left', va='top', fontsize=AXIS_LABEL_FONTSIZE-2, bbox=bbox_props, transform=ax.transData)
    if block3_from != "None" and block3_till != "None":
        print(f"From: {block3_from}, till: {block3_till}")
        blocks = True
        ax.axvline(x=block3_from,linestyle='dashed',color='k',alpha=0.5,linewidth=1)
        ax.axvline(x=block3_till,linestyle='dashed',color='k',alpha=0.5,linewidth=1)
        # ax.axvspan(block3_from,block3_till, alpha=0.1,color='k',hatch='X')
        ax.axvspan(block3_from,block3_till, alpha=0.1,color='k')
        # count spikes for that block
        for i in range(len(translated_peaks)):
            if translated_peaks[i] >= block3_from and translated_peaks[i] <= block3_till:
                block3_spike_ts.append(translated_peaks[i])
                block3_spike_val.append(normalized_signal[peaks][i])
        block3_text = str(len(block3_spike_ts))+" Spikes"
        canvas.fig.text(block3_from+1, block_text_y_pos, block3_text, ha='left', va='top', fontsize=AXIS_LABEL_FONTSIZE-2, bbox=bbox_props, transform=ax.transData)


    my_title = 'Subject: ' + subject
    canvas.fig.suptitle(my_title, fontsize=FIGURE_TITLE_FONT_SIZE)
    ax.set_ylabel('dF/F (%)', fontsize=AXIS_LABEL_FONTSIZE)
    if show_norm_as == "Z-Score":
        ax.set_ylabel('Z-Score', fontsize=AXIS_LABEL_FONTSIZE)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.xaxis.set_visible(False)
    ax.xaxis.set_ticks_position('none')
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.set_ylabel('Peaks', fontsize=AXIS_LABEL_FONTSIZE)
    ax2.set_xlabel('Time (sec)', fontsize=AXIS_LABEL_FONTSIZE)
    ax2.set_yticks([])
    
    ax2.legend(loc='lower left', bbox_to_anchor=(0.85, 0))
    # set the spacing between subplots 
    canvas.fig.subplots_adjust(left=MY_LEFT, 
                    bottom=MY_BOTTOM,  
                    right=0.95,  
                    top=MY_TOP,  
                    wspace=0,  
                    hspace=0) 
       
    # export
    if export == True:
        file_name = file_beginning + "_Spikes.csv"
        if file_beginning != subject:
            file_name = file_beginning +"_"+subject+ "_Spikes.csv"
        dump_file_path = os.path.join(dump_path,file_name)
        val = "dF/F (%)"
        if show_norm_as == "Z-Score":
            val = "Z-Score"
        if blocks == False:
            raw_df = pd.DataFrame({"(sec)":translated_peaks,
                                    val: normalized_signal[peaks]})
        else: # create dataframe with joined blocks that might have different number of rows each
            dfs = [pd.DataFrame({"(sec) Total Spikes":translated_peaks,
                                    val +" Total Spikes" : normalized_signal[peaks]})]
            if block1_from != "None" and block1_till != "None":
                # add nans as missing values to the blocks, so that they all have the same number of rows as total spikes
                while len(block1_spike_ts) < total_peaks:
                    block1_spike_ts.append(np.nan)
                    block1_spike_val.append(np.nan)
                block1_ts_col_name = "(sec) Spikes from "+str(block1_from)+" till "+str(block1_till)
                block1_val_col_name = val + "Spikes from "+str(block1_from)+" till "+str(block1_till)
                df =  pd.DataFrame({
                                    block1_ts_col_name:block1_spike_ts,
                                    block1_val_col_name:block1_spike_val})
                dfs.append(df)
            if block2_from != "None" and block2_till != "None":
                while len(block2_spike_ts) < total_peaks:
                    block2_spike_ts.append(np.nan)
                    block2_spike_val.append(np.nan)
                block2_ts_col_name = "(sec) Spikes from "+str(block2_from)+" till "+str(block2_till)
                block2_val_col_name = val + "Spikes from "+str(block2_from)+" till "+str(block2_till)
                df =  pd.DataFrame({
                                    block2_ts_col_name:block2_spike_ts,
                                    block2_val_col_name:block2_spike_val})
                dfs.append(df)
            if block3_from != "None" and block3_till != "None":
                while len(block3_spike_ts) < total_peaks:
                    block3_spike_ts.append(np.nan)
                    block3_spike_val.append(np.nan)
                block3_ts_col_name = "(sec) Spikes from "+str(block3_from)+" till "+str(block3_till)
                block3_val_col_name = val + "Spikes from "+str(block3_from)+" till "+str(block3_till)
                df =  pd.DataFrame({
                                    block3_ts_col_name:block3_spike_ts,
                                    block3_val_col_name:block3_spike_val})
                dfs.append(df)
            raw_df = pd.concat(dfs,axis=1)
            
        settings_df = get_settings_df(settings_dict)
        settings_df.to_csv(dump_file_path,header=False)            
        with open(dump_file_path,'a') as f:
            f.write("\n")
        raw_df.to_csv(dump_file_path, index=False,mode='a') 
             
    if export == True and save_plots == True:
        plot_file_name = file_beginning + "_Spikes.png"
        plot_file_name2 = file_beginning + "_Spikes.svg"
        if file_beginning != subject:
            plot_file_name = file_beginning +"_"+subject+ "_Spikes.png"
            plot_file_name2 = file_beginning +"_"+subject+ "_Spikes.svg"
        dump_plot_file_path = os.path.join(dump_path,plot_file_name)
        canvas.fig.savefig(dump_plot_file_path)
        # also as svg
        dump_plot_file_path = os.path.join(dump_path,plot_file_name2)
        canvas.fig.savefig(dump_plot_file_path, format='svg', dpi=DPI4SVG)
    else:
        canvas.draw()
    
def plot_peaks_with_event(canvas,subject,data,options,blocks_windows_values,event_name,event2_name,event_data,export,save_plots,group_name,settings_dict,export_loc_data):
    settings_dict[0]["subject"] = subject
    settings_dict[0]["subject_group_name"] = group_name
    show_norm_as = settings_dict[0]["show_norm_as"]
    dump_path,file_beginning = export_loc_data 
    # a list of three tuples. Each tuple is seconds from, seconds till
    time_windows_values = blocks_windows_values
    block1_from = time_windows_values[0][0]
    block1_till = time_windows_values[0][1]
    block2_from = time_windows_values[1][0]
    block2_till = time_windows_values[1][1]
    block3_from = time_windows_values[2][0]
    block3_till = time_windows_values[2][1]
    # check if there were any blocks
    blocks = False
    block1_spike_ts = []
    block1_spike_val = []
    block2_spike_ts = []
    block2_spike_val = []
    block3_spike_ts = []
    block3_spike_val = []
    # get peak params from gui
    height = None
    threshold = None
    distance = None
    prominence = None
    width = None
    wlen = None
    relHeight = 0.5
    plateau = None
    
    # if export was false, check also last parameter from options
    if export == False:
        export = options[-1]
        if export == True: # assume that user wants to plot as well
            save_plots = True
    
    if len(options[0]) > 0 and options[0] !="None":
        try:
            height = float(options[0]) 
        except:
            print("height problem")
    if len(options[1]) > 0 and options[1] !="None":
        try:
            threshold = float(options[1]) 
        except:
            print("threshold problem")
    if len(options[2]) > 0 and options[2] !="None":
        try:
            distance = float(options[2]) 
        except:
            print("distance problem")
    if len(options[3]) > 0 and options[3] !="None":
        try:
            prominence = float(options[3]) 
        except:
            print("prominence problem")
    if len(options[4]) > 0 and options[4] !="None":
        try:
            width = float(options[4]) 
            print("width",width)
        except:
            print("width problem")
    if len(options[5]) > 0 and options[5] !="None":
        try:
            wlen = float(options[5]) 
        except:
            print("wlen problem")
    if len(options[6]) > 0 and options[6] !="None":
        try:
            relHeight = float(options[6]) 
        except:
            print("relHeight problem")
    if len(options[7]) > 0 and options[7] !="None":
        try:
            plateau = float(options[7]) 
        except:
            print("plateau problem")
    print(height,threshold,distance,prominence,width,wlen,relHeight,plateau)
    
    ts = data["ts"]
    total_seconds = ts[-1]-ts[0]
    normalized_signal = data["normalized_signal"]
    # start time from zero
    ts_reset = [i*total_seconds/len(ts) for i in range(len(ts))]
    # distance from seconds to samples
    if distance != None and distance != "None":
        # find how many samples are per second now
        # since time is in second find when one starts
        samples_in_second = 0
        for i in range(len(ts_reset)):
            if ts_reset[i] > 1:
                samples_in_second = len(ts_reset[:i])
                break
        distance = distance*samples_in_second
    #safety check if user entered 0 they probably mean no distance
    if distance == 0:
        distance = None
    peaks, _ = find_peaks(normalized_signal, height= height, threshold = threshold, distance = distance,
                          prominence = prominence, width = width, wlen = wlen, rel_height = relHeight, plateau_size = plateau)
    # adjust peaks results to the existing timestamps
    translated_peaks = [el*total_seconds/len(normalized_signal) for el in peaks]
    total_peaks = len(translated_peaks)
    
    # plot
    # clear previous figure
    canvas.fig.clf()
    ax = canvas.fig.add_subplot(3,1,(1, 2))
    ax2 = canvas.fig.add_subplot(3,1,3,sharex=ax)
    
    ax.plot(ts_reset,normalized_signal,
             color=SIGNAL_COLOR_RGB,
             linewidth=1,
             label = "normalized\nsignal")
#    ax.plot(peaks, normalized_signal[peaks], "xr")
    ax.plot(translated_peaks, normalized_signal[peaks], "xr")
    
    total_label = str(total_peaks)+ " Total Spikes"
    ax2.eventplot(translated_peaks, orientation='horizontal', linelengths=0.5, color=CONTROL_COLOR_RGB,label=total_label)

    # get y limits of the top plot to place the block text under it
    bottom, top = ax.get_ylim()
    block_text_y_pos = bottom-0.6
    # bbox_props = dict(boxstyle='round', facecolor='k', alpha=0.1, hatch='X')
    bbox_props = dict(boxstyle='round', facecolor='k', alpha=0.1)
    # Plot windows if there are any
    if block1_from != "None" and block1_till != "None":
        print(f"From: {block1_from}, till: {block1_till}")
        blocks = True
        ax.axvline(x=block1_from,linestyle='dashed',color='k',alpha=0.5,linewidth=1)
        ax.axvline(x=block1_till,linestyle='dashed',color='k',alpha=0.5,linewidth=1)
        # ax.axvspan(block1_from,block1_till, alpha=0.1,color='k',hatch='X')
        ax.axvspan(block1_from,block1_till, alpha=0.1,color='k')
        # count spikes for that block
        for i in range(len(translated_peaks)):
            if translated_peaks[i] >= block1_from and translated_peaks[i] <= block1_till:
                block1_spike_ts.append(translated_peaks[i])
                block1_spike_val.append(normalized_signal[peaks][i])
        block1_text = str(len(block1_spike_ts))+" Spikes"
        canvas.fig.text(block1_from+1, block_text_y_pos, block1_text, ha='left', va='top', fontsize=AXIS_LABEL_FONTSIZE-2, bbox=bbox_props, transform=ax.transData)
    if block2_from != "None" and block2_till != "None":
        print(f"From: {block2_from}, till: {block2_till}")
        blocks = True
        ax.axvline(x=block2_from,linestyle='dashed',color='k',alpha=0.5,linewidth=1)
        ax.axvline(x=block2_till,linestyle='dashed',color='k',alpha=0.5,linewidth=1)
        # ax.axvspan(block2_from,block2_till, alpha=0.1,color='k',hatch='X')
        ax.axvspan(block2_from,block2_till, alpha=0.1,color='k')
        # count spikes for that block
        for i in range(len(translated_peaks)):
            if translated_peaks[i] >= block2_from and translated_peaks[i] <= block2_till:
                block2_spike_ts.append(translated_peaks[i])
                block2_spike_val.append(normalized_signal[peaks][i])
        block2_text = str(len(block2_spike_ts))+" Spikes"
        canvas.fig.text(block2_from+1, block_text_y_pos, block2_text, ha='left', va='top', fontsize=AXIS_LABEL_FONTSIZE-2, bbox=bbox_props, transform=ax.transData)
    if block3_from != "None" and block3_till != "None":
        print(f"From: {block3_from}, till: {block3_till}")
        blocks = True
        ax.axvline(x=block3_from,linestyle='dashed',color='k',alpha=0.5,linewidth=1)
        ax.axvline(x=block3_till,linestyle='dashed',color='k',alpha=0.5,linewidth=1)
        # ax.axvspan(block3_from,block3_till, alpha=0.1,color='k',hatch='X')
        ax.axvspan(block3_from,block3_till, alpha=0.1,color='k')
        # count spikes for that block
        for i in range(len(translated_peaks)):
            if translated_peaks[i] >= block3_from and translated_peaks[i] <= block3_till:
                block3_spike_ts.append(translated_peaks[i])
                block3_spike_val.append(normalized_signal[peaks][i])
        block3_text = str(len(block3_spike_ts))+" Spikes"
        canvas.fig.text(block3_from+1, block_text_y_pos, block3_text, ha='left', va='top', fontsize=AXIS_LABEL_FONTSIZE-2, bbox=bbox_props, transform=ax.transData)
    
    if len(event_data[0][0]) > 0:
        # first event
        first_evt_onset = event_data[0][0]
        first_evt_onset_translated = [el-ts[0] for el in first_evt_onset]
        first_evt_offset = event_data[0][1]
        first_evt_offset_translated = [el-ts[0] for el in first_evt_offset]
#        print(first_evt_onset)
#        print(first_evt_onset_translated)
        # plot tone onsets as vertical lines
        # first element in tone_on_off is a list with onset times
        ctr = 0
        for i in range(len(first_evt_onset_translated)):
            if first_evt_onset_translated[i] >= ts_reset[0] and first_evt_onset_translated[i] <= ts_reset[-1]:
                if ctr == 0:  # label only first one for the legend
                    ax.axvline(x=first_evt_onset_translated[i],color='k',alpha=0.5,linewidth=1,label=event_name)                   
                else:
                    ax.axvline(x=first_evt_onset_translated[i],color='k',alpha=0.5,linewidth=1)
                ctr+=1
        # if there is more than one event, and blocks were not selected, plot the second one as well
        if (len(event_data) > 1 and len(event_data[1][0]) > 0) and blocks == False:
            ctr = 0
            second_evt_onset = event_data[1][0]
            second_evt_onset_translated = [el-ts[0] for el in second_evt_onset]
            second_evt_offset = event_data[1][1]
            second_evt_offset_translated = [el-ts[0] for el in second_evt_offset]
#            print(second_evt_onset)
#            print(second_evt_onset_translated)
            for i in range(len(second_evt_onset_translated)):
                if second_evt_onset_translated[i] >= ts_reset[0] and second_evt_onset_translated[i] <= ts_reset[-1]:
                    if ctr == 0:  # label only first one for the legend
                        ax.axvline(x=second_evt_onset_translated[i],linestyle='dashed',color='k',alpha=0.5,linewidth=1,label=event2_name)                   
                    else:
                        ax.axvline(x=second_evt_onset_translated[i],linestyle='dashed',color='k',alpha=0.5,linewidth=1)
                    ctr+=1

    
    my_title = 'Subject: ' + subject
    canvas.fig.suptitle(my_title, fontsize=FIGURE_TITLE_FONT_SIZE)
    ax.set_ylabel('dF/F (%)', fontsize=AXIS_LABEL_FONTSIZE)
    if show_norm_as == "Z-Score":
        ax.set_ylabel('Z-Score', fontsize=AXIS_LABEL_FONTSIZE)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.xaxis.set_visible(False)
    ax.xaxis.set_ticks_position('none')
    ax.legend(loc=MY_LEGEND_LOC, bbox_to_anchor=MY_LEGENG_POS)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.set_ylabel('Peaks', fontsize=AXIS_LABEL_FONTSIZE)
    ax2.set_xlabel('Time (sec)', fontsize=AXIS_LABEL_FONTSIZE)
    ax2.set_yticks([])
    
    ax2.legend(loc='lower left', bbox_to_anchor=(0.85, 0))
    # set the spacing between subplots 
    canvas.fig.subplots_adjust(left=MY_LEFT, 
                    bottom=MY_BOTTOM,  
                    right=MY_RIGHT,  
                    top=MY_TOP,  
                    wspace=0,  
                    hspace=0) 
    
    # export
    if export == True:
        if file_beginning != subject:
            file_name = file_beginning + "_" + subject + "_" + event_name+"_Spikes.csv"
        else: 
            file_name = subject + "_" + event_name+"_Spikes.csv"
        if event2_name != "---" and blocks == False:
            file_name = file_beginning + "_" + event_name+"_"+event2_name+"_Spikes.csv"  
        dump_file_path = os.path.join(dump_path,file_name)

        val = "dF/F (%)"
        if show_norm_as == "Z-Score":
            val = "Z-Score"
        if blocks == False:
            raw_df = pd.DataFrame({"(sec)":translated_peaks,
                                    val: normalized_signal[peaks]})
        else: # create dataframe with joined blocks that might have different number of rows each
            dfs = [pd.DataFrame({"(sec) Total Spikes":translated_peaks,
                                    val +" Total Spikes" : normalized_signal[peaks]})]
            if block1_from != "None" and block1_till != "None":
                # add nans as missing values to the blocks, so that they all have the same number of rows as total spikes
                while len(block1_spike_ts) < total_peaks:
                    block1_spike_ts.append(np.nan)
                    block1_spike_val.append(np.nan)
                block1_ts_col_name = "(sec) Spikes from "+str(block1_from)+" till "+str(block1_till)
                block1_val_col_name = val + "Spikes from "+str(block1_from)+" till "+str(block1_till)
                df =  pd.DataFrame({
                                    block1_ts_col_name:block1_spike_ts,
                                    block1_val_col_name:block1_spike_val})
                dfs.append(df)
            if block2_from != "None" and block2_till != "None":
                while len(block2_spike_ts) < total_peaks:
                    block2_spike_ts.append(np.nan)
                    block2_spike_val.append(np.nan)
                block2_ts_col_name = "(sec) Spikes from "+str(block2_from)+" till "+str(block2_till)
                block2_val_col_name = val + "Spikes from "+str(block2_from)+" till "+str(block2_till)
                df =  pd.DataFrame({
                                    block2_ts_col_name:block2_spike_ts,
                                    block2_val_col_name:block2_spike_val})
                dfs.append(df)
            if block3_from != "None" and block3_till != "None":
                while len(block3_spike_ts) < total_peaks:
                    block3_spike_ts.append(np.nan)
                    block3_spike_val.append(np.nan)
                block3_ts_col_name = "(sec) Spikes from "+str(block3_from)+" till "+str(block3_till)
                block3_val_col_name = val + "Spikes from "+str(block3_from)+" till "+str(block3_till)
                df =  pd.DataFrame({
                                    block3_ts_col_name:block3_spike_ts,
                                    block3_val_col_name:block3_spike_val})
                dfs.append(df)
            raw_df = pd.concat(dfs,axis=1)

        settings_df = get_settings_df(settings_dict)
        settings_df.to_csv(dump_file_path,header=False)            
        with open(dump_file_path,'a') as f:
            f.write("\n")
        raw_df.to_csv(dump_file_path, index=False,mode='a') 
        event1_file = file_beginning + "_" + event_name+".csv"
        dump_file_path = os.path.join(dump_path,event1_file)
        # now save event data
        evt1_df = pd.DataFrame({"onset (s)":first_evt_onset_translated,
                                    "offset (s)" : first_evt_offset_translated})
        evt1_df.to_csv(dump_file_path, index=False)
        if (len(event_data) > 1 and len(event_data[1][0]) > 0) and blocks == False:
            event2_file = file_beginning + "_" + event2_name+".csv"
            dump_file_path = os.path.join(dump_path,event2_file)
            evt2_df = pd.DataFrame({"onset (s)":second_evt_onset_translated,
                                       "offset (s)" : second_evt_offset_translated})
            evt2_df.to_csv(dump_file_path, index=False)
             
    if export == True and save_plots == True:
        plot_file_name = file_beginning + "_" + subject + "_" + event_name+"_Spikes.png"
        plot_file_name2 = file_beginning + "_" + subject + "_" + event_name+"_Spikes.svg"
        if file_beginning == subject:
            plot_file_name = file_beginning + "_" + event_name+"_Spikes.png"
            plot_file_name2 = file_beginning + "_" + event_name+"_Spikes.svg"
        if event2_name != "---" and blocks == False:
            plot_file_name = file_beginning + "_" + subject + "_" + event_name+"_"+event2_name+"_Spikes.png" 
            plot_file_name2 = file_beginning + "_" + subject + "_" + event_name+"_"+event2_name+"_Spikes.svg" 
            if file_beginning == subject:
                plot_file_name = file_beginning + "_" + event_name+"_"+event2_name+"_Spikes.png"
                plot_file_name2 = file_beginning + "_" + event_name+"_"+event2_name+"_Spikes.svg"  
        dump_plot_file_path = os.path.join(dump_path,plot_file_name)
        canvas.fig.savefig(dump_plot_file_path)
        # also as svg
        dump_plot_file_path = os.path.join(dump_path,plot_file_name2)
        canvas.fig.savefig(dump_plot_file_path, format='svg', dpi=DPI4SVG)
    else:
        canvas.draw()

# not used now  
def get_batch_normalized(canvas,my_all_normalized,settings_dict,export,export_loc_data):
    show_norm_as = settings_dict[0]["show_norm_as"]
    dump_path,file_beginning = export_loc_data 
    all_normalized = copy.deepcopy(my_all_normalized)
    # get signals and find shortest times
    ts_lengths = []
    all_ts = []
    dfs = []
    # remember all subjects to include in file name
    subjects_end_file = ""
    for el in all_normalized:
        subject,dictionary = el
        subject_string = "_"+subject
        subjects_end_file+=subject_string
        signal = dictionary["normalized_signal"].tolist()
        ts = dictionary["ts"].tolist()
        all_ts.append(ts)
        ts_lengths.append(len(ts))
        df = pd.DataFrame({subject + " dF/F (%)":signal})
        dfs.append(df)
    all_signals_df = pd.concat(dfs, axis=1)
    all_signals_df = all_signals_df.dropna()
    transposed_df = all_signals_df.transpose()
    means_df = transposed_df.mean(axis=0)
    standard_devs_df = transposed_df.std(axis=0)
    # get means and stds as list
    means = means_df.tolist()
    stds = standard_devs_df.tolist()
#    print(all_signals_df)
#    print(transposed_df)
#    print(means_df)
#    print(standard_devs_df)
    idx = ts_lengths.index(min(ts_lengths))
    ts_shortest = all_ts[idx]
    total_seconds = ts_shortest[-1]-ts_shortest[0]
    # start time from zero
    ts_reset = [i*total_seconds/len(ts_shortest) for i in range(len(ts_shortest))]
    # assumes that mean is always calculated on all selected subjects
    se = [stds[i]/np.sqrt(all_signals_df.shape[1]) for i in range(len(stds))]
    positive_std_err_plot = [means[i]+se[i] for i in range(len(se))]
    negative_std_err_plot = [means[i]-se[i] for i in range(len(se))]
    # plot
    # clear previous figure
    canvas.fig.clf()
    ax = canvas.fig.add_subplot(111)
    ax.plot(ts_reset, means,
             color=SIGNAL_COLOR_RGB,
             linewidth=1,
             label = "normalized")
    # plot standard error bands
    ax.fill_between(ts_reset, positive_std_err_plot, negative_std_err_plot,
                      facecolor=SIGNAL_COLOR_RGB, alpha=0.2,label = 'Standard error')
    
    my_title = 'Normalized batch'
    canvas.fig.suptitle(my_title, fontsize=FIGURE_TITLE_FONT_SIZE)
    # hide top and right border
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel('Time (sec)', fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_ylabel('dF/F (%)', fontsize=AXIS_LABEL_FONTSIZE)
    my_pos = (0.94, 1)
    ax.legend(loc=MY_LEGEND_LOC, bbox_to_anchor=my_pos)
    
    # export
    if export == True:
        file_name = file_beginning + subjects_end_file+"normalized.csv"
        dump_file_path = os.path.join(dump_path,file_name)
        time_df = pd.DataFrame({"Time (sec)":ts_reset})
        means_se_df = pd.DataFrame({"Means":means,
                                    "Standard error":se})
        raw_df= pd.concat([time_df,all_signals_df,means_se_df],axis=1)
        settings_df = get_settings_df(settings_dict)
        settings_df.to_csv(dump_file_path,header=False)            
        with open(dump_file_path,'a') as f:
            f.write("\n")
        raw_df.to_csv(dump_file_path, index=False,mode='a') 
        # save plot
        plot_file_name = file_beginning + subjects_end_file+"normalized.png"
        dump_plot_file_path = os.path.join(dump_path,plot_file_name)
        canvas.fig.savefig(dump_plot_file_path)
        # also as svg
        plot_file_name = file_beginning + subjects_end_file+"normalized.svg"
        dump_plot_file_path = os.path.join(dump_path,plot_file_name)
        canvas.fig.savefig(dump_plot_file_path, format='svg', dpi=DPI4SVG)
    else:
        canvas.draw()

# not used now   
def get_batch_normalized_with_event(canvas,my_all_normalized,event_name,event2_name,event_data,settings_dict,export,export_loc_data):
    dump_path,file_beginning = export_loc_data 
    all_normalized = copy.deepcopy(my_all_normalized)
    # get signals and find shortest times
    ts_lengths = []
    all_ts = []
    dfs = []
    # remember all subjects to include in file name
    subjects_end_file = ""
    for el in all_normalized:
        subject,dictionary = el
        subject_string = "_"+subject
        subjects_end_file+=subject_string
        signal = dictionary["normalized_signal"].tolist()
        ts = dictionary["ts"].tolist()
        all_ts.append(ts)
        ts_lengths.append(len(ts))
        df = pd.DataFrame({subject + " dF/F (%)":signal})
        dfs.append(df)
    all_signals_df = pd.concat(dfs, axis=1)
    all_signals_df = all_signals_df.dropna()
    transposed_df = all_signals_df.transpose()
    means_df = transposed_df.mean(axis=0)
    standard_devs_df = transposed_df.std(axis=0)
    # get means and stds as list
    means = means_df.tolist()
    stds = standard_devs_df.tolist()
#    print(all_signals_df)
#    print(transposed_df)
#    print(means_df)
#    print(standard_devs_df)
    idx = ts_lengths.index(min(ts_lengths))
    ts_shortest = all_ts[idx]
    total_seconds = ts_shortest[-1]-ts_shortest[0]
    # start time from zero
    ts_reset = [i*total_seconds/len(ts_shortest) for i in range(len(ts_shortest))]
    # assumes that mean is always calculated on all selected subjects
    se = [stds[i]/np.sqrt(all_signals_df.shape[1]) for i in range(len(stds))]
    positive_std_err_plot = [means[i]+se[i] for i in range(len(se))]
    negative_std_err_plot = [means[i]-se[i] for i in range(len(se))]
    # plot
    # clear previous figure
    canvas.fig.clf()
    ax = canvas.fig.add_subplot(111)
    ax.plot(ts_reset, means,
             color=SIGNAL_COLOR_RGB,
             linewidth=1,
             label = "normalized")
    # plot standard error bands
    ax.fill_between(ts_reset, positive_std_err_plot, negative_std_err_plot,
                      facecolor=SIGNAL_COLOR_RGB, alpha=0.2,label = 'Standard error')
    
    if len(event_data[0][0]) > 0:
        # first event
        first_evt_onset = event_data[0][0]
        first_evt_onset_translated = [el-ts[0] for el in first_evt_onset]
        first_evt_offset = event_data[0][1]
        first_evt_offset_translated = [el-ts[0] for el in first_evt_offset]
#        print(first_evt_onset)
#        print(first_evt_onset_translated)
        # plot tone onsets as vertical lines
        # first element in tone_on_off is a list with onset times
        ctr = 0
        for i in range(len(first_evt_onset_translated)):
            if first_evt_onset_translated[i] >= ts_reset[0] and first_evt_onset_translated[i] <= ts_reset[-1]:
                if ctr == 0:  # label only first one for the legend
                    ax.axvline(x=first_evt_onset_translated[i],color='k',alpha=0.5,linewidth=1,label=event_name)                   
                else:
                    ax.axvline(x=first_evt_onset_translated[i],color='k',alpha=0.5,linewidth=1)
                ctr+=1
        # if there is more than one event, plot the second one as well
        if len(event_data) > 1 and len(event_data[1][0]) > 0:
            ctr = 0
            second_evt_onset = event_data[1][0]
            second_evt_onset_translated = [el-ts[0] for el in second_evt_onset]
            second_evt_offset = event_data[1][1]
            second_evt_offset_translated = [el-ts[0] for el in second_evt_offset]
#            print(second_evt_onset)
#            print(second_evt_onset_translated)
            for i in range(len(second_evt_onset_translated)):
                if second_evt_onset_translated[i] >= ts_reset[0] and second_evt_onset_translated[i] <= ts_reset[-1]:
                    if ctr == 0:  # label only first one for the legend
                        ax.axvline(x=second_evt_onset_translated[i],linestyle='dashed',color='k',alpha=0.5,linewidth=1,label=event2_name)                   
                    else:
                        ax.axvline(x=second_evt_onset_translated[i],linestyle='dashed',color='k',alpha=0.5,linewidth=1)
                    ctr+=1
    
    my_title = 'Normalized batch'
    canvas.fig.suptitle(my_title, fontsize=FIGURE_TITLE_FONT_SIZE)
    # hide top and right border
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel('Time (sec)', fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_ylabel('dF/F (%)', fontsize=AXIS_LABEL_FONTSIZE)
    my_pos = (0.94, 1)
    ax.legend(loc=MY_LEGEND_LOC, bbox_to_anchor=my_pos)
    
    # export
    if export == True:
        file_name = file_beginning + subjects_end_file + "_" + event_name+"_normalized.csv"
        if event2_name != "---":
            file_name = file_beginning + subjects_end_file + "_" + event_name+"_"+event2_name+"_normalized.csv"  
        dump_file_path = os.path.join(dump_path,file_name)
        time_df = pd.DataFrame({"Time (sec)":ts_reset})
        means_se_df = pd.DataFrame({"Means":means,
                                    "Standard error":se})
        raw_df= pd.concat([time_df,all_signals_df,means_se_df],axis=1)
        settings_df = get_settings_df(settings_dict)
        settings_df.to_csv(dump_file_path,header=False)            
        with open(dump_file_path,'a') as f:
            f.write("\n")
        raw_df.to_csv(dump_file_path, index=False,mode='a') 
        event1_file = file_beginning + subjects_end_file + "_" + event_name+"normalized.csv"
        dump_file_path = os.path.join(dump_path,event1_file)
        # now save event data
        evt1_df = pd.DataFrame({"onset (s)":first_evt_onset_translated,
                                    "offset (s)" : first_evt_offset_translated})
        evt1_df.to_csv(dump_file_path, index=False)
        if len(event_data) > 1 and len(event_data[1][0]) > 0:
            event2_file = file_beginning + subjects_end_file + "_" + event2_name+"normalized.csv"
            dump_file_path = os.path.join(dump_path,event2_file)
            evt2_df = pd.DataFrame({"onset (s)":second_evt_onset_translated,
                                       "offset (s)" : second_evt_offset_translated})
            evt2_df.to_csv(dump_file_path, index=False)
        # save plot
        plot_file_name = file_beginning + subjects_end_file + "_" + event_name+"_normalized.png"
        plot_file_name2 = file_beginning + subjects_end_file + "_" + event_name+"_normalized.svg"
        if event2_name != "---":
            plot_file_name = file_beginning + subjects_end_file + "_" + event_name+"_"+event2_name+"_normalized.png" 
            plot_file_name2 = file_beginning + subjects_end_file + "_" + event_name+"_"+event2_name+"_normalized.svg"  
        dump_plot_file_path = os.path.join(dump_path,plot_file_name)
        canvas.fig.savefig(dump_plot_file_path)
        # also as svg
        dump_plot_file_path = os.path.join(dump_path,plot_file_name2)
        canvas.fig.savefig(dump_plot_file_path, format='svg', dpi=DPI4SVG)
    else:
        canvas.draw()
        
# not used
def get_batch_spikes(canvas,options,my_all_normalized,settings_dict,export,export_loc_data):
    dump_path,file_beginning = export_loc_data 
    all_normalized = copy.deepcopy(my_all_normalized)
    # get peak params from gui
     # get peak params from gui
    height = None
    threshold = None
    distance = None
    prominence = None
    width = None
    wlen = None
    relHeight = 0.5
    plateau = None

    
    if len(options[0]) > 0 and options[0] !="None":
        try:
            height = float(options[0]) 
        except:
            print("height problem")
    if len(options[1]) > 0 and options[1] !="None":
        try:
            threshold = float(options[1]) 
        except:
            print("threshold problem")
    if len(options[2]) > 0 and options[2] !="None":
        try:
            distance = float(options[2]) 
        except:
            print("distance problem")
    if len(options[3]) > 0 and options[3] !="None":
        try:
            prominence = float(options[3]) 
        except:
            print("prominence problem")
    if len(options[4]) > 0 and options[4] !="None":
        try:
            width = float(options[4]) 
            print("width",width)
        except:
            print("width problem")
    if len(options[5]) > 0 and options[5] !="None":
        try:
            wlen = float(options[5]) 
        except:
            print("wlen problem")
    if len(options[6]) > 0 and options[6] !="None":
        try:
            relHeight = float(options[6]) 
        except:
            print("relHeight problem")
    if len(options[7]) > 0 and options[7] !="None":
        try:
            plateau = float(options[7]) 
        except:
            print("plateau problem")
        
    # get signals and find shortest times
    ts_lengths = []
    all_ts = []
    dfs = []
    # remember all subjects to include in file name
    subjects_end_file = ""
    for el in all_normalized:
        subject,dictionary = el
        subject_string = "_"+subject
        subjects_end_file+=subject_string
        signal = dictionary["normalized_signal"].tolist()
        ts = dictionary["ts"].tolist()
        all_ts.append(ts)
        ts_lengths.append(len(ts))
        df = pd.DataFrame({subject + "dF/F (%)":signal})
        dfs.append(df)
    all_signals_df = pd.concat(dfs, axis=1)
    all_signals_df = all_signals_df.dropna()
    transposed_df = all_signals_df.transpose()
    means_df = transposed_df.mean(axis=0)
    standard_devs_df = transposed_df.std(axis=0)
    # get means and stds as list
    means = means_df.tolist()
    stds = standard_devs_df.tolist()
#    print(all_signals_df)
#    print(transposed_df)
#    print(means_df)
#    print(standard_devs_df)
    idx = ts_lengths.index(min(ts_lengths))
    ts_shortest = all_ts[idx]
    total_seconds = ts_shortest[-1]-ts_shortest[0]
    # start time from zero
    ts_reset = [i*total_seconds/len(ts_shortest) for i in range(len(ts_shortest))]
    # assumes that mean is always calculated on all selected subjects
    se = [stds[i]/np.sqrt(all_signals_df.shape[1]) for i in range(len(stds))]
    positive_std_err_plot = [means[i]+se[i] for i in range(len(se))]
    negative_std_err_plot = [means[i]-se[i] for i in range(len(se))]
    
    # distance from seconds to samples
    if distance != None and distance != "None":
        # find how many samples are per second now
        # since time is in second find when one starts
        samples_in_second = 0
        for i in range(len(ts_reset)):
            if ts_reset[i] > 1:
                samples_in_second = len(ts_reset[:i])
                break
        distance = distance*samples_in_second
    #safety check if user entered 0 they probably mean no distance
    if distance == 0:
        distance = None
#    print(ts_reset[-1])
    peaks, _ = find_peaks(means, height= height, threshold = threshold, distance = distance,
                          prominence = prominence, width = width, wlen = wlen, rel_height = relHeight, plateau_size = plateau)
#    print("peaks\n",type(peaks[0]))
    # adjust peaks results to the existing timestamps
    translated_peaks = [el*total_seconds/len(means) for el in peaks]
    total_peaks = len(translated_peaks)
#    print(translated_peaks)
    
    # plot
    # clear previous figure
    canvas.fig.clf()
    ax = canvas.fig.add_subplot(3,1,(1, 2))
    ax2 = canvas.fig.add_subplot(3,1,3,sharex=ax)
    
    ax.plot(ts_reset,means,
             color=SIGNAL_COLOR_RGB,
             linewidth=1,
             label = "signal")

    ax.plot(translated_peaks, np.array(means)[peaks], "xr")
    # plot standard error bands
    ax.fill_between(ts_reset, positive_std_err_plot, negative_std_err_plot,
                      facecolor=CONTROL_COLOR_RGB, alpha=0.2,label = 'Standard error')
    
    total_label = str(total_peaks)+ " Spikes Found"
    ax2.eventplot(translated_peaks, orientation='horizontal', linelengths=0.5, color=CONTROL_COLOR_RGB,label = total_label)
    
    my_title = 'Batch Spikes'
    canvas.fig.suptitle(my_title, fontsize=FIGURE_TITLE_FONT_SIZE)
    ax.set_ylabel('dF/F (%)', fontsize=AXIS_LABEL_FONTSIZE)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.xaxis.set_visible(False)
    ax.xaxis.set_ticks_position('none')
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.set_ylabel('Peaks', fontsize=AXIS_LABEL_FONTSIZE)
    ax2.set_xlabel('Time (sec)', fontsize=AXIS_LABEL_FONTSIZE)
    ax2.set_yticks([])
    
    ax2.legend(loc='lower left', bbox_to_anchor=(0.85, 0))
    # set the spacing between subplots 
    canvas.fig.subplots_adjust(left=MY_LEFT, 
                    bottom=MY_BOTTOM,  
                    right=0.95,  
                    top=MY_TOP,  
                    wspace=0,  
                    hspace=0) 
    
    
    # export
    if export == True:
        file_name = file_beginning + subjects_end_file+"_Spikes.csv"
        dump_file_path = os.path.join(dump_path,file_name)
        raw_df= pd.DataFrame({"Time (sec)":translated_peaks,
                                 "Value dF/F (%)" : np.array(means)[peaks]})
        settings_df = get_settings_df(settings_dict)
        settings_df.to_csv(dump_file_path,header=False)            
        with open(dump_file_path,'a') as f:
            f.write("\n")
        raw_df.to_csv(dump_file_path, index=False,mode='a') 
        # save plot
        plot_file_name = file_beginning + subjects_end_file+"_Spikes.png"
        dump_plot_file_path = os.path.join(dump_path,plot_file_name)
        canvas.fig.savefig(dump_plot_file_path)
        # also as svg
        plot_file_name = file_beginning + subjects_end_file+"_Spikes.svg"
        dump_plot_file_path = os.path.join(dump_path,plot_file_name)
        canvas.fig.savefig(dump_plot_file_path, format='svg', dpi=DPI4SVG)
    else:
        canvas.draw()
# not used       
def get_batch_spikes_with_event(canvas,options,my_all_normalized,event_name,event2_name,event_data,settings_dict,export,export_loc_data):
    dump_path,file_beginning = export_loc_data 
    all_normalized = copy.deepcopy(my_all_normalized)
    # get peak params from gui
     # get peak params from gui
    height = None
    threshold = None
    distance = None
    prominence = None
    width = None
    wlen = None
    relHeight = 0.5
    plateau = None

    
    if len(options[0]) > 0 and options[0] !="None":
        try:
            height = float(options[0]) 
        except:
            print("height problem")
    if len(options[1]) > 0 and options[1] !="None":
        try:
            threshold = float(options[1]) 
        except:
            print("threshold problem")
    if len(options[2]) > 0 and options[2] !="None":
        try:
            distance = float(options[2]) 
        except:
            print("distance problem")
    if len(options[3]) > 0 and options[3] !="None":
        try:
            prominence = float(options[3]) 
        except:
            print("prominence problem")
    if len(options[4]) > 0 and options[4] !="None":
        try:
            width = float(options[4]) 
            print("width",width)
        except:
            print("width problem")
    if len(options[5]) > 0 and options[5] !="None":
        try:
            wlen = float(options[5]) 
        except:
            print("wlen problem")
    if len(options[6]) > 0 and options[6] !="None":
        try:
            relHeight = float(options[6]) 
        except:
            print("relHeight problem")
    if len(options[7]) > 0 and options[7] !="None":
        try:
            plateau = float(options[7]) 
        except:
            print("plateau problem")
        
    # get signals and find shortest times
    ts_lengths = []
    all_ts = []
    dfs = []
    # remember all subjects to include in file name
    subjects_end_file = ""
    for el in all_normalized:
        subject,dictionary = el
        subject_string = "_"+subject
        subjects_end_file+=subject_string
        signal = dictionary["normalized_signal"].tolist()
        ts = dictionary["ts"].tolist()
        all_ts.append(ts)
        ts_lengths.append(len(ts))
        df = pd.DataFrame({subject + "dF/F (%)":signal})
        dfs.append(df)
    all_signals_df = pd.concat(dfs, axis=1)
    all_signals_df = all_signals_df.dropna()
    transposed_df = all_signals_df.transpose()
    means_df = transposed_df.mean(axis=0)
    standard_devs_df = transposed_df.std(axis=0)
    # get means and stds as list
    means = means_df.tolist()
    stds = standard_devs_df.tolist()
#    print(all_signals_df)
#    print(transposed_df)
#    print(means_df)
#    print(standard_devs_df)
    idx = ts_lengths.index(min(ts_lengths))
    ts_shortest = all_ts[idx]
    total_seconds = ts_shortest[-1]-ts_shortest[0]
    # start time from zero
    ts_reset = [i*total_seconds/len(ts_shortest) for i in range(len(ts_shortest))]
    # assumes that mean is always calculated on all selected subjects
    se = [stds[i]/np.sqrt(all_signals_df.shape[1]) for i in range(len(stds))]
    positive_std_err_plot = [means[i]+se[i] for i in range(len(se))]
    negative_std_err_plot = [means[i]-se[i] for i in range(len(se))]
    
    # distance from seconds to samples
    if distance != None and distance != "None":
        # find how many samples are per second now
        # since time is in second find when one starts
        samples_in_second = 0
        for i in range(len(ts_reset)):
            if ts_reset[i] > 1:
                samples_in_second = len(ts_reset[:i])
                break
        distance = distance*samples_in_second
    #safety check if user entered 0 they probably mean no distance
    if distance == 0:
        distance = None
#    print(ts_reset[-1])
    peaks, _ = find_peaks(means, height= height, threshold = threshold, distance = distance,
                          prominence = prominence, width = width, wlen = wlen, rel_height = relHeight, plateau_size = plateau)
#    print("peaks\n",type(peaks[0]))
    # adjust peaks results to the existing timestamps
    translated_peaks = [el*total_seconds/len(means) for el in peaks]
    total_peaks = len(translated_peaks)
#    print(translated_peaks)
    
    # plot
    # clear previous figure
    canvas.fig.clf()
    ax = canvas.fig.add_subplot(3,1,(1, 2))
    ax2 = canvas.fig.add_subplot(3,1,3,sharex=ax)
    
    ax.plot(ts_reset,means,
             color=SIGNAL_COLOR_RGB,
             linewidth=1,
             label = "signal")

    ax.plot(translated_peaks, np.array(means)[peaks], "xr")
    # plot standard error bands
    ax.fill_between(ts_reset, positive_std_err_plot, negative_std_err_plot,
                      facecolor=CONTROL_COLOR_RGB, alpha=0.2,label = 'Standard error')
    
    total_label = str(total_peaks)+ " Spikes Found"
    ax2.eventplot(translated_peaks, orientation='horizontal', linelengths=0.5, color=CONTROL_COLOR_RGB,label = total_label)
    
    if len(event_data[0][0]) > 0:
        # first event
        first_evt_onset = event_data[0][0]
        first_evt_onset_translated = [el-ts[0] for el in first_evt_onset]
        first_evt_offset = event_data[0][1]
        first_evt_offset_translated = [el-ts[0] for el in first_evt_offset]
#        print(first_evt_onset)
#        print(first_evt_onset_translated)
        # plot tone onsets as vertical lines
        # first element in tone_on_off is a list with onset times
        ctr = 0
        for i in range(len(first_evt_onset_translated)):
            if first_evt_onset_translated[i] >= ts_reset[0] and first_evt_onset_translated[i] <= ts_reset[-1]:
                if ctr == 0:  # label only first one for the legend
                    ax.axvline(x=first_evt_onset_translated[i],color='k',alpha=0.5,linewidth=1,label=event_name)                   
                else:
                    ax.axvline(x=first_evt_onset_translated[i],color='k',alpha=0.5,linewidth=1)
                ctr+=1
        # if there is more than one event, plot the second one as well
        if len(event_data) > 1 and len(event_data[1][0]) > 0:
            ctr = 0
            second_evt_onset = event_data[1][0]
            second_evt_onset_translated = [el-ts[0] for el in second_evt_onset]
            second_evt_offset = event_data[1][1]
            second_evt_offset_translated = [el-ts[0] for el in second_evt_offset]
#            print(second_evt_onset)
#            print(second_evt_onset_translated)
            for i in range(len(second_evt_onset_translated)):
                if second_evt_onset_translated[i] >= ts_reset[0] and second_evt_onset_translated[i] <= ts_reset[-1]:
                    if ctr == 0:  # label only first one for the legend
                        ax.axvline(x=second_evt_onset_translated[i],linestyle='dashed',color='k',alpha=0.5,linewidth=1,label=event2_name)                   
                    else:
                        ax.axvline(x=second_evt_onset_translated[i],linestyle='dashed',color='k',alpha=0.5,linewidth=1)
                    ctr+=1
    
    my_title = 'Batch Spikes Subjects: ' + " ".join(subjects_end_file.split("_"))
    canvas.fig.suptitle(my_title, fontsize=FIGURE_TITLE_FONT_SIZE)
    ax.set_ylabel('dF/F (%)', fontsize=AXIS_LABEL_FONTSIZE)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.xaxis.set_visible(False)
    ax.xaxis.set_ticks_position('none')
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.set_ylabel('Peaks', fontsize=AXIS_LABEL_FONTSIZE)
    ax2.set_xlabel('Time (sec)', fontsize=AXIS_LABEL_FONTSIZE)
    ax2.set_yticks([])
    
    ax2.legend(loc='lower left', bbox_to_anchor=(0.85, 0))
    # set the spacing between subplots 
    canvas.fig.subplots_adjust(left=MY_LEFT, 
                    bottom=MY_BOTTOM,  
                    right=0.95,  
                    top=MY_TOP,  
                    wspace=0,  
                    hspace=0) 
       
    # export
    if export == True:
        file_name = file_beginning + subjects_end_file+"_Spikes.csv"
        dump_file_path = os.path.join(dump_path,file_name)
        raw_df= pd.DataFrame({"Time (sec)":translated_peaks,
                                 "Value dF/F (%)" : np.array(means)[peaks]})
        settings_df = get_settings_df(settings_dict)
        settings_df.to_csv(dump_file_path,header=False)            
        with open(dump_file_path,'a') as f:
            f.write("\n")
        raw_df.to_csv(dump_file_path, index=False,mode='a') 
        event1_file = file_beginning + subjects_end_file + "_" + event_name+".csv"
        dump_file_path = os.path.join(dump_path,event1_file)
        # now save event data
        evt1_df = pd.DataFrame({"onset (s)":first_evt_onset_translated,
                                    "offset (s)" : first_evt_offset_translated})
        evt1_df.to_csv(dump_file_path, index=False)
        if len(event_data) > 1 and len(event_data[1][0]) > 0:
            event2_file = file_beginning + subjects_end_file + "_" + event2_name+".csv"
            dump_file_path = os.path.join(dump_path,event2_file)
            evt2_df = pd.DataFrame({"onset (s)":second_evt_onset_translated,
                                       "offset (s)" : second_evt_offset_translated})
            evt2_df.to_csv(dump_file_path, index=False)
        # save plot
        plot_file_name = file_beginning + subjects_end_file+"_Spikes.png"
        dump_plot_file_path = os.path.join(dump_path,plot_file_name)
        canvas.fig.savefig(dump_plot_file_path)
        # also as svg
        plot_file_name = file_beginning + subjects_end_file+"_Spikes.svg"
        dump_plot_file_path = os.path.join(dump_path,plot_file_name)
        canvas.fig.savefig(dump_plot_file_path, format='svg', dpi=DPI4SVG)
    else:
        canvas.draw()
        

def get_batch_perievent_normalized(canvas,my_all_dfs,group_names,perievent_options_dict,settings_dict,export_loc_data):
    show_norm_as = settings_dict[0]["show_norm_as"]
    event_name = perievent_options_dict["event_name"] if len(perievent_options_dict["event_name"]) > 0 else perievent_options_dict["event"]
    settings_dict[0]["baseline_from_sec"] =perievent_options_dict["baseline_from"]
    settings_dict[0]["baseline_to_sec"] =perievent_options_dict["baseline_to"]
    settings_dict[0]["auc_pre_from"] =perievent_options_dict["auc_pre_from"]
    settings_dict[0]["auc_pre_to"] =perievent_options_dict["auc_pre_to"]
    settings_dict[0]["auc_post_from"] =perievent_options_dict["auc_post_from"]
    settings_dict[0]["auc_post_to"] =perievent_options_dict["auc_post_to"]
    dump_path,file_beginning = export_loc_data 
    all_dfs = copy.deepcopy(my_all_dfs)
    # contains single dataframe with time
    ts_df = []
     # contains a list of all subject's dataframe(with all normalized trials)
    trials_dfs = []
    # remember all subjects
    subjects_end_file = ""
    # and their group names
    subjects_group_names = ""
    # separate dataframe to show subjects and groups above headers
    group_df = [pd.DataFrame({"headers":["Subject","Group"]})]
    group_trials_dfs = []
    group_idx = 0
    for el in all_dfs:
        subject,df = el
        subject_string = subject+","
        subjects_end_file+=subject_string
        # if first dataframe keep time
        if len(ts_df) == 0:
            ts = df.iloc[:,:1]
            ts_df.append(ts)
        # remove time
        del df['Time (sec)']
        original_cols = list(df.columns)
        new_dfs_list = []
        for i in range(len(original_cols)):
            single_df = pd.DataFrame({i:[subject,group_names[group_idx]]})
            new_dfs_list.append(single_df)
        group_trials_dfs.append(pd.concat(new_dfs_list,axis=1))
        group_idx +=1
        trials_dfs.append(df)
    # add subjects and their group names 
    group_df.extend(group_trials_dfs)
    # create a single dataframe from that
    group_header_df4csv = pd.concat(group_df,axis=1)
    # create headers with group names for standard header
    for name in group_names:
        name2write = name+","
        subjects_group_names+=name2write

    # add subjects to settings
    # exclude last coma
    settings_dict[0]["subject"] = subjects_end_file[:-1]
    settings_dict[0]["subject_group_name"] = subjects_group_names[:-1]
    # create df with all subject info for csv
    ts_df.extend(trials_dfs)
    main_df = pd.concat(ts_df,axis=1)  
    
    # calculate mean across all trials
    means_form_all_trials_df = pd.concat(trials_dfs,axis=1)
    transposed_df = means_form_all_trials_df.transpose()
    mean_by_subject_df = transposed_df.mean(axis=0)
    mean_by_subject = mean_by_subject_df.tolist()
    mean_by_subject_df = pd.DataFrame({"Mean":mean_by_subject})
    standard_devs_df = transposed_df.std(axis=0)
    stds = standard_devs_df.tolist()
    se = [stds[i]/np.sqrt(means_form_all_trials_df.shape[1]) for i in range(len(stds))]
    se_df = pd.DataFrame({"Error":se})
    # add to df
    main_df = pd.concat([main_df,mean_by_subject_df,se_df],axis=1)
    
    positive_std_err_plot = [mean_by_subject[i]+se[i] for i in range(len(se))]
    negative_std_err_plot = [mean_by_subject[i]-se[i] for i in range(len(se))]
    # plot
    # clear previous figure
    canvas.fig.clf()
    ax = canvas.fig.add_subplot(111)
    ts = main_df["Time (sec)"].tolist()
    ax.plot(ts, mean_by_subject,
             color=SIGNAL_COLOR_RGB,
             linewidth=1,
             label = "Average")
    # plot a vertical line at t=0
    ax.axvline(x=0, linewidth=2, color='k', alpha=0.3,label=event_name)
    # plot standard error bands
    ax.fill_between(ts, positive_std_err_plot, negative_std_err_plot,
                      facecolor=SIGNAL_COLOR_RGB, alpha=0.2,label = 'Standard error')
    
    my_title = 'Total Subjects: ' + str(len(all_dfs))
    canvas.fig.suptitle(my_title, fontsize=FIGURE_TITLE_FONT_SIZE)
    # hide top and right border
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel('Time (sec)', fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_ylabel('dF/F (%)', fontsize=AXIS_LABEL_FONTSIZE)
    if show_norm_as == "Z-Score":
        ax.set_ylabel("Z-Score",fontsize=AXIS_LABEL_FONTSIZE)
    ax.legend(loc=MY_LEGEND_LOC, bbox_to_anchor=MY_LEGENG_POS)

    # save csv with data
    file_name = file_beginning +"_perievent_avg.csv"
    dump_file_path = os.path.join(dump_path,file_name)
    settings_df = get_settings_df(settings_dict)
    settings_df.to_csv(dump_file_path,header=False)            
    with open(dump_file_path,'a') as f:
        f.write("\n")
    group_header_df4csv.to_csv(dump_file_path, header=False,index=False,mode='a')
    main_df.to_csv(dump_file_path, index=False,mode='a')
      
    # save plot
    plot_file_name = file_beginning + "_perievent_avg.png"
    dump_plot_file_path = os.path.join(dump_path,plot_file_name)
    canvas.fig.savefig(dump_plot_file_path)
    # also as svg
    # plot_file_name = file_beginning + subjects_end_file+"_perievent_avg.svg"
    plot_file_name = file_beginning +"_perievent_avg.svg"
    dump_plot_file_path = os.path.join(dump_path,plot_file_name)
    canvas.fig.savefig(dump_plot_file_path, format='svg', dpi=DPI4SVG)
    
def get_batch_perievent_zscored(canvas,my_all_dfs,group_names,perievent_options_dict,settings_dict,export_loc_data):
    event_name = perievent_options_dict["event_name"] if len(perievent_options_dict["event_name"]) > 0 else perievent_options_dict["event"]
    settings_dict[0]["baseline_from_sec"] =perievent_options_dict["baseline_from"]
    settings_dict[0]["baseline_to_sec"] =perievent_options_dict["baseline_to"]
    settings_dict[0]["auc_pre_from"] =perievent_options_dict["auc_pre_from"]
    settings_dict[0]["auc_pre_to"] =perievent_options_dict["auc_pre_to"]
    settings_dict[0]["auc_post_from"] =perievent_options_dict["auc_post_from"]
    settings_dict[0]["auc_post_to"] =perievent_options_dict["auc_post_to"]
    dump_path,file_beginning = export_loc_data 
#    print(len(my_all_dfs))
    all_dfs = copy.deepcopy(my_all_dfs)
    # contains single dataframe with time
    ts_df = []
    # contains a list of all subject's dataframe(with all normalized trials)
    zscore_all_dfs = []
    zscores_means_by_subject = []
    # each inner list will contain same trial from different subject
    zscores_by_trials = []
    # remember all subjects to include in file name
    subjects_end_file = ""
    subjects_group_names = ""
    # separate dataframe to show subjects and groups above headers
    group_df = [pd.DataFrame({"headers":["Subject","Group"]})]
    group_trials_dfs = []
    group_idx = 0
    for i in range(len(all_dfs)):
        subject,df = all_dfs[i]
#        print(subject)
#        print(df)
        subject_string = subject + ","
        subjects_end_file+=subject_string
        # get means as list
        means = df['Mean_zscore'].tolist()
        zscores_means_by_subject.append(means)
        del df['Mean_zscore']
        del df['Error']
        # if first dataframe keep time
        if len(ts_df) == 0:
            ts = df.iloc[:,:1]
            ts_df.append(ts)
        # remove time
        del df['Time (sec)']
        original_cols = list(df.columns)
        new_dfs_list = []
        for i in range(len(original_cols)):
            single_df = pd.DataFrame({i:[subject,group_names[group_idx]]})
            new_dfs_list.append(single_df)
        group_trials_dfs.append(pd.concat(new_dfs_list,axis=1))
        group_idx +=1
        # add to means list
        zscore_all_dfs.append(df) 
        # add each trial column to a separate list
        for j in range(len(original_cols)):
            trial_df = df[original_cols[j]]
            try:
                zscores_by_trials[j].append(trial_df)
            except: # handle index out of range exception
                # add empty list for trials
               zscores_by_trials.append([]) 
               zscores_by_trials[j].append(trial_df)
    # add subjects and their group names 
    group_df.extend(group_trials_dfs)
    # create a single dataframe from that
    group_header_df4csv = pd.concat(group_df,axis=1)
    # create headers with group names for standard header
    for name in group_names:
        name2write = name+","
        subjects_group_names+=name2write

    # add subjects to settings
    # exclude last coma
    settings_dict[0]["subject"] = subjects_end_file[:-1]
    settings_dict[0]["subject_group_name"] = subjects_group_names[:-1]
    
    # calculate means
    # join all trials from all subjects to get means ans errors
    df_all_trials = pd.concat(zscore_all_dfs,axis=1)
    transposed_df = df_all_trials.transpose()
    means_df = transposed_df.mean(axis=0)
    means_df = pd.DataFrame({"Mean":means_df.tolist()})
    standard_devs_df = transposed_df.std(axis=0)
    stds = standard_devs_df.tolist()
    se = [stds[i]/np.sqrt(df_all_trials.shape[1]) for i in range(len(stds))]
    se_df = pd.DataFrame({"Error":se})
    
    # calculate means by trial for the heatmap
    zscores_by_trials_list = []
    for i in range(len(zscores_by_trials)):
        single_trial_df = pd.concat(zscores_by_trials[i],axis=1)
        transposed = single_trial_df.transpose()
        means = transposed.mean(axis=0)
        zscores_by_trials_list.append(means.tolist())       
    
    # create main dataframe
    raw_df= pd.concat([ts_df[0],df_all_trials,means_df,se_df],axis=1)

    # get data to plot
    ts = ts_df[0]['Time (sec)'].tolist()
    mean_zscore = means_df['Mean'].tolist()
    positive_std_err_plot = [mean_zscore[i]+se[i] for i in range(len(se))]
    negative_std_err_plot = [mean_zscore[i]-se[i] for i in range(len(se))]
#    
    # plot
    # clear previous figure
    canvas.fig.clf()
    ax = canvas.fig.add_subplot(211)
    ax2 = canvas.fig.add_subplot(212)
    
    # Heat Map based on z score of control fit subtracted signal
    cs = ax.imshow(zscores_by_trials_list, cmap=plt.cm.bwr, interpolation='none', aspect="auto",
                    extent=[-perievent_options_dict['sec_before'], perievent_options_dict['sec_after'], 
                            len(zscores_by_trials_list),0])
    canvas.fig.colorbar(cs, ax=ax,pad=0.01, fraction=0.02)
    # plot the z-score trace for the signal with std error bands
    ax2.plot(ts, mean_zscore, linewidth=2, color=SIGNAL_COLOR_RGB, label='Group Mean')
    ax2.fill_between(ts, positive_std_err_plot, negative_std_err_plot,
                      facecolor=SIGNAL_COLOR_RGB, alpha=0.2,label = 'Standard error')

    ax2.axvline(x=0, linewidth=2, color='k', alpha=0.3, label=event_name)
  
    my_title = 'Total Subjects: ' + str(len(all_dfs)) + "; Event: " + event_name
    canvas.fig.suptitle(my_title, fontsize=FIGURE_TITLE_FONT_SIZE)
    ax.set_title('z-Score Traces',fontsize=TITLE_FONT_SIZE)
    ax.set_ylabel('Trials',fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_xlabel('Seconds from Event Onset',fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_yticks(np.arange(0.5,len(zscores_by_trials_list), 1))
    ax.set_yticklabels(np.arange(1, len(zscores_by_trials_list)+1, 1))
    
    ax2.set_ylabel('z-Score',fontsize=AXIS_LABEL_FONTSIZE)
    ax2.set_xlabel('Seconds',fontsize=AXIS_LABEL_FONTSIZE)
    # hide top and right border
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.legend(loc=MY_LEGEND_LOC, bbox_to_anchor=MY_LEGENG_POS)
    ax2.set_xlim([-perievent_options_dict['sec_before'], perievent_options_dict['sec_after']])
    
    # set the spacing between subplots 
    canvas.fig.subplots_adjust(left=MY_LEFT, 
                    bottom=MY_BOTTOM,  
                    right=MY_RIGHT,  
                    top=MY_TOP,  
                    wspace=MY_WSPACE,  
                    hspace=MY_HSPACE) 
    # save csv with data
    file_name = file_beginning +"_perievent_zscore.csv"
    dump_file_path = os.path.join(dump_path,file_name)
    settings_df = get_settings_df(settings_dict)
    settings_df.to_csv(dump_file_path,header=False)            
    with open(dump_file_path,'a') as f:
        f.write("\n")
    group_header_df4csv.to_csv(dump_file_path, header=False,index=False,mode='a')
    raw_df.to_csv(dump_file_path, index=False,mode='a')
      
    # save plot
    plot_file_name = file_beginning +"_perievent_zscore.png"
    dump_plot_file_path = os.path.join(dump_path,plot_file_name)
    canvas.fig.savefig(dump_plot_file_path)
    # also as svg
    # plot_file_name = file_beginning + subjects_end_file+"_perievent_zscore.svg"
    plot_file_name = file_beginning +"_perievent_zscore.svg"
    dump_plot_file_path = os.path.join(dump_path,plot_file_name)
    canvas.fig.savefig(dump_plot_file_path, format='svg', dpi=DPI4SVG)

    
def get_batch_perievent_zscored_with_trials(canvas,my_all_dfs,group_names,perievent_options_dict,settings_dict,export_loc_data):
    event_name = perievent_options_dict["event_name"] if len(perievent_options_dict["event_name"]) > 0 else perievent_options_dict["event"]
    settings_dict[0]["baseline_from_sec"] =perievent_options_dict["baseline_from"]
    settings_dict[0]["baseline_to_sec"] =perievent_options_dict["baseline_to"]
    settings_dict[0]["auc_pre_from"] =perievent_options_dict["auc_pre_from"]
    settings_dict[0]["auc_pre_to"] =perievent_options_dict["auc_pre_to"]
    settings_dict[0]["auc_post_from"] =perievent_options_dict["auc_post_from"]
    settings_dict[0]["auc_post_to"] =perievent_options_dict["auc_post_to"]
    dump_path,file_beginning = export_loc_data 
#    print(len(my_all_dfs))
    all_dfs = copy.deepcopy(my_all_dfs)
    # contains single dataframe with time
    ts_df = []
    # contains a list of all subject's dataframe(with all normalized trials)
    zscore_all_dfs = []
    zscores_means_by_subject = []
    # each inner list will contain same trial from different subject
    zscores_by_trials = []
    # remember all subjects to include in file name
    subjects_end_file = ""
    subjects_group_names = ""
    # separate dataframe to show subjects and groups above headers
    group_df = [pd.DataFrame({"headers":["Subject","Group"]})]
    group_trials_dfs = []
    group_idx = 0
    for el in all_dfs:
        subject,df = el
#        print(subject)
#        print(df)
        subject_string = subject + ","
        subjects_end_file+=subject_string
        # get means as list
        means = df['Mean_zscore'].tolist()
        zscores_means_by_subject.append(means)
        del df['Mean_zscore']
        del df['Error']
        # if first dataframe keep time
        if len(ts_df) == 0:
            ts = df.iloc[:,:1]
            ts_df.append(ts)
        # remove time
        del df['Time (sec)']
        original_cols = list(df.columns)
        new_dfs_list = []
        for i in range(len(original_cols)):
            single_df = pd.DataFrame({i:[subject,group_names[group_idx]]})
            new_dfs_list.append(single_df)
        group_trials_dfs.append(pd.concat(new_dfs_list,axis=1))
        group_idx +=1
        # add to means list
        zscore_all_dfs.append(df)  
        # add each trial column to a separate list
        for j in range(len(original_cols)):
            trial_df = df[original_cols[j]]
            try:
                zscores_by_trials[j].append(trial_df)
            except: # handle index out of range exception
                # add empty list for trials
               zscores_by_trials.append([]) 
               zscores_by_trials[j].append(trial_df)
    # add subjects and their group names 
    group_df.extend(group_trials_dfs)
    # create a single dataframe from that
    group_header_df4csv = pd.concat(group_df,axis=1)
    # create headers with group names for standard header
    for name in group_names:
        name2write = name+","
        subjects_group_names+=name2write

    # add subjects to settings
    # exclude last coma
    settings_dict[0]["subject"] = subjects_end_file[:-1]
    settings_dict[0]["subject_group_name"] = subjects_group_names[:-1]

    # calculate means
    # join all trials from all subjects to get means ans errors
    df_all_trials = pd.concat(zscore_all_dfs,axis=1)
    transposed_df = df_all_trials.transpose()
    means_df = transposed_df.mean(axis=0)
    means_df = pd.DataFrame({"Mean":means_df.tolist()})
    standard_devs_df = transposed_df.std(axis=0)
    stds = standard_devs_df.tolist()
    se = [stds[i]/np.sqrt(df_all_trials.shape[1]) for i in range(len(stds))]
    se_df = pd.DataFrame({"Error":se})
    
    # calculate means by trial for the heatmap
    zscores_by_trials_list = []
    for i in range(len(zscores_by_trials)):
        single_trial_df = pd.concat(zscores_by_trials[i],axis=1)
        transposed = single_trial_df.transpose()
        means = transposed.mean(axis=0)
        zscores_by_trials_list.append(means.tolist())
    
    # create main dataframe
    raw_df= pd.concat([ts_df[0],df_all_trials,means_df,se_df],axis=1)
    # get data to plot
    ts = ts_df[0]['Time (sec)'].tolist()
    mean_zscore = means_df['Mean'].tolist()

    # plot
    # clear previous figure
    canvas.fig.clf()
    ax = canvas.fig.add_subplot(211)
    ax2 = canvas.fig.add_subplot(212)
    
    # Heat Map based on z score of control fit subtracted signal
    cs = ax.imshow(zscores_by_trials_list, cmap=plt.cm.bwr, interpolation='none', aspect="auto",
                    extent=[-perievent_options_dict['sec_before'], perievent_options_dict['sec_after'], 
                            len(zscores_by_trials_list),0])
    canvas.fig.colorbar(cs, ax=ax,pad=0.01, fraction=0.02)
    # plot the z-score trace for the signal with std error bands
    ax2.plot(ts, mean_zscore, linewidth=2, color=SIGNAL_COLOR_RGB, label='Group Mean')
    for i in range(len(zscores_means_by_subject)):
        if i == 0:
            ax2.plot(ts, zscores_means_by_subject[i], linewidth=0.5, alpha=0.2, color=SIGNAL_COLOR_RGB, label='Subject Mean')
        else:
            ax2.plot(ts, zscores_means_by_subject[i], linewidth=0.5, alpha=0.2, color=SIGNAL_COLOR_RGB)

    ax2.axvline(x=0, linewidth=2, color='k', alpha=0.3, label=event_name)
  
    my_title = 'Total Subjects: ' + str(len(all_dfs)) + "; Event: " + event_name
    canvas.fig.suptitle(my_title, fontsize=FIGURE_TITLE_FONT_SIZE)
    ax.set_title('z-Score Traces',fontsize=TITLE_FONT_SIZE)
    ax.set_ylabel('Trials',fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_xlabel('Seconds from Event Onset',fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_yticks(np.arange(0.5,len(zscores_by_trials_list), 1))
    ax.set_yticklabels(np.arange(1, len(zscores_by_trials_list)+1, 1))
    
    ax2.set_ylabel('z-Score',fontsize=AXIS_LABEL_FONTSIZE)
    ax2.set_xlabel('Seconds',fontsize=AXIS_LABEL_FONTSIZE)
    # hide top and right border
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.legend(loc=MY_LEGEND_LOC, bbox_to_anchor=MY_LEGENG_POS)
    ax2.set_xlim([-perievent_options_dict['sec_before'], perievent_options_dict['sec_after']])
    
    # set the spacing between subplots 
    canvas.fig.subplots_adjust(left=MY_LEFT, 
                    bottom=MY_BOTTOM,  
                    right=MY_RIGHT,  
                    top=MY_TOP,  
                    wspace=MY_WSPACE,  
                    hspace=MY_HSPACE) 
    # save csv with data
    file_name = file_beginning +"_perievent_zscore.csv"
    dump_file_path = os.path.join(dump_path,file_name)
    settings_df = get_settings_df(settings_dict)
    settings_df.to_csv(dump_file_path,header=False)            
    with open(dump_file_path,'a') as f:
        f.write("\n")
    group_header_df4csv.to_csv(dump_file_path, header=False,index=False,mode='a')
    raw_df.to_csv(dump_file_path, index=False,mode='a')
      
    # save plot
    plot_file_name = file_beginning +"_perievent_zscore.png"
    dump_plot_file_path = os.path.join(dump_path,plot_file_name)
    canvas.fig.savefig(dump_plot_file_path)
    # also as svg
    # plot_file_name = file_beginning + subjects_end_file+"_perievent_zscore.svg"
    plot_file_name = file_beginning +"_perievent_zscore.svg"
    dump_plot_file_path = os.path.join(dump_path,plot_file_name)
    canvas.fig.savefig(dump_plot_file_path, format='svg', dpi=DPI4SVG)
    
def get_batch_perievent_auc(canvas,my_all_dfs,group_names,perievent_options_dict,settings_dict,export_loc_data):
    event_name = perievent_options_dict["event_name"] if len(perievent_options_dict["event_name"]) > 0 else perievent_options_dict["event"]
    settings_dict[0]["baseline_from_sec"] =perievent_options_dict["baseline_from"]
    settings_dict[0]["baseline_to_sec"] =perievent_options_dict["baseline_to"]
    settings_dict[0]["auc_pre_from"] =perievent_options_dict["auc_pre_from"]
    settings_dict[0]["auc_pre_to"] =perievent_options_dict["auc_pre_to"]
    settings_dict[0]["auc_post_from"] =perievent_options_dict["auc_post_from"]
    settings_dict[0]["auc_post_to"] =perievent_options_dict["auc_post_to"]
    dump_path,file_beginning = export_loc_data 
    all_dfs = copy.deepcopy(my_all_dfs)
    # contains single dataframe with time
    ts_df = []
    # contains a list of all subject's dataframe(with all normalized trials)
    zscore_all_dfs = []
    # zscores_means_by_subject = []
    zscores_by_trial = [[]]
    all_new_column_names = []
    # remember all subjects to include in file name
    subjects_end_file = ""
    subjects_group_names = ""
    for el in all_dfs:
        subject,df = el
        subject_string = subject + ","
        subjects_end_file+=subject_string
        del df['Mean_zscore']
        del df['Error']
        # if first dataframe keep time
        if len(ts_df) == 0:
            ts = df.iloc[:,:1]
            ts_df.append(ts)
        # remove time
        del df['Time (sec)']
        original_cols = list(df.columns)
        new_columns = []
        trial_idx = 0
        for col in original_cols:
            new_name = subject+"_"+col
            new_columns.append(new_name)
            all_new_column_names.append(new_name)
            # add trial data
            if trial_idx < len(zscores_by_trial):
                trial_zscore = df[col]
                zscores_by_trial[trial_idx].append(trial_zscore)
            else:
                zscores_by_trial.append([])
                trial_zscore = df[col]
                zscores_by_trial[trial_idx].append(trial_zscore)
            trial_idx +=1
            
        # replace old column names with new that contain subject name
        df.columns = new_columns
        # add to means list
        zscore_all_dfs.append(df)    
    
    for name in group_names:
        name2write = name+","
        subjects_group_names+=name2write

    # add subjects to settings
    # exclude last coma
    settings_dict[0]["subject"] = subjects_end_file[:-1]
    settings_dict[0]["subject_group_name"] = subjects_group_names[:-1]
    # calculate means
    # join all trials from all subjects to get means and errors
    df_all_trials = pd.concat(zscore_all_dfs,axis=1)
    # print(df_all_trials)
    transposed_df = df_all_trials.transpose()
    # print(transposed_df)
    means_df = transposed_df.mean(axis=0)
    means_df = pd.DataFrame({"Mean":means_df.tolist()})

    # create data frames by trial
    by_trials_dfs = []
    for i in range(len(zscores_by_trial)):
        single_trial_df = pd.concat(zscores_by_trial[i],axis=1)
        # transpose to get means
        T_single_trial_df = single_trial_df.transpose()
        trial_means = T_single_trial_df.mean(axis=0)
        single_trial_means = pd.DataFrame({i+1:trial_means.tolist()})
        by_trials_dfs.append(single_trial_means)
    # join means by trial into single dataframe
    by_trials_df = pd.concat(by_trials_dfs,axis=1)
    # print(by_trials_df.head(20))   

    mean_zscore = np.asarray(means_df['Mean'].tolist())
    ts = np.asarray(ts_df[0]['Time (sec)'].tolist())
    
    # Quantify changes as an area under the curve
    pre_ind = np.where((np.array(ts)<perievent_options_dict["auc_pre_to"]) & (np.array(ts)>perievent_options_dict["auc_pre_from"]))
   
    AUC_pre = auc(ts[pre_ind], mean_zscore[pre_ind])
    post_ind = np.where((np.array(ts)>perievent_options_dict["auc_post_from"]) & (np.array(ts)<perievent_options_dict["auc_post_to"]))
    AUC_post= auc(ts[post_ind], mean_zscore[post_ind])
    AUC = [AUC_pre, AUC_post]   


    aucs_pre_by_trial =[]
    aucs_post_by_trial = []
    zscore_all = []
    # get each trials zscores
    for col in by_trials_df.columns:
        zscore_all.append(np.asarray(by_trials_df[col].tolist()))
    for i in range(len(zscore_all)):
        pre = auc(ts[pre_ind], zscore_all[i][pre_ind])
        aucs_pre_by_trial.append(pre)
        post = auc(ts[post_ind], zscore_all[i][post_ind])
        aucs_post_by_trial.append(post)

    AUC_pre_err = np.std(np.asarray(aucs_pre_by_trial)/np.sqrt(len(aucs_pre_by_trial)))
    AUC_post_err = np.std(np.asarray(aucs_post_by_trial)/np.sqrt(len(aucs_post_by_trial)))
    # print("errors",AUC_pre_err,AUC_post_err)
        

    # bin auc data in one second bins
    # create indexes
    # print(df_all_trials.head(20))
    col_names = df_all_trials.columns.tolist()
    aucs = [[] for i in range(len(col_names))]
    start = -perievent_options_dict["sec_before"]
    end = start + 1
    while end < perievent_options_dict["sec_after"]:
        ind = np.where((np.array(ts)<end) & (np.array(ts)>start))
        for i in range(len(col_names)):
            trial_zscore = np.asarray(df_all_trials[col_names[i]].tolist())
            # if i == 0:
            #     print(trial_zscore[ind])
            my_auc = auc(ts[ind], trial_zscore[ind])
            aucs[i].append(my_auc)
        start = end
        end = start+1
        
    # plot
    # clear previous figure
    canvas.fig.clf()
    ax = canvas.fig.add_subplot(111)
   
    ax.bar(np.arange(len(AUC)), AUC, color=[.8, .8, .8], align='center', alpha=0.5)
    ax.axhline(y=0,linewidth=0.5,color='k')
    ax.errorbar(np.arange(len(AUC)), AUC, yerr=[AUC_pre_err,AUC_post_err], fmt='none',color='k',capsize=4)
    x1, x2 = 0, 1 # columns indices for labels
    # y, h, col = max(AUC) + 2, 2, 'k'
    y, h, col = max(AUC) + max([AUC_pre_err,AUC_post_err]), 2, 'k'
        
    my_title = 'Total Subjects: ' + str(len(all_dfs))+ "; Event: "+event_name
    canvas.fig.suptitle(my_title, fontsize=FIGURE_TITLE_FONT_SIZE)
    # hide top and right border
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    if min(AUC) < 0:
        ax.set_ylim(min(AUC)-2*h, y+2*h)
    else:
        ax.set_ylim(0-2, y+2*h)
    ax.set_ylabel('AUC',fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_title('Pre vs Post',fontsize=TITLE_FONT_SIZE)
    ax.set_xticks(np.arange(-1, len(AUC)+1))
    ax.set_xticklabels(['', 'PRE', 'POST', ''])
    ax.set_yticks([round(min(AUC), 2),round(max(AUC), 2)])
    
    # save csv with data
    file_name = file_beginning+"_perievent_auc.csv"
    dump_file_path = os.path.join(dump_path,file_name)
    settings_df = get_settings_df(settings_dict)
    settings_df.to_csv(dump_file_path,header=False)            
    with open(dump_file_path,'a') as f:
        f.write("\n")

    # save each trial's auc instead
    raw_df= pd.DataFrame({"Trial":[i+1 for i in range(len(aucs_pre_by_trial))],
                                 "Pre_AUC" : aucs_pre_by_trial,
                                 "Post_AUC":aucs_post_by_trial})
    raw_df.to_csv(dump_file_path, index=False,mode='a')
    
    # export also auc binned in one second bins
    sec_auc_dfs = []
    for i in range(len(aucs)):
        d = pd.DataFrame({col_names[i]:aucs[i]})
        sec_auc_dfs.append(d)
    sec_auc_df = pd.concat(sec_auc_dfs,axis = 1)
    auc_file_name = file_beginning+"_auc_one_second_bins.csv"
    dump_file_path = os.path.join(dump_path,auc_file_name)
    settings_df.to_csv(dump_file_path,header=False)            
    with open(dump_file_path,'a') as f:
        f.write("\n")
    sec_auc_df.to_csv(dump_file_path, index=False,mode='a') 
      
    # save plot
    plot_file_name = file_beginning +"_perievent_auc.png"
    dump_plot_file_path = os.path.join(dump_path,plot_file_name)
    canvas.fig.savefig(dump_plot_file_path)
    # also as svg
    # plot_file_name = file_beginning + subjects_end_file+"_perievent_auc.svg"
    plot_file_name = file_beginning +"_perievent_auc.svg"
    dump_plot_file_path = os.path.join(dump_path,plot_file_name)
    canvas.fig.savefig(dump_plot_file_path, format='svg', dpi=DPI4SVG)
    
def show_polynomial_fitting(canvas, settings_dict,downsampled,signal_name,control_name,subject_name,save_plot,dump_path):
    normalization = settings_dict["normalization"]
    smooth = settings_dict["filter"]
    smooth_window = settings_dict["filter_window"]
    
    # change lists to numpy array for calculations
    ts_arr = np.asarray(downsampled["ts"])
    signal_arr =np.asarray(downsampled["signal"])
    control_arr = np.asarray(downsampled["control"])
    # reset time to start from zero
    total_seconds = ts_arr[-1]-ts_arr[0]
    # start time from zero
    ts_reset = [i*total_seconds/len(ts_arr) for i in range(len(ts_arr))]

    if smooth == True:
        print("Start smoothing",smooth_window)
        a = 1
        b = np.divide(np.ones((smooth_window,)), smooth_window)
        control_arr = filtfilt(b, a, control_arr)
        signal_arr = filtfilt(b, a, signal_arr)
        print("Done smoothing")

    # in order to suggest if user should normalize using modified method
    # check if signals in both channels do not decrease equally
    # filter out signal values that are below or above 2 standard deviations from the signal mean 
    # fit time axis to the 465nm stream 
    mean_signal = np.mean(signal_arr)
    stdev_signal = np.std(signal_arr)
    mean_control = np.mean(control_arr)
    stdev_control = np.std(control_arr)
    indexes_signal = np.where((signal_arr<mean_signal+2*stdev_signal) & (signal_arr>mean_signal-2*stdev_signal))
    selected_signal = signal_arr[indexes_signal]
    selected_ts_signal_arr = ts_arr[indexes_signal]
    indexes_control = np.where((control_arr<mean_control+2*stdev_control) & (control_arr>mean_control-2*stdev_control))
    selected_control = control_arr[indexes_control]
    selected_ts_control_arr = ts_arr[indexes_control]
    # fit time axis to the 465nm stream  
    bls_Ca = np.polynomial.polynomial.Polynomial.fit(selected_ts_signal_arr,selected_signal,1)
    # fit time axis the 405nm stream
    bls_ref = np.polynomial.polynomial.Polynomial.fit(selected_ts_control_arr,selected_control,1)
    #######################
    # Before removing +/-2xstdev
    # bls_Ca = np.polynomial.polynomial.Polynomial.fit(ts_reset,signal_arr,1)
    # # fit time axis the 405nm stream
    # bls_ref = np.polynomial.polynomial.Polynomial.fit(ts_reset,control_arr,1)
    ##################################################################################
    # the below returns first: slope, second: intercept
    print("bls_Ca",bls_Ca.convert().coef[::-1])
    # the below returns first: slope, second: intercept
    print("bls_ref",bls_ref.convert().coef[::-1])
    # put those values in a dictionary
    slope_intercept_dict = {"signal_slope_intercept":bls_Ca.convert().coef[::-1],
                            "control_slope_intercept":bls_ref.convert().coef[::-1]
    }

    if normalization == 'Standard Polynomial Fitting':      
        # Before removing +/-2xstdev             
        # # https://stackoverflow.com/questions/45338872/matlab-polyval-function-with-three-outputs-equivalent-in-python-numpy
        # mu = np.mean(control_arr)
        # std = np.std(control_arr, ddof=0)
        # # Call np.polyfit(), using the shifted and scaled version of control_arr
        # cscaled = np.polynomial.polynomial.Polynomial.fit((control_arr - mu)/std, signal_arr, 1)
        # # Create a poly1d object that can be called
        # # https://numpy.org/doc/stable/reference/routines.polynomials.html
        # pscaled = Polynomial(cscaled.convert().coef)
        # # Inputs to pscaled must be shifted and scaled using mu and std
        # F0 = pscaled((control_arr - mu)/std)    
        #############################################################################################
        # filter out signal values that are below or above 2 standard deviations from the signal mean 
        mean_signal = np.mean(signal_arr)
        stdev_signal = np.std(signal_arr)
        indexes = np.where((signal_arr<mean_signal+2*stdev_signal) & (signal_arr>mean_signal-2*stdev_signal))
        selected_signal = signal_arr[indexes]
        selected_control = control_arr[indexes]
        mu = np.mean(selected_control)
        std = np.std(selected_control)
        cscaled = np.polynomial.polynomial.Polynomial.fit((selected_control - mu)/std, selected_signal, 1)
        pscaled = Polynomial(cscaled.convert().coef)
        F0 = pscaled((control_arr - mu)/std)
        #################################################################################
        # plot
        # clear previous figure
        canvas.fig.clf()
        ax = canvas.fig.add_subplot(111)
        
        ax.plot(ts_reset, signal_arr, linewidth=1, color=SIGNAL_COLOR_RGB, label='F'+signal_name)
        ax.plot(ts_reset, F0, linewidth=1, color='k', label='F0')
        
        my_title = "Check polynomial fitting: " + normalization
        ax.set_title(my_title,fontsize=FIGURE_TITLE_FONT_SIZE)
#        canvas.fig.suptitle(my_title, fontsize=FIGURE_TITLE_FONT_SIZE)
        ax.set_ylabel('mV',fontsize=AXIS_LABEL_FONTSIZE)
        ax.set_xlabel('Time (sec)', fontsize=AXIS_LABEL_FONTSIZE)
        # show only intiger yticks
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        # hide top and right border
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.legend(loc=MY_LEGEND_LOC, bbox_to_anchor=MY_LEGENG_POS)
    
        canvas.draw()
    
    if normalization == 'Modified Polynomial Fitting':           
        F0Ca = polyval(ts_reset,bls_Ca.convert().coef)
        F0Ref = polyval(ts_reset,bls_ref.convert().coef)
        #######################################################
        # plot
        # clear previous figure
        canvas.fig.clf()
        ax = canvas.fig.add_subplot(211)
        ax2 = canvas.fig.add_subplot(212)
        
        ax.plot(ts_reset, signal_arr, linewidth=1, color=SIGNAL_COLOR_RGB, label='F'+signal_name)
        ax.plot(ts_reset, F0Ca, linewidth=1, color='k', label='F0')
        
        ax2.plot(ts_reset, control_arr, linewidth=1, color=CONTROL_COLOR_RGB, label='F'+control_name)
        ax2.plot(ts_reset, F0Ref, linewidth=1, color='k', label='F0')
        
        my_title = "Check polynomial fitting: " + normalization
#        ax.set_title(my_title,fontsize=14)
        canvas.fig.suptitle(my_title, fontsize=FIGURE_TITLE_FONT_SIZE)
        ax.set_ylabel('mV',fontsize=AXIS_LABEL_FONTSIZE)
        ax.set_xlabel('Time (sec)', fontsize=AXIS_LABEL_FONTSIZE)
        # show only intiger yticks
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        # hide top and right border
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.legend(loc=MY_LEGEND_LOC, bbox_to_anchor=MY_LEGENG_POS)
        ax2.set_ylabel('mV',fontsize=AXIS_LABEL_FONTSIZE)
        ax2.set_xlabel('Time (sec)', fontsize=AXIS_LABEL_FONTSIZE)
        # show only intiger yticks
        ax2.yaxis.set_major_locator(MaxNLocator(integer=True))
        # hide top and right border
        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        ax2.legend(loc=MY_LEGEND_LOC, bbox_to_anchor=MY_LEGENG_POS)

    if save_plot == True:
        my_path = dump_path
        file_name = subject_name+"_poly_fitting.png"
        save_path = os.path.join(my_path,file_name)
        canvas.fig.savefig(save_path)
        # also as svg
        file_name = subject_name+"_poly_fitting.svg"
        dump_plot_file_path = os.path.join(dump_path,file_name)
        canvas.fig.savefig(dump_plot_file_path, format='svg', dpi=DPI4SVG)
    else:
        canvas.draw()
    return slope_intercept_dict

def show_polynomial_fitting_custom_baseline(canvas, settings_dict,downsampled,baseline_dict,signal_name,control_name,subject_name,save_plot,dump_path):
    normalization = settings_dict["normalization"]
    smooth = settings_dict["filter"]
    smooth_window = settings_dict["filter_window"]
    
    # change lists to numpy array for calculations
    ts_arr = np.asarray(downsampled["ts"])
    signal_arr =np.asarray(downsampled["signal"])
    control_arr = np.asarray(downsampled["control"])
    print(f"All trace: first: {ts_arr[0]}, last: {ts_arr[-1]}")
    # change baseline lists to numpy array for calculations
    ts_baseline_arr = np.asarray(baseline_dict["ts"])
    signal_baseline_arr =np.asarray(baseline_dict["signal"])
    control_baseline_arr = np.asarray(baseline_dict["control"])
    print(f"Baseline: first: {ts_baseline_arr[0]}, last: {ts_baseline_arr[-1]}")
    # reset time to start from zero
    total_seconds = ts_arr[-1]-ts_arr[0]
    total_seconds_baseline = ts_baseline_arr[-1]-ts_baseline_arr[0]
    # start time from zero
    ts_reset = [i*total_seconds/len(ts_arr) for i in range(len(ts_arr))]

    if smooth == True:
        print("Start smoothing",smooth_window)
        a = 1
        b = np.divide(np.ones((smooth_window,)), smooth_window)
        control_arr = filtfilt(b, a, control_arr)
        signal_arr = filtfilt(b, a, signal_arr)
        control_baseline_arr = filtfilt(b, a, control_baseline_arr)
        signal_baseline_arr = filtfilt(b, a, signal_baseline_arr)
        print("Done smoothing")

    # in order to suggest if user should normalize using modified method
    # check if signals in both channels do not decrease equally
    # filter out signal values that are below or above 2 standard deviations from the signal mean  
    mean_baseline_signal = np.mean(signal_baseline_arr)
    stdev_baseline_signal = np.std(signal_baseline_arr)
    mean_baseline_control = np.mean(control_baseline_arr)
    stdev_baseline_control = np.std(control_baseline_arr)
    indexes_signal = np.where((signal_baseline_arr<mean_baseline_signal+2*stdev_baseline_signal) & (signal_baseline_arr>mean_baseline_signal-2*stdev_baseline_signal))
    selected_signal = signal_baseline_arr[indexes_signal]
    selected_ts_baseline_signal_arr = ts_baseline_arr[indexes_signal]
    indexes_control = np.where((control_baseline_arr<mean_baseline_control+2*stdev_baseline_control) & (control_baseline_arr>mean_baseline_control-2*stdev_baseline_control))
    selected_control = control_baseline_arr[indexes_control]
    selected_ts_baseline_control_arr = ts_baseline_arr[indexes_control]
    new_signal_mean = np.mean(selected_signal)
    new_control_mean = np.mean(selected_control)
    # fit time axis to the 465nm stream  
    bls_Ca = np.polynomial.polynomial.Polynomial.fit(selected_ts_baseline_signal_arr,selected_signal,1)
    bls_ref = np.polynomial.polynomial.Polynomial.fit(selected_ts_baseline_control_arr,selected_control,1)

    # put those values in a dictionary
    slope_intercept_dict = {"signal_slope_intercept":bls_Ca.convert().coef[::-1],
                                "control_slope_intercept":bls_ref.convert().coef[::-1]
        }
    ##########################################################
    # # debug
    # df_before = pd.DataFrame({"Time (sec)": ts_baseline_arr,
    #                             "Signal": signal_baseline_arr,
    #                             "Control": control_baseline_arr})
    # df_before.to_csv("custom_baseline_before.csv", index=False)
    # df_signal_after = pd.DataFrame({"Sample": selected_ts_baseline_signal_arr,
    #                             "Signal": selected_signal})
    # df_control_after = pd.DataFrame({"Sample": selected_ts_baseline_control_arr,
    #                             "Control": selected_control})
    # df_signal_after.to_csv("custom_baseline_signal_filtered.csv", index=False)
    # df_control_after.to_csv("custom_baseline_control_filtered.csv", index=False)
    #################################################################################################
    if normalization == 'Standard Polynomial Fitting':    
        # filter out signal values that are below or above 2 standard deviations from the signal mean  
        mean_baseline_signal = np.mean(signal_baseline_arr)
        stdev_baseline_signal = np.std(signal_baseline_arr)
        indexes = np.where((signal_baseline_arr<mean_baseline_signal+2*stdev_baseline_signal) & (signal_baseline_arr>mean_baseline_signal-2*stdev_baseline_signal))
        selected_signal = signal_baseline_arr[indexes]
        selected_control = control_baseline_arr[indexes]
        mu = np.mean(selected_control)
        std = np.std(selected_control)
        cscaled = np.polynomial.polynomial.Polynomial.fit((selected_control - mu)/std, selected_signal, 1)
        pscaled = Polynomial(cscaled.convert().coef)
        F0 = pscaled((control_arr - mu)/std)
        
        ##########################################################################################################           
        # # https://stackoverflow.com/questions/45338872/matlab-polyval-function-with-three-outputs-equivalent-in-python-numpy
        # mu = np.mean(control_baseline_arr)
        # std = np.std(control_baseline_arr, ddof=0)
        # # Call np.polyfit(), using the shifted and scaled version of control_arr
        # cscaled = np.polynomial.polynomial.Polynomial.fit((control_baseline_arr - mu)/std, signal_baseline_arr, 1)
        # # Create a poly1d object that can be called
        # # https://numpy.org/doc/stable/reference/routines.polynomials.html
        # pscaled = Polynomial(cscaled.convert().coef)
        # # Inputs to pscaled must be shifted and scaled using mu and std
        # F0 = pscaled((control_arr - mu)/std)  
        ###########################################################################################################  
        # plot
        # clear previous figure
        canvas.fig.clf()
        ax = canvas.fig.add_subplot(111)
        
        ax.plot(ts_reset, signal_arr, linewidth=1, color=SIGNAL_COLOR_RGB, label='F'+signal_name)
        ax.plot(ts_reset, F0, linewidth=1, color='k', label='F0')
        
        my_title = "Check polynomial fitting: " + normalization
        ax.set_title(my_title,fontsize=FIGURE_TITLE_FONT_SIZE)
#        canvas.fig.suptitle(my_title, fontsize=FIGURE_TITLE_FONT_SIZE)
        ax.set_ylabel('mV',fontsize=AXIS_LABEL_FONTSIZE)
        ax.set_xlabel('Time (sec)', fontsize=AXIS_LABEL_FONTSIZE)
        # show only intiger yticks
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        # hide top and right border
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.legend(loc=MY_LEGEND_LOC, bbox_to_anchor=MY_LEGENG_POS)
    
        canvas.draw()
    
    if normalization == 'Modified Polynomial Fitting': 
        F0Ca = polyval(ts_arr,bls_Ca.convert().coef)
        F0Ref = polyval(ts_arr,bls_ref.convert().coef)  
        # print(f"signal mean, stdv, len: {mean_baseline_signal}, {stdev_baseline_signal}, {len(selected_signal)}") 
        # print(f"control mean, stdv, len: {mean_baseline_control}, {stdev_baseline_control}, {len(selected_control)}") 
        # print(f"signal poly: {bls_Ca.convert().coef[::-1]}\ncontrol poly: {bls_ref.convert().coef[::-1]}")
###############################################################################
        # F0Ca = polyval(ts_reset,bls_Ca.convert().coef)
        # F0Ref = polyval(ts_reset,bls_ref.convert().coef)
        # plot
        # clear previous figure
        canvas.fig.clf()
        ax = canvas.fig.add_subplot(211)
        ax2 = canvas.fig.add_subplot(212)
        
        ax.plot(ts_reset, signal_arr, linewidth=1, color=SIGNAL_COLOR_RGB, label='F'+signal_name)
        ax.plot(ts_reset, F0Ca, linewidth=1, color='k', label='F0')
        
        ax2.plot(ts_reset, control_arr, linewidth=1, color=CONTROL_COLOR_RGB, label='F'+control_name)
        ax2.plot(ts_reset, F0Ref, linewidth=1, color='k', label='F0')
        
        my_title = "Check polynomial fitting: " + normalization
#        ax.set_title(my_title,fontsize=14)
        canvas.fig.suptitle(my_title, fontsize=FIGURE_TITLE_FONT_SIZE)
        ax.set_ylabel('mV',fontsize=AXIS_LABEL_FONTSIZE)
        ax.set_xlabel('Time (sec)', fontsize=AXIS_LABEL_FONTSIZE)
        # show only intiger yticks
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        # hide top and right border
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.legend(loc=MY_LEGEND_LOC, bbox_to_anchor=MY_LEGENG_POS)
        ax2.set_ylabel('mV',fontsize=AXIS_LABEL_FONTSIZE)
        ax2.set_xlabel('Time (sec)', fontsize=AXIS_LABEL_FONTSIZE)
        # show only intiger yticks
        ax2.yaxis.set_major_locator(MaxNLocator(integer=True))
        # hide top and right border
        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        ax2.legend(loc=MY_LEGEND_LOC, bbox_to_anchor=MY_LEGENG_POS)

    if save_plot == True:
        my_path = dump_path
        file_name = subject_name+"_poly_fitting.png"
        save_path = os.path.join(my_path,file_name)
        canvas.fig.savefig(save_path)
        # also as svg
        file_name = subject_name+"_poly_fitting.svg"
        dump_plot_file_path = os.path.join(dump_path,file_name)
        canvas.fig.savefig(dump_plot_file_path, format='svg', dpi=DPI4SVG)
    else:
        canvas.draw()
    return slope_intercept_dict
    
   
    
    