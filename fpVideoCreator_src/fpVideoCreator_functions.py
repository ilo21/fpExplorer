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
Many analysis ideas come from TDT
https://www.tdt.com/support/python-sdk/offline-analysis-examples/fiber-photometry-epoch-averaging-example/

Many of the analysis ideas are inspired on code provided by the following sources:
Tucker-Davis Technologies ( https://www.tdt.com/support/python-sdk/offline-analysis-examples/fiber-photometry-epoch-averaging-example/)
Dr. David Barker (pMAT; Bruno C.A. et al. 2021, Pharm BioChem Behav 201, https://doi.org/10.1016/j.pbb.2020.173093)
Dr. Patrik Mulholland (Braunsheidel K.M. et al. 2019, J Neurosci 39(46), https://doi.org/10.1523/JNEUROSCI.1674-19.2019)


"""
# https://github.com/LABSN/tdtpy
# https://www.tdt.com/support/python-sdk/offline-analysis-examples/introduction-to-python-tdt-package/
import tdt
import tdt
from tdt.TDTfilter import combine_time
from tdt.TDTfilter import get_valid_ind
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
from matplotlib import animation
from matplotlib.animation import FuncAnimation, FFMpegWriter
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.signal import find_peaks, filtfilt
from scipy.interpolate import interp1d
import os
import copy
import cv2
from tqdm import tqdm
import time
import warnings



#################################
# PLOT GENERAL SETTINGS
##############################

DPI4SVG = 1200

DPI4ANIMATION = 200
# the bigger the number the slower the animation
DELAY = 50 
# parameters for joined visualization
DESIRED_VIDEO_HEIGHT = 600
DESIRED_VIDEO_WIDTH = 1200
DESIRED_ANIMATION_HEIGHT = 300
DESIRED_ANIMATION_WIDTH = DESIRED_VIDEO_WIDTH

VIDEO_EXTENTIONS = ["avi"]
CAM_EVT_NAME = ["Cam1","Cam2"]
################################
# return raw data structure
def get_raw_data(path):
    '''Returns a raw data extracted by tdt from all recording files
    '''
#    raw_data = tdt.read_block("C://Users//ilosz01//OneDrive - Linköpings universitet//FPproject//Data//1//fear_test")
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

def get_cam_events(raw_data):
    '''
    Returns a list of selected cam events from file
    '''
    # tdt StructTypes are dictionaries
    all_epocs_list = [key for key in raw_data.epocs.keys()]
    my_event_list = []
    # iterate over all events that start with Ptr
    for evt in all_epocs_list:
        if evt in CAM_EVT_NAME:
            my_event_list.append(evt)
    return my_event_list

# works for just a single video file in folder
def is_there_a_video(path,cam_event):
    there_is = False
    file_path = None
    files = os.listdir(path)
    for f in files:
        file_splited = f.split('.')
        if len(file_splited) > 1:
            ext = file_splited[-1]
            file_name = file_splited[0]
            file_name_splited = file_name.split("_")
            if len(file_name_splited) > 1:
                file_name_ending = file_name_splited[-1]
                if ext in VIDEO_EXTENTIONS and file_name_ending == cam_event:
                    file_path = os.path.join(path,f)
                    there_is = True
    print(f"Video path: {file_path}")
    return there_is, file_path

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
            print("Problem getting on off event data by tdt.epoc_filter")
            print(event_name_split[0]) # some names end with _ and that gets replaced in data tank
            if event_name_split[0][-1] == "_":
                adjusted_name = event_name_split[0][:-1]+"/"
                print(adjusted_name)
                evt_on_off= tdt.epoc_filter(raw_data, adjusted_name, values=[int(event_name_split[1])])
                on_off = evt_on_off.time_ranges
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
    # print(f"\nEvent {event} onsets: {on_off[0]}")
    return on_off

def get_cam_event_on_off(raw_data,cam_evt):
    ''' Returns onset and offset times of the cam1 event
        First element of this array are onset times, 
        Second el are offset times
    '''
    on_off = [[],[]]
    data = raw_data.epocs[cam_evt]
    onsets = []
    try:
        onsets = data.onset
        offset = data.offset
    except:
        try:
            onsets = data.notes.ts
            offset = []
        except:
            onsets = []
            offset = []
    on_off[0] = onsets
    on_off[1] = offset[:-2] if len(offset)>0 else offset
    print(f"\nTotal frames in cam event: {len(on_off[0])}")
    print(f"First frame (sec): {on_off[0][0]} Last frame (sec): {on_off[0][-1]}")
    minutes = int(int(on_off[0][-1])/60)
    seconds = int(on_off[0][-1])%60
    print('duration (M:S) = ' + str(minutes) + ':' + str(seconds))
    if len(on_off[0]) == 0:
        print("Could not get onsets any of known ways")
    return on_off



def get_frequency(raw_data,chnl_name):
    return raw_data.streams[chnl_name].fs



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

def filter_data_around_event_directly(raw_data,event_name,sec_before,sec_after,settings_dict,signal_name,control_name):
    event_name_split = event_name.split(" ")
    before = -sec_before
    till = sec_before+sec_after+0.1
    trange = [before,till]
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

# plot raw but downsampled perievents
def get_normalized_trials(subject,modified_data,current_trials,event,sec_before,sec_after,settings_dict,signal_name,control_name):
    settings_dict[0]["subject"] = subject
    show_norm_as = settings_dict[0]["show_norm_as"]
    raw_df = None
    event_name = event
    GCaMP_perievent_data = modified_data.streams[signal_name].filtered_downsampled
    control_perievent_data = modified_data.streams[control_name].filtered_downsampled
    if len(current_trials) > 0:
        # only those trials that are selected
        # print(current_trials)
        GCaMP_perievent_data = [GCaMP_perievent_data[trial-1] for trial in current_trials]
        control_perievent_data = [control_perievent_data[trial-1] for trial in current_trials]
    
    N = settings_dict[0]["downsample"]
    # !!! Create times in sec first to get last time value not total samples!!!
    t1From0 = np.linspace(1/N,(sec_before+sec_after),(sec_before+sec_after)*N)
    # print(t1From0)
    ts1 = -sec_before + t1From0
    # print(ts1)
    ts2 = ts1

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
            
            # https://stackoverflow.com/questions/45338872/matlab-polyval-function-with-three-outputs-equivalent-in-python-numpy
            mu = np.mean(x)
            std = np.std(x, ddof=0)
            # Call np.polyfit(), using the shifted and scaled version of control_arr
            cscaled = np.polynomial.polynomial.Polynomial.fit((x - mu)/std, y, 1)
            # Create a poly1d object that can be called
            # https://numpy.org/doc/stable/reference/routines.polynomials.html
            pscaled = Polynomial(cscaled.convert().coef)
            # Inputs to pscaled must be shifted and scaled using mu and std
            F0 = pscaled((x - mu)/std)
        #    print("F0?",F0[:20])
            dffnorm = (y - F0)/F0 * 100
            # find all values of the normalized DF/F that are negative so you can next shift up the curve 
            # to make 0 the mean value for DF/F
            negative = dffnorm[dffnorm<0]
            dff=dffnorm-np.mean(negative)

            if show_norm_as == "Z-Score":
                median_all = np.median(dff)
                mad = stats.median_absolute_deviation(dff)
                dff = (dff - median_all)/mad

            y_dff_all.append(dff)
            
    elif settings_dict[0]["normalization"] == 'Modified Polynomial Fitting':
        for x, y in zip(control_perievent_data, GCaMP_perievent_data):
            x = np.array(x)
            y = np.array(y)
            bls_signal = np.polynomial.polynomial.Polynomial.fit(ts1, y, 1)
            F0_signal = polyval(ts1,bls_signal.convert().coef)
            dFF_signal = (y - F0_signal)/F0_signal *100
            
            bls_control = np.polynomial.polynomial.Polynomial.fit(ts2,x,1)
            F0_control = polyval(ts2,bls_control.convert().coef)
            dFF_control = (x - F0_control)/F0_control *100
            dFFnorm = dFF_signal - dFF_control
            # find all values of the normalized DF/F that are negative so you can next shift up the curve 
            # to make 0 the mean value for DF/F
            negative = dFFnorm[dFFnorm<0]
            dff = dFFnorm-np.mean(negative)

            if show_norm_as == "Z-Score":
                median_all = np.median(dff)
                mad = stats.median_absolute_deviation(dff)
                dff = (dff - median_all)/mad

            y_dff_all.append(dff)
    
    
    # create a list of dataframes to join later
    dfs = []
    time_df = pd.DataFrame({"Time (sec)":ts1})
    dfs.append(time_df)
    if show_norm_as == 'Z-Score':
        norm_type = " (z-score)"
    else:
        norm_type = " (%df/F)"
    # for i in range(total_plots):
    for i in range(len(current_trials)):
        # header_signal = "Trial"+str(i+1)+norm_type
        header_signal = "Trial"+str(current_trials[i])+norm_type
        signal_data = y_dff_all[i]
        df = pd.DataFrame({header_signal:signal_data})
        dfs.append(df)
    raw_df= pd.concat(dfs,axis=1)

    return raw_df


def find_closest_time_pair_indexes(time_stamp,found_index_list=[],times_to_compare=[],time_margin=0.05):
    '''
        time_stamp - time of the event
        found_index_list - what was previously detected but was not a pair of before and after indexes
        times_to_compare - list of frame times (from cam event)
    '''
    if len(found_index_list) > 0 and len(times_to_compare) > 0:
        before = []
        before_idx = []
        after = []
        after_idx = []
        for idx in found_index_list:
            if times_to_compare[idx] < time_stamp:
                before.append(times_to_compare[idx])
                before_idx.append(idx)
            else:
                after.append(times_to_compare[idx])
                after_idx.append(idx)
        try:
            closest_before = max(before)
            closest_before_idx = before_idx[before.index(closest_before)]
            closest_after = min(after)
            closest_after_idx = after_idx[after.index(closest_after)]
            print(f"New Before ts: {closest_before}, new After ts: {closest_after}")
            print(f"Before idx were: {found_index_list}, new idx are: {[closest_before_idx,closest_after_idx]}")
            return [closest_before_idx,closest_after_idx]
        except:
            return []
    elif len(found_index_list) == 0:
        print(f"New margin: {margin}")
        margin = time_margin
        before_and_after_frame_idx = []
        for i in range(len(times_to_compare)):
            if times_to_compare[i] > time_stamp-margin and times_to_compare[i] < time_stamp+margin:
                print(f"Frame ts: {times_to_compare[i]}")
                before_and_after_frame_idx.append(i)
        return before_and_after_frame_idx
    
def shift_indexes_left(input_list,shift_no):  # make sure the shift_no is positive           
    shifted_list = []
    for i in range(len(input_list)):
        el = input_list[i]
        trial_no = el[0]
        before,after = el[1]
        before_updated = before - shift_no
        after_updated = after - shift_no
        if before_updated > 0 and after_updated > 0:
            shifted_list.append((trial_no,[before_updated,after_updated]))
        else:
            print(f"Could not shift trial: {trial_no} data")
            shifted_list.append(el)
    return shifted_list

def shift_indexes_right(input_list,shift_no):             
    shifted_list = []
    for i in range(len(input_list)):
        el = input_list[i]
        trial_no = el[0]
        before,after = el[1]
        before_updated = before + shift_no
        after_updated = after + shift_no        
        shifted_list.append((trial_no,[before_updated,after_updated]))

    return shifted_list
#################################################################################################
# PLOTS

def create_trial_animations(data,trials,event_name,fps,saving_path,subject,cam_event):
    global step_samples
    all_trials_df = data
    which_trials = trials
    sampling_rate = fps
    if not os.path.exists(saving_path):
        os.mkdir(saving_path)
    headers = all_trials_df.columns.tolist()
    # print(f"Headers {headers}")
    # get time and trial values
    time_list = all_trials_df.iloc[:,0].tolist()
    # list of tuples: column idx, trial
    indexes_with_trials = []
    # change trials to column indexes in the dataframe
    for i in range(len(headers)):
        header_split = headers[i].split(" ")
        if header_split[0].startswith("Trial"):
            try:
                # convert all after "Trial" to int
                trial_col = int(header_split[0][5:])
                for trial in which_trials:
                    if trial_col == trial:
                        indexes_with_trials.append((i,trial))
            except:
                pass
                # print("Could not find trial number in the header")
    # for trial in which_trials:
    print()
    print(f"Total trials to animate: {len(indexes_with_trials)}")
    if len(indexes_with_trials) > 0:
        for i in tqdm(range(len(indexes_with_trials)), total= len(indexes_with_trials), mininterval = 3, desc="Saving animations"):
            # print(f"Trial: {trial}")
            trial_col = indexes_with_trials[i][0]
            trial_data = all_trials_df.iloc[:,trial_col].tolist()
            video_duration_seconds = round(time_list[-1] + abs(time_list[0]))
            # For the display animation the duration is going to be frames  interval  1000 (in seconds)
            # how many frames in total video duration
            delay = int(1000/sampling_rate)
            total_frames = int(video_duration_seconds*1000/delay)
            # how many samples per frame
            total_samples = len(trial_data)
            step_samples = int(total_samples/total_frames)
            single_step = step_samples
            ##############################################################################
            fig = plt.figure(figsize=(18, 4))
            ax = fig.add_subplot(111)
            trial_no, unit = headers[trial_col].split(" ")
            fig.suptitle(trial_no, fontsize=14)
            ax.set_xlabel('Time (sec)', fontsize=14)
            ax.set_ylabel(unit[1:-1], fontsize=14)
            fig.tight_layout()
            ax.plot(time_list, trial_data,color='k')
            ax.axvline(x=0,color='k',alpha=0.5,linewidth=1)
            # function that draws each frame of the animation
            def animate(i):
                global step_samples
                if step_samples <= total_samples:
                    plt.cla()
                    x = time_list[:step_samples]
                    y = trial_data[:step_samples]
                    step_samples +=single_step
                    ax.plot(time_list, trial_data, color='k', linewidth=1)
                    ax.plot(x, y, color = 'red', linewidth=1.2)
                    ax.axvline(x=x[-1],color='red',alpha=0.5,linewidth=1)
                    # ax.plot(x[-1], y[-1],marker='D', color = 'red', markersize=8)
                    ax.axvline(x=0,color='k',alpha=0.5,linewidth=1)
                    ax.set_xlabel('Time (sec)', fontsize=14)
                    ax.set_ylabel(unit[1:-1], fontsize=14)
                

            # function to initialize a figure    
            def init():
                ax.set_xlim([time_list[0],time_list[-1]])
                ax.set_ylim([min(trial_data),max(trial_data)])

            # run the animation
            ani = FuncAnimation(fig, animate, frames=total_frames, interval=delay, repeat=False, blit=False)
            # saving to m4 using ffmpeg writer
            writervideo = FFMpegWriter(fps=sampling_rate)
            # create a file name 
            ani_file = subject+"_"+cam_event+"_"+str(indexes_with_trials[i][1])+"_trial_evt_"+event_name+".mp4"
            path2save = os.path.join(saving_path,ani_file)
            ani.save(path2save,writer=writervideo,dpi=DPI4ANIMATION)
        print("Done creating animations\n")
        return True
    else:
        print("Could not read trials dataframe. No animations created")
        return False

def create_trial_videos(video_path,sec_before,sec_after,before_and_after_frameidx_for_all_trials,event_name,saving_path,path4videos,subject,cam_event):
    # create path for saving video images/frames
    head,tail = os.path.split(video_path)
    video_name = tail.split(".")[0]
    path4frames = os.path.join(saving_path,video_name+"_frames")
    fps = None
    frame_count = None
    if not os.path.exists(path4frames):
    #############################################################
    # Save frames first
        print()
        try: 
            os.mkdir(path4frames)
            time.sleep(3)
            # Convert frame to image and save to file
            # Open video file
            video = cv2.VideoCapture(video_path)

            # number of frames in video
            frame_count = int(video.get(cv2.CAP_PROP_FRAME_COUNT))
            fps = video.get(cv2.CAP_PROP_FPS)      # OpenCV v2.x used "CV_CAP_PROP_FPS"
            duration = frame_count/fps
            print(f"Total video frames: {frame_count}")
            print(f"Total video duration: {duration}")
            print(f"FPS: {fps}")
            minutes = int(duration/60)
            seconds = round(duration%60,2)
            print('duration (M:S) = ' + str(minutes) + ':' + str(seconds))
            for i in tqdm(range(frame_count), total= frame_count, mininterval = 3, desc="Saving video frames"):
                ret, frame = video.read()
                if ret:
                    image_path = os.path.join(path4frames,f"{i}.jpg")
                    saved_file = cv2.imwrite(image_path, frame)
                    if saved_file == False:
                        print("\nimwrite function failed to save jpg frame.\nTry changing the saving path to local or with no special characters.\n")
                        return False
            video.release()
        except:
            print("\nCould not create directory for saving video frames\n")
            video.release()
            return False
    else: # if frame images created get the fps 
        # doublecheck that files are there
        if len(os.listdir(path4frames)) == 0:
            print("\nFolder with frame images is empty\n")
            return False
        # Open video file
        video = cv2.VideoCapture(video_path)
        frame_count = int(video.get(cv2.CAP_PROP_FRAME_COUNT))
        fps = video.get(cv2.CAP_PROP_FPS)      # OpenCV v2.x used "CV_CAP_PROP_FPS"
        duration = frame_count/fps
        print(f"Total video frames: {frame_count}")
        print(f"Total video duration: {duration}")
        print(f"FPS: {fps}")
        minutes = int(duration/60)
        seconds = round(duration%60,2)
        print('duration (M:S) = ' + str(minutes) + ':' + str(seconds))
        video.release()
    ##############################################################
    # save fragmants of that video (one per trial)
    print()
    if fps == None:
        print("\nCould not open video file and check fps\n")
        return False
    # how many frames before and after
    total_frames_before = sec_before*fps
    total_frames_after = sec_after*fps
    for i in tqdm(range(len(before_and_after_frameidx_for_all_trials)), total= len(before_and_after_frameidx_for_all_trials), mininterval = 3, desc="Saving trial videos"):
        trial,frame_before_after = before_and_after_frameidx_for_all_trials[i] # it is a tuple trial, [frame_idx_before, frame_idx_after]
        from_idx = int(frame_before_after[0]-total_frames_before)
        till_idx = int(frame_before_after[1]+total_frames_after-1)
        print(f"\nTrial {trial}; from idx: {from_idx}, till idx: {till_idx}")
        # check if the last index frame exists
        if till_idx > frame_count:
            print("\nThe video is shorter than the requested recording data. Try reducing time window around the event\n")
            return False
        # create a list of paths (jpg files should be subsequent numers)
        all_indexes = np.arange(from_idx,till_idx)
        # print(f"Total files: {len(all_indexes)}")
        images4cv2=[]
        for i in range(len(all_indexes)):
            # frame files should have intigers as names to sort them and fing them easily
            file_name = str(all_indexes[i])+".jpg"
            full_file_path = os.path.join(path4frames,file_name)
            try:
                images4cv2.append(cv2.imread(full_file_path))
            except:
                print("File does not exist")

        height,width,layers=images4cv2[0].shape
        # print(f"height,width {height}, {width}")
        # outfile name = same convention as animation ani_file = str(trial)+"_trial_evt_"+event_name+".mp4"
        out_file_name = subject+"_"+cam_event+"_"+str(trial)+"_trial_evt_"+event_name+".avi"
        path2save = os.path.join(path4videos,out_file_name)
        video_out=cv2.VideoWriter(path2save, cv2.VideoWriter_fourcc(*"MJPG"), fps, (width,height))

        for j in range(len(images4cv2)):
            video_out.write(images4cv2[j])

        video_out.release()
        cv2.destroyAllWindows()
    ##############################################################
    # save videos together with animations (one per trial)
    print()
    size = (DESIRED_VIDEO_WIDTH, DESIRED_VIDEO_HEIGHT+DESIRED_ANIMATION_HEIGHT)
    for i in tqdm(range(len(before_and_after_frameidx_for_all_trials)), total= len(before_and_after_frameidx_for_all_trials), mininterval = 3, desc="Saving combined videos"):
        trial,frame_before_after = before_and_after_frameidx_for_all_trials[i] # it is a tuple trial, [frame_idx_before, frame_idx_after]
        # create path to video file
        video_file_name = subject+"_"+cam_event+"_"+str(trial)+"_trial_evt_"+event_name+".avi"
        path2video = os.path.join(path4videos,video_file_name)
        # create path to animation file
        ani_file_name = subject+"_"+cam_event+"_"+str(trial)+"_trial_evt_"+event_name+".mp4"
        path2ani = os.path.join(path4videos,ani_file_name)
        # create path to save output
        out_file_name = subject+"_"+cam_event+"_"+str(trial)+"_trial_evt_"+event_name+"_combined.avi"
        path2out = os.path.join(path4videos,out_file_name)
        # get video
        cap1 = cv2.VideoCapture(path2video)
        # get animation
        cap2 = cv2.VideoCapture(path2ani)
        fourcc = cv2.VideoWriter_fourcc(*'XVID')
        videowriter = cv2.VideoWriter(path2out,fourcc,fps,size)
        # ctr = 0
        while cap1.isOpened() or cap2.isOpened():
            # try:
                # ctr+=1
            # video
            ret,frame1=cap1.read()
            # animation
            ret,frame2=cap2.read()
            if ret:
                    # # Print the shape of frame1 and frame2 to check if they have valid data
                    # print(str(ctr)+" Frame 1 shape:", frame1.shape if frame1 is not None else None)
                    # print(str(ctr)+" Frame 2 shape:", frame2.shape if frame2 is not None else None)
                    # print()
                frame1 = cv2.resize(frame1, (DESIRED_VIDEO_WIDTH,DESIRED_VIDEO_HEIGHT), interpolation = cv2.INTER_AREA)
                frame2 = cv2.resize(frame2, (DESIRED_ANIMATION_WIDTH,DESIRED_ANIMATION_HEIGHT), interpolation = cv2.INTER_AREA)
                together = cv2.vconcat([frame1, frame2])

                videowriter.write(together)
                # cv2.imshow('view' , together)
            else:
                break
            if cv2.waitKey(1) & 0xFF == ord('q'):
                cv2.destroyAllWindows()
                break
            # except Exception as e:
            #     print("Error:", str(e))
        cap1.release()
        cap2.release()
        videowriter.release()
        cv2.destroyAllWindows()
###########################################################################################

    print("Done")
    return True

