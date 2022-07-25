"""
 Copyright (c) 2022 CSAN_LiU

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
Created on Fri 29  4 15:08:47 2022

@author: ilosz01
"""


from PyQt5 import QtGui 
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.QtCore import Qt, pyqtSignal, pyqtSlot, QTimer
import pyqtgraph as pg
from pyqtgraph.dockarea import *
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib.ticker import MaxNLocator
from scipy.signal import find_peaks, filtfilt
import pandas as pd
import numpy as np
from numpy.polynomial import Polynomial
from numpy.polynomial.polynomial import polyval
from scipy.interpolate import interp1d
from scipy import stats
import math
import os
import warnings
import fpExplorer_functions
import fpExplorer
warnings.filterwarnings("ignore")

#icon https://logomakr.com/9OpQfD
# description how to create icons resources
# https://www.learnpyqt.com/courses/packaging-and-distribution/packaging-pyqt5-pyside2-applications-windows-pyinstaller/
import resources



# how low in downsampling can go
MAX_DOWNSAMPLE_PCT = 50
# how high in downsampling can you go
MIN_DOWNSAMPLE_PCT = 0.5
# default_samples to average 1% of original sampling rate
DEFAULT_DOWNSAMPLE_PCT = 10
# change to hardcoded default of 100 Hz
DEFAULT_HZ = 100
MAX_SMOOTH_WINDOW = 100
DEFAULT_SMOOTH_WINDOW = 10
DEFAULT_EXPORT_FOLDER = "_fpExplorerAnalysis"

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
###################################################
# App icon
ICO = ':/icons/app_icon'
# define stylesheet
STYLESHEET = \
            ".QPushButton{\n" \
            + "padding: 7px;\n" \
            + "margin-left:20px;\n" \
            + "margin-right:20px;\n" \
            + "margin-top:5px;\n" \
            + "margin-bottom:5px;\n" \
            + "font-size: 14px;\n" \
            + "font-weight: bold;\n" \
            + "}" 

#####################################
# CLASS FOR SELECT CSV DATA WINDOW  #
#####################################
        
class SelectCvsDataWindow(QMainWindow):
    # pyqt Signal has to be defined on class level
    # and not in init method !!
    close_select_data_sig = pyqtSignal()
    # create a list with information acquired from window
    # and send it as a signal
    got_valid_user_input_sig = pyqtSignal(list)
    def __init__(self,parent_window,user_input):
        super(SelectCvsDataWindow, self).__init__()
        self.setWindowTitle("Select Custom Data")
        self.setWindowIcon(QtGui.QIcon(ICO))
        self.resize(650,250)
        self.parent_window = parent_window
        self.parent_window.app_closed.connect(self.exit_app)

        # list that will contain user input
        self.user_input_list = user_input
        if len(self.user_input_list) > 0:
            self.user_input_dict = self.user_input_list[0]
        else:
            self.user_input_dict = {}
        self.isValid = False

        self.bold_label_stylesheet = "QLabel {font-weight: bold}"
##########################################
        self.main_widget = QWidget()
        self.setCentralWidget(self.main_widget)
        # Create an outer layout
        self.outer_layout = QVBoxLayout()
        # create layout for file selection
        self.form_layout = QFormLayout()
        self.form_layout.setContentsMargins(10,10,10,10)
        self.experiment_name_label = QLabel("Experiment name")
        self.experiment_name_label.setStyleSheet(self.bold_label_stylesheet)
        self.experiment_name_text = QLineEdit("")
        if "selected_experiment" in self.user_input_dict:
            self.experiment_name_text.setText(self.user_input_dict["selected_experiment"])
        self.form_layout.addRow(self.experiment_name_label,self.experiment_name_text)
        self.subject_name_label = QLabel("Subject name")
        self.subject_name_label.setStyleSheet(self.bold_label_stylesheet)
        self.subject_name_text = QLineEdit("")
        if "subject_name" in self.user_input_dict:
            self.subject_name_text.setText(self.user_input_dict["subject_name"])
        self.form_layout.addRow(self.subject_name_label,self.subject_name_text)
        self.form_layout.addRow(QLabel("\n\nSelect csv file with signal and control:"))
        self.select_data_file_btn = QPushButton("Select")
        self.signal_control_file_path_text = QLineEdit("")
        self.form_layout.addRow(self.select_data_file_btn,self.signal_control_file_path_text)
        self.form_layout.addRow(QLabel("Select csv file with the events:"))
        self.select_events_file_btn = QPushButton("Select")
        self.events_file_path_text = QLineEdit("")
        self.form_layout.addRow(self.select_events_file_btn,self.events_file_path_text)

        self.outer_layout.addLayout(self.form_layout)
        # create layout for read button
        self.read_layout = QHBoxLayout()
        self.read_layout.setContentsMargins(10,10,10,10)
        self.read_layout.setAlignment(Qt.AlignRight)
        self.read_btn = QPushButton('Read Data')
        # use stylesheet
        self.read_btn.setStyleSheet(STYLESHEET)
        self.read_layout.addWidget(self.read_btn)
        # nest the inner layout into the outer layout
        self.outer_layout.addLayout(self.read_layout)
        
        # set the window's main layout
        self.main_widget.setLayout(self.outer_layout)

        self.select_data_file_btn.clicked.connect(self.select_data_btn_clicked)
        self.select_events_file_btn.clicked.connect(self.select_events_file_btn_clicked)
        self.read_btn.clicked.connect(self.read_custom_info)

    def select_data_btn_clicked(self):
        # set text field to the value of selected path
        self.signal_control_file_path_text.setText(QFileDialog.getOpenFileName(self,"Select file with raw data")[0])

    def select_events_file_btn_clicked(self):
        # set text field to the value of selected path
        self.events_file_path_text.setText(QFileDialog.getOpenFileName(self,"Select file with raw data")[0])

    def read_custom_info(self):
        # reset old infor if experiment changed
        if "selected_experiment" in self.user_input_dict:
            if self.experiment_name_text.text() != self.user_input_dict["selected_experiment"]:
                self.user_input_dict = {}
        self.user_input_dict["selected_experiment"] = self.experiment_name_text.text()
        self.user_input_dict["subject_name"] = self.subject_name_text.text()
        if len(self.user_input_dict["subject_name"]) > 0:
            self.isValid = True
        else:
            self.show_info_dialog("Subject name cannot be empty")
        if self.isValid == True:
            data_path = self.signal_control_file_path_text.text()
            if os.path.exists(data_path):
                events_path = self.events_file_path_text.text()
                if len(events_path) > 0:
                    if os.path.exists(events_path):
                        if "subject_paths" not in self.user_input_dict:
                            self.user_input_dict["subject_paths"] = {self.user_input_dict["subject_name"]:(data_path,events_path)}
                        else:
                            # if self.user_input_dict["subject_name"] not in self.user_input_dict["subject_paths"]:
                            self.user_input_dict["subject_paths"][self.user_input_dict["subject_name"]] = (data_path,events_path)
                            # else:
                            #     self.isValid = False
                            #     self.show_info_dialog("Subject already added to preview")
                    else:
                        self.isValid = False
                        self.show_info_dialog("Please enter a valid event file path")
                else:
                    if "subject_paths" not in self.user_input_dict:
                        self.user_input_dict["subject_paths"] = {self.user_input_dict["subject_name"]:(data_path,events_path)}
                    else:
                        # if self.user_input_dict["subject_name"] not in self.user_input_dict["subject_paths"]:
                        self.user_input_dict["subject_paths"][self.user_input_dict["subject_name"]] = (data_path,events_path)
                        # else:
                        #     self.isValid = False
                        #     self.show_info_dialog("Subject already added to preview")
            else:
                self.isValid = False
                self.show_info_dialog("Please enter a valid data file path")
        if self.isValid == True:   # check if user input is still valid
            self.got_valid_user_input_sig.emit([self.user_input_dict])

    # show info popup
    def show_info_dialog(self, text):
        msgBox = QMessageBox()
        msgBox.setWindowIcon(QtGui.QIcon(ICO))
        msgBox.setText(text)
        msgBox.setWindowTitle("Info!")
        msgBox.setStandardButtons(QMessageBox.Ok)
        msgBox.exec()

    def exit_app(self):
        self.close()

    # emit signal when window is closed
    def closeEvent(self, event): 
        try: 
            self.close_select_data_sig.emit()
        except:
            print("Main does not exist")

##############################################
# CLASS FOR PREVIEW EVENT-BASED EXPERIMENTS  #
##############################################
            
class PreviewEventBasedWidget(QWidget):
#    close_wait_sig = pyqtSignal()
    close_export_window_signal = pyqtSignal()
    disable_buttons_signal = pyqtSignal()
    enable_buttons_signal = pyqtSignal()
    batch_analysis_complete_signal = pyqtSignal()
    done_batch_processing_sig = pyqtSignal()
    def __init__(self,parent_window,init_params):
        super(PreviewEventBasedWidget, self).__init__()
        # list of init_params to start previewing
        # first element is a dictionary from select data folder window
        # second element is a dict of paths to all data folders
        self.preview_init_params = init_params
        self.parent_window = parent_window
        # list with dictionary with settings from main app window
        self.settings_dict = self.parent_window.settings_dict
        # self.parent_window.start_batch_processing_sig.connect(self.start_batch_analysis)
              
        # options read from options section
        self.options = {"subject":self.preview_init_params[0]["subject_name"],
                        "subject_group_name":"",
                        "event":"---",
                        "event2":"---",
                        "plot_raw":True,
                        "plot_downsampled":False,
                        "plot_normalized":False}

        # remember subject(key):group name(value) for later in a form dictionary
        self.group_names_dict = {self.options["subject"]:self.options["subject_group_name"]}
        
        # a dictionary with subject:extracted raw data
        self.raw_data_dict = {}
        self.events_dict = {}
        # start a widget by reading the subjects data on the list and adding that to dict
        for key,val in self.preview_init_params[0]["subject_paths"].items():
            self.raw_data_dict[key] = self.read_data_csv(val[0])
            if len(val[1])>0:
                self.events_dict[key] = self.read_data_csv(val[1])
        # get last timestamp in seconds
        self.last_raw_ts = self.raw_data_dict[self.options["subject"]].iloc[:,0].iloc[-1]
        # read first subject's frequency and create suggested downsampled rate
        self.current_fs = self.get_frequency(self.options["subject"])
        # # self.suggested_downsample_samples = int(int(self.current_fs)*DEFAULT_DOWNSAMPLE_PCT/100)
        self.suggested_downsample_samples = DEFAULT_HZ
        self.settings_dict[0]['downsample'] = self.suggested_downsample_samples
        # self.suggested_downsample_rate = round(self.current_fs/self.suggested_downsample_samples)
        # self.settings_dict[0]["entered_downsample"] = self.suggested_downsample_rate
        # # print("Current suggested downsample rate:",self.settings_dict[0]['downsample'])
        # keep trimmed data separately with latest trimming settings
        # first element of that list is beginning and end seconds to trim
        # second element is trimmed data dict ts:timestamps, signal:data,control:data
        self.trimmed_raw_data_dict = {}
        
        # create a list of available events for current subject
        self.events_from_current_subject = []
        if self.options["subject"] in self.events_dict:
            self.events_from_current_subject = self.get_events(self.events_dict[self.options["subject"]])
            
        # store downsampled data key (subject): dict("ts","signal","control")
        self.downsampled_dict = {}
        # store normalized data key (subject): dict("ts","normalized signal")
        self.normalized_dict = {}
        # get data of current events on and offset times
        self.event_data = []
        self.perievent_window = None
        # options sent from perievent window
        self.perievent_options_dict = {}
        self.export_window = None
        self.peak_options_window = None
        self.recent_peak_values = []
        # # dictionary of complete (all trials) data around event
        # # key is subject name and value is a dictionary of key: event value: trials
        # self.filtered_data = {}
        # remember previous event from perievent window
        self.current_perievent = None
        # remember selected trials (integers indicating trial number)
        self.current_trials = []
        # remember how many were there
        self.total_current_trials = 0
        # key is subject, value is a dict where key is event and value is list of most recent trials for that event
        self.trials_dict = {}
        
        # initialize export settings
        self.raw_export = False
        self.trimmed_export = False
        self.downsampled_export = False
        self.normalized_export = False
        self.perievent_export = False
        self.spikes_export = False
        self.save_plots = False
        self.export_path = ""
        self.export_begining = ""
        # add that to perievent options
        self.perievent_options_dict["export_path"] = self.export_path
        self.perievent_options_dict["file_beginning"] = self.export_begining
        
        # keep track id save peaks window has to be opened later
        # that is in case user wants to export both: perievents and spikes
        self.open_spikes_later = False
        # if user wanted to run spikes on batch
        self.batch_peaks = False
        # if user wanted to run perievent on batch
        self.batch_perievent = False
        
        self.disable_buttons_signal.connect(self.disable_buttons)
        self.enable_buttons_signal.connect(self.enable_buttons)
        
        self.bold_label_stylesheet = "QLabel {font-weight: bold}"
        self.bold_cb_stylesheet = "QCheckBox {font-weight: bold}"

        # create gui
        self.setupGUI()

    # receives a list with dictionary of signal and control channel names 
    def setupGUI(self):
        self.layout = QVBoxLayout()
        self.layout.setContentsMargins(0,0,0,0)
        self.setLayout(self.layout)
        # area on the left with options
        self.splitter = QSplitter()
        self.splitter.setOrientation(Qt.Horizontal)
        self.layout.addWidget(self.splitter)
############################################
# Options; left side for
        # create widget with options
        self.options_widget = QWidget()
        # create main layout for widget with options
        self.options_main_layout = QVBoxLayout()
        self.options_main_layout.setContentsMargins(10,10,0,10)
        self.options_layout = QFormLayout()
        self.options_layout.setVerticalSpacing(15)
        self.experiment_name_text = QLineEdit(self.preview_init_params[0]["selected_experiment"])
        self.experiment_name_text.setReadOnly(True)
        self.experiment_label = QLabel("Experiment")
        self.experiment_label.setStyleSheet(self.bold_label_stylesheet)
        self.options_layout.addRow(self.experiment_label,self.experiment_name_text)
        # add drop down menu for subjects from the list
        self.subject_comboBox = QComboBox()
        self.subject_comboBox.addItems([self.preview_init_params[0]["subject_name"]])
        # self.subject_comboBox.currentIndexChanged.connect(self.on_subject_change)
        self.subject_label = QLabel("Subject")
        self.subject_label.setStyleSheet(self.bold_label_stylesheet)
        self.options_layout.addRow(self.subject_label, self.subject_comboBox)
        self.subject_group_name = QLineEdit("")
        self.options_layout.addRow("Subject group name", self.subject_group_name)
        self.last_raw_ts = round(self.last_raw_ts,2)
        self.last_raw_ts_text = QLineEdit(str(self.last_raw_ts))
        self.last_raw_ts_text.setReadOnly(True)
        self.options_layout.addRow("Total seconds", self.last_raw_ts_text)
        # adjust trimming when subject changed
        self.trim_beginning_select = QComboBox()
        if len(self.events_from_current_subject) > 0:
            self.trim_beginning_select.addItems(["Trim first seconds", "Trim by event"])
        else:
            self.trim_beginning_select.addItems(["Trim first seconds"])
        self.trim_beginning_select.currentIndexChanged.connect(self.begin_trim_changed)
        self.trim_beginning_sec = QLineEdit("0")
        self.trim_beginning_sec.setValidator(QtGui.QIntValidator())
        self.trim_beginning_event = QComboBox()
        self.trim_beginning_event.addItems(self.events_from_current_subject)
        self.trim_beginning_event.setEnabled(False)
        self.trim_ending_select = QComboBox()
        if len(self.events_from_current_subject) > 0:
            self.trim_ending_select.addItems(["Trim last seconds", "Trim by event"])
        else:
            self.trim_ending_select.addItems(["Trim last seconds"])
        self.trim_ending_select.currentIndexChanged.connect(self.end_trim_changed)
        self.trim_ending_sec = QLineEdit("0")
        self.trim_ending_sec.setValidator(QtGui.QIntValidator())
        self.trim_end_event = QComboBox()
        self.trim_end_event.addItems(self.events_from_current_subject)
        self.trim_end_event.setEnabled(False)
        self.trim_begin_label = QLabel("Trim beginning")
        self.trim_begin_label.setStyleSheet(self.bold_label_stylesheet)
        self.trim_end_label = QLabel("Trim ending")
        self.trim_end_label.setStyleSheet(self.bold_label_stylesheet)
        self.options_layout.addRow(self.trim_begin_label,self.trim_beginning_select)
        self.options_layout.addRow("Seconds",self.trim_beginning_sec)
        self.options_layout.addRow("Event to trim",self.trim_beginning_event)
        self.options_layout.addRow(self.trim_end_label,self.trim_ending_select)
        self.options_layout.addRow("Seconds",self.trim_ending_sec)
        self.options_layout.addRow("Event to trim",self.trim_end_event)
        self.event_from_data_comboBox = QComboBox()
        self.event_from_data_comboBox.addItems(["---",*self.events_from_current_subject])
        self.show_event_label = QLabel("Show Event")
        self.show_event_label.setStyleSheet(self.bold_label_stylesheet)
        self.options_layout.addRow(self.show_event_label,self.event_from_data_comboBox)       
        self.event_name_text = QLineEdit("")
        self.event_name_text.setToolTip("i.e, Tone")
        self.options_layout.addRow("Event name (optional)",self.event_name_text)
        self.event2_from_data_comboBox = QComboBox()
        self.event2_from_data_comboBox.addItems(["---",*self.events_from_current_subject])
        self.showevent2_label = QLabel("Show Event2")
        self.showevent2_label.setStyleSheet(self.bold_label_stylesheet)
        self.options_layout.addRow(self.showevent2_label,self.event2_from_data_comboBox)       
        self.event2_name_text = QLineEdit("")
        self.event2_name_text.setToolTip("i.e, Tone")
        self.options_layout.addRow("Event2 name (optional)",self.event2_name_text)
        self.downsample_cb = QCheckBox("Downsampling")
        self.downsample_cb.setStyleSheet(self.bold_cb_stylesheet)
        self.perform_label = QLabel("Perform:")
        self.options_layout.addRow(self.perform_label,self.downsample_cb)
        self.perform_label.setStyleSheet(self.bold_label_stylesheet)
        self.normalize_cb = QCheckBox("Normalization (+ Downsampling)")
        self.normalize_cb.setStyleSheet(self.bold_cb_stylesheet)
        self.options_layout.addRow("",self.normalize_cb)
        # choose what to show on the plot
        self.raw_plot_cb = QCheckBox("Raw data")
        # self.raw_plot_cb.setChecked(True)
        self.separate_signal_contol_cb = QCheckBox("Separate signal and control")
        self.downsampled_plot_cb = QCheckBox("Downsampled data")
        self.downsampled_plot_cb.setEnabled(False)
        self.normalized_plot_cb = QCheckBox("Normalized data")
        self.normalized_plot_cb.setEnabled(False)
        self.show_on_plot_label = QLabel("Show:")
        self.show_on_plot_label.setStyleSheet(self.bold_label_stylesheet)
        # self.options_layout.addRow(self.show_on_plot_label,self.raw_plot_cb)
        self.options_layout.addRow(self.show_on_plot_label,self.downsampled_plot_cb)
        # self.options_layout.addRow("",self.downsampled_plot_cb)
        self.options_layout.addRow("",self.separate_signal_contol_cb)
        self.options_layout.addRow("",self.normalized_plot_cb)
              
        # nest the inner layout into the outer layout
        self.options_main_layout.addLayout(self.options_layout)
        
        # add apply button
        self.apply_btn = QPushButton("Plot/Apply")
        self.options_main_layout.addWidget(self.apply_btn)
        # add perievent button
        self.perievent_analysis_btn = QPushButton("Perievent analysis")
        if len(self.events_from_current_subject) == 0:
            self.perievent_analysis_btn.setEnabled(False)
        self.options_main_layout.addWidget(self.perievent_analysis_btn)
        # add find peaks button
        self.find_peaks_btn = QPushButton("Spike Detection")
        self.options_main_layout.addWidget(self.find_peaks_btn)
        # add export data button
        self.export_data_btn = QPushButton('Export Data')
        self.options_main_layout.addWidget(self.export_data_btn)  
        # add check polynomial fit button
        self.check_poly_btn = QPushButton("*Check Polynomial Fitting")
        self.check_btn_stylesheet = \
            ".QPushButton{\n" \
            + "padding-left: 5px;\n" \
            + "padding-right: 5px;\n" \
            + "padding-top: 7px;\n"\
            + "padding-bottom: 5px;\n" \
            + "margin-left:20px;\n" \
            + "margin-right:20px;\n" \
            + "margin-top:5px;\n" \
            + "margin-bottom:5px;\n" \
            + "font-size: normal;\n" \
            + "font-weight: normal;\n" \
            + "}" 
        self.check_poly_btn.setStyleSheet(self.check_btn_stylesheet)
        self.options_main_layout.addWidget(self.check_poly_btn)
               
        self.options_widget.setLayout(self.options_main_layout)
        self.splitter.addWidget(self.options_widget)
#####################################################       
# Plots; biggest part of preview (top right)
        # area on the top right for the plots
        self.splitter2 = QSplitter()
        self.splitter2.setOrientation(Qt.Vertical)
        self.splitter.addWidget(self.splitter2)
        # Plot area widget
        self.plot_area_widget = QWidget()
        self.plot_area_widget_layout = QVBoxLayout()
        self.canvas = fpExplorer.MplCanvas(self, width=12, height=8, dpi=100)
        # Create toolbar, passing canvas as first parament, parent (self, the MainWindow) as second.
        self.toolbar = NavigationToolbar(self.canvas, self)
        # add to layout
        self.plot_area_widget_layout.addWidget(self.toolbar)
        self.plot_area_widget_layout.addWidget(self.canvas)
        self.plot_area_widget.setLayout(self.plot_area_widget_layout)
        self.splitter2.addWidget(self.plot_area_widget)
        # create widget for available trials
        self.trials_widget = QWidget()
        self.trials_layout = QHBoxLayout()
        self.trials_widget.setLayout(self.trials_layout)
        self.current_trials_widget = QWidget()
        # self.current_trials_layout = QHBoxLayout()
        # self.current_trials_widget.setLayout(self.current_trials_layout)
        # self.trials_layout.addWidget(self.current_trials_widget)
        self.splitter2.addWidget(self.trials_widget)
        self.trials_button_group = []
        
        # plot initial plot of raw data of the first subject
        self.plot_raw(self.canvas,self.options["subject"],False,False,("",""))
##############################################
# Previous/Next buttons on the bottom        
        # add next/previous buttons at the bottom
        self.preview_buttons_layout = QHBoxLayout()
        self.preview_buttons_layout.setAlignment(Qt.AlignRight)
        self.preview_buttons_widget = QWidget()
        self.include_cb = QCheckBox("Include Subject")
        self.next_btn = QPushButton("Next")
        self.previous_btn = QPushButton("Previous")
        # enable these only if there are more than 1 subjects
        if len(self.preview_init_params[0]["subject_paths"]) < 2:
            self.next_btn.setEnabled(False)
            self.previous_btn.setEnabled(False)
#        self.previous_btn.setFixedWidth(120)
#        self.next_btn.setFixedWidth(120)
        self.preview_buttons_layout.addWidget(self.include_cb)
        self.preview_buttons_layout.addWidget(self.previous_btn)
        self.preview_buttons_layout.addWidget(self.next_btn)
        self.preview_buttons_widget.setLayout(self.preview_buttons_layout)
        self.splitter2.addWidget(self.preview_buttons_widget)
        
#        self.splitter2.setSizes([int(self.height()*0.95), int(self.height()*0.05)])
        
        # connect buttons to fpExplorer_functions
        self.downsample_cb.stateChanged.connect(self.adjust_downsampled_plot_cb)
        self.normalize_cb.stateChanged.connect(self.adjust_normalized_plot_cb)
        self.downsampled_plot_cb.stateChanged.connect(self.adjuct_separate)
        self.raw_plot_cb.stateChanged.connect(self.adjust_separate_control_signal_cb)
        self.separate_signal_contol_cb.stateChanged.connect(self.adjust_downsampled2separate)
        self.apply_btn.clicked.connect(self.apply_btn_clicked)
        self.check_poly_btn.clicked.connect(self.check_poly_btn_clicked)
        self.find_peaks_btn.clicked.connect(self.peaks_btn_clicked)
        self.export_data_btn.clicked.connect(self.export_data_btn_clicked)
        self.perievent_analysis_btn.clicked.connect(self.perievent_analysis_btn_clicked)
        self.next_btn.clicked.connect(self.next_btn_clicked)
        self.previous_btn.clicked.connect(self.previous_btn_clicked)
        # disable run on batch button
        self.parent_window.run_on_batch_btn.setEnabled(False)
  
    @pyqtSlot()   
    def disable_buttons(self):
        self.apply_btn.setEnabled(False)
        self.check_poly_btn.setEnabled(False)
        self.find_peaks_btn.setEnabled(False)
        self.perievent_analysis_btn.setEnabled(False)
        self.export_data_btn.setEnabled(False)
        self.next_btn.setEnabled(False)
        self.previous_btn.setEnabled(False)
        self.parent_window.disable_buttons()
        QApplication.processEvents() 
    
    @pyqtSlot()     
    def enable_buttons(self):
        self.apply_btn.setEnabled(True)
        self.check_poly_btn.setEnabled(True)
        self.find_peaks_btn.setEnabled(True)
        if len(self.events_from_current_subject) > 0:
            self.perievent_analysis_btn.setEnabled(True)
        self.export_data_btn.setEnabled(True)
        self.next_btn.setEnabled(True)
        self.previous_btn.setEnabled(True)
        self.parent_window.enable_all_buttons()
        # all in parent window but run on batch
        self.parent_window.run_on_batch_btn.setEnabled(False)
        QApplication.processEvents() 
            
        
    def begin_trim_changed(self):
        if self.trim_beginning_select.currentText() == "Trim by event":
            self.trim_beginning_sec.setEnabled(False)
            self.trim_beginning_event.setEnabled(True)
        if self.trim_beginning_select.currentText() == "Trim first seconds":
            self.trim_beginning_sec.setEnabled(True)
            self.trim_beginning_event.setEnabled(False)
        
    def end_trim_changed(self):
        if self.trim_ending_select.currentText() == "Trim by event":
            self.trim_ending_sec.setEnabled(False)
            self.trim_end_event.setEnabled(True)
        if self.trim_ending_select.currentText() == "Trim last seconds":
            self.trim_ending_sec.setEnabled(True)
            self.trim_end_event.setEnabled(False)
            
    def reset_export_settings(self):
        self.raw_export = False
        self.trimmed_export = False
        self.downsampled_export = False
        self.normalized_export = False
        self.perievent_export = False
        self.save_plots = False
    
    # function to read options for plotting on the left side of the plot    
    def read_options(self):
        # always check if user wnats to include the current subject
        if self.include_cb.isChecked():
            # add that subject to the list of selected_subjects
            if self.options["subject"] not in self.parent_window.selected_subjects:
                self.parent_window.selected_subjects.append(self.options["subject"])
        # if it is unchecked and used to be in selected, remove from selected
        else:
            if self.options["subject"] in self.parent_window.selected_subjects:
               self.parent_window.selected_subjects.remove(self.options["subject"]) 
        # set current subject
        self.options["subject"] = self.subject_comboBox.currentText()
        self.options["subject_group_name"] = self.subject_group_name.text()
        # update value
        self.group_names_dict[self.options["subject"]] = self.options["subject_group_name"]

        # get trimming settings
        # check what method was selected
        if self.trim_beginning_select.currentText() == "Trim first seconds":
            try:
                trim_beginning = int(self.trim_beginning_sec.text())
            except:
                trim_beginning = 0
                self.trim_beginning_sec.setText("0")
        else: # if trim by event
            trim_begin_event = self.trim_beginning_event.currentText()
            begin_evt_data_onsets = self.get_event_on_off(self.events_dict[self.options["subject"]], trim_begin_event)[0]
            trim_beginning = 0 # assign zero in case there is no such event
            if len(begin_evt_data_onsets) > 0:
                # use the first onset as trim begin
                trim_beginning = math.ceil(begin_evt_data_onsets[0])
            self.trim_beginning_sec.setText(str(trim_beginning))
        if self.trim_ending_select.currentText() == "Trim last seconds":
            try:                
                trim_end = int(self.trim_ending_sec.text())
            except:                
                trim_end = 0
                self.trim_ending_sec.setText("0")
        else: # if trim by event
            trim_end_event = self.trim_end_event.currentText()
            end_evt_data_onsets = self.get_event_on_off(self.events_dict[self.options["subject"]], trim_end_event)[0]
            trim_end = 0 # assign zero in case there is no such event
            if len(end_evt_data_onsets) > 0:
                # use last onset as trim begin
                trim_end = math.ceil(self.last_raw_ts - end_evt_data_onsets[-1])
            self.trim_ending_sec.setText(str(trim_end))
        # if it was not trimmed yet, add to dictionary
        if self.options["subject"] not in self.trimmed_raw_data_dict:
            # key is the subject name, value is a list
            # first element of that list is beginning and end seconds to trim
            # second element is trimmed data dict ts:timestamps, signal:data,control:data
            self.trimmed_raw_data_dict[self.options["subject"]] = [(trim_beginning,trim_end),self.trim_raw_data(self.raw_data_dict[self.options["subject"]],
                                                                                                                      trim_beginning,
                                                                                                                      trim_end
                                                                                                                      )]
        else:   # if already in trimmed check previous trimming settings
            trimmed = self.trimmed_raw_data_dict[self.options["subject"]]
            begin,end = trimmed[0]
            if trim_beginning != begin or trim_end != end:
                # if new trimming params, replace the lists in dict
                self.trimmed_raw_data_dict[self.options["subject"]] = [(trim_beginning,trim_end),self.trim_raw_data(self.raw_data_dict[self.options["subject"]],
                                                                                                                      trim_beginning,
                                                                                                                      trim_end
                                                                                                                      )]
        if len(self.trimmed_raw_data_dict[self.options["subject"]][1]["ts"]) == 0:
            self.show_info_dialog("No data left after trimming.\nTry to trim with different values.")
            return False
        else:
            # set current trimming
            self.options["trim_start"] = self.trim_beginning_sec.text()
            self.options["trim_end"] = self.trim_ending_sec.text()
            self.options["trim_begin_evt"] = self.trim_beginning_event.currentText()
            self.options["trim_end_evt"] = self.trim_end_event.currentText()
            
            
            # get trimmed timestamps and update new session duration
            trimmed_ts = self.trimmed_raw_data_dict[self.options["subject"]][1]["ts"]
            total_trimmed = trimmed_ts[-1]-trimmed_ts[0]
            # start time from zero
            ts_reset = [i*total_trimmed/len(trimmed_ts) for i in range(len(trimmed_ts))]
            last_ts = round(ts_reset[-1],2)
            self.last_raw_ts_text.setText(str(last_ts))
            
            # set event
            self.options["event"] = self.event_from_data_comboBox.currentText()
            self.options["event_name"] = self.event_name_text.text()
            self.options["event2"] = self.event2_from_data_comboBox.currentText()
            self.options["event2_name"] = self.event2_name_text.text()
            # reset the event list to make sure it always has current events
            self.event_data = []
            if self.options["event"] != "---":
                evt = self.get_event_on_off(self.events_dict[self.options["subject"]], self.options["event"])
                if len(evt[0]) == 0:
                    self.show_info_dialog("Some "+self.options["event"]+" event data is missing.")
                self.event_data.append(evt)
            if self.options["event2"] != "---":
                evt = self.get_event_on_off(self.events_dict[self.options["subject"]], self.options["event2"])
                if len(evt[0])==0:
                    self.show_info_dialog("Some "+self.options["event2"]+" event data is missing.")
                self.event_data.append(evt)
            # if downsample was selected
            if self.downsample_cb.isChecked() or self.downsampled_export == True or self.save_plots == True or self.separate_signal_contol_cb.isChecked():
                # add to downdampled dict
                self.downsampled_dict[self.options["subject"]] = fpExplorer_functions.downsample(self.trimmed_raw_data_dict[self.options["subject"]][1],
                                            self.settings_dict[0]["downsample"])
            # if normalize was selected
            if self.normalize_cb.isChecked() or self.normalized_export == True:
                # always downsample first before normalizing
                # downsample
                self.downsampled_dict[self.options["subject"]] = fpExplorer_functions.downsample(self.trimmed_raw_data_dict[self.options["subject"]][1],
                                                                                        self.settings_dict[0]["downsample"])
                # check settings for method to normalize
                if self.settings_dict[0]["normalization"] == "Modified Polynomial Fitting":                   
                    # normalize from downsampled
                    self.normalized_dict[self.options["subject"]] = fpExplorer_functions.normalize_dff(
                                                                                        self.downsampled_dict[self.options["subject"]],
                                                                                        self.settings_dict[0]["show_norm_as"],
                                                                                        self.settings_dict[0]["filter"],
                                                                                        self.settings_dict[0]["filter_window"])
                if self.settings_dict[0]["normalization"] == "Standard Polynomial Fitting":                
                    # normalize from downsampled
                    self.normalized_dict[self.options["subject"]] = fpExplorer_functions.normalize_pMat(
                                                                                        self.downsampled_dict[self.options["subject"]],
                                                                                        self.settings_dict[0]["show_norm_as"],
                                                                                        self.settings_dict[0]["filter"],
                                                                                        self.settings_dict[0]["filter_window"])
            # check what to show on the plot
            self.options["plot_separate"] = True if self.separate_signal_contol_cb.isChecked() else False 
            self.options["plot_downsampled"] = True if self.downsampled_plot_cb.isChecked() else False
            self.options["plot_normalized"] = True if self.normalized_plot_cb.isChecked() else False
            if self.options["plot_separate"] == False and self.options["plot_downsampled"] == False and self.options["plot_normalized"] == False:
                self.raw_plot_cb.setChecked(True)
            self.options["plot_raw"] = True if self.raw_plot_cb.isChecked() else False
            if self.raw_plot_cb.isChecked(): # if it was set as checked the first time, uncheck it in the background(it is not visible anymore)
                self.raw_plot_cb.setChecked(False)
            
            return True

    def clear_include_trials_bar(self):
        new_trials_widget = QWidget()
        new_trials_layout = QHBoxLayout()
        new_trials_widget.setLayout(new_trials_layout)
        # replace with updated widget
        self.trials_layout.replaceWidget(self.current_trials_widget,new_trials_widget)
        self.current_trials_widget.deleteLater()
        self.current_trials_widget = new_trials_widget
            
        
    # function to plot with user selected options       
    def apply_btn_clicked(self):
        self.disable_buttons_signal.emit()
        if self.read_options() == True: 
            wrong_event_order_info = "If you want to show just one event on the plots,\nselect your event as first event.\nAnd leave Event2 empty."
            # before you proceed with ploting check correct order of events if only one event is selected
            if self.options["event"] == "---" and self.options["event2"] != "---":
                self.show_info_dialog(wrong_event_order_info)
            # plot
            # if just raw was checked
            elif (self.options["plot_raw"]==True and self.options["plot_separate"]==False
                and self.options["plot_downsampled"]==False and self.options["plot_normalized"]==False):
                if self.options["event"] == "---":
                    self.plot_trimmed(self.canvas,
                                        self.options["subject"],
                                        self.trimmed_raw_data_dict[self.options["subject"]][1],
                                        self.raw_export,
                                        self.save_plots,
                                        (self.export_path,self.export_begining))
                else:
                    custom_event_name = self.options["event"] if len(self.options["event_name"])==0 else self.options["event_name"]
                    custom_event_name2 = self.options["event2"] if len(self.options["event2_name"])==0 else self.options["event2_name"]
                    self.plot_with_event(self.canvas,
                                        self.options["subject"],
                                        self.trimmed_raw_data_dict[self.options["subject"]][1],
                                        custom_event_name,
                                        custom_event_name2,
                                        self.event_data,
                                        self.raw_export,
                                        self.save_plots,
                                        (self.export_path,self.export_begining))
            # if only downsampled was checked
            elif (self.options["plot_raw"]==False and self.options["plot_separate"]==False and self.options["subject"] in self.downsampled_dict and self.options["plot_downsampled"]==True 
                and self.options["plot_normalized"]==False):
                if self.options["event"] == "---": # no event selected
                    self.plot_downsampled_alone(self.canvas,
                                            self.options,
                                            self.downsampled_dict[self.options["subject"]],
                                            self.downsampled_export,
                                            self.save_plots,
                                            (self.export_path,self.export_begining),
                                            self.settings_dict)
                else:   # with event
                    custom_event_name = self.options["event"] if len(self.options["event_name"])==0 else self.options["event_name"]
                    custom_event_name2 = self.options["event2"] if len(self.options["event2_name"])==0 else self.options["event2_name"]
                    self.plot_downsampled_alone_with_event(self.canvas,
                                            self.options,
                                            self.downsampled_dict[self.options["subject"]],
                                            custom_event_name,
                                            custom_event_name2,
                                            self.event_data,
                                            self.downsampled_export,
                                            self.save_plots,
                                            (self.export_path,self.export_begining),
                                            self.settings_dict)
            # if only normalized was checked
            elif (self.options["plot_raw"]==False and self.options["plot_separate"]==False and self.options["plot_downsampled"]==False
                and self.options["subject"] in self.normalized_dict and self.options["plot_normalized"]==True):
                if self.options["event"] == "---": # no event selected
                    fpExplorer_functions.plot_normalized_alone(self.canvas,
                                                self.options,
                                                self.normalized_dict[self.options["subject"]],
                                                self.normalized_export,
                                                self.save_plots,
                                                (self.export_path,self.export_begining),
                                                self.settings_dict)
                else:   # with event
                    custom_event_name = self.options["event"] if len(self.options["event_name"])==0 else self.options["event_name"]
                    custom_event_name2 = self.options["event2"] if len(self.options["event2_name"])==0 else self.options["event2_name"]
                    fpExplorer_functions.plot_normalized_alone_with_event(self.canvas,
                                            self.options,
                                            self.normalized_dict[self.options["subject"]],
                                            custom_event_name,
                                            custom_event_name2,
                                            self.event_data,
                                            self.normalized_export,
                                            self.save_plots,
                                            (self.export_path,self.export_begining),
                                            self.settings_dict)
            # if downsampled and normalized were selected
            elif (self.options["plot_raw"]==False and self.options["plot_separate"]==False and self.options["subject"] in self.downsampled_dict and self.options["plot_downsampled"]==True
                and self.options["subject"] in self.normalized_dict and self.options["plot_normalized"]==True):
                if self.options["event"] == "---": # no event selected
                    fpExplorer_functions.plot_downsampled_and_normalized_alone(self.canvas,
                                            self.options["subject"],
                                            self.downsampled_dict[self.options["subject"]],
                                            self.normalized_dict[self.options["subject"]],
                                            self.settings_dict[0]["show_norm_as"])
                else:   # with event
                    custom_event_name = self.options["event"] if len(self.options["event_name"])==0 else self.options["event_name"]
                    custom_event_name2 = self.options["event2"] if len(self.options["event2_name"])==0 else self.options["event2_name"]
                    fpExplorer_functions.plot_downsampled_and_normalized_with_event(self.canvas,
                                            self.options["subject"],
                                            self.downsampled_dict[self.options["subject"]],
                                            self.normalized_dict[self.options["subject"]],
                                            self.settings_dict[0]["show_norm_as"],
                                            custom_event_name,
                                            custom_event_name2,
                                            self.event_data)
            # if user wants signal and control separately
            elif (self.options["plot_raw"]==False and self.options["plot_separate"]==True 
                and self.options["plot_downsampled"]==False and self.options["plot_normalized"]==False):
                # if there is downsampled data, plot that, if not raw
                downsampled_to_plot = {}           
                if self.options["subject"] in self.downsampled_dict:
                    downsampled_to_plot = self.downsampled_dict[self.options["subject"]]
                    self.show_info_dialog("Separate signal and control plots\nshow downsampled data.")
                if self.options["event"] == "---": # no event selected
                    fpExplorer_functions.plot_separate_only(self.canvas,
                                            self.options["subject"],
                                            self.trimmed_raw_data_dict[self.options["subject"]][1],
                                            downsampled_to_plot)
                else:   # with event
                    custom_event_name = self.options["event"] if len(self.options["event_name"])==0 else self.options["event_name"]
                    custom_event_name2 = self.options["event2"] if len(self.options["event2_name"])==0 else self.options["event2_name"]
                    fpExplorer_functions.plot_separate_with_event(self.canvas,
                                            self.options["subject"],
                                            self.trimmed_raw_data_dict[self.options["subject"]][1],
                                            downsampled_to_plot,
                                            custom_event_name,
                                            custom_event_name2,
                                            self.event_data)
            # if user wants signal and control separately and specificaly chose downsample
            elif (self.options["plot_raw"]==False and self.options["plot_separate"]==True 
                and self.options["subject"] in self.downsampled_dict and self.options["plot_downsampled"]==True 
                and self.options["plot_normalized"]==False):
                self.show_info_dialog("Separate signal and control plots\nwill show downsampled data.")
                if self.options["event"] == "---": # no event selected
                    fpExplorer_functions.plot_separate_only(self.canvas,
                                            self.options["subject"],
                                            self.trimmed_raw_data_dict[self.options["subject"]][1],
                                            self.downsampled_dict[self.options["subject"]])
                else:   # with event
                    custom_event_name = self.options["event"] if len(self.options["event_name"])==0 else self.options["event_name"]
                    custom_event_name2 = self.options["event2"] if len(self.options["event2_name"])==0 else self.options["event2_name"]
                    fpExplorer_functions.plot_separate_with_event(self.canvas,
                                            self.options["subject"],
                                            self.trimmed_raw_data_dict[self.options["subject"]][1],
                                            self.downsampled_dict[self.options["subject"]],
                                            custom_event_name,
                                            custom_event_name2,
                                            self.event_data)
            # if user wants signal and control separately and normalized
            elif (self.options["plot_raw"]==False and self.options["plot_separate"]==True 
                and self.options["plot_downsampled"]==False 
                and self.options["subject"] in self.normalized_dict and self.options["plot_normalized"]==True):
                # if there is downsampled data, plot that, if not raw
                downsampled_to_plot = {}
                if self.options["subject"] in self.downsampled_dict:
                    downsampled_to_plot = self.downsampled_dict[self.options["subject"]]
                    self.show_info_dialog("Separate signal and control plots\nwill show downsampled data.")
                if self.options["event"] == "---": # no event selected
                    fpExplorer_functions.plot_separate_with_normalized(self.canvas,
                                            self.options["subject"],
                                            self.trimmed_raw_data_dict[self.options["subject"]][1],
                                            downsampled_to_plot,
                                            self.normalized_dict[self.options["subject"]],
                                            self.settings_dict[0]["show_norm_as"])
                else:   # with event
                    custom_event_name = self.options["event"] if len(self.options["event_name"])==0 else self.options["event_name"]
                    custom_event_name2 = self.options["event2"] if len(self.options["event2_name"])==0 else self.options["event2_name"]
                    fpExplorer_functions.plot_separate_with_normalized_with_event(self.canvas,
                                            self.options["subject"],
                                            self.trimmed_raw_data_dict[self.options["subject"]][1],
                                            downsampled_to_plot,
                                            self.normalized_dict[self.options["subject"]],
                                            self.settings_dict[0]["show_norm_as"],
                                            custom_event_name,
                                            custom_event_name2,
                                            self.event_data)
            # if user wants signal and control separately and normalized and specificaly chose downsample
            elif (self.options["plot_raw"]==False and self.options["plot_separate"]==True 
                and self.options["plot_downsampled"]==True and  self.options["subject"] in self.downsampled_dict
                and self.options["subject"] in self.normalized_dict and self.options["plot_normalized"]==True):
                self.show_info_dialog("Separate signal and control plots\nwill show downsampled data.")
                if self.options["event"] == "---": # no event selected
                    fpExplorer_functions.plot_separate_with_normalized(self.canvas,
                                            self.options["subject"],
                                            self.trimmed_raw_data_dict[self.options["subject"]][1],
                                            self.downsampled_dict[self.options["subject"]],
                                            self.normalized_dict[self.options["subject"]],
                                            self.settings_dict[0]["show_norm_as"])
                else:   # with event
                    custom_event_name = self.options["event"] if len(self.options["event_name"])==0 else self.options["event_name"]
                    custom_event_name2 = self.options["event2"] if len(self.options["event2_name"])==0 else self.options["event2_name"]
                    fpExplorer_functions.plot_separate_with_normalized_with_event(self.canvas,
                                            self.options["subject"],
                                            self.trimmed_raw_data_dict[self.options["subject"]][1],
                                            self.downsampled_dict[self.options["subject"]],
                                            self.normalized_dict[self.options["subject"]],
                                            self.settings_dict[0]["show_norm_as"],
                                            custom_event_name,
                                            custom_event_name2,
                                            self.event_data)
        if self.current_perievent != None:
            self.clear_include_trials_bar()
        self.enable_buttons_signal.emit()

    def next_btn_clicked(self):
        # check if user wnats to include the current subject
        if self.include_cb.isChecked():
            # add that subject to the list of selected_subjects
            if self.subject_comboBox.currentText() not in self.parent_window.selected_subjects:
                self.parent_window.selected_subjects.append(self.subject_comboBox.currentText())
        # if it is unchecked and used to be in selected, remove from selected
        else:
            if self.subject_comboBox.currentText() in self.parent_window.selected_subjects:
               self.parent_window.selected_subjects.remove(self.subject_comboBox.currentText()) 
        # get index of current subject
        current_subject_index = self.preview_init_params[0][0]["subject_names"].index(self.subject_comboBox.currentText())
        if current_subject_index < len(self.preview_init_params[0][0]["subject_names"])-1:
            next_index = current_subject_index + 1
        else:
            next_index = 0
        # set current subject to the next subject
        self.subject_comboBox.setCurrentText(self.preview_init_params[0][0]["subject_names"][next_index])
        
        
    def previous_btn_clicked(self):
        # check if user wnats to include the current subject
        if self.include_cb.isChecked():
            # add that subject to the list of selected_subjects
            if self.subject_comboBox.currentText() not in self.parent_window.selected_subjects:
                self.parent_window.selected_subjects.append(self.subject_comboBox.currentText())
        # if it is unchecked and used to be in selected, remove from selected
        else:
            if self.subject_comboBox.currentText() in self.parent_window.selected_subjects:
               self.parent_window.selected_subjects.remove(self.subject_comboBox.currentText()) 
        # get index of current subject
        current_subject_index = self.preview_init_params[0][0]["subject_names"].index(self.subject_comboBox.currentText())
        if current_subject_index - 1 >= 0:
            previous_index = current_subject_index - 1
        else:
            previous_index = len(self.preview_init_params[0][0]["subject_names"])-1
        # set current subject to the previous subject
        self.subject_comboBox.setCurrentText(self.preview_init_params[0][0]["subject_names"][previous_index])
        
    
    # check and set trimming for current subject
    def check_trimming(self):
        self.options["subject"] = self.subject_comboBox.currentText()
        # if subject was not previewed, set trimming to 0
        if self.options["subject"] not in self.trimmed_raw_data_dict:
            self.trim_beginning_sec.setText("0")
            self.trim_ending_sec.setText("0")
        else:
            # set trimming to most recent
            trimmed = self.trimmed_raw_data_dict[self.options["subject"]]
            begin,end = trimmed[0]
            self.trim_beginning_sec.setText(str(begin))
            self.trim_ending_sec.setText(str(end))
    # check and set event for current subject      
    def check_events(self):
        self.options["subject"] = self.subject_comboBox.currentText()
        # if subject was not previewed yet, read data and add to dict
        if self.options["subject"] not in self.raw_data_dict:
            self.get_raw_data(self.options["subject"],
                          self.parent_window.preview_params[1][self.options["subject"]])
            if self.raw_data_dict[self.options["subject"]] == None:
                return
        # check events for new data
        # first check if there are the same events
        previous_events = self.events_from_current_subject
        self.events_from_current_subject = self.get_events(self.events_dict[self.options["subject"]]) 
        same = True
        if len(previous_events) != self.events_from_current_subject:
            same = False
        else: # check if they are the same events
            for el in self.events_from_current_subject:
                if el not in previous_events:
                    same = False
        self.trim_beginning_select.setCurrentText("Trim first seconds")
        self.trim_ending_select.setCurrentText("Trim last seconds")
        # set class variables first:
        # set current trimming
        self.options["trim_start"] = self.trim_beginning_sec.text()
        self.options["trim_end"] = self.trim_ending_sec.text()
        self.options["trim_begin_evt"] = ""
        self.options["trim_end_evt"] = ""
        if same == True:
            self.event_from_data_comboBox.setCurrentText(self.options["event"])
            # set event name also to previous
            if self.options["event"] != "---":
                self.event_name_text.setText(self.options["event_name"])
            self.event2_from_data_comboBox.setCurrentText(self.options["event2"])
            # set event name also to previous
            if self.options["event2"] != "---":
                self.event2_name_text.setText(self.options["event2_name"])
        else: # if events changed
            # check trimming events
            if self.options["trim_begin_evt"] in self.events_from_current_subject:
                # update combo box first
                self.trim_beginning_event.clear()
                self.trim_beginning_event.addItems(self.events_from_current_subject)
                # set to previous value
                self.trim_beginning_event.setCurrentText(self.options["trim_begin_evt"])
            else: # update to new events set trimming back to seconds
                self.trim_beginning_event.clear()
                self.trim_beginning_event.addItems(self.events_from_current_subject)
#                self.trim_beginning_select.setCurrentText("Trim first seconds")
            if self.options["trim_end_evt"] in self.events_from_current_subject:
                # update combo box first
                self.trim_end_event.clear()
                self.trim_end_event.addItems(self.events_from_current_subject)
                # set to previous value
                self.trim_end_event.setCurrentText(self.options["trim_end_evt"])
            else: # update to new events set trimming back to seconds
                self.trim_end_event.clear()
                self.trim_end_event.addItems(self.events_from_current_subject)
#                self.trim_ending_select.setCurrentText("Trim last seconds")
                
            # if previously selected event is in current events, set that filed the same as previously
            if self.options["event"] in self.events_from_current_subject or self.options["event"] == "---":
                # update combo box first
                self.event_from_data_comboBox.clear()
                self.event_from_data_comboBox.addItems(["---",*self.events_from_current_subject])
                # clear any previous event names
                self.event_name_text.setText("")
                # set to previous event
                self.event_from_data_comboBox.setCurrentText(self.options["event"])
                # set event name also to previous
                if self.options["event"] != "---":
                    self.event_name_text.setText(self.options["event_name"])
            else:
                self.event_from_data_comboBox.clear()
                self.event_from_data_comboBox.addItems(["---",*self.events_from_current_subject])
                # clear any previous event names
                self.event_name_text.setText("")
            # now check for event 2
            if self.options["event2"] in self.events_from_current_subject or self.options["event2"] == "---":
                # update combo box first
                self.event2_from_data_comboBox.clear()
                self.event2_from_data_comboBox.addItems(["---",*self.events_from_current_subject])
                # clear any previous event names
                self.event2_name_text.setText("")
                # set to previous event
                self.event2_from_data_comboBox.setCurrentText(self.options["event2"])
                # set event name also to previous
                if self.options["event2"] != "---":
                    self.event2_name_text.setText(self.options["event2_name"])
            else:
                self.event2_from_data_comboBox.clear()
                self.event2_from_data_comboBox.addItems(["---",*self.events_from_current_subject])
                # clear any previous event names
                self.event2_name_text.setText("")
        
    # check trimming and event at the same time        
    def options_pre_check(self):
        # always check if user wnats to include the current subject
        if self.include_cb.isChecked():
            # add that subject to the list of selected_subjects
            if self.options["subject"] not in self.parent_window.selected_subjects:
                self.parent_window.selected_subjects.append(self.options["subject"])
        # if it is unchecked and used to be in selected, remove from selected
        else:
            if self.options["subject"] in self.parent_window.selected_subjects:
               self.parent_window.selected_subjects.remove(self.options["subject"]) 
        self.check_events()
        self.check_trimming()
        # set include to actual status
        if self.subject_comboBox.currentText() in self.parent_window.selected_subjects:
            self.include_cb.setChecked(True)
        else:
            self.include_cb.setChecked(False)
            
    def on_subject_change(self):
        self.options_pre_check()
        if self.raw_data_dict[self.options["subject"]] == None:
            return
        # update last timestamp
        # get last timestamp in seconds
        self.last_raw_ts = fpExplorer_functions.get_last_timestamp(
                                     self.raw_data_dict[self.options["subject"]],
                                     self.preview_init_params[0][0]["signal_name"],
                                     )
        # read new subject's frequency and update suggested downsampled rate (only if it is different)
        new_fs = fpExplorer_functions.get_frequency(self.raw_data_dict[self.preview_init_params[0][0]["subject_names"][0]],self.preview_init_params[0][0]["signal_name"])
        if new_fs != self.current_fs:
            # self.suggested_downsample_samples = int(int(self.current_fs)*DEFAULT_DOWNSAMPLE_PCT/100)
            self.suggested_downsample_samples = DEFAULT_HZ
            self.settings_dict[0]['downsample'] = self.suggested_downsample_samples
            self.suggested_downsample_rate = int(int(self.current_fs)/self.suggested_downsample_samples)
            self.settings[0]["entered_downsample"] = self.suggested_downsample_rate
        if self.subject_comboBox.currentText() in self.group_names_dict:
            self.subject_group_name.setText(self.group_names_dict[self.subject_comboBox.currentText()])
        else:
            self.subject_group_name.setText("")
        # clear long lasting check boxes to start faster
        self.downsample_cb.setChecked(False)
        self.normalize_cb.setChecked(False)
        # plot raw by default
        self.raw_plot_cb.setChecked(True)
        # clear "include trials
        # # create new buttons
        new_trials_widget = QWidget()
        new_trials_layout = QHBoxLayout()
        new_trials_widget.setLayout(new_trials_layout)
        # replace with updated widget
        self.trials_layout.replaceWidget(self.current_trials_widget,new_trials_widget)
        self.current_trials_widget.deleteLater()
        self.current_trials_widget = new_trials_widget
        # self.current_perievent = None
        self.current_trials = []
        self.trials_button_group = []
        self.total_current_trials = 0
        # plot raw, trimmed data with or without event
        self.apply_btn_clicked()
            
    def adjust_downsampled_plot_cb(self):
        if self.downsample_cb.isChecked():
            self.downsampled_plot_cb.setChecked(True)
            self.downsampled_plot_cb.setEnabled(True)
        else:
            self.downsampled_plot_cb.setChecked(False)
            self.downsampled_plot_cb.setEnabled(False)
            
    def adjust_normalized_plot_cb(self):
        if self.normalize_cb.isChecked():
            self.normalized_plot_cb.setChecked(True)
            self.normalized_plot_cb.setEnabled(True)
        else:
            self.normalized_plot_cb.setChecked(False)
            self.normalized_plot_cb.setEnabled(False)
            
    def adjust_separate_control_signal_cb(self):
        if self.raw_plot_cb.isChecked():
            self.separate_signal_contol_cb.setChecked(False)
            
    def adjuct_separate(self):
        if self.downsampled_plot_cb.isChecked():
            self.separate_signal_contol_cb.setChecked(False)
    def adjust_downsampled2separate(self):
        if self.separate_signal_contol_cb.isChecked():
            self.downsampled_plot_cb.setChecked(False)
            self.raw_plot_cb.setChecked(False)
            
    def check_poly_btn_clicked(self):
        self.disable_buttons_signal.emit()
        if self.normalize_cb.isChecked():
            if self.read_options() == True:
                slope_intercept_dict = self.show_polynomial_fitting(self.canvas,
                                                self.settings_dict[0], 
                                                self.downsampled_dict[self.options["subject"]],
                                                self.options["subject"],
                                                False,
                                                "")
                # show the equation of fitted lines
                signal_slope_intercept = slope_intercept_dict["signal_slope_intercept"]
                control_slope_intercept = slope_intercept_dict["control_slope_intercept"]
                signal_b_direction = "+ " if signal_slope_intercept[1] >= 0 else "- "
                control_b_direction = "+ " if control_slope_intercept[1] >= 0 else "- "
                info_text = "\nSignal fitted: y = "+str(round(signal_slope_intercept[0],4))+"x "+signal_b_direction+str(round(signal_slope_intercept[1],2))\
                            +"\nControl fitted: y = "+str(round(control_slope_intercept[0],4))+"x "+control_b_direction+str(round(control_slope_intercept[1],2))+"\n"
                # if they slopes of signal and control have oposite directions, add suggestion
                if signal_slope_intercept[0]*control_slope_intercept[0] < 0:
                    info_text = "!!! Warning !!!\nWe suggest using Modified Polynomial Fitting to normalize this data.\nYou can change normalization method under application's Settings.\n"+info_text
                self.show_info_dialog(info_text)
        else: # prompt user to check normalize first
            self.show_info_dialog("Please check Normalize check box first.")
        self.enable_buttons_signal.emit()
            
    def export_data_btn_clicked(self):
        # pass a list with first element=subject
        self.export_window = fpExplorer.ExportDataWindow(self.parent_window,self,[self.subject_comboBox.currentText(),
                                                    self.preview_init_params,self.export_path,self.export_begining])
        self.export_window.got_export_selection_sig.connect(self.export_options_received)
        self.export_window.show()
        self.disable_buttons_signal.emit()
        
    def perievent_analysis_btn_clicked(self):
        self.disable_buttons_signal.emit()
        # update group name
        self.options["subject_group_name"] = self.subject_group_name.text()
        # update value
        self.group_names_dict[self.options["subject"]] = self.options["subject_group_name"]
        # show window with options to plot perievent
        # pass a list with first element=subject
        # second=events for current subjects
        # third element=raw data for current subject
        # check if on batch and check common events
        common_events = []
        subjects_event_sets = []
        if self.batch_perievent == True:
            for subject in self.parent_window.batch_export_settings_dict["batch_subjects"]:
                if subject not in self.raw_data_dict:
                    self.get_raw_data(subject,self.parent_window.preview_params[1][subject])
                if self.raw_data_dict[subject] != None:
                    evt_list = self.get_events(self.events_dict[subject])
                    subjects_event_sets.append(set(evt_list))
            # find intersection of all subjects events (common events)  
            common_events = list(set.intersection(*subjects_event_sets))  

        self.perievent_window = fpExplorer.PeriEventOptionsWindow(self.parent_window,[self.options["subject"],
                                                                           self.events_from_current_subject,
                                                                           self.raw_data_dict[self.options["subject"]]],
                                                                           self.perievent_options_dict,
                                                                           self.batch_perievent,
                                                                           common_events)
        self.perievent_window.got_peri_event_options_sig.connect(self.get_perievent_options)
        self.perievent_window.show()
            
            

    def read_trials(self):
        checked = []
        for btn in self.trials_button_group:
            # print(btn.text(),btn.isChecked())
            if btn.isChecked():
                checked.append(int(btn.text()))
        return checked
            
    def read_trials_range(self):
        checked = []
        start_trial = abs(int(self.from_trial.text()))
        end_trial = abs(int(self.to_trial.text()))
        if start_trial > end_trial:
            self.show_info_dialog("Last trial cannot be earlier than the first.")
            return checked
        if start_trial < 1 or start_trial >= self.total_current_trials or end_trial < 1 or end_trial > self.total_current_trials:
            self.show_info_dialog("All trials need to be within the available range.")
            return checked
        checked = np.linspace(start_trial,end_trial,(end_trial-start_trial+1),dtype = int)
        return checked


    @pyqtSlot(list) 
    # receives a list with perievent options read from perievent window
    def get_perievent_options(self,perievent_options):
        self.perievent_options_dict = perievent_options[0]
        print("perievent options",self.perievent_options_dict)
        if self.current_perievent == self.perievent_options_dict['event'] and len(self.current_trials)>0: 
            # read the trials radio buttons
            try:
                if self.total_current_trials < 10:
                    print("Less than 10 events")
                    if len(self.read_trials()) > 0:
                        self.current_trials = self.read_trials()
                else:
                    print("More than 9 events")
                    # if user did not enter the right range, don't override previous
                    if len(self.read_trials_range()) > 0:
                        self.current_trials = self.read_trials_range()
                print("Trials read from previous",self.current_trials)
                if self.options["subject"] in self.trials_dict:
                    self.trials_dict[self.options["subject"]][self.current_perievent] = self.current_trials
                else:
                    self.trials_dict[self.options["subject"]] = {self.current_perievent:self.current_trials}
            except:
                if self.total_current_trials < 10:
                    self.trials_button_group = []
                    # create new buttons
                    new_trials_widget = QWidget()
                    new_trials_layout = QHBoxLayout()
                    new_trials_widget.setLayout(new_trials_layout)
                    text_label = QLabel("Include Trials:")
                    new_trials_layout.addWidget(text_label)
                    for i in range(self.total_current_trials):
                            # create btn
                            # add button to a group
                            # add to list
                            # add to layout
                            btn = QCheckBox(str(i+1))
                            if i+1 in self.current_trials:
                                btn.setChecked(True)
                            new_trials_layout.addWidget(btn)
                            self.trials_button_group.append(btn)                        
                    # replace with updated widget
                    self.trials_layout.replaceWidget(self.current_trials_widget,new_trials_widget)
                    self.current_trials_widget.deleteLater()
                    self.current_trials_widget = new_trials_widget
                    print("Trials not read from previous",self.current_trials)
                else: # for more than 9 trials select the range of trials
                    # create new buttons
                    new_trials_widget = QWidget()
                    new_trials_layout = QHBoxLayout()
                    new_trials_layout.setAlignment(Qt.AlignRight)
                    new_trials_widget.setLayout(new_trials_layout)
                    text_label = QLabel("Include From Trial:")
                    new_trials_layout.addWidget(text_label)
                    self.from_trial = QLineEdit(str(1))
                    self.from_trial.setValidator(QtGui.QIntValidator())
                    new_trials_layout.addWidget(self.from_trial)
                    till_trial_label = QLabel("To Trial:")
                    new_trials_layout.addWidget(till_trial_label)
                    self.to_trial = QLineEdit(str(self.total_current_trials))
                    self.to_trial.setValidator(QtGui.QIntValidator())
                    self.to_trial.setToolTip("Max: "+str(self.total_current_trials))
                    new_trials_layout.addWidget(self.to_trial)                   
                    # replace with updated widget
                    self.trials_layout.replaceWidget(self.current_trials_widget,new_trials_widget)
                    self.current_trials_widget.deleteLater()
                    self.current_trials_widget = new_trials_widget
                    print("Trials not read from previous",self.current_trials)

        if self.batch_perievent == True:    
            if "perievent" in self.parent_window.batch_export_settings_dict:
                if self.parent_window.batch_export_settings_dict["perievent"] == True:
                    # collect subject data if possible: list of tuples: subject+each all trials as dataframe
                    all_subjects_peri_normalized_dfs = []
                    # list of tuples: subject+df
                    all_subjects_zscored_dfs = []
                    # check if user wanted data for each subject
                    if self.parent_window.batch_export_settings_dict["export_for_single_subjects"] == True:
                        for i in range(len(self.parent_window.batch_export_settings_dict["batch_subjects"])):
                            subject = self.parent_window.batch_export_settings_dict["batch_subjects"][i]
                            # add the group name to group names dictionary if it was not there or update group names from batch options
                            self.group_names_dict[subject] = self.parent_window.batch_export_settings_dict["batch_subjects_group_names"][i]
                            if subject not in self.raw_data_dict: # if data has not beed read yet
                                self.get_raw_data(subject,
                                                  self.parent_window.preview_params[1][subject])
                                if self.raw_data_dict[subject] != None:
                                    # if subject data has not been trimmed
                                    # if it was not trimmed yet, add to dictionary with zero trimming
                                    self.trimmed_raw_data_dict[subject] = [(0,0),self.trim_raw_data(self.raw_data_dict[self.options["subject"]],
                                                                                                                      0,
                                                                                                                      0
                                                                                                                      )]
                            if subject not in self.trimmed_raw_data_dict and self.raw_data_dict[subject] != None: # if raw data read but not trimmed assume trimming 0
                                # key is the subject name, value is a list
                                # first element of that list is beginning and end seconds to trim
                                # second element is trimmed data dict ts:timestamps, signal:data,control:data
                                self.trimmed_raw_data_dict[subject] = [(0,0),self.trim_raw_data(self.raw_data_dict[self.options["subject"]],
                                                                                                                      0,
                                                                                                                      0
                                                                                                                      )]
                            if self.raw_data_dict[subject] != None:
                                # filter around trimmed
                                data_dict = self.filter_data_around_event(self.raw_data_dict[self.options["subject"]],
                                                    self.events_dict[self.options["subject"]],
                                                    self.perievent_options_dict,
                                                    self.settings_dict)
                                # # if run on batch was selected, select all trials by default
                                # self.current_trials = [i+1 for i in range(len(data.streams[self.preview_init_params[0][0]["signal_name"]].filtered))]
                                # if run on batch was selected, get selected trials
                                self.current_trials = self.perievent_options_dict["trials"]
                                # add to dictionary
                                self.trials_dict[subject] = {self.perievent_options_dict["event"]:self.perievent_options_dict["trials"]}
                                # analyze
                                analyzed_perievent_dict = fpExplorer_functions.analyze_perievent_data(data_dict,
                                                                                    self.current_trials,
                                                                                self.perievent_options_dict,
                                                                                self.settings_dict,
                                                                                "",
                                                                                "")
                                # check if there is any data to plot 
                                if (len(data_dict["signal"]) > 0) and (len(data_dict["control"]) > 0):
                                    # show buttons of available trials
                                    self.total_current_trials = len(data_dict["signal"])
                                    if self.current_perievent == None:
                                        self.current_perievent = self.perievent_options_dict['event']
                                        # add trials to view
                                        self.current_trials_layout = QHBoxLayout()
                                        self.current_trials_widget.setLayout(self.current_trials_layout)
                                        text_label = QLabel("Include Trials:")
                                        self.current_trials_layout.addWidget(text_label)
                                        for i in range(self.total_current_trials):
                                            # create btn
                                            # add button to a group
                                            # add to list
                                            # add to layout
                                            btn = QCheckBox(str(i+1))
                                            btn.setChecked(True)
                                            self.current_trials_layout.addWidget(btn)
                                            self.trials_button_group.append(btn)
                                            # self.current_trials.append(i+1)
                                        self.trials_layout.addWidget(self.current_trials_widget)
                                        print("Before plot raw",self.current_trials)

                                    # check if selected current trials are in subject's data
                                    available = [i+1 for i in range(len(data.streams[self.preview_init_params[0][0]["signal_name"]].filtered))]
                                    all_trials_available = True
                                    for el in self.current_trials:
                                        if el not in available:
                                            all_trials_available = False
                                            break
                                    if all_trials_available == False:
                                        self.show_info_dialog("Some of the selected trials not present\nin subject "+subject+" data.")
                                    else:
                                        # create subfolder with subject name
                                        subject_subfolder = os.path.join(self.parent_window.batch_export_settings_dict["dump_path"],subject)
                                        if not os.path.exists(subject_subfolder):
                                            try:
                                                os.mkdir(subject_subfolder)
                                            except:
                                                self.show_info_dialog("Problem creating subfolder")
                                                # save under main folder then
                                                subject_subfolder = self.parent_window.batch_export_settings_dict["dump_path"]
                                        # plot and save normalized preview
                                        all_trials_df = self.plot_raw_perievents(self.canvas,
                                                  self.options["subject"],
                                                  data_dict,
                                                  self.current_trials,
                                                  self.perievent_options_dict,
                                                  self.settings_dict,
                                                  self.perievent_options_dict["export"],
                                                  self.save_plots,
                                                  self.group_names_dict[self.options["subject"]],
                                                  (self.export_path,self.export_begining)) 
                                        if self.parent_window.batch_export_settings_dict["export_group_data"] == True:
                                            all_subjects_peri_normalized_dfs.append((subject,all_trials_df))
                                        if self.perievent_options_dict["plot_avg"] == True:
                                            fpExplorer_functions.plot_perievent_average_alone(self.canvas,
                                                      self.options["subject"],
                                                      self.current_trials,
                                                      self.perievent_options_dict,
                                                      analyzed_perievent_dict,
                                                      self.perievent_options_dict["export"],
                                                      self.save_plots,
                                                      self.group_names_dict[self.options["subject"]],
                                                      self.settings_dict,
                                                      "",
                                                      "",
                                                      (self.export_path,self.export_begining))
                                        if (self.perievent_options_dict["plot_zscore"] == True or self.perievent_options_dict["plot_zscore_trials"] == True 
                                            or self.perievent_options_dict["plot_auc"] == True): # we need zcsore data for auc
                                            if self.perievent_options_dict["plot_zscore"] == True:
                                                z_score_df = fpExplorer_functions.plot_perievent_zscore_alone(self.canvas,
                                                                            subject,
                                                                            self.perievent_options_dict,
                                                                            analyzed_perievent_dict,
                                                                            self.parent_window.batch_export_settings_dict["export_for_single_subjects"],
                                                                            self.save_plots,
                                                                            self.group_names_dict[subject],
                                                                            self.settings_dict,
                                                                            (subject_subfolder,self.parent_window.batch_export_settings_dict["file_begin"]))
                                            else:
        #                                    if self.perievent_options_dict["plot_zscore_trials"] == True:
                                                z_score_df = fpExplorer_functions.plot_perievent_zscore_with_trials_alone(self.canvas,
                                                                            subject,
                                                                            self.perievent_options_dict,
                                                                            analyzed_perievent_dict,
                                                                            self.parent_window.batch_export_settings_dict["export_for_single_subjects"],
                                                                            self.save_plots,
                                                                            self.group_names_dict[subject],
                                                                            self.settings_dict,
                                                                            (subject_subfolder,self.parent_window.batch_export_settings_dict["file_begin"]))
                                            if self.parent_window.batch_export_settings_dict["export_group_data"] == True: # save for later in order to not repeat tasks
                                                all_subjects_zscored_dfs.append((subject,z_score_df))
                                        if self.perievent_options_dict["plot_auc"] == True:
                                            fpExplorer_functions.plot_perievent_auc_alone(self.canvas,
                                                                        subject,
                                                                        self.perievent_options_dict,
                                                                        analyzed_perievent_dict,
                                                                        self.parent_window.batch_export_settings_dict["export_for_single_subjects"],
                                                                        self.save_plots,
                                                                        self.group_names_dict[subject],
                                                                        self.settings_dict,
                                                                        (subject_subfolder,self.parent_window.batch_export_settings_dict["file_begin"]))
                                else:
                                    self.show_info_dialog("Not enough data that satisfy your request.\nTry changing event or times around the event.")
                                    break
                    # check if user wanted group analysis
                    if self.parent_window.batch_export_settings_dict["export_group_data"] == True:
                        if self.parent_window.batch_export_settings_dict["export_for_single_subjects"] == True:
                            if self.perievent_options_dict["plot_avg"] == True:
                                if len(all_subjects_peri_normalized_dfs) > 0: # if there already is data from single subjects
                                    fpExplorer_functions.get_batch_perievent_normalized(self.canvas,
                                                                                 all_subjects_peri_normalized_dfs,
                                                                                 self.parent_window.batch_export_settings_dict["batch_subjects_group_names"],
                                                                                 self.perievent_options_dict,
                                                                                 self.settings_dict,
                                                                                 (self.parent_window.batch_export_settings_dict["dump_path"],self.parent_window.batch_export_settings_dict["file_begin"]))
                            if self.perievent_options_dict["plot_zscore"] == True:
                                if len(all_subjects_zscored_dfs) > 0: # if there already is data from single subjects
                                    fpExplorer_functions.get_batch_perievent_zscored(self.canvas,
                                                                          all_subjects_zscored_dfs,
                                                                          self.parent_window.batch_export_settings_dict["batch_subjects_group_names"],
                                                                          self.perievent_options_dict,
                                                                          self.settings_dict,
                                                                          (self.parent_window.batch_export_settings_dict["dump_path"],self.parent_window.batch_export_settings_dict["file_begin"]))
                            if self.perievent_options_dict["plot_zscore_trials"] == True:
                                if len(all_subjects_zscored_dfs) > 0: # if there already is data from single subjects
                                    fpExplorer_functions.get_batch_perievent_zscored_with_trials(self.canvas,
                                                                          all_subjects_zscored_dfs,
                                                                          self.parent_window.batch_export_settings_dict["batch_subjects_group_names"],
                                                                          self.perievent_options_dict,
                                                                          self.settings_dict,
                                                                          (self.parent_window.batch_export_settings_dict["dump_path"],self.parent_window.batch_export_settings_dict["file_begin"]))
                            if self.perievent_options_dict["plot_auc"] == True:
                                if len(all_subjects_zscored_dfs) > 0: # if there already is data from single subjects
                                    fpExplorer_functions.get_batch_perievent_auc(self.canvas,
                                                                  all_subjects_zscored_dfs,
                                                                  self.parent_window.batch_export_settings_dict["batch_subjects_group_names"],
                                                                 self.perievent_options_dict,
                                                                 self.settings_dict,
                                                                 (self.parent_window.batch_export_settings_dict["dump_path"],self.parent_window.batch_export_settings_dict["file_begin"]))
                        else: # user did not select single subject data
                            if (self.perievent_options_dict["plot_avg"] == True or self.perievent_options_dict["plot_zscore"] == True 
                                or self.perievent_options_dict["plot_zscore_trials"] == True or self.perievent_options_dict["plot_auc"] == True):
                                # collect subject data if possible: list of tuples: subject+each all trials as dataframe
                                all_subjects_peri_normalized_dfs = []
                                # list of tuples: subject+df
                                all_subjects_zscored_dfs = []
                                for subject in self.parent_window.batch_export_settings_dict["batch_subjects"]:
                                    if subject not in self.raw_data_dict: # if data has not beed read yet
                                        self.get_raw_data(subject,
                                                          self.parent_window.preview_params[1][subject])
                                        if self.raw_data_dict[subject] != None:
                                            # if subject data has not been trimmed
                                            # if it was not trimmed yet, add to dictionary with zero trimming
                                            self.trimmed_raw_data_dict[subject] = [(0,0),self.trim_raw_data(self.raw_data_dict[self.options["subject"]],
                                                                                                                      0,
                                                                                                                      0
                                                                                                                      )]
                                    if self.raw_data_dict[subject] != None:
                                        if subject not in self.trimmed_raw_data_dict: # if raw data read but not trimmed assume trimming 0
                                            # key is the subject name, value is a list
                                            # first element of that list is beginning and end seconds to trim
                                            # second element is trimmed data dict ts:timestamps, signal:data,control:data
                                            self.trimmed_raw_data_dict[subject] = [(0,0),self.trim_raw_data(self.raw_data_dict[self.options["subject"]],
                                                                                                                      0,
                                                                                                                      0
                                                                                                                      )]
                                        # filter around trimmed
                                        data_dict = self.filter_data_around_event(self.raw_data_dict[self.options["subject"]],
                                                    self.events_dict[self.options["subject"]],
                                                    self.perievent_options_dict,
                                                    self.settings_dict)
                                        # if run on batch was selected, select all trials by default
                                        # self.current_trials = [i+1 for i in range(len(data.streams[self.preview_init_params[0][0]["signal_name"]].filtered))]
                                        # if run on batch was selected, get selected trials
                                        self.current_trials = self.perievent_options_dict["trials"]
                                        # add to dictionary
                                        self.trials_dict[subject] = {self.perievent_options_dict["event"]:self.perievent_options_dict["trials"]}
                                        # check if there is any data to plot 
                                        if (len(data_dict["signal"]) > 0) and (len(data_dict["control"]) > 0):
                                            # check if selected current trials are in subject's data
                                            available = [i+1 for i in range(len(data_dict["signal"]))]
                                            all_trials_available = True
                                            for el in self.current_trials:
                                                if el not in available:
                                                    all_trials_available = False
                                                    break
                                            if all_trials_available == False:
                                                self.show_info_dialog("Some of the selected trials not present\nin subject "+subject+" data.")
                                            else:
                                                if (self.perievent_options_dict["plot_zscore"] == True or self.perievent_options_dict["plot_zscore_trials"] == True 
                                                    or self.perievent_options_dict["plot_auc"] == True):
                                                    # analyze
                                                    analyzed_perievent_dict = fpExplorer_functions.analyze_perievent_data(data_dict,
                                                                                                    self.current_trials,
                                                                                                self.perievent_options_dict,
                                                                                                self.settings_dict,
                                                                                                "",
                                                                                                "")

                                                    # show buttons of available trials
                                                    self.total_current_trials = len(data_dict["signal"])
                                                    if self.current_perievent == None:
                                                        self.current_perievent = self.perievent_options_dict['event']
                                                        # add trials to view
                                                        self.current_trials_layout = QHBoxLayout()
                                                        self.current_trials_widget.setLayout(self.current_trials_layout)
                                                        text_label = QLabel("Include Trials:")
                                                        self.current_trials_layout.addWidget(text_label)
                                                        for i in range(len(data_dict["signal"])):
                                                            # create btn
                                                            # add button to a group
                                                            # add to list
                                                            # add to layout
                                                            btn = QCheckBox(str(i+1))
                                                            btn.setChecked(True)
                                                            self.current_trials_layout.addWidget(btn)
                                                            self.trials_button_group.append(btn)
                                                            # self.current_trials.append(i+1)
                                                        self.trials_layout.addWidget(self.current_trials_widget)
                                                        
                                                if self.perievent_options_dict["plot_avg"] == True:
                                                    # perieventnormalized preview
                                                    all_trials_df = self.plot_raw_perievents(self.canvas,
                                                                                            self.options["subject"],
                                                                                            data_dict,
                                                                                            self.current_trials,
                                                                                            self.perievent_options_dict,
                                                                                            self.settings_dict,
                                                                                            self.perievent_options_dict["export"],
                                                                                            self.save_plots,
                                                                                            self.group_names_dict[self.options["subject"]],
                                                                                            (self.export_path,self.export_begining)) 
                                                    all_subjects_peri_normalized_dfs.append((subject,all_trials_df))
                                                if (self.perievent_options_dict["plot_zscore"] == True or self.perievent_options_dict["plot_zscore_trials"] == True 
                                                    or self.perievent_options_dict["plot_auc"] == True):
                                                    if self.perievent_options_dict["plot_zscore"] == True:
                                                        z_score_df = fpExplorer_functions.plot_perievent_zscore_alone(self.canvas,
                                                                                subject,
                                                                                self.perievent_options_dict,
                                                                                analyzed_perievent_dict,
                                                                                self.parent_window.batch_export_settings_dict["export_for_single_subjects"],
                                                                                self.save_plots,
                                                                                self.parent_window.preview_params[1][subject],
                                                                                self.settings_dict,
                                                                                (self.parent_window.batch_export_settings_dict["dump_path"],self.parent_window.batch_export_settings_dict["file_begin"]))
                                                        all_subjects_zscored_dfs.append((subject,z_score_df))
                                                    else:
        #                                            if self.perievent_options_dict["plot_zscore_trials"] == True:
                                                        z_score_df = fpExplorer_functions.plot_perievent_zscore_with_trials_alone(self.canvas,
                                                                                subject,
                                                                                self.perievent_options_dict,
                                                                                analyzed_perievent_dict,
                                                                                self.parent_window.batch_export_settings_dict["export_for_single_subjects"],
                                                                                self.save_plots,
                                                                                self.parent_window.preview_params[1][subject],
                                                                                self.settings_dict,
                                                                                (self.parent_window.batch_export_settings_dict["dump_path"],self.parent_window.batch_export_settings_dict["file_begin"]))
                                                        all_subjects_zscored_dfs.append((subject,z_score_df))
                                        else:
                                            self.show_info_dialog("Not enough data that satisfy your request.\nTry changing event or times around the event.")
                                            self.perievent_window.enable_buttons()
                                            break
                                
                                if self.perievent_options_dict["plot_avg"] == True:
                                    if len(all_subjects_peri_normalized_dfs) > 0:
                                        fpExplorer_functions.get_batch_perievent_normalized(self.canvas,
                                                                                 all_subjects_peri_normalized_dfs,
                                                                                 self.parent_window.batch_export_settings_dict["batch_subjects_group_names"],
                                                                                 self.perievent_options_dict,
                                                                                 self.settings_dict,
                                                                                 (self.parent_window.batch_export_settings_dict["dump_path"],self.parent_window.batch_export_settings_dict["file_begin"]))
                                if len(all_subjects_zscored_dfs) > 0:
                                    if self.perievent_options_dict["plot_zscore"] == True:
                                        fpExplorer_functions.get_batch_perievent_zscored(self.canvas,
                                                                          all_subjects_zscored_dfs,
                                                                          self.parent_window.batch_export_settings_dict["batch_subjects_group_names"],
                                                                          self.perievent_options_dict,
                                                                          self.settings_dict,
                                                                          (self.parent_window.batch_export_settings_dict["dump_path"],self.parent_window.batch_export_settings_dict["file_begin"]))
                                    if self.perievent_options_dict["plot_zscore_trials"] == True:
                                        fpExplorer_functions.get_batch_perievent_zscored_with_trials(self.canvas,
                                                                          all_subjects_zscored_dfs,
                                                                          self.parent_window.batch_export_settings_dict["batch_subjects_group_names"],
                                                                          self.perievent_options_dict,
                                                                          self.settings_dict,
                                                                          (self.parent_window.batch_export_settings_dict["dump_path"],self.parent_window.batch_export_settings_dict["file_begin"]))
                                    if self.perievent_options_dict["plot_auc"] == True:
                                        fpExplorer_functions.get_batch_perievent_auc(self.canvas,
                                                                      all_subjects_zscored_dfs,
                                                                      self.parent_window.batch_export_settings_dict["batch_subjects_group_names"],
                                                                     self.perievent_options_dict,
                                                                     self.settings_dict,
                                                                     (self.parent_window.batch_export_settings_dict["dump_path"],self.parent_window.batch_export_settings_dict["file_begin"]))
                    
                # close popup window
                self.perievent_window.close()
                # self.reset_export_settings()
                # reset batch peaks when done
                self.batch_perievent = False
                # update trials according to batch settings
                self.trials_button_group = []
                # create new buttons
                new_trials_widget = QWidget()
                new_trials_layout = QHBoxLayout()
                new_trials_widget.setLayout(new_trials_layout)
                text_label = QLabel("Include Trials:")
                new_trials_layout.addWidget(text_label)
                for i in range(self.total_current_trials):
                        # create btn
                        # add button to a group
                        # add to list
                        # add to layout
                        btn = QCheckBox(str(i+1))
                        if i+1 in self.current_trials:
                            btn.setChecked(True)
                        new_trials_layout.addWidget(btn)
                        self.trials_button_group.append(btn)                        
                # replace with updated widget
                self.trials_layout.replaceWidget(self.current_trials_widget,new_trials_widget)
                self.current_trials_widget.deleteLater()
                self.current_trials_widget = new_trials_widget
                
                # reset open spikes later
                if self.open_spikes_later == True:
                    self.peaks_btn_clicked()
                # close run on batch window if still open
                else: 
                    self.reset_export_settings()
                    if self.parent_window.run_on_batch_window != None:
                        self.parent_window.run_on_batch_window.close()         
                    
        else: # if not on batch
            if self.options["subject"] not in self.trimmed_raw_data_dict: # if raw data read but not trimmed assume trimming 0
                # key is the subject name, value is a list
                # first element of that list is beginning and end seconds to trim
                # second element is trimmed data dict ts:timestamps, signal:data,control:data
                self.trimmed_raw_data_dict[self.options["subject"]] = [(0,0),self.trim_raw_data(self.raw_data_dict[self.options["subject"]],
                                                                                                                      0,
                                                                                                                      0
                                                                                                                      )]
            # filter around event
            data_dict = self.filter_data_around_event(self.raw_data_dict[self.options["subject"]],
                                                    self.events_dict[self.options["subject"]],
                                                    self.perievent_options_dict,
                                                    self.settings_dict)

            # update export path settings
            if len(self.perievent_options_dict["export_path"]) > 0 and os.path.split(self.perievent_options_dict["export_path"])[1]==DEFAULT_EXPORT_FOLDER:
                # check if folder exists first
                try: 
                    os.mkdir(self.perievent_options_dict["export_path"])
                except:
                    if not os.path.exists(self.perievent_options_dict["export_path"]):
                        self.show_info_dialog("Problem creating subfolder")
                self.export_path = self.perievent_options_dict["export_path"]
                self.export_begining = self.perievent_options_dict["file_beginning"]
            
            # check if there is any data to plot 
            if (len(data_dict["signal"]) > 0) and (len(data_dict["control"]) > 0):

                # show buttons of available trials
                print("How many trials are there?",len(data_dict["signal"]))
                self.total_current_trials = len(data_dict["signal"])
                if self.current_perievent == None:
                    self.current_perievent = self.perievent_options_dict['event']
                    if self.total_current_trials < 10:
                        # add trials to view
                        self.current_trials_layout = QHBoxLayout()
                        self.current_trials_widget.setLayout(self.current_trials_layout)
                        text_label = QLabel("Include Trials:")
                        self.current_trials_layout.addWidget(text_label)
                        for i in range(self.total_current_trials):
                            # create btn
                            # add button to a group
                            # add to list
                            # add to layout
                            btn = QCheckBox(str(i+1))
                            btn.setChecked(True)
                            self.current_trials_layout.addWidget(btn)
                            self.trials_button_group.append(btn)
                            self.current_trials.append(i+1)
                        self.trials_layout.addWidget(self.current_trials_widget)
                    else: # for more than 9 trials select the range of trials
                        # create new buttons
                        self.current_trials_layout = QHBoxLayout()
                        self.current_trials_layout.setAlignment(Qt.AlignRight)
                        self.current_trials_widget.setLayout(self.current_trials_layout)
                        text_label = QLabel("Include From Trial:")
                        self.current_trials_layout.addWidget(text_label)
                        self.from_trial = QLineEdit(str(1))
                        self.from_trial.setValidator(QtGui.QIntValidator())
                        self.current_trials_layout.addWidget(self.from_trial)
                        till_trial_label = QLabel("To Trial:")
                        self.current_trials_layout.addWidget(till_trial_label)
                        self.to_trial = QLineEdit(str(self.total_current_trials))
                        self.to_trial.setValidator(QtGui.QIntValidator())
                        self.to_trial.setToolTip("Max: "+str(self.total_current_trials))
                        self.current_trials_layout.addWidget(self.to_trial)  
                        self.trials_layout.addWidget(self.current_trials_widget)
                        self.current_trials = np.linspace(1,self.total_current_trials,self.total_current_trials,dtype=int)            
                    if self.options["subject"] in self.trials_dict:
                        self.trials_dict[self.options["subject"]][self.current_perievent] = self.current_trials
                    else:
                        self.trials_dict[self.options["subject"]] = {self.current_perievent:self.current_trials}
                elif self.current_perievent != self.perievent_options_dict['event']:  # if new event different from previous       
                    # clear current trials and buttons group
                    self.current_trials = []
                    if self.total_current_trials < 10:
                        self.trials_button_group = []
                        # create new buttons
                        new_trials_widget = QWidget()
                        new_trials_layout = QHBoxLayout()
                        new_trials_widget.setLayout(new_trials_layout)
                        text_label = QLabel("Include Trials:")
                        new_trials_layout.addWidget(text_label)
                        if self.options["subject"] not in self.trials_dict:
                            for i in range(self.total_current_trials):
                                # create btn
                                # add button to a group
                                # add to list
                                # add to layout
                                btn = QCheckBox(str(i+1))
                                btn.setChecked(True)
                                new_trials_layout.addWidget(btn)
                                self.trials_button_group.append(btn)
                                self.current_trials.append(i+1)
                            self.trials_dict[self.options["subject"]] = {self.current_perievent:self.current_trials}
                        else: # subject already had some events trials previewed
                            if self.perievent_options_dict['event'] in self.trials_dict[self.options["subject"]]:
                                self.current_trials = self.trials_dict[self.options["subject"]][self.perievent_options_dict['event']]
                                for i in range(self.total_current_trials):
                                    btn = QCheckBox(str(i+1))
                                    if i+1 in self.current_trials:
                                        btn.setChecked(True)
                                    new_trials_layout.addWidget(btn)
                                    self.trials_button_group.append(btn)
                            else: # add event for that subject
                                for i in range(self.total_current_trials):
                                    # create btn
                                    # add button to a group
                                    # add to list
                                    # add to layout
                                    btn = QCheckBox(str(i+1))
                                    btn.setChecked(True)
                                    new_trials_layout.addWidget(btn)
                                    self.trials_button_group.append(btn)
                                    self.current_trials.append(i+1)
                                self.trials_dict[self.options["subject"]][self.perievent_options_dict['event']] = self.current_trials
                        # replace with updated widget
                        self.trials_layout.replaceWidget(self.current_trials_widget,new_trials_widget)
                        self.current_trials_widget.deleteLater()
                        self.current_trials_widget = new_trials_widget
                    else: # for more than 9 trials select the range of trials
                        # create new buttons
                        new_trials_widget = QWidget()
                        new_trials_layout = QHBoxLayout()
                        new_trials_layout.setAlignment(Qt.AlignRight)
                        new_trials_widget.setLayout(new_trials_layout)
                        text_label = QLabel("Include From Trial:")
                        new_trials_layout.addWidget(text_label)
                        if self.options["subject"] not in self.trials_dict:
                            self.from_trial = QLineEdit(str(1))
                            self.from_trial.setValidator(QtGui.QIntValidator())
                            new_trials_layout.addWidget(self.from_trial)
                            till_trial_label = QLabel("To Trial:")
                            new_trials_layout.addWidget(till_trial_label)
                            self.to_trial = QLineEdit(str(self.total_current_trials))
                            self.to_trial.setValidator(QtGui.QIntValidator())
                            self.to_trial.setToolTip("Max: "+str(self.total_current_trials))
                            new_trials_layout.addWidget(self.to_trial)  
                            self.current_trials = np.linspace(1,self.total_current_trials,self.total_current_trials,dtype=int) 
                            self.trials_dict[self.options["subject"]] = {self.current_perievent:self.current_trials}
                        else: # subject already had some events trials previewed
                            if self.perievent_options_dict['event'] in self.trials_dict[self.options["subject"]]:
                                self.current_trials = self.trials_dict[self.options["subject"]][self.perievent_options_dict['event']]
                                self.from_trial = QLineEdit(str(self.current_trials[0]))
                                self.from_trial.setValidator(QtGui.QIntValidator())
                                new_trials_layout.addWidget(self.from_trial)
                                till_trial_label = QLabel("To Trial:")
                                new_trials_layout.addWidget(till_trial_label)
                                self.to_trial = QLineEdit(str(self.current_trials[-1]))
                                self.to_trial.setValidator(QtGui.QIntValidator())
                                self.to_trial.setToolTip("Max: "+str(self.total_current_trials))
                                new_trials_layout.addWidget(self.to_trial)  
                            else: # add event for that subject
                                self.from_trial = QLineEdit(str(1))
                                self.from_trial.setValidator(QtGui.QIntValidator())
                                new_trials_layout.addWidget(self.from_trial)
                                till_trial_label = QLabel("To Trial:")
                                new_trials_layout.addWidget(till_trial_label)
                                self.to_trial = QLineEdit(str(self.total_current_trials))
                                self.to_trial.setValidator(QtGui.QIntValidator())
                                self.to_trial.setToolTip("Max: "+str(self.total_current_trials))
                                new_trials_layout.addWidget(self.to_trial)  
                                self.current_trials = np.linspace(1,self.total_current_trials,self.total_current_trials,dtype=int)    
                                self.trials_dict[self.options["subject"]][self.perievent_options_dict['event']] = self.current_trials
                                         
                        # replace with updated widget
                        self.trials_layout.replaceWidget(self.current_trials_widget,new_trials_widget)
                        self.current_trials_widget.deleteLater()
                        self.current_trials_widget = new_trials_widget
                    # set new event
                    self.current_perievent = self.perievent_options_dict['event']
                else: # if same event and/or new subject
                    if self.total_current_trials < 10:
                        # create new buttons
                        new_trials_widget = QWidget()
                        new_trials_layout = QHBoxLayout()
                        new_trials_widget.setLayout(new_trials_layout)
                        text_label = QLabel("Include Trials:")
                        new_trials_layout.addWidget(text_label)
                        if len(self.trials_button_group) > 0:
                            for btn in self.trials_button_group:
                                new_trials_layout.addWidget(btn)
                        else:
                            if self.options["subject"] not in self.trials_dict:
                                for i in range(self.total_current_trials):
                                    # create btn
                                    # add button to a group
                                    # add to list
                                    # add to layout
                                    btn = QCheckBox(str(i+1))
                                    btn.setChecked(True)
                                    new_trials_layout.addWidget(btn)
                                    self.trials_button_group.append(btn)
                                    self.current_trials.append(i+1)
                                self.trials_dict[self.options["subject"]] = {self.current_perievent:self.current_trials}
                            else: # subject already had some events trials previewed
                                if self.perievent_options_dict['event'] in self.trials_dict[self.options["subject"]]:
                                    self.current_trials = self.trials_dict[self.options["subject"]][self.perievent_options_dict['event']]
                                    for i in range(self.total_current_trials):
                                        btn = QCheckBox(str(i+1))
                                        if i+1 in self.current_trials:
                                            btn.setChecked(True)
                                        new_trials_layout.addWidget(btn)
                                        self.trials_button_group.append(btn)
                                else: # add event for that subject
                                    for i in range(self.total_current_trials):
                                        # create btn
                                        # add button to a group
                                        # add to list
                                        # add to layout
                                        btn = QCheckBox(str(i+1))
                                        btn.setChecked(True)
                                        new_trials_layout.addWidget(btn)
                                        self.trials_button_group.append(btn)
                                        self.current_trials.append(i+1)
                                    self.trials_dict[self.options["subject"]][self.current_perievent] = self.current_trials
                    else: # more than 9 events
                        # create new buttons
                        new_trials_widget = QWidget()
                        new_trials_layout = QHBoxLayout()
                        new_trials_layout.setAlignment(Qt.AlignRight)
                        new_trials_widget.setLayout(new_trials_layout)
                        text_label = QLabel("Include From Trial:")
                        new_trials_layout.addWidget(text_label)
                        if self.options["subject"] not in self.trials_dict:
                            self.from_trial = QLineEdit(str(1))
                            self.from_trial.setValidator(QtGui.QIntValidator())
                            new_trials_layout.addWidget(self.from_trial)
                            till_trial_label = QLabel("To Trial:")
                            new_trials_layout.addWidget(till_trial_label)
                            self.to_trial = QLineEdit(str(self.total_current_trials))
                            self.to_trial.setValidator(QtGui.QIntValidator())
                            self.to_trial.setToolTip("Max: "+str(self.total_current_trials))
                            new_trials_layout.addWidget(self.to_trial)  
                            self.current_trials = np.linspace(1,self.total_current_trials,self.total_current_trials,dtype=int) 
                            self.trials_dict[self.options["subject"]] = {self.current_perievent:self.current_trials}
                        else: # subject already had some events trials previewed
                            if self.perievent_options_dict['event'] in self.trials_dict[self.options["subject"]]:
                                self.current_trials = self.trials_dict[self.options["subject"]][self.perievent_options_dict['event']]
                                self.from_trial = QLineEdit(str(self.current_trials[0]))
                                self.from_trial.setValidator(QtGui.QIntValidator())
                                new_trials_layout.addWidget(self.from_trial)
                                till_trial_label = QLabel("To Trial:")
                                new_trials_layout.addWidget(till_trial_label)
                                self.to_trial = QLineEdit(str(self.current_trials[-1]))
                                self.to_trial.setValidator(QtGui.QIntValidator())
                                self.to_trial.setToolTip("Max: "+str(self.total_current_trials))
                                new_trials_layout.addWidget(self.to_trial) 
                            else: # add event for that subject
                                self.from_trial = QLineEdit(str(1))
                                self.from_trial.setValidator(QtGui.QIntValidator())
                                new_trials_layout.addWidget(self.from_trial)
                                till_trial_label = QLabel("To Trial:")
                                new_trials_layout.addWidget(till_trial_label)
                                self.to_trial = QLineEdit(str(self.total_current_trials))
                                self.to_trial.setValidator(QtGui.QIntValidator())
                                self.to_trial.setToolTip("Max: "+str(self.total_current_trials))
                                new_trials_layout.addWidget(self.to_trial)  
                                self.current_trials = np.linspace(1,self.total_current_trials,self.total_current_trials,dtype=int)    
                                self.trials_dict[self.options["subject"]][self.perievent_options_dict['event']] = self.current_trials

                    # replace with updated widget
                    self.trials_layout.replaceWidget(self.current_trials_widget,new_trials_widget)
                    self.current_trials_widget.deleteLater()
                    self.current_trials_widget = new_trials_widget
                    print("Trials dictionary")
                    print(self.trials_dict)
            ##########################################

                # if user clicked on preview in perievent options window, show preview
                if self.perievent_options_dict["preview"] == True:
                    if self.total_current_trials > 9:
                        self.show_info_dialog("You can preview max 9 events at the same time.\nBut you can include as many as you want in the analysis.")
                    # plot normalized preview
                    self.plot_raw_perievents(self.canvas,
                                                  self.options["subject"],
                                                  data_dict,
                                                  self.current_trials,
                                                  self.perievent_options_dict,
                                                  self.settings_dict,
                                                  self.perievent_options_dict["export"],
                                                  self.save_plots,
                                                  self.group_names_dict[self.options["subject"]],
                                                  (self.export_path,self.export_begining))  
                # if user clicked on analyze in perievent options window, analyze
                if self.perievent_options_dict["analyze"] == True:
                    analyzed_perievent_dict = fpExplorer_functions.analyze_perievent_data(data_dict,
                                                                                self.current_trials,
                                                                               self.perievent_options_dict,
                                                                               self.settings_dict,
                                                                               "",
                                                                               "")
                    # find out what to show on plot
                    # only average
                    if (self.perievent_options_dict["plot_avg"] == True and self.perievent_options_dict["plot_zscore"] == False
                        and self.perievent_options_dict["plot_zscore_trials"] == False and self.perievent_options_dict["plot_auc"] == False):
                        fpExplorer_functions.plot_perievent_average_alone(self.canvas,
                                                      self.options["subject"],
                                                      self.current_trials,
                                                      self.perievent_options_dict,
                                                      analyzed_perievent_dict,
                                                      self.perievent_options_dict["export"],
                                                      self.save_plots,
                                                      self.group_names_dict[self.options["subject"]],
                                                      self.settings_dict,
                                                      "",
                                                      "",
                                                      (self.export_path,self.export_begining))
                    # only zscore with error
                    elif (self.perievent_options_dict["plot_avg"] == False and self.perievent_options_dict["plot_zscore"] == True
                            and self.perievent_options_dict["plot_zscore_trials"] == False and self.perievent_options_dict["plot_auc"] == False):
                        fpExplorer_functions.plot_perievent_zscore_alone(self.canvas,
                                                      self.options["subject"],
                                                      self.perievent_options_dict,
                                                      analyzed_perievent_dict,
                                                      self.perievent_options_dict["export"],
                                                      self.save_plots,
                                                      self.group_names_dict[self.options["subject"]],
                                                      self.settings_dict,
                                                      (self.export_path,self.export_begining))
                    # only zscore with trials
                    elif (self.perievent_options_dict["plot_avg"] == False and self.perievent_options_dict["plot_zscore"] == False
                            and self.perievent_options_dict["plot_zscore_trials"] == True and self.perievent_options_dict["plot_auc"] == False):
                        fpExplorer_functions.plot_perievent_zscore_with_trials_alone(self.canvas,
                                                      self.options["subject"],
                                                      self.perievent_options_dict,
                                                      analyzed_perievent_dict,
                                                      self.perievent_options_dict["export"],
                                                      self.save_plots,
                                                      self.group_names_dict[self.options["subject"]],
                                                      self.settings_dict,
                                                      (self.export_path,self.export_begining))
                    # only auc
                    elif (self.perievent_options_dict["plot_avg"] == False and self.perievent_options_dict["plot_zscore"] == False
                            and self.perievent_options_dict["plot_zscore_trials"] == False and self.perievent_options_dict["plot_auc"] == True):
                        fpExplorer_functions.plot_perievent_auc_alone(self.canvas,
                                                      self.options["subject"],
                                                      self.perievent_options_dict,
                                                      analyzed_perievent_dict,
                                                      self.perievent_options_dict["export"],
                                                      self.save_plots,
                                                      self.group_names_dict[self.options["subject"]],
                                                      self.settings_dict,
                                                      (self.export_path,self.export_begining))
                    # avg and zscore with error
                    elif (self.perievent_options_dict["plot_avg"] == True and self.perievent_options_dict["plot_zscore"] == True
                            and self.perievent_options_dict["plot_zscore_trials"] == False and self.perievent_options_dict["plot_auc"] == False):
                        fpExplorer_functions.plot_perievent_avg_zscore(self.canvas,
                                                      self.options["subject"],
                                                      self.current_trials,
                                                      self.perievent_options_dict,
                                                      analyzed_perievent_dict)
                    # avg and zscore with trials
                    elif (self.perievent_options_dict["plot_avg"] == True and self.perievent_options_dict["plot_zscore"] == False
                            and self.perievent_options_dict["plot_zscore_trials"] == True and self.perievent_options_dict["plot_auc"] == False):
                        fpExplorer_functions.plot_perievent_avg_zscore_trials(self.canvas,
                                                      self.options["subject"],
                                                      self.current_trials,
                                                      self.perievent_options_dict,
                                                      analyzed_perievent_dict)
                    # avg and auc
                    elif (self.perievent_options_dict["plot_avg"] == True and self.perievent_options_dict["plot_zscore"] == False
                            and self.perievent_options_dict["plot_zscore_trials"] == False and self.perievent_options_dict["plot_auc"] == True):
                        fpExplorer_functions.plot_perievent_avg_auc(self.canvas,
                                                      self.options["subject"],
                                                      self.current_trials,
                                                      self.perievent_options_dict,
                                                      analyzed_perievent_dict)
                    #zscore with error and auc
                    elif (self.perievent_options_dict["plot_avg"] == False and self.perievent_options_dict["plot_zscore"] == True
                            and self.perievent_options_dict["plot_zscore_trials"] == False and self.perievent_options_dict["plot_auc"] == True):
                        fpExplorer_functions.plot_perievent_zscore_auc(self.canvas,
                                                      self.options["subject"],
                                                      self.perievent_options_dict,
                                                      analyzed_perievent_dict)
                    #zscore with trials and auc
                    elif (self.perievent_options_dict["plot_avg"] == False and self.perievent_options_dict["plot_zscore"] == False
                            and self.perievent_options_dict["plot_zscore_trials"] == True and self.perievent_options_dict["plot_auc"] == True):
                        fpExplorer_functions.plot_perievent_zscore_trials_auc(self.canvas,
                                                      self.options["subject"],
                                                      self.perievent_options_dict,
                                                      analyzed_perievent_dict)
                    # all 4 plots (zscore with error)
                    elif (self.perievent_options_dict["plot_avg"] == True and self.perievent_options_dict["plot_zscore"] == True
                            and self.perievent_options_dict["plot_zscore_trials"] == False and self.perievent_options_dict["plot_auc"] == True):
                        fpExplorer_functions.plot_all_perievent(self.canvas,
                                                      self.options["subject"],
                                                      self.current_trials,
                                                      self.perievent_options_dict,
                                                      analyzed_perievent_dict)
                    # all 4 plots (zscore with trials)
                    elif (self.perievent_options_dict["plot_avg"] == True and self.perievent_options_dict["plot_zscore"] == False
                            and self.perievent_options_dict["plot_zscore_trials"] == True and self.perievent_options_dict["plot_auc"] == True):
                        fpExplorer_functions.plot_all_perievent_zscore_trials(self.canvas,
                                                      self.options["subject"],
                                                      self.current_trials,
                                                      self.perievent_options_dict,
                                                      analyzed_perievent_dict)
                # if user clicked on analyze in perievent options window, analyze
                if self.perievent_options_dict["export"] == True:
                    analyzed_perievent_dict = fpExplorer_functions.analyze_perievent_data(data_dict,
                                                                                self.current_trials,
                                                                               self.perievent_options_dict,
                                                                               self.settings_dict,
                                                                               "",
                                                                               "")
                    # save downsampled preview
                    self.plot_raw_perievents(self.canvas,
                                                  self.options["subject"],
                                                  data_dict,
                                                  self.current_trials,
                                                  self.perievent_options_dict,
                                                  self.settings_dict,
                                                  self.perievent_options_dict["export"],
                                                  self.save_plots,
                                                  self.group_names_dict[self.options["subject"]],
                                                  (self.export_path,self.export_begining)) 
                    if self.perievent_options_dict["plot_avg"] == True:
                        fpExplorer_functions.plot_perievent_average_alone(self.canvas,
                                                      self.options["subject"],
                                                      self.current_trials,
                                                      self.perievent_options_dict,
                                                      analyzed_perievent_dict,
                                                      self.perievent_options_dict["export"],
                                                      self.save_plots,
                                                      self.group_names_dict[self.options["subject"]],
                                                      self.settings_dict,
                                                      "",
                                                      "",
                                                      (self.export_path,self.export_begining))
                    if self.perievent_options_dict["plot_zscore"] == True:
                        fpExplorer_functions.plot_perievent_zscore_alone(self.canvas,
                                                      self.options["subject"],
                                                      self.perievent_options_dict,
                                                      analyzed_perievent_dict,
                                                      self.perievent_options_dict["export"],
                                                      self.save_plots,
                                                      self.group_names_dict[self.options["subject"]],
                                                      self.settings_dict,
                                                      (self.export_path,self.export_begining))
                    if self.perievent_options_dict["plot_zscore_trials"] == True:
                        fpExplorer_functions.plot_perievent_zscore_with_trials_alone(self.canvas,
                                                      self.options["subject"],
                                                      self.perievent_options_dict,
                                                      analyzed_perievent_dict,
                                                      self.perievent_options_dict["export"],
                                                      self.save_plots,
                                                      self.group_names_dict[self.options["subject"]],
                                                      self.settings_dict,
                                                      (self.export_path,self.export_begining))
                    if self.perievent_options_dict["plot_auc"] == True:
                        fpExplorer_functions.plot_perievent_auc_alone(self.canvas,
                                                      self.options["subject"],
                                                      self.perievent_options_dict,
                                                      analyzed_perievent_dict,
                                                      self.perievent_options_dict["export"],
                                                      self.save_plots,
                                                      self.group_names_dict[self.options["subject"]],
                                                      self.settings_dict,
                                                      (self.export_path,self.export_begining))
                # close popup window
                self.perievent_window.close()
                self.enable_buttons_signal.emit()
                if self.open_spikes_later == False:
                    self.reset_export_settings()
                    if self.export_window != None:
                        self.export_window.close()
                elif self.open_spikes_later == True:
                    self.peaks_btn_clicked()
            else:   # if there was not enough data
                self.show_info_dialog("Not enough data that satisfy your request.\nTry changing event or times around the event.")
                self.perievent_window.enable_buttons()


       
    @pyqtSlot(list) 
    def export_options_received(self,export_options):
        if len(export_options[0]["dump_path"]) == 0 or os.path.exists(export_options[0]["dump_path"]) == False:
            if (len(export_options[0]["dump_path"]) > 0 and os.path.split(export_options[0]["dump_path"])[1] == DEFAULT_EXPORT_FOLDER):
                try: 
                    os.mkdir(export_options[0]["dump_path"])
                except:
                    if not os.path.exists(export_options[0]["dump_path"]):
                        self.show_info_dialog("Problem creating subfolder")
            else:
                self.show_info_dialog("Please provide a valid export path!")
                # close the old window open new window
                self.export_window.close()
                # pass a list with first element=subject
                self.export_window = fpExplorer.ExportDataWindow(self.parent_window,self,[self.subject_comboBox.currentText(),
                                                            self.preview_init_params,self.export_path,self.export_begining])
                self.export_window.got_export_selection_sig.connect(self.export_options_received)
                self.export_window.show()
                self.disable_buttons_signal.emit()
        if os.path.exists(export_options[0]["dump_path"]) == True:
            # read export options
            self.export_path = export_options[0]["dump_path"]
            self.export_begining = export_options[0]["file_begin"]
            self.raw_export = export_options[0]["raw"]
            self.downsampled_export = export_options[0]["downsampled"]
            self.normalized_export = export_options[0]["normalized"]
            self.perievent_export = export_options[0]["perievent"]
            self.spikes_export = export_options[0]["spikes"]
            self.save_plots = export_options[0]["save_plot"]
            self.export_path = export_options[0]["dump_path"]
            self.export_begining = export_options[0]["file_begin"]
            self.export_data()

                
    # function to plot with user selected options       
    def export_data(self):
        if self.read_options() == True:
            wrong_event_order_info = "If you want to show just one event on the plots,\nselect your event as first event.\nAnd leave Event2 empty."
            # before you proceed with ploting check correct order of events  if only one event is selected
            if self.options["event"] == "---" and self.options["event2"] != "---":
                self.show_info_dialog(wrong_event_order_info)
            # save fittings if save all plots was selected
            if self.save_plots == True:
                self.show_polynomial_fitting(self.canvas,
                                                self.settings_dict[0], 
                                                self.downsampled_dict[self.options["subject"]],
                                                self.options["subject"],
                                                self.save_plots,
                                                self.export_path)
            # raw was checked
            if self.raw_export == True:
                if self.options["event"] == "---":
                    self.plot_raw(self.canvas,self.options["subject"],self.raw_export,self.save_plots,(self.export_path,self.options["subject"]))
                else:
                    custom_event_name = self.options["event"] if len(self.options["event_name"])==0 else self.options["event_name"]
                    custom_event_name2 = self.options["event2"] if len(self.options["event2_name"])==0 else self.options["event2_name"]
                    self.plot_with_event(self.canvas,
                                        self.options["subject"],
                                        self.trimmed_raw_data_dict[self.options["subject"]][1],
                                        custom_event_name,
                                        custom_event_name2,
                                        self.event_data,
                                        self.raw_export,
                                        self.save_plots,
                                        (self.export_path,self.export_begining))
            if self.downsampled_export == True:
                if self.options["event"] == "---": # no event selected
                    self.plot_downsampled_alone(self.canvas,
                                            self.options,
                                            self.downsampled_dict[self.options["subject"]],
                                            self.downsampled_export,
                                            self.save_plots,
                                            (self.export_path,self.export_begining),
                                            self.settings_dict)
                else:   # with event
                    custom_event_name = self.options["event"] if len(self.options["event_name"])==0 else self.options["event_name"]
                    custom_event_name2 = self.options["event2"] if len(self.options["event2_name"])==0 else self.options["event2_name"]
                    self.plot_downsampled_alone_with_event(self.canvas,
                                            self.options,
                                            self.downsampled_dict[self.options["subject"]],
                                            custom_event_name,
                                            custom_event_name2,
                                            self.event_data,
                                            self.downsampled_export,
                                            self.save_plots,
                                            (self.export_path,self.export_begining),
                                            self.settings_dict)
            if self.normalized_export == True:
                if self.options["event"] == "---": # no event selected
                    fpExplorer_functions.plot_normalized_alone(self.canvas,
                                                self.options,
                                                self.normalized_dict[self.options["subject"]],
                                                self.normalized_export,
                                                self.save_plots,
                                                (self.export_path,self.export_begining),
                                                self.settings_dict)
                else:   # with event
                    custom_event_name = self.options["event"] if len(self.options["event_name"])==0 else self.options["event_name"]
                    custom_event_name2 = self.options["event2"] if len(self.options["event2_name"])==0 else self.options["event2_name"]
                    fpExplorer_functions.plot_normalized_alone_with_event(self.canvas,
                                                            self.options,
                                                            self.normalized_dict[self.options["subject"]],
                                                            custom_event_name,
                                                            custom_event_name2,
                                                            self.event_data,
                                                            self.normalized_export,
                                                            self.save_plots,
                                                            (self.export_path,self.export_begining),
                                                            self.settings_dict)
            if self.spikes_export == True or self.perievent_export == True:
                if self.perievent_export == True and self.spikes_export == False:   # show perievent settings window
                    self.open_spikes_later = False
                    # update path and add that to perievent options
                    self.perievent_options_dict["export_path"] = self.export_path
                    self.perievent_options_dict["file_beginning"] = self.export_begining
                    self.perievent_analysis_btn_clicked()  
                elif self.spikes_export == True and self.perievent_export == False:
                    self.peaks_btn_clicked()
                else: # both were selected
                    self.open_spikes_later = True # opens spikes window later
                    # open perievent window first
                    # update path and add that to perievent options
                    self.perievent_options_dict["export_path"] = self.export_path
                    self.perievent_options_dict["file_beginning"] = self.export_begining
                    self.perievent_analysis_btn_clicked()  
            else:
                self.reset_export_settings()
        
            self.close_export_window_signal.emit()
            self.enable_buttons_signal.emit()
        
    def peaks_btn_clicked(self):
        self.disable_buttons_signal.emit()
        self.peak_options_window = fpExplorer.PeakSettingsWindow(self.parent_window,
                                                      self.options["subject"],
                                                      self.recent_peak_values,
                                                      (self.export_path,self.export_begining),
                                                      self.batch_peaks)
        self.peak_options_window.got_peak_options_sig.connect(self.peak_options_received)
        self.peak_options_window.peak_window_closed_sig.connect(self.reset_open_later)
        self.peak_options_window.show()
        
    @pyqtSlot()   
    def reset_open_later(self):
        # in case user closed window that was set to open later
        self.open_spikes_later = False
      
    @pyqtSlot(list) 
    def peak_options_received(self,peak_options):
        # update group name
        self.options["subject_group_name"] = self.subject_group_name.text()
        # update value
        self.group_names_dict[self.options["subject"]] = self.options["subject_group_name"]
        self.recent_peak_values = peak_options[0]
        export_path,file_begin = peak_options[1]
        if self.batch_peaks == True:
            if "spikes" in self.parent_window.batch_export_settings_dict:
                if self.parent_window.batch_export_settings_dict["spikes"] == True:
                    if self.parent_window.batch_export_settings_dict["export_for_single_subjects"] == True:
                        for subject in self.parent_window.batch_export_settings_dict["batch_subjects"]:
                            if self.raw_data_dict[subject] != None:
                                # create subfolder with subject name
                                subject_subfolder = os.path.join(self.parent_window.batch_export_settings_dict["dump_path"],subject)
                                if not os.path.exists(subject_subfolder):
                                    try:
                                        os.mkdir(subject_subfolder)
                                    except:
                                        self.show_info_dialog("Problem creating subfolders")
                                        # save under main folder
                                        subject_subfolder = self.parent_window.batch_export_settings_dict["dump_path"]
                                if self.options["event"] == "---":
                                    try:
                                        fpExplorer_functions.plot_peaks(self.canvas,subject,self.normalized_dict[subject],
                                                            self.recent_peak_values,
                                                            self.spikes_export,
                                                            self.save_plots,
                                                            self.group_names_dict[subject],
                                                            self.settings_dict,
                                                            (subject_subfolder,self.parent_window.batch_export_settings_dict["file_begin"]))
                                    except:
                                        self.show_info_dialog("Could not calculate spikes.\nTry with different parameters.")
                                else:
                                    custom_event_name = self.options["event"] if len(self.options["event_name"])==0 else self.options["event_name"]
                                    custom_event_name2 = self.options["event2"] if len(self.options["event2_name"])==0 else self.options["event2_name"]
                                    try:
                                        fpExplorer_functions.plot_peaks_with_event(self.canvas,subject,self.normalized_dict[subject],
                                                                        self.recent_peak_values,
                                                                        custom_event_name,
                                                                        custom_event_name2,
                                                                        self.event_data,
                                                                        self.spikes_export,
                                                                        self.save_plots,
                                                                        self.group_names_dict[subject],
                                                                        self.settings_dict,
                                                                        (subject_subfolder,self.parent_window.batch_export_settings_dict["file_begin"]))
                                    except:
                                        self.show_info_dialog("Could not calculate spikes.\nTry with different parameters.")
                # close popup window
                self.peak_options_window.close()
                self.reset_export_settings()
                # reset batch peaks when done
                self.batch_peaks = False
                # reset open spikes later
                if self.open_spikes_later == True:
                    self.open_spikes_later = False
                # close run on batch window if still open
                if self.parent_window.run_on_batch_window != None:
                    self.parent_window.run_on_batch_window.close()
        else: # if it is not for batch analysis
            
            # if user entered path or if path was already there, update
            if len(export_path) > 0 or os.path.split(export_path)[1]==DEFAULT_EXPORT_FOLDER:
                try: 
                    os.mkdir(export_path)
                except:
                    if not os.path.exists(export_path):
                        self.show_info_dialog("Problem creating subfolder")
                   
                self.export_path = export_path
                self.export_begining = file_begin
            # if subject was not previewed yet, read data and add to dict
            if self.options["subject"] not in self.raw_data_dict:
                self.get_raw_data(self.options["subject"],
                              self.parent_window.preview_params[1][self.options["subject"]])
            # get trimming settings
            try:
                trim_beginning = int(self.trim_beginning_sec.text())
                trim_end = int(self.trim_ending_sec.text())
            except:
                trim_beginning = 0
                trim_end = 0
                self.show_info_dialog("Wrong trimming values.\nMake sure you entered valid integers")
                self.trim_beginning_sec.setText("0")
                self.trim_ending_sec.setText("0")
            # if it was not trimmed yet, add to dictionary
            if self.options["subject"] not in self.trimmed_raw_data_dict:
                # key is the subject name, value is a list
                # first element of that list is beginning and end seconds to trim
                # second element is trimmed data dict ts:timestamps, signal:data,control:data
                self.trimmed_raw_data_dict[self.options["subject"]] = [(trim_beginning,trim_end),self.trim_raw_data(self.raw_data_dict[self.options["subject"]],
                                                                                                                          trim_beginning,
                                                                                                                          trim_end
                                                                                                                          )]
            else:   # if already in trimmed check previous trimming settings
                trimmed = self.trimmed_raw_data_dict[self.options["subject"]]
                begin,end = trimmed[0]
                if trim_beginning != begin or trim_end != end:
                    # if new trimming params, replace the lists in dict
                    self.trimmed_raw_data_dict[self.options["subject"]] = [(trim_beginning,trim_end),self.trim_raw_data(self.raw_data_dict[self.options["subject"]],
                                                                                                                          trim_beginning,
                                                                                                                          trim_end
                                                                                                                          )]
            # always downsample first before normalizing
            # downsample
            self.downsampled_dict[self.options["subject"]] = fpExplorer_functions.downsample(self.trimmed_raw_data_dict[self.options["subject"]][1],
                                                                                         self.settings_dict[0]["downsample"])
            # check settings for method to normalize
            if self.settings_dict[0]["normalization"] == "Modified Polynomial Fitting":                   
                # normalize from downsampled
                self.normalized_dict[self.options["subject"]] = fpExplorer_functions.normalize_dff(
                                                                                        self.downsampled_dict[self.options["subject"]],
                                                                                        self.settings_dict[0]["show_norm_as"],
                                                                                        self.settings_dict[0]["filter"],
                                                                                        self.settings_dict[0]["filter_window"])
            if self.settings_dict[0]["normalization"] == "Standard Polynomial Fitting":                
                # normalize from downsampled
                self.normalized_dict[self.options["subject"]] = fpExplorer_functions.normalize_pMat(
                                                                                        self.downsampled_dict[self.options["subject"]],
                                                                                        self.settings_dict[0]["show_norm_as"],
                                                                                        self.settings_dict[0]["filter"],
                                                                                        self.settings_dict[0]["filter_window"])
            # set event
            self.options["event"] = self.event_from_data_comboBox.currentText()
            self.options["event_name"] = self.event_name_text.text()
            self.options["event2"] = self.event2_from_data_comboBox.currentText()
            self.options["event2_name"] = self.event2_name_text.text()
            # reset the event list to make sure it always has current events
            self.event_data = []
            if self.options["event"] != "---":
                evt = self.get_event_on_off(self.events_dict[self.options["subject"]], self.options["event"])
                print(evt)
                if len(evt[0]) == 0:
                    self.show_info_dialog("Some "+self.options["event"]+" event data is missing.")
                self.event_data.append(evt)
            if self.options["event2"] != "---":
                evt = self.get_event_on_off(self.events_dict[self.options["subject"]], self.options["event2"])
                if len(evt[0])==0:
                    self.show_info_dialog("Some "+self.options["event2"]+" event data is missing.")
                self.event_data.append(evt)
                
            
            wrong_event_order_info = "If you want to show just one event on the plots,\nselect your event as first event.\nAnd leave Event2 empty."
            # before you proceed with ploting check correct order of events  if only one event is selected
            if self.options["event"] == "---" and self.options["event2"] != "---":
                self.show_info_dialog(wrong_event_order_info)
              
            if self.options["event"] == "---":
                try:
                    fpExplorer_functions.plot_peaks(self.canvas,self.options["subject"],self.normalized_dict[self.options["subject"]],
                                         self.recent_peak_values,
                                         self.spikes_export,
                                         self.save_plots,
                                         self.group_names_dict[self.options["subject"]],
                                         self.settings_dict,
                                         (self.export_path,self.export_begining))
                except:
                    self.show_info_dialog("Could not calculate spikes.\nTry with different parameters.")
            else:
                custom_event_name = self.options["event"] if len(self.options["event_name"])==0 else self.options["event_name"]
                custom_event_name2 = self.options["event2"] if len(self.options["event2_name"])==0 else self.options["event2_name"]
                try:
                    fpExplorer_functions.plot_peaks_with_event(self.canvas,self.options["subject"],self.normalized_dict[self.options["subject"]],
                                                    self.recent_peak_values,
                                                    custom_event_name,
                                                    custom_event_name2,
                                                    self.event_data,
                                                    self.spikes_export,
                                                    self.save_plots,
                                                    self.group_names_dict[self.options["subject"]],
                                                    self.settings_dict,
                                                    (self.export_path,self.export_begining))
                except:
                    self.show_info_dialog("Could not calculate spikes.\nTry with different parameters.")
            # close popup window
            self.peak_options_window.close()
            self.reset_export_settings()
            if self.export_window != None:
                self.export_window.close()
                
#     @pyqtSlot()     
#     def start_batch_analysis(self):
#         # update group names
#         self.group_names_dict = self.parent_window.batch_export_settings_dict["updated_group_names"]
#         # update group name field of the current subject
#         self.subject_group_name.setText(self.group_names_dict[self.subject_comboBox.currentText()])

#         # update current trimming
#         # get trimming settings
#         # check what method was selected
#         if self.parent_window.batch_export_settings_dict["trim_begin_select"] == "Trim first seconds":
#             self.trim_beginning_sec.setText(self.parent_window.batch_export_settings_dict["trim_begin"])
#         else: # if trim by event
#             trim_begin_event = self.parent_window.batch_export_settings_dict["trim_begin_event"]
#             begin_evt_data_onsets = fpExplorer_functions.get_event_on_off(self.raw_data_dict[self.options["subject"]], trim_begin_event)[0]
#             trim_beginning = 0 # assign zero in case there is no such event
#             if len(begin_evt_data_onsets) > 0:
#                 # use the first onset as trim begin
#                 trim_beginning = math.ceil(begin_evt_data_onsets[0])
#             else:
#                 self.show_info_dialog("There was a problem reading event "+trim_begin_event+"\nfor subject "+self.options["subject"]+"\nThis data was not trimmed.")
#             self.trim_beginning_sec.setText(str(trim_beginning))
#         if self.parent_window.batch_export_settings_dict["trim_end_select"] == "Trim last seconds":
#             self.trim_ending_sec.setText(self.parent_window.batch_export_settings_dict["trim_end"])
#         else: # if trim by event
#             trim_end_event = self.parent_window.batch_export_settings_dict["trim_end_event"]
#             end_evt_data_onsets = fpExplorer_functions.get_event_on_off(self.raw_data_dict[self.options["subject"]], trim_end_event)[0]
#             trim_end = 0 # assign zero in case there is no such event
#             if len(end_evt_data_onsets) > 0:
#                 last_raw_ts = fpExplorer_functions.get_last_timestamp(
#                                      self.raw_data_dict[self.options["subject"]],
#                                      self.preview_init_params[0][0]["signal_name"],
#                                      )
#                 # use last onset as trim begin
#                 trim_end = math.ceil(last_raw_ts - end_evt_data_onsets[-1])
#             else:
#                 self.show_info_dialog("There was a problem reading event "+trim_end_event+"\nfor subject "+self.options["subject"]+"\nThis data was not trimmed.")
#             self.trim_ending_sec.setText(str(trim_end))


#         # set include to actual status
#         if self.subject_comboBox.currentText() in self.parent_window.selected_subjects:
#             self.include_cb.setChecked(True)
#         else:
#             self.include_cb.setChecked(False)
#         # disable buttons while processing
#         self.disable_buttons_signal.emit()
#         print("Started batch processing")
#         # new_trim_start = int(self.parent_window.batch_export_settings_dict["trim_begin"])
#         # new_trim_end = int(self.parent_window.batch_export_settings_dict["trim_end"])
#         new_trim_start = int(self.trim_beginning_sec.text())
#         new_trim_end = int(self.trim_ending_sec.text())
        
#         # set event
#         self.options["event"] = self.event_from_data_comboBox.currentText()
#         self.options["event_name"] = self.event_name_text.text()
#         self.options["event2"] = self.event2_from_data_comboBox.currentText()
#         self.options["event2_name"] = self.event2_name_text.text()
#         # reset the event list to make sure it always has current events
#         self.event_data = []
#         if self.options["event"] != "---":
#             evt = fpExplorer_functions.get_event_on_off(self.raw_data_dict[self.options["subject"]], self.options["event"])
#             if len(evt[0]) == 0:
#                 self.show_info_dialog("Some "+self.options["event"]+" event data is missing.")
#             self.event_data.append(evt)
#         if self.options["event2"] != "---":
#             evt = fpExplorer_functions.get_event_on_off(self.raw_data_dict[self.options["subject"]], self.options["event2"])
#             if len(evt[0])==0:
#                 self.show_info_dialog("Some "+self.options["event2"]+" event data is missing.")
#             self.event_data.append(evt)
                   
#         wrong_event_order_info = "If you want to show just one event on the plots,\nselect your event as first event.\nAnd leave Event2 empty."
#         # before you proceed with ploting check correct order of events  if only one event is selected
#         if self.options["event"] == "---" and self.options["event2"] != "---":
#             self.show_info_dialog(wrong_event_order_info)
            
#         # create a popup window that will be on during processing?
#         # disable all buttons?
#         if self.parent_window.batch_export_settings_dict["normalized"] == True or self.parent_window.batch_export_settings_dict["spikes"] == True:
#             # create normalized data for most recent settings
#             for i in range(len(self.parent_window.batch_export_settings_dict["batch_subjects"])):
#                 subject = self.parent_window.batch_export_settings_dict["batch_subjects"][i]
#                 # create separate options dictionary for batch analysis
#                 self.batch_options_dict = {"subject":subject,
#                                             "subject_group_name":self.parent_window.batch_export_settings_dict["batch_subjects_group_names"][i]}
#                 # add the group name to group names dictionary if it was not there or update group names from batch options
#                 self.group_names_dict[subject] = self.parent_window.batch_export_settings_dict["batch_subjects_group_names"][i]
#                 # if subject was not previewed yet, read data and add to dict
#                 if subject not in self.raw_data_dict:
#                     self.get_raw_data(subject, self.parent_window.preview_params[1][subject])
#                 if  self.raw_data_dict[subject] != None:
#                     # trim to batch export settings
#                     # key is the subject name, value is a list
#                     # first element of that list is beginning and end seconds to trim
#                     # second element is trimmed data dict ts:timestamps, signal:data,control:data
#                     self.trimmed_raw_data_dict[subject] = [(new_trim_start,new_trim_end),fpExplorer_functions.trim_raw_data(self.raw_data_dict[subject],
#                                                                                     self.preview_init_params[0][0]["signal_name"],
#                                                                                     self.preview_init_params[0][0]["control_name"],
#                                                                                     new_trim_start,
#                                                                                     new_trim_end
#                                                                                     )]
                    
#                     # always downsample first before normalizing
#                     # downsample
#                     self.downsampled_dict[subject] = fpExplorer_functions.downsample(self.trimmed_raw_data_dict[subject][1],
#                                                                                 self.settings_dict[0]["downsample"])
#                     # check settings for method to normalize
#                     if self.settings_dict[0]["normalization"] == "Modified Polynomial Fitting":                   
#                         # normalize from downsampled
#                         self.normalized_dict[subject] = fpExplorer_functions.normalize_dff(
#                                                                                                 self.downsampled_dict[subject],
#                                                                                                 self.settings_dict[0]["show_norm_as"],
#                                                                                                 self.settings_dict[0]["filter"],
#                                                                                                 self.settings_dict[0]["filter_window"])
#                     if self.settings_dict[0]["normalization"] == "Standard Polynomial Fitting":                
#                         # normalize from downsampled
#                         self.normalized_dict[subject] = fpExplorer_functions.normalize_pMat(
#                                                                                                 self.downsampled_dict[subject],
#                                                                                                 self.settings_dict[0]["show_norm_as"],
#                                                                                                 self.settings_dict[0]["filter"],
#                                                                                                 self.settings_dict[0]["filter_window"])
#                     if self.parent_window.batch_export_settings_dict["export_for_single_subjects"] == True:
#                         # create subfolder with subject name
#                         subject_subfolder = os.path.join(self.parent_window.batch_export_settings_dict["dump_path"],subject)
#                         if not os.path.exists(subject_subfolder):
#                             try:
#                                 os.mkdir(subject_subfolder)
#                             except:
#                                 self.show_info_dialog("Problem creating subfolder")
#                                 # save under main folder
#                                 subject_subfolder = self.parent_window.batch_export_settings_dict["dump_path"]
#                         # save also polynomial fitting by default
#                         fpExplorer_functions.show_polynomial_fitting(self.canvas,
#                                                 self.settings_dict[0], 
#                                                 self.downsampled_dict[subject],
#                                                 self.preview_init_params[0][0]["signal_name"],
#                                                 self.preview_init_params[0][0]["control_name"],
#                                                 subject,
#                                                 True,
#                                                 subject_subfolder)
#                         if self.parent_window.batch_export_settings_dict["normalized"] == True:
#                             if self.options["event"] == "---":
#                                 fpExplorer_functions.plot_normalized_alone(self.canvas,
#                                                         self.batch_options_dict,
#                                                         self.normalized_dict[subject],
#                                                         True,
#                                                         True,
#                                                         (subject_subfolder,self.parent_window.batch_export_settings_dict["file_begin"]),
#                                                         self.settings_dict)
#                             else:
#                                 custom_event_name = self.options["event"] if len(self.options["event_name"])==0 else self.options["event_name"]
#                                 custom_event_name2 = self.options["event2"] if len(self.options["event2_name"])==0 else self.options["event2_name"]
#                                 fpExplorer_functions.plot_normalized_alone_with_event(self.canvas,
#                                                             self.batch_options_dict,
#                                                             self.normalized_dict[subject],
#                                                             custom_event_name,
#                                                             custom_event_name2,
#                                                             self.event_data,
#                                                             True,
#                                                             True,
#                                                             (subject_subfolder,self.parent_window.batch_export_settings_dict["file_begin"]),
#                                                             self.settings_dict)
#                         if self.parent_window.batch_export_settings_dict["spikes"] == True:
#                             self.spikes_export = True
#                             self.save_plots = True
#             ### end for loop for each subject   
#             if self.parent_window.batch_export_settings_dict["export_group_data"] == True:
#                 if self.parent_window.batch_export_settings_dict["spikes"] == True:
#                     # set batch peaks to true
#                     self.batch_peaks = True
#             if self.parent_window.batch_export_settings_dict["spikes"] == True:
#                 # set batch peaks to true
#                 self.batch_peaks = True
#                 # remember to open it later
#                 self.open_spikes_later = True                  
#         if self.parent_window.batch_export_settings_dict["perievent"] == True:
#             self.batch_perievent = True
#             self.save_plots = True
#             self.perievent_analysis_btn_clicked()
#         else:
#             if self.open_spikes_later == False:
#                 print("Done batch processing")
#                 self.done_batch_processing_sig.emit()
#                 self.enable_buttons_signal.emit()
#             elif self.open_spikes_later == True:
#                 self.peaks_btn_clicked()
            
            
    def get_frequency(self,subject):
        data_df = self.raw_data_dict[subject]
        ts = data_df.iloc[:,0].to_numpy()
        interval_float = ts[1]-ts[0]
        fs = 1/interval_float
        return fs

    def read_data_csv(self,path_to_data):
        df = pd.read_csv(path_to_data)
        return df

    def get_events(self,df):
        all_events = df.iloc[:,0].unique()
        events = []
        for evt in all_events:
            # exclude cam or tickfrom events
            if evt.lower().startswith("cam") == False and evt.lower().startswith("tick") == False:
                events.append(evt)
        return events

    def get_event_on_off(self,event_dataframe, event):
        ''' Returns onset and offset times of the tones
            First element of this array are onset times, 
            Second el are offset times
        '''
        on_off = [[],[]]
        # first column should be Event name, second event onset, last event offset
        column_names = list(event_dataframe.columns)
        # select only rows with given event name
        selected_event_df = event_dataframe.loc[event_dataframe[column_names[0]] == event]
        onsets = selected_event_df[column_names[1]].tolist()
        # validate (could be i.e inf)
        valid_onsets = []
        for el in onsets:
            try:
                val = float(el)
                valid_onsets.append(val)
            except:
                pass
        offsets = selected_event_df[column_names[2]].tolist()
        # validate
        valid_offsets = []
        for el in offsets:
            try:
                val = float(el)
                valid_offsets.append(val)
            except:
                pass
        if len(valid_onsets) > 0:
            on_off[0] = valid_onsets
        else:
            print("There are no onsets for the event",event)
        if len(valid_offsets) > 0:
            on_off[1] = valid_offsets
        else:
            print("There are no offsets for the event",event)
        return on_off

    def filter_data_around_event(self,raw_data,events_df,perievent_options_dict,settings_dict):
        filtered = {}
        event_name = perievent_options_dict["event"]
        before = perievent_options_dict["sec_before"]
        till = perievent_options_dict["sec_after"]+0.1
        evt_onsets = self.get_event_on_off(events_df, event_name)[0]
        signal_data = raw_data.iloc[:,1].to_numpy()
        control_data = raw_data.iloc[:,2].to_numpy()
        # timestamps in data do not match exactly timestamps in event file
        evt_onsets_from_data = []
        raw_ts = raw_data.iloc[:,0].to_numpy()
        # round raw to 3 decimal places just in case the original time is not found
        rounded_arr = np.around(raw_ts,decimals=3)
        for el in evt_onsets:
            if el in raw_ts:
                evt_onsets_from_data.append(el)
            else:              
                rounded_onset = round(el,3)
                if rounded_onset in rounded_arr:
                    # pick only the first one
                    result_idx = np.where(rounded_arr == rounded_onset)[0][0]
                    evt_onsets_from_data.append(raw_ts[result_idx])
        if len(evt_onsets_from_data) == 0:
            print("Could not find matching event onset times times")
            return filtered
        
        # chop data into windows around the events
        signals = []
        controls = []
        for evt in evt_onsets_from_data:
            condition = np.greater_equal(raw_ts,evt-before)&np.less(raw_ts,evt+till)
            signal_temp = np.extract(condition,signal_data)
            control_temp = np.extract(condition,control_data)
            signals.append(signal_temp)
            controls.append(control_temp)
        
        # downsample data as well
        N = settings_dict[0]["downsample"] 
        
        control_downsampled = []
        signal_downsampled = []
        # get frequency
        interval_float = raw_ts[1]-raw_ts[0]
        fs = 1/interval_float
        try:
            for lst in controls:
                ts = np.linspace(1, len(lst), len(lst))/fs
                # round last time to full integer
                last = math.floor(ts[-1])
                # make new times till total expected
                ts_adjusted = np.linspace(1/N,last,last*N)
                control_interpolated = interp1d(ts, lst)
                resampled_data = control_interpolated(ts_adjusted)
                control_downsampled.append(resampled_data)
            for lst in signals:
                ts = np.linspace(1, len(lst), len(lst))/fs
                # round last time to full integer
                last = math.floor(ts[-1])
                ts_adjusted = np.linspace(1/N,last,last*N)
                signal_interpolated = interp1d(ts, lst)
                resampled_data = signal_interpolated(ts_adjusted)
                signal_downsampled.append(resampled_data)
            filtered["signal"] = signal_downsampled
            filtered["control"] = control_downsampled   
        except:
            print("Not enough data that satisfy your request. Try changing event or times around event.")
        
        return filtered

    def trim_raw_data(self,data_df,beginning_sec,ending_sec):
        # first column should be time, second signal, last control
        # return numpy array for all GCaMP channel raw data
        signal_data = data_df.iloc[:,1].to_numpy()
        # return numpy arrayfor all control raw data
        control_data = data_df.iloc[:,2].to_numpy()
        ts_raw = data_df.iloc[:,0].to_numpy()
        
        # calculate the first beginning_sec seconds from a signal channel
        t0 = int(beginning_sec * self.current_fs) # int rounds it to the nearest integer
        t1 = int(ending_sec * self.current_fs)
    #    print(t0,t1)
        if t1 == 0:
            ts = ts_raw[t0:]
            signal_trimmed = signal_data[t0:]
            control_trimmed = control_data[t0:]
        else:
            ts = ts_raw[t0:-t1]
            signal_trimmed = signal_data[t0:-t1]
            control_trimmed = control_data[t0:-t1]
    #    print(len(signal_trimmed),len(control_trimmed))
        # create a dictionary
        trimmed_dict = {"ts":ts,"signal":signal_trimmed,"control":control_trimmed}
        return trimmed_dict

    # create settings dataframe for csv
    #  headers   
    def get_settings_df(self,settings_dict):
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

    # show info popup
    def show_info_dialog(self, text):
        msgBox = QMessageBox()
        msgBox.setWindowIcon(QtGui.QIcon(ICO))
        msgBox.setText(text)
        msgBox.setWindowTitle("Info!")
        msgBox.setStandardButtons(QMessageBox.Ok)
        msgBox.exec()

#######################################################################################################################
# PLOTS
    # when app starts plot just raw data  
    def plot_raw(self,canvas,subject,export,save_plots,export_loc_data):
        # a tuple path+file beginning
        dump_path,file_beginning = export_loc_data
        data_df = self.raw_data_dict[subject]
        # return numpy array for all GCaMP channel raw data
        signal_data = data_df.iloc[:,1].to_numpy()
        # return numpy arrayfor all control raw data
        control_data = data_df.iloc[:,2].to_numpy()
        # first column should be time, second signal, last control
        column_names = list(data_df.columns)
        
        ###############################
        # create new timestamps for the data (from zero)
        # ts = create_timestamps(raw_data,signal_name,GCaMP_data)
        times_from_file = data_df.iloc[:,0].to_numpy()
        num_samples = len(times_from_file)
        ts = np.linspace(1, num_samples, num_samples) / self.current_fs

        # plot
        ax = canvas.fig.add_subplot(111)
        ax.cla()  # Clear the canvas
        # plot the lines
        ax.plot(ts, signal_data,
                color=SIGNAL_COLOR_RGB,
                linewidth=1,
                label = column_names[1])
        ax.plot(ts,control_data,
                color=CONTROL_COLOR_RGB,
                linewidth=1,
                label = column_names[2])

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
            sig = column_names[1]
            cont = column_names[2]
            raw_df= pd.DataFrame({"Time (sec)":ts,
                                    sig +' (mV)' : signal_data, 
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

    def plot_trimmed(self,canvas,subject,signal_dict,export,save_plots,export_loc_data):
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
            sig = "Signal"
            cont = "Control"
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

    def plot_with_event(self,canvas,subject,signal_dict,event_name,event2_name,event_data,export,save_plots,export_loc_data):
        # a tuple path+file beginning
        dump_path,file_beginning = export_loc_data
        
        GCaMP_data = signal_dict["signal"]
        control_data = signal_dict["control"]
        ts = signal_dict["ts"]
            
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
            sig = "Signal"
            cont = "Control"
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

    def plot_downsampled_alone(self,canvas,options_dict,downsampled_signal_dict,export,save_plots,export_loc_data,settings_dict):
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
            settings_df = self.get_settings_df(settings_dict)
            sig = "Signal"
            cont = "Control"
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

    def plot_downsampled_alone_with_event(self,canvas,options_dict,downsampled_signal_dict,event_name,event2_name,event_data,export,save_plots,export_loc_data,settings_dict):
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
            settings_df = self.get_settings_df(settings_dict)
            sig = "Signal"
            cont = "Control"
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

    def show_polynomial_fitting(self,canvas, settings_dict,downsampled,subject_name,save_plot,dump_path):
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
        # fit time axis to the 465nm stream 
        bls_Ca = np.polynomial.polynomial.Polynomial.fit(ts_reset,signal_arr,1)
        # fit time axis the 405nm stream
        bls_ref = np.polynomial.polynomial.Polynomial.fit(ts_reset,control_arr,1)
        # the below returns first: slope, second: intercept
        print("bls_Ca",bls_Ca.convert().coef[::-1])
        # the below returns first: slope, second: intercept
        print("bls_ref",bls_ref.convert().coef[::-1])
        # put those values in a dictionary
        slope_intercept_dict = {"signal_slope_intercept":bls_Ca.convert().coef[::-1],
                                "control_slope_intercept":bls_ref.convert().coef[::-1]
        }

        if normalization == 'Standard Polynomial Fitting':                   
            # https://stackoverflow.com/questions/45338872/matlab-polyval-function-with-three-outputs-equivalent-in-python-numpy
            mu = np.mean(control_arr)
            std = np.std(control_arr, ddof=0)
            # Call np.polyfit(), using the shifted and scaled version of control_arr
            cscaled = np.polynomial.polynomial.Polynomial.fit((control_arr - mu)/std, signal_arr, 1)
            # Create a poly1d object that can be called
            # https://numpy.org/doc/stable/reference/routines.polynomials.html
            pscaled = Polynomial(cscaled.convert().coef)
            # Inputs to pscaled must be shifted and scaled using mu and std
            F0 = pscaled((control_arr - mu)/std)    
            # plot
            # clear previous figure
            canvas.fig.clf()
            ax = canvas.fig.add_subplot(111)
            
            ax.plot(ts_reset, signal_arr, linewidth=1, color=SIGNAL_COLOR_RGB, label='F signal')
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
            # plot
            # clear previous figure
            canvas.fig.clf()
            ax = canvas.fig.add_subplot(211)
            ax2 = canvas.fig.add_subplot(212)
            
            ax.plot(ts_reset, signal_arr, linewidth=1, color=SIGNAL_COLOR_RGB, label='F signal')
            ax.plot(ts_reset, F0Ca, linewidth=1, color='k', label='F0')
            
            ax2.plot(ts_reset, control_arr, linewidth=1, color=CONTROL_COLOR_RGB, label='F control')
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

    # plot raw but downsampled perievents
    def plot_raw_perievents(self,canvas,subject,modified_data,current_trials,perievent_options_dict,settings_dict,export,save_plots,group_name,export_loc_data):
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
        # get downsampled data
        GCaMP_perievent_data = modified_data["signal"]
        control_perievent_data = modified_data["control"]
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
        for i in range(total_plots):
            header_signal = "Trial"+str(i+1)+norm_type
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
            settings_df = self.get_settings_df(settings_dict)
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
        

########################### end PreviewEventBasedWidget class   

        
            