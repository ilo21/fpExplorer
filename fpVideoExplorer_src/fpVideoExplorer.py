"""
 Copyright (c) 2023 CSAN_LiU

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
Created on Thu Jan  7 12:08:47 2021

@author: ilosz01
"""

from tkinter import E
import fpVideoExplorer_functions
import os
import shutil
import math
from PyQt5 import QtGui 
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.QtCore import Qt, pyqtSignal, pyqtSlot, QTimer
import pyqtgraph as pg
from pyqtgraph.dockarea import *
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import types
import pandas as pd
import numpy as np
import cv2
import sys
if os.name == "posix": # ubuntu: qt plugin internally used by opencv is not compatible with pyqt5
    if sys.platform != "darwin": # if not mac
        # unset environment variable after importing cv2
        os.environ.pop("QT_QPA_PLATFORM_PLUGIN_PATH")
import copy
import resources
import warnings
warnings.filterwarnings("ignore")

# App icon
ICO = ':/icons/app_icon'
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
# how low in downsampling can go
MAX_DOWNSAMPLE_PCT = 50
# how high in downsampling can you go
MIN_DOWNSAMPLE_PCT = 0.5
# default_samples to average 1% of original sampling rate
DEFAULT_DOWNSAMPLE_PCT = 10
# change to hardcoded default of 100 Hz
DEFAULT_HZ = 100
MAX_SMOOTH_WINDOW = 10000
DEFAULT_SMOOTH_WINDOW = 10
DEFAULT_EXPORT_FOLDER = "_fpVideoExplorerAnalysis"
SHIFT_MAX = 100
VIDEO_DELAY_SEC = 0.2 # the camera starts recording later than fp data collection
###############################
# CLASS FORM MAIN GUI WINDOW  #
###############################
            
class MyMainWidget(QMainWindow):
    # signal that select folder window is open
    select_data_open_signal = pyqtSignal()
    got_channel_names = pyqtSignal()
    # emit signal to close all other windows
    app_closed = pyqtSignal() 
    def __init__(self):
        super(MyMainWidget, self).__init__()
        self.name = "fpVideoExplorer"
        self.app_width = 800
        self.app_height = 200
        self.top_buttons_height = 50
        self.setWindowIcon(QtGui.QIcon(ICO))
        self.setWindowTitle(self.name)
        self.resize(self.app_width,self.app_height)

        self.select_data_window_content = []
        self.select_data_window = None
        # window with general settings
        self.settings_window = None
        # second element is a dict of paths to all data folders
        self.preview_params = []
        # widget with preview
        self.preview_widget = None
        # list with general settings dictionary as first element 
        self.settings_dict = [{"downsample":None,
                              "entered_downsample":None,
                              "normalization": "Standard Polynomial Fitting",
                              "show_norm_as":"dF/F (in % )",
                              "filter":False,
                              "filter_window":DEFAULT_SMOOTH_WINDOW,
                              "subject":""}]
        # window with general settings
        self.settings_window = None

        # use docks from pyqtgraph
        self.area = DockArea()
        self.setCentralWidget(self.area)
        self.top_dock_widget = Dock("Dock1", size=(self.app_width, self.top_buttons_height))
        self.top_dock_widget.hideTitleBar()
        self.area.addDock(self.top_dock_widget,'left')
        self.bottom_dock_widget = Dock("Dock2", size=(self.app_width,self.app_height-self.top_buttons_height))
        self.bottom_dock_widget.hideTitleBar()
        self.area.addDock(self.bottom_dock_widget, 'bottom', self.top_dock_widget)  # place the bottom dock at bottom edge of top dock
        
        # define fields
        self.select_data_btn = QPushButton('Select TDT Data')
        self.settings_btn = QPushButton('Settings')
        self.settings_btn.setEnabled(False)

        # Create a layout to manage the buttons size and position
        self.box_layout = QHBoxLayout()
        # set the spacing around the layout
        self.box_layout.setContentsMargins(10,10,10,10)
        # place widgets to the layout in their proper positions
        self.box_layout.addWidget(self.select_data_btn)
        self.box_layout.addWidget(self.settings_btn)
        # create a widget for dock from pyqtgraph
        self.top_widget = QWidget()
        # add that widget to the dock
        self.top_widget.setLayout(self.box_layout)
        self.top_dock_widget.addWidget(self.top_widget)
        
        
        # connect buttons to fpVideoExplorer_functions
        self.select_data_btn.clicked.connect(self.select_data_btn_clicked)
        self.settings_btn.clicked.connect(self.settings_btn_clicked)
        self.got_channel_names.connect(self.show_preview)
        # use stylesheet
        self.area.setStyleSheet(STYLESHEET)

    # define keypress events
    def keyPressEvent(self,event):
        # if enter is pressed start button clicked
        if event.key() == Qt.Key_Return and self.select_data_btn.isEnabled:
            self.select_data_btn_clicked()

    def select_data_btn_clicked(self):
        self.disable_buttons()
        self.select_data_window = SelectDataWindow(self,self.select_data_window_content)
        self.select_data_window.got_user_input_sig.connect(self.get_select_data_window_user_input)
        self.select_data_window.close_select_data_sig.connect(self.got_select_data_window_closed)
        self.select_data_window.show()
        # emit signal that the window is open
        self.select_data_open_signal.emit()

    def settings_btn_clicked(self):
        self.disable_buttons()
        QApplication.processEvents() 
        self.settings_window = SettingsWindow(self,self.settings_dict)
        self.settings_window.got_user_settings_sig.connect(self.get_settings)
        self.settings_window.show()

    def enable_all_buttons(self):
        self.select_data_btn.setEnabled(True)
        if (len(self.select_data_window_content)>0 or len(self.select_custom_data_window_content)>0):
            self.settings_btn.setEnabled(True)
        QApplication.processEvents() 

    def disable_buttons(self):
        self.select_data_btn.setEnabled(False)
        self.settings_btn.setEnabled(False)
        QApplication.processEvents() 

    @pyqtSlot()
    # receives signal that select folder window is closed
    # enables buttons 
    def got_select_data_window_closed(self):
        self.enable_all_buttons()

    @pyqtSlot(list) 
    # receives a list with options data read from select data window
    def get_select_data_window_user_input(self,user_input):

        self.select_data_window_content = user_input       
        if self.select_data_window_content[0]["subject_experiment"] == True:
            self.batch_paths_dict = fpVideoExplorer_functions.create_list_of_paths(self.select_data_window_content[0]["main_path"],
                                                                   self.select_data_window_content[0]["subject_names"],
                                                                   self.select_data_window_content[0]["selected_experiment"])
            # if path doesn't exist remove it
            keys2remove = []
            remaining = []
            for key in self.batch_paths_dict:
                if not os.path.exists(self.batch_paths_dict[key]):
                    keys2remove.append(key)
                else:
                    remaining.append(key)
            for key in keys2remove:
                del self.batch_paths_dict[key]
            self.select_data_window_content[0]["subject_names"] = remaining
        elif self.select_data_window_content[0]["experiment_subject"] == True:
            valid_subjects,self.batch_paths_dict = fpVideoExplorer_functions.create_list_of_paths_experiment_subjects(self.select_data_window_content[0]["main_path"],
                                                                    self.select_data_window_content[0]["subject_names"],
                                                                    self.select_data_window_content[0]["selected_experiment"])
            # update subject list for that experiment
            self.select_data_window_content[0]["subject_names"] = valid_subjects
            
       
        # ask user for signal and control channels
        data = fpVideoExplorer_functions.get_raw_data(self.batch_paths_dict[self.select_data_window_content[0]["subject_names"][0]])
        if data != None:
            self.signal_control_window = SignalControlWindow(self,fpVideoExplorer_functions.get_channel_names(data))
            self.signal_control_window.got_signal_name_sig.connect(self.got_signal_name_sig)
            self.signal_control_window.show()   
        else:
            self.show_info_dialog("Could not read subject "+self.select_data_window_content[0]["subject_names"][0]+" file.\nTry excluding it from your data folder.") 

    @pyqtSlot(list) 
    # receives a list with dictionary of signal and control channel names
    def got_signal_name_sig(self,channel_names):
        if len(channel_names) > 0:
            self.select_data_window_content[0]["signal_name"] = channel_names[0]["signal_name"]
            self.select_data_window_content[0]["control_name"] = channel_names[0]["control_name"]
            print(self.select_data_window_content)
            # clear also previously selected subjects
            self.selected_subjects = []
            self.got_channel_names.emit()
        else:
            if self.settings_window != None:
                self.settings_window.close()
       
    @pyqtSlot()   
    def show_preview(self):
        self.select_data_window.close()
        self.select_custom_data_window_content = []
        if self.preview_widget != None: # clear preview before creating next one
            self.preview_widget.close()    
            self.preview_widget.deleteLater()
            self.preview_widget = None
        # create list of init_params to start previewing
        self.preview_params = [self.select_data_window_content,self.batch_paths_dict]
        print("\npreview_params from main\n",self.preview_params)
        
        # automatically add preview to the main window
        # if the event based experiment type was selected set preview to PreviewEventBasedWidget
        # if self.select_data_window_content[0]["event_based"] == True:
            # try reading the first data tank on the list 
            # to verify if there are indeed events
        if fpVideoExplorer_functions.check_events(self.batch_paths_dict[self.select_data_window_content[0]["subject_names"][0]]):
            self.enable_all_buttons()
            # show window with options to plot perievent
            self.preview_widget = PreviewEventBasedWidget(self,self.preview_params)
            # add widget to dock
            self.bottom_dock_widget.addWidget(self.preview_widget)
        else:
            self.show_info_dialog("No events found in your data.")

    # add raw data structure and subject to self.raw_data dictionary   
    def get_raw_data(self,subject,my_path):
        self.raw_data_dict[subject] = fpVideoExplorer_functions.get_raw_data(my_path)
        if self.raw_data_dict[subject] == None:
            # remove from combo box
            index = self.subject_comboBox.findText(subject)  # find the index of text
            self.subject_comboBox.removeItem(index)  # remove item from index
            # and remove from available subjects
            self.preview_init_params[0][0]["subject_names"].remove(subject)
            self.show_info_dialog("Problem reading subject's "+subject+" file.\nSubject will not be available for the analysis.")

    @pyqtSlot(list) 
    # receives a single element list with settings from settings window
    def get_settings(self,settings):
        self.settings_dict[0]["downsample"] = settings[0]["downsample"] 
        self.settings_dict[0]["normalization"] = settings[0]["normalization"]
        self.settings_dict[0]["show_norm_as"] = settings[0]["show_norm_as"]
        self.settings_dict[0]["filter"] = settings[0]["filter"]
        print("settings:",self.settings_dict)
        self.enable_all_buttons()

    # show info popup
    def show_info_dialog(self, text):
        msgBox = QMessageBox()
        msgBox.setWindowIcon(QtGui.QIcon(ICO))
        msgBox.setText(text)
        msgBox.setWindowTitle("Info!")
        msgBox.setStandardButtons(QMessageBox.Ok)
        msgBox.exec()

    # # emit signal when app is closed
    # def closeEvent(self, event):  
    #     self.app_closed.emit() 

####### end MyMainWidget class

#################################
# CLASS FOR SELECT DATA WINDOW  #
#################################
        
class SettingsWindow(QMainWindow):
    # pyqt Signal has to be defined on class level
    # and not in init method !!
    # create a list with information acquired from window
    # and send it as a signal
    got_user_settings_sig = pyqtSignal(list)
    def __init__(self,parent_window,settings):
        super(SettingsWindow, self).__init__()
        self.setWindowTitle("Settings")
        self.setWindowIcon(QtGui.QIcon(ICO))
        self.resize(650,250)
        self.parent_window = parent_window
        self.parent_window.app_closed.connect(self.exit_app)
        self.settings = settings
        self.current_fs = 0
        if self.parent_window.preview_widget != 0:
            self.current_fs = self.parent_window.preview_widget.current_fs
            # current_subject = self.parent_window.preview_widget.options["subject"]
            # self.current_fs = fpVideoExplorer_functions.get_frequency(self.parent_window.preview_widget.raw_data_dict[current_subject],self.parent_window.preview_widget.preview_init_params[0][0]["signal_name"])
            print("Frequency:",self.current_fs)

        # calculate min, max and sugested sample rates
        # self.suggested_downsample_samples = round(self.current_fs*DEFAULT_DOWNSAMPLE_PCT/100)
        # self.suggested_downsample_rate = round(self.current_fs/self.suggested_downsample_samples)
        self.min_downsampe_rate = round(self.current_fs*MIN_DOWNSAMPLE_PCT/100)
        self.max_downsample_rate = self.current_fs
        # self.min_samples = round(self.current_fs/self.max_downsample_rate)
        # self.max_samples = round(self.current_fs/self.min_downsampe_rate)
        
        # create widget with settings
        self.settings_main_widget = QWidget()
        self.setCentralWidget(self.settings_main_widget)
        # create main layout for widget with settings
        self.settings_main_layout = QVBoxLayout()
        self.settings_main_layout.setContentsMargins(10,10,10,10)
        self.settings_layout = QFormLayout()
        self.settings_layout.setContentsMargins(10,10,10,10)
        self.settings_layout.setVerticalSpacing(15)
        # self.downsample_text = QComboBox()
        self.downsample_text = QLineEdit(str(self.settings[0]["downsample"]))
        self.downsample_text.setValidator(QtGui.QIntValidator())
        self.downsample_text.setToolTip("Integers from "+str(self.min_downsampe_rate)+" to "+str(round(self.max_downsample_rate)))
        self.downsample_label = QLabel("Downsample to (in Hz)\nOriginal rate: "+str(round(self.current_fs))+" Hz")
        self.settings_layout.addRow(self.downsample_label,self.downsample_text)
        # self.update_rate_btn = QPushButton("Preview new rate")
        # self.after_rate_label = QLabel("After downsampling: "+str(round(self.current_fs/int(self.downsample_text.text()))) +" Hz")
        # self.settings_layout.addRow(self.after_rate_label,self.update_rate_btn)
        self.normalization_method_comboBox = QComboBox()
        self.normalization_method_comboBox.addItems(["Standard Polynomial Fitting","Modified Polynomial Fitting"])
        self.normalization_method_comboBox.setCurrentText(self.settings[0]["normalization"])
        self.settings_layout.addRow("Method of normalization",self.normalization_method_comboBox)
        self.normalization_show_comboBox = QComboBox()
        self.normalization_show_comboBox.addItems(["df/F ( in % )","Z-Score"])
        self.normalization_show_comboBox.setCurrentText(self.settings[0]["show_norm_as"])
        self.settings_layout.addRow("Show normalized data as",self.normalization_show_comboBox)
        self.filter_cb = QCheckBox("Filter")
        self.filter_cb.setChecked(self.settings[0]["filter"])
        self.settings_layout.addRow("Smooth data:",QLabel("The window around each sample (in samples):"))
        self.smooth_text = QLineEdit(str(self.settings[0]["filter_window"]))
        self.smooth_text.setToolTip("From 0 to "+str(MAX_SMOOTH_WINDOW))
        self.smooth_text.setValidator(QtGui.QIntValidator())
        self.settings_layout.addRow(self.filter_cb,self.smooth_text)
        self.settings_main_layout.addLayout(self.settings_layout)
        self.save_settings_layout = QVBoxLayout()
        self.save_settings_layout.setAlignment(Qt.AlignRight)
        self.save_settings_btn = QPushButton("Save settings")
        self.save_settings_layout.addWidget(self.save_settings_btn)
        self.settings_main_layout.addLayout(self.save_settings_layout)
        
        # set the window's main layout
        self.settings_main_widget.setLayout(self.settings_main_layout)
        # # make the hint button with updated sampling rate smaller than regular buttons
        # self.update_rate_btn.setStyleSheet("QPushButton {padding: 2px; margin-left:0px; margin-right:0px; margin-top:0px; margin-bottom:0px; font-size: 10pt; font-weight: normal;}")
        
        # use main stylesheet
        self.setStyleSheet(STYLESHEET)
        
        self.save_settings_btn.clicked.connect(self.save_settings_btn_clicked)
        # self.update_rate_btn.clicked.connect(self.update_rate)

    def save_settings_btn_clicked(self):
        # read user settings
        self.read_user_settings()
        # send values to main app window
        self.got_user_settings_sig[list].emit(self.settings)
        self.close()
        
    
    def read_user_settings(self):
        # reads all fields and sets new values for main app window
        # read downsampling settings
        if int(self.downsample_text.text()) >= self.min_downsampe_rate and int(self.downsample_text.text()) <= self.max_downsample_rate: 
            self.settings[0]["downsample"] = int(self.downsample_text.text())
            # samples = round(self.current_fs/int(self.downsample_text.currentText()))
            # self.settings[0]["downsample"] = samples
            self.settings[0]["entered_downsample"] = round(self.current_fs/int(self.downsample_text.text()))
        else:
            self.show_info_dialog("Downsample was not updated.\nEnter values between "+str(self.min_downsampe_rate)+" and " + str(self.max_downsample_rate))
        # don't loose filter fraq information
        if self.filter_cb.isChecked() == False:
            self.settings[0]["filter_window"] = DEFAULT_SMOOTH_WINDOW
        else:
            try:
                fraq = int(self.smooth_text.text())
                if fraq > 0 and fraq <= MAX_SMOOTH_WINDOW:
                    self.settings[0]["filter_window"] = fraq
                else:
                    self.show_info_dialog("Fraction of the data has to be between 0 and "+str(MAX_SMOOTH_WINDOW))
            except:
                self.show_info_dialog("Fraction of the data has to be between 0 and "+str(MAX_SMOOTH_WINDOW))
        # read normalization method
        self.settings[0]["normalization"] = self.normalization_method_comboBox.currentText()
        self.settings[0]["show_norm_as"] = self.normalization_show_comboBox.currentText()
        self.settings[0]["filter"] = self.filter_cb.isChecked()
        
    # show info popup
    def show_info_dialog(self, text):
        msgBox = QMessageBox()
        msgBox.setWindowIcon(QtGui.QIcon(ICO))
        msgBox.setText(text)
        msgBox.setWindowTitle("Info!")
        msgBox.setStandardButtons(QMessageBox.Ok)
        msgBox.exec()
        
    def closeEvent(self, event): 
        try:
            # on closing, enable back parent's buttons
            self.parent_window.enable_all_buttons()
        except:
            print("Main does not exist")

    def exit_app(self):
        self.close()
                   
############################## end SettingsWindow class

#################################
# CLASS FOR SELECT DATA WINDOW  #
#################################
        
class SelectDataWindow(QMainWindow):
    # pyqt Signal has to be defined on class level
    # and not in init method !!
    close_select_data_sig = pyqtSignal()
    # create a list with information acquired from window
    # and send it as a signal
    got_user_input_sig = pyqtSignal(list)
    def __init__(self,parent_window,user_input):
        super(SelectDataWindow, self).__init__()
        self.setWindowTitle("Select TDT Data")
        self.setWindowIcon(QtGui.QIcon(ICO))
        self.resize(650,250)
        self.parent_window = parent_window
        self.parent_window.app_closed.connect(self.exit_app)

        self.bold_label_stylesheet = "QLabel {font-weight: bold; font-size: 14px;}"
        
        # list that will contain user input
        self.user_input_list = user_input
        self.experiments_window  = None
        self.signal_control_window = None
        if len(self.user_input_list) > 0:
            self.my_data_folder_path = self.user_input_list[0]['main_path']
            self.subject_experiment_folder_structure = self.user_input_list[0]['subject_experiment']
            self.experiment_subject_folder_structure = self.user_input_list[0]['experiment_subject']
        else:
            self.my_data_folder_path = None
            self.subject_experiment_folder_structure = True
            self.experiment_subject_folder_structure = False
            
        
        self.main_widget = QWidget()
        self.setCentralWidget(self.main_widget)
        # Create an outer layout
        self.outer_layout = QVBoxLayout()
        # layout for experiment check boxes
        self.experiment_layout = QHBoxLayout()
        # set the spacing around the layout
        self.experiment_layout.setContentsMargins(10,10,10,10)
        self.experiment_layout.setAlignment(Qt.AlignCenter)

        self.event_based_label = QLabel("Select Your Event Based Experiments")
        self.event_based_label.setStyleSheet(self.bold_label_stylesheet)
        self.experiment_layout.addWidget(self.event_based_label)
        # nest the inner layout into the outer layout
        self.outer_layout.addLayout(self.experiment_layout)
        
        # add layout for selection of folder structure
        self.folder_structure_layout = QHBoxLayout()
        # set the spacing around the layout
        self.folder_structure_layout.setContentsMargins(10,10,10,10)
        self.folder_structure_layout.setAlignment(Qt.AlignCenter)
        self.folder_structure_button_group = QButtonGroup()
        self.folder_structure_button_group.setExclusive(True)
        self.subject_experiment_radio = QRadioButton("Folder Structure: Subject->Experiments")
        self.subject_experiment_radio.setChecked(self.subject_experiment_folder_structure)
        self.experiment_subject_radio = QRadioButton("Folder Structure: Experiment->Subjects")
        self.experiment_subject_radio.setChecked(self.experiment_subject_folder_structure)
        self.folder_structure_button_group.addButton(self.subject_experiment_radio)
        self.folder_structure_button_group.addButton(self.experiment_subject_radio)
        self.folder_structure_layout.addWidget(self.subject_experiment_radio)
        self.folder_structure_layout.addWidget(self.experiment_subject_radio)
        # nest the inner layout into the outer layout
        self.outer_layout.addLayout(self.folder_structure_layout)
        
        
        # create layout for select folder
        self.select_folder_layout = QHBoxLayout()
        self.select_folder_layout.setContentsMargins(10,10,10,10)
        self.select_folder_btn = QPushButton('Select Main Data Folder')

        if len(self.user_input_list) == 0:
            self.select_folder_text = QLineEdit("")
        else:
            self.select_folder_text = QLineEdit(self.my_data_folder_path)
        self.select_folder_text.setReadOnly(False)
        self.select_folder_layout.addWidget(self.select_folder_btn)
        self.select_folder_layout.addWidget(self.select_folder_text)
        # nest the inner layout into the outer layout
        self.outer_layout.addLayout(self.select_folder_layout)
        
        # create layout for read button
        self.submit_btn_layout = QHBoxLayout()
        self.submit_btn_layout.setAlignment(Qt.AlignRight)
        self.submit_btn = QPushButton('Read Data')
        self.submit_btn_layout.addWidget(self.submit_btn)
        # nest the inner layout into the outer layout
        self.outer_layout.addLayout(self.submit_btn_layout)
        
        # set the window's main layout
        self.main_widget.setLayout(self.outer_layout)
        
        # use stylesheet
        self.main_widget.setStyleSheet(STYLESHEET)
        
        
        self.select_folder_btn.clicked.connect(self.select_folder_btn_clicked)
        self.submit_btn.clicked.connect(self.read_data_folder)
        
        
    def select_folder_btn_clicked(self):
        # set text field to the value of selected path
        self.select_folder_text.setText(QFileDialog.getExistingDirectory(self,"Select folder with raw data"))
        
    def read_data_folder(self):
        if len(self.select_folder_text.text()) == 0:
            self.show_info_dialog("You have not selected data folder")
        else:
            self.my_data_folder_path = self.select_folder_text.text()
            self.subject_folder_names = []
            self.experiment_folder_names = []
            try:
                # check which folder structure to use to find experiment and subjects
                if self.subject_experiment_radio.isChecked():                    
                    my_subject_folder_paths = os.listdir(self.my_data_folder_path)                    
                    for el in my_subject_folder_paths:
                        # check if subfolder is a folder not a file
                        subfolder_split = el.split(".")
                        if len(subfolder_split) == 1:
                            self.subject_folder_names.append(el)
                            my_experiment_folders = os.listdir(os.path.join(self.my_data_folder_path,el))
                            for exp in my_experiment_folders:
                                exp_split = exp.split(".")
                                if len(exp_split) == 1:
                                    if exp not in self.experiment_folder_names:
                                        self.experiment_folder_names.append(exp)
                elif self.experiment_subject_radio.isChecked():
                    my_experiment_folder_paths = os.listdir(self.my_data_folder_path)                    
                    for el in my_experiment_folder_paths:
                        # check if subfolder is a folder not a file
                        subfolder_split = el.split(".")
                        if len(subfolder_split) == 1:
                            self.experiment_folder_names.append(el)
                            my_subjects_folders = os.listdir(os.path.join(self.my_data_folder_path,el))
                            for subj in my_subjects_folders:
                                subj_split = subj.split(".")
                                if len(subj_split) == 1:
                                    if subj not in self.subject_folder_names:
                                        self.subject_folder_names.append(subj)
                # proceed only if any subjects or experiments were found
                if (len(self.subject_folder_names) == 0 or len(self.experiment_folder_names) == 0):
                    self.show_info_dialog("Wrong data folder structure.\n(Folder structure:\nSubjects>Experiments)")
                else:
                    # sort lists
                    self.subject_folder_names.sort()
                    self.experiment_folder_names.sort()
                    self.experiments_window = ExperimentsWindow(self.parent_window,self,self.experiment_folder_names)
                    self.experiments_window.got_user_input_sig.connect(self.get_user_experiment)
                    self.experiments_window.show()
            except:
                self.show_info_dialog("Wrong data folder structure.\n(Folder structure:\nSubjects>Experiments)")
                
    @pyqtSlot(str) 
    # receives a string with experiment name from experiments window
    def get_user_experiment(self, exp):
       
        if self.folder_structure_button_group.checkedButton().text() == "Folder Structure: Subject->Experiments":
            self.subject_experiment_folder_structure = True
            self.experiment_subject_folder_structure = False
        elif self.folder_structure_button_group.checkedButton().text() == "Folder Structure: Experiment->Subjects":
            self.subject_experiment_folder_structure = False
            self.experiment_subject_folder_structure = True
                
        data_dict = {"main_path":self.my_data_folder_path,
                     "subject_names":self.subject_folder_names,
                     "selected_experiment":exp,
                     "subject_experiment":self.subject_experiment_folder_structure,
                     "experiment_subject":self.experiment_subject_folder_structure
                }
        # clear previous inputs
        self.user_input_list = [data_dict]
        # send that to main app window
        self.got_user_input_sig[list].emit(self.user_input_list)
#        self.close()
        
    # emit signal when window is closed
    def closeEvent(self, event): 
        try: 
            self.close_select_data_sig.emit()
        except:
            print("Main does not exist")

    def exit_app(self):
        self.close()
        
    # show info popup
    def show_info_dialog(self, text):
        msgBox = QMessageBox()
        msgBox.setWindowIcon(QtGui.QIcon(ICO))
        msgBox.setText(text)
        msgBox.setWindowTitle("Info!")
        msgBox.setStandardButtons(QMessageBox.Ok)
        msgBox.exec()
        
################### end SelectDataWindow class
        
########################################
# CLASS FOR SELECTING EXPERIMENT NAME  #
########################################   
        
class ExperimentsWindow(QMainWindow):
    # pyqt signal to send experiment name
    got_user_input_sig = pyqtSignal(str)
    def __init__(self,main_window,parent_window,user_input):
        super(ExperimentsWindow, self).__init__()
        self.setWindowTitle("Experiments")
        self.setWindowIcon(QtGui.QIcon(ICO))
        self.resize(300,150)
        self.parent_window = parent_window
        self.main_window = main_window
        self.main_window.app_closed.connect(self.exit_app)
        self.experiments = user_input
            
        self.main_widget = QWidget()
        self.setCentralWidget(self.main_widget)
        
        self.main_layout = QVBoxLayout()
        self.main_layout.setContentsMargins(10,10,10,10)
        self.radio_btn_layout = QVBoxLayout()
        self.radio_btn_layout.setContentsMargins(10,10,10,30)
        self.radio_btn_layout.setAlignment(Qt.AlignCenter)
        self.experiment_button_group = QButtonGroup()
        self.experiment_button_group.setExclusive(True)
        text_label = QLabel("What is the name of your experiment?\n")
        self.radio_btn_layout.addWidget(text_label)
        for i in range(len(self.experiments)):
            # create radio btn
            # add button to a group
            # add to list
            # add to layout
            radio_btn = QRadioButton(self.experiments[i])
            if i == 0: # check if this is the first one
                radio_btn.setChecked(True)
            self.radio_btn_layout.addWidget(radio_btn)
            self.experiment_button_group.addButton(radio_btn)
            
        self.main_layout.addLayout(self.radio_btn_layout)
        
        self.submit_layout = QHBoxLayout()
        self.submit_layout.setAlignment(Qt.AlignRight)
        self.submit_btn = QPushButton('Read Data')
        self.submit_layout.addWidget(self.submit_btn)
        self.main_layout.addLayout(self.submit_layout)
        
        
        # set the window's main layout
        self.main_widget.setLayout(self.main_layout)
        
        self.submit_btn.clicked.connect(self.read_experiments)
        
        
    def read_experiments(self):
        if self.experiment_button_group.checkedId() == -1:
            print("No experiment selected")
            self.close()
        else:
            self.got_user_input_sig[str].emit(self.experiment_button_group.checkedButton().text())
#            print(self.experiment_button_group.checkedButton().text())
            self.close()

    def exit_app(self):
        self.close()
            
################### end ExperimentsWindow class   

########################################
# CLASS FOR SELECTING EXPERIMENT NAME  #
########################################   
        
class SignalControlWindow(QMainWindow):
    # pyqt signal to send experiment name
    got_signal_name_sig = pyqtSignal(list)
    def __init__(self,parent_window,channel_names_list):
        super(SignalControlWindow, self).__init__()
        self.setWindowTitle("Select Signal")
        self.setWindowIcon(QtGui.QIcon(ICO))
        self.resize(300,100)
        self.parent_window = parent_window
        self.parent_window.app_closed.connect(self.exit_app)
        self.channel_names = channel_names_list
        self.selected_channel_names = []
            
        self.main_widget = QWidget()
        self.setCentralWidget(self.main_widget)
        
        self.main_layout = QFormLayout()
        self.main_layout.setContentsMargins(10,10,10,10)
        self.signal_name_cb = QComboBox()
        self.signal_name_cb.addItems(self.channel_names)
        self.main_layout.addRow("Signal channel name:",self.signal_name_cb)
        self.control_name_cb = QComboBox()
        self.control_name_cb.addItems(self.channel_names)
        self.main_layout.addRow("Control channel name:",self.control_name_cb)
        self.ok_btn = QPushButton("OK")
        self.main_layout.addRow("",self.ok_btn)
                
        self.main_widget.setLayout(self.main_layout)
        
        self.ok_btn.clicked.connect(self.read_names)
        
    def read_names(self):
        self.ok_btn.setEnabled(False)
        QApplication.processEvents() 
        signal_name = self.signal_name_cb.currentText()
        control_name = self.control_name_cb.currentText()
        if signal_name != control_name:
            temp_dict = {"signal_name":signal_name,
                         "control_name":control_name}
            self.selected_channel_names.append(temp_dict)
        self.close()
        
    def closeEvent(self,evt):
        self.got_signal_name_sig[list].emit(self.selected_channel_names)

    def exit_app(self):
        self.close()
################### end SignalControlWindow class   

##############################################
# CLASS FOR PREVIEW EVENT-BASED EXPERIMENTS  #
##############################################
            
class PreviewEventBasedWidget(QWidget):
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
              
        # options read from options section
        self.options = {"subject":self.preview_init_params[0][0]["subject_names"][0],
                        "event":""
                        }
        
        # a dictionary with subject:extracted raw data
        self.raw_data_dict = {}
        self.cam_data_dict = {}
        # start a widget by reading the first subject data on the list and adding that to dict
        self.get_raw_data(self.preview_init_params[0][0]["subject_names"][0],
                          self.preview_init_params[1][self.preview_init_params[0][0]["subject_names"][0]])
        # read first subject's frequency and create suggested downsampled rate
        self.current_fs = fpVideoExplorer_functions.get_frequency(self.raw_data_dict[self.preview_init_params[0][0]["subject_names"][0]],self.preview_init_params[0][0]["signal_name"])
        # self.suggested_downsample_samples = int(int(self.current_fs)*DEFAULT_DOWNSAMPLE_PCT/100)
        self.suggested_downsample_samples = DEFAULT_HZ
        self.settings_dict[0]['downsample'] = self.suggested_downsample_samples
        self.suggested_downsample_rate = round(self.current_fs/self.suggested_downsample_samples)
        self.settings_dict[0]["entered_downsample"] = self.suggested_downsample_rate
        # print("Current suggested downsample rate:",self.settings_dict[0]['downsample'])
        # keep trimmed data separately with latest trimming settings
        # first element of that list is beginning and end seconds to trim
        # second element is trimmed data dict ts:timestamps, signal:data,control:data
        self.trimmed_raw_data_dict = {}
        
        # create a list of available events for current subject
        self.events_from_current_subject = fpVideoExplorer_functions.get_events(self.raw_data_dict[self.preview_init_params[0][0]["subject_names"][0]])
        self.cam_events_from_current_subject = fpVideoExplorer_functions.get_cam_events(self.raw_data_dict[self.preview_init_params[0][0]["subject_names"][0]])
        print(f"Cam events: {self.cam_events_from_current_subject}")
        if len(self.cam_events_from_current_subject) > 0:
            self.options["cam_event"] = self.cam_events_from_current_subject[0]
        else:
            self.show_info_dialog("No Cam events timing video frames for the first subject.")
        # store downsampled data key (subject): dict("ts","signal","control")
        self.downsampled_dict = {}
        # store normalized data key (subject): dict("ts","normalized signal")
        self.normalized_dict = {}
        self.perievent_window = None
        # options sent from perievent window
        self.perievent_options_dict = {}
        self.export_window = None
        # remember previous event from perievent window
        self.current_perievent = None
        # remember selected trials (integers indicating trial number)
        self.current_trials = []
        # remember how many were there
        self.total_current_trials = 0
        # key is subject, value is a dict where key is event and value is list of most recent trials for that event
        self.trials_dict = {}

        self.select_trials_window = None

        self.current_trials_df = None
        self.current_saving_path = None
        self.video_path = None
        self.data_by_event = None
        self.cam_msg_displayed = 0
        
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
        
        
        self.disable_buttons_signal.connect(self.disable_buttons)
        self.enable_buttons_signal.connect(self.enable_buttons)
        
        self.bold_label_stylesheet = "QLabel {font-weight: bold}"
        self.bold_cb_stylesheet = "QCheckBox {font-weight: bold}"
        self.normal_btn_stylesheet = \
            ".QPushButton{\n" \
            + "padding-left: 5px;\n" \
            + "padding-right: 5px;\n" \
            + "padding-top: 1px;\n"\
            + "padding-bottom: 1px;\n" \
            + "margin-left:0px;\n" \
            + "margin-right:0px;\n" \
            + "margin-top:0px;\n" \
            + "margin-bottom:0px;\n" \
            + "font-size: 6;\n" \
            + "font-weight: bold;\n" \
            + "}" 
        # create gui
        self.setupGUI()

    # receives a list with dictionary of signal and control channel names 
    def setupGUI(self):
        self.layout = QVBoxLayout()
        self.layout.setContentsMargins(5,0,5,0)
        self.setLayout(self.layout)

        # create widget with options
        self.options_widget = QWidget()
        # create main layout for widget with options
        self.options_main_layout = QVBoxLayout()
        # self.options_main_layout.setContentsMargins(10,10,0,10)
        self.options_layout = QFormLayout()
        self.options_layout.setFieldGrowthPolicy(self.options_layout.AllNonFixedFieldsGrow)
        self.options_layout.setVerticalSpacing(15)
        self.experiment_name_text = QLineEdit(self.preview_init_params[0][0]["selected_experiment"])
        self.experiment_name_text.setReadOnly(True)
        self.experiment_label = QLabel("Experiment")
        self.experiment_label.setStyleSheet(self.bold_label_stylesheet)
        self.options_layout.addRow(self.experiment_label,self.experiment_name_text)
        # add drop down menu for subjects from the list
        self.subject_comboBox = QComboBox()
        self.subject_comboBox.addItems(self.preview_init_params[0][0]["subject_names"])
        self.subject_comboBox.currentIndexChanged.connect(self.on_subject_change)
        self.subject_label = QLabel("Subject")
        self.subject_label.setStyleSheet(self.bold_label_stylesheet)
        self.options_layout.addRow(self.subject_label, self.subject_comboBox)
        
        self.cam_event_from_data_comboBox = QComboBox()
        self.cam_event_from_data_comboBox.addItems(self.cam_events_from_current_subject)
        self.cam_event_from_data_comboBox.currentIndexChanged.connect(self.on_cam_evt_change)
        self.show_cam_event_label = QLabel("Select Frames")
        self.show_cam_event_label.setStyleSheet(self.bold_label_stylesheet)
        self.options_layout.addRow(self.show_cam_event_label,self.cam_event_from_data_comboBox)  
        self.event_from_data_comboBox = QComboBox()
        self.event_from_data_comboBox.addItems(self.events_from_current_subject)
        self.show_event_label = QLabel("Select Event")
        self.show_event_label.setStyleSheet(self.bold_label_stylesheet)
        self.options_layout.addRow(self.show_event_label,self.event_from_data_comboBox)       
        self.event_name_text = QLineEdit("")
        self.event_name_text.setToolTip("i.e, Tone")
        self.options_layout.addRow("Event name (optional)",self.event_name_text)
        self.select_window_label = QLabel("Select your window around the event:")
        self.select_window_label.setStyleSheet(self.bold_label_stylesheet)
        self.options_layout.addRow(self.select_window_label)
        self.before_sec_label = QLabel("Time before event (sec):")
        self.before_sec_text = QLineEdit("15")
        self.before_sec_text.setValidator(QtGui.QIntValidator())
        self.options_layout.addRow(self.before_sec_label,self.before_sec_text)
        self.after_sec_label = QLabel("Time after event (sec):")
        self.after_sec_text = QLineEdit("15")
        self.after_sec_text.setValidator(QtGui.QIntValidator())
        self.options_layout.addRow(self.after_sec_label,self.after_sec_text)
        self.shift_label = QLabel("Shift n Video Frames:")
        self.shift_text = QLineEdit("0")
        shift_info = "Between 0 and " +str(SHIFT_MAX)+" frames"
        self.shift_text.setToolTip(shift_info)
        self.shift_text.setValidator(QtGui.QIntValidator())
        self.options_layout.addRow(self.shift_label,self.shift_text)
        self.select_folder_btn = QPushButton("Save to folder")
        self.select_folder_btn.setStyleSheet(self.normal_btn_stylesheet)
        default_saving_path = os.path.join(self.parent_window.preview_params[1][self.options["subject"]],DEFAULT_EXPORT_FOLDER)
        self.path2save_text = QLineEdit(default_saving_path)
        self.options_layout.addRow(self.select_folder_btn,self.path2save_text)
              
        # nest the inner layout into the outer layout
        self.options_main_layout.addLayout(self.options_layout)
        
 
        # add perievent button
        self.perievent_analysis_btn = QPushButton("Save Trial Videos")
        self.options_main_layout.addWidget(self.perievent_analysis_btn)
               
        self.options_widget.setLayout(self.options_main_layout)
        self.layout.addWidget(self.options_widget)

        self.perievent_analysis_btn.clicked.connect(self.perievent_analysis_btn_clicked)
        self.select_folder_btn.clicked.connect(self.select_folder_btn_clicked)

    # show a reminder when cam event is changed
    def on_cam_evt_change(self):
        if self.cam_msg_displayed == 0 and len(self.cam_events_from_current_subject) > 1:
            self.show_info_dialog("Make sure that selected chanel for signal and control are for that recording.\nMake sure to select the events for the same channel.")
            self.cam_msg_displayed +=1
     
    @pyqtSlot()   
    def disable_buttons(self):
        self.perievent_analysis_btn.setEnabled(False)
        self.parent_window.disable_buttons()
        QApplication.processEvents() 
    
    @pyqtSlot()     
    def enable_buttons(self):
        self.perievent_analysis_btn.setEnabled(True)
        self.parent_window.enable_all_buttons()
        QApplication.processEvents() 

    def select_folder_btn_clicked(self):
        # set text field to the value of selected path
        self.path2save_text.setText(QFileDialog.getExistingDirectory(self,"Save to folder"))
            
            
    # function to read options for plotting on the left side of the plot    
    def read_options(self):
        # set current subject
        self.options["subject"] = self.subject_comboBox.currentText()
        shift = int(self.shift_text.text())
        if shift < 0 or shift > SHIFT_MAX:
            self.show_info_dialog("If you need to shift so much, you might have a problem with your data.") 
            return False
        else:
            self.options["shift"] = shift
        # if subject was not previewed yet, read data and add to dict
        if self.options["subject"] not in self.raw_data_dict:
            self.get_raw_data(self.options["subject"],
                          self.parent_window.preview_params[1][self.options["subject"]])
        if self.raw_data_dict[self.options["subject"]] == None:
            self.show_info_dialog("Could not read raw data.")
            return False
        else:
            # slect cam event
            self.options["cam_event"] = self.cam_event_from_data_comboBox.currentText()
            print("\nSelected Cam Event:",self.options["cam_event"])
            self.get_cam_data(self.options["subject"],
                          self.raw_data_dict[self.options["subject"]],self.options["cam_event"])
            # check if there even is a video there
            has_video, self.video_path = fpVideoExplorer_functions.is_there_a_video(self.parent_window.preview_params[1][self.options["subject"]],self.options["cam_event"])
            if has_video == False:
                self.show_info_dialog("There is no video for this subject.\nCheck that the video ends with _"+self.options["cam_event"])
                return False
        try:
            self.options["sec_before"] = abs(int(self.before_sec_text.text()))
            self.options["sec_after"] = abs(int(self.after_sec_text.text()))  
        except:
            self.before_sec_text.setText("0")
            self.after_sec_text.setText("0")
            self.options["sec_before"] = 0
            self.options["sec_after"] = 0
            self.show_info_dialog("Wrong time window around the event.\Try only intigers.")
            return False
        
        # set event
        self.options["event"] = self.event_from_data_comboBox.currentText()
        self.options["event_name"] = self.event_name_text.text()
        self.event_onsets = fpVideoExplorer_functions.get_event_on_off(self.raw_data_dict[self.options["subject"]], self.options["event"])

        # # create trials list 
        # self.current_trials = [i+1 for i in range(len(self.event_onsets[0]))]
        print("\nEvent:",self.options["event"])
        print("This event onsets", self.event_onsets[0])
        if len(self.event_onsets[0]) == 0:
            self.show_info_dialog("Some "+self.options["event"]+" event data is missing.")
            return False
        if os.path.exists(self.path2save_text.text()):
            self.current_saving_path = self.path2save_text.text()
        else:
            try: # if path did not exist, try creating one
                os.mkdir(self.path2save_text.text())
                self.current_saving_path = self.path2save_text.text()
            except:
                self.show_info_dialog("Selected saving path does not exist")
                return False
                                                                  
        return True

 
    # check and set event for current subject      
    def check_events(self):
        self.cam_msg_displayed = 0
        self.options["subject"] = self.subject_comboBox.currentText()
        # if subject was not previewed yet, read data and add to dict
        if self.options["subject"] not in self.raw_data_dict:
            self.get_raw_data(self.options["subject"],
                          self.parent_window.preview_params[1][self.options["subject"]])
            if self.raw_data_dict[self.options["subject"]] == None:
                return
        # check events for new data
        previous_cam_events = self.cam_events_from_current_subject
        self.cam_events_from_current_subject = fpVideoExplorer_functions.get_cam_events(self.raw_data_dict[self.options["subject"]])
        print("\nNew Cam events",self.cam_events_from_current_subject)
        same_cam = True
        if len(previous_cam_events) != self.cam_events_from_current_subject:
            same_cam = False
        else: # check if they are the same events
            for el in self.cam_events_from_current_subject:
                if el not in previous_cam_events:
                    same_cam = False
        if same_cam == True:
            if self.cam_event_from_data_comboBox.setCurrentText() != self.options["cam_event"]:
                self.cam_event_from_data_comboBox.setCurrentText(self.options["cam_event"])
        else: # if events changed
            # if previously selected event is in current events, set that filed the same as previously
            if self.options["cam_event"] in self.cam_events_from_current_subject:
                # update combo box first
                self.cam_event_from_data_comboBox.clear()
                self.cam_event_from_data_comboBox.addItems(self.cam_events_from_current_subject)
                # set to previous event
                self.cam_event_from_data_comboBox.setCurrentText(self.options["cam_event"])
            else:
                self.cam_event_from_data_comboBox.clear()
                self.cam_event_from_data_comboBox.addItems(self.cam_events_from_current_subject)
       
        # first check if there are the same events
        previous_events = self.events_from_current_subject
        self.events_from_current_subject = fpVideoExplorer_functions.get_events(self.raw_data_dict[self.options["subject"]])  
        same = True
        if len(previous_events) != self.events_from_current_subject:
            same = False
        else: # check if they are the same events
            for el in self.events_from_current_subject:
                if el not in previous_events:
                    same = False
        if same == True:
            self.event_from_data_comboBox.setCurrentText(self.options["event"])
            # set event name also to previous
            self.event_name_text.setText(self.options["event_name"])
        else: # if events changed
            # if previously selected event is in current events, set that filed the same as previously
            if self.options["event"] in self.events_from_current_subject:
                # update combo box first
                self.event_from_data_comboBox.clear()
                self.event_from_data_comboBox.addItems(self.events_from_current_subject)
                # clear any previous event names
                self.event_name_text.setText("")
                # set to previous event
                self.event_from_data_comboBox.setCurrentText(self.options["event"])
                # set event name also to previous
                self.event_name_text.setText(self.options["event_name"])
            else:
                self.event_from_data_comboBox.clear()
                self.event_from_data_comboBox.addItems(self.events_from_current_subject)
                # clear any previous event names
                self.event_name_text.setText("")
            
        
    # check trimming and event at the same time        
    def options_pre_check(self):
        self.check_events()
        # update saving path
        saving_path = os.path.join(self.parent_window.preview_params[1][self.options["subject"]],DEFAULT_EXPORT_FOLDER)
        self.path2save_text.setText(saving_path) 
            
    def on_subject_change(self):
        self.options_pre_check()
        if self.raw_data_dict[self.options["subject"]] == None:
            return
        # read new subject's frequency and update suggested downsampled rate (only if it is different)
        new_fs = fpVideoExplorer_functions.get_frequency(self.raw_data_dict[self.preview_init_params[0][0]["subject_names"][0]],self.preview_init_params[0][0]["signal_name"])
        if new_fs != self.current_fs:
            # self.suggested_downsample_samples = int(int(self.current_fs)*DEFAULT_DOWNSAMPLE_PCT/100)
            self.suggested_downsample_samples = DEFAULT_HZ
            self.settings_dict[0]['downsample'] = self.suggested_downsample_samples
            self.suggested_downsample_rate = int(int(self.current_fs)/self.suggested_downsample_samples)
            self.settings[0]["entered_downsample"] = self.suggested_downsample_rate
        
        
            
    def perievent_analysis_btn_clicked(self):
        self.disable_buttons_signal.emit()
        all_data_read = self.read_options()
        if all_data_read:
            self.data_by_event = fpVideoExplorer_functions.filter_data_around_event_directly(self.raw_data_dict[self.options["subject"]],
                                                                            self.options["event"],
                                                                            self.options["sec_before"],
                                                                            self.options["sec_after"],
                                                                            self.settings_dict,
                                                                            self.preview_init_params[0][0]["signal_name"],
                                                                            self.preview_init_params[0][0]["control_name"])
            # show popup with trials selection
            if (len(self.data_by_event.streams[self.preview_init_params[0][0]["signal_name"]].filtered) > 0) and (len(self.data_by_event.streams[self.preview_init_params[0][0]["control_name"]].filtered) > 0):
                trials_no = len(self.data_by_event.streams[self.preview_init_params[0][0]["signal_name"]].filtered)
                self.select_trials_window = SelectTrialsWindow(self,trials_no)  
                self.select_trials_window.got_trials_sig.connect(self.selected_trials_received)
                self.select_trials_window.show()
            else:   # if there was not enough data
                self.show_info_dialog("Not enough data that satisfy your request.\nTry changing event or times around the event.")
                self.enable_buttons_signal.emit()
        else:
            self.enable_buttons_signal.emit()

    @pyqtSlot(list)     
    def selected_trials_received(self,trials):
        print("\nSelected event:", self.options["event"])
        print(f"Selected trials: {trials}")
        # check if user named the event
        event = self.options["event"]
        if len(self.options["event_name"]) > 0:
            event = self.options["event_name"]
        # list of trials
        self.current_trials = trials
        # get normalized trials to preview
        self.current_trials_df = fpVideoExplorer_functions.get_normalized_trials(
                                                    self.options["subject"],
                                                    self.data_by_event,
                                                    self.current_trials,
                                                    event,
                                                    self.options["sec_before"],
                                                    self.options["sec_after"],
                                                    self.settings_dict,
                                                    self.preview_init_params[0][0]["signal_name"],
                                                    self.preview_init_params[0][0]["control_name"]
                                                    )  
        # print(f"Trials df:\n{self.current_trials_df}")
        # look for matching trial in frames
        if self.options["subject"] in self.cam_data_dict:
            current_cam_ons = self.cam_data_dict[self.options["subject"]][0]
            # get mean interval
            diffs = np.diff(current_cam_ons)
            mean_frame_interval = round(np.mean(diffs),3)
            print(f"Mean interval between frames from cam event: {mean_frame_interval}\n")
            video = cv2.VideoCapture(self.video_path)
            fps = video.get(cv2.CAP_PROP_FPS) # check FPS to shift frames later if it is 10fps
            video.release()
            # find frame ind before and after for each trial (for each trial there will be a tuple: trial, list of size=2(before and after))
            all_trials_before_and_after_frame_idx = []
            for i in range(len(self.current_trials)):
                trial_idx = self.current_trials[i]-1 # trials start from 1
                trial_ts = self.event_onsets[0][trial_idx] + current_cam_ons[0]# correct for delay + VIDEO_DELAY_SEC?...First cam event is the offset
                before_and_after_frame_idx = []
                print(f"Trial {self.current_trials[i]}, trial ts: {trial_ts}; before and after from cam evt:")
                for j in range(len(current_cam_ons)):
                    if current_cam_ons[j] > trial_ts-mean_frame_interval and current_cam_ons[j] < trial_ts+mean_frame_interval:
                        print(f"Frame {j} ts: {current_cam_ons[j]}")
                        before_and_after_frame_idx.append(j)
                if len(before_and_after_frame_idx)==2:
                    # print(f"Trial {self.current_trials[i]}, pair of frames found")
                    all_trials_before_and_after_frame_idx.append((self.current_trials[i],before_and_after_frame_idx))
                else:
                    if len(before_and_after_frame_idx) > 2:
                        # include only the closest to the event time
                        updated_before_after = fpVideoExplorer_functions.find_closest_time_pair_indexes(trial_ts,found_index_list=before_and_after_frame_idx,times_to_compare=current_cam_ons)
                        if len(updated_before_after) == 2 :
                            print(f"Narrowed down successfully")
                            all_trials_before_and_after_frame_idx.append((self.current_trials[i],updated_before_after))
                        else: # if a pair was found, try increasing mean frame interval and repeat the above
                            updated_before_after = fpVideoExplorer_functions.find_closest_time_pair_indexes(trial_ts,times_to_compare=current_cam_ons,time_margin=mean_frame_interval+0.03)
                            if len(updated_before_after) == 2:
                                all_trials_before_and_after_frame_idx.append((self.current_trials[i],updated_before_after))
                            elif len(updated_before_after) > 2:
                                # include only the closest to the event time
                                updated_before_after_pair = fpVideoExplorer_functions.find_closest_time_pair_indexes(trial_ts,found_index_list=updated_before_after,times_to_compare=current_cam_ons)
                                if len(updated_before_after_pair) == 2 :
                                    all_trials_before_and_after_frame_idx.append((self.current_trials[i],updated_before_after_pair))
                                else: # give up
                                    self.show_info_dialog("Could not find before and after frame pair for the video")
                                    all_trials_before_and_after_frame_idx = []
                                    break
        if len(all_trials_before_and_after_frame_idx) != 0:
            add_shift = 0
            if fps == 10:
                add_shift = 2 # add additional shift if the video was slower
            shifted = fpVideoExplorer_functions.shift_indexes_right(all_trials_before_and_after_frame_idx,self.options["shift"]+add_shift) 
            print(f"shifted trials before and after: {shifted}")
            ################################################################################################################
            # create animated trials first
            # for now include all trials
            fps = int(round(1/mean_frame_interval))
            # create folder first
            head,tail = os.path.split(self.video_path)
            video_name = tail.split(".")[0]
            path4videos = os.path.join(self.current_saving_path,video_name+"_shift"+str(self.options["shift"])+"_videos")
            # create as many animations as in current trials
            created_animations = fpVideoExplorer_functions.create_trial_animations(self.current_trials_df,self.current_trials,event,fps,path4videos,self.options["subject"],self.options["cam_event"])
            ####################################################################################################################
            # create as many videos as animations
            if created_animations:
            # if True:
                created_trial_videos = fpVideoExplorer_functions.create_trial_videos(self.video_path,
                                                                                self.options["sec_before"],
                                                                                self.options["sec_after"],
                                                                                shifted,
                                                                                event,
                                                                                self.current_saving_path,
                                                                                path4videos,
                                                                                self.options["subject"],
                                                                                self.options["cam_event"])
                if created_trial_videos == False:
                    self.show_info_dialog("Did not create trial videos or combined videos.\nPlease check the terminal for more details.")
            else:
                self.show_info_dialog("No trial animations created. Therefore no videos either")
            #############################################################################################################################
        else:
            self.show_info_dialog("Did not create animations or trial videos")
        self.enable_buttons_signal.emit()               


    # add raw data structure and subject to self.raw_data dictionary   
    def get_raw_data(self,subject,my_path):
        self.raw_data_dict[subject] = fpVideoExplorer_functions.get_raw_data(my_path)
        if self.raw_data_dict[subject] == None:
            # remove from combo box
            index = self.subject_comboBox.findText(subject)  # find the index of text
            self.subject_comboBox.removeItem(index)  # remove item from index
            # and remove from available subjects
            self.preview_init_params[0][0]["subject_names"].remove(subject)
            self.show_info_dialog("Problem reading subject's "+subject+" file.\nSubject will not be available for the analysis.")

    # add cam1 event data and subject to self.raw_data dictionary        
    def get_cam_data(self,subject,data,cam_evt):
        self.cam_data_dict[subject] = fpVideoExplorer_functions.get_cam_event_on_off(data,cam_evt)
        
    # show info popup
    def show_info_dialog(self, text):
        msgBox = QMessageBox()
        msgBox.setWindowIcon(QtGui.QIcon(ICO))
        msgBox.setText(text)
        msgBox.setWindowTitle("Info!")
        msgBox.setStandardButtons(QMessageBox.Ok)
        msgBox.exec()
        

########################### end PreviewEventBasedWidget class  

###############################
# CLASS FOR SELECTING TRIALS  #
###############################   
        
class SelectTrialsWindow(QMainWindow):
    # pyqt signal to send experiment name
    got_trials_sig = pyqtSignal(list)
    def __init__(self,parent_window,trials_list):
        super(SelectTrialsWindow, self).__init__()
        self.setWindowTitle("Select Trials")
        self.setWindowIcon(QtGui.QIcon(ICO))
        self.resize(400,100)
        self.parent_window = parent_window
        # self.parent_window.app_closed.connect(self.exit_app)
        self.channel_names = trials_list
        self.checked_trials = []

        self.bold_label_stylesheet = "QLabel {font-weight: bold}"
            
        self.main_widget = QWidget()
        self.setCentralWidget(self.main_widget)
        
        self.main_layout = QVBoxLayout()
        self.main_layout.setContentsMargins(10,10,10,10)
        self.trials_label = QLabel("Select trials:")
        self.trials_label.setStyleSheet(self.bold_label_stylesheet)
        self.main_layout.addWidget(self.trials_label)

        # update trials according to batch settings
        self.trials_button_group = []
        # create new buttons
        self.trials_layout = QHBoxLayout()
        self.trials_layout.setAlignment(Qt.AlignCenter)
        for i in range(trials_list):
            # create btn
            # add button to a group
            # add to list
            # add to layout
            btn = QCheckBox(str(i+1))
            btn.setChecked(True)
            self.trials_layout.addWidget(btn)
            self.trials_button_group.append(btn)    
        self.main_layout.addLayout(self.trials_layout)
        self.btn_layout = QHBoxLayout()
        self.btn_layout.setAlignment(Qt.AlignRight)
        self.select_btn = QPushButton("Done")
        self.btn_layout.addWidget(self.select_btn)
        self.main_layout.addLayout(self.btn_layout)  
        self.main_widget.setLayout(self.main_layout)
        
        self.select_btn.clicked.connect(self.read_trials)
        
    def read_trials(self):
        for btn in self.trials_button_group:
            # print(btn.text(),btn.isChecked())
            if btn.isChecked():
                self.checked_trials.append(int(btn.text()))
        # print("selected trials",self.checked_trials)
        self.select_btn.setEnabled(False)
        QApplication.processEvents() 
        self.close()
        
    def closeEvent(self,evt):
        self.got_trials_sig[list].emit(self.checked_trials)

    def exit_app(self):
        self.close()
################### end SelectTrialsWindow class  
 

################################################################
#                                                              #
# EXECUTE GUI FROM MAIN                                        #
#                                                              #
################################################################
if __name__ == "__main__":
    # Always start by initializing Qt (only once per application)
    app = QApplication([])
    main_widget = MyMainWidget()
    main_widget.show()
    # # show select data window on top when app is started
    # main_widget.select_data_window.show()
    app.exec_()

    print('Done')