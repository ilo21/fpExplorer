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
Created on Thu Jan  7 12:08:47 2021

@author: ilosz01
"""

import fpExplorer_functions
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
import warnings
warnings.filterwarnings("ignore")

#icon https://logomakr.com/9OpQfD
# description how to create icons resources
# https://www.learnpyqt.com/courses/packaging-and-distribution/packaging-pyqt5-pyside2-applications-windows-pyinstaller/
import resources


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
            
# MAX_DOWNSAMPLE = 1000
# how low in downsampling can go
MAX_DOWNSAMPLE_PCT = 50
# how high in downsampling can you go
MIN_DOWNSAMPLE_PCT = 0.5
# default_samples to average 1% of original sampling rate
DEFAULT_DOWNSAMPLE_PCT = 1
MAX_SMOOTH_FRAQ = 100
DEFAULT_SMOOTH_FRAQ = 1
DEFAULT_EXPORT_FOLDER = "_fpExplorerAnalysis"
###############################
# CLASS FORM MAIN GUI WINDOW  #
###############################
            
class MyMainWidget(QMainWindow):
    # signal that select folder window is open
    select_data_open_signal = pyqtSignal()
    got_channel_names = pyqtSignal()
    start_batch_processing_sig = pyqtSignal()
    # emit signal to close all other windows
    app_closed = pyqtSignal() 
    def __init__(self):
        super(MyMainWidget, self).__init__()
        self.name = "fpExplorer"
        self.app_width = 1200
        self.app_height = 600
        self.top_buttons_height = 50
        self.bottom_buttons_height = 50
        self.setWindowIcon(QtGui.QIcon(ICO))
        self.setWindowTitle(self.name)
        self.resize(self.app_width,self.app_height)
        
        # main apps parameters to remember
        # user input from select data window
        self.select_data_window_content = []
        self.select_data_window = None
        # when app is started creeate select data window to show on top of the app
        if (len(self.select_data_window_content)==0):
            self.select_data_window = SelectDataWindow(self,self.select_data_window_content)
            self.select_data_window.got_user_input_sig.connect(self.get_select_data_window_user_input)
        # list with general settings dictionary as first element 
        self.settings_dict = [{"downsample":None,
                              "entered_downsample":None,
                              "normalization": "Standard Polynomial Fitting",
                              "show_norm_as":"dF/F (in % )",
                              "filter":False,
                              "filter_fraction":DEFAULT_SMOOTH_FRAQ,
                              "subject":""}]
        # window with general settings
        self.settings_window = None
        # dict containing paths to all data tanks in the folder
        # subject:path
        self.batch_paths_dict = {}
        # list of init_params to start previewing
        # first element is a dictionary from select data folder window
        # second element is a dict of paths to all data folders
        self.preview_params = []
        # widget with preview
        self.preview_widget = None
        # a list of subjects selected by the user
        self.selected_subjects = []
        # remember previous run on batch settings
        self.batch_export_settings_dict = {}
        
        self.run_on_batch_window = None
       
        
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
        self.select_data_btn = QPushButton('Select Data')
#        self.preview_single_btn = QPushButton('Preview Single')
#        self.run_on_batch_btn = QPushButton('Run On Batch')
        self.settings_btn = QPushButton('Settings')
        self.run_on_batch_btn = QPushButton('Run on Batch')
        self.open_doc_btn = QPushButton('Documentation')
        
        # disable all buttons except settings when the app starts
        self.settings_btn.setEnabled(False)
        self.run_on_batch_btn.setEnabled(False)
        
        # Create a layout to manage the buttons size and position
        self.box_layout = QHBoxLayout()
        # set the spacing around the layout
        self.box_layout.setContentsMargins(10,10,10,10)
        # place widgets to the layout in their proper positions
        self.box_layout.addWidget(self.select_data_btn)
        self.box_layout.addWidget(self.run_on_batch_btn)
        self.box_layout.addWidget(self.settings_btn)
        self.box_layout.addWidget(self.open_doc_btn)
        # create a widget for dock from pyqtgraph
        self.top_widget = QWidget()
        # add that widget to the dock
        self.top_dock_widget.addWidget(self.top_widget)
        self.top_widget.setLayout(self.box_layout)
        
        # Create a layout for preview area
        self.preview_main_layout = QHBoxLayout()
        
        # connect buttons to fpExplorer_functions
        self.select_data_btn.clicked.connect(self.select_data_btn_clicked)
        self.settings_btn.clicked.connect(self.settings_btn_clicked)
        self.run_on_batch_btn.clicked.connect(self.run_on_batch_btn_clicked)
        self.open_doc_btn.clicked.connect(self.open_doc_btn_clicked)
        # signal that select data window is open
        self.select_data_open_signal.connect(self.got_select_data_window_open)
        self.got_channel_names.connect(self.show_preview)
        # use stylesheet
        self.area.setStyleSheet(STYLESHEET)

        
    def select_data_btn_clicked(self):
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
        
    def run_on_batch_btn_clicked(self):
        if self.preview_widget != None:
            self.preview_widget.disable_buttons_signal.emit()
            # if trimming was already set in batch prefer that option
            # if not check the most recent trimming in preview
            if "trim_begin" not in self.batch_export_settings_dict:
                try:
                    trim_beginning = int(self.preview_widget.trim_beginning_sec.text())
                    trim_end = int(self.preview_widget.trim_ending_sec.text())
                except:
                    trim_beginning = 0
                    trim_end = 0
                # add to settings dictionary
                self.batch_export_settings_dict["trim_begin"] = str(trim_beginning)
                self.batch_export_settings_dict["trim_end"] = str(trim_end)
            
            # pre check include check box
            # always check if user wnats to include the current subject
            if self.preview_widget.include_cb.isChecked():
                # add that subject to the list of selected_subjects
                if self.preview_widget.options["subject"] not in self.selected_subjects:
                    self.selected_subjects.append(self.preview_widget.options["subject"])
            # if it is unchecked and used to be in selected, remove from selected
            else:
                if self.preview_widget.options["subject"] in self.selected_subjects:
                   self.selected_subjects.remove(self.preview_widget.options["subject"]) 
            # set include to actual status
            if self.preview_widget.subject_comboBox.currentText() in self.selected_subjects:
                self.preview_widget.include_cb.setChecked(True)
            else:
                self.preview_widget.include_cb.setChecked(False)
            
        self.run_on_batch_window = RunOnBatchWindow(self,self.preview_widget,[self.batch_paths_dict,self.selected_subjects,self.batch_export_settings_dict])
        self.run_on_batch_window.got_batch_options_sig.connect(self.got_run_on_batch_settings)
        self.run_on_batch_window.show()
        self.disable_buttons()
        
    def open_doc_btn_clicked(self):
        self.disable_buttons()
        if self.preview_widget != None:
            self.preview_widget.disable_buttons()
        # get path for doc file
        doc_folder = os.path.join(os.path.abspath(os.getcwd()),"Documentation")
        doc_path = os.path.join(doc_folder,"docs.pdf")
        # open file dialog to ask for 
        doc_save_folder = QFileDialog.getExistingDirectory(self,"Select folder to save documentation file")
        # copy pdf with documentation in the selected location
        if len(doc_save_folder) > 0:
            # check if it exists
            new_path = os.path.join(doc_save_folder,"docs.pdf")
            if not os.path.isfile(new_path):
                print(doc_save_folder)
                shutil.copy(doc_path, doc_save_folder)
        self.enable_all_buttons()
        if self.preview_widget != None:
            self.preview_widget.enable_buttons()
        
    def enable_all_buttons(self):
        self.select_data_btn.setEnabled(True)
        self.settings_btn.setEnabled(True)
        self.run_on_batch_btn.setEnabled(True)
        self.open_doc_btn.setEnabled(True)
        QApplication.processEvents() 
        
    def disable_buttons(self):
        self.select_data_btn.setEnabled(False)
        self.settings_btn.setEnabled(False)
        self.run_on_batch_btn.setEnabled(False)
        self.open_doc_btn.setEnabled(False)
        QApplication.processEvents() 
        
    @pyqtSlot()   
    def show_preview(self):
        self.select_data_window.close()
        # enable all buttons
        self.enable_all_buttons()
        # create list of init_params to start previewing
        self.preview_params = [self.select_data_window_content,self.batch_paths_dict]
#        print("preview_params from main\n",self.preview_params)
        if self.preview_widget != None: # clear preview before creating next one
            self.preview_widget.close()    
            self.preview_widget.deleteLater()
            self.preview_widget = None
        # automatically add preview to the main window
        # if the event based experiment type was selected set preview to PreviewEventBasedWidget
        if self.select_data_window_content[0]["event_based"] == True:
            # try reading the first data tank on the list 
            # to verify if there are indeed events
            if fpExplorer_functions.check_events(self.batch_paths_dict[self.select_data_window_content[0]["subject_names"][0]]):
                self.preview_widget = PreviewEventBasedWidget(self,self.preview_params)
                self.preview_widget.done_batch_processing_sig.connect(self.close_batch_window)
                # add widget to dock
                self.preview_widget.setLayout(self.preview_main_layout)
                self.bottom_dock_widget.addWidget(self.preview_widget)
            else:
                self.show_info_dialog("No events found in your data.\nTry choosing Whole Trace Analysis in Select Data window.")
            
        elif self.select_data_window_content[0]["continuous"] == True:
#            self.show_info_dialog("Continuous data not supported yet.")
            self.preview_widget = PreviewContinuousWidget(self,self.preview_params)
            self.preview_widget.done_batch_processing_sig.connect(self.close_batch_window)
            # add widget to dock
            self.preview_widget.setLayout(self.preview_main_layout)
            self.bottom_dock_widget.addWidget(self.preview_widget)
        print(self.select_data_window_content)
        
    @pyqtSlot()   
    def close_batch_window(self):
        if self.run_on_batch_window != None:
            self.run_on_batch_window.close()
        
    @pyqtSlot(list) 
    # receives a single element list with settings from settings window
    def get_settings(self,settings):
        self.settings_dict[0]["downsample"] = settings[0]["downsample"] 
        self.settings_dict[0]["normalization"] = settings[0]["normalization"]
        self.settings_dict[0]["show_norm_as"] = settings[0]["show_norm_as"]
        self.settings_dict[0]["filter"] = settings[0]["filter"]
        print("settings:",self.settings_dict)
        self.enable_all_buttons()
#        QApplication.processEvents() 
        
    # receives batch export settings from run on batch window
    @pyqtSlot(list)   
    def got_run_on_batch_settings(self,batch_export_settings):
#        print(batch_export_settings)
        # update list of selected subjects
        self.selected_subjects = batch_export_settings[0]["batch_subjects"]
        # update export settings dictionary
        self.batch_export_settings_dict = batch_export_settings[0]
        print("Batch settings dict",self.batch_export_settings_dict)
#        self.run_on_batch_window.close()
        self.start_batch_processing_sig.emit()
    
    @pyqtSlot(list) 
    # receives a list with options data read from select data window
    def get_select_data_window_user_input(self,user_input):

        self.select_data_window_content = user_input       
        if self.select_data_window_content[0]["subject_experiment"] == True:
            self.batch_paths_dict = fpExplorer_functions.create_list_of_paths(self.select_data_window_content[0]["main_path"],
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
            valid_subjects,self.batch_paths_dict = fpExplorer_functions.create_list_of_paths_experiment_subjects(self.select_data_window_content[0]["main_path"],
                                                                    self.select_data_window_content[0]["subject_names"],
                                                                    self.select_data_window_content[0]["selected_experiment"])
            # update subject list for that experiment
            self.select_data_window_content[0]["subject_names"] = valid_subjects
            
       
        # ask user for signal and control channels
        data = fpExplorer_functions.get_raw_data(self.batch_paths_dict[self.select_data_window_content[0]["subject_names"][0]])
        self.signal_control_window = SignalControlWindow(self,fpExplorer_functions.get_channel_names(data))
        self.signal_control_window.got_signal_name_sig.connect(self.got_signal_name_sig)
        self.signal_control_window.show()    
            
        
    @pyqtSlot() 
    # receives signal that select folder window is open
    # disables buttons 
    def got_select_data_window_open(self):
        self.disable_buttons()
#        QApplication.processEvents()
        
    @pyqtSlot()
    # receives signal that select folder window is closed
    # enables buttons 
    def got_select_data_window_closed(self):
        print("enable")
        self.enable_all_buttons()
#        QApplication.processEvents()
        
    @pyqtSlot(list) 
    # receives a list with dictionary of signal and control channel names
    def got_signal_name_sig(self,channel_names):
        if len(channel_names) > 0:
            self.select_data_window_content[0]["signal_name"] = channel_names[0]["signal_name"]
            self.select_data_window_content[0]["control_name"] = channel_names[0]["control_name"]
#            print(self.select_data_window_content)
            # clear also previous batch export settings
            self.batch_export_settings_dict = {}
            # clear also previously selected subjects
            self.selected_subjects = []
            self.got_channel_names.emit()
        else:
            if self.settings_window != None:
                self.settings_window.close()
        #    self.show_info_dialog("Please double check your channel names") 
        
    # show info popup
    def show_info_dialog(self, text):
        msgBox = QMessageBox()
        msgBox.setWindowIcon(QtGui.QIcon(ICO))
        msgBox.setText(text)
        msgBox.setWindowTitle("Info!")
        msgBox.setStandardButtons(QMessageBox.Ok)
        msgBox.exec()

    # emit signal when app is closed
    def closeEvent(self, event):  
        self.app_closed.emit() 
        
        
################### end MyMainWidget class   
        
class RunOnBatchWindow(QMainWindow):
    
    got_batch_options_sig = pyqtSignal(list)
    def __init__(self,parent_window,prewiew_window,batch_params):
        super(RunOnBatchWindow, self).__init__()
        self.setWindowTitle("Run On Batch")
        self.setWindowIcon(QtGui.QIcon(ICO))
        self.resize(650,350)
        self.parent_window = parent_window
        self.parent_window.app_closed.connect(self.exit_app)
        self.preview_window = prewiew_window
        # subject:subject group name dictionary
        self.group_names_dict = self.preview_window.group_names_dict
        
        self.all_paths_dict = batch_params[0]
        self.pre_selected_subjects = batch_params[1]
        self.batch_params_dict = batch_params[2]
        self.selected_subjects = []
        self.selected_subjects_group_names = []
        self.updated_group_dict = self.group_names_dict
        
#        print(self.all_paths_dict)
        
        self.dump_path = ""
        self.file_begin = ""
        self.trim_begin = ""
        self.trim_end = ""
        self.export_options = {}
        self.valid_export_data = False
        
        self.bold_stylesheet = "QRadioButton, .QCheckBox, QPushButton {font-weight: bold}"
        
        self.events_present = False
        # check if there are any events
        for v in self.all_paths_dict.values():
            if fpExplorer_functions.check_events(v)==True:
                self.events_present = True
            break
        
        # create gui
        self.setupGUI()
        
    def setupGUI(self):
        # use docks from pyqtgraph
        self.export_loc_area = DockArea()
        self.setCentralWidget(self.export_loc_area)
        
        self.export_loc_dock_widget = Dock("Dock1", size=(1, 1))
        self.export_loc_dock_widget.hideTitleBar()
        self.export_loc_dock_widget.setContentsMargins(10,10,10,10)
        self.export_loc_area.addDock(self.export_loc_dock_widget,'left')
        self.trim_dock_widget = Dock("Dock2", size=(1, 1))
        self.trim_dock_widget.hideTitleBar()
        self.trim_dock_widget.setContentsMargins(10,10,10,10)
        self.export_loc_area.addDock(self.trim_dock_widget,'bottom',self.export_loc_dock_widget)
        self.selected_dock_widget = Dock("Dock3", size=(1,1))
        self.selected_dock_widget.hideTitleBar()
        self.selected_dock_widget.setContentsMargins(10,10,10,10)
        self.export_loc_area.addDock(self.selected_dock_widget, 'bottom',self.trim_dock_widget)  # place at the bottom
        self.options_dock_widget = Dock("Dock4", size=(1,100))
        self.options_dock_widget.hideTitleBar()
        self.options_dock_widget.setContentsMargins(10,10,10,10)
        self.export_loc_area.addDock(self.options_dock_widget, 'right',self.selected_dock_widget)  
        
##########################    
        # EXPORT LOCATION
        self.export_loc_layout = QFormLayout()
        self.export_loc_layout.setContentsMargins(20,20,20,20)
        self.export_loc_widget = QWidget()
        self.select_folder_btn = QPushButton(" Select Export Folder ")
        self.select_folder_btn.setStyleSheet(self.bold_stylesheet)
        self.selected_folder_text = QLineEdit("")
        if "dump_path" in self.batch_params_dict:
            self.selected_folder_text.setText(self.batch_params_dict["dump_path"])
        self.export_loc_layout.addRow(self.select_folder_btn,self.selected_folder_text)
        self.suggested_file_beginning_text = QLineEdit("batch_analysis")
        if "file_begin" in self.batch_params_dict:
            self.suggested_file_beginning_text.setText(self.batch_params_dict["file_begin"])
        self.export_loc_layout.addRow("Suggested file name beginning:",self.suggested_file_beginning_text)
        
##########################
        # TRIMMING
        self.trimming_widget = QWidget()
        self.trimming_layout = QFormLayout()
        self.trimming_layout.setContentsMargins(20,20,20,20)
        self.trim_begin_text = QLineEdit("")
        if "trim_begin" in self.batch_params_dict:
            self.trim_begin_text.setText(self.batch_params_dict["trim_begin"])
        self.trimming_layout.addRow("Trim all beginning (seconds)",self.trim_begin_text)
        self.trim_end_text = QLineEdit("")
        if "trim_end" in self.batch_params_dict:
            self.trim_end_text.setText(self.batch_params_dict["trim_end"])
        self.trimming_layout.addRow("Trim all ending (seconds)",self.trim_end_text)
        
##########################
        # AVAILABLE SUBJECTS
        self.selected_layout = QVBoxLayout()
        self.selected_widget = QWidget()
        self.subjets_layout = QFormLayout() 
        self.subjets_layout.addRow(QLabel("Select Subjects:"),QLabel("Enter Group Name:"))
        self.subjets_layout.addRow(QLabel(""),QLabel(""))
        # find all subjects names
        self.all_subject_names = []
        self.all_group_names = []
        for key in self.all_paths_dict:
            self.all_subject_names.append(key)
        self.subjects_button_group = QButtonGroup()
        self.subjects_button_group.setExclusive(False)
        # create checkboxes for all
        for el in self.all_subject_names:
            radio_btn = QRadioButton(el)
            group_name = ""
            if el in self.group_names_dict:
                group_name = self.group_names_dict[el]
            if el in self.pre_selected_subjects:
                radio_btn.setChecked(True)
            subject_group_name = QLineEdit(group_name)
            self.all_group_names.append(subject_group_name)
            self.subjets_layout.addRow(radio_btn,subject_group_name)
            self.subjects_button_group.addButton(radio_btn)
        self.selected_layout.addLayout(self.subjets_layout)
        self.select_all_subjects_btn = QRadioButton("Select All")
        self.subjets_layout.addRow(QLabel(""),QLabel(""))
        self.select_all_subjects_btn.setStyleSheet(self.bold_stylesheet)
        # self.selected_layout.addWidget(QLabel(""))
        self.selected_layout.addWidget(self.select_all_subjects_btn)
        

#############################
        # TYPE OF ANALYSIS
        self.options_layout = QVBoxLayout()
        self.options_widget = QWidget()
        self.select_analysis_label = QLabel("Select Analysis:")
        self.options_layout.addWidget(self.select_analysis_label)
        self.normalized_cb = QCheckBox("Normalized Data (single subjects only)")
        if "normalized" in self.batch_params_dict:
            self.normalized_cb.setChecked(self.batch_params_dict["normalized"])
        self.perievent_cb = QCheckBox("Perievent Data")
        if "perievent" in self.batch_params_dict:
            self.perievent_cb.setChecked(self.batch_params_dict["perievent"])
        self.spike_cb = QCheckBox("Spikes (single subjects only)")
        if "spikes" in self.batch_params_dict:
            self.spike_cb.setChecked(self.batch_params_dict["spikes"])
        self.options_layout.addWidget(self.normalized_cb)
        if self.events_present == True: # show event option only for eventspreview
            self.options_layout.addWidget(self.perievent_cb)
        self.options_layout.addWidget(self.spike_cb)
        self.options_layout.addWidget(QLabel(""))
        self.export_for_single_subjects_cb = QCheckBox("Export Each Subject Data")
        if "export_for_single_subjects" in self.batch_params_dict:
            self.export_for_single_subjects_cb.setChecked(self.batch_params_dict["export_for_single_subjects"])
        self.export_group_data_cb = QCheckBox("Export Group Analysis Data")
        if "export_group_data" in self.batch_params_dict:
            self.export_group_data_cb.setChecked(self.batch_params_dict["export_group_data"])
        self.export_for_single_subjects_cb.setStyleSheet(self.bold_stylesheet)
        self.export_group_data_cb.setStyleSheet(self.bold_stylesheet)
        self.options_layout.addWidget(self.export_for_single_subjects_cb)
        if self.events_present == True: # group analysis only for events
            self.options_layout.addWidget(self.export_group_data_cb)
        self.options_layout.addWidget(QLabel(""))
        self.run_btn = QPushButton("Run and Export")
        self.run_btn.setStyleSheet(STYLESHEET)
        self.options_layout.addWidget(self.run_btn)
       
        
        # add widgets to docks
        self.export_loc_widget.setLayout(self.export_loc_layout)
        self.export_loc_dock_widget.addWidget(self.export_loc_widget)
        
        self.trimming_widget.setLayout(self.trimming_layout)
        self.trim_dock_widget.addWidget(self.trimming_widget)
        
        self.selected_widget.setLayout(self.selected_layout)
        self.selected_dock_widget.addWidget(self.selected_widget)

        self.options_widget.setLayout(self.options_layout)
        self.options_dock_widget.addWidget(self.options_widget)
        
        self.select_folder_btn.clicked.connect(self.select_folder_btn_clicked)
        self.select_all_subjects_btn.toggled.connect(self.select_all_btn_clicked)
        self.run_btn.clicked.connect(self.run_btn_clicked)
        
#        # use stylesheet
#        self.export_loc_area.setStyleSheet(STYLESHEET)
        
    def select_all_btn_clicked(self):
        if self.select_all_subjects_btn.isChecked():
            for item in self.selected_widget.findChildren(QRadioButton):
                # only subject names
                if item.text() in self.all_subject_names:
                    item.setChecked(True)
        else:
            for item in self.selected_widget.findChildren(QRadioButton):
                # only subject names
                if item.text() in self.all_subject_names:
                    item.setChecked(False)
                    
    def select_folder_btn_clicked(self):
        # set text field to the value of selected path
        self.selected_folder_text.setText(QFileDialog.getExistingDirectory(self,"Choose Folder to Export Data"))
        
        
    def read_and_validate(self):
        ok = True
        # get export location
        self.dump_path = self.selected_folder_text.text()
        self.file_begin = self.suggested_file_beginning_text.text()
        # get trimming
        try:
            trim_start = int(self.trim_begin_text.text())
            trim_end = int(self.trim_end_text.text())
        except:
            trim_start = 0
            trim_end = 0
        # get selected subjects
        # check all radiobuttons
        for item in self.selected_widget.findChildren(QRadioButton):
            # only subject names
            if item.text() in self.all_subject_names:
                group_name = ""
                # now fing matching group name
                for i in range(len(self.all_subject_names)):
                    if item.text() == self.all_subject_names[i]:
                        group_name = self.all_group_names[i].text()
                        self.updated_group_dict[item.text()] = group_name
                if item.isChecked():
                    if item.text() not in self.selected_subjects:
                        self.selected_subjects.append(item.text())
                        self.selected_subjects_group_names.append(group_name)
                        # # now fing matching group name
                        # for i in range(len(self.all_subject_names)):
                        #     if item.text() == self.all_subject_names[i]:
                        #         self.selected_subjects_group_names.append(self.all_group_names[i].text())
                        
        # get what to export               
        normalized = self.normalized_cb.isChecked()
        perievent = self.perievent_cb.isChecked()
        spikes = self.spike_cb.isChecked()
        
        export_for_single_subjects = self.export_for_single_subjects_cb.isChecked()
        export_group_data = self.export_group_data_cb.isChecked()
        
        # add export options to dictionary
        self.export_options["normalized"] = normalized
        self.export_options["perievent"] = perievent
        self.export_options["spikes"] = spikes
        self.export_options["export_for_single_subjects"] = export_for_single_subjects
        self.export_options["export_group_data"] = export_group_data
        self.export_options["dump_path"] = self.dump_path
        self.export_options["file_begin"] = self.file_begin
        self.export_options["trim_begin"] = str(trim_start)
        self.export_options["trim_end"] = str(trim_end)
        self.export_options["batch_subjects"] = self.selected_subjects
        self.export_options["batch_subjects_group_names"] = self.selected_subjects_group_names
        self.export_options["updated_group_names"] = self.updated_group_dict
               
        # if path is empty or no subjects were selected set ok to False
        if len(self.dump_path) == 0 or len(self.selected_subjects) == 0:
            ok = False
        if ((normalized==False and perievent==False and spikes==False) 
            or (export_for_single_subjects==False and export_group_data==False)):
            ok = False
                        
        return ok
        
    def run_btn_clicked(self):
        self.select_folder_btn.setEnabled(False)
        self.run_btn.setEnabled(False)
        self.valid_export_data = self.read_and_validate()
        if self.valid_export_data == False:
            self.show_info_dialog("Could not export data.\nPlease review your export options.")
            self.select_folder_btn.setEnabled(True)
            self.run_btn.setEnabled(True)
        else:
            print("valid export settings")
            # send list with a dictionary
            self.got_batch_options_sig[list].emit([self.export_options])
            
    # emit signal when window is closed
    def closeEvent(self, event):  
        try:
            self.parent_window.enable_all_buttons()
            self.preview_window.enable_buttons_signal.emit()
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
        self.setWindowTitle("Select Data")
        self.setWindowIcon(QtGui.QIcon(ICO))
        self.resize(650,250)
        self.parent_window = parent_window
        self.parent_window.app_closed.connect(self.exit_app)
        
        # list that will contain user input
        self.user_input_list = user_input
        self.experiments_window  = None
        self.signal_control_window = None
        if len(self.user_input_list) > 0:
            self.my_data_folder_path = self.user_input_list[0]['main_path']
            self.event_based = self.user_input_list[0]['event_based']
            self.continuous = self.user_input_list[0]['continuous']
            self.subject_experiment_folder_structure = self.user_input_list[0]['subject_experiment']
            self.experiment_subject_folder_structure = self.user_input_list[0]['experiment_subject']
        else:
            self.my_data_folder_path = None
            self.event_based = True
            self.continuous = False
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
        self.experiment_button_group = QButtonGroup()
        self.experiment_button_group.setExclusive(True)
        self.event_based_radio = QRadioButton("Event Based Experiment")
        self.event_based_radio.setChecked(self.event_based)
        self.continuous_radio = QRadioButton("Whole Trace Analysis")
        self.continuous_radio.setChecked(self.continuous)
        self.experiment_button_group.addButton(self.event_based_radio)
        self.experiment_button_group.addButton(self.continuous_radio)
        self.experiment_layout.addWidget(self.event_based_radio)
        self.experiment_layout.addWidget(self.continuous_radio)
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
        # create a dictionary of user input from select data window and experiments window
        if self.experiment_button_group.checkedId() == -1:
            print("No experiment type selected")
            self.close()
        else:
            if self.experiment_button_group.checkedButton().text() == "Event Based Experiment":
                self.event_based = True
                self.continuous = False
            elif self.experiment_button_group.checkedButton().text() == "Whole Trace Analysis":
                self.event_based = False
                self.continuous = True
            if self.folder_structure_button_group.checkedButton().text() == "Folder Structure: Subject->Experiments":
                self.subject_experiment_folder_structure = True
                self.experiment_subject_folder_structure = False
            elif self.folder_structure_button_group.checkedButton().text() == "Folder Structure: Experiment->Subjects":
                self.subject_experiment_folder_structure = False
                self.experiment_subject_folder_structure = True
                
        data_dict = {"main_path":self.my_data_folder_path,
                     "subject_names":self.subject_folder_names,
                     "selected_experiment":exp,
                     "event_based":self.event_based,
                     "continuous":self.continuous,
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
# CLASS FOR PREVIEW CONTINUOUS EXPERIMENTS  #
##############################################
        
class PreviewContinuousWidget(QWidget):
#    close_wait_sig = pyqtSignal()
    close_export_window_signal = pyqtSignal()
    disable_buttons_signal = pyqtSignal()
    enable_buttons_signal = pyqtSignal()
    done_batch_processing_sig = pyqtSignal()
    def __init__(self,parent_window,init_params):
        super(PreviewContinuousWidget, self).__init__()
        # list of init_params to start previewing
        # first element is a dictionary from select data folder window
        # second element is a dict of paths to all data folders
        self.preview_init_params = init_params
        
        self.parent_window = parent_window
        # list with dictionary with settings from main app window
        self.settings_dict = self.parent_window.settings_dict
        self.parent_window.start_batch_processing_sig.connect(self.start_batch_analysis)
              
        # options read from options section
        self.options = {"subject":self.preview_init_params[0][0]["subject_names"][0],
                        "subject_group_name":"",
                        "plot_raw":True,
                        "plot_downsampled":False,
                        "plot_normalized":False}

        # remember subject(key):group name(value) for later in a form dictionary
        self.group_names_dict = {self.options["subject"]:self.options["subject_group_name"]}
        
        # a dictionary with subject:extracted raw data
        self.raw_data_dict = {}
        # start a widget by reading the first subject data on the list and adding that to dict
        self.get_raw_data(self.preview_init_params[0][0]["subject_names"][0],
                          self.preview_init_params[1][self.preview_init_params[0][0]["subject_names"][0]])
        
        # get last timestamp in seconds
        self.last_raw_ts = fpExplorer_functions.get_last_timestamp(
                                     self.raw_data_dict[self.preview_init_params[0][0]["subject_names"][0]],
                                     self.preview_init_params[0][0]["signal_name"],
                                     )
        # read first subject's frequency and create suggested downsampled rate
        self.current_fs = fpExplorer_functions.get_frequency(self.raw_data_dict[self.preview_init_params[0][0]["subject_names"][0]],self.preview_init_params[0][0]["signal_name"])
        self.suggested_downsample_samples = int(int(self.current_fs)*DEFAULT_DOWNSAMPLE_PCT/100)
        self.settings_dict[0]['downsample'] = self.suggested_downsample_samples
        # keep trimmed data separately with latest trimming settings
        # first element of that list is beginning and end seconds to trim
        # second element is trimmed data dict ts:timestamps, signal:data,control:data
        self.trimmed_raw_data_dict = {}
        # store downsampled data key (subject): dict("ts","signal","control")
        self.downsampled_dict = {}
        # store normalized data key (subject): dict("ts","normalized signal")
        self.normalized_dict = {}
        self.export_window = None
        # initialize export settings
        # export settings
        self.raw_export = False
        self.trimmed_export = False
        self.downsampled_export = False
        self.normalized_export = False
        self.perievent_export = False
        self.spikes_export = False
        self.save_plots = False
        self.export_path = ""
        self.export_begining = ""
        
        self.peak_options_window = None
        self.recent_peak_values = []
        # if user wanted to run spikes on batch
        self.batch_peaks = False
        # all normalized datafor current batch analysis
        self.all_normalized = []
        
        self.disable_buttons_signal.connect(self.disable_buttons)
        self.enable_buttons_signal.connect(self.enable_buttons)
        
        self.bold_label_stylesheet = "QLabel {font-weight: bold}"
        self.bold_cb_stylesheet = "QCheckBox {font-weight: bold}"
        
        # create gui
        self.setupGUI()
        
    def setupGUI(self):
        self.layout = QVBoxLayout()
        self.layout.setContentsMargins(0,0,0,0)
        self.setLayout(self.layout)
        # area on the left with options
        self.splitter = QSplitter()
        self.splitter.setOrientation(Qt.Horizontal)
        self.layout.addWidget(self.splitter)
###########################################
# Options; left side for
        # create widget with options
        self.options_widget = QWidget()
        # create main layout for widget with options
        self.options_main_layout = QVBoxLayout()
        self.options_main_layout.setContentsMargins(10,10,0,10)
        self.options_layout = QFormLayout()
        self.options_layout.setVerticalSpacing(15)
        self.experiment_name_text = QLineEdit(self.preview_init_params[0][0]["selected_experiment"])
        self.experiment_name_text.setReadOnly(True)
        self.experiment_label = QLabel("Experiment")
        self.experiment_label.setStyleSheet(self.bold_label_stylesheet)
        self.options_layout.addRow(self.experiment_label,self.experiment_name_text)
        # add drop down menu for subjects from the list
        self.subject_comboBox = QComboBox()
        self.subject_comboBox.addItems(self.preview_init_params[0][0]["subject_names"])
        # adjust trimming when subject changed
        self.subject_comboBox.currentIndexChanged.connect(self.on_subject_change)
        self.subject_label = QLabel("Subject")
        self.subject_label.setStyleSheet(self.bold_label_stylesheet)
        self.options_layout.addRow(self.subject_label, self.subject_comboBox)
        self.subject_group_name = QLineEdit("")
        self.options_layout.addRow("Subject group name", self.subject_group_name)
        self.last_raw_ts = round(self.last_raw_ts,2)
        self.last_raw_ts_text = QLineEdit(str(self.last_raw_ts))
        self.last_raw_ts_text.setReadOnly(True)
        self.options_layout.addRow("Total seconds", self.last_raw_ts_text)
        self.trim_beginning_sec = QLineEdit("0")
        self.trim_beginning_sec.setValidator(QtGui.QIntValidator())
        self.trim_ending_sec = QLineEdit("0")
        self.trim_ending_sec.setValidator(QtGui.QIntValidator())
        self.options_layout.addRow("Trim beginning (seconds)",self.trim_beginning_sec)
        self.options_layout.addRow("Trim ending (seconds)",self.trim_ending_sec)
        self.downsample_cb = QCheckBox("Downsample")
        self.perform_label = QLabel("Perform:")
        self.perform_label.setStyleSheet(self.bold_label_stylesheet)
        self.downsample_cb.setStyleSheet(self.bold_cb_stylesheet)
        self.options_layout.addRow(self.perform_label,self.downsample_cb)
        self.normalize_cb = QCheckBox("Normalize Downsampled")
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
        # self.options_layout.addRow("Show on the plot:",self.raw_plot_cb)
        self.show_label = QLabel("Show:")
        self.show_label.setStyleSheet(self.bold_label_stylesheet)
        self.options_layout.addRow(self.show_label,self.downsampled_plot_cb)
        # self.options_layout.addRow("",self.downsampled_plot_cb)
        self.options_layout.addRow("",self.separate_signal_contol_cb)
        self.options_layout.addRow("",self.normalized_plot_cb)
              
        # nest the inner layout into the outer layout
        self.options_main_layout.addLayout(self.options_layout)
        
        # add apply button
        self.apply_btn = QPushButton("Plot")
        self.options_main_layout.addWidget(self.apply_btn)
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
            + "padding: 7px;\n" \
            + "margin-left:20px;\n" \
            + "margin-right:20px;\n" \
            + "margin-top:10px;\n" \
            + "margin-bottom:10px;\n" \
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
        self.canvas = MplCanvas(self, width=12, height=8, dpi=100)
        # Create toolbar, passing canvas as first parament, parent (self, the MainWindow) as second.
        self.toolbar = NavigationToolbar(self.canvas, self)
        # add to layout
        self.plot_area_widget_layout.addWidget(self.toolbar)
        self.plot_area_widget_layout.addWidget(self.canvas)
        self.plot_area_widget.setLayout(self.plot_area_widget_layout)
        self.splitter2.addWidget(self.plot_area_widget)
        
        # plot initial plot of raw data of the first subject
        self.plot_bare_raw_data(self.preview_init_params[0][0]["subject_names"][0],
                                self.raw_data_dict[self.preview_init_params[0][0]["subject_names"][0]])
        
##############################################
# Previous/Next buttons on the bottom        
        # add next/previous buttons at the bottom
        self.preview_buttons_layout = QHBoxLayout()
        self.preview_buttons_layout.setAlignment(Qt.AlignRight)
        self.preview_buttons_widget = QWidget()
        self.include_cb = QCheckBox("Include Subject")
        self.next_btn = QPushButton("Next")
        self.previous_btn = QPushButton("Previous")
#        self.previous_btn.setFixedWidth(120)
#        self.next_btn.setFixedWidth(120)
        self.preview_buttons_layout.addWidget(self.include_cb)
        self.preview_buttons_layout.addWidget(self.previous_btn)
        self.preview_buttons_layout.addWidget(self.next_btn)
        self.preview_buttons_widget.setLayout(self.preview_buttons_layout)
        self.splitter2.addWidget(self.preview_buttons_widget)
        
        # connect buttons to fpExplorer_functions
        self.downsample_cb.stateChanged.connect(self.adjust_downsampled_plot_cb)
        self.normalize_cb.stateChanged.connect(self.adjust_normalized_plot_cb)
        self.downsampled_plot_cb.stateChanged.connect(self.adjuct_separate)
        self.raw_plot_cb.stateChanged.connect(self.adjust_separate_control_signal_cb)
        self.separate_signal_contol_cb.stateChanged.connect(self.adjust_downsampled2separate)
        self.apply_btn.clicked.connect(self.apply_btn_clicked)
        self.find_peaks_btn.clicked.connect(self.peaks_btn_clicked)
        self.export_data_btn.clicked.connect(self.export_data_btn_clicked)
        self.check_poly_btn.clicked.connect(self.check_poly_btn_clicked)
        self.next_btn.clicked.connect(self.next_btn_clicked)
        self.previous_btn.clicked.connect(self.previous_btn_clicked)
        
    @pyqtSlot()   
    def disable_buttons(self):
        self.apply_btn.setEnabled(False)
        self.find_peaks_btn.setEnabled(False)
        self.export_data_btn.setEnabled(False)
        self.check_poly_btn.setEnabled(False)
        self.next_btn.setEnabled(False)
        self.previous_btn.setEnabled(False)
        self.parent_window.disable_buttons()
        QApplication.processEvents() 
    
    @pyqtSlot()     
    def enable_buttons(self):
        self.apply_btn.setEnabled(True)
        self.find_peaks_btn.setEnabled(True)
        self.export_data_btn.setEnabled(True)
        self.check_poly_btn.setEnabled(True)
        self.next_btn.setEnabled(True)
        self.previous_btn.setEnabled(True)
        self.parent_window.enable_all_buttons()
        QApplication.processEvents() 
        
    # function to read options for plotting on the left side of the plot    
    def read_options(self):
        # set current subject
        self.options["subject"] = self.subject_comboBox.currentText()
        # if subject was not previewed yet, read data and add to dict
        if self.options["subject"] not in self.raw_data_dict:
            self.get_raw_data(self.options["subject"],
                          self.parent_window.preview_params[1][self.options["subject"]])
        self.options["subject_group_name"] = self.subject_group_name.text()
        # update value
        self.group_names_dict[self.options["subject"]] = self.options["subject_group_name"]
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
            self.trimmed_raw_data_dict[self.options["subject"]] = [(trim_beginning,trim_end),fpExplorer_functions.trim_raw_data(self.raw_data_dict[self.options["subject"]],
                                                                                                                      self.preview_init_params[0][0]["signal_name"],
                                                                                                                      self.preview_init_params[0][0]["control_name"],
                                                                                                                      trim_beginning,
                                                                                                                      trim_end
                                                                                                                      )]
        else:   # if already in trimmed check previous trimming settings
            trimmed = self.trimmed_raw_data_dict[self.options["subject"]]
            begin,end = trimmed[0]
            if trim_beginning != begin or trim_end != end:
                # if new trimming params, replace the lists in dict
                self.trimmed_raw_data_dict[self.options["subject"]] = [(trim_beginning,trim_end),fpExplorer_functions.trim_raw_data(self.raw_data_dict[self.options["subject"]],
                                                                                                                      self.preview_init_params[0][0]["signal_name"],
                                                                                                                      self.preview_init_params[0][0]["control_name"],
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
            
            # get trimmed timestamps and update new session duration
            trimmed_ts = self.trimmed_raw_data_dict[self.options["subject"]][1]["ts"]
            total_trimmed = trimmed_ts[-1]-trimmed_ts[0]
            # start time from zero
            ts_reset = [i*total_trimmed/len(trimmed_ts) for i in range(len(trimmed_ts))]
            last_ts = round(ts_reset[-1],2)
            self.last_raw_ts_text.setText(str(last_ts))
            
            # if downsample was selected
            if self.downsample_cb.isChecked() or self.downsampled_export == True or self.save_plots == True:
                # add to downdampled dict
                self.downsampled_dict[self.options["subject"]] = fpExplorer_functions.downsample_tdt(self.trimmed_raw_data_dict[self.options["subject"]][1],
                                            self.settings_dict[0]["downsample"])
            # if normalize was selected
            if self.normalize_cb.isChecked() or self.normalized_export == True:
                # always downsample first before normalizing
                # downsample
                self.downsampled_dict[self.options["subject"]] = fpExplorer_functions.downsample_tdt(self.trimmed_raw_data_dict[self.options["subject"]][1],
                                                                                        self.settings_dict[0]["downsample"])
                # check settings for method to normalize
                if self.settings_dict[0]["normalization"] == "Modified Polynomial Fitting":                   
                    # normalize from downsampled
                    self.normalized_dict[self.options["subject"]] = fpExplorer_functions.normalize_dff(self.raw_data_dict[self.options["subject"]],
                                                                                        self.downsampled_dict[self.options["subject"]],
                                                                                        self.settings_dict[0]["show_norm_as"],
                                                                                        self.settings_dict[0]["filter"],
                                                                                        self.settings_dict[0]["filter_fraction"])
                if self.settings_dict[0]["normalization"] == "Standard Polynomial Fitting":                
                    # normalize from downsampled
                    self.normalized_dict[self.options["subject"]] = fpExplorer_functions.normalize_pMat(self.raw_data_dict[self.options["subject"]],
                                                                                        self.downsampled_dict[self.options["subject"]],
                                                                                        self.settings_dict[0]["show_norm_as"],
                                                                                        self.settings_dict[0]["filter"],
                                                                                        self.settings_dict[0]["filter_fraction"])
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
            
        
    # function to plot with user selected options       
    def apply_btn_clicked(self):
        self.disable_buttons_signal.emit()
        if self.read_options() == True:
            # plot
            # if just raw was checked
            if (self.options["plot_raw"]==True and self.options["plot_separate"]==False
                and self.options["plot_downsampled"]==False and self.options["plot_normalized"]==False):
                fpExplorer_functions.plot_trimmed(self.canvas,
                                        self.options["subject"],
                                        self.trimmed_raw_data_dict[self.options["subject"]][1],
                                        self.raw_export,
                                        self.save_plots,
                                        (self.export_path,self.export_begining),
                                        self.preview_init_params[0][0]["signal_name"],
                                        self.preview_init_params[0][0]["control_name"])
            # if only downsampled was checked
            elif (self.options["plot_raw"]==False and self.options["plot_separate"]==False and self.options["subject"] in self.downsampled_dict and self.options["plot_downsampled"]==True 
                and self.options["plot_normalized"]==False):
                fpExplorer_functions.plot_downsampled_alone(self.canvas,
                                            self.options,
                                            self.downsampled_dict[self.options["subject"]],
                                            self.downsampled_export,
                                            self.save_plots,
                                            (self.export_path,self.export_begining),
                                            self.settings_dict,
                                            self.preview_init_params[0][0]["signal_name"],
                                            self.preview_init_params[0][0]["control_name"])
            # if only normalized was checked
            elif (self.options["plot_raw"]==False and self.options["plot_separate"]==False and self.options["plot_downsampled"]==False
                and self.options["subject"] in self.normalized_dict and self.options["plot_normalized"]==True):
                fpExplorer_functions.plot_normalized_alone(self.canvas,
                                                self.options,
                                                self.normalized_dict[self.options["subject"]],
                                                self.normalized_export,
                                                self.save_plots,
                                                (self.export_path,self.export_begining),
                                                self.settings_dict)
            # if downsampled and normalized were selected
            elif (self.options["plot_raw"]==False and self.options["plot_separate"]==False and self.options["subject"] in self.downsampled_dict and self.options["plot_downsampled"]==True
                and self.options["subject"] in self.normalized_dict and self.options["plot_normalized"]==True):
                fpExplorer_functions.plot_downsampled_and_normalized_alone(self.canvas,
                                            self.options["subject"],
                                            self.downsampled_dict[self.options["subject"]],
                                            self.normalized_dict[self.options["subject"]],
                                            self.settings_dict[0]["show_norm_as"])
            # if user wants signal and control separately
            elif (self.options["plot_raw"]==False and self.options["plot_separate"]==True 
                and self.options["plot_downsampled"]==False and self.options["plot_normalized"]==False):
                # if there is downsampled data, plot that, if not raw
                downsampled_to_plot = {}
                if self.options["subject"] in self.downsampled_dict:
                    downsampled_to_plot = self.downsampled_dict[self.options["subject"]]
                    self.show_info_dialog("Separate signal and control plots\nshow downsampled data.")
                
                fpExplorer_functions.plot_separate_only(self.canvas,
                                            self.options["subject"],
                                            self.trimmed_raw_data_dict[self.options["subject"]][1],
                                            downsampled_to_plot)
            # if user wants signal and control separately and specificaly chose downsample
            elif (self.options["plot_raw"]==False and self.options["plot_separate"]==True 
                and self.options["subject"] in self.downsampled_dict and self.options["plot_downsampled"]==True 
                and self.options["plot_normalized"]==False):
                self.show_info_dialog("Separate signal and control plots\nwill show downsampled data.")
                fpExplorer_functions.plot_separate_only(self.canvas,
                                            self.options["subject"],
                                            self.trimmed_raw_data_dict[self.options["subject"]][1],
                                            self.downsampled_dict[self.options["subject"]])
            # if user wants signal and control separately and normalized
            elif (self.options["plot_raw"]==False and self.options["plot_separate"]==True 
                and self.options["plot_downsampled"]==False 
                and self.options["subject"] in self.normalized_dict and self.options["plot_normalized"]==True):
                # if there is downsampled data, plot that, if not raw
                downsampled_to_plot = {}
                if self.options["subject"] in self.downsampled_dict:
                    downsampled_to_plot = self.downsampled_dict[self.options["subject"]]
                    self.show_info_dialog("Separate signal and control plots\nwill show downsampled data.")
                fpExplorer_functions.plot_separate_with_normalized(self.canvas,
                                            self.options["subject"],
                                            self.trimmed_raw_data_dict[self.options["subject"]][1],
                                            downsampled_to_plot,
                                            self.normalized_dict[self.options["subject"]],
                                            self.settings_dict[0]["show_norm_as"])
            # if user wants signal and control separately and normalized and specificaly chose downsample
            elif (self.options["plot_raw"]==False and self.options["plot_separate"]==True 
                and self.options["plot_downsampled"]==True and  self.options["subject"] in self.downsampled_dict
                and self.options["subject"] in self.normalized_dict and self.options["plot_normalized"]==True):
                self.show_info_dialog("Separate signal and control plots\nwill show downsampled data.")
                fpExplorer_functions.plot_separate_with_normalized(self.canvas,
                                            self.options["subject"],
                                            self.trimmed_raw_data_dict[self.options["subject"]][1],
                                            self.downsampled_dict[self.options["subject"]],
                                            self.normalized_dict[self.options["subject"]],
                                            self.settings_dict[0]["show_norm_as"])
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
        
    def check_poly_btn_clicked(self):
        self.disable_buttons_signal.emit()
        if self.normalize_cb.isChecked():
            if self.read_options() == True:
                fpExplorer_functions.show_polynomial_fitting(self.canvas,
                                                self.settings_dict[0], 
                                                self.downsampled_dict[self.options["subject"]],
                                                self.preview_init_params[0][0]["signal_name"],
                                                self.preview_init_params[0][0]["control_name"],
                                                self.options["subject"],
                                                False,
                                                "")
        else: # prompt user to check normalize first
            self.show_info_dialog("Please check Normalize check box first.")
        self.enable_buttons_signal.emit()
        
    def peaks_btn_clicked(self):
        self.disable_buttons_signal.emit()
        self.peak_options_window = PeakSettingsWindow(self.parent_window,
                                                      self.options["subject"],
                                                      self.recent_peak_values, 
                                                      (self.export_path,self.export_begining),
                                                      self.batch_peaks)
        self.peak_options_window.got_peak_options_sig.connect(self.peak_options_received)
        self.peak_options_window.show()
      
    @pyqtSlot(list) 
    def peak_options_received(self,peak_options):
        # save peak parameters 
        self.recent_peak_values = peak_options[0]
        export_path,file_begin = peak_options[1]
        if self.batch_peaks == True:
            if "spikes" in self.parent_window.batch_export_settings_dict:
                if self.parent_window.batch_export_settings_dict["spikes"] == True:
                    if self.parent_window.batch_export_settings_dict["export_for_single_subjects"] == True:
                        for subject in self.parent_window.batch_export_settings_dict["batch_subjects"]:
                            # create subfolder with subject name
                            subject_subfolder = os.path.join(self.parent_window.batch_export_settings_dict["dump_path"],subject)
                            if not os.path.exists(subject_subfolder):
                                try:
                                    os.mkdir(subject_subfolder)
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
                                except:
                                    self.show_info_dialog("Problem creating subfolder")
                            else: # if subfolder already exists
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
                    if self.parent_window.batch_export_settings_dict["export_group_data"] == True:
                        pass # single subjects only
                        # spikes for averaged signal
#                        try:
                        # fpExplorer_functions.get_batch_spikes(self.canvas,
                        #                                        self.recent_peak_values,
                        #                                        self.all_normalized,
                        #                                        self.settings_dict,
                        #                                        self.parent_window.batch_export_settings_dict["spikes"],
                        #                                        (self.parent_window.batch_export_settings_dict["dump_path"],self.parent_window.batch_export_settings_dict["file_begin"]))
#                        except:
#                            self.show_info_dialog("Could not calculate spikes.\nTry with different parameters.")
                # close popup window
                self.peak_options_window.close()
                self.reset_export_settings()
                # reset batch peaks when done
                self.batch_peaks = False
                if self.parent_window.run_on_batch_window != None:
                    self.parent_window.run_on_batch_window.close()
        else: # if it is not for batch analysis
            # if user entered path or if path was already there, update
            if len(export_path) > 0:
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
                self.trimmed_raw_data_dict[self.options["subject"]] = [(trim_beginning,trim_end),fpExplorer_functions.trim_raw_data(self.raw_data_dict[self.options["subject"]],
                                                                                                                          self.preview_init_params[0][0]["signal_name"],
                                                                                                                          self.preview_init_params[0][0]["control_name"],
                                                                                                                          trim_beginning,
                                                                                                                          trim_end
                                                                                                                          )]
            else:   # if already in trimmed check previous trimming settings
                trimmed = self.trimmed_raw_data_dict[self.options["subject"]]
                begin,end = trimmed[0]
                if trim_beginning != begin or trim_end != end:
                    # if new trimming params, replace the lists in dict
                    self.trimmed_raw_data_dict[self.options["subject"]] = [(trim_beginning,trim_end),fpExplorer_functions.trim_raw_data(self.raw_data_dict[self.options["subject"]],
                                                                                                                          self.preview_init_params[0][0]["signal_name"],
                                                                                                                          self.preview_init_params[0][0]["control_name"],
                                                                                                                          trim_beginning,
                                                                                                                          trim_end
                                                                                                                          )]
            # always downsample first before normalizing
            # downsample
            self.downsampled_dict[self.options["subject"]] = fpExplorer_functions.downsample_tdt(self.trimmed_raw_data_dict[self.options["subject"]][1],
                                                                                         self.settings_dict[0]["downsample"])
            # check settings for method to normalize
            if self.settings_dict[0]["normalization"] == "Modified Polynomial Fitting":                   
                # normalize from downsampled
                self.normalized_dict[self.options["subject"]] = fpExplorer_functions.normalize_dff(self.raw_data_dict[self.options["subject"]],
                                                                                        self.downsampled_dict[self.options["subject"]],
                                                                                        self.settings_dict[0]["show_norm_as"],
                                                                                        self.settings_dict[0]["filter"],
                                                                                        self.settings_dict[0]["filter_fraction"])
            if self.settings_dict[0]["normalization"] == "Standard Polynomial Fitting":                
                # normalize from downsampled
                self.normalized_dict[self.options["subject"]] = fpExplorer_functions.normalize_pMat(self.raw_data_dict[self.options["subject"]],
                                                                                        self.downsampled_dict[self.options["subject"]],
                                                                                        self.settings_dict[0]["show_norm_as"],
                                                                                        self.settings_dict[0]["filter"],
                                                                                        self.settings_dict[0]["filter_fraction"])
    
                
            # close popup window
            self.peak_options_window.close()
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
            self.reset_export_settings()
            self.enable_buttons_signal.emit()
            if self.export_window != None:
                self.close_export_window_signal.emit()
          
        

        
    def export_data_btn_clicked(self):
        # pass a list with first element=subject
        # second element is a list with one dict with init params (i.e. continuous, event_based...)
        self.export_window = ExportDataWindow(self.parent_window,self,[self.options["subject"],
                                                    self.preview_init_params,self.export_path,self.export_begining])
        self.export_window.got_export_selection_sig.connect(self.export_options_received)
        self.export_window.show()
        self.disable_buttons_signal.emit()
        
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
                self.export_window = ExportDataWindow(self,[self.subject_comboBox.currentText(),
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
            self.spikes_export = export_options[0]["spikes"]
            self.save_plots = export_options[0]["save_plot"]
            self.export_data()
        
    # function to plot with user selected options       
    def export_data(self):
        if self.read_options() == True:
            # save fittings if save all plots was selected
            if self.save_plots == True:
                fpExplorer_functions.show_polynomial_fitting(self.canvas,
                                                self.settings_dict[0], 
                                                self.downsampled_dict[self.options["subject"]],
                                                self.preview_init_params[0][0]["signal_name"],
                                                self.preview_init_params[0][0]["control_name"],
                                                self.options["subject"],
                                                self.save_plots,
                                                self.export_path)
            # raw was checked
            if self.raw_export == True:
                self.plot_bare_raw_data(self.options["subject"],self.raw_data_dict[self.options["subject"]])
            if self.downsampled_export == True:           
                fpExplorer_functions.plot_downsampled_alone(self.canvas,
                                            self.options,
                                            self.downsampled_dict[self.options["subject"]],
                                            self.downsampled_export,
                                            self.save_plots,
                                            (self.export_path,self.export_begining),
                                            self.settings_dict,
                                            self.preview_init_params[0][0]["signal_name"],
                                            self.preview_init_params[0][0]["control_name"])
            if self.normalized_export == True:
                fpExplorer_functions.plot_normalized_alone(self.canvas,
                                                self.options,
                                                self.normalized_dict[self.options["subject"]],
                                                self.normalized_export,
                                                self.save_plots,
                                                (self.export_path,self.export_begining),
                                                self.settings_dict)
            if self.spikes_export == True:
                self.peaks_btn_clicked()
            else:   # if perievent export was not selected, reset export settings, otherwise, reset it after perievent export   
                self.reset_export_settings()

                self.close_export_window_signal.emit()
                self.enable_buttons_signal.emit()
        
        
    # add raw data structure and subject to self.raw_data dictionary   
    def get_raw_data(self,subject,my_path):
        self.raw_data_dict[subject] = fpExplorer_functions.get_raw_data(my_path)
        
    # when app starts plot just raw data  
    def plot_bare_raw_data(self,subject,raw_data):
        fpExplorer_functions.plot_raw(self.canvas,
                           subject,
                           raw_data,
                           self.preview_init_params[0][0]["signal_name"],
                           self.preview_init_params[0][0]["control_name"],
                           self.raw_export,
                           self.save_plots,
                           (self.export_path,self.export_begining)
                           )
     
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
        self.check_trimming()
        # set include to actual status
        if self.subject_comboBox.currentText() in self.parent_window.selected_subjects:
            self.include_cb.setChecked(True)
        else:
            self.include_cb.setChecked(False)

            
    def on_subject_change(self):
        self.options_pre_check()
        # update last timestamp
        # get last timestamp in seconds
        self.last_raw_ts = fpExplorer_functions.get_last_timestamp(
                                     self.raw_data_dict[self.options["subject"]],
                                     self.preview_init_params[0][0]["signal_name"],
                                     )
        # read new subject's frequency and update suggested downsampled rate (only if it is different)
        new_fs = fpExplorer_functions.get_frequency(self.raw_data_dict[self.preview_init_params[0][0]["subject_names"][0]],self.preview_init_params[0][0]["signal_name"])
        if new_fs != self.current_fs:
            self.suggested_downsample_samples = int(int(self.current_fs)*DEFAULT_DOWNSAMPLE_PCT/100)
            self.settings_dict[0]['downsample'] = self.suggested_downsample_samples
            self.settings[0]["entered_downsample"] = None
        if self.subject_comboBox.currentText() in self.group_names_dict:
            self.subject_group_name.setText(self.group_names_dict[self.subject_comboBox.currentText()])
        else:
            self.subject_group_name.setText("")
        # clear long lasting check boxes to start faster
        self.downsample_cb.setChecked(False)
        self.normalize_cb.setChecked(False)
        # plot raw by default
        self.raw_plot_cb.setChecked(True)
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
        
    def reset_export_settings(self):
        self.raw_export = False
        self.trimmed_export = False
        self.downsampled_export = False
        self.normalized_export = False
        self.perievent_export = False
        self.save_plots = False
        
    @pyqtSlot()     
    def start_batch_analysis(self):
        # update group names
        self.group_names_dict = self.parent_window.batch_export_settings_dict["updated_group_names"]
        # update group name field of the current subject
        self.subject_group_name.setText(self.group_names_dict[self.subject_comboBox.currentText()])
        # update current trimming
        self.trim_beginning_sec.setText(self.parent_window.batch_export_settings_dict["trim_begin"])
        self.trim_ending_sec.setText(self.parent_window.batch_export_settings_dict["trim_end"])
        # set include to actual status
        if self.subject_comboBox.currentText() in self.parent_window.selected_subjects:
            self.include_cb.setChecked(True)
        else:
            self.include_cb.setChecked(False)
        # disable buttons while processing
        self.disable_buttons_signal.emit()
        print("Started batch processing")
        new_trim_start = int(self.parent_window.batch_export_settings_dict["trim_begin"])
        new_trim_end = int(self.parent_window.batch_export_settings_dict["trim_end"])
        
        
        # create a popup window that will be on during processing?
        # disable all buttons?
        if self.parent_window.batch_export_settings_dict["normalized"] == True or self.parent_window.batch_export_settings_dict["spikes"] == True:
            # create normalized data for most recent settings
            for subject in self.parent_window.batch_export_settings_dict["batch_subjects"]:
                # create separate options dictionary for batch analysis
                self.batch_options_dict = {"subject":subject,
                                            "subject_group_name":self.group_names_dict[subject]}
                # if subject was not previewed yet, read data and add to dict
                if subject not in self.raw_data_dict:
#                    print(self.preview_init_params[1])
                    self.get_raw_data(subject, self.parent_window.preview_params[1][subject])
                    
                # trim to batch export settings
                # key is the subject name, value is a list
                # first element of that list is beginning and end seconds to trim
                # second element is trimmed data dict ts:timestamps, signal:data,control:data
                self.trimmed_raw_data_dict[subject] = [(new_trim_start,new_trim_end),fpExplorer_functions.trim_raw_data(self.raw_data_dict[subject],
                                                                                  self.preview_init_params[0][0]["signal_name"],
                                                                                  self.preview_init_params[0][0]["control_name"],
                                                                                  new_trim_start,
                                                                                  new_trim_end
                                                                                  )]
                
                # always downsample first before normalizing
                # downsample
                self.downsampled_dict[subject] = fpExplorer_functions.downsample_tdt(self.trimmed_raw_data_dict[subject][1],
                                                                             self.settings_dict[0]["downsample"])
                # check settings for method to normalize
                if self.settings_dict[0]["normalization"] == "Modified Polynomial Fitting":                   
                    # normalize from downsampled
                    self.normalized_dict[subject] = fpExplorer_functions.normalize_dff(self.raw_data_dict[subject],
                                                                                            self.downsampled_dict[subject],
                                                                                            self.settings_dict[0]["show_norm_as"],
                                                                                            self.settings_dict[0]["filter"],
                                                                                            self.settings_dict[0]["filter_fraction"])
                if self.settings_dict[0]["normalization"] == "Standard Polynomial Fitting":                
                    # normalize from downsampled
                    self.normalized_dict[subject] = fpExplorer_functions.normalize_pMat(self.raw_data_dict[subject],
                                                                                            self.downsampled_dict[subject],
                                                                                            self.settings_dict[0]["show_norm_as"],
                                                                                            self.settings_dict[0]["filter"],
                                                                                            self.settings_dict[0]["filter_fraction"])
                if self.parent_window.batch_export_settings_dict["export_for_single_subjects"] == True:
                    # create subfolder with subject name
                    subject_subfolder = os.path.join(self.parent_window.batch_export_settings_dict["dump_path"],subject)
                    if not os.path.exists(subject_subfolder):
                        try:
                            os.mkdir(subject_subfolder)
                            # save also polynomial fitting by default
                            fpExplorer_functions.show_polynomial_fitting(self.canvas,
                                                    self.settings_dict[0], 
                                                    self.downsampled_dict[subject],
                                                    self.preview_init_params[0][0]["signal_name"],
                                                    self.preview_init_params[0][0]["control_name"],
                                                    subject,
                                                    True,
                                                    subject_subfolder)
                            
                            if self.parent_window.batch_export_settings_dict["normalized"] == True:
                                fpExplorer_functions.plot_normalized_alone(self.canvas,
                                                            self.batch_options_dict,
                                                            self.normalized_dict[subject],
                                                            True,
                                                            True,
                                                            (subject_subfolder,self.parent_window.batch_export_settings_dict["file_begin"]),
                                                            self.settings_dict)
                        except:
                            self.show_info_dialog("Problem creating subfolder")
                    else: # if subfolder already folder exists
                        fpExplorer_functions.show_polynomial_fitting(self.canvas,
                                                    self.settings_dict[0], 
                                                    self.downsampled_dict[subject],
                                                    self.preview_init_params[0][0]["signal_name"],
                                                    self.preview_init_params[0][0]["control_name"],
                                                    subject,
                                                    True,
                                                    subject_subfolder)
                            
                        if self.parent_window.batch_export_settings_dict["normalized"] == True:
                            fpExplorer_functions.plot_normalized_alone(self.canvas,
                                                            self.batch_options_dict,
                                                            self.normalized_dict[subject],
                                                            True,
                                                            True,
                                                            (subject_subfolder,self.parent_window.batch_export_settings_dict["file_begin"]),
                                                            self.settings_dict)
                    if self.parent_window.batch_export_settings_dict["spikes"] == True:
                        self.spikes_export = True
                        self.save_plots = True
                    
            ### end for loop for each subject   
            if self.parent_window.batch_export_settings_dict["export_group_data"] == True: 
                # average signal
                # create a list of all normalized signals
                # self.all_normalized = [(subject,self.normalized_dict[subject]) for subject in self.parent_window.batch_export_settings_dict["batch_subjects"]]
    #            print(all_normalized)
                if self.parent_window.batch_export_settings_dict["normalized"] == True:
                    pass # only for single subjects
                    # fpExplorer_functions.get_batch_normalized(self.canvas,
                    #                                    self.all_normalized,
                    #                                    self.settings_dict,
                    #                                    self.parent_window.batch_export_settings_dict["normalized"],
                    #                                    (self.parent_window.batch_export_settings_dict["dump_path"],self.parent_window.batch_export_settings_dict["file_begin"]))
                if self.parent_window.batch_export_settings_dict["spikes"] == True:
                    # set batch peaks to true
                    self.batch_peaks = True
            if self.parent_window.batch_export_settings_dict["spikes"] == True:
                # set batch peaks to true
                self.batch_peaks = True
                self.peaks_btn_clicked()
            else:
                print("Done batch processing")
                self.done_batch_processing_sig.emit()
                self.enable_buttons_signal.emit()
                
#        print("Done batch processing")
#        self.done_batch_processing_sig.emit()
#        self.enable_buttons_signal.emit()
        
    # show info popup
    def show_info_dialog(self, text):
        msgBox = QMessageBox()
        msgBox.setWindowIcon(QtGui.QIcon(ICO))
        msgBox.setText(text)
        msgBox.setWindowTitle("Info!")
        msgBox.setStandardButtons(QMessageBox.Ok)
        msgBox.exec()
        
######################## end PreviewContinuousWidget class
        
     
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
        self.parent_window.start_batch_processing_sig.connect(self.start_batch_analysis)
              
        # options read from options section
        self.options = {"subject":self.preview_init_params[0][0]["subject_names"][0],
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
        # start a widget by reading the first subject data on the list and adding that to dict
        self.get_raw_data(self.preview_init_params[0][0]["subject_names"][0],
                          self.preview_init_params[1][self.preview_init_params[0][0]["subject_names"][0]])
        # get last timestamp in seconds
        self.last_raw_ts = fpExplorer_functions.get_last_timestamp(
                                     self.raw_data_dict[self.preview_init_params[0][0]["subject_names"][0]],
                                     self.preview_init_params[0][0]["signal_name"],
                                     )
        # read first subject's frequency and create suggested downsampled rate
        self.current_fs = fpExplorer_functions.get_frequency(self.raw_data_dict[self.preview_init_params[0][0]["subject_names"][0]],self.preview_init_params[0][0]["signal_name"])
        self.suggested_downsample_samples = int(int(self.current_fs)*DEFAULT_DOWNSAMPLE_PCT/100)
        self.settings_dict[0]['downsample'] = self.suggested_downsample_samples
        self.suggested_downsample_rate = int(int(self.current_fs)/self.suggested_downsample_samples)
        self.settings_dict[0]["entered_downsample"] = self.suggested_downsample_rate
        # print("Current suggested downsample rate:",self.settings_dict[0]['downsample'])
        # keep trimmed data separately with latest trimming settings
        # first element of that list is beginning and end seconds to trim
        # second element is trimmed data dict ts:timestamps, signal:data,control:data
        self.trimmed_raw_data_dict = {}
        
        # create a list of available events for current subject
        self.events_from_current_subject = fpExplorer_functions.get_events(self.raw_data_dict[self.preview_init_params[0][0]["subject_names"][0]])
            
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
        # all normalized datafor current batch analysis
        self.all_normalized = []
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
###########################################
# Options; left side for
        # create widget with options
        self.options_widget = QWidget()
        # create main layout for widget with options
        self.options_main_layout = QVBoxLayout()
        self.options_main_layout.setContentsMargins(10,10,0,10)
        self.options_layout = QFormLayout()
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
        self.subject_group_name = QLineEdit("")
        self.options_layout.addRow("Subject group name", self.subject_group_name)
        self.last_raw_ts = round(self.last_raw_ts,2)
        self.last_raw_ts_text = QLineEdit(str(self.last_raw_ts))
        self.last_raw_ts_text.setReadOnly(True)
        self.options_layout.addRow("Total seconds", self.last_raw_ts_text)
        # adjust trimming when subject changed
        self.trim_beginning_select = QComboBox()
        self.trim_beginning_select.addItems(["Trim first seconds", "Trim by event"])
        self.trim_beginning_select.currentIndexChanged.connect(self.begin_trim_changed)
        self.trim_beginning_sec = QLineEdit("0")
        self.trim_beginning_sec.setValidator(QtGui.QIntValidator())
        self.trim_beginning_event = QComboBox()
        self.trim_beginning_event.addItems(self.events_from_current_subject)
        self.trim_beginning_event.setEnabled(False)
        self.trim_ending_select = QComboBox()
        self.trim_ending_select.addItems(["Trim last seconds", "Trim by event"])
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
        self.apply_btn = QPushButton("Plot")
        self.options_main_layout.addWidget(self.apply_btn)
        # add perievent button
        self.perievent_analysis_btn = QPushButton("Perievent analysis")
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
            + "padding: 7px;\n" \
            + "margin-left:10px;\n" \
            + "margin-right:10px;\n" \
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
        self.canvas = MplCanvas(self, width=12, height=8, dpi=100)
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
        self.plot_bare_raw_data(self.preview_init_params[0][0]["subject_names"][0],
                                self.raw_data_dict[self.preview_init_params[0][0]["subject_names"][0]])
##############################################
# Previous/Next buttons on the bottom        
        # add next/previous buttons at the bottom
        self.preview_buttons_layout = QHBoxLayout()
        self.preview_buttons_layout.setAlignment(Qt.AlignRight)
        self.preview_buttons_widget = QWidget()
        self.include_cb = QCheckBox("Include Subject")
        self.next_btn = QPushButton("Next")
        self.previous_btn = QPushButton("Previous")
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
        self.perievent_analysis_btn.setEnabled(True)
        self.export_data_btn.setEnabled(True)
        self.next_btn.setEnabled(True)
        self.previous_btn.setEnabled(True)
        self.parent_window.enable_all_buttons()
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

        # if subject was not previewed yet, read data and add to dict
        if self.options["subject"] not in self.raw_data_dict:
            self.get_raw_data(self.options["subject"],
                          self.parent_window.preview_params[1][self.options["subject"]])
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
            begin_evt_data_onsets = fpExplorer_functions.get_event_on_off(self.raw_data_dict[self.options["subject"]], trim_begin_event)[0]
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
            end_evt_data_onsets = fpExplorer_functions.get_event_on_off(self.raw_data_dict[self.options["subject"]], trim_end_event)[0]
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
            self.trimmed_raw_data_dict[self.options["subject"]] = [(trim_beginning,trim_end),fpExplorer_functions.trim_raw_data(self.raw_data_dict[self.options["subject"]],
                                                                                                                      self.preview_init_params[0][0]["signal_name"],
                                                                                                                      self.preview_init_params[0][0]["control_name"],
                                                                                                                      trim_beginning,
                                                                                                                      trim_end
                                                                                                                      )]
        else:   # if already in trimmed check previous trimming settings
            trimmed = self.trimmed_raw_data_dict[self.options["subject"]]
            begin,end = trimmed[0]
            if trim_beginning != begin or trim_end != end:
                # if new trimming params, replace the lists in dict
                self.trimmed_raw_data_dict[self.options["subject"]] = [(trim_beginning,trim_end),fpExplorer_functions.trim_raw_data(self.raw_data_dict[self.options["subject"]],
                                                                                                                      self.preview_init_params[0][0]["signal_name"],
                                                                                                                      self.preview_init_params[0][0]["control_name"],
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
                evt = fpExplorer_functions.get_event_on_off(self.raw_data_dict[self.options["subject"]], self.options["event"])
                if len(evt[0]) == 0:
                    self.show_info_dialog("Some "+self.options["event"]+" event data is missing.")
                self.event_data.append(evt)
            if self.options["event2"] != "---":
                evt = fpExplorer_functions.get_event_on_off(self.raw_data_dict[self.options["subject"]], self.options["event2"])
                if len(evt[0])==0:
                    self.show_info_dialog("Some "+self.options["event2"]+" event data is missing.")
                self.event_data.append(evt)
            # if downsample was selected
            if self.downsample_cb.isChecked() or self.downsampled_export == True or self.save_plots == True or self.separate_signal_contol_cb.isChecked():
                # add to downdampled dict
                self.downsampled_dict[self.options["subject"]] = fpExplorer_functions.downsample_tdt(self.trimmed_raw_data_dict[self.options["subject"]][1],
                                            self.settings_dict[0]["downsample"])
            # if normalize was selected
            if self.normalize_cb.isChecked() or self.normalized_export == True:
                # always downsample first before normalizing
                # downsample
                self.downsampled_dict[self.options["subject"]] = fpExplorer_functions.downsample_tdt(self.trimmed_raw_data_dict[self.options["subject"]][1],
                                                                                        self.settings_dict[0]["downsample"])
                # check settings for method to normalize
                if self.settings_dict[0]["normalization"] == "Modified Polynomial Fitting":                   
                    # normalize from downsampled
                    self.normalized_dict[self.options["subject"]] = fpExplorer_functions.normalize_dff(self.raw_data_dict[self.options["subject"]],
                                                                                        self.downsampled_dict[self.options["subject"]],
                                                                                        self.settings_dict[0]["show_norm_as"],
                                                                                        self.settings_dict[0]["filter"],
                                                                                        self.settings_dict[0]["filter_fraction"])
                if self.settings_dict[0]["normalization"] == "Standard Polynomial Fitting":                
                    # normalize from downsampled
                    self.normalized_dict[self.options["subject"]] = fpExplorer_functions.normalize_pMat(self.raw_data_dict[self.options["subject"]],
                                                                                        self.downsampled_dict[self.options["subject"]],
                                                                                        self.settings_dict[0]["show_norm_as"],
                                                                                        self.settings_dict[0]["filter"],
                                                                                        self.settings_dict[0]["filter_fraction"])
            # check what to show on the plot
            self.options["plot_separate"] = True if self.separate_signal_contol_cb.isChecked() else False 
            self.options["plot_downsampled"] = True if self.downsampled_plot_cb.isChecked() else False
            self.options["plot_normalized"] = True if self.normalized_plot_cb.isChecked() else False
            if self.options["plot_separate"] == False and self.options["plot_downsampled"] == False and self.options["plot_normalized"] == False:
                self.raw_plot_cb.setChecked(True)
            self.options["plot_raw"] = True if self.raw_plot_cb.isChecked() else False
            if self.raw_plot_cb.isChecked(): # if it was set as checked the first time, uncheck it in the background(it is not visible anymore)
                self.raw_plot_cb.setChecked(False)
            # self.options["plot_separate"] = True if self.separate_signal_contol_cb.isChecked() else False 
            # self.options["plot_downsampled"] = True if self.downsampled_plot_cb.isChecked() else False
            # self.options["plot_normalized"] = True if self.normalized_plot_cb.isChecked() else False
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
                    fpExplorer_functions.plot_trimmed(self.canvas,
                                        self.options["subject"],
                                        self.trimmed_raw_data_dict[self.options["subject"]][1],
                                        self.raw_export,
                                        self.save_plots,
                                        (self.export_path,self.export_begining),
                                        self.preview_init_params[0][0]["signal_name"],
                                        self.preview_init_params[0][0]["control_name"])
                else:
                    custom_event_name = self.options["event"] if len(self.options["event_name"])==0 else self.options["event_name"]
                    custom_event_name2 = self.options["event2"] if len(self.options["event2_name"])==0 else self.options["event2_name"]
                    fpExplorer_functions.plot_with_event(self.canvas,
                                        self.options["subject"],
                                        self.trimmed_raw_data_dict[self.options["subject"]][1],
                                        custom_event_name,
                                        custom_event_name2,
                                        self.event_data,
                                        self.raw_export,
                                        self.save_plots,
                                        (self.export_path,self.export_begining),
                                        self.preview_init_params[0][0]["signal_name"],
                                        self.preview_init_params[0][0]["control_name"])
            # if only downsampled was checked
            elif (self.options["plot_raw"]==False and self.options["plot_separate"]==False and self.options["subject"] in self.downsampled_dict and self.options["plot_downsampled"]==True 
                and self.options["plot_normalized"]==False):
                if self.options["event"] == "---": # no event selected
                    fpExplorer_functions.plot_downsampled_alone(self.canvas,
                                            self.options,
                                            self.downsampled_dict[self.options["subject"]],
                                            self.downsampled_export,
                                            self.save_plots,
                                            (self.export_path,self.export_begining),
                                            self.settings_dict,
                                            self.preview_init_params[0][0]["signal_name"],
                                            self.preview_init_params[0][0]["control_name"])
                else:   # with event
                    custom_event_name = self.options["event"] if len(self.options["event_name"])==0 else self.options["event_name"]
                    custom_event_name2 = self.options["event2"] if len(self.options["event2_name"])==0 else self.options["event2_name"]
                    fpExplorer_functions.plot_downsampled_alone_with_event(self.canvas,
                                            self.options,
                                            self.downsampled_dict[self.options["subject"]],
                                            custom_event_name,
                                            custom_event_name2,
                                            self.event_data,
                                            self.downsampled_export,
                                            self.save_plots,
                                            (self.export_path,self.export_begining),
                                            self.settings_dict,
                                            self.preview_init_params[0][0]["signal_name"],
                                            self.preview_init_params[0][0]["control_name"])
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
        # check events for new data
        # first check if there are the same events
        previous_events = self.events_from_current_subject
        self.events_from_current_subject = fpExplorer_functions.get_events(self.raw_data_dict[self.options["subject"]])  
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
#        self.check_trimming()
#        self.read_options()
        self.check_events()
        self.check_trimming()
        # set include to actual status
        if self.subject_comboBox.currentText() in self.parent_window.selected_subjects:
            self.include_cb.setChecked(True)
        else:
            self.include_cb.setChecked(False)
            
    def on_subject_change(self):
        self.options_pre_check()
        # update last timestamp
        # get last timestamp in seconds
        self.last_raw_ts = fpExplorer_functions.get_last_timestamp(
                                     self.raw_data_dict[self.options["subject"]],
                                     self.preview_init_params[0][0]["signal_name"],
                                     )
        # read new subject's frequency and update suggested downsampled rate (only if it is different)
        new_fs = fpExplorer_functions.get_frequency(self.raw_data_dict[self.preview_init_params[0][0]["subject_names"][0]],self.preview_init_params[0][0]["signal_name"])
        if new_fs != self.current_fs:
            self.suggested_downsample_samples = int(int(self.current_fs)*DEFAULT_DOWNSAMPLE_PCT/100)
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
                fpExplorer_functions.show_polynomial_fitting(self.canvas,
                                                self.settings_dict[0], 
                                                self.downsampled_dict[self.options["subject"]],
                                                self.preview_init_params[0][0]["signal_name"],
                                                self.preview_init_params[0][0]["control_name"],
                                                self.options["subject"],
                                                False,
                                                "")
        else: # prompt user to check normalize first
            self.show_info_dialog("Please check Normalize check box first.")
        self.enable_buttons_signal.emit()
            
    def export_data_btn_clicked(self):
        # pass a list with first element=subject
        self.export_window = ExportDataWindow(self.parent_window,self,[self.subject_comboBox.currentText(),
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
                evt_list = fpExplorer_functions.get_events(self.raw_data_dict[subject])
                subjects_event_sets.append(set(evt_list))
            # find intersection of all subjects events (common events)  
            common_events = list(set.intersection(*subjects_event_sets))    
        self.perievent_window = PeriEventOptionsWindow(self.parent_window,[self.options["subject"],
                                                                           self.events_from_current_subject,
                                                                           self.raw_data_dict[self.options["subject"]]],
                                                                           self.perievent_options_dict,
                                                                           self.batch_perievent,
                                                                           common_events)
        self.perievent_window.got_peri_event_options_sig.connect(self.get_perievent_options)
        self.perievent_window.show()
            
    # add raw data structure and subject to self.raw_data dictionary   
    def get_raw_data(self,subject,my_path):
        self.raw_data_dict[subject] = fpExplorer_functions.get_raw_data(my_path)
    # when app starts plot just raw data  
    def plot_bare_raw_data(self,subject,raw_data):
        fpExplorer_functions.plot_raw(self.canvas,
                           subject,
                           raw_data,
                           self.preview_init_params[0][0]["signal_name"],
                           self.preview_init_params[0][0]["control_name"],
                           self.raw_export,
                           self.save_plots,
                           (self.export_path,self.export_begining)
                           )

    def read_trials(self):
        checked = []
        for btn in self.trials_button_group:
            # print(btn.text(),btn.isChecked())
            if btn.isChecked():
                checked.append(int(btn.text()))
        # print("trials",checked)
        return checked


    @pyqtSlot(list) 
    # receives a list with perievent options read from perievent window
    def get_perievent_options(self,perievent_options):
        self.perievent_options_dict = perievent_options[0]
        print("perievent options",self.perievent_options_dict)
        if self.current_perievent == self.perievent_options_dict['event'] and len(self.current_trials)>0: 
            # read the trials radio buttons
            try:
                self.current_trials = self.read_trials()
            except:
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
                                # if subject data has not been trimmed
                                # if it was not trimmed yet, add to dictionary with zero trimming
                                self.trimmed_raw_data_dict[subject] = [(0,0),fpExplorer_functions.trim_raw_data(self.raw_data_dict[subject],
                                                                                                            self.preview_init_params[0][0]["signal_name"],
                                                                                                            self.preview_init_params[0][0]["control_name"],
                                                                                                            0,
                                                                                                            0
                                                                                                            )]
                            if subject not in self.trimmed_raw_data_dict: # if raw data read but not trimmed assume trimming 0
                                # key is the subject name, value is a list
                                # first element of that list is beginning and end seconds to trim
                                # second element is trimmed data dict ts:timestamps, signal:data,control:data
                                self.trimmed_raw_data_dict[subject] = [(0,0),fpExplorer_functions.trim_raw_data(self.raw_data_dict[subject],
                                                                                                            self.preview_init_params[0][0]["signal_name"],
                                                                                                            self.preview_init_params[0][0]["control_name"],
                                                                                                            0,
                                                                                                            0
                                                                                                            )]
                            # filter around trimmed
                            data = fpExplorer_functions.filter_data_around_event(self.raw_data_dict[subject],self.trimmed_raw_data_dict[subject],
                                                    self.perievent_options_dict,
                                                    self.settings_dict,
                                                    self.preview_init_params[0][0]["signal_name"],
                                                    self.preview_init_params[0][0]["control_name"])
                            # analyze
                            analyzed_perievent_dict = fpExplorer_functions.analyze_perievent_data(data,
                                                                                self.current_trials,
                                                                               self.perievent_options_dict,
                                                                               self.settings_dict,
                                                                               self.preview_init_params[0][0]["signal_name"],
                                                                               self.preview_init_params[0][0]["control_name"])
                            # check if there is any data to plot 
                            if (len(data.streams[self.preview_init_params[0][0]["signal_name"]].filtered) > 0) and (len(data.streams[self.preview_init_params[0][0]["control_name"]].filtered) > 0):

                                # show buttons of available trials
                                self.total_current_trials = len(data.streams[self.preview_init_params[0][0]["signal_name"]].filtered)
                                if self.current_perievent == None:
                                    self.current_perievent = self.perievent_options_dict['event']
                                    # add trials to view
                                    self.current_trials_layout = QHBoxLayout()
                                    self.current_trials_widget.setLayout(self.current_trials_layout)
                                    text_label = QLabel("Include Trials:")
                                    self.current_trials_layout.addWidget(text_label)
                                    for i in range(len(data.streams[self.preview_init_params[0][0]["signal_name"]].filtered)):
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
                                all_trials_df = fpExplorer_functions.plot_raw_perievents(self.canvas,
                                                      subject,
                                                      data,
                                                      self.current_trials,
                                                      self.perievent_options_dict,
                                                      self.settings_dict,
                                                      self.preview_init_params[0][0]["signal_name"],
                                                      self.preview_init_params[0][0]["control_name"],
                                                      self.parent_window.batch_export_settings_dict["export_for_single_subjects"],
                                                      self.save_plots,
                                                      self.group_names_dict[subject],
                                                      (subject_subfolder,self.parent_window.batch_export_settings_dict["file_begin"]))  
                                if self.parent_window.batch_export_settings_dict["export_group_data"] == True:
                                    all_subjects_peri_normalized_dfs.append((subject,all_trials_df))
                                if self.perievent_options_dict["plot_avg"] == True:
                                    fpExplorer_functions.plot_perievent_average_alone(self.canvas,
                                                                  subject,
                                                                  self.perievent_options_dict,
                                                                  analyzed_perievent_dict,
                                                                  self.parent_window.batch_export_settings_dict["export_for_single_subjects"],
                                                                  self.save_plots,
                                                                  self.group_names_dict[subject],
                                                                  self.settings_dict,
                                                                  self.preview_init_params[0][0]["signal_name"],
                                                                  self.preview_init_params[0][0]["control_name"],
                                                                  (subject_subfolder,self.parent_window.batch_export_settings_dict["file_begin"]))
                                if (self.perievent_options_dict["plot_zscore"] == True or self.perievent_options_dict["plot_zscore_trials"] == True 
                                    or self.perievent_options_dict["plot_auc"] == True): # we need zcsore data for auc
                                    if self.perievent_options_dict["plot_zscore"] == True:
                                        z_score_df = fpExplorer_functions.plot_perievent_zscore_alone(self.canvas,
                                                                      subject,
                                                                      data,
                                                                      self.perievent_options_dict,
                                                                      analyzed_perievent_dict,
                                                                      self.preview_init_params[0][0]["signal_name"],
                                                                      self.parent_window.batch_export_settings_dict["export_for_single_subjects"],
                                                                      self.save_plots,
                                                                      self.group_names_dict[subject],
                                                                      self.settings_dict,
                                                                      (subject_subfolder,self.parent_window.batch_export_settings_dict["file_begin"]))
                                    else:
#                                    if self.perievent_options_dict["plot_zscore_trials"] == True:
                                        z_score_df = fpExplorer_functions.plot_perievent_zscore_with_trials_alone(self.canvas,
                                                                      subject,
                                                                      data,
                                                                      self.perievent_options_dict,
                                                                      analyzed_perievent_dict,
                                                                      self.preview_init_params[0][0]["signal_name"],
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
                                        # if subject data has not been trimmed
                                        # if it was not trimmed yet, add to dictionary with zero trimming
                                        self.trimmed_raw_data_dict[subject] = [(0,0),fpExplorer_functions.trim_raw_data(self.raw_data_dict[subject],
                                                                                                                    self.preview_init_params[0][0]["signal_name"],
                                                                                                                    self.preview_init_params[0][0]["control_name"],
                                                                                                                    0,
                                                                                                                    0
                                                                                                                    )]
                                    if subject not in self.trimmed_raw_data_dict: # if raw data read but not trimmed assume trimming 0
                                        # key is the subject name, value is a list
                                        # first element of that list is beginning and end seconds to trim
                                        # second element is trimmed data dict ts:timestamps, signal:data,control:data
                                        self.trimmed_raw_data_dict[subject] = [(0,0),fpExplorer_functions.trim_raw_data(self.raw_data_dict[subject],
                                                                                                                    self.preview_init_params[0][0]["signal_name"],
                                                                                                                    self.preview_init_params[0][0]["control_name"],
                                                                                                                    0,
                                                                                                                    0
                                                                                                                    )]
                                    # filter around trimmed
                                    data = fpExplorer_functions.filter_data_around_event(self.raw_data_dict[subject],self.trimmed_raw_data_dict[subject],
                                                                self.perievent_options_dict,
                                                                self.settings_dict,
                                                                self.preview_init_params[0][0]["signal_name"],
                                                                self.preview_init_params[0][0]["control_name"])
                                    # check if there is any data to plot 
                                    if (len(data.streams[self.preview_init_params[0][0]["signal_name"]].filtered) > 0) and (len(data.streams[self.preview_init_params[0][0]["control_name"]].filtered) > 0):
                                        if (self.perievent_options_dict["plot_zscore"] == True or self.perievent_options_dict["plot_zscore_trials"] == True 
                                            or self.perievent_options_dict["plot_auc"] == True):
                                            # analyze
                                            analyzed_perievent_dict = fpExplorer_functions.analyze_perievent_data(data,
                                                                                            self.current_trials,
                                                                                           self.perievent_options_dict,
                                                                                           self.settings_dict,
                                                                                           self.preview_init_params[0][0]["signal_name"],
                                                                                           self.preview_init_params[0][0]["control_name"])

                                            # show buttons of available trials
                                            self.total_current_trials = len(data.streams[self.preview_init_params[0][0]["signal_name"]].filtered)
                                            if self.current_perievent == None:
                                                self.current_perievent = self.perievent_options_dict['event']
                                                # add trials to view
                                                self.current_trials_layout = QHBoxLayout()
                                                self.current_trials_widget.setLayout(self.current_trials_layout)
                                                text_label = QLabel("Include Trials:")
                                                self.current_trials_layout.addWidget(text_label)
                                                for i in range(len(data.streams[self.preview_init_params[0][0]["signal_name"]].filtered)):
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
                                                
                                        if self.perievent_options_dict["plot_avg"] == True:
                                            # perieventnormalized preview
                                            all_trials_df = fpExplorer_functions.plot_raw_perievents(self.canvas,
                                                                      subject,
                                                                      data,
                                                                      self.current_trials,
                                                                      self.perievent_options_dict,
                                                                      self.settings_dict,
                                                                      self.preview_init_params[0][0]["signal_name"],
                                                                      self.preview_init_params[0][0]["control_name"],
                                                                      self.parent_window.batch_export_settings_dict["export_for_single_subjects"],
                                                                      self.save_plots,
                                                                      self.group_names_dict[subject],
                                                                      (self.parent_window.batch_export_settings_dict["dump_path"],self.parent_window.batch_export_settings_dict["file_begin"]))  
                                            all_subjects_peri_normalized_dfs.append((subject,all_trials_df))
                                        if (self.perievent_options_dict["plot_zscore"] == True or self.perievent_options_dict["plot_zscore_trials"] == True 
                                            or self.perievent_options_dict["plot_auc"] == True):
                                            if self.perievent_options_dict["plot_zscore"] == True:
                                                z_score_df = fpExplorer_functions.plot_perievent_zscore_alone(self.canvas,
                                                                          subject,
                                                                          data,
                                                                          self.perievent_options_dict,
                                                                          analyzed_perievent_dict,
                                                                          self.preview_init_params[0][0]["signal_name"],
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
                                                                          data,
                                                                          self.perievent_options_dict,
                                                                          analyzed_perievent_dict,
                                                                          self.preview_init_params[0][0]["signal_name"],
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
                self.trimmed_raw_data_dict[self.options["subject"]] = [(0,0),fpExplorer_functions.trim_raw_data(self.raw_data_dict[self.options["subject"]],
                                                                                                        self.preview_init_params[0][0]["signal_name"],
                                                                                                        self.preview_init_params[0][0]["control_name"],
                                                                                                        0,
                                                                                                        0
                                                                                                        )]
            # filter around trimmed
            data = fpExplorer_functions.filter_data_around_event(self.raw_data_dict[self.options["subject"]],self.trimmed_raw_data_dict[self.options["subject"]],
                                                    self.perievent_options_dict,
                                                    self.settings_dict,
                                                    self.preview_init_params[0][0]["signal_name"],
                                                    self.preview_init_params[0][0]["control_name"])

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
            if (len(data.streams[self.preview_init_params[0][0]["signal_name"]].filtered) > 0) and (len(data.streams[self.preview_init_params[0][0]["control_name"]].filtered) > 0):

                # show buttons of available trials
                print("How many trials are there?",len(data.streams[self.preview_init_params[0][0]["signal_name"]].filtered))
                self.total_current_trials = len(data.streams[self.preview_init_params[0][0]["signal_name"]].filtered)
                if self.current_perievent == None:
                    self.current_perievent = self.perievent_options_dict['event']
                    # add trials to view
                    self.current_trials_layout = QHBoxLayout()
                    self.current_trials_widget.setLayout(self.current_trials_layout)
                    text_label = QLabel("Include Trials:")
                    self.current_trials_layout.addWidget(text_label)
                    for i in range(len(data.streams[self.preview_init_params[0][0]["signal_name"]].filtered)):
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
                elif self.current_perievent != self.perievent_options_dict['event']:  # if new event different from previous       
                    # clear current trials and buttons group
                    self.current_trials = []
                    self.trials_button_group = []
                    # create new buttons
                    new_trials_widget = QWidget()
                    new_trials_layout = QHBoxLayout()
                    new_trials_widget.setLayout(new_trials_layout)
                    text_label = QLabel("Include Trials:")
                    new_trials_layout.addWidget(text_label)
                    for i in range(len(data.streams[self.preview_init_params[0][0]["signal_name"]].filtered)):
                        # create btn
                        # add button to a group
                        # add to list
                        # add to layout
                        btn = QCheckBox(str(i+1))
                        btn.setChecked(True)
                        new_trials_layout.addWidget(btn)
                        self.trials_button_group.append(btn)
                        self.current_trials.append(i+1)
                    # replace with updated widget
                    self.trials_layout.replaceWidget(self.current_trials_widget,new_trials_widget)
                    self.current_trials_widget.deleteLater()
                    self.current_trials_widget = new_trials_widget

                    # set new event
                    self.current_perievent = self.perievent_options_dict['event']
                else: # if same event and/or new subject
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
                        for i in range(len(data.streams[self.preview_init_params[0][0]["signal_name"]].filtered)):
                            # create btn
                            # add button to a group
                            # add to list
                            # add to layout
                            btn = QCheckBox(str(i+1))
                            btn.setChecked(True)
                            new_trials_layout.addWidget(btn)
                            self.trials_button_group.append(btn)
                            self.current_trials.append(i+1)
                    # replace with updated widget
                    self.trials_layout.replaceWidget(self.current_trials_widget,new_trials_widget)
                    self.current_trials_widget.deleteLater()
                    self.current_trials_widget = new_trials_widget
            ##########################################

                # if user clicked on preview in perievent options window, show preview
                if self.perievent_options_dict["preview"] == True:
                    # plot normalized preview
                    fpExplorer_functions.plot_raw_perievents(self.canvas,
                                                  self.options["subject"],
                                                  data,
                                                  self.current_trials,
                                                  self.perievent_options_dict,
                                                  self.settings_dict,
                                                  self.preview_init_params[0][0]["signal_name"],
                                                  self.preview_init_params[0][0]["control_name"],
                                                  self.perievent_options_dict["export"],
                                                  self.save_plots,
                                                  self.group_names_dict[self.options["subject"]],
                                                  (self.export_path,self.export_begining))  
                # if user clicked on analyze in perievent options window, analyze
                if self.perievent_options_dict["analyze"] == True:
                    analyzed_perievent_dict = fpExplorer_functions.analyze_perievent_data(data,
                                                                                self.current_trials,
                                                                               self.perievent_options_dict,
                                                                               self.settings_dict,self.preview_init_params[0][0]["signal_name"],
                                                                               self.preview_init_params[0][0]["control_name"])
                    # find out what to show on plot
                    # only average
                    if (self.perievent_options_dict["plot_avg"] == True and self.perievent_options_dict["plot_zscore"] == False
                        and self.perievent_options_dict["plot_zscore_trials"] == False and self.perievent_options_dict["plot_auc"] == False):
                        fpExplorer_functions.plot_perievent_average_alone(self.canvas,
                                                      self.options["subject"],
                                                      self.perievent_options_dict,
                                                      analyzed_perievent_dict,
                                                      self.perievent_options_dict["export"],
                                                      self.save_plots,
                                                      self.group_names_dict[self.options["subject"]],
                                                      self.settings_dict,
                                                      self.preview_init_params[0][0]["signal_name"],
                                                      self.preview_init_params[0][0]["control_name"],
                                                      (self.export_path,self.export_begining))
                    # only zscore with error
                    elif (self.perievent_options_dict["plot_avg"] == False and self.perievent_options_dict["plot_zscore"] == True
                            and self.perievent_options_dict["plot_zscore_trials"] == False and self.perievent_options_dict["plot_auc"] == False):
                        fpExplorer_functions.plot_perievent_zscore_alone(self.canvas,
                                                      self.options["subject"],
                                                      data,
                                                      self.perievent_options_dict,
                                                      analyzed_perievent_dict,
                                                      self.preview_init_params[0][0]["signal_name"],
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
                                                      data,
                                                      self.perievent_options_dict,
                                                      analyzed_perievent_dict,
                                                      self.preview_init_params[0][0]["signal_name"],
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
                                                      data,
                                                      self.perievent_options_dict,
                                                      analyzed_perievent_dict,
                                                      self.preview_init_params[0][0]["signal_name"])
                    # avg and zscore with trials
                    elif (self.perievent_options_dict["plot_avg"] == True and self.perievent_options_dict["plot_zscore"] == False
                            and self.perievent_options_dict["plot_zscore_trials"] == True and self.perievent_options_dict["plot_auc"] == False):
                        fpExplorer_functions.plot_perievent_avg_zscore_trials(self.canvas,
                                                      self.options["subject"],
                                                      data,
                                                      self.perievent_options_dict,
                                                      analyzed_perievent_dict,
                                                      self.preview_init_params[0][0]["signal_name"])
                    # avg and auc
                    elif (self.perievent_options_dict["plot_avg"] == True and self.perievent_options_dict["plot_zscore"] == False
                            and self.perievent_options_dict["plot_zscore_trials"] == False and self.perievent_options_dict["plot_auc"] == True):
                        fpExplorer_functions.plot_perievent_avg_auc(self.canvas,
                                                      self.options["subject"],
                                                      self.perievent_options_dict,
                                                      analyzed_perievent_dict)
                    #zscore with error and auc
                    elif (self.perievent_options_dict["plot_avg"] == False and self.perievent_options_dict["plot_zscore"] == True
                            and self.perievent_options_dict["plot_zscore_trials"] == False and self.perievent_options_dict["plot_auc"] == True):
                        fpExplorer_functions.plot_perievent_zscore_auc(self.canvas,
                                                      self.options["subject"],
                                                      data,
                                                      self.perievent_options_dict,
                                                      analyzed_perievent_dict,
                                                      self.preview_init_params[0][0]["signal_name"])
                    #zscore with trials and auc
                    elif (self.perievent_options_dict["plot_avg"] == False and self.perievent_options_dict["plot_zscore"] == False
                            and self.perievent_options_dict["plot_zscore_trials"] == True and self.perievent_options_dict["plot_auc"] == True):
                        fpExplorer_functions.plot_perievent_zscore_trials_auc(self.canvas,
                                                      self.options["subject"],
                                                      data,
                                                      self.perievent_options_dict,
                                                      analyzed_perievent_dict,
                                                      self.preview_init_params[0][0]["signal_name"])
                    # all 4 plots (zscore with error)
                    elif (self.perievent_options_dict["plot_avg"] == True and self.perievent_options_dict["plot_zscore"] == True
                            and self.perievent_options_dict["plot_zscore_trials"] == False and self.perievent_options_dict["plot_auc"] == True):
                        fpExplorer_functions.plot_all_perievent(self.canvas,
                                                      self.options["subject"],
                                                      data,
                                                      self.perievent_options_dict,
                                                      analyzed_perievent_dict,
                                                      self.preview_init_params[0][0]["signal_name"])
                    # all 4 plots (zscore with trials)
                    elif (self.perievent_options_dict["plot_avg"] == True and self.perievent_options_dict["plot_zscore"] == False
                            and self.perievent_options_dict["plot_zscore_trials"] == True and self.perievent_options_dict["plot_auc"] == True):
                        fpExplorer_functions.plot_all_perievent_zscore_trials(self.canvas,
                                                      self.options["subject"],
                                                      data,
                                                      self.perievent_options_dict,
                                                      analyzed_perievent_dict,
                                                      self.preview_init_params[0][0]["signal_name"])
                # if user clicked on analyze in perievent options window, analyze
                if self.perievent_options_dict["export"] == True:
                    analyzed_perievent_dict = fpExplorer_functions.analyze_perievent_data(data,
                                                                                self.current_trials,
                                                                               self.perievent_options_dict,
                                                                               self.settings_dict,self.preview_init_params[0][0]["signal_name"],
                                                                               self.preview_init_params[0][0]["control_name"])
                    # save downsampled preview
                    fpExplorer_functions.plot_raw_perievents(self.canvas,
                                                  self.options["subject"],
                                                  data,
                                                  self.current_trials,
                                                  self.perievent_options_dict,
                                                  self.settings_dict,
                                                  self.preview_init_params[0][0]["signal_name"],
                                                  self.preview_init_params[0][0]["control_name"],
                                                  self.perievent_options_dict["export"],
                                                  self.save_plots,
                                                  self.group_names_dict[self.options["subject"]],
                                                  (self.export_path,self.export_begining))  
                    if self.perievent_options_dict["plot_avg"] == True:
                        fpExplorer_functions.plot_perievent_average_alone(self.canvas,
                                                      self.options["subject"],
                                                      self.perievent_options_dict,
                                                      analyzed_perievent_dict,
                                                      self.perievent_options_dict["export"],
                                                      self.save_plots,
                                                      self.group_names_dict[self.options["subject"]],
                                                      self.settings_dict,
                                                      self.preview_init_params[0][0]["signal_name"],
                                                      self.preview_init_params[0][0]["control_name"],
                                                      (self.export_path,self.export_begining))
                    if self.perievent_options_dict["plot_zscore"] == True:
                        fpExplorer_functions.plot_perievent_zscore_alone(self.canvas,
                                                      self.options["subject"],
                                                      data,
                                                      self.perievent_options_dict,
                                                      analyzed_perievent_dict,
                                                      self.preview_init_params[0][0]["signal_name"],
                                                      self.perievent_options_dict["export"],
                                                      self.save_plots,
                                                      self.group_names_dict[self.options["subject"]],
                                                      self.settings_dict,
                                                      (self.export_path,self.export_begining))
                    if self.perievent_options_dict["plot_zscore_trials"] == True:
                        fpExplorer_functions.plot_perievent_zscore_with_trials_alone(self.canvas,
                                                      self.options["subject"],
                                                      data,
                                                      self.perievent_options_dict,
                                                      analyzed_perievent_dict,
                                                      self.preview_init_params[0][0]["signal_name"],
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
                self.export_window = ExportDataWindow(self.parent_window,self,[self.subject_comboBox.currentText(),
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
                fpExplorer_functions.show_polynomial_fitting(self.canvas,
                                                self.settings_dict[0], 
                                                self.downsampled_dict[self.options["subject"]],
                                                self.preview_init_params[0][0]["signal_name"],
                                                self.preview_init_params[0][0]["control_name"],
                                                self.options["subject"],
                                                self.save_plots,
                                                self.export_path)
            # raw was checked
            if self.raw_export == True:
                if self.options["event"] == "---":
                    self.plot_bare_raw_data(self.options["subject"],self.raw_data_dict[self.options["subject"]])
                else:
                    custom_event_name = self.options["event"] if len(self.options["event_name"])==0 else self.options["event_name"]
                    custom_event_name2 = self.options["event2"] if len(self.options["event2_name"])==0 else self.options["event2_name"]
                    fpExplorer_functions.plot_with_event(self.canvas,
                                        self.options["subject"],
                                        self.raw_data_dict[self.options["subject"]],
                                        custom_event_name,
                                        custom_event_name2,
                                        self.event_data,
                                        self.raw_export,
                                        self.save_plots,
                                        (self.export_path,self.export_begining),
                                        self.preview_init_params[0][0]["signal_name"],
                                        self.preview_init_params[0][0]["control_name"])
            if self.downsampled_export == True:
                if self.options["event"] == "---": # no event selected
                    fpExplorer_functions.plot_downsampled_alone(self.canvas,
                                            self.options,
                                            self.downsampled_dict[self.options["subject"]],
                                            self.downsampled_export,
                                            self.save_plots,
                                            (self.export_path,self.export_begining),
                                            self.settings_dict,
                                            self.preview_init_params[0][0]["signal_name"],
                                            self.preview_init_params[0][0]["control_name"])
                else:   # with event
                    custom_event_name = self.options["event"] if len(self.options["event_name"])==0 else self.options["event_name"]
                    custom_event_name2 = self.options["event2"] if len(self.options["event2_name"])==0 else self.options["event2_name"]
                    fpExplorer_functions.plot_downsampled_alone_with_event(self.canvas,
                                            self.options,
                                            self.downsampled_dict[self.options["subject"]],
                                            custom_event_name,
                                            custom_event_name2,
                                            self.event_data,
                                            self.downsampled_export,
                                            self.save_plots,
                                            (self.export_path,self.export_begining),
                                            self.settings_dict,
                                            self.preview_init_params[0][0]["signal_name"],
                                            self.preview_init_params[0][0]["control_name"])
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
        self.peak_options_window = PeakSettingsWindow(self.parent_window,
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
                    # if self.parent_window.batch_export_settings_dict["export_group_data"] == True:
                    #     pass # single subjects only
                        # # spikes for averaged signal
                        # if self.options["event"] == "---":
                        #     try:
                        #         fpExplorer_functions.get_batch_spikes(self.canvas,
                        #                                        self.recent_peak_values,
                        #                                        self.all_normalized,
                        #                                        self.settings_dict,
                        #                                        self.parent_window.batch_export_settings_dict["spikes"],
                        #                                        (self.parent_window.batch_export_settings_dict["dump_path"],self.parent_window.batch_export_settings_dict["file_begin"]))
                        #     except:
                        #         self.show_info_dialog("Could not calculate spikes.\nTry with different parameters.")
                        # else: # events
                        #     custom_event_name = self.options["event"] if len(self.options["event_name"])==0 else self.options["event_name"]
                        #     custom_event_name2 = self.options["event2"] if len(self.options["event2_name"])==0 else self.options["event2_name"]
                            
                        #     fpExplorer_functions.get_batch_spikes_with_event(self.canvas,
                        #                                           self.recent_peak_values,
                        #                                           self.all_normalized,
                        #                                           custom_event_name,
                        #                                           custom_event_name2,
                        #                                           self.event_data,
                        #                                           self.settings_dict,
                        #                                           self.parent_window.batch_export_settings_dict["spikes"],
                        #                                           (self.parent_window.batch_export_settings_dict["dump_path"],self.parent_window.batch_export_settings_dict["file_begin"]))
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
            if len(export_path) > 0 and os.path.split(export_path)[1]==DEFAULT_EXPORT_FOLDER:
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
                self.trimmed_raw_data_dict[self.options["subject"]] = [(trim_beginning,trim_end),fpExplorer_functions.trim_raw_data(self.raw_data_dict[self.options["subject"]],
                                                                                                                          self.preview_init_params[0][0]["signal_name"],
                                                                                                                          self.preview_init_params[0][0]["control_name"],
                                                                                                                          trim_beginning,
                                                                                                                          trim_end
                                                                                                                          )]
            else:   # if already in trimmed check previous trimming settings
                trimmed = self.trimmed_raw_data_dict[self.options["subject"]]
                begin,end = trimmed[0]
                if trim_beginning != begin or trim_end != end:
                    # if new trimming params, replace the lists in dict
                    self.trimmed_raw_data_dict[self.options["subject"]] = [(trim_beginning,trim_end),fpExplorer_functions.trim_raw_data(self.raw_data_dict[self.options["subject"]],
                                                                                                                          self.preview_init_params[0][0]["signal_name"],
                                                                                                                          self.preview_init_params[0][0]["control_name"],
                                                                                                                          trim_beginning,
                                                                                                                          trim_end
                                                                                                                          )]
            # always downsample first before normalizing
            # downsample
            self.downsampled_dict[self.options["subject"]] = fpExplorer_functions.downsample_tdt(self.trimmed_raw_data_dict[self.options["subject"]][1],
                                                                                         self.settings_dict[0]["downsample"])
            # check settings for method to normalize
            if self.settings_dict[0]["normalization"] == "Modified Polynomial Fitting":                   
                # normalize from downsampled
                self.normalized_dict[self.options["subject"]] = fpExplorer_functions.normalize_dff(self.raw_data_dict[self.options["subject"]],
                                                                                        self.downsampled_dict[self.options["subject"]],
                                                                                        self.settings_dict[0]["show_norm_as"],
                                                                                        self.settings_dict[0]["filter"],
                                                                                        self.settings_dict[0]["filter_fraction"])
            if self.settings_dict[0]["normalization"] == "Standard Polynomial Fitting":                
                # normalize from downsampled
                self.normalized_dict[self.options["subject"]] = fpExplorer_functions.normalize_pMat(self.raw_data_dict[self.options["subject"]],
                                                                                        self.downsampled_dict[self.options["subject"]],
                                                                                        self.settings_dict[0]["show_norm_as"],
                                                                                        self.settings_dict[0]["filter"],
                                                                                        self.settings_dict[0]["filter_fraction"])
            # set event
            self.options["event"] = self.event_from_data_comboBox.currentText()
            self.options["event_name"] = self.event_name_text.text()
            self.options["event2"] = self.event2_from_data_comboBox.currentText()
            self.options["event2_name"] = self.event2_name_text.text()
            # reset the event list to make sure it always has current events
            self.event_data = []
            if self.options["event"] != "---":
                evt = fpExplorer_functions.get_event_on_off(self.raw_data_dict[self.options["subject"]], self.options["event"])
                if len(evt[0]) == 0:
                    self.show_info_dialog("Some "+self.options["event"]+" event data is missing.")
                self.event_data.append(evt)
            if self.options["event2"] != "---":
                evt = fpExplorer_functions.get_event_on_off(self.raw_data_dict[self.options["subject"]], self.options["event2"])
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
                
    @pyqtSlot()     
    def start_batch_analysis(self):
        # update group names
        self.group_names_dict = self.parent_window.batch_export_settings_dict["updated_group_names"]
        # update group name field of the current subject
        self.subject_group_name.setText(self.group_names_dict[self.subject_comboBox.currentText()])
        # update current trimming
        self.trim_beginning_sec.setText(self.parent_window.batch_export_settings_dict["trim_begin"])
        self.trim_ending_sec.setText(self.parent_window.batch_export_settings_dict["trim_end"])
        # set include to actual status
        if self.subject_comboBox.currentText() in self.parent_window.selected_subjects:
            self.include_cb.setChecked(True)
        else:
            self.include_cb.setChecked(False)
        # disable buttons while processing
        self.disable_buttons_signal.emit()
        print("Started batch processing")
        new_trim_start = int(self.parent_window.batch_export_settings_dict["trim_begin"])
        new_trim_end = int(self.parent_window.batch_export_settings_dict["trim_end"])
        
        # set event
        self.options["event"] = self.event_from_data_comboBox.currentText()
        self.options["event_name"] = self.event_name_text.text()
        self.options["event2"] = self.event2_from_data_comboBox.currentText()
        self.options["event2_name"] = self.event2_name_text.text()
        # reset the event list to make sure it always has current events
        self.event_data = []
        if self.options["event"] != "---":
            evt = fpExplorer_functions.get_event_on_off(self.raw_data_dict[self.options["subject"]], self.options["event"])
            if len(evt[0]) == 0:
                self.show_info_dialog("Some "+self.options["event"]+" event data is missing.")
            self.event_data.append(evt)
        if self.options["event2"] != "---":
            evt = fpExplorer_functions.get_event_on_off(self.raw_data_dict[self.options["subject"]], self.options["event2"])
            if len(evt[0])==0:
                self.show_info_dialog("Some "+self.options["event2"]+" event data is missing.")
            self.event_data.append(evt)
                   
        wrong_event_order_info = "If you want to show just one event on the plots,\nselect your event as first event.\nAnd leave Event2 empty."
        # before you proceed with ploting check correct order of events  if only one event is selected
        if self.options["event"] == "---" and self.options["event2"] != "---":
            self.show_info_dialog(wrong_event_order_info)
            
        # create a popup window that will be on during processing?
        # disable all buttons?
        if self.parent_window.batch_export_settings_dict["normalized"] == True or self.parent_window.batch_export_settings_dict["spikes"] == True:
            # create normalized data for most recent settings
            for i in range(len(self.parent_window.batch_export_settings_dict["batch_subjects"])):
                subject = self.parent_window.batch_export_settings_dict["batch_subjects"][i]
                # create separate options dictionary for batch analysis
                self.batch_options_dict = {"subject":subject,
                                            "subject_group_name":self.parent_window.batch_export_settings_dict["batch_subjects_group_names"][i]}
                # add the group name to group names dictionary if it was not there or update group names from batch options
                self.group_names_dict[subject] = self.parent_window.batch_export_settings_dict["batch_subjects_group_names"][i]
                # if subject was not previewed yet, read data and add to dict
                if subject not in self.raw_data_dict:
                    self.get_raw_data(subject, self.parent_window.preview_params[1][subject])
                    
                # trim to batch export settings
                # key is the subject name, value is a list
                # first element of that list is beginning and end seconds to trim
                # second element is trimmed data dict ts:timestamps, signal:data,control:data
                self.trimmed_raw_data_dict[subject] = [(new_trim_start,new_trim_end),fpExplorer_functions.trim_raw_data(self.raw_data_dict[subject],
                                                                                  self.preview_init_params[0][0]["signal_name"],
                                                                                  self.preview_init_params[0][0]["control_name"],
                                                                                  new_trim_start,
                                                                                  new_trim_end
                                                                                  )]
                
                # always downsample first before normalizing
                # downsample
                self.downsampled_dict[subject] = fpExplorer_functions.downsample_tdt(self.trimmed_raw_data_dict[subject][1],
                                                                             self.settings_dict[0]["downsample"])
                # check settings for method to normalize
                if self.settings_dict[0]["normalization"] == "Modified Polynomial Fitting":                   
                    # normalize from downsampled
                    self.normalized_dict[subject] = fpExplorer_functions.normalize_dff(self.raw_data_dict[subject],
                                                                                            self.downsampled_dict[subject],
                                                                                            self.settings_dict[0]["show_norm_as"],
                                                                                            self.settings_dict[0]["filter"],
                                                                                            self.settings_dict[0]["filter_fraction"])
                if self.settings_dict[0]["normalization"] == "Standard Polynomial Fitting":                
                    # normalize from downsampled
                    self.normalized_dict[subject] = fpExplorer_functions.normalize_pMat(self.raw_data_dict[subject],
                                                                                            self.downsampled_dict[subject],
                                                                                            self.settings_dict[0]["show_norm_as"],
                                                                                            self.settings_dict[0]["filter"],
                                                                                            self.settings_dict[0]["filter_fraction"])
                if self.parent_window.batch_export_settings_dict["export_for_single_subjects"] == True:
                    # create subfolder with subject name
                    subject_subfolder = os.path.join(self.parent_window.batch_export_settings_dict["dump_path"],subject)
                    if not os.path.exists(subject_subfolder):
                        try:
                            os.mkdir(subject_subfolder)
                        except:
                            self.show_info_dialog("Problem creating subfolder")
                            # save under main folder
                            subject_subfolder = self.parent_window.batch_export_settings_dict["dump_path"]
                    # save also polynomial fitting by default
                    fpExplorer_functions.show_polynomial_fitting(self.canvas,
                                              self.settings_dict[0], 
                                              self.downsampled_dict[subject],
                                              self.preview_init_params[0][0]["signal_name"],
                                              self.preview_init_params[0][0]["control_name"],
                                              subject,
                                              True,
                                              subject_subfolder)
                    if self.parent_window.batch_export_settings_dict["normalized"] == True:
                        if self.options["event"] == "---":
                            fpExplorer_functions.plot_normalized_alone(self.canvas,
                                                       self.batch_options_dict,
                                                       self.normalized_dict[subject],
                                                       True,
                                                       True,
                                                       (subject_subfolder,self.parent_window.batch_export_settings_dict["file_begin"]),
                                                       self.settings_dict)
                        else:
                            custom_event_name = self.options["event"] if len(self.options["event_name"])==0 else self.options["event_name"]
                            custom_event_name2 = self.options["event2"] if len(self.options["event2_name"])==0 else self.options["event2_name"]
                            fpExplorer_functions.plot_normalized_alone_with_event(self.canvas,
                                                           self.batch_options_dict,
                                                           self.normalized_dict[subject],
                                                           custom_event_name,
                                                           custom_event_name2,
                                                           self.event_data,
                                                           True,
                                                           True,
                                                           (subject_subfolder,self.parent_window.batch_export_settings_dict["file_begin"]),
                                                           self.settings_dict)
                    if self.parent_window.batch_export_settings_dict["spikes"] == True:
                        self.spikes_export = True
                        self.save_plots = True
            ### end for loop for each subject   
            if self.parent_window.batch_export_settings_dict["export_group_data"] == True:
                # average signal
                # create a list of all normalized signals
                self.all_normalized = [(subject,self.normalized_dict[subject]) for subject in self.parent_window.batch_export_settings_dict["batch_subjects"]]
    #            print(all_normalized)
                # if self.parent_window.batch_export_settings_dict["normalized"] == True:
                #     if self.options["event"] == "---":
                #         pass # only for single subjects
                        # fpExplorer_functions.get_batch_normalized(self.canvas,
                        #                                self.all_normalized,
                        #                                self.settings_dict,
                        #                                self.parent_window.batch_export_settings_dict["normalized"],
                        #                                (self.parent_window.batch_export_settings_dict["dump_path"],self.parent_window.batch_export_settings_dict["file_begin"]))
                    # else:
                    #     pass # only for single subjects
                        # custom_event_name = self.options["event"] if len(self.options["event_name"])==0 else self.options["event_name"]
                        # custom_event_name2 = self.options["event2"] if len(self.options["event2_name"])==0 else self.options["event2_name"]
                        # fpExplorer_functions.get_batch_normalized_with_event(self.canvas,
                        #                                           self.all_normalized,
                        #                                           custom_event_name,
                        #                                           custom_event_name2,
                        #                                           self.event_data,
                        #                                           self.settings_dict,
                        #                                           self.parent_window.batch_export_settings_dict["normalized"],
                        #                                           (self.parent_window.batch_export_settings_dict["dump_path"],self.parent_window.batch_export_settings_dict["file_begin"]))
                if self.parent_window.batch_export_settings_dict["spikes"] == True:
                    # set batch peaks to true
                    self.batch_peaks = True
            if self.parent_window.batch_export_settings_dict["spikes"] == True:
                # set batch peaks to true
                self.batch_peaks = True
                # remember to open it later
                self.open_spikes_later = True                  
#            else: # if user didn't want spike analysis
#                print("Done batch processing")
#                self.done_batch_processing_sig.emit()
#                self.enable_buttons_signal.emit()
        if self.parent_window.batch_export_settings_dict["perievent"] == True:
            self.batch_perievent = True
            self.save_plots = True
            self.perievent_analysis_btn_clicked()
        else:
            if self.open_spikes_later == False:
                print("Done batch processing")
                self.done_batch_processing_sig.emit()
                self.enable_buttons_signal.emit()
            elif self.open_spikes_later == True:
                self.peaks_btn_clicked()
            
            
        
        
    # show info popup
    def show_info_dialog(self, text):
        msgBox = QMessageBox()
        msgBox.setWindowIcon(QtGui.QIcon(ICO))
        msgBox.setText(text)
        msgBox.setWindowTitle("Info!")
        msgBox.setStandardButtons(QMessageBox.Ok)
        msgBox.exec()
        

########################### end PreviewEventBasedWidget class     
        
##############################################
# CLASS FOR  ADVANCED PEAK DETECTION WINDOW  #
##############################################
        
class AdvancedPeakSettingsWindow(QMainWindow):
    
    got_advanced_options_sig = pyqtSignal(list)
    got_advanced_options4export_sig = pyqtSignal(list)
    disable_btn_sig = pyqtSignal()
    def __init__(self,parent_window,recent_peak_values,export_loc,batch):
        super(AdvancedPeakSettingsWindow, self).__init__()
        self.setWindowTitle("Advanced Peak Settings")
        self.setWindowIcon(QtGui.QIcon(ICO))
        self.resize(300,300)
        self.parent_window = parent_window
        # [self.height,self.threshold,self.distance,self.prominence,self.width,self.wlen,
        # self.relheight,self.plateau]
        self.recent_vals = recent_peak_values
        self.export_path,self.export_file_begin = export_loc
        self.batch = batch
        
        
        # create widget with settings
        self.settings_main_widget = QWidget()
        self.setCentralWidget(self.settings_main_widget)
        # create main layout for widget with settings
        self.settings_main_layout = QVBoxLayout()
        self.settings_main_layout.setContentsMargins(10,10,10,10)
        self.settings_layout = QFormLayout()
        self.height_text = QLineEdit("None")
        if len(self.recent_vals) > 0:
            self.height_text.setText(str(self.recent_vals[0]))
        self.settings_layout.addRow("Height:",self.height_text)
        self.threshold_text = QLineEdit("None")
        if len(self.recent_vals) > 0:
            self.threshold_text.setText(str(self.recent_vals[1]))
        self.settings_layout.addRow("Threshold:",self.threshold_text)
        self.distance_text = QLineEdit("None")
        if len(self.recent_vals) > 0:
            self.distance_text.setText(str(self.recent_vals[2]))
        self.settings_layout.addRow("Distance (sec):",self.distance_text)
        self.prominence_text = QLineEdit("None")
        if len(self.recent_vals) > 0:
            self.prominence_text.setText(str(self.recent_vals[3]))
        self.settings_layout.addRow("Prominence:",self.prominence_text)
        self.width_text = QLineEdit("None")
        if len(self.recent_vals) > 0:
            self.width_text.setText(str(self.recent_vals[4]))
        self.settings_layout.addRow("Width:",self.width_text)
        self.wlen_text = QLineEdit("None")
        if len(self.recent_vals) > 0:
            self.wlen_text.setText(str(self.recent_vals[5]))
        self.settings_layout.addRow("Wlen:",self.wlen_text)
        self.relheight_text = QLineEdit("None")
        if len(self.recent_vals) > 0:
            self.relheight_text.setText(str(self.recent_vals[6]))
        self.settings_layout.addRow("RelHeight:",self.relheight_text)
        self.plateau_text = QLineEdit("None")
        if len(self.recent_vals) > 0:
            self.plateau_text.setText(str(self.recent_vals[7]))
        self.settings_layout.addRow("PlateauSize:",self.plateau_text)
        self.settings_main_layout.addLayout(self.settings_layout)
        self.save_settings_layout = QVBoxLayout()
        self.save_settings_layout.setAlignment(Qt.AlignRight)
        self.save_settings_btn = QPushButton("Find")
        if self.batch == False:
            self.export_btn = QPushButton("Find and Export")
            self.save_settings_layout.addWidget(self.export_btn)
            self.export_btn.clicked.connect(self.validate_and_export)
        else:
            self.export_batch_btn = QPushButton("Run on Batch and Export")
            self.save_settings_layout.addWidget(self.export_batch_btn)
            self.export_batch_btn.clicked.connect(self.validate_and_export_batch)
        self.save_settings_layout.addWidget(self.save_settings_btn)
        self.settings_main_layout.addLayout(self.save_settings_layout)
        
        # set the window's main layout
        self.settings_main_widget.setLayout(self.settings_main_layout)
        
        # use main stylesheet
        self.setStyleSheet(STYLESHEET)
        
        self.save_settings_btn.clicked.connect(self.validate)

        self.disable_btn_sig.connect(self.disable_btn)
     
    def disable_btn(self):
        self.save_settings_btn.setEnabled(False)
        try:
            self.export_btn.setEnabled(False)
        except:
            self.export_batch_btn.setEnabled(False)
        
    def validate(self):
        self.disable_btn_sig.emit()
        QApplication.processEvents() 
        self.height = self.height_text.text()
        self.threshold = self.threshold_text.text()
        self.distance = self.distance_text.text()
        self.prominence = self.prominence_text.text()
        self.width = self.width_text.text()
        self.wlen = self.wlen_text.text()
        self.relheight = self.relheight_text.text()
        self.plateau = self.plateau_text.text()
#        print([self.height,self.threshold,self.distance,self.prominence,self.width,self.wlen,
#                                        self.relheight,self.plateau])
        # last element for export
        self.got_advanced_options_sig.emit([[self.height,self.threshold,self.distance,self.prominence,self.width,self.wlen,
                                        self.relheight,self.plateau, False], (self.export_path,self.export_file_begin)])
            
    def validate_and_export(self):
        self.disable_btn_sig.emit()
        QApplication.processEvents() 
        self.height = self.height_text.text()
        self.threshold = self.threshold_text.text()
        self.distance = self.distance_text.text()
        self.prominence = self.prominence_text.text()
        self.width = self.width_text.text()
        self.wlen = self.wlen_text.text()
        self.relheight = self.relheight_text.text()
        self.plateau = self.plateau_text.text()
#        print([self.height,self.threshold,self.distance,self.prominence,self.width,self.wlen,
#                                        self.relheight,self.plateau])

        # last element for export
        self.got_advanced_options4export_sig.emit([[self.height,self.threshold,self.distance,self.prominence,self.width,self.wlen,
                                            self.relheight,self.plateau, True], (self.export_path,self.export_file_begin)])
            
    def validate_and_export_batch(self):
        self.disable_btn_sig.emit()
        QApplication.processEvents() 
        self.height = self.height_text.text()
        self.threshold = self.threshold_text.text()
        self.distance = self.distance_text.text()
        self.prominence = self.prominence_text.text()
        self.width = self.width_text.text()
        self.wlen = self.wlen_text.text()
        self.relheight = self.relheight_text.text()
        self.plateau = self.plateau_text.text()
        # send ordered list of parameters for scipy.find_peaks, and last element for export
        # [self.height,self.threshold,self.distance,self.prominence,self.width,self.wlen,
        # self.relheight,self.plateau,export]
        self.options2send = [[self.height,self.threshold,self.distance,self.prominence,self.width,self.wlen,
                                            self.relheight,self.plateau, True], (self.export_path,self.export_file_begin)]
        self.got_advanced_options_sig.emit(self.options2send)

        
        
    # show info popup
    def show_info_dialog(self, text):
        msgBox = QMessageBox()
        msgBox.setWindowIcon(QtGui.QIcon(ICO))
        msgBox.setText(text)
        msgBox.setWindowTitle("Info!")
        msgBox.setStandardButtons(QMessageBox.Ok)
        msgBox.exec()
        
    def closeEvent(self, event):  
        # on close, enable back parent's buttons
        if self.parent_window != None:
            self.parent_window.enable_buttons()
        
        
########################### end AdvancedPeakSettingsWindow class   
            
#####################################
# CLASS FOR PEAK DETECTION WINDOW  #
####################################
        
class PeakSettingsWindow(QMainWindow):
    
    got_peak_options_sig = pyqtSignal(list)
    disable_btns_sig = pyqtSignal()
    peak_window_closed_sig = pyqtSignal()
    def __init__(self,parent_window,subject_name,recent_peak_values,export_loc,batch_peaks):
        super(PeakSettingsWindow, self).__init__()
        self.setWindowTitle("Find Peak Settings")
        self.setWindowIcon(QtGui.QIcon(ICO))
        self.resize(400,300)
        self.parent_window = parent_window
        self.parent_window.app_closed.connect(self.exit_app)
        self.advanced_window = None
        self.options2send = []
        # [self.height,self.threshold,self.distance,self.prominence,self.width,self.wlen,
        # self.relheight,self.plateau]
        self.recent_vals = recent_peak_values
        self.export_path,self.export_file_begin = export_loc
        self.subject_name = subject_name
        self.batch = batch_peaks
        
        # use docks from pyqtgraph
        self.area = DockArea()
        self.setCentralWidget(self.area)
        self.top_dock_widget = Dock("Dock1")
        self.top_dock_widget.hideTitleBar()
        self.area.addDock(self.top_dock_widget,'left')
        self.bottom_dock_widget = Dock("Dock2")
        self.bottom_dock_widget.hideTitleBar()
        self.area.addDock(self.bottom_dock_widget, 'bottom', self.top_dock_widget)  # place the bottom dock at bottom edge of top dock       

        
        
        # create widget with settings
        self.settings_main_widget = QWidget()
#        self.setCentralWidget(self.settings_main_widget)
        # create main layout for widget with settings
        self.settings_main_layout = QVBoxLayout()
        self.settings_main_layout.setContentsMargins(10,10,10,10)
        self.settings_layout = QFormLayout()
        self.prominence_text = QLineEdit("None")
        if len(self.recent_vals) > 0:
            self.prominence_text.setText(str(self.recent_vals[3]))
        self.settings_layout.addRow("Prominence:",self.prominence_text)
        self.distance_text = QLineEdit("None")
        if len(self.recent_vals) > 0:
            self.distance_text.setText(str(self.recent_vals[2]))
        self.settings_layout.addRow("Distance (sec):",self.distance_text)
        
        self.settings_main_layout.addLayout(self.settings_layout)
        self.save_settings_layout = QVBoxLayout()
        self.save_settings_layout.setAlignment(Qt.AlignRight)
        self.save_settings_btn = QPushButton("Find")
        self.advance_settings_btn = QPushButton("Advanced Settings")
        self.save_settings_layout.addWidget(self.save_settings_btn)
        self.save_settings_layout.addWidget(self.advance_settings_btn)
        self.settings_main_layout.addLayout(self.save_settings_layout)
        
        # add that widget to the dock
        self.top_dock_widget.addWidget(self.settings_main_widget)
        self.settings_main_widget.setLayout(self.settings_main_layout)
        
        # create widget for bottom dock
        self.bottom_widget = QWidget()
        
        # create main layout for widget
        self.bottom_layout = QVBoxLayout()
        self.bottom_layout.setContentsMargins(10,10,10,10)
        
        if self.batch == False:
            self.export_loc_layout = QFormLayout()
            self.select_folder_btn = QPushButton(" Select Export Folder ")
            self.selected_folder_text = QLineEdit("")
            # always start with default npath
#            if len(self.export_path) > 0:
#                self.selected_folder_text.setText(self.export_path)
#            else: # create default subfolder
            default_path = os.path.join(self.parent_window.preview_params[1][self.subject_name],DEFAULT_EXPORT_FOLDER)
            self.selected_folder_text.setText(default_path)
            
            self.export_loc_layout.addRow(self.select_folder_btn,self.selected_folder_text)
            self.suggested_file_beginning_text = QLineEdit(self.subject_name)
#            if len(self.export_file_begin) > 0:
#                self.suggested_file_beginning_text.setText(self.export_file_begin)
            self.export_loc_layout.addRow("Suggested file name beginning:",self.suggested_file_beginning_text)
            self.bottom_layout.addLayout(self.export_loc_layout)
            
            self.export_btn = QPushButton("Find and Export")
            self.bottom_layout.addWidget(self.export_btn)
            self.select_folder_btn.clicked.connect(self.select_folder_btn_clicked)
            self.export_btn.clicked.connect(self.validate_and_export)
            self.export_btn.setStyleSheet(STYLESHEET)
        else:
            self.export_batch_btn = QPushButton("Run on Batch and Export")
            self.bottom_layout.addWidget(self.export_batch_btn)
            self.export_batch_btn.clicked.connect(self.validate_and_export_batch)
            self.export_batch_btn.setStyleSheet(STYLESHEET)
            
        self.bottom_dock_widget.addWidget(self.bottom_widget)
        self.bottom_widget.setLayout(self.bottom_layout)

        # use main stylesheet
        self.settings_main_widget.setStyleSheet(STYLESHEET)
        

        
        self.save_settings_btn.clicked.connect(self.validate)
        self.advance_settings_btn.clicked.connect(self.show_advanced_options)
        self.disable_btns_sig.connect(self.disable_buttons)
        
    def show_advanced_options(self):
        self.disable_btns_sig.emit()
        QApplication.processEvents() 
        # read parent window first
        self.distance = self.distance_text.text()
        self.prominence = self.prominence_text.text()
        try:
            self.export_path = self.selected_folder_text.text()
            self.export_file_begin = self.suggested_file_beginning_text.text()
        except:
            pass # if it doesn't exist in this view do nothing
        if len(self.recent_vals) > 0:
            self.recent_vals[2] = self.distance
            self.recent_vals[3] = self.prominence
        else:
            self.recent_vals = ["None","None",self.distance,self.prominence,"None","None",
                                        "None","None", False]
        self.advanced_window  = AdvancedPeakSettingsWindow(self,self.recent_vals,(self.export_path,self.export_file_begin),self.batch)
        self.advanced_window.got_advanced_options_sig.connect(self.advanced_options_received)
        self.advanced_window.got_advanced_options4export_sig.connect(self.advanced_options_received2export)
        self.advanced_window.show()
        
    def validate(self):
        self.disable_btns_sig.emit()
        QApplication.processEvents() 
        self.distance = self.distance_text.text()
        self.prominence = self.prominence_text.text()
#        self.export_path = self.selected_folder_text.text()
#        self.export_file_begin = self.suggested_file_beginning_text.text()
        # send ordered list of parameters for scipy.find_peaks, and last element for export
        # [self.height,self.threshold,self.distance,self.prominence,self.width,self.wlen,
        # self.relheight,self.plateau]
        self.options2send = [["None","None",self.distance,self.prominence,"None","None",
                                        "None","None", False], (self.export_path,self.export_file_begin)]
        self.got_peak_options_sig.emit(self.options2send)
        
    def validate_and_export(self):
        self.disable_btns_sig.emit()
        QApplication.processEvents() 
        self.distance = self.distance_text.text()
        self.prominence = self.prominence_text.text()

        self.export_path = self.selected_folder_text.text()
        self.export_file_begin = self.suggested_file_beginning_text.text()
        if len(self.export_path) > 0:
            # send ordered list of parameters for scipy.find_peaks, and last element for export
            # [self.height,self.threshold,self.distance,self.prominence,self.width,self.wlen,
            # self.relheight,self.plateau]
            self.options2send = [["None","None",self.distance,self.prominence,"None","None",
                                            "None","None", True], (self.export_path,self.export_file_begin)]
            self.got_peak_options_sig.emit(self.options2send)
        else:
            self.show_info_dialog("Please enter valid export path")
            
    def validate_and_export_batch(self):
        self.disable_btns_sig.emit()
        QApplication.processEvents() 
        self.distance = self.distance_text.text()
        self.prominence = self.prominence_text.text()
        # send ordered list of parameters for scipy.find_peaks, and last element for export
        # [self.height,self.threshold,self.distance,self.prominence,self.width,self.wlen,
        # self.relheight,self.plateau]
        self.options2send = [["None","None",self.distance,self.prominence,"None","None",
                                        "None","None", True], (self.export_path,self.export_file_begin)]
        self.got_peak_options_sig.emit(self.options2send)
        
    def select_folder_btn_clicked(self):
        # set text field to the value of selected path
        self.selected_folder_text.setText(QFileDialog.getExistingDirectory(self,"Choose Folder to Export Data"))

    def disable_buttons(self):
        self.save_settings_btn.setEnabled(False)
        try:
            self.export_btn.setEnabled(False)
            self.select_folder_btn.setEnabled(False)
        except:
            self.export_batch_btn.setEnabled(False)
        self.advance_settings_btn.setEnabled(False)

     
    def enable_buttons(self):
        self.save_settings_btn.setEnabled(True)
        try:
            self.export_btn.setEnabled(True)
            self.select_folder_btn.setEnabled(True)
        except:
            self.export_batch_btn.setEnabled(True)
        self.advance_settings_btn.setEnabled(True)
            
    @pyqtSlot(list)
    def advanced_options_received(self,advanced_options):
        self.options2send = advanced_options
        self.got_peak_options_sig.emit(self.options2send)
        
    @pyqtSlot(list)
    # if user clicked on export button check export path
    def advanced_options_received2export(self,advanced_options):
        self.options2send = advanced_options
        export_path = self.options2send[1][0]
        if len(export_path) > 0:
            self.got_peak_options_sig.emit(self.options2send)
        else:
            self.show_info_dialog("Please enter valid export path")
            if self.advanced_window != None:
                self.advanced_window.close()
        
    # show info popup
    def show_info_dialog(self, text):
        msgBox = QMessageBox()
        msgBox.setWindowIcon(QtGui.QIcon(ICO))
        msgBox.setText(text)
        msgBox.setWindowTitle("Info!")
        msgBox.setStandardButtons(QMessageBox.Ok)
        msgBox.exec()
        
    def closeEvent(self, event): 
        self.peak_window_closed_sig.emit()
        # on close, enable back parent's buttons
        if self.parent_window.preview_widget != None:
            self.parent_window.preview_widget.enable_buttons()
        if self.advanced_window != None:
            self.advanced_window.close()

    def exit_app(self):
        self.close()
        
        
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
            current_subject = self.parent_window.preview_widget.options["subject"]
            self.current_fs = fpExplorer_functions.get_frequency(self.parent_window.preview_widget.raw_data_dict[current_subject],self.parent_window.preview_widget.preview_init_params[0][0]["signal_name"])
            # print("Frequency:",self.current_fs)

        self.suggested_downsample_samples = int(int(self.current_fs)*DEFAULT_DOWNSAMPLE_PCT/100)
        self.suggested_downsample_rate = int(int(self.current_fs)/self.suggested_downsample_samples)
        self.min_downsampe_rate = int(int(self.current_fs)*MIN_DOWNSAMPLE_PCT/100)
        self.max_downsample_rate = int(int(self.current_fs)*MAX_DOWNSAMPLE_PCT/100)
        
        # create widget with settings
        self.settings_main_widget = QWidget()
        self.setCentralWidget(self.settings_main_widget)
        # create main layout for widget with settings
        self.settings_main_layout = QVBoxLayout()
        self.settings_main_layout.setContentsMargins(10,10,10,10)
        self.settings_layout = QFormLayout()
        self.settings_layout.setContentsMargins(10,10,10,10)
        self.settings_layout.setVerticalSpacing(15)
        self.downsample_text = QLineEdit(str(int(int(self.current_fs)/self.settings[0]["downsample"])))
        self.downsample_text.setValidator(QtGui.QIntValidator())
        # self.downsample_text.setToolTip("Integers from 2 to "+str(MAX_DOWNSAMPLE))
        self.downsample_text.setToolTip("Between "+str(self.min_downsampe_rate)+" and "+str(self.max_downsample_rate))
        # self.settings_layout.addRow("Downsample (How many samples to average)\nRecommended: 1-2% of sampling frequency",self.downsample_text)
        downsample_label = "Downsample to lower sampling rate.\nOriginal: "+str(int(self.current_fs))+"    Suggested: "+str(self.suggested_downsample_rate)
        self.settings_layout.addRow(downsample_label,self.downsample_text)
        self.normalization_method_comboBox = QComboBox()
        self.normalization_method_comboBox.addItems(["Standard Polynomial Fitting","Modified Polynomial Fitting"])
        self.normalization_method_comboBox.setCurrentText(self.settings[0]["normalization"])
        self.settings_layout.addRow("Method of normalization",self.normalization_method_comboBox)
        self.normalization_show_comboBox = QComboBox()
        self.normalization_show_comboBox.addItems(["df/F ( in % )","Z-Score"])
        self.normalization_show_comboBox.setCurrentText(self.settings[0]["show_norm_as"])
        self.settings_layout.addRow("Show normalized data as",self.normalization_show_comboBox)
        self.filter_cb = QCheckBox("Gaussian filter")
        self.filter_cb.setChecked(self.settings[0]["filter"])
        self.settings_layout.addRow("Smooth data:",QLabel("The fraction (0-100%) of the data used when estimating each y-value:"))
        # The fraction of the data used when estimating each y-value
#        self.smooth_text = QLineEdit(str(DEFAULT_SMOOTH_FRAQ*100))
        self.smooth_text = QLineEdit(str(self.settings[0]["filter_fraction"]))
        self.smooth_text.setToolTip("From 0 to "+str(MAX_SMOOTH_FRAQ))
        self.settings_layout.addRow(self.filter_cb,self.smooth_text)
        self.settings_main_layout.addLayout(self.settings_layout)
        self.save_settings_layout = QVBoxLayout()
        self.save_settings_layout.setAlignment(Qt.AlignRight)
        self.save_settings_btn = QPushButton("Save settings")
        self.save_settings_layout.addWidget(self.save_settings_btn)
        self.settings_main_layout.addLayout(self.save_settings_layout)
        
        # set the window's main layout
        self.settings_main_widget.setLayout(self.settings_main_layout)
        
        # use main stylesheet
        self.setStyleSheet(STYLESHEET)
        
        self.save_settings_btn.clicked.connect(self.save_settings_btn_clicked)
        
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
            samples = int(int(self.current_fs)/int(self.downsample_text.text()))
            self.settings[0]["downsample"] = samples
            self.settings[0]["entered_downsample"] = int(self.downsample_text.text())
        else:
            self.show_info_dialog("Downsample was not updated.\nEnter values between "+str(self.min_downsampe_rate)+" and " + str(self.max_downsample_rate))
        # don't loose filter fraq information
        if self.filter_cb.isChecked() == False:
            self.settings[0]["filter_fraction"] = DEFAULT_SMOOTH_FRAQ
        else:
            try:
                fraq = float(self.smooth_text.text())
                if fraq > 0 and fraq <= MAX_SMOOTH_FRAQ:
                    self.settings[0]["filter_fraction"] = fraq
                else:
                    self.show_info_dialog("Fraction of the data has to be between 0 and "+str(MAX_SMOOTH_FRAQ))
            except:
                self.show_info_dialog("Fraction of the data has to be between 0 and "+str(MAX_SMOOTH_FRAQ))
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
        
#################################################
# CLASS FOR PERIEVENT OPTIONS WINDOW           #
################################################
        
class PeriEventOptionsWindow(QMainWindow):
    # pyqt Signal has to be defined on class level
    # and not in init method !!
    # create a list with information acquired from window
    # and send it as a signal
    got_peri_event_options_sig = pyqtSignal(list)
    disable_buttons_signal = pyqtSignal()
    def __init__(self,parent_window,subject_info,perievent_options,batch,common_events):
        super(PeriEventOptionsWindow, self).__init__()
        self.setWindowTitle("Peri-event options")
        self.setWindowIcon(QtGui.QIcon(ICO))
        self.resize(500,500)
        self.parent_window = parent_window
        self.parent_window.app_closed.connect(self.exit_app)
        # list; first el=subject
        # second=list of possible events
        # third el=raw data for current subject
        self.subject_info = subject_info
        
        # dictionery with user entered options
        self.options_dict = perievent_options
        # reset preview button clicked to false
        # reset analyze button clicked to false
        # reset export button clicked to false
        self.options_dict["preview"] = False
        self.options_dict["analyze"] = False
        self.options_dict["export"] = False
        
        self.batch = batch
        
        self.export_path = self.options_dict["export_path"]
        self.export_file_begin = self.options_dict["file_beginning"]
        
        # use docks from pyqtgraph
        self.area = DockArea()
        self.setCentralWidget(self.area)
        self.top_dock_widget = Dock("Dock1")
        self.top_dock_widget.hideTitleBar()
        self.area.addDock(self.top_dock_widget,'left')
        self.bottom_dock_widget = Dock("Dock2")
        self.bottom_dock_widget.hideTitleBar()
        self.area.addDock(self.bottom_dock_widget, 'bottom', self.top_dock_widget)  # place the bottom dock at bottom edge of top dock       
        self.export_dock_widget = Dock("Dock3")
        self.export_dock_widget.hideTitleBar()
        self.area.addDock(self.export_dock_widget, 'bottom', self.bottom_dock_widget) 
        
        # create widget for top dock
        self.main_widget = QWidget()
        
        # create main layout for widget with settings
        self.main_layout = QVBoxLayout()
        self.main_layout.setContentsMargins(10,10,10,10)
        self.event_layout = QHBoxLayout()
        self.event_layout.setContentsMargins(10,10,10,10)
        self.event_label = QLabel("Select Your Event:")
        self.event_from_file_comboBox = QComboBox()
        if self.batch == True:
            print(common_events)
            self.event_from_file_comboBox.addItems(common_events)
        else:
            self.event_from_file_comboBox.addItems(self.subject_info[1])
        if 'event' in self.options_dict:
            self.event_from_file_comboBox.setCurrentText(self.options_dict['event'])
        self.event_name_label = QLabel( "Event name (optional):")
        self.event_name_text = QLineEdit("")
        if 'event_name' in self.options_dict:
            self.event_name_text.setText(self.options_dict['event_name'])
        self.event_name_text.setToolTip("i.e, Tone")
        self.event_layout.addWidget(self.event_label)
        self.event_layout.addWidget(self.event_from_file_comboBox)
        self.event_layout.addWidget(self.event_name_label)
        self.event_layout.addWidget(self.event_name_text)
        self.main_layout.addLayout(self.event_layout)
        
        self.message_layout = QHBoxLayout()
        self.message_layout.setContentsMargins(10,0,0,0)
        self.message_layout.addWidget(QLabel( "Here you select your data window for all Peri-event analysis!\n"))
        self.main_layout.addLayout(self.message_layout)
        self.event_window_layout = QHBoxLayout()
        self.event_window_layout.setContentsMargins(10,0,10,10)
        self.before_sec_label = QLabel("Time before event (sec):")
        self.before_sec_text = QLineEdit("30")
        if 'sec_before' in self.options_dict:
            self.before_sec_text.setText(str(self.options_dict['sec_before']))
        self.before_sec_text.setValidator(QtGui.QIntValidator())
        self.after_sec_label = QLabel("Time after event (sec):")
        self.after_sec_text = QLineEdit("30")
        if 'sec_after' in self.options_dict:
            self.after_sec_text.setText(str(self.options_dict['sec_after']))
        self.after_sec_text.setValidator(QtGui.QIntValidator())
        self.event_window_layout.addWidget(self.before_sec_label)
        self.event_window_layout.addWidget(self.before_sec_text)
        self.event_window_layout.addWidget(self.after_sec_label)
        self.event_window_layout.addWidget(self.after_sec_text)
        self.main_layout.addLayout(self.event_window_layout)
        
        self.plot_raw_btn = QPushButton("Preview All Events")
        self.main_layout.addWidget(self.plot_raw_btn)
        
        # add that widget to the dock
        self.top_dock_widget.addWidget(self.main_widget)
        self.main_widget.setLayout(self.main_layout)
        
        # create widget for bottom dock
        self.bottom_widget = QWidget()
        
        # create main layout for widget
        self.bottom_layout = QVBoxLayout()
        self.bottom_layout.setContentsMargins(10,10,10,10)
        self.baseline_label = QLabel("\n   BASELINE and AUC (Area Under the Curve):\n\n   "+
                                     "Select time window for baseline. And time windows for AUC-pre and AUC-post.\n   "+
                                     "Time windows for AUC pre and post have to be the same size.\n   "+
                                     "Use negative integers to indicate time before the event. Time is in seconds.\n")
        self.bottom_layout.addWidget(self.baseline_label)
        self.baseline_layout = QGridLayout()
        self.baseline_layout.setVerticalSpacing(15)
        self.baseline_layout.setHorizontalSpacing(15)
        self.baseline_layout.setContentsMargins(10,10,10,10)
        self.baseline_from_label = QLabel("Baseline From:")
        self.baseline_from_text = QLineEdit("-30")
        if 'baseline_from' in self.options_dict:
            self.baseline_from_text.setText(str(self.options_dict["baseline_from"]))
        self.baseline_from_text.setValidator(QtGui.QIntValidator())
        self.baseline_to_label = QLabel("Baseline To:")
        self.baseline_to_text = QLineEdit("-10")
        if 'baseline_to' in self.options_dict:
            self.baseline_to_text.setText(str(self.options_dict["baseline_to"]))
        self.baseline_to_text.setValidator(QtGui.QIntValidator())
        self.baseline_layout.addWidget(self.baseline_from_label,0,0)
        self.baseline_layout.addWidget(self.baseline_from_text,0,1)
        self.baseline_layout.addWidget(self.baseline_to_label,0,2)
        self.baseline_layout.addWidget(self.baseline_to_text,0,3)
        self.auc_pre_from_label = QLabel("AUC-pre From:")     
        self.auc_pre_from_text = QLineEdit("-30")
        if 'auc_pre_from' in self.options_dict:
            self.auc_pre_from_text.setText(str(self.options_dict["auc_pre_from"]))
        self.auc_pre_from_text.setValidator(QtGui.QIntValidator())
        self.auc_pre_to_label = QLabel("AUC-pre To:")
        self.auc_pre_to_text = QLineEdit("0")
        if 'auc_pre_to' in self.options_dict:
            self.auc_pre_to_text.setText(str(self.options_dict["auc_pre_to"]))
        self.auc_pre_to_text.setValidator(QtGui.QIntValidator())
        self.auc_post_from_label = QLabel("AUC-post From:")
        self.auc_post_from_text = QLineEdit("0")
        if 'auc_post_from' in self.options_dict:
            self.auc_post_from_text.setText(str(self.options_dict["auc_post_from"]))
        self.auc_post_from_text.setValidator(QtGui.QIntValidator())
        self.auc_post_to_label = QLabel("AUC-post To:")
        self.auc_post_to_text = QLineEdit("30")
        if 'auc_post_to' in self.options_dict:
            self.auc_post_to_text.setText(str(self.options_dict["auc_post_to"]))
        self.auc_post_to_text.setValidator(QtGui.QIntValidator())
        self.baseline_layout.addWidget(self.auc_pre_from_label,1,0)
        self.baseline_layout.addWidget(self.auc_pre_from_text,1,1)
        self.baseline_layout.addWidget(self.auc_pre_to_label,1,2)
        self.baseline_layout.addWidget(self.auc_pre_to_text,1,3)
        self.baseline_layout.addWidget(self.auc_post_from_label,2,0)
        self.baseline_layout.addWidget(self.auc_post_from_text,2,1)
        self.baseline_layout.addWidget(self.auc_post_to_label,2,2)
        self.baseline_layout.addWidget(self.auc_post_to_text,2,3)      
        self.select_plot_label = QLabel("Show on plot:")
        self.plot_average_cb = QCheckBox("Average")
        self.plot_average_cb.setChecked(True)
        if "plot_avg" in self.options_dict:
            self.plot_average_cb.setChecked(self.options_dict["plot_avg"])
        self.plot_zscore_cb = QCheckBox("Z-Score with error")
        if "plot_zscore" in self.options_dict:
            self.plot_zscore_cb.setChecked(self.options_dict["plot_zscore"])
        self.plot_zscore_with_trials_cb = QCheckBox("Z-Score with trials")
        if "plot_zscore_trials" in self.options_dict:
            self.plot_zscore_with_trials_cb.setChecked(self.options_dict["plot_zscore_trials"])
        self.plot_auc_cb = QCheckBox("Area Under the Curve (AUC)")  
        if "plot_auc" in self.options_dict:
            self.plot_auc_cb.setChecked(self.options_dict["plot_auc"])
        self.baseline_layout.addWidget(self.select_plot_label,3,0)
        self.baseline_layout.addWidget(self.plot_average_cb,3,1)
        self.baseline_layout.addWidget(self.plot_zscore_cb,3,2) 
        self.baseline_layout.addWidget(self.plot_auc_cb,3,3)
        self.baseline_layout.addWidget(self.plot_zscore_with_trials_cb,4,2)
        self.bottom_layout.addLayout(self.baseline_layout)
        
        self.analyze_btn = QPushButton("Analyze")
        if self.batch == False:
            self.bottom_layout.addWidget(self.analyze_btn)
        
        
        # add that widget to the dock
        self.bottom_dock_widget.addWidget(self.bottom_widget)
        self.bottom_widget.setLayout(self.bottom_layout)
        
        # add widget for export destination and button
        self.export_widget = QWidget()
        # create main layout for widget
        self.export_layout = QVBoxLayout()
        self.export_layout.setContentsMargins(10,10,10,10)
        
        if self.batch == False:
            self.export_loc_layout = QFormLayout()
            self.export_loc_layout.setContentsMargins(10,10,10,10)
            self.select_folder_btn = QPushButton(" Select Export Folder ")
            self.selected_folder_text = QLineEdit('')
#            if len(self.export_path) > 0:
#                self.selected_folder_text.setText(self.export_path)
#            else: # create default subfolder
            # always create default subfolder
            default_path = os.path.join(self.parent_window.preview_params[1][self.subject_info[0]],DEFAULT_EXPORT_FOLDER)
            self.selected_folder_text.setText(default_path)
            self.export_loc_layout.addRow(self.select_folder_btn,self.selected_folder_text)
            self.suggested_file_beginning_text = QLineEdit(subject_info[0])
#            if len(self.export_file_begin) > 0:
#                self.suggested_file_beginning_text.setText(self.export_file_begin)
            self.export_loc_layout.addRow("Suggested file name beginning:",self.suggested_file_beginning_text)
            self.export_layout.addLayout(self.export_loc_layout)
            self.export_btn = QPushButton("Export Data")
            self.export_layout.addWidget(self.export_btn)
            
            self.select_folder_btn.clicked.connect(self.select_folder_btn_clicked)
            self.export_btn.clicked.connect(self.export_btn_clicked)
            self.export_btn.setStyleSheet(STYLESHEET)
        else:
            self.export_batch_btn = QPushButton("Run on Batch and Export")
            self.export_layout.addWidget(self.export_batch_btn)
            self.export_batch_btn.clicked.connect(self.analyze_btn_clicked)
            self.export_batch_btn.setStyleSheet(STYLESHEET)
        
        self.export_widget.setLayout(self.export_layout)
        self.export_dock_widget.addWidget(self.export_widget)
        
        # use main stylesheet
        self.main_widget.setStyleSheet(STYLESHEET)
        self.bottom_widget.setStyleSheet(STYLESHEET)
        
        
        self.plot_raw_btn.clicked.connect(self.plot_raw_btn_clicked)
        self.analyze_btn.clicked.connect(self.analyze_btn_clicked)
        self.plot_zscore_cb.stateChanged.connect(self.zscore_toggle)
        self.plot_zscore_with_trials_cb.stateChanged.connect(self.zscore_with_trials_toggle)
        
        
        self.disable_buttons_signal.connect(self.disable_buttons)
        
    def zscore_toggle(self):
        if self.plot_zscore_cb.isChecked():
            self.plot_zscore_with_trials_cb.setChecked(False)
            
    def zscore_with_trials_toggle(self):
        if self.plot_zscore_with_trials_cb.isChecked():
            self.plot_zscore_cb.setChecked(False)
        
    def disable_buttons(self):
        self.plot_raw_btn.setEnabled(False)
        self.analyze_btn.setEnabled(False)
        try:
            self.export_btn.setEnabled(False)
            self.select_folder_btn.setEnabled(False)
        except:
            self.export_batch_btn.setEnabled(False)
        QApplication.processEvents() 
        
    def enable_buttons(self):
        self.plot_raw_btn.setEnabled(True)
        self.analyze_btn.setEnabled(True)
        try:
            self.export_btn.setEnabled(True)
            self.select_folder_btn.setEnabled(True)
        except:
            self.export_batch_btn.setEnabled(True)
        QApplication.processEvents() 
        
    
    def plot_raw_btn_clicked(self):
        self.disable_buttons_signal.emit()
        self.read_options()
        if len(self.options_dict["event"])==0:
            self.show_info_dialog("No events were found to run the analysis")
            self.close()
        # set preview to true
        self.options_dict["preview"] = True
        # send values to main app window
        self.got_peri_event_options_sig[list].emit([self.options_dict])
#        self.close()
        
    def analyze_btn_clicked(self):
        self.disable_buttons_signal.emit()
        self.read_options()
        if self.validate() == True:
            # set analyze to true
            self.options_dict["analyze"] = True
            # send values to main app window
            self.got_peri_event_options_sig[list].emit([self.options_dict])
#            self.close()
    
    def export_btn_clicked(self):
        self.disable_buttons_signal.emit()
        self.read_options()
        if self.validate() == True:
            # check if user defined path
            if len(self.export_path) > 0:
                # set export to true
                self.options_dict["export"] = True
                self.options_dict["export_path"] = self.export_path
                self.options_dict["file_beginning"] = self.export_file_begin
                # send values to main app window
                self.got_peri_event_options_sig[list].emit([self.options_dict])
#                self.close()
            else:
                self.show_info_dialog("Please select export path.")
            
    def validate(self):
        valid = False
        # make sure that for analysis, baseline is correct
        if self.options_dict["baseline_from"] == self.options_dict["baseline_to"]:
            self.show_info_dialog("Baseline window needs to me longer than 0 seconds.")
        elif self.options_dict["baseline_from"] > self.options_dict["baseline_to"]:
            self.show_info_dialog("Start of baseline window needs to be earlier\nthan the end of baseline window.")
        # check if baseline is within perievent window
        if self.options_dict["baseline_from"] < -self.options_dict["sec_before"] or self.options_dict["baseline_to"] > self.options_dict["sec_after"]:
            self.show_info_dialog("Baseline time window needs to be\nwithin perievent time window.")
        # check if auc windows are within range
        if (self.options_dict["auc_pre_from"] < -self.options_dict["sec_before"] or self.options_dict["auc_pre_from"] > self.options_dict["sec_after"] 
            or self.options_dict["auc_post_to"] > self.options_dict["sec_after"] or self.options_dict["auc_post_to"] < -self.options_dict["sec_before"]):
            self.show_info_dialog("Please make sure that your AUC time window\nis within perievent time window.")
        # check id auc windows are equal
        elif abs(self.options_dict["auc_pre_from"]+self.options_dict["auc_pre_to"]) != abs(self.options_dict["auc_post_from"]+self.options_dict["auc_post_to"]):
            self.show_info_dialog("Please make sure the AUC-pre and AUC-post\ntime windows are equal.")
        else:
            valid = True
        # this should always be true, unless there are no events!
        if len(self.options_dict["event"])==0:
            valid = False
            self.show_info_dialog("No events were found to run the analysis")
            self.close()
        return valid
    
    def read_options(self):
        self.options_dict["event"] = self.event_from_file_comboBox.currentText()
        self.options_dict["event_name"] = self.event_name_text.text()
        self.options_dict["plot_avg"] = self.plot_average_cb.isChecked()
        self.options_dict["plot_zscore"] = self.plot_zscore_cb.isChecked()
        self.options_dict["plot_zscore_trials"] = self.plot_zscore_with_trials_cb.isChecked()
        self.options_dict["plot_auc"] = self.plot_auc_cb.isChecked()
        try:
            self.options_dict["sec_before"] = abs(int(self.before_sec_text.text()))
            self.options_dict["sec_after"] = abs(int(self.after_sec_text.text()))
            self.options_dict["baseline_from"] = int(self.baseline_from_text.text())
            self.options_dict["baseline_to"] = int(self.baseline_to_text.text())
            self.options_dict["auc_pre_from"] = int(self.auc_pre_from_text.text())
            self.options_dict["auc_pre_to"] = int(self.auc_pre_to_text.text())
            self.options_dict["auc_post_from"] = int(self.auc_post_from_text.text())
            self.options_dict["auc_post_to"] = int(self.auc_post_to_text.text())
        except:
            self.show_info_dialog("Enter valid intigers for\nseconds before and after event.")
            self.options_dict["sec_before"] = 0
            self.options_dict["sec_after"] = 0
            self.options_dict["baseline_from"] = 0
            self.options_dict["baseline_to"] = 0
        try:
            self.export_path = self.selected_folder_text.text()
            self.export_file_begin = self.suggested_file_beginning_text.text()
        except:
            pass # if don't exist, do nothing
            
    def select_folder_btn_clicked(self):
        # set text field to the value of selected path
        self.selected_folder_text.setText(QFileDialog.getExistingDirectory(self,"Choose Folder to Export Data"))
            
        
    # show info popup
    def show_info_dialog(self, text):
        msgBox = QMessageBox()
        msgBox.setWindowIcon(QtGui.QIcon(ICO))
        msgBox.setText(text)
        msgBox.setWindowTitle("Info!")
        msgBox.setStandardButtons(QMessageBox.Ok)
        msgBox.exec()
        
    def closeEvent(self, event):  
        # on close, enable back parent's buttons
        if self.parent_window.preview_widget != None:
            self.parent_window.preview_widget.enable_buttons()

    def exit_app(self):
        self.close()
        
############################## end PeriEventOptionsWindow    
        
#################################################
# CLASS FOR EXPORT DATA WINDOW                 #
################################################
        
class ExportDataWindow(QMainWindow):
    
    got_export_selection_sig = pyqtSignal(list)
    def __init__(self,main_window,parent_window,params):
        super(ExportDataWindow, self).__init__()
        self.setWindowTitle("Export Data")
        self.setWindowIcon(QtGui.QIcon(ICO))
        self.resize(400,400)
        self.parent_window = parent_window
        self.parent_window.close_export_window_signal.connect(self.close_all)
        self.main_window = main_window
        self.main_window.app_closed.connect(self.exit_app)
        # list; first el=subject
        # second=list of possible events
        # third el=raw data for current subject
        self.subject_info = params[0]
        self.init_dict = params[1]
        self.export_path = params[2]
        self.file_begin = params[3]
        self.event_based = self.init_dict[0][0]["event_based"]
        
        # tbd
        self.export_options_list = []
        
        # create main widget
        self.export_main_widget = QWidget()
        self.setCentralWidget(self.export_main_widget)
        # create main layout
        self.export_main_layout = QVBoxLayout()
        self.export_main_layout.setContentsMargins(10,10,10,10)
        
        self.subject_label = QLabel("SUBJECT "+self.subject_info)
        self.subject_label.setAlignment(Qt.AlignCenter)
        self.big_label_stylesheet = "QLabel {font-size: 16pt}"
        self.subject_label.setStyleSheet(self.big_label_stylesheet)
        self.export_main_layout.addWidget(self.subject_label)
        
        self.selection_layout = QVBoxLayout()
        self.selection_layout.setContentsMargins(20,20,20,20)
        
        self.export_loc_layout = QFormLayout()
        self.export_loc_layout.setContentsMargins(20,20,20,20)
        self.export_loc_widget = QWidget()
        self.select_folder_btn = QPushButton(" Select Export Folder ")
#        self.select_folder_btn.setStyleSheet(self.bold_stylesheet)
        self.selected_folder_text = QLineEdit("")
#        if len(self.export_path) > 0:
#            self.selected_folder_text.setText(self.export_path)
#        else: # create default subfolder
        # always create default subfolder for current subject
        default_path = os.path.join(self.parent_window.parent_window.preview_params[1][self.subject_info],DEFAULT_EXPORT_FOLDER)
        self.selected_folder_text.setText(default_path)
        self.export_loc_layout.addRow(self.select_folder_btn,self.selected_folder_text)
        self.suggested_file_beginning_text = QLineEdit(self.subject_info)
        self.export_loc_layout.addRow("Suggested file name beginning:",self.suggested_file_beginning_text)
        self.export_main_layout.addLayout(self.export_loc_layout)
        
        
        self.select_label = QLabel("Select Data for Export:\n")
        self.selection_layout.addWidget(self.select_label)
        self.raw_cb = QCheckBox("Raw Data")
#        self.trimmed_cb = QCheckBox("Trimmed Data")
        self.downsampled_cb = QCheckBox("Downsampled Trimmed Data")
        self.normalized_cb = QCheckBox("Normalized Data")
        self.perievent_cb = QCheckBox("Perievent Data")
        self.spike_cb = QCheckBox("Spikes")
        self.save_plots_cb = QCheckBox("Save Plots")
        self.bold_stylesheet = "QCheckBox {font-weight: bold}"
        self.save_plots_cb.setChecked(True)
        self.save_plots_cb.setStyleSheet(self.bold_stylesheet)
        self.select_all_cb = QCheckBox("Select All")
        self.select_all_cb.setStyleSheet(self.bold_stylesheet)
        self.selection_layout.addWidget(self.raw_cb)
#        self.selection_layout.addWidget(self.trimmed_cb)
        self.selection_layout.addWidget(self.downsampled_cb)
        self.selection_layout.addWidget(self.normalized_cb)
        if self.event_based == True:
            self.selection_layout.addWidget(self.perievent_cb)
        self.selection_layout.addWidget(self.spike_cb)
        self.selection_layout.addWidget(self.save_plots_cb)  
        self.selection_layout.addWidget(QLabel(""))   
        self.selection_layout.addWidget(self.select_all_cb)  
        self.export_main_layout.addLayout(self.selection_layout)
        
        self.export_btn_layout = QHBoxLayout()
        self.export_btn_layout.setContentsMargins(20,20,20,20)
        self.export_btn_layout.setAlignment(Qt.AlignRight)
        self.export_btn = QPushButton("Export")
        self.export_btn_layout.addWidget(self.export_btn)
        self.export_main_layout.addLayout(self.export_btn_layout)
        
       
        self.export_main_widget.setLayout(self.export_main_layout)
        
        self.select_folder_btn.clicked.connect(self.select_folder_btn_clicked)
        self.select_all_cb.stateChanged.connect(self.select_all_toggle)
        self.export_btn.clicked.connect(self.read_export_options)
        
    def select_all_toggle(self):
        if self.select_all_cb.isChecked():
            self.raw_cb.setChecked(True)
#            self.trimmed_cb.setChecked(True)
            self.downsampled_cb.setChecked(True)
            self.normalized_cb.setChecked(True)
#            self.event_data_cb.setChecked(True)
            self.perievent_cb.setChecked(True)
            self.spike_cb.setChecked(True)
            self.save_plots_cb.setChecked(True)
        else:
            self.raw_cb.setChecked(False)
#            self.trimmed_cb.setChecked(False)
            self.downsampled_cb.setChecked(False)
            self.normalized_cb.setChecked(False)
            self.perievent_cb.setChecked(False)
            self.spike_cb.setChecked(False)
            self.save_plots_cb.setChecked(False)
            
    def select_folder_btn_clicked(self):
        # set text field to the value of selected path
        self.selected_folder_text.setText(QFileDialog.getExistingDirectory(self,"Choose Folder to Export Data"))
            
    
    # read user options and sent to preview        
    def read_export_options(self):
        self.export_btn.setEnabled(False)
        self.select_folder_btn.setEnabled(False)
        QApplication.processEvents() 
        export_options_dict = {}
        export_options_dict["dump_path"] = self.selected_folder_text.text()
        export_options_dict["file_begin"] = self.suggested_file_beginning_text.text()
        export_options_dict["raw"] = self.raw_cb.isChecked()
#        export_options_dict["trimmed"] = self.trimmed_cb.isChecked()
        export_options_dict["downsampled"] = self.downsampled_cb.isChecked()
        export_options_dict["normalized"] = self.normalized_cb.isChecked()
        export_options_dict["perievent"] = self.perievent_cb.isChecked()
        export_options_dict["spikes"] = self.spike_cb.isChecked()
        export_options_dict["save_plot"] = self.save_plots_cb.isChecked()
        self.export_options_list.append(export_options_dict)
        self.got_export_selection_sig[list].emit(self.export_options_list)
        
    @pyqtSlot()   
    def close_all(self):
        self.close()
        
    def closeEvent(self, event):  
        # on close, enable back parent's buttons
        self.parent_window.enable_buttons()

    def exit_app(self):
        self.close()
     
############################## end ExportDataWindow  

       
#################################################
# CLASS FOR MATPLOTLIB CANVAS                  #
################################################
        
class MplCanvas(FigureCanvasQTAgg):

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        super(MplCanvas, self).__init__(self.fig)  
        self.fig.canvas.mpl_connect('scroll_event',self.mousewheel_move)
# scrolling only zooms x axis
# inspired by: https://stackoverflow.com/questions/11551049/matplotlib-plot-zooming-with-scroll-wheel        
    def mousewheel_move(self,event):
        try:
            ax=event.inaxes
            ax._pan_start = types.SimpleNamespace(
                    lim=ax.viewLim.frozen(),
                    trans=ax.transData.frozen(),
                    trans_inverse=ax.transData.inverted().frozen(),
                    bbox=ax.bbox.frozen(),
                    x=event.x,
                    y=event.y)
            if event.button == 'up':
                ax.drag_pan(3, event.key, event.x+10, event.y)
            else: #event.button == 'down':
                ax.drag_pan(3, event.key, event.x-10, event.y)
            fig=ax.get_figure()
            fig.canvas.draw_idle()
        except:
            pass
        
################################################################
#                                                              #
# EXECUTE GUI FROM MAIN                                        #
#                                                              #
################################################################
if __name__ == "__main__":
    # Always start by initializing Qt (only once per application)
    app = QApplication([])
    main_widget = MyMainWidget()
    main_widget.showMaximized()
    # show select data window on top when app is started
    main_widget.select_data_window.show()
    app.exec_()

    print('Done')