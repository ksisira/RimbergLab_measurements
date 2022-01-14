# -*- coding: utf-8 -*-
"""
Created on Fri Dec 24 12:32:23 2021

@author: Sisira

Optimized PID values for measured Y fluctuations and power fluctuations.
"""
"Import modules"

import instrument_classes_module as icm #compiled function module
import minicircuits_usb_devices as mud
import gpib_classes_module as gcm
import serial_classes_module as scm
import numpy as np #Scientific computing module
import os, sys #system functions
import shutil
import datetime as dt#date and time
import subprocess as sp #opening files in notepad
import time#for transient delay
import functions as fn
import pickle as pkl
import measurement_modules as mm
import matplotlib.pyplot as plt
from scipy import stats
import optimize_pid_simulation as ops

"""Make directories"""
current_time=dt.datetime.now()
directory_name=('../{}_{}_{}_{}_pid_optimize'.
                format(current_time.month, current_time.day,current_time.hour,
                       current_time.minute))
os.mkdir(directory_name)
scripts_dir = directory_name+'/scripts'
shutil.copytree(os.getcwd(),scripts_dir)
current_file = __file__
fn.copy_file(current_file, directory_name+'/script')

"""Copy relevant data files"""
prompt_cont = input(('''Copy the following files from relevant data folder to
                     directory. \n
                     1. Parameters.dat \n
                     2. Output_signal/fridge_out_signal.dat \n
                     3. Output_signal/fridge_out_off_resonance.dat \n
                     4. Error signal data.dat (From the ng sweep for same
                     bias point, photons and lock-in values). \n
                     5. quadrature_signal_on_resonance.dat
                     (from the lockin_fft_feedback_on folder after preliminary
                     correction)\n'''+
                         'Continue? (y/exits)'))
folder_name = input('Copy the folder name containing data.')

time_step = 1 #in ns
parameters = input('''Give following values [photons,bias,tc,vna_power,
                                    input att, kappa_int, kappa_ext].''')                  
(res_power, res_Y, out) = ops.optimize_pid(directory_name, parameters, 
folder_name)

                    
                     
