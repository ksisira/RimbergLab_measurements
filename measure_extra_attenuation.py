# -*- coding: utf-8 -*-
"""
Created on January 25 2021 5 30 PM.
Modified Dec 20th, 2021.

@author: Sisira
"""

"""
Get extra attenuation between the spectrum analyzer input measurement and
VNA input power.

Refer lab note (Part II, Page 13).

Measurements:
1. Connect from VNA Port 1 - Cable 2 - VNA Port 2.
2. Connect from VNA Port 1 - Cable 6 - VNA Port 2.
3. Connect from VNA Port 1 - Cable 10 - SP4T 1. 
   SP4T COM - Cable 9 - VNA Port 2.
"""
"Import modules"
import instrument_classes_module as icm #compiled function module
import minicircuits_usb_devices as mud
import numpy as np #Scientific computing module
import os, sys #system functions
import datetime as dt#date and time
import subprocess as sp #opening files in notepad
import matplotlib.pyplot as plt #Plotting module
import pickle as pkl
import functions as fn
from scipy import interpolate

line = fn.lineno

"""Call instruments"""
na=icm.agilent_e5071c(15) #call the network analyzer
switch_in = mud.minicircuits_rc4spdta18()
switch_sa = mud.minicircuits_rc1sp4ta18()

"""Make directories"""
current_time=dt.datetime.now()
directory_name='{}_{}_{}_attenuation'.format(
                current_time.month,
                current_time.day,current_time.year)
os.mkdir(directory_name)
current_file = __file__
fn.copy_file(current_file, directory_name+'/script')

"Parameter file"

power_array = '[-55]'
span = '[400e6]'
center_freq='[5.75e9]'
el_delay='[63.3e-9]'
date='{}_{}_{}'.format(current_time.month,
      current_time.day,current_time.year)

par_array=[]
line0 = line()
par_array.append([line()-line0,"Date of scan",date])
par_array.append([line()-line0,"Channel number",1])
par_array.append([line()-line0,"Number of traces",2])
par_array.append([line()-line0,"Number of vna sweep points",1601])
par_array.append([line()-line0,"VNA Averaging number",50])
par_array.append([line()-line0,"Input power (in dBm)",power_array])
par_array.append([line()-line0,"Span (in Hz)",span])
par_array.append([line()-line0,"Center frequency",center_freq])
par_array.append([line()-line0,"VNA IF bandwidth",10e3])
par_array.append([line()-line0,"Electrical delay",el_delay])

filename_parameter=str(directory_name)+'/Parameters.dat'
np.savetxt(filename_parameter,np.row_stack((par_array)),delimiter=';',fmt='%s')
programName = "notepad.exe"
sp.Popen([programName,filename_parameter]) #opens file in notepad
prompt=input('''Check the notepad file. Continue with the saved parameter\
             values? (y to continue/exits otherwise)''')
if prompt=='y': #prompts to check the input values
    print('Running..')
else:
    sys.exit('Script aborted. Confirm parameter values.')
    
parameters=np.loadtxt(filename_parameter,delimiter=';',usecols=2,dtype=np.str) 
#loads the parameter values

"Parameter values"

parameters = np.array([eval(i) for i in parameters])
[date,ch,num_traces,vna_points,vna_av,in_power_array,vna_span_arr,
 vna_center_freq_arr,
 vna_if_bw,el_delay_arr]=parameters

"""----------------------------------MEASUREMENT ---------------------------"""
na.set_freqs(vna_center_freq_arr[0],vna_span_arr[0],'span',ch)
na.set_power(in_power_array[0],ch)

"Measurement loop"

"""Measurement 1"""
prompt=input('''Connect the configuration for measurement 1.
             (y to continue/exits otherwise)''')
if prompt=='y': #prompts to check the input values
    print('Running..')
else:
    sys.exit('Script aborted. Confirm parameter values.')
def vna_in():
    switch_in.set_switch('A',0)
vna_in()
for inp in in_power_array:
    for freq_no in range(len(vna_center_freq_arr)):
        par_array, freq, logmag_vna, phase = fn.vna_scan(
            vna_center_freq_arr[freq_no], vna_span_arr[freq_no], inp, 
            vna_points, vna_av, 
            el_delay_arr[freq_no], vna_if_bw)
        
            
        filename_str=str(directory_name)+('/vna_line_attenuation')
        filename="%s.dat" %filename_str
        np.savetxt(filename,np.column_stack((freq, logmag_vna, phase)))
            
        fig, ax1=plt.subplots(figsize=(20,10))
        color = 'tab:red'
        ax1.plot(freq*1e-9, logmag_vna, color=color)
        ax1.set_xlabel('Frequency (in GHz)',fontsize = 30)
        ax1.set_ylabel('|S11| (in dB)',fontsize = 30,color=color)
        ax1.tick_params(axis='y', labelcolor=color)
        ax2 = ax1.twinx()
        color = 'tab:blue'
        ax2.set_ylabel('Phase (deg)', color=color,fontsize = 30)
        ax2.plot(freq*1e-9, phase, color=color)
        ax2.tick_params(axis='y', labelcolor=color)
        
        filename="%s.jpeg" %filename_str
        plt.savefig(filename)
        plt.close(fig)

"""Measurement 2"""
prompt=input('''Connect the configuration for measurement 2.
             (y to continue/exits otherwise)''')
if prompt=='y': #prompts to check the input values
    print('Running..')
else:
    sys.exit('Script aborted. Confirm parameter values.')
def phase_mod_sig_in():
    switch_in.set_switch('A',1)
    switch_in.set_switch('B',0)
phase_mod_sig_in()
for inp in in_power_array:
    for freq_no in range(len(vna_center_freq_arr)):
        par_array, freq, logmag_pm, phase = fn.vna_scan(
            vna_center_freq_arr[freq_no], vna_span_arr[freq_no], inp, 
            vna_points, vna_av, 
            el_delay_arr[freq_no], vna_if_bw)
        
            
        filename_str=str(directory_name)+('/pm_line_attenuation')
        filename="%s.dat" %filename_str
        np.savetxt(filename,np.column_stack((freq, logmag_pm, phase)))
            
        fig, ax1=plt.subplots(figsize=(20,10))
        color = 'tab:red'
        ax1.plot(freq*1e-9, logmag_pm, color=color)
        ax1.set_xlabel('Frequency (in GHz)',fontsize = 30)
        ax1.set_ylabel('|S11| (in dB)',fontsize = 30,color=color)
        ax1.tick_params(axis='y', labelcolor=color)
        ax2 = ax1.twinx()
        color = 'tab:blue'
        ax2.set_ylabel('Phase (deg)', color=color,fontsize = 30)
        ax2.plot(freq*1e-9, phase, color=color)
        ax2.tick_params(axis='y', labelcolor=color)
        
        filename="%s.jpeg" %filename_str
        plt.savefig(filename)
        plt.close(fig)            
na.toggle_output(0)
na.set_power(in_power_array[0],ch)

"""Measurement 3"""
prompt=input('''Connect the configuration for measurement 3.
             (y to continue/exits otherwise)''')
if prompt=='y': #prompts to check the input values
    print('Running..')
else:
    sys.exit('Script aborted. Confirm parameter values.')
def sa_sig():
    switch_in.set_switch('A',0)
    switch_in.set_switch('B',1)
    switch_sa.set_switch(1)
sa_sig()
for inp in in_power_array:
    for freq_no in range(len(vna_center_freq_arr)):
        par_array, freq, logmag_sa, phase = fn.vna_scan(
            vna_center_freq_arr[freq_no], vna_span_arr[freq_no], inp, 
            vna_points, vna_av, 
            el_delay_arr[freq_no], vna_if_bw)
        
            
        filename_str=str(directory_name)+('/sa_line_attenuation')
        filename="%s.dat" %filename_str
        np.savetxt(filename,np.column_stack((freq, logmag_sa, phase)))
            
        fig, ax1=plt.subplots(figsize=(20,10))
        color = 'tab:red'
        ax1.plot(freq*1e-9, logmag_sa, color=color)
        ax1.set_xlabel('Frequency (in GHz)',fontsize = 30)
        ax1.set_ylabel('|S11| (in dB)',fontsize = 30,color=color)
        ax1.tick_params(axis='y', labelcolor=color)
        ax2 = ax1.twinx()
        color = 'tab:blue'
        ax2.set_ylabel('Phase (deg)', color=color,fontsize = 30)
        ax2.plot(freq*1e-9, phase, color=color)
        ax2.tick_params(axis='y', labelcolor=color)
        
        filename="%s.jpeg" %filename_str
        plt.savefig(filename)
        plt.close(fig)            
na.toggle_output(0)
na.set_power(in_power_array[0],ch)

"""Interpolate"""
f = interpolate.interp1d(freq, 
                         abs(logmag_pm)-abs(logmag_sa)-abs(logmag_vna),
                         'cubic')
with open(
        directory_name+'/interpolated_vna_sa_attenuation.pkl', 'wb') as k:
    pkl.dump(f, k)
    
with open('interpolated_vna_sa_attenuation.pkl', 'wb') as k:
    pkl.dump(f, k)

"""Get attenuation"""

plt.figure(figsize=(20,10))
plt.scatter(freq, abs(logmag_pm)-abs(logmag_sa)-abs(logmag_vna))
plt.plot(freq, f(freq), color = 'r')
att = np.mean(abs(logmag_pm)-abs(logmag_sa)-abs(logmag_vna))
plt.title('$P_{{SA}}$ = $P_{{VNA}}$ + {} dB'.format(att))    
plt.legend()
plt.xlabel('Frequency (in GHz)',fontsize = 30)
plt.ylabel('|S11| (in dB)',fontsize = 30)
plt.savefig(str(directory_name)+'/attenuation.jpeg')
filename = str(directory_name)+'/attenuation.dat'
np.savetxt(filename, np.column_stack((freq, logmag_vna, logmag_pm,
                                               logmag_sa)))

"""Save results file"""
results_file = directory_name+'/results_attenuation_sa_vna.txt'
file = open(results_file, "a")
file.write('SA detected power = VNA input power + {} dB'.format(att))
file.close()
