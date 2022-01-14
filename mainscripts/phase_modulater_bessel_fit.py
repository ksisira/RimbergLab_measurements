# -*- coding: utf-8 -*-
"""
Created on Sat Feb  8 18:58:07 2020.
Modified Dec 20th, 2021.

@author: Sisira
Get power spectrum of phase modulated signal and bessel fit.

Connections:
    1. Signal generator to Phase modulator(PM) HMC-C010 Input.
    2. Modulation generator to PM control.
    3. PM out to spectrum analyzer.
"""

"""Import modules"""
import functions as fn
import instrument_classes_module as icm
import serial_classes_module as scm
import gpib_classes_module as gcm
import numpy as np #Scientific computing module
import os, sys #system functions
import time #for transient delay
import subprocess as sp #opening files in notepad
import matplotlib.pyplot as plt #Plotting module
import datetime as dt #date module
from scipy import special, stats

"""Conditions"""
mod_freq = '[30e6]'
carrier_power = '[-30, -40]'
carrier_freq = '[5.75e9]'

"""Call instruments"""
sa = icm.keysight_n9000b(24)
modgen = scm.Novatech_409B(4)
sg=gcm.keysight_n5183b(25)

"""Make directories"""
current_time=dt.datetime.now()
directory_name='{}_{}_{}_{}_power_mod_calibrate'.format(current_time.month,
                current_time.day,current_time.hour,current_time.minute)
os.mkdir(directory_name)
current_file = __file__
fn.copy_file(current_file, directory_name+'/script')

"""Parameter file"""
date='{}_{}_{}'.format(current_time.month,
      current_time.day,current_time.year)
line = fn.lineno
par_array=[]
line0 = line()
par_array=[]
par_array.append([line()-line0,"Date of scan",date])
par_array.append([line()-line0,"Carrier frequency (GHz)",carrier_freq])
par_array.append([line()-line0,"Carrier power (in dBm)",carrier_power])
par_array.append([line()-line0,"Modulation frequency (MHz)",mod_freq])
par_array.append([line()-line0,"Modulation start voltage",0])
par_array.append([line()-line0,"Spectrum analyzer central frequency",carrier_freq])
par_array.append([line()-line0,"Spectrum analyzer span",5*eval(mod_freq)[0]])
par_array.append([line()-line0,"Resolution bandwidth",10e3])
par_array.append([line()-line0,"video bandwidth",10e3])
par_array.append([line()-line0,"averaging",50])
par_array.append([line()-line0,"Peak threshold", -90])
par_array.append([line()-line0,"Number of modulation sweep points", 50])
par_array.append([line()-line0,"Modulation amplitude sweep step (Vpp)", 10e-3])
par_array.append([line()-line0,"Reference level (dBm)", 0])

filename_parameter_in=str(directory_name)+'/Parameters.dat'
np.savetxt(filename_parameter_in,np.row_stack((par_array)),
           delimiter=';',fmt='%s')
programName = "notepad.exe"
sp.Popen([programName,filename_parameter_in]) #opens file in notepad
prompt=input('''Check the notepad file. Continue with the saved parameter\
             values? (y to continue/exits otherwise)''')
if prompt=='y': #prompts to check the input values
    print('Running..')
else:
    sys.exit('Script aborted. Confirm parameter values')
    
parameters=np.loadtxt(filename_parameter_in,delimiter=';',usecols=2,
                      dtype = np.str) #loads the parameter values

"""Parameter values"""
parameters = np.array([eval(i) for i in parameters])
[date,carrier_freq,carrier_power,mod_freq,mod_start_vpp,sa_central_freq,
 sa_span,res_bw,video_bw,avg_no,peak_threshold,mod_points,mod_step,
 ref_level]=parameters
 
avg_no =int(avg_no)
mod_points=int(mod_points)

"""Setting the parameters"""
modgen.set_state_freq(0,1,mod_freq[0]*1e-6,mod_freq[0]*1e-6)
modgen.set_state_Vpp(0,1,0,0)
modgen.set_state_phase(0,1, 0.0, 0.0)
sg.set_power(carrier_power[0])
sa.set_timeout(120e3)
sa.set_resolution_bandwidth(res_bw)
sa.set_video_bandwidth(video_bw)
sa.toggle_continuous_sweep('off')
sa.set_averaging(avg_no)
sg.toggle_output(1)

output_peak_freq=[]
output_peak_power=[]
mod_vpp_array=[]
output_power=[]
variables=[]

"""----------------------------Measurements---------------------------------"""
for ghz in range(len(carrier_freq)):
    sg.set_frequency(carrier_freq[ghz])
    sa.set_frequency_center(sa_central_freq[ghz])
    
    for mhz in range(len(mod_freq)):
        modgen.set_state_freq(0,1,mod_freq[mhz]*1e-6,mod_freq[mhz]*1e-6)
        sa.set_frequency_span(sa_span)
        
        for power in carrier_power:
            sg.set_power(power)
            
            for i in range(0,int(mod_points)):
                mod_vpp = mod_start_vpp + i*mod_step
                mod_vpp_array.append(mod_vpp)
                modgen.set_state_Vpp(0,1,mod_vpp,mod_vpp)
                time.sleep(1)
                filename_data=(directory_name+
                               '/power_spectrum_{}GHz_{}MHz_{}dbm_{}Vpp.dat'.
                               format(carrier_freq[ghz]*1e-9,
                                      mod_freq[mhz]*1e-6,power,mod_vpp))
                (par_array, freq, power_out) = fn.spec_scan(
                        sa_central_freq[ghz], sa_span,
                        res_bw, video_bw, avg_no, ref_level, filename_data)
                peak_freq,peak_power=sa.get_peaks(peak_threshold)
                peak_freq=np.asarray(peak_freq)
                peak_freq_sort=[]
                peak_power_sort=[]
                for j in range(-2,3):
                    pos = np.where(peak_freq==carrier_freq[ghz]+
                                   j*mod_freq[mhz])[0]
                    if len(pos)==0:
                        peak_freq_sort.append(np.NaN)
                        peak_power_sort.append(np.NaN)
                    else:
                        peak_power_sort.append(peak_power[pos[0]])
                        peak_freq_sort.append(peak_freq[pos[0]])
                output_peak_freq.append(peak_freq_sort)
                output_peak_power.append(peak_power_sort)
                variables.append([carrier_freq[ghz],mod_freq[mhz],
                                  power,mod_vpp])
                print(i,carrier_freq[ghz]*1e-9,
                      mod_freq[mhz]*1e-6,power,mod_vpp,peak_power_sort)    

modgen.set_state_Vpp(0,1,0,0)
sg.toggle_output(0)

"""Plots and Data"""
mod_vpp_array=np.asarray(mod_vpp_array)
output_peak_freq=np.asarray(output_peak_freq)
output_peak_power=np.asarray(output_peak_power)
output_peak_amp_square = 10**(output_peak_power/10)

l = len(output_peak_freq)
for j in range(0,l,mod_points):
    plt.figure(figsize=(20,10))
    ax=plt.subplot(111)
    for i in range(5):
        deltaf = round((
                output_peak_freq[l-1,i] - output_peak_freq[l-1,2])*1e-6,2)
        plt.scatter(mod_vpp_array[j:j+mod_points],
                    output_peak_power[j:j+mod_points,i], 
                    label = '''$\Delta$f = {} MHz'''.format(deltaf))
        plt.plot(mod_vpp_array[j:j+mod_points],
                 output_peak_power[j:j+mod_points,i])

    plt.legend()
    plt.xlabel('Modulation Vpp (V)')
    plt.ylabel('Sideband power (dbm)')
    filename_plot=(directory_name+
                   '/modulation_logpower_dependence_{}GHz_{}MHz_{}dbm.jpeg'.
                   format(variables[j][0]*1e-9,
                          variables[j][1]*1e-6,variables[j][2]))
    plt.savefig(filename_plot)
    filename = (directory_name+
                '/peak_data_{}GHz_{}MHz_{}dbm.dat'.format(variables[j][0]*1e-9,
                            variables[j][1]*1e-6,variables[j][2]))
    np.savetxt(filename,np.column_stack((mod_vpp_array[j:j+mod_points],
                                         output_peak_power[j:j+mod_points])))

for j in range(0,l,mod_points):
    plt.figure(figsize=(20,10))
    ax=plt.subplot(111)
    norm_factor = output_peak_amp_square[j,2]
    for i in range(5):
        deltaf = round((
                output_peak_freq[l-1,i] - output_peak_freq[l-1,2])*1e-6,2)
        plt.scatter(mod_vpp_array[j:j+mod_points],
                    output_peak_amp_square[j:j+mod_points,i]/norm_factor, 
                    label = '''$\Delta$f = {} MHz'''.format(deltaf))
        plt.plot(mod_vpp_array[j:j+mod_points],
                 output_peak_amp_square[j:j+mod_points,i]/norm_factor)

    plt.legend()
    plt.xlabel('Modulation Vpp (V)')
    plt.ylabel('Sideband amplitude square')
    filename_plot=(directory_name+
                   '/modulation_ampsquare_dependence_{}GHz_{}MHz_{}dbm.jpeg'.
                   format(variables[j][0]*1e-9,
                          variables[j][1]*1e-6,variables[j][2]))
    plt.savefig(filename_plot)


beta = np.linspace(0,2,1000)
zeroth = special.jv(0,beta)**2
first = special.jv(1,beta)**2
second =special.jv(2,beta)**2
diff = abs(zeroth-first)
pos = np.where(diff == np.min(diff))[0][0]
intersect = beta[pos]

fit_parameters=[]
for j in range(0,l,mod_points):
    norm_factor = output_peak_amp_square[j,2]
    data_zeroth = output_peak_amp_square[j:j+mod_points,2]/norm_factor
    data_first = (output_peak_amp_square[j:j+mod_points,1]+
                  output_peak_amp_square[j:j+mod_points,3])/2/norm_factor
    data_diff = abs(data_zeroth-data_first)
    pos_data = np.where(data_diff == np.nanmin(
            data_diff[:int(mod_points/2)]))[0]
    if len(pos_data) == 0:
        fit_parameters.append([np.NaN,np.NaN])
        continue
    mod_vpp_subset = mod_vpp_array[j:j+mod_points]
    slope, intercept, r_value, p_value, std_err = stats.linregress(
            [0,mod_vpp_subset[pos_data[0]]],[0,intersect])
    fit_parameters.append([slope,intercept])
    slope_rev, intercept_rev, r_value, p_value, std_err = stats.linregress(
            [0,intersect],[0,mod_vpp_subset[pos_data[0]]])
    beta_optimal = 1.08
    optimal_voltage = slope_rev*beta_optimal+intercept_rev
    beta_optimal_2 = 1.84
    optimal_voltage_2 = slope_rev*beta_optimal+intercept_rev

fit_parameters=np.asarray(fit_parameters)
for k in range(3):
    plt.figure(figsize=(20,10))
    ax=plt.subplot(111)
    for j in range(len(fit_parameters)):
        modulation_index = fit_parameters[j,0]*mod_vpp_array[
                j*mod_points:(j+1)*mod_points]+fit_parameters[j,1]
        beta = np.linspace(0,np.max(modulation_index),1000)
        fit = special.jv(2-k,beta)**2
        plt.plot(modulation_index,output_peak_amp_square[
                j*mod_points:(j+1)*mod_points,k]/output_peak_amp_square[
                        j*mod_points,2], 
        label ='{}GHz, {} MHz, {} dBm'.format(variables[j*mod_points][0]*1e-9,
                variables[j*mod_points][1]*1e-6,variables[j*mod_points][2]))
        plt.scatter(modulation_index, output_peak_amp_square[
                j*mod_points:(j+1)*mod_points,
                k]/output_peak_amp_square[j*mod_points,2])
        if j==len(fit_parameters)-1:
            np.savetxt(directory_name+'/J_{}_data.dat'.format(abs(k-2)),
                       np.column_stack((modulation_index,
                                        output_peak_amp_square[
                                                j*mod_points:(j+1)*mod_points,
                                                k]/output_peak_amp_square[
                                                        j*mod_points,2])))
        if k !=2:
            plt.plot(modulation_index,output_peak_amp_square[
                    j*mod_points:(j+1)*mod_points,
                    -(k+1)]/output_peak_amp_square[j*mod_points,2])
            plt.scatter(modulation_index,
                        output_peak_amp_square[j*mod_points:(j+1)*mod_points,
                                               -(k+1)]/output_peak_amp_square[
                                                       j*mod_points,2])
            np.savetxt(directory_name+'/J_{}_data.dat'.format(k-2),
                       np.column_stack((modulation_index,
                                        output_peak_amp_square[
                                                j*mod_points:(j+1)*mod_points,
                                                -(k+1)]/output_peak_amp_square[
                                                        j*mod_points,2])))
        
    plt.plot(beta,fit,label='J_{}'.format(2-k))
    plt.legend()
    plt.xlabel('Modulation index')
    plt.ylabel('Normalized amplitude square')
    filename_plot=directory_name+'/bessel_fit_{}th_band.jpeg'.format(2-k)
    plt.savefig(filename_plot)


plt.figure(figsize=(20,10))
ax=plt.subplot(111)
for j in range(len(fit_parameters)):
    modulation_index = fit_parameters[j,0]*mod_vpp_array[
            j*mod_points:(j+1)*mod_points]+fit_parameters[j,1]
    plt.plot(modulation_index,
             mod_vpp_array[j*mod_points:(j+1)*mod_points],
             label ='{}GHz, {} MHz, {} dBm'.format(variables[
                     j*mod_points][0]*1e-9,variables[j*mod_points][1]*1e-6,
            variables[j*mod_points][2]))
    plt.scatter(modulation_index,mod_vpp_array[j*mod_points:(j+1)*mod_points])
plt.legend()
plt.xlabel('Modulation index')
plt.ylabel('Modulation Vpp (V)')
filename_plot=directory_name+'/modulation_index_vs_vpp.jpeg'
plt.savefig(filename_plot)

"""Save results file"""
results_file = directory_name+'/results_final_setting.txt'
file = open(results_file, "a")
file.write("Intersection value (J0=J1) of beta is {}\n".format(intersect))
file.write('''Slope and intercept fit from modulation voltage to beta are\
           {} and {}\n'''.format(slope,intercept))
file.write('''Slope and intercept fit from beta to modulation voltage are\
           {} and {}\n'''.format(slope_rev, intercept_rev))
file.write('Modulation voltage corresponding to beta=1.08 is {} V'.
           format(optimal_voltage))
file.write('Modulation voltage corresponding to beta=1.842 is {} V'.
           format(optimal_voltage_2))
file.close()

np.savetxt(directory_name+'/beta_to_mod.dat',
           np.column_stack((beta,slope_rev*beta+intercept_rev)))