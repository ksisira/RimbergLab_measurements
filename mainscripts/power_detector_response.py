# -*- coding: utf-8 -*-
"""
Updated 9/24/2021
Modified 12/22/2021

@author:Sisira

Returns the power detector SNR off-resonance.
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


"""Conditions"""
bias_point='[[0.0,0.0]]'
bias_point_bg='[[0.5,0.0]]'
photons=1
vna_power = -55.4
twpa_bool=True
temp=24

"""Call instruments"""
na=icm.agilent_e5071c(15) #call the network analyzer
gate_meter=icm.hp_34401a(5)
flux_meter=icm.hp_34401a(3)
awg = icm.tektronix_awg520(10)
switch_in = mud.minicircuits_rc4spdta18()
switch_sa = mud.minicircuits_rc1sp4ta18()
daq=icm.national_instruments_bnc2090()
sg=gcm.keysight_n5183b(25)
modgen=scm.Novatech_409B(4)
sa=icm.keysight_n9000b(24)
bp_ctrl = icm.hp_6612c(7)
lockin=gcm.SRS_844(8)
dig = icm.alazartech_ats9462()
if twpa_bool:
    twpa_pump = icm.hp_83711b(13)
else:
    twpa_pump = None
inst_call = [na, sg, modgen, gate_meter, flux_meter, daq, awg, switch_in,
             switch_sa, bp_ctrl, twpa_pump, sa, lockin, dig]
    
"""Make directories"""
current_time=dt.datetime.now()
directory_name=('../{}_{}_{}_{}_power_detector_output_off_resonance_{}_photon'.
                format(current_time.month, current_time.day,current_time.hour,
                       current_time.minute,
                       photons))
os.mkdir(directory_name)
current_file = __file__
fn.copy_file(current_file, directory_name+'/script')

"Parameter file"

date='{}_{}_{}'.format(current_time.month,
      current_time.day,current_time.year)

line = fn.lineno
par_array=[]
line0 = line()
par_array.append([line()-line0,"Date of scan",date])
par_array.append([line()-line0,"Bias point?", bias_point])
par_array.append([line()-line0,"Bias point for bg scan?", bias_point_bg])
par_array.append([line()-line0,"Temperature(mK)",temp])
par_array.append([line()-line0,"Input photons at cavity?",photons])
par_array.append([line()-line0,"TWPA state?",int(twpa_bool)])
par_array.append([line()-line0,"Gate source? (0 for DAQ and 1 for AWG)",0])
par_array.append([line()-line0,"Voltage divider for gate", 101])
par_array.append([line()-line0,"External Resisitor on flux",10e3])
par_array.append([line()-line0,"TWPA frequency?",6.79e9])
par_array.append([line()-line0,"TWPA power?",-3.16])
par_array.append([line()-line0,"Averaging number VNA bias scan (broad)",100])
par_array.append([line()-line0,"VNA Power (dBm)",vna_power])
par_array.append([line()-line0,"VNA smooting aperture",1.5])
par_array.append([line()-line0,"Electrical delay?",63e-9])
par_array.append([line()-line0,"Signal generator carrier power to J0 diff", 16.78])
par_array.append([line()-line0,"Modulation frequency",30e6])
par_array.append([line()-line0,"Modulation voltage (Vpp)",0.181])
par_array.append([line()-line0,"Bandpass control voltage reference (V)",2.44])

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
[date, bias_point, bias_point_bg, temp, photons, twpa_state, 
 gate_source_no, gate_voltage_divider, flux_res,
 twpa_freq, twpa_power,
 vna_av_bias, vna_power_bias, smooth_param, el_delay,
 sg_carrier_J0_dB, mod_freq, mod_volt,
 bp_volt]=parameters

if gate_source_no == 0:
    gate_source = 'daq'
else:
    raise ValueError('No other gate sources available.')

    
"Bandpass function extract"
with open('interpolated_bandpass_freqvsctrlvolt.pkl', 'rb') as handle:
    bp_fn = pkl.load(handle)

bp_offset = bp_volt-bp_fn([5.75e9])[0]
def bp_volt_fn(bp_freq):
    return bp_offset + bp_fn([bp_freq])[0]

"""SA-VNA attenuation"""
with open('interpolated_vna_sa_attenuation.pkl', 'rb') as handle:
    sa_vna_fn = pkl.load(handle)
    
"""--------------------------Measurement Flow-------------------------------"""
"""Reset"""
if twpa_state:
    twpa_pump.toggle_output(0)
    
"Bias points load"
with open('../interpolated_gate_bias_conversion.pkl', 'rb') as handle:
    gate_fn = pkl.load(handle)    
with open('../interpolated_flux_bias_conversion.pkl', 'rb') as handle:
    flux_fn = pkl.load(handle)

bias = bias_point[0]
flux = flux_fn([bias[0]])[0]
fn.bias_source_set_volt('daq', flux, 0, daq = daq)
gate = gate_fn([bias[1]])[0]
fn.bias_source_set_volt(gate_source, gate, 1, 
                        volt_div = gate_voltage_divider, daq = daq)
(res_freq, res_freq_meas, kappa_int, kappa_ext,
 f_logmag, f_phase) = mm.setup_resonance_scan(
        vna_power, 
        vna_av_bias, el_delay,
        bias, bias_point_bg, temp, smooth_param, 
        directory_name, flux_fn, gate_fn,
        gate_source,
        gate_voltage_divider,
        scan_no=2, inst_call=inst_call)
att_mean_vna_sa = sa_vna_fn(res_freq_meas)
sa_J0_power_est = vna_power-att_mean_vna_sa
carrier_power_est = round(sa_J0_power_est+sg_carrier_J0_dB,2)

"""Set up bandpass"""
bp_ctrl.set_voltage(bp_volt_fn(res_freq_meas))
bp_ctrl.toggle_output(1)

"""TWPA parameters"""
if twpa_state:
    twpa_pump.set_frequency(twpa_freq)
    twpa_pump.set_power(twpa_power)
    twpa_pump.toggle_output(1)

carrier_power = carrier_power_est

"""Input signal measurement"""
in_sig_parameters = [carrier_power, res_freq_meas, mod_freq, mod_volt, 0]
peak_power_J0 = mm.powerspectrum_in(directory_name, in_sig_parameters, 
                                    inst_call)

#"""Correct carrier power"""
#sg_carrier_J0_dB_corrected = carrier_power-peak_power_J0
#carrier_power = round(sa_J0_power_est+sg_carrier_J0_dB_corrected,2)
#in_sig_parameters = [carrier_power, res_freq_meas, mod_freq, mod_volt, 1]
#peak_power_J0 = mm.powerspectrum_in(directory_name, in_sig_parameters, 
#                                    inst_call)

"""VNA scan with TWPA on"""
vna_scan_parameters = [bias, temp,
                       vna_power, smooth_param, el_delay, 
                       res_freq_meas]
mm.res_scan_twpa_on(directory_name, vna_scan_parameters, scan_no = 3, 
                    inst_call = inst_call)

"""Response measurement using spectrum analyzer"""
output_parameters = [carrier_power, res_freq_meas, mod_freq, mod_volt,
                     bias, flux_fn, gate_fn, gate_source, gate_voltage_divider,
                     photons]
mm.output_measure(directory_name, output_parameters, inst_call, 
                  '_on_resonance')

"""Move off-resonance"""
flux = flux_fn([bias[0]])[0]
fn.bias_source_set_volt('daq', flux, 0, daq = daq)

"""Repeat output signal measurement"""
output_parameters = [carrier_power, res_freq_meas, mod_freq, mod_volt,
                     bias, flux_fn, gate_fn, gate_source, gate_voltage_divider,
                     photons]
mm.output_measure(directory_name, output_parameters, inst_call,
                  '_off_resonance')

"""Reset"""
bp_ctrl.set_voltage(0)
bp_ctrl.toggle_output(0)
if gate_source_no == 1:
    awg.toggle_output(1,0)
fn.bias_source_set_volt('daq', 0, 0, daq = daq)
fn.bias_source_set_volt(gate_source, 0, 1, daq = daq)
sg.toggle_output(0)
modgen.set_voltage_Vpp(0,0)