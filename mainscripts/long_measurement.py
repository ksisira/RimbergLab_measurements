# -*- coding: utf-8 -*-
"""
Updated 9/24/2021
Modified 12/24/2021

@author:Sisira

Measure fluctuations with and without feedback.

Steps:
    1. Get bias interoplations.
    2. Get input attenuation and fix input power accordingly.
    3. Get resonance scan for the bias point and photon number.
    4. Get the input signal power spectrum.
    5. Get the output signal power spectrum.
    6. Lock-in scan at off-resonance, on-resonance noise
        and on-resonance signal. Measures input ng and Y.
    7. Error signal scan for gate sweep without feedback.
    8. Lock-in scan at off-resonance, on-resonance noise
        and on-resonance signal for lower time constant.
        Measures input ng and Y.
    9. Optimize and output the PID values.
    10. Continue measurements with and without feedback.

Manually set:
    1. Connect input to phase modulator from signal generator, output to
    switch B COM and control from modulation source.
    2. Connect B0 to A1.
    3. Connect B1 to SP4T1.
    4. Connect C1 to bandpass filter in.
    5. Connect BPF out to directional coupler.
    6. Connect directional coupler one output to SP4T 2.
    7. Connect directional coupler second output to DCOM through power detector.
    6. Connect SP4T COM to spectrum analyzer.
    7. Connect fridge in to A COM.
    8. Connect fridge out to C COM.
    9. Connect D0 to lock-in amplifier.
    10. Connect D1 to SP4T3.
    11. Connect VNA P1 to 4SPDTA0.
    12. Connect VNA P2 to 4SPDTC0.
    
    VNA smoothing as required.
    BP filter set up.
    Switch on power supply for room temp amplifier, SDLVA.
     Digitizer measurements with internal clock, set to measure ng and Y of
    lock-in amplifier.
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
from scipy.fft import fft, fftfreq

"""Flags"""
feed_bool = True
twpa_bool = True
bias_scan = True
att_scan = True

"""Conditions"""
bias_point='[[0.0,0.6]]'
bias_point_bg='[[0.5,0.0]]'
temp=24
photons = 10

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
directory_name=('../{}_{}_{}_{}_long_{}_mK_{}_ng_{}_phi0_{}_photons'.
                format(current_time.month, current_time.day,current_time.hour,
                       current_time.minute,temp,
                       eval(bias_point)[0][1], eval(bias_point)[0][0],
                       photons))
os.mkdir(directory_name)
scripts_dir = directory_name+'/scripts'
shutil.copytree(os.getcwd(),scripts_dir)
current_file = __file__
fn.copy_file(current_file, directory_name+'/script')

"Parameter file"

date='{}_{}_{}'.format(current_time.month,
      current_time.day,current_time.year)
gate_check_array = '[0,5e-3,10e-3,-5e-3,-10e-3]'

line = fn.lineno
par_array=[]
line0 = line()
par_array.append([line()-line0,"Date of scan",date])
par_array.append([line()-line0,"Bias point?", bias_point])
par_array.append([line()-line0,"Bias point for bg scan?", bias_point_bg])
par_array.append([line()-line0,"Temperature(mK)",temp])
par_array.append([line()-line0,"Input photons at cavity?",photons])
par_array.append([line()-line0,"Feedback state?", int(feed_bool)])
par_array.append([line()-line0,"TWPA state?",int(twpa_bool)])
par_array.append([line()-line0,"Gate source? (0 for DAQ and 1 for AWG)",0])
par_array.append([line()-line0,"Voltage divider for gate", 101])
par_array.append([line()-line0,"External resistor on flux",10e3])
par_array.append([line()-line0,"TWPA frequency?",6.79e9])
par_array.append([line()-line0,"TWPA power?",-3.16])
par_array.append([line()-line0,"Averaging number VNA bias scan (broad)",100])
par_array.append([line()-line0,"VNA Power for bias scan? (dBm)",-60])
par_array.append([line()-line0,"VNA smooting aperture",1.5])
par_array.append([line()-line0,"Electrical delay?",63e-9])
par_array.append([line()-line0,"Gate check points for bias scan",gate_check_array])
par_array.append([line()-line0,"Max gate voltage for bias scan (V)", 60e-3])
par_array.append([line()-line0,"Number of sweep points for gate for bias scan",41])
par_array.append([line()-line0,"Max flux voltage for bias scan (V)",0.6])
par_array.append([line()-line0,"Number of sweep points for flux for bias scan",41])
par_array.append([line()-line0,"Signal generator carrier power to J0 diff", 16.78])
par_array.append([line()-line0,"Modulation frequency",30e6])
par_array.append([line()-line0,"Modulation voltage (Vpp)",0.118])
par_array.append([line()-line0,"Lockin sensitivity (nV)",100e6])
par_array.append([line()-line0,"Lock-in time constant for sweep (sec)",300e-6])
par_array.append([line()-line0,"Lock-in time constant for feedback (sec)",300e-6])
par_array.append([line()-line0,"Lockin phase in deg",-5.3])
par_array.append([line()-line0,"Lockin expand",1])
par_array.append([line()-line0,"Lockin filter",6])
par_array.append([line()-line0,"Digitizer sampling frequency", 1e6])
par_array.append([line()-line0,"Acquisitions for error signal", 1])
par_array.append([line()-line0,"Acquisition time for setup (sec)", 2])
par_array.append([line()-line0,"Counts for long measurement", 750])
par_array.append([line()-line0,"Acquisition time per buffer for long (sec)", 10])
par_array.append([line()-line0,"Data save sampling frequency setup", 20])
par_array.append([line()-line0,"Data save sampling frequency sweep", 1e3])
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
[date, bias_point, bias_point_bg, temp, photons, feed_bool, twpa_state, 
 gate_source_no, gate_voltage_divider, flux_res,
 twpa_freq, twpa_power,
 vna_av_bias, vna_power_bias, smooth_param, el_delay,
 gate_check_array, gate_max_bias, gate_sweep_no_bias, 
 flux_max_bias, flux_sweep_no_bias, 
 sg_carrier_J0_dB, mod_freq, mod_volt,
 sens, tc_sweep, tc_feed, phase_ref, expand, filt, clock_rate,
 acq_no, acq_time_setup, counts, acq_time_long,
 sample_rate_setup, sample_rate_sweep,
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
    
fn.bias_source_set_volt(gate_source, 0, 1, 
                        daq = daq) 
    
"""Bias Scan"""
if bias_scan:
    bias_scan_parameters = [temp,
                            gate_source_no, gate_voltage_divider,
                            vna_av_bias, vna_power_bias, 
                            smooth_param, el_delay,
                            gate_check_array, gate_max_bias, 
                            gate_sweep_no_bias,
                            flux_max_bias, flux_sweep_no_bias]
    mm.get_bias_position(directory_name, bias_scan_parameters, inst_call)
else:
    prompt_cont = input(('Copy flux, gate and bg interpolation files.'+
                         'Continue? (y/exits)'))
"Bias points load"
with open(directory_name +
          '/interpolated_gate_bias_conversion.pkl', 'rb') as handle:
    gate_fn = pkl.load(handle)    
with open(directory_name +
          '/interpolated_flux_bias_conversion.pkl', 'rb') as handle:
    flux_fn = pkl.load(handle)
with open('interpolated_kerr_values.pkl', 'rb') as handle: #fn of ng,phi
    kerr_fn = pkl.load(handle)
with open('interpolated_res_freq_values.pkl', 'rb') as handle:
        # fn of ng, phi
        res_fn = pkl.load(handle)

if gate_source_no == 1:
    awg.toggle_output(1,1)

"""Input attenuation and photon number in the cavity"""
bias = bias_point[0]    
kerr_value=kerr_fn(bias[1],bias[0])[0] #Hz

if abs(kerr_value*1e-6) < 0.1:
    prompt_bias = [0.2,0] #This is set assuming bias points is near (0.25,0)
    flux_replace = flux_fn([prompt_bias[0]])[0]
    fn.bias_source_set_volt('daq', flux_replace, 0, daq = daq)
    gate_replace = gate_fn([prompt_bias[1]])[0]
    kerr_value=kerr_fn(prompt_bias[1],prompt_bias[0])[0] #Hz
    fn.bias_source_set_volt(gate_source, gate_replace, 1, 
                            volt_div = gate_voltage_divider, daq = daq)
    time.sleep(2)
    (res_freq, res_freq_meas,
     kappa_int, kappa_ext,
     f_logmag, f_phase)= mm.setup_resonance_scan(vna_power_bias, 
                         vna_av_bias, el_delay,
                         prompt_bias, bias_point_bg, temp, smooth_param, 
                         directory_name, flux_fn, gate_fn, gate_source,
                         gate_voltage_divider,
                         scan_no=11, inst_call=inst_call)
    kerr_bias_point = [prompt_bias]
else:
    flux = flux_fn([bias[0]])[0]
    fn.bias_source_set_volt('daq', flux, 0, daq = daq) 
    gate = gate_fn([bias[1]])[0]
    fn.bias_source_set_volt(gate_source, gate, 1, 
                            volt_div = gate_voltage_divider, daq = daq)
    time.sleep(2.0)
    (res_freq, res_freq_meas,
     kappa_int, kappa_ext,
     f_logmag, f_phase) = mm.setup_resonance_scan(vna_power_bias, 
                         vna_av_bias, el_delay,
                         bias, bias_point_bg, temp, smooth_param, 
                         directory_name, flux_fn, gate_fn,
                         gate_source,
                         gate_voltage_divider,
                         scan_no=1, inst_call=inst_call)
    kerr_bias_point = bias_point

"""Estimate the photon number for 1.5 MHz shift"""
est_photon_no = abs(1.5e6/kerr_value)

"""At (0,0), 3.5 photons correspond to 11 nW roughly. Get the estimate
for other bias points"""
power_max = est_photon_no*11/3.5
if power_max>20:
    power_max=20
    
"""Fix this as the max point for power sweep"""
power_array_nW = np.round(np.linspace(0.3, power_max, 15),3)
  
input_att_parameters = [kerr_bias_point, temp, smooth_param, el_delay, 
                        power_array_nW, res_freq, kerr_value,
                        kappa_int, kappa_ext, f_logmag, f_phase]
if att_scan:
    input_att_dB, input_att, power_array, n_array = mm.input_attenuation_scan(
            directory_name,                                                            
            input_att_parameters, inst_call)
    
    kerr_value=kerr_fn(bias[1],bias[0])[0] #Hz
    if abs(kerr_value*1e-6) < 0.1:
        flux = flux_fn([bias[0]])[0]
        fn.bias_source_set_volt('daq', flux, 0, daq = daq)
        gate = gate_fn([bias[1]])[0]
        fn.bias_source_set_volt(gate_source, gate, 1, 
                                volt_div = gate_voltage_divider, daq = daq)
        time.sleep(2.0)
        (res_freq, res_freq_meas,
         kappa_int, kappa_ext,
         f_logmag, f_phase) = mm.setup_resonance_scan(vna_power_bias, 
                             vna_av_bias, el_delay,
                             bias, bias_point_bg, temp, smooth_param, 
                             directory_name, flux_fn, gate_fn,
                             gate_source,
                             gate_voltage_divider,
                             scan_no=1, inst_call=inst_call)
        h = 6.626e-34
        kappa_tot = kappa_int + kappa_ext
        n_array = 4*kappa_ext*power_array_nW*1e-9/(h*res_freq_meas*
                                                   kappa_tot**2*
                                                   np.pi*2)/abs(input_att)
        filename_resonance=directory_name+'/no_of_photons.dat'
        np.savetxt(filename_resonance, np.column_stack((power_array, 
                                                        power_array_nW,
                                                        n_array)))
        
        fig, ax = plt.subplots(figsize=(20,10))
        ax.spines['bottom'].set_linewidth(5)
        ax.spines['top'].set_linewidth(5)
        ax.spines['left'].set_linewidth(5)
        ax.spines['right'].set_linewidth(5)
        ax.tick_params(axis='x',direction='out', labelsize=40,width=2,length=10)
        ax.tick_params(axis='y',direction='out', labelsize=40,width=2,length=10)
        plt.scatter(power_array, n_array, s = 60)
        plt.xlabel('Input power (dBm)',fontsize = 30)
        plt.ylabel('Average photons',fontsize = 30)
        plt.yscale('log')
        plt.plot(power_array, n_array, 'r', lw = 4)
        filename_plot = directory_name+'/input_power_to_photons_at_res.jpeg'
        plt.title(('Input attenuation = {} dB,\n'.format(round(input_att_dB,2))+
                   '$\Phi_{{ext}}$ = {} $\Phi_0$ and $n_g$ = {} e\n'.format(
                           bias_point[0][0], bias_point[0][1])+
                           'Kerr coefficient = {} MHz'.format(
                                   round(kerr_value*1e-6,3))), 
                           fontsize=30, 
        y=0.87)
        plt.tight_layout()
        plt.savefig(filename_plot)
else:
    prompt_cont = input(('Copy input attenuation file.'+
                         'Continue? (y/exits)'))
    input_att_dB = input('Input attenuation in dB?')
    filename = directory_name+'/no_of_photons.dat'
    att_file = np.loadtxt(filename)
    n_array = att_file[:,2]
    power_array_nW = att_file[:,1]
slope_kerr, intercept_kerr, r_value, p_value, std_err = stats.linregress(
            n_array,power_array_nW)

"""VNA scan for required photon number"""
vna_power_dBm = 10*np.log10((slope_kerr*photons+intercept_kerr)*1e-6)
(res_freq, res_freq_meas, kappa_int, kappa_ext,
 f_logmag, f_phase) = mm.setup_resonance_scan(
        vna_power_dBm, 
        vna_av_bias, el_delay,
        bias, bias_point_bg, temp, smooth_param, 
        directory_name, flux_fn, gate_fn,
        gate_source,
        gate_voltage_divider,
        scan_no=2, inst_call=inst_call)
att_mean_vna_sa = sa_vna_fn(res_freq_meas)
sa_J0_power_est = vna_power_dBm-att_mean_vna_sa
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
vna_scan_parameters = [bias_point, temp,
                       vna_power_dBm, smooth_param, el_delay, 
                       res_freq_meas,
                       bias_point_bg, gate_source, gate_voltage_divider,
                       flux_fn, gate_fn]
mm.res_scan_twpa_on(directory_name, vna_scan_parameters, scan_no = 3, 
                    inst_call = inst_call)

"""Response measurement using spectrum analyzer"""
output_parameters = [carrier_power, res_freq_meas, mod_freq, mod_volt,
                     bias, flux_fn, gate_fn, gate_source, gate_voltage_divider,
                     photons]
mm.output_measure(directory_name, output_parameters, inst_call)

"""Error signal Setup"""
ng_feed_corrected_array=[]
ng_dc_feed_mean = bias[1]
count=0
if feed_bool:
    feed_prompt = input('''Continue to feedback measurement? 
                        Switch on the 1:11 voltage divider. (y/n)''')
    while feed_prompt=='y': #prompts to check the feedback loop to continue.
        count = count+1
        print('Continuing measurements..')
#        bias[1] = ng_dc_feed_mean
        PID_prompt = input(' PID optimized values. Enter [P,I,D,UL,LL]')
        fft_parameters = [photons, bias, temp, mod_volt,
                           sens, tc_feed, phase_ref, expand, filt, 
                           clock_rate, sample_rate_setup, acq_no, 
                           acq_time_setup, res_freq_meas,
                           carrier_power, mod_freq,
                           gate_source_no, gate_voltage_divider, 
                           PID_prompt]
        ng_dc_feed_mean = mm.feedback_on_off(directory_name, 
                                             fft_parameters, inst_call, 
                                             feed_str = '_feedback_on_{}'.
                                             format(count))
        ng_feed_corrected_array.append(ng_dc_feed_mean)
        reset_prompt = input('''Reset feedback. Switch off 1:11 and
                             switch on 1:101 (y)''')
        if reset_prompt!='y':
            print('NOTE!!!!!! FEEDBACK NOT RESET!!!!!!!!!.')
        feed_prompt = input('''Continue to feedback measurement? (y/n)''')
    else:
        count = count+1
#        bias[1] = ng_dc_feed_mean
        fft_parameters = [photons, bias, temp, mod_volt,
                           sens, tc_feed, phase_ref, expand, filt, 
                           clock_rate, sample_rate_setup, acq_no, 
                           acq_time_long, res_freq_meas,
                           carrier_power, mod_freq,
                           gate_source_no, gate_voltage_divider, 
                           PID_prompt,counts]
        mm.feedback_long(directory_name, fft_parameters, inst_call, 
                         feed_str = '_feedback_on_{}'.format(count))
        reset_prompt = input('''Reset feedback. Switch off 1:11 and
                             switch on 1:101 (y)''')
        if reset_prompt!='y':
            print('NOTE!!!!!! FEEDBACK NOT RESET!!!!!!!!!.')
        print('Feedback measurements completed.')
else:
    print('Script written for feedback configuration.')

if twpa_state:
    twpa_pump.toggle_output(0)
    
ng_feed_corrected_array = np.asarray(ng_feed_corrected_array)

"""Get new bias point"""
if feed_bool:
    filename=(directory_name+
              '/Lockin_fft_feedback_on_1/quadrature_signal_feed_dc')
    ng_initial = np.loadtxt(filename+'.dat', encoding = 'utf-16')
    ng_in_mean = np.mean(ng_initial[:,1])
    filename=(directory_name+
              '/Lockin_fft_feedback_on_{}/quadrature_signal_feed_dc_{}'.format(
                      count,counts-1))
    cor_data = np.loadtxt(filename+'.dat', encoding = 'utf-16')
    ng_final_mean = np.mean(cor_data[:,1])
    new_ng = round(bias[1]/ng_in_mean*ng_final_mean,3)
    
"""Check if the resonance moved"""  
fn.vna_connect(switch_in)
time.sleep(1)
(res_freq, res_freq_meas, kappa_int, kappa_ext,
 f_logmag, f_phase) = mm.setup_resonance_scan(
        vna_power_dBm, 
        vna_av_bias, el_delay,
        [bias[0],new_ng], bias_point_bg, temp, smooth_param, 
        directory_name, flux_fn, gate_fn,
        gate_source,
        gate_voltage_divider,
        scan_no=4, inst_call=inst_call)

"""Check if the resonance moved"""  
fn.vna_connect(switch_in)
time.sleep(1)
(res_freq, res_freq_meas, kappa_int, kappa_ext,
 f_logmag, f_phase) = mm.setup_resonance_scan(
        vna_power_dBm, 
        vna_av_bias, el_delay,
        bias, bias_point_bg, temp, smooth_param, 
        directory_name, flux_fn, gate_fn,
        gate_source,
        gate_voltage_divider,
        scan_no=5, inst_call=inst_call)

"""Reset"""
bp_ctrl.set_voltage(0)
bp_ctrl.toggle_output(0)
if gate_source_no == 1:
    awg.toggle_output(1,0)
fn.bias_source_set_volt('daq', 0, 0, daq = daq)
fn.bias_source_set_volt(gate_source, 0, 1, daq = daq)
fn.error_sig_measure(switch_in, switch_sa)
fn.vna_connect(switch_in)

"""Save results file"""
results_file = directory_name+'/results_feedback_fft.txt'
file = open(results_file, "a")
file.write('Kerr value for phi, ng = {} is {} MHz.\n'.format(
        bias,kerr_value*1e-6))
file.write('Input attenuation in dB is {}.\n'.format(input_att_dB))
file.write('VNA power for {} photon is {} dBm.\n'.format(
        photons,vna_power_dBm))
file.write('bg feedback corrected to following values.{}\n'.format(
        ng_feed_corrected_array))
file.close()

"""Plot long measurement"""
ng_measured=[]
for i in range(counts):
    filename=(directory_name+
          '/Lockin_fft_feedback_on_{}/quadrature_signal_feed_dc_{}'.format(
                  count,i))
    data_load = np.loadtxt(filename+'.dat', encoding = 'utf-16')
    ng_measured.append(data_load[:,1])
ng_measured = (np.asarray(ng_measured)).flatten()
plt.figure()
plt.plot(ng_measured[100:])
plt.savefig(directory_name+'/long_ng.jpeg')

Y_measured_off=[]
for i in range(counts):
    filename=(directory_name+
          '/Lockin_fft_feedback_on_{}/quadrature_signal_feed_dc_{}'.format(
                  count,i))
    data_load = np.loadtxt(filename+'.dat', encoding = 'utf-16')
    Y_load = data_load[:,2]
    Y_measured_off.append(Y_load)
Y_measured_off = (np.asarray(Y_measured_off)).flatten()
plt.figure()
l=len(Y_measured_off)
Y_fft_off = fft(np.abs(Y_measured_off[200:]))
plt.loglog(Y_fft_off)
plt.savefig(directory_name+'long_Y_fft.jpeg')
#
#Y_measured_on=[]
#for i in range(counts):
#    filename=('../12_28_12_8_feedback_fft_24_mK_0.6_ng_0.0_phi0_5_photons'+
#          '/Lockin_fft_feedback_on_{}/quadrature_signal_feed_dc_{}'.format(
#                  count,i))
#    data_load = np.loadtxt(filename+'.dat', encoding = 'utf-16')
#    Y_load = data_load[:,2]
#    Y_measured_on.append(Y_load)
#Y_measured_on = (np.asarray(Y_measured_on)).flatten()
#l=len(Y_measured_on)
#Y_fft_on = fft(Y_measured_on[:l//2])
#plt.loglog(Y_fft_on)
#plt.ylim(1,1e2)