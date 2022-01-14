# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 17:31:07 2020

@author: Sisira

Module containing some useful functions for running experiments.
"""
import instrument_classes_module as icm
import serial_classes_module as scm
import matplotlib.pyplot as plt
import minicircuits_usb_devices as mud
import gpib_classes_module as gcm
import usb_classes_module as ucm
import pickle as pkl
import numpy as np
import datetime as dt
import os, sys
import shutil
import inspect
import time
from scipy import interpolate, stats
from lmfit.models import QuadraticModel

mod_gen_address = 4

"""----------------------------------SCANS----------------------------------"""
def vna_scan(center_freq, span, power, sweep_points, av_no, 
             el_delay, IF_bw, address = 15, channel = 1):
    """Takes a scan on the VNA and returns parameter array, and data"""
    na=icm.agilent_e5071c(address)
    par_array={'channel_number' : channel,
               'number_of_traces' : 2,
               'center_freq_Hz' : center_freq,
               'span_Hz' : span,
               'power_dbm' : power,
               'sweep_points' : sweep_points,
               'average_no' : av_no,
               'electrical_delay_sec' : el_delay,
               'IF_bandwidth_Hz' : IF_bw}
    na.set_timeout(60e3)
    na.allocate_traces('1_2')
    na.set_channel(channel)
    na.set_num_traces(2)
    na.set_measurement('S21',channel)
    na.set_sweep_type('lin',channel)
    na.set_sweep_mode('STEP')
    na.set_format('MLOG',channel,1)
    na.set_format('PPH',channel,2)
    na.toggle_averaging(1,channel)
    na.set_trigger_source('man')
    na.toggle_averaging_trigger(1)
    na.toggle_continuous_triggering(1,channel)
    na.set_sweep_points(sweep_points,channel)
    na.set_averaging(av_no,channel)
    na.set_electrical_delay(el_delay,channel,2)
    na.set_if_bandwidth(IF_bw,channel)
    na.set_freqs(center_freq,span,'span',channel)
    na.set_power(power,channel)
    na.toggle_output(1)
    na.trigger(1)
    
    freq=na.get_frequency_data(channel)
    freq=np.fromstring(freq,sep=',')
    for trace in range(1,3):
        na.autoscale(channel,trace)
        na.transfer_data_to_memory(channel,trace)
        tracedata=na.get_trace_data(channel,trace)
        tracedata=np.fromstring(tracedata,sep=',')
        if trace == 1:
            logmag = tracedata
        else:
            phase = tracedata
    
    na.toggle_output(0)
    
    return (par_array, freq, logmag, phase)

def vna_scan_save(center_freq, span, power, sweep_points, av_no, 
               el_delay, IF_bw, 
               filename_str_logmag = None,filename_str_phase = None, 
               address = 15, channel = 1, bg = False,
               f_mag=None, f_phase=None, save = False, plot = False):
    par_array, freq, logmag, phase = vna_scan(
            center_freq, span, power, sweep_points, av_no, 
            el_delay, IF_bw)
    """Takes a scan on the VNA and saves data and figures after subtracting
    background depending on the bool parameters passed."""
    if bg:
        bg_reduced_data_mag = logmag - f_mag(freq)
        bg_reduced_data_mag = np.asarray(bg_reduced_data_mag)
        bg_reduced_data_phase = phase - f_phase(freq)
        bg_reduced_data_phase = np.asarray(bg_reduced_data_phase)
        if save:
            filename_logmag=filename_str_logmag+'.dat'
            np.savetxt(filename_logmag,np.column_stack((freq,logmag,
                                                        bg_reduced_data_mag)))
            filename_phase=filename_str_phase+'.dat'
            np.savetxt(filename_phase,np.column_stack((freq,phase,
                                                       bg_reduced_data_phase)))
        res_amp_pos = np.where(
                bg_reduced_data_mag == min(bg_reduced_data_mag))[0]
        res_amp = bg_reduced_data_mag[res_amp_pos][0]
        center_freq_ref = freq[res_amp_pos][0]
        if plot:
            plt.figure()
            plt.plot((center_freq_ref*1e-9),(res_amp),'ro')
            plt.plot(freq*1e-9,bg_reduced_data_mag)
            plt.xlabel('Frequency (GHz)',fontsize = 20)
            plt.ylabel('|S11|-bg (dB)',fontsize = 20)
            plt.title('Resonance at {} GHz'.format(round(
                    center_freq_ref*1e-9,5)),
                  fontsize=20)
            filename_logmag=filename_str_logmag+'.jpeg'
            plt.savefig(filename_logmag)
            plt.close()
        
            plt.figure()
            color = 'tab:red'
            plt.plot(freq*1e-9,bg_reduced_data_phase,color=color)
            plt.xlabel('Frequency (GHz)',fontsize = 20)
            plt.ylabel('Phase-bg (deg)',color=color,fontsize = 20)
            plt.title('Resonance at {} GHz'.format(round(
                    center_freq_ref*1e-9,5)),
                      fontsize=20)
            filename_phase=filename_str_phase+'.jpeg'
            plt.savefig(filename_phase)
            plt.close()
        
        return (par_array, freq, logmag, phase,
                res_amp_pos, res_amp, center_freq_ref, 
                bg_reduced_data_mag, bg_reduced_data_phase)
    else:
        if save:
            filename_logmag=filename_str_logmag+'.dat'
            np.savetxt(filename_logmag,np.column_stack((freq,logmag)))
            filename_phase=filename_str_phase+'.dat'
            np.savetxt(filename_phase,np.column_stack((freq,phase)))
        res_amp_pos = np.where(logmag == min(logmag))[0]
        res_amp = logmag[res_amp_pos][0]
        center_freq_ref = freq[res_amp_pos][0]
        if plot:
            plt.figure()
            plt.plot((center_freq_ref*1e-9),(res_amp),'ro')
            plt.plot(freq*1e-9,logmag)
            plt.xlabel('Frequency (GHz)',fontsize = 20)
            plt.ylabel('|S11|(dB)',fontsize = 20)
            plt.title('Resonance at {} GHz'.format(round(
                    center_freq_ref*1e-9,5)),
                  fontsize=20)
            filename_logmag=filename_str_logmag+'.jpeg'
            plt.savefig(filename_logmag)
            plt.close()
        
            plt.figure()
            color = 'tab:red'
            plt.plot(freq*1e-9,phase,color=color)
            plt.xlabel('Frequency (GHz)',fontsize = 20)
            plt.ylabel('Phase (deg)',color=color,fontsize = 20)
            plt.title('Resonance at {} GHz'.format(round(
                    center_freq_ref*1e-9,5)),
                      fontsize=20)
            filename_phase=filename_str_phase+'.jpeg'
            plt.savefig(filename_phase)
            plt.close()
        return (par_array, freq, logmag, phase,
                res_amp_pos, res_amp, center_freq_ref)

def spec_scan(center_freq, span, res_bw, video_bw, avg_no, ref_level, 
              filename = None, 
              save = True,
              address = 24, coupling = 'ac', trace = 'aver'):
    """Takes a scan in spectrum analyzer and saves if bool parameter 
    is passed."""
    sa = icm.keysight_n9000b(24)
    par_array={'center_freq_Hz' : center_freq,
               'span_Hz' : span,
               'res_bw' : res_bw,
               'video_bw' : video_bw,
               'avg_no' : avg_no,
               'ref_level' : ref_level}
    sa.set_timeout(120e3)
    sa.set_frequency_span(span)
    sa.set_frequency_center(center_freq)
    sa.set_resolution_bandwidth(res_bw)
    sa.set_video_bandwidth(video_bw)
    sa_points = span/res_bw
    sa.set_sweep_points(sa_points + 1)
    sa.set_ref_level(ref_level)
    sa.toggle_continuous_sweep('off')
    sa.set_averaging(avg_no)
    sa.set_trace_type(trace)
    sa.set_input_coupling(coupling)
    data_string=sa.read_trace()
    data=np.fromstring(data_string, dtype=np.float, sep=',' )
    freq=data[::2]
    power=data[1::2]
    if save:
        np.savetxt(filename,np.column_stack((freq,power)))
    sa.set_input_coupling('ac')
    return (par_array, freq, power)
    
def lockin_scan(sens, tc, phase, srate_int, wait, spoints, save = False,
                filename = None):
    """Lockin collects data according to the parameters and transfers to remote
    interface."""
    par_array={'Sensitivity' : sens,
               'Time constant' : tc,
               'Phase reference' : phase,
               'Sample rate string' : srate_int,
               'Wait time' : wait,
               'Bin points' : spoints}
    lockin=gcm.SRS_844(8)
    lockin.set_timeout(60e3)
    lockin.set_sensitivity_nV(sens)
    lockin.set_time_constant_sec(tc)
    lockin.set_reference_phase_deg(phase)
    lockin.write('SRAT {}'.format(srate_int))
    lockin.write('REST')
    lockin.write('STRT')
    time.sleep(wait+1)
    lockin.write('PAUS')
    last_bin = int(lockin.query('SPTS?'))
    x = lockin.query('TRCA?1,{},{}'.format(int(last_bin-spoints),spoints))
    print('x extracted')
    y = lockin.query('TRCA?2,{},{}'.format(int(last_bin-spoints),spoints))
    print('y extracted')
    x = np.fromstring(x, dtype=np.float,sep=',')
    y = np.fromstring(y, dtype=np.float,sep=',')
    if save:
        np.savetxt(filename, np.column_stack((x[:-1],y[:-1])))
    return (par_array, x[:-1], y[:-1])

def best_gate_response(gate_check_array, source_flux, source_gate, 
                       voltage_divider_gate,
                       center_freq, span, power, 
                       points, av, el_delay, IF_bw, res_amp_best = 0,
                       flux_source_no = 0, gate_source_no = 1,
                       daq = None):
    """Takes vna scans for the corresponding set of gate values and picks the
    gate value with maximum dip in amplitude."""
    gate_save = gate_check_array[0]
    for j in range(len(gate_check_array)):
        flux=0
        bias_source_set_volt(source_flux,flux,flux_source_no, daq = daq)
        gate=gate_check_array[j]
        bias_source_set_volt(source_gate, gate, gate_source_no, 
                             volt_div = voltage_divider_gate, daq = daq)
        print('Gate input is {} V and flux input is {} V'.format(
                gate,flux))
        time.sleep(2.0)
        (par_array, freq, logmag, phase,
         res_amp_pos, res_amp, center_freq_ref) = vna_scan_save(
                 center_freq, span, power, points, av, 
                 el_delay, IF_bw)
        print(res_amp)
        if res_amp < -40:
            gate_save = gate
            break
        else:
            res_amp_best = np.min((res_amp_best,res_amp))
            if res_amp_best == res_amp:
                gate_save = gate
    return gate_save

"""----------------------------SWITCH FUNCTIONS-----------------------------"""
def vna_connect(switch_in):
    """Configuration switches for a vna scan. Connect cable from VNA P1 to
    4SPDTA0 and VNA P2 to 4SPDTC0.
    """
    switch_in.set_switch('A',0)
    switch_in.set_switch('C',0)
    
def error_sig_measure(switch_in, switch_sa,
                      mod_gen_address = mod_gen_address):
    """
    Pound locking setup.
    1. Connect input to phase modulator from signal generator, output to
    switch B COM and control from modulation source.
    2. Connect B0 to A1.
    3. Connect B1 to SP4T1.
    4. Connect C1 to bandpass filter in.
    5. Connect BPF out to power splitter.
    6. Connect power splitter one output to SP4T 2.
    7. Connect power splitter second output to DCOM through power detector.
    6. Connect SP4T COM to spectrum analyzer.
    7. Connect fridge in to A COM.
    8. Connect fridge out to C COM.
    9. Connect D0 to lock-in amplifier.
    10. Connect D1 to SP4T3.
    11. Connect VNA P1 to 4SPDTA0.
    12. Connect VNA P2 to 4SPDTC0.
    """
    sg=gcm.keysight_n5183b(25)
    modgen=scm.Novatech_409B(mod_gen_address)
    sg.toggle_output(0)
    modgen.set_voltage_Vpp(0,0)
    switch_in.set_switch('A',1)
    switch_in.set_switch('B',0)
    switch_in.set_switch('C',1)
    switch_sa.set_switch(2)
    switch_in.set_switch('D',0)
    
def power_detector_out(switch_in, switch_sa,
                       mod_gen_address = mod_gen_address):
    """
    Pound locking setup.
    1. Connect input to phase modulator from signal generator, output to
    switch B COM and control from modulation source.
    2. Connect B0 to A1.
    3. Connect B1 to SP4T1.
    4. Connect C1 to bandpass filter in.
    5. Connect BPF out to power splitter.
    6. Connect power splitter one output to SP4T 2.
    7. Connect power splitter second output to DCOM through power detector.
    6. Connect SP4T COM to spectrum analyzer.
    7. Connect fridge in to A COM.
    8. Connect fridge out to C COM.
    9. Connect D0 to lock-in amplifier.
    10. Connect D1 to SP4T3.
    11. Connect VNA P1 to 4SPDTA0.
    12. Connect VNA P2 to 4SPDTC0.
    """
    sg=gcm.keysight_n5183b(25)
    modgen=scm.Novatech_409B(mod_gen_address)
    sg.toggle_output(0)
    modgen.set_voltage_Vpp(0,0)
    switch_in.set_switch('A',1)
    switch_in.set_switch('B',0)
    switch_in.set_switch('C',1)
    switch_in.set_switch('D',1)
    switch_sa.set_switch(3)
    
def spec_connect_in(switch_in, switch_sa,
                    mod_gen_address = mod_gen_address):
    """
    1. Connect input to phase modulator from signal generator, output to
    switch B COM and control from modulation source.
    2. Connect B0 to A1.
    3. Connect B1 to SP4T1.
    4. Connect C1 to bandpass filter in.
    5. Connect BPF out to power splitter.
    6. Connect power splitter one output to SP4T 2.
    7. Connect power splitter second output to DCOM through power detector.
    6. Connect SP4T COM to spectrum analyzer.
    7. Connect fridge in to A COM.
    8. Connect fridge out to C COM.
    9. Connect D0 to lock-in amplifier.
    10. Connect D1 to SP4T3.
    11. Connect VNA P1 to 4SPDTA0.
    12. Connect VNA P2 to 4SPDTC0.
    """
    sg=gcm.keysight_n5183b(25)
    modgen=scm.Novatech_409B(mod_gen_address)
    sg.toggle_output(0)
    modgen.set_voltage_Vpp(0,0)
    switch_in.set_switch('A',0)
    switch_in.set_switch('B',1)
    switch_in.set_switch('C',0)
    switch_in.set_switch('D',1)
    switch_sa.set_switch(1)
    
def spec_connect_out(switch_in, switch_sa,
                     mod_gen_address = mod_gen_address):
    """
    Pound locking setup. Measure power spectrum of output.
    1. Connect input to phase modulator from signal generator, output to
    switch B COM and control from modulation source.
    2. Connect B0 to A1.
    3. Connect B1 to SP4T1.
    4. Connect C1 to bandpass filter in.
    5. Connect BPF out to power splitter.
    6. Connect power splitter one output to SP4T 2.
    7. Connect power splitter second output to DCOM through power detector.
    6. Connect SP4T COM to spectrum analyzer.
    7. Connect fridge in to A COM.
    8. Connect fridge out to C COM.
    9. Connect D0 to lock-in amplifier.
    10. Connect D1 to SP4T3.
    11. Connect VNA P1 to 4SPDTA0.
    12. Connect VNA P2 to 4SPDTC0.
    """
    sg=gcm.keysight_n5183b(25)
    modgen=scm.Novatech_409B(mod_gen_address)
    sg.toggle_output(0)
    modgen.set_voltage_Vpp(0,0)
    switch_in.set_switch('A',1)
    switch_in.set_switch('B',0)
    switch_in.set_switch('C',1)
    switch_sa.set_switch(2)
    switch_in.set_switch('D',0)

"""-------------------------BIAS SETTING FUNCTIONS--------------------------"""    
def set_bias_dcstate(channel,freq=0.1e-3,amp=0,mod_state=0):
    """Sets the voltage dc state in Keysight 33600A."""
    wavegen = ucm.keysight_33600A(addr=0)
    wavegen.set_waveform(channel,'SQU')
    freq_str = wavegen.set_frequency(channel,freq)
    volt_str = wavegen.set_voltage(channel,2*abs(amp))
    wavegen.write('SOUR{}:PHAS:SYNC'.format(channel))
    if amp < 0:
        wavegen.write('OUTP{}:POL INV'.format(channel))    
    else:
        wavegen.write('OUTP{}:POL NORM'.format(channel))
    pol_str = wavegen.query('OUTP{}:POL?'.format(channel)).split()[0]
#    voltoff_str = wavegen.set_voltage_offset(channel,amp/2)
    if mod_state:
        wavegen.write('SOUR{}:SUM:STAT 1'.format(channel))
    mod_state_str = wavegen.query('SOUR{}:SUM:STAT?'.format(channel))
    return (freq_str,0.5*volt_str,mod_state_str,pol_str)

def awg_pulse(volt_on, volt_off, freq, duty_cycle = 50, 
              volt_divider = 15, channel = 1, 
              address = 10):
    """Continously outputs a pulse from awg with corresponding
    high and low states."""
    if (abs(volt_on) > 75e-3) or (abs(volt_off) > 75e-3):
        sys.exit('Voltage amplitude too high.')
    awg = icm.tektronix_awg520(address)
    awg.toggle_fg_mode(1)
    awg.set_fg_shape(channel,'puls')
    awg.set_fg_frequency(freq)
    awg.set_fg_duty_cycle(channel, duty_cycle)
    volt_pp_in = (volt_on - volt_off)*volt_divider
    awg.set_fg_voltage_pp(channel, volt_pp_in)
    volt_offset_in = (volt_on + volt_off)*volt_divider/2
    awg.set_fg_voltage_offset(channel, volt_offset_in)
    awg.run()
    
def bias_source_set_volt(source, bias_value, source_no, volt_div =1,
                         daq = None):
    """Sets the bias voltage according to the voltage source under use."""
    if source == 'wavegen':
        set_bias_dcstate(source_no, amp=bias_value, mod_state=0)
    elif source == 'daq':
        if daq:
            daq.set_voltage(source_no, bias_value*volt_div)
        else:
            raise ValueError ('No DAQ instrument detected/ called.')
    elif source == 'awg':
        awg_dc(channel = source_no, volt_amp = bias_value, 
               volt_divider = volt_div)
        
def awg_dc(volt_amp, volt_divider = 15, channel = 1, address = 10):
    """Sets the DC bias output voltage in awg."""
    if abs(volt_amp) > 75e-3:
        sys.exit('Voltage amplitude too high.')
    awg = icm.tektronix_awg520(address)
    awg.toggle_fg_mode(1)
    awg.set_fg_shape(channel,'dc')
    volt_offset = volt_amp*volt_divider
    awg.set_fg_voltage_offset(channel, volt_offset)
    awg.run()
    awg.toggle_output(channel, 1)
    
def awg_sine(volt_amp, volt_off, freq, 
              volt_divider = 15, channel = 1, 
              address = 10):
    """Continously outputs a offset sine wave from awg with corresponding
    amplitude and offset."""
    if (abs(volt_amp) > 75e-3) or (abs(volt_off) > 75e-3):
        sys.exit('Voltage amplitude too high.')
    awg = icm.tektronix_awg520(address)
    awg.toggle_fg_mode(1)
    awg.set_fg_shape(channel,'sin')
    awg.set_fg_frequency(freq)
    volt_pp_in = 2*volt_amp*volt_divider
    awg.set_fg_voltage_pp(channel, volt_pp_in)
    volt_offset_in = volt_off*volt_divider
    awg.set_fg_voltage_offset(channel, volt_offset_in)
    awg.run()
    
def create_waveform_awg(volt_array, sample_rate, filename, channel = 1,
                        address = 10):
    awg = icm.tektronix_awg520(address)
    length = len(volt_array)
    marker_array = np.zeros(length)
    awg.set_frequency_reference(channel, 'ext')
    awg.send_waveform(volt_array, marker_array, marker_array,
                      filename, sample_rate)

"""------------------------------DATA ANALYSIS-----------------------------"""        
def interpolate_flux_bias(flux_array, center_freq_array, flux_sweep_no,
                          flux_max, directory_name, add_str = ''):
    """Extracts the flux periodicity in Phi_0."""
    flux_interp_fn = interpolate.interp1d(flux_array,center_freq_array,'cubic')
    flux_half_phi0 = 1
    flux_zero_phi0 = 1
    for i in range(1,flux_sweep_no):
        if (center_freq_array[i]<center_freq_array[i+1] and 
            center_freq_array[i]<center_freq_array[i-1]):
            if abs(flux_half_phi0) > abs(flux_array[i]):
                flux_half_phi0 = flux_array[i]
                flux_min_pos = i
        if (center_freq_array[i] > center_freq_array[i+1] and 
            center_freq_array[i] > center_freq_array[i-1]):
            if abs(flux_zero_phi0) > abs(flux_array[i]):
                flux_zero_phi0 = flux_array[i]
                flux_max_pos = i
    current_to_volt_fn = interpolate.interp1d([flux_array[0],
                                               flux_array[flux_sweep_no]],
                                              [round(-flux_max,3),
                                               round(flux_max,3)])
    
    flux_min_fine_array = np.linspace(flux_array[flux_min_pos-3],
                                      flux_array[flux_min_pos+3],50)
    flux_min_pos = np.where(flux_interp_fn(flux_min_fine_array) == np.min(
            flux_interp_fn(flux_min_fine_array)))[0]
    flux_half_phi0 = flux_min_fine_array[flux_min_pos[0]]
    
    flux_max_fine_array = np.linspace(flux_array[flux_max_pos-3],
                                      flux_array[flux_max_pos+3],50)
    flux_max_pos = np.where(flux_interp_fn(flux_max_fine_array) == np.max(
            flux_interp_fn(flux_max_fine_array)))[0]
    flux_zero_phi0 = flux_max_fine_array[flux_max_pos[0]]
    
    slope_volt, intercept_volt, r_value, p_value, std_err = stats.linregress(
                (0,0.5),(current_to_volt_fn(flux_zero_phi0),
                 current_to_volt_fn(flux_half_phi0)))
    (slope_current, intercept_current, r_value, 
     p_value, std_err) = stats.linregress(
                (flux_zero_phi0,flux_half_phi0),(0,0.5))
    flux_phi0_array = np.linspace(0,1,51)
    flux_volt_array = slope_volt*flux_phi0_array+intercept_volt
    np.savetxt(directory_name+'/flux_control_{}.dat'.format(add_str),
               np.column_stack((flux_phi0_array,flux_volt_array)))
    
    fig=plt.figure(figsize=(20,10))
    plt.scatter(flux_phi0_array,flux_volt_array,s=20)
    plt.plot(flux_phi0_array,flux_volt_array,'r')
    plt.xlabel('Flux (in $phi_0$)', fontsize=20)
    plt.ylabel('Flux input voltage (V)',fontsize=20)
    plt.savefig(directory_name+'/flux_control_{}.jpeg'.format(add_str))
    plt.close(fig)
    
    flux_array = np.asarray(flux_array)
    center_freq_array=np.asarray(center_freq_array)
    fig=plt.figure(figsize=(20,10))
    data_x=slope_current*flux_array+intercept_current
    plt.scatter(data_x,center_freq_array*1e-9,s=20)
    flux_current_array = np.linspace(flux_zero_phi0,flux_half_phi0)
    data_x=np.linspace(0,0.5)
    plt.plot(data_x,flux_interp_fn(flux_current_array)*1e-9,'r')
    plt.xlabel('Flux (in $phi_0$)', fontsize=20)
    plt.ylabel('Resonance in GHz',fontsize=20)
    plt.savefig(directory_name+'/flux_fit_{}.jpeg'.format(add_str))
    plt.close(fig)
    
    fig=plt.figure(figsize=(20,10))
    data_x=current_to_volt_fn(flux_array)
    plt.scatter(data_x,center_freq_array*1e-9,s=20)
    flux_current_array = np.linspace(flux_zero_phi0,flux_half_phi0)
    data_phi0=np.linspace(0,0.5)
    data_x=slope_volt*data_phi0+intercept_volt
    plt.plot(data_x,flux_interp_fn(flux_current_array)*1e-9,'r')
    plt.xlabel('Volt', fontsize=20)
    plt.ylabel('Resonance in GHz',fontsize=20)
    plt.savefig(directory_name+'/flux_fit_volt_{}.jpeg'.format(add_str))
    plt.close(fig)
    
    f = interpolate.interp1d(flux_phi0_array,flux_volt_array,'cubic')
    with open(
            directory_name+'/interpolated_flux_bias_conversion{}.pkl'.format(
                    add_str), 'wb') as k:
        pkl.dump(f, k)
    print(abs(flux_zero_phi0-flux_half_phi0))
    if (abs(flux_zero_phi0-flux_half_phi0)>30e-6 or
        abs(flux_zero_phi0-flux_half_phi0)<20e-6):
        prompt=input('''Flux periodicity ={} uA. Not accurate. Continue? 
                     (y to continue/exits otherwise)'''.format(
                     abs(flux_zero_phi0-flux_half_phi0)))
        if prompt=='y': #prompts to check the input values
            print('Continuing measurements..')
        else:
            raise ValueError('Flux interpolation not accurate.')
    return (slope_volt, intercept_volt)

def bg_from_flux_sweep(freq, center_freq_array, j, points, directory_name,
                       directory_name_broad, add_str = '',
                       string_pass = False, logmag_str = None,
                       phase_str = None, replace = 60):
    """Extracts background from the flux sweep scan by replacing the resonance
    response region with the corresponding values from another bias point.
    Not very accurate."""
    index_min=np.where(center_freq_array==np.min(center_freq_array))[0][0]
    index_max=np.where(center_freq_array==np.max(center_freq_array))[0][0]
    data_mag=np.zeros((2,points))
    data_phase=np.zeros((2,points))
    i=index_min
    if string_pass:
        filename_str_logmag = logmag_str.format(index_min)
        filename_str_phase = phase_str.format(index_min)
    else:
        filename_str_logmag=(str(directory_name_broad)+
                             '/flux_sweep_logmag_broad_'+str(j)+
                             '_th_gate_'+str(i)+'_th_flux')
        filename_str_phase=(str(directory_name_broad)+
                            '/flux_sweep_phase_broad_'+str(j)+
                            '_th_gate_'+str(i)+'_th_flux')
    filename_logmag="%s.dat" %filename_str_logmag
    filename_phase="%s.dat" %filename_str_phase
    data = np.loadtxt(filename_logmag)
    data_mag[0] = data[:,1]
    data = np.loadtxt(filename_phase)
    data_phase[0] = data[:,1]
    i=index_max
    if string_pass:
        filename_str_logmag = logmag_str.format(index_max)
        filename_str_phase = phase_str.format(index_max)
    else:
        filename_str_logmag=(str(directory_name_broad)+
                             '/flux_sweep_logmag_broad_'+str(j)+
                             '_th_gate_'+str(i)+'_th_flux')
        filename_str_phase=(str(directory_name_broad)+
                            '/flux_sweep_phase_broad_'+str(j)+
                            '_th_gate_'+str(i)+'_th_flux')
    filename_logmag="%s.dat" %filename_str_logmag
    filename_phase="%s.dat" %filename_str_phase
    data = np.loadtxt(filename_logmag)
    data_mag[1] = data[:,1]
    data = np.loadtxt(filename_phase)
    data_phase[1] = data[:,1]
    
    for i in range(2):
        if i==0:
            k=index_min
        else:
            k=index_max
        freq_diff = abs(freq-center_freq_array[k])
        resonance = np.where(freq_diff==np.min(freq_diff))[0][0]
        data_mag[i,resonance-replace:resonance+replace] = data_mag[
                (i+1)%2,resonance-replace:resonance+replace]
        data_phase[i,resonance-replace:resonance+replace] = data_phase[
                (i+1)%2,resonance-replace:resonance+replace]
        
    bg_mag = (data_mag[0]+data_mag[1])/2
    bg_phase = (data_phase[0]+data_phase[1])/2
    
    step = 10
    f_logmag = interpolate.interp1d(freq[::step],bg_mag[::step],'linear')
    with open((directory_name+'/interpolated_bg_logmag_{}.pkl'
               .format(add_str)), 'wb') as k:
        pkl.dump(f_logmag, k)
        
    f_phase = interpolate.interp1d(freq[::step],bg_phase[::step],'cubic')
    with open((directory_name+'/interpolated_bg_phase_{}.pkl'
               .format(add_str)), 'wb') as k:
        pkl.dump(f_phase, k)
    
    fig, ax = plt.subplots(figsize=(20, 10))
    plt.plot(freq[::step]*1e-9, f_logmag(freq[::step]), lw=4, color = 'r')
    plt.plot(freq*1e-9,bg_mag)
    plt.xlabel('Frequency (GHz)',fontsize = 30)
    plt.ylabel('|S11| (dB)',fontsize = 30)
    ax.spines['bottom'].set_linewidth(5)
    ax.spines['top'].set_linewidth(5)
    ax.spines['left'].set_linewidth(5)
    ax.spines['right'].set_linewidth(5)
    ax.tick_params(axis='x',direction='out', labelsize=40,width=2,length=10)
    ax.tick_params(axis='y',direction='out', labelsize=40,width=2,length=10)
    filename_plot =str(directory_name)+'/bg_mag{}.jpeg'.format(add_str)
    plt.savefig(filename_plot)
    
    
    fig, ax = plt.subplots(figsize=(20,10))
    plt.plot(freq[::step]*1e-9, f_phase(freq[::step]), lw=4, color = 'r')
    plt.plot(freq*1e-9,bg_phase)
    plt.xlabel('Frequency (GHz)',fontsize = 30)
    plt.ylabel('Phase (deg)',fontsize = 30)
    ax.spines['bottom'].set_linewidth(5)
    ax.spines['top'].set_linewidth(5)
    ax.spines['left'].set_linewidth(5)
    ax.spines['right'].set_linewidth(5)
    ax.tick_params(axis='x',direction='out', labelsize=40,width=2,length=10)
    ax.tick_params(axis='y',direction='out', labelsize=40,width=2,length=10)
    filename_plot = str(directory_name)+'/bg_phase.jpeg'.format(add_str)
    plt.savefig(filename_plot)  
    
    np.savetxt((str(directory_name)+'/bg_freq_logmag_phase{}.dat'.format(
            add_str)),np.column_stack((freq,bg_mag,bg_phase)))
        
    return (bg_mag, bg_phase, f_logmag, f_phase)

def interpolate_gate_bias(flux, gate_array, center_freq_array, 
                          directory_name, 
                          gate_max, add_str = '',
                          fine_points=50):
    """Extracts the gate periodicity in e. Best for gate values with a
    clear concave or convex response."""
    gate_reading_to_in_volt_fn = interpolate.interp1d(
            [gate_array[0], gate_array[-1]],
            [round(-gate_max,3), round(gate_max,3)])
    gate_in_array = gate_reading_to_in_volt_fn(gate_array)
    gate_interp_fn = interpolate.interp1d(gate_in_array,
                                          center_freq_array)
    min_point_pos_array=[]
    max_point_pos_array=[]
    zero_points={}
    
    gate_sweep_no = len(gate_array)
    for j in range(3,gate_sweep_no-2):
        if check_lessthan(center_freq_array[j],
                               center_freq_array[j-3:j+3]):
            min_point_pos_array.append(j) #adds the position of minimum
        if check_greaterthan(center_freq_array[j],
                                  center_freq_array[j-3:j+3]):
            max_point_pos_array.append(j) #adds the position of maximum
    for j in min_point_pos_array:#removes repeating values. selects the middle.
        if (j-1) in min_point_pos_array and (j-2) in min_point_pos_array:
            min_point_pos_array.remove(j)
            min_point_pos_array.remove(j-2)
    for j in min_point_pos_array:#removes repeating values. selects the left.
        if (j-1) in min_point_pos_array:
            min_point_pos_array.remove(j)
            
    for j in max_point_pos_array:#removes repeating values. selects the middle.
        if (j-1) in max_point_pos_array and (j-2) in max_point_pos_array:
            max_point_pos_array.remove(j)
            max_point_pos_array.remove(j-2)
    for j in max_point_pos_array:#removes repeating values. selects the left.
        if (j-1) in max_point_pos_array:
            max_point_pos_array.remove(j)
    
    print(min_point_pos_array)
    print(max_point_pos_array)      
    if flux%1==0: #gate response is convex
        if (len(min_point_pos_array) == 2) or (len(max_point_pos_array) == 2):
            zp_1, zp_2 = min_point_pos_array[:2]
        else:
            for j in max_point_pos_array:
                diff=gate_sweep_no        
                if (center_freq_array[j]-center_freq_array[
                        j+1] > center_freq_array[j]-center_freq_array[j-1]):
                #search left for actual minimum.
                    for k in min_point_pos_array:
                        if k < j and j-k < diff:
                            diff = j-k
                            zero_points['{}'.format(j)]=k
                            #final value will be the k value closest to j.
                else:
                #search right for actual minimum.
                    diff=gate_sweep_no
                    for k in min_point_pos_array:
                        if k > j and k-j < diff:
                            diff=k-j
                            zero_points['{}'.format(j)]=k
                            #final value will be the k value closest to j.
            zero_points_final=[]
            for x in zero_points.values():
                if x not in zero_points_final:
                    zero_points_final.append(x) #removes repetition.
            zp_1, zp_2 = zero_points_final[:2] #zero points are first two.
        
        gate_min_fine_array_1 = np.linspace(gate_in_array[zp_1-2],
                                            gate_in_array[zp_1+2],fine_points)
        gate_min_pos_1 = np.where(gate_interp_fn(
                gate_min_fine_array_1) == np.min(gate_interp_fn(
                        gate_min_fine_array_1)))[0]
        gate_0 = gate_min_fine_array_1[gate_min_pos_1[0]]
        
        gate_min_fine_array_2 = np.linspace(gate_in_array[zp_2-2],
                                            gate_in_array[zp_2+2], fine_points)
        gate_min_pos_2 = np.where(gate_interp_fn(
                gate_min_fine_array_2) == np.min(gate_interp_fn(
                        gate_min_fine_array_2)))[0]
        gate_2 = gate_min_fine_array_2[gate_min_pos_2[0]]
    else: #gate response is concave
        if (len(min_point_pos_array) == 2) or (len(max_point_pos_array) == 2):
            zp_1, zp_2 = max_point_pos_array[:2]
        else:
            for j in min_point_pos_array:
                diff=gate_sweep_no        
                if (center_freq_array[j]-center_freq_array[
                        j+1] < center_freq_array[j]-center_freq_array[j-1]):
                #search left for actual maximum.
                    for k in max_point_pos_array:
                        if k < j and j-k < diff:
                            diff = j-k
                            zero_points['{}'.format(j)]=k
                            #final value will be the k value closest to j.
                else:
                #search right for actual maximum.
                    diff=gate_sweep_no
                    for k in max_point_pos_array:
                        if k > j and k-j < diff:
                            diff=k-j
                            zero_points['{}'.format(j)]=k
                            #final value will be the k value closest to j.
            zero_points_final=[]
            for x in zero_points.values():
                if x not in zero_points_final:
                    zero_points_final.append(x)#removes repetition.
            zp_1,zp_2 = zero_points_final[:2] #zero points are first two.
        
        
        gate_max_fine_array_1 = np.linspace(gate_in_array[zp_1-2],
                                            gate_in_array[zp_1+2],fine_points)
        gate_max_pos_1 = np.where(gate_interp_fn(
                gate_max_fine_array_1) == np.max(gate_interp_fn(
                        gate_max_fine_array_1)))[0]
        gate_0 = gate_max_fine_array_1[gate_max_pos_1[0]]
        
        gate_max_fine_array_2 = np.linspace(gate_in_array[zp_2-3],
                                            gate_in_array[zp_2+3],fine_points)
        gate_max_pos_2 = np.where(gate_interp_fn(
                gate_max_fine_array_2) == np.max(gate_interp_fn(
                        gate_max_fine_array_2)))[0]
        gate_2 = gate_max_fine_array_2[gate_max_pos_2[0]]
                    
    if zp_1 < int(gate_sweep_no/4): 
        #chooses 0,2 such that 0 point is close to 0V.
        period_start=-2
        period_stop=0
    elif zp_2 > int(gate_sweep_no/4*3):
        period_start = 0
        period_stop=2
    else:
        period_start=0
        period_stop=2
    (slope_sample, intercept_sample, 
     r_value, p_value, std_err) = stats.linregress(
             (gate_0, gate_2),(period_start,period_stop))
    print(gate_0-gate_2)

    data_e=slope_sample*gate_in_array+intercept_sample
    fig, ax=plt.subplots(figsize=(12,8))
    plt.scatter(data_e, center_freq_array*1e-9, s=20)
    gate_sample_array = np.linspace(gate_0, gate_2)
    data_x=np.linspace(period_start, period_stop)
    plt.plot(data_x, gate_interp_fn(gate_sample_array)*1e-9,'r')
    plt.xlabel('Gate (in e)', fontsize=20)
    plt.ylabel('Resonance (GHz)',fontsize=20)
    ax.spines['bottom'].set_linewidth(5)
    ax.spines['top'].set_linewidth(5)
    ax.spines['left'].set_linewidth(5)
    ax.spines['right'].set_linewidth(5)
    ax.tick_params(axis='x',direction='out', labelsize=20,width=2,length=10)
    ax.tick_params(axis='y',direction='out', labelsize=20,width=2,length=10)
    plt.savefig(directory_name+'/gate_e_vs_freq{}.jpeg'.format(add_str))

    np.savetxt(directory_name+'/gate_V_vs_e{}.dat'.format(add_str),
               np.column_stack((
                       gate_in_array, data_e)))
    f = interpolate.interp1d(data_e, gate_in_array, 'cubic')
    filename_pkl = (directory_name+
                    '/interpolated_gate_bias_conversion{}.pkl'.format(add_str))
    with open(filename_pkl, 'wb') as k:
        pkl.dump(f, k)
    if (abs(gate_0-gate_2)>55e-3 or
        abs(gate_0-gate_2)<45e-3):
        prompt=input('''Gate periodicity ={} mV. Not accurate. Continue? 
                     (y to continue/exits otherwise)'''.format(
                     abs(gate_0-gate_2)))
        if prompt=='y': #prompts to check the input values
            print('Continuing measurements..')
        else:
            raise ValueError('Gate interpolation not accurate.')
    return (slope_sample,intercept_sample, data_e)

def extrapolate_gate_vs_resonance(ng_array, res_array):
    """Extrapolates the resonance in the quasiparticle region. Square root
    of out.best_values equals resonance."""
    mod = QuadraticModel()
    pars = mod.guess(res_array**2, x = ng_array)
    out = mod.fit(res_array**2, pars, x=ng_array)
    print(out.best_values)
    return out
    
"""------------------------------MISCELLANEOUS------------------------------"""        
def time_string():
    """Returns the current time in a string."""
    current_time=dt.datetime.now()
    time_str ='_{}_{}_{}_{}_'.format(current_time.month,
                current_time.day,current_time.hour,current_time.minute)
    return time_str

def date():
    """Returns the date."""
    current_time=dt.datetime.now()
    date='{}_{}_{}'.format(current_time.month,
      current_time.day,current_time.year)
    return date

def lineno():
    """Returns the current line number in program."""
    return inspect.currentframe().f_back.f_lineno

def check_lessthan(x, array):
    """If x is greater than any of the elements in array, returns False. If x
    is less than all of the elements in array, or equal to any of the
    elements in array, returns True."""
    for j in array:
        if x>j:
            return False
    else:
        return True
    
def check_greaterthan(x,array):
    """If x is less than any of the elements in array, returns False. If x
    is greater than all of the elements in array, or equal to any of the
    elements in array, returns True."""
    for j in array:
        if x<j:
            return False
    else:
        return True
    
def copy_file(source, destination):
    shutil.copy(source, destination) 
