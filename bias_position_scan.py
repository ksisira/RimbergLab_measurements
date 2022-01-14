# -*- coding: utf-8 -*-
"""
Created on March 7 2021

@author: Sisira

Measurement: Resonance data for a flux sweep near n_g = 0 and gate sweep at
given flux points.

1. At flux input  = 0V, picks a good gate point which has the minimum
amplitude. (from -10,-5,0,5,10 mV)
2. Do a flux sweep at the good gate point.
3. Plot 2D.
4. Gets the bias location using interpolation.
5. Extract bg from flux = 0 and flux = 0.5 Phi_0.
6. At the given flux points, do the gate scan.
7. Plot 2D resonance data for the varying gate.
8. Get the interpolation for gate sweep.
9. Save pickle files containing interpolation.

Manually set:
    Smoothing parameter for VNA.
"""
import os
import time
import functions as fn
import numpy as np
import matplotlib.pyplot as plt
import pickle as pkl

    
def get_bias_position(directory_name, parameters, inst_call):
    directory_name_bias_scan = directory_name + '/Bias_scan'
    os.mkdir(directory_name_bias_scan)
    directory_name_broad = directory_name_bias_scan + '/broad'
    os.mkdir(directory_name_broad)
    line = fn.lineno
    [temp, gate_source_no, gate_voltage_divider,
     vna_av_bias, vna_power_bias, smooth_param, el_delay,
     gate_check_array, gate_max_bias, gate_sweep_no_bias,
     flux_max_bias, flux_sweep_no_bias] = parameters

    "Parameter file"

    par_array=[]
    line0 = line()
    par_array.append([line()-line0,"Temperature(mK)",temp])
    par_array.append([line()-line0,"Gate source? (0 for DAQ and 1 for AWG)",gate_source_no])
    par_array.append([line()-line0,"Voltage divider for gate", gate_voltage_divider])
    par_array.append([line()-line0,"Channel number",1])
    par_array.append([line()-line0,"Number of traces",2])
    par_array.append([line()-line0,"Number of VNA sweep points (broad)",721])
    par_array.append([line()-line0,"Averaging number (broad)",vna_av_bias])
    par_array.append([line()-line0,"Power? (in dBm)",vna_power_bias])
    par_array.append([line()-line0,"Span (broad) (in Hz)",300e6])
    par_array.append([line()-line0,"Center frequency (broad)",5.75e9])
    par_array.append([line()-line0,"VNA smooting aperture",smooth_param])
    par_array.append([line()-line0,"Electrical delay?",el_delay])
    par_array.append([line()-line0,"IF Bandwidth (broad)",10e3])
    par_array.append([line()-line0,"Gate check points",gate_check_array])
    par_array.append([line()-line0,"Max gate voltage(in V)",gate_max_bias])
    par_array.append([line()-line0,"Number of sweep points for gate",gate_sweep_no_bias])
    par_array.append([line()-line0,"Max flux voltage (in V)",flux_max_bias])
    par_array.append([line()-line0,"Number of sweep points for flux",flux_sweep_no_bias])
    
    filename_parameter=str(directory_name_bias_scan)+'/Parameters.dat'
    np.savetxt(filename_parameter,np.row_stack((par_array)),delimiter=';',
               fmt='%s')
    bias_parameters=np.loadtxt(filename_parameter,delimiter=';',
                               usecols=2,dtype=np.str) 
    #loads the parameter values

    "Parameter values"

    bias_parameters = np.array([eval(i) for i in bias_parameters])
    [temp, gate_source_no, gate_voltage_divider,
     channel, num_traces, points, av, power, span, center_freq,smooth_p,
     el_delay, IF_bw, gate_check_array,
     gate_max, gate_sweep_no, flux_max, flux_sweep_no]=bias_parameters
     
    [na, sg, modgen, gate_meter, flux_meter, daq, awg, switch_in,
     switch_sa, bp_ctrl, twpa_pump, sa, lockin, dig] = inst_call
 
    channel=int(channel)
    num_traces=int(num_traces)
    flux_sweep_no=int(flux_sweep_no)
    gate_sweep_no=int(gate_sweep_no)
    points=int(points)
    gate_step=2*gate_max/gate_sweep_no
    flux_step=2*flux_max/flux_sweep_no
    gate_source_no = int(gate_source_no)

    """--------------------------------MEASUREMENT -------------------------"""
    fn.vna_connect(switch_in)
    for tr in [1,2]:
        na.toggle_smoothing(0, trace = tr)
    if gate_source_no == 1:
        gate_source = 'awg'
        awg.toggle_output(1,1)
    elif gate_source_no == 0:
        gate_source = 'daq'

    "Gate point pick"
    gate_save = fn.best_gate_response(
            gate_check_array, 'daq', gate_source, gate_voltage_divider,
            center_freq, span, power, points, av, 
            el_delay, IF_bw, daq = daq)
    
    "Flux dependance"

    gate_array=[]
    flux_array=[]
    center_freq_array=[]
    for i in range(flux_sweep_no+1):
        flux=round((-flux_max+(i*flux_step)),3)
        fn.bias_source_set_volt('daq', flux, 0, daq = daq)
#        daq.set_voltage(0, flux)
        for j in range(1):
            gate_sample=gate_save
            print(('Gate input is {} V and flux input is {} V'
                   .format(gate_sample,flux)))
            fn.bias_source_set_volt(gate_source, gate_sample, 1, 
                                    volt_div = gate_voltage_divider,
                                    daq = daq)
            time.sleep(2)
            flux_reading=flux_meter.get_current()
            gate_reading=gate_meter.get_voltage()
#            flux_reading = flux/11e3
#            gate_reading = gate_sample
            gate_array.append(gate_reading)
            flux_array.append(flux_reading)
        
            filename_str_logmag=(str(directory_name_broad)+
                                 '/flux_sweep_logmag_broad_'+str(j)+
                                 '_th_gate_'+str(i)+'_th_flux')
            filename_str_phase=(str(directory_name_broad)+
                                '/flux_sweep_phase_broad_'+str(j)+
                                '_th_gate_'+str(i)+'_th_flux')
            (par_array, freq, logmag, phase,
             res_amp_pos, res_amp, center_freq_ref) = fn.vna_scan_save(
                     center_freq, span, power, points, av, 
                     el_delay, IF_bw,
                     filename_str_logmag,filename_str_phase, 
                     save = True, plot = True)
            if i != 0:
                if center_freq_array[-1] == center_freq_ref:
                    center_freq_ref = center_freq_ref + 0.1*span/points
            center_freq_array.append(center_freq_ref)
    fn.bias_source_set_volt(gate_source, 0, 1, 
                            volt_div = gate_voltage_divider,
                            daq = daq)
#    daq.set_voltage(0,0.0)
    fn.bias_source_set_volt('daq', 0, 0, daq = daq)
    gate_array = np.asarray(gate_array)
    flux_array = np.asarray(flux_array)
    center_freq_array = np.asarray(center_freq_array)
    np.savetxt(directory_name_bias_scan+'/flux_sweep_bias_Volt_Amp_Hz.dat',
               np.column_stack((
            gate_array,flux_array,center_freq_array)))

    "Flux_plot"
    for j in range(1):
        fig, ax = plt.subplots(figsize=(12,8))
        x = flux_array
        x=np.asarray(x)
        logmag_array = np.zeros((flux_sweep_no+1,points))
        for i in range(flux_sweep_no+1):
            filename_str_logmag=(str(directory_name_broad)+
                                 '/flux_sweep_logmag_broad_'+str(j)+
                                 '_th_gate_'+str(i)+'_th_flux')
            filename_logmag="%s.dat" %filename_str_logmag
            data = np.loadtxt(filename_logmag)
            logmag_array[i] = data[:,1]
        plt.pcolormesh(x*1e6,freq*1e-9,logmag_array.T,cmap = plt.cm.coolwarm)
        cbar = plt.colorbar()
        plt.ylabel('$Frequency (GHz)$',fontsize=50)
        plt.xlabel('Flux ($\mu$A)',fontsize=50)
        cbar.set_label('Reflection coefficient',fontsize=40,labelpad=20)
        plt.title(('$V_g$ = {} mV, Input power (RT) = {} dBm, \n'.format(
                round(gate_reading*1e3,2), power) +
                  'Temperature = {} mK'.format(temp)), 
            fontsize = 20)
        ax.spines['bottom'].set_linewidth(5)
        ax.spines['top'].set_linewidth(5)
        ax.spines['left'].set_linewidth(5)
        ax.spines['right'].set_linewidth(5)
        ax.tick_params(axis='x',direction='out', labelsize=40,width=2,length=10)
        ax.tick_params(axis='y',direction='out', labelsize=40,width=2,length=10)
        cbar.ax.tick_params(direction='out', labelsize=40)
        fig.tight_layout()
        plt.savefig(directory_name+'/flux_dependence.jpeg')
        plt.close()

    "Flux bias location"

    slope_volt, intercept_volt = fn.interpolate_flux_bias(
            flux_array, center_freq_array, flux_sweep_no, flux_max, 
            directory_name_bias_scan)

    "Background"

    bg_mag, bg_phase, f_mag, f_phase = fn.bg_from_flux_sweep(
            freq, center_freq_array, j, points, directory_name,
            directory_name_broad)
#    bg_data = np.loadtxt(directory_name+'/bg_freq_logmag_phase.dat')
#    f_mag = bg_data[:,1]
#    f_phase = bg_data[:,2]

    "Gate sweep" 

    for i in range(1):
        flux_array=[]
        gate_array=[]
        center_freq_array=[]
        flux=0
        flux_in = flux*slope_volt+intercept_volt
#        daq.set_voltage(0,flux_in)
        fn.bias_source_set_volt('daq', flux_in, 0, daq = daq)
        for j in range(gate_sweep_no+1):
            gate_sample=(-gate_max+(j*gate_step))
            print('Gate input is {} V and flux input is {} V'.format(
                    gate_sample, flux_in))
            fn.bias_source_set_volt(gate_source, gate_sample, 1, 
                                    volt_div = gate_voltage_divider,
                                    daq = daq)
            time.sleep(2)
            flux_reading=flux_meter.get_current()
            gate_reading=gate_meter.get_voltage()
#            flux_reading = flux/11e3
#            gate_reading = gate_sample
            gate_array.append(gate_reading)
            flux_array.append(flux_reading)
        
            filename_str_logmag=(str(directory_name_broad)+
                                 '/gate_sweep_logmag_broad_'+str(j)+
                                 '_th_gate_'+str(i)+'_th_flux')
            filename_str_phase=(str(directory_name_broad)+
                                '/gate_sweep_phase_broad_'+
                                str(j)+'_th_gate_'+str(i)+'_th_flux')
            (par_array, freq, logmag, phase,
             res_amp_pos, res_amp, center_freq_ref, 
             bg_reduced_data_mag, bg_reduced_data_phase) = fn.vna_scan_save(
                     center_freq, span, power, points, av, 
                     el_delay, IF_bw,
                     filename_str_logmag,filename_str_phase, 
                     save = True, plot = True, bg = True,
                     f_mag = f_mag, f_phase = f_phase)
            center_freq_array.append(center_freq_ref)        
        fn.bias_source_set_volt(gate_source, 0, 1, 
                                volt_div = gate_voltage_divider,
                                daq = daq)
        gate_array = np.asarray(gate_array)
        flux_array = np.asarray(flux_array)
        center_freq_array = np.asarray(center_freq_array)
        np.savetxt((directory_name_bias_scan+
                    '/gate_sweep_bias_Volt_Amp_Hz_{}thflux.dat'.format(i)),
                    np.column_stack((gate_array,flux_array,center_freq_array)))
#    daq.set_voltage(0,0.0)
    fn.bias_source_set_volt('daq', 0, 0, daq = daq)
    if gate_source_no == 1:
        awg.toggle_output(1,0)
    fn.bias_source_set_volt(gate_source, 0, 1, 
                                    volt_div = gate_voltage_divider,
                                    daq = daq)
    
    "Gate plot and interpolate"
    for i in range(1):
        bias_array = np.loadtxt((directory_name_bias_scan+
                                 '/gate_sweep_bias_Volt_Amp_Hz_{}thflux.dat'
                                 .format(i)))
        x = bias_array[:,0]
        x=np.asarray(x)
        center_freq_array = bias_array[:,2]
        logmag_array = np.zeros((gate_sweep_no+1,points))
        for j in range(gate_sweep_no):
            filename_str_logmag=(str(directory_name_broad)+
                                 '/gate_sweep_logmag_broad_'+str(j)+
                                 '_th_gate_'+str(i)+'_th_flux')
            filename_logmag="%s.dat" %filename_str_logmag
            data = np.loadtxt(filename_logmag)
            logmag_array[j] = data[:,1] - bg_mag #background reduced
    
        fig, ax = plt.subplots(figsize=(12,8))
        plt.pcolormesh(x*1e3,freq*1e-9,logmag_array.T,cmap = plt.cm.coolwarm)
        cbar = plt.colorbar()
        plt.ylabel('$Frequency (GHz)$',fontsize=50)
        plt.xlabel('Gate (mV)',fontsize=50)
        cbar.set_label('Reflection coefficient - bg',fontsize=40,labelpad=20)
        plt.title(('Flux = 0 $\Phi_0$, Input power (RT)' +
                   '= {} dBm, Temperature = {} mK'.
                  format(power, temp)), fontsize = 20)
        ax.spines['bottom'].set_linewidth(5)
        ax.spines['top'].set_linewidth(5)
        ax.spines['left'].set_linewidth(5)
        ax.spines['right'].set_linewidth(5)
        ax.tick_params(axis='x',direction='out', labelsize=40,width=2,length=10)
        ax.tick_params(axis='y',direction='out', labelsize=40,width=2,length=10)
        cbar.ax.tick_params(direction='out', labelsize=40)
        ymax = np.max(center_freq_array*1e-9+10e-3)
        ymin = np.min(center_freq_array*1e-9-10e-3)
        plt.ylim((ymin,ymax))
        fig.tight_layout()
        plt.savefig(directory_name+'/gate_dependence_{}th_flux.jpeg'.format(i))
        plt.close()
        flux_sweep_array_phi0 = [0]
        if flux_sweep_array_phi0[i] == 0 or flux_sweep_array_phi0[i] == 0.5:
            (slope_sample,intercept_sample, data_e) = fn.interpolate_gate_bias(
                    flux_sweep_array_phi0[i], x, center_freq_array, 
                    directory_name_bias_scan, gate_max, 
                    add_str = '_{}th_flux'.format(i))

    "Bias points load"
    filename = (directory_name_bias_scan + 
                '/interpolated_gate_bias_conversion_0th_flux.pkl')
    with open(filename, 'rb') as handle:
        gate_fn = pkl.load(handle)
    filename = directory_name + '/interpolated_gate_bias_conversion'
    with open(filename+'.pkl', 'wb') as k:
        pkl.dump(gate_fn, k)

    filename = (directory_name_bias_scan+
                '/interpolated_flux_bias_conversion.pkl')
    with open(filename, 'rb') as handle:
        flux_fn = pkl.load(handle)
    filename = directory_name+'/interpolated_flux_bias_conversion'
    with open(filename+'.pkl', 'wb') as k:
        pkl.dump(flux_fn, k)    