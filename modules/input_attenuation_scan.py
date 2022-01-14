# -*- coding: utf-8 -*-
"""
Created on March 13 2021

@author: Sisira

Finds the input attenuation from VNA to the sample. Does this by checking the
dependance of Kerr renormalized shift to power.

1. For a set of bias points, does a power sweep in a span of 200 MHz.
2. The linear slope of resonant frequency vs power(W) is inversely 
proportional to the input attenuation. Eq (43) in cCPT characterization paper.
"""
"Import modules"

import numpy as np #Scientific computing module
import os #system functions
import matplotlib.pyplot as plt #Plotting module
from scipy import interpolate, stats #Scientific analysis module
import functions as fn
import pickle as pkl

def input_attenuation_scan(directory_name, parameters, inst_call):
    directory_name_att = directory_name + '/Input_attenuation_scan'
    os.mkdir(directory_name_att)
    directory_name_data = directory_name_att + '/data'
    os.mkdir(directory_name_data)
    line = fn.lineno
    [bias_point, temp, smooth_param, el_delay, 
     power_array_nW, center_freq_ref, kerr_value,
     kappa_int, kappa_ext, f_logmag, f_phase] = parameters
    power_array_nW_str = '{}'.format([k for k in power_array_nW])
     
    "Parameter file"

    par_array=[]
    line0 = line()
    par_array.append([line()-line0,"Bias point? ",bias_point])
    par_array.append([line()-line0,"Temperature(mK)",temp])
    par_array.append([line()-line0,"Channel number",1])
    par_array.append([line()-line0,"Number of traces",2])
    par_array.append([line()-line0,"Electrical delay?",el_delay])
    par_array.append([line()-line0,"VNA smooting aperture",smooth_param])
    par_array.append([line()-line0,"Number of sweep points (refined)",601])
    par_array.append([line()-line0,"Averaging number (refined)",100])
    par_array.append([line()-line0,"Span (refined) (in Hz)",30e6])
    par_array.append([line()-line0,"IF Bandwidth (refined)",5e3])
    par_array.append([line()-line0,"Center frequency", center_freq_ref])
    par_array.append([line()-line0,"Kerr value (MHz)", kerr_value])
    par_array.append([line()-line0,"Internal damping rate (Hz)", kappa_int])
    par_array.append([line()-line0,"External damping rate (Hz)", kappa_ext])
    par_array.append([line()-line0,"Power? (in nW)", power_array_nW_str])
    
    filename_parameter=str(directory_name_att)+'/Parameters.dat'
    np.savetxt(filename_parameter,np.row_stack((par_array)),delimiter=';',
               fmt='%s')
    scan_parameters=np.loadtxt(filename_parameter,delimiter=';',
                               usecols=2,dtype=np.str) 
    #loads the parameter values

    "Parameter values"

    scan_parameters = np.array([eval(i) for i in scan_parameters])
    [bias_point, temp, 
     channel, num_traces, el_delay, smooth_param,
     points, av, span, IF_bw, center_freq, kerr_value,
     kappa_int, kappa_ext, power_array_nW] = scan_parameters
     
    [na, sg, modgen, gate_meter, flux_meter, daq, awg, switch_in,
     switch_sa, bp_ctrl, twpa_pump, sa, lockin, dig] = inst_call
 
    channel=int(channel)
    num_traces=int(num_traces)
    points=int(points)
    print(bias_point)
    
    with open('interpolated_kappa_int_values.pkl', 'rb') as handle: #fn of ng,phi
        kappa_int_fn = pkl.load(handle)
    with open('interpolated_kappa_ext_values.pkl', 'rb') as handle: #fn of ng,phi
        kappa_ext_fn = pkl.load(handle)
    
    """--------------------------------MEASUREMENT -------------------------"""
    fn.vna_connect(switch_in)
    for tr in [1,2]:
        na.set_smoothing(smooth_param, channel = tr)
        na.toggle_smoothing(1, trace = tr)
#    if (center_freq_ref*1e-9 < 5.8083) and (center_freq_ref*1e9 > 5.699):
#        kappa_ext = kappa_ext_fn(center_freq_ref)
#        kappa_int = kappa_int_fn(center_freq_ref)
    kappa_tot = kappa_int+kappa_ext

    power_array_nW=np.asarray(power_array_nW)
    power_array = 10*np.log10(power_array_nW*1e-6)
    
    center_freq_array=[]
    for power in power_array:
        filename_str_logmag=(directory_name_data+
                             '/freq_sweep_logmag_ref_{}_dBm'.
                             format(round(power,2)))
        filename_str_phase=(directory_name_data+
                            '/freq_sweep_phase_ref_{}_dBm'.
                            format(round(power,2)))
        (par_array, freq_fine, logmag_fine, phase_fine,
         res_amp_pos_ref, res_amp, center_freq_final, 
         bg_reduced_data_mag, bg_reduced_data_phase) = fn.vna_scan_save(
                 center_freq, span, round(power,2), points, av, 
                 el_delay, IF_bw,
                 filename_str_logmag,filename_str_phase, 
                 save = True, plot = True, bg = True,
                 f_mag = f_logmag, f_phase = f_phase)
        center_freq_array.append(center_freq_final)
    center_freq_array=np.asarray(center_freq_array)
    fig=plt.figure(figsize=(20,10))
    plt.scatter(power_array_nW, center_freq_array*1e-9)
    slope, intercept, r_value, p_value, std_err = stats.linregress(
            power_array_nW, center_freq_array*1e-9)
    plt.xlabel('Power (nW)',fontsize = 30)
    plt.ylabel('Frequency (GHz)',fontsize = 30)
    linear_fit = slope*power_array_nW+intercept
    plt.plot(power_array_nW, linear_fit)
    filename_plot = directory_name_att+'/power_dependance.jpeg'
    plt.ylim((min(linear_fit)-1e-4,max(linear_fit)+1e-4))
    plt.title('Kerr shift = {} MHz'.format(kerr_value*1e-6))
    fig.savefig(filename_plot)
    plt.close(fig)
    
    
    f0=intercept*1e9
    h = 6.626e-34
    input_att = 1e-18*4*kerr_value*kappa_ext/(h*
                                              f0*slope*kappa_tot**2*np.pi*2)
    input_att_dB = np.log10(abs(input_att))*10
    
    n_array = 4*kappa_ext*power_array_nW*1e-9/(h*
                                               f0*kappa_tot**2*
                                               np.pi*2)/abs(input_att)
    fig, ax = plt.subplots(figsize=(20,10))
    ax.spines['bottom'].set_linewidth(5)
    ax.spines['top'].set_linewidth(5)
    ax.spines['left'].set_linewidth(5)
    ax.spines['right'].set_linewidth(5)
    ax.tick_params(axis='x',direction='out', labelsize=40,width=2,length=10)
    ax.tick_params(axis='y',direction='out', labelsize=40,width=2,length=10)
    plt.scatter(n_array, (center_freq_array-center_freq_ref)*1e-6, s = 60)
    plt.xlabel('Average Input photons',fontsize = 30)
    plt.ylabel('Frequency shift from {} GHz(MHz)'.format(
            round(center_freq_ref*1e-9,3)),fontsize = 30)
    plt.plot(n_array, (linear_fit-center_freq_ref*1e-9)*1e3, 'r', lw = 4)
    filename_plot = directory_name+'/photon_no_dependance.jpeg'
#    plt.ylim((min(center_freq_array*1e-9)-1e-4,
#              max(center_freq_array*1e-9)+1e-4))
    plt.title(('Input attenuation = {} dB,\n'.format(round(input_att_dB,2))+
               '$\Phi_{{ext}}$ = {} $\Phi_0$ and $n_g$ = {} e\n'.format(
                       bias_point[0][0], bias_point[0][1])+
                       'Kerr coefficient = {} MHz'.format(
                               round(kerr_value*1e-6,3))), 
                       fontsize=30, 
    y=0.87)
    plt.tight_layout()
    plt.savefig(filename_plot)
    
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
    filename_plot = directory_name+'/input_power_to_photons.jpeg'
    plt.title(('Input attenuation = {} dB,\n'.format(round(input_att_dB,2))+
               '$\Phi_{{ext}}$ = {} $\Phi_0$ and $n_g$ = {} e\n'.format(
                       bias_point[0][0], bias_point[0][1])+
                       'Kerr coefficient = {} MHz'.format(
                               round(kerr_value*1e-6,3))), 
                       fontsize=30, 
    y=0.87)
    plt.tight_layout()
    plt.savefig(filename_plot)

    filename_resonance=directory_name_att+'/resonance_data.dat'
    np.savetxt(filename_resonance, np.column_stack((center_freq_array,
                                                    power_array,
                                                    n_array)))
    filename_resonance=directory_name+'/no_of_photons.dat'
    np.savetxt(filename_resonance, np.column_stack((power_array, 
                                                    power_array_nW,
                                                    n_array)))

    return (input_att_dB, input_att, power_array, n_array)
