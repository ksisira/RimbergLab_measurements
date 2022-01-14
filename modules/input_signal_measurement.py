# -*- coding: utf-8 -*-
"""
Measures the input signal going into the fridge.
"""
import numpy as np #Scientific computing module
import os #system functions
import time
import matplotlib.pyplot as plt #Plotting module
import functions as fn
    
def input_signal_scan(directory_name, parameters, inst_call):
    line = fn.lineno
    [carrier_power, res_freq_meas, mod_freq, mod_volt,
     scan_no] = parameters
    directory_name_pm_tone = (directory_name + 
                              '/Phase_modulated_input_tone_{}'.format(
                                      scan_no))
    os.mkdir(directory_name_pm_tone)
     
    "Parameter file"

    par_array=[]
    line0 = line()
    par_array.append([line()-line0,"Carrier power (dBm)", carrier_power])
    par_array.append([line()-line0,"Modulation frequency",mod_freq])
    par_array.append([line()-line0,"Modulation voltage (Vpp)",mod_volt])
    par_array.append([line()-line0,"Spectrum analyzer center frequency",res_freq_meas])
    par_array.append([line()-line0,"Spectrum analyzer span",10e3])
    par_array.append([line()-line0,"Resolution bandwidth",1])
    par_array.append([line()-line0,"video bandwidth",1])
    par_array.append([line()-line0,"averaging",10])
    par_array.append([line()-line0,"Peak threshold input", -95])
    par_array.append([line()-line0,"Reference level (dBm)", -40])
    
    filename_parameter=str(directory_name_pm_tone)+'/Parameters.dat'
    np.savetxt(filename_parameter,np.row_stack((par_array)),delimiter=';',
               fmt='%s')
    scan_parameters=np.loadtxt(filename_parameter,delimiter=';',
                               usecols=2,dtype=np.str) 
    #loads the parameter values

    "Parameter values"

    scan_parameters = np.array([eval(i) for i in scan_parameters])
    [carrier_power, mod_freq, mod_volt, res_freq_meas,
     span, res_bw, video_bw, avg, peak_threshold, ref_level] = scan_parameters
    
    [na, sg, modgen, gate_meter, flux_meter, daq, awg, switch_in,
     switch_sa, bp_ctrl, twpa_pump, sa, lockin, dig] = inst_call
    
    """--------------------------------MEASUREMENT -------------------------"""
    sg.set_mode()
    sg.set_frequency(res_freq_meas)
    sg.set_power(carrier_power)
    modgen.set_state_freq(0,1,mod_freq*1e-6,mod_freq*1e-6)
    modgen.set_state_Vpp(0,1,mod_volt,mod_volt)
    modgen.set_state_phase(0,1, 0.0, 0.0)
    
    fn.spec_connect_in(switch_in, switch_sa)
    modgen.set_state_Vpp(0,1,mod_volt,mod_volt)
    sg.toggle_output(1)
    time.sleep(1)
    filename_data=(directory_name_pm_tone+'/fridge_in_J0.dat')
    (par_array, freq_J0, power_J0) = fn.spec_scan(res_freq_meas, span,
    res_bw, video_bw, avg, ref_level, filename_data)
    peak_power_0 = round(power_J0[int(span/2)],2)
    
    filename_data=(directory_name_pm_tone+'/fridge_in_Jm1.dat')
    (par_array, freq_Jm1, power_Jm1) = fn.spec_scan(res_freq_meas-mod_freq, 
    span,
    res_bw, video_bw, avg, ref_level, filename_data)
    peak_power_m1 = round(power_Jm1[int(span/2)],2)
    
    filename_data=(directory_name_pm_tone+'/fridge_in_Jp1.dat')
    (par_array, freq_Jp1, power_Jp1) = fn.spec_scan(res_freq_meas+mod_freq, 
    span,
    res_bw, video_bw, avg, ref_level, filename_data)
    peak_power_p1 = round(power_Jp1[int(span/2)],2)
    print(peak_power_0, peak_power_m1, peak_power_p1)
    
    freq_J0 = freq_J0-res_freq_meas
    freq_Jm1 = freq_Jm1-res_freq_meas
    freq_Jp1 = freq_Jp1-res_freq_meas

    fig, (ax2, ax1, ax3) = plt.subplots(1, 3,figsize=(20,10), sharey=True)
    # hide the spines between ax1, ax2 and ax3
    ax1.spines['bottom'].set_linewidth(5)
    ax1.spines['top'].set_linewidth(5)
    ax2.spines['bottom'].set_linewidth(5)
    ax2.spines['top'].set_linewidth(5)
    ax3.spines['bottom'].set_linewidth(5)
    ax3.spines['top'].set_linewidth(5)
    ax2.spines['left'].set_linewidth(5)
    ax3.spines['right'].set_linewidth(5)
    
    ax2.tick_params(axis='x',direction='out', labelsize=40,width=2,length=10)
    ax1.tick_params(axis='x',direction='out', labelsize=40,width=2,length=10)
    ax3.tick_params(axis='x',direction='out', labelsize=40,width=2,length=10)
    ax2.tick_params(axis='y',direction='out', labelsize=40,width=2,length=10)
    ax1.set_xlabel('$\Delta \omega$ (MHz)',fontsize = 30)
    ax2.set_ylabel('Power (dbm)',fontsize = 30)
    ax2.set_xticks([-mod_freq*1e-6])
    ax1.set_xticks([0])
    ax3.set_xticks([mod_freq*1e-6])
    
    ax2.plot((freq_Jm1)*1e-6, power_Jm1,label='$J_{-1}$')
    ax1.plot((freq_J0)*1e-6, power_J0,label='$J_0$')
    ax3.plot((freq_Jp1)*1e-6, power_Jp1,label='$J_{1}$')
    filename_plot_in=(directory_name+'/powerspectrum_in_{}.jpeg'.format(
            scan_no))
    
    ax1.plot((freq_Jm1)*1e-6, power_Jm1,
             label='$J_{-1}$')
    ax3.plot((freq_Jm1)*1e-6, power_Jm1,
             label='$J_{-1}$')
    
    ax2.plot((freq_J0)*1e-6, power_J0,label='$J_0$')
    ax3.plot((freq_J0)*1e-6, power_J0,label='$J_0$')
    
    ax1.plot((freq_Jp1)*1e-6, power_Jp1,label='$J_{1}$')
    ax2.plot((freq_Jp1)*1e-6, power_Jp1,label='$J_{1}$')
    
    ax1.set_xlim(min(freq_J0*1e-6), max(freq_J0*1e-6))
    ax2.set_xlim(min(freq_Jm1*1e-6), max(freq_Jm1*1e-6))
    ax3.set_xlim(min(freq_Jp1*1e-6), max(freq_Jp1*1e-6))
    
    ax2.spines['right'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax3.spines['left'].set_visible(False)
    
    ax2.yaxis.tick_left()
    ax2.tick_params(labelright='off')
    ax1.tick_params(labelright='off')
    ax1.tick_params(labelleft='off')
    ax3.tick_params(labelleft='off')
    plt.tight_layout()
    plt.suptitle(('Carrier power = {} dBm, Modulation frequency = {} MHz\n'.format(
            carrier_power, mod_freq*1e-6)+
    'Modulation index = 1.08, $J_0$ power = {}\n'.format(peak_power_0)+
    r'$J_{{-1}}$ power = {}, $J_1$ power = {}'.format(peak_power_m1,
        peak_power_p1)),fontsize=30, y = 0.9)
    plt.savefig(filename_plot_in)
    
    return peak_power_0

    