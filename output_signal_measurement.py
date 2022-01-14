# -*- coding: utf-8 -*-
"""
Measures the input signal going into the fridge.
"""
import numpy as np #Scientific computing module
import os #system functions
import time
import matplotlib.pyplot as plt #Plotting module
import functions as fn
    
def output_signal_scan(directory_name, parameters, inst_call,add_str=''):
    directory_name_out = directory_name + '/Output_signal'+add_str
    os.mkdir(directory_name_out)
    line = fn.lineno
    [carrier_power, res_freq_meas, mod_freq, mod_volt,
     bias, flux_fn, gate_fn, gate_source, gate_voltage_divider,
     photons] = parameters
     
    "Parameter file"

    par_array=[]
    line0 = line()
    par_array.append([line()-line0,"Carrier power (dBm)", carrier_power])
    par_array.append([line()-line0,"Modulation frequency",mod_freq])
    par_array.append([line()-line0,"Modulation voltage (Vpp)",mod_volt])
    par_array.append([line()-line0,"Spectrum analyzer center frequency",res_freq_meas])
    par_array.append([line()-line0,"Spectrum analyzer span",100e3])
    par_array.append([line()-line0,"Resolution bandwidth",1])
    par_array.append([line()-line0,"video bandwidth",1])
    par_array.append([line()-line0,"averaging",10])
    par_array.append([line()-line0,"Reference level (dBm)", -40])
    
    filename_parameter=str(directory_name_out)+'/Parameters.dat'
    np.savetxt(filename_parameter,np.row_stack((par_array)),delimiter=';',
               fmt='%s')
    scan_parameters=np.loadtxt(filename_parameter,delimiter=';',
                               usecols=2,dtype=np.str) 
    #loads the parameter values

    "Parameter values"

    scan_parameters = np.array([eval(i) for i in scan_parameters])
    [carrier_power, mod_freq, mod_volt, res_freq_meas,
     span, res_bw, video_bw, avg, ref_level] = scan_parameters
    
    [na, sg, modgen, gate_meter, flux_meter, daq, awg, switch_in,
     switch_sa, bp_ctrl, twpa_pump, sa, lockin, dig] = inst_call
    
    def avg_data(freq, power_sa):
        left = np.flip(power_sa[:int(span/2)])
        right = power_sa[int(span/2)+1:]
        avg_data = 10**((left+right)/2/10)
        signal_freq = freq[int(span/2)]
        freq_left = np.flip(signal_freq-freq[:int(span/2)])
        freq_right = freq[int(span/2)+1:]-signal_freq
        return (freq_left, freq_right, avg_data)
    
    """--------------------------------MEASUREMENT -------------------------"""
    fn.spec_connect_out(switch_in, switch_sa)
    
    flux = flux_fn([bias[0]+0.5])[0]
    fn.bias_source_set_volt('daq', flux, 0, daq = daq)
    gate = gate_fn([bias[1]])[0]
    fn.bias_source_set_volt(gate_source, gate, 1, 
                            volt_div = gate_voltage_divider, daq = daq)
    sg.set_mode()
    sg.set_frequency(res_freq_meas)
    sg.set_power(carrier_power)
    modgen.set_state_freq(0, 1,mod_freq*1e-6,mod_freq*1e-6)
    modgen.set_state_Vpp(0, 1,mod_volt,mod_volt)
    modgen.set_state_phase(0, 1 , 0.0, 0.0)
    modgen.set_state_Vpp(0, 1, mod_volt, mod_volt)
    sg.toggle_output(1)
    filename_data=(directory_name_out+'/fridge_out_off_resonance_sig.dat')
    (par_array, freq_off_res_sig, power_off_res_sig) = fn.spec_scan(
            res_freq_meas, 
            span,
            res_bw, video_bw, avg, ref_level, filename_data)
    off_res_sig_avg = avg_data(freq_off_res_sig, power_off_res_sig)
    modgen.set_state_Vpp(0, 1, 0, 0)
    sg.toggle_output(0)
    time.sleep(3)
    filename_data=(directory_name_out+'/fridge_out_off_resonance.dat')
    (par_array, freq_off_res, power_off_res) = fn.spec_scan(res_freq_meas, 
        span,
        res_bw, video_bw, avg, ref_level, filename_data)
    off_res_avg = avg_data(freq_off_res, power_off_res)
    
    flux = flux_fn([bias[0]])[0]
    fn.bias_source_set_volt('daq', flux, 0, daq = daq)
    time.sleep(3)
    filename_data=(directory_name_out+'/fridge_out_noise_on_resonance.dat')
    (par_array, freq_noise_on_res, power_noise_on_res) = fn.spec_scan(
            res_freq_meas, 
            span,
            res_bw, video_bw, avg, ref_level, filename_data)
    noise_on_res_avg = avg_data(freq_noise_on_res, power_noise_on_res)
    
    sg.set_mode()
    sg.set_frequency(res_freq_meas)
    sg.set_power(carrier_power)
    modgen.set_state_freq(0, 1,mod_freq*1e-6,mod_freq*1e-6)
    modgen.set_state_Vpp(0, 1,mod_volt,mod_volt)
    modgen.set_state_phase(0, 1 , 0.0, 0.0)
    modgen.set_state_Vpp(0, 1, mod_volt, mod_volt)
    sg.toggle_output(1)
    time.sleep(1)
    filename_data=(directory_name_out+'/fridge_out_signal.dat')
    (par_array, freq_sig, power_sig) = fn.spec_scan(
            res_freq_meas, 
            span,
            res_bw, video_bw, avg, ref_level, filename_data)
    sig_avg = avg_data(freq_sig, power_sig)
    
    fn.power_detector_out(switch_in, switch_sa)
    modgen.set_state_Vpp(0, 1, mod_volt, mod_volt)
    sg.toggle_output(1)
    time.sleep(1)
    filename_data=(directory_name_out+'/power_detector_out.dat')
    (par_array, freq_pd, power_pd) = fn.spec_scan(
            mod_freq, 
            span,
            res_bw, video_bw, avg, ref_level, filename_data)
    pd_peak = power_pd[int(span/2)]
    
    fig, ax = plt.subplots(figsize=(20,10))
    ax.spines['bottom'].set_linewidth(5)
    ax.spines['top'].set_linewidth(5)
    ax.spines['left'].set_linewidth(5)
    ax.spines['right'].set_linewidth(5)
    ax.tick_params(axis='x',direction='out', labelsize=40,width=2,length=10)
    ax.tick_params(axis='y',direction='out', labelsize=40,width=2,length=10)
    ax.set_xlabel('Frequency (MHz)',fontsize = 30)
    ax.set_ylabel('Power (dbm)',fontsize = 30)
    ax.plot(freq_pd*1e-6, power_pd)
    filename=(directory_name+'/power_detector_out.jpeg')
    plt.suptitle(('Carrier power = {} dBm, Modulation frequency = {} MHz\n'.format(
            carrier_power, mod_freq*1e-6)+
    'Signal measured = {} dBm, Photons = {}\n'.format(pd_peak, photons)+
    '$\Phi_{{ext}}$ = {} $\Phi_0$ and $n_g$ = {} e'.format(
                       bias[0], bias[1])), fontsize=30, y = 0.9)
    plt.savefig(filename)
    filename=(directory_name_out+'/power_detector_out.jpeg')
    plt.savefig(filename)
    
    fig, ax = plt.subplots(figsize=(20,10))
    ax.spines['bottom'].set_linewidth(5)
    ax.spines['top'].set_linewidth(5)
    ax.spines['left'].set_linewidth(5)
    ax.spines['right'].set_linewidth(5)
    ax.tick_params(axis='x',direction='out', labelsize=40,width=2,length=10)
    ax.tick_params(axis='y',direction='out', labelsize=40,width=2,length=10)
    ax.set_xlabel('Frequency (Hz)',fontsize = 30)
    ax.set_ylabel('Power (dBm)',fontsize = 30)
    ax.set_xscale('log')
    ax.set_yscale('log')
    start = 2
    plt.plot(off_res_avg[0][start:], off_res_avg[2][start:], 
             label = 'Off resonance noise', lw =4)
    plt.plot(off_res_sig_avg[0][start:], off_res_sig_avg[2][start:], 
             label = 'Off resonance signal', lw =4)
    plt.plot(noise_on_res_avg[0][start:], noise_on_res_avg[2][start:],
             label = 'On resonance noise', lw =4)
    plt.plot(sig_avg[0][start:], sig_avg[2][start:],
             label = 'On resonance with signal on', lw =4)
    filename=(directory_name+'/output_noise.jpeg')
    plt.legend(fontsize = 20, loc = 'best')
    plt.suptitle(('Carrier power = {} dBm, Modulation frequency = {} MHz\n'.
                  format(
            carrier_power, mod_freq*1e-6)+
    'Resonance frequency = {} GHz, Photons = {}\n'.format(
            round(res_freq_meas*1e-9,5),
                           photons)+
    '$\Phi_{{ext}}$ = {} $\Phi_0$ and $n_g$ = {} e'.format(
                       bias[0], bias[1])), fontsize=30, y = 0.9)
    plt.savefig(filename)    