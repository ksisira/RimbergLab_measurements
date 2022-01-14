# -*- coding: utf-8 -*-
"""
Measure fluctuations in lock-in amplifier
at off-resonance, on-resonance noise and signal on noise - both with and
without feedback.
Also measure the fluctuations using spectrum analyzer.

Manually set:
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
import numpy as np #Scientific computing module
import os #system functions
import time
import matplotlib.pyplot as plt #Plotting module
from scipy.fft import fft, fftfreq
 #Scientific analysis module
import functions as fn
import pickle as pkl
from scipy import interpolate, stats

def lockin_fft(directory_name, parameters, inst_call, feed_str = 'feed_on'):
    directory_name_fft = directory_name + '/Lockin_fft'+feed_str
    os.mkdir(directory_name_fft)
    line = fn.lineno
    [photons, bias, temp, mod_volt,
     sens, tc, phase_ref, expand, filt, 
     clock_rate, sample_rate, acq_no, 
     acq_time, res_freq_meas,
     carrier_power, mod_freq,
     gate_source_no, gate_voltage_divider, 
     PID_prompt] = parameters
     
    "Parameter file"

    par_array=[]
    line0 = line()
    par_array.append([line()-line0,"Modulation voltage (Vpp)",mod_volt])
    par_array.append([line()-line0,"Lockin sensitivity (in nV)",sens])
    par_array.append([line()-line0,"Lockin time constant (in sec)",tc])
    par_array.append([line()-line0,"Lockin start phase in deg",phase_ref])
    par_array.append([line()-line0,"Lockin expand",expand])
    par_array.append([line()-line0,"Lockin filter",filt])
    par_array.append([line()-line0,"Digitizer sampling frequency", clock_rate])
    par_array.append([line()-line0,"Acquisitions",acq_no])
    par_array.append([line()-line0,"Acquisition time per buffer (sec)",acq_time])
    par_array.append([line()-line0,"Data save sampling frequency", sample_rate])
    par_array.append([line()-line0,"Digitizer input range", 8])
    par_array.append([line()-line0,"Digitizer impedance", 1e6])
    par_array.append([line()-line0,"Spectrum analyzer span",100e3])
    par_array.append([line()-line0,"Resolution bandwidth",1])
    par_array.append([line()-line0,"video bandwidth",1])
    par_array.append([line()-line0,"averaging",20])
    par_array.append([line()-line0,"Reference level (dBm)", -40])
    par_array.append([line()-line0,"PID, UL and LL values", PID_prompt])
    par_array.append([line()-line0,"Bias", str(bias)])
    
    filename_parameter=str(directory_name_fft)+'/Parameters.dat'
    np.savetxt(filename_parameter,np.row_stack((par_array)),delimiter=';',
               fmt='%s')
    scan_parameters=np.loadtxt(filename_parameter,delimiter=';',
                               usecols=2,dtype=np.str) 
    #loads the parameter values

    "Parameter values"

    scan_parameters = np.array([eval(i) for i in scan_parameters])
    [mod_volt,
     sens, tc, phase_ref, expand, filt, clock_rate,
     acq_no, acq_time, sample_rate, input_range_V, impedance,
     span, res_bw, video_bw, avg, ref_level, PID_prompt,
     bias] = scan_parameters
     
    [na, sg, modgen, gate_meter, flux_meter, daq, awg, switch_in,
             switch_sa, bp_ctrl, twpa_pump, sa, lockin, dig] = inst_call
    
    impedance = int(impedance)
    acq_no = int(acq_no)
    
    with open(directory_name +
              '/interpolated_flux_bias_conversion.pkl', 'rb') as handle:
        flux_fn = pkl.load(handle)
        
    with open(directory_name +
              '/interpolated_gate_bias_conversion.pkl', 'rb') as handle:
        gate_fn = pkl.load(handle)
        
    if gate_source_no == 1:
        gate_source = 'awg'
        awg.toggle_output(1,1)
    elif gate_source_no == 0:
        gate_source = 'daq'      
        
    def avg_data(freq, power_sa):
        left = np.flip(power_sa[:int(span/2)])
        right = power_sa[int(span/2)+1:]
        avg_data = 10**((left+right)/2/10)
        signal_freq = freq[int(span/2)]
        freq_left = np.flip(signal_freq-freq[:int(span/2)])
        freq_right = freq[int(span/2)+1:]-signal_freq
        return (freq_left, freq_right, avg_data)
    
    """------------------------------MEASUREMENT ---------------------------"""
    dig.configure_board(input_range_V,input_range_V, ref_source = 'int',
                            trigger_timeout_sec = 10e-6,
                            impedance_1 = impedance, impedance_2 = impedance,
                            coupling_1 = 'dc', coupling_2 = 'dc',
                            sample_rate = clock_rate)
    samples = int((clock_rate*acq_time)//32)*32
    
    feed_prompt = input('''Switch on the feedback.
                        (y to continue/exits otherwise)''')
    if feed_prompt=='y': #prompts to check the input values
        print('Running..')
        
        acquire_data = dig.acquire(3, 0, samples, 1, acq_no,
                               input_range_V=input_range_V,
                               clock_rate=clock_rate, 
                               sample_rate=sample_rate)
        stitched_data=dig.stitch_data(acquire_data[0],acquire_data[1],1)
        l = len(stitched_data[0])
        dt = 1/sample_rate
        time_array = np.linspace(0,dt*l,l)
        data_feed_dc = (time_array, stitched_data[0], stitched_data[1])
        filename=(directory_name_fft+'/quadrature_signal_feed_dc')
        np.savetxt(filename+'.dat', np.column_stack((data_feed_dc[0],
                                                     data_feed_dc[1],
                                                     data_feed_dc[2])),
        encoding = 'utf-16')
    
    ng_array = np.linspace(-1,1,1001)
    gate_V_array = gate_fn(ng_array)
    (slope_V_ng, intercept_V_ng, 
     r, p, se) = stats.linregress(gate_V_array, ng_array)
    
    """PLOT 9"""
    
    fig, ax = plt.subplots(figsize=(20,10))
    ax.spines['bottom'].set_linewidth(5)
    ax.spines['top'].set_linewidth(5)
    ax.spines['left'].set_linewidth(5)
    ax.spines['right'].set_linewidth(5)
    ax.tick_params(axis='x',direction='out', labelsize=40,width=2,length=10)
    ax.tick_params(axis='y',direction='out', labelsize=40,width=2,length=10)
    plt.xlabel('Time (sec)',fontsize = 30)
    plt.ylabel('ng (e)',fontsize = 30)
    ng_dc_feed = (slope_V_ng*data_feed_dc[1]/gate_voltage_divider+
                  intercept_V_ng)
    ng_dc_feed_mean = np.mean(ng_dc_feed)
    plt.plot(data_feed_dc[0],
             ng_dc_feed, 
             label = 'feedback on', lw=4)
    plt.title(('Sensitivity = {:E} nV,'.format(sens)+
               'Time constant = {} s, \nFilter slope = {} dB/oct,'.format(
                       tc, filt)+
                       ' Photons = {}\n'.format(photons)+
                       '$\Phi_{{ext}}$ = {} $\Phi_0$ and $n_g$ = {} e'.format(
                       bias[0], bias[1])), fontsize=20, 
    y=0.65)
    plt.legend(fontsize = 20)
    filename = directory_name_fft + '/ng_dc_with_feedback'
    plt.savefig(filename+'.jpeg')
    plt.close()
    
    """PLOT 10"""
    fig, ax = plt.subplots(figsize=(20,10))
    ax.spines['bottom'].set_linewidth(5)
    ax.spines['top'].set_linewidth(5)
    ax.spines['left'].set_linewidth(5)
    ax.spines['right'].set_linewidth(5)
    ax.tick_params(axis='x',direction='out', labelsize=40,width=2,length=10)
    ax.tick_params(axis='y',direction='out', labelsize=40,width=2,length=10)
    plt.xlabel('Time (sec)',fontsize = 30)
    plt.ylabel('Lock-in output Y (V)',fontsize = 30)
    plt.plot(data_feed_dc[0],
             data_feed_dc[2], 
             label = 'feedback on', lw=4)
    plt.title(('Sensitivity = {:E} nV,'.format(sens)+
               'Time constant = {} s, \nFilter slope = {} dB/oct,'.format(
                       tc, filt)+
                       ' Photons = {}\n'.format(photons)+
                       '$\Phi_{{ext}}$ = {} $\Phi_0$ and $n_g$ = {} e'.format(
                       bias[0], bias[1])), fontsize=20, 
    y=0.65)
    plt.legend(fontsize = 20)
    filename = directory_name_fft + '/Y_Quadrature_dc_feedback'
    plt.savefig(filename+'.jpeg')
    plt.close()
    
    time_array = data_feed_dc[0]
    N = len(time_array)
    T = 1/sample_rate
    tf = fftfreq(N, T)[:N//2]
    
    ngf_feed_dc = fft(data_feed_dc[1]/gate_voltage_divider*1e3)
    yf_feed_dc = fft(data_feed_dc[2])
    
    valid_feed_ng_dc = 2/N*np.abs(ngf_feed_dc[:N//2])
    valid_feed_y_dc = 2/N*np.abs(yf_feed_dc[:N//2])
    
    
    
    """PLOT 16"""
    fig, ax1 = plt.subplots(figsize=(20,10))
    ax1.spines['bottom'].set_linewidth(5)
    ax1.spines['top'].set_linewidth(5)
    ax1.spines['left'].set_linewidth(5)
    ax1.spines['right'].set_linewidth(5)
    ax1.tick_params(axis='x',direction='out', labelsize=40,width=2,length=10)
    ax1.tick_params(axis='y',direction='out', labelsize=40,width=2,length=10)
    ax1.loglog(tf[1:], valid_feed_ng_dc[1:],
               label = 'Feedback on', lw=4)
    plt.title(('Sensitivity = {:E} nV,'.format(sens)+
               'Time constant = {} s, \nFilter slope = {} dB/oct,'.format(
                       tc, filt)+
                       ' Photons = {}\n'.format(photons)+
                       '$\Phi_{{ext}}$ = {} $\Phi_0$ and $n_g$ = {}e'.
                       format(
                       bias[0], bias[1])), fontsize=20, 
    y=0.65)
    plt.legend(fontsize = 20)
    ax1.set_xscale('log')
    ax1.set_xlabel('Frequency (Hz)',fontsize = 30)
    ax1.set_ylabel('FFT ng', fontsize = 30)
    filename = directory_name_fft + '/Feedback_FFT_ng_DC'
    plt.savefig(filename+'.jpeg')
    plt.close()
    
    """PLOT 17"""
    fig, ax3 = plt.subplots(figsize=(20,10))
    ax3.spines['bottom'].set_linewidth(5)
    ax3.spines['top'].set_linewidth(5)
    ax3.spines['left'].set_linewidth(5)
    ax3.spines['right'].set_linewidth(5)
    ax3.tick_params(axis='x',direction='out', labelsize=40,width=2,length=10)
    ax3.tick_params(axis='y',direction='out', labelsize=40,width=2,length=10)
    ax3.set_xlabel('Frequency (Hz)',fontsize = 30)
    ax3.set_ylabel('FFT Y', fontsize = 30)
    ax3.set_xscale('log')
    ax3.loglog(tf[1:], valid_feed_y_dc[1:],
               label = 'Feedback on Y', lw=4)
    plt.title(('Sensitivity = {:E} nV,'.format(sens)+
               'Time constant = {} s, \nFilter slope = {} dB/oct,'.format(
                       tc, filt)+
                       ' Photons = {}\n'.format(photons)+
                       '$\Phi_{{ext}}$ = {} $\Phi_0$ and $n_g$ = {}e'.
                       format(
                       bias[0], bias[1])), fontsize=20, 
    y=0.65)
    plt.legend(fontsize = 20)
    filename = directory_name_fft + '/Feedback_FFT_Y_DC'
    plt.savefig(filename+'.jpeg')
    plt.close()
    
    np.savetxt(directory_name_fft +'/fft.dat', np.column_stack((
            tf, valid_feed_ng_dc, valid_feed_y_dc)))
    
    return np.round(ng_dc_feed_mean,3)