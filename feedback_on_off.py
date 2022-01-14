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
    
    filename_data=(directory_name+'/Output_signal/fridge_out_signal.dat')
    output_noise_prev = np.loadtxt(filename_data)
    sig_avg = avg_data(output_noise_prev[:,0], output_noise_prev[:,1])
    
    """------------------------------MEASUREMENT ---------------------------"""
    fn.error_sig_measure(switch_in, switch_sa)
    #Error signal switch setting also sets the spectrum analyzer to measure the
    #output spectrum.
    
    "Take the off-resonance noise by moving the bias point away."
    if bias[0]<0.5:
        flux_off_res = 0.5
    else:
        flux_off_res = 0
    print(flux_off_res)
    flux = flux_fn([flux_off_res])[0]
    fn.bias_source_set_volt('daq', flux, 0, daq = daq)
    gate = gate_fn([bias[1]])[0]
    fn.bias_source_set_volt(gate_source, gate, 1, 
                            volt_div = gate_voltage_divider, daq = daq)
    time.sleep(3)
#    lockin.set_sensitivity_nV(sens)
    lockin.set_time_constant_sec(tc)
#    lockin.set_reference_phase_deg(phase_ref)
#    lockin.set_filter_slope(filt)
    
    #Initially set the offset to zero.
#    lockin.set_expand(1, 1)
#    lockin.set_expand(2, 1)
#    lockin.set_offset('X',0)
#    lockin.set_offset('Y',0)
#    time.sleep(3)
#    #omit last few points as digitizer measures wrong values towards the end.
    points_omit = int(sample_rate*0.1)
#    
#    #Digitizer impedance set to 1 MOhm because the PID is 1MOhm impedance.
#    #Internal reference as we are not worried about phase locking upto 1 kHz.
#    #No trigger set either.
#    dig.configure_board(input_range_V,input_range_V, ref_source = 'int',
#                            trigger_timeout_sec = 10e-6,
#                            impedance_1 = impedance, impedance_2 = impedance,
#                            coupling_1 = 'dc', coupling_2 = 'dc',
#                            sample_rate = clock_rate)
    samples = int((clock_rate*acq_time)//32)*32
#    acquire_data = dig.acquire(3, 0, samples, 1, acq_no,
#                               input_range_V=input_range_V,
#                               clock_rate=clock_rate, 
#                               sample_rate=sample_rate)
#    stitched_data=dig.stitch_data(acquire_data[0],acquire_data[1],1)
#    data_off_res_0_offset = (acquire_data[2], stitched_data[0], stitched_data[1])
#    filename=(directory_name_fft+'/quadrature_noise_off_resonance_0_offset')
#    np.savetxt(filename+'.dat', np.column_stack((data_off_res_0_offset[0],
#                                                 data_off_res_0_offset[1],
#                                                 data_off_res_0_offset[2])),
#    encoding = 'utf-16')
#    
#    #Get the mean data. Fix the y-mean value as offset.
#    ng_noise_mean = np.mean(data_off_res_0_offset[1]
#    [:-points_omit])/gate_voltage_divider
#    y_noise_mean = np.mean(data_off_res_0_offset[2][:-points_omit])
#    ng_noise_err = np.std(data_off_res_0_offset[1]
#    [:-points_omit])/gate_voltage_divider
#    y_noise_err = np.std(data_off_res_0_offset[1][:-points_omit])
#    offset_y = y_noise_mean/10*100
#    print(y_noise_mean, offset_y)
#    lockin.set_offset('Y', -offset_y)
#    lockin.set_expand(1, expand)
#    lockin.set_expand(2, expand)
#    time.sleep(3)
#    
    acquire_data = dig.acquire(3, 0, samples, 1, acq_no,
                               input_range_V=input_range_V,
                               clock_rate=clock_rate, 
                               sample_rate=sample_rate)
    stitched_data=dig.stitch_data(acquire_data[0],acquire_data[1],1)
    data_off_res_1_offset = (acquire_data[2], stitched_data[0], stitched_data[1])
    filename=(directory_name_fft+'/quadrature_noise_off_resonance_1_offset')
    np.savetxt(filename+'.dat', np.column_stack((data_off_res_1_offset[0],
                                                 data_off_res_1_offset[1],
                                                 data_off_res_1_offset[2])),
    encoding = 'utf-16')
    
#    """PLOT 1"""
#    plt.figure()
#    plt.plot(data_off_res_0_offset[0][:-points_omit], 
#             data_off_res_0_offset[2][:-points_omit], label = 'y_raw')
#    plt.plot(data_off_res_1_offset[0][:-points_omit],
#             data_off_res_1_offset[2][:-points_omit], label = 'y_offset')
#    plt.title('Mean ng and Y = {}, {} mV\n and STD = {}, {} mV'.format(
#            round(ng_noise_mean*1e3,2), round(y_noise_mean*1e3,2),
#            round(ng_noise_err*1e3,2), 
#            round(y_noise_err*1e3),2))
#    plt.legend()
#    filename=(directory_name_fft+'/quadrature_offset')
#    plt.savefig(filename+'.jpeg')
#    plt.close()
#    
    "Move the bias point back to resonance and measure the noise."
    flux = flux_fn([bias[0]])[0]
    fn.bias_source_set_volt('daq', flux, 0, daq = daq)
    time.sleep(3)
    
    acquire_data = dig.acquire(3, 0, samples, 1, acq_no,
                               input_range_V=input_range_V,
                               clock_rate=clock_rate, 
                               sample_rate=sample_rate)
    stitched_data=dig.stitch_data(acquire_data[0],acquire_data[1],1)
    data_on_res_0_sig = (acquire_data[2], stitched_data[0], stitched_data[1])
    filename=(directory_name_fft+'/quadrature_noise_on_resonance')
    np.savetxt(filename+'.dat', np.column_stack((data_on_res_0_sig[0],
                                                 data_on_res_0_sig[1],
                                                 data_on_res_0_sig[2])),
    encoding = 'utf-16')
    
    "Switch on the signal and measure fluctuations on resonance."
    sg.set_mode()
    sg.set_frequency(res_freq_meas)
    sg.set_power(carrier_power)
    sg.toggle_output(1)
    modgen.set_state_freq(0,1,mod_freq*1e-6,mod_freq*1e-6)
    modgen.set_state_Vpp(0,1,mod_volt,mod_volt)
    modgen.set_state_phase(0,1,0.0,0.0)
    time.sleep(3)
    
    acquire_data = dig.acquire(3, 0, samples, 1, acq_no,
                               input_range_V=input_range_V,
                               clock_rate=clock_rate, 
                               sample_rate=sample_rate)
    stitched_data=dig.stitch_data(acquire_data[0],acquire_data[1],1)
    data_on_res_1_sig = (acquire_data[2], stitched_data[0], stitched_data[1])
    filename=(directory_name_fft+'/quadrature_signal_on_resonance')
    np.savetxt(filename+'.dat', np.column_stack((data_on_res_1_sig[0],
                                                 data_on_res_1_sig[1],
                                                 data_on_res_1_sig[2])),
    encoding = 'utf-16')
    
    #Measure the fluctuations with signal on resonance using spectrum analyzer.
    filename_data=(directory_name_fft+'/fridge_out_signal_on_resonance.dat')
    (par_array, freq_sig, power_sig) = fn.spec_scan(
            res_freq_meas, 
            span,
            res_bw, video_bw, avg, ref_level, filename_data)
    sig_avg_on_res = avg_data(freq_sig, power_sig)
    
    """PLOT 2"""
    
    fig, ax = plt.subplots(figsize=(20,10))
    ax.spines['bottom'].set_linewidth(5)
    ax.spines['top'].set_linewidth(5)
    ax.spines['left'].set_linewidth(5)
    ax.spines['right'].set_linewidth(5)
    ax.tick_params(axis='x',direction='out', labelsize=40,width=2,length=10)
    ax.tick_params(axis='y',direction='out', labelsize=40,width=2,length=10)
    plt.xlabel('Time (sec)',fontsize = 30)
    plt.ylabel('Lock-in output Y (V)',fontsize = 30)
    plt.plot(data_off_res_1_offset[0][:-points_omit],
             data_off_res_1_offset[2][:-points_omit], 
             label = 'off resonance noise', lw=4)
    plt.plot(data_on_res_0_sig[0][:-points_omit],
             data_on_res_0_sig[2][:-points_omit], 
             label = 'on resonance noise', lw=4)
    plt.plot(data_on_res_1_sig[0][:-points_omit],
             data_on_res_1_sig[2][:-points_omit], 
             label = 'on resonance signal', lw=4)
    plt.title(('Sensitivity = {:E} nV,'.format(sens)+
               'Time constant = {} s, \nFilter slope = {} dB/oct,'.format(
                       tc, filt)+
                       ' Photons = {}\n'.format(photons)+
                       '$\Phi_{{ext}}$ = {} $\Phi_0$ and $n_g$ = {} e'.format(
                       bias[0], bias[1])), fontsize=20, 
    y=0.65)
    plt.legend(fontsize = 20)
    filename = directory_name_fft + '/Y_Quadrature_Setup'
    plt.savefig(filename+'.jpeg')
    plt.close()
    
    
    "Get fluctuations off resonance with signal on."
    flux = flux_fn([flux_off_res])[0]
    fn.bias_source_set_volt('daq', flux, 0, daq = daq)
    time.sleep(3)
    
    acquire_data = dig.acquire(3, 0, samples, 1, acq_no,
                               input_range_V=input_range_V,
                               clock_rate=clock_rate, 
                               sample_rate=sample_rate)
    stitched_data=dig.stitch_data(acquire_data[0],acquire_data[1],1)
    data_off_res_1_sig = (acquire_data[2], stitched_data[0], stitched_data[1])
    filename=(directory_name_fft+'/quadrature_signal_off_resonance')
    np.savetxt(filename+'.dat', np.column_stack((data_off_res_1_sig[0],
                                                 data_off_res_1_sig[1],
                                                 data_off_res_1_sig[2])),
    encoding = 'utf-16')
    
    #Measure the fluctuations with signal off-resonance using spectrum analyzer.
    filename_data=(directory_name_fft+'/fridge_out_signal_off_resonance.dat')
    (par_array, freq_sig, power_sig) = fn.spec_scan(
            res_freq_meas, 
            span,
            res_bw, video_bw, avg, ref_level, filename_data)
    sig_avg_off_res = avg_data(freq_sig, power_sig)
    
    "Get the bias point back to resonance."
    flux = flux_fn([bias[0]])[0]
    fn.bias_source_set_volt('daq', flux, 0, daq = daq)
    time.sleep(3)
    
    feed_prompt = input('''Switch on the feedback.
                        (y to continue/exits otherwise)''')
    if feed_prompt=='y': #prompts to check the input values
        print('Running..')
        
        acquire_data = dig.acquire(3, 0, samples, 1, acq_no,
                               input_range_V=input_range_V,
                               clock_rate=clock_rate, 
                               sample_rate=sample_rate)
        stitched_data=dig.stitch_data(acquire_data[0],acquire_data[1],1)
        data_feed_dc = (acquire_data[2], stitched_data[0], stitched_data[1])
        filename=(directory_name_fft+'/quadrature_signal_feed_dc')
        np.savetxt(filename+'.dat', np.column_stack((data_feed_dc[0],
                                                     data_feed_dc[1],
                                                     data_feed_dc[2])),
        encoding = 'utf-16')
        #Measure power spectrum with feedback.
        filename_data=(directory_name_fft+'/fridge_out_signal_feed_dc.dat')
        (par_array, freq_sig, power_sig) = fn.spec_scan(
                res_freq_meas, 
                span,
                res_bw, video_bw, avg, ref_level, filename_data)
        sig_avg_feed_dc = avg_data(freq_sig, power_sig)   
    
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
    ng_dc_feed = (slope_V_ng*data_feed_dc[1]
    [:-points_omit]/gate_voltage_divider+intercept_V_ng)
    ng_dc_feed_mean = np.mean(ng_dc_feed)
    ng_sig = (slope_V_ng*data_on_res_1_sig[1]
    [:-points_omit]/gate_voltage_divider+intercept_V_ng)
    plt.plot(data_feed_dc[0][:-points_omit],
             ng_dc_feed, 
             label = 'feedback on', lw=4)
    plt.plot(data_on_res_1_sig[0][:-points_omit],
             ng_sig, 
             label = 'feedback off', lw=4)
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
    plt.plot(data_on_res_1_sig[0][:-points_omit],
             data_on_res_1_sig[2][:-points_omit], 
             label = 'feedback off', lw=4)
    plt.plot(data_feed_dc[0][:-points_omit],
             data_feed_dc[2][:-points_omit], 
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
    
    time_array = data_off_res_1_offset[0][:-points_omit]
    N = len(time_array)
    T = 1/sample_rate
    tf = fftfreq(N, T)[:N//2]
    
    ngf_off_res_noise = fft(data_off_res_1_offset[1]
    [:-points_omit]/gate_voltage_divider*1e3)
    yf_off_res_noise = fft(data_off_res_1_offset[2][:-points_omit])
    ngf_on_res_noise = fft(data_on_res_0_sig[1]
    [:-points_omit]/gate_voltage_divider*1e3)
    yf_on_res_noise = fft(data_on_res_0_sig[2][:-points_omit])
    ngf_on_res_sig = fft(data_on_res_1_sig[1]
    [:-points_omit]/gate_voltage_divider*1e3)
    yf_on_res_sig = fft(data_on_res_1_sig[2][:-points_omit])
    ngf_feed_dc = fft(data_feed_dc[1][:-points_omit]/gate_voltage_divider*1e3)
    yf_feed_dc = fft(data_feed_dc[2][:-points_omit])
    ngf_off_res_sig = fft(data_off_res_1_sig[1]
    [:-points_omit]/gate_voltage_divider*1e3)
    yf_off_res_sig = fft(data_off_res_1_sig[2][:-points_omit])
    
    valid_off_res_ng = 2/N*np.abs(ngf_off_res_noise[:N//2])
    valid_on_res_ng = 2/N*np.abs(ngf_on_res_noise[:N//2])
    valid_sig_ng = 2/N*np.abs(ngf_on_res_sig[:N//2])
    valid_off_res_y = 2/N*np.abs(yf_off_res_noise[:N//2])
    valid_on_res_y = 2/N*np.abs(yf_on_res_noise[:N//2])
    valid_sig_y = 2/N*np.abs(yf_on_res_sig[:N//2])
    
    valid_feed_ng_dc = 2/N*np.abs(ngf_feed_dc[:N//2])
    valid_feed_y_dc = 2/N*np.abs(yf_feed_dc[:N//2])
    
    valid_off_res_ng_sig = 2/N*np.abs(ngf_off_res_sig[:N//2])
    valid_off_res_y_sig = 2/N*np.abs(yf_off_res_sig[:N//2])
    
    """PLOT 11"""
    
    fig, ax = plt.subplots(figsize=(20,10))
    ax.spines['bottom'].set_linewidth(5)
    ax.spines['top'].set_linewidth(5)
    ax.spines['left'].set_linewidth(5)
    ax.spines['right'].set_linewidth(5)
    ax.tick_params(axis='x',direction='out', labelsize=40,width=2,length=10)
    ax.tick_params(axis='y',direction='out', labelsize=40,width=2,length=10)
    plt.xlabel('Frequency (Hz)', fontsize = 30)
    plt.ylabel('FFT Y', fontsize = 30)
    plt.loglog(tf[1:], valid_off_res_y[1:],
               label = 'off resonance noise', lw=4)
    plt.loglog(tf[1:], valid_on_res_y[1:],
               label = 'on resonance noise', lw=4)
    plt.loglog(tf[1:], valid_sig_y[1:],
               label = 'on signal', lw=4)
    plt.xscale('log')
    plt.title(('Sensitivity = {:E} nV,'.format(sens)+
               'Time constant = {} s, \nFilter slope = {} dB/oct,'.format(
                       tc, filt)+
                       ' Photons = {}\n'.format(photons)+
                       '$\Phi_{{ext}}$ = {} $\Phi_0$ and $n_g$ = {}e'.
                       format(
                       bias[0], bias[1])), fontsize=20, 
    y=0.65)
    plt.legend(fontsize = 20)
    filename = directory_name_fft + '/Y_Quadrature_Setup_FFT'
    plt.savefig(filename+'.jpeg')
    plt.close()
    
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
    ax1.loglog(tf[1:], valid_sig_ng[1:],
               label = 'Feedback off', lw=4)
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
    ax3.loglog(tf[1:], valid_sig_y[1:],
               label = 'Feedback off Y', lw=4)
    ax3.loglog(tf[1:], valid_feed_y_dc[1:],
               label = 'Feedback on Y', lw=4)
    ax3.loglog(tf[1:], valid_off_res_y_sig[1:],
               label = 'Off resonance', lw=4) 
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
    
    """PLOT 19"""
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
    plt.plot(sig_avg_off_res[0][start:], sig_avg_off_res[2][start:],
             label = 'Off resonance', lw =4)
    plt.plot(sig_avg_on_res[0][start:], sig_avg_on_res[2][start:],
             label = 'Feedback off', lw =4)
    plt.plot(sig_avg[0][start:], sig_avg[2][start:],
             label = 'Beginning noise', lw =4)
    plt.plot(sig_avg_feed_dc[0][start:], sig_avg_feed_dc[2][start:],
             label = 'Feedback on', lw =4)
    plt.legend(fontsize = 20, loc = 'best')
    plt.suptitle(('Carrier power = {} dBm, Modulation frequency = {} MHz\n'.
                  format(
            carrier_power, mod_freq*1e-6)+
    'Resonance frequency = {} GHz, Photons = {}\n'.format(
            round(res_freq_meas*1e-9,5),
                           photons)+
    '$\Phi_{{ext}}$ = {} $\Phi_0$ and $n_g$ = {} e\n'.format(
                       bias[0], bias[1])+
            'PID, UL and LL are {}'.format(PID_prompt)), fontsize=20, y = 0.9)
    filename=(directory_name_fft+'/output_noise_feedback_dc.jpeg')
    plt.savefig(filename)
    
    np.savetxt(directory_name_fft +'/fft.dat', np.column_stack((
            tf, valid_off_res_ng, valid_off_res_y,
            valid_on_res_ng, valid_on_res_y,
            valid_sig_ng, valid_sig_y)))
    
    return np.round(ng_dc_feed_mean,3)