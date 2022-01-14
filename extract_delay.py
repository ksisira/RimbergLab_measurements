# -*- coding: utf-8 -*-
"""
Measure off-resonance fluctuations to get the offset, get the fluctuations
off-resonance, on-resonance noise and signal.
"""
import numpy as np #Scientific computing module
import os #system functions
import time
import matplotlib.pyplot as plt #Plotting module
from scipy.fft import fft, fftfreq
 #Scientific analysis module
import functions as fn
import pickle as pkl

def extract_delay(directory_name, parameters, inst_call):
    directory_name_off = directory_name + '/extract_delay'
    os.mkdir(directory_name_off)
    line = fn.lineno
    [photons, bias, temp, mod_volt,
     sens, tc, phase_ref, expand, filt, 
     clock_rate, sample_rate, acq_no, acq_time, res_freq_meas,
     carrier_power, mod_freq, gate_voltage_divider,manual_dng,
     PID_VD] = parameters
     
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
    par_array.append([line()-line0,"Manual dng", manual_dng])
    par_array.append([line()-line0,"Additional voltage divider", PID_VD])
    
    filename_parameter=str(directory_name_off)+'/Parameters.dat'
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
     manual_dng, PID_VD] = scan_parameters
     
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
    
    """--------------------------------MEASUREMENT -------------------------"""
    fn.error_sig_measure(switch_in, switch_sa)
    
    flux = flux_fn([bias[0]+0.5])[0]
    fn.bias_source_set_volt('daq', flux, 0, daq = daq)  
    lockin.set_sensitivity_nV(sens)
    lockin.set_time_constant_sec(tc)
    lockin.set_reference_phase_deg(phase_ref)
    lockin.set_filter_slope(filt)
    
    lockin.set_expand(1, 1)
    lockin.set_expand(2, 1)
    lockin.set_offset('X',0)
    lockin.set_offset('Y',0)
    time.sleep(3)
    
    points_omit = int(sample_rate*0.1)
    dig.configure_board(input_range_V,input_range_V, ref_source = 'int',
                            trigger_timeout_sec = 10e-6,
                            impedance_1 = impedance, impedance_2 = impedance,
                            coupling_1 = 'dc', coupling_2 = 'dc',
                            sample_rate = clock_rate)
    samples = int((clock_rate*acq_time)//32)*32
    acquire_data = dig.acquire(3, 0, samples, 1, 1,
                               input_range_V=input_range_V,
                               clock_rate=clock_rate, 
                               sample_rate=sample_rate)
    stitched_data=dig.stitch_data(acquire_data[0],acquire_data[1],1)
    data_off_res_0_off = (acquire_data[2], stitched_data[0], stitched_data[1])
    filename=(directory_name_off+'/quadrature_noise_off_resonance_0_offset')
    np.savetxt(filename+'.dat', np.column_stack((data_off_res_0_off[0],
                                                 data_off_res_0_off[1],
                                                 data_off_res_0_off[2])),
    encoding = 'utf-16')
    
    
    ng_noise_mean = np.mean(data_off_res_0_off[1]
    [:-points_omit])/gate_voltage_divider
    y_noise_mean = np.mean(data_off_res_0_off[2][:-points_omit])
    ng_noise_err = np.std(data_off_res_0_off[1]
    [:-points_omit])/gate_voltage_divider
    y_noise_err = np.std(data_off_res_0_off[1][:-points_omit])
    offset_y = y_noise_mean/10*100
    print(y_noise_mean, offset_y)
    lockin.set_offset('Y', -offset_y)
#    print('changed')
    lockin.set_expand(1, expand)
    lockin.set_expand(2, expand)
    time.sleep(3)
    
    acquire_data = dig.acquire(3, 0, samples, 1, 1,
                               input_range_V=input_range_V,
                               clock_rate=clock_rate, 
                               sample_rate=sample_rate)
    stitched_data=dig.stitch_data(acquire_data[0],acquire_data[1],1)
    data_off_res_1_off = (acquire_data[2], stitched_data[0], stitched_data[1])
    filename=(directory_name_off+'/quadrature_noise_off_resonance_1_offset')
    np.savetxt(filename+'.dat', np.column_stack((data_off_res_1_off[0],
                                                 data_off_res_1_off[1],
                                                 data_off_res_1_off[2])),
    encoding = 'utf-16')
    
    
    flux = flux_fn([bias[0]])[0]
    fn.bias_source_set_volt('daq', flux, 0, daq = daq)
    time.sleep(3)
    
    acquire_data = dig.acquire(3, 0, samples, 1, 1,
                               input_range_V=input_range_V,
                               clock_rate=clock_rate, 
                               sample_rate=sample_rate)
    stitched_data=dig.stitch_data(acquire_data[0],acquire_data[1],1)
    data_on_res_0_sig = (acquire_data[2], stitched_data[0], stitched_data[1])
    filename=(directory_name_off+'/quadrature_noise_on_resonance')
    np.savetxt(filename+'.dat', np.column_stack((data_on_res_0_sig[0],
                                                 data_on_res_0_sig[1],
                                                 data_on_res_0_sig[2])),
    encoding = 'utf-16')
    
    sg.set_mode()
    sg.set_frequency(res_freq_meas)
    sg.set_power(carrier_power)
    sg.toggle_output(1)
    modgen.set_state_freq(0,1,mod_freq*1e-6,mod_freq*1e-6)
    modgen.set_state_Vpp(0,1,mod_volt,mod_volt)
    modgen.set_state_phase(0,1,0.0,0.0)
    time.sleep(3)
    
    acquire_data = dig.acquire(3, 0, samples, 1, 1,
                               input_range_V=input_range_V,
                               clock_rate=clock_rate, 
                               sample_rate=sample_rate)
    stitched_data=dig.stitch_data(acquire_data[0],acquire_data[1],1)
    data_on_res_1_sig = (acquire_data[2], stitched_data[0], stitched_data[1])
    filename=(directory_name_off+'/quadrature_signal_on_resonance')
    np.savetxt(filename+'.dat', np.column_stack((data_on_res_1_sig[0],
                                                 data_on_res_1_sig[1],
                                                 data_on_res_1_sig[2])),
    encoding = 'utf-16')
    
    """dng addition"""
    manual_dV = round(gate_fn(manual_dng)*gate_voltage_divider*PID_VD,3)
    prompt = input('Set the PID manual to {} V. Switch on PID connection (y)'.
                   format(manual_dV))
    if prompt!='y':
        print('DELTA NG NOT ACTIVE!!!!!!!')
    acquire_data = dig.acquire(3, 0, samples, 1, 1,
                               input_range_V=input_range_V,
                               clock_rate=clock_rate, 
                               sample_rate=sample_rate)
    stitched_data=dig.stitch_data(acquire_data[0],acquire_data[1],1)
    data_on_res_1_dng = (acquire_data[2], stitched_data[0], stitched_data[1])
    filename=(directory_name_off+'/quadrature_signal_delta_ng')
    np.savetxt(filename+'.dat', np.column_stack((data_on_res_1_dng[0],
                                                 data_on_res_1_dng[1],
                                                 data_on_res_1_dng[2])),
    encoding = 'utf-16')
    
    """dng switching"""
    prompt = input('Get ready to switch off dng while acquisition.')
    if prompt!='y':
        print('DELTA NG SWITCH DID NOT OCCUR!!!!!!!')
    acquire_data = dig.acquire(3, 0, samples, 1, 1,
                               input_range_V=input_range_V,
                               clock_rate=clock_rate, 
                               sample_rate=sample_rate)
    stitched_data=dig.stitch_data(acquire_data[0],acquire_data[1],1)
    data_on_res_1_sw = (acquire_data[2], stitched_data[0], stitched_data[1])
    filename=(directory_name_off+'/quadrature_signal_switch')
    np.savetxt(filename+'.dat', np.column_stack((data_on_res_1_sw[0],
                                                 data_on_res_1_sw[1],
                                                 data_on_res_1_sw[2])),
    encoding = 'utf-16')
    
    plt.figure()
    plt.plot(data_off_res_0_off[0][:-points_omit],
             data_off_res_0_off[2][:-points_omit], label = 'y_raw')
    plt.plot(data_off_res_1_off[0][:-points_omit],
             data_off_res_1_off[2][:-points_omit], label = 'y_offset')
    plt.title('Mean ng and Y = {}, {} mV\n and STD = {}, {} mV'.format(
            round(ng_noise_mean*1e3,2), round(y_noise_mean*1e3,2),
            round(ng_noise_err*1e3,2), 
            round(y_noise_err*1e3),2))
    plt.legend()
    filename=(directory_name_off+'/quadrature_offset')
    plt.savefig(filename+'.jpeg')
    plt.close()
    
    fig, ax = plt.subplots(figsize=(20,10))
    ax.spines['bottom'].set_linewidth(5)
    ax.spines['top'].set_linewidth(5)
    ax.spines['left'].set_linewidth(5)
    ax.spines['right'].set_linewidth(5)
    ax.tick_params(axis='x',direction='out', labelsize=40,width=2,length=10)
    ax.tick_params(axis='y',direction='out', labelsize=40,width=2,length=10)
    plt.xlabel('Time (sec)',fontsize = 30)
    plt.ylabel('ng (mV)',fontsize = 30)
    plt.plot(data_off_res_1_off[0][:-points_omit],
             data_off_res_1_off[1][:-points_omit]/gate_voltage_divider*1e3, 
             label = 'off resonance noise', lw=4)
    plt.plot(data_on_res_0_sig[0][:-points_omit],
             data_on_res_0_sig[1][:-points_omit]/gate_voltage_divider*1e3, 
             label = 'on resonance noise', lw=4)
    plt.plot(data_on_res_1_sig[0][:-points_omit],
             data_on_res_1_sig[1][:-points_omit]/gate_voltage_divider*1e3, 
             label = 'on resonance signal', lw=4)
    plt.xlim(0, tc*10)
    plt.title(('Sensitivity = {:E} nV,'.format(sens)+
               'Time constant = {} s, \nFilter slope = {} dB/oct,'.format(
                       tc, filt)+
                       ' Photons = {}\n'.format(photons)+
                       '$\Phi_{{ext}}$ = {} $\Phi_0$ and $n_g$ = {} e'.format(
                       bias[0], bias[1])), fontsize=30, 
    y=0.65)
    plt.legend(fontsize = 20)
    filename = directory_name + '/ng_Quadrature_Setup'
    plt.savefig(filename+'.jpeg')
    filename = directory_name_off + '/ng_Quadrature_Setup'
    plt.savefig(filename+'.jpeg')
    plt.close()
    
    fig, ax = plt.subplots(figsize=(20,10))
    ax.spines['bottom'].set_linewidth(5)
    ax.spines['top'].set_linewidth(5)
    ax.spines['left'].set_linewidth(5)
    ax.spines['right'].set_linewidth(5)
    ax.tick_params(axis='x',direction='out', labelsize=40,width=2,length=10)
    ax.tick_params(axis='y',direction='out', labelsize=40,width=2,length=10)
    plt.xlabel('Time (sec)',fontsize = 30)
    plt.ylabel('ng (mV)',fontsize = 30)
    plt.plot(data_on_res_1_sig[0][:-points_omit],
             data_on_res_1_sig[1][:-points_omit]/gate_voltage_divider*1e3, 
             label = 'on resonance signal', lw=4)
    plt.plot(data_on_res_1_dng[0][:-points_omit],
             data_on_res_1_dng[1][:-points_omit]/gate_voltage_divider*1e3, 
             label = 'dng signal', lw=4)
    plt.plot(data_on_res_1_sw[0][:-points_omit],
             data_on_res_1_sw[1][:-points_omit]/gate_voltage_divider*1e3, 
             label = 'switched signal', lw=4)
    plt.xlim(0, tc*10)
    plt.title(('Sensitivity = {:E} nV,'.format(sens)+
               'Time constant = {} s, \nFilter slope = {} dB/oct,'.format(
                       tc, filt)+
                       ' Photons = {}\n'.format(photons)+
                       '$\Phi_{{ext}}$ = {} $\Phi_0$ and $n_g$ = {} e'.format(
                       bias[0], bias[1])), fontsize=30, 
    y=0.65)
    plt.legend(fontsize = 20)
    filename = directory_name + '/ng_Quadrature_switch'
    plt.savefig(filename+'.jpeg')
    filename = directory_name_off + '/ng_Quadrature_switch'
    plt.savefig(filename+'.jpeg')
    plt.close()
    
    fig, ax = plt.subplots(figsize=(20,10))
    ax.spines['bottom'].set_linewidth(5)
    ax.spines['top'].set_linewidth(5)
    ax.spines['left'].set_linewidth(5)
    ax.spines['right'].set_linewidth(5)
    ax.tick_params(axis='x',direction='out', labelsize=40,width=2,length=10)
    ax.tick_params(axis='y',direction='out', labelsize=40,width=2,length=10)
    plt.xlabel('Time (sec)',fontsize = 30)
    plt.ylabel('Lock-in output Y (V)',fontsize = 30)
    plt.plot(data_off_res_1_off[0][:-points_omit],
             data_off_res_1_off[2][:-points_omit], 
             label = 'off resonance noise', lw=4)
    plt.plot(data_on_res_0_sig[0][:-points_omit],
             data_on_res_0_sig[2][:-points_omit], 
             label = 'on resonance noise', lw=4)
    plt.plot(data_on_res_1_sig[0][:-points_omit],
             data_on_res_1_sig[2][:-points_omit], 
             label = 'on resonance signal', lw=4)
    plt.xlim(0, tc*10)
    plt.title(('Sensitivity = {:E} nV,'.format(sens)+
               'Time constant = {} s, \nFilter slope = {} dB/oct,'.format(
                       tc, filt)+
                       ' Photons = {}\n'.format(photons)+
                       '$\Phi_{{ext}}$ = {} $\Phi_0$ and $n_g$ = {} e'.format(
                       bias[0], bias[1])), fontsize=30, 
    y=0.65)
    plt.legend(fontsize = 20)
    filename = directory_name + '/Y_Quadrature_Setup'
    plt.savefig(filename+'.jpeg')
    filename = directory_name_off + '/Y_Quadrature_Setup'
    plt.savefig(filename+'.jpeg')
    plt.close()
    
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
             label = 'on resonance signal', lw=4)
    plt.plot(data_on_res_1_dng[0][:-points_omit],
             data_on_res_1_dng[2][:-points_omit], 
             label = 'dng signal', lw=4)
    plt.plot(data_on_res_1_sw[0][:-points_omit],
             data_on_res_1_sw[2][:-points_omit], 
             label = 'switched signal', lw=4)
    plt.xlim(0, tc*10)
    plt.title(('Sensitivity = {:E} nV,'.format(sens)+
               'Time constant = {} s, \nFilter slope = {} dB/oct,'.format(
                       tc, filt)+
                       ' Photons = {}\n'.format(photons)+
                       '$\Phi_{{ext}}$ = {} $\Phi_0$ and $n_g$ = {} e'.format(
                       bias[0], bias[1])), fontsize=30, 
    y=0.65)
    plt.legend(fontsize = 20)
    filename = directory_name + '/Y_Quadrature_switch'
    plt.savefig(filename+'.jpeg')
    filename = directory_name_off + '/Y_Quadrature_switch'
    plt.savefig(filename+'.jpeg')
    plt.close()
    
    time_array = data_off_res_1_off[0][:-points_omit]
    N = len(time_array)
    T = 1/sample_rate
    tf = fftfreq(N, T)[:N//2]
    
    ngf_off_res_noise = fft(data_off_res_1_off[1]
    [:-points_omit]/gate_voltage_divider)
    yf_off_res_noise = fft(data_off_res_1_off[2][:-points_omit])
    ngf_on_res_noise = fft(data_on_res_0_sig[1]
    [:-points_omit]/gate_voltage_divider)
    yf_on_res_noise = fft(data_on_res_0_sig[2][:-points_omit])
    ngf_on_res_sig = fft(data_on_res_1_sig[1]
    [:-points_omit]/gate_voltage_divider)
    yf_on_res_sig = fft(data_on_res_1_sig[2][:-points_omit])
    valid_off_res_ng = 2/N*np.abs(ngf_off_res_noise[:N//2])
    valid_on_res_ng = 2/N*np.abs(ngf_on_res_noise[:N//2])
    valid_sig_ng = 2/N*np.abs(ngf_on_res_sig[:N//2])
    valid_off_res_y = 2/N*np.abs(yf_off_res_noise[:N//2])
    valid_on_res_y = 2/N*np.abs(yf_on_res_noise[:N//2])
    valid_sig_y = 2/N*np.abs(yf_on_res_sig[:N//2])
    
    
    fig, ax = plt.subplots(figsize=(20,10))
    ax.spines['bottom'].set_linewidth(5)
    ax.spines['top'].set_linewidth(5)
    ax.spines['left'].set_linewidth(5)
    ax.spines['right'].set_linewidth(5)
    ax.tick_params(axis='x',direction='out', labelsize=40,width=2,length=10)
    ax.tick_params(axis='y',direction='out', labelsize=40,width=2,length=10)
    plt.xlabel('Frequency (Hz)',fontsize = 30)
    plt.plot(tf[4:], valid_off_res_ng[4:],
             label = 'off resonance noise', lw=4)
    plt.plot(tf[4:], valid_on_res_ng[4:],
             label = 'on resonance noise', lw=4)
    plt.plot(tf[4:], valid_sig_ng[4:],
             label = 'on signal', lw=4)
    plt.xscale('log')
    plt.yscale('log')
    plt.title(('Sensitivity = {:E} nV,'.format(sens)+
               'Time constant = {} s, \nFilter slope = {} dB/oct,'.format(
                       tc, filt)+
                       ' Photons = {}\n'.format(photons)+
                       '$\Phi_{{ext}}$ = {} $\Phi_0$ and $n_g$ = {} e'.format(
                       bias[0], bias[1])), fontsize=30, 
    y=0.65)
    plt.legend(fontsize = 20)
#    ylim_x = np.max([np.max(valid_off_res_x), np.max(valid_on_res_x),
#                   np.max(valid_sig_x)])
#    plt.ylim(top = ylim_x)
    filename = directory_name + '/ng_Quadrature_Setup_FFT'
    plt.savefig(filename+'.jpeg')
    filename = directory_name_off + '/ng_Quadrature_Setup_FFT'
    plt.savefig(filename+'.jpeg')
    plt.close()
    
    fig, ax = plt.subplots(figsize=(20,10))
    ax.spines['bottom'].set_linewidth(5)
    ax.spines['top'].set_linewidth(5)
    ax.spines['left'].set_linewidth(5)
    ax.spines['right'].set_linewidth(5)
    ax.tick_params(axis='x',direction='out', labelsize=40,width=2,length=10)
    ax.tick_params(axis='y',direction='out', labelsize=40,width=2,length=10)
    plt.xlabel('Frequency (Hz)',fontsize = 30)
    plt.plot(tf[4:], valid_off_res_y[4:],
             label = 'off resonance noise', lw=4)
    plt.plot(tf[4:], valid_on_res_y[4:],
             label = 'on resonance noise', lw=4)
    plt.plot(tf[4:], valid_sig_y[4:],
             label = 'on signal', lw=4)
    plt.xscale('log')
    plt.yscale('log')
    plt.title(('Sensitivity = {:E} nV,'.format(sens)+
               'Time constant = {} s, \nFilter slope = {} dB/oct,'.format(
                       tc, filt)+
                       ' Photons = {}\n'.format(photons)+
                       '$\Phi_{{ext}}$ = {} $\Phi_0$ and $n_g$ = {} e'.format(
                       bias[0], bias[1])), fontsize=30, 
    y=0.65)
    plt.legend(fontsize = 20)
#    ylim_y = np.max([np.max(valid_off_res_y), np.max(valid_on_res_y),
#                   np.max(valid_sig_y)])
#    plt.ylim(top = ylim_y)
    filename = directory_name + '/Y_Quadrature_Setup_FFT'
    plt.savefig(filename+'.jpeg')
    filename = directory_name_off + '/Y_Quadrature_Setup_FFT'
    plt.savefig(filename+'.jpeg')
    plt.close()
    
    fn.error_sig_measure(switch_in, switch_sa)
    
    np.savetxt(directory_name_off +'/fft.dat', np.column_stack((
            tf, valid_off_res_ng, valid_off_res_y,
            valid_on_res_ng, valid_on_res_y,
            valid_sig_ng, valid_sig_y)))
    
    return (offset_y)