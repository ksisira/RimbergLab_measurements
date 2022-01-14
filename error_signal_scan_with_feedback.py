# -*- coding: utf-8 -*-
"""
Acquires error signal for varying gate near ng_(0) with or without feedback.
Takes vna scans when configured without feedback to extract resonance.
"""
import numpy as np #Scientific computing module
import os #system functions
import time
import matplotlib.pyplot as plt #Plotting module
 #Scientific analysis module
import functions as fn
import pickle as pkl
from scipy import interpolate

def scan(directory_name, parameters, inst_call, feed_str = ''):
    directory_name_sweep = directory_name + '/Lockin_sweep_ng' + feed_str
    os.mkdir(directory_name_sweep)
    directory_name_vna = directory_name_sweep + '/vna'
    if feed_str=='':
        os.mkdir(directory_name_vna)
    directory_name_lockin = directory_name_sweep + '/lockin'
    os.mkdir(directory_name_lockin)
    line = fn.lineno
    [photons, bias, temp, mod_freq, mod_volt, res_freq_meas,
     sens, tc, phase_ref, expand, filt, 
     clock_rate, sample_rate, acq_no, acq_time,
     bp_volt_fn, ng_stop, ng_start, sweep_points,
     vna_power_dBm, el_delay, smooth_param,
     gate_source_no, gate_voltage_divider, twpa_state] = parameters
     
    [na, sg, modgen, gate_meter, flux_meter, daq, awg, switch_in,
             switch_sa, bp_ctrl, twpa_pump, sa, lockin, dig] = inst_call
     
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
    par_array.append([line()-line0,"Gate sweep start (e)",ng_start])
    par_array.append([line()-line0,"Gate sweep stop (e)",ng_stop])
    par_array.append([line()-line0,"Carrier sweep points",sweep_points])
    par_array.append([line()-line0,"Channel number",1])
    par_array.append([line()-line0,"Number of traces",2])
    par_array.append([line()-line0,"Power? (dBm)",vna_power_dBm])
    par_array.append([line()-line0,"Electrical delay?",el_delay])
    par_array.append([line()-line0,"VNA smooting aperture",smooth_param])
    par_array.append([line()-line0,"Number of sweep points (refined)",601])
    par_array.append([line()-line0,"Averaging number (refined)",30])
    par_array.append([line()-line0,"Span (refined) (Hz)",100e6])
    par_array.append([line()-line0,"IF Bandwidth (refined)",1e3])

    filename_parameter=str(directory_name_sweep)+'/Parameters.dat'
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
     ng_start, ng_stop, sweep_points, channel, num_traces, vna_power_dBm,
     el_delay, smooth_param, vna_points, vna_avg, vna_span,
     IF_bw] = scan_parameters
     
    impedance = int(impedance)
    acq_no = int(acq_no)
    sweep_points = int(sweep_points)
    channel=int(channel)
    num_traces=int(num_traces)
    vna_points=int(vna_points)
    
    if gate_source_no == 1:
        gate_source = 'awg'
        awg.toggle_output(1,1)
    elif gate_source_no == 0:
        gate_source = 'daq'
    with open(directory_name +
          '/interpolated_gate_bias_conversion.pkl', 'rb') as handle:
        gate_fn = pkl.load(handle)
    ng_e_array = np.linspace(ng_start, ng_stop, sweep_points+1)
    
    bg_data = np.loadtxt(directory_name+'/bg_freq_logmag_phase.dat')
    f_logmag = interpolate.interp1d(bg_data[:,0], bg_data[:,1], 'cubic')    
    f_phase = interpolate.interp1d(bg_data[:,0],bg_data[:,2], 'cubic')
    
    """--------------------------------MEASUREMENT -------------------------"""
    carrier_freq_set = res_freq_meas
    sg.set_frequency(carrier_freq_set)
    bp_ctrl.set_voltage(bp_volt_fn(carrier_freq_set))
    gate_array=[]
    freq_array=[]
    ng_array=[]
    y_array=[]
    ng_array_err=[]
    y_array_err=[]
    points_omit = int(sample_rate*0.1)
    for j in range(sweep_points+1):
        print(j, ng_e_array[j])
        gate = gate_fn([ng_e_array[j]])[0]
        fn.bias_source_set_volt(gate_source, gate, 1, 
                                volt_div = gate_voltage_divider, daq = daq)
        time.sleep(2.0)
        gate_reading=gate_meter.get_voltage()
        gate_array.append(gate_reading)
        print('Gate at sample {} mV'.format(gate_reading*1e3))
        
        if feed_str == '':
            """VNA Scan"""
            fn.vna_connect(switch_in)
            if twpa_state:
                twpa_pump.toggle_output(0)
            filename_str_logmag=(str(directory_name_vna)+
                                 '/gate_sweep_logmag_broad_'+str(j)+
                                 '_th_gate')
            filename_str_phase=(str(directory_name_vna)+
                                '/gate_sweep_phase_broad_'+str(j)+
                                '_th_gate')
            (par_array, freq, logmag, phase,
             res_amp_pos, res_amp, center_freq_ref, 
             bg_reduced_data_mag, bg_reduced_data_phase) = fn.vna_scan_save(
                     res_freq_meas, vna_span, vna_power_dBm, 
                     vna_points, vna_avg, 
                     el_delay, IF_bw,
                     filename_str_logmag, filename_str_phase, 
                     save = True, plot = True, bg = True,
                     f_mag = f_logmag, f_phase = f_phase)
            freq_array.append(center_freq_ref)
        
        "Error signal measurement"
        if twpa_state:
            twpa_pump.toggle_output(1)
        time.sleep(1)
        fn.error_sig_measure(switch_in, switch_sa)
        modgen.set_state_Vpp(0,1,mod_volt,mod_volt)
        sg.toggle_output(1)
        time.sleep(1)
        
        filename=(directory_name_lockin+'/quadrature_data_'+str(j)+
                             '_th_gate')
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
        data_dig = (acquire_data[2], stitched_data[0], stitched_data[1])
        np.savetxt(filename+'.dat', np.column_stack((data_dig[0], 
                                                     data_dig[1], 
                                                     data_dig[2])),
                   encoding = 'utf-16')
        
        if feed_str == '':
            if j==0:
                plt.figure()
#                plt.plot(data_dig[0][:-points_omit],
#                         data_dig[1][:-points_omit], label = 'X')
                plt.plot(data_dig[0][:-points_omit],
                         data_dig[2][:-points_omit], label = 'Y')
                plt.legend()
                plt.savefig(directory_name_sweep+'/sample_data.jpeg')
                plt.close()
        else:
            plt.figure()
#            plt.plot(data_dig[0][:-points_omit],
#                     data_dig[1][:-points_omit], label = 'X')
            plt.plot(data_dig[0][:-points_omit],
                     data_dig[2][:-points_omit], label = 'Y')
            plt.legend()
            plt.savefig(directory_name_lockin+'/feedback_Y_{}_th_gate.jpeg'.
                        format(j))
            plt.close()
            
            plt.figure()
#            plt.plot(data_dig[0][:-points_omit],
#                     data_dig[1][:-points_omit], label = 'X')
            plt.plot(data_dig[0][:-points_omit],
                     data_dig[1][:-points_omit]/gate_voltage_divider,
                     label = '$n_g$')
            plt.legend()
            plt.savefig(directory_name_lockin+'/feedback_ng_{}_th_gate.jpeg'.
                        format(j))
            plt.close()
            
        ng_mean = np.mean(data_dig[1][:-points_omit])/gate_voltage_divider
        y_mean = np.mean(data_dig[2][:-points_omit])
        ng_err = np.std(data_dig[1][:-points_omit])/gate_voltage_divider
        y_err = np.std(data_dig[2][:-points_omit])
        ng_array.append(ng_mean)
        y_array.append(y_mean)
        ng_array_err.append(ng_err)
        y_array_err.append(y_err)
        
    gate_array = np.asarray(gate_array)
    ng_array = np.asarray(ng_array)
    y_array = np.asarray(y_array)
    ng_array_err = np.asarray(ng_array_err)
    y_array_err = np.asarray(y_array_err)
    
    filename = directory_name+'/error_signal_data{}.dat'.format(feed_str)
    if feed_str == '':
        freq_array = np.asarray(freq_array)
        np.savetxt(filename,np.column_stack((ng_e_array,ng_array,y_array,
                                             ng_array_err, y_array_err,
                                             gate_array, freq_array)))
    else:
        np.savetxt(filename,np.column_stack((ng_e_array,ng_array,y_array,
                                             ng_array_err, y_array_err,
                                             gate_array)))
    
    fig, ax = plt.subplots(figsize=(20,10))
    ax.spines['bottom'].set_linewidth(5)
    ax.spines['top'].set_linewidth(5)
    ax.spines['left'].set_linewidth(5)
    ax.spines['right'].set_linewidth(5)
    ax.tick_params(axis='x',direction='out', labelsize=40,width=2,length=10)
    ax.tick_params(axis='y',direction='out', labelsize=40,width=2,length=10)
    plt.xlabel('$n_g$ - $n_g^{{(0)}}$ (e)',fontsize = 30)
    plt.ylabel('Corrected $n_g$ (mV)',fontsize = 30)
    plt.scatter(ng_e_array - bias[1], ng_array*1e3, s =40)
    plt.errorbar(ng_e_array - bias[1],ng_array*1e3, yerr = ng_array_err*1e3, 
                 ms=40)
    plt.title(('Sensitivity = {:E} nV,'.format(sens)+
               'Time constant = {} s\n Filter slope = {} dB/oct,'.format(
                       tc, filt)+
                       'Photons = {}\n'.format(photons)+
                       '$\Phi_{{ext}}$ = {} $\Phi_0$ and $n_g^{{(0)}}$ = {} e'.
                       format(
                       bias[0], bias[1])), fontsize=30, 
    y=0.86)
    filename = directory_name+'/Error_signal_ng_sweep_ng{}.jpeg'.format(
            feed_str)
    plt.savefig(filename)
    
    fig, ax = plt.subplots(figsize=(20,10))
    ax.spines['bottom'].set_linewidth(5)
    ax.spines['top'].set_linewidth(5)
    ax.spines['left'].set_linewidth(5)
    ax.spines['right'].set_linewidth(5)
    ax.tick_params(axis='x',direction='out', labelsize=40,width=2,length=10)
    ax.tick_params(axis='y',direction='out', labelsize=40,width=2,length=10)
    plt.xlabel('$n_g$-$n_g^{{(0)}}$ (e)',fontsize = 30)
    plt.ylabel('Lock-in output Y (mV)',fontsize = 30)
    plt.scatter(ng_e_array-bias[1], y_array*1e3, s =40)
    plt.errorbar(ng_e_array - bias[1], y_array*1e3, yerr = y_array_err*1e3, 
                 ms=40)
    plt.title(('Sensitivity = {:E} nV,'.format(sens)+
               'Time constant = {} s\n Filter slope = {} dB/oct,'.format(
                       tc, filt)+
                       'Photons = {}\n'.format(photons)+
                       '$\Phi_{{ext}}$ = {} $\Phi_0$ and $n_g^{{(0)}}$ = {} e'.
                       format(
                       bias[0], bias[1])), fontsize=30, 
    y=0.86)
    filename = directory_name+'/Error_signal_ng_sweep_Y{}.jpeg'.format(
            feed_str)
    plt.savefig(filename)
    
    if feed_str=='':
        fig, ax = plt.subplots(figsize=(20,10))
        ax.spines['bottom'].set_linewidth(5)
        ax.spines['top'].set_linewidth(5)
        ax.spines['left'].set_linewidth(5)
        ax.spines['right'].set_linewidth(5)
        ax.tick_params(axis='x',direction='out', labelsize=40,width=2,length=10)
        ax.tick_params(axis='y',direction='out', labelsize=40,width=2,length=10)
        plt.xlabel('$\omega_c$ - $\omega_0$ (MHz)',fontsize = 30)
        plt.ylabel('Corrected $n_g$ (mV)',fontsize = 30)
        plt.scatter((freq_array-res_freq_meas)*1e-6, ng_array*1e3, s =40)
        plt.title(('Sensitivity = {:E} nV,'.format(sens)+
                   'Time constant = {} s\n Filter slope = {} dB/oct,'.format(
                           tc, filt)+
                       'Photons = {}\n'.format(photons)+
                       '$\Phi_{{ext}}$ = {} $\Phi_0$ and $n_g^{{(0)}}$ = {} e'.
                       format(
                       bias[0], bias[1])), fontsize=30, 
        y=0.86)
        filename = (directory_name+'/Error_signal_ng_sweep_f0vsng{}.jpeg'.
                    format(feed_str))
        plt.savefig(filename)
    
        fig, ax = plt.subplots(figsize=(20,10))
        ax.spines['bottom'].set_linewidth(5)
        ax.spines['top'].set_linewidth(5)
        ax.spines['left'].set_linewidth(5)
        ax.spines['right'].set_linewidth(5)
        ax.tick_params(axis='x',direction='out', labelsize=40,width=2,length=10)
        ax.tick_params(axis='y',direction='out', labelsize=40,width=2,length=10)
        plt.xlabel('$\omega_c$ - $\omega_0$ (MHz)',fontsize = 30)
        plt.ylabel('Lock-in output Y (mV)',fontsize = 30)
        plt.scatter((freq_array-res_freq_meas)*1e-6, y_array*1e3, s =40)
        plt.title(('Sensitivity = {:E} nV,'.format(sens)+
                   'Time constant = {} s\n Filter slope = {} dB/oct,'.format(
                           tc, filt)+
                           'Photons = {}\n'.format(photons)+
                           '$\Phi_{{ext}}$ = {} $\Phi_0$ and $n_g^{{(0)}}$ = {} e'.
                           format(
                           bias[0], bias[1])), fontsize=30, 
        y=0.86)
        filename = directory_name+'/Error_signal_ng_sweep_f0vsY{}.jpeg'.format(
                feed_str)
        plt.savefig(filename)
    
        fig, ax = plt.subplots(figsize=(20,10))
        ax.spines['bottom'].set_linewidth(5)
        ax.spines['top'].set_linewidth(5)
        ax.spines['left'].set_linewidth(5)
        ax.spines['right'].set_linewidth(5)
        ax.tick_params(axis='x',direction='out', labelsize=40,width=2,length=10)
        ax.tick_params(axis='y',direction='out', labelsize=40,width=2,length=10)
        plt.xlabel('$n_g$ (e)',fontsize = 30)
        plt.ylabel('$\omega_0$ (GHz)',fontsize = 30)
        plt.scatter(ng_e_array, freq_array*1e-9, s = 40)
        plt.plot(ng_e_array, freq_array*1e-9, lw = 3, color = 'r')
        plt.title(('Sensitivity = {:E} nV,'.format(sens)+
                   'Time constant = {} s\n Filter slope = {} dB/oct,'.format(
                           tc, filt)+
                           'Photons = {}\n'.format(photons)+
                           '$\Phi_{{ext}}$ = {} $\Phi_0$ and $n_g^{{(0)}}$ = {} e'.
                           format(
                           bias[0], bias[1])), fontsize=30, 
        y=0.86)
        filename = directory_name+'/ng_vs_f0{}.jpeg'.format(feed_str)
        plt.savefig(filename)
    
    
    modgen.set_state_Vpp(0,1,0,0)
    sg.toggle_output(0)
    gate = gate_fn([bias[1]])[0]
    fn.bias_source_set_volt(gate_source, gate, 1, 
                            volt_div = gate_voltage_divider,
                            daq = daq)
    time.sleep(2.0)
    if feed_str=='':
        return (ng_e_array, y_array, gate_array, freq_array)
    else:
        return (ng_e_array, y_array, gate_array)