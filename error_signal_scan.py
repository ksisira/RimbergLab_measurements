# -*- coding: utf-8 -*-
"""
Acquires error signal for varying carrier frequency.
Spectrum analyzer measures the power spectrum output and plots the peak
value.
Digitizer acquires both the quadrature.
If phase_ref zero, output the slope of the parameteric plot.
"""
import numpy as np #Scientific computing module
import os #system functions
import time
import matplotlib.pyplot as plt #Plotting module
 #Scientific analysis module
import functions as fn
from scipy import stats
    
def scan(directory_name, parameters, inst_call):
    directory_name_sweep = directory_name + '/Lockin_sweep'
    os.mkdir(directory_name_sweep)
    directory_name_sa = directory_name_sweep + '/sa'
    os.mkdir(directory_name_sa)
    directory_name_lockin = directory_name_sweep + '/lockin'
    os.mkdir(directory_name_lockin)
    line = fn.lineno
    [photons, bias, temp, mod_freq, mod_volt, res_freq_meas,
     sens, tc, phase_ref, expand, filt, 
     clock_rate, sample_rate, acq_no, acq_time,
     bp_volt_fn, carrier_step, sweep_points] = parameters
     
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
    par_array.append([line()-line0,"Digitizer input range", 4])
    par_array.append([line()-line0,"Digitizer ipedance", 1e6])
    par_array.append([line()-line0,"Carrier step (in Hz)",carrier_step])
    par_array.append([line()-line0,"Carrier sweep points",sweep_points])
    par_array.append([line()-line0,"Spectrum analyzer span",10e3])
    par_array.append([line()-line0,"Resolution bandwidth",1])
    par_array.append([line()-line0,"video bandwidth",1])
    par_array.append([line()-line0,"averaging",2])
    par_array.append([line()-line0,"Reference level (dBm)", -40])
    
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
     carrier_step, sweep_points, sa_span, res_bw,
     video_bw, avg, ref_level] = scan_parameters
     
    [na, sg, modgen, gate_meter, flux_meter, daq, awg, switch_in,
     switch_sa, bp_ctrl, twpa_pump, sa, lockin, dig] = inst_call
     
    impedance = int(impedance)
    acq_no = int(acq_no)
    sweep_points = int(sweep_points)
    
    plt.figure()
    
    """--------------------------------MEASUREMENT -------------------------"""
    carrier_freq_set = res_freq_meas - (
            sweep_points/2)*carrier_step - carrier_step
    freq_array=[]
    x_array=[]
    y_array=[]
    x_array_err=[]
    y_array_err=[]
    output_peak_power =[]
    points_omit = int(sample_rate*0.1)

    for i in range(sweep_points+1):
        print(i)
        carrier_freq_set = carrier_freq_set + carrier_step
        print(round(carrier_freq_set*1e-9,6))
        freq_array.append(carrier_freq_set)
        sg.set_frequency(carrier_freq_set)
        bp_ctrl.set_voltage(bp_volt_fn(carrier_freq_set))
        
        "Output signal measurement"
        fn.spec_connect_out(switch_in, switch_sa)
        modgen.set_state_Vpp(0,1,mod_volt,mod_volt)
        sg.toggle_output(1)
        time.sleep(1)
        
        filename_data=(directory_name_sa+'/fridge_out_broad_{}_GHz.dat'
                       .format(round(carrier_freq_set*1e-9,6)))
        (par_array, freq, power) = fn.spec_scan(carrier_freq_set,
        5*mod_freq,10e3, 10e3, 2, ref_level, filename_data)
        plt.plot(power)
        
        filename_data=(directory_name_sa+'/fridge_out_left_{}_GHz.dat'
                       .format(round(carrier_freq_set*1e-9,6)))
        (par_array, freq, power) = fn.spec_scan(carrier_freq_set-mod_freq,
        sa_span,res_bw, video_bw, avg, ref_level, filename_data)
        filename_data=(directory_name_sa+'/fridge_out_right_{}_GHz.dat'
                       .format(round(carrier_freq_set*1e-9,6)))
        (par_array, freq, power) = fn.spec_scan(carrier_freq_set+mod_freq,
        sa_span,res_bw, video_bw, avg, ref_level, filename_data)
        filename_data=(directory_name_sa+'/fridge_out_center_{}_GHz.dat'
                       .format(round(carrier_freq_set*1e-9,6)))
        (par_array, freq, power) = fn.spec_scan(carrier_freq_set,
        sa_span,res_bw, video_bw, avg, ref_level, filename_data)
        peak_power = power[int(sa_span/2)]
        output_peak_power.append(peak_power)
        
        "Error signal measurement"
        fn.error_sig_measure(switch_in, switch_sa)
        modgen.set_state_Vpp(0,1,mod_volt,mod_volt)
        sg.toggle_output(1)
        time.sleep(1)
        
        filename=(directory_name_lockin+'/quadrature_data_{}_GHz'
                       .format(round(carrier_freq_set*1e-9,6)))
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
        
        if i==0:
            plt.figure()
            plt.plot(data_dig[0][:-points_omit],
                     data_dig[1][:-points_omit], label = 'X')
            plt.plot(data_dig[0][:-points_omit],
                     data_dig[2][:-points_omit], label = 'Y')
            plt.legend()
            plt.savefig(directory_name_sweep+'/sample_data.jpeg')
            plt.close()
            
        x_mean = np.mean(data_dig[1][:-points_omit])
        y_mean = np.mean(data_dig[2][:-points_omit])
        x_err = np.std(data_dig[1][:-points_omit])
        y_err = np.std(data_dig[2][:-points_omit])
        x_array.append(x_mean)
        y_array.append(y_mean)
        x_array_err.append(x_err)
        y_array_err.append(y_err)
        
    plt.savefig(directory_name_lockin+'/sa_output_broad.jpeg')
    plt.close()
        
    freq_array = np.asarray(freq_array)
    x_array = np.asarray(x_array)
    y_array = np.asarray(y_array)
    x_array_err = np.asarray(x_array_err)
    y_array_err = np.asarray(y_array_err)
    output_peak_power = np.asarray(output_peak_power)
    
    fig, ax = plt.subplots(figsize=(20,10))
    ax.spines['bottom'].set_linewidth(5)
    ax.spines['top'].set_linewidth(5)
    ax.spines['left'].set_linewidth(5)
    ax.spines['right'].set_linewidth(5)
    ax.tick_params(axis='x',direction='out', labelsize=40,width=2,length=10)
    ax.tick_params(axis='y',direction='out', labelsize=40,width=2,length=10)
    plt.xlabel('$\omega_c$ - $\omega_0$ (MHz)',fontsize = 30)
    plt.ylabel('Lock-in output X (mV)',fontsize = 30)
    plt.scatter((freq_array-res_freq_meas)*1e-6, x_array*1e3, s =40)
    plt.errorbar((freq_array-res_freq_meas)*1e-6,x_array*1e3, 
                 yerr = x_array_err*1e3, 
                 ms=40)
    plt.title(('Sensitivity = {:E} nV,'.format(sens)+
               'Time constant = {} s\n Filter slope = {} dB/oct,'.format(
                       tc, filt)+
                       'Photons = {}\n'.format(photons)+
                       '$\Phi_{{ext}}$ = {} $\Phi_0$ and $n_g$ = {} e'.format(
                       bias[0], bias[1])), fontsize=30, 
    y=0.86)
    filename = directory_name+'/Error_signal_sweep_X.jpeg'
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
    plt.errorbar((freq_array-res_freq_meas)*1e-6,y_array*1e3, 
                 yerr = y_array_err*1e3, 
                 ms=40)
    plt.axvline(x=0, ymin = np.min(y_array*1e3), ymax = np.max(y_array*1e3),
                color = 'r')
    plt.title(('Sensitivity = {:E} nV,'.format(sens)+
               'Time constant = {} s\n Filter slope = {} dB/oct,'.format(
                       tc, filt)+
                       'Photons = {}\n'.format(photons)+
                       '$\Phi_{{ext}}$ = {} $\Phi_0$ and $n_g$ = {} e'.format(
                       bias[0], bias[1])), fontsize=30, 
    y=0.86)
    filename = directory_name+'/Error_signal_sweep_Y.jpeg'
    plt.savefig(filename)
    plt.close()
    
    fig, ax = plt.subplots(figsize=(20,10))
    ax.spines['bottom'].set_linewidth(5)
    ax.spines['top'].set_linewidth(5)
    ax.spines['left'].set_linewidth(5)
    ax.spines['right'].set_linewidth(5)
    ax.tick_params(axis='x',direction='out', labelsize=40,width=2,length=10)
    ax.tick_params(axis='y',direction='out', labelsize=40,width=2,length=10)
    plt.xlabel('$\omega_c$ - $\omega_0$ (MHz)',fontsize = 30)
    plt.ylabel('Spectrum analyzer power at $\omega_c$ (dBm)',fontsize = 30)
    plt.scatter((freq_array-res_freq_meas)*1e-6, output_peak_power, s =40)
    plt.plot((freq_array-res_freq_meas)*1e-6, output_peak_power)
    plt.title(('Sensitivity = {:E} nV,'.format(sens)+
               'Time constant = {} s\n Filter slope = {} dB/oct,'.format(
                       tc, filt)+
                       'Photons = {}\n'.format(photons)+
                       '$\Phi_{{ext}}$ = {} $\Phi_0$ and $n_g$ = {} e'.format(
                       bias[0], bias[1])), fontsize=30, 
    y=0.86)
    filename = directory_name+'/Error_signal_sweep_SA.jpeg'
    plt.savefig(filename)
    plt.close()
    
    fig, ax = plt.subplots(figsize=(20,10))
    ax.spines['bottom'].set_linewidth(5)
    ax.spines['top'].set_linewidth(5)
    ax.spines['left'].set_linewidth(5)
    ax.spines['right'].set_linewidth(5)
    ax.tick_params(axis='x',direction='out', labelsize=40,width=2,length=10)
    ax.tick_params(axis='y',direction='out', labelsize=40,width=2,length=10)
    plt.xlabel('X (V)',fontsize = 30)
    plt.ylabel('Y (V)',fontsize = 30)
    plt.scatter(x_array, y_array, s = 40)
    linear_fit_points = 10
    mid_point = int(sweep_points/2)
    slope, intercept, r_value, p_value, std_err = stats.linregress(
            x_array[mid_point-linear_fit_points:mid_point+
                    linear_fit_points+1], 
            y_array[mid_point-linear_fit_points:mid_point+
                    linear_fit_points+1])
    angle_fit = slope*x_array[mid_point-linear_fit_points:mid_point+
                              linear_fit_points+1]+intercept
    if phase_ref == 0:
        plt.plot(x_array[mid_point-linear_fit_points:mid_point+
                         linear_fit_points+1], angle_fit, lw=3)
    plt.title(('Sensitivity = {:E} nV,'.format(sens)+
               'Time constant = {} s\n Filter slope = {} dB/oct,'.format(
                       tc, filt)+
                       'Photons = {}\n'.format(photons)+
                       '$\Phi_{{ext}}$ = {} $\Phi_0$ and $n_g$ = {} e'.format(
                       bias[0], bias[1])), fontsize=30, 
    y=0.86)
    filename = directory_name+'/XY.jpeg'
    plt.savefig(filename)
    plt.close()
    
    filename = directory_name+'/error_signal_data.dat'
    np.savetxt(filename,np.column_stack((freq_array,x_array,y_array,
                                         x_array_err, y_array_err, 
                                         output_peak_power)))
    
    modgen.set_state_Vpp(0,1,0,0)
    sg.toggle_output(0)
    return (slope, intercept, np.arctan(slope))
    
    