# -*- coding: utf-8 -*-
"""
Sweeps power around estimated carrier power and marks down the peak power for
an input signal at the required resonance frequency.
Lower value corresponds to the required carrier power.
"""
import numpy as np #Scientific computing module
import os #system functions
import time
import matplotlib.pyplot as plt #Plotting module
 #Scientific analysis module
import functions as fn
from scipy import stats
    
def scan(directory_name, parameters, inst_call):
    directory_name_sweep = directory_name + '/carrier_power_calibration'
    os.mkdir(directory_name_sweep)
    directory_name_data = directory_name_sweep + '/data'
    os.mkdir(directory_name_data)
    line = fn.lineno
    [photons, bias, temp, mod_freq, mod_volt, res_freq_meas,
     bp_volt_fn, power_array, carrier_power_est] = parameters
    power_array_str = '{}'.format([k for k in power_array])
     
    "Parameter file"

    par_array=[]
    line0 = line()
    par_array.append([line()-line0,"Modulation voltage (Vpp)",mod_volt])
    par_array.append([line()-line0,"Spectrum analyzer span",10e3])
    par_array.append([line()-line0,"Resolution bandwidth",1])
    par_array.append([line()-line0,"video bandwidth",1])
    par_array.append([line()-line0,"averaging",10])
    par_array.append([line()-line0,"Reference level (dBm)", -40])
    par_array.append([line()-line0,"Carrier power str", power_array_str])
    
    filename_parameter=str(directory_name_sweep)+'/Parameters.dat'
    np.savetxt(filename_parameter,np.row_stack((par_array)),delimiter=';',
               fmt='%s')
    scan_parameters=np.loadtxt(filename_parameter,delimiter=';',
                               usecols=2,dtype=np.str) 
    #loads the parameter values

    "Parameter values"

    scan_parameters = np.array([eval(i) for i in scan_parameters])
    [mod_volt,
     sa_span, res_bw,
     video_bw, avg, ref_level, power_array] = scan_parameters
     
    [na, sg, modgen, gate_meter, flux_meter, daq, awg, switch_in,
     switch_sa, bp_ctrl, twpa_pump, sa, lockin, dig] = inst_call
    
    sg.set_frequency(res_freq_meas)
    bp_ctrl.set_voltage(bp_volt_fn(res_freq_meas))
    fn.spec_connect_out(switch_in, switch_sa)
    modgen.set_state_Vpp(0,1,mod_volt,mod_volt)
    modgen.set_state_freq(0,1,mod_freq*1e-6,mod_freq*1e-6)
    modgen.set_state_phase(0,1, 0.0, 0.0)
    sg.toggle_output(1)
    time.sleep(1)
    
    """--------------------------------MEASUREMENT -------------------------"""
    output_peak_power =[]
    for power in power_array:
        "Output signal measurement"
        sg.set_power(power) 
        filename_data=(directory_name_data+'/fridge_out_center_{}_dbm.dat'
                       .format(power))
        (par_array, freq, power_data) = fn.spec_scan(res_freq_meas,
        sa_span, res_bw, video_bw, avg, ref_level, filename_data)
        peak_power = power_data[int(sa_span/2)]+carrier_power_est-power
        output_peak_power.append(peak_power)
        print(power, peak_power)
        
    output_peak_power = np.asarray(output_peak_power)
    pos = np.where(output_peak_power == np.min(output_peak_power))[0]
    carrier_power = power_array[pos[0]]
    
    fig, ax = plt.subplots(figsize=(20,10))
    ax.spines['bottom'].set_linewidth(5)
    ax.spines['top'].set_linewidth(5)
    ax.spines['left'].set_linewidth(5)
    ax.spines['right'].set_linewidth(5)
    ax.tick_params(axis='x',direction='out', labelsize=40,width=2,length=10)
    ax.tick_params(axis='y',direction='out', labelsize=40,width=2,length=10)
    plt.xlabel('Carrier power (dBm)',fontsize = 30)
    plt.ylabel('Spectrum analyzer power at $\omega_c$ (dBm)',fontsize = 30)
    plt.scatter(power_array, output_peak_power, s =40)
    plt.plot(power_array, output_peak_power)
    plt.title(('$\Phi_{{ext}}$ = {} $\Phi_0$ and $n_g$ = {} e\n'.format(
                       bias[0], bias[1])+
        'Carrier power for n = {} is {} dBm'.format(
                photons, carrier_power)), fontsize=30, 
    y=0.86)
    filename = directory_name+'/SA_power_response.jpeg'
    plt.savefig(filename)
    plt.close()
    
    filename = directory_name+'/power_data.dat'
    np.savetxt(filename,np.column_stack((power_array, 
                                         output_peak_power)))
    
    modgen.set_state_Vpp(0,1,0,0)
    sg.toggle_output(0)

    return carrier_power
    
    