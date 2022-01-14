# -*- coding: utf-8 -*-
"""
Created on March 7 2021

@author: Sisira

Measurement: Characterize resonance using frequency fluctuations model.

1. Load bias point for bg.
2. Do a resonance scan.
3. Load actual bias point and repeat scan.
3. Background subtract and apply frequency fluctuations model.

Manually set:
    Smoothing parameter for VNA.
"""
import numpy as np #Scientific computing module
import matplotlib.pyplot as plt #Plotting module
import functions as fn
import pickle as pkl
import lmfit
import cmath
import os
import time
from scipy import special, interpolate

def s11_sigma_rot(f,f0,kint,kext,sigma,phi):
    radius = kext/(kint+kext)
    ktot = kint+kext
    s11 = 1-radius*np.sqrt(np.pi/2)*ktot/sigma*special.wofz(
            (ktot*1j-2*(f-f0))/(2*np.sqrt(2)*sigma))
    x = s11.real
    y = s11.imag
    x_1 = (x-1)*np.cos(phi) + y*np.sin(phi) + 1
    y_1 = -(x-1)*np.sin(phi) + y*np.cos(phi)
    return x_1 + y_1*1j
    
class ResonatorModel(lmfit.model.Model):
    def __init__(self, *args, **kwargs):
        super().__init__(s11_sigma_rot, *args, **kwargs)
        
def characterize_resonance(directory_name, parameters, scan_no, inst_call):
    directory_name_res_scan = (directory_name + '/Res_characterize_scan_{}'.
                               format(scan_no))
    os.mkdir(directory_name_res_scan)
    line = fn.lineno
    [bias, bias_point_bg, temp,
     vna_power, smooth_param, el_delay, center_freq,
     flux_fn, gate_fn, gate_source, gate_voltage_divider] = parameters
    
    [na, sg, modgen, gate_meter, flux_meter, daq, awg, switch_in,
     switch_sa, bp_ctrl, twpa_pump, sa, lockin, dig] = inst_call

    "Parameter file"

    par_array=[]
    line0 = line()
    par_array.append([line()-line0,"Bias point? ",bias])
    par_array.append([line()-line0,"Temperature(mK)",temp])
    par_array.append([line()-line0,"Channel number",1])
    par_array.append([line()-line0,"Number of traces",2])
    par_array.append([line()-line0,"Power? (dBm)",vna_power])
    par_array.append([line()-line0,"Electrical delay?",el_delay])
    par_array.append([line()-line0,"VNA smooting aperture",smooth_param])
    par_array.append([line()-line0,"Number of sweep points (refined)",201])
    par_array.append([line()-line0,"Averaging number (refined)",50])
    par_array.append([line()-line0,"Span (refined) (Hz)",30e6])
    par_array.append([line()-line0,"IF Bandwidth (refined)",1e3])
    par_array.append([line()-line0,"Center frequency", center_freq])

    filename_parameter=str(directory_name_res_scan)+'/Parameters.dat'
    np.savetxt(filename_parameter,np.row_stack((par_array)),delimiter=';',
               fmt='%s')
    scan_parameters=np.loadtxt(filename_parameter,delimiter=';',
                               usecols=2,dtype=np.str) 
    #loads the parameter values

    "Parameter values"

    scan_parameters = np.array([eval(i) for i in scan_parameters])
    [bias, temp, 
     channel, num_traces, power, el_delay, smooth_param,
     points, av, span, IF_bw, center_freq] = scan_parameters
 
    channel=int(channel)
    num_traces=int(num_traces)
    points=int(points)
    
    resonator = ResonatorModel()
    
    with open('interpolated_kappa_int_values.pkl', 'rb') as handle: 
        #fn of ng,phi
        kappa_int_fn = pkl.load(handle)
    with open('interpolated_kappa_ext_values.pkl', 'rb') as handle: 
        #fn of ng,phi
        kappa_ext_fn = pkl.load(handle)

    """------------------------------MEASUREMENT ---------------------------"""
    fn.vna_connect(switch_in)
    for tr in [1,2]:
        na.set_smoothing(smooth_param, channel = tr)
        if smooth_param == 0:
            na.toggle_smoothing(0, trace = tr)
        else:
            na.toggle_smoothing(1, trace = tr)
    
    bias_bg = bias_point_bg[0]
    flux = flux_fn([bias_bg[0]])[0]
    fn.bias_source_set_volt('daq', flux, 0, daq = daq)    
    gate = gate_fn([bias_bg[1]])[0]
    fn.bias_source_set_volt(gate_source, gate, 1, 
                            volt_div = gate_voltage_divider, daq = daq)
    time.sleep(2.0)
    
    filename_str_logmag=directory_name_res_scan+ '/logmag_bg'
    filename_str_phase=directory_name_res_scan+ '/phase_bg'
    (par_array, freq, logmag_bg, phase_bg,
     res_amp_pos, res_amp, center_freq_ref) = fn.vna_scan_save(
             center_freq, span, power, points, av, 
             el_delay, IF_bw,
             filename_str_logmag, filename_str_phase, 
             save = True, plot = True, bg = False)
    
    f_logmag = interpolate.interp1d(freq, logmag_bg, 'cubic')    
    f_phase = interpolate.interp1d(freq, phase_bg, 'cubic')
    
    flux = flux_fn([bias[0]])[0]
    fn.bias_source_set_volt('daq', flux, 0, daq = daq)    
    gate = gate_fn([bias[1]])[0]
    fn.bias_source_set_volt(gate_source, gate, 1, 
                            volt_div = gate_voltage_divider, daq = daq)
    time.sleep(2.0)
    
    filename_str_logmag=directory_name_res_scan+ '/logmag'
    filename_str_phase=directory_name_res_scan+ '/phase'
    (par_array, freq, logmag, phase,
     res_amp_pos, res_amp, center_freq_ref, 
     bg_reduced_data_mag, bg_reduced_data_phase) = fn.vna_scan_save(
             center_freq, span, power, points, av, 
             el_delay, IF_bw,
             filename_str_logmag,filename_str_phase, 
             save = True, plot = True, bg = True,
             f_mag = f_logmag, f_phase = f_phase)
    mod = 10**(bg_reduced_data_mag/20)
    theta = bg_reduced_data_phase*np.pi/180
    xy = [cmath.rect(mod[k], theta[k]) for k in range(len(mod))]
    xy = np.asarray(xy)
    f = freq
    
    true_params = resonator.make_params(kint=0.1,kext=1.23,
                           sigma = 0.01, f0 = center_freq_ref*1e-6,
                           phi = 0.0)
    true_params['kint'].set(max=5,min=0.01)
    true_params['kext'].set(min=0.1,max=5)
    true_params['sigma'].set(min=0,max=5)
    true_params['f0'].set(min=min(f*1e-6),max=max(f*1e-6))
    true_params['phi'].set(min=-0.2, max=0.2)
#        true_params['f0'].set(vary=False)
    result = resonator.fit(xy,f=f*1e-6, params=true_params)
    fit_s11 = resonator.eval(params=result.params, f=f*1e-6)
    kint = result.best_values['kint']
    kext = result.best_values['kext']
    f0 = result.best_values['f0']
    sigma = result.best_values['sigma']
    phi = result.best_values['phi']
    
    fig, ((ax1,ax_fit),(ax5,ax3)) =plt.subplots(2,2, figsize=(20,20))
    ax1.plot(freq*1e-9, mod, lw=3)
    ax1.set_xlabel('Frequency (GHz)',fontsize = 30)
    ax1.set_ylabel('|S11|',fontsize = 30)
    ax1.spines['bottom'].set_linewidth(6)
    ax1.spines['top'].set_linewidth(6)
    ax1.spines['left'].set_linewidth(6)
    ax1.spines['right'].set_linewidth(6)
    ax1.tick_params(axis='x',direction='out', labelsize=30,width=2,length=10)
    ax1.tick_params(axis='y',direction='out', labelsize=30,width=2,length=10)
        
    ax3.plot(freq*1e-9, bg_reduced_data_phase)
    ax3.set_xlabel('Frequency (GHz)',fontsize = 30)
    ax3.set_ylabel('Phase (deg)',fontsize = 30)
    ax3.spines['bottom'].set_linewidth(6)
    ax3.spines['top'].set_linewidth(6)
    ax3.spines['left'].set_linewidth(6)
    ax3.spines['right'].set_linewidth(6)
    ax3.tick_params(axis='x',direction='out', labelsize=30,width=2,length=10)
    ax3.tick_params(axis='y',direction='out', labelsize=30,width=2,length=10)
        
        
    ax5.spines['bottom'].set_linewidth(6)
    ax5.spines['top'].set_linewidth(6)
    ax5.spines['left'].set_linewidth(6)
    ax5.spines['right'].set_linewidth(6)
    ax5.tick_params(axis='x',direction='out', labelsize=30,width=2,length=10)
    ax5.tick_params(axis='y',direction='out', labelsize=30,width=2,length=10)
    ax5.scatter(xy.real, xy.imag, s=40)
    ax5.plot(fit_s11.real, fit_s11.imag, color='r', lw=3)
    ax5.axis('equal')
    ax5.set_title('${S}_{11}$ in complex plane',fontsize=30,y=0)
        
    ax_fit.spines['bottom'].set_linewidth(6)
    ax_fit.spines['top'].set_linewidth(6)
    ax_fit.spines['left'].set_linewidth(6)
    ax_fit.spines['right'].set_linewidth(6)
    string = ('$f_0$ = {} GHz\n $\phi$ = {} \n'.format(
            round(f0*1e-3,6), round(phi,2))+
    '$\kappa_{{int}}$ = {} MHz\n $\kappa_{{ext}}$ = {} MHz \n'.format(
            round(kint,3), round(kext,3))+ 
            r'$\sigma$ = {} MHz'.format(round(sigma,2))+'\n'+
            r'$\chi_{{\nu}}^2$ ={}'.format(round(result.redchi,5)))
    print(f0)
    if (f0*1e-3 < 5.8083) and (f0*1e-3 > 5.699):
        add_string = ('\n Published $\kappa_{{int}}$ = {} MHz'.
        format(round(kappa_int_fn(f0*1e6)*1e-6,2))+
            '\n Published $\kappa_{{ext}}$ = {} MHz'.format(
                    round(kappa_ext_fn(f0*1e6)*1e-6,2)))
    else:
        add_string = '\n $f_0$ not in published data.'
    ax_fit.text(0.1,0.3, string+add_string, fontsize=30)
    ax_fit.set_xticks([])
    ax_fit.set_yticks([])
    fig.subplots_adjust(wspace=0.4)
    fig.suptitle('$\Phi_{{ext}}$ = {} $\Phi_0$ and $n_g$ = {} e'.format(
            bias[0], bias[1]), fontsize=30,y=0.91)
    plt.savefig(directory_name+'/freq_fluc_fit_{}_dBm_{}.jpeg'.format(round(
            power,2), scan_no))
    plt.close()
    return (center_freq, round(f0*1e-3,6)*1e9, center_freq_ref,
            round(kint,3)*1e6, 
            round(kext,3)*1e6, f_logmag, f_phase)
    
def setup_resonance_scan(vna_power_bias, vna_av_bias, el_delay,
                         bias, bias_point_bg, temp, smooth_param, 
                         directory_name, flux_fn, gate_fn, gate_source,
                         gate_voltage_divider, scan_no,
                         inst_call):
    res_check = fn.vna_scan_save(5.75e9,300e6, vna_power_bias, 961, 
                                 vna_av_bias,
                                 el_delay, 10e3)
    res_freq = res_check[6]
    vna_scan_parameters = [bias, bias_point_bg, temp,
                           vna_power_bias, smooth_param, el_delay, 
                           res_freq,
                           flux_fn, gate_fn, gate_source, gate_voltage_divider]
    (res_freq, res_freq_meas, res_freq_ref,
     kappa_int, kappa_ext,
     f_logmag, f_phase) = characterize_resonance(directory_name, 
                                                    vna_scan_parameters,
                                                    scan_no=scan_no,
                                                    inst_call=inst_call)
    if (abs(res_freq-res_freq_meas) > 5e6):
        raise ValueError(
                'Measured resonance frequency too far from initial estimate.')
    return (res_freq, res_freq_meas, res_freq_ref, kappa_int, kappa_ext,
            f_logmag, f_phase)

def get_twpa_resonance(directory_name, parameters, scan_no, inst_call):
    directory_name_res_scan = (directory_name + '/Res_characterize_scan_{}'.
                               format(scan_no))
    os.mkdir(directory_name_res_scan)
    line = fn.lineno
    [bias_point, temp,
     vna_power, smooth_param, el_delay, center_freq,
     bias_point_bg, gate_source, gate_voltage_divider,
     flux_fn, gate_fn] = parameters
     
    [na, sg, modgen, gate_meter, flux_meter, daq, awg, switch_in,
             switch_sa, bp_ctrl, twpa_pump, sa, lockin, dig] = inst_call

    "Parameter file"

    par_array=[]
    line0 = line()
    par_array.append([line()-line0,"Bias point? ",bias_point])
    par_array.append([line()-line0,"Temperature(mK)",temp])
    par_array.append([line()-line0,"Channel number",1])
    par_array.append([line()-line0,"Number of traces",2])
    par_array.append([line()-line0,"Power? (dBm)",vna_power])
    par_array.append([line()-line0,"Electrical delay?",el_delay])
    par_array.append([line()-line0,"VNA smooting aperture",smooth_param])
    par_array.append([line()-line0,"Number of sweep points (refined)",201])
    par_array.append([line()-line0,"Averaging number (refined)",30])
    par_array.append([line()-line0,"Span (refined) (Hz)",25e6])
    par_array.append([line()-line0,"IF Bandwidth (refined)",1e3])
    par_array.append([line()-line0,"Center frequency", center_freq])

    filename_parameter=str(directory_name_res_scan)+'/Parameters.dat'
    np.savetxt(filename_parameter,np.row_stack((par_array)),delimiter=';',
               fmt='%s')
    scan_parameters=np.loadtxt(filename_parameter,delimiter=';',
                               usecols=2,dtype=np.str) 
    #loads the parameter values

    "Parameter values"

    scan_parameters = np.array([eval(i) for i in scan_parameters])
    [bias_point, temp, 
     channel, num_traces, power, el_delay, smooth_param,
     points, av, span, IF_bw, center_freq] = scan_parameters
 
    channel=int(channel)
    num_traces=int(num_traces)
    points=int(points)
    
    """--------------------------------MEASUREMENT -------------------------"""
    fn.vna_connect(switch_in)
    
    bias_bg = bias_point_bg[0]
    bias = bias_point[0]
    flux = flux_fn([bias_bg[0]])[0]
    fn.bias_source_set_volt('daq', flux, 0, daq = daq)    
    gate = gate_fn([bias_bg[1]])[0]
    fn.bias_source_set_volt(gate_source, gate, 1, 
                            volt_div = gate_voltage_divider, daq = daq)
    time.sleep(2.0)
    
    filename_str_logmag=directory_name_res_scan+ '/logmag_bg'
    filename_str_phase=directory_name_res_scan+ '/phase_bg'
    (par_array, freq, logmag_bg, phase_bg,
     res_amp_pos, res_amp, center_freq_ref) = fn.vna_scan_save(
             center_freq, span, power, points, av, 
             el_delay, IF_bw,
             filename_str_logmag, filename_str_phase, 
             save = True, plot = True, bg = False)
    
    f_logmag = interpolate.interp1d(freq, logmag_bg, 'cubic')    
    f_phase = interpolate.interp1d(freq, phase_bg, 'cubic')
    print([bias[0]])
    flux = flux_fn([bias[0]])[0]
    fn.bias_source_set_volt('daq', flux, 0, daq = daq)    
    gate = gate_fn([bias[1]])[0]
    fn.bias_source_set_volt(gate_source, gate, 1, 
                            volt_div = gate_voltage_divider, daq = daq)
    time.sleep(2.0)
    
    filename_str_logmag=directory_name_res_scan+ '/logmag'
    filename_str_phase=directory_name_res_scan+ '/phase'
    (par_array, freq, logmag, phase,
     res_amp_pos, res_amp, center_freq_ref,
     bg_reduced_data_mag, bg_reduced_data_phase) = fn.vna_scan_save(
             center_freq, span, power, points, av, 
             el_delay, IF_bw,
             filename_str_logmag,filename_str_phase, 
             save = True, plot = True, bg =True,
             f_mag = f_logmag, f_phase = f_phase)
    
    fig, (ax1,ax3) =plt.subplots(2,1, figsize=(20,20))
    ax1.plot(freq*1e-9, logmag, lw=3)
    ax1.set_xlabel('Frequency (GHz)',fontsize = 30)
    ax1.set_ylabel('|S11|',fontsize = 30)
    ax1.spines['bottom'].set_linewidth(6)
    ax1.spines['top'].set_linewidth(6)
    ax1.spines['left'].set_linewidth(6)
    ax1.spines['right'].set_linewidth(6)
    ax1.tick_params(axis='x',direction='out', labelsize=30,width=2,length=10)
    ax1.tick_params(axis='y',direction='out', labelsize=30,width=2,length=10)
        
    ax3.plot(freq*1e-9, phase, lw=3)
    ax3.set_xlabel('Frequency (GHz)',fontsize = 30)
    ax3.set_ylabel('Phase (deg)',fontsize = 30)
    ax3.spines['bottom'].set_linewidth(6)
    ax3.spines['top'].set_linewidth(6)
    ax3.spines['left'].set_linewidth(6)
    ax3.spines['right'].set_linewidth(6)
    ax3.tick_params(axis='x',direction='out', labelsize=30,width=2,length=10)
    ax3.tick_params(axis='y',direction='out', labelsize=30,width=2,length=10)
    
    fig.suptitle('$\Phi_{{ext}}$ = {} $\Phi_0$ and $n_g$ = {} e'.format(
            bias[0], bias[1]), fontsize=30,y=0.91)
    plt.savefig(directory_name_res_scan+'/resonance_twpa_on.jpeg')
    plt.close()
    return (center_freq, center_freq_ref)
        