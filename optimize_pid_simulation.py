# -*- coding: utf-8 -*-
"""
Optimizes the PID values from the measured Y or power spectral density.

Methods - 

Powell

L-BFGS-B

TNC

SLSQP
"""
import numpy as np
import matplotlib.pyplot as plt
import os, sys
# import functions as fn
from scipy.misc import derivative
from scipy.optimize import minimize
from scipy.fft import fft, fftfreq, ifft
from scipy import interpolate, stats
from lmfit.models import PowerLawModel, ConstantModel
from scipy.signal import blackman

def avg_data(freq, power_sa, span):
    left = np.flip(power_sa[:int(span/2)])
    right = power_sa[int(span/2)+1:]
    avg_data = (left+right)/2
    signal_freq = freq[int(span/2)]
    freq_left = np.flip(signal_freq-freq[:int(span/2)])
    freq_right = freq[int(span/2)+1:]-signal_freq
    return (freq_left, freq_right, avg_data)

def dbm_to_W(value):
    return 10**(value/10)*1e-3

def fft_from_PSD(density_array):
    return np.sqrt(2*density_array)

def ng_sample_meas_from_PSD(fft_array, f_array, ng0):
    fmax = np.max(f_array)
    fft_freq_power=np.concatenate(([0],f_array,np.flip(-f_array)))
    NT = len(fft_freq_power)
    fft_power_both=np.concatenate(([0],NT/2*fft_array,
                                   np.flip(-NT/2*fft_array)))
    random_array=np.exp(np.random.random(len(fft_power_both))*2*np.pi*1j)
    fft_power_both=fft_power_both*random_array
    ng_values = ifft(fft_power_both)
    dt = 1/fmax/2
    time_array = np.linspace(0, NT*dt, NT, endpoint = False)
    return (np.abs(ng_values), time_array)

def ng_tot_array(P, I, D, omega_ar, t_ar, ng_ar, ng_Y_interp, slope,
                 time_factor = 2):
    l = len(t_ar)
    dt = t_ar[1]-t_ar[0]
    ng_tot = np.array([k for k in ng_ar])
    y_ar = np.zeros(l)
    ng_cor_ar = np.zeros(l)
    int_value = 0
    for i in range(l-time_factor):
        y_ar[i] = (ng_Y_interp([ng_tot[i]])[0])
        if i==0:
            int_value = 0
            der_value = 0
        else:
            int_value = int_value + dt*(y_ar[i-1]+y_ar[i])/2
            der_value = (y_ar[i]-y_ar[i-1])/dt
        out_V = -P*(y_ar[i]+I*int_value+D*der_value)/101/11
        ng_cor_ar[i+time_factor] = (slope*out_V)
        ng_tot[i+time_factor] = (ng_tot[i+time_factor]+
              ng_cor_ar[i+time_factor])%2
    return (ng_tot, ng_cor_ar)

def optimize_pid(directory_name, parameters, folder_name):
    directory_name_opt = directory_name + '/optimize_pid'
    os.mkdir(directory_name_opt)
    [photons, bias, tc, vna_power, input_att,
     kappa_int, kappa_ext, delay] = parameters
     
    "Parameter file"

    par_array=[]
    par_array.append([0,"ng0",bias[1]])
    par_array.append([1,"Time constant", tc])
    par_array.append([2,"VNA input power", vna_power])
    par_array.append([3,"Input attenuation to fridge", input_att])
    par_array.append([4,"Delay", delay])
    par_array.append([5,"Folder name", folder_name])

    filename_parameter=str(directory_name_opt)+'/Parameters.dat'
    np.savetxt(filename_parameter,np.row_stack((par_array)),delimiter=';',
               fmt='%s')
    scan_parameters=np.loadtxt(filename_parameter,delimiter=';',
                               usecols=2,dtype=np.str) 
    #loads the parameter values

    "Parameter values"

    scan_parameters = np.array([eval(i) for i in scan_parameters[:-1]])
    [ng0, tc, vna_power, input_att, delay] = scan_parameters
    
    """----------------------------Analysis---------------------------------"""
    
    out_sig_data = np.loadtxt((directory_name+
                               '/Output_signal/fridge_out_signal.dat'))
    out_base_data = np.loadtxt((directory_name+
                                '/Output_signal/fridge_out_off_resonance.dat'))
    ng_Y_data = np.loadtxt(directory_name + '/error_signal_data.dat')
    meas_Y_fluc_sig = np.loadtxt((directory_name + '/Lockin_offset_low_tc' +
                                  '/quadrature_signal_on_resonance.dat'),
    encoding = 'utf-16')
    
    span = len(out_sig_data)-1
    noise = avg_data(out_sig_data[:,0], out_sig_data[:,1], span)
    noise_base = avg_data(out_base_data[:,0], out_base_data[:,1], span)
    noise_W = 10**(noise[2]/10)*1e-3
    noise_base_W = 10**(noise_base[2]/10)*1e-3
    freq = noise[0]
    Sout_sa = noise_W-noise_base_W
    np.savetxt(directory_name_opt+'/power_spectral_density.dat',
               np.column_stack((freq, Sout_sa)))
    
    fig, ax = plt.subplots(figsize=(20,10))
    ax.spines['bottom'].set_linewidth(5)
    ax.spines['top'].set_linewidth(5)
    ax.spines['left'].set_linewidth(5)
    ax.spines['right'].set_linewidth(5)
    ax.tick_params(axis='x',direction='out', labelsize=40,width=2,length=10)
    ax.tick_params(axis='y',direction='out', labelsize=40,width=2,length=10)
    ax.set_xlabel('Frequency (Hz)',fontsize = 30)
    ax.set_ylabel('Power spectral density (W/Hz)',fontsize = 30)
    ax.set_xscale('log')
    ax.set_yscale('log')
    plt.plot(freq[1:], Sout_sa[1:],lw=4)
    filename=(directory_name_opt + '/power_spectral_density.jpeg')
    plt.savefig(filename)
    
    power_in_sample = vna_power-input_att
    p_out_sa =  out_sig_data[int(span/2),1]
    G = dbm_to_W(p_out_sa)/dbm_to_W(power_in_sample)
    kappa_tot = kappa_int+kappa_ext
    
    Sww = (Sout_sa/G/dbm_to_W(power_in_sample)*kappa_tot**2/
       (2*kappa_ext**2)*(freq**2+kappa_tot**2/4))
    np.savetxt(directory_name_opt+'/freq_fluc_density.dat',
               np.column_stack((freq, Sww)))
    step=3
    dng = ng_Y_data[step,0]-ng_Y_data[0,0]
    g_array = np.gradient(ng_Y_data[::step,6], dng)
    np.savetxt(directory_name_opt+'/g_values.dat',
               np.column_stack((ng_Y_data[::step,0], g_array)))
    
    plt.figure()
    plt.plot(ng_Y_data[::step,0], g_array)
    g_interp = interpolate.interp1d(ng_Y_data[::step,0], g_array)
    g=g_interp(ng0)
    if g == np.NaN:
        g = input('g is invalid. Enter value for g (Hz/e).') or 10e6
    Sqq=Sww/g**2
    for i in range(len(Sqq)):
        limit = 1e-11
        if Sqq[i]<limit:
            Sqq[i] = limit
    mod_1 = PowerLawModel()
    mod_2 = ConstantModel()
    mod = mod_1+mod_2
    params = mod_1.make_params()
    params.add('exponent', value=-1, min=-2, max=-1.3, vary=False)
    params += mod_2.guess(Sqq[1:], x = freq[1:])
    params.add('constant', value = 0)
    stop = -1
    start = 1
    out = mod.fit(Sqq[start:stop], params, x = freq[start:stop])
    np.savetxt(directory_name_opt+'/charge_noise_density.dat',
               np.column_stack((freq[start:stop], Sqq[start:stop], 
                                out.best_fit)))
    
    fig, ax = plt.subplots(figsize=(20,10))
    ax.spines['bottom'].set_linewidth(5)
    ax.spines['top'].set_linewidth(5)
    ax.spines['left'].set_linewidth(5)
    ax.spines['right'].set_linewidth(5)
    ax.tick_params(axis='x',direction='out', labelsize=40,width=2,length=10)
    ax.tick_params(axis='y',direction='out', labelsize=40,width=2,length=10)
    ax.set_xlabel('Frequency (Hz)',fontsize = 30)
    ax.set_ylabel('Charge density ($e^2$/Hz)',fontsize = 30)
    ax.set_xscale('log')
    ax.set_yscale('log')         
    plt.plot(freq[1:], Sqq[1:], lw=4)
    plt.plot(freq[start:stop], out.best_fit, lw=4)
    plt.suptitle(('g = {} MHz/e at $n_g$ = {}e\n'.format(
            np.round(g*1e-6,2), bias[1])+
    'Fit amplitude = {:.3g}, exponent (vary {}) = {:.3g} and constant = {:.3g}'.
    format(out.params['amplitude'].value,
           out.params['exponent'].vary,
           out.params['exponent'].value,
           out.params['constant'].value)), fontsize = 20)
    filename=directory_name_opt+'/charge_noise_density_from_SA.jpeg'
    plt.savefig(filename)
    
    time_factor = 2
    dt = delay/time_factor
    omega_max=int(1/dt/2)
    omega_ar = freq[:omega_max]
    Nf = len(omega_ar)
    Sqq_fit = mod.eval(out.params, x=omega_ar)
    NT = 2*Nf+1
    
    Y_measured = meas_Y_fluc_sig[:,2]
    ng_Y_interp = interpolate.interp1d(np.concatenate((
            ng_Y_data[:,0],2-ng_Y_data[:-1,0])),
        np.concatenate((
                ng_Y_data[:,2],ng_Y_data[:-1,2])),
                'cubic')
    slope_mV_ng, intercept_mV_ng, r, p, se = stats.linregress(ng_Y_data[:,5],
                                                              ng_Y_data[:,0])
        
    """Power fit"""
    fft_power_array = fft_from_PSD(Sqq_fit)
    (ng_noise_power, time_ar)=ng_sample_meas_from_PSD(
        fft_power_array, omega_ar, ng0)
    # time_mid_1_pos = np.where(time_ar_long==0.5)[0][0]
    # time_mid_2_pos = time_mid_1_pos+int(len(time_ar_long-1)*domega_fine)
    # time_ar = time_ar_long[time_mid_1_pos:time_mid_2_pos+1]-0.5
    # ng_noise_power = ng_noise_power_long[time_mid_1_pos:time_mid_2_pos+1]
    ng_sample_array_power = ng_noise_power+ng0
    print(time_ar[-1])

    def ng_tot_power_std(var):
        P, I, D = var
        ng_tot_ar = ng_tot_array(P, I, D, omega_ar, 
                                 time_ar, ng_sample_array_power,
                                 ng_Y_interp, slope_mV_ng,
                                 time_factor=time_factor)[0]
        return np.std(ng_tot_ar)
    
    print(np.std(ng_sample_array_power))
    res_power = minimize(ng_tot_power_std, [0.2,5e3,4.7e-5],
                         bounds = ((-10,10), (0,1e6), (0,1e2)),method = 'TNC')
    #'TNC' and 'L-BFGS-B' performing better with random initial conditions
    ng_res_power = ng_tot_array(res_power['x'][0], res_power['x'][1], 
                                res_power['x'][2], omega_ar, 
                                time_ar, ng_sample_array_power, 
                                ng_Y_interp, slope_mV_ng,
                                time_factor=time_factor)
    np.savetxt(directory_name_opt+'/charge_dynamics_power.dat',
               np.column_stack((time_ar, ng_sample_array_power,
                                ng_res_power[0])))
    
    """Y fit"""
    ng_start = ng0-0.1
    ng_stop = ng0+0.1
    ng_start_pos = np.where(abs(
            ng_Y_data[:,0]-ng_start)==np.min(
                    abs(ng_Y_data[:,0]-ng_start)))[0][0]
    ng_stop_pos = np.where(abs(
            ng_Y_data[:,0]-ng_stop)==np.min(
                    abs(ng_Y_data[:,0]-ng_stop)))[0][0]
    Y_max = np.max(ng_Y_data[ng_start_pos:ng_stop_pos,2])
    Y_min = np.min(ng_Y_data[ng_start_pos:ng_stop_pos,2])
    for Y in range(len(Y_measured)):
        if Y_measured[Y] > Y_max:
            Y_measured[Y] = Y_max
        if Y_measured[Y] < Y_min:
            Y_measured[Y] = Y_min
    Y_ng_interp = interpolate.interp1d(
            ng_Y_data[ng_start_pos:ng_stop_pos+1,2],
            ng_Y_data[ng_start_pos:ng_stop_pos+1,0])
    ng_sample_meas = Y_ng_interp(Y_measured)
    ng_sample_interp = interpolate.interp1d(meas_Y_fluc_sig[:,0],
                                            ng_sample_meas)
    ng_sample_array_Y = ng_sample_interp(time_ar)
    ng_noise_Y = ng_sample_array_Y-ng0
    
    def ng_tot_Y_std(var):
        P, I, D = var
        ng_tot_ar = ng_tot_array(P, I, D, omega_ar, 
                                 time_ar, ng_sample_array_Y,
                                 ng_Y_interp, slope_mV_ng)[0]
        return np.std(ng_tot_ar)
    
    print(np.std(ng_sample_array_Y))
    res_Y = minimize(ng_tot_Y_std, [0.2,5e3,4.7e-5],
                         bounds = ((-10,10), (0,1e6), (0,1e2)),method = 'TNC')
    #'TNC' and 'L-BFGS-B' performing better with random initial conditions
    ng_res_Y = ng_tot_array(res_Y['x'][0], res_Y['x'][1], 
                                res_Y['x'][2], omega_ar, 
                                time_ar, ng_sample_array_Y, 
                                ng_Y_interp, slope_mV_ng)
    np.savetxt(directory_name_opt+'/charge_dynamics_Y.dat',
               np.column_stack((time_ar, ng_sample_array_Y,
                                ng_res_Y[0])))        
    
    fig, ax = plt.subplots(figsize=(20,10))
    ax.spines['bottom'].set_linewidth(5)
    ax.spines['top'].set_linewidth(5)
    ax.spines['left'].set_linewidth(5)
    ax.spines['right'].set_linewidth(5)
    ax.tick_params(axis='x',direction='out', labelsize=40,width=2,length=10)
    ax.tick_params(axis='y',direction='out', labelsize=40,width=2,length=10)
    ax.set_xlabel('Time (s)',fontsize = 30)
    ax.set_ylabel('$n_g$(t)',fontsize = 30)
    plt.plot(time_ar, ng_sample_array_power, lw=4, label = 'Feedback off')
    plt.plot(time_ar, ng_res_power[0], lw=4, label = 'Feedback on')
    plt.legend(fontsize=20)
    plt.suptitle(('$\sigma$ before and after feedback {:e}  and {:e} e\n'.
                  format(
            np.std(ng_sample_array_power),
                 res_power.fun)+
                    'Convergence {}, Mean $n_g$ = {}, Set ng = {}\n'.
    format(res_power.success,np.round(np.mean(ng_res_power[0]),3), ng0)+
    'Optimization done from SA data'), fontsize=20)
    filename=(directory_name_opt+'/sample_charge_dynamics_power_fit.jpeg')
    plt.savefig(filename)
    
    fig, ax = plt.subplots(figsize=(20,10))
    ax.spines['bottom'].set_linewidth(5)
    ax.spines['top'].set_linewidth(5)
    ax.spines['left'].set_linewidth(5)
    ax.spines['right'].set_linewidth(5)
    ax.tick_params(axis='x',direction='out', labelsize=40,width=2,length=10)
    ax.tick_params(axis='y',direction='out', labelsize=40,width=2,length=10)
    ax.set_xlabel('Time (s)',fontsize = 30)
    ax.set_ylabel('$n_g$(t)',fontsize = 30)
    plt.plot(time_ar, ng_sample_array_Y, lw=4, label = 'Feedback off')
    plt.plot(time_ar, ng_res_Y[0], lw=4, label = 'Feedback on')
    plt.legend(fontsize=20)
    plt.suptitle(('$\sigma$ before and after feedback {:e}  and {:e} e\n'.format(
            np.std(ng_sample_array_Y),
                 res_Y.fun)+'Convergence {}, Mean $n_g$ = {}, Set ng = {}\n'.
    format(res_Y.success,np.round(np.mean(ng_res_Y[0]),3), ng0)+
    'Optimization done from digitizer data'), fontsize=20)
    filename=(directory_name_opt+'/sample_charge_dynamics_Y_fit.jpeg')
    plt.savefig(filename)
    
    fig, ax = plt.subplots(figsize=(20,10))
    ax.spines['bottom'].set_linewidth(5)
    ax.spines['top'].set_linewidth(5)
    ax.spines['left'].set_linewidth(5)
    ax.spines['right'].set_linewidth(5)
    ax.tick_params(axis='x',direction='out', labelsize=40,width=2,length=10)
    ax.tick_params(axis='y',direction='out', labelsize=40,width=2,length=10)
    ax.set_xlabel('Frequency (Hz)',fontsize = 30)
    ax.set_ylabel('Charge noise ($e^2$/Hz)',fontsize = 30)
    ax.set_xscale('log')
    ax.set_yscale('log')
    NT = len(time_ar)
    dt=time_ar[1]-time_ar[0]
    xf=fftfreq(NT, dt)
    w=blackman(NT)
    yf_off = fft(ng_noise_power*w)
    yf_on = fft(ng_res_power[0]*w)
    plt.plot(freq[1:], Sqq[1:], lw=4,
             label = 'Data without feedback')
    plt.plot(xf[2:NT//2], (np.abs(yf_off[2:NT//2])*2/NT)**2/2, lw=4,
             label = 'Simulated charge noise without feedback')
    plt.plot(xf[2:NT//2], (np.abs(yf_on[2:NT//2])*2/NT)**2/2,'r',lw=4,
             label = 'Simulated charge noise with feedback')
    plt.legend(fontsize=20)
    plt.suptitle(('P = {:.3g}, I = {:.3g} and D = {:.3g}\n'.format(
            res_power['x'][0],res_power['x'][1],
            res_power['x'][2])+
    'Fit amplitude = {:.3g} and exponent = {:.3g}'.format(
            out.params['amplitude'].value,
            out.params['exponent'].value)), fontsize = 20)
    filename=(directory_name_opt+'/charge_noise_density_power_fit.jpeg')
    plt.savefig(filename)
    
#    np.savetxt(directory_name_opt+'/FFT_ng_feedback_power_fit.dat',
#               np.column_stack((xf, 2/NT*np.abs(yf_off)*np.sqrt(domega),
#                                2/NT*np.abs(yf_on)*np.sqrt(domega))))
    
    fig, ax = plt.subplots(figsize=(20,10))
    ax.spines['bottom'].set_linewidth(5)
    ax.spines['top'].set_linewidth(5)
    ax.spines['left'].set_linewidth(5)
    ax.spines['right'].set_linewidth(5)
    ax.tick_params(axis='x',direction='out', labelsize=40,width=2,length=10)
    ax.tick_params(axis='y',direction='out', labelsize=40,width=2,length=10)
    ax.set_xlabel('Frequency (Hz)',fontsize = 30)
    ax.set_ylabel('Charge noise ($e^2$/Hz)',fontsize = 30)
    ax.set_xscale('log')
    ax.set_yscale('log')
    xf=fftfreq(NT, dt)
    yf_off = fft(ng_noise_Y*w)
    yf_on = fft(ng_res_Y[0]*w)
    plt.plot(xf[2:NT//2], (np.abs(yf_off[2:NT//2])*2/NT)**2/2, lw=4,
             label = 'Simulated charge noise without feedback')
    plt.plot(xf[2:NT//2], (np.abs(yf_on[2:NT//2])*2/NT)**2/2,'r',lw=4,
             label = 'Simulated charge noise with feedback')
    plt.legend(fontsize=20)
    plt.suptitle(('P = {:.3g}, I = {:.3g} and D = {:.3g}\n'.format(
            res_Y['x'][0],res_Y['x'][1],
            res_Y['x'][2])), fontsize = 20)
    filename=(directory_name_opt+'/charge_noise_density_Y_fit.jpeg')
    plt.savefig(filename)
#    np.savetxt(directory_name_opt+'/FFT_ng_feedback_Y_fit.dat',
#               np.column_stack((xf, 2/NT*np.abs(yf_off)*np.sqrt(domega),
#                                2/NT*np.abs(yf_on)*np.sqrt(domega))))
    
    fig, ax = plt.subplots(figsize=(20,10))
    ax.spines['bottom'].set_linewidth(5)
    ax.spines['top'].set_linewidth(5)
    ax.spines['left'].set_linewidth(5)
    ax.spines['right'].set_linewidth(5)
    ax.tick_params(axis='x',direction='out', labelsize=40,width=2,length=10)
    ax.tick_params(axis='y',direction='out', labelsize=40,width=2,length=10)
    ax.set_xlabel('Frequency (Hz)',fontsize = 30)
    ax.set_ylabel('Charge noise ($e^2$/Hz)',fontsize = 30)
    ax.set_xscale('log')
    ax.set_yscale('log')
    plt.plot(freq[1:], Sqq[1:], lw=4,
             label = 'Spectrum analyzer measurement')
    plt.plot(xf[2:NT//2], (np.abs(yf_off[2:NT//2])*2/NT)**2/2, lw=4,
             label = 'Digitizer Y measurement')
    plt.legend(fontsize=20)
    filename=(directory_name_opt+'/compare_data.jpeg')
    plt.savefig(filename)

#    
#    plt.figure()
#    plt.plot(time_ar, ng_res[1])
#    plt.plot(time_ar, ng_sample_array-ng0)
#    
#    cor_f = fft(ng_res[1])
#    
#    fig, ax = plt.subplots(figsize=(20,10))
#    ax.spines['bottom'].set_linewidth(5)
#    ax.spines['top'].set_linewidth(5)
#    ax.spines['left'].set_linewidth(5)
#    ax.spines['right'].set_linewidth(5)
#    ax.tick_params(axis='x',direction='out', labelsize=40,width=2,length=10)
#    ax.tick_params(axis='y',direction='out', labelsize=40,width=2,length=10)
#    ax.set_xlabel('Frequency (Hz)',fontsize = 30)
#    ax.set_ylabel('Corrected charge noise ($e^2$/Hz)',fontsize = 30)
#    ax.set_xscale('log')
#    ax.set_yscale('log')
#    plt.plot(xf[2:NT//2], (2.0/NT * np.abs(cor_f[2:NT//2]))**2/domega, lw=4)
#    plt.legend(fontsize=20)
#    plt.suptitle('P = {}, I = {} e3 and D = {} e-5'.format(
#            round(res['x'][0],2), round(res['x'][1]*1e-3,2),
#            round(res['x'][2]*1e5,2)), fontsize = 20)
#    filename=(directory_name_opt+'/corrected_charge_noise_density.jpeg')
#    plt.savefig(filename)
#    
    return (res_power, res_Y, out)