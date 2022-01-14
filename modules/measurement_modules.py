# -*- coding: utf-8 -*-
"""
Created on March 7 2021

@author: Sisira

Import and pass all the experiment modules.
"""

import bias_position_scan as bps
import characterize_resonance as cr
import input_attenuation_scan as ias
import input_signal_measurement as ism
import output_signal_measurement as osm
import lockin_offset_measurement as lom
import error_signal_scan as ess
import error_signal_scan_ng as essng
import carrier_power_calibration as cpc
#import lockin_setup_fluctuations as lsf
import lockin_setup_fluctuations_feed as lsff
import error_signal_scan_with_feedback as essngfeed
import lockin_offset_measurement_feedback as lomfeed
import optimize_pid_simulation as ops
import extract_delay as ed
#import lockin_setup_fluctuations_pid_optimize_dc as lsfpidopt
import feedback_on_off as foo
import feedback_long as fl

def get_bias_position(directory_name, parameters, inst_call):
    bps.get_bias_position(directory_name, parameters, inst_call)
    
def characterize_resonance(directory_name, parameters, scan_no, inst_call):
    return cr.characterize_resonance(directory_name, parameters, scan_no, 
                                     inst_call)
    
def input_attenuation_scan(directory_name, parameters, inst_call):
    return ias.input_attenuation_scan(directory_name, parameters, inst_call)

def setup_resonance_scan(vna_power_bias, vna_av_bias, el_delay,
                         bias, bias_point_bg, temp, smooth_param, 
                         directory_name, flux_fn, gate_fn, gate_source,
                         gate_voltage_divider, scan_no,
                         inst_call):
    return cr.setup_resonance_scan(vna_power_bias, vna_av_bias, el_delay,
                         bias, bias_point_bg, temp, smooth_param, 
                         directory_name, flux_fn, gate_fn, gate_source,
                         gate_voltage_divider, scan_no,
                         inst_call)
    
def powerspectrum_in(directory_name, parameters, inst_call):
    return ism.input_signal_scan(directory_name, parameters, inst_call)

def calibrate_carrier_power(directory_name, parameters, inst_call):
    return cpc.scan(directory_name, parameters, inst_call)

def output_measure(directory_name, parameters, inst_call, add_str=''):
    return osm.output_signal_scan(directory_name, parameters, inst_call,
                                  add_str = add_str)

def lockin_offset(directory_name, parameters, inst_call):
    return lom.lockin_offset(directory_name, parameters, inst_call)

def lockin_offset_feed(directory_name, parameters, inst_call,
                       feed_str = ''):
    return lomfeed.lockin_offset(directory_name, parameters, inst_call,
                                 feed_str=feed_str)

def error_signal_sweep(directory_name, parameters, inst_call):
    return ess.scan(directory_name, parameters, inst_call)

def res_scan_twpa_on(directory_name, parameters, scan_no, inst_call):
    return cr.get_twpa_resonance(directory_name, parameters, scan_no,
                                 inst_call)

def error_signal_sweep_ng(directory_name, parameters, inst_call, 
                          feed_str = ''):
    return essng.scan(directory_name, parameters, inst_call, 
                      feed_str=feed_str)
    
def error_signal_sweep_ng_feed(directory_name, parameters, inst_call, 
                               feed_str = ''):
    return essngfeed.scan(directory_name, parameters, inst_call, 
                          feed_str=feed_str)

#def measure_pulse_fft(directory_name, parameters, inst_call):
#    return lsf.lockin_fft(directory_name, parameters, inst_call)
#
def measure_pulse_fft_feed(directory_name, parameters, inst_call,
                           feed_str = ''):
    return lsff.lockin_fft(directory_name, parameters, inst_call,
                           feed_str= feed_str)
    
def feedback_on_off(directory_name, parameters, inst_call,
                           feed_str = ''):
    return foo.lockin_fft(directory_name, parameters, inst_call,
                          feed_str= feed_str)
    
def feedback_long(directory_name, parameters, inst_call,
                           feed_str = ''):
    return fl.lockin_fft(directory_name, parameters, inst_call,
                          feed_str= feed_str)
    
def optimize_pid(directory_name, parameters):
    return ops.optimize_pid(directory_name, parameters)

def extract_delay(directory_name, parameters, inst_call):
    return ed.extract_delay(directory_name, parameters, inst_call)
    
#def measure_pulse_fft_pid_opt(directory_name, parameters, inst_call,
#                              feed_str = ''):
#    return lsfpidopt.lockin_fft(directory_name, parameters, inst_call,
#                                feed_str= feed_str)
#    
