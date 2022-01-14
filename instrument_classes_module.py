# -*- coding: utf-8 -*-
"""
Created on Wed Nov 14 16:36:20 2018

@author: Ben, Sisira
"""

import numpy as np
import statistics
import struct
import visa
import nidaqmx
import os, time, sys
import ctypes
import tkinter as tk
import atsapi as ats

class gpib_instrument:

    def __init__(self, addr):
        addr_str = 'GPIB0::'+str(addr)+'::INSTR'
        self.instr = visa.ResourceManager().open_resource(addr_str)
    
    def write(self, message, return_output = 0):
        # message = str
        # return_output = 1 to return the instrument's response
        if return_output:
            return self.instr.write(message)
        else:
            self.instr.write(message)
        
    def query(self, message):
        # message = str
        return self.instr.query(message)
    
    def write_raw(self, message, return_output = 0):
        # message = bytes
        if return_output:
            return self.instr.write_raw(message)
        else:
            self.instr.write_raw(message)

    def read_raw(self):
        # message = str
        return self.instr.read_raw()

    def set_timeout(self,timeout):
        # sets timeout (in milliseconds)
        # timeout = float
        self.instr.timeout = timeout

    def wai(self,):
        self.write('*wai')
        
    def opc(self):
        self.query('*opc?')
        
    def rst(self):
        self.write('*rst')
        
    def idn(self):
        return self.query('*idn?')


"""
*******************************************************************************
*******************************************************************************
"""

class keysight_n5183b(gpib_instrument):
    
    def __init__(self, addr):
        # inherit all the methods of gpib_instrument class
        super().__init__(addr)

    def set_frequency(self, freq):
        #freq = float
        message = ':freq '+str(freq)
        self.write(message)
        
    def set_power(self, power, units = 'dbm'):
        # power = float
        # units = str
        #   options: 'dbm', 'mv', 'dBuV', 'dBuVemf', 'uV', 'mVemf', 'uVemf'
        message = ':pow '+str(power)+units
        self.write(message)
        
    def set_phase(self, phase, units = 'rad'):
        # phase = float
        # units = str: 'rad' or 'deg'
        message = ':phas '+str(phase)+units
        self.write(message)
        
    def toggle_output(self, state):
        # turns RF output on or off
        # state = 0,1 (off,on)
        message = ':outp '+str(state)
        self.write(message)
        
    def toggle_modulation(self, state):
        # turns modulation on or off
        # state = 0,1 (off,on)
        message = ':outp:mod '+str(state)
        self.write(message)
        
    def toggle_pulse_mode(self, state):
        # turns pulse mode on or off
        # state = 0,1 (off,on)
        message = ':pulm:stat '+str(state)
        self.write(message)

    def toggle_alc(self, state):
        # turns on and off automatic leveling control
        #   - useful for ultra-narrow-width (UNW) pulse generation
        # state = 0,1 (off,on)
        message = ':pow:alc '+str(state)
        self.write(message)
        
    def set_pulse_source(self, source):
        # sets the source for the pulse modulation tone
        # source = str
        #   options:
        #       'ext'   -external pulse modulates output
        #       'trig'  -internal source, triggered
        #       'frun'  -internal source, free run
        #       'ado'   -internal source, see manual
        #       'doub'  -internal source, see manual
        #       'gate'  -internal source, see manual
        #       'ptr'   -internal source, see manual
        if source == 'ext':
            message = 'pulm:sour '+source
        else:
            message = 'pulm:sour:int '+source
                
        self.write(message)

    def set_pulse_delay(self, delay, units='s'):
        # sets pulse delay
        # delay = float
        # units = str: 's','us','ns'
        message = ':pulm:int:del '+str(delay)+units
        self.write(message)
        
    def set_pulse_width(self, width, units='s'):
        # sets pulse delay
        # width = float
        # units = str: 's','us','ns'
        message = ':pulm:int:pwid '+str(width)+units
        self.write(message)

"""
*******************************************************************************
*******************************************************************************
"""    
    
class agilent_e8257c(gpib_instrument):
    # SIGNAL GENERATOR
    
    def __init__(self, addr):
        # inherit all the methods of gpib_instrument class
        super().__init__(addr)

    def toggle_output(self, state):
        # turns output off or on
        # state = 0,1 (off,on)
        message = ':outp '+str(state)
        self.write(message)
        
    def toggle_modulation(self, state):
        # turns modulation off or on
        # state = 0,1 (off,on)
        message = ':outp:mod '+str(state)
        self.write(message)
    
    def set_frequency(self, freq, units = 'Hz'):
        # set the frequency
        # freq = float , scientific notation OK
        # units = str: 'Hz', 'kHz', 'MHz', 'GHz' all OK
        message = ':freq '+str(freq)+units
        self.write(message)
        
    def set_phase(self, phase, units = 'rad'):
        # phase = float
        # units = str: 'rad' or 'deg'
        message = ':phas '+str(phase)+units
        self.write(message)
        
    def toggle_alc(self, state):
        # turns automatic leveling control off,on
        # state = 0,1 (off,on)
        message = ':pow:alc '+str(state)
        self.write(message)
        
    def set_power(self, power, units = 'dbm'):
        # power = float
        # units = str
        #   options: 'dbm', 'mv', 'dBuV', 'dBuVemf', 'uV', 'mVemf', 'uVemf'
        message = ':pow '+str(power)+units
        self.write(message)
        
    def toggle_pulse_mode(self, state):
        # turns pulse mode on or off
        # state = 0,1 (off,on)
        message = ':pulm:stat '+str(state)
        self.write(message)

    def set_pulse_source(self, source):
        # sets the source for the pulse modulation tone
        # source = str
        #   options:
        #       'ext'   -external pulse modulates output
        #       'squ'   -internal source, square
        #       'trig'  -internal source, triggered
        #       'frun'  -internal source, free run
        #       'doub'  -internal source, see manual
        #       'gate'  -internal source, see manual
        if source == 'ext':
            message = 'pulm:sour '+source
        else:
            message = 'pulm:sour:int '+source
                
        self.write(message)

    def set_pulse_delay(self, delay, units='s'):
        # sets pulse delay
        # delay = float
        # units = str: 's','us','ns'
        message = ':pulm:int:del '+str(delay)+units
        self.write(message)
        
    def set_pulse_width(self, width, units='s'):
        # sets pulse delay
        # width = float
        # units = str: 's','us','ns'
        message = ':pulm:int:pwid '+str(width)+units
        self.write(message)


"""
*******************************************************************************
*******************************************************************************
"""    

class hp_83711b(gpib_instrument):
    # SIGNAL GENERATOR
    
    def __init__(self, addr):
        # inherit all the methods of gpib_instrument class
        super().__init__(addr)

    def set_frequency(self, freq, units = 'Hz'):
        # set the frequency
        # freq = float , scientific notation OK
        # units = str: 'Hz', 'kHz', 'MHz', 'GHz' all OK
        message = 'freq '+str(freq)+units
        self.write(message)

    def set_power(self, power, units = 'dBm'):
        # set the power
        # power = float
        # units = str: 'dBm', 'mV', 'V', etc.
        message = 'pow '+str(power)+units
        self.write(message)

    def toggle_output(self, state):
        # turns output off or on
        # state = 0,1 (off,on)
        message = 'outp '+str(state)
        self.write(message)


"""
*******************************************************************************
*******************************************************************************
"""    

class hp_8648c(gpib_instrument):
    # SIGNAL GENERATOR
    
    def __init__(self, addr):
        super().__init__(addr)

    def set_frequency(self, freq, units='Hz'):
        message = ':freq ' + str(freq) + units
        self.write(message)

    def toggle_output(self, state):
        message = 'outp ' + str(state)
        self.write(message)

    def set_power(self, power, units='dBm'):
        message = 'pow ' + str(power) + units
        self.write(message)

"""
*******************************************************************************
*******************************************************************************
"""    

class hp_34401a(gpib_instrument):
    # MULTIMETER
    
    def __init__(self, addr):
        super().__init__(addr)

    def get_voltage(self, current_type = 'dc'):
        """
        Measures the voltage. Multimeter chooses the best available 
        configuration.
        
        Parameters:
            current_type (str): 'dc' or 'ac'
        
        Returns:
            voltage in V
        """
        message = 'meas:volt:{}?'.format(current_type)
        return float(self.query(message))
    
    def get_current(self):
        message = 'meas:curr:dc?'
        return float(self.query(message))
    
    def set_detector_bandwidth(self,freq_min):
        """
        Sets the detector bandwidth. Useful for low or high frequencies.
        
        Parameters:
            freq_min (int): Minimum value of signal frequency (3, 20 or 200 Hz).
        """
        message = 'det:band {}'.format(freq_min)
        self.write(message)
        
    def read_voltage(self, current_type = 'ac', trig_source = 'imm',
                     sample_count = 1, ac_bw=20):
        """
        Reads the voltage according to the set configuration.
        
        Parameters:
            current_type (str): 'ac' or 'dc' measurement.
            trig_source (str): Choose the trigger from 'bus', 'imm' or 'ext'.
            sample_count (int): Counts per trigger.
        
        Returns:
            volt (tuple): Voltage values corresponding to given counts.
        """
        self.write('conf:volt:{}'.format(current_type))
        if current_type == 'ac':
            self.set_detector_bandwidth(ac_bw)
        self.write('trig:sour {}'.format(trig_source))
        self.write('samp:count {}'.format(sample_count))
        volt = self.query('read?')
        return eval(volt)
    
"""
*******************************************************************************
*******************************************************************************
"""    

class agilent_e3634a(gpib_instrument):
    # DC POWER SUPPLY
    
    def __init__(self, addr):
        super().__init__(addr)

    def toggle_output(self, state):
        message = 'outp ' + str(state)
        self.write(message)

    def apply(self, voltage, current, v_units='V', i_units='A'):
        message = 'appl ' + str(voltage) + v_units + ', ' + str(current) + i_units
        self.write(message)

    def set_voltage(self, voltage, units='V'):
        message = 'volt ' + str(voltage) + units
        self.write(message)

    def set_current(self, current, units='A'):
        message = 'curr ' + str(current) + units
        self.write(message)

    def measure_voltage(self):
        message = 'meas:volt?'
        return self.query(message)

    def measure_current(self):
        message = 'meas:curr?'
        return self.query(message)

"""
*******************************************************************************
*******************************************************************************
"""    

class hp_6612c(gpib_instrument):
    # DC POWER SUPPLY
    
    def __init__(self, addr):
        super().__init__(addr)

    def toggle_output(self, state):
        message = 'outp ' + str(state)
        self.write(message)
    
    def set_voltage(self, volt):
        message = 'volt ' +str(volt)
        self.write(message)
    
    def get_voltage(self):
        message = 'meas:volt?'
        return float(self.query(message))


"""
*******************************************************************************
*******************************************************************************
"""    

class keysight_n9000b(gpib_instrument):
    # SPECTRUM ANALYZER
    
    def __init__(self, addr):
        # inherit all the methods of gpib_instrument class
        super().__init__(addr)
        
    def set_mode(self, mode_str):
        # sets analysis mode
        # mode_str = str
        #   options:
        #       'sa'    - spectrum analysis
        #       'basic' - basic I/Q analysis
        message = 'inst:sel '+mode_str
        self.write(message)
        
    def select_screen(self, screen_str):
        # screen_str = name of screen surrounded by double quotes
        message = ':inst:scr:sel '+screen_str
        self.write(message)
        
    def read_trace(self,trace=1):
        """
        Performs a sweep and returns data
        # format: freq_str, pow_str, freq_str, pow_str, ...
        # I omit the new line character at the end of the trace
        
        Parameters:
            trace (int): Trace number. Default to 1.
        """
        message = 'read:san{}?'.format(trace)
        return self.query(message)[:-1]
    
    def fetch_trace(self):
        # retrieves trace already taken
        # format: freq_str, pow_str, freq_str, pow_str, ...
        # I omit the new line character at the end of the trace
        message = 'fetch:san?'
        return self.query(message)[:-1]
    
    def get_all_peaks(self):
        # get the peaks of a trace in x,y pairs
        # I exclude the final new line character \n
        message = 'trac:math:peak?'
        return self.query(message)[:-1]
    
    def get_peaks(self, threshold,convert=True):
        """
        Get peak from the active trace with the set threshold and zero
        excursion.
        
        Parameters:
            threshold (float): Threshold in dBm for identifying the peak.
            
        Returns:
            freq (list): Frequency of peaks
            power (list): Corresponding power of the frequencies
        """
        self.write('CALC:MARK:PEAK:THR {}'.format(threshold))
        message = 'TRAC:MATH:PEAK?'
        peaks = self.query(message)[:-1]
        if convert:
            peaks = peaks.split(',')
            peaks = [float(ii) for ii in peaks]
        freq = peaks[::2]
        power = peaks[1::2]
        return (freq, power)
    
    def toggle_automatic_resolution_bandwidth(self, state):
        # turn automatic coupling of resolution bw off or on
        # state = 0,1 (off,on)
        message = ':band:auto ' + str(state)
        self.write(message)
        
    def set_resolution_bandwidth(self, bandwidth, units='Hz'):
        # set the resolution bandwidth
        # bandwidth = float: scientific notation OK
        # units = str
        message = ':band ' + str(bandwidth) + units
        self.write(message)
        
    def get_resolution_bandwidth(self):
        # query the resolution bandwidth (in Hz)
        message = ':band?'
        return float(self.query(message))
    
    def set_video_bandwidth(self, bandwidth, units='Hz'):
        # set the video bandwidth
        # bandwidth = float: scientific notation OK
        # units = str
        message = ':band:vid ' + str(bandwidth) + units
        self.write(message)

    def get_video_bandwidth(self):
        # query the video bandwidth (in Hz)
        message = ':band:vid?'
        return float(self.query(message))

    def toggle_automatic_video_bandwidth(self, state):
        # turn automatic coupling of video bw off or on
        # state = 0,1 (off,on)
        message = ':band:vid:auto ' + str(state)
        self.write(message)

    def set_frequency_center(self, frequency, units='Hz'):
        # frequency = float: scientific notation OK
        message = ':freq:cent ' + str(frequency) + units
        self.write(message)

    def set_frequency_span(self, frequency, units='Hz'):
        # frequency = float: scientific notation OK
        message = ':freq:span ' + str(frequency) + units
        self.write(message)

    def set_frequency_start(self, frequency, units='Hz'):
        # frequency = float: scientific notation OK
        message = ':freq:star ' + str(frequency) + units
        self.write(message)

    def set_frequency_stop(self, frequency, units='Hz'):
        # frequency = float: scientific notation OK
        message = ':freq:stop ' + str(frequency) + units
        self.write(message)

    def set_freqs(self, freq1, freq2, interval='span'):
        if interval == 'range':
            self.set_frequency_start(freq1)
            self.set_frequency_stop(freq2)
        else:
            if interval == 'span':
                self.set_frequency_center(freq1)
                self.set_frequency_span(freq2)
            else:
                print('ERROR: Invalid Interval Type!')
                
    def set_averaging(self, counts):
        # counts = int
        # note: may initiate or augment a sweep
        message = ':aver:count ' + str(counts)
        self.write(message)
        
    def set_averaging_type(self, type_str):
        # type_str = str
        #   options: 'log', 'rms', 'scal'
        message = ':aver:type ' + type_str
        self.write(message)
        
    def toggle_continuous_sweep(self, state):
        # on: continuous sweep
        # off: single sweep
        message = ':init:cont ' + str(state)
        self.write(message)

    def restart(self):
        # restart current sweep
        message = ':init:rest'
        self.write(message)
        
    def set_sweep_points(self, num_points):
        # sets the number of frequency points in a sweep
        # max = 20,001
        message = ':swe:poin ' + str(num_points)
        self.write(message)
        
    def abort(self):
        # abort current sweep
        message = ':abor'
        self.write(message)
        
    def set_trace_type(self, type_str):
        # type_str = str
        #   options:
        #       'WRITe'     - clear/write
        #       'AVERage'   - trace average
        #       'MAXHold'
        #       'MINHold'
        message = 'trac:type '+type_str
        self.write(message)
        
    def set_detector_type(self, type_str, trace=1):
        """Sets the detector type for selected trace.
        
        # type_str = str
        #   options:
        #       'AVERage'   - good for noise
        #       'NEGative'
        #       'NORMal'    - works simultaneously for pure tones and noise
        #       'POSitive'  - positive peak, good for measuring pure tones
        #       'SAMPle'    - good for noise
        
        Parameters:
            type_str (str): Detector type as listed above.
            trace (int): Trace number. Selects trace 1 if none given.
        """
        message = ':det:trac{} {}'.format(trace,type_str)
        self.write(message)
        
    def trigger(self, wait=True):
        # force trigger
        message = '*trg'
        self.write(message)
        if wait:
            return self.query('*opc?')

    def set_iq_resolution_bandwidth(self, freq, units='Hz'):
        # set the resolution bandwidth for the IQ spectrum
        message = ':spec:band '+str(freq)+units
        self.write(message)
        
    def set_iq_averaging(self, count):
        # set averaging count for IQ spectrum
        message = ':spec:aver:coun '+str(count)
        self.write(message)
        
    def read_iq_waveform(self, return_time_config = True):
        # return iq waveform and time configuration data
        # the waveform and config will be split by a new line '\n'
        message = 'read:wav0?'
        waveform = self.query(message)
        if return_time_config:
            message2 = 'fetch:wav1?'
            time_config = self.query(message2)[:-1]
            waveform += time_config
        return waveform
    
    def set_ref_level(self, value_dbm):
        """
        Sets the reference value.
        
        Parameters:
            value_dbm (float) : Value of reference in dBm.
        """
        message = 'DISP:WIND:TRAC:Y:RLEV {} dBm'.format(value_dbm)
        self.write(message)
        
    def set_marker_trace(self, marker, trace):
        """
        Places the marker on the trace specified.
        
        Parameters:
            marker (int) : Marker number.
            trace (int) : Trace number.
        """
        message = 'calc:mark{}:trac {}'.format(marker,trace)
        self.write(message)
        
    def set_marker_band_function(self, marker, function):
        """
        Sets the band function for the selected marker.
        
        Parameters:
            marker (int): Marker number.
            function (str): Band function strings 
            Options are 'NOISe','BPOWer','BDENsity','OFF'
        """
        message = 'calc:mark{}:func {}'.format(marker, function)
        self.write(message)
        
    def get_marker_amplitude(self, marker):
        """
        Gets the amplitude of the selected marker.
        
        Parameters:
            marker (int): Marker number.
        Returns:
            amp (float): Amplitude in dBm if band function is off. If the band
            function is set to noise, the amplitude is given in dBm/Hz.
            unit(str): Gives the unit of amplitude - dBm or dBm/Hz.
        """
        message = ':calc:mark{}:y?'.format(marker)
        amp = eval(self.query(message))
        message = 'calc:mark{}:func?'.format(marker)
        function = self.query(message)[:-1]
        if function == 'NOIS' or function == 'BDEN':
            unit = 'dBm/Hz'
        else:
            unit = 'dBm'
        
        return amp, unit
    
    def clear_all_traces(self):
        """
        Clears all the traces.
        """
        message = ':trac:cle:all'
        self.write(message)
        
    def set_marker_freq(self, marker, freq):
        """
        Places the marker at the selected frequency.
        
        Parameters:
            marker(int): Marker number.
            freq (float): Frequency at which the marker is to be placed (in Hz).
        """
        message = ':calc:mark{}:x {}'.format(marker,freq)
        self.write(message)
        
    def toggle_automatic_sweeptime(self,state):
        """
        Toggles the sweep time to automatic or manual.
        
        Parameters:
            state (binary): 0 or 1.
        """
        message = ':SWE:TIME:AUTO {}'.format(state)
        self.write(message)
        
    def set_y_scaleperdivision(self,scale,unit='dB'):
        """
        Sets the scale per division of the amplitude.
        
        Parameters:
            scale (float): Scale per division value.
            unit (str): Unit default set to dB.
        """
        message = 'DISP:WIND:TRAC:Y:PDIV {} {}'.format(scale, unit)
        self.write(message)
        
    def set_ref_level_offset(self,offset):
        """
        Sets the reference level in dB.
        
        Parameters:
            offset (float): Value of offset in dB.
        """
        message = 'DISP:WIND:TRAC:Y:RLEV:OFFS {}'.format(offset)
        self.write(message)
        
    def set_trace_average_hold_number(self,count):
        """
        Sets the number of measurement averages.
        
        Parameters:
            count (int): Hold number.
        """
        message = 'CHP:AVER:COUN {}'.format(count)
        self.write(message)
        
    def set_input_coupling(self, coupling_str):
        message = 'inp:coup ' + coupling_str
        self.write(message)


"""
*******************************************************************************
*******************************************************************************
"""    

class agilent_e4404b(gpib_instrument):
    # SPECTRUM ANALYZER
    # ONLY TESTED FOR 4408B, SHOULD BE THE SAME
    
    def __init__(self, addr):
        # inherit all the methods of gpib_instrument class
        super().__init__(addr)
        
    def abort(self):
        message = ':abor'
        self.write(message)

    def force_trigger(self):
        message = '*trg'
        self.write(message)

    def toggle_coupling(self, state):
        if state:
            cpl_str = 'all'
        else:
            cpl_str = 'none'
        message = ':coup ' + cpl_str
        self.write(message)

    def toggle_continuous_sweep(self, state):
        message = ':init:cont ' + str(state)
        self.write(message)

    def restart(self):
        message = ':init:rest'
        self.write(message)

    def set_input_coupling(self, coupling_str):
        message = ':inp:coup ' + coupling_str
        self.write(message)

    def toggle_averaging(self, state):
        message = ':aver ' + str(state)
        self.write(message)

    def set_averaging(self, counts):
        message = ':aver:count ' + str(counts)
        self.write(message)

    def set_averaging_type(self, type_str):
        message = ':aver:type ' + type_str
        self.write(message)

    def set_resolution_bandwidth(self, bandwidth, units='Hz'):
        message = ':band ' + str(bandwidth) + units
        self.write(message)

    def get_resolution_bandwidth(self):
        message = ':band?'
        return float(self.query(message))

    def toggle_automatic_resolution_bandwidth(self, state):
        message = ':band:auto ' + str(state)
        self.write(message)

    def set_video_bandwidth(self, bandwidth, units='Hz'):
        message = ':band:vid ' + str(bandwidth) + units
        self.write(message)

    def get_video_bandwidth(self):
        message = ':band:vid?'
        return float(self.query(message))

    def toggle_automatic_video_bandwidth(self, state):
        message = ':band:vid:auto ' + str(state)
        self.write(message)

    def set_detection_type(self, type_str):
        message = ':det ' + type_str
        self.write(message)

    def set_frequency_center(self, frequency, units='Hz'):
        message = ':freq:cent ' + str(frequency) + units
        self.write(message)

    def set_frequency_span(self, frequency, units='Hz'):
        message = ':freq:span ' + str(frequency) + units
        self.write(message)

    def set_frequency_start(self, frequency, units='Hz'):
        message = ':freq:star ' + str(frequency) + units
        self.write(message)

    def set_frequency_stop(self, frequency, units='Hz'):
        message = ':freq:stop ' + str(frequency) + units
        self.write(message)

    def set_freqs(self, freq1, freq2, interval='range', channel=None):
        if interval == 'range':
            self.set_frequency_start(freq1)
            self.set_frequency_stop(freq2)
        else:
            if interval == 'span':
                self.set_frequency_center(freq1)
                self.set_frequency_span(freq2)
            else:
                print('ERROR: Invalid Interval Type!')

    def set_sweep_points(self, num_points):
        message = ':swe:poin ' + str(num_points)
        self.write(message)

    def get_trace_data(self, convert=True):
        message = ':trac? trace1'
        trace = self.query(message)
        if convert:
            trace = trace.split(',')
            trace = [float(ii) for ii in trace]
        return trace

"""
*******************************************************************************
*******************************************************************************
"""    

class agilent_e4408b(gpib_instrument):
    # SPECTRUM ANALYZER
    
    def __init__(self, addr):
        # inherit all the methods of gpib_instrument class
        super().__init__(addr)

    def abort(self):
        message = ':abor'
        self.write(message)

    def force_trigger(self):
        message = '*trg'
        self.write(message)

    def toggle_coupling(self, state):
        if state:
            cpl_str = 'all'
        else:
            cpl_str = 'none'
        message = ':coup ' + cpl_str
        self.write(message)

    def toggle_continuous_sweep(self, state):
        message = ':init:cont ' + str(state)
        self.write(message)

    def restart(self):
        message = ':init:rest'
        self.write(message)

    def set_input_coupling(self, coupling_str):
        message = ':inp:coup ' + coupling_str
        self.write(message)

    def toggle_averaging(self, state):
        message = ':aver ' + str(state)
        self.write(message)

    def set_averaging(self, counts):
        message = ':aver:count ' + str(counts)
        self.write(message)

    def set_averaging_type(self, type_str):
        message = ':aver:type ' + type_str
        self.write(message)

    def set_resolution_bandwidth(self, bandwidth, units='Hz'):
        message = ':band ' + str(bandwidth) + units
        self.write(message)

    def get_resolution_bandwidth(self):
        message = ':band?'
        return float(self.query(message))

    def toggle_automatic_resolution_bandwidth(self, state):
        message = ':band:auto ' + str(state)
        self.write(message)

    def set_video_bandwidth(self, bandwidth, units='Hz'):
        message = ':band:vid ' + str(bandwidth) + units
        self.write(message)

    def get_video_bandwidth(self):
        message = ':band:vid?'
        return float(self.query(message))

    def toggle_automatic_video_bandwidth(self, state):
        message = ':band:vid:auto ' + str(state)
        self.write(message)

    def set_detection_type(self, type_str):
        message = ':det ' + type_str
        self.write(message)

    def set_frequency_center(self, frequency, units='Hz'):
        message = ':freq:cent ' + str(frequency) + units
        self.write(message)

    def set_frequency_span(self, frequency, units='Hz'):
        message = ':freq:span ' + str(frequency) + units
        self.write(message)

    def set_frequency_start(self, frequency, units='Hz'):
        message = ':freq:star ' + str(frequency) + units
        self.write(message)

    def set_frequency_stop(self, frequency, units='Hz'):
        message = ':freq:stop ' + str(frequency) + units
        self.write(message)

    def set_freqs(self, freq1, freq2, interval='range', channel=None):
        if interval == 'range':
            self.set_frequency_start(freq1)
            self.set_frequency_stop(freq2)
        else:
            if interval == 'span':
                self.set_frequency_center(freq1)
                self.set_frequency_span(freq2)
            else:
                print('ERROR: Invalid Interval Type!')

    def set_sweep_points(self, num_points):
        message = ':swe:poin ' + str(num_points)
        self.write(message)

    def get_trace_data(self, convert=True):
        message = ':trac? trace1'
        trace = self.query(message)
        if convert:
            trace = trace.split(',')
            trace = [float(ii) for ii in trace]
        return trace

"""
*******************************************************************************
*******************************************************************************
"""    

class agilent_e5071c(gpib_instrument):
    # NETWORK ANALYZER
    
    def __init__(self, addr):
        # inherit all the methods of gpib_instrument class
        super().__init__(addr)
    
    ####################
    # DISPLAY SETTINGS #
    ####################
    
    def allocate_channels(self, alloc_str):
        # splits the display into separate windows
        # channel_str = str
        #   options:
        #       '1'    - one graph in full window
        #       '12'   - two graphs split left/right
        #       '1_2'  - two graphs split top/bottom
        #   see manual for more complicated allocations
        # NOTE: the argument given in the manual has a 'D' in front - this is
        #   handled internally by this method
        message = ':disp:spl D' + alloc_str
        self.write(message)

    def set_channel(self,channel):
        # sets active channel
        # channel = int
        message = ':disp:wind'+str(channel)+':act'
        self.write(message)

    def query_channel(self):
        # queries the currently active channel
        # returns the channel as an int
        message = ':serv:chan:act?'
        channel = self.query(message)
        return int(channel)
        
    def set_num_traces(self, num_traces, channel = None):
        # sets the total number of traces present on channel
        # num_traces = int
        # channel = int: if unspecified, defaults to current channel
        if not channel:
            channel = self.query_channel()
        message = ':calc'+str(channel)+':par:coun '+str(num_traces)
        self.write(message)
        
    def query_num_traces(self, channel = None):
        # queries the total number of traces present on channel
        # if channel is unspecified, defaults to the current channel
        if not channel:
            channel = self.query_channel()
        message = ':calc'+str(channel)+':par:coun?'
        num_traces = self.query(message)
        return int(num_traces)   
    
    def allocate_traces(self, alloc_str, channel = None):
        # splits the layout of traces on a given channel
        # if channel is unspecified, defaults to the current channel
        # alloc_str = str
        #   options:
        #       '1'    - one graph in full window
        #       '12'   - two graphs split left/right
        #       '1_2'  - two graphs split top/bottom
        #   see manual for more complicated allocations
        # NOTE: the argument given in the manual has a 'D' in front - this is
        # handled internally by this method
        if not channel:
            channel = self.query_channel()
        message = ':disp:wind'+str(channel)+':spl D'+alloc_str
        self.write(message)
        
    def set_trace(self, trace, channel = None):
        # sets active trace of a given channel
        # if channel is not given, then the current channel is assumed
        # trace = int
        # channel = int
        if not channel:
            channel = self.query_channel()
        message = ':calc'+str(channel)+':par'+str(trace)+':sel'
        self.write(message)
    
    def query_trace(self, channel = None):
        # queries the active trace on channel
        # channel = int
        # if channel unspecified, defaults to current channel
        # returns trace as an int
        if not channel:
            channel = self.query_channel()
        message = ':serv:chan'+str(channel)+':trac:act?'
        trace = self.query(message)
        return int(trace)
    
    def autoscale(self, channel = None, trace = None):
        # autoscales a given trace on a given channel
        # if channel/trace unspecified, defaults to current
        if not channel:
            channel = self.query_channel()
        if not trace:
            trace = self.query_trace(channel)
        message = ':disp:wind'+str(channel)+':trac'+str(trace)+':y:auto'
        self.write(message)
    
    ###########
    # MARKERS #
    ###########
    # NOTE: MARKER 10 IS THE REFERENCE MARKER
    
    def toggle_marker(self, marker, state, channel = None):
        # turns marker for a given channel on or off
        # if channel unspecified, defaults to current
        # marker = int
        # state = 0,1 (off,on)
        if not channel:
            channel = self.query_channel()
        message = ':calc'+str(channel)+':mark'+str(marker)+' '+str(state)
        self.write(message)
        
    def move_marker(self, marker, freq, channel = None):
        # move a marker for a given channel to a given frequency
        # if channel unspecified, defaults to current
        # marker = int
        # freq = float (scientific notation OK)
        if not channel:
            channel = self.query_channel()
        message = ':calc'+str(channel)+':mark'+str(marker)+':x '+str(freq)
        self.write(message)
        
    def activate_marker(self, marker, channel = None):
        # activate a marker for a given channel
        # if channel unspecified, defaults to current
        if not channel:
            channel = self.query_channel()
        message = ':calc'+str(channel)+':mark'+str(marker)+':act'
        self.write(message)
        
    def marker_search(self, marker, type_str, channel = None, trace = None):
        # executes marker search type for marker on trace on channel
        # if channel/trace unspecified, defaults to current
        # marker = int
        # type_str = str
        #   options:
        #       'min'   - minimum
        #       'max'   - maximum
        #       'peak'  - peak
        #       'lpe'   - peak to the left
        #       'rpe'   - peak to the right
        #       'targ'  - target
        #       'ltar'  - target to the left
        #       'rtar'  - target to the right
        if not channel:
            channel = self.query_channel()
        if not trace:
            trace = self.query_trace()
        message = ':calc'+str(channel)+':trac'+str(trace)+':mark'+str(marker) \
                        +':func' 
        self.write(message+':type '+type_str)
        self.write(message+':exec')
        
    def marker_track(self, marker, type_str, channel = None, trace = None):
        # performs a marker search every sweep according to the type specified
        # if channel/trace unspecified, defaults to current
        # marker = int
        # type_str = str
        #   options:
        #       'min'   - minimum
        #       'max'   - maximum
        #       'peak'  - peak
        #       'lpe'   - peak to the left
        #       'rpe'   - peak to the right
        #       'targ'  - target
        #       'ltar'  - target to the left
        #       'rtar'  - target to the right
        if not channel:
            channel = self.query_channel()
        if not trace:
            trace = self.query_trace()
        message = ':calc'+str(channel)+':trac'+str(trace)+':mark'+str(marker) \
                        +':func' 
        self.write(message+':type '+type_str)
        self.toggle_marker_tracking(marker, 1, channel, trace)
        
    def toggle_marker_tracking(self, marker, state, channel = None, \
                               trace = None):
        # toggles tracking for a given marker
        # marker = int
        # state = 0,1 (off,on)
        if not channel:
            channel = self.query_channel()
        if not trace:
            trace = self.query_trace()
        message = ':calc'+str(channel)+':trac'+str(trace)+':mark'+str(marker) \
                            +':func:trac '+str(state)
        self.write(message)
        
    def toggle_marker_search_range(self, state, channel = None, trace = None):
        # toggles whether to use a manual search range
        # state = 0,1 (off,on)
        if not channel:
            channel = self.query_channel()
        if not trace:
            trace = self.query_trace()
        message = ':calc'+str(channel)+':trac'+str(trace) \
                        +':mark:func:dom '+str(state)
        self.write(message)
    
    def set_marker_search_start(self, freq, channel = None, trace = None):
        # sets the starting frequency of the marker search domain
        if not channel:
            channel = self.query_channel()
        if not trace:
            trace = self.query_trace()
        message = ':calc'+str(channel)+':trac'+str(trace) \
                        +':mark:func:dom:star '+str(freq)
        self.write(message)
            
    def set_marker_search_stop(self, freq, channel = None, trace = None):
        # sets the starting frequency of the marker search domain
        if not channel:
            channel = self.query_channel()
        if not trace:
            trace = self.query_trace()
        message = ':calc'+str(channel)+':trac'+str(trace) \
                        +':mark:func:dom:stop '+str(freq)
        self.write(message)
        
    def toggle_bandwidth_search(self, state, channel = None, trace = None):
        # toggles bandwidth search mode
        # state = 0,1 (off,on)
        if not channel:
            channel = self.query_channel()
        if not trace:
            trace = self.query_trace()
        message = ':calc'+str(channel)+':trac'+str(trace)+':mark:bwid ' \
                    +str(state)
        self.write(message)
        
    def set_bandwidth_threshold(self, marker, threshold, \
                                    channel = None, trace = None):
        # sets the threshold for determining bandwidth
        # threshold = float : +/- 3dB usually
        # marker = int
        if not channel:
            channel = self.query_channel()
        if not trace:
            trace = self.query_trace()
        message = ':calc'+str(channel)+':trac'+str(trace)+':mark'+str(marker) \
                    +':bwid:thr '+str(threshold)
        self.write(message)
    
    def track_resonance(self, marker = 1, channel = None, trace = None):
        # track resonance in the current window, find 3dB points and Q
        if not channel:
            channel = self.query_channel()
        if not trace:
            trace = self.query_trace()
        self.toggle_bandwidth_search(1, channel, trace)
        self.set_bandwidth_threshold(marker, +3, channel, trace)
        self.toggle_marker(marker, 1, channel)
        self.activate_marker(marker, channel)
        self.marker_track(marker, 'min', channel, trace)

    ########################
    # MEASUREMENT SETTINGS #
    ########################
    
    def toggle_output(self, state):
        # turns RF Out off or on
        # state = 0,1 (off,on)
        message = ':outp '+str(state)
        self.write(message)
    
    def set_measurement(self, meas_str, channel = None, trace = None):
        # sets the measurement carried out by a trace on a channel
        # meas_str = str
        #   options:
        #       'S11' - reflection on port 1
        #       'S21' - transmission to port 2 from port 1
        #       'S12' - transmission to port 1 from port 2
        #       'S22' - reflection on port 2
        #   see manual for others
        # channel = int
        # trace = int
        if not channel:
            channel = self.query_channel()
        if not trace:
            trace = self.query_trace()
        message = ':calc'+str(channel)+':par'+str(trace)+':def '+meas_str
        self.write(message)

    def query_measurement(self, channel = None, trace = None):
        # queries the current measurement carried out by a trace on a channel
        if not channel:
            channel = self.query_channel()
        if not trace:
            trace = self.query_trace()
        message = ':calc'+str(channel)+':par'+str(trace)+':def?'
        meas_str = self.query(message)
        return meas_str[:-1]
        
    def set_sweep_type(self, type_str, channel = None):
        # sets the type of sweep on a given channel
        # if channel not specified, defaults to current channel
        # type_str = str
        #   options:
        #       'lin' - linear frequency sweep
        #       'log' - logarithmic frequency sweep
        #       'segm' - segment sweep
        #       'pow' - power sweep
        if not channel:
            channel = self.query_channel()
        message = ':sens'+str(channel)+':swe:type '+type_str
        self.write(message)
        
    def query_sweep_type(self, channel = None):
        # queries type of sweep on channel
        # if channel not specified, defaults to current
        # channel = int
        if not channel:
            channel = self.query_channel()
        message = ':sens'+str(channel)+':swe:type?'
        type_str = self.query(message)
        return type_str[:-1]
    
    def set_sweep_mode(self, mode_str, channel = None):
        # sets sweep mode of a channel
        # if channel unspecified, defaults to current
        # channel = int
        # mode_str = str
        #   options (only caps part necessary):
        #       'STEPped'   - stepped mode
        #       'ANALog'    - swept mode
        #   manual specified a couple more options that don't seem to do
        #   anything for this model network analyzer
        if not channel:
            channel = self.query_channel()
        message = ':sens'+str(channel)+':swe:gen '+mode_str
        self.write(message)

    def query_sweep_mode(self, channel = None):
        # queries sweep mode of a channel
        # if channel unspecified, defaults to current
        # channel = int
        if not channel:
            channel = self.query_channel()
        message = ':sens'+str(channel)+':swe:gen?'
        mode_str = self.query(message)
        return mode_str[:-1]

    def set_sweep_points(self, points, channel = None):
        # sets the number of points in the sweep on the specified channel
        # if channel unspecified, defaults to current
        # points = int
        if not channel:
            channel = self.query_channel()
        message = ':sens'+str(channel)+':swe:poin '+str(points)
        self.write(message)

    def set_frequency_start(self, freq, channel = None):
        # sets the start frequency of the sweep on a channel
        # if channel unspecified, defaults to current
        # freq = float: frequency in Hz
        #               scientific notation also OK
        # channel = int
        if not channel:
            channel = self.query_channel()
        message = ':sens'+str(channel)+':freq:star '+str(freq)
        self.write(message)
        
    def set_frequency_stop(self, freq, channel = None):
        # sets the stop frequency of the sweep on a channel
        # if channel unspecified, defaults to current
        # freq = float: frequency in Hz
        #               scientific notation also OK
        # channel = int
        if not channel:
            channel = self.query_channel()
        message = ':sens'+str(channel)+':freq:stop '+str(freq)
        self.write(message)
        
    def set_frequency_center(self, freq, channel = None):
        # sets the center frequency of the sweep on a channel
        # if channel unspecified, defaults to current
        # freq = float: frequency in Hz
        #               scientific notation also OK
        # channel = int
        if not channel:
            channel = self.query_channel()
        message = ':sens'+str(channel)+':freq:cent '+str(freq)
        self.write(message)
        
    def set_frequency_span(self, freq, channel = None):
        # sets the center frequency of the sweep on a channel
        # if channel unspecified, defaults to current
        # freq = float: frequency in Hz
        #               scientific notation also OK
        # channel = int
        if not channel:
            channel = self.query_channel()
        message = ':sens'+str(channel)+':freq:span '+str(freq)
        self.write(message)
        
    def set_freqs(self, freq1, freq2, interval_type = 'span', channel = None):
        # sets a frequency interval according to the type
        # interval_type = str: 'range' or 'span'
        #   if range:
        #       freq1 = start
        #       freq2 = stop
        #   if span:
        #       freq1 = center
        #       freq2 = span
        # freq1,2 = float
        if not channel:
            channel = self.query_channel()
        if interval_type == 'range':
            self.set_frequency_start(freq1,channel)
            self.set_frequency_stop(freq2,channel)
        elif interval_type == 'span':
            self.set_frequency_center(freq1,channel)
            self.set_frequency_span(freq2,channel)
        else:
            print('ERROR: Invalid Interval Type!')
                
    def set_power(self, power, channel = None):
        # sets the power output of a given channel
        # if channel unspecified, defaults to current
        # power = float: power in dBm
        if not channel:
            channel = self.query_channel()
        message = ':sour'+str(channel)+':pow '+str(power)
        self.write(message)
                
    def set_format(self, format_str, channel = None, trace = None):
        # sets the measurement format of a given trace on a given channel
        # if channel/trace not specified, defaults to current
        # format_str = str
        #   options:
        #       'MLOGarithmic'  - logarithmic magnitude
        #       'PHASe'         - Phase
        #       'UPHase'        - Expanded phase
        #       'SMITh'         - Smith chart: R+jX (resistance/reactance)
        if not channel:
            channel = self.query_channel()
        if not trace:
            trace = self.query_trace(channel)
        message = ':calc'+str(channel)+':trac'+str(trace)+':form '+format_str
        self.write(message)
    
    def query_format(self, channel = None, trace = None):
        # queries the format of a given troce on a given channel
        # if channel/trace not specified, defaults to current
        if not channel:
            channel = self.query_channel()
        if not trace:
            trace = self.query_trace(channel)
        message = ':calc'+str(channel)+':trac'+str(trace)+':form?'
        format_str = self.query(message)
        return format_str[:-1]
    
    def set_electrical_delay(self, delay, channel = None, trace = None):
        # sets the electrical delay for phase measurements
        # if channel/trace unspecified, defaults to current
        # delay = float: time in seconds (scientific notation OK)
        if not channel:
            channel = self.query_channel()
        if not trace:
            trace = self.query_trace(channel)
        message = ':calc'+str(channel)+':trac'+str(trace)+':corr:edel:time ' \
                        + str(delay)
        self.write(message)

    #####################
    # AVERAGING SETINGS #
    #####################
    
    def toggle_averaging(self, state, channel = None):
        # toggle averaging on or off for a given channel
        # if channel unspecified, defaults to current
        # state = int: 0 or 1 (off or on)
        # channel = int
        if not channel:
            channel = self.query_channel()
        message = ':sens'+str(channel)+':aver '+str(state)
        self.write(message)
        
    def set_averaging(self, factor, channel = None):
        # set averaging factor for a given channel
        # if channel unspecified, defaults to current
        # factor = int
        # channel = int
        if not channel:
            channel = self.query_channel()
        message = ':sens'+str(channel)+':aver:coun '+str(factor)
        self.write(message)
        
    def restart_averaging(self, channel = None):
        # restarts averaging on a given channel
        # if channel not specified, defaults to current
        if not channel:
            channel = self.query_channel()
        message = ':sens'+str(channel)+':aver:cle'
        self.write(message)
        
    def toggle_averaging_trigger(self, state):
        # toggles averaging on trigger initiation
        # usage: collect the avg of many traces after a single trigger
        # state = 0,1 (off,on)
        message = ':trig:aver '+str(state)
        self.write(message)
    
    def set_if_bandwidth(self, bandwidth, channel = None):
        # sets the IF bandwidth of a given channel
        #   -can affect noisiness of measurements - see manual
        # if channel unspecified, defaults to current
        # bandwidth = float (frequency in Hz - scientific notation OK)
        # channel = int
        if not channel:
            channel = self.query_channel()
        message = ':sens'+str(channel)+':band '+str(bandwidth)
        self.write(message)
        
    def toggle_smoothing(self, state, trace = None, channel = None):
        # state = 0,1 (off,on)
        if not channel:
            channel = self.query_channel()
        if not trace:
            trace = self.query_trace(channel)
        message = ':calc'+str(channel)+':trac'+str(trace)+':smo '+str(state)
        self.write(message)
        
    def set_smoothing(self, percent, trace = None, channel = None):
        # percent = float, 0.05->25
        if not channel:
            channel = self.query_channel()
        if not trace:
            trace = self.query_trace(channel)
        message = ':calc'+str(channel)+':trac'+str(trace)+':smo:aper ' \
                    +str(percent)
        self.write(message)
    
    ####################
    # TRIGGER SETTINGS #
    ####################

    def toggle_continuous_triggering(self, state, channel = None):
        # turns continuous triggering mode on or off for a given channel
        # if channel unspecified, defaults to current
        # state = 0,1 (off,on)
        if not channel:
            channel = self.query_channel()
        message = ':init'+str(channel)+':cont '+str(state)
        self.write(message)
        
    def set_trigger_source(self, source_str):
        # sets the source of triggers
        # source_str = str
        #   options:
        #       'int' - internal (always on, I think, so no real control)
        #       'ext' - external (back panel)
        #       'man' - manual (use this one for taking data)
        #       'bus' - bus (just responds to '*TRG' SCPI command - useless?)
        message = ':trig:sour '+source_str
        self.write(message)

    def set_trigger_scope(self, scope_str):
        # sets whether all channels or just the active channel are triggered
        # scope_str = str: 'all' or 'act'
        message = ':trig:scop '+scope_str
        self.write(message)

    def trigger(self, wait = True):
        # triggers a single measurement
        # can be used in conjunction with a query of '*OPC?' to see when
        #   the measurement has completed
        # wait = 0,1 or False,True
        #   if 1: waits for operation to complete before returning
        message = ':trig:sing'
        self.write(message)
        if wait:
            return self.query('*OPC?')

    ################
    # INPUT/OUTPUT #
    ################

    def transfer_data_to_memory(self, channel = None, trace = None):
        # transfers data trace to memory trace (so it is not overwritten by
        #   subsequent datasets until calling this function again)
        # if channel/trace unspecified, default to current
        if not channel:
            channel = self.query_channel()
        if not trace:
            trace = self.query_trace(channel)
        message = ':calc'+str(channel)+':trac'+str(trace)+':math:mem'
        self.write(message)
        
    def get_trace_data(self, channel = None, trace = None, \
                               convert = False):
        # gets data from the memory trace on the VNA and returns a local copy
        # if channel/trace unspecified, defaults to current
        # if convert is true, will return the trace data in an array of
        #   floats (or complex numbers in the case of a smith chart)
        # NOTE: DOESN'T NECESSARILY SUPPORT ALL TRACE FORMATS CURRENTLY
        if not channel:
            channel = self.query_channel()
        if not trace:
            trace = self.query_trace(channel)
        message = ':calc'+str(channel)+':trac'+str(trace)+':data:fmem?'
        data_str = self.query(message)

        format_str = self.query_format(channel,trace)
        data_list = data_str.split(',')
        return_list = []
        for ii in range(len(data_list)//2):
            if format_str == 'SMIT':
                resistance = float(data_list[2*ii])
                reactance = float(data_list[((2*ii)+1)])
                if convert:
                    return_list.append(complex(resistance,reactance))
                if not convert:
                    return_list.append(str(complex(resistance,reactance)))
                    return_str = ','.join(return_list)
            else:
                if convert:
                    return_list.append(float(data_list[2*ii]))
                else:
                    return_list.append(data_list[2*ii])
                    return_str = ','.join(return_list)
        if convert:
            return return_list
        else:
            return return_str

            
    def get_frequency_data(self, channel = None, convert = False):
        # gets all the frequency points for a given channel
        # if convert to floats, will return the trace data
        if not channel:
            channel = self.query_channel()
        message = ':sens'+str(channel)+':freq:data?'
        freq_str = self.query(message)
        if convert:
            freq_list = freq_str.split(',')
            freq_list = [float(ii) for ii in freq_list]
            return freq_list
        else:
            return freq_str[:-1]
        
    def get_parameters(self):
        # return all parameters relevant for interpreting trace data
        # power, averaging, electrical delay, format, Sij,
        pass
    
    def get_bandwidth_data(self, channel = None, marker = None, convert = False):
        # gets the bandwidth data corresponsing to the set maximum
        if not channel:
            channel = self.query_channel()
        if not marker:
            marker = 1
        message = ":CALC{}:MARK{}:BWID:DATA?".format(channel,marker)
        bw_str = self.query(message)
        if convert:
            bw_list = bw_str.split(',')
            bw_list = [float(ii) for ii in bw_list]
            return bw_list
        else:
            return bw_str[:-1]
        
       
"""
*******************************************************************************
*******************************************************************************
"""      
    
class tektronix_tds7704b(gpib_instrument):
    # OSCILLOSCOPE
    
    def __init__(self, addr):
        super().__init__(addr)

    def toggle_channel(self, channel, state):
        # turns a channel on or off
        # channel = int: 1-4
        # state = 0,1 (off,on)
        message = 'sel:ch' + str(channel) + ' ' + str(state)
        self.write(message)

    def set_coupling(self, channel, coupling_string):
        # sets the coupling for a given channel
        # channel = int: 1-4
        # coupling_string = str: 'dc' or 'gnd'
        message = 'ch' + str(channel) + ':coup ' + coupling_string
        self.write(message)

    #################
    # VERTICAL AXIS #
    #################
    
    def set_vertical_offset(self, channel, offset):
        # sets the offset of a given channel in volts
        # channel = int: 1-4
        # offset = float: in volts, scientific notation OK
        message = 'ch' + str(channel) + ':offs ' + str(offset)
        self.write(message)

    def set_vertical_position(self, channel, divisions):
        # sets the vertical position in divisions
        # channel = int: 1-4
        # divisions = float
        message = 'ch' + str(channel) + ':pos ' + str(divisions)
        self.write(message)

    def set_vertical_scale(self, channel, scale):
        # sets the vertical scale in volts per division
        # channel = int: 1-4
        # scale = float: volts per division, scientific notation OK
        message = 'ch' + str(channel) + ':sca ' + str(scale)
        self.write(message)
        
    ###################
    # HORIZONTAL AXIS #
    ###################
    
    def toggle_horizontal_delay(self, state):
        # turn trigger delay on or off
        # state = 0,1 (off,on)
        message = 'hor:del:mod '+str(state)
        self.write(message)
        
    def set_horizontal_delay(self, delay):
        # sets the trigger delay in seconds
        # delay = float: in seconds, scientific notation OK
        message = 'hor:del:tim '+str(delay)
        self.write(message)
        
    def set_horizontal_position(self, position):
        # sets the horizontal position of the trigger
        # like delay, but when delay is not toggled
        # probably not too useful for actual data collection
        # position = float: percentage of window
        message = 'hor:pos '+str(position)
        self.write(message)
        
    def set_horizontal_reference(self, ref):
        # sets the reference point of the horizontal window
        #   trigger point + delay falls at this point in the horizontal window
        # ref = float: percentage of window, 0 to 100
        message = 'hor:del:pos '+str(ref)
        self.write(message)
        
    def set_horizontal_samplerate(self, rate):
        # NOT INDEPENDENT: AFFECTED BY RECORD LENGTH AND SCALE
        # sets the sample rate in Hz
        # rate = float: in Hz, scientific notation OK
        message = 'hor:mai:sampler '+str(rate)
        self.write(message)
        
    def set_horizontal_scale(self, scale):
        # NOT INDEPENDENT: AFFECTED BY SAMPLE RATE AND RECORD LENGTH
        # sets horizontal scale in seconds per division
        # scale = float: in seconds/division, scientific notation OK
        message = 'hor:sca '+str(scale)
        self.write(message)
        
    def set_horizontal_record_length(self, samples):
        # NOT INDEPENDENT: AFFECTED BY SAMPLE RATE AND SCALE
        # sets number of samples in a record
        # samples = int
        message = 'hor:reco '+str(samples)
        self.write(message)

    def get_horizontal_record_length(self):
        # query horizontal record length
        message = 'hor:reco?'
        return int(self.query(message))
    

    ########################
    # ACQUISITION SETTINGS #
    ########################
    
    def toggle_single_acquisition(self, state):
        # turn single acquisition off or on (off means continuous acquisition)
        # state = 0,1 (off,on)
        if state:
            mode_str = 'seq'
        else:
            mode_str = 'runst'
        message = 'acq:stopa '+mode_str
        self.write(message)
    
    def set_acquisition_sampling_mode(self, mode_str):
        # sets the type of sampling
        # mode_str = str:
        #   options:
        #       'rt'    - real time
        #       'it'    - interpolated time
        #       'et'    - equivalent time
        message = 'acq:samp '+mode_str
        self.write(message)

    ######################
    # FASTFRAME SETTINGS #
    ######################
    
    def toggle_fastframe(self, state):
        # turn fastframe off or on
        # state = int: 0,1 (off,on)
        message = 'hor:fast:state '+str(state)
        self.write(message)
    
    def set_frame_count(self, count):
        # sets the number of frames to acquire
        # count: int
        message = 'hor:fast:coun '+str(count)
        self.write(message)

    def set_frame_record_length(self, samples):
        # sets the number of samples in each frame
        # samples = int
        message = 'hor:fast:len '+str(samples)
        self.write(message)

    ####################
    # TRIGGER SETTINGS #
    ####################

    def force_trigger(self):
        # force a trigger event
        message = 'trig forc'
        self.write(message)

    def toggle_auto_trigger(self, state):
        # turn auto trigger off or on
        # state = int: 0,1 (off,on)
        if state:
            mode_str = 'auto'
        else:
            mode_str = 'norm'
        message = 'trig:a:mod '+mode_str
        self.write(message)

    def set_trigger_coupling(self, coupling_str):
        # sets the coupling of the main trigger
        # coupling_str = str
        # options: 'ac','dc', see manual for others
        message = 'trig:a:edge:coup '+coupling_str
        self.write(message)

    def set_trigger_slope(self, slope_str):
        # sets whether the trigger event is on a rising or falling edge
        # slope_str = str: 'rise' or 'fall'
        message = 'trig:a:edge:slo '+slope_str
        self.write(message)
        
    def set_trigger_source(self, channel):
        # sets the channel that is the source of the trigger signal
        #   NOTE: might be useful to trigger on AUX IN eventually, but this
        #           isn't currently supported by this method
        # channel = int: 1-4
        message = 'trig:a:edge:sou ch'+str(channel)
        self.write(message)
        
    def set_trigger_holdoff_mode(self, mode_str):
        # sets the type of trigger holdoff
        # mode_str = str: 'time', 'random', or 'auto'
        message = 'trig:a:hold:by '+mode_str
        self.write(message)
        
    def set_trigger_holdoff(self, delay):
        # sets the trigger holdoff in seconds
        # delay = float: scientific notation OK
        message = 'trig:a:hold:tim '+str(delay)
        self.write(message)
        
    def set_trigger_level(self, level=None):
        # sets trigger level in volts
        # level = float: volts, scientific notation OK
        # default behavior (without level passed): set to 50%
        if level:
            message = 'trig:a:lev '+str(level)
        else:
            message = 'trig:a'
        self.write(message)
        
    ################
    # INPUT/OUTPUT #
    ################

    def set_data_encoding(self, encoding_str):
        # sets the format of the output data
        # encoding_str = str
        #   options:
        #       'ascii'
        #       'rib' - a binary data format
        #       see manual for others
        message = 'dat:enc '+encoding_str
        self.write(message)

    def set_data_source(self, channel):
        # sets the channel from which output data will be obtained
        # channel = int: 1-4
        message = 'dat:sou ch'+str(channel)
        self.write(message)

    def get_data(self, channel=None, encoding_str=None):
        # gets data from a given channel
        # if no channel given, assume the current source channel
        # if no encoding given, assume current encoding
        if channel:
            self.set_data_source(channel)
        if encoding_str:
            self.set_data_encoding(encoding_str)
        message = 'curv?'
        return self.query(message)

    def set_data_start(self,sample):
        # sets the first sample to transfer with curv?
        # sample = int: index from 1
        message = 'dat:start '+str(sample)
        self.write(message)
        
    def set_data_stop(self,sample):
        # sets the last sample to transfer with curv?
        # sample = int
        message = 'dat:stop '+str(sample)
        self.write(message)
        
    def get_data_in_chunks(self, chunk_size=None, start=1, stop=None, \
                           channel = None, encoding_str = 'ascii'):
        # get a large data set chunk by chunk and then concatenate
        if channel:
            self.set_data_source(channel)
        self.set_data_encoding(encoding_str)
        if not stop:
            stop = self.get_horizontal_record_length()
        n_samples = stop-start
        if not chunk_size:
            chunk_size = round(n_samples/1000)
            if chunk_size == 0:
                chunk_size = 1
                
        samples_left = n_samples
        data = ''
        first_sample = start
        last_sample = chunk_size
        self.set_data_start(first_sample)
        self.set_data_stop(last_sample)
        while samples_left > 0:
            new_data = self.get_data()
            data += new_data[:-1]+','
            self.wai()
            first_sample += chunk_size
            samples_left -= chunk_size
            if samples_left < chunk_size:
                last_sample = stop
            else:
                last_sample += chunk_size
            self.set_data_start(first_sample)
            self.set_data_stop(last_sample)
            self.wai()
                        
        data = self.convert_ascii_to_volts(data[:-1])
        return data
            
    def get_ymult(self):
        # get factor for scaling ascii to volts
        message = 'wfmo:ymu?'
        return float(self.query(message))
    
    def get_ypos(self):
        # get y position for converting ascii to volts
        message = 'wfmo:yof?'
        return float(self.query(message))
    
    def get_yoff(self):
        # get y offset for converting ascii to volts
        message = 'wfmo:yze?'
        return float(self.query(message))
    
    def convert_ascii_to_volts(self,data_str):
        # convert ascii data to voltages
        data_list = data_str.split(',')
        ymult = self.get_ymult()
        yoff = self.get_yoff()
        ypos = self.get_ypos()
        data_list =[str(((float(item)-ypos)*ymult)+yoff) for item in data_list]
        new_data_str = ','.join(data_list)
        return new_data_str
    
    def get_tscale(self):
        # get the time scale associated with each sample
        message = 'wfmo:xin?'
        return float(self.query(message))
    
    def get_tstart(self):
        # get the starting time of the data
        message = 'wfmo:xzero?'
        return float(self.query(message))
    
    def get_time_axis(self):
        # get the time axis of an acquisition
        # return as csv string
        tscale = self.get_tscale()
        tstart = self.get_tstart()
        records = self.get_horizontal_record_length()
        total_time = tscale*records
        tstop = tstart + total_time
        time_axis = list(np.linspace(tstart,tstop,num=records))
        time_axis = [str(item) for item in time_axis]
        time_str = ','.join(time_axis)
        return time_str
        

"""
*******************************************************************************
*******************************************************************************
"""    

class tektronix_awg7052(gpib_instrument):
    
    def __init__(self,addr):
        super().__init__(addr)

    def force_trigger(self):
        # force a trigger event
        message = '*trg'
        self.write(message)

    def toggle_output(self,channel,state):
        # toggles whether channel is on or off
        # channel   = int: 1 or 2
        # state = 0,1 (off,on)
        message = 'outp'+str(channel)+' '+str(state)
        self.write(message)

    def run(self):
        # initiate output (turns 'run' on)
        message = 'awgc:run'
        self.write(message)
        
    def stop(self):
        # stops output (turns 'run' off)
        message = 'awgc:stop'
        self.write(message)
        
    def toggle_run(self,state):
        # toggles whether the awg is running or not
        # state = 0,1 (off,on)
        if state:
            self.run()
        else:
            self.stop()

    def set_sampling_rate(self,rate):
        # rate = float, scientific notation x.yEz works too, flexible
        message = 'sour:freq '+str(rate)
        self.write(message)

    def set_run_mode(self,mode):
        # sets the run mode, obviously
        # mode = str
        #   -options:
        #       'cont'  -continuous, repeats waveform indefinitely
        #       'trig'  -triggered, outputs one waveform for each trigger
        #       'gat'   -gated, see manual
        #       'seq'   -sequence, see manual
        #       'enh'   -enhanced, see manual
        message = 'awgc:rmode '+mode
        self.write(message)

    def set_frequency_reference(self,source):
        # sets whether to use the internal or an external 10mhz reference
        # source = str: 'int' or 'ext'
        message = 'sour:rosc:sour '+source
        self.write(message)

    #########################
    # WAVEFORM MANIPULATION #
    #########################

    def new_waveform(self, name, size):
        # creates a new (empty) waveform in the current waveform list
        # name = str
        # size = int, number of samples in waveform
        message = 'wlis:wav:new '+'"'+name+'"'+','+str(size)+',REAL'
        self.write(message)
        
    def delete_waveform(self,name):
        # deletes a waveform from the current list
        # name = str
        message = 'wlis:wav:del '+'"'+name+'"'
        self.write(message)
        
    def clear_waveforms(self):
        # deletes all user-defined waveforms from the current list
        message = 'wlis:wav:del ALL'
        self.write(message)
        
    def send_waveform(self, name, w, m1, m2):
        # name = str
        # w = numpy array of floats
        # m1, m2 = numpy arrays containing only 0s and 1s
        n_samples = len(w)
        self.delete_waveform(name)
        self.new_waveform(name,n_samples)
        
        m = ((2**7)*m2) + ((2**6)*m1)
        
        bytes_data = b''
        for ii in range(n_samples):
            bytes_data += struct.pack('fB',w[ii],int(m[ii]))
        
        num_bytes = n_samples*5
        num_bytes = str(num_bytes)
        num_digits = str(len(num_bytes))
        
        num_bytes = num_bytes.encode('ascii')
        num_digits = num_digits.encode('ascii')
        bytes_count = num_digits + num_bytes
        
        bytes_name = name.encode('ascii')
        bytes_samples = str(n_samples)
        bytes_samples = bytes_samples.encode('ascii')
        
        message = b'wlis:wav:data "'+bytes_name+b'",0,'+bytes_samples \
                    +b',#'+bytes_count+bytes_data
        self.write_raw(message)

    def load_waveform(self,channel,name):
        # load a waveform from the current list into channel
        # channel = int: 1 or 2
        # name = str
        message = 'sour'+str(channel)+':wav "'+name+'"'
        self.write(message)

    #######################
    # ANALOG OUT SETTINGS #
    #######################
        
    def set_analog_amplitude(self, channel, amplitude, units='V'):
        # sets the peak to peak voltage amplitude of a given channel
        # channel   = int: 1 or 2
        # amplitude = float
        # units     = str: 'V', 'mV'
        message = 'sour'+str(channel)+':volt '+str(amplitude)+units
        self.write(message)

    ###################
    # MARKER SETTINGS #
    ###################

    def set_marker_low(self, channel, marker, voltage, units='V'):
        # set the low voltage of a given marker on a given channel
        # channel   = int: 1 or 2
        # marker    = int: 1 or 2
        # voltage = float
        # units     = str: 'V', 'mV'
        message = 'sour'+str(channel)+':mark'+str(marker)+':volt:low ' \
                    + str(voltage)+units
        self.write(message)

    def set_marker_high(self, channel, marker, voltage, units='V'):
        # set the low voltage of a given marker on a given channel
        # channel   = int: 1 or 2
        # marker    = int: 1 or 2
        # voltage = float
        # units     = str: 'V', 'mV'
        message = 'sour'+str(channel)+':mark'+str(marker)+':volt:high ' \
                    + str(voltage)+units
        self.write(message)

    def set_marker_delay(self, channel, marker, delay):
        # set the delay for a given channel and marker
        # NOTE: ONLY ACCEPTS DELAY IN PICOSECONDS
        # channel   = int: 1 or 2
        # marker    = int: 1 or 2
        # delay     = float: time in picoseconds
        # units     = str
        message = 'sour' + str(channel) + ':mark' + str(marker) + ':del ' \
                    + str(delay) + 'ps'
        self.write(message)

    #######################
    # FILE SYSTEM METHODS #
    #######################

    def query_cwd(self):
        # return current working directory of awg mass memory
        message = 'mmem:cdir?'
        return self.query(message)

    def mkdir(self, dir_name):
        # makes a new directory in the current working directory
        # dir_name = str
        message = 'mmem:mdir '+'"'+dir_name+'"'
        self.write(message)
        
    def ls(self):
        # query contents of current working directory
        return self.query('mmem:cat?')

    def cd(self,rel_path):
        # change directory relative to current directory
        # rel_path = str
        # '..' moves up one level
        message = 'mmem:cdir '+'"'+rel_path+'"'
        self.write(message)

    def reset_cwd(self):
        message = 'mmem:cdir'
        self.write(message)

    def set_cwd(self,absolute_path):
        # set current working directory of awg mass memory to absolute_path
        # path = str
        self.reset_cwd()
        message = 'mmem:cdir '+'"'+absolute_path+'"'
        self.write(message)

"""
*******************************************************************************
*******************************************************************************
"""

class tektronix_awg520(gpib_instrument):
    
    def __init__(self,addr):
        super().__init__(addr)

    def force_trigger(self):
        # force a trigger event
        message = '*trg'
        self.write(message)

    def toggle_output(self,channel,state):
        # toggles whether channel is on or off
        # channel   = int: 1 or 2
        # state = 0,1 (off,on)
        message = 'outp'+str(channel)+' '+str(state)
        self.write(message)
    
    def set_run_mode(self,mode):
        # sets the run mode, obviously
        # mode = str
        #   -options:
        #       'cont'  -continuous, repeats waveform indefinitely
        #       'trig'  -triggered, outputs one waveform for each trigger
        #       'gat'   -gated, see manual
        #       'enh'   -enhanced, see manual
        message = 'awgc:rmode '+mode
        self.write(message)
        
    def run(self):
        # initiate output (turns 'run' on)
        message = 'awgc:run'
        self.write(message)
        
    def stop(self):
        # stops output (turns 'run' off)
        message = 'awgc:stop'
        self.write(message)
        
    def toggle_run(self,state):
        # toggles whether the awg is running or not
        # state = 0,1 (off,on)
        if state:
            self.run()
        else:
            self.stop()
            
    def set_offset(self, channel, offset, units='V'):
        # set the voltage offset of a given channel
        # channel   = int: 1 or 2
        # offset = float
        # units     = str: 'V', 'mV'
        message = 'sour'+str(channel)+':volt:offs '+str(offset)+units
        self.write(message)
        
    def set_amplitude(self, channel, amplitude, units='V'):
        # sets the peak to peak voltage amplitude of a given channel
        # channel   = int: 1 or 2
        # amplitude = float
        # units     = str: 'V', 'mV'
        message = 'sour'+str(channel)+':volt '+str(amplitude)+units
        self.write(message)

    def set_frequency_reference(self, channel, source):
        # sets the source of the 10MHz reference to use with a channel
        # channel = int: 1 or 2
        # source = str: 'int' or 'ext'
        message = 'sour'+str(channel)+':rosc:sour '+source
        self.write(message)

    def send_waveform(self, w, m1, m2, filename, samplerate):
        # w = array of floats between -1 and 1
        # m1 = array of integers, only 0 and 1 allowed
        # m2 as m1
        # filename = str ending in .wfm
        # samplerate = float up to 1.0E9 (one GS/s)
        
        # need to send a block of bytes data with format
        # cmd + file_counter + header + data_counter + bytes_data + trailer
        
        #converting strings to raw bytes data
        bytes_filename = filename.encode('ascii')
        cmd = b'mmem:data "'+bytes_filename+b'",'
        header = b'MAGIC 1000\r\n'
        
        nsamples = len(w)
        m = m1 + np.multiply(m2, 2)
        #packing waveform and markers into bytes data
        bytes_data = b''
        for ii in range(nsamples):
            bytes_data += struct.pack('<fB',w[ii],int(m[ii]))
        
        num_bytes = len(bytes_data)
        num_digits = len(str(num_bytes))
        
        # converting this info into bytes
        num_digits = str(num_digits)
        num_digits = num_digits.encode('ascii')
        num_bytes = str(num_bytes)
        num_bytes = num_bytes.encode('ascii')
        data_counter = b'#'+num_digits+num_bytes
        
        samplerate_str = '{:.2E}'.format(samplerate)
        samplerate_bytes = samplerate_str.encode('ascii')
        trailer = b'CLOCK '+samplerate_bytes+b'\r\n'
        
        file = header + data_counter + bytes_data + trailer
        
        num_file_bytes = len(file)
        num_file_digits = len(str(num_file_bytes))
        
        num_file_digits = str(num_file_digits)
        num_file_digits = num_file_digits.encode('ascii')
        num_file_bytes = str(num_file_bytes)
        num_file_bytes = num_file_bytes.encode('ascii')
        file_counter = b'#'+num_file_digits+num_file_bytes
        
        message = cmd + file_counter + file
        self.write_raw(message)
        

    def load_waveform(self, channel, filename):
        # loads the waveform at filename into a channel
        # channel   = int: 1 or 2
        # filename  = str
        #   - can be a path, always relative to current working directory
        message = 'sour'+str(channel)+':func:user '+'"'+filename+'"'
        self.write(message)

    ###################
    # MARKER SETTINGS #
    ###################

    def set_marker_low(self, channel, marker, voltage, units='V'):
        # set the low voltage of a given marker on a given channel
        # channel   = int: 1 or 2
        # marker    = int: 1 or 2
        # voltage = float
        # units     = str: 'V', 'mV'
        message = 'sour'+str(channel)+':mark'+str(marker)+':volt:low ' \
                    + str(voltage)+units
        self.write(message)

    def set_marker_high(self, channel, marker, voltage, units='V'):
        # set the low voltage of a given marker on a given channel
        # channel   = int: 1 or 2
        # marker    = int: 1 or 2
        # voltage = float
        # units     = str: 'V', 'mV'
        message = 'sour'+str(channel)+':mark'+str(marker)+':volt:high ' \
                    + str(voltage)+units
        self.write(message)

    def set_marker_delay(self, channel, marker, delay, units='s'):
        # set the delay for a given channel and marker
        # channel   = int: 1 or 2
        # marker    = int: 1 or 2
        # delay     = float
        # units     = str
        message = 'sour' + str(channel) + ':mark' + str(marker) + ':del ' \
                    + str(delay) + units
        self.write(message)

    #######################
    # FILE SYSTEM METHODS #
    #######################

    def set_mass_storage(self, device = 'MAIN'):
        #sets mass storage device
        #device = str
        #options: 'MAIN', 'FLOP', 'NET1', 'NET2', 'NET3'
        message = 'mmem:msis '+'"'+device+'"'
        self.write(message)

    def query_cwd(self):
        # return current working directory of awg mass memory
        message = 'mmem:cdir?'
        return self.query(message)

    def mkdir(self, dir_name):
        # makes a new directory in the current working directory
        # dir_name = str
        message = 'mmem:mdir '+'"'+dir_name+'"'
        self.write(message)
        
    def ls(self):
        # query contents of current working directory
        return self.query('mmem:cat?')

    def cd(self,rel_path):
        # change directory relative to current directory
        # rel_path = str
        # '..' moves up one level
        message = 'mmem:cdir '+'"'+rel_path+'"'
        self.write(message)

    def reset_cwd(self):
        message = 'mmem:cdir'
        self.write(message)

    def set_cwd(self,absolute_path):
        # set current working directory of awg mass memory to absolute_path
        # path = str
        self.reset_cwd()
        message = 'mmem:cdir '+'"'+absolute_path+'"'
        self.write(message)
        
    #######################
    # FG MODE #
    #######################
    def toggle_fg_mode(self, state):
        """
        Sets the function generator mode on or off.
        Parameters:
            state (int): 0 or 1 for on and off.
        """
        message = 'awgc:fg {}'.format(state)
        self.write(message)
        
    def set_fg_frequency(self, freq_Hz):
        """
        Sets the frequency of the funtion waveform.
        Parameters:
            freq_Hz (float): Value of frequency in Hz.
        """
        message = 'awgc:fg:freq {}'.format(freq_Hz)
        self.write(message)
        
    def set_fg_shape(self, channel, shape_str):
        """
        Sets the shape of the funtion waveform.
        Parameters:
            channel (int): Channel number.
            shape_str (str): Standard waveform shape - 
                            'SIN', 'TRI', 'SQU', 'RAMP', 'PULS'
        """
        message = 'awgc:fg{}:func {}'.format(channel, shape_str)
        self.write(message)
        
    def set_fg_duty_cycle(self, channel, value):
        """
        Sets the duty cycle of the pulse waveform on the fg mode.
        Parameters:
            channel (int): Channel number.
            value (float): Value of the duty cycle from 0.1 to 99.9%.
        """
        message = 'awgc:fg{}:puls:dcyc {}'.format(channel, value)
        self.write(message)
        
    def set_fg_voltage_pp(self, channel, value):
        """
        Sets the voltage amplitude of the function waveform.
        (NOTE: This is set for the 50 ohm impedance. Higher impedance load
        will give double the output voltage.)
        Parameters:
            channel (int): Channel number.
            value (float): Value of the peak to peak voltage (V).
                            Range from 20 mV to 2 V only.
        """
        message = 'awgc:fg{}:volt {}'.format(channel, value)
        self.write(message)
        
    def set_fg_voltage_offset(self, channel, value):
        """
        Sets the offset voltage of the function waveform.
        Parameters:
            channel (int): Channel number.
            value (float): Value of the offset voltage (V).
        """
        message = 'awgc:fg{}:volt:offs {}'.format(channel, value)
        self.write(message)
    
        
"""
*******************************************************************************
*******************************************************************************
"""

class national_instruments_bnc2090:

    def __init__(self):
        self.output_voltages = []
        for ii in range(2):
            task = nidaqmx.Task()
            ch_str = 'Dev1/ao' + str(ii)
            task.ao_channels.add_ao_voltage_chan(ch_str)
            self.output_voltages.append(0.0)
            task.write(0.0)
            task.close()

    def set_voltage(self, output_ind, voltage):
        task = nidaqmx.Task()
        ch_str = 'Dev1/ao' + str(output_ind)
        task.ao_channels.add_ao_voltage_chan(ch_str)
        self.output_voltages[output_ind] = voltage
        task.write(voltage)
        task.close()

    def get_voltage(self, input_ind, samples=None):
        task = nidaqmx.Task()
        ch_str = 'Dev1/ai' + str(input_ind)
        task.ai_channels.add_ai_voltage_chan(ch_str)
        if samples:
            result = task.read(samples)
        else:
            result = task.read()
        task.close()
        return result

    def get_mean_voltage(self, input_ind, samples, return_stdev=False):
        voltages = self.get_voltage(input_ind, samples)
        mean = statistics.mean(voltages)
        if return_stdev:
            stdev = statistics.stdev(voltages)
            return (
             mean, stdev)
        else:
            return mean

    def gui(self):
        window = tk.Tk()
        window.geometry('710x450')
        window.title('Analog Outputs')
        window.focus_force()
        
        frame0 = tk.Frame(window,width=800,height=300,bd=3,relief='ridge')
        frame0.grid(row=0,column=0,padx=10,pady=10)
        
        ao0 = tk.DoubleVar()
        ao0.set(0.0)
        
        ao0_coarse = tk.DoubleVar()
        ao0_coarse.set(0.0)
        
        ao0_fine = tk.DoubleVar()
        ao0_fine.set(0.0)
        
        def update_voltage_ao0(*args):
            coarse_voltage = ao0_coarse.get()
            fine_voltage = ao0_fine.get()*0.001
            voltage = round(coarse_voltage+fine_voltage,3)
            ao0.set(voltage)
            self.set_voltage(0,voltage)
        
        def coarse_increment_ao0():
            coarse_val = ao0_coarse.get()+0.1
            ao0_coarse.set(coarse_val)
            update_voltage_ao0()
            
        def coarse_decrement_ao0():
            coarse_val = ao0_coarse.get()-0.1
            ao0_coarse.set(coarse_val)
            update_voltage_ao0()
            
        def fine_increment_ao0():
            fine_val = ao0_fine.get()+1
            ao0_fine.set(fine_val)
            update_voltage_ao0()
            
        def fine_decrement_ao0():
            fine_val = ao0_fine.get()-1
            ao0_fine.set(fine_val)
            update_voltage_ao0()
        
        
        ao0_lbl = tk.Label(frame0, text = 'Analog Out 0', font='helvetica 16 bold')
        ao0_lbl.grid(row=0,column=0, columnspan=2, sticky= 's', padx=10,pady=10)
        
        ao0_box = tk.Entry(frame0, textvariable=ao0,\
                               width=10)
        ao0_box.grid(row=1,column=0, sticky = 'n', pady=10)
        
        box_label = tk.Label(frame0, text = 'Volts')
        box_label.grid(row=1,column=1,sticky='wn',pady=10)
        
        out_coarse_scale = tk.Scale(frame0, from_=-10, to=10, orient='horizontal', \
                             resolution = 0.1, command = update_voltage_ao0, \
                             digits = 3, label = 'V', length = 400, \
                             tickinterval = 2.5, variable = ao0_coarse)
        out_coarse_scale.grid(row=0,column=2,padx=10,pady=10)
        
        out_fine_scale = tk.Scale(frame0, from_=-100, to=100, orient='horizontal', \
                             resolution = 1, command = update_voltage_ao0, \
                             digits = 3, label = 'mV', length = 400, \
                             tickinterval = 25, variable = ao0_fine)
        out_fine_scale.grid(row=1,column=2,padx=10,pady=10)
        
        
        
        coarse_dec_button = tk.Button(frame0,text='-',command=coarse_decrement_ao0,\
                                      font = 'helvetica 16')
        coarse_dec_button.grid(row=0,column=3,sticky='ew',padx=10,pady=10)
        coarse_inc_button = tk.Button(frame0,text='+',command=coarse_increment_ao0,\
                                      font = 'helvetica 16')
        coarse_inc_button.grid(row=0,column=4,sticky='ew',padx=10,pady=10)
        
        
        fine_dec_button = tk.Button(frame0,text='-',command=fine_decrement_ao0,\
                                      font = 'helvetica 16')
        fine_dec_button.grid(row=1,column=3,sticky='ew',padx=10,pady=10)
        fine_inc_button = tk.Button(frame0,text='+',command=fine_increment_ao0,\
                                      font = 'helvetica 16')
        fine_inc_button.grid(row=1,column=4,sticky='ew',padx=10,pady=10)
        
        ao0.trace('w',update_voltage_ao0)
        
        """
        ************************************************
        """
        
        frame1 = tk.Frame(window,width=800,height=300,bd=3,relief='ridge')
        frame1.grid(row=1,column=0,padx=10,pady=10)
        
        ao1 = tk.DoubleVar()
        ao1.set(0.0)
        
        ao1_coarse = tk.DoubleVar()
        ao1_coarse.set(0.0)
        
        ao1_fine = tk.DoubleVar()
        ao1_fine.set(0.0)
        
        def update_voltage_ao1(*args):
            coarse_voltage = ao1_coarse.get()
            fine_voltage = ao1_fine.get()*0.001
            voltage = round(coarse_voltage+fine_voltage,3)
            ao1.set(voltage)
            self.set_voltage(1,voltage)
        
        def coarse_increment_ao1():
            coarse_val = ao1_coarse.get()+0.1
            ao1_coarse.set(coarse_val)
            update_voltage_ao1()
            
        def coarse_decrement_ao1():
            coarse_val = ao1_coarse.get()-0.1
            ao1_coarse.set(coarse_val)
            update_voltage_ao1()
            
        def fine_increment_ao1():
            fine_val = ao1_fine.get()+1
            ao1_fine.set(fine_val)
            update_voltage_ao1()
            
        def fine_decrement_ao1():
            fine_val = ao1_fine.get()-1
            ao1_fine.set(fine_val)
            update_voltage_ao1()
        
        
        ao1_lbl = tk.Label(frame1, text = 'Analog Out 1', font='helvetica 16 bold')
        ao1_lbl.grid(row=3,column=0, columnspan=2, sticky= 's', padx=10,pady=10)
        
        ao1_box = tk.Entry(frame1, textvariable=ao1,\
                               width=10)
        ao1_box.grid(row=4,column=0, sticky = 'n', pady=10)
        
        box_label = tk.Label(frame1, text = 'Volts')
        box_label.grid(row=4,column=1,sticky='wn',pady=10)
        
        out_coarse_scale = tk.Scale(frame1, from_=-10, to=10, orient='horizontal', \
                             resolution = 0.1, command = update_voltage_ao1, \
                             digits = 3, label = 'V', length = 400, \
                             tickinterval = 2.5, variable = ao1_coarse)
        out_coarse_scale.grid(row=3,column=2,padx=10,pady=10)
        
        out_fine_scale = tk.Scale(frame1, from_=-100, to=100, orient='horizontal', \
                             resolution = 1, command = update_voltage_ao1, \
                             digits = 3, label = 'mV', length = 400, \
                             tickinterval = 25, variable = ao1_fine)
        out_fine_scale.grid(row=4,column=2,padx=10,pady=10)
        
        
        
        coarse_dec_button = tk.Button(frame1,text='-',command=coarse_decrement_ao1,\
                                      font = 'helvetica 16')
        coarse_dec_button.grid(row=3,column=3,sticky='ew',padx=10,pady=10)
        coarse_inc_button = tk.Button(frame1,text='+',command=coarse_increment_ao1,\
                                      font = 'helvetica 16')
        coarse_inc_button.grid(row=3,column=4,sticky='ew',padx=10,pady=10)
        
        
        fine_dec_button = tk.Button(frame1,text='-',command=fine_decrement_ao1,\
                                      font = 'helvetica 16')
        fine_dec_button.grid(row=4,column=3,sticky='ew',padx=10,pady=10)
        fine_inc_button = tk.Button(frame1,text='+',command=fine_increment_ao1,\
                                      font = 'helvetica 16')
        fine_inc_button.grid(row=4,column=4,sticky='ew',padx=10,pady=10)
        
        ao1.trace('w',update_voltage_ao1)
        
        window.mainloop()
        
"""
*******************************************************************************
*******************************************************************************
"""

class alazartech_ats9462():    
    def __init__(self, boardId = 1):
        self.board = ats.Board(systemId = 1, boardId = 1)
    
    def set_timebase(self, source, sample_rate = 180e6):
        """
        Sets the timebase to internal, external or external 10 MHz clock.
        Parameters:
            source (str): Selects the source.
                'int'/'fast_ext'/'ext_10'/'slow_ext'
            sample_rate (float): The sample rate in Hz. Note only selected
                values can be used. Check manual.
        """
        if source == 'int':
            source_value = ats.INTERNAL_CLOCK
            sample_rate_value = eval('ats.SAMPLE_RATE_{}MSPS'.format(
                    int(sample_rate*1e-6)))
        elif source == 'ext_10':
            source_value = ats.EXTERNAL_CLOCK_10MHz_REF
            sample_rate_value = eval('ats.SAMPLE_RATE_{}MSPS'.format(
                    int(sample_rate)))
        elif source == 'fast_ext':
            source_value = ats.FAST_EXTERNAL_CLOCK
            sample_rate_value = ats.SAMPLE_RATE_USER_DEF
        elif source == 'slow_ext':
            source_value = ats.SLOW_EXTERNAL_CLOCK
            sample_rate_value = ats.SAMPLE_RATE_USER_DEF
        self.board.setCaptureClock(source_value, sample_rate_value, 
                                   ats.CLOCK_EDGE_RISING,0)
        
    def set_input_params(self, channel, input_range_V,
                         coupling = 'dc',
                         impedance = 1e6):
        """
        Sets the parameters for input - desired input range, termination,
        coupling.
        Parameters:
            channel (int): selects the channel 1/2
            input_range (float): specifies the input range. 
                Available options for 50 Ohm - 0.2,0.4,0.8,2,4
                Available options for 1 MOhm - 0.2,0.4,0.8,2,4,8,16
            coupling (string): AC/DC
            impedance (int): Sets the termination - 50/75/300/1e6
        """
        channel_value = eval('ats.CHANNEL_{}'.format(chr(64+channel)))
        coupling_value = eval('ats.{}_COUPLING'.format(coupling.upper()))
        if impedance == 1e6:
            impedance_value = ats.IMPEDANCE_1M_OHM
        else:
            impedance_value = eval('ats.IMPEDANCE_{}_OHM'.format(impedance))
        if input_range_V < 1:
            input_range_value = eval('ats.INPUT_RANGE_PM_{}_MV'.format(
                    int(input_range_V*1e3)))
        else:
            input_range_value = eval('ats.INPUT_RANGE_PM_{}_V'.format(
                    int(input_range_V)))
        self.board.inputControl(channel_value, coupling_value,
                                input_range_value, impedance_value)
    
    def toggle_bandwidth_filter(self, channel, state = False):
        """
        Toggles the low pass filter that attenuate signal about 20 MHz.
        Parameters:
            channel (int): Selects the channel.
            state (bool): Enables/Disables LPF.
        """
        self.board.setBWLimit(channel, int(state))
        
    def set_trigger_operation(self, source_1, slope_1, fraction_1,
                              source_2 = None, slope_2 = None,
                              fraction_2 = None):
        """
        Configures the two trigger engines. If source_2 is None, only J
        tigger engine is used. If a (J or K) option is required, set the
        parameters for labels _2.
        Parameters:
            source_1 (int/str): Trigger source to use for J (1/2/'ext'/'dis')
                'dis' disables this engine.
            slope_1 (str): Selects the rising/falling edge to activate engine
                J. ('pos'/'neg')
            fraction_1 (float): Sets the trigger level configuration for engine
                J. Value is written as the fraction of the full scale input 
                range. If input range is y mV, required level is x mV, 
                the fraction is set as x/y.
            source_2 (int/str): Trigger source to use for K (1/2/'ext'/'dis')
                'dis' disables this engine. By default, set to None which
                disables this engine.
            slope_2 (str): Selects the rising/falling edge to activate engine
                K. ('pos'/'neg')
            fraction_2 (float): Sets the trigger level configuration for engine
                K. Value is written as the fraction of the full scale input 
                range. If input range is y mV, required level is x mV, 
                the fraction is set as x/y.
        """
        if source_2:
            trigger_op_config = ats.TRIG_ENGINE_OP_J_OR_K
            if source_2 == 'ext':
                source_2_value = ats.TRIG_EXTERNAL
            else:
                source_2_value = eval('ats.TRIG_CHAN_{}'.format(
                        chr(source_2+64)))
            if slope_2 == 'pos':
                slope_2_value = ats.TRIGGER_SLOPE_POSITIVE
            else:
                slope_2_value = ats.TRIGGER_SLOPE_NEGATIVE
            trigger_level_2 = int(128*(1+fraction_2))
        else:
            trigger_op_config = ats.TRIG_ENGINE_OP_J
            source_2_value = ats.TRIG_DISABLE
            slope_2_value = ats.TRIGGER_SLOPE_POSITIVE
            trigger_level_2 = 128
        if source_1 == 'ext':
            source_1_value = ats.TRIG_EXTERNAL
        elif source_1 == 'dis':
            source_1_value = ats.TRIG_DISABLE
        else:
            source_1_value = eval('ats.TRIG_CHAN_{}'.format(chr(source_1+64)))
        if slope_1 == 'pos':
                slope_1_value = ats.TRIGGER_SLOPE_POSITIVE
        else:
            slope_1_value = ats.TRIGGER_SLOPE_NEGATIVE
        trigger_level_1 = int(128*(1+fraction_1))
        self.board.setTriggerOperation(trigger_op_config,
                                       ats.TRIG_ENGINE_J,
                                       source_1_value,
                                       slope_1_value,
                                       trigger_level_1,
                                       ats.TRIG_ENGINE_K,
                                       source_2_value,
                                       slope_2_value,
                                       trigger_level_2)
        
    def set_external_trigger_params(self, coupling = 'ac', 
                                    trigger_range = ats.ETR_1V):
        """
        Sets the external trigger range and coupling.
        Parameters:
            coupling (str): Selects the coupling to AC/DC - 'ac'/'dc'
            range (int): Sets the external trigger range identifier.
                Options are 'ats.ETR_5V','ats.ETR_1V', ETR_TTL, 'ats.ETR_2V5'
        """
        coupling_value = eval('ats.{}_COUPLING'.format(coupling.upper()))
        self.board.setExternalTrigger(coupling_value, trigger_range)
        
    def set_trigger_timeout(self, wait_sec):
        """
        Sets the amount of time the board will wait for a hardware trigger to
        occur before automatically generating a software trigger event. If
        set to zero, the board will wait forever for a trigger event.
        Parameters:
            wait_sec (float): Amount of time to wait in sec.
        """
        self.board.setTriggerTimeOut(int(wait_sec/10e-6))
        
    def set_trigger_delay(self, wait_sec):
        """
        Sets the amount of time the board will wait after it receives a
        trigger before capturing a record.
        Parameters:
            wait_sec (float): Amount of time to wait in sec.
        """
        self.board.setTriggerDelay(int(wait_sec*180e6))
    
    def set_record_size(self, pretrigger, posttrigger):
        """
        Sets the number of pre-trigger and post-trigger samples per record.
        Parameters:
            pretrigger (int): Number of pre-trigger samples.
            posttrigger (int): Number of post-trigger samples.
        """
        self.board.setRecordSize(pretrigger, posttrigger)
        
    def start_capture(self):
        """
        Arm a board to start an acquisition.
        """
        self.board.startCapture()
        
    def force_trigger(self):
        """
        Forces a trigger event.
        """
        self.board.forceTrigger()
        
    def configure_board(self, input_range_V1, input_range_V2, 
                        trigger_source_1 = 'ext', trigger_source_2 = None,
                        trigger_slope_1 = 'pos', trigger_slope_2 = None,
                        trigger_fraction_1 = 0.2, trigger_fraction_2 = None,
                        ref_source = 'fast_ext', sample_rate = None,
                        coupling_1 = 'dc', coupling_2 = 'dc',
                        impedance_1 = 1e6, impedance_2 = 1e6, bw_1 = False,
                        bw_2 = False, trigger_coupling = 'dc', 
                        trigger_range = ats.ETR_1V, delay_sec = 0,
                        trigger_timeout_sec = 10e-6):
        self.set_timebase(ref_source, sample_rate)
        self.set_input_params(1, input_range_V1,
                              coupling_1, impedance_1)
        self.set_input_params(2, input_range_V2,
                              coupling_2, impedance_2)
        self.toggle_bandwidth_filter(1, bw_1)
        self.toggle_bandwidth_filter(2, bw_2)
        self.set_trigger_operation(trigger_source_1, trigger_slope_1, 
                                   trigger_fraction_1, trigger_source_2,
                                   trigger_slope_2, trigger_fraction_2)
        self.set_external_trigger_params(trigger_coupling, trigger_range)
        self.set_trigger_delay(delay_sec)
        self.set_trigger_timeout(trigger_timeout_sec)

    def acquire(self, channels, pretrigger, posttrigger, 
                recordsperbuffer, buffersperacq,
                timeout_sec = 1000, input_range_V=4,
                clock_rate=1e6, sample_rate=1e3):
        """
        Dual port AutoDMA acquisition. Traditional AutoDMA.
        Samples->Record->Buffer->Acquisition.
        Buffer Organization for both channels - R1A, R1B, R2A, R2B, R3A, R3B..
        """
        samplesperrecord = pretrigger + posttrigger
        memorySize_samples, bitsPerSample = self.board.getChannelInfo()
        bytesPerSample = (bitsPerSample.value + 7) // 8
        bytesPerRecord = bytesPerSample * samplesperrecord
        channelcount = int(channels/3)+1
        bytesperbuffer = bytesPerRecord * recordsperbuffer * channelcount
        bufferCount = 4 #Makes 4 buffers available. Use a minimum of 2.
#        buffer_array=[]
        if (channels==1) or (channels==2):
            ch=1
            data_1={}
        else:
            ch=2
            data_1={}
            data_2={}
        dn = int(clock_rate/sample_rate)

        # Allocate DMA buffers
    
        sample_type = ctypes.c_uint8
        if bytesPerSample > 1:
            sample_type = ctypes.c_uint16
    
        buffers = []
        for i in range(bufferCount):
            buffers.append(ats.DMABuffer(self.board.handle, 
                                         sample_type, bytesperbuffer))
        self.set_record_size(pretrigger, posttrigger)
        recordsperacq = recordsperbuffer * buffersperacq
        #Configure the board to make an AutoDMA acquisition.
        self.board.beforeAsyncRead(channels,
                                   -pretrigger,
                                   samplesperrecord,
                                   recordsperbuffer,
                                   recordsperacq,0)
        for buffer in buffers:
            #Make buffers available to be filled by board.
            self.board.postAsyncBuffer(buffer.addr, buffer.size_bytes)
            
        start = time.clock() # Keep track of when acquisition started
        try:
            self.board.startCapture() # Arm the board to begin the acquisition.
            print("Capturing %d buffers. Press <enter> to abort" %
                  buffersperacq)
            buffersCompleted = 0
            bytesTransferred = 0
            while (buffersCompleted < buffersperacq and not
                   ats.enter_pressed()):
                # Wait for the buffer at the head of the list of available
                # buffers to be filled by the board.
                buffer = buffers[buffersCompleted % len(buffers)]
                self.board.waitAsyncBufferComplete(buffer.addr, 
                                                   int(timeout_sec*1e3))
#                buffer_array.append(buffer)
                buffersCompleted += 1
                bytesTransferred += buffer.size_bytes
                value_array = buffer.buffer
                check_start=time.time()
#                volt_array = self.code_to_volt(value_array, input_range_V)
                volt_array = 2*input_range_V*(value_array/(2**16-1)-0.5)
                check_stop=time.time()
                print('Time elapsed is {}'.format(check_stop-check_start))
#                volt_array = value_array
                n_per_chan = len(volt_array)//ch
                data_1[str(buffersCompleted)] = []
                data_1[str(buffersCompleted)].append(volt_array[
                        :n_per_chan:dn])
                if ch==2:
                    data_2[str(buffersCompleted)] = []
                    data_2[str(buffersCompleted)].append(volt_array[
                            n_per_chan:-1:dn])
        
                # Add the buffer to the end of the list of available buffers.
                self.board.postAsyncBuffer(buffer.addr, buffer.size_bytes)
        finally:
            self.board.abortAsyncRead()
        # Compute the total transfer time, and display performance information.
        transferTime_sec = time.clock() - start
        print("Capture completed in %f sec" % transferTime_sec)
        buffersPerSec = 0
        bytesPerSec = 0
        recordsPerSec = 0
        if transferTime_sec > 0:
            buffersPerSec = buffersCompleted / transferTime_sec
            bytesPerSec = bytesTransferred / transferTime_sec
            recordsPerSec = recordsperbuffer*buffersCompleted/transferTime_sec
        print("Captured %d buffers (%f buffers per sec)" %
              (buffersCompleted, buffersPerSec))
        print("Captured %d records (%f records per sec)" %
              (recordsperbuffer * buffersCompleted, recordsPerSec))
        print("Transferred %d bytes (%f bytes per sec)" %
              (bytesTransferred, bytesPerSec))
        total_time=(samplesperrecord*recordsperbuffer*buffersperacq)/clock_rate
        sample_number = int(sample_rate*total_time)
        time_array = np.linspace(0, total_time, sample_number)
        
        return (data_1, data_2, time_array)
        
    def stitch_data(self, data_1, data_2, buffersperacq):
        stitched_data_1 = []
        stitched_data_2 = []
        for i in range(1,buffersperacq+1):
            stitched_data_1.append(data_1[str(i)][0])
            stitched_data_2.append(data_2[str(i)][0])
        stitched_data_1 = np.asarray(stitched_data_1)
        stitched_data_2 = np.asarray(stitched_data_2)
        stitched_data_1 = stitched_data_1.flatten()
        stitched_data_2 = stitched_data_2.flatten()
        return (stitched_data_1, stitched_data_2)
    
    def save_data(self, acquired_data, clock_rate, sample_rate, channels = 1,
                  split = True, input_range_V = 2):
        (buffer_array, samplesperrecord, recordsperbuffer, 
         buffersperacq) = acquired_data
        data = {}
        for c in range(channels):
            data[str(c+1)] = []
        for i in range(buffersperacq):
            value_array = buffer_array[i].buffer
            volt_array = self.code_to_volt(value_array, input_range_V)
            n_per_chan = len(volt_array)//channels
            dn = int(clock_rate/sample_rate)
            for c in range(channels):
                data[str(c+1)].append(volt_array[
                        c*n_per_chan:(c+1)*n_per_chan:dn])
        for c in range(channels):
            data[str(c+1)] = np.asarray(data[str(c+1)])
            if not split:
                data[str(c+1)] = data[str(c+1)].flatten()
                total_time = (samplesperrecord*recordsperbuffer*
                          buffersperacq)/clock_rate
            else:
                total_time = (samplesperrecord*recordsperbuffer)/clock_rate
        sample_number = int(sample_rate*total_time)
        time_array = np.linspace(0, total_time, sample_number)
#        time_array = np.array([k for k in time_array], dtype = 'uint16')
        if channels == 1:
            return (time_array, data['1'])
        else:
            return (time_array, data['1'], data['2'])
        
    def code_to_volt(self, value_array, input_range_V):
        coderange = 65535 #(2**16-1)
        array = np.array([2*
                input_range_V*(k/coderange-0.5) for k in value_array])
        return array
        
                
            