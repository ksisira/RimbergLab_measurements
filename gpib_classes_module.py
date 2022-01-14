# -*- coding: utf-8 -*-
"""
Created on Sun Feb  9 13:18:31 2020

@author: Sisira
Module containing classes for instruments requiring gpib connection.
"""
import numpy as np

from instrument_classes_module import gpib_instrument as gpib

class SRS_844(gpib):
    def __init__(self, addr):
        """
        inherit all the methods of gpib_instrument class from
        instrument_classes_module
        """
        super().__init__(addr)
        
    def get_output_display(self):
        """
        Gets the current display.
        
        Returns:
            str: current display
        """
        
        chone_str=self.query('DDEF?1')
        chone_out=int(chone_str)
        chtwo_str=self.query('DDEF?2')
        chtwo_out=int(chtwo_str)
        chone_disp=['X','Rvolts','Rdbm']
        chtwo_disp=['Y','theta']
        return chone_disp[chone_out]+chtwo_disp[chtwo_out]
    
    def set_output_display(self,display='XY'):
        """
        Sets the display into (X,Y), (Rvolts,theta) or (Rdbm,theta)
        
        Parameters:
            display (str): 'XY','Rvolts','Rdbm'
         
         Returns:
            str: Current display 
        """
        
        if display=='XY':
            self.write('DDEF1,0')
            self.write('DDEF2,0')
            
        elif display=='Rvolts':
            self.write('DDEF1,1')
            self.write('DDEF2,1')
            
        elif display=='Rdbm':
            self.write('DDEF1,2')
            self.write('DDEF2,1')
            
        return self.get_output_display()
    
    def get_reference_phase(self):
        """
        Gets the reference phase in degrees.
        
        Returns:
            float: Phase reference in deg
        """
        
        phase_str=self.query('PHAS?')
        return float(phase_str)
    
    def set_reference_phase_deg(self,phase):
        """
        Sets the reference phase in degrees.
        
        Parameters:
            phase (float): The reference phase in degrees.
        
        Returns:
            float: the reference phase.
        """
        
        self.write('PHAS{}'.format(phase))
        return self.get_reference_phase()
    
    def get_sensitivity(self):
        """
        Gets the sensitivity.
        
        Returns:
            str: Value of sensitivity.
        """
        
        sens_index=self.query('SENS?')
        sens_str=['100 nVrms / -127 dBm', '300 nVrms / -117 dBm',
                  '1 µVrms / -107 dBm','3 µVrms / -97 dBm',
                  '10 µVrms / -87 dBm','30 µVrms / -77 dBm',
                  '100 µVrms / -67 dBm', '300 µVrms / -57 dBm',
                  '1 mVrms / -47 dBm','3 mVrms / -37 dBm',
                  '10 mVrms / -27 dBm','30 mVrms / -17 dBm',
                  '100 mVrms / -7 dBm','300 mVrms / +3 dBm',
                  '1 Vrms / +13 dBm']
        return sens_str[int(sens_index)]
    
    def set_sensitivity_nV(self,sens):
        """
        Sets the sensitivity.
        
        Parameters:
            sens(float): Sensitivity in nV.
            
        Return:
            str: Value of sensitivity.
        """
        
        sens_value=np.array([100,300,1e3,3e3,10e3,30e3,100e3,300e3,1e6,3e6,
                             10e6,30e6,100e6,300e6,1e9])
        sens_index=np.where(sens_value==sens)
        self.write('SENS{}'.format(sens_index[0][0]))
        return self.get_sensitivity()
    
    def get_time_constant(self):
        """
        Gets the time constant.
        
        Returns:
            str: Value of time constant.
        """
        
        tc_index=self.query('OFLT?')
        tc_str=['100 µs', '300 µs', '1 ms', '3 ms', '10 ms', '30 ms', '100 ms',
                '300 ms', '1 s', '3 s', '10 s', '30 s', '100 s', '300 s', 
                '1 ks', '3 ks', '10 ks', '30 ks']
        return tc_str[int(tc_index)]
    
    def set_time_constant_sec(self,tc):
        """
        Sets the time constant.
        
        Parameters:
            tc(float): Time constant in seconds.
            
        Returns:
            str: Value of time constant.
        """
        
        tc_value=np.array([100e-6,300e-6,1e-3,3e-3,10e-3,30e-3,100e-3,
                           300e-3,1,3,10,30,100,300,1e3,3e3,10e3,30e3])
        tc_index=np.where(tc_value==tc)
        self.write('OFLT{}'.format(tc_index[0][0]))
        return self.get_time_constant()
    
    def push_sensitivity(self,act):
        """
        Push sensitivity one up or down.
        
        Parameters:
            act (str): 'up' or 'down'.
            
        Returns:
            Value of sensitivity.
        """
        if act=='up':
            self.write('KEYP 4')
        if act=='down':
            self.write('KEYP 5')
        return self.get_sensitivity()
    
    def push_time_constant(self,act):
        """
        Push time constant one up or down.
        
        Parameters:
            act (str): 'up' or 'down'.
            
        Returns:
            Value of time constant.
        """
        if act=='up':
            self.write('KEYP 0')
        if act=='down':
            self.write('KEYP 1')
        return self.get_time_constant()
    
    def set_expand(self, channel, expand = 1, display = 0):
        """
        Set expand for the output.
        
        Parameters:
            channel (int): Channel number
            display (int): Display for channel. 0 implies X/Y.
            expand (int): Expansion factor 1/10/100
            
        Returns:
            Value of expand.
        """
        value = int(np.log10(expand))
        self.write('DEXP{},{},{}'.format(channel, display, value))
        
    def set_filter_slope(self, slope):
        """
        Sets the filter slope.
        
        Parameters:
            slope (int): Filter slope value - 6/12/18/24.
            
        Returns:
            str: Value of filter slope.
        """
        slope_index=slope//6
        self.write('OFSL{}'.format(slope_index))
        return self.get_filter_slope()
    
    def get_filter_slope(self):
        """
        Gets the filter slope.
        
        Returns:
            str: Value of slope.
        """
        
        slope_int=self.query('OFSL?')
        slope_str=['0', '6 dB/oct', '12 dB/oct', '18 dB/oct', '24 dB/oct']
        return slope_str[int(slope_int)]
    
    def get_offset(self, quantity):
        """
        Gets the offset of channel.
        
        Parameters:
            quantity (str): Quadrature offset to measure - 'X'/'Y' 
        Returns:
            offset_percent (float): Offset in percentage of full scale.
        """
        if quantity == 'X':
            add_str = '1,0'
        elif quantity == 'Y':
            add_str = '2,0'
        offset_percent=self.query('DOFF?{}'.format(add_str))
        return eval(offset_percent)
    
    def set_offset(self, quantity, percent):
        """
        Sets the offset in percentage of full scale.
        
        Parameters:
            quantity (str): Quadrature offset to measure - 'X'/'Y'
            percent (float): Offset in percentage of full scale (from -110% to
                                (110%)).
        Returns:
            offset_percent (float): Offset in percentage of full scale.
        """
        if quantity == 'X':
            add_str = '1,0'
        elif quantity == 'Y':
            add_str = '2,0'
        self.write('DOFF{},{}'.format(add_str,percent))
        offset_percent = self.get_offset(quantity)
        print(offset_percent)
        return offset_percent
        
    
class keysight_n5183b(gpib):
    
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
        
    def sweep_trigger_single(self):
        """
        This command aborts the current sweep, then either arms or arms and 
        starts a single list, depending on the trigger type.
        """
        message = ':TSW'
        self.write(message)
        
    def set_sweep_center_freq(self,freq):
        """
        This command sets the center frequency for a step sweep.
        
        Parameters:
            freq (float): Center frequency for a step sweep in Hz.
        """
        message = ':FREQ:CENT {}'.format(freq)
        self.write(message)
    
    def set_sweep_span(self,span):
        """
        This command sets the length of the frequency range for a step sweep.
        Span setting is symmetrically divided by the selected center frequency.
        
        Parameters:
            span (float): Value of span in Hz.
        """
        
        message = ':FREQ:SPAN {}'.format(span)
        self.write(message)
        
    def set_sweep_dwell_time(self,time):
        """
        This command enables you to set the dwell time for a step sweep.
        The variable <value> is expressed in units of seconds 
        with a 0.001 resolution.
        
        Parameters:
            time (float): Dwell time in seconds.
        """
        
        message = ':SWE:DWEL {}'.format(time)
        self.write(message)
    
    def set_sweep_points(self,points):
        """
        This command defines the number of step sweep points.
        
        Parameters:
            points (int): Number of points in the sweep.
        """
        message = ':SWE:POIN {}'.format(int(points))
        self.write(message)
        
    def set_mode(self,mode = 'CW'):
        """
        This command sets the frequency mode of the signal generator to 
        CW or swept.
        
        Parameters:
            mode (str): Mode for enabling continuous ('CW') or list/step
            sweep ('LIST').
        """
        message = ':FREQ:MODE {}'.format(mode)
        self.write(message)
        
        
class tektronix_tds7704b(gpib):
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
    
    def get_horizontal_scale(self):
        # NOT INDEPENDENT: AFFECTED BY SAMPLE RATE AND RECORD LENGTH
        # gets horizontal scale in seconds per division
        message = 'hor:sca?'
        return self.query(message)
    

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
        
    def acquire(self, state):
        # Starts/stops an acquisition in the selected acquisition mode 
        #(cont/single)
        if state:
            mode_str = 'ON'
        else:
            mode_str = 'OFF'
        message = 'ACQ:STATE {}'.format(mode_str)
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
            
