# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 18:41:18 2020

@author: Sisira
Module containing classes for instruments requiring usb connection.
"""
import visa
import numpy as np

class usb_instrument:

    def __init__(self, addr,model_str):
        addr_str = 'USB{}::{}::INSTR'.format(addr,model_str)
        self.instr = visa.ResourceManager().open_resource(addr_str)
    
    def read(self):
        # message = str
        return self.instr.read()
    
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
        
class keysight_33600A(usb_instrument):
    def __init__(self, addr,model_str='0x0957::0x5707::MY53804562'):
        #2391::22279::MY53804562::0
        # inherit all the methods of gpib_instrument class
        super().__init__(addr,model_str)
        
    def set_waveform(self,channel,wave_str):
        """
        Sets the waveform of source channel.
        
        Parameters:
            channel (int): Source channel
            wave_str (str): 'SIN', 'SQUare', 'TRIangle', 'PULS'
        """
        message = 'SOUR{}:APPL:{}'.format(channel,wave_str)
        self.write(message)
        
    def set_frequency(self,channel,frequency):
        """
        Sets the frequency of source channel.
        
        Parameters:
            channel (int): Source channel
            frequency (float): Value of frequency in Hz
        Returns:
            float: Value of frequency in Hz
        """
        message = 'SOUR{}:FREQ {}'.format(channel,frequency)
        self.write(message)
        freq_str = self.query('SOUR{}:FREQ?'.format(channel))
        return float(freq_str)
    
    def set_phase_deg(self,channel,phase):
        """
        Sets the phase of source channel in degrees.
        
        Parameters:
            channel (int): Source channel
            phase (float): Value of phase in deg.
        Returns:
            float: Value of phase in deg.
        """
        message = 'SOUR{}:PHAS {}'.format(channel,phase)
        self.write(message)
        phase_str = self.query('SOUR{}:phase?'.format(channel))
        return float(phase_str)
    
    def set_voltage(self,channel,voltage):
        """
        Sets the peak to peak voltage of source channel.
        
        Parameters:
            channel (int): Source channel
            voltage (float): Value of voltage in Vpp.
        Returns:
            float: Value of voltage in Vpp.
        """
        message = 'SOUR{}:VOLT {}'.format(channel,voltage)
        self.write(message)
        volt_str = self.query('SOUR{}:VOLT?'.format(channel))
        return float(volt_str)
    
    def set_voltage_offset(self,channel,volt_off):
        """
        Sets the offset voltage of source channel.
        
        Parameters:
            channel (int): Source channel
            volt_off (float): Value of offset voltage in V.
        Returns:
            float: Value of offset voltage in V.
        """
        message = 'SOUR{}:VOLT:OFFS {}'.format(channel,volt_off)
        self.write(message)
        volt_str = self.query('SOUR{}:VOLT:OFFS?'.format(channel))
        return float(volt_str)
    
    def toggle_output(self,channel,state):
        """
        Switches the output to on/off state.
        
        Parameters:
            channel (int): Source channel
            state (bool): State of the output 0/1
        Returns:
            int: Output state on/off
        """
        message = 'OUTP{} {}'.format(channel,state)
        self.write(message)
        out_str = self.query('OUTP{}?'.format(channel))
        return int(out_str)
    
    def toggle_FM(self,channel,state):
        """
        Switches the state of frequency modulation in source channel.
        
        Parameters:
            channel (int): Source channel
            state (bool): State of the frequency modulation.
        Returns:
            int: State of FM.
        """
        message = 'SOUR{}:FM:STAT {}'.format(channel,state)
        self.write(message)
        state_str = self.query('SOUR{}:FM:STAT?'.format(channel))
        return int(state_str)
    
    def set_FM_source(self,channel,source):
        """
        Sets the source of FM in source channel.
        
        Parameters:
            channel (int): Source channel
            source (str): Source of FM - INTernal, EXTernal, CH1, CH2
        Returns:
            str: Source of the FM.
        """
        message = 'SOUR{}:FM:SOUR {}'.format(channel,source)
        self.write(message)
        source_str = self.query('SOUR{}:FM:SOUR?'.format(channel))
        return source_str.split()[0]
    
    def set_FM_deviation(self,channel,freq_dev):
        """
        Sets the frequency deviation of source channel.
        
        Parameters:
            channel (int): Source channel
            freq_dev (float): Value of frequency deviation in Hz.
        Returns:
            float: Value of frequency deviation in Hz.
        """
        message = 'SOUR{}:FM:DEV {}'.format(channel,freq_dev)
        self.write(message)
        freq_str = self.query('SOUR{}:FM:DEV?'.format(channel))
        return float(freq_str)
    
    def set_freq_state(self,freq1,freq2):
        """
        Sets the frequency of both channels.
        
        Parameters:
            freq1: Frequency of first channel.
            freq2: Frequency of second channel.
            
        Returns:
            tuple: Frequency of channel 1 and channel 2 respectively.
        """
        ch_one = self.set_frequency(1,freq1)
        ch_two = self.set_frequency(2,freq2)
        return (ch_one,ch_two)
    
    def set_volt_state(self,volt1,volt2):
        """
        Sets the voltage of both channels.
        
        Parameters:
            volt1: Voltage of first channel in Vpp.
            volt2: Voltage of second channel in Vpp.
            
        Returns:
            tuple: Voltage of channel 1 and channel 2 respectively.
        """
        ch_one = self.set_voltage(1,volt1)
        ch_two = self.set_voltage(2,volt2)
        return (ch_one,ch_two)
    
    def set_offset_state(self,offset1,offset2):
        """
        Sets the voltage offset of both channels.
        
        Parameters:
            offset1: Voltage offset of first channel in V.
            offset2: Voltage offset of second channel in V.
            
        Returns:
            tuple: Voltage offset of channel 1 and channel 2 respectively.
        """
        ch_one = self.set_voltage_offset(1,offset1)
        ch_two = self.set_voltage_offset(2,offset2)
        return (ch_one,ch_two)
    
    def set_phase_state(self,phase1,phase2):
        """
        Sets the phase of both channels.
        
        Parameters:
            phase1: Phase offset of first channel in deg.
            phase2: Phase offset of second channel in deg.
            
        Returns:
            tuple: Phase offset of channel 1 and channel 2 respectively.
        """
        ch_one = self.set_phase_deg(1,phase1)
        ch_two = self.set_phase_deg(2,phase2)
        return (ch_one,ch_two)
    
    def toggle_tracking(self, channel=1, state=1):
        """
        Toggles the tracking of the two sources to create identical outputs.
        
        Parameters:
            channel (int): The other channel tracks this channel.
            state (str): State of the channel to on, off or inverted. Options
            'on' , 'off', or 'inv'
        """
        message = 'sour{}:trac {}'.format(channel,state)
        self.write(message)
        
    def set_input_load(self, channel=1, load='def'):
        """Sets the input load of the channel.
        
        Parameters:
            channel (int): Source channel
            load (str): load value. Generally used 'def' for 50 Ohms
            and 'inf' for high impedance.
        """
        message = 'outp{}:load {}'.format(channel, load)
        self.write(message)
    

    