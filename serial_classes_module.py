# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 14:23:03 2020

@author: Sisira
Module containing classes for instruments requiring serial connection.
"""
import visa
import numpy as np

class serial_instrument:

    def __init__(self, addr):
        addr_str = 'ASRL'+str(addr)+'::INSTR'
        self.instr = visa.ResourceManager().open_resource(addr_str)
        
    def baud_rate(self,bd_rate=None):
        if bd_rate:
            self.instr.baud_rate = bd_rate
        return self.instr.baud_rate
    
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
        
class Novatech_409B(serial_instrument):
    #Connects via RS-232
    #The saved state is echo disabled, if cleared using CLR commamnd, 
    #make sure echo is disabled before running this class.
    #initializes with 
        #baudrate to 19200
        #echo disabled (functions are written accordingly)
    
    def __init__(self, addr):
        # inherit all the methods of serial_instrument class
        super().__init__(addr)
        self.baud_rate(19200)
        self.write('e d')
        self.read()
        self.auto_io()
        self.print_state()
        
#    def set_ref(self, ref):#Not needed if /R option is installed.
#        if ref == 'e':
#            self.write('C e')
#            self.read()
#            self.write(' 1 aa')
#            self.read()
#        elif ref == 'i':
#            self.write('C i')
#            self.read()
        
    def get_state(self):
        self.write('QUE')
        freq=[]
        phase=[]
        amp=[]
        for i in range(4):
            state_str = self.read()
            str_array=state_str.split( )
            freq.append(int(str_array[0],16))
            phase.append(int(str_array[1],16))
            amp.append(int(str_array[2],16))
        self.read()
        freq=np.asarray(freq)
        phase=np.array(phase)
        amp=np.array(amp)
        return freq*1e-7, phase*np.pi/8192, amp/1023
    
    def print_state(self, channel=None):
        freq,phase,amp=self.get_state()
        if channel==None:
            print('Frequency of channels are '+str(freq)+' MHz')
            print('Phase of channel are '+str(phase)+' radians')
            print('Amplitude of channel are '+str(amp)+' Vpp')
        else:
            print('Frequency of channel '+str(channel)+ ' is '
                  +str(freq[channel])+' MHz')
            print('Phase of channel '+str(channel)+ ' is '
                  +str(phase[channel])+' radians')
            print('Amplitude of channel '+str(channel)+ ' is '
                  +str(amp[channel])+' V')
            
    def set_frequency_MHz(self,channel,freq):
        freq_str='F{channel} {freq}'.format(channel=channel,freq=freq)
        self.write(freq_str)
        self.read()
        self.print_state(channel)
        
    def set_voltage_Vpp(self,channel,amp):
        amp_str='V{channel} {amp}'.format(channel=channel,amp=int(amp*1023))
        self.write(amp_str)
        self.read()
        self.print_state(channel)
        
    def set_power_dbm(self,channel,power): #assuming 50 ohms
        #power in dbm = 10log(P/1mW,10) where P=Vpp**2/8/z0 where z0=50 ohms
        amp=np.sqrt(10**(power/10)*1e-3*50*8)
        amp_str='V{channel} {amp}'.format(channel=channel,amp=int(amp*1023))
        self.write(amp_str)
        self.read()
        self.print_state(channel)
        
    def set_phase_rad(self,channel,phase):
        phase_str='P{channel} {phase}'.format(channel=channel,
                    phase=int(phase*8192/np.pi))
        self.write(phase_str)
        self.read()
        self.print_state(channel)
        
    def manual_io_wait(self):
        self.write('I m')
        self.read()
        
    def manual_io_run(self):
         self.write('I p')
         self.read()
         
    def auto_io(self):
        self.write('I a')
        self.read()
        
    def phase_sync(self,state = True):
        if state:
            self.write('M a')
            self.read()
        else:
            self.write('M n')
            self.read()
        
    def save(self):
        self.write('S')
        self.read()
        
    def reset(self):
        self.write('R')
        self.read()
        
    def set_state_freq(self,channel0,channel1,freq0,freq1):
        self.manual_io_wait()
        self.set_frequency_MHz(channel0,freq0)
        self.set_frequency_MHz(channel1,freq1)
        self.phase_sync()
        self.manual_io_run()
        
    def set_state_power(self,channel0,channel1,power0,power1):
        self.manual_io_wait()
        self.set_power_dbm(channel0,power0)
        self.set_power_dbm(channel1,power1)
        self.phase_sync()
        self.manual_io_run()
        
    def set_state_Vpp(self,channel0,channel1,volt0,volt1):
        self.manual_io_wait()
        self.set_voltage_Vpp(channel0,volt0)
        self.set_voltage_Vpp(channel1,volt1)
        self.phase_sync()
        self.manual_io_run()
        
    def set_state_phase(self,channel0,channel1,phase0,phase1):
        self.manual_io_wait()
        self.set_phase_rad(channel0,phase0)
        self.set_phase_rad(channel1,phase1)
        self.manual_io_run()
        
        
        