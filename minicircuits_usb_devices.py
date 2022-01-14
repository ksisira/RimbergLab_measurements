# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 15:45:00 2020

@author: Sisira
Module containing classes for miniciruits switches via USB control
"""

import clr # pythonnet
# Reference the DLL for switch
clr.AddReference('mcl_RF_Switch_Controller_NET45')
# Reference the DLL for attenuator
clr.AddReference('mcl_RUDAT_NET45')

from mcl_RF_Switch_Controller_NET45 import USB_RF_SwitchBox
from mcl_RUDAT_NET45 import USB_RUDAT

class usb_instrument_switch:

    def __init__(self, serial_number):
        self.sw = USB_RF_SwitchBox()  # Create an instance of the switch class
        self.sw.Connect(SN=serial_number)
        Responses = self.sw.Send_SCPI(":SN?", "")       # Read serial number
        print (str(Responses[2]))   # Python interprets the response as a 
        #tuple[function return (0 or 1), command parameter, response parameter]

        Responses = self.sw.Send_SCPI(":MN?", "")            # Read model name
        print (str(Responses[2]))
        
class usb_instrument_attenuator:

    def __init__(self, serial_number):
        self.att = USB_RUDAT()  # Create an instance of the attenuator class
        self.att.Connect(SN=serial_number)
        Responses = self.att.Send_SCPI(":SN?", "")       # Read serial number
        print (str(Responses[1]))   # Python interprets the response as a 
        #tuple[function return (0 or 1), command parameter, response parameter]

        Responses = self.att.Send_SCPI(":MN?", "")            # Read model name
        print (str(Responses[1]))
        
class minicircuits_rc1sp4ta18(usb_instrument_switch):
    def __init__(self, serial_number='02003310056'):
        # inherit all the methods of usb_instrument_switch class
        super().__init__(serial_number)
        
    def set_switch(self,switch_state):
        """
        Sets the switch state
        
        Parameters:
            switch_state (int): switch number 1,2,3,4
        """
        switch_addr = "SP4TA:STATE:{}".format(switch_state)
        self.sw.Send_SCPI(switch_addr,"")      
        # Set switch state (SW A, COM<>4)
        Responses = self.sw.Send_SCPI("SP4TA:STATE?", "")    
        # Read switch state (SP4TA)
        print (str(Responses[2]))
        
class minicircuits_rc4spdta18(usb_instrument_switch):
    def __init__(self, serial_number='02002180110'):
        # inherit all the methods of gpib_instrument class
        super().__init__(serial_number)
        
    def set_switch(self,switch_name,switch_state):
        """
        Sets the switch state of the corresponding switch name
        
        Parameters:
            switch_name (str): switch name A, B, C, D
            switch_state (int): switch number 0, 1
        """
        self.sw.Set_Switch(switch_name,switch_state) 
        # Set switch state (SW Number, COM<>0/1)
        Responses = self.sw.GetSwitchesStatus(0)
        # Read switch state in integer (SPDT)
        return bin(Responses[1])
        
    def set_switch_sync(self,state):
        """
        Sets all the four switches simultaneously.
        
        Parameters:
            state (str): State of all the four switches in binary with LSB=A
            and so on. For eg: 0101 represents D->1, C->2, B->1 and A->2.
        """
        self.sw.Set_SwitchesPort(int(state,2))
        Responses = self.sw.GetSwitchesStatus(0)
        return bin(Responses[1])
        
        
class minicircuits_rcdat800060(usb_instrument_attenuator):
    def __init__(self, item_label='1'):
        # inherit all the methods of usb_instrument_attenuator class
        if item_label == '1':
            serial_number = '11805240010'
        elif item_label == '2':
            serial_number = '11805240007'
        else:
            print('Item label not identified.')
        super().__init__(serial_number)
        self.set_startup_attenuation_mode('L')
        Responses = self.get_attenuation_dB()
        self.store_last_attenuation_value()
        print('Attenuation is {} dB'.format(Responses))
        
    def get_attenuation_dB(self):
        """
        Gets the attenuation.
        """
        Responses = self.att.Send_SCPI(':ATT?','')
        self.store_last_attenuation_value()
        return float(Responses[1])
    
    def set_attenuation_dB(self, attenuation):
        """
        Sets the attenuation.
        
        Parameters:
            attenuation (float): attenuation in dB.
        """
        self.att.Send_SCPI(':SETATT={}'.format(attenuation),'')
        Responses = self.get_attenuation_dB()
        self.store_last_attenuation_value()
        return Responses
    
    def get_startup_attenuation_mode(self):
        """
        Gets the startup attenuation mode.
        (L for last attenuation before device was last powered off,
            F for fixed user defined value,
            N for the factory default value)
        """
        Response_mode = self.att.Send_SCPI(':STARTUPATT:INDICATOR?','')
        if Response_mode[1] == 'F':
            Response_att = self.get_startup_attenuation_value_dB()
        else:
            Response_att = None
        self.store_last_attenuation_value()
        return (Response_mode[1],Response_att)
        
    def set_startup_attenuation_mode(self,mode,value=None):
        """
        Sets the startup attenuation mode.
        
        Parameters:
            mode (str): Chooses the attenuation value defined by the mode. 
            (L for last attenuation before device was last powered off,
            F for fixed user defined value,
            N for the factory default value)
            value (float): Set the startup attenuation value if in fixed mode.
        """
        self.att.Send_SCPI(':STARTUPATT:INDICATOR:{}'.format(mode),'')
        if mode == 'F' and value:
            self.set_startup_attenuation_value_dB(value)
        Responses = self.get_startup_attenuation_mode()
        self.store_last_attenuation_value()
        return Responses
    
    def store_last_attenuation_value(self):
        """
        Stores the last attenuation value. Necessary to run this if the
        attenuator is set to startup with the last known value.
        """
        self.att.Send_SCPI(':LASTATT:STORE:INITIATE','')
        
    def get_startup_attenuation_value_dB(self):
        """
        Gets the startup attenuation value. This command is only applicable
        if the startup mode is set to fixed.
        """
        Responses = self.att.Send_SCPI(':STARTUPATT:VALUE?','')
        self.store_last_attenuation_value()
        return Responses[1]
    
    def set_startup_attenuation_value_dB(self, attenuation):
        """
        Sets the startup attenuation value. This command is only applicable
        if the startup mode is set to fixed.
        
        Parameters:
            attenuation (float): The value of attenuation in dB to be set for
            fixed attenuation mode.
        """
        self.att.Send_SCPI(':STARTUPATT:VALUE:{}'.format(attenuation),'')
        Responses = self.get_startup_attenuation_value_dB()
        self.store_last_attenuation_value()
        return Responses
    
    
        
