import numpy as np
import profileclass
import math

class PinnaclePlanarDose:

    def __init__(self, filename):
        self._array_txt = np.genfromtxt(filename, delimiter = ',')
        self.planar_values = np.array(self._array_txt[1:len(self._array_txt),1:(len(self._array_txt[0])-1)])
        self.ycoord = np.array(self._array_txt[1:len(self._array_txt),0])
        self.xcoord = np.array(self._array_txt[0,1:len(self._array_txt[0])-1])
        
    
    def GetYcentralprofile(self):
        prof_vals = self.planar_values[0:len(self.ycoord),math.trunc(len(self.xcoord)/2)]
        _mat = np.array([self.ycoord,prof_vals])
        return profileclass.Profile(_mat.transpose())
    
    def GetXcentralprofile(self):
        prof_vals = self.planar_values[math.trunc(len(self.ycoord)/2),0:len(self.xcoord)]
        _mat = np.array([self.xcoord,prof_vals])
        return profileclass.Profile(_mat.transpose())
    
    
    
    
