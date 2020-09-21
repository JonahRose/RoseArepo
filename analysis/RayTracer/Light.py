# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 18:23:48 2020

@author: Jonah Rose
"""

import numpy as np

class Light():
    
    def __init__(self, pos, brightness):
        self.brightness = np.array(brightness)
        self.pos = np.array(pos)
        self.color = np.array([1,1,1])
        
        
        return