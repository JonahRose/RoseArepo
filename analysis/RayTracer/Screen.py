# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 18:26:15 2020

@author: Jonah Rose
"""

class Screen():
    
    def __init__(self, pixels, objs, lights, cam):
        self.obj_pos = []
        self.cam_pos = []
        self.lig_pos = []
        
        self.objs = []
        self.lights = []
        
        self.pixels = []
        
        
        return