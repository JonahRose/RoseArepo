# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 18:26:15 2020

@author: Jonah Rose
"""

import numpy as np

class Screen():
    
    def __init__(self, pixels, objs, lights, cam):
        self.obj_pos = []
        self.cam_pos = []
        self.lig_pos = []
        
        self.objs = objs
        self.lights = lights
        self.cam = cam
        
        self.pixels = np.zeros((pixels[0], pixels[1], 3), dtype=int)
        
        return
    
    def get_pos(self):
        for obj in self.objs:
            self.obj_pos.append(obj.pos)
        for lig in self.lights:
            self.lig_pos.append(lig.pos)
        self.cam_pos.append(self.cam.pos)
        return
    
    #if ray intersects object return pixel location
    def hit(self, ray):
        min_point = None
        min_dist = 1e10
        for i,obj in enumerate(self.objs):
            hit_point, hit_dist = obj.intersect(ray)
            if hit_point is not None:                                
                if hit_dist < min_dist:
                    min_dist = hit_dist
                    min_point = hit_point
                    ray.obj_idx = i
        return min_point
                
                
                
        
        