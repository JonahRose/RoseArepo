# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 18:36:22 2020

@author: Jonah Rose
"""
import numpy as np

class Ray():
    
    def __init__(self, start_pos, end_pos, shadow=False):
        if type(start_pos) == type([]):
            start_pos = np.array(start_pos)
        if type(end_pos) == type([]):
            end_pos = np.array(end_pos)
            
        self.start_pos = start_pos
        self.end_pos = end_pos 
        
        self.p1 = end_pos - start_pos
        
        self.obj_idx = False
        
        self.hit = False
        self.shadow = shadow
        return
    
    def __mul__(self, scalar):
        self.p1 = self.p1 * scalar
        self.end_pos = self.start_pos + self.p1
        return self
    
    def __rmul__(self, scalar):
        self.p1 = self.p1 * scalar
        self.end_pos = self.start_pos + self.p1
        return self
    
    def __truediv__(self, scalar):
        p1 = self.p1 / scalar
        end_pos = self.start_pos + p1
        ray = Ray(self.start_pos, end_pos)
        return ray
    
    def __neg__(self):
        ray = Ray(self.end_pos, self.start_pos)
        return ray
    
    def __add__(self, ray):
        return Ray(self.start_pos, ray.start_pos + ray.p1)
    
    def __sub__(self, ray):
        return self + -ray
    
    def __repr__(self):
        return f"{self.start_pos} + t*{self.p1}"
    
    def cross(self, ray):
        return np.cross(self.p1, ray.p1)
    
    def mag(self):
        return np.sqrt(np.dot(self.p1, self.p1))
    
    def dot(self, ray):
        return np.sum(self.p1 * ray.p1)

        
        
    