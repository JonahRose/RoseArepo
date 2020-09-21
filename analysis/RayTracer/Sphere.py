# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 18:25:02 2020

@author: Jonah Rose
"""

from RoseArepo.analysis.RayTracer.Ray import Ray
from RoseArepo.analysis.RayTracer.Light import Light
import numpy as np

class Sphere():
    
    def __init__(self, pos, radius, color):
        self.pos = pos
        self.r = radius
        self.color = np.array(color)
        
        return
    
    #TODO: Doesn't account for hitting back of sphere? maybe?
    def intersect(self, ray):
        if ray.mag() == 0:
            return None, None
        if self.r == 0:
            return None, None
        L = Ray(ray.start_pos, self.pos)
        if L.dot(ray) < 0 and ray.shadow: #bouncing off surface?
            return None, None
        if L.mag() < self.r-1: #-1 for rounding
            raise ValueError("Ray starts inside sphere?")
        t_ca =  ray * (L.dot(ray) / (ray.mag() *ray.mag()))
        d = -L + t_ca
        if d.mag() > self.r:
            return None, None
        t_hc_length = np.sqrt(self.r*self.r - d.mag()*d.mag())
        t_hc = -t_ca / t_ca.mag() * t_hc_length
        P = t_hc.end_pos
        dist = Ray(ray.start_pos, P).mag()
        return P, dist

"""
s = Sphere([25,25,25], 20, color=[255,0,0])
r = Ray(np.array([0,25,25]), np.array([100, 25, 25]))
h,d = s.intersect(r)
l = Light([50,25,25], [1,1,1])
sr = Ray(h, l.pos, shadow=True)
print(sr)
sh,sd = s.intersect(sr)
print(sh)
"""