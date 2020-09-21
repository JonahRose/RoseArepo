# -*- coding: utf-8 -*-
"""
Created on Sun Sep 20 14:47:04 2020

@author: astro
"""

from RoseArepo.analysis.RayTracer.Ray import Ray
import numpy as np

class Plane():
    
    def __init__(self, low_left, dx, dy, color):
        self.low_left = np.array(low_left)
        self.dx = dx #index 0
        self.dy = dy #index 1?
        self.area = dx*dy
        self.color = color
        
        self.norm = Ray(low_left, self.low_left + Ray(low_left, [0,0,dx]).cross(Ray(low_left, [0,dy,0])))
        return
    
    def in_bounds(self, P):
        if P[0] < self.low_left[0] or P[0] > self.low_left[0]+self.dx:
            return False
        if P[1] < self.low_left[1] or P[1] > self.low_left[1]+self.dy:
            return False        
        return True
        
    
    def intersect(self, ray):
        if ray.mag() == 0:
            return None, None
        if self.area == 0:
            return None, None
        ldn = ray.dot(self.norm)
        if np.abs(ldn) < 1e-6: #ray parallel to plane
            return None, None
        num = -1 * (self.norm.p1[0]*ray.start_pos[0] +\
                    self.norm.p1[1]*ray.start_pos[1] +\
                    self.norm.p1[2]*ray.start_pos[2])
        denom = self.norm.p1[0]*ray.p1[0] +\
                self.norm.p1[1]*ray.p1[1] +\
                self.norm.p1[2]*ray.p1[2]
        t = num / denom
        P = ray.start_pos + t*ray.p1
        if not self.in_bounds(P):
            return None, None
        dist = Ray(ray.start_pos, P).mag()
        return P, dist

"""
p = Plane([100,0,0], 100,100, [255,255,255])
print(p.norm)
ray = Ray([25,25,25], [25,25,100])
h, d = p.intersect(ray)
print(h,d)
"""