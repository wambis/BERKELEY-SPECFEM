#!/usr/bin/env python
import numpy as np


F=open("STATIONS.wamba","r") 
ff=open("STATIONS","w")
for l in F:
    net, stn, lat, lon, depth, B =l.split() 
    ll= "%s %7s %9s %9s %9s %9s"%(net,stn,lat,lon,"0.00000","0.00000") 
    ff.write("%s\n"%(ll)) 


