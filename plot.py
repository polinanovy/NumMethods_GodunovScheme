#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

fop = open('RESULT_A_01t')
x = []; rho = []; v = []; p = []
for line in fop.readlines():
    sl = line.split()
    x.append(float(sl[0])); rho.append(float(sl[1])); v.append(float(sl[2])); p.append(float(sl[3]))
    
fop = open('RESULT')
x_new = []; rho_new = []; v_new = []; p_new = []
for line in fop.readlines():
    sl = line.split()
    x_new.append(float(sl[0])); rho_new.append(float(sl[1])); v_new.append(float(sl[2])); p_new.append(float(sl[3]))    

fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(90, 30))
font = {'size'   : 40, 'family' : 'sans-serif'}
mpl.rc('font', **font)
plt.rc('text', usetex=True)

ax[0].set_xlabel('x')
ax[0].set_ylabel(r'$\rho$')
ax[0].axvline(0, color = 'k')
ax[0].plot(x, rho, 'k-')
ax[0].plot(x_new, rho_new, 'r-')
ax[0].grid()

ax[1].set_xlabel('x')
ax[1].set_ylabel('v')
ax[1].axvline(0, color = 'k')
ax[1].plot(x, v, 'k-')
ax[1].plot(x_new, v_new, 'r-')
ax[1].grid()

ax[2].set_xlabel('x')
ax[2].set_ylabel('p')
ax[2].axvline(0, color = 'k')
ax[2].plot(x, p, 'k-')
ax[2].plot(x_new, p_new, 'r-')
ax[2].grid()

ax[1].set_title('N=40, C=0.6, Configuration A')
plt.savefig('plot_A.png',bbox_inches='tight')
plt.clf()
plt.cla()
plt.close()