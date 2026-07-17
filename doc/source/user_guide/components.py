#!/usr/bin/env python
"""components.py
Schematic from Sutterley and Velicogna, Remote Sensing (2019)
"""

import numpy as np
import scipy.interpolate
import matplotlib
import matplotlib.colors
import matplotlib.patches
import matplotlib.pyplot as plt

# adjust fonts
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Lato']
matplotlib.rcParams['font.monospace'] = ['DejaVu Sans Mono']

# create an initial rough estimate of the crustal surface
domain_length = 350000
x = np.linspace(0, domain_length, 1300).astype(np.float64)
crust = np.zeros_like(x)
crust[:140] = 250.0
crust[140:190] = np.linspace(250, 225, 50)
crust[190:210] = 225.0
crust[210:260] = np.linspace(225, 250, 50)
crust[260:500] = 250.0
crust[500:750] = 250.0 - 2e-8 * (x[500:750] - x[500]) ** 2
crust[750:850] = np.linspace(crust[749], 10, 100)
crust[850:1025] = 10 + 2e-8 * (x[850:1025] - x[850]) ** 2
crust[1025:] = crust[1024]

# create interpolated line with smooth transitions
SPL = scipy.interpolate.UnivariateSpline(x, crust, k=2)
XI = np.arange(0, domain_length, 100)
NI = len(XI)
C1 = SPL(XI)

# create sea surface
S1 = 160.0 * np.ones((NI)) + 2.0 * np.sin(np.arange(NI) * np.pi / 180.0 * 5)
(ii,) = np.nonzero(C1 >= np.max(S1) - 5)
S1[ii] = np.nan

# create ice surface
ice_surface = -0.0025 * (np.arange(NI) / 1.5 - NI / 3.625) ** 2 + 500.0
ice_surface[:974] = 250.0
ice_surface[1900:2100] = np.linspace(ice_surface[1899], 210, 200)
ice_surface[2100:] = 210.0
# create smooth spline of ice surface
SPL1 = scipy.interpolate.UnivariateSpline(XI, ice_surface, k=2)
I1 = SPL1(XI)
I1[:965] = np.nan
I1[2010:] = np.nan

# initialize output figure
f1, ax1 = plt.subplots(num=1, figsize=(12, 4), facecolor='#fcfcfc')

# add snowflakes
fontdict = dict(family='monospace', color='gray', weight='normal', size=22)
ax1.text(
    135000,
    530,
    '\u2744',
    fontdict=fontdict,
    ha='center',
    va='center',
    alpha=0.65,
)
ax1.text(
    145000,
    540,
    '\u2744',
    fontdict=fontdict,
    ha='center',
    va='center',
    alpha=0.65,
)
ax1.text(
    155000,
    530,
    '\u2744',
    fontdict=fontdict,
    ha='center',
    va='center',
    alpha=0.65,
)
# add trees
for x1 in (5000, 10000, 15000, 20000):
    color = matplotlib.colors.colorConverter.to_rgba('forestgreen', alpha=0.15)
    ax1.fill(
        [x1 - 2000, x1 + 2000, x1, x1 - 2000],
        [258, 258, 298, 258],
        facecolor=color,
        edgecolor='forestgreen',
    )
    ax1.plot([x1, x1], [250, 256], lw=2.0, color='forestgreen')
    ax1.plot([x1, x1], [258, 298], lw=0.5, color='forestgreen')
    for y in range(0, 30, 5):
        ax1.plot(
            [
                x1 - (2000.0 / 40.0) * (30 - y),
                x1,
                x1 + (2000.0 / 40.0) * (30 - y),
            ],
            [263.0 + y + 2.5, 263.0 + y, 263.0 + y + 2.5],
            lw=0.5,
            color='forestgreen',
        )
# add icebergs
facecolor = matplotlib.colors.colorConverter.to_rgba('gray', alpha=0.15)
edgecolor = matplotlib.colors.colorConverter.to_rgba('gray', alpha=0.65)
for icex in (210000, 220000, 230000):
    ax1.fill(
        [icex, icex + 5000, icex + 5000, icex, icex],
        [140, 140, 175, 175, 140],
        facecolor=facecolor,
        edgecolor=edgecolor,
        lw=1.5,
    )
# add glacier
ax1.fill_between(XI, C1, y2=I1, facecolor=facecolor, edgecolor=edgecolor, lw=2)

# add lake surface
ax1.plot(XI[360:715], 250.0 * np.ones((355)), color='#1da2d8', lw=2)
# add crust
ax1.plot(XI, C1, color='black', lw=2)
# add sea surface
ax1.plot(XI, S1, color='#1da2d8', lw=2)
# fill between surfaces
ax1.fill_between(
    XI[965:1050], C1[965:1050], y2=I1[965:1050], alpha=0.35, color='steelblue'
)
ax1.fill_between(XI, C1, y2=S1, alpha=0.35, color='#1da2d8')
ax1.fill_between(XI, C1, alpha=0.25, color='#835C3B')
ax1.fill_between(
    XI[300:710],
    C1[300:710],
    y2=250.0 * np.ones((410)),
    alpha=0.35,
    color='#1da2d8',
)
# add groundwater
color = matplotlib.colors.colorConverter.to_rgba('#1da2d8', alpha=0.35)
gw = matplotlib.patches.Ellipse(
    (22000, 190), 27500, 22.5, linewidth=2, facecolor=color, edgecolor='#1da2d8'
)
ax1.add_patch(gw)

# add labels
ax1.text(
    275000, 525, 'Atmospheric\nPressure', ha='center', va='center', fontsize=16
)
ax1.annotate(
    '',
    xy=(252500, 495),
    xytext=(252500, 555),
    arrowprops=dict(arrowstyle='<->'),
)
ax1.annotate(
    '',
    xy=(297500, 495),
    xytext=(297500, 555),
    arrowprops=dict(arrowstyle='<->'),
)
ax1.text(
    45000,
    300,
    'Terrestrial\nWater Storage',
    ha='center',
    va='center',
    fontsize=16,
)
ax1.text(
    22000, 140, 'Groundwater\nStorage', ha='center', va='center', fontsize=16
)
# t = ax1.text(33000,160,'Groundwater Storage',ha='center',va='center',fontsize=16,
#     bbox=dict(boxstyle="round",ec='#1da2d8',fc=color,lw=2))
# bb = t.get_bbox_patch()
# bb.set_boxstyle("round", rounding_size=1.0, pad=0.6)
ax1.text(33000, 40, 'Earthquakes', ha='center', va='center', fontsize=16)
ax1.annotate(
    '', xy=(3000, 40), xytext=(14000, 40), arrowprops=dict(arrowstyle='<-')
)
ax1.annotate(
    '', xy=(52000, 40), xytext=(63000, 40), arrowprops=dict(arrowstyle='->')
)
ax1.text(
    126500, 205, 'Elastic\nResponse', ha='center', va='center', fontsize=16
)
ax1.annotate(
    '',
    xy=(107500, 175),
    xytext=(107500, 235),
    arrowprops=dict(arrowstyle='<->'),
)
ax1.annotate(
    '',
    xy=(145500, 175),
    xytext=(145500, 235),
    arrowprops=dict(arrowstyle='<->'),
)
ax1.text(
    126500,
    85,
    'Glacial\nIsostatic\nAdjustment',
    ha='center',
    va='center',
    fontsize=16,
)
ax1.annotate(
    '',
    xy=(94500, 30),
    xytext=(107500, 80),
    arrowprops=dict(arrowstyle='<-', connectionstyle='arc3,rad=-0.3'),
)
ax1.annotate(
    '',
    xy=(158500, 30),
    xytext=(145500, 80),
    arrowprops=dict(arrowstyle='<-', connectionstyle='arc3,rad=0.3'),
)
# ax1.text(90000,530,'Accumulation',ha='left',va='center',fontsize=16)
# ax1.text(103000,350,'Surface\nAblation',ha='right',va='center',fontsize=16)
# ax1.text(220000,220,'Dynamic\nAblation',ha='center',va='center',fontsize=16)
# ax1.text(145000,370,'Glacier',color='gray',ha='center',va='center',fontsize=16)
ax1.text(
    145000, 370, 'Glacier\nMass Balance', ha='center', va='center', fontsize=16
)
ax1.text(346000, 25, 'Solid Earth Tides', ha='right', va='center', fontsize=16)
ax1.text(
    346000,
    103,
    'Tidal and Non-Tidal\nOceanic Processes',
    ha='right',
    va='center',
    fontsize=16,
)
ax1.text(
    257000, 103, 'Sea Level\nChange', ha='center', va='center', fontsize=16
)

# set x and y limits
ax1.set_xlim(XI[0], XI[-1])
ax1.set_ylim(0, 570)
# no axis
ax1.axis('off')

# adjust subplot within figure
f1.subplots_adjust(left=0.01, right=0.99, bottom=0.01, top=0.99, wspace=0.05)
plt.show()
