#region modules

import numpy as np
import matplotlib.pyplot as plt
import math
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn
from scipy import stats
from matplotlib.offsetbox import AnchoredText
import matplotlib.dates as mdates
pd.plotting.register_matplotlib_converters(explicit=True)
# seaborn.set()
np.seterr(all='ignore')

#endregion


class snow(object):
    _dir = r'D:\DRIVE\TUBITAK\Snow_Module'
    _data = "Snow.csv"
    def __init__(self):
        self._working_directory = None
        self.Data_file = None
        self.df = None
        self.Tmin = None
        self.Tmean = None
        self.Tmax = None
        self.J = []
        self.t = None
        self.Lat = 46.0
        self.Date = None
        self.PEcalulated = None
        self.stataionel = 500
        self.zoneeleation = 1300
        self.tlapsrate = -0.0065
        self.plapsrate = 0.0004
        self.degreeday = 3.74
        self.snowpackinertia = 0.25
        self.snowstorage = None
        self.Pupdated = None
        self.actualmelt = None
        self.LiqP = None
        self.Date = None
        self.width = 15  # Figure width
        self.height = 10  # Figure height
    @property

    def process_path(self):
        return self._working_directory

    @process_path.setter
    def process_path(self, value):
        self._working_directory = value
        pass

    def DataRead(self):
        self.df = pd.read_csv(self.Data_file, sep=',', parse_dates=[0], header=0)
        self.df = self.df.set_index('Date')
        # get date index df.index.to_pydatetime()

    def InitData(self):
        self.Timn = self.df.Tmin
        self.Tmean = self.df.Tmean
        self.Tmax = self.df.Tmax
        self.P = self.df.P
        self.df['dayofyear'] = pd.DatetimeIndex(self.df.index.to_pydatetime()).day
        # self.J = np.array(self.df.dayofyear[:])
        self.PE = self.df.PE
        self.n = self.df.__len__()
        self.Pupdated = np.zeros(self.n)
        self.PEcalulated = np.zeros(self.n)
        self.Date = self.df.index.to_pydatetime()
        # pd.dt.
    # def Initzone(self):

    def Juliandate(self):
        Date = self.df.index.to_pydatetime()
        for i, date in enumerate(Date):
            self.J.append((Date[i].timetuple()).tm_yday)
        self.J = np.array(self.J)

    def Pecalc(self):
        self.Juliandate()
        for i in range(self.n):
            teta = 0.4093*math.sin(self.J[i]/58.1-1.405)
            cosGz = max(0.001,math.cos(self.Lat/57.3-teta))
            Gz = math.acos(cosGz)
            cosOM = max(-1,min(1-cosGz/math.cos(self.Lat/57.3)/math.cos(teta),1))
            OM = math.acos(cosOM)
            Eta = 1+math.cos(self.J[i]/58.1)/30
            cosPz = cosGz+math.cos(self.Lat/57.3)*math.cos(teta)*(math.sin(OM)/OM-1)
            Rad = 446*OM*Eta*cosPz
            self.PEcalulated[i] = max(0,Rad*(self.Tmean[i]+5)/28.5/100)

    def Zone(self):
        self.snowstorage = np.zeros(self.n)
        Tminz = self.Timn + (self.zoneeleation - self.stataionel)*self.tlapsrate
        Tmeanz = self.Tmean + (self.zoneeleation - self.stataionel)*self.tlapsrate
        Tmaxz = self.Tmax + (self.zoneeleation - self.stataionel)*self.tlapsrate
        Pretotalz = self.P*math.exp((self.zoneeleation-self.stataionel)*self.plapsrate)
        # persnow = (D3<=0,1,IF(B3>=0,0,1-D3/(D3-B3)))
        persnow = np.where(Tmaxz < 0.1,1,np.where(Tminz >= 0,0, 1 -Tmaxz/(Tmaxz-Tminz) ))
        self.LiqP =Pretotalz*(1-persnow)
        SolP = Pretotalz*persnow
        prepp = np.average(SolP)*0.9*365.25
        snowbeforemelt = np.zeros(len(SolP))
        snowbeforemelt[0] = 0
        snowpacktemp =np.zeros(len(SolP))
        snowpacktemp[0] = 0
        potmelt = np.zeros(len(SolP))
        potmelt[0] = np.where(snowpacktemp[0] == 0 , min(snowbeforemelt[0],max(0,self.degreeday*Tmeanz[0])),0)
        moderatingfactor = np.zeros(len(SolP))
        moderatingfactor[0] = np.where(snowbeforemelt[0] < prepp, snowbeforemelt[0]/prepp,1)
        self.actualmelt = np.zeros(len(SolP))
        self.actualmelt[0] = (0.9*moderatingfactor[0]+0.1)*potmelt[0]
        self.snowstorage[0] = - self.actualmelt[0] + snowbeforemelt[0]

        for i in range(1,self.n):
            snowbeforemelt[i] = SolP[i] + self.snowstorage[i-1]
            snowpacktemp[i] = min(0,self.snowpackinertia*snowpacktemp[i-1]+(1-self.snowpackinertia)*Tmeanz[i])
            potmelt[i] = np.where(snowpacktemp[i] == 0, min(snowbeforemelt[i], max(0, self.degreeday * Tmeanz[i])), 0)
            moderatingfactor[i] = np.where(snowbeforemelt[i] < prepp, snowbeforemelt[i] / prepp, 1)
            self.actualmelt[i] = (0.9 * moderatingfactor[i] + 0.1) * potmelt[i]
            self.snowstorage[i] = - self.actualmelt[i] + snowbeforemelt[i]

    def updateP(self):
        self.Pupdated = self.actualmelt + self.LiqP

    def updatedf(self):
        self.df['Julian_Date'] = self.J
        self.df['Snow_Storage'] = self.snowstorage
        self.df['Pupdated'] = self.Pupdated
        self.df['PECalculated'] = self.PEcalulated

    def draw(self):
        # plt.bar(a.df['P'], 'b--',a.df['Pupdated'], 'r-')
        # plt.bar(a.Date, 'b--',a.df['Pupdated'], 'r-')

        fig, ax1 = plt.subplots(figsize=(self.width, self.height))
        ax1.set_title('CemaNeige Snow Model Result', fontsize=16, fontweight='bold', style='italic')
        color = 'tab:red'
        ax1.set_xlabel('Date', style='italic', fontweight='bold', fontsize=14)
        ax1.set_ylabel('Precipitation , mm', color=color, style='italic', fontweight='bold', fontsize=14)
        ax1.set_ylim(0, max(self.df.Pupdated) * 1.1)
        ax1.tick_params(axis='x', labelrotation=45)
        ax1.bar(self.Date, self.P)
        ax1.bar(self.Date, self.Pupdated)
        # plt.bar(a.Date,a.P)
        # plt.bar(a.Date,a.Pupdated)
        plt.gca().legend(('Precipitation', 'Precipitation + Snow Melting'))
        fig.tight_layout()
        plt.show()

    def drawstorage(self):
        fig, ax1 = plt.subplots(figsize=(self.width, self.height))
        color = 'tab:red'
        ax1.set_title('CemaNeige Snow Model Result', fontsize=16, fontweight='bold', style='italic')
        ax1.set_xlabel('Date', style='italic', fontweight='bold', fontsize=14)
        ax1.set_ylabel('Snow Storage ', color=color, style='italic', fontweight='bold', fontsize=14)
        ax1.tick_params(axis='x', labelrotation=45)
        plt.plot(self.df['Snow_Storage'])
        fig.tight_layout()
        plt.show()

    def drawpe(self):
        fig, ax1 = plt.subplots(figsize=(self.width, self.height))
        ax1.set_title('CemaNeige Snow Model Result', fontsize=16, fontweight='bold', style='italic')
        ax1.set_xlabel('Date', style='italic', fontweight='bold', fontsize=14)
        color = 'tab:red'
        ax1.set_ylabel('Potentional Evaporation, mm', color=color, style='italic', fontweight='bold', fontsize=14)
        ax1.tick_params(axis='x', labelrotation=45)
        plt.plot(self.df['PECalculated'])
        fig.tight_layout()
        plt.show()

# Initilize object
snow = snow()
# Process path
snow.process_path = r'D:\DRIVE\TUBITAK\Snow_Module'
# Data file
snow.Data_file = os.path.join(snow.process_path, "Snow.csv")
# Calculate POT

snow.DataRead()
snow.InitData()
snow.Pecalc()
snow.Zone()
snow.updateP()
snow.updatedf()
snow.draw()
snow.drawpe()
snow.drawstorage()