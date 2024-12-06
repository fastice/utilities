# -*- coding: utf-8 -*-
"""
Created on Mon Oct  9 12:42:03 2017

@author: ian
"""
from utilities.myerror import myerror
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from datetime import datetime
'''
    class defintion for an object based on a shape file point, polygon, or
    polyline
'''


class shpplot:
    """ \ngeoimage - object for scalar or velocity PS data + geodat data """

    def __init__(self, xdom=None, ydom=None, xplot=None, y1plot=None,
                 y2plot=None, y3plot=None,
                 e1plot=None, e2plot=None, e3plot=None, dates=None,
                 plotType=None, verbose=True, name=None, labels=None):
        # points where data are from
        self.plotType = None
        self.xdom, self.ydom = [], []
        self.xplot, self.y1plot, self.y2plot, self.y3plot = [], [], [], []
        self.e1plot, self.e2plot, self.e3plot = [], [], []
        self.dates, self.labels = [], []
        self.name = ''
        # domain types - need to expand list
        validTypes = ['point', 'linestring']
        if plotType is not None:
            if plotType.lower() in validTypes:
                self.plotType = plotType.lower()
        # area (domain) to plot, could be a point, polygon, or line
        if xdom is not None:
            self.xdom = xdom
        if ydom is not None:
            self.ydom = ydom
        # x, y data to plot
        if xplot is not None:
            self.xplot = xplot
        if y1plot is not None:
            self.y1plot = y1plot
        if y2plot is not None:
            self.y2plot = y2plot
        if y3plot is not None:
            self.y3plot = y3plot
        if e1plot is not None:
            self.e1plot = e1plot
        if e2plot is not None:
            self.e2plot = e2plot
        if e3plot is not None:
            self.e3plot = e3plot
        if dates is not None:
            self.dates = dates
        if labels is not None:
            self.labels = labels
        if name is not None:
            self.name = name
        self.verbose = verbose

        if self.verbose:
            print('Type ', self.plotType)

    def xydom(self):
        '''simple return domain '''
        return self.xdom, self.ydom

    def setName(self, name):
        ''' set shape name'''
        self.name = name

    def appendPoint(self, newDate, newY, label):
        ''' append a point to the list for plotting - allow for up to 3
        fields'''
        self.xplot.append(newDate)
        self.labels.append(label)
        if len(newY) <= 0 or len(newY) > 3:
            myerror('shpplot.appendPoint: Invalide number of y points '+newY)
        self.y1plot.append(newY[0])
        if len(newY) >= 2:
            self.y2plot.append(newY[1])
        if len(newY) == 3:
            self.y3plot.append(newY[2])

    def appendPointError(self, newY):
        '''append errors for a point - allow for up to 3 filelds '''
        if len(newY) <= 0 or len(newY) > 3:
            myerror(f'shpplot.appendPointError: Invalid # of y points {newY}')
        self.e1plot.append(newY[0])
        if len(newY) >= 2:
            self.e2plot.append(newY[1])
        if len(newY) == 3:
            self.e3plot.append(newY[2])

    def appendProfile(self, newDate, y1, label, y2=None, y3=None):
        ''' append a point to the list for plotting - allow for up to 3
        fields'''
        self.dates.append(newDate)
        self.labels.append(label)
        #
        self.y1plot.append(y1)
        if len(y1) != len(self.xplot):
            myerror('appendProfile length of x and y data do not match')
        if y2 is not None:
            self.y2plot.append(y2)
        if y3 is not None:
            self.y3plot.append(y3)

    def appendProfileError(self, y1, y2=None, y3=None):
        #
        self.e1plot.append(y1)
        if len(y1) != len(self.xplot):
            myerror('appendProfileError length of x and y data do not match')
        if y2 is not None:
            self.e2plot.append(y2)
        if y3 is not None:
            self.e3plot.append(y3)

    def plotParamHandle(self, markerStyle, nUnique):
        ''' Unpack markers for label/name '''
        # This is the case where a single marker serves all - just replicate
        if type(markerStyle) != list:
            markerStyle = [markerStyle]*nUnique
        else:
            # Just repeat sequence at least enough to make nUnique
            markerStyle = markerStyle*nUnique
            # trim to get exact number needed
            markerStyle = markerStyle[0:nUnique]
            print('msu', markerStyle)
        return markerStyle

    def firstElement(self, x):
        '''If passed a list, just return first element as scalar'''
        if type(x) == list:
            return x[0]
        return x

    def timeSeriesTrend(self, ax=None, color=None, splitTrendByLabel=False,
                        showPlot=True, order=None, timeRange=None):
        #
        if color is None:
            color = 'r'
        pColor = self.firstElement(color)
        y1 = np.array(self.y1plot)
        iGood = np.isfinite(y1)
        ytmp = y1[iGood]
        xtmp = np.array(self.xplot)[iGood]
        ltmp = np.array(self.labels)[iGood]
        # filter out unique labels if byLabel
        if splitTrendByLabel:
            uniqueLabels = np.unique(ltmp)
        else:
            uniqueLabels = [None]
        # override plot order with command line spec
        if order is None:
            uniqueLabels = uniqueLabels[order]
        print(uniqueLabels)
        data = list(zip(xtmp, ytmp, ltmp))
        # sort by label
        # print(data)
        dataSorted = sorted(data, key=lambda x: x[0])
        # print(dataSorted)
        # loop over uniqueLabels - not not byLabel it will go together
        for uniqueLabel in uniqueLabels:
            # unpack sorted data
            x1, y1, labels1 = zip(*dataSorted)
            x1, y1, labels1 = np.array(x1), np.array(y1), np.array(labels1)
            dmid = x1[int(len(x1)/2)]
            dt = np.array([])
            x1use, y1use = np.array([]), np.array([])
            for xel, yel, curLabel in zip(x1, y1, labels1):
                useMe = False
                if uniqueLabel is None:
                    useMe = True
                elif curLabel == uniqueLabel:
                    useMe = True
                if useMe:
                    dt = np.append(dt, (xel-dmid).days)
                    y1use = np.append(y1use, yel)
                    x1use = np.append(x1use, xel)
            if len(y1use) >= 2:
                print(dt)
                print(y1use)
                slope, intercept, r_value, p_value, std_err = \
                    stats.linregress(dt, y1use)
                # fit=np.polyfit(dt, y1, 1)
                # yfit=fit[1] + fit[0]*dt
                yfit = intercept + slope*dt
                print('stats ', r_value, r_value**2, p_value, std_err)
                # print('Fit ', fit)
                if ax is None and showPlot:
                    plt.plot_date(x1use, yfit, '-', color=pColor)
                elif showPlot:
                    ax.plot_date(x1use, yfit, '-', color=pColor,
                                 label='$\mathregular{r^2}=$' +
                                 str(round(r_value**2, 2)) + ',p=' +
                                 str(round(p_value, 2))+',s=' +
                                 str(round(slope*365., 1)))
        if timeRange is not None:
            ax.set_xlim(timeRange)
        return x1use, yfit

    def timeSeries(self, color='r', yIndex=1, markerSize=8, ax=None,
                   markerStyle='o', lineStyle='', lineWidth=1, useLabel=True,
                   deTrend=False, timeRange=None, order=None, yRange=None):
        ''' Plot time series data'''
        # modify here later to plot other
        try:
            y = [self.y1plot, self.y2plot, self.y3plot][yIndex-1]
        except Exception:
            print('Invalid yIndex {0:d}'.format(yIndex))
        # zip,  x,y, and labels
        y1 = np.array(y)
        # cull out the good points
        iGood = np.isfinite(y1)
        ytmp = y1[iGood]
        xtmp = np.array(self.xplot)[iGood]
        ltmp = np.array(self.labels)[iGood]
        data = list(zip(xtmp, ytmp, ltmp))
        # sort by label
        dataSorted = sorted(data, key=lambda x: x[0])
        # unpack sorted data
        x1, y1, labels1 = zip(*dataSorted)
        x1, y1, labels1 = np.array(x1), np.array(y1), np.array(labels1)
        if useLabel:
            # build list and count of unique labels
            uniqueLabels = np.unique(self.labels)
            nUnique = len(uniqueLabels)
            # override plot order with command line spec
            if order is not None:
                uniqueLabels = uniqueLabels[order]
            # print(uniqueLabels)

            # get markStyle, markerSizes
            markerStyles = self.plotParamHandle(markerStyle, nUnique)
            markerSizes = self.plotParamHandle(markerSize, nUnique)
            lineStyles = self.plotParamHandle(lineStyle, nUnique)
            lineWidths = self.plotParamHandle(lineWidth, nUnique)
            colors = self.plotParamHandle(color, nUnique)
            for lab, marker, mSize, lStyle, lWidth, pColor in zip(uniqueLabels,
                                                                  markerStyles,
                                                                  markerSizes,
                                                                  lineStyles,
                                                                  lineWidths,
                                                                  colors):
                ii = np.where(lab == labels1)
                # if detrend, compute fit and remove
                if deTrend:
                    xt, yt = self.timeSeriesTrend(showPlot=False)
                    y1[ii] -= yt
                if ax is None:
                    p = plt.plot_date(x1[ii], y1[ii], color=pColor,
                                      marker=marker, markersize=mSize,
                                      linestyle=lStyle, linewidth=lWidth,
                                      label=self.name+'-'+lab)
                    myAx = plt.gca()
                else:
                    myAx = ax
                    p = ax.plot_date(x1[ii], y1[ii], color=pColor,
                                     marker=marker, markersize=mSize,
                                     linestyle=lStyle, linewidth=lWidth,
                                     label=self.name+'-'+lab)
                    for dd, ss in zip(x1[ii], y1[ii]):
                        print(f'{dd.strftime("%y %m %d")} {ss}')
                    print('----')
                if timeRange is not None:
                    myAx.set_xlim(timeRange)
                if yRange is not None:
                    myAx.set_ylim(yRange)
        #
        # No label case, just plot the whole series
        else:
            # in theory this should not be needed, but wil work if user
            # specifies a list
            marker = self.firstElement(markerStyle)
            mSize = self.firstElement(markerSize)
            lStyle = self.firstElement(lineStyle)
            lWidth = self.firstElement(lineWidth)
            pColor = self.firstElement(color)
            # if detrend, compute fit and remove
            if deTrend:
                xt, yt = self.timeSeriesTrend(showPlot=False)
                y1 -= yt
            if ax is None:
                p = plt.plot_date(x1, y1, color=pColor, marker=marker,
                                  markersize=mSize, linestyle=lStyle,
                                  linewidth=lWidth, label=self.name)
            else:
                p = ax.plot_date(x1, y1, color=pColor, marker=marker,
                                 markersize=mSize, linestyle=lStyle,
                                 linewidth=lWidth, label=self.name)
            ax.set_xlim(timeRange)
        return p

    def timeSeriesError(self, color=None, markerSize=None, ax=None,
                        markerStyle=None):
        # modify here later to plot other
        y = self.y1plot
        e = self.e1plot
        if ax is None:
            for xp, yp, ep in zip(self.xplot, y, e):
                plt.plot_date([xp, xp], [yp-ep, yp+ep], 'k-')
        else:
            for xp, yp, ep in zip(self.xplot, y, e):
                ax.plot_date([xp, xp], [yp-ep, yp+ep], 'k-')
        return

    def profError(self, index, color=None, markerSize=None, ax=None,
                  markerStyle=None):
        # modify here later to plot other
        y = self.y1plot[index]
        e = self.e1plot[index]
        x = self.xplot

        if ax is None:
            for xp, yp, ep in zip(x, y, e):
                # exit()
                plt.plot([xp, xp], [yp-ep, yp+ep], 'k-')
        else:
            for xp, yp, ep in zip(x, y, e):
                ax.plot([xp, xp], [yp-ep, yp+ep], 'k-')
        return

    def profPlot(self, index, color=None, markerSize=None, ax=None,
                 markerStyle=None, lineStyle=None, lineWidth=None,
                 dateFormat=None, useLabel=None):
        print(lineStyle)
        if dateFormat is None:
            dateFormat = "%Y-%m-%d"
        if lineStyle is None:
            lineStyle = '-'
        if lineWidth is None:
            lineWidth = 1
        if color is None:
            color = 'r'
        if markerStyle is None:
            markerStyle = ''
        if markerSize is None:
            markerSize = 12
        lab = ''

        if useLabel is None:
            if self.dates[index] is not None:
                lab = self.dates[index].strftime(dateFormat)
        else:
            lab = self.labels[index]
            if self.dates[index] is not None:
                lab += ':'+self.dates[index].strftime(dateFormat)
        if ax is None:
            print('Label ', self.labels)
            p = plt.plot(self.xplot, self.y1plot[index], color=color,
                         marker=markerStyle, markersize=markerSize,
                         linestyle=lineStyle, linewidth=lineWidth, label=lab)
        else:
            print('Label ', self.labels)
            p = ax.plot(self.xplot, self.y1plot[index], color=color,
                        marker=markerStyle, markersize=markerSize,
                        linestyle=lineStyle, linewidth=lineWidth, label=lab)
        return p

    def readTextTimeSeriesData(self, file):
        #
        # open file
        fpIn = open(file, 'r')
        # init defaults
        tsFormat = None
        label = ''
        ncol = 1
        # this works only for points (at present anyway)
        self.plotType = 'point'
        # loop over lines in file
        for line in fpIn:
            # check for header info
            if 'dateformat' in line:
                tsFormat = line.split('=')[-1].strip()
            elif 'ndatacol' in line:
                ncol = int(line.split('=')[-1])
            elif 'name' in line:
                self.name = line.split('=')[-1].strip()
            elif 'label' in line:
                label = line.split('=')[-1].strip()
            elif ';' in line or line.isspace():
                # comment
                continue
            #
            # process the data
            elif tsFormat is not None:
                # process date
                if 1:
                    dateString = line.split()[0]
                    myDate = datetime.strptime(dateString, tsFormat)
                else:
                    myerror('Invalid date format -- ' + tsFormat +
                            'for date' + dateString + 'from ' + line)
                # process data
                newY = [float(x) for x in line.split()[1:]]
                if len(newY) > 3 or len(newY) < 1 or len(newY) != ncol:
                    myerror('invalid number of columns\n'+line)
                # append data
                self.appendPoint(myDate, newY, label)
            else:
                myerror('error parsing time series\n'+line)
        return
