#!/usr/bin/python
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from matplotlib.ticker import FormatStrFormatter
import numpy as np

class MplCanvas(FigureCanvasQTAgg):

    def __init__(self, parent=None, width=7, height=10, dpi=100, fig=None):
        """

        :type fig: Figure
        """
        if fig:
            self.fig = fig.axes.figure
        else:
            self.fig = Figure(figsize=(width, height), dpi=dpi)
            self.axes = self.fig.add_subplot(111)
            self.noDataText = self.fig.text(0.4, 0.5, "No Data", bbox={'facecolor':'white','alpha':1,'edgecolor':'none','pad':1})
        super(MplCanvas, self).__init__(self.fig)
        self.colorList = ["red", "blue", "green", "black", "purple"]

    def updateFig(self, rho, yvals, s=10, legend=None, keepLims=None, xFormatter: FormatStrFormatter = None,
                  yFormatter: FormatStrFormatter=None, color="black", *args, **kwargs):
        try:
            for txt in self.fig.texts:
                txt.set_visible(False)
        except:
            pass
        if keepLims:
            xlim = self.axes.get_xlim()
            ylim = self.axes.get_ylim()
        self.axes.cla()

        self.axes.set_xticklabels(self.axes.get_xticks(), size=8)
        self.axes.set_yticklabels(self.axes.get_yticks(), size=8)

        if xFormatter:
            self.axes.xaxis.set_major_formatter(xFormatter)
        else:
            self.axes.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))

        if yFormatter:
            self.axes.yaxis.set_major_formatter(yFormatter)
        else:
            self.axes.yaxis.set_major_formatter(FormatStrFormatter('%.2E'))

        if type(yvals) == list:
            for n, yval in enumerate(yvals):
                self.axes.scatter(rho, yval, s=s, color=self.colorList[n])
        elif type(yvals) == np.ndarray:
            if kwargs.get("marker"):
                self.axes.scatter(rho, yvals, s=s, color=color, marker=kwargs.get("marker"))
            else:
                self.axes.scatter(rho, yvals, s=s, color=color)
        elif yvals == None:
            pass

        if legend:
            self.fig.legend(legend, fontsize=8)



        # for c in fig.axes.get_children():
        #      if isinstance(c, matplotlib.collections.PathCollection):
        #          offsets = np.array(c.get_offsets()).T
        #          self.axes.scatter(offsets[0], offsets[1], s=3)
        #self.fig.tight_layout(pad=.9)
        if keepLims:
            self.axes.set_xlim(xlim)
            self.axes.set_ylim(ylim)
        self.fig.canvas.draw()
        self.draw_idle()

    def add_scatter(self, x, y, color, s=10, *args, **kwargs):
        self.axes.scatter(x, y, color=color, s=s, marker=kwargs.get("marker"))
        self.fig.canvas.draw()
        self.draw_idle()

    def updateXMin(self, val):
        xmax = self.axes.get_xlim()[1]
        self.axes.set_xlim(val, xmax)
        self.draw_idle()

    def updateXMax(self, val):
        xmin = self.axes.get_xlim()[0]
        self.axes.set_xlim(xmin, val)
        self.draw_idle()

    def updateYMin(self, val):
        ymax = self.axes.get_ylim()[1]
        self.axes.set_ylim(val, ymax)
        self.draw_idle()

    def updateYMax(self, val):
        ymin = self.axes.get_ylim()[0]
        self.axes.set_ylim(ymin, val)
        self.draw_idle()

    def _saveFig(self, f):
        self.fig.savefig(f, dpi=300, format="png")