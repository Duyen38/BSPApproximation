import matplotlib.pyplot as plt
import matplotlib.widgets as widgets
import numpy as np

# Defining the cursor
class SnaptoCursor():
    def __init__(self, ax, x, y, y_point):
        self.ax = ax        
        self.lx = ax.axhline(color='k', alpha=0.2)  # the horiz line
        self.ly = ax.axvline(color='k', alpha=0.2)  # the vert line
        self.marker = ax.annotate("^", xy=(0,0), xytext=(10,15),textcoords="offset points",
                    bbox=dict(boxstyle="round", fc="w"),
                    arrowprops=dict(arrowstyle="->", ec="k", lw=1))
        
        self.marker.set_visible(False)
        self.x = x
        self.y = y
        self.y_point = y_point
        self.txt = ax.text(0.7, 0.9, '')

    def mouse_move(self, event):
        if not event.inaxes: 
            return
        x, y = event.xdata, event.ydata
        indx = min(np.searchsorted(self.x, x), len(self.x) - 1)
        x = self.x[indx]
        y = self.y[indx]
        delta_y = self.y[indx] - self.y_point[indx]
        # update the line positions
        self.lx.set_ydata(y)
        self.ly.set_xdata(x)
        self.marker.xy = (x,y)
        self.txt ='(%1.2f, %1.2f)\ndelta y=%1.2f' % (x, y, delta_y)
        self.marker.set_text(self.txt)
        self.marker.set_visible(True)
        self.ax.figure.canvas.draw()
        
# x = np.arange(0.0, 1.0, 0.01)
# y = np.sin(2*2*np.pi*t)
# fig, ax = plt.subplots()

# #cursor = Cursor(ax)
# cursor = SnaptoCursor(ax, x, y)
# cid =  plt.connect('motion_notify_event', cursor.mouse_move)

# ax.plot(x, y,)
# plt.axis([0, 1, -1, 1])
# plt.show()
