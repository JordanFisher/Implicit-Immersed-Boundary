from pylab import rc, gca, grid

class Parameters:
    def __init__(self):
        self.FontSize = 18
        self.TickFontSize = 15
        self.Grid = True
        self.GridWidth = 1

    def DefaultGrid(self):
        if self.Grid:
            grid(True,linewidth=self.GridWidth)
        
def StartLatex():
    rc('text', usetex=True)
    rc('font', family='serif')

def SetTickFont(FontSize):
    ax = gca()
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(FontSize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(FontSize)
