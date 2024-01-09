
import numpy as np
from Qt.QtWidgets import QVBoxLayout

from .analysis import AnalysisWidget
from .plot import SinglePlotWidget

class QuickViewWidget(AnalysisWidget):
    
    def __init__(self, session, callback=None, parent=None, stretch=0):
        # Add buttons to bottom layout
        super().__init__(parent=parent, session=session, title="Quick view")
        self.callback = callback
        
        layout = QVBoxLayout(self.widget)
        layout.setContentsMargins(0, 0, 0, 0)  

        plotter = SinglePlotWidget(
            xlabel=r"Wavelength [$\mathrm{\AA}$]",
            ylabel="Flux",
            figsize=(100, 3), 
            toolbar=True, 
            toolbar_left_right=True
        )
        current_index = 0
        
        line, = plotter.ax.plot([], [])
        
        S = len(session.spectra)
        
        def update_canvas(index):
            nonlocal current_index
            if index < 0 or index > (S - 1):
                return
            wavelength, flux, ivar, meta = session.spectra[index]
            line.set_data(wavelength, flux)
            plotter.ax.set_xlim(wavelength[[0, -1]])
            ylims = (np.nanmin(flux), np.nanmax(flux))
            edge = 0.05 * np.ptp(ylims)
            ylims = [ylims[0] - edge, ylims[1] + edge]
            # limit to 0
            ylims[0] = max(ylims[0], 0)
            
            plotter.ax.set_ylim(ylims)
            plotter.figure.tight_layout()
            plotter.canvas.draw()
            current_index = index
            plotter.action_left.setEnabled(current_index > 0)
            plotter.action_right.setEnabled(current_index < (S - 1))
            plotter.canvas.setFocus()
            # set a new nav stack so that when the user pans around on THIS order, the 'home'
            # button will take them to the original view for THIS order, not all orders
            plotter.reset_current_as_home()
                        
        plotter.action_left.triggered.connect(lambda: update_canvas(current_index - 1))
        plotter.action_right.triggered.connect(lambda: update_canvas(current_index + 1))
        
        update_canvas(current_index)

        layout.addWidget(plotter)                
        
        '''
        self.button_do_not_rectify = PushButton("No continuum normalization")
        self.button_rectify = PrimaryPushButton("Continuum normalize all spectra")
        
        self.button_do_not_rectify.clicked.connect(self.button_do_not_rectify_clicked)
        self.button_rectify.clicked.connect(self.button_rectify_clicked)
        self.bottomLayout.addWidget(self.button_do_not_rectify)
        self.bottomLayout.addStretch(1)
        self.bottomLayout.addWidget(self.button_rectify)
        '''
        
        self.widget.setParent(self.card)
        self.widget.show()
        
        return None
    
    def button_do_not_rectify_clicked(self):
        if self.callback is not None:
            self.callback()
        
    def button_rectify_clicked(self):
        if self.callback is not None:
            self.callback()
            


    
'''
if __name__ == '__main__':
    # enable dpi scale
    import sys
    from Qt import QtCore
    QtCore.QApplication.setHighDpiScaleFactorRoundingPolicy(QtCore.Qt.HighDpiScaleFactorRoundingPolicy.PassThrough)
    QtCore.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling)
    QtCore.QApplication.setAttribute(QtCore.Qt.AA_UseHighDpiPixmaps)

    app = QtCore.QApplication(sys.argv)
    w = QuickViewWidget()
    w.resize(600, 600)
    w.show()
    app.exec_()    
'''    