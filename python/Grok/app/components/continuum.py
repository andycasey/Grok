
# coding:utf-8

import numpy as np
from Qt.QtWidgets import QVBoxLayout
from qfluentwidgets import PrimaryPushButton, PushButton

from .analysis import AnalysisWidget
from .plot import SinglePlotWidget


class ContinuumRectificationWidget(AnalysisWidget):
    
    def __init__(self, session, callback=None, parent=None, stretch=0):
        # Add buttons to bottom layout
        super().__init__(parent=parent, session=session, title="Continuum Rectification")
        self.callback = callback
        
        layout = QVBoxLayout(self.widget)
        layout.setContentsMargins(0, 0, 0, 0)  
        
        figure = SinglePlotWidget(
            xlabel=r"Wavelength [$\mathrm{\AA}$]",
            ylabel="Rectified flux",
            figsize=(100, 3), 
            toolbar=True, 
            toolbar_left_right=True
        )
        current_index = 0
        
        line, = figure.ax.plot([], [])
        
        S = len(session.spectra)
        
        def update_canvas(index):
            nonlocal current_index
            if index < 0 or index > (S - 1):
                return
            wavelength, flux, ivar, meta = session.spectra[index]
            line.set_data(wavelength, flux / np.nanmedian(flux))
            figure.ax.set_xlim(wavelength[[0, -1]])
            figure.ax.set_ylim(0, 1.2)
            figure.canvas.draw()
            current_index = index
            figure.action_left.setEnabled(current_index > 0)
            figure.action_right.setEnabled(current_index < (S - 1))
            figure.canvas.setFocus()
            
        
        figure.action_left.triggered.connect(lambda: update_canvas(current_index - 1))
        figure.action_right.triggered.connect(lambda: update_canvas(current_index + 1))
        
        update_canvas(current_index)

        layout.addWidget(figure)                
        
        self.button_do_not_rectify = PushButton("No continuum normalization")
        self.button_rectify = PrimaryPushButton("Continuum normalize all spectra")
        
        self.button_do_not_rectify.clicked.connect(self.button_do_not_rectify_clicked)
        self.button_rectify.clicked.connect(self.button_rectify_clicked)
        self.bottomLayout.addWidget(self.button_do_not_rectify)
        self.bottomLayout.addStretch(1)
        self.bottomLayout.addWidget(self.button_rectify)

        self.widget.setParent(self.card)
        self.widget.show()
        
        return None
    
    def button_do_not_rectify_clicked(self):
        if self.callback is not None:
            self.callback()
        
    def button_rectify_clicked(self):
        if self.callback is not None:
            self.callback()
            
