
# coding:utf-8

import numpy as np
from Qt import QtCore
from Qt.QtGui import QDoubleValidator
from Qt.QtWidgets import (QHBoxLayout, QLabel, QVBoxLayout)
from qfluentwidgets import (ComboBox, LineEdit, PrimaryPushButton, PushButton)

from ..view.gallery_interface import SeparatorWidget

from .analysis import AnalysisWidget
from .plot import SinglePlotWidget


class RadialVelocityWidget(AnalysisWidget):
    
    
    def __init__(self, session, callback=None, parent=None, stretch=0):
        
        # Add buttons to bottom layout
        super().__init__(parent=parent, session=session, title="Radial Velocity")
        self.session = session
        self.callback = callback
        '''
        layout = QVBoxLayout(self.widget)
        layout.addWidget(QLabel("Test"))
        
        
        self.button_do_not_rectify = PushButton("No continuum normalization")
        self.button_rectify = PrimaryPushButton("Rectify")
        self.bottomLayout.addWidget(self.button_do_not_rectify)
        self.bottomLayout.addStretch(1)
        self.bottomLayout.addWidget(self.button_rectify)
        '''
        
        # Get the order closest to the wavelength range.
        index = session._get_closest_spectrum_index(8500)                
        wavelength, flux, *_ = session.spectra[index]
        
        
        layout = QVBoxLayout(self.widget)
        layout.setContentsMargins(0, 0, 0, 0)  
        
        #index = 28
        plot = SinglePlotWidget(
            x=wavelength,
            y=flux / np.nanmedian(flux),            
            xlabel=r"Wavelength [$\mathrm{\AA}$]",
            ylabel="Rectified flux",
            figsize=(100, 2), toolbar=False, toolbar_left_right=False)
        plot.ax.set_xlim(wavelength[[0, -1]])
        
        layout.addWidget(plot)

        
        
        rv_opt_layout = QHBoxLayout()
        rv_opt_layout.setSpacing(30)
        
        lhs_layout = QVBoxLayout()
        middle_layout = QVBoxLayout()
        rhs_layout = QVBoxLayout()
        
        
        wavelength_layout = QHBoxLayout()
        wavelength_layout.addWidget(QLabel("Wavelength range:"))
        wavelength_combo = ComboBox()
        wavelength_combo.addItems([
            "8000 - 9000 Å"
        ])
        wavelength_layout.addWidget(wavelength_combo)


        template_layout = QHBoxLayout()
        template_combo = ComboBox()
        template_combo.setPlaceholderText("Auto")
        template_combo.addItems([
            "Auto",
            "Select from file",
            "Synthesise spectrum"
        ])
        for combo in (wavelength_combo, template_combo):
            combo.setFixedWidth(180)
        template_layout.addWidget(QLabel("Template:"))
        template_layout.addWidget(template_combo)

        lhs_layout.addLayout(wavelength_layout)
        lhs_layout.addLayout(template_layout)

    
        # Middle layout
        continuum_layout = QHBoxLayout()
        continuum_layout.setAlignment(QtCore.Qt.AlignTop)
        continuum_layout.addWidget(QLabel("Continuum:"))
        continuum_combo = ComboBox()
        continuum_combo.setPlaceholderText("Sinusoids")
        continuum_combo.addItems([
            "Sinusoids",
            "Polynomial",
            "Spline",
        ])
        continuum_layout.addWidget(continuum_combo)        
        middle_layout.addLayout(continuum_layout)
        

        rv_opt_layout.addLayout(lhs_layout)
        rv_opt_layout.addLayout(middle_layout)
        rv_opt_layout.addStretch(1)
        rv_opt_layout.addLayout(rhs_layout)    
        
        layout.addLayout(rv_opt_layout)
                
        self.button_no_shift = PushButton("Spectrum is already at rest")
        self.button_no_shift.clicked.connect(self.on_no_shift_button_pushed)
        self.button_measure = PushButton("Measure")
        self.button_measure.clicked.connect(self.on_shift_button_clicked)
        self.button_measure_and_shift = PrimaryPushButton("Measure and shift")
        self.button_measure_and_shift.clicked.connect(self.on_measure_and_shift_button_clicked)
        self.bottomLayout.addWidget(self.button_no_shift)
        self.bottomLayout.addStretch(1)
        self.bottomLayout.addWidget(QLabel("Relative velocity:"))
        
        self.v_rel = LineEdit(self)
        self.v_rel.setValidator(QDoubleValidator(-1000, 1000, 3, self))
        self.v_rel.setText("0.0")
        self.v_rel.setFixedWidth(80)
        self.v_rel.setAlignment(QtCore.Qt.AlignHCenter)
        
        self.v_rel.textEdited.connect(self.on_v_rel_changed)
        
        self.bottomLayout.addWidget(self.v_rel)
        self.bottomLayout.addWidget(QLabel("km/s"))
        self.bottomLayout.addWidget(SeparatorWidget())
        self.bottomLayout.addWidget(self.button_measure)
        self.bottomLayout.addWidget(self.button_measure_and_shift)
        
        self.button_measure.setVisible(False)
        #measure_layout = QHBoxLayout()
        #measure_layout.addStretch(1)
        #measure_layout.addWidget(PushButton("Measure relative velocity"), 0, QtCore.Qt.AlignRight)
        #layout.addLayout(measure_layout)
                        
        return None    
    
    def on_no_shift_button_pushed(self):
        self.v_rel.setText("0.0")
        if self.callback is not None:
            self.callback()
            
    def on_measure_and_shift_button_clicked(self):
        print("on measure clicked")
        if self.callback is not None:
            self.callback()
            
    def on_shift_button_clicked(self):
        print("on shift clicked")
        if self.callback is not None:
            self.callback()
    
    def on_v_rel_changed(self):
        
        state = self.v_rel.validator().validate(self.v_rel.text(), 0)[0]
        if state == QDoubleValidator.Acceptable:    
            self.button_measure.setVisible(True)
            self.button_measure_and_shift.setText("Shift")
        else:
            print("warning")
        print("edited")