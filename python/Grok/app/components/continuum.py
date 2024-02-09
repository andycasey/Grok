
# coding:utf-8

import numpy as np
from Qt import QtCore
from Qt.QtWidgets import QVBoxLayout, QHBoxLayout, QGridLayout, QLabel, QSlider, QSpacerItem
from qfluentwidgets import PrimaryPushButton, PushButton, StrongBodyLabel, ComboBox, LineEdit, CheckBox, SpinBox, SwitchButton, HollowHandleStyle, Slider, SwitchButton

from .analysis import AnalysisWidget
from .plot import BasePlotWidget

from Grok.specutils import continuum


class ContinuumRectificationWidget(AnalysisWidget):
    
    def __init__(self, session, callback=None, parent=None, stretch=0):
        # Add buttons to bottom layout
        super().__init__(parent=parent, session=session, title="Continuum Rectification")
        self.callback = callback
        
        layout = QVBoxLayout(self.widget)
        layout.setContentsMargins(0, 0, 0, 0)  
        
        self.figure = BasePlotWidget(
            figsize=(100, 5), 
            toolbar=True, 
            toolbar_left_right=True
        )
        self.figure.axes = ax_rect, ax_orig = self.figure.figure.subplots(2, 1, gridspec_kw={"height_ratios": [2, 1]})
        ax_rect.set_xticks([])
        ax_orig.set_xlabel("Wavelength [Å]")
        ax_orig.set_ylabel("Flux")
        ax_rect.set_ylabel("Rectified flux")
        #self.figure.tight_layout()
        self.figure.canvas.draw()
        
    
        self.current_index = 0
        
        self._line, = ax_orig.plot([], [])
        self._line_continuum, = ax_orig.plot([], [], c="tab:red")
        self._line_rect, = ax_rect.plot([], [])
        self._line_H, = ax_rect.plot([], [], c="tab:green")
        

        layout.addWidget(self.figure)          
        
        control_layout = QHBoxLayout()
        # use a grid layout
        grid_layout = QGridLayout()
        
        # create a combo box for the continuum type
        #Lengthgrid_layout.addWidget(QLabel(""))
            
        self.continuum_degree_label = QLabel("Degree:")
        self.continuum_degree = SpinBox(self)
        self.continuum_degree.setRange(0, 100)
        self.continuum_degree.setMaximumWidth(120)
                
        grid_layout.addWidget(QLabel("Sines-and-cosine options"), 0, 0, 1, 4, QtCore.Qt.AlignLeft)
                
        grid_layout.addWidget(self.continuum_degree_label, 1, 0)
        grid_layout.addWidget(self.continuum_degree, 1, 1)
        
        self.continuum_length_scale = LineEdit()
        self.continuum_length_scale.setMaximumWidth(60)
        self.continuum_length_scale_label = QLabel("Length scale:")
        grid_layout.addWidget(self.continuum_length_scale_label, 1, 2)
        grid_layout.addWidget(self.continuum_length_scale, 1, 3, alignment=QtCore.Qt.AlignRight)
        
        
        col_reg = 5

        spacer = QLabel(" ")
        spacer.setMinimumWidth(30)
        grid_layout.addWidget(spacer, 0, 4)
        
        alpha_max_width = 60
        self.alpha_W = LineEdit()
        self.alpha_W.setMaximumWidth(alpha_max_width)
        self.alpha_H = LineEdit()
        self.alpha_H.setMaximumWidth(alpha_max_width)
        grid_layout.addWidget(QLabel("Regularization"), 0, col_reg, 1, 5, QtCore.Qt.AlignLeft)
        
        grid_layout.addWidget(QLabel("On W:"), 1, col_reg)
        grid_layout.addWidget(self.alpha_W, 1, col_reg + 1, alignment=QtCore.Qt.AlignRight)
        grid_layout.addWidget(QLabel("On H:"), 1, col_reg + 2)
        grid_layout.addWidget(self.alpha_H, 1, col_reg + 3, alignment=QtCore.Qt.AlignRight)
        
        
        self.l1_l2_ratio = Slider(QtCore.Qt.Horizontal, self)
        self.l1_l2_ratio.setRange(0, 100)
        self.l1_l2_ratio.setTickInterval(25)
        self.l1_l2_ratio.setMaximumWidth(120)
        l1_l2_layout = QHBoxLayout()
        l1_l2_layout.addWidget(QLabel("L1"))
        l1_l2_layout.addWidget(self.l1_l2_ratio)
        l1_l2_layout.addWidget(QLabel("L2"))
        l1_l2_layout.addStretch(1)
        grid_layout.addLayout(l1_l2_layout, 1, col_reg + 4, alignment=QtCore.Qt.AlignLeft)#3, 1, 2, QtCore.Qt.AlignHCenter)
        
        control_layout.addLayout(grid_layout)
        
        layout.addLayout(control_layout)
        
        self.button_do_not_rectify = PushButton("No continuum normalization")
        self.button_rectify = PrimaryPushButton("Continuum normalize all spectra")
        
        self.button_do_not_rectify.clicked.connect(self.button_do_not_rectify_clicked)
        self.button_rectify.clicked.connect(self.button_rectify_clicked)
        self.bottomLayout.addWidget(self.button_do_not_rectify)
        self.bottomLayout.addStretch(1)
        self.bottomLayout.addWidget(self.button_rectify)

        self.widget.setParent(self.card)
        self.widget.show()

        
        
        self.figure.action_left.triggered.connect(lambda: self.update_canvas(self.current_index - 1))
        self.figure.action_right.triggered.connect(lambda: self.update_canvas(self.current_index + 1))
        
        self.update_canvas(self.current_index)
        
        self.alpha_W.textChanged.connect(self.redraw)
        self.alpha_H.textChanged.connect(self.redraw)
        self.l1_l2_ratio.valueChanged.connect(self.redraw)
        self.continuum_length_scale.textChanged.connect(self.redraw)
        #self.continuum_type.currentTextChanged.connect(self.hide_show_continuum_controls)
        self.continuum_degree.valueChanged.connect(self.redraw)
        
        
        #self.hide_show_continuum_controls()
        
        return None
    
    def redraw(self, *args):
        print(f"redraw {args}")
        self.update_canvas(self.current_index)
    
    def update_canvas(self, index):
        #nonlocal current_index
        S = len(self.session.spectra)
        
        if index < 0 or index > (S - 1):
            return
        wavelength, flux, ivar, meta = self.session.spectra[index]
        
        try:
            L = float(self.continuum_length_scale.text())
        except:
            L = 2 * np.ptp(wavelength)
        print(f"L is {L}")
        try:
            alpha_W = float(self.alpha_W.text())
            alpha_H = float(self.alpha_H.text())
            l1_ratio = self.l1_l2_ratio.value() / 100
        except:
            import logging
            logging.exception("wha")
            return None
        else:
            print(alpha_W, alpha_H, l1_ratio)
            (W, H, theta, continuum_, chi2, dof) = continuum.fit_nmf_sinusoids(
                wavelength, 
                flux, 
                ivar, 
                int(self.continuum_degree.value()), 
                L,
                alpha_W=alpha_W,
                alpha_H=alpha_H,
                l1_ratio=l1_ratio
            )
            
            self._line.set_data(wavelength, flux)
        self._line_rect.set_data(wavelength, flux / continuum_)
        self._line_H.set_data(wavelength, 1 - (W @ H).flatten())
        self._line_continuum.set_data(wavelength, continuum_)
        
        for ax in self.figure.axes:
            ax.set_xlim(wavelength[[0, -1]])
                    
        ylims = (np.nanmin(flux), np.nanmax(flux))
        edge = 0.10 * np.ptp(ylims)
        ylims = [ylims[0] - edge, ylims[1] + edge]
        # limit to 0
        ylims[0] = max(ylims[0], 0)
        self.figure.axes[0].set_ylim(0, 1.2)
        self.figure.axes[1].set_ylim(ylims)
        self.figure.canvas.draw()
        self.current_index = index
        self.figure.action_left.setEnabled(self.current_index > 0)
        self.figure.action_right.setEnabled(self.current_index < (S - 1))
        #self.figure.canvas.setFocus()
        
        
    def button_do_not_rectify_clicked(self):
        
        if self.callback is not None:
            self.callback()
        
    def button_rectify_clicked(self):
        if self.callback is not None:
            self.callback()
            
