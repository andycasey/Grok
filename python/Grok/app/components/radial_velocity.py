
# coding:utf-8
import os
import numpy as np
from Qt import QtCore
from Qt.QtGui import QDoubleValidator
from Qt.QtWidgets import (QHBoxLayout, QLabel, QVBoxLayout, QFileDialog)
from qfluentwidgets import (ComboBox, LineEdit, PrimaryPushButton, PushButton)

from ..common.config import cfg, qconfig
from ..view.gallery_interface import SeparatorWidget

from .analysis import AnalysisWidget
from .plot import SinglePlotWidget

from specutils import Spectrum1D, continuum, rv

def _parse_wavelength_range(text):
    return tuple(map(float, map(str.strip, text.strip("Å").split("-"))))

class RadialVelocityWidget(AnalysisWidget):
    
    
    def __init__(self, session, callback=None, parent=None, stretch=0):
        
        # Add buttons to bottom layout
        super().__init__(parent=parent, session=session, title="Radial Velocity")
        self.session = session
        if callback is None:
            callback = lambda: None
        self.callback = callback
        
        layout = QVBoxLayout(self.widget)
        layout.setContentsMargins(0, 0, 0, 0)  
    
        wavelength_range = _parse_wavelength_range(cfg.get(cfg.rv_default_wavelength_range))
    
        self.plot = SinglePlotWidget(
            xlabel=r"Wavelength [$\mathrm{\AA}$]",
            ylabel="Rectified flux",
            figsize=(100, 2), 
            toolbar=True, 
            toolbar_left_right=False
        )
        self.plot.ax.set_xlim(wavelength_range)
        self.plot.ax.set_ylim(0, 1.2)
        self._plot_observed, = self.plot.ax.plot([], [], c="tab:blue")
        self._plot_template, = self.plot.ax.plot([], [], c="tab:red")
                
        layout.addWidget(self.plot)
        
        rv_opt_layout = QHBoxLayout()
        rv_opt_layout.setSpacing(30)
        
        lhs_layout = QVBoxLayout()
        middle_layout = QVBoxLayout()
        #rhs_layout = QVBoxLayout()
        
        wavelength_layout = QHBoxLayout()
        wavelength_layout.addWidget(QLabel("Wavelength range:"))
        self.combo_wavelength_range = ComboBox()
        self.combo_wavelength_range.addItems(cfg.rv_default_wavelength_range.options)
        self.combo_wavelength_range.setCurrentText(cfg.get(cfg.rv_default_wavelength_range))
        wavelength_layout.addWidget(self.combo_wavelength_range)


        template_layout = QHBoxLayout()
        template_combo = ComboBox()
        template_combo.addItems(cfg.rv_default_template_type.options)
        template_combo.setCurrentText(cfg.get(cfg.rv_default_template_type))

        for combo in (self.combo_wavelength_range, template_combo):
            combo.setFixedWidth(180)

            
        template_layout.addWidget(QLabel("Template:"))
        template_layout.addWidget(template_combo)
        
        button_select_template = PushButton("Select")

        self.template_path_label = QLabel("")
        
        self.set_template_path(cfg.get(cfg.rv_default_template_path))
        
        template_path_layout = QHBoxLayout()
        template_path_layout.addSpacing(20)
        template_path_layout.addWidget(self.template_path_label)
        template_path_layout.addStretch(1)
        template_path_layout.addWidget(button_select_template)



        lhs_layout.addLayout(wavelength_layout)
        lhs_layout.addLayout(template_layout)
        lhs_layout.addLayout(template_path_layout)
        
        #lhs_layout.addLayout(template_path_layout)

    
        # Middle layout
        continuum_layout = QHBoxLayout()
        continuum_layout.setAlignment(QtCore.Qt.AlignTop)
        continuum_layout.addWidget(QLabel("Continuum:"))
        self.combo_continuum_type = ComboBox()
        self.combo_continuum_type.addItems(cfg.rv_default_continuum_type.options)
        self.combo_continuum_type.setCurrentText(cfg.get(cfg.rv_default_continuum_type))        
        continuum_layout.addWidget(self.combo_continuum_type)        
        
        degree_layout = QHBoxLayout()
        degree_layout.addSpacing(20)
        degree_layout.addWidget(QLabel("Degree:"))
        degree_layout.addStretch(1)
        self.combo_continuum_degree = ComboBox()
        self.combo_continuum_degree.addItems(tuple(map(str, range(10))))
        self.combo_continuum_degree.setCurrentText(cfg.get(cfg.rv_default_continuum_degree))
        degree_layout.addWidget(self.combo_continuum_degree)
        
        middle_layout.addLayout(continuum_layout)
        middle_layout.addLayout(degree_layout)
        middle_layout.addStretch(1)
        

        rv_opt_layout.addLayout(lhs_layout)
        rv_opt_layout.addLayout(middle_layout)
        
        layout.addLayout(rv_opt_layout)
                
        self.button_no_shift = PushButton("Spectrum is already at rest")
        self.button_no_shift.clicked.connect(self.on_no_shift_button_pushed)
        self.button_measure = PushButton("Measure")
        self.button_measure.clicked.connect(self.on_measure_clicked)
        self.button_measure_and_shift = PrimaryPushButton("Measure and shift")
        self.button_measure_and_shift.clicked.connect(self.on_measure_and_shift_button_clicked)
        self.bottomLayout.addWidget(self.button_no_shift)
        self.bottomLayout.addStretch(1)
        self.bottomLayout.addWidget(QLabel("Relative velocity:"))
        
        self.v_rel = LineEdit(self)
        self.v_rel.setValidator(QDoubleValidator(-1000, 1000, 3, self))
        self.v_rel.setFixedWidth(80)
        self.v_rel.setAlignment(QtCore.Qt.AlignHCenter)
        
        self.v_rel.textEdited.connect(self.on_v_rel_changed)
        
        self.bottomLayout.addWidget(self.v_rel)
        self.bottomLayout.addWidget(QLabel("km/s"))
        self.bottomLayout.addWidget(SeparatorWidget())
        self.bottomLayout.addWidget(self.button_measure)
        self.bottomLayout.addWidget(self.button_measure_and_shift)
        
        self.button_measure.setVisible(False)
        self._separate_measure_and_shift = False
        
        self.combo_wavelength_range.currentTextChanged.connect(self.on_wavelength_range_changed)
        self.combo_continuum_degree.currentTextChanged.connect(self.fit_continuum_and_update_current_spectrum_plot)
        self.combo_continuum_type.currentTextChanged.connect(self.fit_continuum_and_update_current_spectrum_plot)
        button_select_template.clicked.connect(self.on_select_template_clicked)
        
        # Load the spectrum
        self.current_spectrum = {}

        index = self.get_nearest_index(wavelength_range)
        self.get_spectrum(index)
        
        # get current v_rel from session        
        v_rel = self.session.get_v_rel(index)
        self.v_rel.setText(f"{v_rel:.1f}")
        self.fit_continuum()
        self.update_current_spectrum_plot()
        return None
    
    def get_nearest_index(self, wavelength_range):
        return self.session._get_closest_spectrum_index(np.mean(wavelength_range))

    
    def get_spectrum(self, index):
        # Get the order closest to the wavelength range.
                
        # Load the spectrum.
        (wavelength, flux, ivar, *_) = self.session.get_spectrum(index, rest_frame=False, rectified=False)

        self.current_spectrum.update(
            wavelength=wavelength,
            flux=flux,
            ivar=ivar,
        )

    def update_current_spectrum_plot(self):   
        try:
            v_rel = float(self.v_rel.text())     
        except:
            v_rel = 0
        wavelength = self.current_spectrum["wavelength"]
        flux = self.current_spectrum["flux"]
        continuum = self.current_spectrum["continuum"]
        
        self._plot_observed.set_data(
            wavelength * (1 - v_rel / 299792.458),
            flux / continuum
        )
        self.plot.canvas.draw()
        
    
    def fit_continuum(self):
        
        deg = int(self.combo_continuum_degree.currentText().strip())
        wl, flux, ivar = (self.current_spectrum["wavelength"], self.current_spectrum["flux"], self.current_spectrum["ivar"])
        c, theta = continuum.fit_polynomial(wl, flux, ivar, deg=deg)
        self.current_spectrum["continuum"] = c
    
    
    def fit_continuum_and_update_current_spectrum_plot(self):
        self.fit_continuum()
        self.update_current_spectrum_plot()
        
    
    
    

    def on_wavelength_range_changed(self):        
        wavelength_range = _parse_wavelength_range(self.combo_wavelength_range.currentText())        
        self.plot.ax.set_xlim(wavelength_range)
        self.plot.reset_current_as_home()
        
        index = self.get_nearest_index(wavelength_range)
        
        self.get_spectrum(index)
        self.fit_continuum()
        self.update_current_spectrum_plot()
        
    
    def set_template_path(self, path, max_label_length=50):
        if not path: 
            self.template_path = ""
            return None
            
        self.template_spectrum = Spectrum1D.read(path)                
        self.template_path = path
        self._plot_template.set_data(self.template_spectrum.wavelength, self.template_spectrum.flux)

        # Show a shortened version of the path
        if len(path) > max_label_length:
            basename = os.path.basename(path)
            directory = os.path.dirname(path)
            D = max_label_length - len(basename) - 3
            label = f"{directory[:D]}.../{basename}"
        else:
            label = path
            
        self.template_path_label.setText(label)
        self.plot.canvas.draw()
        if not cfg.get(cfg.rv_default_template_path):
            qconfig.set(
                cfg.rv_default_template_path,
                path
            )
        
        
    def on_select_template_clicked(self):
        try:
            directory = os.path.dirname(cfg.get(cfg.rv_default_template_path))
        except:
            directory = ""
            
        path, _ = QFileDialog.getOpenFileName(
            caption="Select template spectrum", 
            directory=directory, 
            filter="*.fits"
        )
        if path:
            self.set_template_path(path)

        
    def set_v_rel(self):
        try:
            v_rel = float(self.v_rel.text())
        except:
            v_rel = 0.0
        finally:
            self.session.set_v_rel(v_rel)
            self.v_rel.setText(f"{v_rel:.1f}")
        return None
        
    def on_no_shift_button_pushed(self):
        self.v_rel.setText("0")
        self.update_current_spectrum_plot()
        self.set_v_rel()
        self.callback()
            
    def on_measure_and_shift_button_clicked(self):        
        self.set_v_rel()
        if not self._separate_measure_and_shift:
            self.on_measure_clicked()
        self.update_current_spectrum_plot()
        self.callback()
            
    def on_measure_clicked(self):        
        # TODO: Put this to Session object?
        v_rel, *_ = rv.measure_relative_velocity(
            self.current_spectrum["wavelength"],
            self.current_spectrum["flux"] / self.current_spectrum["continuum"],
            self.template_spectrum.wavelength, self.template_spectrum.flux            
        )
        self.v_rel.setText(f"{v_rel:.1f}")
        self.update_current_spectrum_plot()
    
    def on_v_rel_changed(self):       
        state = self.v_rel.validator().validate(self.v_rel.text(), 0)[0]
        if state == QDoubleValidator.Acceptable:
            self.update_current_spectrum_plot()
            if not self._separate_measure_and_shift:    
                self.button_measure.setVisible(True)
                self.button_measure_and_shift.setText("Shift")
                self._separate_measure_and_shift = True
        
    