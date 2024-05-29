
import numpy as np
from Qt.QtWidgets import QVBoxLayout
from Qt import QtCore

from Qt.QtWidgets import QVBoxLayout, QHBoxLayout, QGridLayout, QLabel, QSlider, QSpacerItem
from Qt.QtWidgets import QSizePolicy
from qfluentwidgets import PrimaryPushButton, PushButton, StrongBodyLabel, ComboBox, LineEdit, CheckBox, SpinBox, SwitchButton, HollowHandleStyle, Slider, SwitchButton

from .analysis import AnalysisWidget
from .plot import SinglePlotWidget
from .plot import BasePlotWidget

import h5py as h5
from Grok.models.basis import NMFSpectralBasis, FourierBasis
from Grok.models.stellar_spectrum import LinearStellarSpectrumModel


class LinearSpectralModelWidget(AnalysisWidget):
    
    def __init__(self, session, callback=None, parent=None, stretch=0):

        # Add buttons to bottom layout
        super().__init__(parent=parent, session=session, title="Quick fit")
        self.callback = callback
        
        layout = QVBoxLayout(self.widget)
        layout.setContentsMargins(0, 0, 0, 0)  
        
        self.figure = BasePlotWidget(
            figsize=(100, 5), 
            toolbar=True, 
            toolbar_left_right=True,            
        )
        
        self.figure.toolbar_layout.addStretch(1)            
        self.action_fit_current_order = PushButton("Fit")        
        self.action_fit_all_orders = PrimaryPushButton("Fit all")
        for button in (self.action_fit_current_order, self.action_fit_all_orders):
            button.setMinimumWidth(100)
            self.figure.toolbar_layout.addWidget(button)
    
        
        self.figure.axes = ax_orig, ax_rect = self.figure.figure.subplots(2, 1, gridspec_kw={"height_ratios": [2, 1]})#, sharex=True)
        ax_orig.set_xticks([])
        ax_rect.set_xlabel("Wavelength [Å]")
        ax_orig.set_ylabel("Flux")
        ax_rect.set_ylabel("Rectified flux")
        
        

        # hide the spines between ax and ax2
        ax_orig.spines.bottom.set_visible(False)
        ax_rect.spines.top.set_visible(False)
        ax_orig.xaxis.tick_top()
        ax_orig.tick_params(labeltop=False)  # don't put tick labels at the top
        ax_rect.xaxis.tick_bottom()

        d = .5  # proportion of vertical to horizontal extent of the slanted line
        color = '#ffffff'
        kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12,
                    linestyle="none", color=color, mec=color, mew=1, clip_on=False)
        ax_orig.plot([0, 1], [0, 0], transform=ax_orig.transAxes, **kwargs)
        ax_rect.plot([0, 1], [1, 1], transform=ax_rect.transAxes, **kwargs)        


        self.figure.figure.subplots_adjust(hspace=0.01, wspace=0.01)
        
        self.figure.canvas.draw()

        self.current_index = 0
        
        self._plot_data, = ax_orig.plot([], [], c="#ffffff")
        self._plot_model, = ax_orig.plot([], [], c="tab:red")
        self._plot_continuum, = ax_orig.plot([], [], c="tab:blue")
        
        self._plot_rectified_data, = ax_rect.plot([], [])
        self._plot_rectified_model, = ax_rect.plot([], [], c="tab:red")
        
        layout.addWidget(self.figure)    
        
        control_layout = QHBoxLayout()
        # use a grid layout
        grid_layout = QGridLayout()
        
        # create a combo box for the continuum type
        #Lengthgrid_layout.addWidget(QLabel(""))
            
        self.continuum_degree_label = QLabel("Degree:")
        self.continuum_degree = SpinBox(self)
        self.continuum_degree.setRange(0, 100)
        self.continuum_degree.setValue(10) # TODO: move to default options
        self.continuum_degree.setMaximumWidth(120)
                
        grid_layout.addWidget(QLabel("Fourier basis options"), 0, 0, 1, 4, QtCore.Qt.AlignLeft)
                
        grid_layout.addWidget(self.continuum_degree_label, 1, 0)
        grid_layout.addWidget(self.continuum_degree, 1, 1)
        
        self.continuum_length_scale = LineEdit()
        self.continuum_length_scale.setText("1000") # TODO: move to default options
        self.continuum_length_scale.setMaximumWidth(60)
        
        grid_layout.addWidget(QLabel("Absorption basis options"), 2, 0, 1, 4, QtCore.Qt.AlignLeft)
        grid_layout.addWidget(QLabel("Spectral resolution:"), 3, 0)
        self.Ro = LineEdit()
        self.Ro.setText("50,000")
        self.Ro.setMaximumWidth(100)
        grid_layout.addWidget(self.Ro, 3, 1)
        
        grid_layout.addWidget(QLabel("Stellar relative velocity:"), 4, 0)
        self.stellar_v_rel = LineEdit() 
        # set stellar_v_rel initially to zero
        self.stellar_v_rel.setText("0")        
        self.stellar_v_rel.setMaximumWidth(60)
        
        grid_layout.addWidget(self.stellar_v_rel, 4, 1)
        grid_layout.addWidget(QLabel("Telluric relative velocity:"), 4, 2)
        self.telluric_v_rel = LineEdit()
        self.telluric_v_rel.setText("0")
        self.telluric_v_rel.setMaximumWidth(60)
        
        
        #grid_layout.addWidget(self.action_fit_current_order, 0, 3)
        grid_layout.addWidget(self.telluric_v_rel, 4, 3)
        
        
        spacer = QLabel(" ")
        spacer.setMinimumWidth(30)        
        grid_layout.addWidget(spacer, 0, 4)
                
        control_layout.addLayout(grid_layout)
        layout.addLayout(control_layout)
        
        self.action_fit_current_order.clicked.connect(self.fit_current_order)
        self.action_fit_all_orders.clicked.connect(self.fit_all_orders) 
        
        self.figure.action_left.triggered.connect(lambda: self.update_canvas(self.current_index - 1))
        self.figure.action_right.triggered.connect(lambda: self.update_canvas(self.current_index + 1))
        
        self.update_canvas(self.current_index)        
        
        self.widget.setParent(self.card)
        self.widget.show()
        return None
    
    def fit_all_orders(self):
        print("OK")
        with h5.File("Grok/bosz-highres-optical-basis-vectors.h5") as fp:
            stellar_basis = NMFSpectralBasis(
                fp["vacuum_wavelength"][:], 
                fp["stellar_basis_vectors"][:].T,
                v_rel=float(self.stellar_v_rel.text()),
                Ro=float(self.Ro.text().replace(",", "")),
                Ri=85_000
            )
                        
        model = LinearStellarSpectrumModel(
            stellar_basis=stellar_basis,
            continuum_basis=FourierBasis(int(self.continuum_degree.text()))
        )

        print("fitting all orders")
        wavelength, flux, ivar, meta = self.session.get_spectral_order(self.current_index)
    
        self._fitted_result = λ, order_index, y, *_, continuum = model.fit(self.session.spectra)                
        print("done")
        print(model.θ[:32])
        self.update_canvas(self.current_index, update_axis_limits=False)
        return None
    
    def fit_current_order(self):

        
        with h5.File("Grok/bosz-highres-optical-basis-vectors.h5") as fp:
            stellar_basis = NMFSpectralBasis(
                fp["vacuum_wavelength"][:], 
                fp["stellar_basis_vectors"][:].T,
                v_rel=float(self.stellar_v_rel.text()),
                Ro=float(self.Ro.text().replace(",", "")),
                Ri=85_000
            )
                        
        model = LinearStellarSpectrumModel(
            stellar_basis=stellar_basis,
            continuum_basis=FourierBasis(int(self.continuum_degree.text()))
        )

        wavelength, flux, ivar, meta = self.session.get_spectral_order(self.current_index)
    
        λ, order_index, y, *_, continuum = model._fit(wavelength, flux, ivar)
        self._fitted_result = (λ, self.current_index * np.ones_like(order_index), y, *_, continuum)
            
        self.update_canvas(self.current_index, update_axis_limits=False)
        return None

                
    
    def update_canvas(self, index, update_axis_limits=True):
        S = self.session.n_orders
        
        if index < 0 or index > (S - 1):
            return
        wavelength, flux, ivar, meta = self.session.get_spectral_order(index)

            
        self._plot_data.set_data(wavelength, flux)
        try:
            λ, order_index, y, *_, continuum = self._fitted_result
        except AttributeError:
            for item in (self._plot_model, self._plot_rectified_model, self._plot_rectified_data, self._plot_continuum):
                item.set_data([], [])
        else:
            order_mask = (order_index == index)
            
            λm = λ[order_mask]
            cm = continuum[order_mask]
            ym = y[order_mask]
            self._plot_model.set_data(λm, ym)
            self._plot_continuum.set_data(λm, cm)
            
            try:
                self._plot_rectified_data.set_data(λm, flux / cm)
                self._plot_rectified_model.set_data(λm, ym / cm)
            except ValueError:
                for item in (self._plot_rectified_data, self._plot_rectified_model):
                    item.set_data([], [])
            
        self.current_index = index
        self.figure.action_left.setEnabled(self.current_index > 0)
        self.figure.action_right.setEnabled(self.current_index < (S - 1))
            
            
        #self._line_rect.set_data(wavelength, flux / continuum_)
        #self._line_H.set_data(wavelength, 1 - (W @ H).flatten())
        #self._line_continuum.set_data(wavelength, continuum_)
        
        if update_axis_limits:                
            for ax in self.figure.axes:
                ax.set_xlim(wavelength[[0, -1]])
                        
            ylims = (np.nanmin(flux), np.nanmax(flux))
            edge = 0.10 * np.ptp(ylims)
            ylims = [ylims[0] - edge, ylims[1] + edge]
            # limit to 0
            break_point = 1.2
            ylims[0] = max(ylims[0], break_point) # 
            self.figure.axes[1].set_ylim(0, break_point)
            self.figure.axes[0].set_ylim(ylims)
            
        self.figure.canvas.draw()
        #self.figure.canvas.setFocus()
        