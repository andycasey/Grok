
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
                
        def update_canvas(index):
            nonlocal current_index
            if index < 0 or index > (session.n_orders - 1):
                return
            
            λ, flux, *_ = session.get_spectral_order(index)
            
            line.set_data(λ, flux)
            plotter.ax.set_xlim(λ[[0, -1]])
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
            plotter.action_right.setEnabled(current_index < (session.n_orders - 1))
            plotter.canvas.setFocus()
            plotter.reset_current_as_home()
                        
        plotter.action_left.triggered.connect(lambda: update_canvas(current_index - 1))
        plotter.action_right.triggered.connect(lambda: update_canvas(current_index + 1))
        
        update_canvas(current_index)

        layout.addWidget(plotter)                
        
        self.widget.setParent(self.card)
        self.widget.show()
        
        return None
    
    def button_do_not_rectify_clicked(self):
        if self.callback is not None:
            self.callback()
        
    def button_rectify_clicked(self):
        if self.callback is not None:
            self.callback()
            
