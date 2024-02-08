

# coding:utf-8

# coding:utf-8
import time
from typing import List

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from astropy.io.fits import getval
from PyQt5 import QtCore
# import requried things
from PyQt5.QtCore import (QEvent, QItemSelectionModel, QModelIndex, QObject,
                          Qt, QUrl, QVariant, pyqtSignal, pyqtSlot, QThread)
from PyQt5.QtGui import (QColor, QDesktopServices, QDoubleValidator,
                         QIcon, QIntValidator, QPainter, QPalette, QPen)
from PyQt5.QtWidgets import (QApplication, QComboBox, QCompleter,
                             QFileDialog, QFrame, QGridLayout, QHBoxLayout,
                             QHeaderView, QItemDelegate, QLabel,
                             QListWidget, QListWidgetItem, QSizePolicy,
                             QStackedWidget, QStyle, QStyleOptionViewItem,
                             QTableView, QTableWidget, QTableWidgetItem,
                             QVBoxLayout, QWidget)
from qfluentwidgets import (Action, BodyLabel, BreadcrumbBar, CaptionLabel,
                            CardWidget, CheckableMenu, CheckBox, ComboBox,
                            CommandBar, EditableComboBox, FlowLayout, FluentIcon)
from qfluentwidgets import FluentIcon as FIF
from qfluentwidgets import (Flyout, IconWidget, InfoBarIcon, LineEdit,
                            MenuAnimationType, MenuIndicatorType,
                            MenuItemDelegate, Pivot, PrimaryDropDownPushButton,
                            PrimaryPushButton, PrimarySplitPushButton,
                            PushButton, RoundMenu, ScrollArea, SearchLineEdit,
                            SmoothScrollArea, SmoothScrollDelegate,
                            StrongBodyLabel, TabBar, TabCloseButtonDisplayMode,
                            TableItemDelegate, TableView, TableWidget,
                            TextWrap, Theme, TitleLabel, ToolButton,
                            ToolTipFilter, TransparentDropDownPushButton,
                            applyThemeColor, isDarkTheme, qrouter,
                            setCustomStyleSheet, setTheme, toggleTheme)

from ..common.config import cfg




import numpy as np
from Qt import QtCore
from Qt.QtGui import QDoubleValidator
from Qt.QtWidgets import (QHBoxLayout, QLabel, QVBoxLayout)
from qfluentwidgets import (ComboBox, LineEdit, PrimaryPushButton, PushButton)


from session import Session

from .analysis import AnalysisWidget
from .plot import BasePlotWidget, SinglePlotWidget
from astropy.table import Table


class TransitionsTableModel(QtCore.QAbstractTableModel):

    COLUMNS = [
        ("", "use"),
        ("λ", "wl"),
        ("Species", "species"),
        ("χ", "E_lower"),
        ("log(gf)", "log_gf"),
        ("Γ", "gamma_rad"),
        ("Deg", "poly_deg"),
        ("Tol", "wl_tol"),
        ("Profile", "profile"),
        ("EW", "equivalent_width"),
        ("A(X)", "abundance"),
    ]

    def __init__(self, data, **kwargs):
        super().__init__(**kwargs)
        
        # convert to a recarray 
        if data is None:
            self._row_count = 0
            self._data = np.recarray(
                [],
                dtype=[
                    ("use", bool),
                    ("wl", float),
                    ("species", str),
                    ("E_lower", float),
                    ("log_gf", float),
                    ("gamma_rad", float),
                    ("poly_deg", int),
                    ("wl_tol", float),
                    ("profile", str),
                    ("equivalent_width", float),
                    ("abundance", float),
                ]
            )          
        else:
            self._data = data
            self._row_count = len(data)
        '''
        self._row_count = len(linelist)
        if linelist:            
            rows = [
                (True, 1e8 * l.wl, f"{l.species}", l.E_lower, l.log_gf, l.gamma_rad, -1, 0, "Gaussian", np.nan, np.nan)
                for l in linelist
            ]
            self._data = np.rec.fromrecords(rows, names=[c for h, c in self.COLUMNS])
        else:
        '''
            
        return None
    

    def rowCount(self, parent):
        return self._row_count

    def columnCount(self, parent):
        return len(self.COLUMNS)

    def headerData(self, col, orientation, role):
        if orientation == QtCore.Qt.Horizontal \
        and role == QtCore.Qt.DisplayRole:
            return self.COLUMNS[col][0]
        return None


    def flags(self, index):
        #if not index.isValid(): return
        flags = QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable
        if index.column() == 0:
            flags |= QtCore.Qt.ItemIsUserCheckable 
        elif index.column() == 8:
            flags |= QtCore.Qt.ItemIsEditable
        return flags


    def data(self, index, role):
        if not index.isValid():
            return None
        
        if index.column() == 0:
            if role == QtCore.Qt.CheckStateRole:
                return int(self._data[self.COLUMNS[index.column()][1]][index.row()])
            else:
                return None

        elif role != QtCore.Qt.DisplayRole:
            return None

        return QVariant(f"{self._data[self.COLUMNS[index.column()][1]][index.row()]}")
    
    
    def setData(self, index, value, role=QtCore.Qt.DisplayRole):
        header, column = self.COLUMNS[index.column()]
        if column != "use":
            return False
        
        '''
        old_value = self.data(index, role)
        
        # Enable tri-state for first column, where the order is :
        # false, true, upper limit
        print(f"setting {index} {role} to {value} from {old_value}")
        if old_value == 2:
            value = 1
        elif old_value == 1:
            value = 0
                        
        self._data[index.row()][index.column()] = value
        self.dataChanged.emit(index, index)
        return value
        '''
        self._data["use"][index.row()] = value
        return value
                

    def sort(self, column, order):
        #self.emit(QtCore.SIGNAL("layoutAboutToBeChanged()"))
        _, column = self.COLUMNS[column]        
        try:
            self._data.sort(order=column) # sort in place                        
        except:            
            None
        else:
            if order == QtCore.Qt.DescendingOrder:
                self._data = self._data[::-1]
                            
        self.dataChanged.emit(
            self.createIndex(0, 0),
            self.createIndex(self.rowCount(0), self.columnCount(0))
        )
        #self.emit(QtCore.SIGNAL("layoutChanged()"))



class ComboDelegate(QItemDelegate):
    editorItems=['Gaussian', 'Voight','Lorentzian']
    #height = 25
    #width = 200
    def createEditor(self, parent, option, index):
        print(f"creating editor with parent {parent} {option} {index}")
        editor = TransparentDropDownPushButton(parent)

        self.menu = RoundMenu()
        self.menu.addAction(Action('Gaussian'))
        self.menu.addAction(Action('Voight'))
        self.menu.addAction(Action('Lorentzian'))
        
        #editor = ComboBox(parent)
        # editor.addItems(self.editorItems)
        # editor.setEditable(True)
        editor.setMenu(self.menu)
        #editor.currentIndexChanged.connect(self.currentIndexChanged)
        
        #editor.setCurrentIndex(-1)#self.editorItems.index(index.data()))
        #editor.addItems(self.editorItems)
        return editor

    '''
    def setEditorData(self,editor,index):
        z = 0
        for item in self.editorItems:
            #ai = QListWidgetItem(item)
            editor.addItem(item)
            if item == index.data():
                editor.setCurrentItem(editor.item(z))
            z += 1
        #editor.setGeometry(0,index.row()*self.height,self.width,self.height*len(self.editorItems))
    '''
    def setModelData(self, editor, model, index):
        print(self, editor, model, index)
        return None
        editorIndex=editor.currentIndex()
        text=self.editorItems[editorIndex]
        model.setData(index, text)
        # print '\t\t\t ...setModelData() 1', text

    @pyqtSlot()
    def currentIndexChanged(self): 
        self.commitData.emit(self.sender())
            


        
class ExcitationIonizationBalanceWidget(BasePlotWidget):    
    
    def __init__(
        self, 
        figsize=(30, 8),
        parent=None, 
        size_policy=(QSizePolicy.Expanding, QSizePolicy.Fixed),
        resize_interval=100
    ):
        super().__init__(toolbar=True, toolbar_left_right=False, toolbar_position="top", parent=parent, figsize=figsize, size_policy=size_policy, resize_interval=resize_interval)
        
        self.axes = self.canvas.figure.subplots(3, 1, gridspec_kw={"height_ratios": [1, 1, 1]})
        ax_profile, ax_excitation, ax_rew = self.axes
        ax_excitation.set_xlabel(r"$E_\mathrm{lower}$ [eV]")
        ax_excitation.set_ylabel(r"A(X)")
        
        ax_rew.set_xlabel(r"$\log_{10}(EW/\lambda)$")
        for ax in (ax_excitation, ax_rew):
            ax.set_ylabel(r"A(X)")
        
        ax_profile.set_xlabel(r"Wavelength $\lambda$ [$\mathrm{\AA}$]")
        ax_profile.set_ylabel(r"Rectified flux")
        self.figure.tight_layout()
        self.figure.canvas.draw()
        
        return None
    




class StellarParametersCOGWidget(AnalysisWidget):
    
    signalStatus = QtCore.Signal(str)
    
    def __init__(self, session=None, callback=None, parent=None):
        super().__init__(parent=parent, session=session, title="Stellar Parameters")
        self.parent = parent
        layout = QHBoxLayout(self.widget)
        layout.setContentsMargins(0, 0, 0, 0)    
        
        self.tableView = TableView(self)        

        self.table_model = TransitionsTableModel([])
        self.tableView.setModel(self.table_model)
        #self.tableView.scrollDelagate = SmoothScrollDelegate(self.tableView)
        self.tableView.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        
        #self.tableView.setItemDelegate(CustomTableItemDelegate(self.tableView))
        self.tableView.setItemDelegateForColumn(8, ComboDelegate(self))
        # self.tableView.setSelectRightClickedRow(True)

        # enable border
        self.tableView.setBorderVisible(True)
        self.tableView.setBorderRadius(8)

        self.tableView.setWordWrap(False)

        self.tableView.verticalHeader().hide()
        self.tableView.horizontalHeader().setSectionResizeMode(QHeaderView.Interactive)        
        self.tableView.horizontalHeader().setSectionResizeMode(0, QHeaderView.Fixed)
        self.tableView.setColumnWidth(0, 45)        
        self.tableView.resizeColumnsToContents()
        
        self.tableView.horizontalHeader().setSectionsMovable(True)
        self.tableView.setSortingEnabled(True)
        self.tableView.setDragEnabled(True) 
        self.tableView.setTabKeyNavigation(True)
        self.tableView.selectionModel().selectionChanged.connect(self.line_selection_changed)
        
        self.tableView.setMinimumSize(400, 300)
        self.tableView.setSelectionMode(QTableView.SingleSelection)


        button_atomic_line_list = PushButton("Select atomic line list")
        button_atomic_line_list.clicked.connect(self.select_atomic_line_list)
        transitions_button_layout = QHBoxLayout()
        transitions_button_layout.addWidget(button_atomic_line_list)
        transitions_button_layout.addStretch(1)
        
        left_layout = QVBoxLayout()
        left_layout.addLayout(transitions_button_layout)
        left_layout.addWidget(self.tableView)

        
        self.line_edit_teff = LineEdit(self)
        self.line_edit_teff.setValidator(QIntValidator(2_500, 10_000, self))
        self.line_edit_teff.setText(f"{cfg.get(cfg.initial_teff):,}")
        self.line_edit_teff.setFixedWidth(80)
        self.line_edit_teff.setAlignment(Qt.AlignHCenter)

        self.line_edit_logg = LineEdit(self)
        self.line_edit_logg.setValidator(QDoubleValidator(0, 5, 3, self))
        self.line_edit_logg.setText("4.5")
        self.line_edit_logg.setFixedWidth(80)
        self.line_edit_logg.setAlignment(Qt.AlignHCenter)        
        
        self.line_edit_feh = LineEdit(self)
        self.line_edit_feh.setValidator(QDoubleValidator(-5, 1, 3, self))
        self.line_edit_feh.setText("0")
        self.line_edit_feh.setFixedWidth(80)
        self.line_edit_feh.setAlignment(Qt.AlignHCenter)        

        self.line_edit_v_micro = LineEdit(self)
        self.line_edit_v_micro.setValidator(QDoubleValidator(0, 10, 3, self))
        self.line_edit_v_micro.setText("1")
        self.line_edit_v_micro.setFixedWidth(80)
        self.line_edit_v_micro.setAlignment(Qt.AlignHCenter)
                
        param_grid_layout = QGridLayout()
        param_grid_layout.setColumnStretch(0, 1)
        param_grid_layout.addWidget(QLabel("Effective temperature:"), 0, 0)
        param_grid_layout.addWidget(self.line_edit_teff, 0, 2)
        param_grid_layout.addWidget(QLabel("K"), 0, 3)
        
        param_grid_layout.addWidget(QLabel("Surface gravity:"), 1, 0)
        param_grid_layout.addWidget(self.line_edit_logg, 1, 2)
        param_grid_layout.addWidget(QLabel("dex"), 1, 3)
        
        param_grid_layout.addWidget(QLabel("Metallicity ([M/H]):"), 2, 0)
        param_grid_layout.addWidget(self.line_edit_feh, 2, 2)
        param_grid_layout.addWidget(QLabel("dex"), 2, 3)
        
        param_grid_layout.addWidget(QLabel("Microturbulence:"), 3, 0)
        param_grid_layout.addWidget(self.line_edit_v_micro, 3, 2)
        param_grid_layout.addWidget(QLabel("km/s"), 3, 3)
        
        
        param_button_layout = QHBoxLayout()
        param_button_layout.setContentsMargins(0, 10, 0, 10)
        param_button_layout.addStretch(1)
        
        push_button_solve = PrimarySplitPushButton("Solve stellar parameters")
                
        self.actions_hold_stellar_parameter_constant = [
            Action("Hold effective temperature constant"),
            Action("Hold surface gravity constant"),
            Action("Hold metallicity constant"),
            Action("Hold microturbulence constant"),
        ]
        for action in self.actions_hold_stellar_parameter_constant:
            action.setCheckable(True)
        
        menu = CheckableMenu(parent=self, indicatorType=MenuIndicatorType.CHECK)
        menu.addActions(self.actions_hold_stellar_parameter_constant)
                
        push_button_solve.setFlyout(menu)
        push_button_solve.clicked.connect(self.solve_stellar_parameters)
        
        push_button_compute = PushButton("Compute line abundances")
        push_button_compute.clicked.connect(self.compute_line_abundances)
        
        param_button_layout.addWidget(push_button_compute)
        param_button_layout.addWidget(push_button_solve)
                
        left_layout.addLayout(param_grid_layout)
        left_layout.addLayout(param_button_layout)
        left_layout.addStretch(1)
        
        layout.addLayout(left_layout)
        plot_layout = QVBoxLayout()
        plot_layout.setContentsMargins(0, 30, 0, 0)
        self.excitation_ionization_plot = ExcitationIonizationBalanceWidget(parent=self)
        
        self._excitation_scatter = self.excitation_ionization_plot.axes[1].scatter(
            [], [],
            edgecolor="k",
            lw=0.5,
        )
        self._excitation_mean = self.excitation_ionization_plot.axes[1].plot(
            [-100, 100],
            [np.nan, np.nan],
            lw=0.5, ls=":"
        )
        self._rew_scatter = self.excitation_ionization_plot.axes[2].scatter(
            [], [], edgecolor="k", lw=0.5, c="#FFFFFF", 
        )
        self._rew_mean = self.excitation_ionization_plot.axes[2].plot(
            [-100, 100],
            [np.nan, np.nan],
            lw=0.5, ls=":", c="#FFFFFF", # TODO: don't set color here, use from theme
        )        
        
        plot_layout.addWidget(self.excitation_ionization_plot)
        layout.addLayout(plot_layout)
        
        
        self.bottomLayout.addWidget(PushButton("Compute abundances with these stellar parameters"))
        self.bottomLayout.addStretch(1)

                
        self.widget.setParent(self.card)
        self.widget.show()
        
        return None
    
    
    def solve_stellar_parameters(self):
        
        Teff0 = float(self.line_edit_teff.text().replace(",", ""))
        logg0 = float(self.line_edit_logg.text())
        metallicity0 = float(self.line_edit_feh.text())
        vmic0 = float(self.line_edit_v_micro.text())
        
        
        self.session.cog_read_linelist(
            "/Users/andycasey/research/Grok/python/Grok/Melendez_et_al_Fe.moog.sorted", 
            format="moog",
        )
        
        self.session.cog_interpolate_atmosphere(Teff0, logg0, metallicity0, metallicity0)
        
        ews = list(np.loadtxt("/Users/andycasey/research/Grok/python/Grok/Melendez_et_al_Fe.moog.sorted", usecols=(5, ), skiprows=1))
    
        E_lower = [l.E_lower for l in self.session.cog_linelist]
        REW = [np.log10(ew/(l.wl * 1e11)) for ew, l in zip(ews, self.session.cog_linelist)]
            
        # set colors of points based on species
        self._excitation_scatter.set_offsets(np.array([E_lower, np.nan * np.ones_like(E_lower)]).T)
        self._rew_scatter.set_offsets(np.array([REW, np.nan * np.ones_like(E_lower)]).T) 
        
        vmin, vmax = (0, 1)
        n = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
        cmap = mpl.colors.LinearSegmentedColormap.from_list(
            "",
            ["tab:blue", "tab:orange"]#, "tab:red", "tab:green"]
        )
        m = mpl.cm.ScalarMappable(norm=n, cmap=cmap)
        
        colors = m.to_rgba([l.species.charge for l in self.session.cog_linelist])
        for scat in (self._excitation_scatter, self._rew_scatter):
            scat.set_facecolor(colors)
            scat.set_clim(vmin=vmin, vmax=vmax)
        
        self.excitation_ionization_plot.axes[1].set_xlim(get_limits(E_lower))
        self.excitation_ionization_plot.axes[2].set_xlim(get_limits(REW))
        self.excitation_ionization_plot.canvas.draw()

        self._ylim_ptp = None
                
        
        kwds = dict(callback=self.update_excitation_ionization_axes)
        
        # Move this logic to the session call.
        
        max_step_sizes = [1000.0, 1.0, 0.3, 0.5]
        any_held_constant = False
        for i, action in enumerate(self.actions_hold_stellar_parameter_constant):
            if action.isChecked():
                max_step_sizes[i] = 0
                any_held_constant = True
                
        if any_held_constant:
            kwds["max_step_sizes"] = max_step_sizes
        
        foo = self.session.cog_solve_stellar_parameters(
            ews,
            Teff0=Teff0,
            logg0=logg0,
            vmic0=vmic0,
            metallicity0=metallicity0,
            **kwds
        )
        print(f"foo: {foo}")
        
        
    def update_stellar_parameter_line_edits(self, teff, logg, vmic, m_h):
        self.line_edit_teff.setText(f"{teff:,.0f}")
        self.line_edit_logg.setText(f"{logg:.3f}")
        self.line_edit_v_micro.setText(f"{vmic:.3f}")
        self.line_edit_feh.setText(f"{m_h:.3f}")
                
    
        
    def update_excitation_ionization_axes(self, params, residuals, line_abundances):
        print(params, residuals)#, line_abundances)
        
        teff, logg, vmic, m_h = params
        self.line_edit_teff.setText(f"{teff:,.0f}")
        self.line_edit_logg.setText(f"{logg:.3f}")
        self.line_edit_v_micro.setText(f"{vmic:.3f}")
        self.line_edit_feh.setText(f"{m_h:.3f}")
        
        ylim = get_limits(line_abundances)
        try:
            self._ylim_ptp
        except:
            self._ylim_ptp = ylim
        else:

            if self._ylim_ptp is not None:
                ylim_mean = np.mean(ylim)
                edge = 0.5 * np.ptp(self._ylim_ptp)
                ylim = (ylim_mean - edge, ylim_mean + edge)
            
        for ax in self.excitation_ionization_plot.axes[1:]:
            ax.set_ylim(ylim)     
        
        mean = np.nanmean(line_abundances)
        self._excitation_mean[0].set_ydata([mean, mean])
        self._rew_mean[0].set_ydata([mean, mean])

        '''
        # Crude animation (should do this when QKorgProcess is in use)
        if np.any(np.isfinite(self.line_abundances)):
            print("STEPPING")
            step = 10
            dx = (np.array(line_abundances) - np.array(self.line_abundances)) / step
            for i in range(step):
                y = self.line_abundances + i * dx
                self._excitation_scatter.set_offsets(np.array([E_lower, y]).T)
                self._rew_scatter.set_offsets(np.array([REW, y]).T)
                self.excitation_ionization_plot.canvas.draw()
                QApplication.processEvents()
            from time import sleep
            sleep(3)                    
        else:
        ''' 
        x = self._excitation_scatter.get_offsets().T[0]
        self._excitation_scatter.set_offsets(np.array([x, line_abundances]).T)
        x = self._rew_scatter.get_offsets().T[0]
        self._rew_scatter.set_offsets(np.array([x, line_abundances]).T)
        #(np.array([REW, line_abundances]).T)
        self.excitation_ionization_plot.canvas.draw()
        QApplication.processEvents()                    

    
    def get_current_stellar_parameters(self):    
        Teff = float(self.line_edit_teff.text().replace(",", ""))
        logg = float(self.line_edit_logg.text())
        metallicity = float(self.line_edit_feh.text())
        vmic = float(self.line_edit_v_micro.text())   
        return (Teff, logg, metallicity, vmic)

    
    def compute_line_abundances(self):
    
        Teff0, logg0, metallicity0, vmic0 = self.get_current_stellar_parameters()
        
        self.session.cog_interpolate_atmosphere(Teff0, logg0, metallicity0, metallicity0)
        
        ews = list(np.loadtxt("/Users/andycasey/research/Grok/python/Grok/Melendez_et_al_Fe.moog.sorted", usecols=(5, ), skiprows=1))
        abundances = self.session.cog_ews_to_abundances(ews, vmic=vmic0)
        
        self.update_excitation_ionization_axes((Teff0, logg0, vmic0, metallicity0), None, abundances)
            
        # update figure
        
    
    
    def resizeEvent(self, e):
        super().resizeEvent(e)
    
    def select_atomic_line_list(self):
        # TODO: We should handle the format better by using a custom dialog, but this is a quick fix.
        
        filenames, selected_filter = QFileDialog.getOpenFileNames(
            self, 
            caption="Select atomic line list", 
            directory="", 
            # TODO: get filter order from app config
            filter="MOOG (*);;VALD (*);;Kurucz (*);;TurboSpectrum (*)"
        )
        if filenames:            
            linelist = self.session.cog_read_linelist(
                filenames[0], 
                format=selected_filter.split("(")[0].lower().strip()
            )
            
            # Create recarray view for the table view 
            rows = []
            for line in linelist:
                rows.append((1, line.wl * 1e8, f"{line.species}", line.E_lower, line.log_gf, line.gamma_rad, -1, 0, "Gaussian", np.nan, np.nan))

            data = np.rec.fromrecords(
                rows,
                names=[c for h, c in TransitionsTableModel.COLUMNS]
            )
            self.table_model = TransitionsTableModel(data)
            self.tableView.setModel(self.table_model)
            self.tableView.resizeColumnsToContents()

            
    def line_current_changed(self):
        """Force the first column (checkbox) to always be in focus."""
        model = self.tableView.selectionModel()
        current_index = model.currentIndex()
        if current_index.column() > 0:
            model.setCurrentIndex(
                self.table_model.index(current_index.row(), 0), 
                QItemSelectionModel.Select
            )
        return None
    
    def line_selection_changed(self):        
        try:
            index = self.tableView.selectionModel().selectedRows()[0]
        except IndexError:
            return None
        else:
            print(f"selection changed {index}")
        
        

def get_limits(data, percent=0.05):
    min_v, max_v = (np.nanmin(data), np.nanmax(data))
    ptp = (max_v - min_v)
    offset = percent * ptp
    return (min_v - offset, max_v + offset)        
