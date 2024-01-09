# coding:utf-8
import time
from typing import List

import matplotlib.pyplot as plt
import numpy as np
from astropy.io.fits import getval
from PyQt5 import QtCore
# import requried things
from PyQt5.QtCore import (QEvent, QItemSelectionModel, QModelIndex, QObject,
                          Qt, QUrl, QVariant, pyqtSignal, pyqtSlot)
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
from session import Session

from ..common.config import (ACRONYM, EXAMPLE_URL, FEEDBACK_URL,
                             FILENAME_SUFFIX, HELP_URL, cfg)
from ..common.icon import Icon
from ..common.style_sheet import StyleSheet
from ..components.plot import (ExcitationIonizationBalanceWidget,
                               SinglePlotWidget)
from ..components.tool_bar import ToolBar
from ..components.analysis import AnalysisWidget
from ..components.quick_view import QuickViewWidget
from ..view.gallery_interface import SeparatorWidget


class ExampleCard2(QWidget):
    """ Example card """

    def __init__(self, title, widget: QWidget, leftButtonText=None, rightButtonText=None, stretch=0, parent=None):
        super().__init__(parent=parent)
        self.widget = widget
        self.stretch = stretch

        self.titleLabel = StrongBodyLabel(title, self)
        self.card = QFrame(self)

        self.sourceWidget = QFrame(self.card)
        
        self.leftButtonText = leftButtonText
        self.rightButtonText = rightButtonText
        #self.sourcePathLabel = BodyLabel(
        #    self.tr('Source code'), 
        #    self.sourceWidget
        #)
        #self.linkIcon = IconWidget(FluentIcon.LINK, self.sourceWidget)

        self.vBoxLayout = QVBoxLayout(self)
        self.cardLayout = QVBoxLayout(self.card)
        self.topLayout = QHBoxLayout()
        self.bottomLayout = QHBoxLayout(self.sourceWidget)

        self.__initWidget()

    def __initWidget(self):
        #self.linkIcon.setFixedSize(16, 16)
        self.__initLayout()

        #self.sourceWidget.setCursor(Qt.PointingHandCursor)
        self.sourceWidget.installEventFilter(self)

        self.card.setObjectName('card')
        self.sourceWidget.setObjectName('sourceWidget')

    def __initLayout(self):
        self.vBoxLayout.setSizeConstraint(QVBoxLayout.SetMinimumSize)
        self.cardLayout.setSizeConstraint(QVBoxLayout.SetMinimumSize)
        self.topLayout.setSizeConstraint(QHBoxLayout.SetMinimumSize)

        self.vBoxLayout.setSpacing(12)
        self.vBoxLayout.setContentsMargins(0, 0, 0, 0)
        self.topLayout.setContentsMargins(12, 12, 12, 12)
        self.bottomLayout.setContentsMargins(18, 18, 18, 18)
        self.cardLayout.setContentsMargins(0, 0, 0, 0)

        self.vBoxLayout.addWidget(self.titleLabel, 0, Qt.AlignTop)
        self.vBoxLayout.addWidget(self.card, 0, Qt.AlignTop)
        self.vBoxLayout.setAlignment(Qt.AlignTop)

        self.cardLayout.setSpacing(0)
        self.cardLayout.setAlignment(Qt.AlignTop)
        self.cardLayout.addLayout(self.topLayout, 0)
        self.cardLayout.addWidget(self.sourceWidget, 0, Qt.AlignBottom)

        self.widget.setParent(self.card)
        self.topLayout.addWidget(self.widget)
        if self.stretch == 0:
            self.topLayout.addStretch(1)

        self.widget.show()

        #self.bottomLayout.addWidget(self.sourcePathLabel, 0, Qt.AlignLeft)
        #self.bottomLayout.addStretch(1)
        #self.bottomLayout.addWidget(self.linkIcon, 0, Qt.AlignRight)
        
        if self.leftButtonText is not None:
            leftButton = PushButton(self.leftButtonText)
            #leftButton.clicked.connect(self.enable_norm)
            self.bottomLayout.addWidget(leftButton, 0, Qt.AlignLeft)
            
        self.bottomLayout.addStretch(1)
        if self.rightButtonText is not None:
            rightButton = PrimaryPushButton(self.rightButtonText, self)
            #rightButton.clicked.connect(self.enable_norm)
            
            self.bottomLayout.addWidget(rightButton, 0, Qt.AlignRight)

        self.bottomLayout.setAlignment(Qt.AlignLeft | Qt.AlignVCenter)

    #def enable_norm(self):
    #    self.parent.vBoxLayout.addWidget(self.parent.continuum_card)
                

    #def eventFilter(self, obj, e):
    #    if obj is self.sourceWidget:
    #        if e.type() == QEvent.MouseButtonRelease:
    #            QDesktopServices.openUrl(QUrl(self.sourcePath))

    #        return super().eventFilter(obj, e)
    
#from style_sheet import StyleSheet


class KeyPressFilter(QObject):

    def eventFilter(self, widget, event):
        try:
            print(f"ok: {event.type()} {event.key()} {event.text()}")
        except:
            None
        '''
        if event.type() == QEvent.KeyPress:
            text = event.text()
            if event.modifiers():
                text = event.keyCombination().key().name.decode(encoding="utf-8")
            print(f"{event} {event.type}: {text}")
        '''            
        return False
    



class TopAndBottomAnalysisWidget(QWidget):
    
    def __init__(self, session, title=None, parent=None, stretch=0):
        super().__init__(parent=parent)
        self.session = session
        self.stretch = stretch
        self.widget = QWidget(self)

        self.titleLabel = StrongBodyLabel(title, self)
        self.card = QFrame(self)

        self.sourceWidget = QFrame(self.card)
        
        self.vBoxLayout = QVBoxLayout(self)
        self.cardLayout = QVBoxLayout(self.card)
        self.topLayout = QHBoxLayout()
        self.bottomLayout = QHBoxLayout(self.sourceWidget)

        self.__initWidget()    
        
    def __initWidget(self):
        #self.linkIcon.setFixedSize(16, 16)
        self.__initLayout()

        #self.sourceWidget.setCursor(Qt.PointingHandCursor)
        self.sourceWidget.installEventFilter(self)

        self.card.setObjectName('card')
        self.sourceWidget.setObjectName('sourceWidget')

    def __initLayout(self):
        self.vBoxLayout.setSizeConstraint(QVBoxLayout.SetMinimumSize)
        self.cardLayout.setSizeConstraint(QVBoxLayout.SetMinimumSize)
        self.topLayout.setSizeConstraint(QHBoxLayout.SetMinimumSize)

        self.vBoxLayout.setSpacing(12)
        self.vBoxLayout.setContentsMargins(0, 0, 0, 0)
        self.topLayout.setContentsMargins(12, 12, 12, 12)
        self.bottomLayout.setContentsMargins(18, 18, 18, 18)
        self.cardLayout.setContentsMargins(0, 0, 0, 0)

        self.vBoxLayout.addWidget(self.titleLabel, 0, Qt.AlignTop)
        self.vBoxLayout.addWidget(self.card, 0, Qt.AlignTop)
        self.vBoxLayout.setAlignment(Qt.AlignTop)

        self.cardLayout.setSpacing(0)
        self.cardLayout.setAlignment(Qt.AlignTop)
        self.cardLayout.addLayout(self.topLayout, 0)
        self.cardLayout.addWidget(self.sourceWidget, 0, Qt.AlignBottom)

        self.topLayout.addWidget(self.widget)
        if self.stretch == 0:
            self.topLayout.addStretch(1)


        '''
        #self.bottomLayout.addWidget(self.sourcePathLabel, 0, Qt.AlignLeft)
        #self.bottomLayout.addStretch(1)
        #self.bottomLayout.addWidget(self.linkIcon, 0, Qt.AlignRight)
        
        if self.leftButtonText is not None:
            leftButton = PushButton(self.leftButtonText)
            #leftButton.clicked.connect(self.enable_norm)
            self.bottomLayout.addWidget(leftButton, 0, Qt.AlignLeft)
            
        self.bottomLayout.addStretch(1)
        if self.rightButtonText is not None:
            rightButton = PrimaryPushButton(self.rightButtonText, self)
            #rightButton.clicked.connect(self.enable_norm)
            
            self.bottomLayout.addWidget(rightButton, 0, Qt.AlignRight)

        '''
        self.bottomLayout.setAlignment(Qt.AlignLeft | Qt.AlignVCenter)
        
        self.widget.setParent(self.card)
        self.widget.show()
        
        
    def eventFilter(self, widget, event):
        try:
            print(f"ok: {event.type()} {event.key()} {event.text()}")
        except:
            None
        '''
        if event.type() == QEvent.KeyPress:
            text = event.text()
            if event.modifiers():
                text = event.keyCombination().key().name.decode(encoding="utf-8")
            print(f"{event} {event.type}: {text}")
        '''            
        return False             
    





class CustomTableItemDelegate(TableItemDelegate):
    """ Custom table item delegate """

    def initStyleOption(self, option: QStyleOptionViewItem, index: QModelIndex):
        super().initStyleOption(option, index)
        if index.column() != 1:
            return

        if isDarkTheme():
            option.palette.setColor(QPalette.Text, Qt.white)
            option.palette.setColor(QPalette.HighlightedText, Qt.white)
        else:
            option.palette.setColor(QPalette.Text, Qt.red)
            option.palette.setColor(QPalette.HighlightedText, Qt.red)    

class TransitionsTableModel(QtCore.QAbstractTableModel):

    header = ["", "λ", "Species", "χ", "log(gf)", "Γ", "Deg", "Tol", "Profile", "EW", "A(X)"]
    attrs = (
        "is_acceptable",
        "_repr_wavelength", 
        "_repr_element", 
        "chi",
        "loggf",
        "C6",
        "poly",
        "wl_tol",
        "option",
        "equivalent_width",
        "A(X)",
    )

    def __init__(self, parent, data, *args):
        super().__init__(parent, *args)
        self._data = data

    def rowCount(self, parent):
        return len(self._data)

    def columnCount(self, parent):
        return len(self.header)

    def data(self, index, role):
        if not index.isValid():
            return None

        value = self._data[index.row()][index.column()]

        if index.column() == 0:
            if role == QtCore.Qt.CheckStateRole:
                return value #QtCore.Qt.Checked if value else QtCore.Qt.Unchecked
            else:
                return None

        elif role != QtCore.Qt.DisplayRole:
            return None

        return QVariant(value)
    

    def headerData(self, col, orientation, role):
        if orientation == QtCore.Qt.Horizontal \
        and role == QtCore.Qt.DisplayRole:
            return self.header[col]
        return None


    def setData(self, index, value, role=QtCore.Qt.DisplayRole):
        print(f"setData {index} {value}")
        attr = self.attrs[index.column()]
        if attr != "is_acceptable":
            return False

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

    def sort(self, column, order):
        #self.emit(QtCore.SIGNAL("layoutAboutToBeChanged()"))
        print(column, order)
        self._data = sorted(self._data, key=lambda sm: sm[column])
        
        if order == QtCore.Qt.DescendingOrder:
            self._data.reverse()

        self.dataChanged.emit(self.createIndex(0, 0),
            self.createIndex(self.rowCount(0), self.columnCount(0)))
        #self.emit(QtCore.SIGNAL("layoutChanged()"))

    def flags(self, index):
        #if not index.isValid(): return
        flags = QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable
        if index.column() == 0:
            flags |= QtCore.Qt.ItemIsUserCheckable 
        elif index.column() == 8:
            flags |= QtCore.Qt.ItemIsEditable
        return flags


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
            

class StellarParametersWidget(AnalysisWidget):
    
    def __init__(self, session, callback=None, parent=None):
        super().__init__(parent=parent, session=session, title="Stellar Parameters")
        layout = QVBoxLayout(self.widget)
        layout.setContentsMargins(0, 0, 0, 0)
        
        self.primaryButton = PrimaryPushButton("Measure all equivalent widths")
        top_layout = QHBoxLayout()
        top_layout.addStretch(1)
        top_layout.addWidget(self.primaryButton)
        layout.addLayout(top_layout)        
        
        lr_layout = QHBoxLayout()
        left_layout = QVBoxLayout()
        # stellar params
        left_layout.setContentsMargins(10, 10, 10, 10)
        left_layout.setSpacing(10)
        
        line_edit_teff = LineEdit(self)
        line_edit_teff.setValidator(QIntValidator(2_500, 10_000, self))
        line_edit_teff.setText(f"{cfg.get(cfg.initial_teff):,}")
        line_edit_teff.setFixedWidth(80)
        line_edit_teff.setAlignment(Qt.AlignHCenter)

        line_edit_logg = LineEdit(self)
        line_edit_logg.setValidator(QDoubleValidator(0, 5, 3, self))
        line_edit_logg.setText("4.5")
        line_edit_logg.setFixedWidth(80)
        line_edit_logg.setAlignment(Qt.AlignHCenter)        
        
        line_edit_feh = LineEdit(self)
        line_edit_feh.setValidator(QDoubleValidator(-5, 1, 3, self))
        line_edit_feh.setText("0")
        line_edit_feh.setFixedWidth(80)
        line_edit_feh.setAlignment(Qt.AlignHCenter)        

        line_edit_v_micro = LineEdit(self)
        line_edit_v_micro.setValidator(QDoubleValidator(0, 10, 3, self))
        line_edit_v_micro.setText("1")
        line_edit_v_micro.setFixedWidth(80)
        line_edit_v_micro.setAlignment(Qt.AlignHCenter)
                
        param_layout = QGridLayout()
        param_layout.setColumnStretch(0, 1)
        param_layout.addWidget(QLabel("Effective temperature:"), 0, 0)
        param_layout.addWidget(line_edit_teff, 0, 2)
        param_layout.addWidget(QLabel("K"), 0, 3)
        
        param_layout.addWidget(QLabel("Surface gravity:"), 1, 0)
        param_layout.addWidget(line_edit_logg, 1, 2)
        param_layout.addWidget(QLabel("dex"), 1, 3)
        
        param_layout.addWidget(QLabel("Metallicity ([M/H]):"), 2, 0)
        param_layout.addWidget(line_edit_feh, 2, 2)
        param_layout.addWidget(QLabel("dex"), 2, 3)
        
        param_layout.addWidget(QLabel("Microturbulence:"), 3, 0)
        param_layout.addWidget(line_edit_v_micro, 3, 2)
        param_layout.addWidget(QLabel("km/s"), 3, 3)
        
        left_layout.addLayout(param_layout)
        
        
        
        #right_layout = QVBoxLayout()
        #right_layout.addWidget(ExcitationIonizationBalanceWidget(parent=self), 0, Qt.AlignTop)
        
        
        self.tableView = TableView(self)        
        #header = ["", "Wavelength", "Species", "χ", "log(gf)", "Γ", "Deg", "Tol", "Profile", "Automask", "EW", "σ(EW)", "A(X)", "σ(A(X))"]
        self.table_model = TransitionsTableModel(self.tableView, [
            [0, 5160, "Fe I", 4.22, -1.01, 0.003, 3, 0.1, "Gaussian", "0.1 ± 0.1", "7.5 ± 0.1"],
            [1, 5172, "Fe II", 4.1, -3, None, -1, 0.1, "Voight", "0.1 ± 0.1", "7.5 ± 0.1"],         
        ])
        self.tableView.setModel(self.table_model)
        self.tableView.scrollDelagate = SmoothScrollDelegate(self.tableView)
        self.tableView.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        
        
        #    attrs = ("is_acceptable", "_repr_wavelength", "_repr_element", 
        #    "equivalent_width")
        
        #self.tableView.setFixedSize(400, 800)

        # NOTE: use custom item delegate
        #self.tableView.setItemDelegate(CustomTableItemDelegate(self.tableView))
        self.tableView.setItemDelegateForColumn(8, ComboDelegate(self))



        # select row on right-click
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
        #self.tableView.selectionModel().currentChanged.connect(self.line_current_changed)
        self.tableView.selectionModel().selectionChanged.connect(self.line_selection_changed)
        
        #        self.tableView.selectionModel().
        self.tableView.setMinimumSize(800, 200)
        self.tableView.setSelectionMode(QTableView.SingleSelection)
        #self.tableView.setSelectionBehavior(..)

        self.setStyleSheet("Demo{background: rgb(255, 255, 255)} ")
        #layout.setContentsMargins(50, 30, 50, 30)
        lr_layout.addLayout(left_layout)
        lr_layout.addStretch(1)
        #lr_layout.addLayout(right_layout)
        layout.addWidget(self.tableView)
        layout.addWidget(ExcitationIonizationBalanceWidget(parent=self))
        layout.addLayout(lr_layout)
                
        
        #layout.addWidget(QLabel(" " g gglshd))  
        

        #self.primaryButton = PrimaryPushButton("Measure all equivalent widths")
        #self.bottomLayout.addStretch(1)
        #self.bottomLayout.addWidget(self.primaryButton)
        self.bottomLayout.addWidget(PushButton("Compute abundances with these stellar parameters"))
        self.bottomLayout.addStretch(1)
        push_button_solve = PrimarySplitPushButton("Solve stellar parameters")
        
        #menu = RoundMenu()
        menu = CheckableMenu(parent=self, indicatorType=MenuIndicatorType.CHECK)
        menu.addAction(Action("Hold effective temperature constant"))
        menu.addAction(Action("Hold surface gravity constant"))
        menu.addAction(Action("Hold metallicity constant"))
        menu.addAction(Action("Hold microturbulence constant"))
        for i in range(3):            
            menu.actions()[i].setCheckable(True)
                
        push_button_solve.setFlyout(menu)
        
        self.bottomLayout.addWidget(push_button_solve)
        
        self.widget.setParent(self.card)
        self.widget.show()
        
        return None
    
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
        


class FloatLineEdit(LineEdit):

    valueChanged = pyqtSignal(str)

    def __init__(self, value, lower, upper, decimals, parent=None):
        super().__init__(parent)
        self.setText(str(value))
        #self.setFixedSize(136, 33)
        self.setClearButtonEnabled(True)
        self.setValidator(QDoubleValidator(lower, upper, decimals, self))

        self.textEdited.connect(self._onTextEdited)

    def _onTextEdited(self, text):
        """ text edited slot """
        
        state = self.validator().validate(text, 0)[0]
        if state == QDoubleValidator.Acceptable:
            self.valueChanged.emit(text)
            print(f"ok {text}")
        else:
            print("nope")
            
        return False

    

        


class SessionInterface(ScrollArea):
    """ Dialog interface """

    def __init__(self, session, parent=None):
        super().__init__(parent=parent)
        
        self.session = session
        
        self.view = QWidget(self)
        self.toolBar = ToolBar(session, parent=self)
        self.vBoxLayout = QVBoxLayout(self.view)

        self.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.setViewportMargins(0, self.toolBar.height(), 0, 0)
        self.setWidget(self.view)
        self.setWidgetResizable(True)

        self.vBoxLayout.setSpacing(30)
        self.vBoxLayout.setAlignment(Qt.AlignTop)
        self.vBoxLayout.setContentsMargins(36, 20, 36, 36)

        self.view.setObjectName('view')
        StyleSheet.ANALYSIS_INTERFACE.apply(self)
        
        self.setObjectName('fooInterface')

        initial_visibility = False

        def continuum_callback():
            self.cog_eqw = StellarParametersWidget(session, self)
            # The hide/show makes the widget appear less 'jumpy'
            self.cog_eqw.setVisible(False)
            self.vBoxLayout.addWidget(self.cog_eqw, 0, Qt.AlignTop)
            self.cog_eqw.setVisible(True)
            

        # Radial velocity stuff
        #self.continuum = ContinuumRectificationWidget(session, callback=continuum_callback, parent=self)
        #self.continuum.setVisible(initial_visibility)
        
        
        def rv_callback():
            self.continuum.setVisible(True)
        
        qv = QuickViewWidget(session, parent=self)
        self.vBoxLayout.addWidget(qv)
        
        #card = RadialVelocityWidget(session, callback=rv_callback, parent=self)
        #self.vBoxLayout.addWidget(card, 0, Qt.AlignTop)
        
        # Continuum normalization

        #self.vBoxLayout.addWidget(self.continuum, 0, Qt.AlignTop)
        
        
        # add stellar parameter analysis? --> differential, etc.

    def add_analysis_widget(self, widget_class):#, from_action=None):
        widget = widget_class(self.session)
        
        # not sure that this is the right thing to do here, but without the
        # visible False/True and the processEvents, it doesn't scroll to the new widget
        widget.setVisible(False)
        self.vBoxLayout.addWidget(widget, 0, Qt.AlignTop)
        widget.setVisible(True)
                
        QApplication.processEvents()
        
        w = self.vBoxLayout.itemAt(self.vBoxLayout.count() - 1).widget()        
        self.verticalScrollBar().setValue(w.y())
        #if from_action is not None:
        #    from_action.setIcon(FluentIcon.DATE_TIME.qicon())
        

    def resizeEvent(self, e):
        super().resizeEvent(e)
        self.toolBar.resize(self.width(), self.toolBar.height())
        
        


class SessionTabsInterface(QWidget):
    """ Session tabs interface """

    def __init__(self, parent=None):
        super().__init__(parent=parent)
        self.parent = parent
        self.tabCount = 0

        self.tabBar = TabBar(self)
        self.tabBar.setMovable(True)
        self.stackedWidget = QStackedWidget(self)
        
        self.tabView = QWidget(self)

        self.hBoxLayout = QHBoxLayout(self)
        self.vBoxLayout = QVBoxLayout(self.tabView)

        # add items to pivot
        self.tabBar.setTabMaximumWidth(200)

        self.hBoxLayout.addWidget(self.tabView, 1)
        self.hBoxLayout.setContentsMargins(0, 0, 0, 0)

        self.vBoxLayout.addWidget(self.tabBar)
        self.vBoxLayout.addWidget(self.stackedWidget)
        self.vBoxLayout.setContentsMargins(0, 0, 0, 0)
        
        self.connectSignalToSlot()

        self.tabBar.setCloseButtonDisplayMode(TabCloseButtonDisplayMode.ON_HOVER)

        session = Session([
            "/Users/andycasey/Downloads/hd122563blue_multi.fits",
            "/Users/andycasey/Downloads/hd122563red_multi.fits"
        ])

        self.myInterface = SessionInterface(session, self)        
        self.addMySubInterface(self.myInterface, 'myInterface', 'HD 122563')

        qrouter.setDefaultRouteKey(
            self.stackedWidget, self.myInterface.objectName())
        
        
        
    def connectSignalToSlot(self):
        self.tabBar.tabAddRequested.connect(self.addTab)
        self.tabBar.tabCloseRequested.connect(self.removeTab)

        self.stackedWidget.currentChanged.connect(self.onCurrentIndexChanged)


    def addSubInterface(self, widget: QLabel, objectName, text, icon):
        widget.setObjectName(objectName)
        widget.setAlignment(Qt.AlignTop | Qt.AlignLeft)
        self.stackedWidget.addWidget(widget)
        self.tabBar.addTab(
            routeKey=objectName,
            text=text,
            icon=icon,
            onClick=lambda: self.stackedWidget.setCurrentWidget(widget)
        )

    def addMySubInterface(self, thing, objectName, text):
        thing.setObjectName(objectName)
        thing.setAlignment(Qt.AlignTop | Qt.AlignLeft)

        self.stackedWidget.addWidget(thing)
        self.tabBar.addTab(
            routeKey=objectName,
            text=text,
            icon=None,
            onClick=lambda: self.stackedWidget.setCurrentWidget(thing)            
        )


    def onCurrentIndexChanged(self, index):
        widget = self.stackedWidget.widget(index)
        if not widget:
            return

        self.tabBar.setCurrentTab(widget.objectName())
        qrouter.push(self.stackedWidget, widget.objectName())


    def addTab(self):
        # load new file to select
        filenames, _ = QFileDialog.getOpenFileNames(
            self, 
            caption="Select files", 
            directory="", 
            filter=f"FITS (*.fits *.fit *.fits.gz *.fit.gz);;CSV (*.csv);;ASCII (*.txt);;{ACRONYM} sessions (*.{FILENAME_SUFFIX});;All files (*)"
        )
        if filenames:            
            
            session = Session(filenames)
            
            # Get a suggested name.
            try:
                name = getval(filenames[0], "OBJECT", 0)
            except:
                name = f"Untitled-{self.tabCount}"
            
            # load the file, parse a name from it.        
            widget = SessionInterface(session, self)

            self.addMySubInterface(widget, name, name)        
            self.tabCount += 1
            
            # Set the current session to the newest one
            self.stackedWidget.setCurrentWidget(widget)
            self.tabBar.setCurrentTab(widget.objectName())

            # Make sure the analysis tab is in view
            interface, = self.parent.parent.findChildren(AnalysisInterface)
            self.parent.parent.stackedWidget.setCurrentWidget(interface, False)
            


    def removeTab(self, index):
        # ask to save before quit
        item = self.tabBar.tabItem(index)
        print(index, item, item.routeKey)
        widget = self.findChild(QLabel, item.routeKey())

        self.stackedWidget.removeWidget(widget)
        self.tabBar.removeTab(index)
        self.tabCount -= 1
        try:
            widget.deleteLater()
        except:
            print("cannot delete later")


class AnalysisInterface(ScrollArea):
    """ Analysis Interface """

    def __init__(self, parent=None):
        super().__init__(parent=parent)
        self.parent = parent
        self.view = SessionTabsInterface(self)        
        self.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.setViewportMargins(0, 0, 0, 0)
        self.setWidget(self.view)
        self.setWidgetResizable(True)
        self.view.setObjectName('view')
        #StyleSheet.GALLERY_INTERFACE.apply(self)
        self.setObjectName('AnalysisInterface')

