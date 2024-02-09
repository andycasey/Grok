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
from Grok.session import Session

from ..common.config import (ACRONYM, EXAMPLE_URL, FEEDBACK_URL,
                             FILENAME_SUFFIX, HELP_URL, cfg)
from ..common.icon import Icon
from ..common.style_sheet import StyleSheet
from ..components.tool_bar import ToolBar
from ..components.analysis import AnalysisWidget
from ..components.quick_view import QuickViewWidget
from ..view.gallery_interface import SeparatorWidget




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
    






'''
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
'''
    

        


class SessionInterface(ScrollArea):
    """ Dialog interface """

    def __init__(self, session, parent=None):
        super().__init__(parent=parent)
        self.parent = parent
        
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
        
        qv = QuickViewWidget(session, parent=self)
        self.vBoxLayout.addWidget(qv)
        

    def add_analysis_widget(self, widget_class):
        widget = widget_class(self.session, parent=self)
        
        # not sure that this is the right thing to do here, but without the
        # visible False/True and the processEvents, it doesn't scroll to the new widget
        widget.setVisible(False)
        self.vBoxLayout.addWidget(widget, 0, Qt.AlignTop)
        widget.setVisible(True)
                
        QApplication.processEvents()
        
        w = self.vBoxLayout.itemAt(self.vBoxLayout.count() - 1).widget()        
        self.verticalScrollBar().setValue(w.y())
        

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
 
        """
        session = Session(
            [
                "/Users/andycasey/Downloads/hd122563blue_multi.fits",
                "/Users/andycasey/Downloads/hd122563red_multi.fits"
            ],
            #synthesis=self.parent.parent.korg_process
        )

        self.myInterface = SessionInterface(session, self)        
        self.addMySubInterface(self.myInterface, 'myInterface', '18 Sco')
        #HD 122563')
        """
        #qrouter.setDefaultRouteKey(
        #    self.stackedWidget, 
        #    self.myInterface.objectName()
        #)
        
        
        
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
            
            session = Session(
                filenames, 
                #synthesis=self.parent.parent.korg_process
            )
            
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

