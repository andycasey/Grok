# coding: utf-8
from Qt.QtCore import QUrl, QSize, QObject, QEvent, QProcess
from Qt.QtGui import QIcon, QDesktopServices, QKeyEvent
from Qt.QtWidgets import QApplication

from qfluentwidgets import (Action, NavigationAvatarWidget, NavigationItemPosition, MessageBox, FluentWindow,
                            SplashScreen, FolderListDialog, NavigationWidget)
from qfluentwidgets import FluentIcon as FIF

from .home_interface import HomeInterface
from .setting_interface import SettingInterface
from .analysis_interface import AnalysisInterface

from ..common.config import NAME, SUPPORT_URL, cfg
from ..common.icon import Icon
from ..common.signal_bus import signalBus
from ..common.translator import Translator
from ..common import resource # removing this causes pixmap is a null pixmap warnings

from Grok.utils import suppress_stderr

from astropy.table import Table

class KeyPressFilter(QObject):

    def eventFilter(self, widget, event):
        try:
            print(f"__debug: {widget} {event.type()} {event.key()} {event.text()}")
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
    
#from Grok.synthesis.korg import QKorgProcess

class MainWindow(FluentWindow):

    def __init__(self, filenames=None):        
        super().__init__()        
        self.initWindow()        

        #self.korg_process = QKorgProcess()
        
        # create sub interface
        self.homeInterface = HomeInterface(self)
        self.analysisInterface = AnalysisInterface(self)
        self.settingInterface = SettingInterface(self)        

        # enable acrylic effect
        self.navigationInterface.setAcrylicEnabled(True)        

        self.connectSignalToSlot()        

        # add items to navigation interface
        self.initNavigation()        
        self.splashScreen.finish()        
                
        self.installEventFilter(KeyPressFilter(parent=self))           

        if filenames:
            print(f"adding with filenames: {filenames}")
            self.analysisInterface.view.createTab(filenames)        


    def connectSignalToSlot(self):
        signalBus.micaEnableChanged.connect(self.setMicaEffectEnabled)
        signalBus.switchToSampleCard.connect(self.switchToSample)
        #signalBus.supportSignal.connect(self.onSupport)

    def initNavigation(self):
        # add navigation items
        t = Translator()
        self.addSubInterface(self.homeInterface, FIF.HOME, self.tr('Home'))
        self.addSubInterface(self.analysisInterface, FIF.APPLICATION, self.tr('Analysis'))
        self.navigationInterface.addSeparator()
        
        # add an action for new analysis
        self.navigationInterface.addItem(
            "newAnalysis",
            FIF.ADD, 
            "New analysis",
            onClick=self.analysisInterface.view.addTab,
            tooltip="Start new analysis",
            selectable=False
        )
        
        # add custom widget to bottom
        self.navigationInterface.addWidget(
            routeKey='avatar',
            widget=NavigationAvatarWidget(cfg.get(cfg.userName), ':/gallery/images/shoko.png'),
            onClick=self.onAvatarClicked,
            position=NavigationItemPosition.BOTTOM
        )
        self.addSubInterface(
            self.settingInterface, FIF.SETTING, self.tr('Settings'), NavigationItemPosition.BOTTOM)



    def initWindow(self):
        self.resize(900, 1024)
        self.setMinimumWidth(780)        
        self.setWindowIcon(QIcon(':/gallery/images/logo.png'))   
        with suppress_stderr():
            # Suppress the following warning:
            # > qt.qpa.fonts: Populating font family aliases took 338 ms. Replace uses of missing font family "Segoe UI" with one that exists to avoid this cost. 
            self.setWindowTitle(NAME)        
        print("B")

        self.setMicaEffectEnabled(cfg.get(cfg.micaEnabled))        

        # create splash screen        
        self.splashScreen = SplashScreen(self.windowIcon(), self)
        self.splashScreen.setIconSize(QSize(106, 106))
        self.splashScreen.raise_()        

        desktop = QApplication.desktop().availableGeometry()
        w, h = desktop.width(), desktop.height()
        self.move(w//2 - self.width()//2, h//2 - self.height()//2)
        self.show()        
        QApplication.processEvents()


    def onAvatarClicked(self):
        pass


    def resizeEvent(self, e):
        super().resizeEvent(e)
        if hasattr(self, 'splashScreen'):
            self.splashScreen.resize(self.size())

    def switchToSample(self, routeKey, index):
        """ switch to sample """
        interfaces = self.findChildren(HomeInterface)
        for w in interfaces:
            if w.objectName() == routeKey:
                self.stackedWidget.setCurrentWidget(w, False)
                w.scrollToCard(index)
