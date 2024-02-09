# coding:utf-8
import os
from qfluentwidgets import (SettingCardGroup, SwitchSettingCard, FolderListSettingCard,
                            OptionsSettingCard, PushSettingCard, SettingCard,
                            HyperlinkCard, PrimaryPushSettingCard, ScrollArea,
                            ComboBoxSettingCard, ExpandLayout, Theme, CustomColorSettingCard,
                            setTheme, setThemeColor, RangeSettingCard, isDarkTheme)
from qfluentwidgets import FluentIcon as FIF
from qfluentwidgets import InfoBar
from PyQt5.QtCore import Qt, pyqtSignal, QUrl, QStandardPaths
from PyQt5.QtGui import QDesktopServices
from PyQt5.QtWidgets import QWidget, QLabel, QFileDialog

from ..common.config import cfg, HELP_URL, FEEDBACK_URL, AUTHOR, VERSION, YEAR, isWin11
from ..common.signal_bus import signalBus
from ..common.style_sheet import StyleSheet

# import what we need
from typing import Union

from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QColor, QIcon, QPainter
from PyQt5.QtWidgets import (QFrame, QHBoxLayout, QLabel, QToolButton,
                             QVBoxLayout, QPushButton)

from qfluentwidgets import LineEdit

from ..common.config import qconfig, ConfigItem
from ..common.icon import FluentIconBase

class TextSettingCard(SettingCard):
    
    textEdited = pyqtSignal(str)
    
    def __init__(self, icon: Union[str, QIcon, FluentIconBase], title, content=None,
                 configItem: ConfigItem = None, parent=None):
        """
        Parameters
        ----------
        icon: str | QIcon | FluentIconBase
            the icon to be drawn

        title: str
            the title of card

        content: str
            the content of card

        configItem: ConfigItem
            configuration item operated by the card

        parent: QWidget
            parent widget
        """
        super().__init__(icon, title, content, parent)
        self.configItem = configItem
        self.lineEdit = LineEdit(self)
        
        if configItem:
            self.setValue(qconfig.get(configItem))
            configItem.valueChanged.connect(self.setValue)

        # add switch button to layout
        self.hBoxLayout.addWidget(self.lineEdit, 0, Qt.AlignRight)
        self.hBoxLayout.addSpacing(16)

        self.lineEdit.textEdited.connect(self.__onTextEdited)

    def __onTextEdited(self):
        """ switch button checked state changed slot """
        value = self.lineEdit.text()
        print(f"onTextEdited: {value}")
        self.setValue(value)
        self.textEdited.emit(value)


    def setValue(self, value, **kwargs):
        print(f"setValue: {value} {kwargs}")
        if self.configItem:
            try:            
                qconfig.set(self.configItem, int(value)) # HACK
            except:
                # TODO
                None

        self.lineEdit.setText(str(value))
        



class SettingInterface(ScrollArea):
    """ Setting interface """

    def __init__(self, parent=None):
        super().__init__(parent=parent)
        self.scrollWidget = QWidget()
        self.expandLayout = ExpandLayout(self.scrollWidget)

        # setting label
        self.settingLabel = QLabel(self.tr("Settings"), self)

        # music folders
        self.musicInThisPCGroup = SettingCardGroup(
            self.tr("Music on this PC"), self.scrollWidget)
        self.musicFolderCard = FolderListSettingCard(
            cfg.musicFolders,
            self.tr("Local music library"),
            directory=QStandardPaths.writableLocation(
                QStandardPaths.MusicLocation),
            parent=self.musicInThisPCGroup
        )


        # personalization
        self.personalGroup = SettingCardGroup(
            self.tr('Personalization'), self.scrollWidget)
        
        # RV
        self.radialVelocityParametersGroup = SettingCardGroup("Radial velocity", self.scrollWidget)
        self.rv_default_wavelength_range = ComboBoxSettingCard(
            cfg.rv_default_wavelength_range,
            FIF.LANGUAGE,
            self.tr('Default wavelength range'),
            self.tr('Set your preferred wavelength range for radial velocity measurements'),
            texts=cfg.rv_default_wavelength_range.options,
            parent=self.radialVelocityParametersGroup
        )

        '''
        self.rv_default_template_type = OptionsSettingCard(
            cfg.rv_default_template_type,
            FIF.BRUSH,
            self.tr('Default template type'),
            self.tr("Default template type"),
            texts=[
                self.tr("Auto"),
                self.tr("Select from file"),
                self.tr("Synthesize spectrum"),
                self.tr('Use system setting')                
            ],
            parent=self.radialVelocityParametersGroup
        )
        '''        
        self.rv_default_template_type = ComboBoxSettingCard(
            cfg.rv_default_template_type,
            FIF.BRUSH,
            'Default template type',
            "Set the default rest-frame template type for radial velocity measurements",
            texts=cfg.rv_default_template_type.options,
            parent=self.radialVelocityParametersGroup
        )        
        self.rv_default_template_path = PushSettingCard(
            "Select template spectrum",
            FIF.DOWNLOAD,
            "Radial velocity template path",
            cfg.get(cfg.rv_default_template_path),
            self.radialVelocityParametersGroup
        )        
        
        
        
        self.rv_default_continuum_type = ComboBoxSettingCard(
            cfg.rv_default_continuum_type,
            FIF.BRUSH,
            "Default continuum function",
            "Set the default continuum function for radial velocity measurements",
            texts=cfg.rv_default_continuum_type.options,
            parent=self.radialVelocityParametersGroup
        )
        self.rv_default_continuum_degree = ComboBoxSettingCard(
            cfg.rv_default_continuum_degree,
            FIF.BRUSH,
            "Default polynomial degree for continuum",
            "Set the default polynomial degree for continuum",
            texts=cfg.rv_default_continuum_degree.options,
            parent=self.radialVelocityParametersGroup
        )

        # initial stellar params
        self.initialStellarParametersGroup = SettingCardGroup("Curve of growth analysis", self.scrollWidget)
        self.initialTeffCard = TextSettingCard(
            FIF.TRANSPARENT,
            "Initial effective temperature",
            content="Used for fitting",
            configItem=cfg.initial_teff,
            parent=self.initialStellarParametersGroup
        )
        
        '''
        self.userName = SettingCard(
            FIF.PEOPLE,
            self.tr('User name'),
            #self.tr('Display user name in the navigation bar'),
            cfg.userName,
            self.personalGroup
        )
        '''
        
        
        self.micaCard = SwitchSettingCard(
            FIF.TRANSPARENT,
            self.tr('Mica effect'),
            self.tr('Apply semi transparent to windows and surfaces'),
            cfg.micaEnabled,
            self.personalGroup
        )
        self.themeCard = OptionsSettingCard(
            cfg.themeMode,
            FIF.BRUSH,
            self.tr('Application theme'),
            self.tr("Change the appearance of your application"),
            texts=[
                self.tr('Light'), self.tr('Dark'),
                self.tr('Use system setting')
            ],
            parent=self.personalGroup
        )
        self.themeColorCard = CustomColorSettingCard(
            cfg.themeColor,
            FIF.PALETTE,
            self.tr('Theme color'),
            self.tr('Change the theme color of you application'),
            self.personalGroup
        )
        self.zoomCard = OptionsSettingCard(
            cfg.dpiScale,
            FIF.ZOOM,
            self.tr("Interface zoom"),
            self.tr("Change the size of widgets and fonts"),
            texts=[
                "100%", "125%", "150%", "175%", "200%",
                self.tr("Use system setting")
            ],
            parent=self.personalGroup
        )
        self.languageCard = ComboBoxSettingCard(
            cfg.language,
            FIF.LANGUAGE,
            self.tr('Language'),
            self.tr('Set your preferred language for UI'),
            texts=['简体中文', '繁體中文', 'English', self.tr('Use system setting')],
            parent=self.personalGroup
        )

        # material
        self.materialGroup = SettingCardGroup(
            self.tr('Material'), self.scrollWidget)
        self.blurRadiusCard = RangeSettingCard(
            cfg.blurRadius,
            FIF.ALBUM,
            self.tr('Acrylic blur radius'),
            self.tr('The greater the radius, the more blurred the image'),
            self.materialGroup
        )

        # update software
        self.updateSoftwareGroup = SettingCardGroup(
            self.tr("Software update"), self.scrollWidget)
        self.updateOnStartUpCard = SwitchSettingCard(
            FIF.UPDATE,
            self.tr('Check for updates when the application starts'),
            self.tr('The new version will be more stable and have more features'),
            configItem=cfg.checkUpdateAtStartUp,
            parent=self.updateSoftwareGroup
        )

        # application
        self.aboutGroup = SettingCardGroup(self.tr('About'), self.scrollWidget)
        self.helpCard = HyperlinkCard(
            HELP_URL,
            self.tr('Open help page'),
            FIF.HELP,
            self.tr('Help'),
            self.tr(
                'Discover new features and learn useful tips about PyQt-Fluent-Widgets'),
            self.aboutGroup
        )
        self.feedbackCard = PrimaryPushSettingCard(
            self.tr('Provide feedback'),
            FIF.FEEDBACK,
            self.tr('Provide feedback'),
            self.tr('Help us improve PyQt-Fluent-Widgets by providing feedback'),
            self.aboutGroup
        )
        self.aboutCard = PrimaryPushSettingCard(
            self.tr('Check update'),
            FIF.INFO,
            self.tr('About'),
            '© ' + self.tr('Copyright') + f" {YEAR}, {AUTHOR}. " +
            self.tr('Version') + " " + VERSION,
            self.aboutGroup
        )

        self.__initWidget()

    def __initWidget(self):
        self.resize(1000, 800)
        self.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.setViewportMargins(0, 80, 0, 20)
        self.setWidget(self.scrollWidget)
        self.setWidgetResizable(True)
        self.setObjectName('settingInterface')

        # initialize style sheet
        self.scrollWidget.setObjectName('scrollWidget')
        self.settingLabel.setObjectName('settingLabel')
        StyleSheet.SETTING_INTERFACE.apply(self)

        self.micaCard.setEnabled(isWin11())

        # initialize layout
        self.__initLayout()
        self.__connectSignalToSlot()

    def __initLayout(self):
        self.settingLabel.move(36, 30)

        # add cards to group
        self.musicInThisPCGroup.addSettingCard(self.musicFolderCard)
        #self.musicInThisPCGroup.addSettingCard(self.downloadFolderCard)

        self.radialVelocityParametersGroup.addSettingCard(self.rv_default_wavelength_range)
        self.radialVelocityParametersGroup.addSettingCard(self.rv_default_template_type)
        self.radialVelocityParametersGroup.addSettingCard(self.rv_default_template_path)
        self.radialVelocityParametersGroup.addSettingCard(self.rv_default_continuum_type)
        self.radialVelocityParametersGroup.addSettingCard(self.rv_default_continuum_degree)

        self.initialStellarParametersGroup.addSettingCard(self.initialTeffCard)

        self.personalGroup.addSettingCard(self.micaCard)
        self.personalGroup.addSettingCard(self.themeCard)
        self.personalGroup.addSettingCard(self.themeColorCard)
        self.personalGroup.addSettingCard(self.zoomCard)
        self.personalGroup.addSettingCard(self.languageCard)

        self.materialGroup.addSettingCard(self.blurRadiusCard)

        self.updateSoftwareGroup.addSettingCard(self.updateOnStartUpCard)

        self.aboutGroup.addSettingCard(self.helpCard)
        self.aboutGroup.addSettingCard(self.feedbackCard)
        self.aboutGroup.addSettingCard(self.aboutCard)

        # add setting card group to layout
        self.expandLayout.setSpacing(28)
        self.expandLayout.setContentsMargins(36, 10, 36, 0)
        self.expandLayout.addWidget(self.radialVelocityParametersGroup)
        self.expandLayout.addWidget(self.musicInThisPCGroup)
        self.expandLayout.addWidget(self.initialStellarParametersGroup)
        self.expandLayout.addWidget(self.personalGroup)
        self.expandLayout.addWidget(self.materialGroup)
        self.expandLayout.addWidget(self.updateSoftwareGroup)
        self.expandLayout.addWidget(self.aboutGroup)

    def __showRestartTooltip(self):
        """ show restart tooltip """
        InfoBar.success(
            self.tr('Updated successfully'),
            self.tr('Configuration takes effect after restart'),
            duration=1500,
            parent=self
        )



    '''
    def __onDownloadFolderCardClicked(self):
        """ download folder card clicked slot """
        folder = QFileDialog.getExistingDirectory(
            self, self.tr("Choose folder"), "./")
        if not folder or cfg.get(cfg.downloadFolder) == folder:
            return

        cfg.set(cfg.downloadFolder, folder)
        self.downloadFolderCard.setContent(folder)
    '''

    def on_rv_default_template_path_clicked(self):
        try:
            directory = os.path.dirname(cfg.get(cfg.rv_default_template_path))
        except:
            directory = ""
            
        path, _ = QFileDialog.getOpenFileName(
            caption="Select template spectrum", 
            directory=directory, 
            filter="*.fits"
        )
        if not path: return

        cfg.set(cfg.rv_default_template_path, path)
        self.rv_default_template_path.setContent(path)
                
        
    def __connectSignalToSlot(self):
        """ connect signal to slot """
        cfg.appRestartSig.connect(self.__showRestartTooltip)


        # personalization
        self.rv_default_template_path.clicked.connect(self.on_rv_default_template_path_clicked)
        self.themeCard.optionChanged.connect(lambda ci: setTheme(cfg.get(ci)))
        self.themeColorCard.colorChanged.connect(setThemeColor)
        self.micaCard.checkedChanged.connect(signalBus.micaEnableChanged)

        # about
        self.feedbackCard.clicked.connect(
            lambda: QDesktopServices.openUrl(QUrl(FEEDBACK_URL)))
