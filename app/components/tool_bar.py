from PyQt5.QtCore import Qt, QUrl
from PyQt5.QtGui import QColor, QDesktopServices
from PyQt5.QtWidgets import QHBoxLayout, QVBoxLayout, QWidget, QAction, QSizePolicy
from qfluentwidgets import CaptionLabel, FluentIcon, PushButton, TitleLabel, DropDownToolButton, RoundMenu

from ..common.config import EXAMPLE_URL, HELP_URL
from ..view.gallery_interface import SeparatorWidget

from astropy.io.fits import getval

class ToolBar(QWidget):
    """ Tool bar """

    def __init__(self, session, title=None, subtitle=None, parent=None):
        super().__init__(parent=parent)
        
        if title is None:
            try:
                title = getval(session.input_paths[0], "OBJECT", 0)
            except:
                title = "Untitled"
                        
        self.titleLabel = TitleLabel(title, self)
    
        if subtitle is None:        
            try:
                ra = getval(session.input_paths[0], "RA", 0)
                dec = getval(session.input_paths[0], "DEC", 0)
            except:
                subtitle = ""
            else:
                subtitle = f"{ra} {dec}"

        self.subtitleLabel = CaptionLabel(subtitle, self)


        self.saveButton = PushButton(
            self.tr('Save'), 
            self, 
            FluentIcon.SAVE
        )
        # add 'save as'
        self.saveAsButton = PushButton(
            self.tr('Save as'), 
            self, 
            FluentIcon.SAVE_AS
        )
        # add 'export'
        self.exportButton = PushButton(
            self.tr('Export tables'),
            self,
            FluentIcon.DOCUMENT
        )
        self.labelButton = PushButton(
            "Label",
            self,
            FluentIcon.TAG
        )
        
        self.searchMenu = RoundMenu(parent=self)
        self.searchMenu.addAction(QAction("Search SIMBAD", self))
        self.searchMenu.addAction(QAction("Search Gaia data archive", self))
        self.searchMenu.addAction(QAction("Search NASA/ADS", self))

        self.dropDownSearchButton = DropDownToolButton(FluentIcon.SEARCH, self)
        self.dropDownSearchButton.setMenu(self.searchMenu)
                


        self.vBoxLayout = QVBoxLayout(self)
        self.buttonLayout = QHBoxLayout(self)

        self.__initWidget()
        

    def __initWidget(self):
        self.setFixedHeight(138)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.vBoxLayout.setSpacing(0)
        self.vBoxLayout.setContentsMargins(36, 22, 36, 12)
        self.vBoxLayout.addWidget(self.titleLabel)
        self.vBoxLayout.addSpacing(4)
        self.vBoxLayout.addWidget(self.subtitleLabel)
        self.vBoxLayout.addSpacing(4)
        self.vBoxLayout.addLayout(self.buttonLayout, 1)
        self.vBoxLayout.setAlignment(Qt.AlignTop)

        self.buttonLayout.setSpacing(4)
        self.buttonLayout.setContentsMargins(0, 0, 0, 0)
        self.buttonLayout.addWidget(self.saveButton, 0, Qt.AlignLeft)
        self.buttonLayout.addWidget(self.saveAsButton, 0, Qt.AlignLeft)
        self.buttonLayout.addWidget(self.labelButton, 0, Qt.AlignLeft)
        self.buttonLayout.addWidget(self.exportButton, 0, Qt.AlignLeft)
        # add separator
        self.buttonLayout.addWidget(SeparatorWidget(), 0, Qt.AlignLeft)
        self.buttonLayout.addStretch(1)    
        self.buttonLayout.addWidget(self.dropDownSearchButton)#, alignment=Qt.AlignRight)
        #self.buttonLayout.setAlignment(Qt.AlignVCenter)


        self.subtitleLabel.setTextColor(QColor(96, 96, 96), QColor(216, 216, 216))
