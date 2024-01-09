from Qt import QtCore
from Qt.QtGui import QColor
from Qt.QtWidgets import QHBoxLayout, QVBoxLayout, QWidget, QAction, QSizePolicy
from qfluentwidgets import (
    CaptionLabel, FluentIcon, PushButton, TitleLabel, DropDownToolButton, RoundMenu,
    PrimaryDropDownPushButton, DropDownPushButton
)

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
        self.exportMenu = RoundMenu(parent=self)
        self.exportMenu.addActions([
            QAction("Spectra", self),
            QAction("Tables", self),
            QAction("Figures", self),
        ])
        self.exportButton = DropDownPushButton(FluentIcon.DOWNLOAD, "Export", self)
        self.exportButton.setMenu(self.exportMenu)
        
        self.codeButton = PushButton(
            "Script",
            self,
            FluentIcon.CODE
        )
        
        self.searchMenu = RoundMenu(parent=self)
        self.searchMenu.addActions([
            QAction("Search CDS Portal", self),
            QAction("Search SIMBAD", self),
            QAction("Search Gaia archive", self),
            QAction("Search NASA/ADS", self),
        ])

        self.dropDownSearchButton = DropDownToolButton(FluentIcon.SEARCH, self)
        self.dropDownSearchButton.setMenu(self.searchMenu)
                
        
        self.analysisMenu = RoundMenu(parent=self)
        
        self.stellarParameterMenu = RoundMenu("Stellar parameters", parent=self)
        self.stellarParameterMenu.addActions([
            QAction("Curve-of-growth", self),
            QAction("Differential analysis", self),
            QAction("Spectral fitting", self)
        ])
        
        self.analysisMenu.addAction(QAction("Radial velocity", self))
        self.analysisMenu.addAction(QAction("Continuum", self))
        self.analysisMenu.addMenu(self.stellarParameterMenu)
        self.analysisMenu.addAction(QAction("Synthesis", self))
            
        self.addAnalysis = PrimaryDropDownPushButton(FluentIcon.ADD, "Add analysis", self)
        self.addAnalysis.setMenu(self.analysisMenu)


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
        self.vBoxLayout.setAlignment(QtCore.Qt.AlignTop)

        self.buttonLayout.setSpacing(4)
        self.buttonLayout.setContentsMargins(0, 0, 0, 0)
        self.buttonLayout.addWidget(self.dropDownSearchButton)
        self.buttonLayout.addWidget(SeparatorWidget())        
        self.buttonLayout.addWidget(self.saveButton)
        self.buttonLayout.addWidget(self.saveAsButton)
        self.buttonLayout.addWidget(self.codeButton)
        self.buttonLayout.addWidget(self.exportButton)
        self.buttonLayout.addStretch(1)    
        self.buttonLayout.addWidget(self.addAnalysis, 0, QtCore.Qt.AlignRight)

        self.subtitleLabel.setTextColor(QColor(96, 96, 96), QColor(216, 216, 216))
