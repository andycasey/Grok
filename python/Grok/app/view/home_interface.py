# coding:utf-8
import os
from PyQt5.QtCore import Qt, QRectF
from PyQt5.QtGui import QPixmap, QPainter, QColor, QBrush, QPainterPath, QLinearGradient
from PyQt5.QtWidgets import QWidget, QVBoxLayout, QLabel, QHBoxLayout

from qfluentwidgets import ScrollArea, isDarkTheme, FluentIcon, PushButton
from ..common.config import cfg, NAME, HELP_URL, REPO_URL, EXAMPLE_URL, FEEDBACK_URL
from ..common.icon import Icon, FluentIconBase
from ..components.link_card import LinkCardView
from ..components.sample_card import SampleCardView
from ..common.style_sheet import StyleSheet
from ..components.updates_card import UpdatesCardView

class BannerWidget(QWidget):
    """ Banner widget """

    def __init__(self, parent=None):
        super().__init__(parent=parent)
        self.setFixedHeight(336)

        self.vBoxLayout = QVBoxLayout(self)
        self.galleryLabel = QLabel(NAME, self)
        self.banner = QPixmap(':/gallery/images/header1.png')
        self.linkCardView = LinkCardView(self)

        self.galleryLabel.setObjectName('galleryLabel')

        self.vBoxLayout.setSpacing(0)
        self.vBoxLayout.setContentsMargins(0, 20, 0, 0)
        self.vBoxLayout.addWidget(self.galleryLabel)
        self.vBoxLayout.addWidget(self.linkCardView, 1, Qt.AlignBottom)
        self.vBoxLayout.setAlignment(Qt.AlignLeft | Qt.AlignTop)

        self.linkCardView.addCard(
            ':/gallery/images/logo.png',
            self.tr('Getting started'),
            self.tr('A guide to analyzing high-resolution stellar spectra.'),
            HELP_URL
        )

        self.linkCardView.addCard(
            FluentIcon.GITHUB,
            self.tr('GitHub repository'),
            self.tr(
                'The source code repository announces updates and potential issues.'),
            REPO_URL
        )

        self.linkCardView.addCard(
            FluentIcon.CODE,
            self.tr('Code samples'),
            self.tr(
                'Automate your analysis.'),
            EXAMPLE_URL
        )

        self.linkCardView.addCard(
            FluentIcon.FEEDBACK,
            self.tr('Send feedback'),
            self.tr(f'Provide feedback and help us improve {NAME}.'),
            FEEDBACK_URL
        )

    def paintEvent(self, e):
        super().paintEvent(e)
        painter = QPainter(self)
        painter.setRenderHints(
            QPainter.SmoothPixmapTransform | QPainter.Antialiasing)
        painter.setPen(Qt.NoPen)

        path = QPainterPath()
        path.setFillRule(Qt.WindingFill)
        w, h = self.width(), self.height()
        path.addRoundedRect(QRectF(0, 0, w, h), 10, 10)
        path.addRect(QRectF(0, h-50, 50, 50))
        path.addRect(QRectF(w-50, 0, 50, 50))
        path.addRect(QRectF(w-50, h-50, 50, 50))
        path = path.simplified()

        # init linear gradient effect
        gradient = QLinearGradient(0, 0, 0, h)

        # draw background color
        if not isDarkTheme():
            gradient.setColorAt(0, QColor(207, 216, 228, 255))
            gradient.setColorAt(1, QColor(207, 216, 228, 0))
        else:
            gradient.setColorAt(0, QColor(0, 0, 0, 255))
            gradient.setColorAt(1, QColor(0, 0, 0, 0))
            
        painter.fillPath(path, QBrush(gradient))

        # draw banner image
        pixmap = self.banner.scaled(
            self.size(), transformMode=Qt.SmoothTransformation)
        painter.fillPath(path, QBrush(pixmap))


class HomeInterface(ScrollArea):
    """ Home interface """

    def __init__(self, parent=None):
        super().__init__(parent=parent)
        self.parent = parent
        self.banner = BannerWidget(self)
        self.view = QWidget(self)
        self.vBoxLayout = QVBoxLayout(self.view)

        self.__initWidget()
        
        button_layout = QHBoxLayout()
        button_layout.setContentsMargins(40, 0, 0, 0)
        self.button_pre_compile_korg = PushButton("Pre-compile Korg")
        self.button_pre_compile_korg.clicked.connect(self.pre_compile_korg)
        
        button_layout.addWidget(self.button_pre_compile_korg)
        button_layout.addStretch(1)
        self.vBoxLayout.addLayout(button_layout)
        
        self.show_korg_updates()
    
    def pre_compile_korg(self):
        
        '''
        self.cog_line_list = self.parent.korg_process.read_linelist(
            "/Users/andycasey/Downloads/linelist_mm.txt", 
            format="moog",
        )
        print(f"got cog line list: {self.cog_line_list}")
        print(f"OK DONE")
        
        def my_callback(ll):
            print(f"got a ll with {len(ll)} entries: {ll[0]}")
        
        self.parent.korg_process.async_read_linelist(
            "/Users/andycasey/Downloads/linelist_mm.txt", 
            format="moog",
            callback=my_callback            
        )
        '''
        
        '''
        
        ll = self.parent.korg_process.read_linelist(
            "/Users/andycasey/research/Grok/python/Grok/Melendez_et_al_Fe.moog", 
            format="moog",
        )  
        A_X = self.parent.korg_process.format_A_X(0.0, 0.0)
        atm = self.parent.korg_process.interpolate_marcs(5777, 4.4)
        import numpy as np
        ews = list(np.loadtxt("/Users/andycasey/research/Grok/python/Grok/Melendez_et_al_Fe.moog", usecols=(5, ), skiprows=1))
        foo = self.parent.korg_process.ews_to_abundances(atm, ll, A_X, ews)
        '''
        #self.parent.korg_process.format_A_X(0, 0))
        #print(f"got foo: {foo}")
    

    def __initWidget(self):
        self.view.setObjectName('view')
        self.setObjectName('homeInterface')
        StyleSheet.HOME_INTERFACE.apply(self)

        self.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.setWidget(self.view)
        self.setWidgetResizable(True)

        self.vBoxLayout.setContentsMargins(0, 0, 0, 36)
        self.vBoxLayout.setSpacing(40)
        self.vBoxLayout.addWidget(self.banner)
        self.vBoxLayout.setAlignment(Qt.AlignTop)



    def show_korg_updates(self):
        
        RecentUpdatesView = UpdatesCardView(
            self.tr("Recent updates"), 
            self.view
        )        
        RecentUpdatesView.addSampleCard(
            icon=f"{os.path.dirname(__file__)}/../resource/images/icons/Korg_black.svg",
            title="Korg v0.27.1 released",
            content=self.tr(
                """
                Merged pull requests:

- comments, type hinting (#242) (@ajwheeler)
- tweak parsing of turbospectrum linelists so that isotopes can be spec… (#243) (@ajwheeler)
- properly support [alpha/H] as a fitting parameter (#244) (@ajwheeler)"""
                ),
            routeKey="basicInputInterface",
            index=0
            
        )
        self.vBoxLayout.addWidget(RecentUpdatesView)
        RecentUpdatesView.setCursor(Qt.PointingHandCursor)

        
