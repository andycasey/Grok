from Qt import QtCore
from Qt.QtWidgets import QWidget, QVBoxLayout, QHBoxLayout, QFrame

from qfluentwidgets import StrongBodyLabel

class AnalysisWidget(QWidget):
    
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

        self.vBoxLayout.addWidget(self.titleLabel, 0, QtCore.Qt.AlignTop)
        self.vBoxLayout.addWidget(self.card, 0, QtCore.Qt.AlignTop)
        self.vBoxLayout.setAlignment(QtCore.Qt.AlignTop)

        self.cardLayout.setSpacing(0)
        self.cardLayout.setAlignment(QtCore.Qt.AlignTop)
        self.cardLayout.addLayout(self.topLayout, 0)
        self.cardLayout.addWidget(self.sourceWidget, 0, QtCore.Qt.AlignBottom)

        self.topLayout.addWidget(self.widget)
        if self.stretch == 0:
            self.topLayout.addStretch(1)


        '''
        #self.bottomLayout.addWidget(self.sourcePathLabel, 0, QtCore.Qt.AlignLeft)
        #self.bottomLayout.addStretch(1)
        #self.bottomLayout.addWidget(self.linkIcon, 0, QtCore.Qt.AlignRight)
        
        if self.leftButtonText is not None:
            leftButton = PushButton(self.leftButtonText)
            #leftButton.clicked.connect(self.enable_norm)
            self.bottomLayout.addWidget(leftButton, 0, QtCore.Qt.AlignLeft)
            
        self.bottomLayout.addStretch(1)
        if self.rightButtonText is not None:
            rightButton = PrimaryPushButton(self.rightButtonText, self)
            #rightButton.clicked.connect(self.enable_norm)
            
            self.bottomLayout.addWidget(rightButton, 0, QtCore.Qt.AlignRight)

        '''
        self.bottomLayout.setAlignment(QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter)
        
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