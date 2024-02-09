# coding:utf-8
import os
import sys

from Qt.QtCore import Qt, QTranslator
from Qt.QtWidgets import QApplication
from Grok.utils import suppress

with suppress():
    # No advertising, please.
    from qfluentwidgets import FluentTranslator

from Grok.app.common.config import cfg
from Grok.app.view.main_window import MainWindow


def main():
        
    # enable dpi scale
    if cfg.get(cfg.dpiScale) == "Auto":
        QApplication.setHighDpiScaleFactorRoundingPolicy(
            Qt.HighDpiScaleFactorRoundingPolicy.PassThrough)
        QApplication.setAttribute(Qt.AA_EnableHighDpiScaling)
    else:
        os.environ["QT_ENABLE_HIGHDPI_SCALING"] = "0"
        os.environ["QT_SCALE_FACTOR"] = str(cfg.get(cfg.dpiScale))

    QApplication.setAttribute(Qt.AA_UseHighDpiPixmaps)

    # create application
    app = QApplication(sys.argv)
    app.setAttribute(Qt.AA_DontCreateNativeWidgetSiblings)

    # internationalization
    locale = cfg.get(cfg.language).value
    translator = FluentTranslator(locale)
    galleryTranslator = QTranslator()
    galleryTranslator.load(locale, "gallery", ".", ":/gallery/i18n")

    app.installTranslator(translator)
    app.installTranslator(galleryTranslator)

    # create main window
    w = MainWindow(sys.argv[1:])
    w.show()

    app.exec_()


if __name__ == "__main__":
    main()