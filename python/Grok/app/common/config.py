# coding:utf-8
import sys
import os
from enum import Enum

#from Qt.Qtcore import QLocale
from PyQt5.QtCore import QLocale
from pathlib import Path
from qfluentwidgets import (qconfig, QConfig, ConfigItem, OptionsConfigItem, BoolValidator,
                            OptionsValidator, RangeConfigItem, RangeValidator, ConfigValidator,
                            FolderListValidator, Theme, FolderValidator, ConfigSerializer, __version__)


class Language(Enum):
    """ Language enumeration """

    CHINESE_SIMPLIFIED = QLocale(QLocale.Chinese, QLocale.China)
    CHINESE_TRADITIONAL = QLocale(QLocale.Chinese, QLocale.HongKong)
    ENGLISH = QLocale(QLocale.English)
    AUTO = QLocale()


class LanguageSerializer(ConfigSerializer):
    """ Language serializer """

    def serialize(self, language):
        return language.value.name() if language != Language.AUTO else "Auto"

    def deserialize(self, value: str):
        return Language(QLocale(value)) if value != "Auto" else Language.AUTO


def isWin11():
    return sys.platform == 'win32' and sys.getwindowsversion().build >= 22000

class OptionalPathValidator(ConfigValidator):
    
    def validate(self, value):
        return value is None or Path(value).exists()
    
    def correct(self, value):
        if value is None:
            return ""
        path = Path(value)
        return f"{path.absolute()}"
    

class Config(QConfig):
    """ Config of application """

    # folders
    musicFolders = ConfigItem(
        "Folders", "LocalMusic", [], FolderListValidator())
    #downloadFolder = ConfigItem(
    #    "Folders", "Download", "app/download", FolderValidator())


    # main window
    micaEnabled = ConfigItem("MainWindow", "MicaEnabled", isWin11(), BoolValidator())
    dpiScale = OptionsConfigItem(
        "MainWindow", "DpiScale", "Auto", OptionsValidator([1, 1.25, 1.5, 1.75, 2, "Auto"]), restart=True)
    language = OptionsConfigItem(
        "MainWindow", "Language", Language.AUTO, OptionsValidator(Language), LanguageSerializer(), restart=True)

    # Material
    blurRadius  = RangeConfigItem("Material", "AcrylicBlurRadius", 15, RangeValidator(0, 40))

    # software update
    checkUpdateAtStartUp = ConfigItem("Update", "CheckUpdateAtStartUp", True, BoolValidator())

    userName = ConfigItem("Personalization", "User name", os.getlogin())

    # Radial velocity
    
    rv_default_continuum_type = OptionsConfigItem(
        "Radial Velocity",
        "rv_default_continuum_type",
        "Polynomial",
        OptionsValidator([
        #    "Sinusoids",
            "Polynomial",
        #    "Spline"
        ])
    )
    
    rv_default_wavelength_range = OptionsConfigItem(
        "Radial Velocity",
        "rv_default_wavelength_range",
        "8400 - 8800 Å",
        OptionsValidator([
            "8400 - 8800 Å",
            "5100 - 5200 Å",
        ])
    )
    rv_default_template_type = OptionsConfigItem(
        "Radial Velocity", 
        "rv_default_template_type",
        "From file",
        OptionsValidator([
            #"Auto", 
            "From file", 
            #"Synthesize spectrum"
        ]),
    )
    rv_default_template_path = ConfigItem(
        "Radial Velocity", "rv_default_template_path", None, OptionalPathValidator()
    )
    
    rv_default_continuum_degree = OptionsConfigItem(
        "Radial Velocity",
        "rv_default_continuum_degree",
        "2",
        OptionsValidator(list(map(str, range(0, 10))))
    )


    
    
    initial_teff = ConfigItem("Initial Stellar Parameters", "Initial effective temperature", 5777, RangeValidator(2500, 10_000))
    


YEAR = 2023
NAME = "Grok"
ACRONYM = "Grok"
FILENAME_SUFFIX = "grok"
AUTHOR = f"{NAME} Collaboration"
VERSION = __version__
HELP_URL = f"https://github.com/andycasey/{NAME}"
REPO_URL = HELP_URL
EXAMPLE_URL = HELP_URL
FEEDBACK_URL = f"https://github.com/andycasey/{NAME}/issues"
RELEASE_URL = f"https://github.com/andycasey/{NAME}/releases/latest"
SUPPORT_URL = FEEDBACK_URL


cfg = Config()
cfg.themeMode.value = Theme.AUTO
qconfig.load('app/config/config.json', cfg)
