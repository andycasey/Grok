


def get_qt_backend_name():

    from qfluentwidgets import QObject
    return QObject.__mro__[0].__module__.split(".")[0]


