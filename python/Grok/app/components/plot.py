
import threading
import numpy as np
from collections import namedtuple
from weakref import WeakKeyDictionary
from PyQt5.QtCore import Qt, QSize, QTimer, pyqtSignal, QObject, QEvent
from PyQt5.QtGui import QResizeEvent
from PyQt5.QtWidgets import QWidget, QStackedWidget, QVBoxLayout, QLabel, QHBoxLayout, QFrame, QSizePolicy

from PyQt5.QtWidgets import QApplication, QFrame, QVBoxLayout, QLabel, QWidget, QHBoxLayout
from qfluentwidgets import (FluentIcon, IconWidget, FlowLayout, isDarkTheme,
                            Theme, applyThemeColor, SmoothScrollArea, SearchLineEdit, StrongBodyLabel,
                            BodyLabel, CommandBar)
from qfluentwidgets import RoundMenu, Action
from time import time

# Remember to import matplotlib after Qt.
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import FigureCanvas

# TODO: Import qconfig from qfluentwidgets and when the theme changes, update the matplotlib style file
plt.style.use({
    "font.family": "sans-serif",
    #"font.sans-serif": ["Segoe UI"],
    "axes.labelsize": "10",
    "lines.linewidth": "2",
    "lines.markersize": "10",
    "xtick.labelsize": "10",
    "ytick.labelsize": "10",

    #"figure.facecolor": "#272727f2", #-> main background color
    #"figure.facecolor": "#202020",
    "axes.facecolor": "#292929", # rgba(41,41,41,255)            
    "figure.facecolor": "#292929", # rgba(41,41,41,255)
    "axes.edgecolor": "white",
    "axes.labelcolor": "white",
    "xtick.color": "white",
    "xtick.labelcolor": "white",
    "ytick.color": "white",
    "ytick.labelcolor": "white"
})        
        
        
class ExcitationIonizationBalanceWidget(QWidget):    
    
    def __init__(
        self, 
        x=None,
        y=None,
        xlabel=None,
        ylabel=None,
        figsize=(8, 6),
        parent=None, 
        size_policy=(QSizePolicy.Expanding, QSizePolicy.Fixed),
        resize_interval=50
    ):
        super().__init__(parent)
        self.parent = parent
        self.resize_interval = resize_interval
        self.figure = Figure(figsize=figsize)
        self.canvas = FigureCanvas(self.figure)
        self.canvas.setParent(self)
        self.canvas.setFocusPolicy(Qt.WheelFocus)
        self.canvas.setFocus()
        self.canvas.setSizePolicy(*size_policy)
        
        self.canvas.mpl_connect("figure_enter_event", self._focus)
                    
        self.axes = self.canvas.figure.subplots(3, 1)
        self.figure.tight_layout()
        self.figure.canvas.draw()
        
        self.layout = QVBoxLayout(self)
        self.layout.addWidget(self.canvas)

        #self.layout.addStretch(1)    
        self.installEventFilter(self)
        return None
    
    def _focus(self, event):
        """ Set the focus of the canvas. """
        self.canvas.setFocus()    
        
        
        
    
    def resizeEvent(self, e):
        # Matplotlib wants to redraw the canvas while we are resizing, makes it all yucky
        try:
            self.resizeTimer
        except:
            
            self.resizeTimer = QTimer(self)
            self.resizeTimer.setSingleShot(True)
            self.resizeTimer.timeout.connect(self.afterResizeEvent)
        finally:
            self.resizeTimer.start(self.resize_interval)
                
        return None        

    def afterResizeEvent(self):        
        self.figure.tight_layout()
        self.figure.canvas.draw()
        try:
            self.resizeTimer.stop()
            del self.resizeTimer
        except:
            None
        return None
    
    def eventFilter(self, widget, event):
        try:
            print(f"plot widget {widget} {event.type()} {event.key()} {event.text()}")
        except:
            None
            
        if event.type() == 51:
            if event.key() == Qt.Key_Left:
                try:
                    self.page_left.trigger()
                except:
                    return False
                else:
                    return True
            elif event.key() == Qt.Key_Right:
                try:
                    self.page_right.trigger()
                except:
                    return False
                else:
                    return True
                
                                

        '''
        if event.type() == QEvent.KeyPress:
            text = event.text()
            if event.modifiers():
                text = event.keyCombination().key().name.decode(encoding="utf-8")
            print(f"{event} {event.type}: {text}")
        '''            
        return False
        

from matplotlib.backend_bases import NavigationToolbar2, _Mode, MouseButton, tools, cbook
    
class SinglePlotWidget(QWidget):    
    
    def __init__(
        self, 
        x=None,
        y=None,
        xlabel=None,
        ylabel=None,
        figsize=(8, 2),
        parent=None, 
        size_policy=(QSizePolicy.Expanding, QSizePolicy.Fixed),
        toolbar=False,
        toolbar_left_right=True,
        resize_interval=50
    ):
        super().__init__(parent)
        self.parent = parent
        self.resize_interval = resize_interval
        self.figure = Figure(figsize=figsize)
        self.canvas = FigureCanvas(self.figure)
        self.canvas.setParent(self)
        self.canvas.setFocusPolicy(Qt.WheelFocus)
        self.canvas.setFocus()
        self.canvas.setSizePolicy(*size_policy)
        
        self.canvas.mpl_connect("figure_enter_event", self._focus)
                    
        self.ax = self.canvas.figure.subplots()
        if x is not None and y is not None:
            self.ax.plot(x, y)
        self.ax.set_xlabel(xlabel)
        self.ax.set_ylabel(ylabel)
        self.figure.tight_layout()
        self.figure.canvas.draw()
        
        self.layout = QVBoxLayout(self)
        self.layout.addWidget(self.canvas)
                        
        if toolbar:
            toolbar = CommandBar(parent=self)
            toolbar.setFocusPolicy(Qt.NoFocus)
            
            self.action_home = Action(FluentIcon.HOME, "Home", self)
            self.action_zoom = Action(FluentIcon.ZOOM, "Zoom", self)
            self.action_zoom.setCheckable(True)
            self.action_pan = Action(FluentIcon.MOVE, "Pan", self)
            self.action_pan.setCheckable(True)
            self.action_home.triggered.connect(self.home)
            self.action_zoom.triggered.connect(self.zoom)
            self.action_pan.triggered.connect(self.pan)

            toolbar.addAction(self.action_home)
            toolbar.addAction(self.action_zoom)
            toolbar.addAction(self.action_pan)
        
            if toolbar_left_right:
                self.page_left = Action(FluentIcon.PAGE_LEFT, "Left", self)
                self.page_right = Action(FluentIcon.PAGE_RIGHT, "Right", self)
                toolbar.addSeparator()
                toolbar.addAction(self.page_left)
                toolbar.addAction(self.page_right)
            self.layout.addWidget(toolbar)
        
            # necessary to enable panning/zoom    
            self._nav_stack = cbook._Stack()
            # listen for keyboard events
            self._key_press = self.canvas.mpl_connect("key_press_event", self._key_press_handler)
            self._id_press = self.canvas.mpl_connect("button_press_event", self._zoom_pan_handler)
            self._id_release = self.canvas.mpl_connect("button_release_event", self._zoom_pan_handler)
            self._id_drag = self.canvas.mpl_connect("motion_notify_event", self.mouse_move)
            self._pan_info = None
            self._zoom_info = None
            self._last_cursor = tools.Cursors.POINTER
            self.mode = _Mode.NONE 
            return None
        
        else:
            if toolbar_left_right:
                # Only left/right, set to align in middle
                toolbar = CommandBar(parent=self)
                toolbar.setFocusPolicy(Qt.NoFocus)
                self.page_left = Action(FluentIcon.PAGE_LEFT, "Left", self)
                self.page_right = Action(FluentIcon.PAGE_RIGHT, "Right", self)
                toolbar.addAction(self.page_left)
                toolbar.addAction(self.page_right)
                self.layout.addWidget(toolbar, alignment=Qt.AlignCenter)           
        
        #self.layout.addStretch(1)    
        self.installEventFilter(self)

    def reset_current_as_home(self):
        self._nav_stack = cbook._Stack()
        self.push_current()        

    def _key_press_handler(self, event):
        if event.key in ("left", "right"):
            try:
                dict(left=self.page_left, right=self.page_right)[event.key].trigger()
            except:
                None

    def _zoom_pan_handler(self, event):
        if self.mode == _Mode.PAN:
            if event.name == "button_press_event":
                self.press_pan(event)
            elif event.name == "button_release_event":
                self.release_pan(event)
        if self.mode == _Mode.ZOOM:
            if event.name == "button_press_event":
                self.press_zoom(event)
            elif event.name == "button_release_event":
                self.release_zoom(event)

    _PanInfo = namedtuple("_PanInfo", "button axes cid")

    def home(self, *args):
        self._nav_stack.home()
        #self.set_history_buttons()
        self._update_view()   
        self.figure.tight_layout()
        self.figure.canvas.draw()             

    def set_history_buttons(self):
        pass

    def _update_view(self):
        nav_info = self._nav_stack()
        if nav_info is None:
            return
        # Retrieve all items at once to avoid any risk of GC deleting an Axes
        # while in the middle of the loop below.
        items = list(nav_info.items())
        for ax, (view, (pos_orig, pos_active)) in items:
            ax._set_view(view)
            # Restore both the original and modified positions
            
            # TODO: for some reason, setting the pos_orig causes the figure to be smaller each time            
            #ax._set_position(pos_orig, 'original')
            ax._set_position(pos_active, 'active')
        self.canvas.draw_idle()        

    def press_pan(self, event):
        """Callback for mouse button press in pan/zoom mode."""
        if (event.button not in [MouseButton.LEFT, MouseButton.RIGHT]
                or event.x is None or event.y is None):
            return
        axes = [a for a in self.canvas.figure.get_axes()
                if a.in_axes(event) and a.get_navigate() and a.can_pan()]

        if not axes:
            return
        if self._nav_stack() is None:
            self.push_current()  # set the home button to this view
        for ax in axes:
            ax.start_pan(event.x, event.y, event.button)
        self.canvas.mpl_disconnect(self._id_drag)
        id_drag = self.canvas.mpl_connect("motion_notify_event", self.drag_pan)
        self._pan_info = self._PanInfo(
            button=event.button, axes=axes, cid=id_drag)

    def release_pan(self, event):
        """Callback for mouse button release in pan/zoom mode."""
        if self._pan_info is None:
            return
        self.canvas.mpl_disconnect(self._pan_info.cid)
        self._id_drag = self.canvas.mpl_connect(
            'motion_notify_event', self.mouse_move)
        for ax in self._pan_info.axes:
            ax.end_pan()
        self.canvas.draw_idle()
        self._pan_info = None
        self.push_current()        
    
    def drag_pan(self, event):
        """Callback for dragging in pan/zoom mode."""
        for ax in self._pan_info.axes:
            # Using the recorded button at the press is safer than the current
            # button, as multiple buttons can get pressed during motion.
            ax.drag_pan(self._pan_info.button, event.key, event.x, event.y)
        self.canvas.draw_idle()    
    
    def set_message(self, message):
        pass
    
    def mouse_move(self, event):
        self._update_cursor(event)
        #self.set_message(self._mouse_event_to_message(event))    

    def push_current(self):
        """Push the current view limits and position onto the stack."""
        self._nav_stack.push(
            WeakKeyDictionary(
                {ax: (ax._get_view(),
                      # Store both the original and modified positions.
                      (ax.get_position(True).frozen(),
                       ax.get_position().frozen()))
                 for ax in self.canvas.figure.axes}))
        self.set_history_buttons()

    def _update_cursor(self, event):
        """
        Update the cursor after a mouse move event or a tool (de)activation.
        """
        if self.mode and event.inaxes and event.inaxes.get_navigate():
            if (self.mode == _Mode.ZOOM
                    and self._last_cursor != tools.Cursors.SELECT_REGION):
                self.canvas.set_cursor(tools.Cursors.SELECT_REGION)
                self._last_cursor = tools.Cursors.SELECT_REGION
            elif (self.mode == _Mode.PAN
                  and self._last_cursor != tools.Cursors.MOVE):
                self.canvas.set_cursor(tools.Cursors.MOVE)
                self._last_cursor = tools.Cursors.MOVE
        elif self._last_cursor != tools.Cursors.POINTER:
            self.canvas.set_cursor(tools.Cursors.POINTER)
            self._last_cursor = tools.Cursors.POINTER
                
    def pan(self):
        if not self.canvas.widgetlock.available(self):
            self.set_message("pan unavailable")
            return
        if self.mode == _Mode.PAN:
            self.mode = _Mode.NONE
            self.canvas.widgetlock.release(self)
        else:
            self.mode = _Mode.PAN
            self.action_zoom.setChecked(False)
            self.canvas.widgetlock(self)
            
        for a in self.canvas.figure.get_axes():
            a.set_navigate_mode(self.mode._navigate_mode)        
    
    def zoom(self, *args):
        if not self.canvas.widgetlock.available(self):
            self.set_message("zoom unavailable")
            return
        """Toggle zoom to rect mode."""
        if self.mode == _Mode.ZOOM:
            self.mode = _Mode.NONE
            self.canvas.widgetlock.release(self)
        else:
            self.mode = _Mode.ZOOM
            self.action_pan.setChecked(False)
            self.canvas.widgetlock(self)
        for a in self.canvas.figure.get_axes():
            a.set_navigate_mode(self.mode._navigate_mode)    
    
    _ZoomInfo = namedtuple("_ZoomInfo", "direction start_xy axes cid cbar")
    
    def press_zoom(self, event):
        """Callback for mouse button press in zoom to rect mode."""
        if (event.button not in [MouseButton.LEFT, MouseButton.RIGHT]
                or event.x is None or event.y is None):
            return
        axes = [a for a in self.canvas.figure.get_axes()
                if a.in_axes(event) and a.get_navigate() and a.can_zoom()]
        if not axes:
            return
        if self._nav_stack() is None:
            self.push_current()  # set the home button to this view
        id_zoom = self.canvas.mpl_connect(
            "motion_notify_event", self.drag_zoom)
        # A colorbar is one-dimensional, so we extend the zoom rectangle out
        # to the edge of the Axes bbox in the other dimension. To do that we
        # store the orientation of the colorbar for later.
        if hasattr(axes[0], "_colorbar"):
            cbar = axes[0]._colorbar.orientation
        else:
            cbar = None
        self._zoom_info = self._ZoomInfo(
            direction="in" if event.button == 1 else "out",
            start_xy=(event.x, event.y), axes=axes, cid=id_zoom, cbar=cbar)    

    def release_zoom(self, event):
        """Callback for mouse button release in zoom to rect mode."""
        if self._zoom_info is None:
            return

        # We don't check the event button here, so that zooms can be cancelled
        # by (pressing and) releasing another mouse button.
        self.canvas.mpl_disconnect(self._zoom_info.cid)
        self.remove_rubberband()

        start_x, start_y = self._zoom_info.start_xy
        key = event.key
        # Force the key on colorbars to ignore the zoom-cancel on the
        # short-axis side
        if self._zoom_info.cbar == "horizontal":
            key = "x"
        elif self._zoom_info.cbar == "vertical":
            key = "y"
        # Ignore single clicks: 5 pixels is a threshold that allows the user to
        # "cancel" a zoom action by zooming by less than 5 pixels.
        if ((abs(event.x - start_x) < 5 and key != "y") or
                (abs(event.y - start_y) < 5 and key != "x")):
            self.canvas.draw_idle()
            self._zoom_info = None
            return

        for i, ax in enumerate(self._zoom_info.axes):
            # Detect whether this Axes is twinned with an earlier Axes in the
            # list of zoomed Axes, to avoid double zooming.
            twinx = any(ax.get_shared_x_axes().joined(ax, prev)
                        for prev in self._zoom_info.axes[:i])
            twiny = any(ax.get_shared_y_axes().joined(ax, prev)
                        for prev in self._zoom_info.axes[:i])
            ax._set_view_from_bbox(
                (start_x, start_y, event.x, event.y),
                self._zoom_info.direction, key, twinx, twiny)

        self.canvas.draw_idle()
        self._zoom_info = None
        self.push_current()

    def drag_zoom(self, event):
        """Callback for dragging in zoom mode."""
        start_xy = self._zoom_info.start_xy
        ax = self._zoom_info.axes[0]
        (x1, y1), (x2, y2) = np.clip(
            [start_xy, [event.x, event.y]], ax.bbox.min, ax.bbox.max)
        key = event.key
        # Force the key on colorbars to extend the short-axis bbox
        if self._zoom_info.cbar == "horizontal":
            key = "x"
        elif self._zoom_info.cbar == "vertical":
            key = "y"
        if key == "x":
            y1, y2 = ax.bbox.intervaly
        elif key == "y":
            x1, x2 = ax.bbox.intervalx

        self.draw_rubberband(event, x1, y1, x2, y2)
        
    def draw_rubberband(self, event, x0, y0, x1, y1):
        height = self.canvas.figure.bbox.height
        y1 = height - y1
        y0 = height - y0
        rect = [int(val) for val in (x0, y0, x1 - x0, y1 - y0)]
        self.canvas.drawRectangle(rect)

    def remove_rubberband(self):
        self.canvas.drawRectangle(None)
        
    def _focus(self, event):
        """ Set the focus of the canvas. """
        self.canvas.setFocus()    
        
    
    def resizeEvent(self, e):
        # Matplotlib wants to redraw the canvas while we are resizing, makes it all yucky
        try:
            self.resizeTimer
        except:
            
            self.resizeTimer = QTimer(self)
            self.resizeTimer.setSingleShot(True)
            self.resizeTimer.timeout.connect(self.afterResizeEvent)
        finally:
            self.resizeTimer.start(self.resize_interval)
                
        return None        

    def afterResizeEvent(self):        
        self.figure.tight_layout()
        self.figure.canvas.draw()
        try:
            self.resizeTimer.stop()
            del self.resizeTimer
        except:
            None
        return None

        
    def eventFilter(self, widget, event):
        try:
            print(f"plot widget {widget} {event.type()} {event.key()} {event.text()}")
        except:
            None
            
                    
        if event.type() == 51:
            if event.key() == Qt.Key_Left:
                try:
                    self.page_left.trigger()
                except:
                    return False
                else:
                    return True
            elif event.key() == Qt.Key_Right:
                try:
                    self.page_right.trigger()
                except:
                    return False
                else:
                    return True
                
                                

        '''
        if event.type() == QEvent.KeyPress:
            text = event.text()
            if event.modifiers():
                text = event.keyCombination().key().name.decode(encoding="utf-8")
            print(f"{event} {event.type}: {text}")
        '''            
        return False
        
    
if __name__ == '__main__':
    # enable dpi scale
    from PyQt5.QtCore import QModelIndex, Qt, QRect, QSize
    import sys
    QApplication.setHighDpiScaleFactorRoundingPolicy(Qt.HighDpiScaleFactorRoundingPolicy.PassThrough)
    QApplication.setAttribute(Qt.AA_EnableHighDpiScaling)
    QApplication.setAttribute(Qt.AA_UseHighDpiPixmaps)

    app = QApplication(sys.argv)
    w = SinglePlotWidget()
    w.resize(600, 600)
    w.show()
    app.exec_()    