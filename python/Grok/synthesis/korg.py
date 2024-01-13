
import os
import subprocess
import threading
import queue
import logging
import sys
import asyncio
import re
from random import choices
from string import ascii_lowercase
from time import sleep

try:
    from Qt.QtCore import QProcess
except:
    logging.exception(f"Could not import Qt.QtCore.QProcess: QKorgProcess will not work.")
    QProcess = object

from dataclasses import dataclass
from typing import Union, Sequence

@dataclass(frozen=True)
class Species:
    formula: str
    charge: int
    
    def __repr__(self):
        return f"{self.formula} {'I' * self.charge}"

@dataclass(frozen=True)
class Line:
    wl: float
    log_gf: float
    species: Species
    E_lower: float
    gamma_rad: float
    gamma_stark: float
    vdW: Union[float, Sequence[float]]

    def __repr__(self):
        return f"{self.__class__.__name__}({round(self.wl*1e8, 6)} Å, log gf = {round(self.log_gf, 2)})"

class LineList(list):
    def __init__(self, julia_variable_name=None):
        super(LineList, self).__init__()
        self.julia_variable_name = julia_variable_name
        return None



class BaseKorg:
    pass


class Korg(BaseKorg):
    
    def __init__(self):
        from juliacall import Main as jl
        jl.seval("using Korg")
        self.jl = jl
        return None

LINE_PATTERN = re.compile(
    "\((?P<wl>[-\d\.]+(e-?\d+)?), (?P<log_gf>[-\d\.]+(e-?\d+)?), (?P<species_repr>\w+\sI+), (?P<E_lower>[-\d\.]+(e[-+]\d+)?), (?P<gamma_rad>[-\d\.]+(e-?\d+)?), (?P<gamma_stark>[-\d\.]+(e-?\d+)?), (?P<vdW>[-\d\.]+(e-?\d+)?)\)"
)

class BaseKorgProcess(BaseKorg):
    pass
    
    

class QKorgProcess(QProcess, BaseKorgProcess):
        
    def __init__(self, **kwargs):
        super(QKorgProcess, self).__init__(**kwargs)
        # get current file path
        path = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            "korg_process.jl"
        )
        self.finished.connect(self._handle_finished)
        self.stateChanged.connect(self._handle_state_changed)
        self.readyReadStandardOutput.connect(self._handle_stdout)
        self.readyReadStandardError.connect(self._handle_stderr)
        # TODO: Get the julia executable path from the application config,
        #       or from the juliacall bridge
        self.start("julia", [path])
        self.variable_names = set()
        self.stdout = []
        self.stderr = []
        self.callbacks = {}
        return None    

    def _handle_finished(self):
        print(f"{self} finished")
        
    def _handle_state_changed(self, state):
        print(f"{self} state changed: {state}")
        
    def _handle_stderr(self):
        data = self.readAllStandardError()
        stderr = bytes(data).decode("utf8")
        self.stderr.append(stderr)
        print(f"stderr: {stderr}")
    
    def _handle_stdout(self):
        data = self.readAllStandardOutput()
        stdout = bytes(data).decode("utf8")
        print(f"stdout: {stdout}")        
        self.stdout.append(stdout)
        try:
            command, response = self.stdout[-1].split(">>>", 1)
        except:
            None
        else:
            r = self.callbacks.pop(command, None)
            if r is not None:
                callback, variable_name = r
                callback(response, variable_name)
        
    
    def assign_variable_name(self, variable_name=None, prefix="", k=12, **kwargs):
        if variable_name is not None:
            return variable_name
        while True:
            variable_name = "".join(choices(ascii_lowercase, k=12))
            if variable_name not in self.variable_names:
                self.variable_names.add(variable_name)
                return variable_name


    def eval(self, command, variable_name=None, **kwargs):
        N = len(self.stdout)
        if variable_name is None:
            variable_name = self.assign_variable_name(**kwargs)
        command = f"{variable_name} = {command}"
        self.write((f"{command}\n").encode("utf8"))
        self.waitForReadyRead() # critically important to have this
        while True:
            try:
                for line in self.stdout[N:]:
                    if line.startswith(f"{command}>>>"):
                        command, response = self.stdout[-1].split(">>>", 1)
                        return (response, variable_name)
            except:
                continue
    
    
    def async_eval(self, command, callback, variable_name=None, **kwargs):
        if variable_name is None:
            variable_name = self.assign_variable_name(**kwargs)
        
        self.callbacks[command] = (callback, variable_name)
        self.write((f"{command}\n").encode("utf8"))
        return None
    
    def air_to_vacuum(self, v) -> float:
        response, _ = self.eval(f"Korg.air_to_vacuum({v})")
        return float(response)
    
    
    def read_linelist(self, path, format="vald"):
        response, variable_name = self.eval(f"Korg.read_linelist(\"{path}\", format=\"{format}\")")
        return _parse_linelist(response, variable_name)

    def async_read_linelist(self, path, format="vald", callback=None):
        def outer_callback(response, variable_name):
            return callback(_parse_linelist(response, variable_name))
        self.async_eval(f"Korg.read_linelist(\"{path}\", format=\"{format}\")", outer_callback)
        return None


def _parse_linelist(response, variable_name):
    linelist = LineList(julia_variable_name=variable_name)
    for match in re.finditer(LINE_PATTERN, response):
        linelist.append(
            Line(
                wl=float(match.group("wl")),
                log_gf=float(match.group("log_gf")),
                species=Species(
                    formula=match.group("species_repr").split()[0],
                    charge=match.group("species_repr").count("I"),
                ),
                E_lower=float(match.group("E_lower")),
                gamma_rad=float(match.group("gamma_rad")),
                gamma_stark=float(match.group("gamma_stark")),
                vdW=float(match.group("vdW")),
            )
        )
    return linelist            
    



def enqueue_stdout(out, queue, callback):
    for line in iter(out.readline, b''):
        queue.put(line)
        if len(line) > 0:
            callback()
    out.close()
    
def enqueue_stderr(out, queue, callback):
    for line in iter(out.readline, b''):
        queue.put(line)
        if len(line) > 0:
            callback()
    out.close()    





class KorgProcess(BaseKorgProcess):
    
    def __init__(self):
        path = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            "korg_process.jl"
        )
        
        self.process = subprocess.Popen(
            ["julia", path],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            stdin=subprocess.PIPE,
            bufsize=1,
            encoding="utf-8",
        )
        # This non-blocking doesn't seem to work. https://stackoverflow.com/questions/375427/a-non-blocking-read-on-a-subprocess-pipe-in-python/4896288#4896288
        #os.set_blocking(self.process.stdout.fileno(), False)
        self.stdout_queue = queue.Queue()
        self.stdout_thread = threading.Thread(target=enqueue_stdout, args=(self.process.stdout, self.stdout_queue, self.handle_stdout))
        self.stdout_thread.daemon = True
        self.stdout_thread.start()

        self.stderr_queue = queue.Queue()

        self.stderr_thread = threading.Thread(target=enqueue_stderr, args=(self.process.stderr, self.stderr_queue, self.handle_stderr))
        self.stderr_thread.daemon = True
        self.stderr_thread.start()
        
        self.callbacks = {}
        self.variable_names = set()
        self.stdout = []
        self.stderr = []
        self.ready = False
        return None
    
    def assign_variable_name(self, variable_name=None, prefix="", k=12, **kwargs):
        if variable_name is not None:
            return variable_name
        while True:
            variable_name = "".join(choices(ascii_lowercase, k=12))
            if variable_name not in self.variable_names:
                self.variable_names.add(variable_name)
                return variable_name
        
            
        
    
    def precompile(self):
        return self.eval(f"precompile()")

    def handle_stdout(self):
        stdout = self.stdout_queue.get()
        self.stdout.append(stdout)
        try:
            command, response = self.stdout[-1].split(">>>", 1)
        except:
            None
        else:
            callback = self.callbacks.pop(command, None)
            if callback is not None:
                callback(response)


    def handle_stderr(self):
        stderr = self.stdout_queue.get()
        self.stderr.append(stderr)
        print(f"stderr: {stderr}")

    def send(self, msg):
        self.ready = False
        self.process.stdin.write(msg)
        self.process.stdin.flush()
        return None

    def eval(self, command, timeout=30):
        N = len(self.stdout)
        self.send(f"{command}\n")
        while True:
            try:
                for line in self.stdout[N:]:
                    if line.startswith(f"{command}>>>"):
                        command, response = self.stdout[-1].split(">>>", 1)
                        return response
            except:
                continue
    
    
    def async_eval(self, command, callback):
        self.callbacks[command] = callback
        self.send(f"{command}\n")
        return None
    
    def air_to_vacuum(self, v):
        return float(self.eval(f"Korg.air_to_vacuum({v})"))
    
    def read_linelist(self, path, format="vald", **kwargs):        
        variable_name = self.assign_variable_name(**kwargs)
        command = f"{variable_name} = Korg.read_linelist(\"{path}\", format=\"{format}\")"
        
        linelist = LineList(julia_variable_name=variable_name)
        for match in re.finditer(self.LINE_PATTERN, self.eval(command)):
            linelist.append(
                Line(
                    wl=float(match.group("wl")),
                    log_gf=float(match.group("log_gf")),
                    species=Species(
                        formula=match.group("species_repr").split()[0],
                        charge=match.group("species_repr").count("I"),
                    ),
                    E_lower=float(match.group("E_lower")),
                    gamma_rad=float(match.group("gamma_rad")),
                    gamma_stark=float(match.group("gamma_stark")),
                    vdW=float(match.group("vdW")),
                )
            )
        return linelist
    
    def format_A_X(self, metals_H, alpha_H, **kwargs):
        variable_name = self.assign_variable_name(**kwargs)
        command = f"{variable_name} = Korg.format_A_X({metals_H}, {alpha_H})"
        return self.eval(command)

    def interpolate_marcs(self, Teff, logg, A_X, **kwargs):
        variable_name = self.assign_variable_name(**kwargs)
        