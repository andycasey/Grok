
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


class JuliaReference:
    def __init__(self, value, variable_name=None):
        self.value = value
        self.__julia_variable_name__ = variable_name
        
    def __repr__(self):
        return f"{self.value}"

def ajr(o):
    """Represent objects as a Julia variable name, if it has one."""
    try:
        return o.__julia_variable_name__
    except:
        return o


class LineList(list):
    def __init__(self, julia_variable_name=None):
        super(LineList, self).__init__()
        self.__julia_variable_name__ = julia_variable_name
        return None





class BaseKorg:
    pass


class Korg(BaseKorg):
    
    def __init__(self):
        from juliacall import Main as jl
        jl.seval("using Pkg")
        jl.seval('Pkg.develop(path="/Users/andycasey/research/Grok/Korg.jl")')
        jl.seval("using Korg")
        self.jl = jl
        return None

    def format_A_X(self, default_metals_H, default_alpha_H, **kwargs):
        return self.jl.Korg.format_A_X(default_metals_H, default_alpha_H, **kwargs)
    
    def read_linelist(self, path, format="vald"):
        return self.jl.Korg.read_linelist(path, format=format)
    
    def interpolate_marcs(self, Teff, logg, A_X):
        return self.jl.Korg.interpolate_marcs(Teff, logg, A_X)
    
    def ews_to_abundances(self, atm, ll, A_X, ews, **kwargs):
        return self.jl.Korg.Fit.ews_to_abundances(atm, ll, A_X, ews, **kwargs)
    

    def ews_to_stellar_parameters(self, ll, EWs, Teff0=5000.0, logg0=3.5, vmic0=1.0, metallicity0=0.0, **kwargs):
        return self.jl.Korg.Fit.ews_to_stellar_parameters(
            ll, 
            EWs, 
            Teff0,
            logg0,
            vmic0,
            metallicity0,
            **kwargs
        )

    


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
        #self.setProcessChannelMode(QProcess.MergedChannels)
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
        #self.precompile() 
        return None
    
    def precompile(self):
        self.async_eval("precompile()", lambda *_: None)

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
        print(f"eval: {command}")
        self.write((f"{command}\n").encode("utf8"))
        #self.waitForBytesWritten()
        self.waitForReadyRead() # critically important to have this
        while True:
            try:
                print(f"Waiting on {command}. Was {N} lines, now {len(self.stdout)} lines")
                for line in self.stdout[N:]:
                    print(f" -> {line}")
                    if line.startswith(f"{command}>>>"):
                        _command, response = line.split(">>>", 1)
                        assert _command == command
                        return (response, variable_name)
            except:
                self.waitForReadyRead()

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
    
    
    def read_linelist(self, path, format):
        response, variable_name = self.eval(f"Korg.read_linelist(\"{path}\", format=\"{format}\")")
        return _parse_linelist(response, variable_name)

    def async_read_linelist(self, path, format, callback=None):
        def outer_callback(response, variable_name):
            return callback(_parse_linelist(response, variable_name))
        self.async_eval(f"Korg.read_linelist(\"{path}\", format=\"{format}\")", outer_callback)
        return None

    def format_A_X(self, metals_H, alpha_H):
        response, _ = self.eval(f"Korg.format_A_X({metals_H}, {alpha_H})")
        return list(map(float, response.strip()[1:-1].split(",")))
        
    '''
    def interpolate_marcs(self, Teff, logg, A_X):
        def callback(response, variable_name):
            print(f"got a response: {response}")
        self.async_eval(f"Korg.interpolate_marcs({Teff}, {logg}, {A_X})", callback=callback)
        #print(f"marcs: {response}")
        #return variable_name
    '''
    
    def interpolate_marcs(self, Teff, logg):
        response, variable_name = self.eval(f"Korg.interpolate_marcs({Teff}, {logg})")
        return JuliaReference(response, variable_name)
    
    
    def ews_to_abundances(self, atm, ll, A_X, ews):
        response, _ = self.eval(f"Korg.ews_to_abundances({ajr(atm)}, {ajr(ll)}, {ajr(A_X)}, {ajr(ews)})")
        return response
    

        

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
    

