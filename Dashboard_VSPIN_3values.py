# -*- coding: utf-8 -*-
"""
Created on Mon Sep  1 15:11:40 2025

@author: Mikel
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import time
import threading
import numpy as np
import random

import matplotlib
matplotlib.use("QtAgg")  # o "TkAgg" si no tienes Qt
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Circle

from labjack import ljm
from opcua import Client as OPCClient

# -------------------- CONFIG --------------------
SIMULAR_RESULTADOS = False
DEBUG = False
IP_T4 = "192.168.106.141"       # LabJack T4
REFRESH_HZ = 100                # refresco UI (Hz)

# LJTCS: 5.9Ω y ganancia x20
R_SHUNT = 5.9
GAIN = 20.0
V_TO_MA = 1000.0 / (R_SHUNT * GAIN)  # ≈ 8.474576271
AIN_RANGE = 2.5

# MAPEO DE VARIABLES
# AIN4 -> RPM
# AIN5 -> DIR
# AIN6 -> MOD
AIN_NAMES = ["AIN4", "AIN5", "AIN6"]

# Escalado radial usando AngleValVSpin
# AngleValVSpin = 0.25 * mod_mA - 1.0
# Si mod_mA va de 4 a 20 mA -> AngleValVSpin va de 0 a 4
ANGLE_VAL_MIN = 0.0
ANGLE_VAL_MAX = 4.0
R_MIN = 0.15
R_MAX = 1.00

# Reconexión
LJ_RETRY_SEC  = 2.0
OPC_RETRY_SEC = 2.0
QUIET_RETRIES = True
PRINT_EVERY_N_FAILS = 10

# OPC UA
opc_url  = "opc.tcp://192.168.106.151:4840"
opc_user = "OpcUaOperator"
opc_pass = "kuka"

# ----------------- UTILIDADES -------------------
def clamp(x, lo, hi):
    return lo if x < lo else hi if x > hi else x

def map_linear(x, a0, a1, b0, b1):
    if a1 == a0:
        return b0
    t = (x - a0) / (a1 - a0)
    return b0 + t * (b1 - b0)

def volts_to_mA(v, clip=False):
    mA = float(v) * V_TO_MA
    if clip:
        mA = clamp(mA, 4.0, 20.0)
    return mA

def compute_vspin_values(v_rpm, v_dir, v_mod, simulate=False):
    if simulate:
        rpm_mA = random.uniform(4.0, 20.0)
        dir_mA = random.uniform(4.0, 20.0)
        mod_mA = random.uniform(4.0, 20.0)
    else:
        rpm_mA = volts_to_mA(v_rpm, clip=True)
        dir_mA = volts_to_mA(v_dir, clip=True)
        mod_mA = volts_to_mA(v_mod, clip=True)

    rpm = max(0.0, 7500.0 * rpm_mA - 30000.0)
    angle_dir = (22.5 * dir_mA - 90.0 + 360.0) % 360.0
    angle_val = 0.25 * mod_mA - 1.0

    return rpm, angle_dir, angle_val, rpm_mA, dir_mA, mod_mA

# ----------------- LABJACK I/O ------------------
class LabJackReader:
    def __init__(self, ip):
        self.ip = ip
        self.handle = None
        self.connected = False
        self.lock = threading.Lock()
        self._last_attempt = 0.0
        self._fails = 0
        self._last_state = None

    def _log_state_change(self, now_connected):
        if self._last_state is None or self._last_state != now_connected:
            self._last_state = now_connected
            if now_connected:
                print("LabJack conectado")
            else:
                print("LabJack desconectado")

    def open(self):
        try:
            self.handle = ljm.openS("T4", "ETHERNET", self.ip)
            info = ljm.getHandleInfo(self.handle)
            if DEBUG:
                print(f"LJ conectado: T{info[0]} via {info[1]} (serial {info[2]})")
            self.configure()
            self.connected = True
            self._fails = 0
            self._log_state_change(True)
        except Exception as e:
            self.connected = False
            self._log_state_change(False)
            self._fails += 1
            if not QUIET_RETRIES and (self._fails % PRINT_EVERY_N_FAILS == 0):
                print(f"No se pudo conectar al LabJack ({e})")

    def ensure_connected(self):
        now = time.monotonic()
        if self.connected:
            return True
        if now - self._last_attempt < LJ_RETRY_SEC:
            return False
        self._last_attempt = now
        self.open()
        return self.connected

    def configure(self):
        # Habilitar AIN4, AIN5, AIN6 como analógicos
        mask = (1 << 4) | (1 << 5) | (1 << 6)
        ljm.eWriteName(self.handle, "DIO_ANALOG_ENABLE", mask)

        for ch in (4, 5, 6):
            ljm.eWriteName(self.handle, f"AIN{ch}_RANGE", AIN_RANGE)
            ljm.eWriteName(self.handle, f"AIN{ch}_RESOLUTION_INDEX", 2)

    def read_triplet(self):
        with self.lock:
            if not self.ensure_connected():
                raise RuntimeError("LabJack no conectado")
            try:
                v_rpm, v_dir, v_mod = ljm.eReadNames(self.handle, 3, AIN_NAMES)
                return v_rpm, v_dir, v_mod
            except Exception as e:
                try:
                    if self.handle is not None:
                        ljm.close(self.handle)
                except Exception:
                    pass
                self.handle = None
                self.connected = False
                self._log_state_change(False)
                if DEBUG:
                    print(f"Lectura LJ fallida: {e}")
                raise

    def close(self):
        with self.lock:
            try:
                if self.handle is not None:
                    ljm.close(self.handle)
            finally:
                self.handle = None
                self.connected = False

# ----------------- OPC UA -----------------------
class OPCConnector:
    def __init__(self, url, user, password):
        self.url = url
        self.user = user
        self.password = password
        self.client = None
        self.connected = False
        self._last_attempt = 0.0
        self.lock = threading.Lock()
        self._fails = 0
        self._last_state = None

    def _log_state_change(self, now_connected):
        if self._last_state is None or self._last_state != now_connected:
            self._last_state = now_connected
            print("OPC UA conectado" if now_connected else "OPC UA desconectado")

    def open(self):
        with self.lock:
            try:
                self.client = OPCClient(self.url, timeout=2.0)
                self.client.set_user(self.user)
                self.client.set_password(self.password)
                self.client.connect()
                self.connected = True
                self._fails = 0
                self._log_state_change(True)
            except Exception as e:
                self.connected = False
                self._log_state_change(False)
                self._fails += 1
                if not QUIET_RETRIES and (self._fails % PRINT_EVERY_N_FAILS == 0):
                    print(f"No se pudo conectar OPC UA ({e})")
                try:
                    if self.client:
                        self.client.disconnect()
                except Exception:
                    pass
                self.client = None

    def ensure_connected(self):
        now = time.monotonic()
        if self.connected:
            return True
        if now - self._last_attempt < OPC_RETRY_SEC:
            return False
        self._last_attempt = now
        self.open()
        return self.connected

    def close(self):
        with self.lock:
            try:
                if self.client is not None:
                    self.client.disconnect()
            finally:
                self.client = None
                self.connected = False

# ----------------- UI POLAR ---------------------
class PolarUI:
    def __init__(self, lj_reader: LabJackReader, opc: OPCConnector):
        self.lj = lj_reader
        self.opc = opc
        self.dt = 1.0 / max(1, REFRESH_HZ)

        # Ventana
        self.fig = plt.figure(figsize=(12, 9), num="Visualizador Vspin")
        try:
            self.fig.canvas.manager.set_window_title("Vspin")
        except Exception:
            pass

        # --------- CABECERA ----------
        self.fig.text(0.03, 0.975, "Labjack T4:", fontsize=12, ha="left", va="center", zorder=5)
        self.txt_lj = self.fig.text(0.28, 0.975, "Desconectado", fontsize=12, ha="right", va="center", zorder=5)
        self.ax_lj = self.fig.add_axes([0.30, 0.958, 0.035, 0.035])
        self.ax_lj.set_axis_off()
        self.ax_lj.set_aspect('equal')
        self.circ_lj = Circle((0.5, 0.5), 0.38, transform=self.ax_lj.transAxes, color="red")
        self.ax_lj.add_patch(self.circ_lj)

        self.fig.text(0.03, 0.945, "Device Connector:", fontsize=12, ha="left", va="center", zorder=5)
        self.txt_opc = self.fig.text(0.28, 0.945, "Desconectado", fontsize=12, ha="right", va="center", zorder=5)
        self.ax_opc = self.fig.add_axes([0.30, 0.928, 0.035, 0.035])
        self.ax_opc.set_axis_off()
        self.ax_opc.set_aspect('equal')
        self.circ_opc = Circle((0.5, 0.5), 0.38, transform=self.ax_opc.transAxes, color="red")
        self.ax_opc.add_patch(self.circ_opc)

        # --------- GRÁFICO POLAR ----------
        self.ax_polar = self.fig.add_axes([0.07, 0.06, 0.86, 0.80], projection="polar")
        self.ax_polar.set_title("COMPENSACION RADIAL (en tiempo real)", pad=6)
        self.ax_polar.set_rlim(0, 1.1)
        self.ax_polar.set_rticks([])
        self.ax_polar.set_theta_zero_location("E")
        self.ax_polar.set_theta_direction(1)

        for ang_deg in (0, 90, 180, 270):
            ang = np.deg2rad(ang_deg)
            self.ax_polar.plot([ang, ang], [0, 1.05], ls="--", lw=1, alpha=0.5)

        theta_full = np.linspace(0, 2 * np.pi, 361)
        self.ax_polar.plot(theta_full, np.full_like(theta_full, R_MIN), ls=":", lw=1, alpha=0.4)

        self.vector, = self.ax_polar.plot([0, 0], [0, 1.0], lw=3, color="blue")
        self.head, = self.ax_polar.plot([0], [0], marker="o", ms=6, color="blue")

        # Texto debug y RPM en dos líneas separadas
        self.text_debug = self.ax_polar.text(
            0.62, 0.05, "", transform=self.ax_polar.transAxes,
            ha="left", va="bottom", fontsize=11
        )
        self.text_rpm = self.ax_polar.text(
            0.62, 0.015, "", transform=self.ax_polar.transAxes,
            ha="left", va="bottom", fontsize=11, fontweight="bold"
        )

        self._last_opc_check = 0.0
        self.ani = FuncAnimation(
            self.fig, self.update, interval=int(self.dt * 1000),
            blit=False, cache_frame_data=False
        )

    def _set_lj_status(self, ok: bool):
        self.circ_lj.set_color("green" if ok else "red")
        self.txt_lj.set_text("Conectado" if ok else "Desconectado")

    def _set_opc_status(self, ok: bool):
        self.circ_opc.set_color("green" if ok else "red")
        self.txt_opc.set_text("Conectado" if ok else "Desconectado")

    def update(self, _frame):
        try:
            v_rpm, v_dir, v_mod = self.lj.read_triplet()

            rpm, angle_dir_deg, angle_val, rpm_mA, dir_mA, mod_mA = compute_vspin_values(
                v_rpm, v_dir, v_mod, SIMULAR_RESULTADOS
            )

            self._set_lj_status(True)

            theta = np.deg2rad(angle_dir_deg)

            # El radio ahora sale de AngleValVSpin, no directamente de mod_mA
            angle_val_clamped = clamp(angle_val, ANGLE_VAL_MIN, ANGLE_VAL_MAX)
            r = map_linear(angle_val_clamped, ANGLE_VAL_MIN, ANGLE_VAL_MAX, R_MIN, R_MAX)

            self.vector.set_data([theta, theta], [0, r])
            self.head.set_data([theta], [r])

            self.text_debug.set_text(
                f"AngleDir={angle_dir_deg:.2f}°, AngleVal={angle_val:.2f}, "
                f"mod={mod_mA:.2f} mA, r={r:.2f}"
            )
            self.text_rpm.set_text(f"RPM: {rpm:.1f}")

            if DEBUG:
                print(
                    f"AngleDir={angle_dir_deg:.2f}°, AngleVal={angle_val:.2f}, "
                    f"mod={mod_mA:.2f} mA, r={r:.2f}, RPM={rpm:.1f}"
                )

        except Exception:
            self._set_lj_status(False)
            self.vector.set_data([np.nan, np.nan], [np.nan, np.nan])
            self.head.set_data([np.nan], [np.nan])
            self.text_debug.set_text("—")
            self.text_rpm.set_text("RPM: —")

        now = time.monotonic()
        if now - self._last_opc_check >= OPC_RETRY_SEC:
            self._last_opc_check = now
            ok = self.opc.ensure_connected()
            self._set_opc_status(bool(ok))

        return self.vector,

    def show(self):
        try:
            plt.show(block=True)
        except TypeError:
            plt.show()

# ----------------- MAIN -------------------------
def main():
    lj = LabJackReader(IP_T4)
    opc = OPCConnector(opc_url, opc_user, opc_pass)
    ui = PolarUI(lj, opc)

    try:
        ui.show()
    except KeyboardInterrupt:
        print("\nInterrumpido por el usuario.")
    finally:
        lj.close()
        opc.close()
        print("Cierres realizados correctamente.")

if __name__ == "__main__":
    main()
