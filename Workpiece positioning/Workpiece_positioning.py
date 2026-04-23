# -*- coding: utf-8 -*-
"""
Created on Fri Sep  5 10:45:35 2025

@author: CAMacabado
"""

import tkinter as tk
from tkinter import ttk
import threading
import time
from threading import Lock

import numpy as np
import math

import os, re
from tkinter import filedialog

# import config

# ======================= CONFIGURACIÓN ========================================
# Montaje Vspin
longitud_pinza = 111.6  # [mm]
longitud_hta = 25.6     # [mm]
longitud_contact = 3.3  # [mm]
longitud_pivotamiento = longitud_pinza + longitud_hta - longitud_contact

sentido_horario = True   # direccion de la trayectoria

# LABJACK T4
IP_T4 = "192.168.106.141"        # IP fija del LABJACK T4
analogicas = [4, 5]              # canales analógicos habilitados (FIO/EIO)
names_lj = [f"AIN{ch}" for ch in analogicas]

# KUKA DEVICE CONNECTOR (OPC UA)
opc_url = "opc.tcp://192.168.106.151:4840"
opc_user = "OpcUaOperator"
opc_pass = "kuka"
# ............................................................................. 
'Variables OPC UA'
variables = {
    "X": "ns=8;s=ns=7%3Bi=5004??krlvar:/System/R1#$POS_ACT.CartesianCoordinates.X",
    "Y": "ns=8;s=ns=7%3Bi=5004??krlvar:/System/R1#$POS_ACT.CartesianCoordinates.Y",
    "A": "ns=8;s=ns=7%3Bi=5004??krlvar:/System/R1#$POS_ACT.Orientation.A",
}

DEBUG = False

# ------------------  FILTRO DE CAPTURA  --------------------------------------
MIN_ANGLEVAL_FOR_POINT = 0.20   # solo se graba si AngleValVSpin >= este valor
# MIN_ANGLEVAL_FOR_POINT = 0.0   # DEMO


# ======================= CONFIG REFERENCIAS ==================================
# Distancia máxima entre puntos consecutivos al mostrar la trayectoria de referencia
REF_TOL_MM = 0.10                       # en [mm]
# Rotación antihoraria aplicada a los puntos de referencia cargados
REF_ROT_DEG     = 0.0                   # grados [º]
REF_ROT_ORIGIN  = (218.711, 2.994)      # (x0, y0) del centro de rotación
# Desplazamiento adicional
REF_SHIFT_X     = 0.0                   # en [mm]
REF_SHIFT_Y     = 0.0                   # en [mm]

# -------------- Parámetros del BEST-FIT --------------------------------------
BESTFIT_MAX_ROT_DEG   = 10.0   # |θ| ≤ 10°
BESTFIT_TRIM          = 0.80   # usa el 80% de inliers con menor error
BESTFIT_COARSE_STEP   = 0.5    # paso del barrido angular (°)
BESTFIT_ICP_ITERS     = 35
BESTFIT_ICP_TOL       = 1e-6
BESTFIT_NN_GATE_MM    = 80.0   # puerta de distancia para NN (mm)
BESTFIT_MAX_SAMPLES   = 4000   # límite de puntos para el barrido (rendimiento)
BESTFIT_ICP_SAMPLES   = 3000   # muestreo para ICP (rendimiento)
# -----------------------------------------------------------------------------

# =================== Dependencias hardware  ==================================
try:
    from opcua import Client as OPCClient
except Exception as e:
    OPCClient = None
    print("[AVISO] No se pudo importar opcua.Client:", e)

try:
    from labjack import ljm
except Exception as e:
    ljm = None
    print("[AVISO] No se pudo importar labjack.ljm:", e)


# ======================= UTILIDADES ==========================================
def configurar_canales_analogicos(handle, canales):
    """
    Configura como analógicos los FIO/EIO indicados en 'canales'
    usando el registro DIO_ANALOG_ENABLE.
    """
    if ljm is None:
        raise RuntimeError("LJM no disponible.")
    mascara = sum(1 << ch for ch in canales)
    ljm.eWriteName(handle, "DIO_ANALOG_ENABLE", mascara)
    print(f"[LJ] Configurados como analógicos: {canales}, máscara = {mascara} (decimal)")

def volts_to_mA(v, clip=True):
    # LJTCS: 5.9Ω y ganancia x20
    R_SHUNT = 5.9
    GAIN = 20.0
    mA = v * 1000.0 / (R_SHUNT * GAIN)
    if clip:
        mA = max(4.0, min(20.0, mA))
    return mA

# def calcular_pto_contacto(datos, L):
#     modulo = datos["AngleValVSpin"]
#     direccion = datos["AngleDirVSpin"]
#     x = datos["X"]
#     y = datos["Y"]
#     A = datos["A"]
#     angulo = A - 90 + direccion
#     x_contacto = x + L * np.sin(np.radians(modulo)) * np.cos(np.radians(angulo))
#     y_contacto = y + L * np.sin(np.radians(modulo)) * np.sin(np.radians(angulo))
#     return np.array([x_contacto, y_contacto])

def normalizar_angulo_deg(a):
    return np.mod(a, 360.0)

def getAngle(inicio, final):
    """
    Angulo de trayectoria en grados [0,360)
    """
    ang = math.degrees(math.atan2(final[1] - inicio[1], final[0] - inicio[0]))
    return ang + 360 if ang < 0 else ang

def rebaba_normal_mm(L, ang_cabezal_deg, ang_trayect_deg,
                     dir_vspin_deg, mod_vspin_deg):
    """
    MISMO cálculo que calcular_rebaba, pero devolviendo rebaba con signo opcional.
    Si quieres EXACTAMENTE lo del script viejo, luego hacemos abs().
    """
    ang_cabezal = normalizar_angulo_deg(ang_cabezal_deg)
    A_vspin = dir_vspin_deg + ang_cabezal - 90.0
    A_vspin = normalizar_angulo_deg(A_vspin)

    dif = A_vspin - ang_trayect_deg
    alpha = mod_vspin_deg * np.sin(np.radians(dif))

    rebaba = L * np.tan(np.radians(alpha))
    rebaba = np.nan_to_num(rebaba, nan=0.0)
    return rebaba  # puede ser +/-


def calcular_pto_contacto_normal(datos, L, angulo_trayectoria_deg, sentido_horario):
    """
    Replica el esquema del script de palpado:
    1) rebaba escalar [mm] usando rebaba_normal_mm (== calcular_rebaba)
    2) desplazar X,Y a lo largo de la normal geométrica.
    """
    # --- Lectura de datos ---
    x = float(datos["X"])
    y = float(datos["Y"])
    A = float(datos["A"])

    dir_vspin = float(datos["AngleDirVSpin"])
    mod_vspin = float(datos["AngleValVSpin"])

    # 1) Rebaba en mm (como en el script viejo)
    rebaba = rebaba_normal_mm(
        L,
        ang_cabezal_deg=A,
        ang_trayect_deg=angulo_trayectoria_deg,
        dir_vspin_deg=dir_vspin,
        mod_vspin_deg=mod_vspin
    )

    # El script viejo hace abs() y luego juega con el signo global
    # (direccion_antihoraria). Aquí asumimos rebaba siempre >= 0:
    d = abs(rebaba)

    # 2) Tangente de trayectoria (de pos_previa -> pos_actual)
    th = np.radians(angulo_trayectoria_deg)
    tx = np.cos(th)
    ty = np.sin(th)

    # 3) Normales (izquierda / derecha)
    n_izq = np.array([-ty,  tx])
    n_der = np.array([ ty, -tx])

    # 4) Normal "hacia fuera" según sentido de recorrido
    #    - sentido_horario = True  -> exterior = izquierda
    #    - sentido_horario = False -> exterior = derecha
    if sentido_horario:
        n_out = n_izq
    else:
        n_out = n_der

    # 5) Punto de contacto
    x_contacto = x + d * n_out[0]
    y_contacto = y + d * n_out[1]

    return np.array([x_contacto, y_contacto])

def fit_circle_kasa(pts):
    """
    Ajuste de círculo (método de Kasa) en 2D.
    pts: iterable de (x_mm, y_mm)
    Devuelve: (cx, cy, r)
    """
    P = np.asarray(pts, dtype=float)
    if P.shape[0] < 3:
        raise ValueError("Se necesitan al menos 3 puntos para ajustar un círculo.")
    x = P[:, 0]
    y = P[:, 1]
    A = np.c_[2*x, 2*y, np.ones_like(x)]
    b = x**2 + y**2
    sol, *_ = np.linalg.lstsq(A, b, rcond=None)
    cx, cy, c0 = sol
    r2 = c0 + cx*cx + cy*cy
    r = math.sqrt(max(r2, 0.0))
    return float(cx), float(cy), float(r)

def fit_line_tls(pts):
    """
    Ajuste de línea por mínimos cuadrados ortogonales.
    Devuelve (x0, y0, vx, vy) donde (x0,y0) es el centroide y (vx,vy) dirección unitaria.
    """
    P = np.asarray(pts, dtype=float)
    if P.shape[0] < 2:
        raise ValueError("Se necesitan al menos 2 puntos para ajustar una línea.")
    mu = P.mean(axis=0)
    Q = P - mu
    # vector propio principal (SVD)
    _, _, vh = np.linalg.svd(Q, full_matrices=False)
    vx, vy = vh[0, 0], vh[0, 1]
    nrm = math.hypot(vx, vy) or 1.0
    return float(mu[0]), float(mu[1]), float(vx/nrm), float(vy/nrm)

# ....................... Cargar referencia ...................................
def get_piece_zone(file_path):
    return os.path.splitext(os.path.basename(file_path))[0]

def read_commands_from_file(file_path):
    if not file_path:
        return []
    with open(file_path, "r", encoding="utf-8", errors="ignore") as f:
        return [line.strip() for line in f.readlines()]

def parse_line(line, cmd_type):
    # X 123.45, Y -67.8    (tolerante a espacios)
    pattern = r'X\s*(-?\d+\.?\d*)\s*,\s*Y\s*(-?\d+\.?\d*)'
    if cmd_type == 'LIN':
        m = re.search(pattern, line)
        return (float(m.group(1)), float(m.group(2))) if m else None
    elif cmd_type == 'CIRC':
        ms = re.findall(pattern, line)
        return [(float(x), float(y)) for x, y in ms] if ms else []

def parse_LIN(line, _cmd_type):
    pattern = r'X\s*(-?\d+\.?\d*)\s*,\s*Y\s*(-?\d+\.?\d*)\s*,\s*Z\s*(-?\d+\.?\d*)\s*,\s*A\s*(-?\d+\.?\d*)\s*,\s*B\s*(-?\d+\.?\d*)\s*,\s*C\s*(-?\d+\.?\d*)'
    m = re.search(pattern, line)
    if not m:
        # fallback: al menos X/Y
        xy = parse_line(line, 'LIN')
        if xy:
            return {'X': xy[0], 'Y': xy[1], 'Z': 0.0, 'A': 0.0, 'B': 0.0, 'C': 0.0}
        raise ValueError("No se pudieron parsear X,Y,Z,A,B,C de una LIN")
    x, y, z, a, b, c = (float(m.group(i)) for i in range(1, 7))
    return {'X': x, 'Y': y, 'Z': z, 'A': a, 'B': b, 'C': c}

def calculate_circle_center(p1, p2, p3):
    A = np.array([[p1[0]-p2[0], p1[1]-p2[1]],
                  [p1[0]-p3[0], p1[1]-p3[1]]])
    B = np.array([(p1[0]**2 - p2[0]**2 + p1[1]**2 - p2[1]**2)/2.0,
                  (p1[0]**2 - p3[0]**2 + p1[1]**2 - p3[1]**2)/2.0])
    center = np.linalg.solve(A, B)
    return center

def calculate_angle(point, center):
    return np.arctan2(point[1] - center[1], point[0] - center[0])

def arc_length(radius, start_angle, end_angle, cambio_arco, sentido):
    if start_angle * end_angle < 0 and end_angle < 0 and sentido == 'Antihorario':
        end_angle += 2*np.pi
    if cambio_arco and start_angle > end_angle:
        end_angle += 2*np.pi
    length = radius * (end_angle - start_angle)
    return abs(length), end_angle

def calcular_recta(p1, p2, desired_distance):
    d = np.hypot(p2[0]-p1[0], p2[1]-p1[1])
    if d <= desired_distance:
        return [p1, p2]
    n = int(np.ceil(d / desired_distance))
    return [(p1[0] + j*(p2[0]-p1[0])/n, p1[1] + j*(p2[1]-p1[1])/n) for j in range(n+1)]

def sentido_recorrido(a0, a1, a2):
    a0 = a0 % (2*np.pi); a1 = a1 % (2*np.pi); a2 = a2 % (2*np.pi)
    dif  = (a2 - a0) % (2*np.pi)
    dif1 = (a1 - a0) % (2*np.pi)
    dif2 = (a2 - a1) % (2*np.pi)
    sentido_completo = "Antihorario" if dif < np.pi else "Horario"
    sentido_interm   = "Antihorario" if (dif1 < np.pi and dif2 < np.pi) else "Horario"
    cambio_sentido = (sentido_completo != sentido_interm)
    return cambio_sentido, sentido_completo

def divide_arc3(a0, a1, a2, R, C, desired_distance):
    cambiar, sentido = sentido_recorrido(a0, a1, a2)
    L, a2 = arc_length(R, a0, a2, cambiar, sentido)
    N = int(max(1, np.ceil(L / desired_distance)))
    if N % 2 != 0:  # par
        N += 1
    angs = np.linspace(a0, a2, N+1)
    return [(C[0] + R*np.cos(ang), C[1] + R*np.sin(ang)) for ang in angs], N

def calcular_curva(p1, p2, p3, desired_distance):
    C = calculate_circle_center(p1, p2, p3)
    C = (C[0], C[1])
    R = np.hypot(p1[0]-C[0], p1[1]-C[1])
    a0, a1, a2 = (calculate_angle(p, C) for p in (p1, p2, p3))
    pts, _ = divide_arc3(a0, a1, a2, R, C, desired_distance)
    return pts

def generate_segments(commands, desired_distance):
    all_pts = []
    segmentos = []
    params = {}
    lista_alturas = []
    # 1) recolecta vértices y etiqueta tipo segmento
    for cmd in commands:
        cmd_type = cmd.split()[0]
        if cmd_type == 'LIN':
            params = parse_LIN(cmd, cmd_type)
            all_pts.append((params['X'], params['Y']))
            lista_alturas.append(params['Z'])
            segmentos.append(1)  # recta
        elif cmd_type == 'CIRC':
            pts = parse_line(cmd, cmd_type) or []
            all_pts.extend(pts)
            lista_alturas.append(params.get('Z', 0.0))
            segmentos.append(2)  # arco

    # 2) expandir en puntos (distancia <= desired_distance)
    expanded = []
    ii = 0
    acum = 0
    for seg in segmentos:
        acum += seg
        if ii > 0:
            if seg == 1:
                expanded.extend(calcular_recta(all_pts[acum-2], all_pts[acum-1], desired_distance))
            else:
                expanded.extend(calcular_curva(all_pts[acum-3], all_pts[acum-2], all_pts[acum-1], desired_distance))
        ii += 1

    return expanded, None, params, lista_alturas


# ======================= WIDGETS ==============================================
class StatusWidget(ttk.Frame):
    """
    Muestra:
      [●] NOMBRE  -  Estado: Conectado/Desconectado
    con el círculo en verde/rojo según el estado.
    """
    def __init__(self, parent, name: str, connected: bool = False, *args, **kwargs):
        super().__init__(parent, *args, **kwargs)

        self.name = name
        self._connected = connected

        self.dot = tk.Canvas(self, width=14, height=14, highlightthickness=0, bd=0)
        self.dot.grid(row=0, column=0, padx=(0, 6))

        self.lbl_name = ttk.Label(self, text=self.name, font=("Segoe UI", 10, "bold"))
        self.lbl_name.grid(row=0, column=1, sticky="w")

        self.lbl_state = ttk.Label(self, text="", font=("Segoe UI", 10))
        self.lbl_state.grid(row=0, column=2, padx=(8, 0), sticky="w")

        style = ttk.Style()
        frame_bg = style.lookup("TFrame", "background")
        if not frame_bg:
            frame_bg = self.master.cget("background") if isinstance(self.master, tk.Widget) else "white"
        self.dot.configure(bg=frame_bg)

        self._dot_id = None
        self.set_connected(self._connected)
        self.grid_columnconfigure(2, weight=1)

    def set_connected(self, connected: bool):
        self._connected = connected
        color = "#28a745" if connected else "#dc3545"  # verde/rojo
        text = "Conectado" if connected else "Desconectado"
        self.dot.delete("all")
        self._dot_id = self.dot.create_oval(2, 2, 12, 12, fill=color, outline=color)
        self.lbl_state.configure(text=f"Estado: {text}")

    def toggle(self):
        self.set_connected(not self._connected)

# ======================= BACKEND CONEXIONES ===================================
class ConnectionManager:
    """
    Gestiona conexiones a KDC (OPC UA) y LabJack T4.
    Incluye un worker de salud en hilo dedicado (~100 ms) con timeouts cortos (50 ms),
    debouncing ligero y reintentos en background.
    """
    def __init__(self):
        self.opc_client = None
        self.lj_handle  = None
        self._opc_io = Lock()
        self._lj_io  = Lock()

        # Estado y debouncing
        self._raw_opc_ok = False
        self._raw_lj_ok  = False
        self._opc_fail = 0          # contador de fallos consecutivos OPC
        self._lj_fail  = 0          # contador de fallos consecutivos LabJack
        self._fail_threshold = 1    # nº de fallos antes de marcar rojo

        # Reintentos
        self._last_opc_retry = 0.0
        self._last_lj_retry  = 0.0
        self._retry_interval = 1.0  # segundos entre reintentos automáticos

        # Sincronización
        self._lock = Lock()
        self._stop = False

        # Worker parámetros
        self._probe_period = 0.10   # 100 ms
        self._opc_timeout  = 0.05   # 50 ms
        self._lj_timeout   = 0.05

        self._worker = None

    # ---------- Conexiones explícitas ----------
    def connect_opc(self, timeout=5):
        if OPCClient is None:
            raise RuntimeError("opcua.Client no disponible.")
        client = OPCClient(opc_url)
        client.set_user(opc_user)
        client.set_password(opc_pass)

        ok = [False]; exc = [None]
        def _do_connect():
            try:
                client.connect()
                ok[0] = True
            except Exception as e:
                exc[0] = e
        t = threading.Thread(target=_do_connect, daemon=True)
        t.start()
        t.join(timeout=timeout)

        if not ok[0]:
            if t.is_alive():
                raise TimeoutError(f"Timeout conectando a OPC UA en {timeout}s")
            else:
                raise exc[0] if exc[0] else RuntimeError("Fallo desconocido conectando a OPC UA")

        self.opc_client = client
        return True

    def connect_labjack(self):
        if ljm is None:
            raise RuntimeError("labjack.ljm no disponible.")
        handle = ljm.openS("T4", "ETHERNET", IP_T4)  # ajusta a "ANY","ANY" si conviene
        info = ljm.getHandleInfo(handle)
        print(f"[LJ] Conectado: T{info[0]} via {info[1]} (serial {info[2]})")
        configurar_canales_analogicos(handle, analogicas)
        # Configuración de resolución
        ljm.eWriteName(handle, "AIN4_RESOLUTION_INDEX", 2)
        ljm.eWriteName(handle, "AIN5_RESOLUTION_INDEX", 2)
        self.lj_handle = handle
        return True

    def close_all(self):
        if self.opc_client is not None:
            try:
                self.opc_client.disconnect()
                print("[OPC] Desconectado")
            except Exception as e:
                print("[OPC] Error al desconectar:", e)
            self.opc_client = None

        if self.lj_handle is not None and ljm is not None:
            try:
                ljm.close(self.lj_handle)
                print("[LJ] Handle cerrado")
            except Exception as e:
                print("[LJ] Error al cerrar:", e)
            self.lj_handle = None

    # ---------- Utilidades con timeout ultra-corto ----------
    def _call_with_timeout(self, fn, timeout_s):
        res = {"ok": False, "val": None}
        def _runner():
            try:
                res["val"] = fn()
                res["ok"] = True
            except Exception as e:
                res["val"] = e
                res["ok"] = False
        t = threading.Thread(target=_runner, daemon=True)
        t.start()
        t.join(timeout_s)
        if t.is_alive():
            return False, TimeoutError(f"op timeout {timeout_s}s")
        return res["ok"], res["val"]

    # ---------- Sondeo de salud (no bloqueante) ----------
    def _probe_opc_once(self):
        if self.opc_client is None:
            return False
        try:
            with self._opc_io:
                node = self.opc_client.get_node("i=2258")
                ok, _ = self._call_with_timeout(lambda: node.get_value(), timeout_s=self._opc_timeout)
            return bool(ok)
        except Exception:
            return False
    
    def _probe_lj_once(self):
        if self.lj_handle is None or ljm is None:
            return False
        with self._lj_io:
            ok, _ = self._call_with_timeout(
                lambda: ljm.eReadName(self.lj_handle, "SERIAL_NUMBER"),  # más rápido que eReadNames
                timeout_s=self._lj_timeout
            )
        return bool(ok)

    # ---------- Worker de salud + reintentos ----------
    def _try_reconnect_opc(self):
        try:
            self.connect_opc(timeout=2)
            with self._lock:
                self._opc_fail = 0
                self._raw_opc_ok = True
        except Exception as e:
            print("[MON] OPC sigue caído:", e)

    def _try_reconnect_lj(self):
        try:
            self.connect_labjack()
            with self._lock:
                self._lj_fail = 0
                self._raw_lj_ok = True
        except Exception as e:
            print("[MON] LabJack sigue caído:", e)

    def _health_loop(self):
        while not self._stop:
            raw_opc = self._probe_opc_once()
            raw_lj  = self._probe_lj_once()

            with self._lock:
                self._opc_fail = 0 if raw_opc else (self._opc_fail + 1)
                self._lj_fail  = 0 if raw_lj  else (self._lj_fail  + 1)
                self._raw_opc_ok = (self._opc_fail < self._fail_threshold)
                self._raw_lj_ok  = (self._lj_fail  < self._fail_threshold)

                now = time.time()
                if not self._raw_opc_ok and (now - self._last_opc_retry) >= self._retry_interval:
                    self._last_opc_retry = now
                    threading.Thread(target=self._try_reconnect_opc, daemon=True).start()

                if not self._raw_lj_ok and (now - self._last_lj_retry) >= self._retry_interval:
                    self._last_lj_retry = now
                    threading.Thread(target=self._try_reconnect_lj, daemon=True).start()

            time.sleep(self._probe_period)

    # ---------- API para la UI ----------
    def get_status(self):
        with self._lock:
            return self._raw_opc_ok, self._raw_lj_ok

    def start_monitor(self, period_s: float = 0.10):
        self._probe_period = max(0.05, period_s)  # límite inferior 50 ms
        if self._worker is None or not self._worker.is_alive():
            self._stop = False
            self._worker = threading.Thread(target=self._health_loop, daemon=True)
            self._worker.start()

    def stop(self):
        self._stop = True

# ======================= APP (GUI) ============================================
class App(tk.Tk):
    def __init__(self, version_software):
        super().__init__()
        self.version_software = version_software
        self.title(f"Panel de Control v{self.version_software}")
        self.geometry("1200x650")
        self.minsize(800, 500)

        # Tema ttk (opcional)
        try:
            self.style = ttk.Style(self)
            if "vista" in self.style.theme_names():
                self.style.theme_use("vista")
        except Exception:
            pass
        
        # Estado de grabación
        self._recording = False
        self._record_th = None
        self._rec_lock = Lock()
        self.last_sample = None   # último paquete leído
        
        # --- Estado de puntos de referencia cargados ---
        self._ref_points_raw = []
        self._ref_point_ids  = []
        
        # Layout principal
        self.grid_rowconfigure(0, weight=0)  # topbar
        self.grid_rowconfigure(1, weight=1)  # visualización
        self.grid_columnconfigure(0, weight=1)
        
        # ---------- Estado / almacenamiento de "objetos de medición" ---------
        self._measure_objects = {}   # name -> dict(type, color, params, points)
        self._meas_seq = 1
        self._meas_palette = ["#2f9e44", "#0d6efd", "#d63384", "#f59f00", "#6f42c1", "#20c997", "#e03131"]

        # Barra superior
        topbar = ttk.Frame(self, padding=(12, 10))
        topbar.grid(row=0, column=0, sticky="ew")
        # 3 columnas: [estados] [espaciador] [controles derecha]
        topbar.grid_columnconfigure(0, weight=0)
        topbar.grid_columnconfigure(1, weight=1)   # <- empuja todo lo demás a la derecha
        topbar.grid_columnconfigure(2, weight=0)

        # Izquierda: estados
        left_box = ttk.Frame(topbar)
        left_box.grid(row=0, column=0, sticky="w", padx=(0, 8))

        self.status_kdc = StatusWidget(left_box, "KDC", connected=False)
        self.status_kdc.grid(row=0, column=0, padx=(0, 16))

        self.status_vspin = StatusWidget(left_box, "Vspin1000", connected=False)
        self.status_vspin.grid(row=0, column=1)
        
        # Derecha: cluster de controles
        right_cluster = ttk.Frame(topbar)
        right_cluster.grid(row=0, column=2, sticky="e")
        
        # Botón CARGAR (columna 0)
        self.btn_cargar = ttk.Button(right_cluster, text="CARGAR", command=self._on_load_reference)
        self.btn_cargar.grid(row=0, column=0, padx=(0, 12))
        
        # Botón GRABAR (columna 1)
        self.btn_grabar = ttk.Button(right_cluster, text="GRABAR", command=self.on_grabar_toggle)
        self.btn_grabar.grid(row=0, column=1, padx=(0, 12))
        
        # Controles de medición (en medio del cluster, columna 2)
        meas_box = ttk.Frame(right_cluster)
        meas_box.grid(row=0, column=2, padx=(0, 12))
        ttk.Label(meas_box, text="Objeto:").grid(row=0, column=0, padx=(0,4)) # Texto "Objeto" con desplegable: Circulo/Linea
        self.cbo_meas_type = ttk.Combobox(meas_box, width=8, state="readonly",
                                          values=["Círculo", "Línea"])
        self.cbo_meas_type.current(0)
        self.cbo_meas_type.grid(row=0, column=1, padx=(0,6))
        # Nombre del objeto
        self.txt_meas_name = ttk.Entry(meas_box, width=12)
        self.txt_meas_name.insert(0, f"obj{self._meas_seq}")
        self.txt_meas_name.grid(row=0, column=2, padx=(0,6))
        # Añadir objeto
        self.btn_meas_add = ttk.Button(meas_box, text="Añadir", command=self.on_add_measure_obj)
        self.btn_meas_add.grid(row=0, column=3)
        
        # Acciones a la derecha del todo: [+] y [DEFINIR SC] (columna 3)
        actions_right = ttk.Frame(right_cluster)
        actions_right.grid(row=0, column=3)
        
        self.btn_calc = ttk.Button(actions_right, text="+", width=3, command=self._open_calc_dialog)
        self.btn_calc.grid(row=0, column=0, padx=(0, 8))
        
        self.btn_sc = ttk.Button(actions_right, text="DEFINIR SC", command=self._open_sc_dialog)
        self.btn_sc.grid(row=0, column=1)
        
        self.btn_fit = ttk.Button(actions_right, text="BEST-FIT", command=self._open_fit_dialog)
        self.btn_fit.grid(row=0, column=2)
        self.btn_fit.state(["disabled"])  # arranca deshabilitado
        
        # -------------- Mantener parametros del calculo del SC ---------------
        self._coordsys_state = {
            "origin_key": None,     # nombre elegido en el combo (p.ej. "Centro: C1" o "PI3")
            "mode": "2pts",         # "2pts" (dos puntos) o "ref_line" (alineado a línea)
            "ref_key": None,        # nombre del punto (si 2pts) o de la línea (si ref_line)
            "invert": False         # invertir dirección de la línea de referencia
        }
        
        self._coordsys_params = None  # dict con "origin":(x,y), "A":ánguloº, "ux","uy"...
        # ---------------------------------------------------------------------
        # Secuencia para nombrar puntos (intersecciones)
        self._point_seq = 1
        # Sistema de coordenadas activo (origen y ángulo del eje X en grados)
        self._coordsys = None  # dict: {"x": ox, "y": oy, "A": ang_deg}
        
        # Deshabilitar UI de medición durante la grabación
        self._set_measure_ui_enabled(False)

        # Zona de visualización
        viz_container = ttk.Frame(self, padding=(12, 0, 12, 12))
        viz_container.grid(row=1, column=0, sticky="nsew")
        viz_container.grid_rowconfigure(0, weight=1)
        viz_container.grid_columnconfigure(0, weight=1)

        self.visual_frame = ttk.Frame(viz_container, padding=6, borderwidth=1, relief="solid")
        self.visual_frame.grid(row=0, column=0, sticky="nsew")

        # Canvas pantalla
        self.visual_canvas = tk.Canvas(self.visual_frame, background="#ffffff") # fondo blanco
        # self.visual_canvas = tk.Canvas(self.visual_frame, background="#111317") # gris oscuro
        self.visual_canvas.pack(fill="both", expand=True)
        
        # ---- Config de trazado (mm -> px) ----
        self._plot_scale = 2.0      # píxeles por mm (ajústalo)
        self._plot_margin = 40      # margen interior para ejes/etiquetas
        self._plot_point_r = 2      # radio punto
        self._plot_max_points = 5000
        self._points_logical = []   # historial lógico (mm)
        self._points_ids = []       # ids de círculos en canvas
        self._origin_px = (0, 0)    # se calcula en _redraw_grid_axes
        self._grid_minor_mm = 10
        self._grid_major_mm = 50
        
        # ---- Vista / autoscale ----
        self._view_center_mm = (0.0, 0.0)   # centro lógico mostrado en el centro del canvas
        self._autoscale_enabled = True
        self._autoscale_scheduled = False
        self._autoscale_interval_ms = 200   # máx 5 Hz de autoscale
        
        # Navegación
        self._nav_enabled = False
        self._dragging = False
        self._drag_px_start = (0, 0)
        self._view_center_mm_start = (0.0, 0.0)
        
        # --- Selección (marquee) ---
        self._selected_indices = set()     # índices dentro de _points_logical
        self._sel_rect_id = None
        self._sel_active = False
        self._sel_start_px = (0, 0)
        self._sel_mode = "set"  # "set" | "add" | "remove"
        c = self.visual_canvas
        # Shift = seleccionar (reemplazar)
        c.bind("<Shift-ButtonPress-1>",  self._on_sel_start_set)
        c.bind("<Shift-B1-Motion>",      self._on_sel_drag)
        c.bind("<Shift-ButtonRelease-1>",self._on_sel_end)
        # Ctrl = añadir
        c.bind("<Control-ButtonPress-1>",  self._on_sel_start_add)
        c.bind("<Control-B1-Motion>",      self._on_sel_drag)
        c.bind("<Control-ButtonRelease-1>",self._on_sel_end)
        # Alt/Option = quitar
        c.bind("<Alt-ButtonPress-1>",    self._on_sel_start_remove)
        c.bind("<Alt-B1-Motion>",        self._on_sel_drag)
        c.bind("<Alt-ButtonRelease-1>",  self._on_sel_end)
        
        # --- Hint de selección (3 líneas) ---
        self._hint_lines = [
            "Shift + arrastrar con el botón izquierdo para seleccionar puntos",
            "Ctrl + arrastrar para añadir puntos",
            "Alt/Option + arrastrar para deseleccionar puntos",
        ]
        self._hint_visible = False
        
        # Enlaza eventos de navegación (se activan solo cuando _nav_enabled=True)
        self._bind_navigation_events()
        self._enable_navigation(show_hint=False) # Activa navegación desde el inicio, pero SIN hint
        
        # Redibuja rejilla/ejes al redimensionar
        self.visual_canvas.bind("<Configure>", lambda e: self._redraw_grid_axes())
        self.after(0, self._redraw_grid_axes)
        
        # --- Overlay del ajuste de círculo ---
        # Actualizar tamaño y posicion al desplazarse por la pantalla
        self._fit_circle_params = None   # (cx_mm, cy_mm, r_mm)
        
        # --- BEST FIT con referencia ---
        self._bf_preview_pts = None  # puntos del preview BEST-FIT (violeta)
        # Resultado del BEST-FIT (para overlay)
        self._bestfit_result = None

        # Backend conexiones
        self.conn = ConnectionManager()

        # Intento de conexión al abrir (no bloquea UI)
        self.after(100, self._connect_on_start)

        # Arranca monitor de salud en background (100 ms)
        self.conn.start_monitor(period_s=0.10)

        # Tick UI muy ligero (solo lee flags y pinta)
        self._ui_period_ms = 100
        self.after(self._ui_period_ms, self._ui_tick)

        # Cierre limpio
        self.protocol("WM_DELETE_WINDOW", self._on_close)
    

    # --------- Acciones de botones -------------------------------------------
    
    def _se2_from_R_t(self, R, t):
        T = np.eye(3); T[:2,:2] = np.asarray(R).reshape(2,2)
        T[:2,2] = np.asarray(t).reshape(2,)
        return T
    
    def _se2_inv(self, T):
        Ri = T[:2,:2].T
        ti = -Ri @ T[:2,2]
        Ti = np.eye(3); Ti[:2,:2] = Ri; Ti[:2,2] = ti
        return Ti
    
    def _compose_ui_se2(self):
        """SE(2) de la transformación que aplicas a la referencia para visualizar (rot+shift)."""
        th = math.radians(REF_ROT_DEG)
        c, s = math.cos(th), math.sin(th)
        R = np.array([[c, -s], [s, c]])
        ox, oy = REF_ROT_ORIGIN
        # Rotación alrededor de (ox,oy)
        t_rot = np.array([ox, oy]) - R @ np.array([ox, oy])
        # + desplazamiento en XY
        t = t_rot + np.array([REF_SHIFT_X, REF_SHIFT_Y])
        T = np.eye(3); T[:2,:2] = R; T[:2,2] = t
        return T

    def _open_fit_dialog(self):
        if self._recording:
            self._notify("Detén la grabación para ejecutar BEST-FIT.")
            return
        if len(self._points_logical) < 3 or not getattr(self, "_ref_points_raw", None):
            self._notify("Necesitas puntos grabados y una referencia cargada.")
            return
    
        dlg = tk.Toplevel(self)
        dlg.title("BEST-FIT (rigid 2D)")
        dlg.transient(self)
        dlg.resizable(False, False)
        dlg.geometry("420x180")
    
        frm = ttk.Frame(dlg, padding=12)
        frm.grid(row=0, column=0, sticky="nsew")
    
        # Campos de salida
        var_dx = tk.StringVar(value="—")
        var_dy = tk.StringVar(value="—")
        var_th = tk.StringVar(value="—")
        var_pv = tk.StringVar(value="—")
        var_rm = tk.StringVar(value="—")
    
        r = 0
        ttk.Label(frm, text="ΔX (mm):").grid(row=r, column=0, sticky="w"); ttk.Label(frm, textvariable=var_dx).grid(row=r, column=1, sticky="w"); r+=1
        ttk.Label(frm, text="ΔY (mm):").grid(row=r, column=0, sticky="w"); ttk.Label(frm, textvariable=var_dy).grid(row=r, column=1, sticky="w"); r+=1
        ttk.Label(frm, text="Rotación θ (°):").grid(row=r, column=0, sticky="w"); ttk.Label(frm, textvariable=var_th).grid(row=r, column=1, sticky="w"); r+=1
        ttk.Label(frm, text="Origen rotación:").grid(row=r, column=0, sticky="w"); ttk.Label(frm, textvariable=var_pv).grid(row=r, column=1, sticky="w"); r+=1
        ttk.Label(frm, text="RMSE (mm):").grid(row=r, column=0, sticky="w"); ttk.Label(frm, textvariable=var_rm).grid(row=r, column=1, sticky="w"); r+=1
    
        frm_btn = ttk.Frame(frm)
        frm_btn.grid(row=r, column=0, columnspan=2, sticky="e", pady=(12,0))
    
        def do_compute():
            res = self._compute_bestfit()
            if not res:
                self._notify("No se pudo calcular BEST-FIT.")
                return None
        
            # -------------------- MOSTRAR A->B (fit directo) --------------------
            R, t = res["R"], res["t"]
            var_dx.set(f"{float(t[0]):.3f}")
            var_dy.set(f"{float(t[1]):.3f}")
            var_th.set(f"{res['theta_deg']:.3f}")
            var_rm.set(f"{res['rmse']:.3f}")
            if res.get("origin_meas") is not None:
                omx, omy = res["origin_meas"]
                var_pv.set(f"({omx:.3f}, {omy:.3f})")
            else:
                var_pv.set("—")
        
            # -------------------- CALCULAR Tba = (A->B)^-1 ----------------------
            Ri = R.T
            ti = -Ri @ t
            theta_ba = -float(res['theta_deg'])
            # El "origen" natural de la inversa es el centro de B (referencia)
            obx, oby = res.get("origin_ref", (0.0, 0.0))
        
            # Guarda todo dentro del resultado para el overlay global
            self._bestfit_result = {
                # A->B (como hasta ahora)
                "R": R, "t": t,
                "tx": float(t[0]), "ty": float(t[1]),
                "theta_deg": float(res["theta_deg"]),
                "origin": res.get("origin_meas", (0.0, 0.0)),
                "rmse": float(res["rmse"]),
                # B->A (Tba) – corrección de BASE mostrada
                "ba_dx": float(ti[0]),
                "ba_dy": float(ti[1]),
                "ba_theta_deg": float(theta_ba),
                "ba_origin": (float(obx), float(oby)),
            }
        
            # Si tu BASE real está en la referencia "cruda" (sin REF_ROT/SHIFT),
            # descomenta lo siguiente para incorporar tu T_ui:
            # T_ui = self._compose_ui_se2()  # ver helper más abajo
            # Tba_se2 = self._se2_inv(self._se2_from_R_t(R, t)) @ T_ui
            # self._bestfit_result["ba_dx"] = float(Tba_se2[0,2])
            # self._bestfit_result["ba_dy"] = float(Tba_se2[1,2])
            # self._bestfit_result["ba_theta_deg"] = math.degrees(math.atan2(Tba_se2[1,0], Tba_se2[0,0]))
            # self._bestfit_result["ba_origin"] = self._bestfit_result["ba_origin"]  # opcional
        
            return self._bestfit_result
        
        def do_preview():
            res = do_compute()
            if not res:
                return
            # Pinta los puntos medidos transformados A->B (violeta) para ver encaje
            pts_violet = self._apply_rigid_transform(self._points_logical, res["R"], res["t"])
            self._set_bestfit_preview(pts_violet)
            self._maybe_center_on_preview(pts_violet)
            # Refresca overlay con A->B y Tba
            # self._show_bestfit_overlay()
        
        def do_apply():
            res = do_compute()
            if not res:
                return
            # Mantén la superposición violeta y muestra overlay con Tba
            pts_violet = self._apply_rigid_transform(self._points_logical, res["R"], res["t"])
            self._set_bestfit_preview(pts_violet)
            self._show_bestfit_overlay()
            dlg.destroy()
        
        # .....................................................................    
        btn_prev  = ttk.Button(frm_btn, text="Preview", command=do_preview)
        btn_apply = ttk.Button(frm_btn, text="Aplicar", command=do_apply)
        btn_close = ttk.Button(frm_btn, text="Cerrar", command=dlg.destroy)
    
        btn_prev.grid(row=0, column=0, padx=(0,8))
        btn_apply.grid(row=0, column=1, padx=(0,8))
        btn_close.grid(row=0, column=2)
    
        dlg.grab_set()
        dlg.wait_window()
    

    def _apply_rigid_transform(self, pts, R, t):
        """Aplica x' = R @ x + t a una lista de puntos [(x,y),...]."""
        if pts is None:
            return []
        # R puede ser np.ndarray 2x2; t puede ser (2,) o (2,1)
        t = np.asarray(t).reshape(2,)
        R = np.asarray(R).reshape(2,2)
        out = []
        for (x, y) in pts:
            xp, yp = R @ np.array([x, y]) + t
            out.append((float(xp), float(yp)))
        return out
    
    def _set_bestfit_preview(self, pts):
        """Guarda puntos lógicos transformados y los pinta."""
        self._bf_preview_pts = list(pts) if pts else []
        self._redraw_bestfit_preview()
    
    def _clear_bestfit_preview(self):
        """Elimina la previsualización BEST-FIT."""
        self._bf_preview_pts = None
        try:
            self.visual_canvas.delete("bfprev")
        except Exception:
            pass
    
    def _redraw_bestfit_preview(self):
        """Repinta los puntos de previsualización con la vista actual."""
        c = self.visual_canvas
        c.delete("bfprev")
        if not self._bf_preview_pts:
            return
        r = max(1, self._plot_point_r)
        # color violeta para distinguirlos
        for (x_mm, y_mm) in self._bf_preview_pts:
            px, py = self._to_canvas(x_mm, y_mm)
            c.create_oval(px - r, py - r, px + r, py + r,
                          fill="#845ef7", outline="", tags=("bfprev",))

    def _fit_to_points(self, pts, padding_mm=10.0, min_px_per_mm=0.2, max_px_per_mm=20.0):
        if not pts:
            return
        xs = [p[0] for p in pts]; ys = [p[1] for p in pts]
        minx, maxx = min(xs), max(xs); miny, maxy = min(ys), max(ys)
        cx = (minx + maxx)/2.0; cy = (miny + maxy)/2.0
        rangex = max(maxx - minx, 1e-6) + 2*padding_mm
        rangey = max(maxy - miny, 1e-6) + 2*padding_mm
        w = max(1, self.visual_canvas.winfo_width()  - 2*self._plot_margin)
        h = max(1, self.visual_canvas.winfo_height() - 2*self._plot_margin)
        s = min(w / rangex, h / rangey)
        s = max(min_px_per_mm, min(s, max_px_per_mm))
        self._view_center_mm = (cx, cy)
        self._plot_scale = s
        self._apply_view_change()
    
    def _maybe_center_on_preview(self, pts):
        if not pts:
            return
        xmin, xmax, ymin, ymax = self._viewport_logical_bbox()
        inside = sum(1 for x,y in pts if xmin <= x <= xmax and ymin <= y <= ymax)
        # Si menos del 20% del preview está dentro de la vista, encaja a preview
        if inside < max(5, int(0.2 * len(pts))):
            self._fit_to_points(pts, padding_mm=10.0)


    def _on_load_reference(self):
        path = filedialog.askopenfilename(
            title="Selecciona archivo de trayectorias (.src o .txt)",
            filetypes=[("KUKA / Texto", "*.src *.txt"), ("Todos", "*.*")]
        )
        if not path:
            self._notify("Carga cancelada.")
            return
        try:
            raw = read_commands_from_file(path)
            commands = [ln for ln in raw if ln.strip().startswith(("LIN", "CIRC"))]
            if not commands:
                self._notify("El archivo no contiene comandos LIN/CIRC.")
                return
            pts, _seg, _params, _zs = generate_segments(commands, REF_TOL_MM)
            self._ref_points_raw = pts or []                 # <--- IMPORTANTE
            self._redraw_reference_points()
            self._notify(f"{len(self._ref_points_raw)} puntos de referencia cargados.")
        except Exception as e:
            self._notify(f"Error al cargar: {e}")

    
    def on_grabar_toggle(self):
        """Arranca/para la grabación. Solo permitido si ambos dispositivos están conectados."""
        opc_ok, lj_ok = self.conn.get_status()
        if not self._recording:
            if not (opc_ok and lj_ok):
                self._notify("Conecta KDC y Vspin1000")
                return
            # Inicia grabación
            self._recording = True
            self.btn_grabar.configure(text="DETENER")
            self._notify("Grabando...")
            # self._clear_points()                    #  (opcional) limpia puntos anteriores, mantiene ejes/cuadrícula
            self._autoscale_enabled = True
            self._disable_navigation()              # desactiva navegación al grabar
            self._set_measure_ui_enabled(False)     # desactiva UI de medición al grabar
            # self._clear_measure_objects()           # (opcional) limpiar objetos previos
            self._record_th = threading.Thread(target=self._record_loop, daemon=True)
            self._record_th.start()
            try: self.btn_sc.state(["disabled"])    # deshabilita SC
            except: pass
        else:
            # Detiene grabación
            self._recording = False
            self.btn_grabar.configure(text="GRABAR")
            self._notify("Grabación detenida (arrastrar=pan, rueda=zoom, Ctrl+0=encajar)")
            self._autoscale_enabled = False
            self._enable_navigation(show_hint=True)  # activa navegación y muestra hint
            self._set_measure_ui_enabled(True)       # activa UI de medición
            try: self.btn_sc.state(["!disabled"])    # habilita SC
            except: pass
    
    def _transform_ref_points(self, pts):
        """Aplica rotación CCW (REF_ROT_DEG) alrededor de REF_ROT_ORIGIN y luego desplazamiento (REF_SHIFT_X/Y)."""
        if not pts:
            return []
        ox, oy = REF_ROT_ORIGIN
        th = math.radians(REF_ROT_DEG)
        ct, st = math.cos(th), math.sin(th)
        dx, dy = REF_SHIFT_X, REF_SHIFT_Y
        out = []
        for (x, y) in pts:
            # rotar alrededor de (ox, oy)
            xr = ox + ct*(x - ox) - st*(y - oy)
            yr = oy + st*(x - ox) + ct*(y - oy)
            # desplazar
            out.append((xr + dx, yr + dy))
        return out

    def _redraw_reference_points(self):
        """Redibuja los puntos de referencia (en gris) según la vista actual y la transformación configurada."""
        c = self.visual_canvas
        c.delete("ref")
        if not self._ref_points_raw:
            return
        pts = self._transform_ref_points(self._ref_points_raw)
        r = max(1, self._plot_point_r)  # radio similar al de puntos normales
        for (x_mm, y_mm) in pts:
            px, py = self._to_canvas(x_mm, y_mm)
            c.create_oval(px - r, py - r, px + r, py + r,
                          fill="#6c757d", outline="", tags=("ref",))


    def _clear_points(self):
        """Borra solo puntos (mantiene ejes/cuadrícula) al inicio de la grabacion"""
        self.visual_canvas.delete("points")
        self._points_ids.clear()
        self._points_logical.clear()
        
    def _clear_measure_objects(self):
        self._measure_objects.clear()
        self.visual_canvas.delete("meas")
        
    def _to_canvas(self, x_mm, y_mm):
        """
        Convertir mm ↔ píxeles 
        Usa el centro lógico de la vista (self._view_center_mm) para pintar """
        ox, oy = self._origin_px
        s = self._plot_scale
        cx_mm, cy_mm = getattr(self, "_view_center_mm", (0.0, 0.0))
        return ox + (x_mm - cx_mm) * s, oy - (y_mm - cy_mm) * s
    
    def _from_canvas(self, x_px, y_px):
        """ Convertir píxeles ↔ mm """
        ox, oy = self._origin_px
        s = max(1e-9, self._plot_scale)
        cx_mm, cy_mm = getattr(self, "_view_center_mm", (0.0, 0.0))
        x_mm = cx_mm + (x_px - ox) / s
        y_mm = cy_mm + (oy - y_px) / s
        return x_mm, y_mm
    
    def _choose_grid_steps(self, target_px=32):
        """
        Devuelve (minor_mm, major_mm) para que la separación entre líneas menores
        sea ≈ target_px (1–2–5 × 10^n). Las mayores = 5 * menores.
        """
        s = max(1e-9, self._plot_scale)  # px/mm
        target_mm = target_px / s
        if target_mm <= 0:
            return 10.0, 50.0
        exp = math.floor(math.log10(target_mm))
        base = 10 ** exp
        for k in (1, 2, 5, 10):
            minor = base * k
            if minor >= target_mm:
                major = minor * 5
                return minor, major
        return base * 10, base * 50

    def _redraw_grid_axes(self):
        """Cuadrícula y ejes en función de la vista actual (no fuerza (0,0) en pantalla)."""
        c = self.visual_canvas
        w = max(1, c.winfo_width())
        h = max(1, c.winfo_height())
        c.delete("grid"); c.delete("axes")   # deja puntos y HUD
    
        # Centro de pantalla en píxeles
        self._origin_px = (w // 2, h // 2)
        ox, oy = self._origin_px
        m = getattr(self, "_plot_margin", 40)
    
        # Paso dinámico agradable
        minor_mm, major_mm = self._choose_grid_steps(target_px=32)
        s = max(1e-9, self._plot_scale)  # px/mm
    
        # Coordenadas lógicas visibles en los bordes útil
        cx_mm, cy_mm = getattr(self, "_view_center_mm", (0.0, 0.0))
        x_left_mm   = cx_mm + (m - ox) / s
        x_right_mm  = cx_mm + (w - m - ox) / s
        y_top_mm    = cy_mm + (oy - m) / s
        y_bottom_mm = cy_mm + (oy - (h - m)) / s
    
        # Funciones para primer múltiplo visible
        def first_multiple_geq(a, step):
            return math.ceil(a / step) * step
    
        # ---- Líneas menores (claritas) ----
        # Verticales
        x_mm = first_multiple_geq(x_left_mm, minor_mm)
        while x_mm <= x_right_mm:
            x_px, _ = self._to_canvas(x_mm, cy_mm)
            c.create_line(x_px, m, x_px, h - m, fill="#f1f3f5", tags="grid")
            x_mm += minor_mm
        # Horizontales
        y_mm = first_multiple_geq(y_bottom_mm, minor_mm)
        while y_mm <= y_top_mm:
            _, y_px = self._to_canvas(cx_mm, y_mm)
            c.create_line(m, y_px, w - m, y_px, fill="#f1f3f5", tags="grid")
            y_mm += minor_mm
    
        # ---- Líneas mayores (un poco más marcadas) ----
        # Verticales
        x_mm = first_multiple_geq(x_left_mm, major_mm)
        while x_mm <= x_right_mm:
            x_px, _ = self._to_canvas(x_mm, cy_mm)
            c.create_line(x_px, m, x_px, h - m, fill="#e5e8eb", width=1, tags="grid")
            x_mm += major_mm
        # Horizontales
        y_mm = first_multiple_geq(y_bottom_mm, major_mm)
        while y_mm <= y_top_mm:
            _, y_px = self._to_canvas(cx_mm, y_mm)
            c.create_line(m, y_px, w - m, y_px, fill="#e5e8eb", width=1, tags="grid")
            y_mm += major_mm
    
        # ---- Ejes orientativos centrados en la vista ----
        c.create_line(m, oy, w - m, oy, width=2, fill="#343a40", tags="axes")  # X (viewport)
        c.create_line(ox, m, ox, h - m, width=2, fill="#343a40", tags="axes")  # Y (viewport)
        
        # Reposiciona el mensaje (instrucciones de seleccion) si está activo
        self._position_selection_hint()

    def _draw_point(self, x_mm, y_mm, selected=False):
        cx, cy = self._to_canvas(x_mm, y_mm)
        r = self._plot_point_r + (1 if selected else 0)
        color = "#d63384" if selected else "#0d6efd" # magenta (seleccion) azul (resto)
        pid = self.visual_canvas.create_oval(cx - r, cy - r, cx + r, cy + r,
                                             fill=color, outline="", tags="points")
        self._points_ids.append(pid)
        
    def _add_point(self, x_mm, y_mm):
        """Añade al historial y dibuja (con recorte suave). Agenda autoscale."""
        self._points_logical.append((x_mm, y_mm))
        self._draw_point(x_mm, y_mm, selected=False)
        if len(self._points_ids) > self._plot_max_points:
            old = self._points_ids.pop(0)
            try: self.visual_canvas.delete(old)
            except Exception: pass
        # autoscale con 'debounce'
        if self._autoscale_enabled and not self._autoscale_scheduled:
            self._autoscale_scheduled = True
            self.after(self._autoscale_interval_ms, self._autoscale_tick)
    
    def _repaint_points(self):
        self.visual_canvas.delete("points")
        self._points_ids.clear()
        for i, (x_mm, y_mm) in enumerate(self._points_logical):
            self._draw_point(x_mm, y_mm, selected=(i in self._selected_indices))
    
    def _autoscale_tick(self):
        self._autoscale_scheduled = False
        self._autoscale()
    
    def _autoscale(self, padding_mm=10.0, min_px_per_mm=0.2, max_px_per_mm=20.0):
        """Ajusta escala/centro para que entren todos los puntos con un margen."""
        if not self._points_logical:
            return
        # bbox lógico
        xs = [p[0] for p in self._points_logical]
        ys = [p[1] for p in self._points_logical]
        minx, maxx = min(xs), max(xs)
        miny, maxy = min(ys), max(ys)
    
        # centro lógico y rangos con padding
        cx = (minx + maxx) / 2.0
        cy = (miny + maxy) / 2.0
        rangex = max(maxx - minx, 1e-6) + 2 * padding_mm
        rangey = max(maxy - miny, 1e-6) + 2 * padding_mm
    
        # área útil en px
        w = max(1, self.visual_canvas.winfo_width()  - 2 * self._plot_margin)
        h = max(1, self.visual_canvas.winfo_height() - 2 * self._plot_margin)
    
        # escala que cabe en ambos ejes
        s = min(w / rangex, h / rangey)
        s = max(min_px_per_mm, min(s, max_px_per_mm))  # clamp
    
        # aplica
        self._view_center_mm = (cx, cy)
        self._plot_scale = s
    
        # redibuja grid/ejes y puntos (y circulo ajuste)
        self._redraw_grid_axes()
        self._redraw_reference_points() # coord referencia pieza
        self._repaint_points()
        self._redraw_fit_circle()
        self._redraw_measure_objects()
        self._redraw_coordsys() # sistema de coordenadas
        self._show_bestfit_overlay() 
        self._redraw_bestfit_preview()
        
        
    # ----------Métodos de soporte para mediciones ----------------------------
    def _set_measure_ui_enabled(self, enabled: bool):
        state = "normal" if enabled else "disabled"
        for w in (getattr(self, "cbo_meas_type", None),
                  getattr(self, "txt_meas_name", None),
                  getattr(self, "btn_meas_add", None)):
            try:
                if isinstance(w, ttk.Combobox):
                    w.state([("!disabled" if enabled else "disabled")])
                elif w is not None:
                    w.configure(state=state)
            except Exception:
                pass
    
    def _next_color(self):
        c = self._meas_palette[(self._meas_seq - 1) % len(self._meas_palette)]
        self._meas_seq += 1
        return c
    
    def _viewport_logical_bbox(self):
        """BBox lógico visible (teniendo en cuenta márgenes)."""
        c = self.visual_canvas
        w = max(1, c.winfo_width()); h = max(1, c.winfo_height())
        m = self._plot_margin
        ox, oy = self._origin_px
        s = max(1e-9, self._plot_scale)
        cx_mm, cy_mm = self._view_center_mm
        x_left   = cx_mm + (m - ox) / s
        x_right  = cx_mm + (w - m - ox) / s
        y_top    = cy_mm + (oy - m) / s
        y_bottom = cy_mm + (oy - (h - m)) / s
        return x_left, x_right, y_bottom, y_top  # (xmin, xmax, ymin, ymax)
    
    def _clip_line_to_rect(self, x0, y0, vx, vy, xmin, xmax, ymin, ymax):
        """Intersección de línea infinita p(t)=p0+t*v con el rectángulo lógico visible."""
        pts = []
        eps = 1e-12
        # Fronteras x = const
        if abs(vx) > eps:
            for x in (xmin, xmax):
                t = (x - x0) / vx
                y = y0 + vy * t
                if ymin - eps <= y <= ymax + eps:
                    pts.append((x, y))
        # Fronteras y = const
        if abs(vy) > eps:
            for y in (ymin, ymax):
                t = (y - y0) / vy
                x = x0 + vx * t
                if xmin - eps <= x <= xmax + eps:
                    pts.append((x, y))
        # De-dup y deja 2 extremos
        uniq = []
        for p in pts:
            if all((abs(p[0]-q[0]) > 1e-9 or abs(p[1]-q[1]) > 1e-9) for q in uniq):
                uniq.append(p)
        if len(uniq) >= 2:
            # coge los dos más alejados
            a, b, bestd = None, None, -1
            for i in range(len(uniq)):
                for j in range(i+1, len(uniq)):
                    d = (uniq[i][0]-uniq[j][0])**2 + (uniq[i][1]-uniq[j][1])**2
                    if d > bestd: bestd, a, b = d, uniq[i], uniq[j]
            return a, b
        # fallback: traza un segmento corto alrededor del centroide
        L = max(xmax-xmin, ymax-ymin) * 1.2
        return (x0 - vx*L, y0 - vy*L), (x0 + vx*L, y0 + vy*L)
    
    def _unique_name(self, base: str) -> str:
        name = base
        k = 2
        while name in self._measure_objects:
            name = f"{base}_{k}"
            k += 1
        return name

    def _redraw_measure_objects(self):
        c = self.visual_canvas
        c.delete("meas")
        if not self._measure_objects:
            return
    
        xmin, xmax, ymin, ymax = self._viewport_logical_bbox()
    
        for name, obj in self._measure_objects.items():
            typ = obj["type"]; color = obj["color"]; prm = obj["params"]
    
            if typ == "circle":
                cx, cy, r = prm["cx"], prm["cy"], prm["r"]
                cxp, cyp = self._to_canvas(cx, cy)
                rp = max(1.0, r * max(1e-9, self._plot_scale))
                c.create_oval(cxp - rp, cyp - rp, cxp + rp, cyp + rp,
                              outline=color, width=2, tags=("meas", f"meas:{name}"))
                c.create_oval(cxp-3, cyp-3, cxp+3, cyp+3,
                              fill=color, outline="", tags=("meas", f"meas:{name}"))
                c.create_text(cxp + rp + 8, cyp, anchor="w",
                              text=f"{name}: C=({cx:.3f},{cy:.3f}) R={r:.3f}",
                              fill=color, font=("Segoe UI", 10, "bold"),
                              tags=("meas", f"meas:{name}"))
    
            elif typ == "line":
                x0, y0, vx, vy = prm["x0"], prm["y0"], prm["vx"], prm["vy"]
                (xa, ya), (xb, yb) = self._clip_line_to_rect(x0, y0, vx, vy, xmin, xmax, ymin, ymax)
                xap, yap = self._to_canvas(xa, ya)
                xbp, ybp = self._to_canvas(xb, yb)
                c.create_line(xap, yap, xbp, ybp, fill=color, width=2,
                              tags=("meas", f"meas:{name}"))
                if abs(xb-xa) > 1e-9:
                    m = (yb-ya)/(xb-xa); b = ya - m*xa
                    label = f"{name}: y={m:.4f}x+{b:.4f}"
                else:
                    label = f"{name}: x={xa:.3f}"
                c.create_text(xap + 8, yap - 8, anchor="nw",
                              text=label, fill=color, font=("Segoe UI", 10, "bold"),
                              tags=("meas", f"meas:{name}"))
            
            elif typ == "point":
                x, y = prm["x"], prm["y"]
                xp, yp = self._to_canvas(x, y)
                c.create_oval(xp-3, yp-3, xp+3, yp+3,
                              fill=color, outline="", tags=("meas", f"meas:{name}"))
                c.create_text(xp+6, yp-6, anchor="nw",
                              text=f"{name}: ({x:.3f}, {y:.3f})",
                              fill=color, font=("Segoe UI", 10, "bold"),
                              tags=("meas", f"meas:{name}"))

                
    # Crear un objeto de medición desde la selección (botón “Añadir”)
    def on_add_measure_obj(self):
        if not self._nav_enabled:
            self._notify("Finaliza la grabación para crear objetos de medición.")
            return
    
        pts = self.get_selected_points()
        if not pts:
            self._notify("Selecciona puntos (Shift/Ctrl/Alt) para crear el objeto.")
            return
    
        typ = (self.cbo_meas_type.get() or "Círculo").lower()
        name = (self.txt_meas_name.get() or "").strip()
        if not name:
            name = f"obj{self._meas_seq}"
        if name in self._measure_objects:
            # hace único el nombre
            i = 2
            base = name
            while f"{base}_{i}" in self._measure_objects:
                i += 1
            name = f"{base}_{i}"
    
        color = self._next_color()
    
        try:
            if typ.startswith("círculo") or typ.startswith("circulo"):
                if len(pts) < 3:
                    self._notify("Círculo: mínimo 3 puntos.")
                    return
                cx, cy, r = fit_circle_kasa(pts)
                params = {"cx": cx, "cy": cy, "r": r}
                self._measure_objects[name] = {"type": "circle", "color": color,
                                               "params": params, "points": pts}
                self._notify(f"Círculo '{name}' creado.")
    
            elif typ.startswith("línea") or typ.startswith("linea"):
                if len(pts) == 2:
                    (x0, y0), (x1, y1) = pts[0], pts[1]
                    vx, vy = (x1-x0), (y1-y0)
                    nrm = math.hypot(vx, vy) or 1.0
                    vx, vy = vx/nrm, vy/nrm
                    x0, y0 = (x0+x1)/2.0, (y0+y1)/2.0  # centro del segmento
                else:
                    x0, y0, vx, vy = fit_line_tls(pts)
                params = {"x0": x0, "y0": y0, "vx": vx, "vy": vy}
                self._measure_objects[name] = {"type": "line", "color": color,
                                               "params": params, "points": pts}
                self._notify(f"Línea '{name}' creada.")
            else:
                self._notify("Tipo no reconocido (usa 'Círculo' o 'Línea').")
                return
    
        except Exception as e:
            self._notify(f"Error creando objeto: {e}")
            return
    
        # pinta y actualiza nombre sugerido siguiente
        self._redraw_measure_objects()
        self.txt_meas_name.delete(0, "end")
        self.txt_meas_name.insert(0, f"obj{self._meas_seq}")
        
        # Deja la selección vacía para crear el siguiente objeto desde cero
        self._apply_selection(set())
    
    def _fit_all(self):
        """Encaja todos los puntos (una sola vez) y congela navegación en esa vista."""
        if not self._points_logical:
            return
        # self._autoscale_enabled = False #Si quieres que después de encajar se congele el auto-escale
        self._autoscale()  # usa tu autoscale existente
    
    # -------------------- Desplazamiento raton -------------------------------
    def _bind_navigation_events(self):
        c = self.visual_canvas
        # Pan con botón izquierdo
        c.bind("<ButtonPress-1>", self._on_pan_start)
        c.bind("<B1-Motion>", self._on_pan_move)
        c.bind("<ButtonRelease-1>", self._on_pan_end)
        # Zoom con rueda (Windows/Mac)
        c.bind("<MouseWheel>", self._on_wheel)
        # Zoom con rueda (Linux X11)
        c.bind("<Button-4>", self._on_wheel)
        c.bind("<Button-5>", self._on_wheel)
        # Atajo encajar todo
        self.bind_all("<Control-Key-0>", lambda e: self._fit_all())
        self.bind_all("<Control-Key-KP_0>", lambda e: self._fit_all()) # (Opcional) también en el keypad:
    
    def _enable_navigation(self, show_hint: bool = False):
        self._nav_enabled = True
        try: 
            self.visual_canvas.configure(cursor="fleur")
            self.visual_canvas.focus_set()
        except Exception: 
            pass
        if show_hint:
            self._show_selection_hint(True)
    
    def _disable_navigation(self):
        self._nav_enabled = False
        self.visual_canvas.configure(cursor="")
        self._show_selection_hint(False)
    
    def _apply_view_change(self):
        """Reaplica vista: redibuja cuadrícula y re-pinta puntos y objetos de ajuste"""
        self._redraw_grid_axes()
        self._redraw_reference_points()  # Coordenadas referencia pieza
        self._repaint_points()
        self._redraw_fit_circle()
        self._redraw_measure_objects()
        self._redraw_coordsys() # Sistema de coordenadas
        self._show_bestfit_overlay() 
        self._redraw_bestfit_preview()
        
        
    def _on_pan_start(self, event):
        # Si hay una selección en curso (marquee), no inicies pan
        if getattr(self, "_sel_active", False):
            return
        # Navegación debe estar habilitada
        if not self._nav_enabled:
            return
    
        self._dragging = True
        self._drag_px_start = (event.x, event.y)
        self._view_center_mm_start = self._view_center_mm
    
    def _on_pan_move(self, event):
        """
        Arrastrar con botón izquierdo → desplaza (pan).
        """
        if not (self._nav_enabled and self._dragging):
            return
        x0, y0 = self._drag_px_start
        dx = event.x - x0
        dy = event.y - y0
        s = max(1e-9, self._plot_scale)  # px/mm
        cx0, cy0 = self._view_center_mm_start
        # Pan: mover contenido con el ratón
        cx = cx0 - dx / s
        cy = cy0 + dy / s
        self._view_center_mm = (cx, cy)
        self._apply_view_change()
    
    def _on_pan_end(self, event):
        self._dragging = False
    
    def _on_wheel(self, event):
        """
        ZOOM pantalla
        (Rueda del ratón → zoom hacia el puntero)
        Si quieres invertir el sentido del zoom, cambia factor = 1.1 ... por 0.9 y al revés.
        """
        if not self._nav_enabled:
            return
    
        # Dirección de zoom
        direction = 0
        if hasattr(event, "delta") and event.delta != 0:
            direction = 1 if event.delta > 0 else -1       # Windows/Mac
        elif hasattr(event, "num"):
            direction = 1 if event.num == 4 else -1        # Linux: 4=up, 5=down
        if direction == 0:
            return
    
        # Factor y límites
        factor = 1.1 if direction > 0 else (1/1.1)
        s_old = self._plot_scale
        s_new = max(0.1, min(50.0, s_old * factor))        # clamp px/mm
        if abs(s_new - s_old) < 1e-9:
            return
    
        # Zoom hacia el puntero: mantenemos el punto bajo el cursor estable
        ox, oy = self._origin_px
        cx_mm, cy_mm = self._view_center_mm
        # coords lógicas del cursor antes del zoom
        mmx = cx_mm + (event.x - ox) / s_old
        mmy = cy_mm + (oy - event.y) / s_old
    
        # nuevo centro para mantener el mismo punto bajo el cursor
        cx_new = mmx - (event.x - ox) / s_new
        cy_new = mmy - (oy - event.y) / s_new
    
        self._plot_scale = s_new
        self._view_center_mm = (cx_new, cy_new)
        self._apply_view_change()
    
    # -------------------- Texto instrucciones seleccion ----------------------
    def _show_selection_hint(self, show=True):
        c = self.visual_canvas
        c.delete("hint")
        self._hint_visible = bool(show)
        if not show:
            return
    
        w = max(1, c.winfo_width())
        x = w - 12
        y = 12
        text = "\n".join(self._hint_lines)
    
        txt_id = c.create_text(
            x, y, anchor="ne",
            text=text,
            fill="#0b7285",
            font=("Segoe UI", 10, "bold"),
            tags=("hint", "hint_text")
        )
        x1, y1, x2, y2 = c.bbox(txt_id)
        pad_x, pad_y = 10, 8
        bg_id = c.create_rectangle(
            x1 - pad_x, y1 - pad_y, x2 + pad_x, y2 + pad_y,
            fill="#e3fafc", outline="#0b7285", width=1,
            tags=("hint", "hint_bg")
        )
        c.tag_lower(bg_id, txt_id)

    def _position_selection_hint(self):
        """Reposiciona el hint tras cambios de tamaño/redibujos (si está visible)."""
        if not self._hint_visible:
            return
        # re-crear para recolocar (más simple que mover bbox)
        self._show_selection_hint(True)
        
    # -------------------- Puntos seleccion -----------------------------------
    def _sel_start(self, mode, event):
        """Inicio común de la selección con rectángulo (set/add/remove)."""
        if not getattr(self, "_nav_enabled", False):
            return "break"
        self._sel_mode = mode
        self._sel_active = True
        self._sel_start_px = (event.x, event.y)
    
        # Crea rectángulo con color según modo
        colors = {"set": "#495057", "add": "#2b8a3e", "remove": "#c92a2a"}
        if self._sel_rect_id is not None:
            try:
                self.visual_canvas.delete(self._sel_rect_id)
            except Exception:
                pass
        self._sel_rect_id = self.visual_canvas.create_rectangle(
            event.x, event.y, event.x, event.y,
            outline=colors.get(mode, "#495057"), dash=(3, 2),
            width=1, fill="", tags=("selrect",)
        )
        return "break"

    def _on_sel_start_set(self, event):
        self._sel_start("set", event)
        return "break"
    
    def _on_sel_start_add(self, event):
        self._sel_start("add", event)
        return "break"
    
    def _on_sel_start_remove(self, event):
        self._sel_start("remove", event)
        return "break"
    
    def _on_sel_drag(self, event):
        if not self._sel_active:
            return "break"
        x0, y0 = self._sel_start_px
        self.visual_canvas.coords(self._sel_rect_id, x0, y0, event.x, event.y)
        return "break"
    
    def _on_sel_end(self, event):
        if not self._sel_active:
            return "break"
        self._sel_active = False
    
        # --- (resto de tu código tal cual) ---
        x0, y0 = self._sel_start_px
        x1, y1 = event.x, event.y
        x_min, x_max = (x0, x1) if x0 <= x1 else (x1, x0)
        y_min, y_max = (y0, y1) if y0 <= y1 else (y1, y0)
    
        x_mm_min, y_mm_max = self._from_canvas(x_min, y_min)
        x_mm_max, y_mm_min = self._from_canvas(x_max, y_max)
    
        hits = set()
        for i, (x_mm, y_mm) in enumerate(self._points_logical):
            if (x_mm_min <= x_mm <= x_mm_max) and (y_mm_min <= y_mm <= y_mm_max):
                hits.add(i)
    
        if self._sel_mode == "set":
            new_sel = hits
        elif self._sel_mode == "add":
            new_sel = self._selected_indices | hits
        elif self._sel_mode == "remove":
            new_sel = self._selected_indices - hits
        else:
            new_sel = self._selected_indices
    
        self._apply_selection(new_sel)
    
        try:
            if self._sel_rect_id is not None:
                self.visual_canvas.delete(self._sel_rect_id)
        finally:
            self._sel_rect_id = None
    
        self._notify(f"{len(self._selected_indices)} puntos seleccionados")
        return "break"
    
    def _apply_selection(self, indices: set):
        self._selected_indices = set(indices)
        self._repaint_points()
        
    # -------------------- Obtener puntos seleccionados -----------------------
    def get_selected_points(self):
        """Lista de (x_mm, y_mm) seleccionados en orden creciente de índice."""
        return [self._points_logical[i] for i in sorted(self._selected_indices)]
    
    # -------------------- BOTON CALCULAR -------------------------------------
    def _open_calc_dialog(self):
        if not self._nav_enabled:
            self._notify("Finaliza la grabación para usar las operaciones de cálculo.")
            return
    
        dlg = tk.Toplevel(self)
        dlg.title("Añadir construcciones adicionales")
        dlg.transient(self)
        dlg.resizable(False, False)
        dlg.geometry("320x240")  # tamaño ventana
    
        op = tk.StringVar(value="intersect")
    
        frm_top = ttk.Frame(dlg, padding=12)
        frm_top.grid(row=0, column=0, sticky="nsew")
    
        # Radios
        r1 = ttk.Radiobutton(
            frm_top, text="Punto de intersección (2 líneas)",
            value="intersect", variable=op)
        r2 = ttk.Radiobutton(
            frm_top, text="Línea (2 puntos)", value="centers", variable=op)
        r1.grid(row=0, column=0, sticky="w")
        r2.grid(row=1, column=0, sticky="w", pady=(4, 0))
    
        # Contenedor dinámico
        frm_dyn = ttk.Frame(frm_top)
        frm_dyn.grid(row=2, column=0, sticky="ew", pady=(10, 0))
    
        # --- Panel 'intersect' ---
        frm_i = ttk.Frame(frm_dyn)
        ttk.Label(frm_i, text="Línea 1:").grid(row=0, column=0, padx=(0, 6))
        cbo_l1 = ttk.Combobox(frm_i, state="readonly", width=24)
        cbo_l1.grid(row=0, column=1)
        ttk.Label(frm_i, text="Línea 2:").grid(row=1, column=0, padx=(0, 6), pady=(6, 0))
        cbo_l2 = ttk.Combobox(frm_i, state="readonly", width=24)
        cbo_l2.grid(row=1, column=1, pady=(6, 0))
    
        # --- Panel 'centers' (ahora usa puntos: centros + intersecciones) ---
        frm_c = ttk.Frame(frm_dyn)
        ttk.Label(frm_c, text="Punto 1:").grid(row=0, column=0, padx=(0, 6))
        cbo_p1 = ttk.Combobox(frm_c, state="readonly", width=28)
        cbo_p1.grid(row=0, column=1)
        ttk.Label(frm_c, text="Punto 2:").grid(row=1, column=0, padx=(0, 6), pady=(6, 0))
        cbo_p2 = ttk.Combobox(frm_c, state="readonly", width=28)
        cbo_p2.grid(row=1, column=1, pady=(6, 0))
    
        # Botonera
        frm_btn = ttk.Frame(frm_top)
        frm_btn.grid(row=3, column=0, sticky="e", pady=(12, 0))
        btn_ok = ttk.Button(frm_btn, text="Calcular")
        btn_close = ttk.Button(frm_btn, text="Cerrar", command=dlg.destroy)
        btn_ok.grid(row=0, column=0, padx=(0, 8))
        btn_close.grid(row=0, column=1)
    
        # ---------- utilidades -----------------------------------------------
        def _point_choices():
            """
            Devuelve lista de etiquetas de puntos disponibles:
            - '<nombre_círculo> (centro)' para cada círculo
            - '<nombre_punto>' para cada punto de intersección guardado
            """
            items = []
            for n, o in self._measure_objects.items():
                if o["type"] == "circle":
                    items.append(f"{n} (centro)")
                elif o["type"] == "point":
                    items.append(n)
            items.sort()
            return items
    
        def _resolve_point(label):
            """
            Convierte la etiqueta seleccionada en coordenadas (x, y).
            """
            if not label:
                return None
            if label.endswith("(centro)"):
                cname = label[:-8].strip()  # quita ' (centro)'
                obj = self._measure_objects.get(cname)
                if obj and obj["type"] == "circle":
                    prm = obj["params"]
                    return (prm["cx"], prm["cy"])
            else:
                obj = self._measure_objects.get(label)
                if obj and obj["type"] == "point":
                    prm = obj["params"]
                    return (prm["x"], prm["y"])
            return None
    
        def _unique_name(base):
            name = base
            i = 2
            while name in self._measure_objects:
                name = f"{base}_{i}"
                i += 1
            return name
    
        # Rellenar combos y alternar vistas
        def refresh_lists():
            lines = sorted([n for n, o in self._measure_objects.items() if o["type"] == "line"])
            pts   = _point_choices()
    
            cbo_l1["values"] = lines; cbo_l2["values"] = lines
            cbo_p1["values"] = pts;   cbo_p2["values"] = pts
    
            if lines:
                cbo_l1.current(0)
                if len(lines) > 1: cbo_l2.current(1)
            if pts:
                cbo_p1.current(0)
                if len(pts) > 1: cbo_p2.current(1)
    
        def show_panel(*_):
            for w in (frm_i, frm_c):
                w.grid_forget()
            if op.get() == "intersect":
                frm_i.grid(row=0, column=0, sticky="w")
            else:
                frm_c.grid(row=0, column=0, sticky="w")
    
        def do_calc():
            mode = op.get()
            if mode == "intersect":
                n1 = cbo_l1.get().strip()
                n2 = cbo_l2.get().strip()
                if not n1 or not n2 or n1 == n2:
                    self._notify("Selecciona dos líneas distintas.")
                    return
                res = self.intersection_of_lines(n1, n2)  # debe guardar el punto y redibujar
                if res is not None:
                    dlg.destroy()
            else:
                s1 = cbo_p1.get().strip()
                s2 = cbo_p2.get().strip()
                if not s1 or not s2 or s1 == s2:
                    self._notify("Selecciona dos puntos distintos.")
                    return
    
                p1 = _resolve_point(s1)
                p2 = _resolve_point(s2)
                if not p1 or not p2:
                    self._notify("No se pudieron resolver las coordenadas de los puntos.")
                    return
    
                # Crear línea a partir de dos puntos
                (x1, y1), (x2, y2) = p1, p2
                vx, vy = (x2 - x1), (y2 - y1)
                nrm = math.hypot(vx, vy) or 1.0
                vx, vy = vx / nrm, vy / nrm
                x0, y0 = (x1 + x2) / 2.0, (y1 + y2) / 2.0
    
                new_name = _unique_name(f"obj{self._meas_seq}")
                color = self._next_color()
                self._measure_objects[new_name] = {
                    "type": "line",
                    "color": color,
                    "params": {"x0": x0, "y0": y0, "vx": vx, "vy": vy},
                    "points": [(x1, y1), (x2, y2)]
                }
                self._redraw_measure_objects()
                self._notify(f"Línea '{new_name}' creada.")
                dlg.destroy()
    
        btn_ok.configure(command=do_calc)
        op.trace_add("write", show_panel)
    
        refresh_lists()
        show_panel()
        dlg.grab_set()   # modal
        dlg.wait_window()

    # --------- Actualizar visualizacion circulo calculado --------------------
    def _redraw_fit_circle(self):
        """Redibuja el círculo ajustado (si existe) usando la vista actual."""
        c = self.visual_canvas
        c.delete("fit")  # solo borra el overlay del círculo
    
        if not self._fit_circle_params:
            return
    
        cx_mm, cy_mm, r_mm = self._fit_circle_params
        cx_px, cy_px = self._to_canvas(cx_mm, cy_mm)
        r_px = max(1e-9, self._plot_scale) * r_mm
    
        # Círculo (verde)
        c.create_oval(cx_px - r_px, cy_px - r_px, cx_px + r_px, cy_px + r_px,
                      outline="#2f9e44", width=2, tags=("fit",))
    
        # Centro (punto verde)
        c.create_oval(cx_px - 3, cy_px - 3, cx_px + 3, cy_px + 3,
                      fill="#2f9e44", outline="", tags=("fit",))
    
        # Etiqueta con coordenadas y radio
        label = f"C=({cx_mm:.3f}, {cy_mm:.3f})  R={r_mm:.3f}"
        c.create_text(cx_px + r_px + 8, cy_px, anchor="w",
                      text=label, fill="#2f9e44",
                      font=("Segoe UI", 11, "bold"),
                      tags=("fit",))
    
    # ------------ Operaciones de calculo con OBJETOS DE MEDICION -------------
    def intersection_of_lines(self, nameA: str, nameB: str):
        LA = self._measure_objects.get(nameA)
        LB = self._measure_objects.get(nameB)
        if not LA or not LB or LA["type"] != "line" or LB["type"] != "line":
            self._notify("Ambos objetos deben ser líneas existentes.")
            return None
    
        a = LA["params"]; b = LB["params"]
        x0, y0, vx1, vy1 = a["x0"], a["y0"], a["vx"], a["vy"]
        x1, y1, vx2, vy2 = b["x0"], b["y0"], b["vx"], b["vy"]
    
        den = vx1*vy2 - vy1*vx2
        if abs(den) < 1e-12:
            self._notify("Líneas paralelas o casi paralelas.")
            return None
    
        t = ((x1 - x0)*vy2 - (y1 - y0)*vx2) / den
        xi = x0 + vx1*t
        yi = y0 + vy1*t
    
        # guarda como objeto 'point' con nombre único y color de paleta
        name = self._unique_name("PI")  # p.ej. "PI", "PI_2", ...
        color = self._next_color()
        self._measure_objects[name] = {
            "type": "point",
            "color": color,
            "params": {"x": xi, "y": yi},
            "points": []
        }
    
        self._redraw_measure_objects()
        self._notify(f"Punto de intersección '{name}' creado.")
        return (xi, yi)
        
    # --- Utilidades para listar puntos y líneas (para los desplegables SC) ---
    def _list_sc_points(self):
        """
        Devuelve nombres de puntos disponibles para el SC: 
            - centros de círculos
            - intersecciones guardadas
        """
        centers = [f"Centro: {n}" for n,o in self._measure_objects.items() if o["type"] == "circle"]
        points  = [n for n,o in self._measure_objects.items() if o["type"] == "point"]
        return sorted(centers + points)
        
    
    def _list_line_names(self):
        return sorted([n for n,o in self._measure_objects.items() if o["type"] == "line"])
    
    def _resolve_named_point(self, key: str):
        """Convierte un nombre del combo a coordenadas (x,y)."""
        key = (key or "").strip()
        # prefijo "Centro: " para círculos y, si no lo tiene, busca un objeto point
        if key.startswith("Centro:"):
            cname = key.split(":",1)[1].strip()
            obj = self._measure_objects.get(cname)
            if not obj or obj["type"] != "circle":
                raise ValueError(f"Centro no válido: {cname}")
            return obj["params"]["cx"], obj["params"]["cy"]
        else:
            obj = self._measure_objects.get(key)
            if not obj or obj["type"] != "point":
                raise ValueError(f"Punto no válido: {key}")
            return obj["params"]["x"], obj["params"]["y"]

    # --------------- Ventana “Definir SC” y lógica ---------------------------
    def _open_sc_dialog(self):
        if not self._nav_enabled:
            self._notify("Finaliza la grabación para definir el SC.")
            return
    
        dlg = tk.Toplevel(self)
        dlg.title("Definir sistema de coordenadas (XY)")
        dlg.transient(self)
        dlg.resizable(False, False) 
        dlg.geometry("480x200")  # tamaño ventana
    
        # --- Vars inicializadas desde el último estado guardado ---
        st = self._coordsys_state
        var_mode   = tk.StringVar(value=st.get("mode", "2pts"))  # "2pts" | "ref_line"
        var_origin = tk.StringVar(value=st.get("origin_key") or "")
        var_refpt  = tk.StringVar(value=st.get("ref_key") or "")  # usado si 2pts
        var_line   = tk.StringVar(value=st.get("ref_key") or "")  # usado si ref_line
        var_invert = tk.BooleanVar(value=bool(st.get("invert", False)))
    
        frm = ttk.Frame(dlg, padding=12)
        frm.grid(row=0, column=0, sticky="nsew")
    
        # --- Origen ---
        row = 0
        ttk.Label(frm, text="Origen:").grid(row=row, column=0, sticky="w", padx=(0,6))
        cbo_origin = ttk.Combobox(frm, state="readonly", width=28, textvariable=var_origin)
        cbo_origin.grid(row=row, column=1, sticky="w")
        row += 1
    
        # --- Modo ---
        ttk.Label(frm, text="Alinear eje X por:").grid(row=row, column=0, sticky="w", padx=(0,6), pady=(8,0))
        frm_mode = ttk.Frame(frm)
        frm_mode.grid(row=row, column=1, sticky="w", pady=(8,0))
        ttk.Radiobutton(frm_mode, text="Dos puntos", value="2pts",     variable=var_mode).grid(row=0, column=0, padx=(0,12))
        ttk.Radiobutton(frm_mode, text="Línea de referencia", value="ref_line", variable=var_mode).grid(row=0, column=1)
        row += 1
    
        # --- Panel dinámico ---
        frm_dyn = ttk.Frame(frm)
        frm_dyn.grid(row=row, column=0, columnspan=2, sticky="ew", pady=(8,0))
        row += 1
    
        # Contenido para "2pts"
        frm_2pts = ttk.Frame(frm_dyn)
        ttk.Label(frm_2pts, text="Punto de referencia:").grid(row=0, column=0, padx=(0,6))
        cbo_refpt = ttk.Combobox(frm_2pts, state="readonly", width=28, textvariable=var_refpt)
        cbo_refpt.grid(row=0, column=1, sticky="w")
    
        # Contenido para "ref_line"
        frm_rl = ttk.Frame(frm_dyn)
        ttk.Label(frm_rl, text="Línea de referencia:").grid(row=0, column=0, padx=(0,6))
        cbo_line = ttk.Combobox(frm_rl, state="readonly", width=28, textvariable=var_line)
        cbo_line.grid(row=0, column=1, sticky="w")
        chk_inv = ttk.Checkbutton(frm_rl, text="Dirección contraria", variable=var_invert)
        chk_inv.grid(row=0, column=2, padx=(10,0))
    
        # --- Botonera ---
        frm_btn = ttk.Frame(frm)
        frm_btn.grid(row=row, column=0, columnspan=2, sticky="e", pady=(12,0))
        btn_preview = ttk.Button(frm_btn, text="Preview")
        btn_close   = ttk.Button(frm_btn, text="Cerrar")
        btn_preview.grid(row=0, column=0, padx=(0,8))
        btn_close.grid(row=0, column=1)

        # --- Relleno y selección por defecto (respetando estado anterior) ---
        def refresh_lists(*_):
            points = self._list_sc_points()
            lines  = self._list_line_names()
    
            cbo_origin["values"] = points
            # origen recordado o primer punto
            if var_origin.get() in points:
                cbo_origin.set(var_origin.get())
            elif points:
                cbo_origin.current(0); var_origin.set(cbo_origin.get())
    
            # panel visible
            for w in (frm_2pts, frm_rl):
                w.grid_forget()
    
            if var_mode.get() == "2pts":
                frm_2pts.grid(row=0, column=0, sticky="w")
                # ref points: todos menos el origen
                pts2 = [p for p in points if p != var_origin.get()]
                cbo_refpt["values"] = pts2
                if var_refpt.get() in pts2:
                    cbo_refpt.set(var_refpt.get())
                elif pts2:
                    cbo_refpt.current(0); var_refpt.set(cbo_refpt.get())
            else:
                frm_rl.grid(row=0, column=0, sticky="w")
                cbo_line["values"] = lines
                if var_line.get() in lines:
                    cbo_line.set(var_line.get())
                elif lines:
                    cbo_line.current(0); var_line.set(cbo_line.get())
    
        def _apply_from_ui(save_state: bool):
            """Calcula/actualiza el SC desde los controles. Si save_state=True, persiste selección."""
            try:
                Ox, Oy = self._resolve_named_point(var_origin.get())
            except Exception as e:
                self._notify(str(e)); return False
    
            if var_mode.get() == "2pts":
                try:
                    Px, Py = self._resolve_named_point(var_refpt.get())
                except Exception as e:
                    self._notify(str(e)); return False
                vx, vy = (Px - Ox), (Py - Oy)
            else:
                name = var_line.get()
                L = self._measure_objects.get(name)
                if not L or L["type"] != "line":
                    self._notify("Selecciona una línea válida."); return False
                vx, vy = L["params"]["vx"], L["params"]["vy"]
                if var_invert.get():
                    vx, vy = -vx, -vy
    
            n = math.hypot(vx, vy)
            if n < 1e-12:
                self._notify("Vector X nulo."); return False
            ux, uy = vx/n, vy/n
            ang = math.degrees(math.atan2(uy, ux))
    
            # guarda parámetros y dibuja
            self._coordsys_params = {"origin": (Ox, Oy), "ux": ux, "uy": uy, "A": ang}
            self._redraw_coordsys()
    
            if save_state:
                self._coordsys_state.update({
                    "origin_key": var_origin.get(),
                    "mode": var_mode.get(),
                    "ref_key": (var_refpt.get() if var_mode.get()=="2pts" else var_line.get()),
                    "invert": bool(var_invert.get())
                })
            return True
    
        btn_preview.configure(command=lambda: _apply_from_ui(save_state=True))
        btn_close.configure(command=lambda: (_apply_from_ui(save_state=True), dlg.destroy()))
    
        # Eventos que alteran listas
        var_mode.trace_add("write", refresh_lists)
        var_origin.trace_add("write", refresh_lists)
    
        refresh_lists()
        dlg.grab_set()
        dlg.wait_window()

    # --------- Dibujo del Sistema de Coordenadas (overlay “coordsys”) --------
    def _redraw_coordsys(self):
        """Dibuja el sistema de coordenadas definido (origen + eje X con flecha) y etiqueta abajo."""
        c = self.visual_canvas
        c.delete("coordsys")        # limpia todo lo del SC
        p = getattr(self, "_coordsys_params", None)
        if not p:
            return
    
        Ox, Oy = p["origin"]
        ux, uy = p["ux"], p["uy"]
        A_deg  = p["A"]
    
        # Origen y punta del eje X (longitud en mm calculada según viewport)
        w = max(1, c.winfo_width()); h = max(1, c.winfo_height())
        m = self._plot_margin
        usable_px = max(1, min(w, h) - 2*m)
        L_mm = max(40.0, usable_px / max(1e-9, self._plot_scale) * 0.35)  # ~35% del lado corto
    
        x_tip = Ox + ux * L_mm
        y_tip = Oy + uy * L_mm
    
        # A canvas
        oxp, oyp   = self._to_canvas(Ox, Oy)
        xtp, ytp   = self._to_canvas(x_tip, y_tip)
    
        # Eje X con flecha
        c.create_line(oxp, oyp, xtp, ytp, width=2, arrow="last", tags=("coordsys",))
        # Marca del origen (cruz pequeña)
        c.create_line(oxp-6, oyp, oxp+6, oyp, width=2, tags=("coordsys",))
        c.create_line(oxp, oyp-6, oxp, oyp+6, width=2, tags=("coordsys",))
        c.create_text(oxp+8, oyp-8, anchor="nw", text="O", font=("Segoe UI", 10, "bold"),
                      tags=("coordsys",))
    
        # HUD abajo (anclado al viewport)
        c.create_text(12, h-16, anchor="sw",
                      text=f"SC: O=({Ox:.3f}, {Oy:.3f})  A={A_deg:.3f}°",
                      font=("Segoe UI", 10, "bold"),
                      tags=("coordsys",))

    # --------- Bucle de grabación --------------------------------------------
    def _record_loop(self):
        """
        Mientras self._recording sea True:
        - Lee variables OPC UA (dict 'variables')
        - Lee LabJack (names_lj)
        - Calcula AngleDirVSpin y AngleValVSpin
        - Calcula punto de contacto proyectando la deflexión sobre la normal
          a la trayectoria.
        Guarda el último paquete en self.last_sample.
        """
        client = self.conn.opc_client
        handle = self.conn.lj_handle
        if client is None or handle is None:
            self._notify("Conexiones no disponibles")
            self._recording = False
            self.after(0, lambda: self.btn_grabar.configure(text="GRABAR"))
            return
    
        period = 0.02  # 20 Hz
        pos_previa = None
        MIN_MOV_MM2 = 1e-4  # umbral de movimiento al cuadrado (p.ej. 0.01 mm^2)
    
        while self._recording:
            paquete = {}
    
            # --- KDC (OPC UA) ---
            for nombre, nodo_variable in variables.items():
                try:
                    with self.conn._opc_io:
                        nodo = client.get_node(nodo_variable)
                        paquete[nombre] = nodo.get_value()
                except Exception as e:
                    if DEBUG:
                        print(f"Error al leer {nombre}: {e}")
                    paquete[nombre] = "ERROR"
    
            # --- LabJack T4 ---
            if ljm is not None:
                try:
                    with self.conn._lj_io:
                        valores_lj = ljm.eReadNames(handle, len(names_lj), names_lj)
                    dir_mA = volts_to_mA(valores_lj[0], clip=True)
                    mod_mA = volts_to_mA(valores_lj[1], clip=True)
                    paquete["AngleDirVSpin"] = (22.5 * dir_mA - 90.0 + 360.0) % 360.0
                    paquete["AngleValVSpin"] = 0.25 * mod_mA - 1.0
                except Exception as e:
                    if DEBUG:
                        print(f"Error lectura LabJack: {e}")
                    paquete["AngleDirVSpin"] = "ERROR"
                    paquete["AngleValVSpin"] = "ERROR"
            else:
                paquete["AngleDirVSpin"] = "ERROR"
                paquete["AngleValVSpin"] = "ERROR"
    
            # --------- Posicion actual ---------------------------------------
            try:
                val = float(paquete.get("AngleValVSpin", -1e9))
                x_act = float(paquete.get("X", "nan"))
                y_act = float(paquete.get("Y", "nan"))
            except Exception:
                val = -1e9
                x_act = y_act = float("nan")
        
            if np.isfinite(x_act) and np.isfinite(y_act):
                pos_actual = np.array([x_act, y_act], dtype=float)
            else:
                print("Posicion actual no registrada.")
            # --------- Punto de contacto -------------------------------------
            if pos_previa is not None:
                dx = pos_actual[0] - pos_previa[0]
                dy = pos_actual[1] - pos_previa[1]
    
                # if dx*dx + dy*dy > MIN_MOV_MM2 and val >= MIN_ANGLEVAL_FOR_POINT:
                if dx*dx + dy*dy > MIN_MOV_MM2:
                    # MISMO ángulo de trayectoria que en el script viejo
                    ang_tray_deg = getAngle(pos_previa, pos_actual)
    
                    try:
                        contact_point = calcular_pto_contacto_normal(
                            paquete,
                            longitud_pivotamiento,
                            ang_tray_deg,
                            sentido_horario  # True si recorres la pieza horario
                        )
                        x_pieza, y_pieza = contact_point
                        # self.after(0, lambda x=x_pieza, y=y_pieza: self._add_point(x, y))
                        self.after(0, lambda x=x_act, y=y_act: self._add_point(x, y))
                    except Exception as e:
                        if DEBUG:
                            print("Error calculando punto de contacto:", e)
    
            pos_previa = pos_actual.copy()
    
            # .............................................................
            # Guarda último paquete
            with self._rec_lock:
                self.last_sample = paquete
    
            # HUD ligero
            self._render_last_sample_light()
    
            time.sleep(period)


    def _render_last_sample_light(self):
        """Muestra los datos con 3 decimales en un HUD legible."""
        try:
            with self._rec_lock:
                sample = dict(self.last_sample) if self.last_sample else None
            if not sample:
                return
    
            def fmt(v):
                try:    return f"{float(v):.3f}"
                except: return str(v)
    
            x_txt  = fmt(sample.get("X", "—"))
            y_txt  = fmt(sample.get("Y", "—"))
            A_txt  = fmt(sample.get("A", "—"))
            dir_txt = fmt(sample.get("AngleDirVSpin", "—"))
            val_txt = fmt(sample.get("AngleValVSpin", "—"))
    
            lines = [
                f"X: {x_txt} mm",
                f"Y: {y_txt} mm",
                f"A: {A_txt}°",
                f"Dirección: {dir_txt}°",
                f"Módulo: {val_txt}"
            ]
    
            c = self.visual_canvas
            c.delete("placeholder")  # si aún estaba el texto inicial
            c.delete("hud")          # borra HUD anterior (no puntos, no ejes)
            txt_id = c.create_text(
                12, 12, anchor="nw",
                text="\n".join(lines),
                fill="#212529",
                font=("Segoe UI", 12, "bold"),
                tags=("hud", "hud_text")
            )
            # Caja de fondo (tipo "transparente" con stipple)
            x1, y1, x2, y2 = c.bbox(txt_id)
            pad_x, pad_y = 8, 6
            bg_id = c.create_rectangle(
                x1 - pad_x, y1 - pad_y, x2 + pad_x, y2 + pad_y,
                fill="#ffffff", outline="#adb5bd", width=1,
                stipple="gray25",  # efecto translúcido
                tags=("hud", "hud_bg")
            )
            c.tag_lower(bg_id, txt_id)   # fondo detrás del texto
        except Exception:
            if DEBUG:
                import traceback; traceback.print_exc()
                
    # ------------------ BEST-FIT ---------------------------------------------
    # devuelve centroides (medido y referencia)
    def _rigid_fit_2d(self, P, Q):
        """
        Kabsch 2D: encuentra R,t que lleva P->Q (sin escala).
        Devuelve:
          {
            "theta_deg": float, "tx": float, "ty": float,
            "origin_meas": (px,py),   # centro del conjunto medido (A)
            "origin_ref":  (qx,qy),   # centro del conjunto referencia (B)
            "rmse": float,
            "R": R(2x2),
            "t": t(2,)
          }
        """
        P = np.asarray(P, float)
        Q = np.asarray(Q, float)
        if len(P) < 3 or len(Q) < 3:
            raise ValueError("Se necesitan ≥3 puntos en cada conjunto.")
        p0 = P.mean(axis=0)
        q0 = Q.mean(axis=0)
        Pc = P - p0
        Qc = Q - q0
        H = Pc.T @ Qc
        U, S, Vt = np.linalg.svd(H)
        R = Vt.T @ U.T
        if np.linalg.det(R) < 0:
            Vt[-1, :] *= -1
            R = Vt.T @ U.T
        t = q0 - R @ p0
        theta = math.degrees(math.atan2(R[1,0], R[0,0]))
        P_fit = (R @ P.T).T + t
        err = np.sqrt(np.mean(np.sum((P_fit - Q)**2, axis=1)))
        return {
            "theta_deg": float(theta),
            "tx": float(t[0]), "ty": float(t[1]),
            "origin_meas": (float(p0[0]), float(p0[1])),
            "origin_ref":  (float(q0[0]), float(q0[1])),
            "rmse": float(err),
            "R": R, "t": t
        }
       
    
    # --- Helpers para ICP -------------------------------------------------
    def _rot_from_deg(self, deg: float) -> np.ndarray:
        th = math.radians(deg); c, s = math.cos(th), math.sin(th)
        return np.array([[c, -s],[s, c]], float)
    
    def _kabsch_2d(self, P, Q):
        P = np.asarray(P, float); Q = np.asarray(Q, float)
        p0 = P.mean(axis=0); q0 = Q.mean(axis=0)
        Pc = P - p0; Qc = Q - q0
        H = Pc.T @ Qc
        U, S, Vt = np.linalg.svd(H)
        R = Vt.T @ U.T
        if np.linalg.det(R) < 0:
            Vt[-1, :] *= -1
            R = Vt.T @ U.T
        t = q0 - R @ p0
        Pf = (R @ P.T).T + t
        rmse = float(np.sqrt(np.mean(np.sum((Pf - Q)**2, axis=1))))
        th = float(math.degrees(math.atan2(R[1,0], R[0,0])))
        return R, t, th, rmse
    
    def _nn_indices_gated(self, P, Q, gate_mm=BESTFIT_NN_GATE_MM, chunk=2000):
        """Vecinos más cercanos con puerta de distancia (dentro de gate)."""
        P = np.asarray(P, float); Q = np.asarray(Q, float)
        gate2 = gate_mm * gate_mm
        n = len(P)
        idx = np.full(n, -1, dtype=int)
        ok  = np.zeros(n, dtype=bool)
        for i in range(0, n, chunk):
            X = P[i:i+chunk]
            d2 = ((X[:,None,:] - Q[None,:,:])**2).sum(axis=2)
            j  = np.argmin(d2, axis=1)
            m  = d2[np.arange(len(j)), j] <= gate2
            idx[i:i+chunk] = j
            ok[i:i+chunk]  = m
        return idx, ok
    
    def _uniform_down(self, P, m):
        if len(P) <= m: return np.asarray(P, float)
        idx = np.linspace(0, len(P)-1, m).round().astype(int)
        return np.asarray(P, float)[idx]
    
    def _coarse_init(self, A, B,
                     theta_max=BESTFIT_MAX_ROT_DEG,
                     step=BESTFIT_COARSE_STEP,
                     gate=BESTFIT_NN_GATE_MM,
                     trim=BESTFIT_TRIM):
        """
        Barrido angular en [-theta_max, +theta_max] sin depender del orden:
        para cada θ, rota A sobre su centroide, traduce por la mediana de
        (B_nn - A_rot) y evalúa un coste robusto (trimmed).
        Devuelve R, t y θ en grados.
        """
        A = np.asarray(A, float); B = np.asarray(B, float)
        if len(A) > BESTFIT_MAX_SAMPLES: A = self._uniform_down(A, BESTFIT_MAX_SAMPLES)
        if len(B) > BESTFIT_MAX_SAMPLES: B = self._uniform_down(B, BESTFIT_MAX_SAMPLES)
    
        ca = A.mean(axis=0)
        best = None
        ths = np.arange(-theta_max, theta_max + 1e-12, step)
        for th in ths:
            R = self._rot_from_deg(th)
            Arot = (A - ca) @ R.T  # rotar alrededor del centroide de A
            j, ok = self._nn_indices_gated(Arot, B, gate)
            if ok.sum() < 3:
                continue
            resid_vec = B[j[ok]] - Arot[ok]
            t_local = np.median(resid_vec, axis=0)  # robusto
            # coste robusto con trimmed
            r = (Arot[ok] + t_local) - B[j[ok]]
            d2 = np.sum(r*r, axis=1)
            k  = max(3, int(trim*len(d2)))
            cost = float(np.mean(np.partition(d2, k-1)[:k]))
            # convertir a forma x' = R @ x + t_total
            t_total = t_local - R @ ca
            if best is None or cost < best[0]:
                best = (cost, R, t_total, th)
        if best is None:
            # fallback grosero: centrar por centroides sin rotación
            R = self._rot_from_deg(0.0)
            t = B.mean(axis=0) - A.mean(axis=0)
            return R, t, 0.0
        _, R, t, th = best
        return R, t, th
    
    def _compute_bestfit(self,
                         trim=BESTFIT_TRIM,
                         gate=BESTFIT_NN_GATE_MM,
                         iters=BESTFIT_ICP_ITERS,
                         tol=BESTFIT_ICP_TOL):
        """
        Alinea TODOS los puntos grabados (self._points_logical) a la referencia
        transformada (self._transform_ref_points(self._ref_points_raw)):
          - inicialización por barrido angular (no depende del orden)
          - ICP recortado con puerta de distancia
          - giro global restringido a ±BESTFIT_MAX_ROT_DEG
        Devuelve dict {R, t, theta_deg, rmse, pivot}
        """
        if not self._points_logical or not getattr(self, "_ref_points_raw", None):
            return None
        if len(self._points_logical) < 3 or len(self._ref_points_raw) < 3:
            return None
    
        # Conjuntos completos (NO se ignoran “líneas sueltas”)
        A_full = np.asarray(self._points_logical, float)
        B_full = np.asarray(self._transform_ref_points(self._ref_points_raw), float)
    
        # Inicialización robusta (barrido angular)
        R, t, _ = self._coarse_init(A_full, B_full)
    
        # Muestreos para ICP (rendimiento)
        A = self._uniform_down(A_full, BESTFIT_ICP_SAMPLES)
        B = self._uniform_down(B_full, BESTFIT_ICP_SAMPLES)
    
        prev_err = None
        for _ in range(iters):
            A_hat = (A @ R.T) + t
            j, ok = self._nn_indices_gated(A_hat, B, gate)
            sel = np.where(ok)[0]
            if len(sel) < 3:
                break
    
            # Recorte robusto
            d2 = np.sum((A_hat[sel] - B[j[sel]])**2, axis=1)
            k  = max(3, int(trim*len(sel)))
            keep = sel[np.argsort(d2)[:k]]
    
            # Ajuste global sobre los pares (originales, no transformados)
            R_fit, t_fit, th_fit, _ = self._kabsch_2d(A[keep], B[j[keep]])
            # Limitar ángulo global
            th_fit = max(-BESTFIT_MAX_ROT_DEG, min(BESTFIT_MAX_ROT_DEG, th_fit))
            R = self._rot_from_deg(th_fit)
            # t consistente con el R clipeado (usando centroides de los pares)
            p0 = A[keep].mean(axis=0); q0 = B[j[keep]].mean(axis=0)
            t  = q0 - R @ p0
    
            # Error actual (RMSE sobre inliers recortados)
            A_hat = (A @ R.T) + t
            d2 = np.sum((A_hat[keep] - B[j[keep]])**2, axis=1)
            err = float(np.sqrt(np.mean(d2)))
            if prev_err is not None and abs(prev_err - err) < tol:
                break
            prev_err = err
    
        theta_deg = float(math.degrees(math.atan2(R[1,0], R[0,0])))
        # Origen de rotación (pivot): (I - R) p = t
        I = np.eye(2)
        pivot = None
        try:
            if abs(theta_deg) > 1e-9:
                p = np.linalg.solve(I - R, t)
                pivot = (float(p[0]), float(p[1]))
        except Exception:
            pivot = None
    
        return {"R": R, "t": t, "theta_deg": theta_deg, "rmse": float(prev_err if prev_err is not None else 0.0), "pivot": pivot}
    
    def _show_bestfit_overlay(self):
        """Pinta un recuadro con el BEST-FIT (A->B) y la corrección de BASE Tba (B->A)."""
        c = self.visual_canvas
        c.delete("bestfit")
        bf = getattr(self, "_bestfit_result", None)
        if not bf:
            return
    
        # A->B (fit directo sobre lo que se muestra)
        dx  = bf.get("tx", 0.0)
        dy  = bf.get("ty", 0.0)
        th  = bf.get("theta_deg", 0.0)
        omx, omy = bf.get("origin", (0.0, 0.0))
    
        # B->A (Tba) – corrección que aplicarías a la BASE actual
        dx_ba = bf.get("ba_dx", 0.0)
        dy_ba = bf.get("ba_dy", 0.0)
        th_ba = bf.get("ba_theta_deg", 0.0)
        obx, oby = bf.get("ba_origin", (0.0, 0.0))
    
        # text = (f"BEST-FIT (A→B):  ΔX={dx:.3f}  ΔY={dy:.3f}  θ={th:.3f}°  "
        #         f"Origen=({omx:.3f}, {omy:.3f})\n"
        #         f"Tba (B→A, aplicar a BASE):  ΔX={dx_ba:.3f}  ΔY={dy_ba:.3f}  θ={th_ba:.3f}°  "
        #         f"Origen=({obx:.3f}, {oby:.3f})")
        
        text = (f"Δ_BASE:  [X={dx_ba:.3f},  Y={dy_ba:.3f},  A={th_ba:.3f}]")
    
        w = max(1, c.winfo_width()); h = max(1, c.winfo_height())
        tid = c.create_text(12, h - 12, anchor="sw", text=text,
                            fill="#212529", font=("Segoe UI", 11, "bold"),
                            tags=("bestfit", "bestfit_text"))
        x1, y1, x2, y2 = c.bbox(tid)
        bg = c.create_rectangle(x1 - 8, y1 - 6, x2 + 8, y2 + 6,
                                fill="#f8f9fa", outline="#adb5bd", width=1,
                                tags=("bestfit", "bestfit_bg"))
        c.tag_lower(bg, tid)
    
    # --------- Conexiones iniciales ----------
    def _connect_on_start(self):
        threading.Thread(target=self._try_connect_opc, daemon=True).start()
        threading.Thread(target=self._try_connect_labjack, daemon=True).start()

    def _try_connect_opc(self):
        self._notify("Conectando KDC...")
        try:
            self.conn.connect_opc(timeout=4)
            print("[OPC] Conectado")
        except Exception as e:
            print("[OPC] No se pudo conectar:", e)

    def _try_connect_labjack(self):
        self._notify("Conectando LabJack...")
        try:
            self.conn.connect_labjack()
        except Exception as e:
            print("[LJ] No se pudo conectar:", e)

    # --------- UI tick (lee flags y pinta) ----------
    def _ui_tick(self):
        """
        BOTON "Grabar" habilitado sólo si ambos conectados o si ya estamos grabando (para poder parar)
        BOTON "Definir SC" habilitado siempre que NO estés grabando)
        """

        opc_ok, lj_ok = self.conn.get_status()
        self.status_kdc.set_connected(bool(opc_ok))
        self.status_vspin.set_connected(bool(lj_ok))
    
        can_record = bool(opc_ok and lj_ok)
        if not self._recording:
            self.btn_grabar.state(["!disabled"] if can_record else ["disabled"])
            # SC habilitado aunque no haya conexión (trabaja con datos ya capturados)
            try:
                self.btn_sc.state(["!disabled"])
            except Exception:
                pass
            # .................................................................
            # BEST-FIT: necesita puntos grabados y referencia cargada
            can_fit = (len(self._points_logical) >= 3) and bool(getattr(self, "_ref_points_raw", None))
            try:
                self.btn_fit.state(["!disabled"] if can_fit else ["disabled"])
            except Exception:
                pass
        else:
            # Siempre permitir detener
            self.btn_grabar.state(["!disabled"])
            # SC deshabilitado mientras grabas
            try:
                self.btn_sc.state(["disabled"])
                self.btn_fit.state(["disabled"])
            except Exception:
                pass
            
        # Habilitar BEST-FIT si hay puntos grabados (≥3)
        try:
            has_points = len(self._points_logical) >= 3
            self.btn_fit.state(["!disabled"] if has_points else ["disabled"])
        except Exception:
            pass

    
        self.after(self._ui_period_ms, self._ui_tick)

    # --------- Utilidades visuales ----------
    def _notify(self, text: str):
        self.visual_canvas.delete("overlay_text")
        w = max(1, self.visual_canvas.winfo_width())
        h = max(1, self.visual_canvas.winfo_height())
        self.visual_canvas.create_text(
            w // 2, h // 2,
            text=text,
            fill="#111",
            font=("Segoe UI", 24, "bold"),
            tags="overlay_text"
        )
        self.after(900, lambda: self.visual_canvas.delete("overlay_text"))

    # --------- Cierre limpio ----------
    def _on_close(self):
        try:
            self._recording = False
            self.conn.stop()
            self._notify("Cerrando conexiones...")
            self.conn.close_all()
        finally:
            self.destroy()

# ======================= MAIN =================================================
if __name__ == "__main__":
    version_software = 2.5
    app = App(version_software)
    app.mainloop()
