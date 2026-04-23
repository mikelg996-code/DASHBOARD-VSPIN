# -*- coding: utf-8 -*-
"""
Created on Mon Nov 10 15:22:52 2025

@author: Mikel
"""

import numpy as np
import pandas as pd

# ========== Kabsch/Procrustes 2D ==========
def rigid_transform_2d(XY_ref, XY_obs, *, grados=True):
    """
    Estima la rotación + traslación que lleva XY_ref -> XY_obs (correspondencia por índice).
    Devuelve delta_x, delta_y, delta_ang (CCW) y un dict con R, t, rmse, n.
    """
    P = np.asarray(XY_ref, dtype=float)   # teórico
    Q = np.asarray(XY_obs, dtype=float)   # real
    if P.shape != Q.shape or P.ndim != 2 or P.shape[1] != 2:
        raise ValueError("Entradas deben ser (N,2) y con el mismo N.")
    if P.shape[0] < 2:
        raise ValueError("Se requieren al menos 2 puntos.")

    cP = P.mean(axis=0); cQ = Q.mean(axis=0)
    Pc = P - cP; Qc = Q - cQ

    H = Pc.T @ Qc
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T

    t = cQ - R @ cP
    ang_rad = np.arctan2(R[1, 0], R[0, 0])
    ang = np.degrees(ang_rad) if grados else ang_rad

    P_pred = (R @ P.T).T + t
    rmse = float(np.sqrt(np.mean(np.sum((Q - P_pred) ** 2, axis=1))))
    return float(t[0]), float(t[1]), float(ang), {"R": R, "t": t, "rmse": rmse, "n": int(P.shape[0]), "pred": P_pred}

# ========== Lectura desde Excel ==========
def cargar_desde_excel(ruta="COORDENADAS.xlsx", hoja=0):
    """
    Lee el layout de tu captura:
      Col B,C -> REAL X,Y
      Col E,F -> TEORICO X,Y
    Soporta encabezados duplicados (X, Y, X.1, Y.1) o selección por posición.
    """
    df = pd.read_excel(ruta, sheet_name=hoja, header=0, usecols="A:F")

    cols = list(df.columns)

    # Intento 1: nombres duplicados X, Y y X.1, Y.1
    if {"X", "Y"}.issubset(cols) and {"X.1", "Y.1"}.issubset(cols):
        real_xy = df[["X", "Y"]].to_numpy(dtype=float)
        teo_xy  = df[["X.1", "Y.1"]].to_numpy(dtype=float)
    # Intento 2: por posición (B,C) y (E,F)
    else:
        if df.shape[1] < 6:
            raise ValueError("Se esperaban 6 columnas (A:F) según la captura.")
        real_xy = df.iloc[:, [1, 2]].to_numpy(dtype=float)
        teo_xy  = df.iloc[:, [4, 5]].to_numpy(dtype=float)

    # Filtra filas con NaN en cualquiera de los cuatro valores
    mask = np.all(np.isfinite(np.hstack([real_xy, teo_xy])), axis=1)
    real_xy = real_xy[mask]
    teo_xy  = teo_xy[mask]

    if len(real_xy) < 2:
        raise ValueError("Quedan menos de 2 puntos válidos tras filtrar NaN.")

    return teo_xy, real_xy  # (teórico, real)

# ========== Script principal ==========
if __name__ == "__main__":
    ruta_excel = "COORDENADAS.xlsx"   # ajusta si está en otra carpeta
    hoja = 0                          # o el nombre de la hoja

    XY_teo, XY_real = cargar_desde_excel(ruta_excel, hoja)
    dx, dy, dang, info = rigid_transform_2d(XY_teo, XY_real, grados=True)

    print("=== Transformación TEÓRICO → REAL (correspondencia por índice) ===")
    print(f"delta_x = {dx:.3f} mm")
    print(f"delta_y = {dy:.3f} mm")
    print(f"delta_ang = {dang:.3f} ° (CCW)")

    print("\nMatriz R y vector t:")
    print(info["R"])
    print(info["t"])

    # Diagnóstico de error
    pred = info["pred"]
    err = np.linalg.norm(XY_real - pred, axis=1)
    print("\nErrores punto a punto (mm):")
    for i, e in enumerate(err, 1):
        print(f"  Punto {i}: {e:.3f}")
    print(f"\nRMSE = {info['rmse']:.3f} mm  (n = {info['n']})")

    # Nota: Para aplicar la inversa (REAL → TEÓRICO): X_teo ≈ R.T @ (X_real - t)
