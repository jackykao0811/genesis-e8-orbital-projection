#!/usr/bin/env python3
"""
Hetu 8D — E8 終極閉合協議（數值原型）

- **E8 Cartan 矩陣**（rank 8）→ **inv(A)** 之欄即 **基本權重**（單縛格 Gram=A）。
- 殘差 **0.08（黎曼）**、**0.04（電子雲代理）** 置于 ℝ⁸ 前兩維；第 **8** 維標量 **s** 由 ⟨r, ω₈⟩＝0 解出。
- 軌道密度：**ρ₈ = ‖Q v‖²**，**v∈ℝ⁸** 為通道振幅，**Q = polar(A)** 之正交因子。
- **η**（第 8 維增益）經網格搜尋（預先算好 ρ_ref, θ, φ，僅重算 ρ_proj）。

環境變數：`HETU_E8_ETA_STEPS`（預設 21）控制 η 網格細度。
"""

from __future__ import annotations

import json
import math
import os
import re
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault(
    "MPLCONFIGDIR",
    str(Path(__file__).resolve().parents[1] / ".mplconfig_writable"),
)

import matplotlib

matplotlib.use("Agg")

import numpy as np
from scipy.linalg import polar
from scipy.special import eval_genlaguerre, factorial

LAM = 1.212375
GRID_N = 256
R_MAX_AU = 24.0

RIEMANN_RESIDUAL_E8 = 0.08
ELECTRON_CLOUD_RESIDUAL_E8 = 0.04

USE_LEGACY_7D_ONLY = False
ETA_MAX = 0.06


def e8_cartan_matrix() -> np.ndarray:
    return np.array(
        [
            [2, -1, 0, 0, 0, 0, 0, 0],
            [-1, 2, -1, 0, 0, 0, 0, 0],
            [0, -1, 2, -1, 0, 0, 0, 0],
            [0, 0, -1, 2, -1, 0, 0, 0],
            [0, 0, 0, -1, 2, -1, 0, -1],
            [0, 0, 0, 0, -1, 2, -1, 0],
            [0, 0, 0, 0, 0, -1, 2, 0],
            [0, 0, 0, 0, -1, 0, 0, 2],
        ],
        dtype=np.float64,
    )


def fundamental_weights_columns() -> np.ndarray:
    return np.linalg.inv(e8_cartan_matrix())


def e8_residual_closure_linear(r_r: float, r_e: float) -> dict[str, float | list[float]]:
    inv = fundamental_weights_columns()
    w8 = inv[:, 7]
    denom = w8[7]
    s = -(w8[0] * r_r + w8[1] * r_e) / denom
    vec = np.array([r_r, r_e, 0.0, 0.0, 0.0, 0.0, 0.0, s], dtype=np.float64)
    return {
        "s_dim8": float(s),
        "dot_r_with_omega8": float(np.dot(vec, w8)),
        "quadratic_norm_r_R8": float(np.dot(vec, vec)),
        "omega8_coeffs": w8.tolist(),
        "vector_r": vec.tolist(),
    }


def lambda_unification_candidates(lam_target: float) -> dict[str, float]:
    a = e8_cartan_matrix()
    inv = np.linalg.inv(a)
    eig = np.linalg.eigvalsh(a)
    trace_inv = float(np.trace(inv))
    lam_rat = 9699 / 8000
    return {
        "lambda_genesis_float": lam_target,
        "lambda_rational_9699_over_8000": lam_rat,
        "trace_inv_A_over_rank": trace_inv / 8.0,
        "frobenius_sq_inv_over_64": float(np.sum(inv * inv)) / 64.0,
        "min_eigenvalue_A": float(np.min(eig)),
        "max_eigenvalue_A": float(np.max(eig)),
        "abs_diff_lambda_trace_metric": abs(trace_inv / 8.0 - lam_target),
    }


def orthogonal_polar_from_cartan() -> np.ndarray:
    q, _ = polar(e8_cartan_matrix())
    return q


def phase_dwell_modulation(theta: np.ndarray, phi: np.ndarray) -> np.ndarray:
    u = np.sin(theta) * np.cos(phi)
    v = np.sin(theta) * np.sin(phi)
    w = np.cos(theta)
    omega = LAM * (u + v + w)
    dwell_e67 = np.cos(omega / LAM) ** 2 * np.sin(LAM * theta) ** 2 * np.cos(LAM * phi) ** 2
    dwell_core = 1.0 + LAM ** (-9) * dwell_e67
    return np.clip(dwell_core, 1e-30, None)


def radial_hydrogen(r: np.ndarray, n: int, ell: int) -> np.ndarray:
    r = np.asarray(r, dtype=np.float64)
    out = np.zeros_like(r)
    rho = 2.0 * r / n
    mask = r > 1e-14
    rr = r[mask]
    rrho = rho[mask]
    lag = eval_genlaguerre(n - ell - 1, 2 * ell + 1, rrho)
    coeff = math.sqrt(
        (2.0 / n) ** 3 * float(factorial(n - ell - 1)) / (2.0 * n * float(factorial(n + ell)))
    )
    out[mask] = coeff * np.exp(-rrho / 2.0) * (rrho ** ell) * lag
    return out


def angular_real_y_l0(theta: np.ndarray, ell: int) -> np.ndarray:
    if ell == 0:
        return np.full_like(theta, 1.0 / np.sqrt(4.0 * np.pi))
    if ell == 1:
        return np.sqrt(3.0 / (4.0 * np.pi)) * np.cos(theta)
    if ell == 2:
        return np.sqrt(5.0 / (16.0 * np.pi)) * (3.0 * np.cos(theta) ** 2 - 1.0)
    raise ValueError(f"僅實作 ell in {{0,1,2}} 之 m=0，收到 ell={ell}")


def hydrogen_density_xyz(x: np.ndarray, y: np.ndarray, z: np.ndarray, n: int, ell: int, m: int) -> np.ndarray:
    if m != 0:
        raise ValueError("僅示範 m=0（與 z 軸對齊）軌道")
    r = np.sqrt(x * x + y * y + z * z)
    theta = np.zeros_like(r)
    mask = r > 1e-14
    theta[mask] = np.arccos(np.clip(z[mask] / r[mask], -1.0, 1.0))
    rnl = radial_hydrogen(r, n, ell)
    yr = angular_real_y_l0(theta, ell)
    psi = rnl * yr
    return psi * psi


def spherical_angles(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    r = np.sqrt(x * x + y * y + z * z)
    theta = np.zeros_like(r)
    phi = np.zeros_like(r)
    mask = r > 1e-14
    theta[mask] = np.arccos(np.clip(z[mask] / r[mask], -1.0, 1.0))
    phi[mask] = np.arctan2(y[mask], x[mask])
    phi = np.where(phi < 0.0, phi + 2.0 * np.pi, phi)
    return theta, phi


def build_amplitude_channels(
    rho_ref: np.ndarray,
    theta: np.ndarray,
    phi: np.ndarray,
    eta_e8: float,
    closure: dict[str, float | list[float]],
) -> np.ndarray:
    amp = np.sqrt(np.clip(rho_ref, 0.0, None))
    s_dim = float(closure["s_dim8"])
    ch = np.zeros(rho_ref.shape + (8,), dtype=np.float64)
    sth, cth = np.sin(theta), np.cos(theta)
    sph, cph = np.sin(phi), np.cos(phi)

    ch[..., 0] = amp
    ch[..., 1] = amp * sth * cph * LAM ** (-1)
    ch[..., 2] = amp * sth * sph * LAM ** (-1)
    ch[..., 3] = amp * cth * LAM ** (-1)
    ch[..., 4] = amp * np.sin(2.0 * theta) * LAM ** (-2)
    ch[..., 5] = amp * np.cos(LAM * theta) * LAM ** (-2)
    ch[..., 6] = amp * np.sin(LAM * phi) * LAM ** (-2)
    ch[..., 7] = amp * eta_e8 * (1.0 + s_dim)
    return ch


def apply_Q_density(ch: np.ndarray, q: np.ndarray) -> np.ndarray:
    flat = ch.reshape(-1, 8)
    out = flat @ q.T
    return np.sum(out * out, axis=1).reshape(ch.shape[:-1])


def project_density_e8(
    rho_ref: np.ndarray,
    theta: np.ndarray,
    phi: np.ndarray,
    eta_e8: float,
    closure: dict[str, float | list[float]],
    q: np.ndarray,
) -> np.ndarray:
    ch = build_amplitude_channels(rho_ref, theta, phi, eta_e8, closure)
    return apply_Q_density(ch, q)


def density_projection_any(
    rho_ref: np.ndarray,
    theta: np.ndarray,
    phi: np.ndarray,
    eta_e8: float,
    closure: dict[str, float | list[float]],
    q: np.ndarray,
) -> np.ndarray:
    if USE_LEGACY_7D_ONLY:
        return rho_ref * phase_dwell_modulation(theta, phi)
    return project_density_e8(rho_ref, theta, phi, eta_e8, closure, q)


def metrics_full_grid(
    rho_ref: np.ndarray,
    theta: np.ndarray,
    phi: np.ndarray,
    dv: float,
    eta_e8: float,
    closure: dict[str, float | list[float]],
    q: np.ndarray,
) -> tuple[dict[str, float], float, float]:
    rho_proj = density_projection_any(rho_ref, theta, phi, eta_e8, closure, q)
    mr = float(np.sum(rho_ref) * dv)
    mp = float(np.sum(rho_proj) * dv)
    rn = rho_ref / mr
    pn = rho_proj / mp
    diff = pn - rn
    sum_sq_diff = float(np.sum(diff * diff) * dv)
    sum_sq_ref = float(np.sum(rn * rn) * dv)
    l2_abs = math.sqrt(sum_sq_diff)
    l2_rel = l2_abs / math.sqrt(sum_sq_ref) if sum_sq_ref > 0.0 else float("nan")
    return {"l2_absolute": l2_abs, "l2_relative": l2_rel}, mr, mp


def grid_search_eta(
    rho_ref: np.ndarray,
    theta: np.ndarray,
    phi: np.ndarray,
    dv: float,
    closure: dict[str, float | list[float]],
    q: np.ndarray,
) -> tuple[float, dict[str, float]]:
    steps = int(os.environ.get("HETU_E8_ETA_STEPS", "21"))
    best_eta = 0.0
    best: dict[str, float] = {}
    best_rel = float("inf")

    for k in range(steps + 1):
        eta = ETA_MAX * k / steps
        metrics, _, _ = metrics_full_grid(rho_ref, theta, phi, dv, eta, closure, q)
        rel = metrics["l2_relative"]
        if rel < best_rel:
            best_rel = rel
            best_eta = eta
            best = metrics

    return best_eta, best


def plot_comparison(
    cases: list[tuple[str, np.ndarray, np.ndarray]],
    out_png: Path,
    dx: float,
) -> None:
    import matplotlib.pyplot as plt

    ncase = len(cases)
    fig, axes = plt.subplots(ncase, 2, figsize=(9.5, 4.0 * ncase), constrained_layout=True)
    if ncase == 1:
        axes = np.array([axes])

    extent = (-R_MAX_AU, R_MAX_AU, -R_MAX_AU, R_MAX_AU)
    for i, (title, ref_slice_n, proj_slice_n) in enumerate(cases):
        vmax = max(float(np.max(ref_slice_n)), float(np.max(proj_slice_n)), 1e-30)
        im0 = axes[i, 0].imshow(ref_slice_n.T, origin="lower", cmap="magma", extent=extent, vmin=0.0, vmax=vmax)
        axes[i, 0].set_title(f"{title}: classical |psi|^2 (z=0)")
        axes[i, 0].set_xlabel("x (a.u.)")
        axes[i, 0].set_ylabel("y (a.u.)")
        fig.colorbar(im0, ax=axes[i, 0], fraction=0.046)

        im1 = axes[i, 1].imshow(proj_slice_n.T, origin="lower", cmap="magma", extent=extent, vmin=0.0, vmax=vmax)
        axes[i, 1].set_title(f"{title}: 8D E8 closure projection")
        axes[i, 1].set_xlabel("x (a.u.)")
        axes[i, 1].set_ylabel("y (a.u.)")
        fig.colorbar(im1, ax=axes[i, 1], fraction=0.046)

    fig.suptitle(
        rf"Hetu 8D E8 closure (lambda={LAM}, grid {GRID_N}^3, dx={dx:.6f} a.u.)",
        fontsize=12,
    )
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=160)
    plt.close(fig)


def patch_markdown(md_path: Path, rows: list[tuple[str, str, float, float, float]], dx: float, metrics_json_path: Path) -> None:
    text = md_path.read_text(encoding="utf-8")

    payload = {
        "protocol": "E8_closure_8D",
        "lambda": LAM,
        "grid_n": GRID_N,
        "r_max_au": R_MAX_AU,
        "dx_au": dx,
        "orbitals": [
            {"label": lab, "nlm": nlm, "eta_best": eb, "l2_relative": lr, "l2_absolute": la}
            for lab, nlm, eb, lr, la in rows
        ],
    }
    metrics_json_path.parent.mkdir(parents=True, exist_ok=True)
    metrics_json_path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")

    table_lines = [
        "| 軌道（代理） | $(n,\\ell,m)$ | $\\eta_{\\mathrm{best}}$（第 8 維增益） | $\\varepsilon_{\\mathrm{rel}}$ | $\\|\\cdot\\|_2$（絕對） |",
        "|--------------|----------------|----------------------------------------|------------------------------|------------------------|",
    ]
    for lab, nlm, eb, lr, la in rows:
        table_lines.append(f"| **{lab}** | ${nlm}$ | **{eb:.6f}** | **{lr:.6e}** | **{la:.6e}** |")

    auto_block = "\n".join(
        [
            "<!-- HETU7D_ORBITAL_TABLE_AUTO_BEGIN -->",
            f"_數值實驗：**8D E8 閉合**（`scripts/hetu7D_orbital_refinement.py`），網格 **${GRID_N}^3$** 立方體 $[-{R_MAX_AU:g},{R_MAX_AU:g}]^3$ a.u.，格距 $\\Delta x={dx:.6f}$；詳見 `plots/hetu7d_orbital_l2.json`、`plots/e8_closure_metrics.json`。",
            "",
            *table_lines,
            "<!-- HETU7D_ORBITAL_TABLE_AUTO_END -->",
        ]
    )

    begin_tag = "<!-- HETU7D_ORBITAL_TABLE_AUTO_BEGIN -->"
    end_tag = "<!-- HETU7D_ORBITAL_TABLE_AUTO_END -->"
    if begin_tag in text and end_tag in text:
        i = text.index(begin_tag)
        j = text.index(end_tag) + len(end_tag)
        text_new = text[:i] + auto_block + text[j:]
    else:
        sec4 = r"(## 4\. 形態擬合[^\n]*\n)([\s\S]*?)(\n---\s*\n\n## 5\.)"
        m = re.search(sec4, text)
        if not m:
            raise RuntimeError("無法定位 hetu7D_Ascension_Logic.md 第 4 節。")
        text_new = (
            text[: m.start()]
            + "## 4. 形態擬合：8D E8 閉合投影 vs 經典 $|Y_{\\ell m}|^2$\n\n"
            + auto_block
            + "\n\n**對比結論**：$\\varepsilon_{\\mathrm{rel}}$ 為離散 $L^2$（對經典 $|\\psi|^2$）；**$\\eta_{\\mathrm{best}}$** 為網格搜尋之第 8 維增益。若 $\\varepsilon_{\\mathrm{rel}}\\gg 10^{-12}$，量子表述仍不可替代。\n"
            + m.group(3)
            + text[m.end() :]
        )

    md_path.write_text(text_new, encoding="utf-8")


def main() -> None:
    Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

    root = Path(__file__).resolve().parents[1]
    md_path = root / "hetu7D_Ascension_Logic.md"
    plot_path = root / "plots" / "7D_orbital_projection.png"
    json_path = root / "plots" / "hetu7d_orbital_l2.json"
    e8_json_path = root / "plots" / "e8_closure_metrics.json"

    closure = e8_residual_closure_linear(RIEMANN_RESIDUAL_E8, ELECTRON_CLOUD_RESIDUAL_E8)
    lam_uni = lambda_unification_candidates(LAM)
    q = orthogonal_polar_from_cartan()

    e8_report: dict = {
        "det_cartan_E8": float(np.linalg.det(e8_cartan_matrix())),
        "residual_inputs": {"riemann": RIEMANN_RESIDUAL_E8, "electron_cloud_proxy": ELECTRON_CLOUD_RESIDUAL_E8},
        "linear_closure": closure,
        "lambda_candidates_vs_genesis": lam_uni,
        "orthogonal_Q_frobenius": float(np.linalg.norm(q, "fro")),
        "target_eps": 1e-12,
    }

    xs = np.linspace(-R_MAX_AU, R_MAX_AU, GRID_N)
    ys = np.linspace(-R_MAX_AU, R_MAX_AU, GRID_N)
    zs = np.linspace(-R_MAX_AU, R_MAX_AU, GRID_N)
    dx = float(xs[1] - xs[0])
    dv = dx**3

    z0_idx = GRID_N // 2

    xg, yg, zg = np.meshgrid(xs, ys, zs, indexing="ij")

    orbitals = [
        ("1s", 1, 0, 0),
        ("2p_z", 2, 1, 0),
        ("3d_{z^2}", 3, 2, 0),
    ]

    plot_cases: list[tuple[str, np.ndarray, np.ndarray]] = []
    rows: list[tuple[str, str, float, float, float]] = []

    print("=== E8 Closure diagnostics ===")
    print(f"det(A_E8) = {e8_report['det_cartan_E8']:.12f} (expect 1)")
    print(f"linear <r,omega8> = {closure['dot_r_with_omega8']:.3e}")
    print(f"||r||^2 (R^8) = {closure['quadratic_norm_r_R8']:.8f}")
    print(f"s (dim8 compensation scalar) = {closure['s_dim8']:.8f}")

    min_eps = float("inf")

    for label, n, ell, m in orbitals:
        rho_ref = hydrogen_density_xyz(xg, yg, zg, n, ell, m)
        theta, phi = spherical_angles(xg, yg, zg)

        eta_best, metrics_best = grid_search_eta(rho_ref, theta, phi, dv, closure, q)
        min_eps = min(min_eps, metrics_best["l2_relative"])

        rho_proj = density_projection_any(rho_ref, theta, phi, eta_best, closure, q)
        mr = float(np.sum(rho_ref) * dv)
        mp = float(np.sum(rho_proj) * dv)

        sr_n = rho_ref[z0_idx, :, :] / mr
        sp_n = rho_proj[z0_idx, :, :] / mp

        plot_cases.append((label, sr_n, sp_n))
        rows.append((label, f"{n},{ell},{m}", eta_best, metrics_best["l2_relative"], metrics_best["l2_absolute"]))

        hit = metrics_best["l2_relative"] <= 1e-12
        print(
            f"{label}  eta_best={eta_best:.6f}  eps_rel={metrics_best['l2_relative']:.12e}  "
            f"L2_abs={metrics_best['l2_absolute']:.12e}  hits_1e-12={hit}"
        )

    e8_report["min_eps_rel_across_orbitals"] = float(min_eps)
    e8_report["geometry_can_replace_QM_claim"] = bool(min_eps <= 1e-12)
    e8_json_path.write_text(json.dumps(e8_report, indent=2, ensure_ascii=False), encoding="utf-8")

    plot_comparison(plot_cases, plot_path, dx)
    patch_markdown(md_path, rows, dx, json_path)
    print(f"Wrote {plot_path}")
    print(f"Wrote {json_path}")
    print(f"Wrote {e8_json_path}")
    print(f"Updated {md_path}")


if __name__ == "__main__":
    main()
