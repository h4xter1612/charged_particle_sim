from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def main():
    # --------------------------------------------------
    # Paths
    # --------------------------------------------------

    script_dir = Path(__file__).resolve().parent
    project_dir = script_dir.parent

    csv_path = project_dir / "results" / "polarization.csv"
    figures_dir = project_dir / "results" / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    # --------------------------------------------------
    # Load numerical trajectory
    # --------------------------------------------------

    data = np.genfromtxt(
        csv_path,
        delimiter=",",
        names=True,
    )

    t = data["t"]
    x = data["x"]
    y = data["y"]
    z = data["z"]
    vx = data["vx"]
    vy = data["vy"]
    vz = data["vz"]
    kinetic_energy = data["kinetic_energy"]

    # --------------------------------------------------
    # Benchmark parameters
    #
    # E(t) = E0 cos(omega_E t) x-hat
    # B    = B0 z-hat
    #
    # Initial conditions:
    # r(0) = (0, 0, 0)
    # v(0) = (0, -0.2, 0.2)
    # --------------------------------------------------

    q = 1.0
    mass = 1.0

    E0 = 0.2
    B0 = 1.0
    omega_E = 0.05

    x0 = 0.0
    y0 = 0.0
    z0 = 0.0

    vx0 = 0.0
    vy0 = -0.2
    vz0 = 0.2

    q_over_m = q / mass
    omega_c = q_over_m * B0

    if np.isclose(abs(omega_c), omega_E):
        raise ValueError(
            "The analytical solution used here assumes omega_E != |Omega_c|."
        )

    # --------------------------------------------------
    # Exact full-orbit solution
    #
    # The perpendicular equations are:
    #
    # dvx/dt = (q/m) E0 cos(omega_E t) + Omega_c vy
    # dvy/dt = -Omega_c vx
    #
    # Eliminating vx gives a forced harmonic oscillator:
    #
    # d2(vy)/dt2 + Omega_c^2 vy
    #     = -Omega_c (q/m) E0 cos(omega_E t)
    #
    # The exact solution contains:
    #   1. a small cyclotron-frequency homogeneous term, and
    #   2. the slow response at omega_E.
    # --------------------------------------------------

    slow_vy_amplitude = (
        -omega_c
        * q_over_m
        * E0
        / (omega_c**2 - omega_E**2)
    )

    gyro_vy_amplitude = vy0 - slow_vy_amplitude

    cos_c = np.cos(omega_c * t)
    sin_c = np.sin(omega_c * t)
    cos_e = np.cos(omega_E * t)
    sin_e = np.sin(omega_E * t)

    vy_exact = (
        gyro_vy_amplitude * cos_c
        - vx0 * sin_c
        + slow_vy_amplitude * cos_e
    )

    vx_exact = (
        gyro_vy_amplitude * sin_c
        + vx0 * cos_c
        + slow_vy_amplitude
        * (omega_E / omega_c)
        * sin_e
    )

    x_exact = (
        x0
        + gyro_vy_amplitude
        * (1.0 - cos_c)
        / omega_c
        + vx0 * sin_c / omega_c
        + slow_vy_amplitude
        * (1.0 - cos_e)
        / omega_c
    )

    y_exact = (
        y0
        + gyro_vy_amplitude
        * sin_c
        / omega_c
        + vx0
        * (cos_c - 1.0)
        / omega_c
        + slow_vy_amplitude
        * sin_e
        / omega_E
    )

    z_exact = z0 + vz0 * t
    vz_exact = np.full_like(t, vz0)

    kinetic_energy_exact = 0.5 * mass * (
        vx_exact**2
        + vy_exact**2
        + vz_exact**2
    )

    # --------------------------------------------------
    # First-order drift theory
    #
    # For omega_E / Omega_c << 1:
    #
    # v_E = E x B / B^2
    #
    # v_E,y = -(E0/B0) cos(omega_E t)
    #
    # and the polarization drift is:
    #
    # v_p = (m / q B^2) dE_perp/dt
    #
    # v_p,x = -(m E0 omega_E / q B0^2)
    #         sin(omega_E t)
    # --------------------------------------------------

    v_exb_y_theory = (
        -E0 / B0
        * cos_e
    )

    v_pol_x_theory = (
        -mass
        * E0
        * omega_E
        / (q * B0**2)
        * sin_e
    )

    temporal_ordering = (
        omega_E / abs(omega_c)
    )

    exb_amplitude_first_order = E0 / B0
    polarization_amplitude_first_order = (
        mass
        * E0
        * omega_E
        / (abs(q) * B0**2)
    )

    # The exact slow response differs from first-order
    # drift theory by the factor 1 / [1 - (omega_E/Omega_c)^2].

    slow_response_factor = (
        1.0
        / (1.0 - temporal_ordering**2)
    )

    exb_amplitude_exact_slow = abs(
        slow_vy_amplitude
    )

    polarization_amplitude_exact_slow = abs(
        slow_vy_amplitude
        * omega_E
        / omega_c
    )

    # --------------------------------------------------
    # Numerical error against the exact full orbit
    # --------------------------------------------------

    max_x_error = np.max(np.abs(x - x_exact))
    max_y_error = np.max(np.abs(y - y_exact))
    max_z_error = np.max(np.abs(z - z_exact))

    max_vx_error = np.max(np.abs(vx - vx_exact))
    max_vy_error = np.max(np.abs(vy - vy_exact))
    max_vz_error = np.max(np.abs(vz - vz_exact))

    max_energy_error = np.max(
        np.abs(
            kinetic_energy
            - kinetic_energy_exact
        )
    )

    # --------------------------------------------------
    # Work-energy consistency
    #
    # The magnetic field performs no work. Therefore:
    #
    # dK/dt = q E dot v = q Ex vx
    # --------------------------------------------------

    Ex = E0 * cos_e
    electric_power = q * Ex * vx

    dt_values = np.diff(t)
    work_increments = (
        0.5
        * (electric_power[:-1] + electric_power[1:])
        * dt_values
    )

    cumulative_work = np.concatenate((
        [0.0],
        np.cumsum(work_increments),
    ))

    energy_from_work = (
        kinetic_energy[0]
        + cumulative_work
    )

    max_work_energy_error = np.max(
        np.abs(
            kinetic_energy
            - energy_from_work
        )
    )

    # --------------------------------------------------
    # Console diagnostics
    # --------------------------------------------------

    print("\n=== Polarization drift analysis ===")
    print(f"CSV file: {csv_path}")
    print(f"Number of saved states: {len(t)}")

    print("\nTemporal ordering:")
    print(
        f"  cyclotron frequency Omega_c = "
        f"{omega_c:.8f}"
    )
    print(
        f"  electric-field frequency omega_E = "
        f"{omega_E:.8f}"
    )
    print(
        f"  omega_E / Omega_c = "
        f"{temporal_ordering:.6e}"
    )

    print("\nFirst-order drift amplitudes:")
    print(
        f"  max |v_E| = "
        f"{exb_amplitude_first_order:.8f}"
    )
    print(
        f"  max |v_p| = "
        f"{polarization_amplitude_first_order:.8f}"
    )

    print("\nExact slow-response amplitudes:")
    print(
        f"  max |v_y,slow| = "
        f"{exb_amplitude_exact_slow:.8f}"
    )
    print(
        f"  max |v_x,slow| = "
        f"{polarization_amplitude_exact_slow:.8f}"
    )
    print(
        f"  finite-frequency correction factor = "
        f"{slow_response_factor:.8f}"
    )

    print("\nExact full-orbit comparison:")
    print(
        f"  max |x_num - x_exact| = "
        f"{max_x_error:.6e}"
    )
    print(
        f"  max |y_num - y_exact| = "
        f"{max_y_error:.6e}"
    )
    print(
        f"  max |z_num - z_exact| = "
        f"{max_z_error:.6e}"
    )
    print(
        f"  max |vx_num - vx_exact| = "
        f"{max_vx_error:.6e}"
    )
    print(
        f"  max |vy_num - vy_exact| = "
        f"{max_vy_error:.6e}"
    )
    print(
        f"  max |vz_num - vz_exact| = "
        f"{max_vz_error:.6e}"
    )
    print(
        f"  max |K_num - K_exact| = "
        f"{max_energy_error:.6e}"
    )

    print("\nWork-energy consistency:")
    print(
        f"  max |K_num - (K0 + W_E)| = "
        f"{max_work_energy_error:.6e}"
    )

    # --------------------------------------------------
    # Figure 1: numerical vs exact perpendicular orbit
    # --------------------------------------------------

    plt.figure(figsize=(8, 7))
    plt.plot(x, y, label="Numerical orbit")
    plt.plot(
        x_exact,
        y_exact,
        "--",
        label="Exact full-orbit solution",
    )

    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("Polarization benchmark: perpendicular motion")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    plt.savefig(
        figures_dir / "polarization_xy_motion.png",
        dpi=200,
    )

    # --------------------------------------------------
    # Figure 2: polarization velocity
    #
    # The numerical full orbit contains a small residual
    # cyclotron oscillation because the initial velocity
    # matches first-order drift theory rather than the
    # exact finite-frequency slow solution.
    # --------------------------------------------------

    slow_vx_exact = (
        slow_vy_amplitude
        * (omega_E / omega_c)
        * sin_e
    )

    plt.figure(figsize=(8, 5))
    plt.plot(
        t,
        vx,
        label="Numerical full-orbit vx",
    )
    plt.plot(
        t,
        slow_vx_exact,
        "--",
        label="Exact slow response",
    )
    plt.plot(
        t,
        v_pol_x_theory,
        ":",
        label="First-order polarization drift",
    )

    plt.xlabel("t")
    plt.ylabel("Velocity in x")
    plt.title("Polarization drift: x velocity")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    plt.savefig(
        figures_dir / "polarization_velocity_x.png",
        dpi=200,
    )

    # --------------------------------------------------
    # Figure 3: time-dependent E x B response
    # --------------------------------------------------

    slow_vy_exact = slow_vy_amplitude * cos_e

    plt.figure(figsize=(8, 5))
    plt.plot(
        t,
        vy,
        label="Numerical full-orbit vy",
    )
    plt.plot(
        t,
        slow_vy_exact,
        "--",
        label="Exact slow response",
    )
    plt.plot(
        t,
        v_exb_y_theory,
        ":",
        label="First-order E x B drift",
    )

    plt.xlabel("t")
    plt.ylabel("Velocity in y")
    plt.title("Time-dependent E x B response")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    plt.savefig(
        figures_dir / "polarization_exb_velocity_y.png",
        dpi=200,
    )

    # --------------------------------------------------
    # Figure 4: position components
    # --------------------------------------------------

    plt.figure(figsize=(8, 6))
    plt.plot(t, x, label="x numerical")
    plt.plot(t, x_exact, "--", label="x exact")
    plt.plot(t, y, label="y numerical")
    plt.plot(t, y_exact, "--", label="y exact")

    plt.xlabel("t")
    plt.ylabel("Position")
    plt.title("Polarization benchmark: position components")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    plt.savefig(
        figures_dir / "polarization_position.png",
        dpi=200,
    )

    # --------------------------------------------------
    # Figure 5: work-energy theorem
    # --------------------------------------------------

    plt.figure(figsize=(8, 5))
    plt.plot(
        t,
        kinetic_energy,
        label="Numerical kinetic energy",
    )
    plt.plot(
        t,
        kinetic_energy_exact,
        "--",
        label="Exact kinetic energy",
    )
    plt.plot(
        t,
        energy_from_work,
        ":",
        label="K0 + electric-field work",
    )

    plt.xlabel("t")
    plt.ylabel("Energy")
    plt.title("Polarization benchmark: energy consistency")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    plt.savefig(
        figures_dir / "polarization_energy.png",
        dpi=200,
    )

    # --------------------------------------------------
    # Figure 6: three-dimensional trajectory
    # --------------------------------------------------

    fig = plt.figure(figsize=(8, 7))
    ax = fig.add_subplot(111, projection="3d")

    ax.plot(
        x,
        y,
        z,
        label="Numerical trajectory",
    )
    ax.plot(
        x_exact,
        y_exact,
        z_exact,
        "--",
        label="Exact full-orbit solution",
    )

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.set_title("Polarization benchmark: 3D trajectory")
    ax.legend()
    plt.tight_layout()

    plt.savefig(
        figures_dir / "polarization_3d_trajectory.png",
        dpi=200,
    )

    plt.show()


if __name__ == "__main__":
    main()
