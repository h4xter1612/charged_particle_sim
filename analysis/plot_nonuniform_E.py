from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def main():
    # --------------------------------------------------
    # Paths
    # --------------------------------------------------

    script_dir = Path(__file__).resolve().parent
    project_dir = script_dir.parent

    csv_path = project_dir / "results" / "nonuniform_E.csv"
    figures_dir = project_dir / "results" / "figures"

    figures_dir.mkdir(parents=True, exist_ok=True)

    # --------------------------------------------------
    # Load numerical trajectory
    # --------------------------------------------------

    data = np.genfromtxt(
        csv_path,
        delimiter=",",
        names=True
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
    # Magnetic field:
    #
    # B = B0 z-hat
    #
    # Electric field:
    #
    # E(x) = E0 cos(k x) x-hat
    #
    # Particle:
    #
    # q = 1
    # m = 1
    #
    # Initial conditions:
    #
    # r(0) = (0, 0, 0)
    # v(0) = (1, 0, 0.2)
    #
    # For the unperturbed magnetic orbit:
    #
    # rho_L = 1
    #
    # and the initial guiding-center x coordinate is
    # X_gc = 0.
    # --------------------------------------------------

    q = 1.0
    mass = 1.0

    B0 = 1.0
    E0 = 0.02
    wavenumber = 0.4

    # Signed cyclotron frequency:
    #
    # Omega_c = q B / m

    omega_c = q * B0 / mass

    gyroperiod = (
        2.0 * np.pi
        / abs(omega_c)
    )

    # --------------------------------------------------
    # Electric field sampled by the actual particle
    # --------------------------------------------------

    Ex = (
        E0
        * np.cos(wavenumber * x)
    )

    # --------------------------------------------------
    # Guiding-center coordinates
    #
    # For a uniform magnetic field along +z:
    #
    # R = r + (v x b) / Omega_c
    #
    # Therefore:
    #
    # X_gc = x + vy / Omega_c
    # Y_gc = y - vx / Omega_c
    #
    # For E = Ex(x) x-hat, X_gc is exactly constant:
    #
    # dX_gc/dt = 0
    #
    # while:
    #
    # dY_gc/dt = -Ex(x) / B0
    #
    # The rapid oscillation in Ex(x) is produced by the
    # finite Larmor radius. Its gyro-average determines
    # the slow drift.
    # --------------------------------------------------

    gc_x = (
        x + vy / omega_c
    )

    gc_y = (
        y - vx / omega_c
    )

    # --------------------------------------------------
    # Initial Larmor radius and FLR ordering
    #
    # rho_L = v_perp / |Omega_c|
    #
    # The expansion parameter for this problem is:
    #
    # k rho_L << 1
    # --------------------------------------------------

    v_perp_initial = np.sqrt(
        vx[0]**2
        + vy[0]**2
    )

    larmor_radius = (
        v_perp_initial
        / abs(omega_c)
    )

    krho = (
        wavenumber
        * larmor_radius
    )

    # --------------------------------------------------
    # Zeroth-order local E x B prediction
    #
    # If the particle had zero Larmor radius, it would
    # simply experience the electric field evaluated at
    # the guiding center:
    #
    # v_E = E x B / B^2
    #
    # For E along +x and B along +z:
    #
    # v_E,y = -Ex / B0
    # --------------------------------------------------

    gc_x_initial = gc_x[0]

    local_exb_drift = (
        -E0
        / B0
        * np.cos(
            wavenumber
            * gc_x_initial
        )
    )

    # --------------------------------------------------
    # Second-order finite-Larmor-radius correction
    #
    # For k rho_L << 1:
    #
    # <v_E> =
    #
    #   (1 - k^2 rho_L^2 / 4)
    #
    #   * [E(R_gc) x B] / B^2
    #
    # This is the standard small-Larmor-radius
    # approximation for the nonuniform-E benchmark.
    # --------------------------------------------------

    flr_factor_second_order = (
        1.0
        - krho**2 / 4.0
    )

    flr_drift_second_order = (
        local_exb_drift
        * flr_factor_second_order
    )

    # --------------------------------------------------
    # Direct gyro-average of the sinusoidal field
    #
    # For an unperturbed circular orbit:
    #
    # x(theta) =
    #     X_gc + rho_L sin(theta)
    #
    # Therefore we can directly average:
    #
    # <E> =
    #     (1 / 2pi)
    #     integral E[x(theta)] dtheta
    #
    # This is equivalent to the Bessel-function result:
    #
    # <E> =
    #     E0 cos(k X_gc) J0(k rho_L)
    #
    # Here we evaluate the gyro-average numerically so
    # that no SciPy dependency is required.
    # --------------------------------------------------

    gyroangle = np.linspace(
        0.0,
        2.0 * np.pi,
        100000,
        endpoint=False
    )

    orbit_x_theory = (
        gc_x_initial
        + larmor_radius
        * np.sin(gyroangle)
    )

    gyroaveraged_E = np.mean(
        E0
        * np.cos(
            wavenumber
            * orbit_x_theory
        )
    )

    flr_drift_gyroaverage = (
        -gyroaveraged_E
        / B0
    )

    flr_factor_gyroaverage = (
        gyroaveraged_E
        / (
            E0
            * np.cos(
                wavenumber
                * gc_x_initial
            )
        )
    )

    # --------------------------------------------------
    # Numerical drift measurement
    #
    # The guiding center contains a small periodic FLR
    # oscillation superimposed on a nearly constant
    # average drift.
    #
    # A linear fit to Y_gc(t) measures the slow
    # numerical drift velocity.
    # --------------------------------------------------

    drift_fit = np.polyfit(
        t,
        gc_y,
        1
    )

    drift_velocity_numerical = (
        drift_fit[0]
    )

    gc_y_intercept = (
        drift_fit[1]
    )

    gc_y_fit = (
        drift_velocity_numerical * t
        + gc_y_intercept
    )

    # Theoretical guiding-center trajectories.

    gc_y_local_theory = (
        gc_y[0]
        + local_exb_drift * t
    )

    gc_y_flr_second_order = (
        gc_y[0]
        + flr_drift_second_order * t
    )

    gc_y_flr_gyroaverage = (
        gc_y[0]
        + flr_drift_gyroaverage * t
    )

    # --------------------------------------------------
    # Relative errors of the different drift models
    # --------------------------------------------------

    local_exb_error = abs(
        (
            drift_velocity_numerical
            - local_exb_drift
        )
        / local_exb_drift
    )

    second_order_flr_error = abs(
        (
            drift_velocity_numerical
            - flr_drift_second_order
        )
        / flr_drift_second_order
    )

    gyroaverage_flr_error = abs(
        (
            drift_velocity_numerical
            - flr_drift_gyroaverage
        )
        / flr_drift_gyroaverage
    )

    # --------------------------------------------------
    # Invariance of X_gc
    #
    # With Ey = 0 and uniform B:
    #
    # dX_gc/dt = 0
    #
    # This provides an independent numerical check.
    # --------------------------------------------------

    max_gc_x_deviation = np.max(
        np.abs(
            gc_x - gc_x[0]
        )
    )

    # --------------------------------------------------
    # Energy conservation
    #
    # This electric field is electrostatic:
    #
    # E_x(x) = E0 cos(kx)
    #
    # with electrostatic potential:
    #
    # phi(x) =
    #     -(E0/k) sin(kx)
    #
    # because:
    #
    # E = -grad(phi)
    #
    # Therefore the total particle energy:
    #
    # H = K + q phi
    #
    # should be conserved.
    # --------------------------------------------------

    electric_potential = (
        -E0
        / wavenumber
        * np.sin(
            wavenumber * x
        )
    )

    total_energy = (
        kinetic_energy
        + q * electric_potential
    )

    initial_total_energy = (
        total_energy[0]
    )

    relative_total_energy_error = (
        total_energy
        - initial_total_energy
    ) / abs(initial_total_energy)

    max_relative_total_energy_error = np.max(
        np.abs(
            relative_total_energy_error
        )
    )

    # --------------------------------------------------
    # Console diagnostics
    # --------------------------------------------------

    print("\n=== Nonuniform-E drift analysis ===")
    print(f"CSV file: {csv_path}")
    print(f"Number of saved states: {len(t)}")

    print("\nFinite-Larmor-radius ordering:")
    print(
        f"  gyroperiod T_c = "
        f"{gyroperiod:.8f}"
    )
    print(
        f"  Larmor radius rho_L = "
        f"{larmor_radius:.8f}"
    )
    print(
        f"  k rho_L = "
        f"{krho:.8f}"
    )

    print("\nFLR correction factors:")
    print(
        f"  second-order factor = "
        f"{flr_factor_second_order:.8f}"
    )
    print(
        f"  direct gyro-average factor = "
        f"{flr_factor_gyroaverage:.8f}"
    )

    print("\nDrift predictions:")
    print(
        f"  local E x B, no FLR = "
        f"{local_exb_drift:.8f}"
    )
    print(
        f"  second-order FLR = "
        f"{flr_drift_second_order:.8f}"
    )
    print(
        f"  direct gyro-average = "
        f"{flr_drift_gyroaverage:.8f}"
    )

    print("\nNumerical measurement:")
    print(
        f"  numerical drift = "
        f"{drift_velocity_numerical:.8f}"
    )

    print("\nRelative errors:")
    print(
        f"  local E x B error = "
        f"{local_exb_error:.6e}"
    )
    print(
        f"  second-order FLR error = "
        f"{second_order_flr_error:.6e}"
    )
    print(
        f"  gyro-average FLR error = "
        f"{gyroaverage_flr_error:.6e}"
    )

    print("\nGuiding-center invariant:")
    print(
        f"  max |X_gc - X_gc(0)| = "
        f"{max_gc_x_deviation:.6e}"
    )

    print("\nTotal-energy conservation:")
    print(
        f"  max relative total-energy error = "
        f"{max_relative_total_energy_error:.6e}"
    )

    # --------------------------------------------------
    # Figure 1:
    # Particle orbit and guiding-center trajectory
    # --------------------------------------------------

    plt.figure(figsize=(8, 7))

    plt.plot(
        x,
        y,
        label="Particle orbit"
    )

    plt.plot(
        gc_x,
        gc_y,
        "--",
        label="Guiding center"
    )

    plt.xlabel("x")
    plt.ylabel("y")

    plt.title(
        "Nonuniform E: particle orbit and guiding-center drift"
    )

    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    plt.savefig(
        figures_dir / "nonuniform_E_xy_orbit.png",
        dpi=200
    )

    # --------------------------------------------------
    # Figure 2:
    # Guiding-center drift comparison
    #
    # This is the main result:
    #
    # numerical full orbit
    # vs.
    # local E x B
    # vs.
    # second-order FLR theory
    # vs.
    # direct gyro-average
    # --------------------------------------------------

    plt.figure(figsize=(8, 6))

    plt.plot(
        t,
        gc_y,
        label="Numerical guiding center"
    )

    plt.plot(
        t,
        gc_y_fit,
        "--",
        label=(
            "Numerical fit: "
            f"v = {drift_velocity_numerical:.6f}"
        )
    )

    plt.plot(
        t,
        gc_y_local_theory,
        ":",
        label=(
            "Local E x B: "
            f"v = {local_exb_drift:.6f}"
        )
    )

    plt.plot(
        t,
        gc_y_flr_second_order,
        "-.",
        label=(
            "Second-order FLR: "
            f"v = {flr_drift_second_order:.6f}"
        )
    )

    plt.xlabel("t")
    plt.ylabel("Guiding-center y")

    plt.title(
        "Nonuniform E: finite-Larmor-radius correction"
    )

    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    plt.savefig(
        figures_dir / "nonuniform_E_guiding_center.png",
        dpi=200
    )

    # --------------------------------------------------
    # Figure 3:
    # Electric field sampled during the gyro-orbit
    #
    # The particle does not see only E(X_gc).
    # It samples weaker and stronger field values around
    # its finite-radius orbit.
    # --------------------------------------------------

    plt.figure(figsize=(8, 5))

    plt.plot(
        t,
        Ex
    )

    plt.axhline(
        gyroaveraged_E,
        linestyle="--",
        label=(
            "Theoretical gyro-averaged E = "
            f"{gyroaveraged_E:.6f}"
        )
    )

    plt.xlabel("t")
    plt.ylabel("E_x")

    plt.title(
        "Nonuniform E: electric field sampled by the particle"
    )

    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    plt.savefig(
        figures_dir / "nonuniform_E_sampled_field.png",
        dpi=200
    )

    # --------------------------------------------------
    # Figure 4:
    # Total energy conservation
    #
    # K alone can vary because the electrostatic field
    # performs work. K + q phi should remain constant.
    # --------------------------------------------------

    plt.figure(figsize=(8, 5))

    step = 20  # Plot only 1 out of every 20 points

    plt.plot(
        t[::step],
        total_energy[::step],
        linewidth=1.2
    )

    plt.xlabel("t")
    plt.ylabel("K + q phi")

    plt.title(
        "Nonuniform E: total-energy conservation"
    )

    plt.grid(True)
    plt.tight_layout()

    plt.savefig(
        figures_dir / "nonuniform_E_total_energy.png",
        dpi=200
    )

    # --------------------------------------------------
    # Figure 5:
    # Three-dimensional trajectory
    # --------------------------------------------------

    fig = plt.figure(figsize=(8, 7))

    ax = fig.add_subplot(
        111,
        projection="3d"
    )

    ax.plot(
        x,
        y,
        z,
        label="Particle trajectory"
    )

    ax.plot(
        gc_x,
        gc_y,
        z,
        "--",
        label="Guiding center"
    )

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")

    ax.set_title(
        "Nonuniform E: 3D trajectory"
    )

    ax.legend()
    plt.tight_layout()

    plt.savefig(
        figures_dir / "nonuniform_E_3d_trajectory.png",
        dpi=200
    )

    plt.show()


if __name__ == "__main__":
    main()
