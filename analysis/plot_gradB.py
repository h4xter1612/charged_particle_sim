from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def main():
    # --------------------------------------------------
    # Paths
    # --------------------------------------------------

    script_dir = Path(__file__).resolve().parent
    project_dir = script_dir.parent

    csv_path = project_dir / "results" / "gradB.csv"
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
    # B(x) = B0 * (1 + gradient * x) z-hat
    #
    # with:
    #
    # B0       = 1
    # gradient = 0.05
    #
    # Particle:
    #
    # q = 1
    # m = 1
    #
    # Initial velocity:
    #
    # v(0) = (1, 0, 0.2)
    # --------------------------------------------------

    q = 1.0
    mass = 1.0

    B0 = 1.0
    gradient = 0.05

    # Spatial derivative of the magnetic-field magnitude:
    #
    # dB/dx = B0 * gradient

    dB_dx = B0 * gradient

    # --------------------------------------------------
    # Local magnetic field along the full particle orbit
    # --------------------------------------------------

    Bz = B0 * (1.0 + gradient * x)

    if np.any(Bz <= 0.0):
        raise ValueError(
            "The benchmark requires Bz > 0 along the entire trajectory."
        )

    # Local signed cyclotron frequency:
    #
    # Omega(x) = q B(x) / m

    omega_local = q * Bz / mass

    # --------------------------------------------------
    # Approximate guiding-center position
    #
    # For a magnetic field along +z:
    #
    # R = r + (v x b) / Omega
    #
    # Therefore:
    #
    # R_x = x + vy / Omega
    # R_y = y - vx / Omega
    #
    # In a weakly nonuniform field this is a
    # first-order guiding-center approximation.
    # --------------------------------------------------

    guiding_center_x = x + vy / omega_local
    guiding_center_y = y - vx / omega_local

    # --------------------------------------------------
    # Magnetic moment
    #
    # mu = m v_perp^2 / (2 B)
    #
    # In a slowly varying magnetic field, mu should be
    # approximately conserved:
    #
    # rho_L / L_B << 1
    #
    # This gives us an additional physical validation
    # of the simulation and of the guiding-center regime.
    # --------------------------------------------------

    v_perp_squared = vx**2 + vy**2

    magnetic_moment = (
        mass * v_perp_squared
        / (2.0 * Bz)
    )

    mu_initial = magnetic_moment[0]

    relative_mu_change = (
        magnetic_moment - mu_initial
    ) / mu_initial

    max_relative_mu_change = np.max(
        np.abs(relative_mu_change)
    )

    # --------------------------------------------------
    # Theoretical grad-B drift
    #
    # General expression:
    #
    # v_gradB =
    #     mu / (q B^2) * (B x grad B)
    #
    # For:
    #
    # B = B z-hat
    # grad B = (dB/dx) x-hat
    #
    # we have:
    #
    # B x grad B
    #     = B * (dB/dx) y-hat
    #
    # therefore:
    #
    # v_gradB,y =
    #     mu / (q B) * dB/dx
    #
    # We evaluate the leading-order prediction using
    # the initial guiding-center quantities.
    # --------------------------------------------------

    B_gc_initial = (
        B0
        * (
            1.0
            + gradient * guiding_center_x[0]
        )
    )

    v_perp_initial_squared = (
        vx[0]**2
        + vy[0]**2
    )

    mu_theory = (
        mass
        * v_perp_initial_squared
        / (2.0 * B_gc_initial)
    )

    drift_velocity_theory = (
        mu_theory
        * dB_dx
        / (q * B_gc_initial)
    )

    # --------------------------------------------------
    # Numerical drift measurement
    #
    # Fit the guiding-center motion:
    #
    # R_y(t) = v_drift * t + constant
    #
    # The slope gives the numerical drift velocity.
    # --------------------------------------------------

    drift_fit = np.polyfit(
        t,
        guiding_center_y,
        1
    )

    drift_velocity_numerical = drift_fit[0]
    guiding_center_intercept = drift_fit[1]

    guiding_center_y_fit = (
        drift_velocity_numerical * t
        + guiding_center_intercept
    )

    # Theoretical straight-line guiding-center motion
    # starting from the same initial guiding center.

    guiding_center_y_theory = (
        guiding_center_y[0]
        + drift_velocity_theory * t
    )

    relative_drift_error = abs(
        (
            drift_velocity_numerical
            - drift_velocity_theory
        )
        / drift_velocity_theory
    )

    # --------------------------------------------------
    # Energy conservation
    #
    # There is no electric field.
    #
    # Therefore:
    #
    # dK/dt = q E dot v = 0
    #
    # and the magnetic field performs no work.
    # --------------------------------------------------

    initial_energy = kinetic_energy[0]

    relative_energy_error = (
        kinetic_energy - initial_energy
    ) / initial_energy

    max_relative_energy_error = np.max(
        np.abs(relative_energy_error)
    )

    # --------------------------------------------------
    # Characteristic scale separation
    #
    # At the initial guiding center:
    #
    # rho_L = v_perp / Omega
    #
    # L_B = B / |grad B|
    #
    # epsilon = rho_L / L_B
    #
    # Guiding-center theory assumes epsilon << 1.
    # --------------------------------------------------

    v_perp_initial = np.sqrt(
        v_perp_initial_squared
    )

    omega_initial = (
        abs(q)
        * B_gc_initial
        / mass
    )

    larmor_radius_initial = (
        v_perp_initial
        / omega_initial
    )

    magnetic_scale_length = (
        B_gc_initial
        / abs(dB_dx)
    )

    epsilon = (
        larmor_radius_initial
        / magnetic_scale_length
    )

    # --------------------------------------------------
    # Console diagnostics
    # --------------------------------------------------

    print("\n=== Grad-B drift analysis ===")
    print(f"CSV file: {csv_path}")
    print(f"Number of saved states: {len(t)}")

    print("\nGuiding-center ordering:")
    print(
        f"  initial Larmor radius rho_L = "
        f"{larmor_radius_initial:.8f}"
    )
    print(
        f"  magnetic scale length L_B  = "
        f"{magnetic_scale_length:.8f}"
    )
    print(
        f"  epsilon = rho_L / L_B      = "
        f"{epsilon:.6e}"
    )

    print("\nTheoretical prediction:")
    print(
        f"  mu_0 = {mu_theory:.8f}"
    )
    print(
        f"  v_gradB,y = "
        f"{drift_velocity_theory:.8f}"
    )

    print("\nNumerical guiding-center measurement:")
    print(
        f"  v_drift,y = "
        f"{drift_velocity_numerical:.8f}"
    )
    print(
        f"  relative drift error = "
        f"{relative_drift_error:.6e}"
    )

    print("\nAdiabatic invariant:")
    print(
        f"  max relative variation of mu = "
        f"{max_relative_mu_change:.6e}"
    )

    print("\nEnergy conservation:")
    print(
        f"  max relative energy error = "
        f"{max_relative_energy_error:.6e}"
    )

    # --------------------------------------------------
    # Figure 1:
    # Particle orbit and approximate guiding center
    # --------------------------------------------------

    plt.figure(figsize=(8, 7))

    plt.plot(
        x,
        y,
        label="Particle orbit"
    )

    plt.plot(
        guiding_center_x,
        guiding_center_y,
        "--",
        label="Approximate guiding center"
    )

    plt.xlabel("x")
    plt.ylabel("y")

    plt.title(
        "Grad-B drift: particle orbit and guiding center"
    )

    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    plt.savefig(
        figures_dir / "gradB_xy_orbit.png",
        dpi=200
    )

    # --------------------------------------------------
    # Figure 2:
    # Guiding-center drift
    #
    # The numerical linear fit is compared with the
    # first-order guiding-center prediction.
    # --------------------------------------------------

    plt.figure(figsize=(8, 5))

    plt.plot(
        t,
        guiding_center_y,
        label="Numerical guiding center"
    )

    plt.plot(
        t,
        guiding_center_y_fit,
        "--",
        label=(
            "Numerical linear fit: "
            f"v = {drift_velocity_numerical:.6f}"
        )
    )

    plt.plot(
        t,
        guiding_center_y_theory,
        ":",
        label=(
            "First-order theory: "
            f"v = {drift_velocity_theory:.6f}"
        )
    )

    plt.xlabel("t")
    plt.ylabel("Guiding-center y")

    plt.title(
        "Grad-B drift: guiding-center motion"
    )

    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    plt.savefig(
        figures_dir / "gradB_guiding_center.png",
        dpi=200
    )

    # --------------------------------------------------
    # Figure 3:
    # Magnetic moment
    #
    # mu should be approximately conserved rather than
    # exactly constant because guiding-center theory is
    # an asymptotic approximation.
    # --------------------------------------------------

    plt.figure(figsize=(8, 5))

    plt.plot(
        t,
        magnetic_moment
    )

    plt.xlabel("t")
    plt.ylabel("Magnetic moment mu")

    plt.title(
        "Grad-B drift: magnetic moment"
    )

    plt.grid(True)
    plt.tight_layout()

    plt.savefig(
        figures_dir / "gradB_magnetic_moment.png",
        dpi=200
    )

    # --------------------------------------------------
    # Figure 4:
    # Kinetic energy
    #
    # A static magnetic field performs no work, so
    # kinetic energy should remain numerically constant.
    # --------------------------------------------------

    plt.figure(figsize=(8, 5))

    plt.plot(
        t,
        kinetic_energy
    )

    plt.xlabel("t")
    plt.ylabel("Kinetic energy")

    plt.title(
        "Grad-B drift: kinetic energy conservation"
    )

    plt.grid(True)
    plt.tight_layout()

    plt.savefig(
        figures_dir / "gradB_kinetic_energy.png",
        dpi=200
    )

    # --------------------------------------------------
    # Figure 5:
    # Three-dimensional particle trajectory
    # --------------------------------------------------

    fig = plt.figure(figsize=(8, 7))

    ax = fig.add_subplot(
        111,
        projection="3d"
    )

    ax.plot(
        x,
        y,
        z
    )

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")

    ax.set_title(
        "Grad-B drift: 3D trajectory"
    )

    plt.tight_layout()

    plt.savefig(
        figures_dir / "gradB_3d_trajectory.png",
        dpi=200
    )

    plt.show()


if __name__ == "__main__":
    main()
