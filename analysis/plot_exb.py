from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def main():
    # --------------------------------------------------
    # Paths
    # --------------------------------------------------

    script_dir = Path(__file__).resolve().parent
    project_dir = script_dir.parent

    csv_path = project_dir / "results" / "exb.csv"
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
    # Uniform crossed fields:
    #
    # E = E0 x-hat
    # B = B0 z-hat
    #
    # Normalized initial conditions:
    #
    # r(0) = (0, 0, 0)
    # v(0) = (1, 0, 0.2)
    # q = 1
    # m = 1
    # --------------------------------------------------

    q = 1.0
    mass = 1.0

    E0 = 0.1
    B0 = 1.0

    x0 = 0.0
    y0 = 0.0
    z0 = 0.0

    vx0 = 1.0
    vy0 = 0.0
    vz0 = 0.2

    # --------------------------------------------------
    # Analytical E x B drift
    #
    # v_E = (E x B) / B^2
    #
    # E = E0 x-hat
    # B = B0 z-hat
    #
    # x-hat x z-hat = -y-hat
    #
    # Therefore:
    #
    # v_E,y = -E0 / B0
    # --------------------------------------------------

    drift_velocity_theory = -E0 / B0

    # Signed cyclotron frequency.
    #
    # Keeping the sign is useful when writing the
    # analytical particle orbit.
    omega = q * B0 / mass

    # --------------------------------------------------
    # Exact analytical particle solution
    #
    # Separate the velocity into:
    #
    # v = v_E + w
    #
    # where w performs pure gyromotion around the
    # drifting guiding center.
    #
    # Since:
    #
    # v_E = (0, -E0/B0, 0)
    #
    # define:
    #
    # u_y = v_y - v_E,y
    # --------------------------------------------------

    uy0 = vy0 - drift_velocity_theory

    cos_wt = np.cos(omega * t)
    sin_wt = np.sin(omega * t)

    vx_exact = (
        vx0 * cos_wt
        + uy0 * sin_wt
    )

    uy_exact = (
        uy0 * cos_wt
        - vx0 * sin_wt
    )

    vy_exact = (
        uy_exact
        + drift_velocity_theory
    )

    vz_exact = np.full_like(t, vz0)

    # Integrating the analytical velocity gives
    # the exact particle position.

    x_exact = (
        x0
        + vx0 * sin_wt / omega
        + uy0 * (1.0 - cos_wt) / omega
    )

    y_exact = (
        y0
        + uy0 * sin_wt / omega
        + vx0 * (cos_wt - 1.0) / omega
        + drift_velocity_theory * t
    )

    z_exact = z0 + vz0 * t

    # --------------------------------------------------
    # Analytical kinetic energy
    #
    # Unlike the pure magnetic-field benchmark,
    # kinetic energy is not necessarily constant
    # instantaneously because:
    #
    # dK/dt = q E dot v
    #
    # Therefore we compare the numerical energy
    # with the exact time-dependent solution rather
    # than treating every change in K as numerical error.
    # --------------------------------------------------

    kinetic_energy_exact = 0.5 * mass * (
        vx_exact**2
        + vy_exact**2
        + vz_exact**2
    )

    # --------------------------------------------------
    # Numerical guiding center
    #
    # First subtract the E x B drift from the
    # instantaneous particle velocity:
    #
    # w = v - v_E
    #
    # For B = B0 z-hat, the gyrocenter is:
    #
    # R = r + (w x z-hat) / omega
    #
    # Therefore:
    #
    # R_x = x + w_y / omega
    # R_y = y - w_x / omega
    #
    # The guiding center should move linearly
    # in the -y direction at v_E.
    # --------------------------------------------------

    wx = vx
    wy = vy - drift_velocity_theory

    guiding_center_x = x + wy / omega
    guiding_center_y = y - wx / omega

    # --------------------------------------------------
    # Measure drift velocity from the numerical
    # guiding-center trajectory.
    #
    # Fit:
    #
    # R_y(t) = v_drift * t + constant
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

    relative_drift_error = abs(
        (
            drift_velocity_numerical
            - drift_velocity_theory
        )
        / drift_velocity_theory
    )

    # --------------------------------------------------
    # Numerical error against exact orbit
    # --------------------------------------------------

    max_x_error = np.max(
        np.abs(x - x_exact)
    )

    max_y_error = np.max(
        np.abs(y - y_exact)
    )

    max_z_error = np.max(
        np.abs(z - z_exact)
    )

    max_vx_error = np.max(
        np.abs(vx - vx_exact)
    )

    max_vy_error = np.max(
        np.abs(vy - vy_exact)
    )

    max_vz_error = np.max(
        np.abs(vz - vz_exact)
    )

    max_energy_error = np.max(
        np.abs(
            kinetic_energy
            - kinetic_energy_exact
        )
    )

    # --------------------------------------------------
    # Console diagnostics
    # --------------------------------------------------

    print("\n=== E x B drift analysis ===")
    print(f"CSV file: {csv_path}")
    print(f"Number of saved states: {len(t)}")

    print("\nTheoretical prediction:")
    print(
        f"  v_E,y = {drift_velocity_theory:.8f}"
    )

    print("\nNumerical guiding-center measurement:")
    print(
        f"  v_drift,y = {drift_velocity_numerical:.8f}"
    )
    print(
        f"  relative drift error = "
        f"{relative_drift_error:.6e}"
    )

    print("\nMaximum absolute trajectory errors:")
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

    print("\nMaximum absolute velocity errors:")
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

    print("\nEnergy comparison:")
    print(
        f"  max |K_num - K_exact| = "
        f"{max_energy_error:.6e}"
    )

    # --------------------------------------------------
    # Figure 1:
    # Particle orbit and guiding-center drift
    # --------------------------------------------------

    plt.figure(figsize=(7, 7))

    plt.plot(
        x,
        y,
        label="Numerical particle orbit"
    )

    plt.plot(
        guiding_center_x,
        guiding_center_y,
        "--",
        label="Numerical guiding center"
    )

    plt.xlabel("x")
    plt.ylabel("y")

    plt.title(
        "E x B drift: particle orbit and guiding center"
    )

    plt.axis("equal")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    plt.savefig(
        figures_dir / "exb_xy_orbit.png",
        dpi=200
    )

    # --------------------------------------------------
    # Figure 2:
    # Numerical vs analytical particle orbit
    # --------------------------------------------------

    plt.figure(figsize=(7, 7))

    plt.plot(
        x,
        y,
        label="Numerical"
    )

    plt.plot(
        x_exact,
        y_exact,
        "--",
        label="Analytical"
    )

    plt.xlabel("x")
    plt.ylabel("y")

    plt.title(
        "E x B drift: numerical vs analytical orbit"
    )

    plt.axis("equal")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    plt.savefig(
        figures_dir / "exb_numerical_vs_analytical.png",
        dpi=200
    )

    # --------------------------------------------------
    # Figure 3:
    # Guiding-center displacement
    #
    # A straight line confirms a constant
    # E x B drift velocity.
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
            "Linear fit: "
            f"v = {drift_velocity_numerical:.6f}"
        )
    )

    plt.xlabel("t")
    plt.ylabel("Guiding-center y")

    plt.title(
        "E x B drift: guiding-center motion"
    )

    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    plt.savefig(
        figures_dir / "exb_guiding_center.png",
        dpi=200
    )

    # --------------------------------------------------
    # Figure 4:
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
        label="Numerical trajectory"
    )

    ax.plot(
        x_exact,
        y_exact,
        z_exact,
        "--",
        label="Analytical trajectory"
    )

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")

    ax.set_title(
        "E x B drift: 3D trajectory"
    )

    ax.legend()
    plt.tight_layout()

    plt.savefig(
        figures_dir / "exb_3d_trajectory.png",
        dpi=200
    )

    # --------------------------------------------------
    # Figure 5:
    # Kinetic energy
    #
    # In crossed E and B fields, the magnetic field
    # still performs no work, but the electric field can:
    #
    # dK/dt = q E dot v
    #
    # Therefore K(t) need not remain constant.
    # --------------------------------------------------

    plt.figure(figsize=(8, 5))

    plt.plot(
        t,
        kinetic_energy,
        label="Numerical"
    )

    plt.plot(
        t,
        kinetic_energy_exact,
        "--",
        label="Analytical"
    )

    plt.xlabel("t")
    plt.ylabel("Kinetic energy")

    plt.title(
        "E x B drift: kinetic energy"
    )

    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    plt.savefig(
        figures_dir / "exb_kinetic_energy.png",
        dpi=200
    )

    plt.show()


if __name__ == "__main__":
    main()
