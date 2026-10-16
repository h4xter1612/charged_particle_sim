from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def main():
    # --------------------------------------------------
    # Paths
    # --------------------------------------------------

    script_dir = Path(__file__).resolve().parent
    project_dir = script_dir.parent

    csv_path = project_dir / "results" / "curvature.csv"
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

    position = np.column_stack((x, y, z))
    velocity = np.column_stack((vx, vy, vz))

    # --------------------------------------------------
    # Benchmark parameters
    #
    # Magnetic field:
    #
    # B = B0 * e_phi
    #
    # where:
    #
    # e_phi = (-y/R, x/R, 0)
    #
    # The magnetic-field magnitude is constant:
    #
    # |B| = B0
    #
    # but the field lines are circular and therefore
    # have finite curvature.
    #
    # Initial conditions:
    #
    # r(0) = (10, 0, 0)
    # v(0) = (0.2, 1.0, 0)
    #
    # At (10, 0, 0):
    #
    # B points in +y,
    # therefore:
    #
    # v_parallel = 1.0
    # v_perp     = 0.2
    # --------------------------------------------------

    q = 1.0
    mass = 1.0
    B0 = 1.0

    # Signed cyclotron frequency.
    #
    # Since |B| is constant in this benchmark,
    # Omega is also constant.

    omega = q * B0 / mass

    # --------------------------------------------------
    # Local cylindrical geometry
    # --------------------------------------------------

    radius = np.sqrt(
        x**2 + y**2
    )

    if np.any(radius <= 0.0):
        raise ValueError(
            "Curvature benchmark encountered R <= 0."
        )

    # Unit vector along the magnetic field:
    #
    # b = e_phi
    #
    # In Cartesian coordinates:
    #
    # e_phi = (-y/R, x/R, 0)

    b_hat = np.column_stack((
        -y / radius,
         x / radius,
         np.zeros_like(radius)
    ))

    # --------------------------------------------------
    # Parallel and perpendicular velocity
    #
    # v_parallel = v dot b
    #
    # v_perp = v - v_parallel b
    # --------------------------------------------------

    v_parallel = np.sum(
        velocity * b_hat,
        axis=1
    )

    v_parallel_vector = (
        v_parallel[:, np.newaxis]
        * b_hat
    )

    v_perpendicular = (
        velocity
        - v_parallel_vector
    )

    v_perpendicular_magnitude = np.linalg.norm(
        v_perpendicular,
        axis=1
    )

    # --------------------------------------------------
    # Approximate guiding-center position
    #
    # To first order:
    #
    # R_gc = r + (v_perp x b) / Omega
    #
    # For a uniform straight field this expression is
    # exact for the gyrocenter.
    #
    # Here the magnetic-field direction varies in space,
    # so it is used as a first-order guiding-center
    # approximation.
    # --------------------------------------------------

    gyroradius_vector = (
        np.cross(
            v_perpendicular,
            b_hat
        )
        / omega
    )

    guiding_center = (
        position
        + gyroradius_vector
    )

    gc_x = guiding_center[:, 0]
    gc_y = guiding_center[:, 1]
    gc_z = guiding_center[:, 2]

    gc_radius = np.sqrt(
        gc_x**2 + gc_y**2
    )

    # --------------------------------------------------
    # Guiding-center ordering
    #
    # Curvature drift theory assumes that the Larmor
    # radius is much smaller than the curvature radius:
    #
    # epsilon = rho_L / R_c << 1
    # --------------------------------------------------

    v_perp_initial = (
        v_perpendicular_magnitude[0]
    )

    larmor_radius_initial = (
        v_perp_initial
        / abs(omega)
    )

    curvature_radius_initial = (
        gc_radius[0]
    )

    epsilon = (
        larmor_radius_initial
        / curvature_radius_initial
    )

    # --------------------------------------------------
    # Theoretical curvature drift
    #
    # General expression:
    #
    # v_c =
    #     (m v_parallel^2 / q B)
    #     * (b x kappa)
    #
    # For circular field lines:
    #
    # kappa = -e_R / R
    #
    # and:
    #
    # e_phi x (-e_R) = +e_z
    #
    # Therefore:
    #
    # v_c,z =
    #     m v_parallel^2 / (q B R)
    #
    # We evaluate the leading-order prediction using
    # the initial guiding-center quantities.
    # --------------------------------------------------

    v_parallel_initial = (
        v_parallel[0]
    )

    drift_velocity_theory = (
        mass
        * v_parallel_initial**2
        / (
            q
            * B0
            * curvature_radius_initial
        )
    )

    # --------------------------------------------------
    # Numerical curvature drift
    #
    # The guiding center should move approximately
    # linearly in z:
    #
    # Z_gc(t) =
    #     Z_gc(0) + v_c,z * t
    #
    # A linear fit to the numerical guiding center
    # gives the measured curvature drift velocity.
    # --------------------------------------------------

    drift_fit = np.polyfit(
        t,
        gc_z,
        1
    )

    drift_velocity_numerical = (
        drift_fit[0]
    )

    gc_z_intercept = (
        drift_fit[1]
    )

    gc_z_fit = (
        drift_velocity_numerical * t
        + gc_z_intercept
    )

    gc_z_theory = (
        gc_z[0]
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

    initial_energy = (
        kinetic_energy[0]
    )

    relative_energy_error = (
        kinetic_energy - initial_energy
    ) / initial_energy

    max_relative_energy_error = np.max(
        np.abs(relative_energy_error)
    )

    # --------------------------------------------------
    # Additional diagnostics
    #
    # The magnetic-field magnitude is constant, but
    # because the field direction changes along the
    # orbit, the instantaneous decomposition into
    # parallel and perpendicular velocities can contain
    # finite-Larmor-radius oscillations.
    # --------------------------------------------------

    mean_gc_radius = np.mean(
        gc_radius
    )

    relative_gc_radius_variation = (
        gc_radius - mean_gc_radius
    ) / mean_gc_radius

    max_relative_gc_radius_variation = np.max(
        np.abs(relative_gc_radius_variation)
    )

    # --------------------------------------------------
    # Console diagnostics
    # --------------------------------------------------

    print("\n=== Curvature drift analysis ===")
    print(f"CSV file: {csv_path}")
    print(f"Number of saved states: {len(t)}")

    print("\nGuiding-center ordering:")
    print(
        f"  initial Larmor radius rho_L = "
        f"{larmor_radius_initial:.8f}"
    )
    print(
        f"  initial curvature radius R_c = "
        f"{curvature_radius_initial:.8f}"
    )
    print(
        f"  epsilon = rho_L / R_c = "
        f"{epsilon:.6e}"
    )

    print("\nInitial velocity decomposition:")
    print(
        f"  v_parallel = "
        f"{v_parallel_initial:.8f}"
    )
    print(
        f"  v_perp     = "
        f"{v_perp_initial:.8f}"
    )

    print("\nTheoretical prediction:")
    print(
        f"  v_curvature,z = "
        f"{drift_velocity_theory:.8f}"
    )

    print("\nNumerical guiding-center measurement:")
    print(
        f"  v_drift,z = "
        f"{drift_velocity_numerical:.8f}"
    )
    print(
        f"  relative drift error = "
        f"{relative_drift_error:.6e}"
    )

    print("\nGuiding-center geometry:")
    print(
        f"  mean guiding-center radius = "
        f"{mean_gc_radius:.8f}"
    )
    print(
        f"  max relative radius variation = "
        f"{max_relative_gc_radius_variation:.6e}"
    )

    print("\nEnergy conservation:")
    print(
        f"  max relative energy error = "
        f"{max_relative_energy_error:.6e}"
    )

    # --------------------------------------------------
    # Figure 1:
    # xy projection
    #
    # The particle follows a curved magnetic-field line
    # while gyrating around the guiding center.
    # --------------------------------------------------

    plt.figure(figsize=(7, 7))

    plt.plot(
        x,
        y,
        label="Particle orbit"
    )

    plt.plot(
        gc_x,
        gc_y,
        "--",
        label="Approximate guiding center"
    )

    plt.xlabel("x")
    plt.ylabel("y")

    plt.title(
        "Curvature drift: xy projection"
    )

    plt.axis("equal")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    plt.savefig(
        figures_dir / "curvature_xy_orbit.png",
        dpi=200
    )

    # --------------------------------------------------
    # Figure 2:
    # Guiding-center drift in z
    #
    # The slope of Z_gc(t) is the numerical curvature
    # drift velocity.
    # --------------------------------------------------

    plt.figure(figsize=(8, 5))

    plt.plot(
        t,
        gc_z,
        label="Numerical guiding center"
    )

    plt.plot(
        t,
        gc_z_fit,
        "--",
        label=(
            "Numerical linear fit: "
            f"v = {drift_velocity_numerical:.6f}"
        )
    )

    plt.plot(
        t,
        gc_z_theory,
        ":",
        label=(
            "First-order theory: "
            f"v = {drift_velocity_theory:.6f}"
        )
    )

    plt.xlabel("t")
    plt.ylabel("Guiding-center z")

    plt.title(
        "Curvature drift: guiding-center motion"
    )

    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    plt.savefig(
        figures_dir / "curvature_guiding_center.png",
        dpi=200
    )

    # --------------------------------------------------
    # Figure 3:
    # Three-dimensional trajectory
    #
    # The guiding center travels around the circular
    # field geometry while simultaneously drifting
    # in the +z direction.
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
        gc_z,
        "--",
        label="Approximate guiding center"
    )

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")

    ax.set_title(
        "Curvature drift: 3D trajectory"
    )

    ax.legend()
    plt.tight_layout()

    plt.savefig(
        figures_dir / "curvature_3d_trajectory.png",
        dpi=200
    )

    # --------------------------------------------------
    # Figure 4:
    # Kinetic energy
    #
    # A static magnetic field performs no work.
    # Therefore the total kinetic energy should remain
    # constant up to numerical integration error.
    # --------------------------------------------------

    plt.figure(figsize=(8, 5))

    plt.plot(
        t,
        kinetic_energy
    )

    plt.xlabel("t")
    plt.ylabel("Kinetic energy")

    plt.title(
        "Curvature drift: kinetic energy conservation"
    )

    plt.grid(True)
    plt.tight_layout()

    plt.savefig(
        figures_dir / "curvature_kinetic_energy.png",
        dpi=200
    )

    # --------------------------------------------------
    # Figure 5:
    # Parallel velocity
    #
    # This helps visualize finite-Larmor-radius
    # oscillations in the local velocity decomposition.
    # --------------------------------------------------

    plt.figure(figsize=(8, 5))

    plt.plot(
        t,
        v_parallel
    )

    plt.xlabel("t")
    plt.ylabel("Parallel velocity")

    plt.title(
        "Curvature drift: local parallel velocity"
    )

    plt.grid(True)
    plt.tight_layout()

    plt.savefig(
        figures_dir / "curvature_parallel_velocity.png",
        dpi=200
    )

    plt.show()


if __name__ == "__main__":
    main()
