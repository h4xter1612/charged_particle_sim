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
    # Electric field:
    #
    # E(t) = E0 cos(omega_E t) x-hat
    #
    # Magnetic field:
    #
    # B = B0 z-hat
    #
    # Particle:
    #
    # q = 1
    # m = 1
    #
    # Initial velocity:
    #
    # v(0) = (0, -0.2, 0.2)
    #
    # At t = 0, the perpendicular velocity is chosen
    # to match the initial E x B drift:
    #
    # v_E(0) = -0.2 y-hat
    # --------------------------------------------------

    q = 1.0
    mass = 1.0

    E0 = 0.2
    B0 = 1.0
    omega_E = 0.05

    # Signed cyclotron frequency:
    #
    # Omega_c = q B / m

    omega_c = q * B0 / mass

    # --------------------------------------------------
    # Time-dependent electric field
    # --------------------------------------------------

    Ex = (
        E0
        * np.cos(omega_E * t)
    )

    dEx_dt = (
        -E0
        * omega_E
        * np.sin(omega_E * t)
    )

    # --------------------------------------------------
    # Theoretical E x B drift
    #
    # v_E = E x B / B^2
    #
    # For:
    #
    # E = Ex x-hat
    # B = B0 z-hat
    #
    # x-hat x z-hat = -y-hat
    #
    # therefore:
    #
    # v_E,y = -Ex / B0
    # --------------------------------------------------

    v_exb_y_theory = (
        -Ex / B0
    )

    # --------------------------------------------------
    # Theoretical polarization drift
    #
    # For a uniform magnetic field and a slowly
    # time-varying perpendicular electric field:
    #
    # v_p =
    #     m / (q B^2) * dE_perp/dt
    #
    # Here dE/dt points entirely along x:
    #
    # v_p,x =
    #     m / (q B0^2) * dEx/dt
    # --------------------------------------------------

    v_pol_x_theory = (
        mass
        / (q * B0**2)
        * dEx_dt
    )

    # --------------------------------------------------
    # Adiabatic ordering
    #
    # Polarization drift is naturally small compared
    # with E x B drift when:
    #
    # omega_E / Omega_c << 1
    # --------------------------------------------------

    temporal_ordering = (
        omega_E / abs(omega_c)
    )

    exb_velocity_amplitude = (
        E0 / B0
    )

    polarization_velocity_amplitude = (
        mass
        * E0
        * omega_E
        / (abs(q) * B0**2)
    )

    # --------------------------------------------------
    # Numerical guiding center
    #
    # First remove the instantaneous E x B motion from
    # the particle velocity:
    #
    # w = v - v_E
    #
    # For B = B0 z-hat:
    #
    # R = r + (w x b) / Omega_c
    #
    # Since v_E only has a y component:
    #
    # w_x = vx
    # w_y = vy - v_E,y
    #
    # and therefore:
    #
    # R_x = x + w_y / Omega_c
    # R_y = y - w_x / Omega_c
    #
    # For this spatially uniform E(t), this transformed
    # center cleanly separates the gyromotion from the
    # slow E x B and polarization motions.
    # --------------------------------------------------

    wx = vx
    wy = vy - v_exb_y_theory

    gc_x = (
        x + wy / omega_c
    )

    gc_y = (
        y - wx / omega_c
    )

    gc_z = z

    # --------------------------------------------------
    # Theoretical guiding-center trajectory
    #
    # From:
    #
    # dR_x/dt = v_p,x
    #
    # integrating gives:
    #
    # R_x(t) =
    # R_x(0)
    # + m/(q B^2) [E(t) - E(0)]
    #
    # Also:
    #
    # dR_y/dt = v_E,y
    #
    # which gives:
    #
    # R_y(t) =
    # R_y(0)
    # - E0/(B0 omega_E) sin(omega_E t)
    # --------------------------------------------------

    gc_x_theory = (
        gc_x[0]
        + mass
        / (q * B0**2)
        * (Ex - Ex[0])
    )

    gc_y_theory = (
        gc_y[0]
        - E0
        / (B0 * omega_E)
        * np.sin(omega_E * t)
    )

    # Parallel motion is unaffected by the perpendicular
    # electric and magnetic fields.

    gc_z_theory = (
        gc_z[0]
        + vz[0] * t
    )

    # --------------------------------------------------
    # Numerical guiding-center velocities
    #
    # Differentiate the reconstructed guiding-center
    # trajectory using second-order finite differences.
    # --------------------------------------------------

    gc_vx_numerical = np.gradient(
        gc_x,
        t,
        edge_order=2
    )

    gc_vy_numerical = np.gradient(
        gc_y,
        t,
        edge_order=2
    )

    # --------------------------------------------------
    # Guiding-center position errors
    # --------------------------------------------------

    max_gc_x_error = np.max(
        np.abs(
            gc_x - gc_x_theory
        )
    )

    max_gc_y_error = np.max(
        np.abs(
            gc_y - gc_y_theory
        )
    )

    # --------------------------------------------------
    # Drift-velocity errors
    # --------------------------------------------------

    max_pol_velocity_error = np.max(
        np.abs(
            gc_vx_numerical
            - v_pol_x_theory
        )
    )

    max_exb_velocity_error = np.max(
        np.abs(
            gc_vy_numerical
            - v_exb_y_theory
        )
    )

    # --------------------------------------------------
    # Energy / work consistency
    #
    # Because E != 0, kinetic energy is not expected
    # to remain constant.
    #
    # The work-energy equation is:
    #
    # dK/dt = q E dot v
    #
    # Here:
    #
    # dK/dt = q Ex vx
    #
    # Numerically integrate the electric power and
    # compare it with the simulated kinetic energy.
    # --------------------------------------------------

    electric_power = (
        q * Ex * vx
    )

    dt_values = np.diff(t)

    work_increments = (
        0.5
        * (
            electric_power[:-1]
            + electric_power[1:]
        )
        * dt_values
    )

    cumulative_work = np.concatenate((
        [0.0],
        np.cumsum(work_increments)
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

    print("\nTheoretical drift amplitudes:")
    print(
        f"  max |v_E| = "
        f"{exb_velocity_amplitude:.8f}"
    )
    print(
        f"  max |v_p| = "
        f"{polarization_velocity_amplitude:.8f}"
    )

    print("\nGuiding-center position comparison:")
    print(
        f"  max |R_x,num - R_x,theory| = "
        f"{max_gc_x_error:.6e}"
    )
    print(
        f"  max |R_y,num - R_y,theory| = "
        f"{max_gc_y_error:.6e}"
    )

    print("\nGuiding-center velocity comparison:")
    print(
        f"  max |v_p,x,num - v_p,x,theory| = "
        f"{max_pol_velocity_error:.6e}"
    )
    print(
        f"  max |v_E,y,num - v_E,y,theory| = "
        f"{max_exb_velocity_error:.6e}"
    )

    print("\nWork-energy consistency:")
    print(
        f"  max |K_num - (K0 + W_E)| = "
        f"{max_work_energy_error:.6e}"
    )

    # --------------------------------------------------
    # Figure 1:
    # Particle trajectory and guiding center in xy
    #
    # The fast particle gyromotion is superimposed on
    # the much slower motion of the guiding center.
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
        label="Numerical guiding center"
    )

    plt.plot(
        gc_x_theory,
        gc_y_theory,
        ":",
        label="Drift theory"
    )

    plt.xlabel("x")
    plt.ylabel("y")

    plt.title(
        "Polarization drift: particle and guiding-center motion"
    )

    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    plt.savefig(
        figures_dir / "polarization_xy_motion.png",
        dpi=200
    )

    # --------------------------------------------------
    # Figure 2:
    # Polarization motion in x
    #
    # The x component of the guiding-center velocity
    # is the polarization drift.
    # --------------------------------------------------

    plt.figure(figsize=(8, 5))

    plt.plot(
        t,
        gc_vx_numerical,
        label="Numerical guiding-center velocity"
    )

    plt.plot(
        t,
        v_pol_x_theory,
        "--",
        label="Polarization-drift theory"
    )

    plt.xlabel("t")
    plt.ylabel("Velocity in x")

    plt.title(
        "Polarization drift velocity"
    )

    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    plt.savefig(
        figures_dir / "polarization_velocity_x.png",
        dpi=200
    )

    # --------------------------------------------------
    # Figure 3:
    # Time-dependent E x B drift in y
    # --------------------------------------------------

    plt.figure(figsize=(8, 5))

    plt.plot(
        t,
        gc_vy_numerical,
        label="Numerical guiding-center velocity"
    )

    plt.plot(
        t,
        v_exb_y_theory,
        "--",
        label="E x B theory"
    )

    plt.xlabel("t")
    plt.ylabel("Velocity in y")

    plt.title(
        "Time-dependent E x B drift"
    )

    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    plt.savefig(
        figures_dir / "polarization_exb_velocity_y.png",
        dpi=200
    )

    # --------------------------------------------------
    # Figure 4:
    # Guiding-center coordinates
    #
    # R_y has the much larger E x B motion.
    # R_x contains the smaller polarization response.
    # --------------------------------------------------

    plt.figure(figsize=(8, 6))

    plt.plot(
        t,
        gc_x,
        label="R_x numerical"
    )

    plt.plot(
        t,
        gc_x_theory,
        "--",
        label="R_x theory"
    )

    plt.plot(
        t,
        gc_y,
        label="R_y numerical"
    )

    plt.plot(
        t,
        gc_y_theory,
        "--",
        label="R_y theory"
    )

    plt.xlabel("t")
    plt.ylabel("Guiding-center position")

    plt.title(
        "Polarization drift: guiding-center coordinates"
    )

    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    plt.savefig(
        figures_dir / "polarization_guiding_center.png",
        dpi=200
    )

    # --------------------------------------------------
    # Figure 5:
    # Work-energy theorem
    #
    # Magnetic fields perform no work, but the
    # time-dependent electric field exchanges energy
    # with the particle.
    # --------------------------------------------------

    plt.figure(figsize=(8, 5))

    plt.plot(
        t,
        kinetic_energy,
        label="Numerical kinetic energy"
    )

    plt.plot(
        t,
        energy_from_work,
        "--",
        label="K0 + electric-field work"
    )

    plt.xlabel("t")
    plt.ylabel("Energy")

    plt.title(
        "Polarization benchmark: work-energy consistency"
    )

    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    plt.savefig(
        figures_dir / "polarization_energy.png",
        dpi=200
    )

    # --------------------------------------------------
    # Figure 6:
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
        gc_z,
        "--",
        label="Numerical guiding center"
    )

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")

    ax.set_title(
        "Polarization drift: 3D trajectory"
    )

    ax.legend()
    plt.tight_layout()

    plt.savefig(
        figures_dir / "polarization_3d_trajectory.png",
        dpi=200
    )

    plt.show()


if __name__ == "__main__":
    main()
