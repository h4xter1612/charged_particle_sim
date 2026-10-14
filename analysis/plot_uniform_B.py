from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt


def main():
    # -----------------------------------------
    # Paths
    # -----------------------------------------
    script_dir = Path(__file__).resolve().parent
    project_dir = script_dir.parent
    csv_path = project_dir / "results" / "uniform_B.csv"
    figures_dir = project_dir / "results" / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    # -----------------------------------------
    # Load data
    # -----------------------------------------
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

    # -----------------------------------------
    # Analytical solution for the benchmark
    #
    # Normalized case:
    # q = 1, m = 1, B = (0,0,1)
    # r(0) = (0,0,0)
    # v(0) = (1,0,0.2)
    #
    # Expected solution:
    # x(t)  = sin(t)
    # y(t)  = cos(t) - 1
    # z(t)  = 0.2 t
    # vx(t) = cos(t)
    # vy(t) = -sin(t)
    # vz(t) = 0.2
    # -----------------------------------------
    x_exact = np.sin(t)
    y_exact = np.cos(t) - 1.0
    z_exact = 0.2 * t

    vx_exact = np.cos(t)
    vy_exact = -np.sin(t)
    vz_exact = np.full_like(t, 0.2)

    kinetic_energy_exact = np.full_like(t, 0.52)

    # -----------------------------------------
    # Error analysis
    # -----------------------------------------
    x_error = x - x_exact
    y_error = y - y_exact
    z_error = z - z_exact

    vx_error = vx - vx_exact
    vy_error = vy - vy_exact
    vz_error = vz - vz_exact

    energy_error = kinetic_energy - kinetic_energy_exact

    max_abs_x_error = np.max(np.abs(x_error))
    max_abs_y_error = np.max(np.abs(y_error))
    max_abs_z_error = np.max(np.abs(z_error))

    max_abs_vx_error = np.max(np.abs(vx_error))
    max_abs_vy_error = np.max(np.abs(vy_error))
    max_abs_vz_error = np.max(np.abs(vz_error))

    max_abs_energy_error = np.max(np.abs(energy_error))
    max_relative_energy_error = np.max(
        np.abs(energy_error) / np.abs(kinetic_energy_exact)
    )

    print("\n=== Uniform-B analysis ===")
    print(f"CSV file: {csv_path}")
    print(f"Number of saved states: {len(t)}")
    print()
    print("Maximum absolute trajectory errors:")
    print(f"  max |x_num - x_exact|  = {max_abs_x_error:.6e}")
    print(f"  max |y_num - y_exact|  = {max_abs_y_error:.6e}")
    print(f"  max |z_num - z_exact|  = {max_abs_z_error:.6e}")
    print()
    print("Maximum absolute velocity errors:")
    print(f"  max |vx_num - vx_exact| = {max_abs_vx_error:.6e}")
    print(f"  max |vy_num - vy_exact| = {max_abs_vy_error:.6e}")
    print(f"  max |vz_num - vz_exact| = {max_abs_vz_error:.6e}")
    print()
    print("Energy errors:")
    print(f"  max absolute energy error = {max_abs_energy_error:.6e}")
    print(f"  max relative energy error = {max_relative_energy_error:.6e}")

    # -----------------------------------------
    # Figure 1: xy orbit
    # -----------------------------------------
    plt.figure(figsize=(6, 6))
    plt.plot(x, y, label="Numerical orbit")
    plt.plot(x_exact, y_exact, "--", label="Analytical orbit")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("Uniform-B benchmark: xy orbit")
    plt.axis("equal")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(figures_dir / "uniform_B_xy_orbit.png", dpi=200)

    # -----------------------------------------
    # Figure 2: 3D trajectory
    # -----------------------------------------
    fig = plt.figure(figsize=(7, 6))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot(x, y, z, label="Numerical trajectory")
    ax.plot(x_exact, y_exact, z_exact, "--", label="Analytical trajectory")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.set_title("Uniform-B benchmark: 3D trajectory")
    ax.legend()
    plt.tight_layout()
    plt.savefig(figures_dir / "uniform_B_3d_trajectory.png", dpi=200)

    # -----------------------------------------
    # Figure 3: kinetic energy vs time
    # -----------------------------------------
    plt.figure(figsize=(7, 5))
    plt.plot(t, kinetic_energy, label="Numerical kinetic energy")
    plt.plot(t, kinetic_energy_exact, "--", label="Analytical kinetic energy")
    plt.xlabel("t")
    plt.ylabel("K")
    plt.title("Uniform-B benchmark: kinetic energy")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(figures_dir / "uniform_B_kinetic_energy.png", dpi=200)

    # -----------------------------------------
    # Figure 4: x(t), y(t), z(t)
    # -----------------------------------------
    plt.figure(figsize=(8, 6))
    plt.plot(t, x, label="x numerical")
    plt.plot(t, x_exact, "--", label="x analytical")
    plt.plot(t, y, label="y numerical")
    plt.plot(t, y_exact, "--", label="y analytical")
    plt.plot(t, z, label="z numerical")
    plt.plot(t, z_exact, "--", label="z analytical")
    plt.xlabel("t")
    plt.ylabel("Position components")
    plt.title("Uniform-B benchmark: position vs time")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(figures_dir / "uniform_B_position_vs_time.png", dpi=200)

    # -----------------------------------------
    # Figure 5: error vs time
    # -----------------------------------------
    plt.figure(figsize=(8, 6))
    plt.plot(t, x_error, label="x error")
    plt.plot(t, y_error, label="y error")
    plt.plot(t, z_error, label="z error")
    plt.xlabel("t")
    plt.ylabel("Numerical - analytical")
    plt.title("Uniform-B benchmark: trajectory errors")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(figures_dir / "uniform_B_trajectory_errors.png", dpi=200)

    plt.show()


if __name__ == "__main__":
    main()
