import pyvista as pv
import pandas as pd
import numpy as np
import struct
import sys
import os

def load_trajectories(bin_path):
    particles = []
    if not os.path.exists(bin_path):
        print(f"!!! Файл не найден: {bin_path}")
        return particles
    with open(bin_path, "rb") as f:
        while True:
            header = f.read(4)
            if not header: break
            n_points = struct.unpack("I", header)[0] # Исправлено: берем элемент 0
            points = np.fromfile(f, dtype=np.float64, count=n_points * 3).reshape(-1, 3)
            particles.append(points)
    return particles

def visualize(exp_id):
    ws = "workspace/"
    plotter = pv.Plotter(title=f"Synchrotron Tracker - Exp {exp_id}")
    plotter.set_background("#0a0a0a")

    try:
        # Загрузка магнитов
        exp_df = pd.read_csv(ws + "tables/experiments.csv", skipinitialspace=True)
        lat_id = int(exp_df[exp_df['exp_id'] == exp_id]['lattice_id'].iloc[0])
        lattice_df = pd.read_csv(ws + "tables/lattice_configs.csv", skipinitialspace=True)
        hw_df = pd.read_csv(ws + "tables/hardware_library.csv", skipinitialspace=True)
        
        current_lat = lattice_df[lattice_df['lattice_id'] == lat_id]
        for _, row in current_lat.iterrows():
            dev_id = str(row['device_type_id']).strip()
            dev_info = hw_df[hw_df['device_type_id'].str.strip() == dev_id]
            if not dev_info.empty:
                stl_path = os.path.join(ws, "geometry", str(dev_info['stl_path'].iloc[0]).strip())
                if os.path.exists(stl_path):
                    mesh = pv.read(stl_path)
                    
                    mesh.scale(0.001, inplace=True) # Масштабируем из мм в метры
                    
                    mesh.rotate_z(np.degrees(row['yaw']), inplace=True)
                    mesh.rotate_y(np.degrees(row['pitch']), inplace=True)
                    mesh.rotate_x(np.degrees(row['roll']), inplace=True)
                    mesh.translate([row['x'], row['y'], row['z']], inplace=True)
                    plotter.add_mesh(mesh, color="silver", opacity=0.2, show_edges=True)
    except Exception as e:
        print(f"Ошибка загрузки геометрии: {e}")

    # Загрузка траекторий
    bin_path = f"{ws}results/trajectories_exp{exp_id}.bin"
    particles = load_trajectories(bin_path)
    if not particles:
        print("ВНИМАНИЕ: Траектории пусты.")
    else:
        colors = ["magenta", "cyan", "lime", "yellow", "red"]
        for i, p in enumerate(particles):
            if len(p) < 2: continue
            poly = pv.PolyData(p)
            cells = np.full((1, len(p) + 1), len(p), dtype=np.int_)
            cells[0, 1:] = np.arange(len(p))
            poly.lines = cells
            plotter.add_mesh(poly, color=colors[i % len(colors)], line_width=4)
            print(f"  [OK] Траектория {i+1}: {len(p)} точек.")

    plotter.show_bounds(grid='back', location='outer', xtitle='X (m)', ytitle='Y (m)', ztitle='Z (m)', color='white')
    plotter.add_axes()
    plotter.reset_camera()
    plotter.show()

if __name__ == "__main__":
    visualize(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
