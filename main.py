import tkinter as tk
from tkinter import ttk
from vehicleGeometry import vehicle
from RK4sim import RK4withLayeredAtmosphere

def run_eoms():
    try:
        mass = float(mass_var.get())
        coneDia = float(coned_var.get())
        noseRa = float(noser_var.get())
        deltc = float(conea_var.get())
        banka = float(banka_var.get())

        print(f"Mass: {mass}, Cone diameter: {coneDia}, Nose radius: {noseRa}, Cone angle: {deltc}, Bank angle: {banka}")

        v0 = float(vel0_var.get())
        gam0 = float(gam0_var.get())
        h0 = float(h0_var.get())

        print(f"Velocity0: {v0}, Gam0: {gam0}, Altitude0: {h0}")

        # Create vehicle object
        stardust = vehicle(mass, coneDia, noseRa, deltc, banka)

        # Run simulation
        toIntegrate = RK4withLayeredAtmosphere(stardust, v0, gam0, h0)
        simulation_result = toIntegrate.runSim()

        print("Simulation Result:", simulation_result)

    except ValueError as e:
        print("Enter valid number")

def populate_entry_corridor():
    entry_rdo.delete(0, tk.END)
    entry_rdo.insert(0, "Hello, Tkinter!")


# Create the main window
root = tk.Tk()
root.title("AE 6355: Planetary Entry, Descent, and Landing")

# Size window (w x h)
root.geometry("770x580")

# Create notebook
notebook = ttk.Notebook(root)
notebook.pack(padx=10, pady=10, fill=tk.BOTH, expand=True)

# Create tabs
tab1 = ttk.Frame(notebook)
tab2 = ttk.Frame(notebook)
tab3 = ttk.Frame(notebook)
notebook.add(tab1, text='Planar EoMs')
notebook.add(tab2, text='Entry Corridor')
notebook.add(tab3, text='Tab 3')

# Content left, tab 1
heading_label1 = tk.Label(tab1, text="Enter Vehicle Parameters", font=("Arial", 14, "bold"))
heading_label1.grid(row=0, column=0, columnspan=2, padx=10, pady=10)

label_mass = tk.Label(tab1, text="Mass:")
label_mass.grid(row=1, column=0, padx=5, pady=10, sticky=tk.W)
mass_var = tk.StringVar(value="46")
entry_mass = tk.Entry(tab1, textvariable=mass_var, width=12)
entry_mass.grid(row=1, column=1, padx=5, pady=10)
label_kg = tk.Label(tab1, text="kilograms")
label_kg.grid(row=1, column=2, padx=5, pady=10, sticky=tk.W)

label_coned = tk.Label(tab1, text="Cone diameter:")
label_coned.grid(row=2, column=0, padx=5, pady=10, sticky=tk.W)
coned_var = tk.StringVar(value="0.8128")
entry_coned = tk.Entry(tab1, textvariable=coned_var, width=12)
entry_coned.grid(row=2, column=1, padx=5, pady=10)
label_m = tk.Label(tab1, text="meters")
label_m.grid(row=2, column=2, padx=5, pady=10, sticky=tk.W)

label_noser = tk.Label(tab1, text="Nose radius:")
label_noser.grid(row=3, column=0, padx=5, pady=10, sticky=tk.W)
noser_var = tk.StringVar(value="0.2202")
entry_noser = tk.Entry(tab1, textvariable=noser_var, width=12)
entry_noser.grid(row=3, column=1, padx=5, pady=10)
label_m = tk.Label(tab1, text="meters")
label_m.grid(row=3, column=2, padx=5, pady=10, sticky=tk.W)

label_conea = tk.Label(tab1, text="Cone angle:")
label_conea.grid(row=4, column=0, padx=5, pady=10, sticky=tk.W)
conea_var = tk.StringVar(value="60")
entry_conea = tk.Entry(tab1, textvariable=conea_var, width=12)
entry_conea.grid(row=4, column=1, padx=5, pady=10)
label_d = tk.Label(tab1, text="degrees")
label_d.grid(row=4, column=2, padx=5, pady=10, sticky=tk.W)

label_banka = tk.Label(tab1, text="Bank angle:")
label_banka.grid(row=5, column=0, padx=5, pady=10, sticky=tk.W)
banka_var = tk.StringVar(value="0")
entry_banka = tk.Entry(tab1, textvariable=banka_var, width=12)
entry_banka.grid(row=5, column=1, padx=5, pady=10)
label_d = tk.Label(tab1, text="degrees")
label_d.grid(row=5, column=2, padx=5, pady=10, sticky=tk.W)

# Vertical line, tab 1
separator = ttk.Separator(tab1, orient="vertical")
separator.grid(row=0, column=3, rowspan=6, sticky="ns", padx=5)

# Content right, tab 1
heading_label2 = tk.Label(tab1, text="Enter Initial Flight Conditions", font=("Arial", 14, "bold"))
heading_label2.grid(row=0, column=4, columnspan=2, padx=5, pady=10)

label_vel0 = tk.Label(tab1, text="Initial velocity:")
label_vel0.grid(row=1, column=4, padx=5, pady=10, sticky=tk.W)
vel0_var = tk.StringVar(value="11067.63087")
entry_vel0 = tk.Entry(tab1, textvariable=vel0_var, width=12)
entry_vel0.grid(row=1, column=5, padx=5, pady=10)
label_ms = tk.Label(tab1, text="m/s")
label_ms.grid(row=1, column=6, padx=5, pady=10, sticky=tk.W)

label_gam0 = tk.Label(tab1, text="Initial flight path angle:")
label_gam0.grid(row=2, column=4, padx=5, pady=10, sticky=tk.W)
gam0_var = tk.StringVar(value="-10")
entry_gam0 = tk.Entry(tab1, textvariable=gam0_var, width=12)
entry_gam0.grid(row=2, column=5, padx=5, pady=10)
label_d = tk.Label(tab1, text="degrees")
label_d.grid(row=2, column=6, padx=5, pady=10, sticky=tk.W)

label_h0 = tk.Label(tab1, text="Initial altitude:")
label_h0.grid(row=3, column=4, padx=5, pady=10, sticky=tk.W)
h0_var = tk.StringVar(value="120000")
entry_h0 = tk.Entry(tab1, textvariable=h0_var, width=12)
entry_h0.grid(row=3, column=5, padx=5, pady=10)
label_m = tk.Label(tab1, text="meters")
label_m.grid(row=3, column=6, padx=5, pady=10, sticky=tk.W)

# Run button, tab 1
button_run = tk.Button(tab1, text="Run Simulation", command=run_eoms)
button_run.grid(row=6, column=2, columnspan=3, padx=10, pady=10, sticky=tk.NSEW)

# Content left, tab 2
heading_label3 = tk.Label(tab2, text="Enter Initial Orbit", font=("Arial", 14, "bold"))
heading_label3.grid(row=0, column=0, columnspan=2, padx=10, pady=10)

label_e = tk.Label(tab2, text="Eccenticity:")
label_e.grid(row=1, column=0, padx=5, pady=10, sticky=tk.W)
e_var = tk.StringVar(value="0.1")
entry_e = tk.Entry(tab2, textvariable=e_var, width=12)
entry_e.grid(row=1, column=1, padx=5, pady=10)

label_pa = tk.Label(tab2, text="Periapsis Altitude:")
label_pa.grid(row=2, column=0, padx=5, pady=10, sticky=tk.W)
pa_var = tk.StringVar(value="400000")
entry_pa = tk.Entry(tab2, textvariable=pa_var, width=12)
entry_pa.grid(row=2, column=1, padx=5, pady=10)
label_meters = tk.Label(tab2, text="meters")
label_meters.grid(row=2, column=2, padx=5, pady=10, sticky=tk.W)

heading_label3 = tk.Label(tab2, text="Enter Deorbit Parameters", font=("Arial", 14, "bold"))
heading_label3.grid(row=3, column=0, columnspan=2, padx=10, pady=10)

label_td = tk.Label(tab2, text="Departure true anomaly:") # theta_D
label_td.grid(row=4, column=0, padx=5, pady=10, sticky=tk.W)
td_var = tk.StringVar(value="180")
entry_td = tk.Entry(tab2, textvariable=td_var, width=12)
entry_td.grid(row=4, column=1, padx=5, pady=10)
label_deg = tk.Label(tab2, text="degrees")
label_deg.grid(row=4, column=2, padx=5, pady=10, sticky=tk.W)

label_nmax = tk.Label(tab2, text="Maximum deceleration:")
label_nmax.grid(row=5, column=0, padx=5, pady=10, sticky=tk.W)
nmax_var = tk.StringVar(value="10")
entry_nmax = tk.Entry(tab2, textvariable=nmax_var, width=12)
entry_nmax.grid(row=5, column=1, padx=5, pady=10)
label_gs = tk.Label(tab2, text="g's")
label_gs.grid(row=5, column=2, padx=5, pady=10, sticky=tk.W)

# Vertical line, tab 2
separator = ttk.Separator(tab2, orient="vertical")
separator.grid(row=0, column=3, rowspan=10, sticky="ns", padx=5)

# Content right, tab 2
heading_label2 = tk.Label(tab2, text="OUTPUT: Overshoot", font=("Arial", 14, "bold"))
heading_label2.grid(row=0, column=4, columnspan=2, padx=5, pady=10)

label_rdo = tk.Label(tab2, text="Deorbit position:")
label_rdo.grid(row=1, column=4, padx=5, pady=10, sticky=tk.W)
rdo_var = tk.StringVar(value="")
entry_rdo = tk.Entry(tab2, textvariable=rdo_var, width=12)
entry_rdo.grid(row=1, column=5, padx=5, pady=10)
label_ms = tk.Label(tab2, text="meters")
label_ms.grid(row=1, column=6, padx=5, pady=10, sticky=tk.W)

label_dvo = tk.Label(tab2, text="Change in velocity:")
label_dvo.grid(row=2, column=4, padx=5, pady=10, sticky=tk.W)
dvo_var = tk.StringVar(value="")
entry_dvo = tk.Entry(tab2, textvariable=dvo_var, width=12)
entry_dvo.grid(row=2, column=5, padx=5, pady=10)
label_ms = tk.Label(tab2, text="m/s")
label_ms.grid(row=2, column=6, padx=5, pady=10, sticky=tk.W)

label_rvo = tk.Label(tab2, text="Re-entry velocity:")
label_rvo.grid(row=3, column=4, padx=5, pady=10, sticky=tk.W)
rvo_var = tk.StringVar(value="")
entry_rvo = tk.Entry(tab2, textvariable=rvo_var, width=12)
entry_rvo.grid(row=3, column=5, padx=5, pady=10)
label_ms = tk.Label(tab2, text="m/s")
label_ms.grid(row=3, column=6, padx=5, pady=10, sticky=tk.W)

label_rgo = tk.Label(tab2, text="Re-entry gamma:")
label_rgo.grid(row=4, column=4, padx=5, pady=10, sticky=tk.W)
rgo_var = tk.StringVar(value="")
entry_rgo = tk.Entry(tab2, textvariable=rgo_var, width=12)
entry_rgo.grid(row=4, column=5, padx=5, pady=10)
label_deg = tk.Label(tab2, text="degrees")
label_deg.grid(row=4, column=6, padx=5, pady=10, sticky=tk.W)

heading_label2 = tk.Label(tab2, text="OUTPUT: Undershoot", font=("Arial", 14, "bold"))
heading_label2.grid(row=5, column=4, columnspan=2, padx=5, pady=10)

label_rdu = tk.Label(tab2, text="Deorbit position:")
label_rdu.grid(row=6, column=4, padx=5, pady=10, sticky=tk.W)
rdu_var = tk.StringVar(value="")
entry_rdu = tk.Entry(tab2, textvariable=rdu_var, width=12)
entry_rdu.grid(row=6, column=5, padx=5, pady=10)
label_ms = tk.Label(tab2, text="meters")
label_ms.grid(row=6, column=6, padx=5, pady=10, sticky=tk.W)

label_dvu = tk.Label(tab2, text="Change in velocity:")
label_dvu.grid(row=7, column=4, padx=5, pady=10, sticky=tk.W)
dvu_var = tk.StringVar(value="")
entry_dvu = tk.Entry(tab2, textvariable=dvu_var, width=12)
entry_dvu.grid(row=7, column=5, padx=5, pady=10)
label_ms = tk.Label(tab2, text="m/s")
label_ms.grid(row=7, column=6, padx=5, pady=10, sticky=tk.W)

label_rvu = tk.Label(tab2, text="Re-entry velocity:")
label_rvu.grid(row=8, column=4, padx=5, pady=10, sticky=tk.W)
rvu_var = tk.StringVar(value="")
entry_rvu = tk.Entry(tab2, textvariable=rvu_var, width=12)
entry_rvu.grid(row=8, column=5, padx=5, pady=10)
label_ms = tk.Label(tab2, text="m/s")
label_ms.grid(row=8, column=6, padx=5, pady=10, sticky=tk.W)

label_rgu = tk.Label(tab2, text="Re-entry gamma:")
label_rgu.grid(row=9, column=4, padx=5, pady=10, sticky=tk.W)
rgu_var = tk.StringVar(value="")
entry_rgu = tk.Entry(tab2, textvariable=rgu_var, width=12)
entry_rgu.grid(row=9, column=5, padx=5, pady=10)
label_deg = tk.Label(tab2, text="degrees")
label_deg.grid(row=9, column=6, padx=5, pady=10, sticky=tk.W)

# Run button
button_run = tk.Button(tab2, text="Populate Entry Corridor", command=populate_entry_corridor)
button_run.grid(row=7, column=0, columnspan=2, padx=10, pady=10, sticky=tk.NSEW)


# Content for Tab 3
label3 = tk.Label(tab3, text="Hello, Tab 3!")
label3.pack(padx=10, pady=10)

# Run event loop
root.mainloop()
