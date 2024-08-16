# Quantum Mechanics, Path integral in potential (interaction)

import numpy as np
from matplotlib.figure import Figure
import matplotlib.animation as animation
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
import tkinter as tk
from tkinter import ttk
from mpl_toolkits.mplot3d import proj3d


def update_potential():
    global qtm0a, qtm0b
    global theta_potential_by_a, rot_potential_by_a
    global theta_potential_by_b, rot_potential_by_b
    theta_potential_by_a_buffer = theta_potential_by_a * 0.
    theta_potential_by_b_buffer = theta_potential_by_b * 0.
    # Quantum A
    for i in range(len(x)):
        dx_squared = (x[i] - x) ** 2
        theta_potential_by_a_buffer = theta_potential_by_a_buffer + dx_squared * np.abs(qtm0a[i])
    theta_potential_by_a = theta_potential_by_a_buffer * potential_coefficient * sign_potential
    rot_potential_by_a = np.sin(theta_potential_by_a) + np.cos(theta_potential_by_a) * 1j
    plt_potential_by_a.set_ydata(theta_potential_by_a)
    # Quantum B
    for j in range(len(x)):
        dx_squared = (x[j] - x) ** 2
        theta_potential_by_b_buffer = theta_potential_by_b_buffer + dx_squared * np.abs(qtm0b[j])
    theta_potential_by_b = theta_potential_by_b_buffer * potential_coefficient * sign_potential
    rot_potential_by_b = np.sin(theta_potential_by_b) + np.cos(theta_potential_by_b) * 1j
    plt_potential_by_b.set_ydata(theta_potential_by_b)


def apply_path_integral():
    global qtm0a
    global ya, za, plt_qtm0a
    global qtm0b
    global yb, zb, plt_qtm0b
    qtm0a_buffer = qtm0a * 0.
    qtm0b_buffer = qtm0b * 0.
    update_potential()
    # Quantum A
    for i in range(len(x)):
        dx_squared = (x[i] - x) ** 2
        if t != 0.:
            theta = np.mod(2. * np.pi * mass * dx_squared / (2. * h * t), (2. * np.pi))
        else:
            theta = 0.
            print('Error: divided by zero')
        rot = np.sin(theta) + np.cos(theta) * 1j
        qtm0a_buffer = qtm0a_buffer + qtm0a[i] * rot
    qtm0a_buffer = qtm0a_buffer * rot_potential_by_b
    qtm0a = qtm0a_buffer / np.sum(np.abs(qtm0a_buffer))
    ya = qtm0a.imag
    za = qtm0a.real
    plt_qtm0a.set_xdata(x)
    plt_qtm0a.set_ydata(ya)
    plt_qtm0a.set_3d_properties(za)
    # Quantum B
    for j in range(len(x)):
        dx_squared = (x[j] - x) ** 2
        if t != 0.:
            theta = np.mod(2. * np.pi * mass * dx_squared / (2. * h * t), (2. * np.pi))
        else:
            theta = 0.
            print('Error: divided by zero')
        rot = np.sin(theta) + np.cos(theta) * 1j
        qtm0b_buffer = qtm0b_buffer + qtm0b[j] * rot
    qtm0b_buffer = qtm0b_buffer * rot_potential_by_a
    qtm0b = qtm0b_buffer / np.sum(np.abs(qtm0b_buffer))
    yb = qtm0b.imag
    zb = qtm0b.real
    plt_qtm0b.set_xdata(x)
    plt_qtm0b.set_ydata(yb)
    plt_qtm0b.set_3d_properties(zb)
    # Potential
    update_potential()


def update_qtm0a():
    global gaussian0a, amplitude0a
    global z0a, y0a, y0a, qtm0a
    global ya, za, plt_qtm0a
    gaussian0a = 1 / (np.sqrt(2 * np.pi) * sigma0a) * np.exp(- (x - mu0a) ** 2 / (2 * sigma0a ** 2))
    amplitude0a = gaussian0a / np.sum(np.abs(gaussian0a))
    z0a = np.sin(k0a * x) * amplitude0a
    y0a = np.cos(k0a * x) * 1j * amplitude0a
    qtm0a = z0a + y0a
    ya = qtm0a.imag
    za = qtm0a.real
    plt_qtm0a.set_xdata(x)
    plt_qtm0a.set_ydata(ya)
    plt_qtm0a.set_3d_properties(za)


def update_qtm0b():
    global gaussian0b, amplitude0b
    global z0b, y0b, y0b, qtm0b
    global yb, zb, plt_qtm0b
    gaussian0b = 1 / (np.sqrt(2 * np.pi) * sigma0b) * np.exp(- (x - mu0b) ** 2 / (2 * sigma0b ** 2))
    amplitude0b = gaussian0b / np.sum(np.abs(gaussian0b))
    z0b = np.sin(k0b * x) * amplitude0b
    y0b = np.cos(k0b * x) * 1j * amplitude0b
    qtm0b = z0b + y0b
    yb = qtm0b.imag
    zb = qtm0b.real
    plt_qtm0b.set_xdata(x)
    plt_qtm0b.set_ydata(yb)
    plt_qtm0b.set_3d_properties(zb)


def set_potential_sign(value):
    global sign_potential
    sign_potential = float(value)
    update_potential()


def set_potential_coefficient(value):
    global potential_coefficient
    potential_coefficient = 1. / float(value)
    update_potential()


def set_k0a(value):
    global k0a
    k0a = float(value)
    update_qtm0a()


def set_x0a(value):
    global mu0a
    mu0a = float(value)
    update_qtm0a()


def set_sigma0a(value):
    global sigma0a
    sigma0a = float(value)
    update_qtm0a()


def set_k0b(value):
    global k0b
    k0b = float(value)
    update_qtm0b()


def set_x0b(value):
    global mu0b
    mu0b = float(value)
    update_qtm0b()


def set_sigma0b(value):
    global sigma0b
    sigma0b = float(value)
    update_qtm0b()


def set_tm(value):
    global tm, t
    tm = float(value)
    t = tm * np.power(10., te)


def set_te(value):
    global te, t
    te = float(value)
    t = tm * np.power(10., te)


# Animation control
def step():
    global cnt
    apply_path_integral()
    cnt += 1


def reset():
    global is_play, cnt
    global cnt_step
    global k0a, sigma0a, mu0a
    global k0b, sigma0b, mu0b
    is_play = False
    cnt = 0
    # Quantum A
    k0a = k0a_init
    sigma0a = sigma0a_init
    mu0a = mu0a_init
    update_qtm0a()
    # Quantum B
    k0b = k0b_init
    sigma0b = sigma0b_init
    mu0b = mu0b_init
    update_qtm0b()


def switch():
    global is_play
    if is_play:
        is_play = False
    else:
        is_play = True


def update(f):
    global cnt
    txt_step.set_text("Delta t * " + str(cnt))
    if is_play:
        cnt += 1
        step()


# Global variables
# Animation control
cnt = 0
is_play = False
cnt_step = 1

# Data structure
range_x = 1000
x = np.arange(- range_x, range_x)

range_quantum_yz = 0.02
range_potential = 100

# Quantum parameter
mass = 9.1093837015 * 1.0E-31   # Unit: kg
tm, te = 1., 7
t = tm * np.power(10., te)  # Unit: s
h = 6.62607015 * 1.0E-34    # Unit: m2 kg / s

pc = 10000
potential_coefficient = 1 / pc

sign_potential = 0.

# Quantum
# Quantum 0a
k0a_init = 0.0
sigma0a_init = 20.
mu0a_init = - range_x / 2.

k0a = k0a_init
sigma0a = sigma0a_init
mu0a = mu0a_init
gaussian0a = 1 / (np.sqrt(2 * np.pi) * sigma0a) * np.exp(- (x - mu0a) ** 2 / (2 * sigma0a ** 2))
amplitude0a = gaussian0a / np.sum(np.abs(gaussian0a))
z0a = np.sin(k0a * x) * amplitude0a
y0a = np.cos(k0a * x) * 1j * amplitude0a
qtm0a = z0a + y0a

# Quantum 0b
k0b_init = - 0.0
sigma0b_init = 20.
mu0b_init = range_x / 2.

k0b = k0b_init
sigma0b = sigma0b_init
mu0b = mu0b_init
gaussian0b = 1 / (np.sqrt(2 * np.pi) * sigma0b) * np.exp(- (x - mu0b) ** 2 / (2 * sigma0b ** 2))
amplitude0b = gaussian0b / np.sum(np.abs(gaussian0b))
z0b = np.sin(k0b * x) * amplitude0b
y0b = np.cos(k0b * x) * 1j * amplitude0b
qtm0b = z0b + y0b

# Potential
theta_potential_by_a = 0. * x
rot_potential_by_a = np.sin(theta_potential_by_a) + np.cos(theta_potential_by_a) * 1j

theta_potential_by_b = 0. * x
rot_potential_by_b = np.sin(theta_potential_by_b) + np.cos(theta_potential_by_b) * 1j

# Generate figure and axes
title_tk = "Quantum Mechanics, Path integral (interaction)"
title_ax0 = "Quantum"
title_ax1 = "Potential:"

x_min0 = - range_x
x_max0 = range_x
y_min0 = - range_quantum_yz
y_max0 = range_quantum_yz
z_min0 = - range_quantum_yz
z_max0 = range_quantum_yz

x_min1 = - range_x
x_max1 = range_x
y_min1 = - range_potential
y_max1 = range_potential

fig = Figure()
ax0 = fig.add_subplot(121, projection='3d')
ax0.set_box_aspect((6, 2, 2))
ax0.grid()
ax0.set_title(title_ax0)
ax0.set_xlabel('x')
ax0.set_ylabel('y')
ax0.set_zlabel('z')
ax0.set_xlim(x_min0, x_max0)
ax0.set_ylim(y_min0, y_max0)
ax0.set_zlim(z_min0, z_max0)

ax1 = fig.add_subplot(122)
ax1.grid()
ax1.set_title(title_ax1)
ax1.set_xlabel('x')
ax1.set_ylabel('Additional phase')
ax1.set_xlim(x_min1, x_max1)
ax1.set_ylim(y_min1, y_max1)

# Text items
txt_step = ax0.text2D(x_min0, y_max0, "Delta t * " + str(cnt))
xz, yz, _ = proj3d.proj_transform(x_min0, y_max0, z_max0, ax0.get_proj())
txt_step.set_position((xz, yz))

# Plot items
# Quantum 0a
ya = qtm0a.imag
za = qtm0a.real
plt_qtm0a, = ax0.plot(x, ya, za, color='blue', linewidth=0.5, label='Quantum A')

# Quantum 0b
yb = qtm0b.imag
zb = qtm0b.real
plt_qtm0b, = ax0.plot(x, yb, zb, color='red', linewidth=0.5, label='Quantum B')

# Potential
plt_potential_by_a, = ax1.plot(x, theta_potential_by_a, color='blue',
                               label='Sum of delta x **2 * abs(quantum A(x))')
plt_potential_by_b, = ax1.plot(x, theta_potential_by_b, color='red',
                               label='Sum of delta x **2 * abs(quantum B(x))')
update_potential()

# Legend
ax0.legend(loc='lower right')
ax1.legend(loc='lower right')

# Embed in Tkinter
root = tk.Tk()
root.title(title_tk)
canvas = FigureCanvasTkAgg(fig, root)
canvas.get_tk_widget().pack(expand=True, fill='both')

toolbar = NavigationToolbar2Tk(canvas, root)
canvas.get_tk_widget().pack()

# Animation
frm_anim = ttk.Labelframe(root, relief='ridge', text='Animation', labelanchor='n')
frm_anim.pack(side='left', fill=tk.Y)
btn_play = tk.Button(frm_anim, text="Play/Pause", command=switch)
btn_play.pack(side='left')
btn_step = tk.Button(frm_anim, text="Step", command=step)
btn_step.pack(side='left')
btn_reset = tk.Button(frm_anim, text="Reset", command=reset)
btn_reset.pack(side='left')

# Quantum A
frm_qtm_a = ttk.Labelframe(root, relief="ridge", text="Quantum A (initial)", labelanchor="n")
frm_qtm_a.pack(side='left', fill=tk.Y)

lbl_k0a = tk.Label(frm_qtm_a, text="k:")
lbl_k0a.pack(side='left')
var_k0a = tk.StringVar(root)
var_k0a.set(str(k0a))
spn_k0a = tk.Spinbox(
    frm_qtm_a, textvariable=var_k0a, format="%.2f", from_=-2., to=2., increment=0.01,
    command=lambda: set_k0a(var_k0a.get()), width=6
    )
spn_k0a.pack(side='left')

lbl_x0a = tk.Label(frm_qtm_a, text="x:")
lbl_x0a.pack(side='left')
var_x0a = tk.StringVar(root)
var_x0a.set(str(mu0a))
spn_x0a = tk.Spinbox(
    frm_qtm_a, textvariable=var_x0a, format="%.1f", from_=x_min0, to=x_max0, increment=100,
    command=lambda: set_x0a(var_x0a.get()), width=6
    )
spn_x0a.pack(side='left')

lbl_s0a = tk.Label(frm_qtm_a, text="sigma:")
lbl_s0a.pack(side='left')
var_s0a = tk.StringVar(root)
var_s0a.set(str(sigma0a))
spn_s0a = tk.Spinbox(
    frm_qtm_a, textvariable=var_s0a, format="%.1f", from_=-600, to=600, increment=10,
    command=lambda: set_sigma0a(var_s0a.get()), width=6
    )
spn_s0a.pack(side='left')

# Quantum B
frm_qtm_b = ttk.Labelframe(root, relief="ridge", text="Quantum B (initial)", labelanchor="n")
frm_qtm_b.pack(side='left', fill=tk.Y)

lbl_k0b = tk.Label(frm_qtm_b, text="k:")
lbl_k0b.pack(side='left')
var_k0b = tk.StringVar(root)
var_k0b.set(str(k0b))
spn_k0b = tk.Spinbox(
    frm_qtm_b, textvariable=var_k0b, format="%.2f", from_=-2., to=2., increment=0.01,
    command=lambda: set_k0b(var_k0b.get()), width=6
    )
spn_k0b.pack(side='left')

lbl_x0b = tk.Label(frm_qtm_b, text="x:")
lbl_x0b.pack(side='left')
var_x0b = tk.StringVar(root)
var_x0b.set(str(mu0b))
spn_x0b = tk.Spinbox(
    frm_qtm_b, textvariable=var_x0b, format="%.1f", from_=x_min0, to=x_max0, increment=100,
    command=lambda: set_x0b(var_x0b.get()), width=6
    )
spn_x0b.pack(side='left')

lbl_s0b = tk.Label(frm_qtm_b, text="sigma:")
lbl_s0b.pack(side='left')
var_s0b = tk.StringVar(root)
var_s0b.set(str(sigma0b))
spn_s0b = tk.Spinbox(
    frm_qtm_b, textvariable=var_s0b, format="%.1f", from_=-600, to=600, increment=10,
    command=lambda: set_sigma0b(var_s0b.get()), width=6
    )
spn_s0b.pack(side='left')

# Delta t
frm_t = ttk.Labelframe(root, relief="ridge", text="Delta t", labelanchor="n")
frm_t.pack(side='left', fill=tk.Y)

lbl_tm = tk.Label(frm_t, text="Mantissa:")
lbl_tm.pack(side='left')
var_tm = tk.StringVar(root)
var_tm.set(str(tm))
spn_tm = tk.Spinbox(
    frm_t, textvariable=var_tm, format="%.1f", from_=-100., to=100., increment=0.1,
    command=lambda: set_tm(var_tm.get()), width=6
    )
spn_tm.pack(side='left')

lbl_te = tk.Label(frm_t, text="Exponent:")
lbl_te.pack(side='left')
var_te = tk.StringVar(root)
var_te.set(str(te))
spn_te = tk.Spinbox(
    frm_t, textvariable=var_te, format="%.1f", from_=-10., to=10., increment=0.1,
    command=lambda: set_te(var_te.get()), width=6
    )
spn_te.pack(side='left')

# Potential coefficient
frm_pc = ttk.Labelframe(root, relief="ridge", text="Potential coefficient", labelanchor="n")
frm_pc.pack(side='left', fill=tk.Y)

lbl_sg = tk.Label(frm_pc, text="sign")
lbl_sg.pack(side='left')
var_sg = tk.StringVar(root)
var_sg.set(str(sign_potential))
spn_sg = tk.Spinbox(
    frm_pc, textvariable=var_sg, format="%.1f", from_=-1., to=1., increment=1,
    command=lambda: set_potential_sign(var_sg.get()), width=4
    )
spn_sg.pack(side='left')

lbl_pc = tk.Label(frm_pc, text="1/")
lbl_pc.pack(side='left')
var_pc = tk.StringVar(root)
var_pc.set(str(pc))
spn_pc = tk.Spinbox(
    frm_pc, textvariable=var_pc, format="%.1f", from_=-100000., to=100000., increment=100,
    command=lambda: set_potential_coefficient(var_pc.get()), width=6
    )
spn_pc.pack(side='left')

# main loop
anim = animation.FuncAnimation(fig, update, interval=100, save_count=100)
root.mainloop()
