"""Graphical interface for the path-length / radiation-interception model.

A small Tkinter application that collects the model inputs, lets the user
toggle optional components (canopy row geometry, scattering / absorbed
fraction, diffuse radiation), runs :mod:`pathlengthdistribution`, and then
plots the path-length distribution alongside the scalar interception
quantities.

Run it with the project's virtual environment (the interpreter that has the
model's dependencies plus Tk)::

    venv/bin/python pathlength_ui.py

The compute layer (:func:`run_model`) has no Tk dependency, so it can be
imported and exercised headlessly (see ``tests/test_ui_smoke.py``).
"""

import os
import sys
import threading
import traceback

import numpy as np

# Make the top-level model module importable regardless of CWD (mirrors
# tests/conftest.py).
REPO_ROOT = os.path.dirname(os.path.abspath(__file__))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

import pathlengthdistribution as pld  # noqa: E402


# ---------------------------------------------------------------------------
# Compute layer (no Tk dependency)
# ---------------------------------------------------------------------------

# Shapes offered in the UI, mapped to the model's shape string.  In the model
# 'sphere' is merely an alias for 'ellipsoid' (a sphere is an ellipsoid with
# equal diameters), so only 'ellipsoid' is offered to avoid a misleading
# duplicate.  The mesh option is 'polymesh' (which also accepts the 'cone'
# alias); it needs a PLY.
SHAPE_CHOICES = ('ellipsoid', 'cylinder', 'prism', 'polymesh')
MESH_SHAPES = ('polymesh', 'cone', 'mesh')


def _needs_ply(shape):
    return str(shape).strip().lower() in MESH_SHAPES


def run_model(params):
    """Run the model for a plain parameter dict and return a results dict.

    ``params`` keys (all angles in degrees):

    Geometry / beam (always used):
        shape, scale_x, scale_y, scale_z, zenith, azimuth, nrays, plyfile,
        bins, normalize.

    Optional component flags and their extra inputs:
        include_direct   -> Gtheta, LAD
        include_canopy   -> sr, sp, phi  (phi defaults to azimuth if None)
        include_scatter  -> rho_l, tau_l, rho_s, Q0  (implies canopy)
        include_diffuse  -> diffuse_level ('crown'|'canopy'), n_zenith

    Returns a dict with the raw ``path_length`` array, its ``bins``/
    ``normalize`` settings, and whichever scalar components were requested.
    Raises ``ValueError``/``Exception`` from the model on bad input; callers
    are expected to catch and present these.
    """
    shape = params['shape']
    sx = params['scale_x']
    sy = params['scale_y']
    sz = params['scale_z']
    zenith = params['zenith']
    azimuth = params['azimuth']
    nrays = params['nrays']
    plyfile = params.get('plyfile', '') or ''
    bins = params.get('bins', 10)
    normalize = params.get('normalize', True)

    if _needs_ply(shape):
        if not plyfile:
            raise ValueError(
                "A PLY file is required for the mesh (polymesh/cone) shape.")
        if not os.path.exists(plyfile):
            raise ValueError("PLY file does not exist: {}".format(plyfile))

    include_scatter = params.get('include_scatter', False)
    include_diffuse = params.get('include_diffuse', False)

    Gtheta = params.get('Gtheta', 0.5)
    LAD = params.get('LAD', 1.0)

    # Canopy row geometry is always supplied (sr, sp, phi); phi defaults to the
    # beam azimuth.  The crown- and canopy-level interception fractions are the
    # model's headline output and are always computed.
    sr = params['sr']
    sp = params['sp']
    phi = params.get('phi', None)
    if phi is None:
        phi = azimuth

    results = {
        'shape': shape,
        'bins': bins,
        'normalize': normalize,
    }

    # --- path-length distribution (always) -------------------------------
    path_length, projected_area = pld.pathlengths(
        shape, sx, sy, sz, zenith, azimuth, nrays,
        plyfile=plyfile, degrees=True)

    results['path_length'] = np.asarray(path_length, dtype=float)
    results['projected_area'] = float(projected_area)
    results['n_intersecting'] = int(results['path_length'].size)
    if results['n_intersecting'] > 0:
        results['mean_path'] = float(results['path_length'].mean())
        results['max_path'] = float(results['path_length'].max())
    else:
        results['mean_path'] = 0.0
        results['max_path'] = 0.0

    # --- crown interception components (always) ---------------------------
    results['S_theta'] = float(pld.silhouette_area(
        shape, sx, sy, sz, zenith, nrays=nrays, plyfile=plyfile,
        degrees=True))
    results['S_zero'] = float(pld.silhouette_area(
        shape, sx, sy, sz, 0.0, nrays=nrays, plyfile=plyfile,
        degrees=True))
    results['P_leaf'] = float(pld.crown_interception(
        Gtheta, LAD, shape, sx, sy, sz, zenith, azimuth, nrays,
        plyfile=plyfile, degrees=True))
    results['sunlit_fraction'] = float(pld.crown_sunlit_fraction(
        Gtheta, LAD, shape, sx, sy, sz, zenith, azimuth, nrays,
        plyfile=plyfile, degrees=True))

    # --- canopy binomial interception (always) ----------------------------
    results['P_canopy'] = float(pld.canopy_interception(
        Gtheta, LAD, shape, sx, sy, sz, zenith, azimuth, nrays, sr, sp,
        phi=phi, plyfile=plyfile, degrees=True))

    # --- scattering / absorbed fraction (three-mode) ----------------------
    if include_scatter:
        rho_l = params['rho_l']
        tau_l = params['tau_l']
        rho_s = params['rho_s']
        Q0 = params.get('Q0', 1.0)
        results['zeta'] = float(1.0 - rho_l - tau_l)
        # NOTE: absorbed_fraction takes (LAD, Gtheta, ...) -- order reversed
        # relative to every other model function.
        results['absorbed'] = float(pld.absorbed_fraction(
            LAD, Gtheta, shape, sx, sy, sz, zenith, azimuth, nrays, sr, sp,
            rho_l, tau_l, rho_s, Q0=Q0, phi=phi, plyfile=plyfile,
            degrees=True))
        results['Q0'] = float(Q0)

    # --- diffuse (hemispherically-integrated) -----------------------------
    if include_diffuse:
        level = params.get('diffuse_level', 'canopy')
        n_zenith = int(params.get('n_zenith', 18))
        kwargs = {}
        if level == 'canopy':
            kwargs['sr'] = params['sr']
            kwargs['sp'] = params['sp']
        results['P_diffuse'] = float(pld.diffuse_interception(
            Gtheta, LAD, shape, sx, sy, sz, nrays, level=level,
            n_zenith=n_zenith, plyfile=plyfile, degrees=True, **kwargs))
        results['diffuse_level'] = level

    return results


def warm_up_jit():
    """Trigger numba JIT compilation for the analytic shapes (mirrors the
    test-suite warm-up).  Cheap tiny calls; safe to run on a background thread.
    Mesh shapes compile lazily on first real use to avoid needing a PLY here.
    """
    for shape in ('ellipsoid', 'cylinder', 'prism'):
        try:
            pld.pathlengths(shape, 10.0, 10.0, 10.0, 0.1, 0.0, 16)
        except Exception:  # pragma: no cover - warm-up must never crash the UI
            pass


# ---------------------------------------------------------------------------
# Tkinter application
# ---------------------------------------------------------------------------

def _build_app():
    """Construct and return the Tk root running the application.

    Imported lazily inside ``main`` so that importing this module for the
    compute layer (or the smoke test) never touches Tk or matplotlib.
    """
    import tkinter as tk
    from tkinter import ttk, filedialog, messagebox

    import matplotlib
    matplotlib.use('TkAgg')
    from matplotlib.figure import Figure
    from matplotlib.patches import Ellipse
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

    root = tk.Tk()
    root.title('Path-length distribution & radiation interception')
    root.minsize(760, 480)

    # ---- variables ------------------------------------------------------
    v = {
        # Default to a valid (non-overlapping) orchard-like crown: the lateral
        # extents (3 x 5) fit within the default spacing (sp=4, sr=6).
        'shape': tk.StringVar(value='ellipsoid'),
        'scale_x': tk.StringVar(value='3.0'),
        'scale_y': tk.StringVar(value='5.0'),
        'scale_z': tk.StringVar(value='6.0'),
        'nrays': tk.StringVar(value='5000'),
        'plyfile': tk.StringVar(value=os.path.join(REPO_ROOT, 'PLY', 'sphere.ply')),
        'zenith': tk.StringVar(value='30.0'),
        'azimuth': tk.StringVar(value='0.0'),
        'bins': tk.StringVar(value='15'),
        'normalize': tk.BooleanVar(value=True),
        'Gtheta': tk.StringVar(value='0.5'),
        'LAD': tk.StringVar(value='1.0'),
        'sr': tk.StringVar(value='6.0'),
        'sp': tk.StringVar(value='4.0'),
        'phi': tk.StringVar(value=''),
        'rho_l': tk.StringVar(value='0.09'),
        'tau_l': tk.StringVar(value='0.04'),
        'rho_s': tk.StringVar(value='0.18'),
        'Q0': tk.StringVar(value='1.0'),
        'diffuse_level': tk.StringVar(value='canopy'),
        'n_zenith': tk.StringVar(value='18'),
        'include_scatter': tk.BooleanVar(value=False),
        'include_diffuse': tk.BooleanVar(value=False),
    }

    main = ttk.Frame(root, padding=10)
    main.grid(row=0, column=0, sticky='nsew')
    root.columnconfigure(0, weight=1)
    root.rowconfigure(0, weight=1)

    # The input column can grow taller than the window once several optional
    # sections are shown, so it lives inside a scrollable canvas.  ``inputs`` is
    # the frame the sections are gridded into; ``input_canvas`` scrolls it.
    input_canvas = tk.Canvas(main, highlightthickness=0, width=250, height=560)
    input_canvas.grid(row=0, column=0, sticky='ns')
    input_scroll = ttk.Scrollbar(main, orient='vertical',
                                 command=input_canvas.yview)
    input_scroll.grid(row=0, column=1, sticky='ns', padx=(0, 12))
    input_canvas.configure(yscrollcommand=input_scroll.set)

    inputs = ttk.Frame(input_canvas)
    _inputs_window = input_canvas.create_window((0, 0), window=inputs,
                                                anchor='nw')

    # Keep the scroll region matching the inner frame's true height, and make
    # the inner frame track the canvas width so sections fill it.
    inputs.bind('<Configure>',
                lambda e: input_canvas.configure(
                    scrollregion=input_canvas.bbox('all')))
    input_canvas.bind(
        '<Configure>',
        lambda e: input_canvas.itemconfigure(_inputs_window, width=e.width))

    def _on_mousewheel(event):
        # macOS/Windows deliver small delta values; scroll by whole units.
        input_canvas.yview_scroll(int(-1 * (event.delta or 0)), 'units')

    input_canvas.bind_all('<MouseWheel>', _on_mousewheel)

    outputs = ttk.Frame(main)
    outputs.grid(row=0, column=2, sticky='nsew')
    main.columnconfigure(2, weight=1)
    main.rowconfigure(0, weight=1)

    def add_row(parent, label, var, r, width=12):
        ttk.Label(parent, text=label).grid(row=r, column=0, sticky='w', pady=1)
        e = ttk.Entry(parent, textvariable=var, width=width)
        e.grid(row=r, column=1, sticky='w', pady=1)
        return e

    # ---- geometry section ----------------------------------------------
    geo = ttk.LabelFrame(inputs, text='Crown Geometry', padding=8)
    geo.grid(row=0, column=0, sticky='ew', pady=(0, 6))
    ttk.Label(geo, text='crown shape').grid(row=0, column=0, sticky='w')
    shape_cb = ttk.Combobox(geo, textvariable=v['shape'], values=SHAPE_CHOICES,
                            width=10, state='readonly')
    shape_cb.grid(row=0, column=1, sticky='w')

    # The three scale entries mean different things per shape, so their labels
    # are set by update_scale_labels() below rather than fixed text.
    scale_label_x = ttk.Label(geo)
    scale_label_x.grid(row=1, column=0, sticky='w', pady=1)
    ttk.Entry(geo, textvariable=v['scale_x'], width=12).grid(
        row=1, column=1, sticky='w', pady=1)
    scale_label_y = ttk.Label(geo)
    scale_label_y.grid(row=2, column=0, sticky='w', pady=1)
    ttk.Entry(geo, textvariable=v['scale_y'], width=12).grid(
        row=2, column=1, sticky='w', pady=1)
    scale_label_z = ttk.Label(geo)
    scale_label_z.grid(row=3, column=0, sticky='w', pady=1)
    ttk.Entry(geo, textvariable=v['scale_z'], width=12).grid(
        row=3, column=1, sticky='w', pady=1)
    add_row(geo, 'nrays', v['nrays'], 4)

    # Per-shape names for (scale_x, scale_y, scale_z).
    _SCALE_LABELS = {
        'ellipsoid': ('x-diameter (m)', 'y-diameter (m)', 'height (m)'),
        'sphere': ('x-diameter (m)', 'y-diameter (m)', 'height (m)'),
        'cylinder': ('x-diameter (m)', 'y-diameter (m)', 'height (m)'),
        'prism': ('x-width (m)', 'y-width (m)', 'height (m)'),
        'polymesh': ('x-scale (×)', 'y-scale (×)', 'z-scale (×)'),
        'cone': ('x-scale (×)', 'y-scale (×)', 'z-scale (×)'),
    }

    def update_scale_labels(*_):
        key = v['shape'].get().strip().lower()
        lx, ly, lz = _SCALE_LABELS.get(
            key, ('scale_x (m)', 'scale_y (m)', 'scale_z (m)'))
        scale_label_x.configure(text=lx)
        scale_label_y.configure(text=ly)
        scale_label_z.configure(text=lz)

    v['shape'].trace_add('write', update_scale_labels)

    ply_label = ttk.Label(geo, text='PLY file')
    ply_entry = ttk.Entry(geo, textvariable=v['plyfile'], width=24)
    ply_button = ttk.Button(geo, text='Browse…')

    def browse_ply():
        path = filedialog.askopenfilename(
            title='Select a PLY mesh',
            initialdir=os.path.join(REPO_ROOT, 'PLY'),
            filetypes=[('PLY mesh', '*.ply'), ('All files', '*.*')])
        if path:
            v['plyfile'].set(path)

    ply_button.configure(command=browse_ply)

    def update_ply_visibility(*_):
        if _needs_ply(v['shape'].get()):
            ply_label.grid(row=5, column=0, sticky='w', pady=(4, 0))
            ply_entry.grid(row=6, column=0, columnspan=2, sticky='w')
            ply_button.grid(row=7, column=0, sticky='w', pady=(2, 0))
        else:
            ply_label.grid_remove()
            ply_entry.grid_remove()
            ply_button.grid_remove()

    v['shape'].trace_add('write', update_ply_visibility)

    # ---- beam section ---------------------------------------------------
    beam = ttk.LabelFrame(inputs, text='Sun direction', padding=8)
    beam.grid(row=1, column=0, sticky='ew', pady=6)
    add_row(beam, 'zenith (deg)', v['zenith'], 0)
    add_row(beam, 'azimuth (deg)', v['azimuth'], 1)

    # ---- distribution options ------------------------------------------
    dist = ttk.LabelFrame(inputs, text='Distribution discretization', padding=8)
    dist.grid(row=2, column=0, sticky='ew', pady=6)
    add_row(dist, 'bins', v['bins'], 0)
    ttk.Checkbutton(dist, text='normalize (density)',
                    variable=v['normalize']).grid(row=1, column=0,
                                                  columnspan=2, sticky='w')

    # ---- leaf params ----------------------------------------------------
    leaf = ttk.LabelFrame(inputs, text='Leaf parameters', padding=8)
    leaf.grid(row=3, column=0, sticky='ew', pady=6)
    add_row(leaf, 'Gtheta (Ross G)', v['Gtheta'], 0)
    add_row(leaf, 'LAD (m2/m3)', v['LAD'], 1)

    # ---- canopy geometry (always shown) --------------------------------
    # Row/plant spacing are core model inputs; the canopy interception
    # fraction is computed on every run, so these are not gated behind a
    # checkbox.
    canopy = ttk.LabelFrame(inputs, text='Canopy configuration', padding=8)
    canopy.grid(row=4, column=0, sticky='ew', pady=6)
    add_row(canopy, 'sr row spacing (m)', v['sr'], 0)
    add_row(canopy, 'sp plant spacing (m)', v['sp'], 1)
    add_row(canopy, 'phi row orient (deg)', v['phi'], 2)
    ttk.Label(canopy, text='(phi blank = azimuth)').grid(
        row=3, column=0, columnspan=2, sticky='w')

    # ---- optional component checkboxes ---------------------------------
    comps = ttk.LabelFrame(inputs, text='Optional components', padding=8)
    comps.grid(row=5, column=0, sticky='ew', pady=6)

    ttk.Checkbutton(comps, text='Scattering / absorbed fraction',
                    variable=v['include_scatter']).grid(row=3, column=0,
                                                        sticky='w')
    scatter_fields = ttk.Frame(comps)
    add_row(scatter_fields, 'rho_l leaf refl', v['rho_l'], 0)
    add_row(scatter_fields, 'tau_l leaf trans', v['tau_l'], 1)
    add_row(scatter_fields, 'rho_s soil refl', v['rho_s'], 2)
    add_row(scatter_fields, 'Q0 incident flux', v['Q0'], 3)
    zeta_label = ttk.Label(scatter_fields, text='zeta = 1 - rho_l - tau_l')
    zeta_label.grid(row=4, column=0, columnspan=2, sticky='w')

    def update_zeta(*_):
        try:
            z = 1.0 - float(v['rho_l'].get()) - float(v['tau_l'].get())
            zeta_label.configure(text='zeta = {:.3f}'.format(z))
        except ValueError:
            zeta_label.configure(text='zeta = 1 - rho_l - tau_l')

    v['rho_l'].trace_add('write', update_zeta)
    v['tau_l'].trace_add('write', update_zeta)

    ttk.Checkbutton(comps, text='Diffuse radiation',
                    variable=v['include_diffuse']).grid(row=5, column=0,
                                                        sticky='w')
    diffuse_fields = ttk.Frame(comps)
    ttk.Label(diffuse_fields, text='level').grid(row=0, column=0, sticky='w')
    ttk.Combobox(diffuse_fields, textvariable=v['diffuse_level'],
                 values=('canopy', 'crown'), width=8,
                 state='readonly').grid(row=0, column=1, sticky='w')
    add_row(diffuse_fields, 'n_zenith', v['n_zenith'], 1)

    # Reveal/hide the optional field groups when their checkbox is toggled.
    def update_component_visibility(*_):
        if v['include_scatter'].get():
            scatter_fields.grid(row=4, column=0, sticky='w', padx=(18, 0))
        else:
            scatter_fields.grid_remove()

        if v['include_diffuse'].get():
            diffuse_fields.grid(row=6, column=0, sticky='w', padx=(18, 0))
        else:
            diffuse_fields.grid_remove()

    for key in ('include_scatter', 'include_diffuse'):
        v[key].trace_add('write', update_component_visibility)

    # ---- run button + status -------------------------------------------
    run_btn = ttk.Button(inputs, text='Run model')
    run_btn.grid(row=6, column=0, sticky='ew', pady=(8, 2))
    status = ttk.Label(inputs, text='Ready.')
    status.grid(row=7, column=0, sticky='w')

    # ---- outputs: headline banner + plot + scalar readout ---------------
    # The headline result is the fraction of incident light the canopy
    # intercepts; give it a large, plain-language banner so it does not get
    # lost in the detailed readout below.
    banner = tk.Frame(outputs, bg='#e8f0fe', bd=1, relief='solid')
    banner.grid(row=0, column=0, sticky='ew', pady=(0, 8))
    outputs.columnconfigure(0, weight=1)

    headline_var = tk.StringVar(value='—')
    headline_caption = tk.StringVar(
        value='Fraction of incident light intercepted by the canopy')
    headline_num = tk.Label(banner, textvariable=headline_var, bg='#e8f0fe',
                            fg='#174ea6', font=('Helvetica', 30, 'bold'))
    headline_num.pack(pady=(8, 0))
    tk.Label(banner, textvariable=headline_caption, bg='#e8f0fe',
             fg='#3c4043', font=('Helvetica', 11)).pack(pady=(0, 8))

    # ---- canopy configuration schematic (live top-down preview) ---------
    # A top-down view of a few crowns illustrating crown size, spacing, and row
    # orientation, with a compass.  Updated live from the input fields (not
    # only on Run) so users can see the layout as they type.
    config_fig = Figure(figsize=(5.2, 2.6), dpi=100)
    config_ax = config_fig.add_subplot(111)
    config_canvas = FigureCanvasTkAgg(config_fig, master=outputs)
    config_canvas.get_tk_widget().grid(row=1, column=0, sticky='nsew')

    fig = Figure(figsize=(5.2, 3.0), dpi=100)
    ax = fig.add_subplot(111)
    ax.set_xlabel('path length (m)')
    ax.set_ylabel('probability density')
    canvas = FigureCanvasTkAgg(fig, master=outputs)
    canvas.get_tk_widget().grid(row=2, column=0, sticky='nsew')

    ttk.Label(outputs, text='Detailed results').grid(
        row=3, column=0, sticky='w', pady=(8, 0))
    readout = tk.Text(outputs, width=44, height=13, wrap='word')
    readout.grid(row=4, column=0, sticky='nsew', pady=(2, 0))
    readout.configure(state='disabled')
    # The two plots take the extra height; the readout grows a little too.
    outputs.rowconfigure(1, weight=2)   # config schematic
    outputs.rowconfigure(2, weight=3)   # distribution plot
    outputs.rowconfigure(4, weight=1)   # readout

    # ---- parameter gathering + validation ------------------------------
    def _f(name, label):
        raw = v[name].get().strip()
        try:
            return float(raw)
        except ValueError:
            raise ValueError("'{}' must be a number (got '{}').".format(
                label, raw))

    def gather_params():
        params = {
            'shape': v['shape'].get(),
            'scale_x': _f('scale_x', 'scale_x'),
            'scale_y': _f('scale_y', 'scale_y'),
            'scale_z': _f('scale_z', 'scale_z'),
            'zenith': _f('zenith', 'zenith'),
            'azimuth': _f('azimuth', 'azimuth'),
            'nrays': int(_f('nrays', 'nrays')),
            'plyfile': v['plyfile'].get().strip(),
            'bins': int(_f('bins', 'bins')),
            'normalize': bool(v['normalize'].get()),
            'Gtheta': _f('Gtheta', 'Gtheta'),
            'LAD': _f('LAD', 'LAD'),
            # Canopy geometry is always supplied.
            'sr': _f('sr', 'sr'),
            'sp': _f('sp', 'sp'),
            'phi': (float(v['phi'].get().strip())
                    if v['phi'].get().strip() else None),
            'include_scatter': bool(v['include_scatter'].get()),
            'include_diffuse': bool(v['include_diffuse'].get()),
        }
        if params['include_scatter']:
            params['rho_l'] = _f('rho_l', 'rho_l')
            params['tau_l'] = _f('tau_l', 'tau_l')
            params['rho_s'] = _f('rho_s', 'rho_s')
            params['Q0'] = _f('Q0', 'Q0')
        if params['include_diffuse']:
            params['diffuse_level'] = v['diffuse_level'].get()
            params['n_zenith'] = int(_f('n_zenith', 'n_zenith'))
        return params

    # ---- rendering ------------------------------------------------------
    def render(results):
        ax.clear()
        pl = results['path_length']
        if pl.size:
            ax.hist(pl, bins=results['bins'], density=results['normalize'],
                    color='#4c72b0', edgecolor='white', alpha=0.85)
            ax.axvline(results['mean_path'], color='#c44e52', linestyle='--',
                       label='mean {:.2f} m'.format(results['mean_path']))
            ax.axvline(results['max_path'], color='#55a868', linestyle=':',
                       label='max {:.2f} m'.format(results['max_path']))
            ax.legend(fontsize=8)
        else:
            ax.text(0.5, 0.5, 'No rays intersected the shape',
                    ha='center', va='center', transform=ax.transAxes)
        ax.set_xlabel('path length (m)')
        ax.set_ylabel('probability density' if results['normalize'] else 'count')
        ax.set_title('Path-length distribution ({})'.format(results['shape']))
        fig.tight_layout()
        canvas.draw()

        # Headline banner: canopy interception is the primary result.  When
        # scattering is enabled the absorbed fraction is the more complete
        # answer, so promote that instead.
        if 'absorbed' in results:
            headline_var.set('{:.0%}'.format(results['absorbed']))
            headline_caption.set(
                'Fraction of incident light absorbed by the canopy')
        else:
            headline_var.set('{:.0%}'.format(results['P_canopy']))
            headline_caption.set(
                'Fraction of incident light intercepted by the canopy')

        lines = [
            'Interception & absorption',
            '  Light intercepted by canopy : {:.1%}'.format(results['P_canopy']),
            '  Intercepted per crown       : {:.1%}'.format(results['P_leaf']),
            '  Leaf area sunlit            : {:.1%}'.format(
                results['sunlit_fraction']),
        ]
        if 'absorbed' in results:
            lines += [
                '  Light absorbed by canopy    : {:.1%}'.format(
                    results['absorbed']),
                '  Leaf absorptivity (zeta)    : {:.3f}'.format(results['zeta']),
            ]
        if 'P_diffuse' in results:
            lines.append(
                '  Diffuse ({:>6}) intercepted : {:.1%}'.format(
                    results['diffuse_level'], results['P_diffuse']))
        lines += [
            '',
            'Path lengths',
            '  Rays intersecting : {}'.format(results['n_intersecting']),
            '  Mean path length  : {:.3f} m'.format(results['mean_path']),
            '  Max path length   : {:.3f} m'.format(results['max_path']),
            '  Projected area    : {:.3f} m^2'.format(results['projected_area']),
            '',
            'Silhouette area',
            '  S(theta)          : {:.3f} m^2'.format(results['S_theta']),
            '  S(0)              : {:.3f} m^2'.format(results['S_zero']),
        ]

        readout.configure(state='normal')
        readout.delete('1.0', 'end')
        readout.insert('1.0', '\n'.join(lines))
        readout.configure(state='disabled')

    # ---- live canopy-configuration schematic ---------------------------
    def _try_float(name, default):
        try:
            return float(v[name].get().strip())
        except (ValueError, AttributeError):
            return default

    def draw_config(*_):
        """Top-down schematic of crown size, spacing, and row orientation.

        Drawn in an absolute compass frame (North up, East right; azimuths
        measured clockwise from North).  The sun arrow points from the sun's
        azimuth toward the canopy.  ``phi`` is the sun azimuth *relative to the
        rows* (per the model), so the rows are drawn at compass bearing
        ``azimuth - phi``.  Along a row, plants are spaced by ``sp`` and crowns
        have along-row diameter ``scale_x``; rows are spaced ``sr`` apart and
        crowns have across-row diameter ``scale_y``.  Crowns may not overlap, so
        ``scale_x`` must not exceed ``sp`` nor ``scale_y`` exceed ``sr``; a
        violation is flagged in red.
        """
        import numpy as _np

        sx = _try_float('scale_x', 0.0)
        sy = _try_float('scale_y', 0.0)
        sp = _try_float('sp', 0.0)
        sr = _try_float('sr', 0.0)
        azimuth = _try_float('azimuth', 0.0)
        # phi blank -> use the sun azimuth (same default the model uses).
        phi_raw = v['phi'].get().strip()
        phi = azimuth if phi_raw == '' else _try_float('phi', 0.0)

        config_ax.clear()

        if sp <= 0 or sr <= 0:
            config_ax.text(0.5, 0.5, 'Enter positive sp and sr\nto preview layout',
                           ha='center', va='center',
                           transform=config_ax.transAxes, fontsize=9)
            config_ax.set_axis_off()
            config_fig.tight_layout()
            config_canvas.draw()
            return

        # Compass bearing (deg clockwise from North) of the row direction.
        row_bearing = azimuth - phi
        beta = _np.radians(row_bearing)
        # Unit vectors: along the row, and across the row (rotate +90deg).  In
        # screen coords N is +y, E is +x, so a bearing b -> (sin b, cos b).
        along = _np.array([_np.sin(beta), _np.cos(beta)])          # sp / scale_x
        across = _np.array([_np.sin(beta + _np.pi / 2),
                            _np.cos(beta + _np.pi / 2)])            # sr / scale_y

        ncols, nrows = 3, 3
        overlap = (sx > sp + 1e-9) or (sy > sr + 1e-9)
        edge = '#c0392b' if overlap else '#2e7d32'
        face = '#f2b8b0' if overlap else '#b7dfb0'

        # Crown ellipse: matplotlib Ellipse angle is CCW from +x for its width
        # axis.  We want the width axis (scale_x) along the row bearing.
        ell_angle = 90.0 - row_bearing  # convert bearing (CW from N) to CCW-from-x
        centers = []
        for i in range(ncols):
            for j in range(nrows):
                c = ((i - (ncols - 1) / 2.0) * sp) * along \
                    + ((j - (nrows - 1) / 2.0) * sr) * across
                centers.append(c)
                config_ax.add_patch(Ellipse(
                    (c[0], c[1]), width=sx, height=sy, angle=ell_angle,
                    facecolor=face, edgecolor=edge, linewidth=1.4,
                    alpha=0.85, zorder=2))

        centers = _np.array(centers)

        # Spacing annotations along the two lattice directions, offset to a side.
        # sp arrow between two in-row neighbours, shifted one row across.
        base_sp = -1.9 * across * sr
        config_ax.annotate('', xy=tuple(base_sp + 0.5 * sp * along),
                           xytext=tuple(base_sp - 0.5 * sp * along),
                           arrowprops=dict(arrowstyle='<->', color='#555'))
        # Offset the label further out along the same perpendicular so it does
        # not sit on the arrow or the crowns.
        config_ax.text(*(base_sp - 0.35 * across * sr), 'sp={:g} m'.format(sp),
                       ha='center', va='center', fontsize=8, color='#555')
        # sr arrow between two rows, shifted one plant along.
        base_sr = -1.9 * along * sp
        config_ax.annotate('', xy=tuple(base_sr + 0.5 * sr * across),
                           xytext=tuple(base_sr - 0.5 * sr * across),
                           arrowprops=dict(arrowstyle='<->', color='#555'))
        config_ax.text(*(base_sr - 0.35 * along * sp), 'sr={:g} m'.format(sr),
                       ha='center', va='center', fontsize=8, color='#555')

        # View extent covers the crown block plus the annotation offsets.
        span = max(centers[:, 0].max() - centers[:, 0].min(),
                   centers[:, 1].max() - centers[:, 1].min())
        ext = 0.5 * span + max(sx, sy, sp, sr) * 2.4
        config_ax.set_xlim(-ext, ext)
        config_ax.set_ylim(-ext, ext)
        config_ax.set_aspect('equal')
        config_ax.set_xticks([])
        config_ax.set_yticks([])

        # Sun-direction arrow: azimuth is the compass bearing the sunlight comes
        # FROM; draw the beam pointing toward the canopy centre.
        sun_dir = _np.array([_np.sin(_np.radians(azimuth)),
                             _np.cos(_np.radians(azimuth))])
        L = 0.9 * ext
        config_ax.annotate(
            '', xy=tuple(0.35 * L * sun_dir), xytext=tuple(L * sun_dir),
            arrowprops=dict(arrowstyle='-|>', color='#e8a600', lw=2.2))
        config_ax.text(*(L * sun_dir), ' sun (az={:g}°, φ={:g}°)'.format(
            azimuth, phi), fontsize=8, color='#b57e00', va='center',
            ha='left' if sun_dir[0] >= 0 else 'right')

        # Compass rose, bottom-right corner (N up, E right).
        cxk, cyk = 0.80 * ext, -0.80 * ext
        r = 0.13 * ext
        config_ax.annotate('', xy=(cxk, cyk + r), xytext=(cxk, cyk - r),
                           arrowprops=dict(arrowstyle='-|>', color='#333', lw=1.2))
        config_ax.annotate('', xy=(cxk + r, cyk), xytext=(cxk - r, cyk),
                           arrowprops=dict(arrowstyle='-', color='#ccc', lw=0.8))
        for lbl, dx, dy in (('N', 0, 1), ('S', 0, -1), ('E', 1, 0), ('W', -1, 0)):
            config_ax.text(cxk + dx * r * 1.3, cyk + dy * r * 1.3, lbl,
                           ha='center', va='center', fontsize=7, color='#333')

        title = 'Canopy layout (top-down)'
        if overlap:
            title += '  —  crowns overlap!'
        config_ax.set_title(title, fontsize=10,
                            color='#c0392b' if overlap else '#000')
        if overlap:
            config_ax.text(
                0.5, 0.01,
                'Crown lateral size exceeds spacing: reduce scale_x≤sp and/or '
                'scale_y≤sr',
                ha='center', va='bottom', transform=config_ax.transAxes,
                fontsize=8, color='#c0392b')

        config_fig.tight_layout()
        config_canvas.draw()

    # Redraw the schematic live whenever a relevant input changes.
    for _key in ('scale_x', 'scale_y', 'sr', 'sp', 'phi', 'azimuth'):
        v[_key].trace_add('write', draw_config)

    # ---- threaded run ---------------------------------------------------
    def on_run():
        try:
            params = gather_params()
        except ValueError as exc:
            messagebox.showerror('Invalid input', str(exc))
            return

        run_btn.configure(state='disabled')
        status.configure(text='Running…')

        def worker():
            try:
                results = run_model(params)
            except Exception as exc:  # surface model errors to the UI thread
                tb = traceback.format_exc()
                root.after(0, lambda: on_error(exc, tb))
                return
            root.after(0, lambda: on_done(results))

        threading.Thread(target=worker, daemon=True).start()

    def on_done(results):
        render(results)
        status.configure(text='Done.')
        run_btn.configure(state='normal')

    def on_error(exc, tb):
        status.configure(text='Error.')
        run_btn.configure(state='normal')
        messagebox.showerror('Model error', str(exc))
        sys.stderr.write(tb)

    run_btn.configure(command=on_run)

    # initial visibility state
    update_ply_visibility()
    update_scale_labels()
    update_component_visibility()
    update_zeta()
    draw_config()

    return root


def _bring_to_front(root):
    """Force the window to draw and come to the foreground.

    On macOS the Tk window otherwise opens blank (the initial draw isn't
    flushed until an event arrives) and behind other applications.  Sizing the
    window from its requested geometry, forcing an idle update, and toggling
    ``-topmost`` reliably works around both.
    """
    root.update_idletasks()
    # Open a little larger than the packed size so the plot/readout have room
    # and it is obvious the window can be dragged bigger.  The window stays
    # fully resizable; grid weights let the plot grow to fill any extra space.
    w = max(root.winfo_reqwidth() + 120, 900)
    h = max(root.winfo_reqheight() + 60, 680)
    root.geometry('{}x{}'.format(w, h))
    # Explicitly keep the window user-resizable.  (Do NOT set a custom macOS
    # MacWindowStyle mask here -- on some Tk builds that drops the resize
    # affordance so no resize cursor appears on the window edges.)
    root.resizable(True, True)
    root.deiconify()
    root.lift()
    root.attributes('-topmost', True)
    root.after(200, lambda: root.attributes('-topmost', False))
    root.focus_force()


def main():
    root = _build_app()
    # Warm up the numba kernels off the UI thread so the first Run is snappy.
    threading.Thread(target=warm_up_jit, daemon=True).start()
    _bring_to_front(root)
    root.mainloop()


if __name__ == '__main__':
    main()
