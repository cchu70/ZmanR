"""
02_simulator.py  –  Zman-seq interactive trajectory simulator (v5)

Two optional trajectories:
  - Trajectory 1: A → B → C  (purple gradient, always shown)
  - Trajectory 2: X → Y → Z  (blue gradient, toggle on/off)

When both are active, each cell's timeline row is split:
  traj 1 rectangles above the circulation line, traj 2 below.
"""

import dash
from dash import dcc, html, Input, Output
from collections import defaultdict
import plotly.graph_objects as go
import numpy as np

# ── defaults ──────────────────────────────────────────────────────────────────

DEFAULT_N            = 20
DEFAULT_LABEL_TIMES  = [-36, -24, -12]
DEFAULT_LABEL_COLORS = ["#115e59", "#0d9488", "#5eead4"]  # dark → mid → light teal
STATE_COLORS         = {"A": "#c4b5fd", "B": "#7c3aed", "C": "#2e1065"}   # purple gradient
STATE_COLORS_2       = {"X": "#fed7aa", "Y": "#f97316", "Z": "#7c2d12"}   # orange gradient
UNLABELED_COLOR      = "#aaaaaa"
FILL_ALPHA           = 0.70

# ── colour helpers ─────────────────────────────────────────────────────────────

def hex_to_rgba(hex_color: str, alpha: float = FILL_ALPHA) -> str:
    h = hex_color.lstrip("#")
    r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
    return f"rgba({r},{g},{b},{alpha})"


def text_color_for(hex_color: str) -> str:
    h = hex_color.lstrip("#")
    r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
    return "white" if 0.299*r + 0.587*g + 0.114*b < 140 else "#333"


# ── simulation ─────────────────────────────────────────────────────────────────

def simulate_cells(n: int, min_circ: float = 12.0, max_circ: float = 24.0, seed: int = 42):
    rng = np.random.default_rng(seed)
    hi  = max(max_circ, min_circ + 1)
    cells = []
    for _ in range(n):
        te     = float(rng.uniform(-48, -1))
        offset = float(rng.uniform(min_circ, hi))
        ce     = max(te - offset, -48.0)
        cells.append((ce, te))
    cells.sort(key=lambda x: x[1])
    return cells


def get_label(ce, te, label_times):
    valid = [t for t in label_times if ce <= t <= te]
    return max(valid) if valid else None


def get_final_state(te, t_b, t_c_abs):
    h = -te
    if h < t_b:      return "A"
    if h < t_c_abs:  return "B"
    return "C"


def get_final_state_2(te, t_xy, t_yz_abs):
    h = -te
    if h < t_xy:      return "X"
    if h < t_yz_abs:  return "Y"
    return "Z"


def draw_rect_traces(cells, labels, cell_color_fn, state_names, transitions,
                     y_lo_off, y_hi_off, state_color_map=None, y_positions=None):
    """
    Build batched fill-rect traces (one dict per unique color) for one trajectory.
    Returns (rect_segs, state_annots).
      y_lo_off / y_hi_off : offsets from cell's integer y position.
      state_color_map     : if provided, each state's rect uses state_color_map[sname]
                            instead of the label-driven cell_color_fn(label).
      y_positions         : explicit y-center per cell; defaults to 1..N.
    """
    rect_segs    = defaultdict(lambda: {"x": [], "y": []})
    state_annots = []
    for i, ((ce, te), label) in enumerate(zip(cells, labels)):
        y    = y_positions[i] if y_positions is not None else i + 1
        y_lo = y + y_lo_off
        y_hi = y + y_hi_off
        for j, sname in enumerate(state_names):
            col = state_color_map[sname] if state_color_map else cell_color_fn(label)
            txt = text_color_for(col)
            a_s = te + transitions[j]
            a_e = (te + transitions[j + 1]) if j + 1 < len(transitions) else 0.0
            a_s = max(a_s, te)
            a_e = min(a_e, 0.0)
            if a_s >= a_e:
                continue
            rect_segs[col]["x"] += [a_s, a_e, a_e, a_s, a_s, None]
            rect_segs[col]["y"] += [y_lo, y_lo, y_hi, y_hi, y_lo, None]
            if a_e - a_s > 0.5:
                state_annots.append(dict(
                    x=(a_s + a_e) / 2, y=(y_lo + y_hi) / 2,
                    text=f"<b>{sname}</b>", showarrow=False,
                    font=dict(size=8 if a_e - a_s < 3 else 9, color=txt),
                    xanchor="center", yanchor="middle",
                ))
    return rect_segs, state_annots


# ── app layout ─────────────────────────────────────────────────────────────────

app = dash.Dash(__name__)

_ls = {"fontFamily": "sans-serif", "fontSize": "13px"}
_sm = {i: str(i) for i in range(0, 49, 6)}
_sg = {"height": "220px"}

app.layout = html.Div([
    html.H2("Zman-seq Trajectory Simulator",
            style={"textAlign": "center", "fontFamily": "sans-serif", "marginBottom": "8px"}),

    # ── control panel ──────────────────────────────────────────────────────
    html.Div([

        # trajectory sliders (both trajectories)
        html.Div([
            html.B("Trajectory 1  (A → B → C)", style={**_ls, "color": "#7c3aed"}),
            html.Div([
                html.Label("A → B:", style={**_ls, "width": "90px", "display": "inline-block"}),
                dcc.Slider(id="t-b", min=1, max=47, step=1, value=12, marks=_sm,
                           tooltip={"placement": "bottom", "always_visible": True}),
            ], style={"marginTop": "8px"}),
            html.Div([
                html.Label("B duration:", style={**_ls, "width": "90px", "display": "inline-block"}),
                dcc.Slider(id="t-c", min=1, max=48, step=1, value=12, marks=_sm,
                           tooltip={"placement": "bottom", "always_visible": True}),
            ], style={"marginTop": "12px"}),

            html.Hr(style={"margin": "14px 0 10px", "borderColor": "#ddd"}),

            html.B("Trajectory 2  (X → Y → Z)", style={**_ls, "color": "#3b82f6"}),
            html.Div([
                html.Label("X → Y:", style={**_ls, "width": "90px", "display": "inline-block"}),
                dcc.Slider(id="t-xy", min=1, max=47, step=1, value=6, marks=_sm,
                           tooltip={"placement": "bottom", "always_visible": True}),
            ], style={"marginTop": "8px"}),
            html.Div([
                html.Label("Y duration:", style={**_ls, "width": "90px", "display": "inline-block"}),
                dcc.Slider(id="t-yz", min=1, max=48, step=1, value=18, marks=_sm,
                           tooltip={"placement": "bottom", "always_visible": True}),
            ], style={"marginTop": "12px"}),
        ], style={"flex": "2", "paddingRight": "28px"}),

        # n cells + circulation sliders
        html.Div([
            html.B("Number of cells", style=_ls),
            dcc.Slider(id="n-cells", min=10, max=200, step=10, value=DEFAULT_N,
                       marks={i: str(i) for i in range(0, 201, 50)},
                       tooltip={"placement": "bottom", "always_visible": True}),
            html.B("Min. circulation time (hr)", style={**_ls, "display": "block", "marginTop": "14px"}),
            dcc.Slider(id="min-circ", min=1, max=47, step=1, value=12,
                       marks={i: str(i) for i in range(0, 48, 6)},
                       tooltip={"placement": "bottom", "always_visible": True}),
            html.B("Max. circulation time (hr)", style={**_ls, "display": "block", "marginTop": "14px"}),
            dcc.Slider(id="max-circ", min=2, max=48, step=1, value=24,
                       marks={i: str(i) for i in range(0, 49, 6)},
                       tooltip={"placement": "bottom", "always_visible": True}),
        ], style={"flex": "1", "paddingRight": "28px"}),

        # display options
        html.Div([
            html.B("Display options", style=_ls),
            dcc.Checklist(
                id="show-traj2",
                options=[{"label": " Show trajectory 2 (X→Y→Z)", "value": "show"}],
                value=[],
                style={**_ls, "marginTop": "10px"},
            ),
            dcc.Checklist(
                id="show-unlabeled",
                options=[{"label": " Include unlabeled cells", "value": "show"}],
                value=["show"],
                style={**_ls, "marginTop": "8px"},
            ),
            html.B("CDF line style", style={**_ls, "display": "block", "marginTop": "14px"}),
            dcc.RadioItems(
                id="cdf-style",
                options=[
                    {"label": " Direct (diagonal)", "value": "linear"},
                    {"label": " Step (hv)",          "value": "hv"},
                ],
                value="linear",
                style={**_ls, "marginTop": "6px"},
                labelStyle={"display": "block"},
            ),
        ], style={"flex": "0 0 220px", "paddingRight": "20px"}),

        # label time points + colors
        html.Div([
            html.B("Labeling time points", style=_ls),
            *[
                html.Div([
                    html.Label(f"Label {i+1}:", style={**_ls, "width": "55px", "display": "inline-block"}),
                    dcc.Input(id=f"lt{i+1}", type="number", value=DEFAULT_LABEL_TIMES[i],
                              style={"width": "62px", "marginRight": "6px"}),
                    html.Label("hr  color:", style={**_ls, "marginRight": "4px"}),
                    dcc.Input(id=f"lc{i+1}", type="text", value=DEFAULT_LABEL_COLORS[i],
                              debounce=True, placeholder="#rrggbb",
                              style={"width": "80px", "fontFamily": "monospace", "fontSize": "12px"}),
                ], style={"marginTop": "6px"})
                for i in range(3)
            ],
        ], style={"flex": "1.2"}),

    ], style={
        "display": "flex", "alignItems": "flex-start",
        "padding": "14px 18px", "border": "1px solid #ddd",
        "borderRadius": "8px", "background": "#fafafa",
        "marginBottom": "12px",
    }),

    # ── plot area ──────────────────────────────────────────────────────────
    html.Div([
        html.Div(dcc.Graph(id="main-plot"), style={"flex": "0 0 55%", "minWidth": 0}),
        html.Div([
            html.Div([
                html.Div(dcc.Graph(id="cnt-bin-plot",   style=_sg), style={"flex": "1"}),
                html.Div(dcc.Graph(id="cnt-state-plot", style=_sg), style={"flex": "1"}),
            ], style={"display": "flex", "gap": "4px"}),
            html.Div([
                html.Div(dcc.Graph(id="prop-bin-plot",   style=_sg), style={"flex": "1"}),
                html.Div(dcc.Graph(id="prop-state-plot", style=_sg), style={"flex": "1"}),
            ], style={"display": "flex", "gap": "4px", "marginTop": "4px"}),
            dcc.Graph(id="cdf-plot", style={"height": "280px", "marginTop": "4px"}),
        ], style={"flex": "0 0 45%", "minWidth": 0}),
    ], style={"display": "flex", "gap": "6px"}),

], style={"maxWidth": "1800px", "margin": "0 auto", "padding": "10px"})


# ── callback ───────────────────────────────────────────────────────────────────

@app.callback(
    Output("main-plot",       "figure"),
    Output("cnt-bin-plot",    "figure"),
    Output("cnt-state-plot",  "figure"),
    Output("prop-bin-plot",   "figure"),
    Output("prop-state-plot", "figure"),
    Output("cdf-plot",        "figure"),
    Input("t-b",          "value"),
    Input("t-c",          "value"),
    Input("t-xy",         "value"),
    Input("t-yz",         "value"),
    Input("n-cells",      "value"),
    Input("lt1", "value"), Input("lt2", "value"), Input("lt3", "value"),
    Input("lc1", "value"), Input("lc2", "value"), Input("lc3", "value"),
    Input("show-unlabeled", "value"),
    Input("show-traj2",     "value"),
    Input("min-circ",     "value"),
    Input("max-circ",     "value"),
    Input("cdf-style",    "value"),
)
def update(t_b, t_c, t_xy, t_yz, n_cells,
           lt1, lt2, lt3, lc1, lc2, lc3,
           show_unlabeled, show_traj2, min_circ, max_circ, cdf_style):

    # ── parse inputs ───────────────────────────────────────────────────────
    t_b      = int(t_b    or 12)
    t_c_dur  = int(t_c    or 12)
    t_xy_val = int(t_xy   or 6)
    t_yz_dur = int(t_yz   or 18)
    n_cells  = int(n_cells or DEFAULT_N)

    t_c_abs  = t_b + t_c_dur
    t_yz_abs = t_xy_val + t_yz_dur

    show_t2 = "show" in (show_traj2 or [])

    raw_colors     = [lc1 or DEFAULT_LABEL_COLORS[0],
                      lc2 or DEFAULT_LABEL_COLORS[1],
                      lc3 or DEFAULT_LABEL_COLORS[2]]
    label_pairs    = [(float(t), c) for t, c in zip([lt1, lt2, lt3], raw_colors) if t is not None]
    label_times    = [p[0] for p in label_pairs]
    label_col_list = [p[1] for p in label_pairs]
    lt_to_color    = dict(label_pairs)

    def cell_color(label):
        return lt_to_color.get(label, UNLABELED_COLOR)

    # ── simulate ───────────────────────────────────────────────────────────
    _circ = (float(min_circ or 12), float(max_circ or 24))
    cells  = simulate_cells(n_cells, *_circ, seed=42)
    labels = [get_label(ce, te, label_times) for ce, te in cells]
    states = [get_final_state(te, t_b, t_c_abs) for _, te in cells]

    state_names   = ["A", "B", "C"]
    transitions_1 = [0, t_b, t_c_abs]

    state_names_2 = ["X", "Y", "Z"]
    transitions_2 = [0, t_xy_val, t_yz_abs]

    if show_t2:
        # Independent population for XYZ — different seed
        cells2  = simulate_cells(n_cells, *_circ, seed=99)
        labels2 = [get_label(ce, te, label_times) for ce, te in cells2]
        states2 = [get_final_state_2(te, t_xy_val, t_yz_abs) for _, te in cells2]

        # Merge and sort all cells by tumor-entry time for the main plot
        combined = sorted(
            [("ABC", cell, lbl, st) for cell, lbl, st in zip(cells,  labels,  states )] +
            [("XYZ", cell, lbl, st) for cell, lbl, st in zip(cells2, labels2, states2)],
            key=lambda x: x[1][1],          # sort by te (ascending = most time in tumor first)
        )
        abc_ypos = [i + 1 for i, (pop, _, _, _) in enumerate(combined) if pop == "ABC"]
        xyz_ypos = [i + 1 for i, (pop, _, _, _) in enumerate(combined) if pop == "XYZ"]
        total_rows = len(combined)
        _abc_n = _xyz_n = 0
        y_tick_labels = []
        for pop, _, _, _ in combined:
            if pop == "ABC":
                _abc_n += 1; y_tick_labels.append(f"A-{_abc_n:02d}")
            else:
                _xyz_n += 1; y_tick_labels.append(f"X-{_xyz_n:02d}")

    # ── main trajectory figure ─────────────────────────────────────────────
    fig = go.Figure()

    # legend: label colors
    for lt, lc in zip(label_times, label_col_list):
        fig.add_trace(go.Scatter(x=[None], y=[None], mode="lines+markers",
                                 line=dict(color=lc, width=2), marker=dict(size=8, color=lc),
                                 name=f"Labeled {int(lt)}hr", showlegend=True))
    fig.add_trace(go.Scatter(x=[None], y=[None], mode="lines",
                             line=dict(color=UNLABELED_COLOR, width=1.5, dash="dash"),
                             name="Unlabeled / pre-label", showlegend=True))

    # legend: traj2 trajectory label (no separate state swatches — rects use label colors)
    if show_t2:
        fig.add_trace(go.Scatter(x=[None], y=[None], mode="lines",
                                 line=dict(color=STATE_COLORS_2["Y"], width=2),
                                 name="XYZ trajectory", showlegend=True))

    # ── helper: batch circulation lines for one population ────────────────
    def collect_circ_lines(pop_cells, pop_labels, ypos):
        gx, gy = [], []
        cs = defaultdict(lambda: {"x": [], "y": []})
        for (ce, te), lbl, y in zip(pop_cells, pop_labels, ypos):
            col      = cell_color(lbl)
            gray_end = lbl if (lbl is not None and ce < lbl) else te
            if ce < gray_end:
                gx += [ce, gray_end, None]
                gy += [y, y, None]
            if lbl is not None and lbl <= te:
                cs[col]["x"] += [lbl, te, None]
                cs[col]["y"] += [y, y, None]
        return gx, gy, cs

    # Assign y positions: when traj2 on, merge-sort; otherwise sequential
    if show_t2:
        abc_y = abc_ypos
        xyz_y = xyz_ypos
        total_rows = total_rows       # set above in simulate block
    else:
        abc_y      = list(range(1, n_cells + 1))
        xyz_y      = []
        total_rows = n_cells

    gx1, gy1, cs1 = collect_circ_lines(cells, labels, abc_y)
    gray_x, gray_y = gx1, gy1
    col_segs = cs1

    if show_t2:
        gx2, gy2, cs2 = collect_circ_lines(cells2, labels2, xyz_y)
        gray_x += gx2
        gray_y += gy2
        for col, seg in cs2.items():
            col_segs[col]["x"] += seg["x"]
            col_segs[col]["y"] += seg["y"]

    if gray_x:
        fig.add_trace(go.Scatter(x=gray_x, y=gray_y, mode="lines",
                                 line=dict(color="#cccccc", width=1.5, dash="dash"),
                                 showlegend=False, hoverinfo="skip"))
    for col, seg in col_segs.items():
        fig.add_trace(go.Scatter(x=seg["x"], y=seg["y"], mode="lines",
                                 line=dict(color=col, width=2.5, dash="dash"),
                                 showlegend=False, hoverinfo="skip"))

    # Entry circles
    fig.add_trace(go.Scatter(x=[ce for ce, _ in cells], y=abc_y,
                             mode="markers",
                             marker=dict(symbol="circle", size=6, color="#cccccc",
                                         line=dict(color="#888", width=1.5)),
                             showlegend=False, hoverinfo="skip"))
    if show_t2:
        fig.add_trace(go.Scatter(x=[ce for ce, _ in cells2], y=xyz_y,
                                 mode="markers",
                                 marker=dict(symbol="circle", size=6, color="#cccccc",
                                             line=dict(color="#888", width=1.5)),
                                 showlegend=False, hoverinfo="skip"))

    # Rects
    all_annots = []
    rect_segs_1, annots_1 = draw_rect_traces(
        cells, labels, cell_color, state_names, transitions_1, -0.20, +0.20,
        y_positions=abc_y)
    all_annots += annots_1
    for col, seg in rect_segs_1.items():
        fig.add_trace(go.Scatter(x=seg["x"], y=seg["y"], mode="lines",
                                 fill="toself", fillcolor=hex_to_rgba(col, FILL_ALPHA),
                                 line=dict(color=col, width=2),
                                 showlegend=False, hoverinfo="skip"))

    if show_t2:
        rect_segs_2, annots_2 = draw_rect_traces(
            cells2, labels2, cell_color, state_names_2, transitions_2, -0.20, +0.20,
            y_positions=xyz_y)
        all_annots += annots_2
        for col, seg in rect_segs_2.items():
            fig.add_trace(go.Scatter(x=seg["x"], y=seg["y"], mode="lines",
                                     fill="toself", fillcolor=hex_to_rgba(col, FILL_ALPHA),
                                     line=dict(color=col, width=2),
                                     showlegend=False, hoverinfo="skip"))

    traj1_label = f"A (0–{t_b}hr) → B ({t_b}–{t_c_abs}hr) → C ({t_c_abs}hr+)"
    traj2_label = (f"  |  X (0–{t_xy_val}hr) → Y ({t_xy_val}–{t_yz_abs}hr) → Z ({t_yz_abs}hr+)"
                   if show_t2 else "")

    plot_h = max(400, total_rows * 16 + 100)

    if show_t2:
        y_ticks  = list(range(1, total_rows + 1))
        y_labels = y_tick_labels
        y_range  = [0.4, total_rows + 0.6]
    else:
        y_ticks  = list(range(1, n_cells + 1))
        y_labels = [f"Cell {i:02d}" for i in range(1, n_cells + 1)]
        y_range  = [0.4, n_cells + 0.6]

    fig.update_layout(
        title=dict(text=traj1_label + traj2_label, x=0.5, font=dict(size=12)),
        xaxis=dict(title="Absolute Time (hours)", range=[-50, 3],
                   tickvals=list(range(-48, 1, 6)), gridcolor="#eeeeee", zeroline=False),
        yaxis=dict(title="Cell", range=y_range,
                   tickvals=y_ticks, ticktext=y_labels,
                   gridcolor="#eeeeee"),
        plot_bgcolor="white", paper_bgcolor="white",
        legend=dict(title="Labeling", x=1.01, y=1, xanchor="left",
                    bgcolor="rgba(255,255,255,0.8)", bordercolor="#ccc", borderwidth=1),
        margin=dict(l=90, r=130, t=50, b=50),
        height=plot_h,
    )
    for lt, lc in zip(label_times, label_col_list):
        fig.add_vline(x=lt, line_dash="dot", line_color=lc, line_width=1.2,
                      annotation_text=f"{int(lt)}hr", annotation_position="top",
                      annotation_font=dict(color=lc, size=10))
    fig.add_vline(x=0, line_color="black", line_width=1.5,
                  annotation_text="Harvest", annotation_position="top",
                  annotation_font=dict(size=10))
    for ann in all_annots:
        fig.add_annotation(**ann)

    # ── shared bin/count data ──────────────────────────────────────────────
    include_unlabeled = "show" in (show_unlabeled or [])
    all_bins   = label_times + ([None] if include_unlabeled else [])
    bin_names  = [f"{int(lt)}hr" for lt in label_times] + (["Unlabeled"] if include_unlabeled else [])
    bin_colors = label_col_list + ([UNLABELED_COLOR]     if include_unlabeled else [])

    counts = {b: {s: 0 for s in state_names} for b in all_bins}
    for label, state in zip(labels, states):
        if label in counts:
            counts[label][state] += 1

    counts2 = {b: {s: 0 for s in state_names_2} for b in all_bins}
    if show_t2:
        for label, state2 in zip(labels2, states2):
            if label in counts2:
                counts2[label][state2] += 1

    _base_layout = dict(
        plot_bgcolor="white", paper_bgcolor="white",
        margin=dict(l=45, r=20, t=35, b=35),
        legend=dict(x=1.01, y=1, xanchor="left", font=dict(size=10),
                    bgcolor="rgba(0,0,0,0)", borderwidth=0),
    )

    # ── count by time bin (all states stacked together) ────────────────────
    cbin   = go.Figure()
    base_c = [0] * len(all_bins)
    for s, sc in STATE_COLORS.items():
        y = [counts[b][s] for b in all_bins]
        cbin.add_trace(go.Bar(name=f"State {s}", x=bin_names, y=y, base=base_c[:],
                              marker_color=sc))
        base_c = [b + v for b, v in zip(base_c, y)]
    if show_t2:
        for s, sc in STATE_COLORS_2.items():
            y2 = [counts2[b][s] for b in all_bins]
            cbin.add_trace(go.Bar(name=f"State {s}", x=bin_names, y=y2, base=base_c[:],
                                  marker_color=sc))
            base_c = [b + v for b, v in zip(base_c, y2)]
    cbin.update_layout(barmode="overlay", title="Count by time bin",
                       yaxis_title="# cells", xaxis_title="Time bin",
                       yaxis=dict(gridcolor="#eeeeee"),
                       **_base_layout)

    # ── count by state (grouped, XYZ bars use blue gradient) ──────────────
    all_state_names = state_names + (state_names_2 if show_t2 else [])
    cstate = go.Figure()
    for b, bname, bc in zip(all_bins, bin_names, bin_colors):
        y_vals = [counts[b][s] for s in state_names]
        if show_t2:
            y_vals += [counts2[b][s] for s in state_names_2]
        cstate.add_trace(go.Bar(name=bname, x=all_state_names, y=y_vals,
                                marker_color=bc, legendgroup=f"bn{bname}"))
    cstate.update_layout(barmode="group", title="Count by state",
                         yaxis_title="# cells", xaxis_title="State",
                         yaxis=dict(gridcolor="#eeeeee"),
                         **_base_layout)

    # ── proportion by time bin (all states stacked, combined denominator) ──
    pbin  = go.Figure()
    if show_t2:
        tot_b = [sum(counts[b].values()) + sum(counts2[b].values()) for b in all_bins]
    else:
        tot_b = [sum(counts[b].values()) for b in all_bins]
    base_p = [0.0] * len(all_bins)
    for s, sc in STATE_COLORS.items():
        prop = [counts[b][s] / t if t > 0 else 0 for b, t in zip(all_bins, tot_b)]
        pbin.add_trace(go.Bar(name=f"State {s}", x=bin_names, y=prop, base=base_p[:],
                               marker_color=sc))
        base_p = [b + p for b, p in zip(base_p, prop)]
    if show_t2:
        for s, sc in STATE_COLORS_2.items():
            prop2 = [counts2[b][s] / t if t > 0 else 0 for b, t in zip(all_bins, tot_b)]
            pbin.add_trace(go.Bar(name=f"State {s}", x=bin_names, y=prop2, base=base_p[:],
                                   marker_color=sc))
            base_p = [b + p for b, p in zip(base_p, prop2)]
    pbin.update_layout(
        barmode="overlay",
        title="Proportion by time bin",
        yaxis_title="Proportion", xaxis_title="Time bin",
        yaxis=dict(gridcolor="#eeeeee", range=[0, 1.05]),
        **_base_layout,
    )

    # ── proportion by state (stacked) ─────────────────────────────────────
    tot_s  = [sum(counts[b][s]  for b in all_bins) for s in state_names]
    tot_s2 = [sum(counts2[b][s] for b in all_bins) for s in state_names_2]
    base   = [0.0] * len(all_state_names)
    pstate = go.Figure()
    for b, bname, bc in zip(all_bins, bin_names, bin_colors):
        prop = [counts[b][s] / t if t > 0 else 0 for s, t in zip(state_names, tot_s)]
        if show_t2:
            prop += [counts2[b][s] / t if t > 0 else 0 for s, t in zip(state_names_2, tot_s2)]
        pstate.add_trace(go.Bar(name=bname, x=all_state_names, y=prop, base=base[:],
                                 marker_color=bc, legendgroup=f"bn{bname}"))
        base = [b + p for b, p in zip(base, prop)]
    pstate.update_layout(barmode="overlay", title="Proportion by state",
                         yaxis_title="Proportion", xaxis_title="State",
                         yaxis=dict(gridcolor="#eeeeee", range=[0, 1.05]),
                         **_base_layout)

    # ── CDF ────────────────────────────────────────────────────────────────
    sorted_lts_rev = sorted(label_times, reverse=True)
    cdf_fig = go.Figure()

    def add_cdf_traces(lbl_list, st_list, color_map):
        for s, sc in color_map.items():
            cell_labels_s = [lbl for lbl, st in zip(lbl_list, st_list) if st == s]
            n_labeled     = sum(1 for lbl in cell_labels_s if lbl is not None)
            if n_labeled == 0 or not sorted_lts_rev:
                continue
            cum, cum_props = 0, []
            for lt in sorted_lts_rev:
                cum += sum(1 for lbl in cell_labels_s if lbl == lt)
                cum_props.append(cum / n_labeled)
            cum_props[-1] = 1.0
            x_cdf = [0.0] + sorted_lts_rev
            y_cdf = [0.0] + cum_props
            if cdf_style == "hv":
                auc = sum(y_cdf[i] * abs(x_cdf[i+1] - x_cdf[i]) for i in range(len(x_cdf) - 1))
            else:
                auc = sum((y_cdf[i] + y_cdf[i+1]) / 2 * abs(x_cdf[i+1] - x_cdf[i])
                          for i in range(len(x_cdf) - 1))
            cdf_fig.add_trace(go.Scatter(
                x=x_cdf, y=y_cdf, mode="lines",
                line=dict(color=sc, width=2.5, shape=cdf_style),
                name=f"State {s}  (AUC = {np.abs(x_cdf[-1]) - auc:.1f})",
            ))

    add_cdf_traces(labels, states,  STATE_COLORS)
    if show_t2:
        add_cdf_traces(labels, states2, STATE_COLORS_2)

    x_last = sorted_lts_rev[-1] if sorted_lts_rev else -36
    cdf_fig.add_trace(go.Scatter(
        x=[0.0, x_last], y=[0.0, 1.0], mode="lines",
        line=dict(color="black", width=1.5, dash="dash"),
        name="Uniform reference",
    ))
    x_min = x_last - 3
    cdf_fig.update_layout(
        title="CDF: cumulative labeled fraction by state (labeled cells only)",
        xaxis=dict(title="Absolute time (hr)", range=[3, x_min],
                   tickvals=[0] + sorted_lts_rev, gridcolor="#eeeeee"),
        yaxis=dict(title="Cumul. fraction", range=[0, 1.05], gridcolor="#eeeeee"),
        plot_bgcolor="white", paper_bgcolor="white",
        legend=dict(x=1.01, y=1, xanchor="left", font=dict(size=10),
                    bgcolor="rgba(0,0,0,0)", borderwidth=0),
        margin=dict(l=55, r=20, t=40, b=45),
    )
    for lt, lc in zip(label_times, label_col_list):
        cdf_fig.add_vline(x=lt, line_dash="dot", line_color=lc, line_width=1.2,
                          annotation_text=f"{int(lt)}hr", annotation_position="top",
                          annotation_font=dict(color=lc, size=10))
    cdf_fig.add_vline(x=0, line_color="black", line_width=1.5)

    return fig, cbin, cstate, pbin, pstate, cdf_fig


if __name__ == "__main__":
    app.run(debug=True, port=8050)
