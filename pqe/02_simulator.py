"""
02_simulator.py  –  Zman-seq interactive trajectory simulator (v4)

Cells have:
  - circ_entry : when the cell enters blood (random offset before tumor entry)
  - tumor_entry: when it enters the TME (uniform across -48 to -1 hr)

Harvest is at t=0. Cells are colored by the latest dye-injection time point
they were in circulation for. That color drives the dashed line and the fill/
border of the tumor state rectangles.
"""

import dash
from dash import dcc, html, Input, Output
from plotly.subplots import make_subplots
from collections import defaultdict
import plotly.graph_objects as go
import numpy as np

# ── defaults ──────────────────────────────────────────────────────────────────

DEFAULT_N            = 20
DEFAULT_LABEL_TIMES  = [-36, -24, -12]
DEFAULT_LABEL_COLORS = ["#115e59", "#0d9488", "#5eead4"]  # dark → mid → light teal (-36, -24, -12)
STATE_COLORS         = {"A": "#c4b5fd", "B": "#7c3aed", "C": "#2e1065"}  # light→dark purple
UNLABELED_COLOR      = "#aaaaaa"
FILL_ALPHA           = 0.70

# ── colour helpers ────────────────────────────────────────────────────────────

def hex_to_rgba(hex_color: str, alpha: float = FILL_ALPHA) -> str:
    h = hex_color.lstrip("#")
    r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
    return f"rgba({r},{g},{b},{alpha})"


def text_color_for(hex_color: str) -> str:
    h = hex_color.lstrip("#")
    r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
    return "white" if 0.299*r + 0.587*g + 0.114*b < 140 else "#333"

# ── simulation ────────────────────────────────────────────────────────────────

def simulate_cells(n: int, min_circ: float = 12.0, max_circ: float = 24.0):
    rng = np.random.default_rng(42)
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

# ── app layout ────────────────────────────────────────────────────────────────

app = dash.Dash(__name__)

_ls = {"fontFamily": "sans-serif", "fontSize": "13px"}
_sm = {i: str(i) for i in range(0, 49, 6)}

# shared small-graph style
_sg = {"height": "220px"}

app.layout = html.Div([
    html.H2("Zman-seq Trajectory Simulator",
            style={"textAlign": "center", "fontFamily": "sans-serif", "marginBottom": "8px"}),

    # ── control panel ─────────────────────────────────────────────────────
    html.Div([

        # trajectory transitions
        html.Div([
            html.B("Trajectory transitions (hr after tumor entry)", style=_ls),
            html.Div([
                html.Label("State A → B:", style={**_ls, "width": "110px", "display": "inline-block"}),
                dcc.Slider(id="t-b", min=1, max=47, step=1, value=12, marks=_sm,
                           tooltip={"placement": "bottom", "always_visible": True}),
            ], style={"marginTop": "8px"}),
            html.Div([
                html.Label("State B duration:", style={**_ls, "width": "110px", "display": "inline-block"}),
                dcc.Slider(id="t-c", min=1, max=48, step=1, value=12, marks=_sm,
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

        # bar chart + CDF options
        html.Div([
            html.B("Bar chart options", style=_ls),
            dcc.Checklist(
                id="show-unlabeled",
                options=[{"label": " Include unlabeled cells", "value": "show"}],
                value=["show"],
                style={**_ls, "marginTop": "10px"},
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
        ], style={"flex": "0 0 200px", "paddingRight": "20px"}),

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

    # ── plot area ─────────────────────────────────────────────────────────
    html.Div([

        # main trajectory
        html.Div(dcc.Graph(id="main-plot"), style={"flex": "0 0 55%", "minWidth": 0}),

        # right panel: 2×2 bar grid + CDF below
        html.Div([
            # row 1: count charts
            html.Div([
                html.Div(dcc.Graph(id="cnt-bin-plot",   style=_sg), style={"flex": "1"}),
                html.Div(dcc.Graph(id="cnt-state-plot", style=_sg), style={"flex": "1"}),
            ], style={"display": "flex", "gap": "4px"}),
            # row 2: proportion charts
            html.Div([
                html.Div(dcc.Graph(id="prop-bin-plot",   style=_sg), style={"flex": "1"}),
                html.Div(dcc.Graph(id="prop-state-plot", style=_sg), style={"flex": "1"}),
            ], style={"display": "flex", "gap": "4px", "marginTop": "4px"}),
            # row 3: CDF
            dcc.Graph(id="cdf-plot", style={"height": "280px", "marginTop": "4px"}),
        ], style={"flex": "0 0 45%", "minWidth": 0}),

    ], style={"display": "flex", "gap": "6px"}),

], style={"maxWidth": "1800px", "margin": "0 auto", "padding": "10px"})


# ── callback ──────────────────────────────────────────────────────────────────

@app.callback(
    Output("main-plot",       "figure"),
    Output("cnt-bin-plot",    "figure"),
    Output("cnt-state-plot",  "figure"),
    Output("prop-bin-plot",   "figure"),
    Output("prop-state-plot", "figure"),
    Output("cdf-plot",        "figure"),
    Input("t-b",          "value"),
    Input("t-c",          "value"),
    Input("n-cells",      "value"),
    Input("lt1", "value"), Input("lt2", "value"), Input("lt3", "value"),
    Input("lc1", "value"), Input("lc2", "value"), Input("lc3", "value"),
    Input("show-unlabeled", "value"),
    Input("min-circ",     "value"),
    Input("max-circ",     "value"),
    Input("cdf-style",    "value"),
)
def update(t_b, t_c, n_cells, lt1, lt2, lt3, lc1, lc2, lc3,
           show_unlabeled, min_circ, max_circ, cdf_style):

    # ── parse inputs ──────────────────────────────────────────────────────
    t_b     = int(t_b     or 12)
    t_c_dur = int(t_c     or 12)
    n_cells = int(n_cells or DEFAULT_N)
    t_c_abs = t_b + t_c_dur

    raw_colors = [lc1 or DEFAULT_LABEL_COLORS[0],
                  lc2 or DEFAULT_LABEL_COLORS[1],
                  lc3 or DEFAULT_LABEL_COLORS[2]]
    label_pairs    = [(float(t), c) for t, c in zip([lt1, lt2, lt3], raw_colors) if t is not None]
    label_times    = [p[0] for p in label_pairs]
    label_col_list = [p[1] for p in label_pairs]
    lt_to_color    = dict(label_pairs)

    def cell_color(label):
        return lt_to_color.get(label, UNLABELED_COLOR)

    # ── simulate ──────────────────────────────────────────────────────────
    cells  = simulate_cells(n_cells, float(min_circ or 12), float(max_circ or 24))
    labels = [get_label(ce, te, label_times) for ce, te in cells]
    states = [get_final_state(te, t_b, t_c_abs) for _, te in cells]

    state_names = ["A", "B", "C"]
    transitions = [0, t_b, t_c_abs]

    # ── main trajectory figure ────────────────────────────────────────────
    fig   = go.Figure()
    rect_h = 0.20

    # legend dummies
    for lt, lc in zip(label_times, label_col_list):
        fig.add_trace(go.Scatter(x=[None], y=[None], mode="lines+markers",
                                 line=dict(color=lc, width=2), marker=dict(size=8, color=lc),
                                 name=f"Labeled {int(lt)}hr", showlegend=True))
    fig.add_trace(go.Scatter(x=[None], y=[None], mode="lines",
                             line=dict(color=UNLABELED_COLOR, width=1.5, dash="dash"),
                             name="Unlabeled / pre-label", showlegend=True))

    # dashed circulation lines (batched, added first → behind rectangles)
    gray_x, gray_y = [], []
    col_segs = defaultdict(lambda: {"x": [], "y": []})
    for i, ((ce, te), label) in enumerate(zip(cells, labels)):
        y        = i + 1
        col      = cell_color(label)
        gray_end = label if (label is not None and ce < label) else te
        if ce < gray_end:
            gray_x += [ce, gray_end, None]
            gray_y += [y, y, None]
        if label is not None and label <= te:
            col_segs[col]["x"] += [label, te, None]
            col_segs[col]["y"] += [y, y, None]

    if gray_x:
        fig.add_trace(go.Scatter(x=gray_x, y=gray_y, mode="lines",
                                 line=dict(color="#cccccc", width=1.5, dash="dash"),
                                 showlegend=False, hoverinfo="skip"))
    for col, seg in col_segs.items():
        fig.add_trace(go.Scatter(x=seg["x"], y=seg["y"], mode="lines",
                                 line=dict(color=col, width=2.5, dash="dash"),
                                 showlegend=False, hoverinfo="skip"))

    # circulation entry circles
    fig.add_trace(go.Scatter(x=[ce for ce, _ in cells], y=list(range(1, n_cells+1)),
                             mode="markers",
                             marker=dict(symbol="circle", size=6, color="#cccccc",
                                         line=dict(color="#888", width=1.5)),
                             showlegend=False, hoverinfo="skip"))

    # tumor rectangles (after dashes → on top)
    rect_segs   = defaultdict(lambda: {"x": [], "y": []})
    state_annots = []
    for i, ((ce, te), label) in enumerate(zip(cells, labels)):
        y   = i + 1
        col = cell_color(label)
        txt = text_color_for(col)
        for j, sname in enumerate(state_names):
            a_s = te + transitions[j]
            a_e = (te + transitions[j+1]) if j+1 < len(transitions) else 0.0
            a_s = max(a_s, te);  a_e = min(a_e, 0.0)
            if a_s >= a_e:
                continue
            rect_segs[col]["x"] += [a_s, a_e, a_e, a_s, a_s, None]
            rect_segs[col]["y"] += [y-rect_h, y-rect_h, y+rect_h, y+rect_h, y-rect_h, None]
            if a_e - a_s > 0.5:
                state_annots.append(dict(
                    x=(a_s+a_e)/2, y=y, text=f"<b>{sname}</b>", showarrow=False,
                    font=dict(size=9 if a_e-a_s < 3 else 10, color=txt),
                    xanchor="center", yanchor="middle"))

    for col, seg in rect_segs.items():
        fig.add_trace(go.Scatter(x=seg["x"], y=seg["y"], mode="lines",
                                 fill="toself", fillcolor=hex_to_rgba(col, FILL_ALPHA),
                                 line=dict(color=col, width=2),
                                 showlegend=False, hoverinfo="skip"))

    plot_h = max(400, n_cells * 16 + 100)
    fig.update_layout(
        title=dict(text=f"A (0–{t_b}hr) → B ({t_b}–{t_c_abs}hr) → C ({t_c_abs}hr+) after tumor entry",
                   x=0.5, font=dict(size=12)),
        xaxis=dict(title="Absolute Time (hours)", range=[-50, 3],
                   tickvals=list(range(-48, 1, 6)), gridcolor="#eeeeee", zeroline=False),
        yaxis=dict(title="Cell", range=[0.4, n_cells+0.6],
                   tickvals=list(range(1, n_cells+1)),
                   ticktext=[f"Cell {i:02d}" for i in range(1, n_cells+1)],
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
    for ann in state_annots:
        fig.add_annotation(**ann)

    # ── shared bin/count data ─────────────────────────────────────────────
    include_unlabeled = "show" in (show_unlabeled or [])
    all_bins   = label_times + ([None] if include_unlabeled else [])
    bin_names  = [f"{int(lt)}hr" for lt in label_times] + (["Unlabeled"] if include_unlabeled else [])
    bin_colors = label_col_list + ([UNLABELED_COLOR]     if include_unlabeled else [])

    counts = {b: {s: 0 for s in state_names} for b in all_bins}
    for label, state in zip(labels, states):
        if label in counts:
            counts[label][state] += 1

    _base_layout = dict(
        plot_bgcolor="white", paper_bgcolor="white",
        margin=dict(l=45, r=20, t=35, b=35),
        legend=dict(x=1.01, y=1, xanchor="left", font=dict(size=10),
                    bgcolor="rgba(0,0,0,0)", borderwidth=0),
    )

    # ── count by time bin (grouped) ───────────────────────────────────────
    cbin = go.Figure()
    for s, sc in STATE_COLORS.items():
        cbin.add_trace(go.Bar(name=f"State {s}",
                              x=bin_names, y=[counts[b][s] for b in all_bins],
                              marker_color=sc, legendgroup=f"st{s}"))
    cbin.update_layout(barmode="group", title="Count by time bin",
                       yaxis_title="# cells", xaxis_title="Time bin",
                       yaxis=dict(gridcolor="#eeeeee"),
                       **_base_layout)

    # ── count by state (grouped) ──────────────────────────────────────────
    cstate = go.Figure()
    for b, bname, bc in zip(all_bins, bin_names, bin_colors):
        cstate.add_trace(go.Bar(name=bname,
                                x=state_names, y=[counts[b][s] for s in state_names],
                                marker_color=bc, legendgroup=f"bn{bname}"))
    cstate.update_layout(barmode="group", title="Count by state",
                         yaxis_title="# cells", xaxis_title="State",
                         yaxis=dict(gridcolor="#eeeeee"),
                         **_base_layout)

    # ── proportion by time bin (stacked, manual base) ─────────────────────
    pbin  = go.Figure()
    tot_b = [sum(counts[b].values()) for b in all_bins]
    base  = [0.0] * len(all_bins)
    for s, sc in STATE_COLORS.items():
        prop = [counts[b][s] / t if t > 0 else 0 for b, t in zip(all_bins, tot_b)]
        pbin.add_trace(go.Bar(name=f"State {s}", x=bin_names, y=prop, base=base[:],
                               marker_color=sc, legendgroup=f"st{s}"))
        base = [b + p for b, p in zip(base, prop)]
    pbin.update_layout(barmode="overlay", title="Proportion by time bin",
                       yaxis_title="Proportion", xaxis_title="Time bin",
                       yaxis=dict(gridcolor="#eeeeee", range=[0, 1.05]),
                       **_base_layout)

    # ── proportion by state (stacked, manual base) ────────────────────────
    pstate  = go.Figure()
    tot_s   = [sum(counts[b][s] for b in all_bins) for s in state_names]
    base    = [0.0] * len(state_names)
    for b, bname, bc in zip(all_bins, bin_names, bin_colors):
        prop = [counts[b][s] / t if t > 0 else 0 for s, t in zip(state_names, tot_s)]
        pstate.add_trace(go.Bar(name=bname, x=state_names, y=prop, base=base[:],
                                 marker_color=bc, legendgroup=f"bn{bname}"))
        base = [b + p for b, p in zip(base, prop)]
    pstate.update_layout(barmode="overlay", title="Proportion by state",
                         yaxis_title="Proportion", xaxis_title="State",
                         yaxis=dict(gridcolor="#eeeeee", range=[0, 1.05]),
                         **_base_layout)

    # ── CDF figure ────────────────────────────────────────────────────────
    # x-axis runs from 0 → most-recent label → ... → earliest label (reversed time)
    # Denominator: only LABELED cells in each state (unlabeled excluded)
    # Cumulation: most-recent label first (e.g. -12, then -24, then -36)
    sorted_lts_rev = sorted(label_times, reverse=True)   # e.g. [-12, -24, -36]

    cdf_fig = go.Figure()
    for s, sc in STATE_COLORS.items():
        cell_labels_s = [lbl for lbl, st in zip(labels, states) if st == s]
        n_labeled     = sum(1 for lbl in cell_labels_s if lbl is not None)
        if n_labeled == 0 or not sorted_lts_rev:
            continue

        # Cumulative proportion, most-recent label first; last point forced to 1.0
        cum, cum_props = 0, []
        for lt in sorted_lts_rev:
            cum += sum(1 for lbl in cell_labels_s if lbl == lt)
            cum_props.append(cum / n_labeled)
        cum_props[-1] = 1.0   # ensure the last point is exactly 100 %

        # Points: (0, 0), (-12, p1), (-24, p1+p2), (-36, 1.0)
        x_cdf = [0.0] + sorted_lts_rev
        y_cdf = [0.0] + cum_props

        # AUC: trapezoid rule for diagonal lines, left-endpoint rectangles for hv steps
        if cdf_style == "hv":
            auc = sum(
                y_cdf[i] * abs(x_cdf[i+1] - x_cdf[i])
                for i in range(len(x_cdf) - 1)
            )
        else:
            auc = sum(
                (y_cdf[i] + y_cdf[i+1]) / 2 * abs(x_cdf[i+1] - x_cdf[i])
                for i in range(len(x_cdf) - 1)
            )

        line_shape = cdf_style  # "linear" or "hv"
        cdf_fig.add_trace(go.Scatter(
            x=x_cdf, y=y_cdf, mode="lines",
            line=dict(color=sc, width=2.5, shape=line_shape),
            name=f"State {s}  (AUC = {np.abs(x_cdf[-1]) - auc:.1f})",
        ))

    # Diagonal reference line: (0, 0) → (earliest label time, 1.0)
    x_last = sorted_lts_rev[-1] if sorted_lts_rev else -36
    cdf_fig.add_trace(go.Scatter(
        x=[0.0, x_last], y=[0.0, 1.0], mode="lines",
        line=dict(color="black", width=1.5, dash="dash"),
        name="Uniform reference",
    ))

    # x-axis: 0 on left, most-negative on right
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
