
# Shared figure library for the glioma radiation-signature manuscript figures.
import os, json, csv, textwrap
import numpy as np, pandas as pd
import matplotlib as mpl, matplotlib.pyplot as plt
from PIL import Image

_AP = json.load(open('handoff/artifact_paths.json'))
ART, PATHS = _AP['ART'], _AP['PATHS']
def P(name): return PATHS[name]
def A(name): return ART[name]

# ---- fixed palette (consistent across all figures) ----
HIGH = '#D55E00'   # high risk  (warm, vermilion)
LOW  = '#0072B2'   # low risk   (cool, blue)
COH  = {'TCGA': '#009E73', 'CGGA693': '#CC79A7', 'CGGA325': '#56B4E9'}
COH_LABEL = {'TCGA': 'TCGA-GBMLGG', 'CGGA693': 'CGGA-693', 'CGGA325': 'CGGA-325'}
SUB  = {'RS1': '#4393C3', 'RS2': '#F4A582', 'RS3': '#B2182B'}
RT_ARM = {'RT': '#000000', 'no_RT': '#8C8C8C'}
GREY = '#7F7F7F'; DARK = '#333333'
POS = '#B2182B'; NEG = '#2166AC'   # signed coefficients / HR>1 vs <1 (CVD-safe red/blue)
SEQ_CMAP = 'viridis'; DIV_CMAP = 'RdBu_r'

SIZE = (3.5, 3.1)          # one panel size for every figure
RECT = [0.21, 0.17, 0.74, 0.74]  # one axes rectangle per standard panel
DPI = 600
MANIFEST = []

def setup():
    apply_figure_style(sizes=(8, 7, 6))
    apply_sci_style(base=8, dpi=300)
    mpl.rcParams['savefig.bbox'] = None
    mpl.rcParams['axes.unicode_minus'] = False

def new_fig(rect=None, size=None):
    fig = plt.figure(figsize=size or SIZE)
    ax = fig.add_axes(rect or RECT)
    return fig, ax

def wide_size(ncols=2, gap_px=90):
    """Figure size spanning ncols standard panels + gaps (so the assembler never rescales)."""
    return (SIZE[0]*ncols + gap_px/DPI*(ncols-1), SIZE[1])

def label(ax, letter, dx=None, dy=1.04):
    # place relative to the figure so that labels sit at the same spot on every panel
    fig = ax.figure
    fig.text(0.015, 0.975, str(letter).upper(), fontsize=14, fontweight='bold', family='sans-serif', ha='left', va='top')

def finish(fig, figid, letter, sources, ax_for_label=None, allow_cross=True):
    """Gate (English-only + overlaps) then save FigN<L>.png (600 dpi) + .pdf on a fixed canvas."""
    label(ax_for_label or fig.axes[0], letter)
    fig.canvas.draw()
    res = check_figure(fig)
    assert not res['overlaps'], (figid, letter, res['overlaps'])
    os.makedirs('figures/panels', exist_ok=True)
    stem = f'figures/panels/{figid}{letter.upper()}'
    with mpl.rc_context({'savefig.bbox': None}):
        fig.savefig(stem + '.png', dpi=DPI, facecolor='white')
        fig.savefig(stem + '.pdf', facecolor='white')
    plt.close(fig)
    for s in (sources if isinstance(sources, (list, tuple)) else [sources]):
        MANIFEST.append({'figure': figid, 'panel': letter.upper(), 'file': stem + '.png',
                         'data_file': s, 'artifact_id': ART.get(s, '')})
    return res

def assemble(figid, rows, gap=90):
    """Strict grid assembler: every panel pasted at native pixel size, left-to-right, rows top-to-bottom.
    Standard panels share one size by construction; wide panels are rendered at exactly k*cell + (k-1)*gap."""
    ims = {L: Image.open(f'figures/panels/{figid}{L}.png').convert('RGB') for r in rows for L in r}
    cw = min(im.width for im in ims.values()); ch = min(im.height for im in ims.values())
    for L, im in ims.items():
        k = round((im.width + gap)/(cw + gap)); assert abs(im.width - (k*cw + (k-1)*gap)) <= 2, (L, im.size, cw)
        assert im.height == ch, (L, im.size, ch)
    row_w = [sum(ims[L].width for L in r) + gap*(len(r)-1) for r in rows]
    W = max(row_w); H = len(rows)*ch + gap*(len(rows)-1)
    canvas = Image.new('RGB', (W, H), 'white'); y = 0
    for r in rows:
        x = 0
        for L in r:
            canvas.paste(ims[L], (x, y)); x += ims[L].width + gap
        y += ch + gap
    os.makedirs('figures', exist_ok=True)
    canvas.save(f'figures/{figid}.png', dpi=(DPI, DPI)); canvas.save(f'figures/{figid}.pdf', 'PDF', resolution=float(DPI))
    return canvas.size

def write_manifest():
    pd.DataFrame(MANIFEST).drop_duplicates().to_csv('figure_manifest.csv', index=False)

# ---- KM helper (step curves from pre-computed survival tables) ----
def km_panel(fig, curves, colors, styles=None, labels=None, xmax_years=None, risk_times=(0,2,4,6,8,10),
             ylabel='Overall survival', annot=None, legend_loc='upper right', outer=None):
    """curves: dict name -> DataFrame(time_days, surv, n_risk). Draws step curves in the upper
    part of RECT and an n-at-risk table in the lower part. Returns (ax_main, ax_tab)."""
    x0, y0, w, h = outer or RECT
    n = len(curves)
    tab_h = 0.055 + 0.045*n
    ax = fig.add_axes([x0, y0 + tab_h + 0.06, w, h - tab_h - 0.06])
    axt = fig.add_axes([x0, y0, w, tab_h]); axt.axis('off')
    tmax = 0
    for i, (name, df) in enumerate(curves.items()):
        t = df['time_days'].values/365.25; s = df['surv'].values
        ax.step(t, s, where='post', color=colors[name], lw=1.3,
                ls=(styles or {}).get(name, '-'), label=(labels or {}).get(name, name))
        tmax = max(tmax, t.max())
    xmax = xmax_years or float(np.ceil(tmax))
    rts = [r for r in risk_times if r <= xmax]
    ax.set_xlim(0, xmax); ax.set_ylim(0, 1.02); ax.set_xticks(rts)
    ax.set_ylabel(ylabel); ax.set_xlabel('Time (years)')
    ax.legend(loc=legend_loc, fontsize=6, handlelength=1.6)
    if annot:
        ax.text(0.03, 0.06, annot, transform=ax.transAxes, fontsize=6, va='bottom', ha='left')
    # n at risk
    axt.set_xlim(0, xmax); axt.set_ylim(-0.5, n - 0.5 + 0.9)
    axt.text(-0.02*xmax, n - 0.5 + 0.45, 'No. at risk', fontsize=5.5, ha='right', va='center', color=DARK)
    for i, (name, df) in enumerate(curves.items()):
        row = n - 1 - i
        axt.text(-0.02*xmax, row, (labels or {}).get(name, name), fontsize=5.5, ha='right', va='center', color=colors[name])
        tt = df['time_days'].values/365.25; nr = df['n_risk'].values
        for r in rts:
            idx = np.searchsorted(tt, r, side='right') - 1
            val = int(nr[max(idx, 0)]) if r <= tt.max() else 0
            axt.text(r, row, str(val), fontsize=5.5, ha='center', va='center')
    return ax, axt

def fmt_hr(hr, lo, hi, p):
    return f'HR {hr:.2f} ({lo:.2f}-{hi:.2f})\n' + fmt_p(p).replace('$', '')

def pstr(p):
    """Plain-text P string (no mathtext)."""
    if p < 1e-4: 
        e = int(np.floor(np.log10(p))); m = p/10**e
        return f'P = {m:.1f}x10^{e}'
    return f'P = {p:.3f}' if p >= 0.001 else f'P = {p:.1e}'
