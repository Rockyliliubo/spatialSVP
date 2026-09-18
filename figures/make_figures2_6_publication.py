from __future__ import annotations

import math
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import FancyBboxPatch, Rectangle, Circle
from matplotlib.lines import Line2D
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / 'analysis' / 'results'
OUT = Path(__file__).resolve().parent
OUT.mkdir(exist_ok=True)
SOURCE_DATA = OUT / 'source_data'
SOURCE_DATA.mkdir(exist_ok=True)

mpl.rcParams.update({
    'font.family': 'Arial',
    'font.size': 8.0,
    'axes.titlesize': 9.2,
    'axes.labelsize': 8.2,
    'xtick.labelsize': 7.2,
    'ytick.labelsize': 7.2,
    'legend.fontsize': 7.2,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
    'svg.fonttype': 'none',
    'axes.linewidth': 0.7,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.04,
})

COL = {
    'ink': '#24323D',
    'muted': '#657682',
    'grid': '#DDE5EA',
    'line': '#C8D2D9',
    'panel': '#F8FAFB',
    'white': '#FFFFFF',
    'grey': '#7B8790',
    'grey_mid': '#A9B4BB',
    'grey_light': '#EEF2F4',
    'absolute': '#D55E00',
    'absolute_light': '#F7E4D8',
    'competitive': '#009E73',
    'competitive_light': '#DDF3EC',
    'reference': '#0072B2',
    'reference_light': '#DDECF5',
    'pathway': '#E69F00',
    'pathway_light': '#F9EDCF',
    'pc': '#CC79A7',
    'pc_light': '#F2E0EC',
    'lib': '#009E73',
    'lib_light': '#DDF3EC',
    'primary': '#1F4E79',
    'loss': '#D55E00',
    'gain': '#0072B2',
    'retain': '#009E73',
    'breast': '#0072B2',
    'mouse': '#E69F00',
    'lymph': '#009E73',
}

SEQ_BLUE = LinearSegmentedColormap.from_list('seqblue', ['#F7F7F7', '#D8E3F0', '#86AAD3', '#0072B2'])


def mm(x: float) -> float:
    return x / 25.4


def style_axis(ax, grid_axis='x'):
    ax.set_facecolor(COL['white'])
    for s in ['top', 'right']:
        ax.spines[s].set_visible(False)
    ax.spines['left'].set_color(COL['line'])
    ax.spines['bottom'].set_color(COL['line'])
    ax.tick_params(color=COL['line'], labelcolor=COL['ink'], width=0.65, length=3)
    if grid_axis:
        ax.grid(True, axis=grid_axis, color=COL['grid'], linewidth=0.55, zorder=0)
    ax.set_axisbelow(True)


def panel_label(ax, letter: str, title: str, subtitle: str | None = None, y: float = 1.035, sub_y: float = 0.99):
    ax.text(-0.075, y, letter, transform=ax.transAxes, ha='left', va='bottom',
            fontsize=11.5, fontweight='bold', color=COL['ink'])
    ax.text(0.025, y, title, transform=ax.transAxes, ha='left', va='bottom',
            fontsize=9.5, fontweight='bold', color=COL['ink'])
    if subtitle:
        ax.text(0.025, sub_y, subtitle, transform=ax.transAxes, ha='left', va='bottom',
                fontsize=7.0, color=COL['muted'])


def rounded_box(ax, xy, width, height, facecolor, edgecolor=None, lw=0.8, radius=0.02, transform=None, zorder=1):
    if transform is None:
        transform = ax.transAxes
    p = FancyBboxPatch(xy, width, height,
                       boxstyle=f'round,pad=0.008,rounding_size={radius}',
                       transform=transform, facecolor=facecolor,
                       edgecolor=edgecolor or 'none', linewidth=lw, zorder=zorder)
    ax.add_patch(p)
    return p


def save_all(fig, base: Path, png_dpi=300, tiff_dpi=600):
    base.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(base.with_suffix('.pdf'), facecolor='white')
    fig.savefig(base.with_suffix('.svg'), facecolor='white')
    fig.savefig(base.with_suffix('.png'), dpi=png_dpi, facecolor='white')
    tmp = base.with_name(base.name + '_tmp600.png')
    fig.savefig(tmp, dpi=tiff_dpi, facecolor='white')
    im = Image.open(tmp).convert('RGB')
    im.save(base.with_name(base.name + '_600dpi.tiff'), compression='tiff_lzw', dpi=(tiff_dpi, tiff_dpi))
    tmp.unlink(missing_ok=True)


def wilson(k: int, n: int, z: float = 1.959963984540054) -> Tuple[float, float]:
    p = k / n
    den = 1 + z*z/n
    center = (p + z*z/(2*n)) / den
    half = z * math.sqrt((p*(1-p)/n) + z*z/(4*n*n)) / den
    return center-half, center+half


# -----------------------------------------------------------------------------
# Figure 2 data
# -----------------------------------------------------------------------------

def build_fig2_data() -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    labels = {
        'SPARKX_ACAT': 'SPARK-X / ACAT',
        'STpathway_dCov': 'STpathway dCov',
        'SPARKX_fgsea': 'SPARK-X / fgsea',
        'spatialSVP_permutation_none': 'spatialSVP permutation (total)',
        'spatialSVP_permutation_pcs': 'spatialSVP permutation (10 PC)',
        'spatialSVP_matched_none': 'spatialSVP matched (total)',
        'spatialSVP_matched_pcs': 'spatialSVP matched (10 PC)',
    }
    endpoint = {k: ('absolute' if ('permutation' in k or k in {'SPARKX_ACAT', 'STpathway_dCov'}) else 'competitive')
                for k in labels}
    tier = {k: ('10 PC' if k.endswith('_pcs') else ('total' if k.startswith('spatialSVP') else 'external'))
            for k in labels}
    bench = pd.read_csv(RESULTS / 'method_benchmark' / 'benchmark_summary.csv')
    a0 = bench.loc[bench.scenario.eq('global_null')].copy()
    a = pd.DataFrame({
        'method': a0.method.map(labels), 'endpoint': a0.method.map(endpoint),
        'tier': a0.method.map(tier),
        'rejections': np.rint(a0.any_rejection * a0.n_replicates).astype(int),
        'replicates': a0.n_replicates.astype(int), 'rate': a0.any_rejection,
        'lower': a0.any_rejection_lo, 'upper': a0.any_rejection_hi,
    })
    scenarios = {'coherent_weak': 'coherent weak', 'coherent_strong': 'coherent strong',
                 'opposing_members': 'opposing members'}
    keep = ['SPARKX_ACAT', 'STpathway_dCov', 'SPARKX_fgsea',
            'spatialSVP_matched_none', 'spatialSVP_matched_pcs']
    b0 = bench.loc[bench.scenario.isin(scenarios) & bench.method.isin(keep)].copy()
    b = pd.DataFrame({
        'method': b0.method.map(labels), 'endpoint': b0.method.map(endpoint),
        'scenario': b0.scenario.map(scenarios),
        'detected': np.rint(b0.planted_detection * b0.n_replicates * 3).astype(int),
        'planted': (b0.n_replicates * 3).astype(int), 'rate': b0.planted_detection,
        'lower': b0.planted_detection_lo, 'upper': b0.planted_detection_hi,
    })
    c0 = pd.read_csv(RESULTS / 'direct_method_comparison' / 'GSM8452847_Pt-1A' / 'call_counts.csv')
    c = pd.DataFrame({'method': c0.method.map(labels), 'endpoint': c0.method.map(endpoint),
                      'called': c0.n_BH_lt_0_05.astype(int)})
    if a.method.isna().any() or b.method.isna().any() or c.method.isna().any():
        raise ValueError('Unmapped benchmark method label')
    return a, b, c


def marker_for(method: str, tier: str):
    if 'SPARK-X' in method and 'ACAT' in method:
        return 'D'
    if 'fgsea' in method:
        return '^'
    if 'STpathway' in method:
        return 's'
    return 'o'


def draw_fig2(a: pd.DataFrame, b: pd.DataFrame, c: pd.DataFrame):
    fig = plt.figure(figsize=(mm(180), mm(210)), facecolor='white')
    gs = fig.add_gridspec(3, 1, height_ratios=[1.00, 1.34, 0.96],
                          left=0.145, right=0.985, top=0.94, bottom=0.06, hspace=0.72)

    # A calibration
    ax = fig.add_subplot(gs[0, 0])
    panel_label(ax, 'A', 'Global-null calibration',
                'Replicates with at least one BH-adjusted rejection; Wilson 95% CI', y=1.16, sub_y=1.085)
    order = a.method.tolist()[::-1]
    y = np.arange(len(order))
    ax.axvspan(0, 0.05, color=COL['grey_light'], zorder=0)
    ax.axvline(0.05, color=COL['grey'], ls=(0, (4, 3)), lw=1.0, zorder=1)
    for yi, method in zip(y, order):
        r = a[a.method == method].iloc[0]
        color = COL['absolute'] if r.endpoint == 'absolute' else COL['competitive']
        marker = marker_for(method, r.tier)
        mfc = color if r.tier != '10 PC' else 'white'
        ax.errorbar(r.rate, yi, xerr=[[r.rate-r.lower], [r.upper-r.rate]], fmt=marker,
                    ms=5.2, mfc=mfc, mec=color, mew=1.0, ecolor=color,
                    elinewidth=1.0, capsize=2.4, zorder=3)
        ax.text(min(r.upper+0.006, 0.128), yi, f'{100*r.rate:.0f}%', va='center', ha='left',
                fontsize=6.6, color=COL['ink'])
    ax.set_yticks(y, order)
    for lab in ax.get_yticklabels():
        if lab.get_text() == 'spatialSVP matched (total)':
            lab.set_fontweight('bold'); lab.set_color(COL['competitive'])
    ax.set_xlim(-0.002, 0.135)
    ax.set_xticks([0, 0.05, 0.10], ['0%', '5%', '10%'])
    ax.set_xlabel('Replicates with any rejection')
    style_axis(ax, 'x')
    ax.text(0.995, 1.01, 'dashed line: 5%', transform=ax.transAxes,
            ha='right', va='bottom', fontsize=6.0, color=COL['muted'])

    # B power small multiples
    methods_b = ['SPARK-X / ACAT', 'STpathway dCov', 'SPARK-X / fgsea',
                 'spatialSVP matched (total)', 'spatialSVP matched (10 PC)']
    scenarios = ['coherent weak', 'coherent strong', 'opposing members']
    sub = gs[1, 0].subgridspec(1, 3, wspace=0.15)
    axs = []
    for j, sc in enumerate(scenarios):
        axb = fig.add_subplot(sub[0, j])
        axs.append(axb)
        if j == 0:
            panel_label(axb, 'B', 'Detection of planted pathway programs',
                        'Three signal architectures; 300 planted sets per scenario', y=1.26, sub_y=1.18)
        axb.set_title(sc, fontsize=8.0, fontweight='bold', color=COL['ink'], pad=4,
                      bbox=dict(facecolor=COL['grey_light'], edgecolor='none', boxstyle='round,pad=0.22'))
        ypos = np.arange(len(methods_b))[::-1]
        for yi, method in zip(ypos, methods_b):
            r = b[(b.method == method) & (b.scenario == sc)].iloc[0]
            color = COL['absolute'] if r.endpoint == 'absolute' else COL['competitive']
            marker = marker_for(method, '10 PC' if '10 PC' in method else ('total' if 'spatialSVP' in method else 'external'))
            mfc = 'white' if '10 PC' in method else color
            axb.errorbar(r.rate*100, yi,
                         xerr=[[max(0.0,(r.rate-r.lower)*100)], [max(0.0,(r.upper-r.rate)*100)]],
                         fmt=marker, ms=5.1, mfc=mfc, mec=color, mew=1.0, ecolor=color,
                         elinewidth=0.9, capsize=2.0, zorder=3)
            val = 100*r.rate
            txt = f'{val:.1f}%' if abs(val-round(val))>0.01 else f'{val:.0f}%'
            if val >= 92:
                axb.text(val-2.0, yi+0.16, txt, ha='right', va='bottom', fontsize=5.5, color=COL['ink'])
            else:
                axb.text(val+2.2, yi+0.02, txt, ha='left', va='center', fontsize=5.5, color=COL['ink'])
        axb.set_xlim(-5, 110); axb.set_ylim(-0.75, 4.25)
        axb.set_xticks([0,25,50,75,100], ['0','25','50','75','100'])
        axb.set_xlabel('Detected (%)' if j==1 else '')
        axb.set_yticks(ypos, methods_b if j==0 else ['']*len(methods_b))
        if j==0:
            for lab in axb.get_yticklabels():
                if lab.get_text() == 'spatialSVP matched (total)':
                    lab.set_fontweight('bold'); lab.set_color(COL['competitive'])
        style_axis(axb, 'x')
    # Compact callout below the lowest row, avoiding all data points
    axo = axs[2]
    rounded_box(axo, (0.04, 0.01), 0.60, 0.105, COL['pathway_light'], edgecolor='#E9D49B', lw=0.7, radius=0.02, zorder=1)
    axo.text(0.34, 0.062, 'Mean-score cancellation:\nprimary test detects 1/300',
             transform=axo.transAxes, ha='center', va='center', fontsize=6.1, color=COL['ink'])

    # C Pt-1A calls
    axc = fig.add_subplot(gs[2, 0])
    panel_label(axc, 'C', 'Direct comparison on Pt-1A',
                'Same 13,708 genes, 2,746 spots, and 50 Hallmark pathways', y=1.15, sub_y=1.075)
    order_c = a.method.tolist()[::-1]
    ypos = np.arange(len(order_c))
    for yi, method in zip(ypos, order_c):
        r = c[c.method == method].iloc[0]
        color = COL['absolute'] if r.endpoint == 'absolute' else COL['competitive']
        alpha = 1.0 if method == 'spatialSVP matched (total)' else 0.82
        edge = COL['ink'] if method == 'spatialSVP matched (total)' else 'none'
        axc.barh(yi, r.called, color=color, alpha=alpha, height=0.62,
                 edgecolor=edge, linewidth=0.8, zorder=2)
        axc.text(r.called + 0.6, yi, f'{int(r.called)}', va='center', ha='left', fontsize=6.9,
                 fontweight='bold' if method == 'spatialSVP matched (total)' else 'normal', color=COL['ink'])
    axc.set_yticks(ypos, order_c)
    for lab in axc.get_yticklabels():
        if lab.get_text() == 'spatialSVP matched (total)':
            lab.set_fontweight('bold'); lab.set_color(COL['competitive'])
    axc.set_xlim(0, 53); axc.set_xticks([0,10,20,30,40,50])
    axc.set_xlabel('Hallmark pathways called (of 50)')
    style_axis(axc, 'x')
    axc.text(0.995, 0.02, 'Different endpoints - not an accuracy ranking', transform=axc.transAxes,
             ha='right', va='bottom', fontsize=6.6, color=COL['muted'], fontstyle='italic')

    save_all(fig, OUT/'fig2_method_benchmark')
    plt.close(fig)


# -----------------------------------------------------------------------------
# Figure 3 heatmap extraction and plotting
# -----------------------------------------------------------------------------

PATHWAYS_TOP_TO_BOTTOM = [
    'TNFA_SIGNALING_VIA_NFKB','MYOGENESIS','MTORC1_SIGNALING','INTERFERON_GAMMA_RESPONSE',
    'INTERFERON_ALPHA_RESPONSE','HYPOXIA','G2M_CHECKPOINT','EPITHELIAL_MESENCHYMAL_TRANSITION',
    'E2F_TARGETS','COAGULATION','APICAL_JUNCTION','ALLOGRAFT_REJECTION','INFLAMMATORY_RESPONSE',
    'MITOTIC_SPINDLE','IL6_JAK_STAT3_SIGNALING','ANGIOGENESIS','MYC_TARGETS_V1','KRAS_SIGNALING_UP',
    'COMPLEMENT','UV_RESPONSE_DN','PANCREAS_BETA_CELLS','GLYCOLYSIS','APOPTOSIS','OXIDATIVE_PHOSPHORYLATION',
    'CHOLESTEROL_HOMEOSTASIS','ESTROGEN_RESPONSE_LATE','ESTROGEN_RESPONSE_EARLY','TGF_BETA_SIGNALING',
    'MYC_TARGETS_V2','KRAS_SIGNALING_DN','IL2_STAT5_SIGNALING','XENOBIOTIC_METABOLISM',
    'FATTY_ACID_METABOLISM','BILE_ACID_METABOLISM','P53_PATHWAY','UV_RESPONSE_UP',
    'UNFOLDED_PROTEIN_RESPONSE','ADIPOGENESIS','PEROXISOME','PROTEIN_SECRETION','SPERMATOGENESIS',
    'PI3K_AKT_MTOR_SIGNALING','NOTCH_SIGNALING','APICAL_SURFACE','ANDROGEN_RESPONSE',
    'HEDGEHOG_SIGNALING','DNA_REPAIR','HEME_METABOLISM','REACTIVE_OXYGEN_SPECIES_PATHWAY',
    'WNT_BETA_CATENIN_SIGNALING'
]

DISPLAY = {
    'TNFA_SIGNALING_VIA_NFKB':'TNFα signaling via NF-κB',
    'MTORC1_SIGNALING':'mTORC1 signaling',
    'INTERFERON_GAMMA_RESPONSE':'Interferon-γ response',
    'INTERFERON_ALPHA_RESPONSE':'Interferon-α response',
    'G2M_CHECKPOINT':'G2/M checkpoint',
    'EPITHELIAL_MESENCHYMAL_TRANSITION':'Epithelial–mesenchymal transition',
    'IL6_JAK_STAT3_SIGNALING':'IL-6/JAK/STAT3 signaling',
    'KRAS_SIGNALING_UP':'KRAS signaling up',
    'KRAS_SIGNALING_DN':'KRAS signaling down',
    'UV_RESPONSE_DN':'UV response down',
    'UV_RESPONSE_UP':'UV response up',
    'PANCREAS_BETA_CELLS':'Pancreas β cells',
    'OXIDATIVE_PHOSPHORYLATION':'Oxidative phosphorylation',
    'TGF_BETA_SIGNALING':'TGF-β signaling',
    'IL2_STAT5_SIGNALING':'IL-2/STAT5 signaling',
    'P53_PATHWAY':'p53 pathway',
    'PI3K_AKT_MTOR_SIGNALING':'PI3K–AKT–mTOR signaling',
    'REACTIVE_OXYGEN_SPECIES_PATHWAY':'Reactive oxygen species pathway',
    'WNT_BETA_CATENIN_SIGNALING':'WNT/β-catenin signaling',
}

def display_pathway(x: str) -> str:
    if x in DISPLAY:
        return DISPLAY[x]
    return x.replace('_',' ').lower().capitalize().replace('Myc','MYC').replace('E2f','E2F').replace('Dna','DNA')


def load_fig3_matrix() -> pd.DataFrame:
    rec = pd.read_csv(RESULTS / 'cohort_residualization_v1' / 'recurrence_by_cohort.csv')
    rec['pathway'] = rec.pathway.str.replace('^HALLMARK_', '', regex=True)
    tier_names = {'none': 'total', 'libsize': 'library-size residual', 'pcs': '10-PC residual'}
    rec['estimand'] = rec.tier.map(tier_names)
    cols = pd.MultiIndex.from_tuples([
        ('total','GSE274557'), ('total','GSE282302'),
        ('library-size residual','GSE274557'), ('library-size residual','GSE282302'),
        ('10-PC residual','GSE274557'), ('10-PC residual','GSE282302')
    ], names=['estimand','cohort'])
    matrix = rec.pivot(index='pathway', columns=['estimand', 'cohort'], values='positive')
    matrix = matrix.reindex(index=PATHWAYS_TOP_TO_BOTTOM, columns=cols)
    if matrix.shape != (50, 6) or matrix.isna().any().any():
        raise ValueError(f'Unexpected recurrence matrix: {matrix.shape}, missing={int(matrix.isna().sum().sum())}')
    return matrix


def draw_fig3(df: pd.DataFrame):
    fig = plt.figure(figsize=(mm(180), mm(225)), facecolor='white')

    # A: heatmap (large left margin reserved for pathway names)
    ax = fig.add_axes([0.285, 0.255, 0.650, 0.685])
    panel_label(ax, 'A', 'Patient-positive fractions across 50 Hallmark pathways',
                'Outlined row-pairs meet the fixed ≥80% recurrence rule in both cohorts', y=1.075, sub_y=1.025)
    data = df.values
    im = ax.imshow(data, aspect='auto', cmap=SEQ_BLUE, vmin=0, vmax=1, interpolation='nearest')
    for x in np.arange(-0.5, 6, 1):
        ax.axvline(x, color='white', lw=0.8)
    for y in np.arange(-0.5, 50, 1):
        ax.axhline(y, color='white', lw=0.35)
    ax.axvline(1.5, color=COL['ink'], lw=1.4)
    ax.axvline(3.5, color=COL['ink'], lw=1.4)
    ax.set_yticks(np.arange(50), [display_pathway(x) for x in df.index])
    ax.tick_params(axis='y', length=0, pad=3, labelsize=5.85)
    xt = ['Pei\nGSE274557\nn = 13', 'Lyubetskaya\nGSE282302\nn = 39'] * 3
    ax.set_xticks(np.arange(6), xt)
    ax.tick_params(axis='x', length=0, pad=4, labelsize=6.2)

    headers = [('Total organization', COL['primary'], 0),
               ('Library-size residual', COL['lib'], 2),
               ('10-PC residual', COL['pc'], 4)]
    for title, color, start_col in headers:
        rec = (data[:,start_col] >= 0.8) & (data[:,start_col+1] >= 0.8)
        ax.add_patch(Rectangle((start_col-0.5, -3.95), 2, 2.15, transform=ax.transData,
                               facecolor=color, edgecolor='none', clip_on=False, alpha=0.97))
        ax.text(start_col+0.5, -3.15, title, ha='center', va='center', fontsize=7.0,
                color='white', fontweight='bold', clip_on=False)
        ax.text(start_col+0.5, -2.25, f'{int(rec.sum())}/50 recurrent', ha='center', va='center', fontsize=6.2,
                color='white', clip_on=False)
        for i in np.where(rec)[0]:
            ax.add_patch(Rectangle((start_col-0.49, i-0.48), 1.98, 0.96, fill=False,
                                   edgecolor=color, linewidth=1.35, zorder=5))
    ax.set_xlim(-0.5,5.5); ax.set_ylim(49.5,-4.25)
    for s in ax.spines.values(): s.set_visible(False)

    cax = fig.add_axes([0.947, 0.405, 0.012, 0.31])
    cb = fig.colorbar(im, cax=cax)
    cb.set_ticks([0,0.25,0.5,0.75,1.0])
    cb.ax.tick_params(labelsize=6, length=2, color=COL['line'])
    cb.outline.set_edgecolor(COL['line'])
    cb.set_label('Patient-positive fraction', fontsize=6.5, color=COL['ink'])

    # B: full-width summary
    axb = fig.add_axes([0.05, 0.025, 0.90, 0.165])
    axb.set_xlim(0,1); axb.set_ylim(0,1); axb.axis('off')
    panel_label(axb, 'B', 'Recurrence summary and overlap',
                'Residual tiers are sensitivity analyses, not substitutes for the primary total-map estimand', y=0.96, sub_y=0.86)

    rounded_box(axb, (0.01, 0.12), 0.15, 0.72, COL['grey_light'], edgecolor=COL['line'], radius=0.025)
    axb.text(0.085, 0.67, 'TOTAL', ha='center', va='center', fontsize=7.0, fontweight='bold', color=COL['primary'])
    axb.text(0.085, 0.44, '0', ha='center', va='center', fontsize=25, fontweight='bold', color=COL['primary'])
    axb.text(0.085, 0.22, 'recurrent pathways', ha='center', va='center', fontsize=6.3, color=COL['muted'])

    # Two-set overlap for residual estimands
    axb.add_patch(Circle((0.31,0.47), 0.14, transform=axb.transAxes,
                         facecolor=COL['lib_light'], edgecolor=COL['lib'], lw=1.4))
    axb.add_patch(Circle((0.42,0.47), 0.14, transform=axb.transAxes,
                         facecolor=COL['pc_light'], edgecolor=COL['pc'], lw=1.4, alpha=0.88))
    axb.text(0.27,0.47,'4',ha='center',va='center',fontsize=16,fontweight='bold',color=COL['lib'])
    axb.text(0.365,0.47,'10',ha='center',va='center',fontsize=18,fontweight='bold',color=COL['ink'])
    axb.text(0.46,0.47,'1',ha='center',va='center',fontsize=16,fontweight='bold',color=COL['pc'])
    axb.text(0.255,0.76,'Library-size residual\n14 pathways',ha='center',va='center',fontsize=6.4,color=COL['lib'],fontweight='bold')
    axb.text(0.475,0.76,'10-PC residual\n11 pathways',ha='center',va='center',fontsize=6.4,color=COL['pc'],fontweight='bold')

    shared = ['Allograft rejection','Coagulation','E2F targets','EMT','G2/M checkpoint',
              'Hypoxia','Inflammatory response','Interferon-α response','Interferon-γ response','TNFα/NF-κB']
    rounded_box(axb, (0.57, 0.12), 0.42, 0.72, COL['panel'], edgecolor=COL['line'], radius=0.025)
    axb.text(0.78,0.72,'10 pathways shared by both residual tiers',ha='center',va='center',
             fontsize=6.8,fontweight='bold',color=COL['ink'])
    left='\n'.join('• '+x for x in shared[:5]); right='\n'.join('• '+x for x in shared[5:])
    axb.text(0.60,0.60,left,ha='left',va='top',fontsize=5.8,color=COL['ink'],linespacing=1.20)
    axb.text(0.79,0.60,right,ha='left',va='top',fontsize=5.8,color=COL['ink'],linespacing=1.20)
    axb.text(0.5, 0.005,
             'Key qualification: recurrence after residualization is recurrence of residual spatial organization, not of the unadjusted activity map.',
             ha='center', va='bottom', fontsize=6.2, color=COL['muted'], fontstyle='italic')

    save_all(fig, OUT/'fig3_cohort_estimands')
    plt.close(fig)


# -----------------------------------------------------------------------------
# Figure 4 exact scatter extraction and plotting
# -----------------------------------------------------------------------------

def load_fig4_points() -> pd.DataFrame:
    data = pd.read_csv(RESULTS / 'composition_primary_none_v1' / 'results_long.csv')
    specifications = [
        ('full_section', 'marker', 'Marker modules'),
        ('full_section', 'label_transfer', 'Label transfer'),
        ('RCTD_intersection', 'RCTD', 'RCTD'),
    ]
    frames = []
    for analysis, adjusted_arm, label in specifications:
        raw = data.loc[(data.analysis == analysis) & (data.arm == 'raw'),
                       ['section_id', 'pathway', 'p_matched']].rename(columns={'p_matched': 'p_raw'})
        adjusted = data.loc[(data.analysis == analysis) & (data.arm == adjusted_arm),
                            ['section_id', 'pathway', 'p_matched']].rename(columns={'p_matched': 'p_adjusted'})
        paired = raw.merge(adjusted, on=['section_id', 'pathway'], validate='one_to_one')
        paired['adjustment'] = label
        frames.append(paired)
    paired = pd.concat(frames, ignore_index=True)
    paired['neglog10_p_raw'] = -np.log10(paired.p_raw.clip(lower=np.finfo(float).tiny))
    paired['neglog10_p_adjusted'] = -np.log10(paired.p_adjusted.clip(lower=np.finfo(float).tiny))
    if len(paired) != 900:
        raise ValueError(f'Expected 900 paired section-pathway results, found {len(paired)}')
    return paired[['adjustment', 'neglog10_p_raw', 'neglog10_p_adjusted']]


def draw_fig4(points: pd.DataFrame):
    raw_summary = pd.read_csv(RESULTS / 'composition_primary_none_v1' / 'adjustment_summary.csv')
    label_map = {'marker_modules': 'Marker modules', 'label_transfer': 'Label transfer', 'RCTD': 'RCTD'}
    raw_summary['adjustment'] = raw_summary.adjustment.map(label_map)
    summary = raw_summary.rename(columns={
        'retained_raw_calls': 'retained', 'lost_raw_calls': 'lost',
        'gained_calls': 'gained', 'adjusted_calls': 'adjusted_total',
    })[['adjustment', 'retained', 'lost', 'gained', 'adjusted_total']]
    summary['comparisons'] = raw_summary.comparisons.astype(int)
    fig = plt.figure(figsize=(mm(180), mm(142)), facecolor='white')
    gs = fig.add_gridspec(2, 3, height_ratios=[2.0, 1.05], left=0.09, right=0.985,
                          top=0.84, bottom=0.12, hspace=0.52, wspace=0.18)

    names=['Marker modules','Label transfer','RCTD']
    for j,name in enumerate(names):
        ax=fig.add_subplot(gs[0,j])
        if j==0:
            panel_label(ax,'A','Matched-null p-value sensitivity',
                        'Above diagonal = stronger after adjustment', y=1.38, sub_y=1.28)
        sub=points[points.adjustment==name].copy()
        delta=sub.neglog10_p_adjusted-sub.neglog10_p_raw
        colors=np.where(delta>0.15,COL['competitive'],np.where(delta<-0.15,COL['absolute'],COL['grey_mid']))
        ax.scatter(sub.neglog10_p_raw,sub.neglog10_p_adjusted,s=9,c=colors,alpha=0.68,
                   edgecolors='white',linewidths=0.18,zorder=3)
        ax.plot([0,3.35],[0,3.35],ls=(0,(4,3)),lw=0.9,color=COL['grey'],zorder=2)
        ax.set_xlim(0,3.35); ax.set_ylim(0,3.35)
        ax.set_xticks([0,1,2,3]); ax.set_yticks([0,1,2,3])
        ax.set_xlabel(r'$-\log_{10}(p_{matched})$ raw')
        if j==0: ax.set_ylabel(r'$-\log_{10}(p_{matched})$ composition-adjusted')
        else: ax.set_yticklabels([])
        style_axis(ax, 'both')
        ax.set_title(name, fontsize=8.3, fontweight='bold', color=COL['ink'], pad=7)
        ss=summary[summary.adjustment==name].iloc[0]
        ax.text(0.03,0.95,f'{int(ss.comparisons)} paired tests',transform=ax.transAxes,
                ha='left',va='top',fontsize=6.2,color=COL['muted'],
                bbox=dict(boxstyle='round,pad=0.24',facecolor=COL['panel'],edgecolor=COL['line'],linewidth=0.6))
    handles=[Line2D([0],[0],marker='o',ls='none',mfc=COL['competitive'],mec='none',label='stronger after adjustment'),
             Line2D([0],[0],marker='o',ls='none',mfc=COL['absolute'],mec='none',label='weaker after adjustment'),
             Line2D([0],[0],marker='o',ls='none',mfc=COL['grey_mid'],mec='none',label='similar')]
    fig.legend(handles=handles,loc='upper right',bbox_to_anchor=(0.985,0.925),frameon=False,
               ncol=3,handletextpad=0.3,columnspacing=0.9)

    axb=fig.add_subplot(gs[1,:])
    panel_label(axb,'B','BH-significant call reclassification',
                'Raw calls = 82; adjusted totals combine retained and gained calls', y=1.18, sub_y=1.09)
    x=np.arange(3); w=0.28
    axb.bar(x-w/1.8,summary.retained,width=w,color=COL['retain'],label='retained',zorder=3)
    axb.bar(x-w/1.8,summary.lost,bottom=summary.retained,width=w,color=COL['loss'],label='lost after adjustment',zorder=3)
    axb.bar(x+w/1.8,summary.retained,width=w,color=COL['retain'],zorder=3)
    axb.bar(x+w/1.8,summary.gained,bottom=summary.retained,width=w,color=COL['gain'],label='gained after adjustment',zorder=3)
    for i,row in summary.iterrows():
        axb.text(i-w/1.8,84.0,'RAW\n82',ha='center',va='bottom',fontsize=6.0,color=COL['muted'])
        axb.text(i+w/1.8,row.adjusted_total+1.8,f'ADJ\n{row.adjusted_total}',ha='center',va='bottom',fontsize=6.0,color=COL['muted'])
        axb.text(i-w/1.8,row.retained/2,f'{row.retained}',ha='center',va='center',fontsize=6.2,color='white',fontweight='bold')
        axb.text(i-w/1.8,row.retained+row.lost/2,f'{row.lost}',ha='center',va='center',fontsize=6.2,color='white',fontweight='bold')
        axb.text(i+w/1.8,row.retained/2,f'{row.retained}',ha='center',va='center',fontsize=6.2,color='white',fontweight='bold')
        axb.text(i+w/1.8,row.retained+row.gained/2,f'{row.gained}',ha='center',va='center',fontsize=6.2,color='white',fontweight='bold')
    axb.set_xticks(x,summary.adjustment)
    axb.set_ylim(0,100)
    axb.set_ylabel('Significant section-pathway calls')
    style_axis(axb,'y')
    axb.legend(loc='upper right',bbox_to_anchor=(1.0,1.18),frameon=False,ncol=3,
               handletextpad=0.4,columnspacing=1.0)
    axb.text(0.5,-0.28,
             'Composition adjustment changes both the named pathway map and its matched reference maps; losses and gains are therefore both expected.',
             transform=axb.transAxes,ha='center',va='top',fontsize=6.3,color=COL['muted'],fontstyle='italic')

    save_all(fig, OUT/'fig4_composition_primary')
    summary.to_csv(SOURCE_DATA/'fig4_call_reclassification.csv',index=False)
    plt.close(fig)


# -----------------------------------------------------------------------------
# Figure 5: signal removed by residualization
# -----------------------------------------------------------------------------

def load_fig5_data() -> pd.DataFrame:
    data = pd.read_csv(RESULTS / 'residualization_diagnostic' / 'paired_summary.csv')
    expected = {'coherent_weak', 'coherent_strong'}
    if set(data.scenario) != expected or len(data) != 10:
        raise ValueError('Unexpected residualization diagnostic table')
    return data


def draw_fig5(data: pd.DataFrame):
    fig, axes = plt.subplots(1, 2, figsize=(mm(180), mm(82)), facecolor='white')
    fig.subplots_adjust(left=0.085, right=0.985, top=0.79, bottom=0.22, wspace=0.26)
    tiers = ['none', 'libsize', 'pcs1', 'pcs3', 'pcs10']
    labels = ['Total', 'Library\nsize', '1 PC', '3 PCs', '10 PCs']
    scenarios = [('coherent_weak', 'Weak', COL['reference_light']),
                 ('coherent_strong', 'Strong', COL['absolute'])]
    x = np.arange(len(tiers)); width = 0.36
    for ax, metric, title, ylabel, letter in [
        (axes[0], 'signal_R2_removed', 'Planted-pattern variance removed', 'Variance removed (%)', 'A'),
        (axes[1], 'detection', 'Planted pathways detected', 'Detection (%)', 'B'),
    ]:
        panel_label(ax, letter, title, 'Paired 20-replicate mechanism experiment', y=1.24, sub_y=1.14)
        for j, (scenario, label, color) in enumerate(scenarios):
            values = data.loc[data.scenario.eq(scenario)].set_index('tier').loc[tiers, metric].to_numpy() * 100
            bars = ax.bar(x + (j - 0.5) * width, values, width=width, color=color,
                          edgecolor='white', linewidth=0.5, label=label, zorder=3)
            for bar, value in zip(bars, values):
                if value < 99.5 or metric == 'signal_R2_removed':
                    ax.text(bar.get_x() + bar.get_width()/2, value + 2.2, f'{value:.0f}',
                            ha='center', va='bottom', fontsize=5.8, color=COL['ink'])
        ax.set_xticks(x, labels)
        ax.set_ylim(0, 108); ax.set_yticks([0, 25, 50, 75, 100])
        ax.set_ylabel(ylabel)
        style_axis(ax, 'y')
    axes[0].legend(loc='upper left', bbox_to_anchor=(0.0, 1.04), frameon=False, ncol=2)
    fig.text(0.5, 0.055,
             'Principal-component regression can remove the planted spatial pattern and reduce detection.',
             ha='center', va='bottom', fontsize=6.5, color=COL['muted'], fontstyle='italic')
    save_all(fig, OUT / 'fig5_residualization')
    plt.close(fig)


# -----------------------------------------------------------------------------
# Figure 6: DLPFC donor replication
# -----------------------------------------------------------------------------

def load_fig6_data() -> pd.DataFrame:
    data = pd.read_csv(RESULTS / 'cross_tissue_dlpfc_donors_v1' / 'donor_summary.csv')
    if len(data) != 9 or set(data.donor) != {'Br5292', 'Br5595', 'Br8100'}:
        raise ValueError('Unexpected DLPFC donor summary')
    return data


def draw_fig6(data: pd.DataFrame):
    fig = plt.figure(figsize=(mm(180), mm(132)), facecolor='white')
    outer = fig.add_gridspec(2, 1, height_ratios=[1.25, 0.92], left=0.085, right=0.985,
                             top=0.89, bottom=0.12, hspace=0.52)
    donors = ['Br5292', 'Br5595', 'Br8100']
    tiers = ['none', 'libsize', 'pcs']
    tier_labels = ['Total', 'Library-size\nresidual', '10-PC\nresidual']
    top = outer[0].subgridspec(1, 3, wspace=0.22)
    for j, donor in enumerate(donors):
        ax = fig.add_subplot(top[0, j])
        if j == 0:
            panel_label(ax, 'A', 'Absolute-null calls',
                        'Percentage significant under Westfall–Young correction', y=1.34, sub_y=1.23)
        sub = data.loc[data.donor.eq(donor)].set_index('tier').loc[tiers]
        x = np.arange(3); width = 0.36
        hallmark = 100 * sub.absolute_hallmark_WY.to_numpy() / sub.hallmark_sets.to_numpy()
        random = 100 * sub.absolute_random_WY.to_numpy() / sub.random_sets.to_numpy()
        ax.bar(x - width/2, hallmark, width, color=COL['reference'], label='Hallmark', zorder=3)
        ax.bar(x + width/2, random, width, color=COL['pathway'], label='Random sets', zorder=3)
        ax.set_title(donor, fontsize=8.5, fontweight='bold', color=COL['ink'], pad=5)
        ax.set_xticks(x, tier_labels)
        ax.set_ylim(0, 106); ax.set_yticks([0, 25, 50, 75, 100])
        if j == 0:
            ax.set_ylabel('Significant sets (%)')
            ax.legend(loc='upper right', frameon=False)
        else:
            ax.set_yticklabels([])
        style_axis(ax, 'y')

    axb = fig.add_subplot(outer[1, 0])
    panel_label(axb, 'B', 'Matched-null Hallmark calls',
                'Significant pathways of 50 tested', y=1.26, sub_y=1.15)
    colors = {'Br5292': COL['reference'], 'Br5595': COL['competitive'], 'Br8100': COL['pc']}
    markers = {'Br5292': 'o', 'Br5595': 's', 'Br8100': '^'}
    x = np.arange(3)
    for donor in donors:
        values = data.loc[data.donor.eq(donor)].set_index('tier').loc[tiers, 'relative_hallmark_BH'].to_numpy()
        axb.plot(x, values, marker=markers[donor], color=colors[donor], linewidth=1.8,
                 markersize=5.5, label=donor, zorder=3)
        offsets = {'Br5292': (0, 8), 'Br5595': (-7, -13), 'Br8100': (7, 7)}
        for xi, value in zip(x, values):
            axb.annotate(str(int(value)), (xi, value), xytext=offsets[donor], textcoords='offset points',
                         ha='center', va='bottom', fontsize=6.0, color=colors[donor])
    axb.set_xticks(x, tier_labels)
    axb.set_ylim(-1, 34); axb.set_yticks([0, 10, 20, 30])
    axb.set_ylabel('Hallmark pathways called')
    style_axis(axb, 'y')
    axb.legend(loc='upper right', frameon=False, ncol=3)
    fig.text(0.5, 0.035,
             'The total-map contrast replicates in three independently selected human donors.',
             ha='center', va='bottom', fontsize=6.5, color=COL['muted'], fontstyle='italic')
    save_all(fig, OUT / 'fig6_dlpfc_donors')
    plt.close(fig)


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main():
    # Figure 2
    a,b,c=build_fig2_data()
    a.to_csv(SOURCE_DATA/'fig2_global_null.csv',index=False)
    b.to_csv(SOURCE_DATA/'fig2_planted_signal_detection.csv',index=False)
    c.to_csv(SOURCE_DATA/'fig2_Pt1A_calls.csv',index=False)
    draw_fig2(a,b,c)

    # Figure 3
    heat=load_fig3_matrix()
    heat.to_csv(SOURCE_DATA/'fig3_patient_positive_fraction_matrix.csv')
    draw_fig3(heat)

    # Figure 4
    pts=load_fig4_points()
    pts.to_csv(SOURCE_DATA/'fig4_composition_scatter_points.csv',index=False)
    draw_fig4(pts)

    # Figure 5
    residual=load_fig5_data()
    residual.to_csv(SOURCE_DATA/'fig5_residualization_data.csv',index=False)
    draw_fig5(residual)

    # Figure 6
    brain=load_fig6_data()
    brain.to_csv(SOURCE_DATA/'fig6_dlpfc_donor_data.csv',index=False)
    draw_fig6(brain)

    print(f'Created optimized figures in {OUT}')

if __name__=='__main__':
    main()
