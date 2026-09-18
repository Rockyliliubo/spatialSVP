import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from pathlib import Path as FilePath
from matplotlib.patches import FancyBboxPatch, Rectangle, Circle, FancyArrowPatch, Ellipse, Arc, PathPatch
from matplotlib.path import Path
from matplotlib.colors import LinearSegmentedColormap
from PIL import Image

mpl.rcParams.update({
    'font.family': 'Arial',
    'font.size': 15,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
    'svg.fonttype': 'none',
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.06,
})

COL = {
    'ink': '#24323D',
    'muted': '#667681',
    'line': '#C4CED5',
    'panel': '#F8FAFB',
    'neutral': '#E8EFF3',
    'pathway': '#E69F00',
    'pathway_dark': '#B87800',
    'reference': '#0072B2',
    'reference_light': '#56B4E9',
    'competitive': '#009E73',
    'competitive_light': '#DDF3EC',
    'absolute': '#D55E00',
    'absolute_light': '#F8E6DD',
    'white': '#FFFFFF',
    'grey_light': '#EEF2F4',
}

score_cmap = LinearSegmentedColormap.from_list(
    'blue_white_orange', [COL['reference'], '#F7F8F7', COL['pathway']]
)


def rounded(ax, xy, w, h, fc='white', ec=None, lw=1.2, radius=0.02, alpha=1, zorder=0):
    p = FancyBboxPatch(
        xy, w, h, transform=ax.transAxes,
        boxstyle=f'round,pad=0.006,rounding_size={radius}',
        facecolor=fc, edgecolor=ec or 'none', linewidth=lw, alpha=alpha,
        zorder=zorder
    )
    ax.add_patch(p)
    return p


def arrow(ax, p1, p2, color=None, lw=1.7, mutation=14, style='-|>', zorder=5):
    a = FancyArrowPatch(
        p1, p2, transform=ax.transAxes, arrowstyle=style,
        mutation_scale=mutation, linewidth=lw,
        color=color or COL['muted'], zorder=zorder
    )
    ax.add_patch(a)
    return a


def panel_background(ax):
    ax.set_xlim(0,1); ax.set_ylim(0,1); ax.axis('off')
    rounded(ax, (0.004,0.004),0.992,0.992,fc=COL['panel'],ec=COL['line'],lw=1.0,radius=0.018,zorder=-20)


def panel_header(ax, letter, title):
    ax.text(0.018,0.976,letter,transform=ax.transAxes,ha='left',va='top',fontsize=21,fontweight='bold',color=COL['ink'])
    ax.text(0.078,0.976,title,transform=ax.transAxes,ha='left',va='top',fontsize=16.8,fontweight='bold',color=COL['ink'])


def draw_check(ax, x, y, r=0.017):
    ax.add_patch(Circle((x,y),r,transform=ax.transAxes,facecolor=COL['competitive'],edgecolor='none'))
    # two strokes
    ax.plot([x-r*0.45,x-r*0.08],[y,y-r*0.36],transform=ax.transAxes,color='white',lw=1.8,solid_capstyle='round')
    ax.plot([x-r*0.08,x+r*0.48],[y-r*0.36,y+r*0.42],transform=ax.transAxes,color='white',lw=1.8,solid_capstyle='round')


def make_tissue_points(seed=1,nx=15,ny=11):
    rng=np.random.default_rng(seed)
    pts=[]
    for j in range(ny):
        for i in range(nx):
            x=i+0.5*(j%2); y=j*0.88
            cx,cy=(nx-1)/2,(ny-1)*0.88/2
            ell=((x-cx)/(nx*0.48))**2+((y-cy)/(ny*0.46))**2
            notch=((x-(cx+3.8))/2.6)**2+((y-(cy+2.4))/1.7)**2
            lower=y<(0.10+0.055*(x-cx)**2)
            if ell<1 and notch>1 and not lower:
                pts.append((x,y))
    xy=np.array(pts,float)
    xy[:,0]=(xy[:,0]-xy[:,0].min())/(xy[:,0].max()-xy[:,0].min())
    xy[:,1]=(xy[:,1]-xy[:,1].min())/(xy[:,1].max()-xy[:,1].min())
    xy += rng.normal(0,0.004,xy.shape)
    return xy

TISSUE_XY=make_tissue_points()


def tissue_values(kind='pathway',seed=1):
    x,y=TISSUE_XY[:,0],TISSUE_XY[:,1]
    rng=np.random.default_rng(seed)
    if kind=='pathway':
        z=1.15*np.exp(-((x-0.72)**2/(2*0.17**2)+(y-0.64)**2/(2*0.20**2))) \
          -0.85*np.exp(-((x-0.24)**2/(2*0.20**2)+(y-0.30)**2/(2*0.22**2)))+0.25*(x-y)
    elif kind=='random':
        z=0.95*np.sin(2.4*np.pi*(x*0.62+y*0.22))+0.50*np.cos(2.0*np.pi*y)+0.18*x
    elif kind=='random2':
        z=1.10*np.exp(-((x-0.35)**2/(2*0.22**2)+(y-0.72)**2/(2*0.18**2)))-0.65*(x-0.58)+0.12*np.sin(5*y)
    elif kind=='random3':
        z=0.70*np.cos(2.2*np.pi*x)+0.60*np.sin(1.5*np.pi*(x+y))
    else:
        z=rng.normal(0,1,len(x))
    z += rng.normal(0,0.08,len(x))
    return (z-z.mean())/z.std()


def add_tissue_map(fig,parent_ax,bbox,kind='pathway',border=None,shuffle=False,seed=1,dotsize=26,title=None,title_color=None):
    pos=parent_ax.get_position()
    fx=pos.x0+bbox[0]*pos.width; fy=pos.y0+bbox[1]*pos.height
    fw=bbox[2]*pos.width; fh=bbox[3]*pos.height
    ax=fig.add_axes([fx,fy,fw,fh])
    vals=tissue_values(kind,seed).copy()
    if shuffle:
        np.random.default_rng(seed+99).shuffle(vals)
    ax.scatter(TISSUE_XY[:,0],TISSUE_XY[:,1],c=vals,cmap=score_cmap,vmin=-2,vmax=2,
               s=dotsize,edgecolors='white',linewidths=0.35)
    ax.set_xlim(-0.05,1.05); ax.set_ylim(-0.05,1.05); ax.set_aspect('equal')
    ax.set_xticks([]); ax.set_yticks([]); ax.set_facecolor('white')
    for sp in ax.spines.values():
        sp.set_linewidth(1.0); sp.set_color(border or COL['line'])
    if title:
        ax.set_title(title,fontsize=11.5,fontweight='bold',color=title_color or COL['ink'],pad=2)
    return ax


def globe_icon(ax,center=(0.075,0.735),radius=0.048):
    cx,cy=center
    ax.add_patch(Circle((cx,cy),radius,transform=ax.transAxes,facecolor='white',edgecolor=COL['reference'],lw=1.3))
    for wf in [0.75,1.35]:
        ax.add_patch(Ellipse((cx,cy),radius*wf,radius*2,transform=ax.transAxes,facecolor='none',edgecolor=COL['reference_light'],lw=0.9))
    for hf in [0.75,1.45]:
        ax.add_patch(Ellipse((cx,cy),radius*2,radius*hf,transform=ax.transAxes,facecolor='none',edgecolor=COL['reference_light'],lw=0.9))
    # simple land masses, deliberately schematic
    ax.add_patch(Ellipse((cx-0.014,cy+0.010),radius*0.55,radius*0.33,angle=18,transform=ax.transAxes,facecolor=COL['competitive'],edgecolor='none',alpha=0.78))
    ax.add_patch(Ellipse((cx+0.020,cy-0.010),radius*0.36,radius*0.48,angle=-15,transform=ax.transAxes,facecolor=COL['competitive'],edgecolor='none',alpha=0.78))

def matrix_icon(ax,xy=(0.17,0.66),w=0.135,h=0.15):
    x,y=xy
    rounded(ax,(x,y),w,h,fc='white',ec=COL['line'],lw=1.0,radius=0.010)
    vals=np.array([[.2,.7,.4,.9,.3],[.8,.3,.5,.2,.7],[.4,.6,.9,.5,.2],[.3,.8,.2,.6,.9],[.7,.4,.8,.3,.5]])
    n=5
    for i in range(n):
        for j in range(n):
            ax.add_patch(Rectangle((x+0.012+j*(w-0.024)/n,y+0.012+(n-1-i)*(h-0.024)/n),
                                   (w-0.029)/n,(h-0.029)/n,transform=ax.transAxes,
                                   facecolor=score_cmap(vals[i,j]),edgecolor='white',lw=0.2))


fig=plt.figure(figsize=(16,9.6),facecolor='white')
gs=fig.add_gridspec(2,2,left=0.025,right=0.975,top=0.97,bottom=0.035,wspace=0.025,hspace=0.04)
axA,axB,axC,axD=[fig.add_subplot(gs[i,j]) for i in range(2) for j in range(2)]
for ax in [axA,axB,axC,axD]: panel_background(ax)

# ---------------- A ----------------
panel_header(axA,'A','Structured tissue makes many score maps spatial')

globe_icon(axA)
axA.text(0.075,0.665,'tissue\ngeography',transform=axA.transAxes,ha='center',va='top',fontsize=12.2,fontweight='bold',color=COL['reference'])

matrix_icon(axA,(0.165,0.675),0.14,0.145)
axA.text(0.235,0.642,'spatial counts',transform=axA.transAxes,ha='center',va='top',fontsize=12.2,color=COL['ink'])

rounded(axA,(0.36,0.675),0.14,0.145,fc='white',ec=COL['line'],lw=1.0,radius=0.010)
coords=np.array([[.39,.72],[.42,.78],[.45,.71],[.47,.80],[.49,.74],[.42,.82],[.46,.76],[.40,.75]])
axA.scatter(coords[:,0],coords[:,1],s=55,c=COL['reference_light'],edgecolors='white',lw=0.5,transform=axA.transAxes,zorder=3)
axA.text(0.43,0.642,'spot coordinates',transform=axA.transAxes,ha='center',va='top',fontsize=12.2,color=COL['ink'])

rounded(axA,(0.555,0.675),0.14,0.145,fc=COL['neutral'],ec=COL['line'],lw=1.0,radius=0.010)
for cx,cy,lab,cc in [(0.590,0.755,'g1',COL['pathway']),(0.625,0.785,'g2',COL['pathway_dark']),(0.660,0.750,'g3',COL['pathway'])]:
    axA.add_patch(Circle((cx,cy),0.029,transform=axA.transAxes,facecolor=cc,edgecolor='white',lw=0.7))
    axA.text(cx,cy,lab,transform=axA.transAxes,ha='center',va='center',fontsize=9.8,color='white',fontweight='bold')
axA.text(0.625,0.713,'+',transform=axA.transAxes,ha='center',va='center',fontsize=16,color=COL['muted'],fontweight='bold')
axA.text(0.625,0.642,'pathway scoring',transform=axA.transAxes,ha='center',va='top',fontsize=12.2,color=COL['ink'])

rounded(axA,(0.765,0.69),0.16,0.11,fc='white',ec=COL['line'],lw=1.0,radius=0.012)
axA.text(0.845,0.745,'activity maps',transform=axA.transAxes,ha='center',va='center',fontsize=13.2,fontweight='bold',color=COL['ink'])

for p1,p2 in [((0.125,0.748),(0.16,0.748)),((0.306,0.748),(0.355,0.748)),((0.502,0.748),(0.55,0.748)),((0.698,0.748),(0.76,0.748))]:
    arrow(axA,p1,p2,lw=1.4,mutation=12)

add_tissue_map(fig,axA,[0.12,0.235,0.29,0.31],kind='pathway',border=COL['pathway'],seed=3,dotsize=25,title='Named pathway map',title_color=COL['pathway_dark'])
add_tissue_map(fig,axA,[0.59,0.235,0.29,0.31],kind='random',border=COL['reference'],seed=7,dotsize=25,title='Matched random-set map',title_color=COL['reference'])
arrow(axA,(0.425,0.39),(0.575,0.39),lw=1.5,mutation=12,style='<->')
axA.text(0.50,0.43,'both can be spatial',transform=axA.transAxes,ha='center',va='center',fontsize=12.2,color=COL['muted'])

rounded(axA,(0.16,0.055),0.68,0.10,fc='white',ec=COL['line'],lw=1.0,radius=0.018)
axA.text(0.50,0.104,r'$\bf{Spatially\ patterned}\ \neq\ \bf{pathway\!-\!specific}$',transform=axA.transAxes,ha='center',va='center',fontsize=19,color=COL['ink'])

# ---------------- B ----------------
panel_header(axB,'B','Absolute and competitive nulls ask different questions')
rounded(axB,(0.07,0.845),0.86,0.065,fc=COL['neutral'],ec='none',radius=0.018)
axB.text(0.50,0.878,'Same observed score map and same spatial statistic',transform=axB.transAxes,ha='center',va='center',fontsize=12.8,color=COL['ink'])
axB.text(0.50,0.838,'Moran / covariance / Gaussian process',transform=axB.transAxes,ha='center',va='center',fontsize=11.2,color=COL['muted'])

# left card
rounded(axB,(0.045,0.105),0.43,0.69,fc='white',ec=COL['absolute'],lw=1.4,radius=0.018)
axB.text(0.260,0.754,'ABSOLUTE TEST',transform=axB.transAxes,ha='center',va='center',fontsize=15.7,fontweight='bold',color=COL['absolute'])
axB.text(0.260,0.706,'permute spot labels within section',transform=axB.transAxes,ha='center',va='center',fontsize=11.8,color=COL['muted'])
axB.text(0.185,0.647,'observed',transform=axB.transAxes,ha='center',va='center',fontsize=11.0,fontweight='bold',color=COL['pathway_dark'])
axB.text(0.405,0.647,'permuted',transform=axB.transAxes,ha='center',va='center',fontsize=11.0,fontweight='bold',color=COL['muted'])
add_tissue_map(fig,axB,[0.10,0.43,0.21,0.19],kind='pathway',border=COL['pathway'],seed=3,dotsize=14)
arrow(axB,(0.315,0.525),(0.365,0.525),color=COL['absolute'],lw=1.7,mutation=13)
add_tissue_map(fig,axB,[0.37,0.52,0.09,0.09],kind='pathway',shuffle=True,seed=2,dotsize=4.8)
add_tissue_map(fig,axB,[0.37,0.40,0.09,0.09],kind='pathway',shuffle=True,seed=8,dotsize=4.8)
axB.text(0.415,0.385,'...',transform=axB.transAxes,ha='center',va='center',fontsize=16,color=COL['muted'])
rounded(axB,(0.09,0.205),0.34,0.095,fc=COL['absolute_light'],ec='none',radius=0.018)
axB.text(0.260,0.252,'Spatially associated?',transform=axB.transAxes,ha='center',va='center',fontsize=15.0,fontweight='bold',color=COL['absolute'])
axB.text(0.260,0.150,'Answer: any spatial association',transform=axB.transAxes,ha='center',va='center',fontsize=12.3,color=COL['ink'])

# right card
rounded(axB,(0.525,0.105),0.43,0.69,fc='white',ec=COL['competitive'],lw=1.4,radius=0.018)
axB.text(0.740,0.754,'COMPETITIVE MATCHED TEST',transform=axB.transAxes,ha='center',va='center',fontsize=15.4,fontweight='bold',color=COL['competitive'])
axB.text(0.740,0.706,'keep coordinates; change the gene set',transform=axB.transAxes,ha='center',va='center',fontsize=11.8,color=COL['muted'])
axB.text(0.665,0.647,'observed',transform=axB.transAxes,ha='center',va='center',fontsize=11.0,fontweight='bold',color=COL['pathway_dark'])
axB.text(0.885,0.647,'matched maps',transform=axB.transAxes,ha='center',va='center',fontsize=11.0,fontweight='bold',color=COL['reference'])
add_tissue_map(fig,axB,[0.58,0.43,0.21,0.19],kind='pathway',border=COL['pathway'],seed=3,dotsize=14)
arrow(axB,(0.795,0.525),(0.835,0.525),color=COL['competitive'],lw=1.7,mutation=13)
add_tissue_map(fig,axB,[0.84,0.545,0.08,0.08],kind='random',seed=5,dotsize=4.2)
add_tissue_map(fig,axB,[0.84,0.445,0.08,0.08],kind='random2',seed=8,dotsize=4.2)
add_tissue_map(fig,axB,[0.84,0.345,0.08,0.08],kind='random3',seed=11,dotsize=4.2)
axB.text(0.88,0.320,'...  x 1,999',transform=axB.transAxes,ha='center',va='center',fontsize=11.5,color=COL['reference'],fontweight='bold')
rounded(axB,(0.57,0.205),0.34,0.095,fc=COL['competitive_light'],ec='none',radius=0.018)
axB.text(0.740,0.252,'Exceptional among matched sets?',transform=axB.transAxes,ha='center',va='center',fontsize=14.0,fontweight='bold',color=COL['competitive'])
axB.text(0.740,0.150,'Answer: pathway-specific exceptionality',transform=axB.transAxes,ha='center',va='center',fontsize=12.0,color=COL['ink'])
axB.text(0.50,0.050,'Different scientific questions - not an accuracy ranking',transform=axB.transAxes,ha='center',va='center',fontsize=12.0,color=COL['muted'],fontstyle='italic')

# ---------------- C ----------------
panel_header(axC,'C','Each pathway receives a matched reference class')
rounded(axC,(0.045,0.22),0.38,0.60,fc='white',ec=COL['pathway'],lw=1.4,radius=0.018)
axC.text(0.235,0.772,'Target pathway (k genes)',transform=axC.transAxes,ha='center',va='center',fontsize=15.6,fontweight='bold',color=COL['pathway_dark'])
axC.text(0.235,0.720,'mean-expression deciles',transform=axC.transAxes,ha='center',va='center',fontsize=11.8,color=COL['muted'])
counts=np.array([2,3,4,5,6,8,10,8,6,4])
maxc=counts.max(); x0=0.085; y0=0.33; w=0.30; h=0.30; bw=w/13
for i,c in enumerate(counts):
    bx=x0+i*(w/10); bh=h*c/maxc
    axC.add_patch(Rectangle((bx,y0),bw,bh,transform=axC.transAxes,facecolor=COL['pathway'],edgecolor='white',lw=0.3))
    axC.text(bx+bw/2,y0-0.032,str(i+1),transform=axC.transAxes,ha='center',va='top',fontsize=9.8,color=COL['muted'])
axC.text(x0+w/2,y0-0.084,'expression decile',transform=axC.transAxes,ha='center',va='top',fontsize=11.6,color=COL['ink'])
axC.text(0.058,y0+h/2,'genes',transform=axC.transAxes,ha='center',va='center',rotation=90,fontsize=11.4,color=COL['ink'])

arrow(axC,(0.43,0.52),(0.55,0.52),color=COL['competitive'],lw=2.0,mutation=16)
rounded(axC,(0.438,0.565),0.106,0.072,fc=COL['competitive_light'],ec='none',radius=0.016)
axC.text(0.491,0.601,'sample within\neach decile',transform=axC.transAxes,ha='center',va='center',fontsize=9.6,fontweight='bold',color=COL['competitive'])

rounded(axC,(0.555,0.22),0.40,0.60,fc='white',ec=COL['reference'],lw=1.4,radius=0.018)
axC.text(0.755,0.772,'Matched gene-set references',transform=axC.transAxes,ha='center',va='center',fontsize=14.0,fontweight='bold',color=COL['reference'])
axC.text(0.755,0.720,'same size + same decile counts',transform=axC.transAxes,ha='center',va='center',fontsize=11.8,color=COL['muted'])
for r in range(3):
    yy=0.60-r*0.14
    for i,c in enumerate(counts):
        bx=0.605+i*(0.30/10); bh=0.070*c/maxc
        axC.add_patch(Rectangle((bx,yy),0.30/13,bh,transform=axC.transAxes,facecolor=COL['reference'],edgecolor='white',lw=0.2,alpha=0.92-0.16*r))
    axC.text(0.92,yy+0.030,f'ref {r+1}',transform=axC.transAxes,ha='right',va='center',fontsize=9.8,color=COL['muted'])
axC.text(0.755,0.245,'...  x 1,999 per pathway',transform=axC.transAxes,ha='center',va='center',fontsize=13.2,fontweight='bold',color=COL['reference'])

items=[('same set size',0.22),('same decile profile',0.50),('original members excluded',0.78)]
for text,x in items:
    draw_check(axC,x,0.125,0.015)
    axC.text(x,0.085,text,transform=axC.transAxes,ha='center',va='center',fontsize=10.2,color=COL['ink'])
axC.text(0.755,0.285,'weighted sets retain their signed-weight multiset',transform=axC.transAxes,ha='center',va='center',fontsize=9.8,color=COL['muted'])
axC.text(0.50,0.025,'Not explicitly matched on within-set correlation or gene-network structure',transform=axC.transAxes,ha='center',va='bottom',fontsize=10.5,color=COL['muted'],fontstyle='italic')

# ---------------- D ----------------
panel_header(axD,'D','Exceptionality is judged within the matched distribution')
rounded(axD,(0.045,0.44),0.91,0.38,fc='white',ec=COL['line'],lw=1.0,radius=0.018)
pos=axD.get_position()
dax=fig.add_axes([pos.x0+0.10*pos.width,pos.y0+0.515*pos.height,0.43*pos.width,0.20*pos.height])
x=np.linspace(0,1,400); y=np.exp(-0.5*((x-0.48)/0.15)**2); y=y/y.max(); obs=0.86
dax.fill_between(x,0,y,color=COL['reference_light'],alpha=0.55)
dax.plot(x,y,color=COL['reference'],lw=2)
dax.axvline(obs,color=COL['pathway'],lw=3)
dax.fill_between(x,0,y,where=x>=obs,color=COL['pathway'],alpha=0.20)
dax.text(0.31,0.89,'matched references',ha='center',va='top',fontsize=10.8,color=COL['reference'],fontweight='bold')
dax.text(obs-0.018,0.84,'observed',ha='right',va='top',fontsize=10.2,color=COL['pathway_dark'],fontweight='bold')
dax.set_xlim(0,1); dax.set_ylim(0,1.05); dax.set_yticks([]); dax.set_xlabel('spatial statistic',fontsize=10.8,labelpad=2,color=COL['ink'])
dax.tick_params(axis='x',labelsize=8.8,colors=COL['muted'])
for s in ['top','right','left']: dax.spines[s].set_visible(False)
dax.spines['bottom'].set_color(COL['line'])

rounded(axD,(0.59,0.505),0.33,0.23,fc=COL['competitive_light'],ec='none',radius=0.022)
axD.text(0.755,0.672,r'$p=\frac{1+b}{1+n}$',transform=axD.transAxes,ha='center',va='center',fontsize=20,color=COL['competitive'])
axD.text(0.755,0.604,'b: references at least as extreme',transform=axD.transAxes,ha='center',va='center',fontsize=10.0,color=COL['ink'])
axD.text(0.755,0.558,'n = 1,999',transform=axD.transAxes,ha='center',va='center',fontsize=10.0,color=COL['ink'])
axD.text(0.755,0.475,'extreme tail = exceptional map',transform=axD.transAxes,ha='center',va='center',fontsize=11.0,fontweight='bold',color=COL['competitive'])

axD.text(0.50,0.385,'Pt-1A pancreatic tumor: the reference changes the conclusion',transform=axD.transAxes,ha='center',va='center',fontsize=13.6,fontweight='bold',color=COL['ink'])

rounded(axD,(0.055,0.10),0.36,0.235,fc=COL['absolute_light'],ec=COL['absolute'],lw=1.2,radius=0.018)
axD.text(0.235,0.302,'Absolute endpoint',transform=axD.transAxes,ha='center',va='center',fontsize=14.0,fontweight='bold',color=COL['absolute'])
axD.text(0.155,0.222,'50/50',transform=axD.transAxes,ha='center',va='center',fontsize=22,fontweight='bold',color=COL['pathway_dark'])
axD.text(0.155,0.175,'Hallmark',transform=axD.transAxes,ha='center',va='center',fontsize=10.5,color=COL['ink'])
axD.text(0.315,0.222,'100/100',transform=axD.transAxes,ha='center',va='center',fontsize=22,fontweight='bold',color=COL['reference'])
axD.text(0.315,0.175,'random sets',transform=axD.transAxes,ha='center',va='center',fontsize=10.5,color=COL['ink'])
axD.text(0.235,0.125,'spatial association is pervasive',transform=axD.transAxes,ha='center',va='center',fontsize=11.4,fontweight='bold',color=COL['absolute'])

arrow(axD,(0.43,0.22),(0.56,0.22),color=COL['muted'],lw=2.0,mutation=16)
axD.text(0.495,0.268,'change the\nscientific reference',transform=axD.transAxes,ha='center',va='center',fontsize=10.9,color=COL['muted'])

rounded(axD,(0.575,0.10),0.36,0.235,fc=COL['competitive_light'],ec=COL['competitive'],lw=1.2,radius=0.018)
axD.text(0.755,0.302,'Competitive matched endpoint',transform=axD.transAxes,ha='center',va='center',fontsize=13.5,fontweight='bold',color=COL['competitive'])
axD.text(0.755,0.220,'7/50',transform=axD.transAxes,ha='center',va='center',fontsize=27,fontweight='bold',color=COL['pathway_dark'])
axD.text(0.755,0.174,'Hallmark pathways',transform=axD.transAxes,ha='center',va='center',fontsize=10.8,color=COL['ink'])
axD.text(0.755,0.125,'only a subset is exceptional',transform=axD.transAxes,ha='center',va='center',fontsize=11.4,fontweight='bold',color=COL['competitive'])

rounded(axD,(0.08,0.025),0.84,0.045,fc=COL['grey_light'],ec='none',radius=0.012)
axD.text(0.50,0.047,'Random sets are diagnostic comparators, not known negatives; counts reflect different nulls, not accuracy.',transform=axD.transAxes,ha='center',va='center',fontsize=10.2,color=COL['muted'])

out_dir = FilePath(__file__).resolve().parent
out_dir.mkdir(parents=True, exist_ok=True)
base = out_dir / 'fig1_framework'
fig.savefig(base.with_suffix('.png'),dpi=300,facecolor='white')
fig.savefig(base.with_suffix('.pdf'),facecolor='white')
fig.savefig(base.with_suffix('.svg'),facecolor='white')
tmp = out_dir / 'fig1_framework_tmp600.png'
fig.savefig(tmp,dpi=600,facecolor='white')
plt.close(fig)
img=Image.open(tmp).convert('RGB')
img.save(out_dir / 'fig1_framework_600dpi.tiff', dpi=(600,600), compression='tiff_lzw')
tmp.unlink(missing_ok=True)
print('saved',base)
