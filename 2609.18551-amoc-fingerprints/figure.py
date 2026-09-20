import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import hadisst, fingerprints as fp, ersst
from ersst import trend

tr = np.load('na_trend.npy'); ax4 = np.load('na_axes.npy')
nb = np.load('null_boxes.npy')
h = hadisst.HadISST(); e, f = fp.build(h)

fig = plt.figure(figsize=(13.5, 4.3))
plt.subplots_adjust(left=.05, right=.985, top=.88, bottom=.14, wspace=.28)

# (a) trend map
ax = fig.add_subplot(1,3,1)
lo = ((ax4[2]+180)%360)-180; hi = ((ax4[3]+180)%360)-180
ext = [-70, 0, ax4[0], ax4[1]]
m = ax.imshow(tr, origin='lower', extent=ext, aspect='auto', cmap='RdBu_r',
              vmin=-0.8, vmax=0.8, interpolation='nearest')
ax.contour(np.linspace(ext[0],ext[1],tr.shape[1]), np.linspace(ext[2],ext[3],tr.shape[0]),
           np.nan_to_num(tr), levels=[0], colors='k', linewidths=.9)
ax.add_patch(plt.Rectangle((-55,46), 35, 15, fill=False, ec='k', lw=2.2))
ax.plot(-42.5, 53.5, 'k*', ms=13, mec='w', mew=.8)
ax.set_title('(a) HadISST Nov–Mar SST trend, 1871–2024\nblack box = the SST$_{SG}$ fingerprint region',
             fontsize=9.5, loc='left')
ax.set_xlabel('longitude (°E)', fontsize=8); ax.set_ylabel('latitude (°N)', fontsize=8)
ax.tick_params(labelsize=8)
cb = plt.colorbar(m, ax=ax, pad=.02); cb.set_label('°C / century', fontsize=8); cb.ax.tick_params(labelsize=7)
ax.text(-68, 32, 'star = coldest cell, −0.67\nbox mean = −0.06', fontsize=8,
        bbox=dict(fc='w', ec='.6', alpha=.9, pad=2.5))

# (b) decomposition
ax = fig.add_subplot(1,3,2)
y = f['years']; k = (y>=1871)&(y<=2024)
def sm(v, n=15):
    return np.convolve(np.asarray(v,float), np.ones(n)/n, 'valid')
yy = y[k][7:-7]
for lab, key, c in [('subpolar box (46–61°N, 55–20°W)','SG','#1b5e9c'),
                    ('global mean SST','GLOBAL','#b8860b')]:
    v = np.asarray(f[key][k], float); v = v - v[:30].mean()
    ax.plot(yy, sm(v), color=c, lw=2, label=lab)
v = np.asarray(f['SG-G'][k], float); v = v - v[:30].mean()
ax.plot(yy, sm(v), color='#992222', lw=2, label='SST$_{SG-G}$ = box − global')
ax.axhline(0, color='.6', lw=.7)
for key,c in [('SG','#1b5e9c'),('GLOBAL','#b8860b'),('SG-G','#992222')]:
    b,se = trend(y[k], np.asarray(f[key][k],float))
    ax.text(.03, {'SG':.13,'GLOBAL':.07,'SG-G':.01}[key],
            f"{b:+.3f} ± {se:.3f} °C/century", transform=ax.transAxes, color=c, fontsize=8.5)
ax.set_title('(b) the fingerprint is a difference —\nwhich term carries its trend?', fontsize=9.5, loc='left')
ax.set_xlabel('year', fontsize=8); ax.set_ylabel('°C, 15-yr mean, anomaly', fontsize=8)
ax.legend(fontsize=7.5, loc='upper left', framealpha=.9); ax.tick_params(labelsize=8)

# (c) control
ax = fig.add_subplot(1,3,3)
t = nb[:,0]
ax.hist(t, bins=45, color='.75', ec='.45', lw=.4)
ax.axvline(-0.546, color='#992222', lw=2.2)
ax.text(-0.53, ax.get_ylim()[1]*.86, ' real SST$_{SG-G}$\n −0.546', color='#992222', fontsize=8.5)
ax.axvline(np.median(t), color='k', lw=1, ls='--')
ax.text(np.median(t)+.02, ax.get_ylim()[1]*.55, 'median\n−0.048', fontsize=8)
ax.set_title('(c) the control: the same index for 868 ocean boxes\nof identical size — only 0.8% are this negative',
             fontsize=9.5, loc='left')
ax.set_xlabel('box mean − global mean, trend °C/century', fontsize=8)
ax.set_ylabel('number of boxes', fontsize=8); ax.tick_params(labelsize=8)

plt.savefig('warming_hole_decomposition.png', dpi=155)
print('wrote warming_hole_decomposition.png')
