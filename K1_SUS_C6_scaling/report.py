"""Summaries and static scientific plots; incomplete measurements stay incomplete."""
import csv
import json
from pathlib import Path
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def generate(root):
    root=Path(root)
    data=[json.loads(p.read_text()) for p in sorted(root.glob('Q*/result.json'),
                                                   key=lambda p:int(p.parent.name[1:]))]
    if not data:return
    rows=[]
    for r in data:
        s=r.get('block_stats',{});be=r.get('break_even',{})
        rows.append(dict(Q=r['Q'],status=r['status'],stage=r.get('limit_stage',r.get('stage')),
            build_s=r.get('hm_build_seconds'),direct_op_s=r.get('direct_operator_seconds'),
            hm_op_s=r.get('hm_operator_seconds'),field_s=r.get('field_integral_seconds'),
            speedup=r.get('operator_speedup'),rank_mean=s.get('rank_mean'),
            rank_max=s.get('rank_max'),low_rank_blocks=s.get('low_rank_blocks'),
            direct_blocks=s.get('direct_blocks'),direct_pairs=s.get('covered_direct_pairs'),
            direct_pair_fraction=s.get('direct_pair_fraction'),factor_entries=s.get('factor_complex_entries'),
            stats_complete=r.get('block_stats_complete',False),
            peak_rss_mib=max(r.get('peak_rss_mib',0),r.get('supervisor_peak_rss_mib',0)),
            break_even_real=be.get('real'),non_loss_order=be.get('non_loss_order'),
            strict_profit_order=be.get('strict_profit_order'),j_error=r.get('errors',{}).get('j'),
            j_area_error=r.get('errors',{}).get('j_area_weighted'),e_error=r.get('errors',{}).get('e')))
    with (root/'summary.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    (root/'summary.json').write_text(json.dumps(data,indent=2)+'\n')
    plt.rcParams.update({'font.size':10,'axes.grid':True,'grid.alpha':.22})
    fig,ax=plt.subplots(figsize=(8,5),layout='constrained')
    fits={}
    for key,label in [('direct_op_s','Direct C6 operator'),('hm_op_s','H operator'),
                      ('field_s','189-point field integral'),('build_s','H build + cache + near plan')]:
        pts=[(r['Q'],r[key]) for r in rows if r[key] is not None and r[key]>0]
        if pts:
            q,t=np.array(pts).T;ax.loglog(q,t,'o-',label=label)
            if len(pts)>=3 and max(q)/min(q)>=4:
                slope=float(np.polyfit(np.log(q),np.log(t),1)[0])
                fits[key]=dict(empirical_exponent=slope,points=len(q),q_min=int(min(q)),q_max=int(max(q)))
    censored=[(r['Q'],r.get('elapsed_stage_seconds')) for r in data
              if r['status']=='timeout' and r.get('limit_stage')=='hm_build']
    if censored:
        q,t=np.array(censored).T;ax.scatter(q,t,marker='^',facecolors='none',edgecolors='red',s=80,label='Build timed out (lower bound)')
    ax.set(xlabel='Panels per sector Q',ylabel='Wall time [s]',title='C6 scaling: operator and observation measured separately')
    if ax.get_legend_handles_labels()[0]:ax.legend(fontsize=8)
    fig.savefig(root/'timing_scaling.png',dpi=180);plt.close(fig)
    (root/'empirical_fits.json').write_text(json.dumps(fits,indent=2)+'\n')
    kd_edges=np.array([0,5,10,20,40,80,160,np.inf])
    rd_edges=np.array([0,.5,1,2,4,8,16,np.inf])
    groups={};points={}
    for r in data:
        path=root/f"Q{r['Q']}"/'blocks.csv'
        if not path.exists():continue
        points[r['Q']]=[]
        rng=np.random.default_rng(r['Q']);seen=0
        with path.open() as f:
            for b in csv.DictReader(f):
                try:rank=int(b['rank']);kd=float(b['kD']);rd=float(b['R_over_D'])
                except (ValueError,TypeError):continue
                ki=int(np.clip(np.searchsorted(kd_edges,kd,side='right')-1,0,len(kd_edges)-2))
                ri=int(np.clip(np.searchsorted(rd_edges,rd,side='right')-1,0,len(rd_edges)-2))
                key=(r['Q'],ki,ri)
                g=groups.setdefault(key,dict(total=0,low_rank=0,direct=0,rank_sum=0,rank_max=0))
                g['total']+=1
                if rank:
                    g['low_rank']+=1;g['rank_sum']+=rank;g['rank_max']=max(g['rank_max'],rank)
                    seen+=1
                    point=(kd,rd,rank)
                    if len(points[r['Q']])<5000:points[r['Q']].append(point)
                    else:
                        ix=int(rng.integers(seen))
                        if ix<5000:points[r['Q']][ix]=point
                else:g['direct']+=1
    with (root/'rank_by_geometry.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['Q','kD_low','kD_high','R_over_D_low','R_over_D_high',
            'terminal_blocks','low_rank_blocks','direct_blocks','mean_accepted_rank','max_rank'])
        for (q,ki,ri),g in sorted(groups.items()):
            w.writerow([q,kd_edges[ki],kd_edges[ki+1],rd_edges[ri],rd_edges[ri+1],
                        g['total'],g['low_rank'],g['direct'],g['rank_sum']/g['low_rank'] if g['low_rank'] else '',g['rank_max']])
    fig,axes=plt.subplots(1,2,figsize=(10,4),layout='constrained')
    for r in data:
        pts=points.get(r['Q'],[])
        if not pts:continue
        v=np.array(pts)
        label=f"Q={r['Q']}"+(' (partial)' if not r.get('block_stats_complete') else '')
        axes[0].scatter(v[:,0],v[:,2],s=8,alpha=.35,label=label)
        axes[1].scatter(v[:,1],v[:,2],s=8,alpha=.35,label=label)
    axes[0].set(xlabel='k max(D_target, D_source)',ylabel='Accepted ACA rank',xscale='log')
    axes[1].set(xlabel='Center distance / max cluster diameter',ylabel='Accepted ACA rank',xscale='log')
    if axes[0].get_legend_handles_labels()[0]:axes[0].legend(fontsize=7)
    fig.suptitle('Rank vs electrical size and relative separation (accepted blocks)')
    fig.savefig(root/'rank_geometry.png',dpi=180);plt.close(fig)
    lines=['# C6 H-matrix scaling results','',
           'Times are medians of completed samples. Unfinished stages are not extrapolated.',
           'Direct and H operator times include their DFT/inverse DFT but exclude geometry setup and observation fields.',
           'Build includes assembly, cache saving, near-plan setup, and block-statistics instrumentation.',
           'Rank histograms count retained terminal blocks; direct includes near blocks and failed compression.',
           'Partial block statistics are traversal prefixes and are NOT estimates of the whole operator.','',
           '| Q | Status | Direct op (s) | H op (s) | Field (s) | Build (s) | Peak RSS (MiB) | Break-even orders |',
           '|---:|---|---:|---:|---:|---:|---:|---:|']
    def fmt(v):return '—' if v is None else f'{v:.5g}'
    for r in rows:
        lines.append(f"| {r['Q']} | {r['status']} | {fmt(r['direct_op_s'])} | {fmt(r['hm_op_s'])} | {fmt(r['field_s'])} | {fmt(r['build_s'])} | {fmt(r['peak_rss_mib'])} | {fmt(r['break_even_real'])} |")
    lines+=['','Break-even is T_build / (T_direct-op - T_H-op). A nonpositive denominator has no finite break-even.',
            'The real-valued threshold is cost equality; summary.csv also reports the first non-loss and first strictly profitable integer order count.',
            'The count refers to operator applications, excluding initialization (order 0). Shared field time cancels.','',
            '## Rank histogram','']
    for r in data:
        s=r.get('block_stats',{})
        if not s:continue
        total=s['terminal_blocks']
        lines.append(f"Q={r['Q']} ({'complete' if r.get('block_stats_complete') else 'PARTIAL'}; {total} retained blocks):")
        lines.append('')
        for k,v in s['rank_histogram'].items():lines.append(f"- rank {k}: {100*v/max(total,1):.2f}% ({v})")
        lines.append(f"- direct (near + fallback): {100*s['direct_blocks']/max(total,1):.2f}% ({s['direct_blocks']})")
        lines.append('')
    lines+=['## Empirical slopes','',json.dumps(fits,indent=2),'',
            'Fits use only completed measurements with at least three Q values spanning at least 4x.',
            'Finite-range slopes do not establish asymptotic complexity or production accuracy.','',
            '![Timing](timing_scaling.png)','![Rank](rank_geometry.png)']
    (root/'REPORT.md').write_text('\n'.join(lines)+'\n')

if __name__=='__main__':generate(sys.argv[1])
