#!/usr/bin/env python3
"""
DetailedBalanceChecker – Python visualisation
Usage: python3 show_report.py <report.json>
Saves a PNG and prints its path to stdout.
"""

import sys, json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from collections import defaultdict, deque

GREEN = '#1a7a1a'
RED   = '#8b0000'

# ── main ─────────────────────────────────────────────────────────────────────

def main():
    if len(sys.argv) < 2:
        print("Usage: show_report.py <report.json>", file=sys.stderr)
        sys.exit(1)

    with open(sys.argv[1]) as f:
        data = json.load(f)

    name      = data["name"]
    passed    = data["pass"]
    kl        = data["kl"]
    states    = data["allStates"]
    sim_freq  = data["simFreq"]
    boltzmann = data["boltzmann"]
    db_pairs  = data["dbPairs"]
    mat_raw   = data["matrix"]
    tree_raw  = data["treeData"]
    alg_code  = data.get("algCode", "")
    n         = len(states)

    fig = plt.figure(figsize=(18, 13))
    fig.patch.set_facecolor('#f4f6f9')

    # Banner
    banner_col = GREEN if passed else RED
    verdict    = '✓  DETAILED BALANCE:  PASS' if passed else '✗  DETAILED BALANCE:  FAIL'
    fig.text(0.5, 0.975,
             f'{name}     {verdict}',
             ha='center', va='top', fontsize=13, fontweight='bold', color='white',
             bbox=dict(boxstyle='round,pad=0.45', facecolor=banner_col, alpha=1.0))

    # Grid: 3 rows × 3 cols
    gs = gridspec.GridSpec(3, 3, figure=fig,
                           height_ratios=[3.5, 4.5, 3],
                           hspace=0.45, wspace=0.30,
                           top=0.91, bottom=0.04, left=0.04, right=0.97)

    ax_freq  = fig.add_subplot(gs[0, :2])
    ax_code  = fig.add_subplot(gs[0, 2])
    n_trees  = min(n, 3)
    ax_trees = [fig.add_subplot(gs[1, i]) for i in range(n_trees)]
    for i in range(n_trees, 3):
        fig.add_subplot(gs[1, i]).axis('off')
    ax_db    = fig.add_subplot(gs[2, :2])
    ax_mat   = fig.add_subplot(gs[2, 2])

    draw_freq_chart(ax_freq, states, sim_freq, boltzmann, kl)
    draw_code_panel(ax_code, alg_code)
    for i in range(n_trees):
        s     = states[i]
        paths = tree_raw.get(str(s), [])
        draw_tree(ax_trees[i], s, paths, states)
    draw_db_table(ax_db, db_pairs)
    draw_matrix(ax_mat, states, mat_raw)

    out_path = sys.argv[1].replace('.json', '.png')
    fig.savefig(out_path, dpi=130, bbox_inches='tight',
                facecolor=fig.get_facecolor())
    print(out_path)


# ── Frequency bar chart ───────────────────────────────────────────────────────

def draw_freq_chart(ax, states, sim_freq, boltzmann, kl):
    x     = np.arange(len(states))
    w     = 0.35
    ax.bar(x - w/2, sim_freq,  w, label='Simulated',  color='#3474cc', alpha=0.85, zorder=3)
    ax.bar(x + w/2, boltzmann, w, label='Boltzmann',   color='#cc3333', alpha=0.85, zorder=3)
    ax.set_xticks(x)
    ax.set_xticklabels([f'State {s}' for s in states], fontsize=10)
    ax.set_ylabel('Probability', fontsize=10)
    ax.set_title('Simulated vs Boltzmann State Frequencies', fontsize=11, fontweight='bold')
    ax.legend(fontsize=9)
    ax.set_ylim(0, max(max(sim_freq), max(boltzmann)) * 1.28)
    ax.grid(axis='y', alpha=0.35, zorder=0)
    ax.set_facecolor('#fafbff')

    kl_col = GREEN if kl < 0.02 else RED
    kl_str = '✓ Consistent with Boltzmann' if kl < 0.02 else '✗ Significant deviation'
    ax.text(0.97, 0.95, f'KL = {kl:.5f}   {kl_str}',
            transform=ax.transAxes, ha='right', va='top',
            fontsize=9, color=kl_col, fontweight='bold')


# ── Algorithm source code panel ───────────────────────────────────────────────

def draw_code_panel(ax, code):
    ax.axis('off')
    ax.set_title('Algorithm Under Test', fontsize=10, fontweight='bold', pad=4)
    lines = code.split('\n')
    MAX   = 26
    if len(lines) > MAX:
        lines = lines[:MAX] + [f'  … ({len(lines) - MAX} more lines)']
    ax.text(0.03, 0.97, '\n'.join(lines),
            transform=ax.transAxes, fontsize=6.0, fontfamily='monospace',
            verticalalignment='top',
            bbox=dict(boxstyle='round,pad=0.5',
                      facecolor='#fffff0', edgecolor='#ccccaa', alpha=0.95))


# ── Decision tree ─────────────────────────────────────────────────────────────

def _bkey(bits):
    return 'root' if not bits else 'b' + ''.join(map(str, bits))

def _oid(bits, i):
    return _bkey(bits) + f'o{i}'

def build_tree_graph(state, paths):
    nodes    = {}
    edges    = []
    edge_lbl = {}

    nodes['root'] = dict(label=f'S={state}', color='#3a6bb5', tc='white', sq=False)

    all_pfx  = set()
    leaf_pfx = set()
    for path in paths:
        bits = tuple(path['bits'])
        leaf_pfx.add(bits)
        for k in range(len(bits)):
            all_pfx.add(bits[:k+1])

    internal_pfx = all_pfx - leaf_pfx

    for pfx in internal_pfx:
        nodes[_bkey(pfx)] = dict(label='?', color='#666', tc='white', sq=False)
    for pfx in leaf_pfx:
        nodes[_bkey(pfx)] = dict(label='↓', color='#aaa', tc='white', sq=False)

    for pfx in all_pfx:
        parent = _bkey(pfx[:-1])
        child  = _bkey(pfx)
        edges.append((parent, child))
        edge_lbl[(parent, child)] = str(pfx[-1])

    for path in paths:
        bits    = tuple(path['bits'])
        leaf_id = _bkey(bits)
        for i, out in enumerate(path['outcomes']):
            oid   = _oid(bits, i)
            dest  = out['state']
            prob  = out['probStr']
            if len(prob) > 14:
                prob = prob[:11] + '…'
            color = '#2e9e50' if dest != state else '#d47b10'
            nodes[oid] = dict(label=f'→{dest}\n{prob}', color=color, tc='white', sq=True)
            edges.append((leaf_id, oid))

    return nodes, edges, edge_lbl


def bfs_layout(edges, root):
    children = defaultdict(list)
    for u, v in edges:
        children[u].append(v)
    levels  = defaultdict(list)
    visited = {root}
    q       = deque([(root, 0)])
    while q:
        node, lv = q.popleft()
        levels[lv].append(node)
        for ch in children[node]:
            if ch not in visited:
                visited.add(ch)
                q.append((ch, lv + 1))
    max_lv = max(levels) if levels else 0
    pos    = {}
    for lv, lvnodes in levels.items():
        n = len(lvnodes)
        for i, nd in enumerate(lvnodes):
            pos[nd] = ((i + 0.5) / n,
                       1.0 - lv / (max_lv + 1) if max_lv else 0.8)
    return pos


def draw_tree(ax, state, paths, all_states):
    ax.axis('off')
    ax.set_title(f'State {state}', fontsize=9, fontweight='bold', pad=2)
    ax.set_facecolor('#f7f8fc')

    if not paths:
        ax.text(0.5, 0.5, '(no paths)', ha='center', va='center',
                transform=ax.transAxes, color='gray', fontsize=9)
        return

    nodes, edges, edge_lbl = build_tree_graph(state, paths)
    pos = bfs_layout(edges, 'root')

    ax.set_xlim(-0.08, 1.08)
    ax.set_ylim(-0.15, 1.12)

    # Edges
    for u, v in edges:
        if u not in pos or v not in pos:
            continue
        x0, y0 = pos[u]
        x1, y1 = pos[v]
        ax.annotate('', xy=(x1, y1), xytext=(x0, y0),
                    arrowprops=dict(arrowstyle='-|>', color='#888',
                                   lw=0.9, mutation_scale=8), zorder=1)
        lbl = edge_lbl.get((u, v))
        if lbl is not None:
            ax.text((x0+x1)/2 + 0.04, (y0+y1)/2, lbl,
                    fontsize=7, color='#334499', fontweight='bold',
                    ha='left', va='center', zorder=3)

    # Nodes
    for nid, info in nodes.items():
        if nid not in pos:
            continue
        x, y = pos[nid]
        if info['sq']:
            rect = plt.Rectangle((x - 0.07, y - 0.045), 0.14, 0.09,
                                  facecolor=info['color'], edgecolor='white',
                                  lw=0.8, zorder=4, transform=ax.transData)
            ax.add_patch(rect)
            ax.text(x, y, info['label'], ha='center', va='center',
                    fontsize=5.5, color=info['tc'], fontweight='bold', zorder=5)
        else:
            circ = plt.Circle((x, y), 0.045,
                               facecolor=info['color'], edgecolor='white',
                               lw=0.8, zorder=4)
            ax.add_patch(circ)
            ax.text(x, y, info['label'], ha='center', va='center',
                    fontsize=6, color=info['tc'], fontweight='bold', zorder=5)


# ── Detailed balance pair table ───────────────────────────────────────────────

def draw_db_table(ax, db_pairs):
    ax.axis('off')
    ax.set_title('Detailed Balance Check', fontsize=10, fontweight='bold', pad=4)
    if not db_pairs:
        ax.text(0.5, 0.5, '(no pairs)', ha='center', va='center',
                transform=ax.transAxes)
        return

    rows    = []
    rcolors = []
    for p in db_pairs:
        ok = p['pass']
        rows.append([str(p['i']), str(p['j']),
                     'T(i→j)·π(i) = T(j→i)·π(j) ?' ,
                     '✓ PASS' if ok else '✗ FAIL'])
        bg = '#d4f0d4' if ok else '#f0d4d4'
        rcolors.append(['#fff', '#fff', bg, bg])

    tbl = ax.table(cellText=rows,
                   colLabels=['i', 'j', 'Detailed balance', 'Result'],
                   cellColours=rcolors, cellLoc='center', loc='center')
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(10)
    tbl.scale(1.0, 1.6)
    for (r, c), cell in tbl.get_celld().items():
        if r == 0:
            cell.set_facecolor('#c8d8f0')
            cell.set_text_props(fontweight='bold')
        if c == 3 and r > 0:
            cell.set_text_props(fontweight='bold')


# ── Transition matrix table ───────────────────────────────────────────────────

def draw_matrix(ax, states, mat_raw):
    ax.axis('off')
    ax.set_title('Transition Matrix  T[i→j]', fontsize=10, fontweight='bold', pad=4)

    n     = len(states)
    cells = [['' for _ in range(n)] for _ in range(n)]
    for entry in mat_raw:
        fi = states.index(entry['from'])
        ti = states.index(entry['to'])
        s  = entry['str']
        cells[fi][ti] = (s[:10] + '…') if len(s) > 12 else s

    tbl = ax.table(cellText=cells,
                   rowLabels=[f'{s}→' for s in states],
                   colLabels=[f'→{s}' for s in states],
                   cellLoc='center', loc='center')
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(8)
    tbl.scale(1.0, 1.6)
    for (r, c), cell in tbl.get_celld().items():
        if r == 0 or c == -1:
            cell.set_facecolor('#c8d8f0')
            cell.set_text_props(fontweight='bold', fontsize=8)
        elif r > 0 and c >= 0 and r - 1 == c:
            cell.set_facecolor('#fffde0')


if __name__ == '__main__':
    main()
