#! /usr/bin/env python3

import matplotlib
matplotlib.use('Agg')
import pandas as pd
import numpy as np
import sys
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import gridspec
from scipy.stats.mstats import mquantiles
from scipy.stats import beta
from scipy.stats import linregress

def qq(data,ax):
    xmax = 0
    ymax = 0
    alpha = 0.9
    color = '#000000'
    n_quantiles = 100

    q_pos = np.concatenate([
        np.arange(99.) / len(data),
        np.logspace(-np.log10(len(data)) + 2, 0, n_quantiles)
    ])

    q_data = mquantiles(data, prob=q_pos, alphap=0, betap=1, limit=(0, 1))
    q_th = q_pos.copy()
    q_err = np.zeros([len(q_pos), 2])
    for i in range(0, len(q_pos)):
        q_err[i, :] = q_err[i, :] = beta.interval(
            alpha,
            len(data) * q_pos[i],
            len(data) - len(data) * q_pos[i])

    q_err[i, q_err[i, :] < 0] = 1e-15
    slope, intercept, r_value, p_value, std_err = linregress(q_th, q_data)
    xmax = np.max([xmax, -np.log10(q_th[1])])
    ymax = np.max([ymax, -np.log10(q_data[0])])

    ax.plot(
        -np.log10(q_th[n_quantiles - 1:]),
        -np.log10(q_data[n_quantiles - 1:]),
        '-',
        color=color)
    ax.plot(
        -np.log10(q_th[:n_quantiles]),
        -np.log10(q_data[:n_quantiles]),
        '.',
        color=color)
    ax.plot([0, xmax], [0, xmax], '--k',color='#f42e30')
    ax.fill_between(
        -np.log10(q_th),
        -np.log10(q_err[:, 0]),
        -np.log10(q_err[:, 1]),
        color=color,
        alpha=0.1,
    )

def manhattan(df,ax):
    df[p] = -np.log10(df[p])
    df = df.sort_values(chr_name)
    df_grouped = df.groupby((chr_name))

    colors = ['#e4508f','#556fb5',]
    x_labels = []
    x_labels_pos = []
    end = 1000
    for num, (name, group) in enumerate(df_grouped):
        group[bp] = group[bp] + end
        end = group[bp].max() + 1000
        ax.scatter(group[bp], group[p],c=colors[num % len(colors)],s=3)
        x_labels.append(name)
        x_labels_pos.append(group[bp].mean())
    ax.axhline(y=-np.log10(5e-8), color='#2222FF', linestyle='-')
    ax.set_xticks(x_labels_pos)
    ax.set_xticklabels(x_labels)

if len(sys.argv)<=1:
    sys.argv="manhattan_qq_bgenie.py $base".split()

path = '{}.txt.gz'.format(sys.argv[1])
p,bp,chr_name = '-log10p','pos','chr'
df = pd.read_csv(path,sep=' ')
df = df[df[p].notnull()]
df[chr_name] = df[chr_name].astype(int)
df[p] = [10**(-x) for x in df[p].values]

figure_tile = sys.argv[1]
fig = plt.figure(figsize=(24, 6))
gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1])
ax0 = plt.subplot(gs[0])
ax1 = plt.subplot(gs[1])
qq(df[p], ax1)
manhattan(df,ax0)
ax0.set_xlim(left=-3e7,right=2.9e9)
ax0.set_xlabel('Chromosome', fontsize=14)
ax0.set_ylabel('-log10 P', fontsize=14)
ax1.set_ylabel('Observed -log10 P', fontsize=14)
ax1.set_xlabel('Expected -log10 P', fontsize=14)
ax0.spines['right'].set_visible(False)
ax0.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
fig.suptitle(figure_tile, fontsize=20)
fig.tight_layout()
fig.savefig('{}_Manhattan_QQ.png'.format(figure_tile), dpi=300)
