#! /usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright © 2018 Frederike Duembgen <frederike.duembgen@gmail.com>

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font', **{'family': 'DejaVu Sans', 'sans-serif': ['Helvetica']})
rc('text', usetex=True)

"""
plots.py: Plots for ICASSP paper on localization
"""

cmap = plt.get_cmap('Greys')


def create_plot():
    size = (4.5, 4.5)
    pos = [0.1, 0.15, 0.8, 0.8]  # left, bottom, width, height
    fig = plt.figure(figsize=size)
    ax = fig.add_subplot(111)
    plt.grid('on')
    plt.ylabel('RMSE')
    ax.set_position(pos)
    return fig, ax


def plot_against_distance(dict_methods, chosen_eps, epsilons, sigmas, saveas, title, legend=False):
    chosen_sig = np.arange(len(sigmas))
    colors = [cmap((j+1)/len(chosen_eps)) for j in range(len(chosen_eps))]
    fig, ax = create_plot()
    fig.set_size_inches(5, 4.8)
    for i, eps in enumerate(chosen_eps):
        for m in dict_methods.keys():
            label = m if i == 0 else None
            rmses = dict_methods[m]['rmses']
            ls = dict_methods[m]['linestyle']
            ms = dict_methods[m]['marker']
            plt.plot(sigmas[chosen_sig], rmses[chosen_sig, eps], color=colors[i],
                     label=label, linestyle=ls, marker=ms, fillstyle='none')
        #plt.plot(sigmas[chosen_sig],rmses[chosen_sig,eps], color=colors[i],
            #label='${}={:1.2f}$'.format(noise_label, epsilons[eps]), linestyle=ls, marker=ms,
            #fillstyle='none')
    angle = epsilons[eps]
    plt.title(title.format(angle, 180*angle/np.pi))
    plt.xlabel('$\sigma_d$[-]')
    #ax.xaxis.set_label_coords(0.94, -0.025)
    plt.tight_layout()
    if (legend):
        plt.legend(loc='upper left')
    plt.ylim([0, 0.8])
    fig.savefig(saveas)  # ,bbox_extra_artists=(lgd,),bbx_inches='tight')


def plot_against_angles(dict_methods, chosen_sig, sigmas, epsilons, saveas, title, legend=False, gaussian=False):
    chosen_eps = range(len(epsilons))
    colors = [cmap((j+1)/len(chosen_sig)) for j in range(len(chosen_sig))]
    fig, ax = create_plot()
    fig.set_size_inches(5, 5)

    def tick_function(X):
        V = X * 180 / np.pi
        return ["%.1f" % z for z in V]

    if gaussian:
        plot = ax.plot
    else:
        plot = ax.semilogx

    for i, sig in enumerate(chosen_sig):
        for m in dict_methods.keys():
            label = m if i == 0 else None
            rmses = dict_methods[m]['rmses']
            ls = dict_methods[m]['linestyle']
            ms = dict_methods[m]['marker']
            plot(epsilons[chosen_eps], rmses[sig, chosen_eps], linestyle=ls, label=label,
                 marker=ms,  color=colors[i], fillstyle='none')

    #plt.xlim([3,102])
    ax.set_ylim([0, 0.4])
    ax.set_yticks([0, 0.1, 0.2, 0.3, 0.4])
    if legend:
        ax.legend(loc='upper left')
    ax.set_xlabel('$\sigma_\\alpha$[rad]')

    ax_deg = ax.twiny()
    new_tick_locations = np.array([0, 0.2, 0.4])
    ax_deg.set_xlim(ax.get_xlim())
    ax_deg.set_xticks(new_tick_locations)
    ax_deg.set_xticklabels(tick_function(new_tick_locations))
    ax_deg.set_xlabel('$\sigma_\\alpha [^\circ]$')

    # adjust label and title positions
    #ax_deg.xaxis.set_label_coords(0.55, 1.08) #deg1
    ax_deg.xaxis.set_label_coords(0.5, 1.1)  # deg2

    #plt.title(title.format(sigmas[sig]), y=1.12) #deg1
    plt.title(title.format(sigmas[sig]), y=1.15)  # deg2

    #ax.xaxis.set_label_coords(0.55, -0.05) #deg1
    #ax.xaxis.set_label_coords(0.94, -0.025) #deg2

    plt.tight_layout()
    plt.savefig(saveas)


def plot_seaborn(dict_methods, options, method, folder='', matrix=None, figsize=None, ylabel=None, **kwargs):
    import pandas as pd
    import seaborn as sns
    if matrix is None:
        matrix = dict_methods[method]['rmses']

    rhos = np.round(np.linspace(
        options['min_rho'], options['max_rho'], options['n_rhos']), 2)
    rhos_ext = ['{} ({}$^\circ$)'.format(
        r, np.round(180*r/np.pi, 1)) for r in rhos]
    sigmas = np.round(np.linspace(
        options['min_sigma'], options['max_sigma'], options['n_sigma']), 2)
    data = pd.DataFrame(matrix, columns=rhos_ext, index=sigmas)

    f, ax = plt.subplots(figsize=figsize)

    n_ticklabels = 9 if ylabel else 0
    sns.heatmap(data, **kwargs, annot=True,  # fmt="2.2f",
                linewidths=.5, ax=ax,
                xticklabels=10, yticklabels=n_ticklabels)
    if ylabel:
        plt.ylabel('$\sigma_d$')
    plt.xlabel('$\sigma_\\alpha$')
    ax.invert_yaxis()
    title = method 
    plt.title(title)
    method = method.replace(' ', '_')
    plt.savefig('{}/heatmap_{}.eps'.format(folder, method), transparent=True)
