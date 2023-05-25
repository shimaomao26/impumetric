# Date 20230314, Copyright at Shi

import numpy as np
import pandas as pd 

import matplotlib.pyplot as plt
import seaborn as sns

import matplotlib as mpl


# plot setting
plt.rcParams['figure.dpi'] = 600
plt.rcParams['savefig.dpi'] = 600

sns.set(rc={"figure.dpi":600, 'savefig.dpi':600})
sns.set_context('notebook')
sns.set_style("ticks")

plt.rcParams['svg.fonttype'] = 'none'

plt.rcParams['pdf.use14corefonts'] = True





#  MARE
def plot_MARE(df, marecol, Rsqcol, color_seq=None, domain_=None, filename="MARE"):
    """plot MARE against Rsq. 

    Args:
        df (_type_): table of the metrics (not grouped)
        marecol (_type_): mare column name, or a list of columns to plot them at once
        Rsqcol (_type_): Rsq column name, or a list. If it is a list, should be in the same order with marecol
        color_seq: seq of the color, in the same seq with marecol, default: built-in seaborn color
        domain_: seq of the legend, in the same order with marecol, default: marecol
        filename: if None, not save to disk. Otherwise, save to the file with this name (suffix not needed)
    """
    # filter Rsq bins with low variant counts
    Rsq_bins = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 
                0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 3]
    # median of each bin, to plot the dots
    Rsq_bins_x = [0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 
                0.525, 0.575, 0.625, 0.675, 0.725, 0.775, 0.825, 0.875, 0.925, 0.975, 3]

    # Rsq_bins = [i/100 for i in range(0,101)]
    # Rsq_bins_x = [i/100 + 0.005 for i in range(0,101)]

    if color_seq:
        color_seq = color_seq
    else:
        # 8 colors same to the seaborn default. To avoid version updating, use the hard encoded
        color_seq = ['#1f77b4', '#ff7f0e', '#2ca02c', 
                '#d62728', '#9467bd', '#8c564b', 
                '#e377c2', '#7f7f7f']
    if domain_:
        domain_ = domain_
    else:
        # if the legend is not specified, use the colname of mare
        domain_ = marecol


    # function for grouping mare into Rsq bins
    def group_mare(Rsq_bins2, df2, Rsqcol2, marecol2):
        """grouping MARE into bins

        Args:
            Rsq_bins2 (list): Rsq bin intervals
            df2 (df): mare df
            Rsqcol2 (str): rsq column
            marecol2 (str): mare column

        Returns:
            _type_: list of mean MARE
        """
        tmp = []
        for i in range(len(Rsq_bins2)-1):
            # interval
            left = Rsq_bins2[i]
            right = Rsq_bins2[i+1]

            mare_ = df2.loc[df2[(df2[Rsqcol2] >= left) & (df2[Rsqcol2] < right)].index,
                        marecol2]
            if len(mare_) > 50:
                tmp.append(mare_.mean()   )
            else:
                tmp.append(np.nan   )
        return tmp

    plot_df = pd.DataFrame()
    if type(marecol) == str:
        # marecol, Rsqcol, domain_ are str
        plot_df[domain_] = group_mare(Rsq_bins2=Rsq_bins, df2=df, Rsqcol2=Rsqcol, marecol2=marecol)
    elif type(marecol) == list:
        # marecol, Rsqcol, domain_ are list
        for _ in range(len(marecol)):
            plot_df[domain_[_]   ] = group_mare(Rsq_bins2=Rsq_bins, df2=df, Rsqcol2=Rsqcol[_], marecol2=marecol[_])
    else:
        raise TypeError("marecol and Rsqcol should be a column name or list of column names")
    plot_df.index = Rsq_bins_x[:-1]

    # plot
    fig = plt.figure(figsize=(5,5))
    ax = fig.add_axes([.2,.2,.75,.75])

    # expected mare
    x=np.linspace(0.025,0.975,99)
    y=x-x**2

    for _ in range(plot_df.shape[1] ):
        #plt.plot(plot_df.index, plot_df[plot_df.columns[_]  ], alpha=0.6, c=color_seq[ _%8  ], lw=0.9, ls="--")
        plt.plot(plot_df.index, plot_df[plot_df.columns[_]  ], alpha=0.6, color=color_seq[ _%8  ], lw=0.9, ls="--")
        plt.scatter(plot_df.index, plot_df[plot_df.columns[_] ], label=plot_df.columns[_], alpha=0.9, 
        #c="", edgecolors=color_seq[  _%8  ])
        c="none", edgecolors=color_seq[  _%8  ]) # new plt version

    plt.plot(x,y, label="Rsq=EmpRsq", c="black", ls="--")
    plt.xlabel("Rsq")
    plt.ylabel("Mean MARE")
    plt.legend(fontsize="xx-small")

    #plt.xlim(0.5,1)
    #plt.ylim(0,0.25)

    if filename:
        plt.savefig(filename+".png", dpi=600)
        plt.savefig(filename+".pdf", dpi=600)
        plt.savefig(filename+".svg", dpi=600)

    return fig
    # plt.savefig("MARE.png", dpi=600)
    # plt.savefig("MARE.svg", dpi=600)
    # plt.savefig("MARE.pdf", dpi=600)




##########
# Beta_imp
##########
def plot_Beta(df, betacol, EmpRsqcol, color_seq=None, domain_=None, filename="Beta"):
    """plot Beta against EmpRsq

    Args:
        df (_type_): table of the metrics (not grouped)
        betacol (_type_): Beta column name, or a list of columns to plot them at once
        EmpRsqcol (_type_): EmpRsq column name, or a list. If it is a listm should be in the same order with Betacol
        color_seq (_type_, optional): seq of the color, in the same seq with betacol. Defaults to None (built-in seaborn color).
        domain_ (_type_, optional): seq of the legend, in the same order with betacol. Defaults to betacol.
        filename (str, optional): if None, not save to disk. Otherwise, save to the file with this name (syffix not needed). 
        Defaults to "Beta".
    """
    # filter EmpRsq bins with low variant counts

    EmpRsq_bins = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 
                0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 3]

    EmpRsq_bins_x = [0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 
                0.525, 0.575, 0.625, 0.675, 0.725, 0.775, 0.825, 0.875, 0.925, 0.975, 3]

    # Rsq_bins = [i/100 for i in range(0,101)]
    # Rsq_bins_x = [i/100 + 0.005 for i in range(0,101)]

    if color_seq:
        color_seq = color_seq
    else:
        # 8 colors same to the seaborn default. To avoid version updating, use the hard encoded
        color_seq = ['#1f77b4', '#ff7f0e', '#2ca02c', 
                '#d62728', '#9467bd', '#8c564b', 
                '#e377c2', '#7f7f7f']
    if domain_:
        domain_ = domain_
    else:
        # if the legend is not specified, use the colname of mare
        domain_ = betacol

    # function for grouping beta into EmpRsq bins
    def group_beta(EmpRsq_bins2, df2, EmpRsqcol2, betacol2):
        """grouping EmpRsq into bins

        Args:
            Rsq_bins2 (list):EmpRsq bin intervals
            df2 (df): beta df
            EmpRsqcol2 (str): emprsq column
            betacol2 (str): beta column

        Returns:
            _type_: list of mean Beta
        """
        tmp = []
        for i in range(len(EmpRsq_bins2)-1):
            # interval
            left = EmpRsq_bins2[i]
            right = EmpRsq_bins2[i+1]

            beta_ = df2.loc[df2[(df2[EmpRsqcol2] >= left) & (df2[EmpRsqcol2] < right)].index,
                        betacol2]
            if len(beta_) > 50:
                tmp.append(beta_.mean()   )
            else:
                tmp.append(np.nan   )
        return tmp


    plot_df = pd.DataFrame()
    if type(betacol) == str:
        # marecol, Rsqcol, domain_ are str
        plot_df[domain_] = group_beta(EmpRsq_bins2=EmpRsq_bins, df2=df, EmpRsqcol2=EmpRsqcol, betacol2=betacol)
    elif type(betacol) == list:
        # marecol, Rsqcol, domain_ are list
        for _ in range(len(betacol)):
            plot_df[domain_[_]   ] = group_beta(EmpRsq_bins2=EmpRsq_bins, df2=df, EmpRsqcol2=EmpRsqcol[_], betacol2=betacol[_])
    else:
        raise TypeError("beta and EmpRsqcol should be a column name or list of column names")
    plot_df.index = EmpRsq_bins_x[:-1]

    ######
    # plot
    fig = plt.figure(figsize=(5,5))
    ax = fig.add_axes([.2,.2,.75,.75])

    x=np.linspace(0.025,0.975,99)
    y=x

    for _ in range(plot_df.shape[1] ):
        #plt.plot(plot_df.index, plot_df[plot_df.columns[_]  ], alpha=0.6, c=color_seq[ _%8  ], lw=0.9, ls="--")
        plt.plot(plot_df.index, plot_df[plot_df.columns[_]  ], alpha=0.6, color=color_seq[ _%8  ], lw=0.9, ls="--")
        plt.scatter(plot_df.index, plot_df[plot_df.columns[_] ], label=plot_df.columns[_], alpha=0.9, 
        #color="", edgecolors=color_seq[  _%8  ])
        color="none", edgecolors=color_seq[  _%8  ]) # new plt version

    plt.plot(x,y, label="Rsq=EmpRsq", c="black", ls="--")
    plt.xlabel("EmpRsq")
    plt.ylabel("Mean Beta_imp")
    plt.legend(fontsize="xx-small")

    #plt.xlim(0.4,1)
    #plt.ylim(0.4,1)

    if filename:
        plt.savefig(filename+".png", dpi=600)
        plt.savefig(filename+".pdf", dpi=600)
        plt.savefig(filename+".svg", dpi=600)
    # plt.savefig("Beta.png", dpi=600)
    # plt.savefig("Beta.svg", dpi=600)
    # plt.savefig("Beta.pdf", dpi=600)

    return fig



########################################
# Comparision between Rsq and Dosage-Rsq
########################################
def plot_scatter(infofile, xcol="EmpRsq", ycol="LooRsq", xname=None, yname=None, filename="Deviation", contour=True):
    """plot Rsq against EmpRsq (deviation)

    Args:
        infofile (_type_): path to the metrics (could directly parsing the Minimac4 info file)
        xcol (_type_): EmpRsq col
        xname: label of the x-axis, default to the xcol
        ycol (_type_): Rsq col
        yname: label of the y-axis, default to the ycol
        filename (str, optional): if None, not save to disk. Otherwise, save to the file with this name (syffix not needed). 
        contour (bool, optional): if True, adding contour lines to show the intensity. NOTE: sometimes matplotlib may fail 
        the kde plot with the default parameter. This exception is not handled here. Defaults to True.
    """

    info=pd.read_table(infofile, na_values="-").dropna() # read Minimac4 info file
    x=info[xcol]
    y=info[ycol]

    g = sns.JointGrid(height=5, 
                      ratio=6, x=x,y=y, xlim=(0,1), ylim=(0,1))

    #sns.jointplot(x=x, y=y, kind='hist', joint_kws={"bins":100})
    #g.plot(sns.kdeplot, sns.histplot)

    if contour:
        # sns.kdeplot(x, y, gridsize=70, levels=10, fill=False, clip=((0,1),(0,1)),
        sns.kdeplot(x=x, y=y, gridsize=70, levels=10, fill=False, clip=((0,1),(0,1)), # new sns version
                    ax=g.ax_joint, color="#1f77b4" , cut=1, linewidths=2, bw_adjust=0.8, 
                    bw_method="scott", cumulative=False, thres=0.2)

    # sns.scatterplot(x, y,  ec="#7f7f7f", fc="#7f7f7f", s=1, linewidth=0, ax=g.ax_joint)
    sns.scatterplot(x=x, y=y,  ec="#7f7f7f", fc="#7f7f7f", s=1, linewidth=0, ax=g.ax_joint) # new sns version

    g.ax_joint.plot([0,1], [0,1], c="#e45756", lw=2, ls="-", label="y=x")

    #g.plot_marginals(sns.kdeplot)

    sns.histplot(x=x, fill=True, bins=50, linewidth=0, ax=g.ax_marg_x, color="#1f77b4")
    #sns.kdeplot(y=y, linewidth=2, ax=g.ax_marg_y)

    sns.histplot(y=y, fill=True, bins=50, linewidth=0, ax=g.ax_marg_y, color="#1f77b4")

    if filename:
        plt.savefig(filename+".png", dpi=600)
        plt.savefig(filename+".svg", dpi=600)
        plt.savefig(filename+".pdf", dpi=600)

    return g



##########################################################
# Comparision between the imputed dosage and true genotype
##########################################################
def plot_swarm_dip(ds_, xcol, ycol, filename="swarm"):
    """swarm and violin plot of imputed-dosage (diploid)

    Args:
        ds_ (_type_): table of imputed-dosage and true genotype
        xcol: colname of true genotype
        ycol: colname of imputed-dosage
        filename (str, optional): if None, not save to disk. Otherwise, save to the file with this name (syffix not needed). 
        Defaults to "swarm".
    """
    fig = plt.figure(figsize=(6,5))
    #ax = fig.add_axes([.2,.2,.75,.75])

    #sns.violinplot(data=ds_, x=0, y="t", color="#54a24b", inner="box")
    sns.swarmplot(data=ds_, x=xcol, y=ycol, color="black", edgecolors="black", size=3)
    sns.violinplot(data=ds_, x=xcol, y=ycol, scale="width", inner=None, color="white", cut=0, zorder=50)

    slope, intercept=np.polyfit(ds_[xcol], ds_[ycol], deg=1)
    reg0=intercept
    reg1=intercept+slope
    reg2=intercept+slope*2
    plt.plot([0, 1, 2], [reg0, reg1, reg2], zorder=100, c="r")

    plt.ylim(-0.1, 2.1)
    plt.xlim(-0.5, 2.5)
    plt.xticks([0, 1, 2], [0, 1, 2])

    plt.xlabel("True Genotype")
    plt.ylabel("Imputed Dosage")

    if filename:
        plt.savefig(filename+".png", dpi=600)
        plt.savefig(filename+".svg", dpi=600)
        plt.savefig(filename+".pdf", dpi=600)

    return fig

def plot_swarm_hap(ds_, xcol, ycol, filename="swarm"):
    """swarm and violin plot of imputed-dosage (hiploid)

    Args:
        ds_ (_type_): table of imputed-dosage and true allele
        xcol: colname of true allele
        ycol: colname of imputed-dosage
        filename (str, optional): if None, not save to disk. Otherwise, save to the file with this name (syffix not needed). 
        Defaults to "swarm".
    """
    fig = plt.figure(figsize=(4,5))
    #ax = fig.add_axes([.2,.2,.75,.75])

    #sns.violinplot(data=ds_, x=0, y="t", color="#54a24b", inner="box")
    sns.swarmplot(data=ds_, x=xcol, y=ycol, color="black", edgecolors="black", size=3)
    sns.violinplot(data=ds_, x=xcol, y=ycol, scale="width", inner=None, color="white", cut=0, zorder=50)

    slope, intercept=np.polyfit(ds_[xcol], ds_[ycol], deg=1)
    reg0=intercept
    reg1=intercept+slope
    plt.plot([0, 1], [reg0, reg1], zorder=100, c="r")

    plt.ylim(-0.1, 1.1)
    plt.xlim(-0.5, 1.5)
    plt.xticks([0, 1], [0, 1])

    plt.xlabel("True Allele")
    plt.ylabel("Allelic Dosage")

    if filename:
        plt.savefig(filename+".png", dpi=600)
        plt.savefig(filename+".svg", dpi=600)
        plt.savefig(filename+".pdf", dpi=600)

    return fig



##############
# Beta contour
##############
def beta_contour(rsq=None, emprsq=None, color_seq=None, domain_=None, add_label=True, filename="Beta_imp"):
    """plot the heatmap and contour lines of beta. If no rsq and emprsq added, just plot the contour lines. 
    Otherwise, plot the variants as circles (plot variants after getting the background figure is recommended, 
    since different shapes, colors, etc, could be added).

    Args:
        rsq (list): a list of scatters to add. Default to None (adding scatters to the return)
        emprsq (list): a list of scatters to add. Default to None (adding scatters to the return)
        color_seq (list, optional): color seq of the scatters. Defaults to None.
        domain_ (list, optional): legend of the scatters. Defaults to None.
        add_label (bool, optional): if true, will add labels to the edge of contour lines. Defaults to True.

    Returns:
        _type_: fig. Suggest to directly add scatters to fig
    """

    # built-in Beta values
    vhwe = 1000 # just a placeholder for the denominator
    # find the Rsq and EmpRsq coordinates for each beta
    beta = np.arange(101)/100
    df = pd.DataFrame(columns=beta)
    #df.columns = beta
    # for beta from 0 to 1, get all possible SSres
    for _ in beta:
        df[_] = np.linspace(0, vhwe*(1-_*_), num=101)
    # MARE: df/vhwe
    # Rsq (Equation 2-7)
    df_Rsq = df/vhwe + df.columns**2
    # Var(y)
    df_vdosage = df_Rsq*vhwe
    # EmpRsq (Equation 2-6)
    df_rsq = (df_Rsq - df/vhwe)*vhwe*(1/df_vdosage)
    ####
    #END#

    fig=plt.figure(figsize=(6,6), dpi=600)
    ax = fig.add_axes([0.2,0.2, 0.6,0.6])

    # colormesh
    df3 = pd.DataFrame(np.zeros([1001, 1001])  )
    df3.columns=np.arange(1001)/1000
    df3.index=1-np.arange(1001)/1000
    beta_mesh = df3.apply(lambda x : np.sqrt(x.name*x.index) )

    my_cmap = sns.color_palette("Blues", as_cmap=True)
    # cmap = mpl.cm.cool
    norm = mpl.colors.Normalize(vmin=0, vmax=1)
    ax.pcolormesh(beta_mesh.columns,beta_mesh.index,beta_mesh.values, alpha=0.9, cmap=my_cmap,vmin=0, vmax=1, shading='gouraud')

    # Contour
    #fig = plt.figure(figsize=(3,3))
    #pal=sns.color_palette("Blues" , n_colors=101)

    # dashed line
    for i in range(0,101,5):
        #plt.plot(df_rsq.iloc[:,i], df_Rsq.iloc[:,i], color=pal[i], lw=1)
        plt.plot(df_rsq.iloc[:,i], df_Rsq.iloc[:,i], color="black", lw=1, ls="--", alpha=0.8)
        
    # contour label
    if add_label:
        for i in range(10,91,10):
            plt.text(0.94, i*i/10000, str(i/100), fontsize="x-small", c="red", fontweight="bold")

    # Diagnal line
    plt.plot([0,1], [0,1], c="red", lw=1)

    if rsq:
        if color_seq:
            color_seq = color_seq
        else:
            # 8 colors same to the seaborn default. To avoid version updating, use the hard encoded
            color_seq = ['#1f77b4', '#ff7f0e', '#2ca02c', 
                    '#d62728', '#9467bd', '#8c564b', 
                    '#e377c2', '#7f7f7f']

        for _ in range(len(rsq) ):
            # make pseudo label
            plt.scatter(rsq[_], emprsq[_], label=domain_[_], c=color_seq[_])

            # plot real data point
            plt.scatter(rsq[_], emprsq[_], edgecolor="black", c=color_seq[_], zorder=100)
            # plot the 2nd variant
            # plt.scatter(0.839, 0.862, edgecolor="black",  c="#54a24b", zorder=100, marker="D")

            # make pseudo label. This is only used for plotting two variants on the same plot
            #plt.scatter(0.340683, 0.751399, label="rs142572000", edgecolor="black", c="white", zorder=-100)
            #plt.scatter(0.839, 0.862, label="rs671", edgecolor="black", c="white", zorder=-100, marker="D")

            plt.legend(loc='lower center', fontsize="x-small", title_fontsize="x-small")

    plt.xlabel("EmpRsq")
    plt.ylabel("Rsq")

    # Color key
    cb_ax = fig.add_axes([.81,.5,.01,.3])
    plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=my_cmap),cax=cb_ax, label="Beta_imp")

    if filename:
        plt.savefig(filename+".png", dpi=600)
        plt.savefig(filename+".svg", dpi=600)
        plt.savefig(filename+".pdf", dpi=600)

    return fig



##############
# MARE contour
##############
def mare_contour(rsq=None, emprsq=None, color_seq=None, domain_=None, add_label=True, filename="MARE"):
    """plot contour lines of MARE. If no rsq and emprsq added, just plot the contour lines. 
    Otherwise, plot the variants as circles (plot variants after getting the background figure is recommended, 
    since different shapes, colors, etc, could be added).

    Args:
        rsq (list): a list of scatters to add. Default to None (adding scatters to the return)
        emprsq (list): a list of scatters to add. Default to None (adding scatters to the return)
        color_seq (list, optional): color seq of the scatters. Defaults to None.
        domain_ (list, optional): legend of the scatters. Defaults to None.
        add_label (bool, optional): if true, will add labels to the edge of contour lines. Defaults to True.

    Returns:
        _type_: fig. Suggest to directly add scatters to fig
    """

    # built-in MARE values
    vhwe = 1000 # placeholder
    mare = np.arange(101)/100
    df2 = pd.DataFrame(columns=mare)
    # range of SSres
    SSres = np.linspace(0, vhwe, num=101)
    df2.columns = SSres
    # for SSres from (0. to vhwe), get all possible beta
    for _ in SSres:
        df2[_] = np.linspace(0, np.sqrt(1-_/vhwe), num=101)
    df2_Rsq = df2.columns/vhwe + df2**2
    df2_vdosage = df2_Rsq*vhwe
    df2_rsq = (df2_Rsq - SSres/vhwe)*vhwe*(1/df2_vdosage)
    #####
    #END#


    fig=plt.figure(figsize=(6,6), dpi=600)
    ax = fig.add_axes([0.2,0.2, 0.6,0.6])

    # colormesh
    df3 = pd.DataFrame(np.zeros([1001, 1001])  )
    df3.columns=np.arange(1001)/1000
    df3.index=1-np.arange(1001)/1000
    mare_mesh = df3.apply(lambda x : (1-x.name)*x.index)

    my_cmap = sns.color_palette("Blues", as_cmap=True)
    # cmap = mpl.cm.cool
    norm = mpl.colors.Normalize(vmin=0, vmax=1)
    ax.pcolormesh(mare_mesh.columns,mare_mesh.index,mare_mesh.values, alpha=0.9, cmap=my_cmap,vmin=0, vmax=1, shading='gouraud')

    # Contour
    #fig = plt.figure(figsize=(3,3))
    #pal=sns.color_palette("Blues" , n_colors=101)

    # dashed line
    for i in range(0,101,5):
        #plt.plot(df2_rsq.iloc[:,i], df2_Rsq.iloc[:,i], color=pal[i], lw=1)
        plt.plot(df2_rsq.iloc[:,i], df2_Rsq.iloc[:,i], color="black", lw=1, ls="--", alpha=0.8)
        
    # contour label
    if add_label:
        for i in range(10,91,10):
            plt.text(0, i/100, str(i/100), fontsize="x-small", c="red", fontweight="bold")

    # Diagnal line
    plt.plot([0,1], [0,1], c="red", lw=1)

    if rsq:
        if color_seq:
            color_seq = color_seq
        else:
            # 8 colors same to the seaborn default. To avoid version updating, use the hard encoded
            color_seq = ['#1f77b4', '#ff7f0e', '#2ca02c', 
                    '#d62728', '#9467bd', '#8c564b', 
                    '#e377c2', '#7f7f7f']

        for _ in range(len(rsq) ):
            # make pseudo label
            plt.scatter(rsq[_], emprsq[_], label=domain_[_], c=color_seq[_])

            # plot real data point
            plt.scatter(rsq[_], emprsq[_], edgecolor="black", c=color_seq[_], zorder=100)
            # plot the 2nd variant
            # plt.scatter(0.839, 0.862, edgecolor="black",  c="#54a24b", zorder=100, marker="D")

            # make pseudo label. This is only used for plotting two variants on the same plot
            #plt.scatter(0.340683, 0.751399, label="rs142572000", edgecolor="black", c="white", zorder=-100)
            #plt.scatter(0.839, 0.862, label="rs671", edgecolor="black", c="white", zorder=-100, marker="D")

            plt.legend(loc='lower center', fontsize="x-small", title_fontsize="x-small")

    plt.xlabel("EmpRsq")
    plt.ylabel("Rsq")

    # Color key
    cb_ax = fig.add_axes([.81,.5,.01,.3])
    plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=my_cmap),cax=cb_ax, label="MARE")

    if filename:
        plt.savefig(filename+".png", dpi=600)
        plt.savefig(filename+".svg", dpi=600)
        plt.savefig(filename+".pdf", dpi=600)

    return fig







