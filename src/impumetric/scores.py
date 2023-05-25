# Date 20230314, Copyright at Shi

import numpy as np
import pandas as pd

# calculate scores from tabular data

# Duplicate index (variant ID), duplicate sample ID will cause trouble

# a general encoding of missing values will not cause trouble


##################
# Helper functions
##################

# split the vcf file to GT, GP, DS, HDS
# I suggest to do it in bcftools and not in python.

def reshape(geno_file, sample_file, variant_file=None, to_memory=True, to_file=True, filename=None, verbose=False, **kwargs):
    """_summary_

    Args:
        geno_file (_type_): a genotype table, with the first 4 cols: chr, bp, ref, alt (no header) (read from disk)
        sample_file (_type_): the sample order identical to that in the geno file. Geno file and 
        variant_file (_type_, optional): Optional. If provided, will subset the variants. Defaults to None. 
            This work could also be done in bcftools (just subseting a dataframe is usually faster that extracting variants, however, 
            more disk space and memory required).
        to_memory (bool, optional): return an object. Defaults to True.
        to_file (bool, optional): write to a file. Recommended. `to_memory` and `to_file` could work together. Defaults to True.
        filename (_type_, optional): the suffix to the output file. If true, will use `filename.GT.txt.gz`. Otherwise, 
        will add `.GT.txt.gz` as a suffix to the `geno_file` and output to the current work folder. Defaults to None.
        verbose: print more information

    Raises:
        NotImplementedError: _description_

    Returns:
        _type_: return 3 tables: GT, DS, and GP table
    """
    # reshape the matrix to idv * ids
    # geno: genotype matrix
    # sample: sample list
    print("### Reshape ###")
    if verbose:
        print("Now it is reshape the geno file. It will require huge memory to process a large geno file using pandas")
        print("thus using bcftools for large file is recommended.")
    
    # NOTE: geno file is from bcftools query, so it has no header
    geno = pd.read_table(geno_file,header=None)
    if verbose:
        print("geno file should be: the first four columns should be chr, bp, ref, alt, then followed by")
        print("the individual fields. Genotypes are in GT,DS,GP, this order should not be changed")
    
    # sample should be a txt and each line contains one ID
    if verbose:
        print("sample order should be identical to that in the geno file")
    sample = []
    with open(sample_file,'r') as f:
        for line in f:
            sample.append(line.strip('\n'))
    
    # to extract the specific variants.
    # variant ID should be in the order of chr:pos:ref:alt
    variant = []
    if variant_file:
        with open(variant_file,'r') as f:
            for line in f:
                variant.append(line.strip('\n'))
        raise NotImplementedError("Currently subset of variants is not supported")
    
    # idv
    geno.index = geno[0].astype("str")+":"+geno[1].astype("str")+":"+geno[2]+":"+geno[3]
    if len(variant) > 0:
        geno = geno.loc[variant, :]
    
    # gt -> 4, 7, 10, ...
    gt = geno[range(4,geno.shape[1],3)].copy().set_axis(sample, axis=1)
    # 0|1 -> 0+1
    gt = gt.applymap(lambda x: int(x[0])+int(x[-1]))
    
    # ds -> 5, 8, 11, ...
    ds = geno[range(5,geno.shape[1],3)].copy().set_axis(sample, axis=1)
    
    # gp -> 6, 9, 12, ...
    gp = geno[range(6,geno.shape[1],3)].copy().set_axis(sample, axis=1)
    # to tuple
    gp2 = gp # default: to triplet str
    if kwargs.get("to_tuple"):
        gp2 = gp.applymap(lambda x: tuple([float(i) for i in x.split(",")]))
    if kwargs.get("to_str"):
        gp2 = gp
    # expand to df
    if kwargs.get("to_expand"):
        gp2 = pd.DataFrame(index=gp.index)
        for i in sample:
            _ = gp[i].str.split(",", expand=True).add_prefix(i+"_")
            gp2 = pd.concat([gp2, _], axis=1)    
    
    if to_file:
        if filename:
            gt.to_csv(filename+".GT.txt.gz", sep="\t", na_rep="NA")
            ds.to_csv(filename+".DS.txt.gz", sep="\t", na_rep="NA")
            gp2.to_csv(filename+".GP.txt.gz", sep="\t", na_rep="NA")
        else:
            gt.to_csv(geno_file.split("/")[-1]+".GT.txt.gz", sep="\t", na_rep="NA")
            ds.to_csv(geno_file.split("/")[-1]+".DS.txt.gz", sep="\t", na_rep="NA")
            gp2.to_csv(geno_file.split("/")[-1]+".GP.txt.gz", sep="\t", na_rep="NA")
    
    if to_memory:
        return gt,ds,gp2


# Split the triplet GP to 3 tables

def gp_to_3(gp, to_memory=True, to_file=True, filename=None, **kwargs):
    """_summary_

    Args:
        gp (_type_): GP file produced in the workflow (not a file path)
        to_memory (bool, optional): return object. Defaults to True.
        to_file (bool, optional): write to file. Recommended. Defaults to True.
        filename (_type_, optional): If provided, `filename.GP.AA.txt.gz.` will be wrote to the work folder. Defaults to None.

    Raises:
        NotImplementedError: _description_

    Returns:
        _type_: return 3 tables, AA, aa, and Aa (hetero)
    """
    print("### split GP triplets to 3 files ###")

    AA = gp.applymap(lambda x: x.split(",")[0]).astype("float16")
    aa = gp.applymap(lambda x: x.split(",")[2]).astype("float16")
    Aa = gp.applymap(lambda x: x.split(",")[1]).astype("float16")
    
    if to_file:
        if filename:
            AA.to_csv(filename+".GP.AA.txt.gz", sep="\t", na_rep="NA")
            aa.to_csv(filename+".GP.aa.txt.gz", sep="\t", na_rep="NA")
            Aa.to_csv(filename+".GP.Aa.txt.gz", sep="\t", na_rep="NA")
        else:
            raise NotImplementedError("file name is required to save the GP to AA, aa, and Aa files")
    
    if to_memory:
        return AA,aa,Aa


# Get HDS from DS and GP
# NOTE: Due to the rounding error, the converting will be strongly inflated
# Especially for the values near 0 or 1

def gp_to_hds_2(ds, aa, to_memory=True, to_file=True, filename=None, **kwargs):
    """_summary_

    Args:
        ds (_type_): dosage. range: [0-2]
        aa (_type_): probability of aa genotype. range [0-1]
        to_memory (bool, optional): retuen object. Defaults to True.
        to_file (bool, optional): write to file. Defaults to True.
        filename (_type_, optional): write `filename.HDS.hds1.repair.txt.gz` to the work folder. Must be provided. Defaults to None.

    Raises:
        NotImplementedError: _description_

    Returns:
        _type_: return 2 tables, hds1 and hds2
    """
    print("### Recover HDS from DS and GP ###")
    print("Now it is trying to convert the DS and GP to HDS")
    print("Because the float contains only 3 digits in the vcf file, this function will suffer from great rounding error!!!")
    print("This function is for testing purpose, DO NOT use it for production!!!")
    # the gp_to_hds will create many HDS ≈ 0.5, for heterozygotes, due to rounding error
    
    DS = ds.values
    aa = aa.values
    colnames = ds.columns
    rownames = ds.index
    
    if "AA" not in kwargs:
        print("only aa is supplied, because of the rounding error, the sqrt(*) may < 0, suggest to use AA + aa")
        print("however, the program will keep running with DS + aa")
    
    else:
        AA = kwargs["AA"].values
    
    #import math
    hds1 = pd.DataFrame(DS/2 + np.sqrt((DS/2)**2 - aa), index=rownames, columns=colnames)
    hds2 = pd.DataFrame(DS/2 - np.sqrt((DS/2)**2 - aa), index=rownames, columns=colnames)
    # in case nagative value in sqrt, NA is produced

    if "AA" in kwargs:
        hds3 = pd.DataFrame((2-DS)/2 + np.sqrt(((2-DS)/2)**2 - AA), index=rownames, columns=colnames)
        hds4 = pd.DataFrame((2-DS)/2 - np.sqrt(((2-DS)/2)**2 - AA), index=rownames, columns=colnames)
        # in case nagative value in sqrt, NA is produced

        # fill NA in hds1 hds2 with hds3 hds4 (with alt allele DS flip)
        hds1.update(1 - hds3)
        hds2.update(1 - hds4)
    
    # even using hds3 and hds4, they would be still some NA...
    # example: DS 1.948, GP 0.001 0.050 0.949, HDS 0.978 0.970
    # in this case couldn't get additional information, HDS = DS/2
    ds00 = ds[(ds>1.9)|(ds<0.1)]
    
    hds1.update(ds00/2)
    hds2.update(ds00/2)
    
    # if all above functions do not work, just set hds1 = ds (ds < 1) or 1 (ds > 1)
    hds1.update(ds)
    hds2.update(hds1 - 1)
    
    hds1.where(hds1 <= 1, 1, inplace=True) 
    hds2.where(hds1 >= 0, 0, inplace=True) 
        
    if to_file:
        if filename:
            hds1.to_csv(filename+".HDS.hds1.repair.txt.gz", sep="\t", na_rep="NA")
            hds2.to_csv(filename+".HDS.hds2.repair.txt.gz", sep="\t", na_rep="NA")
        else:
            raise NotImplementedError("file name is required to save the HDS1 and HDS2 files")

    if to_memory:
        return hds1, hds2


################
################
# Main functions
################
################

# HWE

def hwe_p(verbose=False, **kwargs):
    """_summary_
    HWE test according to https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2581961/
    Input:
        kwargs:
            gp: parsing the GP matrix in triplets
            AA, aa, and AA: three tables of GP

    Args:
        verbose (bool, optional): if ture, return not only the p value. Defaults to False.

    Raises:
        NotImplementedError: _description_

    Returns:
        _type_: HWE p value
    """
    print("### Hardy-Weinberge Equilibrium ###")

    if "gp" in kwargs:
        n_AA = gp.applymap(lambda x: x.split(",")[0]).astype("float16").sum(axis=1)
        n_aa = gp.applymap(lambda x: x.split(",")[2]).astype("float16").sum(axis=1)
        n_Aa = gp.applymap(lambda x: x.split(",")[1]).astype("float16").sum(axis=1)
        n = gp.shape[1]
    
    elif "AA" in kwargs and "aa" in kwargs and "Aa" in kwargs:
        n_AA = kwargs["AA"].sum(axis=1)
        n_aa = kwargs["aa"].sum(axis=1)
        n_Aa = kwargs["Aa"].sum(axis=1)
        n = kwargs["AA"].shape[1]
    
    else:
        raise NotImplementedError("HWE: only a df contains ',' separated triplates or 3 dfs of AA, aa, Aa genotypes are accepted")
        
    u = 4*n_AA*n_aa - n_Aa**2
    b = (2*n_AA+n_Aa)*(2*n_aa+n_Aa)

    chisq = n*(u/b)**2

    from scipy import stats
    p = 1 - stats.chi2.cdf(chisq , 1)  
    
    if verbose:
        return p, u, b, chisq
    else:
        return p




# MaCH rsq (based on diploid dosage)

def rsq_MaCH(ds, verbose=False):
    """calculate Rsq using DS. Some format, like pgen, stores DS only. Hence only could be used to 
    calculate MaCH Rsq, not Minimac Rsq. See https://onlinelibrary.wiley.com/doi/10.1002/gepi.20533

    Args:
        ds (_type_): DS
        verbose (bool, optional): return Rsq w/wo other metrics. Defaults to False.

    Returns:
        _type_: an array of MaCH Rsq
    """
    print("### MaCH Rsq ###")
    if verbose:
        print("Now it is calculating MaCH Rsq from DS. It is slightly different from Minimac Rsq calculated from HDS.")
    
    freq = ds.apply(np.mean, axis=1)/2
    vg = ds.apply(np.var, axis=1)
    var_hwe = 2*freq*(1-freq)

    rsq_ds = vg/var_hwe

    if verbose:
        return rsq_ds, freq, vg, var_hwe

    else:
        return rsq_ds


# MaCH rsq (based on GP)

def rsq_MaCH_GP(verbose=False, **kwargs):
    """calculate Rsq from DS, an alternative definition of MaCH Rsq, 
    https://onlinelibrary.wiley.com/doi/10.1002/gepi.20533. Generally, this is not necessary.
    Require three input:
        kwargs:
            ds: DS
            gp: the triplet GP
            aa: GP of the aa (or AA)
            Aa: GP of the Aa (hetero)

    Args:
        verbose (bool, optional): return Rsq w/wo other metrics. Defaults to False.

    Raises:
        NotImplementedError: _description_
        NotImplementedError: _description_

    Returns:
        _type_: an array of Rsq
    """
    print("### MaCH Rsq ###")
    if verbose:
        print("Now it is calculating MaCH Rsq from the alternative definition. Usually this is not necessary at all.")

    if "ds" in kwargs:
        vg = kwargs["ds"].apply(np.var, axis=1)
    else:
        raise NotImplementedError("ds in essential") 
        
    if "gp" in kwargs:
        #n_AA = gp.applymap(lambda x: x.split(",")[0]).astype("float16").sum(axis=1)
        n_aa = kwargs["gp"].applymap(lambda x: x.split(",")[2]).astype("float16").sum(axis=1)
        n_Aa = kwargs["gp"].applymap(lambda x: x.split(",")[1]).astype("float16").sum(axis=1)
        n = kwargs["gp"].shape[1]
    
    elif "aa" in kwargs and "Aa" in kwargs:
        #n_AA = AA.sum(axis=1)
        n_aa = kwargs["aa"].sum(axis=1)
        n_Aa = kwargs["Aa"].sum(axis=1)
        n = kwargs["aa"].shape[1]
        
    else:
        raise NotImplementedError("only a df contains ',' separated triplates or 2 dfs of aa, Aa genotypes are accepted")    
    
    var_gp = (4*n_aa+n_Aa)/n - ((2*n_aa+n_Aa)/n)**2
    
    rsq_gp = vg/var_gp
    
    if verbose:
        return rsq_gp, vg, var_gp

    else:
        return rsq_gp


# minimac rsq (based on HDS)

def rsq_minimac(verbose=False, **kwargs):
    """calculate Rsq (the same as Minimac)
    Accept kwargs as input,
        hds1 and hds2: the two splited table of HDS (generated using the helper function above)
        hds: (HDS table, already splited and concatenated)

    Args:
        verbose (bool, optional): return Rsq w/wo other metrics. Defaults to False.

    Raises:
        NotImplementedError: _description_

    Returns:
        _type_: an array of Rsq
    """
    print("### Minimac Rsq ###")
    if verbose:
        print("Now it is calculate Minimac Rsq. It should be equal with the value outputted by Minimac, or the INFO score")
        print("except some rounding error")

    if "hds1" in kwargs and "hds2" in kwargs:
        hds = pd.concat([kwargs["hds1"], kwargs["hds2"]], axis=1) #Prevent duplicate index, use verify_integrity option.
    elif "hds" in kwargs:
        hds = kwargs["hds"]
    else:
        raise NotImplementedError("2 dfs of hds1, hds2, or 1 df or hds are accepted")    
    
    freq = hds.apply(np.mean, axis=1)
    vg = hds.apply(np.var, axis=1)
    var_hwe = freq*(1-freq)

    rsq_hds = vg/var_hwe

    if verbose:
        return rsq_hds.tolist(), freq, vg, var_hwe

    else:
        return rsq_hds.tolist()



# INFO (based on diploid dosage)

def INFO(verbose=False, **kwargs):
    """calculate the INFO score based on GP
    See https://www.mv.helsinki.fi/home/mjxpirin/GWAS_course/material/GWAS10.html 
    Input: DS, GP of aa

    Args:
        verbose (bool, optional): return info score w/wo the individual-level, and per-individual, per-variant matrix. 
        Warning: For large dataset, it will take long time and large disk space to write this tabel in plain text.
        Defaults to False.

    Raises:
        NotImplementedError: _description_
        NotImplementedError: _description_

    Returns:
        _type_: an array of INFO
    """
    print("### IMPUTE INFO score ###")

    if "ds" in kwargs:
        freq = kwargs["ds"].apply(np.mean, axis=1)/2
        #n = kwargs["ds"].shape[1]
    else:
        raise NotImplementedError("ds in essential") 
    
    w = 2*freq*(1-freq) #HW genotype variance
    
    if "gp" in kwargs:
        raise NotImplementedError("gp is not implemented, use aa")  
        
    elif "aa" in kwargs:
        # v_i = 2*p_i(aa) + ds_i - ds_i**2
        v_i = 2*kwargs["aa"]+kwargs["ds"]-kwargs["ds"].applymap(lambda x: x**2)
        
    # v = sum(v_i)
    #v = 2*kwargs["aa"].sum(axis=1)+kwargs["ds"].sum(axis=1)-kwargs["ds"].applymap(lambda x: x**2).sum(axis=1)
    v = v_i.mean(axis=1)
    info = 1-v/w
    
    if verbose:
        print("Now it is calculating INFO for per variant and per individual")
        print("It could be summarized along variants or individuals")
        print("But it is not widely-accepted. Be sure what you want to do!!!")
        # per variant, per individual v/w, may used to find the high-info and low-info regions?
        v_i_w = v_i.divide(w, axis=0)
        # individual-level mean info (I didn't see any usage by others)
        info_2 = 1-v_i_w.mean(axis=0)
    
    if verbose:
        return info, v, w, v_i, v_i_w, info_2
    else:
        return info


# Dosage r^2

def cor_emp(input_ds1,input_ds2):
    """calculate the empirical correlation

    Args:
        input_ds1 (_type_): DS of the first dataset, e.g., the imputed
        input_ds2 (_type_): DS of the second dataset, e.g., the WGS

    Returns:
        _type_: an array of the correlation (Warning: not the squared correlation)
    """
    print("### correlation r ###")
    print("Now it is calculating the Pearson correlation")
    #print("NOTE: it will return r, not r^2, you could just square the r to get Dosage-Rsq (as many people do)")
    #print("However, the negative correlation should be interpreted as r^2 ≤ 0 (worse than no correlation)")
    #print("If you are using the impumetric parsing function, negative r will be converted to r^2 = 0")
    
    #rsq_ds = np.corrcoef(input_ds1, input_ds2)**2
    #return rsq_ds
    
    cor_list = []

    for i in np.arange(len(input_ds1)):
        # to allow different index (eg hg19 and hg38, use iloc)
        a=input_ds1.iloc[i,:]
        b=input_ds2.iloc[i,:]

        # because vcf may contain missing, handle missing values
        cor_list.append(np.corrcoef(a[ ~a.isna().values & ~b.isna().values ],   b[ ~a.isna().values & ~b.isna().values ])[0][1])
        
    return cor_list


# heterozygosity

def heterozygosity(ds, gp_Aa):
    """calculate the expected heterozygosity (based on HWE) and the observed
    This function is based on the GP of Aa, rather than the hard call genotypes.

    Args:
        ds (_type_): DS
        gp_Aa (_type_): GP of AA. This is not the hard call GT.

    Returns:
        _type_: two arrays of the observed and expected heterozygosity
    """
    print("### heterozygosity ###")

    # expected and observed heterozygosity
    h_obs = gp_Aa.mean(axis=1)
    
    freq = ds.mean(axis=1)/2
    h_exp = 2*freq*(1-freq)
    
    return h_obs, h_exp



# Concordance and IQS

def IQS(gt, verbose=False, **kwargs):
    """calculate concordance and IQS. For IQS, refer to https://pubmed.ncbi.nlm.nih.gov/20300623/
    NOTE: this function is to calculate the concordance base on the GP, not hard call GT.

    Input:
        kwargs:
            ds: DS (optional, now deprecated)
            AA: GP of AA
            Aa: GP of Aa
            aa: GP of aa
    Args:
        gt (_type_): matrix of true genotype
        verbose (bool, optional): True: return IQS + concordance/IQS matrix, there matrix will be returned in data frame, 
        thus require large memory.
        False: return IQS. Defaults to False.

    Raises:
        NotImplementedError: _description_
        NotImplementedError: _description_

    Returns:
        _type_: array of IQS, matrix of concordance
    """
    # the concordance of 3 genotypes
    # IQS: concordance adjusted by marginal allele frequency
    # marginal allele frequency 1: AF in the imputed dataset
    # marginal allele frequency 2: AF in the ground-truth set
    print("### Concordance and IQS ###")
    if verbose:
        print("Now it is calculating concordance and IQS.")
        print("NOTE: The concordance and IQS here are developed for GP.")
        print("In this function, if the GP of alt/alt is 0.2, and the true GT is alt/alt")
        print("Then the concordance will be 0.2 (20 percent agreement), not 0 (if using hard call GT)")
        print("Thus it is different from https://pubmed.ncbi.nlm.nih.gov/20300623/")
    
    # marginal DS
    if "ds" in kwargs:
        #freq = kwargs["ds"].mean(axis=1)/2
        n_sample = kwargs["ds"].shape[1]
        n_variant = kwargs["ds"].shape[0]
        gt.index = ds.index
    elif "Aa" in kwargs:
        n_sample = kwargs["Aa"].shape[1]
        n_variant = kwargs["Aa"].shape[0]    
        gt.index = Aa.index   
    else:
        raise NotImplementedError("a table of DS or 3 tables GP is essential") 
    
    # marginal 3 gt
    # 1: AA
    # 2: Aa
    # 3: aa

    # marginal count and frequency of 3 TRUE genotypes
    # NOTE: gt has missing
    # .count --> will count FLASE as 1; .sum --> will treat FALSE as 0 and True as 1
    # .count --> will treat NA as 0; .sum --> will treat True as the 1 but the numberic value as number
    n_1 = gt[gt == 0].count(axis=1)
    n_2 = gt[gt == 1].count(axis=1)
    n_3 = gt[gt == 2].count(axis=1) 
    n = n_1 + n_2 + n_3
    # marginal count of Non-NA
    n_4 = n_sample - gt.isnull().sum(axis=1)
    # if verbose:
    #     print("max n: {} (should equal to the samples size)".format(n.max()))
    #     print("minimum n: {} (should equal to the minimal non-missing sample size)".format(n.min()))
    #     print("max n_3: {} (max count of the alt allele)".format(n_3.max()))
    #     print("minimum n_3: {} (min count of the alt allele)".format(n_3.min()))  
    assert (n == n_4).all(), "n != n_4, maybe some genotypes are non of {0,1,2,NA}" # the summed 0, 1, 2 == sample - NA
    
    # marginal 3 gp
    if "gp" in kwargs:
        raise NotImplementedError("gp is not implemented, use AA, Aa, and aa")  
        
    elif "aa" in kwargs and "Aa" in kwargs and "AA" in kwargs:
        
        # remove missing data, to keep the marginal count consistent
        AA = kwargs["AA"][~gt.isnull()]
        Aa = kwargs["Aa"][~gt.isnull()]
        aa = kwargs["aa"][~gt.isnull()]
        
        # marginal count and frequency of 3 imputed genotypes (sum the probability)
        n1_ = AA.sum(axis = 1)
        n2_ = Aa.sum(axis = 1)
        n3_ = aa.sum(axis = 1)
        
        n_ = n1_ + n2_ + n3_

        # concordance
        # 9 tables of concordance betweeen gt and gp      
        # gt=2
        m33 = aa[gt == 2]
        m23 = Aa[gt == 2]
        m13 = AA[gt == 2]
        
        # gt=1
        m32 = aa[gt == 1]
        m22 = Aa[gt == 1]
        m12 = AA[gt == 1]
        
        # gt=0
        m31 = aa[gt == 0]
        m21 = Aa[gt == 0]
        m11 = AA[gt == 0]
        
        # concordance matrix (not normalized by non-missing genotypes)
        # because the .fillna, missing genotypes were filled with 0
        concor = m33.fillna(0) + m22.fillna(0) + m11.fillna(0)
        concor = concor[~gt.isnull()]
        # chance agreement
        # p_c is fixed by the marginal AF of WGS and imputation
        p_c = (n_1*n1_ + n_2*n2_ + n_3*n3_ )/(n*n_)
        # concordance subtract the chance agreement
        concor0 = concor.subtract( p_c, axis=0).divide(1-p_c, axis=0)
        # summarize along variant or individual
        iqs = concor0.mean(axis=1)
        iqs_i = concor0.mean(axis=0)
        
        # Deprecated
        # if verbose == 2:
        #     # concordance of 3 genotypes
        #     p33 = ( n_3*n3_ )/(n*n_)
        #     concor33 = m33.subtract( p33, axis=0).divide(1-p33, axis=0)
        #     iqs33 = concor33.mean(axis=1)
        #     iqs_i33 = concor33.mean(axis=0)
            
        #     p22 = ( n_2*n2_ )/(n*n_)
        #     concor22 = m22.subtract( p22, axis=0).divide(1-p22, axis=0)
        #     iqs22 = concor22.mean(axis=1)
        #     iqs_i22 = concor22.mean(axis=0)
            
        #     p11 = ( n_1*n1_ )/(n*n_)
        #     concor11 = m11.subtract( p11, axis=0).divide(1-p11, axis=0)
        #     iqs11 = concor11.mean(axis=1)
        #     iqs_i11 = concor11.mean(axis=0)
        
        
    #if verbose == 2:
    #    return iqs,iqs_i,concor, p_c,concor0,m33, m22, m11,   n_3, n_2, n, n_1, n3_, n2_, n1_, n_,   concor33,iqs33,iqs_i33,concor22,iqs22,iqs_i22,concor11,iqs11,iqs_i11
    #elif verbose == 1:
    if verbose:
        # return variant IQS, individual IQS, concordance matrix, IQS matrix
        return iqs,iqs_i,    concor, concor0  #, p_c,concor0,   m33, m22, m11
    else:
        return iqs, iqs_i


# Within- and between-group variance

def get_var(ds, gt, verbose=True):
    """getting the between-group and within-group variance.
    Most values are returned as SS (sum of squares)

    Args:
        ds (_type_): the imputed DS
        gt (_type_): the true GT
        verbose (bool, optional): _description_. Defaults to True.

    Returns:
        _type_: array with the length of n_variants
            tvar, total variance 
            bvar2, SS between aa (GT=2) and mean (AAF) 
            bvar1, SS between Aa (GT=1) and mean (AAF) 
            bvar0, SS between AA (GT=0) and mean (AAF) 
            ivar2, SS of aa (GT=2) 
            ivar1, SS of aa (GT=1) 
            ivar0, SS of aa (GT=0) 
            n2, count of samples with GT == 2 
            n1, count of samples with GT == 1 
            n0, count of samples with GT == 0 
            n, count of samples with GT == 2, 1, or 0 (missing is not counted) 
            u2, mean imputed dosage of GT==2 
            u1, mean imputed dosage of GT==1 
            u0, mean imputed dosage of GT==0 
            u, mean imputed dosage of GT==2, 1, or 0 (missing is not counted, so it will be different from AAF) 
            u_imp, AAF (sites missing in GT but are imputed will contribute to this value)
    """
    print("Now it is calculating within- and between-group variance.")
    if verbose:
        print("This is a latent function, for checking the data in depth.")
        print("For general use, MARE and Beta_imp give a better summary of within- and between-group variance.")
        # n * var = sum(n_i * var_i) + sum(n_i * (u_i - u)**2)
    
    gt.index = ds.index
    
    u = ds[~gt.isnull()].mean(axis=1) # only nonmissing data
    u_imp = ds.mean(axis=1) # mean of imputed dosage (0-2)
    
    m2 = ds[gt == 2]
    m1 = ds[gt == 1]
    m0 = ds[gt == 0]
    
    u2 = m2.mean(axis = 1) # mean of gt = 2
    u1 = m1.mean(axis = 1)
    u0 = m0.mean(axis = 1)
    
    v2 = m2.var(axis = 1) # var of gt = 2
    v1 = m1.var(axis = 1)
    v0 = m0.var(axis = 1)
    
    n2 = gt[gt == 2].count(axis=1) # n of gt = 2
    n1 = gt[gt == 1].count(axis=1)
    n0 = gt[gt == 0].count(axis=1)
    n = n2 + n1 + n0 # non-missing 
    
    # between group
    bvar2 = n2*(u2-u)**2
    bvar1 = n1*(u1-u)**2
    bvar0 = n0*(u0-u)**2
    
    # within group
    ivar2 = v2*n2 
    ivar1 = v1*n1 
    ivar0 = v0*n0
    
    # total
    tvar = (bvar2.fillna(0) + bvar1.fillna(0) + bvar0.fillna(0) + \
            ivar2.fillna(0) + ivar1.fillna(0) + ivar0.fillna(0))/n
    
    if verbose:
        # tvar, bvar, true genotype count, mean dosage of each true genotype 
        return tvar, bvar2, bvar1, bvar0, ivar2, ivar1, ivar0, n2, n1, n0, n, u2, u1, u0, u, u_imp
    else:
        return tvar, bvar2, bvar1, bvar0, ivar2, ivar1, ivar0


def skew_kurtosis(ds, gt):
    """calculate the skewness and kurtosis for each GT group. DS according to the missing sites in GT will be masked.

    Args:
        ds (_type_): the imputed DS
        gt (_type_): the true GT

    return: array
        s2, skewness of GT==2 
        s2, skewness of GT==1 
        s0, skewness of GT==0 
        k2, kurtosis of GT==2 
        k1, kurtosis of GT==1 
        k0, kurtosis of GT==0
    """
    print("Now it is calculating skewness and kurtosis.")
    if verbose:
        print("This is a latent function, for checking the data in depth.")
        print("For general use, MARE and Beta_imp give a better summary of the imputed-dosage distribution.")

    from scipy.stats import skew
    #from scipy.stats import skewtest
    from scipy.stats import kurtosis
    #from scipy.stats import kurtosistest
    
    gt.index = ds.index
    
    m2 = ds[gt == 2]
    m1 = ds[gt == 1]
    m0 = ds[gt == 0]
    
    s2=skew(m2,nan_policy='omit', axis=1)
    s1=skew(m1,nan_policy='omit', axis=1)
    s0=skew(m0,nan_policy='omit', axis=1)
    
    # skewtest is not valid with less than 8 samples
#     s2p=skewtest(m2, nan_policy='omit', axis=1)
#     s1p=skewtest(m1, nan_policy='omit', axis=1)
#     s0p=skewtest(m0, nan_policy='omit', axis=1)   
    
    k2=kurtosis(m2, nan_policy='omit', axis=1)
    k1=kurtosis(m1, nan_policy='omit', axis=1)
    k0=kurtosis(m0, nan_policy='omit', axis=1)
    
#     k2p=kurtosistest(m2, nan_policy='omit', axis=1)
#     k1p=kurtosistest(m1, nan_policy='omit', axis=1)
#     k0p=kurtosistest(m0, nan_policy='omit', axis=1)
    
    return s2,s1,s0,k2,k1,k0,     #s2p, s1p,s0p,k2p,k1p,k0p



# MARE and Beta_imp from diploid data

def mare_beta_dip(ds, true_gt, verbose=False):
    """calculate MARE and Beta_imp from diploid data. DS according to the missing sites in GT will be masked. 
    In case the np.polyfit fail, NA will be returned.

    Args:
        ds (_type_): imputed dosage
        true_gt (_type_): true WGS genotype

    Returns:
        _type_: list
            beta, Beta_imp 
            mare, MARE 
            ssres, residual sum of squares. For advanced users, ssres for each GT group could be obtained by 
            modifying the source code.
    """
    print("### MARE and Beta_imp from diploid data ###")

    beta = []
    ssres = []
    mare = []
    
    vid = ds.index

    # suppress RankWarning: Polyfit may be poorly conditioned, https://numpy.org/doc/stable/reference/generated/numpy.polyfit.html
    import warnings
    warnings.simplefilter('ignore', np.RankWarning)
    # in very few times, it may report the below error:
    # Intel MKL ERROR: Parameter 4 was incorrect on entry to DGELSD. 

    for i in range(len(ds)):
        try:
            lreg = np.polyfit( true_gt.iloc[i, :][~true_gt.iloc[i, :].isna()], 
                                ds.iloc[i, :][~true_gt.iloc[i, :].isna()], deg=1)
            
            intercept = lreg[1]
            slop = lreg[0]

            # mean of each true genotype
            a = intercept
            b = intercept + slop
            c = b + slop


            dose0 = ds.iloc[i, :][true_gt.iloc[i, :] == 0]
            dose1 = ds.iloc[i, :][true_gt.iloc[i, :] == 1]
            dose2 = ds.iloc[i, :][true_gt.iloc[i, :] == 2]

            ssres0 = np.var(dose0)*len(dose0) + len(dose0)* ( (np.mean(dose0) - a)**2  )
            ssres1 = np.var(dose1)*len(dose1) + len(dose1)* ( (np.mean(dose1) - b)**2  )
            ssres2 = np.var(dose2)*len(dose2) + len(dose2)* ( (np.mean(dose2) - c)**2  )

            # handle nan, otherwise if any of ssres0,1,2 is nan, mare_ will be nan
            ssres0 = np.nan_to_num(ssres0)
            ssres1 = np.nan_to_num(ssres1)
            ssres2 = np.nan_to_num(ssres2)

            ssres_ = ssres0 + ssres1 + ssres2

            mare_ = ssres_/(len(dose0) + len(dose1) + len(dose2))/2/(ds.iloc[i, :].mean()/2)/(1 - ds.iloc[i, :].mean()/2)

            beta.append(slop)
            ssres.append(ssres_)
            mare.append(mare_)
        
        # more ways could be used to fit the regression line, however
        # in real work, polyfit rarely fail for non-rare variants, and mare is not developed for rare variants 
        except:
            beta.append(np.nan)
            ssres.append(np.nan)
            mare.append(np.nan)

    if verbose:
        return vid, beta, mare, ssres
    else:
        return vid, beta, mare



# MARE and Beta_imp from haploid data

def mare_beta_hap(loo_ds):
    """calculate MARE and Beta_imp using the loodosage outputed by Minimac4. Since this file is generally small, 
    currently this implementation does not require the pre-processing of vcf file.

    The imputation server may perform some additional QC, so this file may have a different length to that of the info file 
    and array data. Variant ID is returned for checking.

    Args:
        loo_ds (_type_): path to the loodosage vcf file.

    Returns:
        _type_: list
            vid: variant ID list
            beta: Beta_imp
            mare: MARE
    """
    print("### MARE and Beta_imp from haploid data ###")
    # because loodosage file is usually a small vcf, there is no need of further spliting the vcf

    true_gt0=pd.read_table(loo_ds,  
                       comment='#', header=None).iloc[:,9:]

    vid = pd.read_table(loo_ds,  
                       comment='#', header=None).iloc[:,2].tolist()

    #true_gt0 = loo_ds
    gt = pd.concat([true_gt0.applymap(lambda x: x.split(":")[0].split("|")[0]), 
                     true_gt0.applymap(lambda x: x.split(":")[0].split("|")[1])], axis=1)

    #lds = pd.read_table(fname1+i+fname2, comment='#', header=None).iloc[:,9:]
    lds = pd.concat([true_gt0.applymap(lambda x: x.split(":")[1].split("|")[0]), 
                    true_gt0.applymap(lambda x: x.split(":")[1].split("|")[1])], axis=1)
    # the above should be two tables with identical dimension. And no missing data in it.

    # get beta from two tables of true genotype and Loo dosage.
    beta = []
    mare= []

    for i in range(len(gt)):
        real = gt.iloc[i, :].astype("int").values
        dosage = lds.iloc[i, :].astype("float").values
            
            # removing the WGS missing sites. Loo imputation does not contain missing.
            #dosage = dosage[~np.isnan(real)]
            #real = real[~np.isnan(real)]
            
            # old version: if the dosage or real gt was 0, set beta to 0
    #         if dosage.sum() == 0:
    #             beta.append(0)
    #         elif real.sum() == 0:
    #             beta.append(0)
    #         else:
    #             beta.append(np.polyfit(real, dosage, deg=1)[0])
                
            # new version: if the dosage or real gt was 0, set beta to na, because beta was meaningless
            # in this case, mare is still the "residual" variance. since no regression, no residual, set mare also to NA
        if dosage.sum() == 0:
            beta.append(np.nan)
            mare.append(np.nan)
        elif real.sum() == 0:
            beta.append(np.nan)
            mare.append(np.nan)
        elif dosage.sum() == gt.shape[1]:
            beta.append(np.nan)
            mare.append(np.nan)
        elif real.sum() == gt.shape[1]:
            beta.append(np.nan)
            mare.append(np.nan)
        else:
            # in a two variable case, the regression line will always pass the mean, no need to do polyfit
            real = gt.iloc[i, :].astype("int").values
            dosage = lds.iloc[i, :].astype("float").values

            ssres1 = np.var(dosage[real==1])*(real==1).sum()
            ssres0 = np.var(dosage[real==0])*(real==0).sum()

            mare.append( (ssres1 + ssres0)/len(real)/np.mean(dosage)/(1-np.mean(dosage))  )

            #beta.append(np.polyfit(real, dosage, deg=1)[0])
            # new version: for rare variants, polyfit may fail. set them to NA
            try:
                beta.append(np.polyfit(real, dosage, deg=1)[0])
            except:
                beta.append(np.nan)

    return vid, beta, mare

