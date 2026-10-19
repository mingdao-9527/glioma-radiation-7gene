# ExtraCohorts track — pipeline code (run in conda env "glioma": pandas, lifelines 0.30, scikit-survival 0.27)
# Data access: GEO FTP series matrices + GPL annot; GSE108474 CEL files -> RMA in R (see rma_rembrandt.R below)
import requests, gzip, io, re, os, pandas as pd, numpy as np
from urllib3.util.retry import Retry
from requests.adapters import HTTPAdapter
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
from sksurv.metrics import concordance_index_censored as cic, cumulative_dynamic_auc
from sksurv.util import Surv
S=requests.Session(); S.mount('https://',HTTPAdapter(max_retries=Retry(total=6,backoff_factor=1.5,status_forcelist=[500,502,503,504])))
COEF={'SMC4':0.1904986020928168,'BMP2':-0.2444336196813158,'PSMC2':0.1771619980076448,'HOXC10':0.122362764863625,'HIST1H2BH':0.1279315511231319,'E2F7':0.2051678013916627,'MMS19':-0.0840639121833086}
GENES=list(COEF); coef=pd.Series(COEF)
def fetch(url,dest):
    if os.path.exists(dest) and os.path.getsize(dest)>0: return dest
    r=S.get(url,timeout=1800,stream=True); r.raise_for_status()
    with open(dest+'.tmp','wb') as f:
        for ch in r.iter_content(1<<20): f.write(ch)
    os.replace(dest+'.tmp',dest); return dest
def parse_matrix(path):
    meta={}; rows=[]; intable=False
    with gzip.open(path,'rt',errors='replace') as f:
        for line in f:
            if line.startswith('!series_matrix_table_begin'): intable=True; continue
            if line.startswith('!series_matrix_table_end'): break
            if intable: rows.append(line.rstrip('\n').split('\t')); continue
            if line.startswith('!Sample_'):
                k,_,v=line.partition('\t'); vals=[x.strip().strip('"') for x in v.rstrip('\n').split('\t')]
                meta.setdefault(k[1:],[]).append(vals)
    expr=pd.DataFrame(rows[1:],columns=rows[0]).set_index(rows[0][0]); expr=expr.apply(pd.to_numeric,errors='coerce'); expr.columns=[c.strip('"') for c in expr.columns]; expr.index=[i.strip('"') for i in expr.index]
    gsm=meta['Sample_geo_accession'][0]; ch=pd.DataFrame(index=gsm)
    for i,vals in enumerate(meta.get('Sample_characteristics_ch1',[])): ch[f'char{i}']=vals
    ch['title']=meta['Sample_title'][0]; ch['platform']=meta['Sample_platform_id'][0]
    return meta,ch,expr

def collapse(expr, ann, genes, log2=False):
    """max-mean probe per gene; returns genes x samples and probe-choice table"""
    X=expr.copy()
    if log2: X=np.log2(X.clip(lower=1))
    out={}; choice=[]
    for g in genes:
        exact=ann[ann==g].index; multi=ann[ann.fillna('').str.split(' /// ').apply(lambda L: g in L)].index
        probes=[p for p in (exact if len(exact) else multi) if p in X.index]
        if not probes: choice.append({'gene':g,'probe':None,'n_probes':0}); continue
        means=X.loc[probes].mean(1); best=means.idxmax()
        out[g]=X.loc[best]; choice.append({'gene':g,'probe':best,'n_probes':len(probes),'probe_mean':means[best]})
    return pd.DataFrame(out).T, pd.DataFrame(choice)

def score(G):  # G genes x samples (log2), z within cohort
    z=G.sub(G.mean(1),axis=0).div(G.std(1,ddof=1),axis=0)
    return (z.mul(coef.reindex(z.index),axis=0)).sum(0), z

def analyse(name,df,T,E,label,pop,adj_vars):
    """df has risk, T, E and covariates; returns metrics dict, km rows, cox rows"""
    df=df.dropna(subset=['risk',T,E]).copy(); df[E]=df[E].astype(int); df=df[df[T]>0]
    r=df.risk.values; df['risk_z']=(r-r.mean())/r.std(ddof=1)
    med=np.median(r); df['grp']=np.where(r>med,'high','low')
    met={'cohort':name,'population':pop,'label':label,'endpoint':E,'n':len(df),'events':int(df[E].sum()),'median_cut':med,'n_high':int((df.grp=='high').sum()),
         'C_index':cidx(r,df[E],df[T])}
    lr=logrank_test(df[T][df.grp=='high'],df[T][df.grp=='low'],df[E][df.grp=='high'],df[E][df.grp=='low']); met['logrank_p']=lr.p_value
    c1=CoxPHFitter().fit(df.assign(high=(df.grp=='high').astype(int))[[T,E,'high']],T,E); s=c1.summary.loc['high']
    met.update({'HR_high_vs_low':s['exp(coef)'],'HR_lo':s['exp(coef) lower 95%'],'HR_hi':s['exp(coef) upper 95%']})
    c2=CoxPHFitter().fit(df[[T,E,'risk_z']],T,E); s=c2.summary.loc['risk_z']
    met.update({'HR_perSD_unadj':s['exp(coef)'],'HR_perSD_lo':s['exp(coef) lower 95%'],'HR_perSD_hi':s['exp(coef) upper 95%'],'p_perSD':s['p']})
    for gname in ['low','high']:
        kmf=KaplanMeierFitter().fit(df[T][df.grp==gname],df[E][df.grp==gname]); met[f'median_{E}_days_{gname}']=kmf.median_survival_time_
    # time AUC (IPCW from the same cohort)
    y=Surv.from_arrays(df[E].astype(bool).values,df[T].values); tmax=df[T].max()
    for yr in (1,3,5):
        t=365.25*yr
        if t>=tmax*0.98: met[f'AUC_{yr}y']=np.nan; continue
        try: a,_=cumulative_dynamic_auc(y,y,r,[t]); met[f'AUC_{yr}y']=float(a[0])
        except Exception as e: met[f'AUC_{yr}y']=np.nan
    # adjusted Cox
    cox_rows=[]
    adj=[v for v in adj_vars if v in df.columns and df[v].notna().sum()>=0.5*len(df) and df[v].nunique()>1]
    if adj:
        d2=df[[T,E,'risk_z']+adj].dropna()
        try:
            c3=CoxPHFitter(penalizer=0.001).fit(d2,T,E)
            for v in c3.summary.index:
                s=c3.summary.loc[v]; cox_rows.append({'cohort':name,'population':pop,'label':label,'endpoint':E,'n':len(d2),'events':int(d2[E].sum()),'variable':v,'HR':s['exp(coef)'],'HR_lo':s['exp(coef) lower 95%'],'HR_hi':s['exp(coef) upper 95%'],'p':s['p'],'adjusted_for':'+'.join(adj)})
            met['HR_perSD_adj']=c3.summary.loc['risk_z','exp(coef)']; met['p_perSD_adj']=c3.summary.loc['risk_z','p']; met['adjusted_for']='+'.join(adj); met['n_adj']=len(d2)
        except Exception as e: met['adjusted_for']='FAILED: '+str(e)[:60]
    else: met['adjusted_for']='none available'
    km=km_rows(df,label,T,E,'grp',{'cohort':name,'population':pop,'endpoint':E})
    return met,km,cox_rows,df

def prep_cov(cl):
    d=pd.DataFrame(index=cl.index)
    d['age']=pd.to_numeric(cl.get('age'),errors='coerce')
    d['grade_num']=cl.grade.map({'G1':1,'G2':2,'G3':3,'G4':4}) if 'grade' in cl else np.nan
    d['IDH_wt']=cl.IDH_status.map({'Wildtype':1,'Mutant':0}) if 'IDH_status' in cl else np.nan
    d['RT']=cl.RT.map({'yes':1,'no':0}) if 'RT' in cl else np.nan
    d['sex_male']=cl.sex.map({'Male':1,'Female':0}) if 'sex' in cl else np.nan
    return d

def km_rows(df,label,T,E,grpcol,extra=None):
    rows=[]
    for gname,gdf in df.groupby(grpcol):
        kmf=KaplanMeierFitter().fit(gdf[T],gdf[E]); sf=kmf.survival_function_; ci=kmf.confidence_interval_
        for t,s in sf.itertuples():
            rows.append({'label':label,'group':gname,'time_days':t,'surv':s,'ci_low':ci.loc[t].iloc[0],'ci_high':ci.loc[t].iloc[1],'n_risk':int((gdf[T]>=t).sum()),**(extra or {})})
    return rows

def cidx(risk,E,T): return cic(E.astype(bool).values,T.astype(float).values,np.asarray(risk,float))[0]

def strat(name,df,T,E,col,val,lab):
    d=df[(df[col]==val)].dropna(subset=['risk',T,E]); d=d[d[T]>0]
    if len(d)<20 or d[E].sum()<8: return
    r=d.risk.values; z=(r-r.mean())/r.std(ddof=1); c=CoxPHFitter().fit(pd.DataFrame({'T':d[T].values,'E':d[E].astype(int).values,'z':z}),'T','E'); s=c.summary.loc['z']
    med=np.median(r); hi=r>med; lr=logrank_test(d[T][hi],d[T][~hi],d[E][hi],d[E][~hi])
    STR.append({'cohort':name,'stratum':lab,'n':len(d),'events':int(d[E].sum()),'C_index':cidx(r,d[E],d[T]),'HR_perSD':s['exp(coef)'],'HR_lo':s['exp(coef) lower 95%'],'HR_hi':s['exp(coef) upper 95%'],'p':s['p'],'logrank_p_median_split':lr.p_value})

def write_cohort(name,key):
    risk,z,cl=SCORES[name]
    cl2=cl.copy(); cl2.index.name='sample'; cl2.to_csv(f'extra/out/extra_{key}_clinical.csv')
    sc=pd.DataFrame({'risk':risk}); sc['risk_group_median']=np.where(risk>risk.median(),'high','low')
    G=COH[key.split('_')[0] if key.startswith('CGGA301') else key]['G'][cl.index]
    for g in genes: sc[f'{g}_log2']=G.loc[g]; sc[f'{g}_z']=z.loc[g]
    sc.index.name='sample'; sc.to_csv(f'extra/out/extra_{key}_scores.csv'); return len(sc)

def run_cohort(name,G,cl,E,T,eligible_mask,label_prefix,adj_vars=('age','grade_num','IDH_wt')):
    global KM,COX
    cl=cl[eligible_mask].copy(); G=G[cl.index]
    risk,z=score(G); SCORES[name]=(risk,z,cl)
    cov=prep_cov(cl); df=pd.concat([cl[[T,E]],cov],axis=1); df['risk']=risk
    for pop,mask in [('all',np.ones(len(df),bool)),('RT_treated',(df.RT==1).values),('RT_untreated',(df.RT==0).values)]:
        if mask.sum()<20 or df[mask][E].sum()<5: continue
        adj=[a for a in adj_vars if not (pop!='all' and a=='RT')]
        met,km,cox,_=analyse(name,df[mask],T,E,f'{label_prefix} ({pop})',pop,list(adj)+(['RT'] if pop=='all' and df.RT.notna().mean()>=0.8 else []))
        METS.append(met); KM+=km; COX+=cox
    # interaction
    if df.RT.notna().sum()>0 and (df.RT==0).sum()>=10 and (df.RT==1).sum()>=10:
        d=df.dropna(subset=['risk',T,E,'RT']).copy(); d=d[d[T]>0]; d['risk_z']=(d.risk-d.risk.mean())/d.risk.std(ddof=1); d['RTxrisk']=d.RT*d.risk_z
        adj=[v for v in ('age','grade_num','IDH_wt') if d[v].notna().sum()>=0.5*len(d) and d[v].nunique()>1]
        for lab,cols in [('unadjusted',['risk_z','RT','RTxrisk']),('adjusted',['risk_z','RT','RTxrisk']+adj)]:
            d2=d[[T,E]+cols].dropna()
            try:
                c=CoxPHFitter(penalizer=0.001).fit(d2,T,E); s=c.summary
                INTER.append({'cohort':name,'endpoint':E,'model':lab,'n':len(d2),'events':int(d2[E].sum()),'n_RT':int(d2.RT.sum()),'n_noRT':int((d2.RT==0).sum()),
                              'HR_risk_perSD_noRT':s.loc['risk_z','exp(coef)'],'p_risk_noRT':s.loc['risk_z','p'],'HR_RT_at_mean_risk':s.loc['RT','exp(coef)'],'p_RT':s.loc['RT','p'],
                              'HR_interaction':s.loc['RTxrisk','exp(coef)'],'HR_int_lo':s.loc['RTxrisk','exp(coef) lower 95%'],'HR_int_hi':s.loc['RTxrisk','exp(coef) upper 95%'],'p_interaction':s.loc['RTxrisk','p'],'adjusted_for':'+'.join(adj) if lab=='adjusted' else 'none'})
            except Exception as e: INTER.append({'cohort':name,'endpoint':E,'model':lab,'error':str(e)[:80]})
        # RT effect within risk groups (median split) + stratified KM (4 groups)
        med=np.median(d.risk); d['grp']=np.where(d.risk>med,'high','low')
        for gname in ['low','high']:
            dd=d[d.grp==gname]; 
            try:
                c=CoxPHFitter().fit(dd[[T,E,'RT']],T,E); s=c.summary.loc['RT']
                INTER.append({'cohort':name,'endpoint':E,'model':f'RT effect within {gname}-risk','n':len(dd),'events':int(dd[E].sum()),'n_RT':int(dd.RT.sum()),'n_noRT':int((dd.RT==0).sum()),'HR_RT_at_mean_risk':s['exp(coef)'],'HR_int_lo':s['exp(coef) lower 95%'],'HR_int_hi':s['exp(coef) upper 95%'],'p_RT':s['p']})
            except Exception as e: INTER.append({'cohort':name,'endpoint':E,'model':f'RT effect within {gname}-risk','error':str(e)[:80]})
        d['grp4']=d.grp+'_'+d.RT.map({1:'RT',0:'noRT'})
        KM.extend(km_rows(d,f'{label_prefix} RT-stratified',T,E,'grp4',{'cohort':name,'population':'RT_stratified','endpoint':E}))
    return df


# ---- cohort construction (see README §4 for clinical harmonisation rules) and calls:
# run_cohort('GSE107850_EORTC22033',G,clin,'PFS','PFS.time',mask,'EORTC 22033 PFS')  # RT = arm (RT vs TMZ)
# run_cohort('GSE108474_REMBRANDT',G_rma,clin,'OS','OS.time',mask,'REMBRANDT OS')
# run_cohort('GSE43378',...); run_cohort('GSE4412',...); run_cohort('CGGA301_primary_all',...); run_cohort('CGGA301_primary_nonoverlap',...)

# ---- rma_rembrandt.R (conda env affy_rma: r-base 4.4.3, bioconductor-affy 1.84.0, bioconductor-hgu133plus2cdf)
# library(affy); library(hgu133plus2cdf)
# files <- list.files('extra/raw/cel', pattern='CEL.gz$', full.names=TRUE)   # 550 CEL.gz from GSE108474_RAW.tar
# ab <- ReadAffy(filenames=files, compress=TRUE); eset <- rma(ab); ex <- exprs(eset)
# colnames(ex) <- sub('_.*','', basename(colnames(ex)))   # GSM ids
# data.table::fwrite(data.frame(ID_REF=rownames(ex), ex, check.names=FALSE), 'GSE108474_RMA.csv.gz')
