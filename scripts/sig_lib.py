import pandas as pd, numpy as np, os, time, warnings
warnings.filterwarnings('ignore')
from scipy import stats
from statsmodels.stats.multitest import multipletests
from sksurv.util import Surv
from sksurv.linear_model import CoxPHSurvivalAnalysis, CoxnetSurvivalAnalysis
from sksurv.ensemble import RandomSurvivalForest, GradientBoostingSurvivalAnalysis, ComponentwiseGradientBoostingSurvivalAnalysis
from sksurv.svm import FastSurvivalSVM
from sksurv.metrics import concordance_index_censored as cic, cumulative_dynamic_auc
from sklearn.cross_decomposition import PLSRegression
from sklearn.decomposition import PCA
from sklearn.model_selection import KFold
from sklearn.inspection import permutation_importance
from lifelines import KaplanMeierFitter, CoxPHFitter, NelsonAalenFitter
from lifelines.statistics import logrank_test

P='/Users/laowang/.claude-science/orgs/4321deaf-044c-4bf9-be37-0197b924b0df/artifacts/proj_39d29643281e/'
PATHS=dict(
 tcga_expr=P+'fa4f9bc9-82e3-4fff-8a57-a93296d84899/v7f029025_tcga_expr.parquet',
 tcga_clin=P+'5a8cd171-ca02-410c-8e30-4742b6625e71/v50e1bb7a_tcga_clin.parquet',
 c693_expr=P+'fe02517d-acff-475c-b4cb-ba1c33975958/v126a9f6c_cgga693_expr.parquet',
 c693_clin=P+'37c9414a-faf7-4dca-8720-7937b9aca5e2/v61b83afa_cgga693_clin.parquet',
 c325_expr=P+'8b278b26-3f93-4a34-a088-486f9c4c8a8c/v14827af8_cgga325_expr.parquet',
 c325_clin=P+'3831ce55-ac0c-48ce-ac28-647a73bbcfdf/vef125773_cgga325_clin.parquet',
 shared=P+'6b9fd85a-de4d-4c0b-a770-f44e93f9a5be/v60cc31dd_shared_genes.txt',
 pool=P+'eff6c189-7079-4eb7-85d0-abdb292506c4/v2efdc5af_radiation_gene_pool.csv',
 cand=P+'9f56299d-5f1d-49ed-ab9e-894818e2ac53/v1a7eaa9b_pool_candidates.csv',
 kim=P+'9beb867f-780e-44aa-9f6f-6b986cd12b3f/vde1861f3_kim2012_31gene_signature.csv')

def load_cohorts():
    shared=[l.strip() for l in open(PATHS['shared']) if l.strip()]
    C={}
    for name,ek,ck in [('TCGA','tcga_expr','tcga_clin'),('CGGA693','c693_expr','c693_clin'),('CGGA325','c325_expr','c325_clin')]:
        e=pd.read_parquet(PATHS[ek]); c=pd.read_parquet(PATHS[ck])
        idc='sample' if 'sample' in c.columns else 'CGGA_ID'
        c=c.set_index(idc)
        c=c[c.primary_flag.astype(bool) & c.analysis_eligible_OS.astype(bool)]
        c=c[c.index.isin(e.columns)]
        e=e.loc[shared, c.index].astype(float)
        z=e.sub(e.mean(1),axis=0).div(e.std(1).replace(0,np.nan),axis=0).fillna(0)
        C[name]=dict(expr=e, z=z, clin=c)
    return C, shared

def ysurv(S): return Surv.from_arrays(S['OS'].astype(bool).values, S['OS.time'].astype(float).values)
def cindex(risk,S): return cic(S['OS'].astype(bool).values, S['OS.time'].astype(float).values, np.asarray(risk,float))[0]

def fast_unicox(Xdf, time, event, iters=8):
    t=np.asarray(time,float); d=np.asarray(event,float)
    o=np.argsort(-t, kind='stable'); t=t[o]; d=d[o]; X=np.asarray(Xdf,float)[o]
    n=len(t); last=np.zeros(n,int); i=0
    while i<n:
        j=i
        while j+1<n and t[j+1]==t[i]: j+=1
        last[i:j+1]=j; i=j+1
    beta=np.zeros(X.shape[1]); ev=d==1
    for _ in range(iters):
        w=np.exp(X*beta); S0=np.cumsum(w,0)[last]; S1=np.cumsum(w*X,0)[last]; S2=np.cumsum(w*X*X,0)[last]
        m=S1/S0; U=(X[ev]-m[ev]).sum(0); I=(S2/S0-m*m)[ev].sum(0)
        step=U/np.where(I>1e-12,I,1e-12); step=np.clip(step,-2,2); beta=beta+step
    w=np.exp(X*beta); S0=np.cumsum(w,0)[last]; S1=np.cumsum(w*X,0)[last]; S2=np.cumsum(w*X*X,0)[last]; m=S1/S0
    I=(S2/S0-m*m)[ev].sum(0); se=1/np.sqrt(np.where(I>1e-12,I,1e-12)); z=beta/se; p=2*stats.norm.sf(np.abs(z))
    return pd.DataFrame({'coef':beta,'HR':np.exp(beta),'se':se,'z':z,'p':p}, index=Xdf.columns)

class Wrap:
    def __init__(self, pred_fn, feats): self.pred_fn=pred_fn; self.feats=feats
    def predict(self, X): return np.asarray(self.pred_fn(X[self.feats]),float)

# ---------- selectors ----------
def sel_unicox(X,S,p=0.01,top=None):
    u=fast_unicox(X,S['OS.time'],S['OS']); g=u.index[u.p<p].tolist()
    if top: g=u.loc[g].sort_values('p').index[:top].tolist()
    return g
def _coxnet_cv(X,S,l1,n_alphas=50,seed=0):
    y=ysurv(S); m=CoxnetSurvivalAnalysis(l1_ratio=l1,n_alphas=n_alphas,alpha_min_ratio=0.05,max_iter=20000).fit(X,y)
    alphas=m.alphas_; kf=KFold(5,shuffle=True,random_state=seed); sc=np.zeros(len(alphas))
    for tr,te in kf.split(X):
        mm=CoxnetSurvivalAnalysis(l1_ratio=l1,alphas=alphas,max_iter=20000).fit(X.iloc[tr],y[tr])
        for i,a in enumerate(alphas):
            try: sc[i]+=cindex(mm.predict(X.iloc[te],alpha=a),S.iloc[te])
            except Exception: pass
    return alphas[int(np.argmax(sc))]
def sel_lasso(X,S):
    a=_coxnet_cv(X,S,1.0); m=CoxnetSurvivalAnalysis(l1_ratio=1.0,alphas=[a],max_iter=20000).fit(X,ysurv(S))
    return X.columns[m.coef_[:,0]!=0].tolist()
def sel_rsf(X,S,top=30):
    m=RandomSurvivalForest(n_estimators=500,min_samples_leaf=10,max_features='sqrt',n_jobs=8,random_state=0).fit(X,ysurv(S))
    r=permutation_importance(m,X,ysurv(S),n_repeats=5,random_state=0,n_jobs=8)
    imp=pd.Series(r.importances_mean,index=X.columns).sort_values(ascending=False)
    return imp.index[:top].tolist()
def sel_boost(X,S):
    m=ComponentwiseGradientBoostingSurvivalAnalysis(n_estimators=300,learning_rate=0.1,random_state=0).fit(X,ysurv(S))
    coef=pd.Series(m.coef_[1:],index=X.columns); return coef.index[coef!=0].tolist()
def _step_generic(X,S,direction,top=25,max_k=15):
    cand=sel_unicox(X,S,p=1,top=top)
    df=X[cand].copy(); df['T']=S['OS.time'].values; df['E']=S['OS'].values
    def aic(gs):
        if not gs: return np.inf
        c=CoxPHFitter(penalizer=0.001).fit(df[gs+['T','E']],'T','E'); return -2*c.log_likelihood_+2*len(gs)
    if direction in ('forward','both'):
        sel=[]; cur=np.inf
        while len(sel)<max_k:
            best=None
            for g in cand:
                if g in sel: continue
                try:a=aic(sel+[g])
                except Exception: continue
                if a<cur-1e-6 and (best is None or a<best[1]): best=(g,a)
            if best is None: break
            sel.append(best[0]); cur=best[1]
            if direction=='both' and len(sel)>2:
                improved=True
                while improved and len(sel)>1:
                    improved=False
                    for g in list(sel):
                        a=aic([s for s in sel if s!=g])
                        if a<cur-1e-6: sel.remove(g); cur=a; improved=True; break
        return sel
    sel=list(cand); cur=aic(sel)
    while len(sel)>1:
        best=None
        for g in sel:
            try:a=aic([s for s in sel if s!=g])
            except Exception: continue
            if a<cur-1e-6 and (best is None or a<best[1]): best=(g,a)
        if best is None: break
        sel.remove(best[0]); cur=best[1]
    return sel

# ---------- models ----------
def _coxnet_model(l1):
    def f(X,S):
        a=_coxnet_cv(X,S,l1); m=CoxnetSurvivalAnalysis(l1_ratio=l1,alphas=[a],max_iter=20000).fit(X,ysurv(S))
        w=Wrap(lambda Z:m.predict(Z),list(X.columns)); w.coef_=pd.Series(m.coef_[:,0],index=X.columns); w.alpha=a; return w
    return f
m_lasso=_coxnet_model(1.0); m_ridge=_coxnet_model(0.01)
def m_rsf(X,S):
    m=RandomSurvivalForest(n_estimators=500,min_samples_leaf=10,max_features='sqrt',n_jobs=8,random_state=0).fit(X,ysurv(S)); return Wrap(lambda Z:m.predict(Z),list(X.columns))
def m_gbm(X,S):
    m=GradientBoostingSurvivalAnalysis(n_estimators=300,learning_rate=0.05,max_depth=2,subsample=0.8,random_state=0).fit(X,ysurv(S)); return Wrap(lambda Z:m.predict(Z),list(X.columns))
def m_svm(X,S):
    m=FastSurvivalSVM(alpha=1.0,rank_ratio=1.0,max_iter=200,random_state=0).fit(X,ysurv(S)); return Wrap(lambda Z:m.predict(Z),list(X.columns))
def m_coxboost(X,S):
    m=ComponentwiseGradientBoostingSurvivalAnalysis(n_estimators=300,learning_rate=0.1,random_state=0).fit(X,ysurv(S))
    w=Wrap(lambda Z:m.predict(Z),list(X.columns)); w.coef_=pd.Series(m.coef_[1:],index=X.columns); return w
def m_plscox(X,S,ncomp=3):
    naf=NelsonAalenFitter().fit(S['OS.time'],S['OS']); H=naf.cumulative_hazard_at_times(S['OS.time']).values
    mart=S['OS'].values-H; k=min(ncomp,X.shape[1]); pls=PLSRegression(n_components=k,scale=False).fit(X.values,mart)
    T=pls.transform(X.values); cox=CoxPHSurvivalAnalysis(alpha=0.01).fit(T,ysurv(S))
    return Wrap(lambda Z:cox.predict(pls.transform(Z.values)),list(X.columns))
def m_superpc(X,S,thr_q=0.5,ncomp=1):
    u=fast_unicox(X,S['OS.time'],S['OS']); keep=u.index[u.z.abs()>=u.z.abs().quantile(thr_q)].tolist()
    if len(keep)<2: keep=list(X.columns)
    pca=PCA(n_components=min(ncomp,len(keep))).fit(X[keep].values); T=pca.transform(X[keep].values)
    cox=CoxPHSurvivalAnalysis(alpha=0.01).fit(T,ysurv(S))
    return Wrap(lambda Z:cox.predict(pca.transform(Z[keep].values)),list(X.columns))
def make_step_model(direction):
    def f(X,S):
        g=_step_generic(X,S,direction); g=g if len(g)>=1 else list(X.columns)
        c=CoxPHSurvivalAnalysis(alpha=0.01).fit(X[g],ysurv(S)); w=Wrap(lambda Z:c.predict(Z[g]),list(X.columns)); w.coef_=pd.Series(c.coef_,index=g); return w
    return f
MODELS={'Lasso':m_lasso,'Ridge':m_ridge, **{f'Enet[a={a/10:.1f}]':_coxnet_model(a/10) for a in range(1,10)},
        **{f'StepCox[{d}]':make_step_model(d) for d in ['forward','backward','both']},
        'survival-SVM':m_svm,'CoxBoost':m_coxboost,'plsRcox':m_plscox,'SuperPC':m_superpc,'GBM':m_gbm,'RSF':m_rsf}
SELECTORS={'RSF':sel_rsf,'Lasso':sel_lasso,'CoxBoost':sel_boost, **{f'StepCox[{d}]':(lambda d: lambda X,S:_step_generic(X,S,d))(d) for d in ['forward','backward','both']}}
STEP=[f'StepCox[{d}]' for d in ['forward','backward','both']]
GRID=[(None,m) for m in MODELS]
for s in ['RSF','Lasso','CoxBoost']: GRID+=[(s,m) for m in MODELS if m!=s]
for s in STEP: GRID+=[(s,m) for m in ['RSF','Lasso','Ridge','CoxBoost','plsRcox','SuperPC','GBM','survival-SVM']]
assert len(GRID)==101

def run_grid(Xt,St,ext,grid=GRID,log=None):
    sel_cache={}; rows=[]; fitted={}
    for s,m in grid:
        t0=time.time(); name=(f'{s}+{m}' if s else m)
        try:
            if s is None: feats=list(Xt.columns)
            else:
                if s not in sel_cache: sel_cache[s]=SELECTORS[s](Xt,St)
                feats=sel_cache[s]
            if len(feats)<2: raise ValueError('selector returned <2 features')
            mdl=MODELS[m](Xt[feats],St)
            nz=int((mdl.coef_!=0).sum()) if hasattr(mdl,'coef_') else len(feats)
            row={'combo':name,'selector':s or 'none','model':m,'n_feat_in':len(feats),'n_feat_nonzero':nz,'C_train':cindex(mdl.predict(Xt[feats]),St)}
            for n,X,S in ext: row[f'C_{n}']=cindex(mdl.predict(X[feats]),S)
            row['C_ext_mean']=np.mean([row[f'C_{n}'] for n,_,_ in ext]); row['status']='ok'; row['sec']=round(time.time()-t0,1)
            fitted[name]=(mdl,feats)
        except Exception as e:
            row={'combo':name,'selector':s or 'none','model':m,'status':f'failed: {type(e).__name__}: {str(e)[:80]}'}
        rows.append(row)
        if log:
            with open(log,'a') as fh: fh.write(f"{name}\t{row.get('C_ext_mean')}\t{row['status']}\n")
    return pd.DataFrame(rows), fitted, sel_cache

# ---------- validation ----------
def km_split(risk,S,label):
    risk=np.asarray(risk,float); med=np.median(risk); grp=np.where(risk>med,'high','low')
    df=pd.DataFrame({'T':S['OS.time'].values,'E':S['OS'].values,'g':grp,'risk':risk})
    lr=logrank_test(df['T'][df.g=='high'],df['T'][df.g=='low'],df.E[df.g=='high'],df.E[df.g=='low'])
    cph=CoxPHFitter().fit(df.assign(high=(df.g=='high').astype(int))[['T','E','high']],'T','E')
    cph2=CoxPHFitter().fit(df.assign(z=(df.risk-df.risk.mean())/df.risk.std())[['T','E','z']],'T','E')
    met={'label':label,'n':len(df),'events':int(df.E.sum()),'median_cut':med,'n_high':int((df.g=='high').sum()),
         'logrank_p':lr.p_value,'HR_high_vs_low':cph.hazard_ratios_['high'],'HR_lo':np.exp(cph.confidence_intervals_.iloc[0,0]),'HR_hi':np.exp(cph.confidence_intervals_.iloc[0,1]),
         'HR_perSD_unadj':cph2.hazard_ratios_['z'],'HR_perSD_lo':np.exp(cph2.confidence_intervals_.iloc[0,0]),'HR_perSD_hi':np.exp(cph2.confidence_intervals_.iloc[0,1]),'p_perSD':cph2.summary.loc['z','p'],
         'C_index':cindex(risk,S)}
    rows=[]
    for g in ['low','high']:
        kmf=KaplanMeierFitter().fit(df['T'][df.g==g],df.E[df.g==g]); sf=kmf.survival_function_; ci=kmf.confidence_interval_
        met[f'median_OS_days_{g}']=kmf.median_survival_time_
        for t,s in sf.itertuples():
            rows.append({'label':label,'group':g,'time_days':t,'surv':s,'ci_low':ci.loc[t].iloc[0],'ci_high':ci.loc[t].iloc[1],'n_risk':int((df['T'][df.g==g]>=t).sum())})
    return met, pd.DataFrame(rows), grp

def time_auc(risk_tr,S_tr,risk,S,years=(1,3,5)):
    ytr=ysurv(S_tr); yte=ysurv(S); out=[]; tmax=S['OS.time'].max()
    for y in years:
        t=365.25*y
        if t>=tmax*0.98: out.append({'year':y,'AUC':np.nan,'note':'beyond follow-up'}); continue
        try: a,_=cumulative_dynamic_auc(ytr,yte,np.asarray(risk,float),[t]); out.append({'year':y,'AUC':float(a[0])})
        except Exception as e: out.append({'year':y,'AUC':np.nan,'note':str(e)[:60]})
    return out

def calibration_decile(risk,S,label,t=3*365.25,nbins=10):
    df=pd.DataFrame({'T':S['OS.time'].values,'E':S['OS'].values,'risk':np.asarray(risk,float)})
    c=CoxPHFitter().fit(df,'T','E')
    df['pred']=c.predict_survival_function(df[['risk']],times=[t]).T.iloc[:,0].values
    df['dec']=pd.qcut(df.risk.rank(method='first'),nbins,labels=False)+1
    rows=[]
    for d,g in df.groupby('dec'):
        kmf=KaplanMeierFitter().fit(g['T'],g.E)
        obs=float(kmf.survival_function_at_times(t).iloc[0]); ci=kmf.confidence_interval_survival_function_
        idx=ci.index[ci.index<=t]; lo,hi=(ci.loc[idx[-1]].values if len(idx) else (np.nan,np.nan))
        rows.append({'label':label,'time_years':t/365.25,'decile':d,'n':len(g),'events':int(g.E.sum()),'pred_surv':g.pred.mean(),'obs_surv':obs,'obs_lo':lo,'obs_hi':hi,'mean_risk':g.risk.mean()})
    return pd.DataFrame(rows)

def dca(risk,S,label,t=3*365.25,covars=None,model_name='signature'):
    df=pd.DataFrame({'T':S['OS.time'].values,'E':S['OS'].values,'risk':np.asarray(risk,float)}); cols=['risk']
    if covars is not None:
        cv=covars.reset_index(drop=True); df=pd.concat([df,cv],axis=1); cols+=list(cv.columns); df=df.dropna()
    c=CoxPHFitter(penalizer=0.001).fit(df[['T','E']+cols],'T','E')
    p=1-c.predict_survival_function(df[cols],times=[t]).T.iloc[:,0].values
    kmall=KaplanMeierFitter().fit(df['T'],df.E); prev=1-float(kmall.survival_function_at_times(t).iloc[0])
    rows=[]
    for thr in np.arange(0.01,0.90,0.01):
        pos=p>=thr
        if pos.sum()<5: nb=0.0
        else:
            km=KaplanMeierFitter().fit(df['T'][pos],df.E[pos]); ev=1-float(km.survival_function_at_times(t).iloc[0])
            frac=pos.mean(); nb=ev*frac-(1-ev)*frac*thr/(1-thr)
        nb_all=prev-(1-prev)*thr/(1-thr)
        rows.append({'label':label,'model':model_name,'time_years':t/365.25,'threshold':round(thr,2),'NB_model':nb,'NB_all':nb_all,'NB_none':0.0})
    return pd.DataFrame(rows)

def prep_covars(clin):
    d=pd.DataFrame(index=clin.index)
    d['age']=clin.age.astype(float)
    d['grade_num']=clin.grade.map({'G2':2,'G3':3,'G4':4})
    d['IDH_wt']=clin.IDH_status.map({'Wildtype':1,'Mutant':0})
    d['codel']=clin.codel_1p19q.map({'Codel':1,'Non-codel':0})
    d['MGMT_meth']=clin.MGMT_status.map({'Methylated':1,'Unmethylated':0})
    d['RT']=clin.radiation_therapy.map({'yes':1,'no':0})
    d['chemo']=clin['chemotherapy_TMZ'].map({'yes':1,'no':0}) if 'chemotherapy_TMZ' in clin else np.nan
    d['sex_male']=clin.sex.map({'Male':1,'Female':0})
    return d

def mv_cox(risk,S,cov,label,adj,penalizer=0.0):
    r=np.asarray(risk,float)
    df=pd.concat([pd.DataFrame({'T':S['OS.time'].values,'E':S['OS'].values,'risk_z':(r-r.mean())/r.std()},index=S.index),cov.loc[S.index,adj]],axis=1).dropna()
    c=CoxPHFitter(penalizer=penalizer).fit(df,'T','E'); s=c.summary; rows=[]
    for v in s.index:
        rows.append({'label':label,'n':len(df),'events':int(df.E.sum()),'variable':v,'HR':s.loc[v,'exp(coef)'],'HR_lo':s.loc[v,'exp(coef) lower 95%'],'HR_hi':s.loc[v,'exp(coef) upper 95%'],'p':s.loc[v,'p'],'adjusted_for':'+'.join(adj)})
    return pd.DataFrame(rows), c
