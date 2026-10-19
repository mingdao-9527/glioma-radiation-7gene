"""IDH-mutant G2-3 track — consolidated analysis script (reproduces the cells run in session; see README_idhmut.md).
Run inside the 'glioma' conda env. Steps: population -> 7-gene validation -> RT x risk interaction (Cox, IPTW, pooled)
-> 101-combination grid (sig_lib.run_grid, RSF patched to n_jobs=1) -> comparators & random-gene null -> new16 validation."""
import importlib.util, pandas as pd, numpy as np, warnings, os, time, statsmodels.api as sm
warnings.filterwarnings('ignore')
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
from sksurv.linear_model import CoxPHSurvivalAnalysis
from sksurv.ensemble import RandomSurvivalForest
from sklearn.inspection import permutation_importance
base='/Users/laowang/.claude-science/orgs/4321deaf-044c-4bf9-be37-0197b924b0df/artifacts/proj_39d29643281e/'
spec=importlib.util.spec_from_file_location('sig_lib', base+'a3c7b604-ac9e-4584-9df5-760d2149cced/v4460e772_sig_lib.py'); sl=importlib.util.module_from_spec(spec); spec.loader.exec_module(sl)
C, shared = sl.load_cohorts()
coef=pd.read_csv(base+'f64ebd7c-b360-4294-9b2b-df6373050193/v4d3b0618_signature_coefficients.csv').set_index('gene')['coef']
risk_all=pd.read_csv(base+'93fac55d-f9e4-48af-a7e3-09b610a51a31/vb8eae608_risk_scores_all_cohorts.csv').set_index('sample')
pc=pd.read_csv(base+'9f56299d-5f1d-49ed-ab9e-894818e2ac53/v7d2e230d_pool_candidates.csv')
os.makedirs('idhmut',exist_ok=True)
# ---- 1. population + two standardisations of the 7-gene score
pop={}
for k,d in C.items():
    c=d['clin']; s=c[(c.IDH_status=='Mutant')&(c.grade.isin(['G2','G3']))].copy()
    s['risk_cohortz']=risk_all.loc[s.index,'risk']
    e=d['expr'].loc[coef.index, s.index]; z=e.sub(e.mean(1),axis=0).div(e.std(1),axis=0); s['risk_subz']=(z.mul(coef,axis=0)).sum(0)
    pop[k]=s
def cov(s):
    d=pd.DataFrame(index=s.index); d['age']=s.age.astype(float); d['grade3']=(s.grade=='G3').astype(int)
    d['codel']=s.codel_1p19q.map({'Codel':1,'Non-codel':0}); d['MGMT_meth']=s.MGMT_status.map({'Methylated':1,'Unmethylated':0})
    d['RT']=s.radiation_therapy.map({'yes':1,'no':0}); d['TMZ']=s['chemotherapy_TMZ'].map({'yes':1,'no':0}) if 'chemotherapy_TMZ' in s else np.nan
    d['sex_male']=(s.sex=='Male').astype(int); return d
COV={k:cov(s) for k,s in pop.items()}
# ---- 2. validation metrics: sl.km_split (within-subgroup median), sl.time_auc(3,5,10 y), sl.mv_cox with adj sets
ADJ=[('age+grade+codel',['age','grade3','codel']),('age+grade+codel+MGMT',['age','grade3','codel','MGMT_meth']),('age+grade+codel+MGMT+TMZ',['age','grade3','codel','MGMT_meth','TMZ']),('age+grade+codel+RT',['age','grade3','codel','RT']),('age+grade+codel+MGMT+TMZ+RT',['age','grade3','codel','MGMT_meth','TMZ','RT'])]
# (loop over cohort x standardisation x subset all/RT_yes/RT_no; TCGA additionally PFI and DSS) — see session cell; outputs idhmut_validation_metrics.csv, idhmut_km_data.csv, idhmut_adjusted_cox.csv, idhmut_time_auc.csv
# ---- 3. RT x risk interaction
def build(k,s,endpoint='OS'):
    df=pd.DataFrame({'T':s[endpoint+'.time'].astype(float),'E':s[endpoint].astype(float)},index=s.index)
    r=s.risk_cohortz; df['risk_z']=(r-r.mean())/r.std(); df['risk_high']=(r>r.median()).astype(int)
    df=pd.concat([df,COV[k][['age','grade3','codel','RT','MGMT_meth','TMZ','sex_male']]],axis=1); df['cohort']=k; return df
def fit_cox(df,cols,strata=None,weights=None,robust=False):
    d=df[['T','E']+cols+([strata] if strata else [])+([weights] if weights else [])].dropna()
    return CoxPHFitter().fit(d,'T','E',strata=strata,weights_col=weights,robust=robust), d
# per cohort: fit_cox(d,[rv,'RT','RTxrisk']+adj) for rv in risk_z/risk_high, adj in none / age+grade+codel
# IPTW: ps = Logit(RT ~ age+grade3+codel+risk_z+sex_male); w = pRT/ps (treated) or (1-pRT)/(1-ps); clip at 99th pct; fit_cox(..., weights='w', robust=True)
# within-stratum RT effect: fit_cox(d[d.risk_high==g],['RT']+adj); KM per arm -> idhmut_rt_km_data.csv
# pooled: fit_cox(pd.concat(all cohorts),[rv,'RT','RTxrisk']+adj, strata='cohort') -> idhmut_pooled_cox.csv
# ---- 4. de novo grid
def sel_rsf1(X,S,top=30):
    m=RandomSurvivalForest(n_estimators=500,min_samples_leaf=10,max_features='sqrt',n_jobs=1,random_state=0).fit(X,sl.ysurv(S))
    r=permutation_importance(m,X,sl.ysurv(S),n_repeats=5,random_state=0,n_jobs=1)
    return pd.Series(r.importances_mean,index=X.columns).sort_values(ascending=False).index[:top].tolist()
def m_rsf1(X,S):
    m=RandomSurvivalForest(n_estimators=500,min_samples_leaf=10,max_features='sqrt',n_jobs=1,random_state=0).fit(X,sl.ysurv(S)); return sl.Wrap(lambda Z:m.predict(Z),list(X.columns))
sl.SELECTORS['RSF']=sel_rsf1; sl.MODELS['RSF']=m_rsf1
sym=dict(zip(pc.gene,pc.cohort_symbol)); adj=pc[pc.adj_concordant_ge2].gene.tolist()
genes245=[sym[g] for g in adj]; genes64=[sym[g] for g in pc[(~pc.grade_proxy)&(~pc.IDH_proxy)&pc.adj_concordant_ge2].gene]
Z={}
for k,s in pop.items():
    e=C[k]['expr'].loc[:, s.index]; Z[k]=e.sub(e.mean(1),axis=0).div(e.std(1).replace(0,np.nan),axis=0).fillna(0).T
rt=lambda k: pop[k].index[pop[k].radiation_therapy=='yes']
TR,E1,E2=rt('TCGA'),rt('CGGA693'),rt('CGGA325')
St=pop['TCGA'].loc[TR,['OS','OS.time']].astype(float); S1=pop['CGGA693'].loc[E1,['OS','OS.time']].astype(float); S2=pop['CGGA325'].loc[E2,['OS','OS.time']].astype(float)
results={}
for lab,gs in [('all245',genes245),('nonproxy64',genes64)]:
    lb,fitted,selc=sl.run_grid(Z['TCGA'].loc[TR,gs],St,[('CGGA693',Z['CGGA693'].loc[E1,gs],S1),('CGGA325',Z['CGGA325'].loc[E2,gs],S2)])
    lb['feature_set']=lab; results[lab]=(lb,fitted,selc)
LB=pd.concat([results['all245'][0],results['nonproxy64'][0]],ignore_index=True)
best=LB[LB.feature_set=='all245'].sort_values('C_ext_mean',ascending=False).iloc[0]; mdl,feats=results['all245'][1][best.combo]
cf=mdl.coef_[mdl.coef_!=0]; newgenes=cf.index.tolist()
ns={k:(Z[k][newgenes]*cf.values).sum(1) for k in pop}   # linear predictor (= mdl.predict up to a constant)
# ---- 5. comparators: 7-gene (cohort z / subgroup z), clinical age+grade3+codel Cox fitted on TCGA RT IDHmut, clinical+signature, paired bootstrap (1000)
# ---- 6. random-gene null: 1000 random sets (size 16 and 7) from `shared`, and 1000 of size 16 from genes245; ridge CoxPHSurvivalAnalysis(alpha=1.0) on TCGA RT IDHmut; mean external C
# ---- 7. new16 validation across all/RT/noRT subsets and RT x new16 interaction (CGGA out-of-sample; TCGA RT arm in-sample)
