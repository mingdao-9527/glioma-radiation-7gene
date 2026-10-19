import pandas as pd, numpy as np, warnings; warnings.filterwarnings('ignore')
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines.utils import concordance_index
A='/Users/laowang/.claude-science/orgs/4321deaf-044c-4bf9-be37-0197b924b0df/artifacts/proj_39d29643281e'
risk=pd.read_csv(f'{A}/93fac55d-f9e4-48af-a7e3-09b610a51a31/vb8eae608_risk_scores_all_cohorts.csv')
# ---------- CGGA multivariate (no CNV available) ----------
rows=[]
def cox(df,adj,label,pop,T,E):
    r=df.risk.values; df=df.assign(risk_z=(r-np.nanmean(r))/np.nanstd(r,ddof=1))
    X=df[[T,E,'risk_z']+adj].dropna(); X=X.loc[:,X.nunique()>1]
    c=CoxPHFitter().fit(X,T,E); s=c.summary
    return [dict(label=label,population=pop,endpoint=E,n=len(X),events=int(X[E].sum()),variable=v,HR=s.loc[v,'exp(coef)'],HR_lo=s.loc[v,'exp(coef) lower 95%'],HR_hi=s.loc[v,'exp(coef) upper 95%'],p=s.loc[v,'p'],adjusted_for='+'.join([a for a in adj if a in X.columns])) for v in s.index], c, X
for coh in ['CGGA693','CGGA325']:
    d=risk[risk.cohort==coh].copy(); d['T']=d['OS.time'].astype(float); d['E']=d['OS'].astype(float); d=d.dropna(subset=['T','E']); d=d[d['T']>0]
    d['grade_num']=d.grade.map({'G2':2,'WHO II':2,'G3':3,'WHO III':3,'G4':4,'WHO IV':4}); d['IDH_wt']=d.IDH_status.map({'Wildtype':1,'Mutant':0}); d['codel']=d.codel_1p19q.map({'Codel':1,'Non-codel':0})
    d['MGMT_meth']=d.MGMT_status.map({'Methylated':1,'Unmethylated':0}); d['RT01']=d.RT.map({'yes':1,'no':0}); d['TMZ']=d.chemotherapy_TMZ.map({'yes':1,'no':0}); d['sex_male']=d.sex.map({'Male':1,'Female':0})
    print(coh, d.grade.unique(), d.chemotherapy_TMZ.value_counts(dropna=False).to_dict())
    for pop,sub,adj in [('all_primary',d,['age','sex_male','grade_num','IDH_wt','codel','MGMT_meth','RT01','TMZ']),('RT',d[d.RT=='yes'],['age','sex_male','grade_num','IDH_wt','codel','MGMT_meth','TMZ']),
                        ('IDHmut_noncodel_astro',d[(d.IDH_status=='Mutant')&(d.codel_1p19q=='Non-codel')],['age','sex_male','grade_num','MGMT_meth','RT01','TMZ']),('IDHwt',d[d.IDH_status=='Wildtype'],['age','sex_male','grade_num','MGMT_meth','RT01','TMZ'])]:
        try: r,_,_=cox(sub,adj,'CGGA_existing_covariates_plusTMZ (no CNV/TERTp available)',f'{coh}_{pop}','T','E'); rows+=r
        except Exception as e: print('FAIL',coh,pop,e)
cg=pd.DataFrame(rows); cg.to_csv('molgrade/multivariate_cox_cgga_tmz.csv',index=False)
print(cg[cg.variable=='risk_z'][['population','n','events','HR','HR_lo','HR_hi','p']].round(4).to_string())
# ---------- TCGA secondary endpoints ----------
w=pd.read_csv('molgrade/tcga_who2021.csv'); trt=pd.read_csv('molgrade/tcga_treatment_from_gdc.csv')
w['patient12']=w['sample'].str[:12]; w=w.merge(trt[['patient','TMZ_any']],left_on='patient12',right_on='patient',how='left'); w['TMZ_any']=w.TMZ_any.fillna('unknown')
w['sex_male']=w.sex.map({'Male':1,'Female':0}); w['MGMT_meth']=w.MGMT_status.map({'Methylated':1,'Unmethylated':0}); w['RT01']=w.RT.map({'yes':1,'no':0})
w['TMZ_yes']=(w.TMZ_any=='yes').astype(int); w['TMZ_unknown']=(w.TMZ_any=='unknown').astype(int)
w['grade_num']=w.grade.map({'G2':2,'G3':3,'G4':4}); w['IDH_wt']=w.IDH_status.map({'Wildtype':1,'Mutant':0}); w['codel']=w.codel_1p19q.map({'Codel':1,'Non-codel':0})
w['cls_astro']=(w.who2021_class=='Astrocytoma, IDH-mutant').astype(int); w['cls_gbm']=(w.who2021_class=='Glioblastoma, IDH-wildtype').astype(int); w['cls_nec']=(w.who2021_class=='IDH-wildtype lower-grade, NEC').astype(int)
w['classified']=w.who2021_class.isin(['Astrocytoma, IDH-mutant','Glioblastoma, IDH-wildtype','IDH-wildtype lower-grade, NEC','Oligodendroglioma, IDH-mutant & 1p/19q-codeleted'])
w['who_grade']=w.who2021_grade.astype(float); w['CDKN2A']=w.CDKN2AB_homdel; w['plus7minus10']=w.chr7gain_chr10loss; w['EGFRamp']=w.EGFR_amp
sec=[]; km=[]
for ep in ['OS','PFI','DSS']:
    T,E=f'{ep}.time',ep
    for pop,sub in [('all_primary',w[w.classified]),('RT',w[w.classified&(w.RT=='yes')]),('IDHmut_astrocytoma',w[w.who2021_class=='Astrocytoma, IDH-mutant']),('IDHwt',w[w.IDH_status=='Wildtype']),('oligodendroglioma',w[w.who2021_class=='Oligodendroglioma, IDH-mutant & 1p/19q-codeleted'])]:
        sub=sub.dropna(subset=[T,E]); sub=sub[sub[T]>0]
        strat=pop in ('IDHmut_astrocytoma','IDHwt','oligodendroglioma')
        specs={'univariable':[],'adj_original(age,sex,IDH,codel,hist grade,MGMT,RT,TMZ)':['age','sex_male','grade_num','IDH_wt','codel','MGMT_meth','RT01','TMZ_yes','TMZ_unknown'],
               'adj_WHO2021(age,sex,class,integrated grade,CDKN2A,+7/-10,EGFR,MGMT,RT,TMZ)':['age','sex_male','cls_astro','cls_gbm','cls_nec','who_grade','CDKN2A','plus7minus10','EGFRamp','MGMT_meth','RT01','TMZ_yes','TMZ_unknown']}
        for lab,adj in specs.items():
            a=[x for x in adj if not (strat and x in ['cls_astro','cls_gbm','cls_nec','IDH_wt','codel'])]
            if pop=='RT': a=[x for x in a if x!='RT01']
            if sub[E].sum()<10: continue
            try:
                r,c,X=cox(sub,a,lab,f'TCGA_{pop}',T,E)
                rz=[x for x in r if x['variable']=='risk_z'][0]
                # C-index of signature alone (univariable) within the fitted sample
                rz['C_index_signature_alone']=concordance_index(X[T],-X.risk_z,X[E]); rz['C_index_model']=concordance_index(X[T],-c.predict_partial_hazard(X),X[E])
                sec.append(rz)
            except Exception as e: print('FAIL',ep,pop,lab,e)
        # KM by risk group (all-primary median; RT-median for RT pop)
        grp='risk_group_RTmedian' if pop=='RT' else 'risk_group_allmedian'
        g=sub.dropna(subset=[grp])
        if g[E].sum()>=10 and g[grp].nunique()==2:
            hi=g[g[grp]=='high']; lo=g[g[grp]=='low']
            lr=logrank_test(hi[T],lo[T],hi[E],lo[E])
            for name,gg in [('high',hi),('low',lo)]:
                kmf=KaplanMeierFitter().fit(gg[T]/365.25,gg[E])
                med=kmf.median_survival_time_
                s=kmf.survival_function_.reset_index(); s.columns=['time_years','survival']
                ci=kmf.confidence_interval_.reset_index(); s['ci_lo']=ci.iloc[:,1].values; s['ci_hi']=ci.iloc[:,2].values
                s['endpoint']=ep; s['population']=f'TCGA_{pop}'; s['risk_group']=name; s['n']=len(gg); s['events']=int(gg[E].sum()); s['median_years']=med; s['logrank_p']=lr.p_value; s['group_definition']=grp
                km.append(s)
sec=pd.DataFrame(sec); sec.to_csv('molgrade/secondary_endpoints.csv',index=False)
pd.concat(km).to_csv('molgrade/km_secondary_endpoints.csv',index=False)
pd.set_option('display.width',250)
print(sec[['endpoint','population','label','n','events','HR','HR_lo','HR_hi','p','C_index_signature_alone','C_index_model']].round(4).to_string())
kk=pd.concat(km).groupby(['endpoint','population','risk_group']).agg(n=('n','first'),events=('events','first'),median_years=('median_years','first'),logrank_p=('logrank_p','first')).reset_index(); print(kk.round(4).to_string())
