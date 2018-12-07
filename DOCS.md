# gtx
Genetics ToolboX R package
# Mendelian Randomization Analysis

## Introduction
Mendelian Randomization analysis consists of several steps:
1. Extract the instruments (SNPs) related to an exposure e.g., a risk factor such as smoking.
2. (Optional) Prune the instruments to ensure that they are independent (automated for the datasets with rsIDs, otherwise has to be completed by the user).
3. Extract the corresponding SNPs from the chosen outcome(s) e.g., a heart disease.
4. Format both datasets so that they can be used within the [MR-Base platform](http://www.mrbase.org), and ensure that all required variables were supplied (if a custom file is used).
5. Run MR analysis and follow-up analyses (such as the heterogeneity, pleiotropy and single SNP tests).

## GWAS search
This method extracts a GWAS meta-information such as its analysis ID and samplesize using a keyword/keyphrase e.g.,
```R
info<-extract_study_info('body mass')
```

## Instrument extraction (exposure)
The exposure instruments can be extracted by either using the chosen GWASs (_analysis_ or _analysis+entity_); or searching by a keyword/keyphrase such as 'body mass'. In the latest case, all studies which have a specified keyphrase will be considered for instrument extraction. If all studies should be looked up then:
```R
exp<-extract_exposure()
```

### Instrument extraction from the chosen GWASs
The instruments can be extracted from the chosen GWASs if a dataframe __df__ of either one column (_analysis_) or two columns (_analysis_, _entity_) is provided such as
```R
exp<-extract_exposure(analyses=df)
```
This method returns by default chromosome positions as unique SNP identifiers. If rsIDs should be returned as well then:
```R
exp<-extract_exposure(analyses=df,rsid=TRUE)
```
The default p-value threshold is _5e-08_. If it has to be changed then:
```R
exp<-extract_exposure(p=5e-09,analyses=df)
```
### Instrument extraction using a keyword/keyphrase
Note that all above options are applicable in this scenario as well such as:
```R
exp<-extract_exposure(p=5e-09,str='body mass',rsid=TRUE)
```
or with default settings:
```R
exp<-extract_exposure(str='body mass')
```
