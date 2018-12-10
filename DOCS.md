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

## Instrument (exposure) extraction
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
## Outcome extraction
The instruments chosen for the exposure then have to be found in the outcome GWAS. There are two options:
1. The exposure instruments come from the external file (e.g. the paper's supplementary table) or database (e.g. MR-Base).
2. The exposure instruments are extracted from the internal database.

__Note!__ It is advisable to use the internal database (upload all new data to the internal database) if the large-scale MR analysis with >100 instruments has to be conducted. The first option is better to be used for brief lookups only as it might lead to the time-consuming queries and long waiting times.

### Outcome extraction using the external input
The external input (e.g. a file) has to be read into the R dataframe, e.g. __snps__, and should have either 2 columns with _position_ (1st) and _choromosome_ (2nd), or one column with rsIDs (the name of the column for rsIDs should be either "rsid" or "snp"). To ensure that the data structure used is a dataframe:
```R
snps<-as.data.frame(snps)
```
__Note!__ By default, the method considers extracting SNPs identified by their chromosome positions (if rsIDs are used instead, then choose __search__ = "rsid"). The list of outcome studies to be considered should be specified as a dataframe with _analysis_ (1st column) and possibly _entity_ (2nd column), i.e. __df.out__. Hence, if snps are chromosome positions, then:
```R
out<-extract_outcome(snp_list=snps,analyses_out=df.out)
```
If snps are rsIDs, then:
```R
out<-extract_outcome(snp_list=snps,analyses_out=df.out,search="rsid")
```
__Note!__ By default, the method returns the dataframe with SNPs uniquely identified by their chromosome positions but not their rsIDs. If rsIDs should be returned, then:
```R
out<-extract_outcome(snp_list=snps,analyses_out=df.out,rsid=TRUE)
```
### Outcome extraction using the internal database
In this scenario, there is no need to specify a SNP list but the studies to be extracted from the database for both the exposure and outcome as in the previous examples, that is the dataframe __df.exp__ for the exposure and __df.out__ for the outcome with 2 columns of _analysis_ and _entity_ (if applicable) where _entity_ could be an empty column. If rsIDs should be returned, then __rsid__=TRUE. Hence, the method would be:
```R
out<-extract_outcome(analyses_exp=df.exp,analyses_out=df.out,rsid=TRUE)
```
__Note!__ If __rsid__=TRUE, the response could be slower; that is choose this option if needed.

## Data formatting
This method formats all datasets which were extracted from the internal database automatically. It is advisable to format any external data by using *format_data* method from the [TwoSampleMR](https://mrcieu.github.io/TwoSampleMR/) package.
```R
exposure_data<-format_for_mr(data=df) #  for the exposure (default)
outcome_data<-format_for_mr(data=df,type="outcome") # for the outcome
```
If rsIDs has to be used as unique SNP identifiers instead of chromosome positions & alleles, then:
```R
exposure_data<-format_for_mr(data=df,rsid=TRUE) #  for the exposure (default)
outcome_data<-format_for_mr(data=df,type="outcome",rsid=TRUE) # for the outcome
```
## The univariable MR and follow-up analyses
This method runs the MR analysis and follow-up analyses such as the heterogeneity and pleiotropy tests as well as a single SNP MR analysis (all implemented as separate methods in [TwoSampleMR](https://mrcieu.github.io/TwoSampleMR/) but don't require a connection to [MR-Base](http://www.mrbase.org)), using any unique SNP identifier that is it could be chromosome positions & alleles (pos:chrom_ref_alt) or rsIDs or potentially anything else.

__Note!__ It is important to ensure that both the exposure and outcome use same unique SNP identifiers before applying this method. If alleles are used as a part of SNP identification, then they should be aligned in the same way in both the exposure and outcome at the formatting stage.

__Note!__ If rsIDs are used as SNP identifiers instead of any other options, then the exposure instruments will be LD-pruned automatically. In all other cases, the user has to ensure that the instruments are independent before using them with this method.
```R
res<-run_mr_gsk(exposure_data,outcome_data) # default, runs MR and all other analyses
res<-run_mr_gsk(exposure_data,outcome_data,rsid=TRUE) # runs MR and all other analyses, using rsIDs for SNPs
```
It is possible to 'force' this method to skip LD-pruning when rsIDs are used (for example, if this is too time-consuming or there are connection issues with MR-Base):
```R
res<-run_mr_gsk(exposure_data,outcome_data,rsid=TRUE,clump=FALSE)
```
__Note!__ If __clump__=TRUE, then it is necessary to establish a connection to [MR-Base](http://www.mrbase.org). If you have never used the TwoSampleMR package, then the 'mrbase.auth' token has to be generated first. Please, generate this token on your personal machine rather than a server by installing TwoSampleMR and running __available_outcomes()__ in R. You will be prompted to login into your Gmail account and after you have logged in, then the token will be generated and appear in the working directory. This token can then be copied to the server or anywhere else and will authorise your access to MR-Base.

It is also possible to skip any follow-up analyses, where __hetero__ and __pleio__ stand for the heterogeneity and pleiotropy tests, while __single__ stands for the single SNP MR analysis:
```R
res<-run_mr_gsk(exposure_data,outcome_data,hetero=FALSE,pleio=FALSE,single=FALSE)
```
The output of this method is the list of dataframes __res__:
1. MR results
2. Heterogeneity analysis
3. Pleiotropy analysis
4. Single SNP MR analysis
5. Harmonisation data

For example, if you would like to access the MR results, please, choose __res[[1]]__.


