# TESS-ExoClass (TEC)
TESS exoplanet detection filter, crossmatch, classifier, and ranker.

### Problem:
Your transiting planet detection pipeline results in a huge pile of detections.  Now What?  Do I really have to look at all these manually!?

### Answer:
No!  Algorithmic filtering and classifying is available with TESS-ExoClass (TEC).

### Authors:
Christopher J Burke (MIT), Jeffrey Coughlin (SETI/NASA Ames), Susan Mullally (STScI), Fergal Mullally (SETI), Jessie Christiansen (NExSci)

### Results:
Sector Number | Detection Source | Result Tables | Summary PDFs
------------- | ---------------- | ------------- | ------------
3 | NASA Ames SPOC TCE | [Tables](https://www.dropbox.com/sh/lakbewewji4ieap/AABT-7AquQRV4u1RwBwhvz_Ea?dl=0) | [PDFs](https://www.dropbox.com/sh/pd0lsf5u5uk7qil/AABc2UrZcdI75VExtkDmkLqpa?dl=0)
4 | NASA Ames SPOC TCE | [Tables](https://www.dropbox.com/sh/crkk010kl0rjp9z/AACfSft-gQckJZFTd-4F9_Yta?dl=0) | [PDFs](https://www.dropbox.com/sh/gq7shwdw2pi1xjv/AADQCc1H9yIrd99mwPdADynoa?dl=0)
1-3 | NASA Ames SPOC TCE | [Tables](https://www.dropbox.com/sh/im6pg6uepelt0dw/AADm7VH_PYKriETXda080yIRa?dl=0) | [PDFs](https://www.dropbox.com/sh/sf2qiw708z6fgvb/AAALmpgUqP6DrTf2VtlYgF7Fa?dl=0)
5 | NASA Ames SPOC TCE | [Tables](https://www.dropbox.com/sh/0ejf3rc0fneizq7/AAAzti3PoskfWLlAcTrMbHzza?dl=0) | [PDFs](https://www.dropbox.com/sh/7fvfusjsa00gua1/AAC--PPZHcPfrSrp3HaJYq5ja?dl=0)
6 | NASA Ames SPOC TCE | [Tables](https://www.dropbox.com/sh/hqvx611hfwvqvlx/AAAs60DnbHbw7Q-mWabj2nRRa?dl=0) | [PDFs](https://www.dropbox.com/sh/wuvqhal09wgeqo4/AACGgNO9_LFdTLlc68mM1f4Aa?dl=0)

### Description:
TEC gets you from your pile of signal detections to your best transiting planet candidate science faster.  TEC builds off the highly successful Kepler Robovetter classification process ([Coughlin et al. 2016](http://adsabs.harvard.edu/abs/2016ApJS..224...12C), [Thompson et al. 2018](http://adsabs.harvard.edu/abs/2018ApJS..235...38T)) as well as the [DAVE](http://keplertcert.seti.org/DAVE/) K2 vetting tools ([Kostov et al. 2019](http://adsabs.harvard.edu/abs/2019arXiv190107459K)).  TEC is actually neither of those more mature efforts (yet).  TEC is actually the foundation underlying those classifiers.  In other words TEC is the database of attributes/metrics that generically feed any kind of machine learning or classification process.  For instance, the Kepler Robovetter examined a database of ~50 metrics in order for it to carry out its decision tree based classifications.  TEC is currently focused on developing the database of attributes/metrics that can provide more sophisticated classifications using machine learning techniques.  TEC currently uses a crude decision tree like classifier that is tuned manually in order to match the efficacy of the current manual vetting effort to build the MIT TOI alerts.

### What does TEC enable:
- Observing the best candidates: TEC provides a filtered and ranked list of the most viable detections to get you to the 'gem' planet candidates quickly.  In practice, TEC has demonstrated that the number of viable detections can be reduced by a factor of six from the original detections.
- Machine learning attributes and training sets: All machine learning algorithms begin with a database of attributes/metrics.  You can use TEC to enhance the current set of attributes that you have developed or you can use the currently crude classifications from TEC to enhance your training set.
- Statistical Validation: Do you want to statistically validate your lowest signal to noise ratio (SNR) detections?  Then you need to quantify the reliability of your detections and classification.  [Mullally et al. (2018)](http://adsabs.harvard.edu/abs/2018AJ....155..210M) and [Burke et al. (2019)](http://adsabs.harvard.edu/abs/2019arXiv190100506B) recently have shown that important Kepler statistical validations (such as Kepler-452b and Kepler-186f) are not as statistically significant as previously thought without taking into account the possibility that they could be false alarm/systematic detections.  TEC enables automating the classification thus one can provide TEC examples of systematics in order to demonstrate whether reliability is an issue for the lowest SNR detections.
- Planet Occurrence Rates: Are you one of the brave souls that wants to calculate unbiased planet occurrence rates based on TESS data?  You need automated classifications in order to quantify how this step influences the recovered planets.
- Your Pipeline Here: I have my own planet search pipeline, and my own large pile of detections.  Help!  Without too much work providing a flux-time series and detection ephemeris you can run TEC for your own pipeline.

### STATUS:
TEC is a very, very much work in progress from a coding and workflow standpoint.  However, TEC is producing filtered and ranked list of candidates from the TESS produced NASA Ames SPOC Threshold Crossing Events (TCEs) along with summary documents for each detection that TEC deems worthy of attention.  So, even if you ignore the codes, TEC currently can provide a high fidelity list of the most promising planet candidates for observing follow-up, machine learning training sets, or benchmark to your own efforts.  The TOI alerts (focused on providing the most reliable candidates for the TESS prime mission) are a heavily filtered set of candidates.  TEC provides a more thorough search to lower SNR levels and classifications for all TCEs from the NASA Ames pipeline.

From an attribute/metric standpoint most of the flux-level metrics from the Kepler Robovetter are available with TEC.  The notable exceptions are LPP ([Thompson et al. 2015](http://adsabs.harvard.edu/abs/2015ApJ...812...46T)) and Marshall ([Mullally et al. 2016](http://adsabs.harvard.edu/abs/2016PASP..128g4502M)).  Those missing metrics do have publicly available algorithms under the [DAVE](https://github.com/barentsen/dave) project and will potentially be implemented in future TEC updates.  Pixel-level centroid analysis does not exist in TEC yet, and won't for sometime.  There is an effort to independently calculate in-transit difference images, but that is under developement.  TEC relies on the centroid analysis from the NASA Ames pipeline to warn about possible centroid shifts during transit.  The DAVE project is potentially gearing up for more detailed pixel level analysis that I encourage those interested to follow the DAVE project.

### TEC Outputs:
TEC currently consists of 5 primary outputs that are of most interest to the community.  The columns and format of the files are described in doc/TESS-ExoClass_File_Guide.txt
1. Triage list - Based upon a limited set of metrics that are simple, efficient, and robustly deal with 'junky' candidates, TEC does a triage filter on the NASA Ames TCEs.  This is a hard filter in that TCEs that do not pass the TEC triage filter do not get a TEC summary report and do not have the more detailed metrics calculated.  From a classification standpoint, the TCEs that do not pass this filter should be considered instrumental and/or stellar variability in nature.  The TCEs that pass the triage are considered viable transiting planet candidates and/or detached stellar eclipsing binaries.  Contact binaries that are not detached (detached EBs have a flat stable normalized flux time series in between the eclipses) tend to not pass the triage filter stage.
2. Tier 1 list - List of most promising planet candidates that pass the triage filter and all the additional metrics indicate that the signals are consistent with a planet candidate on the target.  In addition, the Tier 1 list is ranked by a heuristic value that gives higher rank to planet candidates that have brighter hosts, smaller radii, lower insolation flux, and etc...(quantities most suitable for follow up).  In addition, the Tier 1 list provides a flag as to whether the planet candidate is previously known from the confirmed planet list hosted by the NASA Exoplanet Archive or a previously known TOI.  End users (especially observers making plans for a telescope) are still encouraged to perform due diligence and manually examine the summary report pdf in order to confirm that the planet candidate is viable.
3. Tier 2 list - List of planet candidates that pass the triage filter, but one or more metrics indicate an issue with the signal that raises concern the signal is not a planet candidate or the signal is not on the primary target of interest (evidence for a centroid offset).  The Tier 2 list is a similarly ranked list as Tier 1, but contains a flag string to indicate which metric has raised concern.  End users are encouraged to examine the summary report pdf in order to provide a judgement as to whether the metric of concern is significant and calculated in a valid manner.  Some of the metrics can have false negatives and be good clean planet candidates.
4. Tier 3 list - List of probable stellar eclipsing binaries.  In particular, this list contains signals with eclipses and a significant secondary that is not due to a planet.
5. Summary report - Users of the [Kepler TCERT Vetting Reports](https://exoplanetarchive.ipac.caltech.edu/docs/KSCI-19105-001.pdf) from the Kepler Science Office will feel right at home with the TEC summary reports.  The TEC summary report repeats the SPOC DV summary and includes the additional information calculated in TEC such as the very helpful [model shift test](https://github.com/JeffLCoughlin/Model-Shift), SPOC centroid information, and a prototype for an independent difference image calculation.

### Documentation
The doc folder has documentation that describes the columns for the TEC outputs, a semi-thorough step-by-step guide for running TEC on the NASA Ames SPOC TCEs along with hints for implementing TEC codes for your own pipeline results.

### Code
TEC code is in the code folder.  It is primarily python.  See doc/tess-exoclass_opsnotes_s04.txt for the guide for running the code for the NASA Ames SPOC TCEs with data products available at MAST.

 
