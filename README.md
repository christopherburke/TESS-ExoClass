# TESS-ExoClass (TEC)
TESS exoplanet detection filter, crossmatch, classifier, and ranker.

### Problem:
Your transiting planet detection pipeline results in a huge pile of detections.  Now What?  Do I really have to look at all these manually!?

### Answer:
No!  Algorithmic filtering and classifying is available with TESS-ExoClass (TEC).

### Authors:
Christopher J Burke (MIT), Jeffrey Coughlin (SETI/NASA Ames), Susan Mullally (STScI), Fergal Mullally (SETI), Jessie Christiansen (NExSci)

### Description:
TEC gets you from your pile of signal detections to your best transiting planet candidate science faster.  TEC builds off the highly successful Kepler Robovetter classification process ( [Coughlin et al. 2016](http://adsabs.harvard.edu/abs/2016ApJS..224...12C), [Thompson et al. 2018](http://adsabs.harvard.edu/abs/2018ApJS..235...38T) ) as well as the [DAVE](http://keplertcert.seti.org/DAVE/) K2 vetting tools ( [Kostov et al. 2019](http://adsabs.harvard.edu/abs/2019arXiv190107459K) ).  TEC is actually neither of those more mature efforts (yet).  TEC is actually the foundation underlying those classifiers.  In other words TEC is the database of attributes/metrics that generically feed any kind of machine learning or classification process.  For instance, the Kepler Robovetter examined a database of ~50 metrics in order for it to carry out its decision tree based classifications.  TEC is currently focused on developing the database of attributes/metrics that can provide more sophisticated classifications using machine learning techniques.  TEC currently uses a crude decision tree like classifier that is tuned manually in order to match the efficacy of the current manual vetting effort to build the MIT TOI alerts.
