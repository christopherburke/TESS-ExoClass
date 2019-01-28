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

### What does TEC enable:
- Observing the best candidates: TEC provides a filtered and ranked list of the most viable detections to get you to the 'gem' planet candidates quickly.  In practice, TEC has demonstrated that the number of viable detections can be reduced by a factor of six from the original detections.
- Machine learning attributes and training sets: All machine learning algorithms begin with a database of attributes/metrics.  You can use TEC to enhance the current set of attributes that you have developed or you can use the currently crude classifications from TEC to enhance your training set.
- Statistical Validation: Do you want to statistically validate your lowest signal to noise ratio (SNR) detections?  Then you need to quantify the reliability of your detections and classification.  [Mullally et al. (2018)](http://adsabs.harvard.edu/abs/2018AJ....155..210M) and [Burke et al. (2019)](http://adsabs.harvard.edu/abs/2019arXiv190100506B) recently have shown that important Kepler statistical validations (such as Kepler-452b and Kepler-186f) are not as statistically significant as previously thought without taking into account the possibility that they could be false alarm/systematic detections.  TEC enables automating the classification thus one can provide TEC examples of systematics in order to demonstrate whether reliability is an issue for the lowest SNR detections.
- Planet Occurrence Rates: Are you one of the brave souls that wants to calculate unbiased planet occurrence rates based on TESS data?  You need automated classifications in order to quantify how this step influences the recovered planets.
- Your Pipeline Here: I have my own planet search pipeline, and my own large pile of detections.  Help!  Without too much work providing a flux-time series and detection ephemeris you can run TEC for your own pipeline.
