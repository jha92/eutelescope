## Example for a EUTelescope analysis using the GBL processors for evaluating a scatter (SUT) with material budget imaging (MBI).


Detailed information on how to use EUTELESCOPE software is available at

http://eutelescope.desy.de



Source the build_env.sh from the EUTel main directory to set the correct environment:

`source $EUTELESCOPE/build_env.sh`



The analysis we want to perform is controlled by a config file (config.cfg), a csv-table (runlist.csv) and steering file templates (*.xml), all located in the working directory.

You can set the variables for your working directory here as well as input parameters for the different processors.

One dataset is stored on AFS to try out this example (run number 110).



You should be able to run the analysis now. For this, execute the following commands one after another and check for each of the processors to successfully finish.

```
jobsub -c config.cfg -csv runlist.csv -g converter 110
jobsub -c config.cfg -csv runlist.csv -g clustering 110
jobsub -c config.cfg -csv runlist.csv -g hitmaker 110
jobsub -c config.cfg -csv runlist.csv -g aligngbl 110
jobsub -c config.cfg -csv runlist.csv -g aligngbl2 110
jobsub -c config.cfg -csv runlist.csv -g aligngbl3 110
jobsub -c config.cfg -csv runlist.csv -g gblkinkestimator 110
```


At every step a ROOT file is created, containing a set of histograms, which you can find at output/histograms/ after completion.

For the last reconstruction step (```gblkinkestimator```), one has to manually include ```SUT.xml``` (between upstream & downstream telescope planes) into the aligned gear-file 
(called ```*_pre_gbl1_gbl2_gbl3.xml```) and name it as ```*_aligned.xml```.
