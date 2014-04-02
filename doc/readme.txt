More info can be found on the wiki
http://www.beast2.org/wiki/index.php/Path_Sampling

To run a path sampling/stepping stone analysis from a GUI, run from the command line:

java -cp $full_path_to_mode_selection_package/lib/MODEL_SELECTION.addon.jar:$full_path_to_beast.jar/beast.jar beast.app.tools.PathSampler

or, if you want to run paired path sampling analysis

java -cp $full_path_to_mode_selection_package/lib/MODEL_SELECTION.addon.jar:$full_path_to_beast.jar/beast.jar beast.app.tools.PairedPathSampler

Note: you need to give the full paths to the jar files, since these are used to set up processes for the steps in the analysis. 
Hint: On Mac or Linux, you can go to the package directory and use `pwd`/lib/MODEL_SELECTION.addon.jar to shorten the first bit.

