package modelselection.app.tools;

import beastfx.app.tools.Application;

//command line interface to PathSampleAnalyser
public class PathSampleAnalyser {
		
	public static void main(final String[] args) throws Exception {
		new Application(new modelselection.inference.PathSampleAnalyser(), "Path Smple Analyser", args);
	}
}
