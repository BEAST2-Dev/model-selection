package test.simstudy;

// Requirements: the file analysis.xml must be in the directory
// from where this script is run.

import beast.base.evolution.tree.Tree;
import beast.base.inference.Logger;
import beast.base.inference.parameter.RealParameter;
import beast.base.parser.NexusParser;
import beast.base.parser.XMLParserException;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.branchratemodel.StrictClockModel;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.substitutionmodel.Frequencies;
import beast.base.evolution.substitutionmodel.HKY;
import beastfx.app.seqgen.MergeDataWith;
import beastfx.app.seqgen.SequenceSimulator;

import java.io.File;
import java.io.IOException;
import java.util.List;

import beagle.BeagleFlag;

public class GenerateXMLs {
	// This script performs a simulation study for the
	// bModelTest model by sampling from the prior, simulate sequences and
	// running ananlyses to recover the model.

	static int N = 100;




	static void process(String dir, String template) throws IllegalArgumentException, IllegalAccessException, IOException, XMLParserException {

		if (!(new File(dir).exists())) {
			new File(dir).mkdirs();
		}
		System.err.print("Processing " + dir);
		
	for (int i = 0; i < N; i++) {

		// set up model to draw samples from
		// set up model to draw samples from
		Sequence A=new Sequence(); A.initByName("taxon","sequence01","value","?");
		Sequence B=new Sequence(); B.initByName("taxon","sequence02","value","?");
		Sequence C=new Sequence(); C.initByName("taxon","sequence03","value","?");
		Sequence D=new Sequence(); D.initByName("taxon","sequence04","value","?");
		Sequence E=new Sequence(); E.initByName("taxon","sequence05","value","?");
		Sequence F=new Sequence(); F.initByName("taxon","sequence06","value","?");
		Sequence G=new Sequence(); G.initByName("taxon","sequence07","value","?");
		Sequence I=new Sequence(); I.initByName("taxon","sequence08","value","?");
		Sequence H=new Sequence(); H.initByName("taxon","sequence09","value","?");
		Sequence J=new Sequence(); J.initByName("taxon","sequence10","value","?");
		
		Sequence A1=new Sequence(); A1.initByName("taxon","sequence11","value","?");
		Sequence B1=new Sequence(); B1.initByName("taxon","sequence12","value","?");
		Sequence C1=new Sequence(); C1.initByName("taxon","sequence13","value","?");
		Sequence D1=new Sequence(); D1.initByName("taxon","sequence14","value","?");
		Sequence E1=new Sequence(); E1.initByName("taxon","sequence15","value","?");
		Sequence F1=new Sequence(); F1.initByName("taxon","sequence16","value","?");
		Sequence G1=new Sequence(); G1.initByName("taxon","sequence17","value","?");
		Sequence I1=new Sequence(); I1.initByName("taxon","sequence18","value","?");
		Sequence H1=new Sequence(); H1.initByName("taxon","sequence19","value","?");
		Sequence J1=new Sequence(); J1.initByName("taxon","sequence20","value","?");


		Sequence A2=new Sequence(); A2.initByName("taxon","sequence21","value","?");
		Sequence B2=new Sequence(); B2.initByName("taxon","sequence22","value","?");
		Sequence C2=new Sequence(); C2.initByName("taxon","sequence23","value","?");
		Sequence D2=new Sequence(); D2.initByName("taxon","sequence24","value","?");
		Sequence E2=new Sequence(); E2.initByName("taxon","sequence25","value","?");
		Sequence F2=new Sequence(); F2.initByName("taxon","sequence26","value","?");
		Sequence G2=new Sequence(); G2.initByName("taxon","sequence27","value","?");
		Sequence I2=new Sequence(); I2.initByName("taxon","sequence28","value","?");
		Sequence H2=new Sequence(); H2.initByName("taxon","sequence29","value","?");
		Sequence J2=new Sequence(); J2.initByName("taxon","sequence30","value","?");

		Sequence A3=new Sequence(); A3.initByName("taxon","sequence31","value","?");
		Sequence B3=new Sequence(); B3.initByName("taxon","sequence32","value","?");
		Sequence C3=new Sequence(); C3.initByName("taxon","sequence33","value","?");
		Sequence D3=new Sequence(); D3.initByName("taxon","sequence34","value","?");
		Sequence E3=new Sequence(); E3.initByName("taxon","sequence35","value","?");
		Sequence F3=new Sequence(); F3.initByName("taxon","sequence36","value","?");
		Sequence G3=new Sequence(); G3.initByName("taxon","sequence37","value","?");
		Sequence I3=new Sequence(); I3.initByName("taxon","sequence38","value","?");
		Sequence H3=new Sequence(); H3.initByName("taxon","sequence39","value","?");
		Sequence J3=new Sequence(); J3.initByName("taxon","sequence40","value","?");

		Alignment data = new Alignment();
		data.initByName("sequence", A, "sequence", B, "sequence", C, "sequence", D, "sequence", E,
				"sequence", F, "sequence", G, "sequence", H, "sequence", I, "sequence", J,
				"sequence", A1, "sequence", B1, "sequence", C1, "sequence", D1, "sequence", E1,
				"sequence", F1, "sequence", G1, "sequence", H1, "sequence", I1, "sequence", J1,
				"sequence", A2, "sequence", B2, "sequence", C2, "sequence", D2, "sequence", E2,
				"sequence", F2, "sequence", G2, "sequence", H2, "sequence", I2, "sequence", J2,				
				"sequence", A3, "sequence", B3, "sequence", C3, "sequence", D3, "sequence", E3,
				"sequence", F3, "sequence", G3, "sequence", H3, "sequence", I3, "sequence", J3
		);
		Tree tree = trees.get(i);

		RealParameter freqs=new RealParameter("0.25 0.25 0.25 0.25") ;
		Frequencies f = new Frequencies();
		f.initByName("frequencies",freqs);
		
		HKY hky = new beast.base.evolution.substitutionmodel.HKY();
		hky.initByName("frequencies", f, 
			"kappa", "1.0"
		);
		StrictClockModel clockmodel = new beast.base.evolution.branchratemodel.StrictClockModel();
		clockmodel.initByName("clock.rate","1.0");


		// change gammaCategoryCount=1 for generating without gamma rate categories
		SiteModel sitemodel = new SiteModel();
		sitemodel.initByName("gammaCategoryCount", 1, "substModel", hky, "shape", "1.0", "proportionInvariant", "0.0");
		MergeDataWith mergewith = new MergeDataWith();
		mergewith.initByName("template", template, "output", dir + "/analysis-out" + i + ".xml");
		SequenceSimulator sim = new SequenceSimulator();
		sim.initByName("data", data, "tree", tree, "sequencelength", 2000, "outputFileName", 
				"gammaShapeSequence.xml", "siteModel", sitemodel, "branchRateModel", clockmodel, 
				"merge", mergewith);
		// produce gammaShapeSequence.xml and merge with analysis.xml to get analysis-out.xml
		sim.run();
		System.err.print('.');
	  }
	System.err.println();
	}


	static List<Tree> trees;

	static String wdir = "/Users/remco/workspace/model-selection/src/test/simstudy/data";

			public static void main(String[] args) throws IOException, IllegalArgumentException, IllegalAccessException, XMLParserException {
		
		Logger.FILE_MODE = beast.base.inference.Logger.LogFileMode.overwrite;

		// set up flags for BEAGLE -- YMMV
		long beagleFlags = BeagleFlag.VECTOR_SSE.getMask() | BeagleFlag.PROCESSOR_CPU.getMask();
		System.setProperty("beagle.preferred.flags", Long.toString(beagleFlags));
		

		for (String in : new String[] {"cc", "ec0.01", "ec0.025", "ec0.05", "ec0.1"}) {
			NexusParser parser = new beast.base.parser.NexusParser();
			File fin = new File(wdir + "/input/" + in + ".trees");
			parser.parseFile(fin);
			trees = parser.trees;
			trees.remove(0); // burn-in
			if (trees.size() != N) {
				throw new IllegalArgumentException("number of trees in tree set != " + N + " but " + trees.size());
			}
	
			process(wdir + "/" + in + "/cc/", wdir + "/template/cc.xml");
			process(wdir + "/" + in + "/ec/", wdir + "/template/ec.xml");
		
		}

		System.err.println("Done");

	}

}
