package beast.gss;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Set;

import beast.app.util.Application;
import beast.app.util.LogFile;
import beast.app.util.OutFile;
import beast.app.util.TreeFile;
import beast.app.util.XMLFile;
import beast.core.BEASTInterface;
import beast.core.BEASTObject;
import beast.core.Description;
import beast.core.Distribution;
import beast.core.Function;
import beast.core.Input;
import beast.core.Runnable;
import beast.core.State;
import beast.core.StateNode;
import beast.core.StateNodeInitialiser;
import beast.core.parameter.RealParameter;
import beast.core.Input.Validate;
import beast.core.Logger;
import beast.core.Logger.LOGMODE;
import beast.core.util.CompoundDistribution;
import beast.core.util.Log;
import beast.evolution.tree.RandomTree;
import beast.evolution.tree.SimpleRandomTree;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeDistribution;
import beast.evolution.tree.TreeInterface;
import beast.evolution.tree.TreeWithMetaDataLogger;
import beast.gss.distributions.GSSTreeDistribution;
import beast.gss.distributions.GSSTreeDistribution.BranchLengthDistribution;
import beast.gss.distributions.KernelDensityEstimatorDistribution;
import beast.gss.distributions.MultivariateKDEDistribution;
import beast.gss.distributions.NormalKDEDistribution;
import beast.math.distributions.MRCAPrior;
import beast.core.MCMC;
import beast.util.JSONProducer;
import beast.util.TreeParser;
import beast.util.XMLParser;
import beast.util.XMLProducer;

@Description("Convert MCMC analysis to GSS analysis and optionally run the analysis")
abstract public class MCMC2IS extends Runnable {
	public Input<XMLFile> model1Input = new Input<>("xml",
			"file name of BEAST XML file containing the model for which to create a GSS XML file for",
			new XMLFile("examples/normalTest-1XXX.xml"), Validate.REQUIRED);
	public Input<Integer> traceBurninInput = new Input<>("traceBurnin",
			"percentage of the log file to disregard as burn-in (default 10)", 10);
	public Input<OutFile> outputInput = new Input<>("output", "where to save the file", new OutFile("beast.xml"));


	public Input<MCMC> mcmcInput = new Input<>("mcmc", "MCMC analysis used to specify model and operations in each of the particles", Validate.REQUIRED);
	public Input<Long> chainLengthInput = new Input<>("chainLength", "number of sample to run a chain for a single step", 100000l);
	public Input<Integer> preBurnInInput = new Input<>("preBurnin", "number of samples that are discarded for the first step, but not the others", 100000);
	public Input<Boolean> doNotRun = new Input<>("doNotRun", "Set up all files but do not run analysis if true. " +
			"This can be useful for setting up an analysis on a cluster", false);
	

		
	java.util.Map<String, TraceLog> traceLogs = new LinkedHashMap<>();

	@Override
	public void initAndValidate() {
	}

	@Override
	public void run() throws Exception {
		XMLParser parser = new XMLParser();
		MCMC mcmc = (MCMC) parser.parseFile(model1Input.get());

		processTraceLogs(mcmc);

		Runnable gss = newInstance(mcmc);
//		gss.alphaInput = new Input<Double>("alpha", "alpha parameter of Beta(alpha,1) distribution used to space out steps, default 0.3" +
//				"If alpha <= 0, uniform intervals are used.", 0.3);
		gss.setInputValue("chainLength",chainLengthInput.get());
		gss.setInputValue("preBurnin",preBurnInInput.get());
		try {
			gss.setInputValue("doNotRun", doNotRun.get());
		} catch (IllegalArgumentException e) {
			// ignore
		}


		CompoundDistribution prior = getPrior(mcmc);
		Distribution samplingDistribution = getAltPrior(prior, mcmc.startStateInput.get());
		gss.setInputValue("samplingDistribution",samplingDistribution);
		setUpInitialisers(mcmc.initialisersInput.get());

		// save
		// String xml = new XMLProducer().toXML(mcmc.get(), );
		String spec = null;
		OutFile file = outputInput.get();
		if (file.getPath().toLowerCase().endsWith(".json")) {
			spec = toJSON(gss);
		} else {
			spec = toXML(gss);
		}
		FileWriter outfile = new FileWriter(file);
		outfile.write(spec);
		outfile.close();
		
		
		// run?
		if (!doNotRun.get()) {
			gss.initAndValidate();
			gss.run();
		}
		Log.warning("Done");
	} // save

	abstract protected Runnable newInstance(MCMC mcmc);		

	private void setUpInitialisers(List<StateNodeInitialiser> initialiser) {
		for (int i = 0; i < initialiser.size(); i++) {
			StateNodeInitialiser init = initialiser.get(i);
			if (init instanceof RandomTree || init instanceof SimpleRandomTree) {				
				// grab last tree from log file
				TreeParser newInit = new TreeParser();
				Tree tree = (Tree) ((BEASTInterface)init).getInput("initial").get();
				String newick = null;
				for (BEASTInterface o : tree.getOutputs()) {
					if (o instanceof GSSTreeDistribution) {
						newick = ((GSSTreeDistribution)o).getLastTre().getRoot().toNewick();
					}
				}
				newick = newick.replaceAll("\\[[^\\[\\]]+\\]", "");
				if (newick != null) {
					newInit.initByName("newick", newick, 
							"IsLabelledNewick", true, 
							"initial", tree,
							"taxa",  ((BEASTInterface)init).getInput("taxa").get());
					newInit.m_taxonset.set(null);
					newInit.setID(((BEASTInterface)init).getID());
					initialiser.set(i, newInit);
				}
			}
		}
	}

	/** create map with log labels to trace logs **/
	private void processTraceLogs(MCMC mcmc) {
		for (Logger l : mcmc.loggersInput.get()) {
			if (!l.modeInput.get().equals(LOGMODE.tree) && l.fileNameInput.get() != null
					&& !l.fileNameInput.get().equals("")) {
				TraceLog tracelog = new TraceLog();
				tracelog.burnInPercentageInput.setValue(traceBurninInput.get(), tracelog);
				String fileName = l.fileNameInput.get();
				LogFile file = new LogFile(fileName);
				if (!file.exists()) {
					fileName = model1Input.get().getParent() + File.separatorChar + file.getName();
					file = new LogFile(fileName);
					if (!file.exists()) {
						throw new IllegalArgumentException("Cannot find log file " + fileName);
					}
				}
				tracelog.traceFileInput.setValue(file, tracelog);
				tracelog.initAndValidate();

				List<String> labels = tracelog.getLabels();
				for (String label : labels) {
					traceLogs.put(label, tracelog);
				}

			}
		}

	}

	private String toJSON(BEASTObject gss) {
		Set<BEASTInterface> beastObjects = new HashSet<>();
		String json = new JSONProducer().toJSON(gss, beastObjects);

		return json + "\n";
	}

	public String toXML(BEASTObject gss) {
		Set<BEASTInterface> beastObjects = new HashSet<>();
		String xml = new XMLProducer().toXML(gss, beastObjects);
		return xml + "\n";
	}

	private CompoundDistribution getAltPrior(CompoundDistribution prior, State state) {
		List<StateNode> stateNodes = new ArrayList<>();
		stateNodes.addAll(state.stateNodeInput.get());
		
		List<Distribution> altPrior = new ArrayList<>();
		for (Distribution d : prior.pDistributions.get()) {
			if (d instanceof TreeDistribution) {
				Distribution altTreeDist = getAltTreeDist((TreeDistribution) d);
				altPrior.add(altTreeDist);
				
				TreeInterface tree = ((TreeDistribution) d).treeInput.get();
				if (tree == null) {
					tree = ((TreeDistribution) d).treeIntervalsInput.get().treeInput.get();
				}
				stateNodes.remove(tree);
			} else if (d instanceof MRCAPrior) {
				altPrior.add(d);
			} else if (d instanceof beast.math.distributions.Prior) {
				beast.math.distributions.Prior p = (beast.math.distributions.Prior) d;
				Distribution altPriorDist = getAltPriorDist(p.m_x.get(), p.getID());
				altPrior.add(altPriorDist);
				Object o = ((beast.math.distributions.Prior) d).m_x.get();
				stateNodes.remove(o);
			} else {
				throw new IllegalArgumentException(
						"Don't know how to handle distributio " + d.getID() + " of type " + d.getClass().getName());
			}

		}
		for (int i = stateNodes.size() - 1; i >= 0; i--) {
			StateNode s = stateNodes.get(i);
			if (s instanceof RealParameter) {
				Distribution altPriorDist = getAltPriorDist(s, s.getID() + "Prior");
				altPrior.add(altPriorDist);
				stateNodes.remove(s);
			}
		}
		for (StateNode s : stateNodes) {
			Log.warning("No working distribution for " + s.getID());
		}
		
		
		
		
		CompoundDistribution cd = new CompoundDistribution();
		cd.initByName("distribution", altPrior);
		cd.setID("GSSPrior");
		return cd;
	}

	private Distribution getAltTreeDist(TreeDistribution d) {
		// attempt to find associated tree log
		String fileName = null;
		Logger l = null;
		Tree tree = (Tree) d.getInput("tree").get();
		for (BEASTInterface o : tree.getOutputs()) {
			if (o instanceof Logger) {
				l = (Logger) o;
				if (l.modeInput.get().equals(LOGMODE.tree)) {
					fileName = l.fileNameInput.get();
					break;
				}
			} else if (o instanceof TreeWithMetaDataLogger) {
				for (BEASTInterface o2 : o.getOutputs()) {
					if (o2 instanceof Logger) {
						l = (Logger) o2;
						if (l.modeInput.get().equals(LOGMODE.tree)) {
							fileName = l.fileNameInput.get();
							break;
						}
					}
				}
			}
		}
		if (fileName.contains("$(tree)")) {
			String treeName = "tree";
			if (l != null) {
				for (final BEASTInterface logger : l.loggersInput.get()) {
					if (logger instanceof BEASTObject) {
						final String id = ((BEASTObject) logger).getID();
						if (id.indexOf(".t:") > 0) {
							treeName = id.substring(id.indexOf(".t:") + 3);
						}
					}
				}
			}
			fileName = fileName.replace("$(tree)", treeName);
		}
		TreeFile file = new TreeFile(fileName);
		if (!file.exists()) {
			fileName = model1Input.get().getParent() + File.separatorChar + file.getName();
			file = new TreeFile(fileName);
			if (!file.exists()) {
				throw new IllegalArgumentException("Cannot find tree file " + fileName);
			}
		}

		// create GSS tree distribution
		GSSTreeDistribution ccDistr = new GSSTreeDistribution();
		ccDistr.initByName("treefile", file, 
				"tree", d.treeInput.get(),
				"burnin", traceBurninInput.get(), 
				"useGammaForBranchLengths", BranchLengthDistribution.useGamma);
		ccDistr.setID(d.getID()+ ".gss");
		return ccDistr;
	}

	private Distribution getAltPriorDist(Function f, String priorID) { //beast.math.distributions.Prior d) {
		//Function f = 
		String id = ((BEASTInterface) f).getID();
		String shortid = id.contains(".") ? id.substring(0, id.lastIndexOf('.')) : id;
		String shortid2 = id.contains(":") ? id.substring(0, id.lastIndexOf(':')-1) + id.substring(id.lastIndexOf(':') + 1) : id;
		String shortid3 = id.contains(".") ? id.substring(0, id.lastIndexOf('.') + 1) : id;
		TraceLog tracelog = traceLogs.get(shortid);
		String label = shortid;
		if (tracelog == null) {
			tracelog = traceLogs.get(id);
			label = id;
		}
		if (tracelog == null) {
			label = id;
			id += "1";
			tracelog = traceLogs.get(id);
		}
		if (tracelog == null) {
			label = shortid2;
			tracelog = traceLogs.get(label);
		}
		if (tracelog == null) {
			label = shortid2;
			id = label + "1";
			tracelog = traceLogs.get(id);
		}
		if (tracelog == null) {
			label = shortid3;
			tracelog = traceLogs.get(label);
		}
		if (tracelog == null) {
			label = shortid3;
			id = label + "1";
			tracelog = traceLogs.get(id);
		}
		if (tracelog == null) {
			Log.warning.println("Did not find entry " + id + " in any log.");
		}

		Distribution altDist = null;

		int dim = f.getDimension();
		if (dim == 1) {
			altDist = new NormalKDEDistribution(tracelog, label, f);
			if (f instanceof RealParameter) {
				RealParameter p = (RealParameter) f;
				Double mean = tracelog.getMean(label);
				Log.warning("Set value " + p.getID() + " to " + mean);
				p.valuesInput.get().set(0, mean);
			}
		} else {
			KernelDensityEstimatorDistribution[] multivariateKDE = new KernelDensityEstimatorDistribution[dim];
			for (int i = 0; i < dim; i++) {
				multivariateKDE[i] = new NormalKDEDistribution(tracelog, label + (i+1), null);
				if (f instanceof RealParameter) {
					RealParameter p = (RealParameter) f;
					Double mean = tracelog.getMean(label + (i+1));
					Log.warning("Set value " + p.getID() + "[" + (i+1) + "] to " + mean);
					if (i < p.valuesInput.get().size()) {
						p.valuesInput.get().set(i, mean);
					} else {
						p.valuesInput.get().add(mean);
					}
				}
			}
			altDist = new MultivariateKDEDistribution(multivariateKDE, f);
		}
		// EmpiricalDist(d.m_x.get(), trace);

		altDist.setID(priorID + ".gss");
		return altDist;
	}

	private CompoundDistribution getPrior(MCMC mcmc) {
		if (!(mcmc.posteriorInput.get() instanceof CompoundDistribution)) {
			throw new IllegalArgumentException(
					"The XML file does not seem to contain an MCMC analysis with posterior CompoundDistribution");
		}
		CompoundDistribution posterior = (CompoundDistribution) mcmc.posteriorInput.get();
		for (Distribution d : posterior.pDistributions.get()) {
			if (d.getID().equals("prior")) {
				if (!(d instanceof CompoundDistribution)) {
					throw new IllegalArgumentException(
							"The XML file does not seem to contain an MCMC analysis with prior CompoundDistribution");
				}
				return (CompoundDistribution) d;
			}
		}
		throw new IllegalArgumentException(
				"The XML file does not seem to contain an MCMC analysis with a 'prior' in the posterior");
	} // getPrior

}
