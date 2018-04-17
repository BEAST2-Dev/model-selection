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
import beast.core.parameter.RealParameter;
import beast.core.Input.Validate;
import beast.core.Logger;
import beast.core.Logger.LOGMODE;
import beast.core.util.CompoundDistribution;
import beast.core.util.Log;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeDistribution;
import beast.evolution.tree.TreeWithMetaDataLogger;
import beast.gss.distributions.GSSTreeDistribution;
import beast.gss.distributions.GSSTreeDistribution.BranchLengthDistribution;
import beast.gss.distributions.KernelDensityEstimatorDistribution;
import beast.gss.distributions.MultivariateKDEDistribution;
import beast.gss.distributions.NormalKDEDistribution;
import beast.math.distributions.MRCAPrior;
import beast.core.MCMC;
import beast.util.JSONProducer;
import beast.util.XMLParser;
import beast.util.XMLProducer;

@Description("MCMC-step class based on GSS sample")
public class MCMC2GSS extends Runnable {
	public Input<XMLFile> model1Input = new Input<>("xml",
			"file name of BEAST XML file containing the model for which to create a GSS XML file for",
			new XMLFile("examples/normalTest-1XXX.xml"), Validate.REQUIRED);
	public Input<Integer> traceBurninInput = new Input<>("traceBurnin",
			"percentage of the log file to disregard as burn-in (default 10)", 10);
	public Input<OutFile> outputInput = new Input<>("output", "where to save the file", new OutFile("beast.xml"));

	
	public Input<Double> alphaInput = new Input<>("alpha", "alpha parameter of Beta(alpha,1) distribution used to space out steps, default 0.3" +
			"If alpha <= 0, uniform intervals are used.", 0.3);
	public Input<Integer> stepsInput = new Input<>("nrOfSteps", "the number of steps to use, default 8", 8);
	public Input<String> rootDirInput = new Input<>("rootdir", "root directory for storing particle states and log files (default /tmp)", "/tmp");
	public Input<MCMC> mcmcInput = new Input<>("mcmc", "MCMC analysis used to specify model and operations in each of the particles", Validate.REQUIRED);
	public Input<Long> chainLengthInput = new Input<>("chainLength", "number of sample to run a chain for a single step", 100000l);
	public Input<Integer> burnInPercentageInput = new Input<>("burnInPercentage", "burn-In Percentage used for analysing log files", 50);
	public Input<Integer> preBurnInInput = new Input<>("preBurnin", "number of samples that are discarded for the first step, but not the others", 100000);
	public Input<String> m_sScriptInput = new Input<>("value", "script for launching a job. " +
			"$(dir) is replaced by the directory associated with the particle " +
			"$(java.class.path) is replaced by a java class path used to launch this application " +
			"$(java.library.path) is replaced by a java library path used to launch this application " +
			"$(seed) is replaced by a random number seed that differs with every launch " +
			"$(host) is replaced by a host from the list of hosts", Validate.REQUIRED);
	public Input<String> m_sHostsInput = new Input<>("hosts", "comma separated list of hosts. " +
			"If there are k hosts in the list, for particle i the term $(host) in the script will be replaced " +
			"by the (i modulo k) host in the list. " +
			"Note that whitespace is removed");
	public Input<Boolean> doNotRun = new Input<>("doNotRun", "Set up all files but do not run analysis if true. " +
			"This can be useful for setting up an analysis on a cluster", false);
	
	public Input<Boolean> deleteOldLogsInput = new Input<Boolean>("deleteOldLogs", "delete existing log files from root dir", false);

		
	java.util.Map<String, TraceLog> traceLogs = new LinkedHashMap<>();

	@Override
	public void initAndValidate() {
	}

	@Override
	public void run() throws Exception {
		XMLParser parser = new XMLParser();
		MCMC mcmc = (MCMC) parser.parseFile(model1Input.get());

		processTraceLogs(mcmc);

		GeneralisedSteppingStone gss = new GeneralisedSteppingStone();
		gss.alphaInput = new Input<Double>("alpha", "alpha parameter of Beta(alpha,1) distribution used to space out steps, default 0.3" +
				"If alpha <= 0, uniform intervals are used.", 0.3);
		gss.stepsInput.setValue(stepsInput.get(), gss);
		gss.rootDirInput.setValue(rootDirInput.get(), gss);
		gss.chainLengthInput.setValue(chainLengthInput.get(), gss);
		gss.burnInPercentageInput.setValue(burnInPercentageInput.get(), gss);
		gss.preBurnInInput.setValue(preBurnInInput.get(), gss);
		gss.m_sScriptInput.setValue(m_sScriptInput.get(), gss);
		gss.m_sHostsInput.setValue(m_sHostsInput.get(), gss);
		gss.doNotRun.setValue(doNotRun.get(), gss);
		gss.deleteOldLogsInput.setValue(deleteOldLogsInput.get(), gss);

		gss.mcmcInput.setValue(mcmc, gss);

		CompoundDistribution prior = getPrior(mcmc);
		Distribution samplingDistribution = getAltPrior(prior);
		gss.samplingDistributionInput.setValue(samplingDistribution, gss);

		String required = getRequiredAttribute();
		// save
		// String xml = new XMLProducer().toXML(mcmc.get(), );
		String spec = null;
		OutFile file = outputInput.get();
		if (file.getPath().toLowerCase().endsWith(".json")) {
			spec = toJSON(gss, required);
		} else {
			spec = toXML(gss, required);
		}
		FileWriter outfile = new FileWriter(file);
		outfile.write(spec);
		outfile.close();
		
		
		// run?
		if (!doNotRun.get()) {
			gss.initAndValidate();
			gss.run();
		}
	} // save

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

	private String getRequiredAttribute() throws IOException {
		BufferedReader fin = new BufferedReader(new FileReader(model1Input.get()));
		StringBuffer buf = new StringBuffer();
		String str = null;
		while (fin.ready()) {
			str = fin.readLine();
			buf.append(str);
			buf.append('\n');
		}
		fin.close();
		String xml = buf.toString();
		if (xml.matches("required=['\"]")) {
			int start = xml.indexOf("required=") + 10;
			int end = xml.indexOf("'", start);
			int end2 = xml.indexOf("\"", start);
			if (end > 0 && end2 > 0) {
				end = Math.min(end, end2);
			} else if (end < 0) {
				end = end2;
			}
			String required = xml.substring(start, end);
			return required;
		}
		return "";
	}

	private String toJSON(GeneralisedSteppingStone gss, String required) {
		Set<BEASTInterface> beastObjects = new HashSet<>();
		String json = new JSONProducer().toJSON(gss, beastObjects);

		json = json.replaceFirst("\\{", "{ required:\"" + required + "\", ");
		return json + "\n";
	}

	public String toXML(BEASTObject gss, String required) {
		Set<BEASTInterface> beastObjects = new HashSet<>();
		String xml = new XMLProducer().toXML(gss, beastObjects);

		xml = xml.replaceFirst("<beast ", "<beast required='" + required + "' ");
		return xml + "\n";
	}

	private CompoundDistribution getAltPrior(CompoundDistribution prior) {
		List<Distribution> altPrior = new ArrayList<>();
		for (Distribution d : prior.pDistributions.get()) {
			if (d instanceof TreeDistribution) {
				Distribution altTreeDist = getAltTreeDist((TreeDistribution) d);
				altPrior.add(altTreeDist);
			} else if (d instanceof MRCAPrior) {
				altPrior.add(d);
			} else if (d instanceof beast.math.distributions.Prior) {
				Distribution altPriorDist = getAltPriorDist((beast.math.distributions.Prior) d);
				altPrior.add(altPriorDist);
			} else {
				throw new IllegalArgumentException(
						"Don't know how to handle distributio " + d.getID() + " of type " + d.getClass().getName());
			}

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
				"useGammaForBranchLengths", BranchLengthDistribution.none);
		ccDistr.setID(d.getID()+ ".gss");
		return ccDistr;
	}

	private Distribution getAltPriorDist(beast.math.distributions.Prior d) {
		Function f = d.m_x.get();
		String id = ((BEASTInterface) f).getID();
		String shortid = id.contains(".") ? id.substring(0, id.lastIndexOf('.')) : id;
		String shortid2 = id.contains(":") ? id.substring(0, id.lastIndexOf(':')-1) + id.substring(id.lastIndexOf(':') + 1) : id;
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

		altDist.setID(d.getID()+ ".gss");
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

	public static void main(String[] args) throws Exception {
		new Application(new MCMC2GSS(), "MCMC2GSS", args);
	}
}
