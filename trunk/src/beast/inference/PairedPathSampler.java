package beast.inference;

import java.text.DecimalFormat;
import java.util.List;
import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import beast.app.BeastMCMC;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Logger;
import beast.core.MCMC;
import beast.core.Operator;
import beast.core.BEASTObject;
import beast.core.StateNode;
import beast.util.XMLParser;
import beast.util.XMLProducer;

@Description("Calculate marginal likelihood through path/stepping stone sampling for comparing two models. "
		+ "Perform multiple steps and calculate estimate."
		+ "Uses multiple threads if specified as command line option to BEAST. "
		+ "This uses the operator schedule of the first model.")
public class PairedPathSampler extends PathSampler {
	public Input<String> model1Input = new Input<String>(
			"model1",
			"file name of BEAST XML file containing the first model that needs to be compared",
			Validate.REQUIRED);
	public Input<String> model2Input = new Input<String>("model2",
			"file name of second model that needs to be compared",
			Validate.REQUIRED);
	
	public enum Scheme {
		sigmoid, uniform
	}
	public Input<Scheme> schemInput = new Input<Scheme>("scheme" ,"sampling scheme, one of " + Scheme.values(), Scheme.sigmoid);
	

	MCMC model1;
	MCMC model2;
	Set<String> mergedSet;
	
	public PairedPathSampler() {
		mcmcInput.setRule(Validate.OPTIONAL);
		m_sScriptInput.setRule(Validate.OPTIONAL);
	}

	@Override
	public void initAndValidate() throws Exception {

		XMLParser parser1 = new XMLParser();
		Object o = parser1.parseFile(new File(model1Input.get()));
		if (!(o instanceof MCMC)) {
			throw new Exception("The model in " + model1Input.get()
					+ " does not appear to be an MCMC analysis.");
		}
		model1 = (MCMC) o;

		XMLParser parser2 = new XMLParser();
		o = parser2.parseFile(new File(model2Input.get()));
		if (!(o instanceof MCMC)) {
			throw new Exception("The model in " + model2Input.get()
					+ " does not appear to be an MCMC analysis.");
		}
		model2 = (MCMC) o;
		
		mergedSet = new HashSet<String>();

		mergeModel2IntoModel1();

		// merge operators
		for (Operator operator : model2.operatorsInput.get()) {
			if (!mergedSet.contains(operator.getID())) {
				model1.operatorsInput.setValue(operator, model1);
			}
		}

		// merge states
		for (StateNode stateNode : model2.startStateInput.get().stateNodeInput
				.get()) {
			if (!mergedSet.contains(stateNode.getID())) {
				model1.startStateInput.get().stateNodeInput
						.setValue(stateNode, model1);
			}
		}

		generateStepFiles();
	}

	private void generateStepFiles() throws Exception {
		// grab info from inputs
		m_sScript = m_sScriptInput.get();
		if (m_sHostsInput.get() != null) {
			m_sHosts = m_sHostsInput.get().split(",");
			// remove whitespace
			for (int i = 0; i < m_sHosts.length; i++) {
				m_sHosts[i] = m_sHosts[i].replaceAll("\\s", "");
			}
		}

		m_nSteps = stepsInput.get();
		if (m_nSteps <= 1) {
			throw new Exception("number of steps should be at least 2");
		}
		burnInPercentage = burnInPercentageInput.get();
		if (burnInPercentage < 0 || burnInPercentage >= 100) {
			throw new Exception("burnInPercentage should be between 0 and 100");
		}
		int preBurnIn = preBurnInInput.get();

		// root directory sanity checks
		File rootDir = new File(rootDirInput.get());
		if (!rootDir.exists()) {
			throw new Exception("Directory " + rootDirInput.get()
					+ " does not exist.");
		}
		if (!rootDir.isDirectory()) {
			throw new Exception(rootDirInput.get() + " is not a directory.");
		}

		// initialise MCMC
		PairedPathSamplingStep step = new PairedPathSamplingStep();
		for (Input<?> input : model1.listInputs()) {
			try {
				if (input.get() instanceof List) {
					for (Object o : (List<?>) input.get()) {
						step.setInputValue(input.getName(), o);
					}
				} else {
					step.setInputValue(input.getName(), input.get());
				}
			} catch (Exception e) {
				// TODO: handle exception
			}
		}
		step.posterior2Input.setValue(model2.posteriorInput.get(), step);

		int chainLength = chainLengthInput.get();
		// set up chain length for a single step
		step.burnInInput.setValue(0, step);
		step.chainLengthInput.setValue(chainLength, step);

		// add posterior logger
		Logger logger = new Logger();
		logger.initByName("fileName", LIKELIHOOD_LOG_FILE, 
				"log", step.posteriorInput.get(),
				"logEvery", chainLength / 1000);
		step.loggersInput.setValue(logger, step);

		// set up directories with beast.xml files in each of them
		String sFormat = "";
		for (int i = m_nSteps; i > 0; i /= 10) {
			sFormat += "#";
		}
		formatter = new DecimalFormat(sFormat);

		XMLProducer producer = new XMLProducer();

		PrintStream[] cmdFiles = new PrintStream[BeastMCMC.m_nThreads];
		for (int i = 0; i < BeastMCMC.m_nThreads; i++) {
			FileOutputStream outStream = (beast.app.util.Utils.isWindows() ? new FileOutputStream(
					rootDirInput.get() + "/run" + i + ".bat")
					: new FileOutputStream(rootDirInput.get() + "/run" + i
							+ ".sh"));
			cmdFiles[i] = new PrintStream(outStream);
		}

		for (int i = 0; i < m_nSteps; i++) {
			if (i < BeastMCMC.m_nThreads) {
				step.burnInInput.setValue(preBurnIn, step);
			} else {
				step.burnInInput.setValue(0, step);
			}
			// create XML for a single step
			double beta = nextBeta(i, m_nSteps, alphaInput.get());
//					betaDistribution != null ? betaDistribution
//					.inverseCumulativeProbability((i + 0.0) / (m_nSteps - 1))
//					: (i + 0.0) / (m_nSteps - 1);
			step.setInputValue("beta", beta);
			String sXML = producer.toXML(step);
			File stepDir = new File(getStepDir(i));
			if (!stepDir.exists() && !stepDir.mkdir()) {
				throw new Exception("Failed to make directory "
						+ stepDir.getName());
			}
			stepDir.setWritable(true, false);
			FileOutputStream xmlFile = new FileOutputStream(
					stepDir.getAbsoluteFile() + "/beast.xml");
			PrintStream out = new PrintStream(xmlFile);
			out.print(sXML);
			out.close();

			String cmd = getCommand(stepDir.getAbsolutePath(), i);
			FileOutputStream cmdFile = (beast.app.util.Utils.isWindows() ? new FileOutputStream(
					stepDir.getAbsoluteFile() + "/run.bat")
					: new FileOutputStream(stepDir.getAbsoluteFile()
							+ "/run.sh"));
			PrintStream out2 = new PrintStream(cmdFile);
			out2.print(cmd);
			out2.close();

			cmdFile = (beast.app.util.Utils.isWindows() ? new FileOutputStream(
					stepDir.getAbsoluteFile() + "/resume.bat")
					: new FileOutputStream(stepDir.getAbsoluteFile()
							+ "/resume.sh"));
			cmd = cmd.replace("-overwrite", "-resume");
			out2 = new PrintStream(cmdFile);
			out2.print(cmd);
			out2.close();
			// TODO: probably more efficient to group cmdFiles in block of
			// #steps/#threads
			// instead of skipping #threads steps every time.
			if (i >= BeastMCMC.m_nThreads) {
				String copyCmd = (beast.app.util.Utils.isWindows() ? "copy "
						+ getStepDir(i - BeastMCMC.m_nThreads)
						+ "\\beast.xml.state " + getStepDir(i) : "cp "
						+ getStepDir(i - BeastMCMC.m_nThreads)
						+ "/beast.xml.state " + getStepDir(i));
				cmdFiles[i % BeastMCMC.m_nThreads].print(copyCmd);
			}
			cmdFiles[i % BeastMCMC.m_nThreads].print(cmd);
			File script = new File(stepDir.getAbsoluteFile()
					+ (beast.app.util.Utils.isWindows() ? "/run.bat"
							: "/run.sh"));
			script.setExecutable(true);
		}
		for (int k = 0; k < BeastMCMC.m_nThreads; k++) {
			cmdFiles[k].close();
		}
	} // initAndValidate


 	/** sigmoid shaped steps **/
	double nextBeta(int step, int pathSteps, Double alpha) {
        if (step == 0) {
            return 1.0;
        } else if (step == pathSteps) {
            return 0.0;
        } else if (step > pathSteps) {
            return -1.0;
        } else {
            double xvalue = ((pathSteps - step)/((double)pathSteps)) - 0.5;
            return Math.exp(alpha*xvalue)/(Math.exp(alpha*xvalue) + Math.exp(-alpha*xvalue));
        }
    }
	

	/**
	 * replace all objects in model2 with those in model1 if they have the same
	 * functionality.
	 */
	private void mergeModel2IntoModel1() throws Exception {
		// collect objects from model 1 and model 2
		Map<String, BEASTObject> objects1 = new HashMap<String, BEASTObject>();
		collectObjects((BEASTObject) model1, objects1);
		Map<String, BEASTObject> objects2 = new HashMap<String, BEASTObject>();
		collectObjects((BEASTObject) model2, objects2);

		// merge those objects in model 2 that have are of the same class
		// and have the same inputs as in model 1
		boolean progress = true;
		while (progress) {
			// iterate over graph multiple times till no more merging takes
			// place
			progress = false;
			Set<String> objectSet = new HashSet<String>(); 
			objectSet.addAll(objects2.keySet());
			for (String id2 : objectSet) {
				BEASTObject plugin1 = objects1.get(id2);
				BEASTObject plugin2 = objects2.get(id2);
				if (plugin1 != null) {
					if (haveCommonInputs(plugin1, plugin2)) {
						mergePlugins(plugin1, plugin2);
						objects2.remove(id2);
						progress = true;
					}
				}
			}
		}

		// ensure IDs are unique
		for (String id2 : objects2.keySet()) {
			if (objects1.keySet().contains(id2)) {
				int i = 2;
				while (objects1.keySet().contains(id2 + i)) {
					i++;
				}
				BEASTObject plugin = objects2.get(id2);
				plugin.setID(id2 + i);
			}
		}
	}

	/** replace plugin2 of model2 by plugin1 **/
	private void mergePlugins(BEASTObject plugin1, BEASTObject plugin2) throws Exception {
		
		System.err.println("Merging " + plugin1.getID());
		mergedSet.add(plugin1.getID());
		
		Set<BEASTObject> outputSet = new HashSet<BEASTObject>();
		outputSet.addAll(plugin2.outputs);
		for (BEASTObject output : outputSet) {
			boolean found = false;
			for (Input<?> input : output.listInputs()) {
				if (input.get() != null) {
					if (input.get() instanceof List) {
						List list = (List<?>) input.get();
						if (list.contains(plugin2)) {
							list.remove(plugin2);
							list.add(plugin1);
							found = true;
						}
					} else if (input.get().equals(plugin2)) {
						input.setValue(plugin1, output);
						found = true;
					}
				}
			}
			if (!found) {
				throw new Exception(
						"Programming error: could not find plugin2 "
								+ plugin2.getID() + " in output "
								+ output.getID());
			}
		}
	}

	/** check whether plugin1 and plugin2 share inputs **/
	private boolean haveCommonInputs(BEASTObject plugin1, BEASTObject plugin2)
			throws IllegalArgumentException, IllegalAccessException, Exception {
		if (!plugin1.getClass().equals(plugin2.getClass())) {
			return false;
		}
		
		for (Input<?> input1 : plugin1.listInputs()) {
			Input<?> input2 = plugin2.getInput(input1.getName());
			if (input1.get() == null || input2.get() == null) {
				if (input1.get() != null || input2.get() != null) {
					// one input is null, the other is not
					return false;
				}
				// both inputs are null, so we are still fine
			} else if (input1.get() instanceof List) {
				List<?> list1 = (List<?>) input1.get();
				List<?> list2 = (List<?>) input2.get();
				if (list1.size() != list2.size()) {
					return false;
				}
				for (Object o : list1) {
					if (!list2.contains(o)) {
						return false;
					}
				}
			} else if (input1.get() instanceof String) {
				String str1 = (String) input1.get();
				String str2 = (String) input2.get();
				if (!str1.trim().equals(str2.trim())) {
					return false;
				}
			} else // it is a primitive or plugin
			if (!input1.get().equals(input2.get())) {
				return false;
			}
		}
		return true;
	}

	/**
	 * Put all objects of a model in a map by ID If no ID is specified, an ID is
	 * generated.
	 */
	private void collectObjects(BEASTObject plugin, Map<String, BEASTObject> objects)
			throws IllegalArgumentException, IllegalAccessException {
		for (BEASTObject plugin2 : plugin.listActivePlugins()) {
			if (plugin2.getID() == null) {
				String id = plugin2.getClass().getName();
				if (id.indexOf('.') >= 0) {
					id = id.substring(id.lastIndexOf('.') + 1);
				}
				int i = 0;
				while (objects.keySet().contains(id + i)) {
					i++;
				}
				plugin2.setID(id+i);
			}
			objects.put(plugin2.getID(), plugin2);
			collectObjects(plugin2, objects);
		}
	}
}
