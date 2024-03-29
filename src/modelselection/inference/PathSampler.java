package modelselection.inference;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.PrintStream;
import java.io.UnsupportedEncodingException;
import java.net.URISyntaxException;
import java.text.DecimalFormat;
import java.util.List;
import java.util.concurrent.CountDownLatch;

import org.apache.commons.math.distribution.BetaDistribution;
import org.apache.commons.math.distribution.BetaDistributionImpl;

import beastfx.app.util.Utils;
import beast.base.core.Description;
import beast.base.inference.Distribution;
import beast.base.core.Input;
import beast.base.inference.Logger;
import beast.base.inference.MCMC;
import beast.base.core.Input.Validate;
import beast.base.inference.CompoundDistribution;
import beast.base.core.Log;
import beast.base.core.ProgramStatus;
import beast.base.util.Randomizer;
import beast.base.parser.XMLProducer;



@Description("Calculate marginal likelihood through path/stepping stone sampling. " +
		"Perform multiple steps and calculate estimate." +
		"Uses multiple threads if specified as command line option to BEAST.")
public class PathSampler extends beast.base.inference.Runnable {
	public static String LIKELIHOOD_LOG_FILE = "likelihood.log";

	public Input<Double> alphaInput = new Input<Double>("alpha", "alpha parameter of Beta(alpha,1) distribution used to space out steps, default 0.3" +
			"If alpha <= 0, uniform intervals are used.", 0.3);
	public Input<Integer> stepsInput = new Input<Integer>("nrOfSteps", "the number of steps to use, default 8", 8);
	public Input<String> rootDirInput = new Input<String>("rootdir", "root directory for storing particle states and log files (default /tmp)", "/tmp");
	public Input<MCMC> mcmcInput = new Input<MCMC>("mcmc", "MCMC analysis used to specify model and operations in each of the particles", Validate.REQUIRED);
	public Input<Long> chainLengthInput = new Input<>("chainLength", "number of sample to run a chain for a single step", 100000L);
	public Input<Integer> burnInPercentageInput = new Input<Integer>("burnInPercentage", "burn-In Percentage used for analysing log files", 50);
	public Input<Integer> preBurnInInput = new Input<>("preBurnin", "number of samples that are discarded for the first step, but not the others", 100000);
	public Input<String> m_sScriptInput = new Input<String>("value", "script for launching a job. " +
			"$(dir) is replaced by the directory associated with the particle " +
			"$(java.class.path) is replaced by a java class path used to launch this application " +
			"$(java.library.path) is replaced by a java library path used to launch this application " +
			"$(java) is replaced by a path to the java executable included in the BEAST release " +
			"$(seed) is replaced by a random number seed that differs with every launch " +
			"$(host) is replaced by a host from the list of hosts", Validate.REQUIRED);
	public Input<String> m_sHostsInput = new Input<String>("hosts", "comma separated list of hosts. " +
			"If there are k hosts in the list, for particle i the term $(host) in the script will be replaced " +
			"by the (i modulo k) host in the list. " +
			"Note that whitespace is removed. " +
			"This can be useful for example when starting BEAST remotely through ssh.");
	public Input<Boolean> doNotRun = new Input<Boolean>("doNotRun", "Set up all files but do not run analysis if true. " +
			"This can be useful for setting up an analysis on a cluster", false);
	
	public Input<Boolean> deleteOldLogsInpuyt = new Input<Boolean>("deleteOldLogs", "delete existing log files from root dir", false);
	
	public Input<Boolean> posterior2priorInput = new Input<Boolean>("posterior2prior", "whether to do steps from posterior to prior or the other way around. "
			+ "Going from posterior to prior is biased towards over estimates, while from prior to posterior the ML estimate "
			+ "is biased towards under estimates.", true);
	
	int m_nSteps;
	String [] m_sHosts;
	String m_sScript;
	int burnInPercentage;

	CountDownLatch m_nCountDown;

    final static String fileSep = System.getProperty("file.separator");

	DecimalFormat formatter;
	String getStepDir(int iParticle) {
		if (rootDirInput.get().equals("")) {
			Log.warning("WARNING: empty rootdir found. This means the top level of the file system "
					+ "will be used for the step directories. Consider specifying the rootdir attribute "
					+ "with a non-empty string.");
		}
		File f = new File(rootDirInput.get());
		return f.getAbsolutePath() + "/step" + formatter.format(iParticle);
	}


	@Override
	public void initAndValidate() {
	}
	
	@Override
	public void run() throws Exception {
		if (rootDirInput.get().startsWith("\\\\")) {
			Log.warning.println("WARNING: the rootdir starts with \\\\ suggesting this refers to a network drive in Windows. "
					+ "This probably does not work in Windows: use a local drive, or use map the share to a drive and use the "
					+ "corresponding drive letter instead.");
		}
		// grab info from inputs
		m_sScript = m_sScriptInput.get();
		if (m_sScript == null) {
			m_sScript = "cd $(dir)\n" +
					"$(java) -cp $(java.class.path) beast.pkgmgmt.launcher.BeastLauncher $(resume/overwrite) -java -seed $(seed) beast.xml\n";
		}
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
			// try to create directory
			if (!rootDir.mkdirs()) {
				throw new Exception("Directory " + rootDirInput.get() + " does not exist and could not be created.");
			} else {
				Log.warning.println("Created directory " + rootDir.getAbsolutePath());
			}
		}
		if (!rootDir.isDirectory()) {
			throw new Exception(rootDirInput.get() + " is not a directory.");
		}
		
		// initialise MCMC
		MCMC mcmc = mcmcInput.get();

		if (!mcmc.getClass().equals(MCMC.class)) {
			System.out.println("WARNING: class is not beast.base.inference.MCMC, which may result in unexpected behavior ");
		}
		PathSamplingStep step = new PathSamplingStep();
		for (Input<?> input : mcmc.listInputs()) {
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
		mcmc = step;
		
		long chainLength = chainLengthInput.get();
		// set up chain length for a single step
		mcmc.burnInInput.setValue(0, mcmc);
		mcmc.chainLengthInput.setValue(chainLength, mcmc);
		
		// add posterior logger
		Logger logger = new Logger();
		Distribution likelihood = extractLikelihood(mcmc); 
		logger.initByName("fileName", LIKELIHOOD_LOG_FILE, "log", likelihood, "logEvery", (int)(chainLength/1000));
		mcmc.loggersInput.setValue(logger, mcmc);

		// set up directories with beast.xml files in each of them
		String sFormat = "";
		for (int i = m_nSteps; i > 0; i /= 10) {
			sFormat += "#";
		}
		formatter = new DecimalFormat(sFormat);
		
		XMLProducer producer = new XMLProducer();
		BetaDistribution betaDistribution = null;
		if (alphaInput.get() > 0){
			betaDistribution = new BetaDistributionImpl(alphaInput.get(), 1.0);
		}
		
		
		PrintStream [] cmdFiles = new PrintStream[ProgramStatus.m_nThreads];
    	for (int i = 0; i < ProgramStatus.m_nThreads; i++) {
    		FileOutputStream outStream = (Utils.isWindows()?
    					new FileOutputStream(rootDirInput.get() + "/run" + i +".bat"):
    					new FileOutputStream(rootDirInput.get() + "/run" + i +".sh"));
    		 cmdFiles[i] = new PrintStream(outStream);
    	}

		
		
		for (int i = 0; i < m_nSteps; i++) {
			if (i < ProgramStatus.m_nThreads) {
				mcmc.burnInInput.setValue(preBurnIn, mcmc);
			} else {
				mcmc.burnInInput.setValue(0, mcmc);
			}
			// create XML for a single step
			double beta = betaDistribution != null ?
					(posterior2priorInput.get() ?
					betaDistribution.inverseCumulativeProbability((m_nSteps - 1.0 - i)/ (m_nSteps - 1)):
					betaDistribution.inverseCumulativeProbability((i+0.0)/ (m_nSteps - 1))
					):(m_nSteps - 1.0 - i)/ (m_nSteps - 1);
			step.setInputValue("beta", beta);
			String sXML = producer.toXML(step);
			File stepDir = new File(getStepDir(i));
			if (!stepDir.exists() && !stepDir.mkdir()) {
				throw new Exception("Failed to make directory " + stepDir.getName());
			}
			stepDir.setWritable(true, false);
        	FileOutputStream xmlFile = new FileOutputStream(stepDir.getAbsoluteFile() + "/beast.xml");
        	PrintStream out = new PrintStream(xmlFile);
            out.print(sXML);
			out.close();
			
			String cmd = getCommand(stepDir.getAbsolutePath(), i);
        	FileOutputStream cmdFile = 
        			(Utils.isWindows()?
        					new FileOutputStream(stepDir.getAbsoluteFile() + "/run.bat"):
        					new FileOutputStream(stepDir.getAbsoluteFile() + "/run.sh"));
        	PrintStream out2 = new PrintStream(cmdFile);
            out2.println(cmd);
			out2.close();

        	cmdFile = 
        			(Utils.isWindows()?
        					new FileOutputStream(stepDir.getAbsoluteFile() + "/resume.bat"):
        					new FileOutputStream(stepDir.getAbsoluteFile() + "/resume.sh"));
        	cmd = cmd.replace("-overwrite", "-resume");
        	out2 = new PrintStream(cmdFile);
            out2.println(cmd);
			out2.close();
//TODO: probably more efficient to group cmdFiles in block of #steps/#threads
//instead of skipping #threads steps every time.
			if (i >= ProgramStatus.m_nThreads) {
				String copyCmd = (Utils.isWindows()
						? getCopyCmd(i)
						: "cp " + getStepDir(i - ProgramStatus.m_nThreads) + "/beast.xml.state " + getStepDir(i)
							);
				cmdFiles[i % ProgramStatus.m_nThreads].println(copyCmd);				
			}
			if (i / ProgramStatus.m_nThreads == 0) {
				cmd = cmd.replace("-resume", "-overwrite");
			}
			cmdFiles[i % ProgramStatus.m_nThreads].println(cmd);
			File script = new File(stepDir.getAbsoluteFile() + 
					(Utils.isWindows()? "/run.bat": "/run.sh"));
			script.setExecutable(true);
		}
		
		
		
    	for (int k = 0; k < ProgramStatus.m_nThreads; k++) {
    		cmdFiles[k].close();
    	}

    	doRuns();
		
	} // run
	
	private String getCopyCmd(int i) {
		String cmd = "copy \"";
		cmd += getStepDir(i - ProgramStatus.m_nThreads);
		cmd += "\\beast.xml.state\" \"";
		cmd += getStepDir(i);
		cmd += "\"";
		StringBuilder buf = new StringBuilder();
		for (char a : cmd.toCharArray()) {
			if (a == '/') {a = '\\';}
			buf.append(a);
		}
		cmd = buf.toString();
		return cmd;
	}


	public static Distribution extractLikelihood(MCMC mcmc) {
		Distribution posterior = mcmc.posteriorInput.get();
		// expect compound distribution with likelihood and prior
		if (!(posterior instanceof CompoundDistribution)) {
			throw new IllegalArgumentException("Expected posterior being a CompoundDistribution");
		}
		CompoundDistribution d = (CompoundDistribution) posterior;
		List<Distribution> list = d.pDistributions.get();
		if (list.size() < 2) {
            throw new IllegalArgumentException("Expected one likelihood and at least one prior distribution.");
		}

		Distribution[] pDists = new Distribution[list.size()];
		int nextPriorIndex = 1;
        for (Distribution pDist: list) {
            final String distID = pDist.getID().toLowerCase();
            if (distID.startsWith("likelihood")) {
                if (pDists[0] == null) pDists[0] = pDist; 
                else throw new IllegalArgumentException("Expected only one likelihood distribution.");
            } else {
                pDists[nextPriorIndex] = pDist;
                nextPriorIndex++;
            }
        }

        return pDists[0];
	}

	String getCommand(String sStepDir, int iStep) {
		sStepDir = sStepDir.replace("\\", "\\\\");
		String sCommand = m_sScript.replaceAll("\\$\\(dir\\)", "\"" + sStepDir + "\"");
		//while (sCommand.matches("$(seed)")) {
			sCommand = sCommand.replaceAll("\\$\\(seed\\)", Math.abs(Randomizer.nextInt())+"");
		//}
		String javaInJrePath = null;
		try {
			javaInJrePath = getJavaInJrePath();
		} catch (UnsupportedEncodingException | URISyntaxException e) {
			e.printStackTrace();
		}	
		sCommand = sCommand.replaceAll("\\$\\(java\\)",  javaInJrePath != null ? "\"" + javaInJrePath + "\"" : "java");
		sCommand = sCommand.replaceAll("\\$\\(java.library.path\\)",  "\"" + sanitise(System.getProperty("java.library.path")) + "\"");
//		sCommand = sCommand.replaceAll("\\$\\(java.class.path\\)", "\"" + sanitise(System.getProperty("java.class.path")) + "\"");
		sCommand = sCommand.replaceAll("\\$\\(java.class.path\\)", "\"" + sanitise(getLauncherJarPath()) + "\"");
		sCommand = sCommand.replaceAll("beastfx.app.beast.BeastMain", "beast.pkgmgmt.launcher.BeastLauncher");
		sCommand = sCommand.replaceAll("beast.app.beastapp.BeastMain", "beast.pkgmgmt.launcher.BeastLauncher");
		if (m_sHosts != null) {
			sCommand = sCommand.replaceAll("\\$\\(host\\)", m_sHosts[iStep % m_sHosts.length]);
		}
		if (iStep < ProgramStatus.m_nThreads) {
			sCommand = sCommand.replaceAll("\\$\\(resume/overwrite\\)", "-overwrite");
		} else {
			sCommand = sCommand.replaceAll("\\$\\(resume/overwrite\\)", "-resume");
		}		
				
		return sCommand;
	}

	
	private String getJavaInJrePath() throws URISyntaxException, UnsupportedEncodingException {
		File launcherJarFile = new File(getLauncherJarPath());
		String javaInJrePath = launcherJarFile.getParentFile().getParentFile().getAbsolutePath() + "/jre/bin/java";
				
		if (Utils.isWindows()) {
			javaInJrePath += ".exe";
		}
		if (new File(javaInJrePath).exists()) {
			return javaInJrePath;
		}
		return null;
	}


	private String getLauncherJarPath() {
		String property = System.getProperty("java.class.path");
		String pathSeparator = System.getProperty("path.separator");
		String [] paths = property.split(pathSeparator);
		for (String path : paths) {
			if (path.toLowerCase().endsWith("launcher.jar")) {
				return path;
			}
		}
		return null;
	}


	private String sanitise(String property) {
		// make absolute paths from relative paths
		String pathSeparator = System.getProperty("path.separator");
		String [] paths = property.split(pathSeparator);
		StringBuilder b = new StringBuilder();
		for (String path : paths) {
			File f = new File(path);
			b.append(f.getAbsolutePath());
			b.append(pathSeparator);
		}
		// chop off last pathSeparator
		property = b.substring(0, b.length() - 1);
		
		// sanitise for windows
		if (Utils.isWindows()) {
			String cwd = System.getProperty("user.dir");
			cwd = cwd.replace("\\", "/");
			property = property.replaceAll(";\\.", ";" +  cwd + ".");
			property = property.replace("\\", "/");
		}
		return property;
	}


	class StepThread extends Thread {
		int stepNr;
		
		StepThread(int stepNr) {
			this.stepNr = stepNr;
		}
		
		@Override
		public void run() {
			try {
				System.err.println("Starting step " + stepNr);
				File stepDir = new File(getStepDir(stepNr));
				if (!stepDir.exists()) {
					throw new Exception("Failed to find directory " + stepDir.getName());
				}
	        	String cmd = 
        			(Utils.isWindows()?
        					stepDir.getAbsoluteFile() + "/run.bat":
        					stepDir.getAbsoluteFile() + "/run.sh");
	        	
				ProcessBuilder pb = new ProcessBuilder(cmd);
				pb.redirectErrorStream(true); // merge stdout and stderr
				Process p = pb.start();
//	        	
//				Process p = Runtime.getRuntime().exec(cmd);
				BufferedReader pout = new BufferedReader((new InputStreamReader(p.getInputStream())));
				String line;
				while ((line = pout.readLine()) != null) {
					//System.out.println(line);
				}
				pout.close();
				
				p.waitFor();
			} catch (Exception e) {
				e.printStackTrace();
			}
			System.err.println("Finished step " + stepNr);
			m_nCountDown.countDown();
		}
	}
	
    public void doRuns() throws Exception {
    	if (doNotRun.get()) {
    		printDoNotRunMessage();
    		return;
    	}
    	long startTime = System.currentTimeMillis();

		for (int i = 0; i < m_nSteps; i++) {
	    	if (ProgramStatus.m_nThreads > 1) {
	    		int nSteps = Math.min(ProgramStatus.m_nThreads, m_nSteps - i);
	    		m_nCountDown = new CountDownLatch(nSteps);
	    		for (int j = 0; j < nSteps; j++) {
	    			new StepThread(i).start();
	    			i++;
	    		}
	    		m_nCountDown.await();	    		
    			for (int j = 0; j < ProgramStatus.m_nThreads && i+j < m_nSteps; j++) {
    				copyStateFile(i-1, i + j);	    			
    				checkLogFiles(i + j);
    			}
	    		i--;
	    	} else {
				File stepDir = new File(getStepDir(i));
				if (!stepDir.exists()) {
					throw new Exception("Failed to find directory " + stepDir.getName());
				}
	        	String cmd = 
	        			(Utils.isWindows()?
	        					stepDir.getAbsoluteFile() + "\\run.bat":
	        					stepDir.getAbsoluteFile() + "/run.sh");
	        	
	        	
//				try {
	    		if (i > 0) {
	    			copyStateFile(i-1, i);
	    		}
				checkLogFiles(i);
					
					System.out.println(cmd);
					
					ProcessBuilder pb = new ProcessBuilder(cmd);
					pb.redirectErrorStream(true); // merge stdout and stderr
					Process p = pb.start();
					
					BufferedReader pout = new BufferedReader((new InputStreamReader(p.getInputStream())));
					String line;
					while ((line = pout.readLine()) != null) {
						System.out.println(line);
					}
					pout.close();
					p.waitFor();
//				} catch (Exception e) {
//					e.printStackTrace();
//				}
	    	}
		}
    	long endTime = System.currentTimeMillis();

    	analyse();

		System.out.println("\n\nTotal wall time: " + (endTime-startTime)/1000 + " seconds\nDone");
    } // run;	


	void analyse() throws Exception {
    	PathSampleAnalyser analyser = new PathSampleAnalyser();
    	double marginalL = analyser.estimateMarginalLikelihood(m_nSteps, alphaInput.get(), rootDirInput.get(), burnInPercentage);
		System.out.println("marginal L estimate = " + marginalL);
	}


	void printDoNotRunMessage() {
		System.out.println("batch files can be found in " + rootDirInput.get());
		System.out.println("Run these and then run"); 
		System.out.println("java beast.inference.PathSampleAnalyser " + m_nSteps + " " + alphaInput.get() + " " + rootDirInput.get() + " " + burnInPercentage);
	}


	/** check for log files in directory for step i **/
	private void checkLogFiles(int i) throws Exception {
		File stepDir = new File(getStepDir(i));
		
		// remove any existing likglihood.log file
		File logFile = new File(stepDir.getPath() + fileSep + "likelihood.log");
		if (logFile.exists()) {
			if (deleteOldLogsInpuyt.get()) {
				System.err.println("WARNING: deleting file " + logFile.getPath());
				logFile.delete();
			} else {
				throw new Exception("Found old log file " + logFile.getPath() + " and will not overwrite (unless deleteOldLogs flag is set to true)");
			}
		}
		
		// process other log and tree files
		for (File file : stepDir.listFiles()) {
			if (file.getPath().endsWith(".log") || 
					file.getPath().endsWith(".trees")) {
				if (deleteOldLogsInpuyt.get()) {
				System.err.println("WARNING: deleting file " + file.getPath());
					file.delete();
				} else {
					throw new Exception("Found old log file " + file.getPath() + " and will not overwrite (unless deleteOldLogs flag is set to true)");
				}
			}
		}
	}

	/** copy beast.xml.state file from previous directory **/
	private void copyStateFile(int iFrom, int iTo) throws Exception {
		File prevStepDir = new File(getStepDir(iFrom));
		File stepDir = new File(getStepDir(iTo));
        InputStream in = new FileInputStream(new File(prevStepDir.getPath() + fileSep + "beast.xml.state"));
        OutputStream out = new FileOutputStream(new File(stepDir.getPath() + fileSep + "beast.xml.state"));
        byte[] buf = new byte[1024];
        int len;
        while ((len = in.read(buf)) > 0) {
            out.write(buf, 0, len);
        }
        in.close();
        out.close();
	}

	@Override
	public boolean hasPartitions() {
		return false;
	}
} // PathSampler
