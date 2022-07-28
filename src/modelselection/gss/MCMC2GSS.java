package modelselection.gss;

import beastfx.app.tools.Application;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.Runnable;
import beast.base.core.Input.Validate;
import beast.base.inference.MCMC;

@Description("Convert MCMC analysis to GSS analysis and optionally run the analysis")
public class MCMC2GSS extends MCMC2IS {
	
	public Input<Double> alphaInput = new Input<>("alpha", "alpha parameter of Beta(alpha,1) distribution used to space out steps, default 0.3" +
			"If alpha <= 0, uniform intervals are used.", 0.3);
	public Input<Integer> stepsInput = new Input<>("nrOfSteps", "the number of steps to use, default 8", 8);
	public Input<String> rootDirInput = new Input<>("rootdir", "root directory for storing particle states and log files (default /tmp)", "/tmp");
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
	public Input<Integer> burnInPercentageInput = new Input<>("burnInPercentage", "burn-In Percentage used for analysing log files", 50);
	public Input<Boolean> deleteOldLogsInput = new Input<Boolean>("deleteOldLogs", "delete existing log files from root dir", false);

	@Override
	protected Runnable newInstance(MCMC mcmc) {		
		GeneralisedSteppingStone gss = new GeneralisedSteppingStone();
		gss.setInputValue("nrOfSteps",stepsInput.get());
		gss.setInputValue("rootdir",rootDirInput.get());
		gss.setInputValue("value",m_sScriptInput.get());
		gss.setInputValue("hosts",m_sHostsInput.get());
		gss.setInputValue("burnInPercentage",burnInPercentageInput.get());
		gss.setInputValue("deleteOldLogs",deleteOldLogsInput.get());
		gss.setInputValue("mcmc", mcmc);
		return gss;
	}
	
	public static void main(String[] args) throws Exception {
		new Application(new MCMC2GSS(), "MCMC2GSS", args);
	}

}
