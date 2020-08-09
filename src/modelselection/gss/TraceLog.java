package modelselection.gss;

import java.io.IOException;
import java.util.List;

import beast.app.util.LogFile;
import beast.core.BEASTInterface;
import beast.core.BEASTObject;
import beast.core.Description;
import beast.core.Input;
import beast.util.LogAnalyser;
import modelselection.cpo.BEASTRunAnalyser;

@Description("BEAST object wrapper for LogAnalyser")
public class TraceLog extends BEASTObject implements BEASTInterface {
	final public Input<LogFile> traceFileInput = new Input<>("logFile","input file containing tree log.");
	final public Input<Integer> burnInPercentageInput = new Input<Integer>("burnInPercentage", "burn-In Percentage used for analysing log files", 50);
	
	LogAnalyser tracelog;

	
	@Override
	public void initAndValidate() {
		if (tracelog != null) {
			return;
		}
		int burnInPercentage = burnInPercentageInput.get();
		if (burnInPercentage < 0 || burnInPercentage >= 100) {
			throw new IllegalArgumentException("burnInPercentage should be between 0 and 100");
		}

		try {
			tracelog = BEASTRunAnalyser.getTraceLog(null, traceFileInput.get(), burnInPercentage);
		} catch (IOException e) {
			e.printStackTrace();
			throw new IllegalArgumentException(e.getMessage());
		}	
	}

	
	public Double [] getTrace(String label) {
		if (tracelog == null) {
			initAndValidate();
		}
		return tracelog.getTrace(label);
	}


	public List<String> getLabels() {
		if (tracelog == null) {
			initAndValidate();
		}
		return tracelog.getLabels();
	}
	
	public Double getMean(String label) {
		double sum = 0;
		Double [] data = tracelog.getTrace(label);
		for (Double d : data) {
			sum += d;
		}
		return sum / data.length;
	}
}
