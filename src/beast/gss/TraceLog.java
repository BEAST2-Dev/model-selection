package beast.gss;

import java.io.IOException;

import beast.app.util.LogFile;
import beast.core.BEASTInterface;
import beast.core.BEASTObject;
import beast.core.Input;
import beast.cpo.BEASTRunAnalyser;
import beast.util.LogAnalyser;

public class TraceLog extends BEASTObject implements BEASTInterface {
	final public Input<LogFile> traceFileInput = new Input<>("logFile","input file containing tree log.");
	final public Input<Integer> burnInPercentageInput = new Input<Integer>("burnInPercentage", "burn-In Percentage used for analysing log files", 50);
	
	LogAnalyser tracelog;
	
	@Override
	public void initAndValidate() {
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
		return tracelog.getTrace(label);
	}
}
