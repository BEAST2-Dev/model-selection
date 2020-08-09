package modelselection.inference;

import java.io.PrintStream;

import beast.core.BEASTObject;
import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.Loggable;
import beast.core.Input.Validate;

@Description("logs difference of two posteriors for Paired Path sampling")
public class DiffLogger extends BEASTObject implements Loggable {
	public Input<Distribution> post1Input = new Input<>("p1", "posterior of first model", Validate.REQUIRED);
	public Input<Distribution> post2Input = new Input<>("p2", "posterior of first model", Validate.REQUIRED);

	Distribution p1, p2;
	
	public DiffLogger() {}
	public DiffLogger(Distribution p1, Distribution p2) {
		initByName("p1", p1, "p2", p2);
	}
	
	@Override
	public void initAndValidate() {
		p1 = post1Input.get();
		p2 = post2Input.get();
	}

	@Override
	public void init(PrintStream out) {
		out.append("diff-posterior\t");
	}


	@Override
	public void log(long nSample, PrintStream out) {
		out.append(p1.getCurrentLogP() - p2.getCurrentLogP() + "\t");
	}


	@Override
	public void close(PrintStream out) {
		// nothing to do
	}

}
