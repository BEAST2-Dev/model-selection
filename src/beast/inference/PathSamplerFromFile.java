package beast.inference;

import java.io.File;

import beast.core.Description;
import beast.core.Input;
import beast.core.MCMC;
import beast.core.Input.Validate;
import beast.util.XMLParser;

@Description("Path sampler that takes a BEAST MCMC specification from an external file")
public class PathSamplerFromFile extends beast.inference.PathSampler {
		public Input<File> model1Input = new Input<File>(
				"model1",
				"file name of BEAST XML file containing the model for which to run the path sampler",
				new File("examples/normalTest-1XXX.xml"),
				Validate.REQUIRED);
		
		@Override
		public void initAndValidate() throws Exception {
			if (!model1Input.get().exists()) {
				return;
			}
			XMLParser parser1 = new XMLParser();
			Object o = parser1.parseFile(model1Input.get());
			if (!(o instanceof MCMC)) {
				throw new Exception("The model in " + model1Input.get()
						+ " does not appear to be an MCMC analysis.");
			}
			mcmcInput.setValue(o, this);
			super.initAndValidate();
		}
	}
