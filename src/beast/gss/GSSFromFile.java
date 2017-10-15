package beast.gss;

import java.io.IOException;

import javax.xml.parsers.ParserConfigurationException;

import org.xml.sax.SAXException;

import beast.app.util.Application;
import beast.app.util.XMLFile;
import beast.core.Description;
import beast.core.Input;
import beast.core.MCMC;
import beast.core.Input.Validate;
import beast.util.XMLParser;
import beast.util.XMLParserException;

@Description("Generalised stepping stone that takes a BEAST MCMC specification from an external file")
public class GSSFromFile extends GeneralisedSteppingStone {
		public Input<XMLFile> model1Input = new Input<>(
				"model1",
				"file name of BEAST XML file containing the model for which to run the path sampler",
				new XMLFile("examples/normalTest-1XXX.xml"),
				Validate.REQUIRED);
		
		@Override
		public void initAndValidate() {
			if (!model1Input.get().exists()) {
				return;
			}
			XMLParser parser1 = new XMLParser();
			Object o;
			try {
				o = parser1.parseFile(model1Input.get());
			} catch (SAXException | IOException | ParserConfigurationException | XMLParserException e) {
				throw new IllegalArgumentException(e);
			}
			if (!(o instanceof MCMC)) {
				throw new IllegalArgumentException("The model in " + model1Input.get()
						+ " does not appear to be an MCMC analysis.");
			}
			mcmcInput.setValue(o, this);
			super.initAndValidate();
		}
	
		public static void main(String[] args) throws Exception {
			new Application(new GSSFromFile(), "GSS from file", args);
		}
}
