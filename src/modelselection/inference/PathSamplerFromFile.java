package modelselection.inference;


import java.io.IOException;

import javax.xml.parsers.ParserConfigurationException;

import org.xml.sax.SAXException;

import beastfx.app.util.XMLFile;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.MCMC;
import beast.base.core.Input.Validate;
import beast.base.parser.XMLParser;
import beast.base.parser.XMLParserException;

@Description("Path sampler that takes a BEAST MCMC specification from an external file")
public class PathSamplerFromFile extends PathSampler {
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
	}
