package beast.app.tools;

import beast.app.util.Application;

// command line interface to PairedPathSampler
public class PairedPathSampler {
	
	public static void main(final String[] args) throws Exception {
    	Application main = null;
        try {
    		beast.inference.PairedPathSampler sampler = new beast.inference.PairedPathSampler();
                main = new Application(sampler);
                main.parseArgs(args, false);
                sampler.initAndValidate();
                sampler.run();
        } catch (Exception e) {
                System.out.println(e.getMessage());
                if (main != null) {
                	System.out.println(main.getUsage());
                }
        }
    }

}
