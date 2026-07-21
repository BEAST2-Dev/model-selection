package modelselection.gss;

import beast.base.inference.MCMC;
import beast.base.inference.Runnable;

/**
 * A module-info.java file cannot contain an abstract class
 */
public abstract class MCMC2Abstract extends Runnable {
    protected abstract Runnable newInstance(MCMC mcmc);
}
