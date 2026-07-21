open module beast.modelselection {
    requires beast.pkgmgmt;
    requires beast.base;
    requires java.xml;
    requires static beast.fx;
    requires static javafx.controls;

    requires org.apache.commons.statistics.distribution;

    exports modelselection.app.tools;
    exports modelselection.core;
    exports modelselection.cpo;
    exports modelselection.gss;
    exports modelselection.gss.coalescent;
    exports modelselection.gss.distribution;
    exports modelselection.inference;

    // Tell the module system to consume this abstract type
    uses modelselection.gss.MCMC2Abstract;

    provides beast.base.core.BEASTInterface with
        modelselection.core.CPOLogger,
        modelselection.cpo.BEASTRunAnalyser,
        modelselection.cpo.CPOAnalyser,
        modelselection.gss.GSSFromFile,
        modelselection.gss.GeneralisedSteppingStone,
        modelselection.gss.GeneralisedSteppingStoneStep,
        modelselection.gss.MCMC2GSS,
        modelselection.gss.MCMC2IS,
        modelselection.gss.TraceLog,
        modelselection.gss.TreeFromTreeSetFileInitialiser,
        modelselection.gss.coalescent.ExponentialProductPosteriorMeansLikelihood,
        modelselection.gss.distribution.GSSTreeDistribution,
        modelselection.gss.distribution.LogTransformedNormalKDEDistribution,
        modelselection.gss.distribution.LogitTransformedNormalKDEDistribution,
        modelselection.gss.distribution.MultivariateKDEDistribution,
        modelselection.gss.distribution.NormalKDEDistribution,
        modelselection.inference.AICMAnalyser,
        modelselection.inference.DiffLogger,
        modelselection.inference.PairedPathSampleAnalyser,
        modelselection.inference.PairedPathSampler,
        modelselection.inference.PairedPathSamplingStep,
        modelselection.inference.PathSampleAnalyser,
        modelselection.inference.PathSampler,
        modelselection.inference.PathSamplerFromFile,
        modelselection.inference.PathSamplingStep;
}
