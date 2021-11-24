package gedi.metagene;

import gedi.util.StringUtils;
import gedi.util.dynamic.impl.ArrayDynamicObject;
import gedi.util.dynamic.impl.BooleanDynamicObject;
import gedi.util.dynamic.impl.DoubleDynamicObject;
import gedi.util.functions.EI;
import gedi.util.program.GediProgram;
import gedi.util.program.GediProgramContext;
import gedi.util.r.RRunner;

import java.io.File;
import java.io.IOException;

public class MetagenePlotter extends GediProgram {
    public MetagenePlotter(TiSSMetagenePlotterParameterSet params) {
        addInput(params.prefix);
        addInput(params.binSizes);
        addInput(params.binNames);
        addInput(params.plot);
        addInput(params.noHeatmaps);
        addInput(params.plotYMax);
        addInput(params.namesFile);
        addInput(params.voidout);

        addOutput(params.voidout2);
    }

    @Override
    public String execute(GediProgramContext context) throws Exception {
        String prefix = getParameter(0);
        String binsizes = getParameter(1);
        String inputBinNames = getParameter(2);
        boolean plot = getParameter(3);
        boolean noHeatmaps = getParameter(4);
        double yMaxValue = getParameter(5);
        String namesFile = getParameter(6);

        String[] names = null;
        if (namesFile != null && !namesFile.isEmpty()) {
            names = EI.lines(namesFile).toArray(String.class);
        }

        String[] binNames = StringUtils.split(inputBinNames, ",");
        int[] bins = EI.wrap(StringUtils.split(binsizes, ",")).mapToInt(Integer::parseInt).toIntArray();
        Integer[] binSizes = new Integer[bins.length];
        for (int i = 0; i < bins.length; i++) {
            binSizes[i] = bins[i];
        }

        plotDataSummed(prefix, plot, yMaxValue);
        plotDataIndividual(prefix, plot, binNames, binSizes, noHeatmaps, names);
        return null;
    }

    private void plotDataSummed(String prefix, boolean plot, double yMaxValue) throws IOException {
        new File(prefix + "/plots/summed").mkdirs();
        RRunner r = new RRunner(prefix+"/plotSummed.R");
        r.set("prefixIn", prefix + "/data/summed");
        r.set("prefixOut", prefix + "/plots/summed");
        r.set("yMaxValue", new DoubleDynamicObject(yMaxValue));
        r.addSource(getClass().getResourceAsStream("/resources/metagene/plotMetageneSummed.R"));
        if (plot) {
            r.run(false);
        } else {
            r.dontrun(false);
        }
    }

    private void plotDataIndividual(String prefix, boolean plot, String[] binNames, Integer[] binSizes, boolean noHeatmaps, String[] names) throws IOException {
        RRunner r = new RRunner(prefix+"/plotIndividual.R");
        new File(prefix + "/plots/individual").mkdirs();
        r.set("prefixIn", prefix + "/data/individual");
        r.set("prefixOut", prefix + "/plots/individual");
        r.set("binNames", new ArrayDynamicObject(binNames));
        r.set("binSizes", new ArrayDynamicObject(binSizes));
        if (names != null) {
            r.set("names", new ArrayDynamicObject(names));
        } else {
            r.set("names", new ArrayDynamicObject(new String[] {}));
        }
        r.addSource(getClass().getResourceAsStream("/resources/metagene/plotMetageneIndividual.R"));
        if (plot && !noHeatmaps) {
            r.run(false);
        } else {
            r.dontrun(false);
        }
    }
}
