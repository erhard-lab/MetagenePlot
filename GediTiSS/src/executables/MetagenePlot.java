package executables;

import gedi.metagene.MetageneDataExtractor;
import gedi.metagene.MetagenePlotter;
import gedi.metagene.TiSSMetagenePlotterParameterSet;
import gedi.util.program.CommandLineHandler;
import gedi.util.program.GediProgram;

public class MetagenePlot {
    public static void main(String[] args) {
        TiSSMetagenePlotterParameterSet params = new TiSSMetagenePlotterParameterSet();

        GediProgram tiss = GediProgram.create("MetagenPlot",
                new MetageneDataExtractor(params),
                new MetagenePlotter(params));

        GediProgram.run(tiss, params.paramFile, new CommandLineHandler("MetagenPlot",
                "Plots metagene plots in heatmap and graph form", args));
    }
}
