package gedi.metagene;

import gedi.core.data.reads.AlignedReadsData;
import gedi.core.data.reads.ReadCountMode;
import gedi.core.reference.Strandness;
import gedi.core.region.GenomicRegionStorage;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.util.FileUtils;
import gedi.util.StringUtils;
import gedi.util.functions.EI;
import gedi.util.io.text.LineOrientedFile;
import gedi.util.io.text.LineWriter;
import gedi.util.program.GediProgram;
import gedi.util.program.GediProgramContext;
import gedi.util.userInteraction.progress.ConsoleProgress;
import gedi.utils.ReadType;
import gedi.utils.WriterUtils;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;

public class MetageneDataExtractor extends GediProgram {
    public MetageneDataExtractor(TiSSMetagenePlotterParameterSet params) {
        addInput(params.prefix);
        addInput(params.readData);
        addInput(params.strandness);
        addInput(params.binSizes);
        addInput(params.binNames);
        addInput(params.locationFile);
        addInput(params.mergeConditions);
        addInput(params.normalizationMode);
        addInput(params.totalReadCounts);
        addInput(params.readType);
        addInput(params.readCountMode);
        addInput(params.perBinNormMode);
        addInput(params.useLog);

        addOutput(params.voidout);
    }

    @Override
    public String execute(GediProgramContext context) throws Exception {
        String prefix = getParameter(0);
        GenomicRegionStorage<AlignedReadsData> cit = getParameter(1);
        Strandness strandness = getParameter(2);
        String binsizes = getParameter(3);
        String inputBinNames = getParameter(4);
        String locationFile = getParameter(5);
        String toMerge = getParameter(6);
        NormalizationMode normalizationMode = getParameter(7);
        List<Double> totalsIn = getParameters(8);
        ReadType readType = getParameter(9);
        ReadCountMode readCountMode = getParameter(10);
        BinNormalizationMode perBinNormMode = getParameter(11);
        boolean useLog = getParameter(12);

        List<MetageneLocation[]> entries = readLocationFile(locationFile);

        double[] totals = new double[totalsIn.size()];
        for (int i = 0 ; i< totalsIn.size(); i++) {
            totals[i] = totalsIn.get(i);
        }

        int[][] mergeConditions = EI.wrap(StringUtils.split(toMerge, "/")).map(m -> EI.wrap(StringUtils.split(m, ",")).mapToInt(Integer::parseInt).toIntArray()).toArray(new int[0][]);
        String[] conditions = cit.getMetaDataConditions();
        String[] binNames = StringUtils.split(inputBinNames, ",");
        int[] bins = EI.wrap(StringUtils.split(binsizes, ",")).mapToInt(Integer::parseInt).toIntArray();

        String[] mergedConditions = null;
        if (mergeConditions != null && mergeConditions.length != 0) {
            mergedConditions = new String[mergeConditions.length];
            for (int i = 0; i < mergeConditions.length; i++) {
                StringBuilder sb = new StringBuilder();
                sb.append(conditions[mergeConditions[i][0]]);
                for (int j = 1; j < mergeConditions[i].length; j++) {
                    sb.append("_").append(conditions[mergeConditions[i][j]]);
                }
                mergedConditions[i] = sb.toString();
            }
        }

        MetageneBinManager binManager = new MetageneBinManager(mergedConditions == null ? conditions : mergedConditions, binNames, bins);
        EI.wrap(entries).progress(new ConsoleProgress(System.err), entries.size(), d -> "Processing...").forEachRemaining(rgr -> {
            MetageneBinner metageneBinner = new MetageneBinner(rgr, bins);
            metageneBinner.extractReadCounts(cit, strandness, totals, readType, readCountMode);
            if (mergeConditions != null && mergeConditions.length != 0) {
                metageneBinner.mergeReadCounts(mergeConditions);
            }
            metageneBinner.normalizeReadCount(normalizationMode);
            metageneBinner.mergeReadcountsIntoBins(perBinNormMode);
            binManager.addMetageneBin(metageneBinner);
        });

        writeOutData(prefix + "/data", mergedConditions == null ? conditions : mergedConditions, binManager, bins, binNames);

        LineWriter w = new LineOrientedFile(getOutputFile(0).getPath()).write();
        w.writeLine("Number of genes used: " + entries.size());
        w.close();
        return null;
    }

    private void writeOutData(String prefix, String[] conditions, MetageneBinManager binManager, int[] bins, String[] binNames) throws IOException {
        LineWriter[] individualDataWriter = WriterUtils.createWriters(prefix + "/individual/", conditions, ".tsv");
        LineWriter[] summedDataWriter = WriterUtils.createWriters(prefix + "/summed/", conditions, ".tsv");

        WriterUtils.writeLineToAll(summedDataWriter, "Position\tType\tValue");

        for (int c = 0; c < individualDataWriter.length; c++) {
            final int writerIndex = c;
            binManager.iterateIndividualEntries(c).forEachRemaining(entry -> {
                writeIndividualEntry(individualDataWriter[writerIndex], entry);
            });
        }

        for (int c = 0; c < summedDataWriter.length; c++) {
            int index = 1;
            for (int bin = 0; bin < bins.length; bin++) {
                double[] sum = binManager.getMeanData(c, bin);
                for (double value : sum) {
                    summedDataWriter[c].writeLine(index + "\t" + binNames[bin] + "\t" + value);
                    index++;
                }
            }
        }

        WriterUtils.closeWriters(individualDataWriter);
        WriterUtils.closeWriters(summedDataWriter);
    }

    private static void writeIndividualEntry(LineWriter writer, double[][] entry) {
        StringBuilder line = new StringBuilder();
        line.append(entry[0][0]);
        boolean first = true;
        for (double[] bin : entry) {
            for (int i = 0; i < bin.length; i++) {
                if (first) {
                    first = false;
                    continue;
                }
                line.append("\t").append(bin[i]);
            }
        }
        try {
            writer.writeLine(line.toString());
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static List<MetageneLocation[]> readLocationFile(String file) throws IOException {
        List<MetageneLocation[]> lst = new ArrayList<>();
        EI.lines(file).forEachRemaining(l -> {
            String[] split = StringUtils.split(l, "\t");
            MetageneLocation[] rgr = new MetageneLocation[split.length];
            for (int i = 0; i < split.length; i++) {
                rgr[i] = MetageneLocation.parse(split[i]);
            }
            lst.add(rgr);
        });
        return lst;
    }
}
