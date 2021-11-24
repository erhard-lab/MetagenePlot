package gedi.metagene;

import gedi.core.data.reads.AlignedReadsData;
import gedi.core.data.reads.ReadCountMode;
import gedi.core.reference.Strandness;
import gedi.core.region.GenomicRegionStorage;
import gedi.util.program.GediParameter;
import gedi.util.program.GediParameterSet;
import gedi.util.program.parametertypes.*;
import gedi.utils.ReadType;

import java.io.File;

public class TiSSMetagenePlotterParameterSet extends GediParameterSet {
    public GediParameter<String> prefix = new GediParameter<String>(this,"prefix", "The prefix used for all output files", false, new StringParameterType());
    public GediParameter<File> paramFile = new GediParameter<File>(this,"${prefix}.param", "File containing the parameters used to call TiSS", false, new FileParameterType());
    public GediParameter<Strandness> strandness = new GediParameter<Strandness>(this,"strandness", "Strandness of the data", false, new EnumParameterType<>(Strandness.class), Strandness.Sense, true);
    public GediParameter<GenomicRegionStorage<AlignedReadsData>> readData = new GediParameter<GenomicRegionStorage<AlignedReadsData>>(this,"readData", "Read data in CIT-format.", false, new StorageParameterType<AlignedReadsData>());

    public GediParameter<String> binSizes = new GediParameter<String>(this,"binSizes", "The sizes of each bin comma separated (e.g.: 40,60,20 this creates three bins with a total of 120 bins)", false, new StringParameterType());
    public GediParameter<String> binNames = new GediParameter<String>(this,"binNames", "Comma separated string containing names for the bins. Need to be the same size (i.e. same amounts of commas as binSizes)", false, new StringParameterType());
    public GediParameter<String> mergeConditions = new GediParameter<String>(this,"mergeConds", "Conditions to merge. Write as 0,2/1,3/4,5 to merge conditions 0 and 2, 1 and 3 and so on.", false, new StringParameterType(), "", true);
    public GediParameter<String> locationFile = new GediParameter<String>(this,"locFile", "File containing the locations to check for. The number of locations per line need to be the same as binSizes and binNames", false, new StringParameterType());
    public GediParameter<Boolean> plot = new GediParameter<Boolean>(this,"plot", "Shall the R-plotting-script be run (needs RServ to be active)", false, new BooleanParameterType(), true);
    public GediParameter<Boolean> useLog = new GediParameter<Boolean>(this,"useLog", "If values should be log10", false, new BooleanParameterType(), true);
    public GediParameter<Boolean> noHeatmaps = new GediParameter<Boolean>(this,"noHeatmaps", "If plotting is enabled, heatmaps will be skipped (needs RServ to be active)", false, new BooleanParameterType(), true);
    public GediParameter<NormalizationMode> normalizationMode = new GediParameter<NormalizationMode>(this,"normMode", "Normalization mode used for normalization (Total=sum over all readcounts, Maximum=Max readcount)", false, new EnumParameterType<>(NormalizationMode.class), NormalizationMode.TOTAL, true);
    public GediParameter<BinNormalizationMode> perBinNormMode = new GediParameter<BinNormalizationMode>(this,"binNormMode", "Normalization mode used per bin", false, new EnumParameterType<>(BinNormalizationMode.class), BinNormalizationMode.MEAN, true);
    public GediParameter<Double> totalReadCounts = new GediParameter<Double>(this,"totalRC", "Total read counts. Either nothing or same length as all conditions", true, new DoubleParameterType(), true);
    public GediParameter<Double> plotYMax = new GediParameter<Double>(this,"plotYMax", "The max value for the y-axis", false, new DoubleParameterType(), -1.d,true);
    public GediParameter<ReadType> readType = new GediParameter<ReadType>(this,"readType", "Density, 3'- or 5'-end", false, new EnumParameterType<>(ReadType.class), ReadType.DENSITY,true);
    public GediParameter<ReadCountMode> readCountMode = new GediParameter<ReadCountMode>(this, "readCountMode", "mode of how to count reads", false, new GediParameterType<ReadCountMode>() {
        @Override
        public ReadCountMode parse(String s) {
            return ReadCountMode.valueOf(s);
        }

        @Override
        public Class<ReadCountMode> getType() {
            return ReadCountMode.class;
        }
    }, ReadCountMode.Weight, true);
    public GediParameter<String> namesFile = new GediParameter<String>(this,"names", "File containing the names for the genes in the locations file. Need to be in the same order!", false, new StringParameterType(), null, true);

    public GediParameter<File> voidout = new GediParameter<File>(this, "${prefix}.dataExtractionSuccessful.txt", "Only needed to make results directory overridable", false, new FileParameterType());
    public GediParameter<File> voidout2 = new GediParameter<File>(this, "${prefix}.void.txt", "Only needed to make results directory overridable", false, new FileParameterType());
}
