package gedi.metagene;

import gedi.core.data.reads.AlignedReadsData;
import gedi.core.data.reads.ReadCountMode;
import gedi.core.reference.Strandness;
import gedi.core.region.GenomicRegionStorage;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.util.ArrayUtils;
import gedi.util.functions.EI;
import gedi.utils.ReadType;
import gedi.utils.TiSSUtils;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class MetageneBinner {
    // The dimension is the number of bins. Equal to dimensions of binSizes.
    private MetageneLocation[] locations;
    private int[] binSizes;
    // first dimension = locationBin; second dimension = condition; third dimension = read count
    private double[][][] readCounts;
    private double[][][] binnedReadCounts;
    private Random random = new Random();

    public MetageneBinner(MetageneLocation[] locations, int[] binSizes) {
        if (locations.length != binSizes.length) {
            throw new IllegalArgumentException("The locations array and binSizes array need to be of the same length! location: " + locations.length + ", binSizes: " + binSizes.length);
        }
        this.locations = locations;
        this.binSizes = binSizes;
    }

    public void extractReadCounts(GenomicRegionStorage<AlignedReadsData> storage, Strandness strandness, double[] totals, ReadType readType, ReadCountMode readCountMode) {
        readCounts = new double[locations.length][][];
        int[] conditions = EI.seq(0,storage.getMetaDataConditions().length).toIntArray();
        for (int i = 0; i < locations.length; i++) {
            if (readType == ReadType.DENSITY) {
//                double[][] test = new double[conditions.length][locations[i].getRegion().getTotalLength()];
//                for (int t = 0; t < test.length; t++) {
//                    double missingChance = random.nextDouble();
//                    for (int tt = 0; tt < test[t].length; tt++) {
//                        test[t][tt] = random.nextDouble() * 100;
//                        if (random.nextDouble() < missingChance) {
//                            test[t][tt] = 0;
//                        }
//                    }
//                }
//                readCounts[i] = test;
                readCounts[i] = TiSSUtils.extractReadDensitiesForConditionsFromSingleFile(storage, conditions, locations[i].getReference(), strandness, locations[i].getRegion(), readCountMode);
            } else {
                readCounts[i] = TiSSUtils.extractReadCountsForConditionsFromSingleFile(storage, conditions, locations[i].getReference(), strandness, locations[i].getRegion(), readType, readCountMode);
            }
//            if (locations[i].getReference().isMinus()) {
            if (locations[i].isReverse()) {
                for (double[] counts : readCounts[i]) {
                    ArrayUtils.reverse(counts);
                }
            }
        }
        if (totals != null && totals.length > 0) {
            for (int bin = 0; bin < readCounts.length; bin++) {
                for (int cond = 0; cond < readCounts[bin].length; cond++) {
                    for (int i = 0; i < readCounts[bin][cond].length; i++) {
                        readCounts[bin][cond][i] /= totals[cond];
                    }
                }
            }
        }
    }

    public void mergeReadCounts(int[][] mergeConditions) {
        double[][][] mergedReadCounts = new double[readCounts.length][mergeConditions.length][];
        for (int bin = 0; bin < readCounts.length; bin++) {
            for (int cond = 0; cond < mergeConditions.length; cond++) {
                mergedReadCounts[bin][cond] = new double[readCounts[bin][cond].length];
                for (int merge : mergeConditions[cond]) {
                    for (int i = 0; i < readCounts[bin][merge].length; i++) {
                        mergedReadCounts[bin][cond][i] += readCounts[bin][merge][i];
                    }
                }
            }
        }
        readCounts = mergedReadCounts;
    }

    public void normalizeReadCount(NormalizationMode normalizationMode) {
        if (readCounts == null || readCounts.length != locations.length) {
            throw new IllegalStateException("readcounts are null or of zero length. Did you miss calling extractReadCounts?");
        }
        for (int condition = 0; condition < readCounts[0].length; condition++) {
            double normalizer = 0;
            for (int location = 0; location < locations.length; location++) {
                for (int i = 0; i < readCounts[location][condition].length; i++) {
                    switch (normalizationMode) {
                        case TOTAL:
                            normalizer += readCounts[location][condition][i];
                            break;
                        case MAXIMUM:
                            normalizer = Math.max(readCounts[location][condition][i], normalizer);
                            break;
                        case NONE:
                            normalizer = 1;
                            break;
                    }
                }
            }
            if (normalizer == 0) {
                continue;
            }
            for (int location = 0; location < locations.length; location++) {
                for (int i = 0; i < readCounts[location][condition].length; i++) {
                    readCounts[location][condition][i] /= normalizer;
                }
            }
        }
    }

    public void mergeReadcountsIntoBins(BinNormalizationMode perBinNormMode) {
        if (readCounts == null || readCounts.length != locations.length) {
            throw new IllegalStateException("readcounts are null or of zero length. Did you miss calling extractReadCounts?");
        }
        binnedReadCounts = new double[locations.length][readCounts[0].length][];
        for (int bin = 0; bin < locations.length; bin++) {
            double basepairsPerBin = locations[bin].getRegion().getTotalLength() / (double) binSizes[bin];
            // In case our region is smaller than our bin-size, we set it all to the mean value
            // TODO: I should write a test for this and check whether this is a good way of handling these cases!
            if (basepairsPerBin < 1) {
                for (int condition = 0; condition < readCounts[bin].length; condition++) {
                    binnedReadCounts[bin][condition] = new double[binSizes[bin]];
                    double avg = ArrayUtils.sum(readCounts[bin][condition])/readCounts[bin][condition].length;
                    for (int i = 0; i < binnedReadCounts[bin][condition].length; i++) {
                        binnedReadCounts[bin][condition][i] = avg;
                    }
                }
                continue;
            }
            for (int condition = 0; condition < readCounts[bin].length; condition++) {
                double nextBinStop = basepairsPerBin;
                int binIndex = 0;
                binnedReadCounts[bin][condition] = new double[binSizes[bin]];
                double readCountSum = 0;
                int readsInBin = 0;
                double maxValue = -1;
                for (int i = 0; i < readCounts[bin][condition].length; i++) {
                    if (i >= Math.ceil(nextBinStop)) {
                        binnedReadCounts[bin][condition][binIndex] = getBinValue(readCountSum, readsInBin, maxValue, perBinNormMode);
                        readCountSum = 0.0;
                        binIndex++;
                        nextBinStop = basepairsPerBin * (binIndex+1);
                        readsInBin = 0;
                        maxValue = -1;
                    }
                    readCountSum += readCounts[bin][condition][i];
                    maxValue = readCounts[bin][condition][i] > maxValue ? readCounts[bin][condition][i] : maxValue;
                    readsInBin++;
                }
                if (readCountSum > 0) {
                    binnedReadCounts[bin][condition][binIndex] = getBinValue(readCountSum, readsInBin, maxValue, perBinNormMode);
                }
            }
        }
    }

    public double[] getBinnedReadCounts(int condition, int location) {
        return binnedReadCounts[location][condition];
    }

    public double[] getReadCounts(int condition, int location) {
        return readCounts[location][condition];
    }

    private double getBinValue(double readCountSum, int readsInBin, double maxValue, BinNormalizationMode normMode) {
        switch (normMode) {
            case MEAN:
                return readCountSum / readsInBin;
            case MAX:
                return maxValue;
            case SUM:
                return readCountSum;
        }
        throw new IllegalArgumentException();
    }
}
