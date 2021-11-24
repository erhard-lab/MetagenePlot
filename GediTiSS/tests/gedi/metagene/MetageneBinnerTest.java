package gedi.metagene;

import gedi.centeredDiskIntervalTree.CenteredDiskIntervalTreeStorage;
import gedi.core.data.reads.AlignedReadsData;
import gedi.core.reference.Strandness;
import gedi.core.region.GenomicRegionStorage;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.util.ArrayUtils;
import gedi.utils.ArrayUtils2;
import org.junit.Assert;
import org.junit.jupiter.api.Test;

import java.lang.reflect.Field;

class MetageneBinnerTest {

    @Test
    void binning() throws NoSuchFieldException, IllegalAccessException {
        // SETUP
        MetageneLocation rgr1 = MetageneLocation.parse("1+:10-20");
        MetageneLocation rgr2 = MetageneLocation.parse("1+:30-35");
        MetageneLocation rgr3 = MetageneLocation.parse("1+:40-60");
        MetageneLocation rgr4 = MetageneLocation.parse("1+:70-75");
        MetageneLocation rgr5 = MetageneLocation.parse("1+:75-80");
        MetageneLocation rgr6 = MetageneLocation.parse("1+:80-81");
        MetageneLocation[] rgrs = new MetageneLocation[]{rgr1,rgr2,rgr3,rgr4,rgr5,rgr6};
        int[] binSized = new int[]{5, 3, 1, 5, 6, 1};
        MetageneBinner metageneBinner = new MetageneBinner(rgrs, binSized);
        // bin - cond - rc
        double[] readcounts1 = new double[] {10,9,8,7,6,5,4,3,2,1};
        double[] readcounts2 = new double[] {3,6,9,12,15};
        double[] readcounts3 = new double[] {2,3,4,6,8,1,2,3,4,52,2,12,13,13,12,4,9,0,0,2};
        double[] readcounts4 = new double[] {0,0,2,2,20};
        double[] readcounts5 = new double[] {0,0,2,2,20};
        double[] readcounts6 = new double[] {2.6};
        double[][][] readCounts = new double[6][1][];
        readCounts[0][0] = readcounts1;
        readCounts[1][0] = readcounts2;
        readCounts[2][0] = readcounts3;
        readCounts[3][0] = readcounts4;
        readCounts[4][0] = readcounts5;
        readCounts[5][0] = readcounts6;
        Field readCountField = metageneBinner.getClass().getDeclaredField("readCounts");
        readCountField.setAccessible(true);
        readCountField.set(metageneBinner, readCounts);

        // EXPECTATIONS
        double[] bin1 = new double[] {9.5, 7.5, 5.5, 3.5, 1.5};
        double[] bin2 = new double[] {4.5, 10.5, 15};
        double[] bin3 = new double[] {7.6};
        double[] bin4 = new double[] {0, 0, 2, 2, 20};
        double[] bin5 = new double[] {4.8, 4.8, 4.8, 4.8, 4.8, 4.8};
        double[] bin6 = new double[] {2.6};

        // EXECUTION
        metageneBinner.mergeReadcountsIntoBins(BinNormalizationMode.MEAN);
        double[] binnedReadCounts1 = metageneBinner.getBinnedReadCounts(0, 0);
        double[] binnedReadCounts2 = metageneBinner.getBinnedReadCounts(0, 1);
        double[] binnedReadCounts3 = metageneBinner.getBinnedReadCounts(0, 2);
        double[] binnedReadCounts4 = metageneBinner.getBinnedReadCounts(0, 3);
        double[] binnedReadCounts5 = metageneBinner.getBinnedReadCounts(0, 4);
        double[] binnedReadCounts6 = metageneBinner.getBinnedReadCounts(0, 5);

        // TESTING
        Assert.assertArrayEquals(bin1, binnedReadCounts1, 0.0001);
        Assert.assertArrayEquals(bin2, binnedReadCounts2, 0.0001);
        Assert.assertArrayEquals(bin3, binnedReadCounts3, 0.0001);
        Assert.assertArrayEquals(bin4, binnedReadCounts4, 0.0001);
        Assert.assertArrayEquals(bin5, binnedReadCounts5, 0.0001);
        Assert.assertArrayEquals(bin6, binnedReadCounts6, 0.0001);
    }

    @Test
    void normalizing() throws NoSuchFieldException, IllegalAccessException {
        // SETUP
        MetageneLocation rgr1 = MetageneLocation.parse("1+:10-20");
        MetageneLocation rgr2 = MetageneLocation.parse("1+:30-35");
        MetageneLocation rgr3 = MetageneLocation.parse("1+:40-60");
        MetageneLocation rgr4 = MetageneLocation.parse("1+:70-75");
        MetageneLocation[] rgrs = new MetageneLocation[]{rgr1,rgr2,rgr3,rgr4};
        int[] binSized = new int[]{5, 3, 1, 5};
        MetageneBinner metageneBinner = new MetageneBinner(rgrs, binSized);
        // bin - cond - rc
        double[] readcounts1 = new double[] {10,9,8,7,6,5,4,3,2,1};
        double[] readcounts2 = new double[] {3,6,9,12,15};
        double[] readcounts3 = new double[] {2,3,4,6,8,1,2,3,4,52,2,12,13,13,12,4,9,0,0,2};
        double[] readcounts4 = new double[] {0,0,2,2,20};
        double[][][] readCounts = new double[4][1][];
        readCounts[0][0] = readcounts1;
        readCounts[1][0] = readcounts2;
        readCounts[2][0] = readcounts3;
        readCounts[3][0] = readcounts4;
        Field readCountField = metageneBinner.getClass().getDeclaredField("readCounts");
        readCountField.setAccessible(true);
        readCountField.set(metageneBinner, readCounts);

        // EXPECTATIONS
        double totalReadCount = ArrayUtils.sum(readcounts1) + ArrayUtils.sum(readcounts2) + ArrayUtils.sum(readcounts3) + ArrayUtils.sum(readcounts4);
        double[] bin1 = ArrayUtils2.divide(readcounts1, totalReadCount);
        double[] bin2 = ArrayUtils2.divide(readcounts2, totalReadCount);
        double[] bin3 = ArrayUtils2.divide(readcounts3, totalReadCount);
        double[] bin4 = ArrayUtils2.divide(readcounts4, totalReadCount);

        // EXECUTION
        metageneBinner.normalizeReadCount(NormalizationMode.TOTAL);
        double[] normalizedReadCounts1 = metageneBinner.getReadCounts(0, 0);
        double[] normalizedReadCounts2 = metageneBinner.getReadCounts(0, 1);
        double[] normalizedReadCounts3 = metageneBinner.getReadCounts(0, 2);
        double[] normalizedReadCounts4 = metageneBinner.getReadCounts(0, 3);

        // TESTING
        Assert.assertArrayEquals(bin1, normalizedReadCounts1, 0.0001);
        Assert.assertArrayEquals(bin2, normalizedReadCounts2, 0.0001);
        Assert.assertArrayEquals(bin3, normalizedReadCounts3, 0.0001);
        Assert.assertArrayEquals(bin4, normalizedReadCounts4, 0.0001);
    }
}