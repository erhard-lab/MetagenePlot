package gedi.metagene;

import gedi.centeredDiskIntervalTree.CenteredDiskIntervalTreeStorage;
import gedi.core.data.reads.AlignedReadsData;
import gedi.core.data.reads.ReadCountMode;
import gedi.core.reference.Strandness;
import gedi.core.region.GenomicRegionStorage;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.utils.ReadType;
import org.junit.Assert;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.util.concurrent.atomic.AtomicReference;

class MetageneBinManagerTest {
    private MetageneBinManager binManager;

    @BeforeEach
    void setUp() {
        String cit = this.getClass().getResource("/resources/metageneBinner/mNETseqSmall.cit").getPath();
        GenomicRegionStorage<AlignedReadsData> storage = CenteredDiskIntervalTreeStorage.load(cit);
        MetageneLocation rgr1 = MetageneLocation.parse("1+:11694-12000");
        MetageneLocation rgr2 = MetageneLocation.parse("1+:12000-12308");
        MetageneLocation[] rgrs = new MetageneLocation[]{rgr1, rgr2};
        int[] binSized = new int[]{10, 20};

        MetageneBinner metageneBinner1 = new MetageneBinner(rgrs, binSized);
        MetageneBinner metageneBinner2 = new MetageneBinner(rgrs, binSized);
        metageneBinner1.extractReadCounts(storage, Strandness.Sense, null, ReadType.DENSITY, ReadCountMode.Weight);
        metageneBinner2.extractReadCounts(storage, Strandness.Sense, null, ReadType.DENSITY, ReadCountMode.Weight);
        metageneBinner1.normalizeReadCount(NormalizationMode.TOTAL);
        metageneBinner2.normalizeReadCount(NormalizationMode.TOTAL);
        metageneBinner1.mergeReadcountsIntoBins(BinNormalizationMode.MEAN);
        metageneBinner2.mergeReadcountsIntoBins(BinNormalizationMode.MEAN);

        binManager = new MetageneBinManager(storage.getMetaDataConditions(), new String[] {"bin1", "bin2"}, binSized);
        binManager.addMetageneBin(metageneBinner1);
        binManager.addMetageneBin(metageneBinner2);
    }

    @Test
    void iterateIndividualEntries() {
        AtomicReference<Double> bin0Sum = new AtomicReference<>(0.0);
        AtomicReference<Double> bin1Sum = new AtomicReference<>(0.0);
        AtomicReference<Double> bin2Sum = new AtomicReference<>(0.0);
        double[][][] cond9entry = new double[1][][];
        binManager.iterateIndividualEntries(9).forEachRemaining(e -> {
            bin0Sum.updateAndGet(v -> v + e[0][0]);
            bin1Sum.updateAndGet(v -> v + e[0][1]);
            bin2Sum.updateAndGet(v -> v + e[0][2]);
            cond9entry[0] = e;
        });

        Assert.assertEquals(0.00263*2, bin0Sum.get(), 0.00002);
        Assert.assertEquals(0.0*2, bin1Sum.get(), 0.00002);
        Assert.assertEquals(0.00454*2, bin2Sum.get(), 0.00002);
        Assert.assertEquals(2, cond9entry[0].length);
        Assert.assertEquals(10, cond9entry[0][0].length);
        Assert.assertEquals(20, cond9entry[0][1].length);
    }

    @Test
    void iterateIndividualMergedEntries() {
    }

    @Test
    void getSummedData() {
        double[] cond9sum = binManager.getSummedData(9, 0);
        Assert.assertEquals(0.00263*2, cond9sum[0], 0.00002);
        Assert.assertEquals(0.0*2, cond9sum[1], 0.00002);
        Assert.assertEquals(0.00454*2, cond9sum[2], 0.00002);
    }

    @Test
    void getSummedMergedData() {
    }
}