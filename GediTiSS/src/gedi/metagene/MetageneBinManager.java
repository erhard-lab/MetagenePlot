package gedi.metagene;

import gedi.util.ArrayUtils;
import gedi.util.functions.EI;
import gedi.utils.ArrayUtils2;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public class MetageneBinManager {
    private String[] conditionNames;
    private String[] binNames;
    private int[] binSizes;
    private List<double[][][]> individialData;
    // condition - bin - values
    private double[][][] summedUpData;

    public MetageneBinManager(String[] conditionNames, String[] binNames, int[] binSizes) {
        this.conditionNames = conditionNames;
        this.binNames = binNames;
        this.binSizes = binSizes;
        individialData = new ArrayList<>();
        summedUpData = createSingleData();
    }

    // Returns an iterator containing a 2D array holding the bins in the first dimension and the specific values
    // in the second.
    public Iterator<double[][]> iterateIndividualEntries(int condition) {
        return EI.wrap(individialData).map(m -> m[condition]);
    }

    public Iterator<double[][]> iterateIndividualMergedEntries(int[] conditions) {
        return EI.wrap(individialData).map(m -> {
            double[][] merged = new double[binSizes.length][];
            for (int c = 0; c < conditions.length; c++) {
                for (int b = 0; b < binSizes.length; b++) {
                    if (merged[b] == null) {
                        merged[b] = m[c][b].clone();
                    } else {
                        merged[b] = ArrayUtils.add(merged[b], m[c][b]);
                    }
                }
            }
            return merged;
        });
    }

    public double[] getSummedData(int condition, int bin) {
        return summedUpData[condition][bin];
    }

    public double[] getMeanData(int condition, int bin) {
        return ArrayUtils2.divide(summedUpData[condition][bin], individialData.size());
    }

    public void addMetageneBin(MetageneBinner bin) {
        double[][][] data = createSingleData();
        for (int i = 0; i < conditionNames.length; i++) {
            for (int j = 0; j < binSizes.length; j++) {
                data[i][j] = bin.getBinnedReadCounts(i, j);
                summedUpData[i][j] = ArrayUtils.add(summedUpData[i][j], data[i][j]);
            }
        }
        individialData.add(data);
    }

    private double[][][] createSingleData() {
        double[][][] data = new double[conditionNames.length][binSizes.length][];
        for (int i = 0; i < data.length; i++) {
            for (int j = 0; j < data[i].length; j++) {
                data[i][j] = new double[binSizes[j]];
            }
        }
        return data;
    }
}
