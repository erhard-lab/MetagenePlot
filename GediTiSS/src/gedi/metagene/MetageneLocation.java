package gedi.metagene;

import gedi.core.reference.Chromosome;
import gedi.core.reference.ReferenceSequence;
import gedi.core.region.GenomicRegion;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.core.region.ReferenceGenomicRegion;
import gedi.util.StringUtils;

public class MetageneLocation {
    private ReferenceGenomicRegion location;
    boolean reverse;

    private MetageneLocation() {

    }

    private MetageneLocation(ReferenceGenomicRegion location, boolean reverse) {
        this.location = location;
        this.reverse = reverse;
    }

    public static MetageneLocation parse(String locationString) {
        String[] split = StringUtils.split(locationString, ":");
        if (split.length != 2) {
            throw new IllegalArgumentException("Could not parse location. Does not separate ref and reg with double-dot: " + locationString);
        }
        String[] regionString = StringUtils.split(split[1], "|");
        if (regionString.length == 0) {
            throw new IllegalArgumentException("Could not parse location. No region present: " + locationString);
        }
        ReferenceGenomicRegion rgr = ImmutableReferenceGenomicRegion.parse(locationString);

        if (rgr.getRegion() == null) {
            throw new IllegalArgumentException("Could not parse the provided region: " + locationString);
        }

        if (rgr.getRegion().getTotalLength() != 0) {
            return new MetageneLocation(rgr, false);
        }

        StringBuilder newRegionString = new StringBuilder();
        for (int i = regionString.length-1; i >= 0; i--) {
            String[] reg = StringUtils.split(regionString[i], "-");
            if (reg.length != 2) {
                throw new IllegalArgumentException("Could not parse location. One or more region parts do not consist of one start and one end:" + locationString);
            }
            newRegionString.append(reg[1]).append("-").append(reg[0]);
            if (i != 0) {
                newRegionString.append("|");
            }
        }

        rgr = ImmutableReferenceGenomicRegion.parse(split[0] + ":" + newRegionString.toString());
        if (rgr.getRegion().getTotalLength() == 0) {
            throw new IllegalArgumentException("Could not parse location. Region part is neither forward nor backward: " + locationString);
        }

        return new MetageneLocation(rgr, true);
    }

    public ReferenceGenomicRegion getLocation() {
        return location;
    }

    public boolean isReverse() {
        return reverse;
    }

    @Override
    public String toString() {
        return location.toLocationString() + "_isReverse:" + (isReverse() ? "true" : "false");
    }

    public ReferenceSequence getReference() {
        return location.getReference();
    }

    public GenomicRegion getRegion() {
        return location.getRegion();
    }
}
