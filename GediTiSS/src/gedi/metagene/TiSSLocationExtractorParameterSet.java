package gedi.metagene;

import gedi.core.genomic.Genomic;
import gedi.util.program.GediParameter;
import gedi.util.program.GediParameterSet;
import gedi.util.program.parametertypes.FileParameterType;
import gedi.util.program.parametertypes.GenomicParameterType;
import gedi.util.program.parametertypes.IntParameterType;
import gedi.util.program.parametertypes.StringParameterType;

import java.io.File;

public class TiSSLocationExtractorParameterSet extends GediParameterSet {
    public GediParameter<String> prefix = new GediParameter<String>(this,"prefix", "The prefix used for all output files", false, new StringParameterType());
    public GediParameter<File> paramFile = new GediParameter<File>(this,"${prefix}.param", "File containing the parameters used to call TiSS", false, new FileParameterType());
    public GediParameter<Genomic> genomic = new GediParameter<Genomic>(this,"genomic", "The indexed GEDI genome.", true, new GenomicParameterType());

    public GediParameter<Integer> tssUpstream = new GediParameter<Integer>(this,"tssUpstream", "Basepairs upstream of the TiSS", false, new IntParameterType());
    public GediParameter<Integer> tssDownstream = new GediParameter<Integer>(this,"tssDownstream", "Basepairs downstream of the TiSS", false, new IntParameterType());
    public GediParameter<Integer> ttsUpstream = new GediParameter<Integer>(this,"ttsUpstream", "Basepairs upstream of the TTS", false, new IntParameterType());
    public GediParameter<Integer> ttsDownstream = new GediParameter<Integer>(this,"ttsDownstream", "Basepairs downstream of the TTS", false, new IntParameterType());

    public GediParameter<File> locationsOut = new GediParameter<File>(this, "${prefix}.tsv", "The final locations for MetagenePlotter", false, new FileParameterType());
    public GediParameter<File> voidout = new GediParameter<File>(this, "${prefix}.void.txt", "Only needed to make results directory overridable", false, new FileParameterType());
}
