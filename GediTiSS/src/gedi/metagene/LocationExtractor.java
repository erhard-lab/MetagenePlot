package gedi.metagene;

import gedi.core.genomic.Genomic;
import gedi.core.region.ArrayGenomicRegion;
import gedi.util.io.text.LineOrientedFile;
import gedi.util.io.text.LineWriter;
import gedi.util.program.GediProgram;
import gedi.util.program.GediProgramContext;

import java.io.IOException;

public class LocationExtractor extends GediProgram {
    public LocationExtractor(TiSSLocationExtractorParameterSet params) {
        addInput(params.prefix);
        addInput(params.genomic);
        addInput(params.tssUpstream);
        addInput(params.tssDownstream);
        addInput(params.ttsUpstream);
        addInput(params.ttsDownstream);

        addOutput(params.locationsOut);
        addOutput(params.voidout);
    }

    @Override
    public String execute(GediProgramContext context) throws Exception {
        String prefix = getParameter(0);
        Genomic genomic = getParameter(1);
        int tssUpstream = getParameter(2);
        int tssDownstream = getParameter(3);
        int ttsUpstream = getParameter(4);
        int ttsDownstream = getParameter(5);

        String locationsOut = getOutputFile(0).getPath();

        extractLocations(genomic, tssUpstream, tssDownstream, ttsUpstream, ttsDownstream, locationsOut);
        return null;
    }

    public static void extractLocations(Genomic genomic, int tssUpstream, int tssDownstream, int ttsUpstream, int ttsDownstream, String locationsOut) throws IOException {
        LineWriter writer = new LineOrientedFile(locationsOut).write();
        genomic.getTranscripts().ei().forEachRemaining(t -> {
            String upstreamtss = t.getReference().toPlusMinusString() + ":";
            String tss = t.getReference().toPlusMinusString() + ":";
            String downstreamtss = t.getReference().toPlusMinusString() + ":";
            String middle = t.getReference().toPlusMinusString() + ":";
            String upstreamtts = t.getReference().toPlusMinusString() + ":";
            String tts = t.getReference().toPlusMinusString() + ":";
            String downstreamtts = t.getReference().toPlusMinusString() + ":";
            if (t.getRegion().getTotalLength() < tssDownstream+ttsUpstream+4) {
                return;
            }
            if (t.getReference().isPlus()) {
                upstreamtss += new ArrayGenomicRegion(t.getRegion().getStart()-tssUpstream, t.getRegion().getStart()).toString2();
                tss += new ArrayGenomicRegion(t.getRegion().getStart(), t.getRegion().getStart()+1).toString2();
                downstreamtss += t.getRegion().map(new ArrayGenomicRegion(1, tssDownstream+1)).toString2();
                middle += t.getRegion().map(new ArrayGenomicRegion(tssDownstream+1, t.getRegion().getTotalLength()-(ttsUpstream+1))).toString2();
                upstreamtts += t.getRegion().map(new ArrayGenomicRegion(t.getRegion().getTotalLength()-(ttsUpstream+1), t.getRegion().getTotalLength()-1)).toString2();
                tts += new ArrayGenomicRegion(t.getRegion().getEnd()-1, t.getRegion().getEnd());
                downstreamtts += new ArrayGenomicRegion(t.getRegion().getEnd(), t.getRegion().getEnd()+ttsDownstream);
            } else {
                downstreamtts += new ArrayGenomicRegion(t.getRegion().getStart()-ttsDownstream, t.getRegion().getStart()).toString2();
                tts += new ArrayGenomicRegion(t.getRegion().getStart(), t.getRegion().getStart()+1).toString2();
                upstreamtts += t.getRegion().map(new ArrayGenomicRegion(1, ttsUpstream+1)).toString2();
                middle += t.getRegion().map(new ArrayGenomicRegion(ttsUpstream+1, t.getRegion().getTotalLength()-(tssDownstream+1))).toString2();
                downstreamtss += t.getRegion().map(new ArrayGenomicRegion(t.getRegion().getTotalLength()-(tssDownstream+1), t.getRegion().getTotalLength()-1)).toString2();
                tss += new ArrayGenomicRegion(t.getRegion().getEnd()-1, t.getRegion().getEnd());
                upstreamtss += new ArrayGenomicRegion(t.getRegion().getEnd(), t.getRegion().getEnd()+tssUpstream);
            }

            writer.writeLine2(upstreamtss + "\t" + tss + "\t" + downstreamtss + "\t" + middle + "\t" + upstreamtts + "\t" + tts + "\t" + downstreamtts);
        });
        writer.close();
    }
}
