# MetagenePlot
A tool to plot metagene plots as well as metagene heatmaps for the gedi toolkit

# Quick guide
To run `MetagenePlot` download the latest release and execute the following console command (with appropriate parameters): `java -cp MetagenePlot.jar -prefix <Prefix> -binSizes <binSizes> -binNames <binNames> -readData <reads>.cit -locFile <locafile>.tsv <Options>`

An example call could look like this:
```
java -cp MetagenePlot.jar -prefix outputFolder -readData reads.cit -strandness SENSE -binSizes "100,10,500" -binNames "upstream,InterestingPosition,Downstream" -locFile locations.tsv -plot -readType DENSITY -normMode TOTAL -binNormMode MAX -readCountMode Weight -mergeConds "0,1/2/3"
```

The following options are available:
```
General:
 -prefix <prefix>                   The prefix used for all output files
 -binSizes <binSizes>               The sizes of each bin comma separated (e.g.: 40,60,20 this creates three bins with a total of 120 bins)
 -binNames <binNames>               Comma separated string containing names for the bins. Need to be the same size (i.e. same amounts of commas as binSizes)

MetageneDataExtractor:
 -readData <readData>               Read data in CIT-format.
 -locFile <locFile>                 File containing the locations to check for. The number of locations per line need to be the same as binSizes and binNames
[ -strandness <strandness>          Strandness of the data (default: Sense) Choices are : AutoDetect,Sense,Antisense,Unspecific]
[ -mergeConds <mergeConds>          Conditions to merge. Write as 0,2/1,3/4,5 to merge conditions 0 and 2, 1 and 3 and so on. (default: )]
[ -normMode <normMode>              Normalization mode used for normalization (Total=sum over all readcounts, Maximum=Max readcount) (default: TOTAL) Choices are : TOTAL,MAXIMUM,NONE]
[ -totalRC <totalRC>                Total read counts. Either nothing or same length as all conditions]
[ -readType <readType>              Density, 3'- or 5'-end (default: DENSITY) Choices are : FIVE_PRIME,DENSITY,THREE_PRIME]
[ -readCountMode <readCountMode>    mode of how to count reads (default: Weight)]
[ -binNormMode <binNormMode>        Normalization mode used per bin (default: MEAN) Choices are : MEAN,MAX,SUM]
[ -useLog                           If values should be log10]

MetagenePlotter:
[ -plot                             Shall the R-plotting-script be run (needs RServ to be active)]
[ -noHeatmaps                       If plotting is enabled, heatmaps will be skipped (needs RServ to be active)]
[ -plotYMax <plotYMax>              The max value for the y-axis (default: -1.0)]
[ -names <names>                    File containing the names for the genes in the locations file. Need to be in the same order!]
```

# Locations file (important!)
